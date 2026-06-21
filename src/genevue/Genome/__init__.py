import mmap
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord
from rich.progress import track
from rich.table import Table

from genevue import SeqFileNotFoundError, GFF3FileNotFoundError, FormatNotSuitableError
from genevue import console
from genevue.Utils import header_detector


class Genome:
    def __init__(self, fna_file_path: str | Path, gff3_file_path: str | Path) -> None:
        console.log("Packup Genome DNA and SimpleGFF3 file.")
        init_start = time.time()
        with console.status("Initialing..."):
            # Parsing fasta file.
            fna_file_path = Path(fna_file_path).resolve()
            if not fna_file_path.exists():
                raise SeqFileNotFoundError(fna_file_path)

            # If the sequences file is in gzip format?
            # That's not suitable for random access.
            # we need raise a Error to tell user using the `bgzip-recompressor`.
            if header_detector(fna_file_path) == "gzip":
                raise FormatNotSuitableError(
                    f"{fna_file_path} is in gzip format, which is not suitable for random access. Use "
                    f"bgzip format file, or use `bgzip-recompressor` tool to recompress and parse again."
                )

            self.genome_seqs = SeqIO.index(fna_file_path, "fasta")
            console.log(
                f"Genome sequence indexing finished in {time.time() - init_start:3f} seconds."
            )

            # Read and format original GFF3.
            parse_gff3_start = time.time()
            gff3_file_path = Path(gff3_file_path).resolve()
            if not gff3_file_path.exists():
                raise GFF3FileNotFoundError(gff3_file_path)
            self.simp_gff, self.genelist, self.gene_txs, self.tx_gene, self.tx_cds = (
                _build_gff3_indices(gff3_file_path)
            )
        console.log(
            f"SimpleGFF3 processing finished in {time.time() - parse_gff3_start:3f} seconds."
        )
        console.log(f"Initial finished in {time.time() - init_start:3f} seconds.")

    # ------------------------------------------------------------
    # Core methods
    # ------------------------------------------------------------
    def get_cds(self, tx_id, description: str = "") -> SeqRecord:
        """
        Return spliced CDS sequence for transcript tx_id.
        Handles strand and phase.
        """
        if tx_id not in self.tx_cds:
            raise KeyError(f"Transcript {tx_id} has no CDS features")

        #           0    1      2    3       4      5
        # list[[tx_id, chr, start, end, strand, phase], ......]
        segs = self.tx_cds[tx_id]

        # pre-test
        chrom, strand = segs[0][1], segs[0][4]
        if segs[0][1] not in self.genome_seqs.keys():
            raise KeyError(f"Chromosome {chrom} not found in fasta")

        record_seq = self.genome_seqs[chrom]

        if record_seq is None:
            raise KeyError(f"No chromosome sequence for {chrom}")

        # extract start
        pieces = []
        for _, _, start, end, _, _ in segs:
            # GFF is 1-based inclusive; Python is 0-based half-open
            frag = record_seq[start - 1 : end].seq
            pieces.append(frag)

        cds_seq = Seq("").join(pieces)

        # handle negative strand
        if strand == "-":
            cds_seq = cds_seq.reverse_complement()

        cds = SeqRecord(seq=cds_seq, id=tx_id, name=tx_id, description=description)

        return cds

    def get_pep(self, tx_id, to_stop=True, strict=False) -> SeqRecord:
        """
        Translate CDS to peptide using Biopython.
        strict=True uses cds=True for strict checks.
        """
        cds_seq = self.get_cds(tx_id)
        return self.translate(cds_seq, to_stop, strict)

    @staticmethod
    def translate(
        cds_seq: SeqRecord,
        to_stop=True,
        strict=False,
        description: Optional[str] = None,
    ) -> SeqRecord:
        # If not divisible by 3, you can choose to trim or not.
        # Here is a gentle approach: trim tail to multiple of 3 when not strict.
        if not strict:
            # Defence for not standard CDS length.
            trim = len(cds_seq) % 3
            if trim:
                cds_seq = cds_seq[:-trim]

        try:
            pep = cds_seq.translate(to_stop=to_stop, cds=strict)
            pep.id = cds_seq.id
            pep.name = cds_seq.name
            pep.description = (
                description if description is not None else cds_seq.description
            )
        except Exception as e:
            raise ValueError(f"Translation failed for {cds_seq.id}: {e}")
        return pep

    # ------------------------------------------------------------
    # Longest CDS/Prot per gene
    # ------------------------------------------------------------
    def get_longest(self, gene_id, to_stop=True, strict=False):
        """
        Return (best_tx_id, best_cds_seq, best_pep_seq)
        best determined by CDS nucleotide length (after phase trimming/splicing)
        """
        if gene_id not in self.gene_txs:
            raise KeyError(f"Gene {gene_id} has no transcripts")

        best = None
        for tx_id in self.gene_txs[gene_id]:
            if tx_id not in self.tx_cds:
                continue
            cds = self.get_cds(tx_id)
            if best is None or len(cds) > len(best[1]):
                best = (tx_id, cds)

        if best is None:
            raise ValueError(f"Gene {gene_id} has no CDS")

        tx_id, cds = best
        pep = self.get_pep(tx_id, to_stop=to_stop, strict=strict)
        return tx_id, cds, pep

    def iter_longest_pep(self, to_stop=True, strict=False):
        for gene_id in self.gene_txs:
            try:
                tx_id, cds, pep = self.get_longest(
                    gene_id, to_stop=to_stop, strict=strict
                )
                yield gene_id, tx_id, pep
            except Exception:
                continue

    def get_cds_pep_by_gene(
        self,
        out_cds_path: Optional[str] = None,
        out_pep_path: Optional[str] = None,
        genelist: Optional[Iterable[str]] = None,
        to_stop: bool = True,
        strict: bool = False,
        warp: int = 70,
    ):
        console.log("Extract longest CDS and Protein sequences.")
        extract_start = time.time()
        cds_ls, pep_ls = [], []
        genelist = self.genelist if genelist is None else genelist
        for gene_id in track(genelist, "Extracting..."):
            tx_id, cds, pep = self.get_longest(gene_id, to_stop, strict)
            cds.id, cds.name, cds.description = gene_id, gene_id, f"tx_id={tx_id}"
            pep.id, pep.name, pep.description = gene_id, gene_id, f"tx_id={tx_id}"
            cds_ls.append(cds)
            pep_ls.append(pep)
        if (out_cds_path is not None) and (out_pep_path is not None):
            cds_handle = FastaWriter(out_cds_path, wrap=warp)
            pep_handle = FastaWriter(out_pep_path, wrap=warp)
            cds_handle.write_file(cds_ls)
            pep_handle.write_file(pep_ls)
        else:
            console.log("NO specified the cds & pep output file path. skip writting.")

        console.log(
            f"Extract finished in {time.time() - extract_start:3f} seconds. "
            f"{len(cds_ls)} CDS sequences, "
            f"{len(pep_ls)} Protein sequences exported."
        )

        return cds_ls, pep_ls

    def get_cds_pep_by_tx(
        self,
        out_cds_path: Optional[str] = None,
        out_pep_path: Optional[str] = None,
        genelist: Optional[Iterable[str]] = None,
        to_stop: bool = True,
        strict: bool = False,
        warp: int = 70,
    ):
        console.log("Extract CDS and Protein sequences.")
        extract_start = time.time()
        cds_ls, pep_ls = [], []
        genelist = self.genelist if genelist is None else genelist
        for gene_id in track(genelist, "Extracting..."):
            for tx_id in self.gene_txs[gene_id]:
                cds = self.get_cds(tx_id, description=f"gene_id={gene_id}")
                pep = self.translate(cds, to_stop, strict)
                cds_ls.append(cds)
                pep_ls.append(pep)
        if (out_cds_path is not None) and (out_pep_path is not None):
            cds_handle = FastaWriter(out_cds_path, wrap=warp)
            pep_handle = FastaWriter(out_pep_path, wrap=warp)
            cds_handle.write_file(cds_ls)
            pep_handle.write_file(pep_ls)
        else:
            console.log("No specified the cds & pep output file path. skip writting.")

        console.log(
            f"Extract finished in {time.time() - extract_start:3f} seconds. "
            f"{len(cds_ls)} CDS sequences, "
            f"{len(pep_ls)} Protein sequences exported."
        )

        return cds_ls, pep_ls


class GenomeLight:
    """
    A lighten class to process cds, pep and gff3 file.

    It's similar to class `Genome`, but you don't need to use translate tool, such as `gffread`, or `translate` from
    `Bio`, so it has much lesser pressure of your hard disk (for I/O and genome sequence storaging usage.).

    In addition, thanks to the addition of header checking, it natively supports the .gzip file format, regardless of
    whether the content is in Fasta or GFF3 format.

    There have 2 situation for your cds & pep data:

    1. These sequences was named by genes' name, not transcripts' name. In this time, we will only make a simple
    `pd.DataFrame` (`simp_gff`), which like this:

    ```
                gene_id	tx_id	chr	start	end	strand	phase
        0	gene1	rna1	1	371878	371957	+	0
    ```

    2. These sequences was named by transcripts' name. since as one gene can associates with many transcripts, we still
    have to make a index for genes and transcripts, like in `Genome`, and note the gene name as descriptions for every
    records.

    You can change that mode easily by parameter `mode`.
    """

    def __init__(
        self,
        cds_file_path: str | Path,
        pep_file_path: str | Path,
        gff3_file_path: str | Path,
    ) -> None:
        console.log("Packup CDS, PEP and GFF3 file.")
        init_start = time.time()
        with console.status("Initialing..."):
            cds_file_path = Path(cds_file_path).resolve()
            if not cds_file_path.exists():
                raise FileNotFoundError(cds_file_path)
            self.cds = SeqIO.index(cds_file_path, "fasta")

            pep_file_path = Path(pep_file_path).resolve()
            if not pep_file_path:
                raise FileNotFoundError(pep_file_path)
            self.pep = SeqIO.index(pep_file_path, "fasta")

            console.log(
                f"Genome sequence indexing finished in {time.time() - init_start:3f} seconds."
            )

            # Read and format original SimpleGFF3.
            parse_gff3_start = time.time()
            gff3_file_path = Path(gff3_file_path).resolve()
            if not gff3_file_path.exists():
                raise FileNotFoundError(gff3_file_path)

            self.simp_gff, self.genelist, self.gene_txs, self.tx_gene, self.tx_cds = (
                _build_gff3_indices(gff3_file_path)
            )

            self.chr = self.simp_gff["chr"].drop_duplicates().to_list()

            console.log(
                f"GFF3 processing finished in {time.time() - parse_gff3_start:3f} seconds."
            )
            console.log(f"Initial finished in {time.time() - init_start:3f} seconds.")


def _build_gff3_indices(gff3_file_path: Path):
    # line_re: Only parse which column 3 has gene/mRNA/transcript/CDS feature.
    # id_re: Extract the ID=(...) in attribule column.
    # parent_re: Extract the Parent=(..) in attribule column.
    line_re = re.compile(
        rb"^(?P<chr>[^\t#][^\t]*)\t[^\t]*\t(gene|mRNA|transcript|CDS)\t(?P<start>\d+)\t(?P<end>\d+)\t[^\t]*\t("
        rb"?P<strand>[+-])\t[^\t]*\t(?P<attrs>[^\n]*)$",
        re.M,
    )
    id_re = re.compile(rb"(?:^|;)ID=([^;]+)")
    parent_re = re.compile(rb"Parent=([^;]+)")

    # Prepare result collector.
    _df = []
    tx_gene = {}
    gene_txs = defaultdict(list)
    tx_cds = defaultdict(list)

    # Load gff3 file as memory-mapped file.
    with open(gff3_file_path, "rb") as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

        for m in line_re.finditer(mm):
            line = m.group(0)

            # Unpack line
            chro, _, feature_type, start, end, _, strand, phase, attrs = line.split(
                b"\t"
            )

            # Check id and parent (may return None in gene feature line)
            id_m = id_re.search(attrs)
            parent_m = parent_re.search(attrs)

            if feature_type == b"gene" and id_m:
                # 1	ensembl	gene	1	201	.	+	.	ID=gene;...
                gene_id = id_m.group(1).decode("ascii", "ignore")
                # chr, gene_id, start, end, strand
                _df.append(
                    {
                        "chr": chro.decode("ascii", "ignore"),
                        "gene_id": gene_id,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand.decode("ascii"),
                    }
                )

            elif feature_type in (b"mRNA", b"transcript") and id_m and parent_m:
                # 1	ensembl	mRNA	1	201	.	+	.	ID=transcript;Parent=gene;...
                tx_id = id_m.group(1).decode("ascii", "ignore")
                gene_id = parent_m.group(1).decode("ascii", "ignore")
                tx_gene[tx_id] = gene_id
                gene_txs[gene_id].append(tx_id)

            elif feature_type == b"CDS" and parent_m:
                # 1	ensembl	CDS	1	201	.	+	0	ID=cds;Parent=transcript;...
                tx_id = parent_m.group(1).decode("ascii", "ignore")
                tx_cds[tx_id].append(
                    [
                        tx_id,
                        chro.decode("ascii", "ignore"),
                        int(start),
                        int(end),
                        strand.decode("ascii", "ignore"),
                        phase.decode("ascii", "ignore"),
                    ]
                )

    simp_gff = (
        pd.merge(
            left=pd.Series(
                data=list(set(tx_gene[i] for i in tx_cds.keys() if i in tx_gene)),
                name="gene_id",
            ),
            right=pd.DataFrame(data=_df),
            on="gene_id",
            how="inner",
        )
        .sort_values(["chr", "start"])
        .reset_index(drop=True)
    )
    simp_gff["order"] = simp_gff.groupby("chr").cumcount() + 1
    simp_gff = simp_gff.sort_values(["chr", "start"]).reset_index(drop=True)
    genelist = simp_gff["gene_id"].to_list()

    table = Table(
        title="Results Summary",
    )
    table.add_column("Entries Type", style="cyan", no_wrap=True)
    table.add_column("Quantity", style="green", justify="right")
    table.add_row("Protein Coding Gene", f"{len(genelist)}")
    table.add_row("Transcriptable Gene", f"{len(gene_txs.keys())}")
    table.add_row("Protein Coding Transcript", f"{len(tx_cds.keys())}")
    table.add_row("All Transcript", f"{len(tx_gene.keys())}")
    table.add_row("CDS", f"{sum(len(sub) for sub in tx_cds.values())}")
    console.print(table)

    return simp_gff, genelist, dict(gene_txs), dict(tx_gene), dict(tx_cds)
