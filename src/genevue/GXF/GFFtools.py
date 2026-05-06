#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>
"""
A series of GFF-handle tools.

"""

import re
from collections import defaultdict
from pathlib import Path

import mmap
import pandas as pd
from rich.table import Table

from genevue import GFF3FileNotFoundError
from genevue import console
from genevue.GXF import feature_line


class SimpleGFF3:
    def __init__(self, gff3_file_path: str | Path):
        if isinstance(gff3_file_path, str):
            gff3_file_path = Path(gff3_file_path)

        if not gff3_file_path.resolve().exists():
            raise GFF3FileNotFoundError(gff3_file_path.name)

        self.gff3_file_path = gff3_file_path

        self.bed = pd.DataFrame()
        self.genelist = []
        self.tx_gene = {}
        self.gene_txs = defaultdict(list)
        self.tx_cds = defaultdict(list)

        self._gff3_indices_simple()

    def _gff3_indices_simple(self):
        # line_re: Only parse which column 3 has gene/mRNA/transcript/CDS feature.
        # id_re: Extract the ID=(...) in attribule column.
        # parent_re: Extract the Parent=(..) in attribule column.
        id_re = re.compile(rb"(?:^|;)ID=([^;]+)")
        parent_re = re.compile(rb"Parent=([^;]+)")

        # Prepare result collector.
        _df = []

        # Load gff3 file as memory-mapped file.
        with open(self.gff3_file_path, "rb") as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

            for m in feature_line.finditer(mm):
                # Unpack line
                chro, feature_type, start, end, strand, phase, attrs = m.groups()

                # Check id and parent (may return None in gene feature line)
                id_m = id_re.search(attrs)
                parent_m = parent_re.search(attrs)

                if feature_type == b"gene" and id_m:
                    # 1    ensembl    gene    1    201    .    +    .    ID=gene;...
                    gene_id = id_m.group(1).decode("ascii", "ignore")
                    # chr, gene_id, start, end, strand
                    _df.append(
                        {
                            "chr": str(chro.decode("ascii", "ignore")),
                            "gene_id": gene_id,
                            "start": int(start),
                            "end": int(end),
                            "strand": strand.decode("ascii"),
                        }
                    )

                elif feature_type in (b"mRNA", b"transcript") and id_m and parent_m:
                    # 1 ensembl mRNA 1 201 . + . ID=transcript;Parent=gene;...
                    tx_id = id_m.group(1).decode("ascii", "ignore")
                    gene_id = parent_m.group(1).decode("ascii", "ignore")
                    self.tx_gene[tx_id] = gene_id
                    self.gene_txs[gene_id].append(tx_id)

                elif feature_type == b"CDS" and id_m and parent_m:
                    # 1 ensembl CDS 1 201 . + 0 ID=cds;Parent=transcript;...
                    cds_id = id_m.group(1).decode("ascii", "ignore")
                    tx_id = parent_m.group(1).decode("ascii", "ignore")
                    self.tx_cds[tx_id].append(
                        [
                            cds_id,
                            str(chro.decode("ascii", "ignore")),
                            int(start),
                            int(end),
                            strand.decode("ascii", "ignore"),
                            phase.decode("ascii", "ignore"),
                        ]
                    )

        self.bed = (
            pd.merge(
                left=pd.Series(
                    data=list(set(self.tx_gene[i] for i in self.tx_cds.keys())),
                    name="gene_id",
                ),
                right=pd.DataFrame(data=_df),
                on="gene_id",
                how="inner",
            )
            .sort_values(["chr", "start"])
            .reset_index(drop=True)
        )
        self.bed["order"] = self.bed.groupby("chr").cumcount() + 1
        self.genelist = self.bed["gene_id"].to_list()

        table = Table(
            title="Results Summary",
        )
        table.add_column("Entries Type", style="cyan", no_wrap=True)
        table.add_column("Quantity", style="green", justify="right")
        table.add_row("Protein Coding Gene", f"{len(self.genelist)}")
        table.add_row("Transcriptable Gene", f"{len(self.gene_txs.keys())}")

        table.add_row("Protein Coding Transcript", f"{len(self.tx_cds.keys())}")
        table.add_row("All Transcript", f"{len(self.tx_gene.keys())}")
        table.add_row("CDS", f"{sum(len(sub) for sub in self.tx_cds.values())}")
        console.print(table)


class GFF3:
    def __init__(self, gff3_file_path: str | Path):
        if isinstance(gff3_file_path, str):
            gff3_file_path = Path(gff3_file_path)

        if not gff3_file_path.resolve().exists():
            raise GFF3FileNotFoundError(gff3_file_path.name)

        self.gff3_file_path = gff3_file_path

        self.bed = pd.DataFrame()
        self.genelist = []
        self.tx_gene = {}
        self.gene_txs = defaultdict(list)
        self.tx_cds = defaultdict(list)

        self._gff3_indices_simple()

    def _gff3_indices_simple(self):
        # line_re: Only parse which column 3 has gene/mRNA/transcript/CDS feature.
        # id_re: Extract the ID=(...) in attribule column.
        # parent_re: Extract the Parent=(..) in attribule column.
        id_re = re.compile(rb"(?:^|;)ID=([^;]+)")
        parent_re = re.compile(rb"Parent=([^;]+)")

        # Prepare result collector.
        _df = []

        # Load gff3 file as memory-mapped file.
        with open(self.gff3_file_path, "rb") as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

            for m in feature_line.finditer(mm):
                # Unpack line
                chro, feature_type, start, end, strand, phase, attrs = m.groups()

                # Check id and parent (may return None in gene feature line)
                id_m = id_re.search(attrs)
                parent_m = parent_re.search(attrs)

                if feature_type == b"gene" and id_m:
                    # 1    ensembl    gene    1    201    .    +    .    ID=gene;...
                    gene_id = id_m.group(1).decode("ascii", "ignore")
                    # chr, gene_id, start, end, strand
                    _df.append(
                        {
                            "chr": str(chro.decode("ascii", "ignore")),
                            "gene_id": gene_id,
                            "start": int(start),
                            "end": int(end),
                            "strand": strand.decode("ascii"),
                        }
                    )

                elif feature_type in (b"mRNA", b"transcript") and id_m and parent_m:
                    # 1 ensembl mRNA 1 201 . + . ID=transcript;Parent=gene;...
                    tx_id = id_m.group(1).decode("ascii", "ignore")
                    gene_id = parent_m.group(1).decode("ascii", "ignore")
                    self.tx_gene[tx_id] = gene_id
                    self.gene_txs[gene_id].append(tx_id)

                elif feature_type == b"CDS" and id_m and parent_m:
                    # 1 ensembl CDS 1 201 . + 0 ID=cds;Parent=transcript;...
                    cds_id = id_m.group(1).decode("ascii", "ignore")
                    tx_id = parent_m.group(1).decode("ascii", "ignore")
                    self.tx_cds[tx_id].append(
                        [
                            cds_id,
                            str(chro.decode("ascii", "ignore")),
                            int(start),
                            int(end),
                            strand.decode("ascii", "ignore"),
                            phase.decode("ascii", "ignore"),
                        ]
                    )

        self.bed = (
            pd.merge(
                left=pd.Series(
                    data=list(set(self.tx_gene[i] for i in self.tx_cds.keys())),
                    name="gene_id",
                ),
                right=pd.DataFrame(data=_df),
                on="gene_id",
                how="inner",
            )
            .sort_values(["chr", "start"])
            .reset_index(drop=True)
        )
        self.bed["order"] = self.bed.groupby("chr").cumcount() + 1
        self.genelist = self.bed["gene_id"].to_list()

        table = Table(
            title="Results Summary",
        )
        table.add_column("Entries Type", style="cyan", no_wrap=True)
        table.add_column("Quantity", style="green", justify="right")
        table.add_row("Protein Coding Gene", f"{len(self.genelist)}")
        table.add_row("Transcriptable Gene", f"{len(self.gene_txs.keys())}")

        table.add_row("Protein Coding Transcript", f"{len(self.tx_cds.keys())}")
        table.add_row("All Transcript", f"{len(self.tx_gene.keys())}")
        table.add_row("CDS", f"{sum(len(sub) for sub in self.tx_cds.values())}")
        console.print(table)


def build_bed():
    """
    A reformer to transform the GFF3 file to BED.
    """
