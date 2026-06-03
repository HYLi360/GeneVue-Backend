#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>
"""
A series of GFF3-handle tools.

"""

import re
from collections import defaultdict
from pathlib import Path

import mmap
from typing import Optional, List, Dict

import pandas as pd
from rich.table import Table

from genevue import GFF3FileNotFoundError
from genevue.GXF import feature_line

from genevue import setup_rich_logger, console

logger = setup_rich_logger(__name__, console)
_SOURCE = "GVFIX"


class GenomicFeature:
    def __init__(
        self,
        seqid: str,
        source: str,
        type: str,
        start: str,
        end: str,
        score: str,
        strand: str,
        phase: str,
        attrs: str,
    ) -> None:
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = int(start)
        self.end = int(end)
        self.score = float(score) if score != "." else 0
        self.strand = strand
        self.phase = phase
        self.attrs = attrs.strip()

    def __str__(self):
        return f"""{self.id}
seqid:   {self.seqid}
source:  {self.source}
type:    {self.type}
start:   {self.start}
end:     {self.end}
score:   {self.score}
strand:  {self.strand}
phase:   {self.phase}
attrs:   {self.attr_dict}
"""

    def __lt__(self, other):
        if self.type != other.type:
            return self.type < other.type
        else:
            return self.start < other.start

    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)

    def __eq__(self, other):
        return (self.type == other.type) and (self.start == other.start)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __gt__(self, other):
        if self.type != other.type:
            return self.type > other.type
        else:
            return self.start > other.start

    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)

    @property
    def attr_dict(self) -> dict[str, str]:
        d = {}
        for attr in self.attrs.split(";"):
            ls = attr.split("=")
            d[ls[0]] = ls[1]
        return d

    @property
    def parent(self) -> str:
        return self.attr_dict.get("Parent", "")

    @property
    def id(self) -> str:
        return self.attr_dict.get("ID", "")


class GenomicFeatureTree:
    def __init__(self, gff3_path: Path):
        self.gff3_path = gff3_path
        self.features: Dict[str, GenomicFeature] = {}
        self.parent2featureid = defaultdict(List[str])
        self.root_features_ls: List[str] = []
        self._parse()

    def _parse(self):
        with open(self.gff3_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                linels = line.strip().split(maxsplit=8)
                feature = GenomicFeature(*linels)
                self.add_features(feature)

    def add_features(self, feature: GenomicFeature):
        fid = feature.id
        if self.features.get(feature.id) is not None:
            counter = 1
            while True:
                fid = f"{feature.id}-{counter}"
                if self.features.get(fid, None) is None:
                    self.features[fid] = feature
                    break
                counter += 1
        else:
            self.features[fid] = feature
        if feature.parent is not None:
            self.parent2featureid[feature.parent].append(fid)
        else:
            self.root_features_ls.append(fid)

    def children_features(self, feature_name: str) -> List:
        if feature_name not in self.parent2featureid.keys():
            return []
        else:
            return self.parent2featureid[feature_name]

    def iter_subfeatures(self, feature_name: str) -> List:
        if feature_name not in self.features:
            return []
        else:
            mid = [feature_name]
            res = []
            while mid:
                for i in mid:
                    if self.children_features(i):
                        mid.extend(self.children_features(i))
                    else:
                        res.append(i)
                    mid.remove(i)
            return res


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

                elif feature_type == b"mRNA" and id_m and parent_m:
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
                    data=list(
                        set(
                            self.tx_gene[i]
                            for i in self.tx_cds.keys()
                            if i in self.tx_gene
                        )
                    ),
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

    def _gff3_indices(self):
        """
        a more completed parser for gff3.
        """
