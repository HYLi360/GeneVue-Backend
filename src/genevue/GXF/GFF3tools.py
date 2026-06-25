"""
src/genevue/plot/CollScatPlot.py

(C) 2026 HYLi360. All rights reserved.

see LICENSE in LICENSE
see side-package LICENSEs (if used) in LICENSE_OF_SIDE_PACKAGES

--------------------
A series of GFF3-handle tools. Mainly powered by Polars, a high-perference DataFrame library.
"""

import gzip
import re
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple
from urllib.parse import unquote as url_unquote

import polars as pl
import rich
from rich.table import Table

from genevue import FormatNotSuitableError, console, setup_rich_logger
from genevue.Utils.FileSystem import check_filetype

logger = setup_rich_logger(__name__, console)
_SOURCE = "GVFIX"


class BlazingGFF3:
    def __init__(self, bgff3_path: Optional[Path] = None):
        self.data = pl.DataFrame()
        self.gff3_path = Path()
        self.bgff3_path = Path()

        if bgff3_path is not None:
            self.read(bgff3_path)

    def build(self, gff3_path: Path) -> None:
        self.gff3_path = gff3_path

        if check_filetype(gff3_path, "gz"):
            gff3 = gzip.open(gff3_path, "rt")
        elif check_filetype(gff3_path, "gff3"):
            gff3 = open(gff3_path)
        else:
            raise FormatNotSuitableError

        data = []

        # main build
        for line in gff3:
            if line.startswith("#"):
                continue

            linels = line.strip().split(maxsplit=8)
            feature = dict(
                zip(
                    [
                        "seqid",
                        "source",
                        "type",
                        "start",
                        "end",
                        "length",
                        "score",
                        "strand",
                        "strand_readable",
                        "phase",
                    ],
                    [
                        linels[0],
                        linels[1],
                        linels[2] if linels[2] != "transcript" else "mRNA",
                        int(linels[3]),
                        int(linels[4]),
                        int(linels[4]) - int(linels[3]) + 1,
                        float(linels[5]) if linels[5] != "." else 0,
                        {"+": 1, "-": -1}.get(linels[6], 0),
                        linels[6],
                        int(linels[7]) if linels[7] != "." else 0,
                    ],
                )
            )

            _attr_string = linels[8]
            if _attr_string[-1] == ";":
                _attr_string = _attr_string[:-1]

            for attr in _attr_string.split(";"):
                ls = attr.split("=", maxsplit=1)
                feature[ls[0]] = url_unquote(ls[1])

            data.append(feature)

        gff3.close()

        self.data = pl.DataFrame(data=data)

        if self.data.filter(pl.col("type") == "mRNA").shape[0] == 0:
            # CDS -> gene
            cds_gene_lookup = self.data.filter(pl.col("type") == "CDS").select(
                pl.col("ID").alias("CDS"), pl.col("Parent").alias("gene")
            )

            # combine
            self.data = self.data.join(
                cds_gene_lookup, left_on="ID", right_on="CDS", how="left"
            )

        else:
            # Search backward to find the corresponding mRNA/CDS -> gene
            mrna_gene_lookup = self.data.filter(pl.col("type") == "mRNA").select(
                pl.col("ID").alias("_mrna_id"), pl.col("Parent").alias("_gene_of_mrna")
            )
            self.data = self.data.join(
                mrna_gene_lookup, left_on="Parent", right_on="_mrna_id", how="left"
            )
            self.data = self.data.with_columns(
                pl.when(pl.col("type") == "gene")
                .then(pl.col("ID"))
                .when(pl.col("type") == "mRNA")
                .then(pl.col("Parent"))
                .otherwise(pl.col("_gene_of_mrna"))
                .alias("gene")
            ).drop("_gene_of_mrna")

    def read(self, bgff3_path: Path) -> None:
        self.bgff3_path = bgff3_path
        self.data = pl.read_parquet(bgff3_path)

    def brief_report(self):
        seqids_df = self.data.group_by("seqid").agg(pl.count().alias("count"))
        rich.print(f"Seqids count:                {seqids_df.shape[0]}")

        features_count = self.data.shape[0]
        small_count = seqids_df.filter(pl.col("count") < (features_count // 20)).shape[
            0
        ]
        rich.print(
            f"Small seqs count\n(features less than {features_count // 20}): {small_count}"
        )
        rich.print(
            f"Genes count:                 {self.data.filter(pl.col('type') == 'gene').shape[0]}"
        )
        rich.print(
            f"mRNAs count:                 {self.data.filter(pl.col('type') == 'mRNA').shape[0]}"
        )
        rich.print(
            f"CDSs count:                  {self.data.filter(pl.col('type') == 'CDS').shape[0]}"
        )
        rich.print(f"Total eatures count:         {features_count}")

        rich.print("Source and type of features:")
        source_and_type_df = (
            self.data.group_by(["source", "type"])
            .agg(pl.count().alias("counts"))
            .sort("counts", descending=True)
        )
        source_and_type_ls = list(source_and_type_df.iter_rows())
        source_and_type_table = Table()
        source_and_type_table.add_column("Source", no_wrap=True)
        source_and_type_table.add_column("Type", no_wrap=True)
        source_and_type_table.add_column("Counts", no_wrap=True)
        source_and_type_table.add_column("Percent", justify="right", no_wrap=True)
        for _source, _type, _counts in source_and_type_ls:
            source_and_type_table.add_row(
                _source, _type, str(_counts), f"{_counts/features_count*100: 3.2f}%"
            )
        rich.print(source_and_type_table)

    # === Search ===
    def search_text(
        self, text: List[str] | str, columns: Optional[List[str]] = None
    ) -> pl.DataFrame:
        """
        Search text(s) which appears in bgff3 certainly column(s), or all columns.
        Notes that it always do `or` operation if you input at least 2 conditions (contains columns, texts, or both).
        """
        # select all columns if you not specified
        if columns is None:
            columns = [
                col for col in self.data.columns if self.data[col].dtype == pl.Utf8
            ]

        # force `text` turn into `list` if only get a string
        textls = [text] if isinstance(text, str) else text

        # build search conditions
        conditions = []
        for col in columns:
            for text in textls:
                conditions.append(pl.col(col).str.contains(re.escape(text)))

        query = conditions[0]
        for cond in conditions[1:]:
            query = query | cond

        # search
        return self.data.filter(query)

    def search_exact(
        self, value: Any, columns: Optional[List[str]] = None
    ) -> pl.DataFrame:
        # select all columns if you not specified
        if columns is None:
            # only scan the columns which have the same data type of value
            columns = [
                col for col in self.data.columns if self.data[col].dtype == type(value)
            ]
        else:
            columns = [
                col
                for col in self.data.columns
                if (self.data[col].dtype == type(value)) and (col in columns)
            ]

        conditions = []
        for col in columns:
            conditions.append(pl.col(col).eq(value))

        query = conditions[0]

        for cond in conditions[1:]:
            query = query | cond

        return self.data.filter(query)

    def search_regex(self, regex: str):
        """
        Regex enhanced searching.
        """
        pass

    # === Range Search ===
    def search_region(self, seqid: str, start: int = 0, end: int = -1):
        if end == -1:
            return self.data.filter(pl.col("seqid") == seqid)

        return self.data.filter(
            (pl.col("seqid") == seqid)
            & (pl.col("start") >= start)
            & (pl.col("end") <= end)
        )

    # === get anything you want ===
    def get_feature(self, feature_id: str) -> List[Dict]:
        features = self.data.filter(pl.col("ID") == feature_id)

        if features.shape[0] == 0:
            return []

        res = []

        for feature in features.iter_rows(named=True):
            res.append({k: w for k, w in feature.items() if w is not None})

        return res

    def get_subfeatures(
        self,
        feature_id: str | List[str],
        recursive: bool = False,
        containing_self: bool = False,
    ) -> pl.DataFrame:
        """
        Get children features of this feature.

        Return a empty pl.DataFrame if
          - feature name doesn't appear in gff3, or
          - this feature has no children (a leaf feature).

        Args:
            feature_id
            recursive
            containing_self
        """

        if not recursive:
            return self.data.filter(pl.col("Parent").is_in(feature_id))

        res: List[pl.DataFrame] = []
        if isinstance(feature_id, str):
            current_layer = [feature_id]
        else:
            current_layer = feature_id
        visited = set()

        if containing_self:
            res.append(self.data.filter(pl.col("ID") == feature_id))

        while current_layer:
            children = self.data.filter(pl.col("Parent").is_in(current_layer))

            if len(children) == 0:
                break

            res.append(children)

            next_layer = children.select(pl.col("ID")).to_series().unique().to_list()
            current_layer = [_id for _id in next_layer if _id not in visited]
            visited.update(current_layer)

        if res:
            return pl.concat(res)

        return pl.DataFrame()

    def get_parent_feature(
        self,
        feature_id: str | List[str],
        max_depth: int = 100,
        until: str = "",
    ) -> pl.DataFrame:
        if isinstance(feature_id, str):
            feature_id = [feature_id]

        result_ids = set(feature_id)
        current_ids = set(feature_id)

        for depth in range(max_depth):
            parents = (
                self.data.filter(pl.col("ID").is_in(list(current_ids)))
                .select("Parent")
                .drop_nulls()["Parent"]
                .to_list()
            )

            if not parents:
                break

            current_ids = set(parents)
            result_ids.update(current_ids)

            if until:
                types = (
                    self.data.filter(pl.col("ID").is_in(list(current_ids)))
                    .select("type")["type"]
                    .to_list()
                )
                if until in types:
                    break

        return self.data.filter(
            pl.col("ID")
            .is_in(list(result_ids))
            .and_((pl.col("type") == until).or_(pl.col("Parent").is_null()))
        )

    def get_gene_by_protein_id(self, protein_id: str):
        if "protein_id" in self.data.columns:
            items = (
                self.data.filter(pl.col("protein_id") == protein_id)["ID"]
                .unique()
                .to_list()
            )
            return self.get_parent_feature(items, until="gene")
        return pl.DataFrame()

    def get_parent_id(self, feature_id: str) -> List[str]:
        df = self.get_parent_feature(feature_id)

        if df.shape[0] == 0:
            return []
        else:
            return df.select(pl.col("ID")).unique().to_series().to_list()

    def get_full_hierarchy_dict(
        self, feature_id: str, containing_self: bool = True
    ) -> List[Dict]:
        subfeatures = self.get_subfeatures(
            feature_id, recursive=True, containing_self=containing_self
        )

        if subfeatures.shape[0] == 0:
            return []

        res = []

        for feature in subfeatures.iter_rows(named=True):
            res.append({k: w for k, w in feature.items() if w is not None})

        return res

    def get_longest_transcript_per_gene(
        self, mode: Literal["mRNA", "CDS"] = "mRNA"
    ) -> pl.DataFrame:
        """
        For each gene, find the transcript (mRNA) with the longest total CDS length.

        GFF3 hierarchy: gene -> mRNA -> CDS.
        Total CDS length per mRNA = sum of all its CDS feature lengths.
        If the 'protein_id' column exists, prefer it over the mRNA ID.

        Returns:
            DataFrame with columns: gene, best_id, [best_protein]
        """
        mrna = self.data.filter(pl.col("type") == "mRNA")
        cds = self.data.filter(pl.col("type") == "CDS")

        _expr = [
            pl.col("ID")
            .sort_by(pl.col("length"), descending=True)
            .first()
            .alias("best_id")
        ]

        if "protein_id" in self.data.columns:
            _expr.append(
                pl.col("protein_id")
                .sort_by("length", descending=True)
                .first()
                .alias("best_protein"),
            )

        if mrna.is_empty() or (mode == "CDS"):
            cds_and_length = cds.group_by(pl.col("ID")).agg(pl.col("length").sum())
            cds = cds.join(cds_and_length, left_on="ID", right_on="ID", how="left")
            return cds.group_by(pl.col("gene")).agg(_expr)
        else:
            mrna_and_length = cds.group_by(pl.col("Parent")).agg(pl.col("length").sum())
            mrna = mrna.join(
                mrna_and_length, left_on="ID", right_on="Parent", how="left"
            )
            # get protein id
            if "protein_id" in self.data.columns:
                mrna = mrna.drop("protein_id").join(
                    cds.select([pl.col("Parent"), pl.col("protein_id")]),
                    left_on="ID",
                    right_on="Parent",
                    how="left",
                )
            return mrna.group_by(pl.col("gene")).agg(_expr)

    # === Aggregative analyse ===
    def get_gene_density(self, bin_size: float | int = 100_000) -> pl.DataFrame:
        """
        Get the gene density from bgff3 (genes per bin_size), useful for chromosome plotting.
        Result is a 3-column DataFrame:
        pl.DataFrame(
            {
                "seqid": [..., ..., ...,],
                "start": [0, 100_000, 200_000,],
                "end":   [100_000, 200_000, 300_000,]
                "genes": [1_111, 2_222, 3_333,],
                "density": [0.011_11, 0.022_22, 0.033_33,],
            }
        )
        """
        genes = self.data.filter(pl.col("type") == "gene")

        return (
            genes.with_columns((pl.col("start") // bin_size).alias("steps"))
            .group_by(["seqid", "steps"])
            .agg(
                [
                    pl.col("start").min().alias("start"),
                    pl.col("end").max().alias("end"),
                    pl.count().alias("genes"),
                ]
            )
            .sort(by=["seqid", "steps"])
        )

    # === Modifing ===
    def add_feature(self):
        pass

    def modify_feature(self):
        pass

    # === Export ===
    def to_parquet(self, parquet_path: Path):
        """recommend."""
        self.data.to_pandas().to_parquet(parquet_path)

    def to_bed(
        self,
        bed_path: Optional[Path] = None,
        mode: Literal["gene", "mRNA", "CDS", "protein"] = "gene",
        chrname_filter: Optional[List[str]] = None,
    ) -> pl.DataFrame:
        if mode == "gene":
            if not chrname_filter:
                chrname_filter = self.data["seqid"].unique().to_list()
            df = self.data.filter(
                (pl.col("type") == "gene") & pl.col("seqid").is_in(chrname_filter)
            ).select(
                pl.col("seqid").alias("chrom"),
                pl.col("start").alias("chromStart"),
                pl.col("end").alias("chromEnd"),
                pl.col("ID").alias("name"),
                pl.col("score"),
                pl.col("strand_readable").alias("strand"),
            )
        elif mode == "mRNA":
            df = self.get_longest_transcript_per_gene(mode="mRNA")
            df = df.join(
                self.to_bed(mode="gene"),
                left_on="gene",
                right_on="name",
                how="right",
            ).select(
                pl.col("chrom"),
                pl.col("chromStart"),
                pl.col("chromEnd"),
                pl.col("best_id").alias("name"),
                pl.col("score"),
                pl.col("strand"),
            )
        elif mode == "CDS":
            df = self.get_longest_transcript_per_gene(mode="CDS")
            df = df.join(
                self.to_bed(mode="gene"),
                left_on="gene",
                right_on="name",
                how="right",
            ).select(
                pl.col("chrom"),
                pl.col("chromStart"),
                pl.col("chromEnd"),
                pl.col("best_id").alias("name"),
                pl.col("score"),
                pl.col("strand"),
            )
        else:
            df = self.get_longest_transcript_per_gene(mode="mRNA")
            df = df.join(
                self.to_bed(mode="gene"),
                left_on="gene",
                right_on="name",
                how="right",
            ).select(
                pl.col("chrom"),
                pl.col("chromStart"),
                pl.col("chromEnd"),
                pl.col("best_protein").alias("name"),
                pl.col("score"),
                pl.col("strand"),
            )

        if bed_path is not None:
            df.to_pandas().to_csv(bed_path, sep="\t", index=False, header=False)

        return df

    @property
    def sources(self) -> List[str]:
        if self.data.shape[0] == 0:
            return []
        else:
            return self.data.select("source").to_series().unique().to_list()

    @property
    def types(self) -> List[str]:
        if self.data.shape[0] == 0:
            return []
        else:
            return self.data.select("type").to_series().unique().to_list()

    @property
    def sources_and_types(self) -> List[Tuple[str, str]]:
        if self.data.shape[0] == 0:
            return []
        else:
            return list(self.data.select(["source", "type"]).unique().iter_rows())
