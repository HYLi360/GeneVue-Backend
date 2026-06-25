#  Copyright (c) 2026 HYLi360. All rights reserved.
#
#  see LICENSE in /LICENSE
#  see side-package LICENSEs (if used) in /LICENSE_OF_SIDE_PACKAGES

import gzip
import re
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Self, Tuple
from urllib.parse import unquote as url_unquote

import polars as pl
import rich
from pandas import DataFrame as pandasDF
from polars import DataFrame as PolarsDF
from rich.table import Table

from genevue import FormatNotSuitableError, console, setup_rich_logger
from genevue.GXF.GFF3tools import GFF3OPS
from genevue.Utils.FileSystem import check_filetype


class GFF3COL(Enum):
    SEQID = "seqid"
    SOURCE = "source"
    TYPE = "type"
    START = "start"
    END = "end"
    LENGTH = "length"
    SCORE = "score"
    STRAND = "strand_readable"
    PHASE = "phase"

    @property
    def expr(self) -> pl.Expr:
        return pl.col(self.value)


class GFF3:
    @property
    def seqid(self):
        return GFF3COL.SEQID

    @property
    def source(self):
        return GFF3COL.SOURCE

    @property
    def type(self):
        return GFF3COL.TYPE

    @property
    def start(self):
        return GFF3COL.START

    @property
    def end(self):
        return GFF3COL.END

    @property
    def length(self):
        return GFF3COL.LENGTH

    @property
    def score(self):
        return GFF3COL.SCORE

    @property
    def strand(self):
        return GFF3COL.STRAND

    def __init__(self):
        self.db: PolarsDF = PolarsDF()

    def load_from_file(self, gff3_path: str | Path) -> None:
        """
        Read a GFF3 file, try to flatten its attribute column, and return it
        as a Polars DataFrame.

        If the file is in Parquet format, only read it with no additional
        operations.

        Keyword Args:
            gff3_path:
                The path of GFF3, which is a plain text (gff3), gzip, bgzip,
                or Apache Parquet file.

        Raises:
            FormatNotAcceptableError:
                GFF3 is not in the format described above.

        Warnings:
            EmptyDFWarning:
                Returns an empty PolarsDF.

        See Also:
            genevue.Formats.GFF3
        """
        if check_filetype(Path(gff3_path), "gz"):
            gff3 = gzip.open(gff3_path, "rt")
        elif check_filetype(Path(gff3_path), "gff3"):
            gff3 = open(gff3_path)
        elif check_filetype(Path(gff3_path), "parquet"):
            self.db = pl.read_parquet(Path(gff3_path))
            return
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

        db = PolarsDF(data=data)

        if db.filter(pl.col("type") == "mRNA").shape[0] == 0:
            # CDS -> gene, no txs
            cds_gene_lookup = db.filter(pl.col("type") == "CDS").select(
                pl.col("ID").alias("CDS"), pl.col("Parent").alias("gene")
            )

            # combine
            db = db.join(cds_gene_lookup, left_on="ID", right_on="CDS", how="left")

        else:
            # with txs
            # Search backward to find the corresponding mRNA/CDS -> gene
            mrna_gene_lookup = db.filter(pl.col("type") == "mRNA").select(
                pl.col("ID").alias("_mrna_id"), pl.col("Parent").alias("_gene_of_mrna")
            )
            db = db.join(
                mrna_gene_lookup, left_on="Parent", right_on="_mrna_id", how="left"
            )
            db = db.with_columns(
                pl.when(pl.col("type") == "gene")
                .then(pl.col("ID"))
                .when(pl.col("type") == "mRNA")
                .then(pl.col("Parent"))
                .otherwise(pl.col("_gene_of_mrna"))
                .alias("gene")
            ).drop("_gene_of_mrna")

        self.db = db

    def brief_report(self) -> None:
        """
        Print a simple report of GFF3 Database.
        """
        seqids_df = self.db.group_by("seqid").agg(pl.count().alias("count"))
        console.print(f"Seqids count:                {seqids_df.shape[0]}")

        features_count = self.db.shape[0]
        small_count = seqids_df.filter(pl.col("count") < (features_count // 20)).shape[
            0
        ]
        console.print(
            "Small seqs count\n"
            f"(features less than {features_count // 20}): {small_count}"
        )
        console.print(
            "Genes count:                 "
            f"{self.db.filter(pl.col('type') == 'gene').shape[0]}"
        )
        console.print(
            "mRNAs count:                 "
            f"{self.db.filter(pl.col('type') == 'mRNA').shape[0]}"
        )
        console.print(
            "CDSs count:                  "
            f"{self.db.filter(pl.col('type') == 'CDS').shape[0]}"
        )
        console.print(f"Total features count:         {features_count}")

        console.print("Source and type of features:")
        source_and_type_df = (
            self.db.group_by(["source", "type"])
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
                _source, _type, str(_counts), f"{_counts / features_count * 100: 3.2f}%"
            )
        console.print(source_and_type_table)

    def search_by_exact_val(self, search_query) -> Self:
        pass

    def look_back_until(self):
        pass

    def get_longest(
        self, target: Literal["mRNA", "CDS", "protein"] = "mRNA"
    ) -> PolarsDF:
        """
        For each gene, find the transcript (mRNA), total CDS or protein
        with the longest CDS length.

        Keyword Args:
            target:
                Features used for filtering the longest CDS.

        Returns:
            DataFrame with columns: gene, best_id
        """
        _expr = [
            pl.col("ID")
            .sort_by(pl.col("length"), descending=True)
            .first()
            .alias("best_id")
        ]

        # if protein_id appears in gff3 db?
        # that is common in gff3 from NCBI
        if "protein_id" in self.db.columns:
            _expr.append(
                pl.col("protein_id")
                .sort_by("length", descending=True)
                .first()
                .alias("best_protein"),
            )

        cds = self.db.select(self.type.expr == "cds")
        mrna = self.db.select(self.type.expr == "mRNA")

        if cds.is_empty() or (target == "CDS"):
            # No mRNA feature, or targets on CDS
            # Use "retrospective_until" to get gene feature in any cases
            # TODO: ???
            cds_and_length = cds.group_by(pl.col("ID")).agg(pl.col("length").sum())
            cds = cds.join(cds_and_length, left_on="ID", right_on="ID", how="left")
            return cds.group_by(pl.col("gene")).agg(_expr)
        else:
            # Have mRNA feature(s)
            mrna_and_length = cds.group_by(pl.col("Parent")).agg(pl.col("length").sum())
            mrna = mrna.join(
                mrna_and_length, left_on="ID", right_on="Parent", how="left"
            )
            # get protein id
            if "protein_id" in self.db.columns:
                mrna = mrna.drop("protein_id").join(
                    cds.select([pl.col("Parent"), pl.col("protein_id")]),
                    left_on="ID",
                    right_on="Parent",
                    how="left",
                )
            return mrna.group_by(pl.col("gene")).agg(_expr)

    def get_chrom_length(self) -> Dict[str, int]:
        res = self.db.group_by(self.seqid.expr).agg(self.end.expr.max())
        return dict(zip(res[self.seqid.value].to_list(), res[self.end.value].to_list()))

    def get_gene_density(self, bin_size: int = 500_000) -> Dict[str, Dict[int, int]]:
        """
        Aggregate statistics on the distribution of the number of genes
        across all chromosomes within an equal interval.

        Keyword Args:
            bin_size: interval size (in 'bp').

        Returns:
            A dictionary in format of

            >>> {
            >>>    'chrom1': {500_000: 1, 1_000_000: 2, 1_377_831: 3},
            >>>    'chrom2': {347_113: 4},
            >>> }

            which

            - chrom1, chrom2 stand for name of chromosomes.

            - 5_000_000: 1 stands for 'there is a gene in interval of [0, 500_000]'.
              for the first, the interval always starts from 0 bp; for the last, the
              key always stands for the length of chromosome; and for other interval,
              the interval always starts from the key ahead it.
        """
        max_chr_len = max(self.get_chrom_length().values())

        res_df = (
            self.db.select(self.type.expr == "gene")
            .group_by(self.seqid.expr)
            .agg(
                self.end.expr.cut(range(0, max_chr_len, bin_size))
                .value_counts()
                .alias("countls")
            )
        )

        res = dict(zip(res_df[self.seqid.value], res_df["countls"]))

        pass

    def to_parquet(self, out_path: str | Path):
        self.db.to_pandas().to_parquet(out_path)


g = GFF3()
print(g.source.value)
