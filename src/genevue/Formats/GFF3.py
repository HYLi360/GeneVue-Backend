#  Copyright (c) 2026 HYLi360. All rights reserved.
#
#  see LICENSE in /LICENSE
#  see side-package LICENSEs (if used) in /LICENSE_OF_SIDE_PACKAGES

import gzip
from enum import Enum
from pathlib import Path
from tempfile import TemporaryFile
from typing import Dict, List, Literal, Optional, Self, Tuple
from urllib.parse import unquote as url_unquote

import polars as pl
from polars import DataFrame as PolarsDF
from rich.table import Table

from genevue import (
    BaseWarning4GeneVue,
    FormatNotSuitableError,
    console,
    setup_rich_logger,
)
from genevue.Utils.FileSystem import check_filetype

logger = setup_rich_logger(__name__, console)


class GFF3DBLinesTooLargeWarning(BaseWarning4GeneVue):
    def __init__(self, lines):
        self.message = (
            f"This GFF3 database has {lines} lines! that may run out your RAM. "
            "Consider using 'GFF3().sink_to_parquet()' before read to avoid."
        )

    def __str__(self) -> str:
        return f"{self.message}"


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

    def load_from_file(self, gff3_path: str | Path, batch_size: int = 500_000) -> None:
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

        data: List[Dict] = []
        db: PolarsDF = PolarsDF()

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

            if len(data) > batch_size:
                db = pl.concat([db, PolarsDF(data=data)], how="diagonal")
                data = []

        db = pl.concat([db, PolarsDF(data=data)], how="diagonal")
        data = []
        gff3.close()

        if db.filter(pl.col("type") == "mRNA").shape[0] == 0:
            # CDS -> gene, no txs
            # Build ID -> gene mapping as a lightweight dict to avoid
            # materializing an intermediate DataFrame copy from join.
            cds_to_gene = dict(
                db.filter(pl.col("type") == "CDS")
                .select("ID", pl.col("Parent").alias("gene"))
                .iter_rows()
            )
            db = db.with_columns(
                pl.col("ID").replace_strict(cds_to_gene, default=None).alias("gene")
            )

        else:
            # with txs
            # Build mRNA.ID → mRNA.Parent(=gene) mapping as a dict
            # instead of a join lookup to keep peak memory at ~1× db.
            mrna_to_gene = dict(
                db.filter(pl.col("type") == "mRNA").select("ID", "Parent").iter_rows()
            )
            db = db.with_columns(
                pl.when(pl.col("type") == "gene")
                .then(pl.col("ID"))
                .when(pl.col("type") == "mRNA")
                .then(pl.col("Parent"))
                .otherwise(pl.col("Parent").replace_strict(mrna_to_gene, default=None))
                .alias("gene")
            )

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
        # protein target rely on certain column
        if (target == "protein") and ("protein_id" not in self.db.columns):
            return PolarsDF()

        if "protein_id" in self.db.columns:
            _expr = (
                pl.col("protein_id")
                .sort_by("length", descending=True)
                .first()
                .alias("best_id")
            )
        else:
            _expr = (
                pl.col("ID")
                .sort_by(pl.col("length"), descending=True)
                .first()
                .alias("best_id")
            )

        cds = self.db.filter(self.type.expr == "CDS")
        mrna = self.db.filter(self.type.expr == "mRNA")

        if cds.is_empty() or (target == "CDS"):
            # No mRNA feature, or targets on CDS
            cds_and_length = dict(
                cds.group_by(pl.col("ID")).agg(pl.col("length").sum()).iter_rows()
            )
            cds = cds.with_columns(
                pl.col("ID").replace_strict(cds_and_length, 0).alias("length")
            )
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

    def get_subfeatures(
        self,
        feature_id: str | List[str],
        recursive: bool = False,
        containing_self: bool = False,
    ) -> PolarsDF:
        """
        Get children features of this feature.

        Return an empty pl.DataFrame if
          - feature name doesn't appear in gff3, or
          - this feature has no children (a leaf feature).

        Args:
            feature_id
            recursive
            containing_self
        """
        if isinstance(feature_id, str):
            feature_id: List[str] = [feature_id]
        else:
            feature_id = list(dict.fromkeys(feature_id))

        if not recursive:
            return self.db.filter(pl.col("Parent").is_in(feature_id))

        res: List[PolarsDF] = []

        if containing_self:
            res.append(self.db.filter(pl.col("ID").is_in(feature_id)))

        visited = set()

        while feature_id:
            children = self.db.filter(pl.col("Parent").is_in(feature_id))

            if len(children) == 0:
                break

            res.append(children)

            feature_id = [_id for _id in children["ID"].unique() if _id not in visited]
            visited.update(feature_id)

        if res:
            return pl.concat(res)

        return PolarsDF()

    def get_parent_feature(
        self,
        feature_id: str | List[str],
        level: int = -1,
        until: str = "",
        containing_self: bool = False,
    ) -> PolarsDF:
        if isinstance(feature_id, str):
            feature_id: List[str] = [feature_id]
        else:
            feature_id = list(dict.fromkeys(feature_id))

        res: List[PolarsDF] = []

        if containing_self:
            res.append(self.db.filter(pl.col("ID").is_in(feature_id)))

        visited = set()

        if level < 0:
            level = 0

        _expr = pl.col("ID").is_in(feature_id)

        if until:
            _expr = _expr and pl.col("Type") == until

        while level + 1:
            parents = self.db.filter(
                pl.col("ID").is_in(
                    self.db.filter(_expr)["Parent"].drop_nulls().unique()
                )
            )

            if len(parents) == 0:
                break

            res.append(parents)

            feature_id = [_id for _id in parents["ID"].unique() if _id not in visited]
            visited.update(feature_id)

            level -= 1

        if res:
            return pl.concat(res)

        return PolarsDF()

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

    def to_bed(
        self,
        bed_path: str | Path = ".",
        mode: Literal["gene", "mRNA", "CDS", "protein"] = "gene",
        chrname_filter: Optional[List[str]] = None,
    ) -> PolarsDF:
        if chrname_filter is None:
            chrname_filter = self.db["seqid"].unique().to_list()

        df = self.db.filter(
            (pl.col("type") == "gene") & pl.col("seqid").is_in(chrname_filter)
        ).select(
            pl.col("seqid").alias("chrom"),
            pl.col("start").alias("chromStart"),
            pl.col("end").alias("chromEnd"),
            pl.col("ID").alias("name"),
            pl.col("score"),
            pl.col("strand_readable").alias("strand"),
        )

        match mode:
            case "gene":
                pass
            case "mRNA":
                mrna_ids = dict(self.get_longest(target="mRNA").iter_rows())
                df = df.with_columns(
                    pl.col("name").replace_strict(mrna_ids, default="")
                ).filter(pl.col("name") != "")
            case "CDS":
                cds_ids = dict(self.get_longest(target="CDS").iter_rows())
                df = df.with_columns(
                    pl.col("name").replace_strict(cds_ids, default="")
                ).filter(pl.col("name") != "")
            case "protein":
                protein_ids = dict(self.get_longest(target="protein").iter_rows())
                df = df.with_columns(
                    pl.col("name").replace_strict(protein_ids, default="")
                ).filter(pl.col("name") != "")

        if str(bed_path) != ".":
            df.write_csv(bed_path, separator="\t", include_header=False)

        return df

    @property
    def db_sources(self) -> List[str]:
        if self.db.shape[0] == 0:
            return []
        else:
            return self.db.select("source").to_series().unique().to_list()

    @property
    def db_types(self) -> List[str]:
        if self.db.shape[0] == 0:
            return []
        else:
            return self.db.select("type").to_series().unique().to_list()

    @property
    def db_sources_and_types(self) -> List[Tuple[str, str]]:
        if self.db.shape[0] == 0:
            return []
        else:
            return list(self.db.select(["source", "type"]).unique().iter_rows())

    @property
    def db_seqids(self) -> List[str]:
        return self.db["seqid"].unique().sort().to_list()
