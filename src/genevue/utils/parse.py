#   Copyright (C) 2025-2026, HYLi360.
#   Free software distributed under the terms of the GNU GPL-3.0 license,
#   and comes with ABSOLUTELY NO WARRANTY.
#   See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

import tempfile
from functools import singledispatch
from pathlib import Path
from typing import Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rich.text import Text

from genevue import console
from genevue.utils import pairwise_re, mcscan_info_line, yn00_res_re
import json


@singledispatch
def record_filter(
    records: Path | list[SeqRecord] | dict[str, SeqRecord], record_id_list: list[str]
) -> dict[str, SeqRecord]: ...


@record_filter.register
def _(records: Path, record_id_list: list) -> dict[str, SeqRecord]:
    if len(record_id_list) == 0:
        raise RuntimeError("record_id_list is empty")
    # resolve and check exists
    path = records.resolve()
    if not path.exists():
        raise FileNotFoundError(f"Can not find sequence file: {str(records)}")
    # read the file
    records = {record.id: record for record in SeqIO.parse(str(path), "fasta")}
    # filt
    new_records = {}
    for record_id in record_id_list:
        try:
            new_records[record_id] = records[record_id]
        except KeyError:
            continue
    return new_records


@record_filter.register
def _(records: list, record_id_list: list) -> dict[str, SeqRecord]:
    if len(record_id_list) == 0:
        raise RuntimeError("record_id_list is empty")
    # expand the list
    records = {record.id: record for record in records}
    # filt
    new_records = {}
    for record_id in record_id_list:
        try:
            new_records[record_id] = records[record_id]
        except KeyError:
            continue
    return new_records


@record_filter.register
def _(records: dict, record_id_list: list) -> dict[str, SeqRecord]:
    if len(record_id_list) == 0:
        raise RuntimeError("record_id_list is empty")
    # filt
    new_records = {}
    for record_id in record_id_list:
        try:
            new_records[record_id] = records[record_id]
        except KeyError:
            continue
    return new_records


def blast6reader(
    blast_6_result_path: str | Path,
    lo_bitscore: Optional[float] = None,
    lo_evalue: Optional[float] = None,
    hi_bitscore: Optional[float] = None,
    hi_evalue: Optional[float] = None,
    genels1: Optional[list] = None,
    genels2: Optional[list] = None,
    reverse: bool = False,
) -> pd.DataFrame:
    blast = pd.read_csv(
        blast_6_result_path,
        sep="\t",
        header=None,
        comment="#",
        dtype={
            0: str,
            1: str,
            2: float,
            3: int,
            4: int,
            5: int,
            6: int,
            7: int,
            8: int,
            9: int,
            10: float,
            11: float,
        },
    )
    if reverse:
        blast[0], blast[1] = blast[1], blast[0]
        genels1, genels2 = genels2, genels1

    blast.columns = [
        "gene1",
        "gene2",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]

    if lo_bitscore is not None:
        blast = blast[blast["bitscore"] >= lo_bitscore]

    if lo_evalue is not None:
        blast = blast[blast["evalue"] >= lo_evalue]

    if hi_bitscore is not None:
        blast = blast[blast["bitscore"] <= hi_bitscore]

    if hi_evalue is not None:
        blast = blast[blast["evalue"] <= hi_evalue]

    if genels1 is not None:
        blast = blast[blast["gene1"].isin(genels1)]

    if genels2 is not None:
        blast = blast[blast["gene2"].isin(genels2)]

    return blast.sort_values(["gene1", "gene2"]).reset_index(drop=True)


def coll_res_reader(mcscan_res_path: str) -> pd.DataFrame:
    with (
        open(mcscan_res_path, "r") as mcscan_res,
        tempfile.NamedTemporaryFile("w+") as tmpf,
    ):
        bkinfo = []
        for line in mcscan_res:
            if mcscan_info_line.match(line):
                bkinfo = list(mcscan_info_line.match(line).groups())
                continue
            else:
                n = line.strip().split()
                n.extend(bkinfo)
                # n has 12 columns
                #     0      1     2      3         4     5     6      7 8    9   10              11
                # gene1 order1 gene2 order2 direction bkidx score pvalue N chr1 chr2 total_direction
                # return dataframe:
                # bkidx total_direction score pvalue gene1 chr1 order1 gene2 chr2 order2 direction
                # which is
                # [5, 11, 6, 7, 0, 9, 1, 2, 10, 3, 4]
                tmpf.write(
                    "\t".join(
                        [
                            n[5],
                            n[11],
                            n[6],
                            n[7],
                            n[0],
                            n[9],
                            n[1],
                            n[2],
                            n[10],
                            n[3],
                            n[4],
                        ]
                    )
                    + "\n"
                )
        tmpf.seek(0)
        return pd.read_csv(
            tmpf.name,
            sep="\t",
            header=None,
            names=[
                "bkidx",
                "total_direction",
                "score",
                "pvalue",
                "gene1",
                "chr1",
                "order1",
                "gene2",
                "chr2",
                "order2",
                "direction",
            ],
            dtype={
                "bkidx": int,
                "total_direction": str,
                "score": int,
                "pvalue": float,
                "gene1": str,
                "chr1": str,
                "order1": int,
                "gene2": str,
                "chr2": str,
                "order2": int,
                "direction": int,
            },
        )


@singledispatch
def sequence_cutter(seq: Seq | SeqRecord, start: int, end: int) -> Seq | SeqRecord: ...


@sequence_cutter.register
def _sequence_cutter_base(seq: Seq, start: int, end: int) -> Seq:
    if (start <= 0) or (end <= 0):
        raise ValueError(
            f"Start or end position should bigger than 0, but received {start} and {end}"
        )
    else:
        console.log(f"Base sequence length: {len(seq)}.")
        if start < end:
            console.log(
                f"Normal cutting mode. cutting range: {min(len(seq), start)} ~ {min(len(seq), end)}."
            )
            return seq[(min(len(seq), start) - 1) : min(len(seq), end)]
        elif start == end:
            console.log(
                f"Point cutting mode. Input range: {start} ~ {end}. Cutting point: {min(len(seq), start)}."
            )
            return Seq(seq[(min(len(seq), start) - 1)])
        else:
            raise ValueError(
                f"End should bigger than Start, but received {start} and {end}"
            )


@sequence_cutter.register
def _(seq: SeqRecord, start: int, end: int):
    seq, sid, name, description, dbxrefs = (
        seq.seq,
        seq.id,
        seq.name,
        seq.description,
        seq.dbxrefs,
    )
    new_seq = _sequence_cutter_base(seq, start, end)
    new_description = " ".join([f"{start}-{end}", description])
    return SeqRecord(new_seq, sid, name, new_description, dbxrefs)


def codeml_pairwise(rst_path: str) -> dict:
    res_name = ["N", "S", "dN", "dS", "omega", "t"]

    with open(rst_path, encoding="utf-8", errors="replace") as rst:
        for line in rst.readlines():
            if pairwise_re.match(line):
                return dict(
                    zip(res_name, [float(i) for i in pairwise_re.match(line).groups()])
                )
    raise ValueError("internal error: falled to parse result of pairwise")


def yn00_result(res_path: str):
    res_name = ["S", "N", "t", "kappa", "omega", "dN", "dNSE", "dS", "dSSE"]

    with open(res_path, encoding="utf-8", errors="replace") as res:
        for line in res.readlines():
            line = line.strip()
            if yn00_res_re.match(line):
                return dict(
                    zip(res_name, [float(i) for i in yn00_res_re.match(line).groups()])
                )
    raise ValueError("internal error: falled to parse result of yn00")


def is_seq_equal(record1: SeqRecord, record2: SeqRecord) -> int:
    if record1.seq != record2.seq:
        console.log(f"{record1.id} and {record2.id} have different sequence data.")
        return 2
    else:
        if record1.id != record2.id:
            console.log(
                f"{record1.id} and {record2.id} have the same sequence data, but have different id."
            )
            return 1
        else:
            console.log(
                f"{record1.id} and {record2.id}... are actually the same record."
            )
            return 0


def species_italic_name(species_name: str) -> Text:
    species_name_list: list[str] = species_name.strip().split()
    res: Text = Text()
    label_first = True
    for name in species_name_list:
        if not label_first:
            res.append(Text(" "))
        if name not in ["subsp.", "ssp.", "var.", "subvar.", "f.", "subf"]:
            res.append(Text(name, style="italic"))
        else:
            res.append(Text(name))
        label_first = False

    return res


def asm_data_report_jsonl_parser(jsonl_path: Path):
    # Return format (for every accession):
    # [
    #   {
    #     accession: ...,
    #     source: ...,
    #     lineage:
    #     {
    #       common_name: ...,
    #       species: ...,
    #       genus: ...,
    #       family: ...,
    #       order: ...,
    #     },
    #     asminfo:
    #     {
    #       name: ...,
    #       level: ...,
    #       type ...,
    #       atgccount: ...,
    #       gccount: ...,
    #       size: ...,
    #       ungappedsize: ...,
    #       contig_count: ...,
    #       contig_l50: ...,
    #       contig_n50: ...,
    #       scaffold_count: ...,
    #       scaffold_l50: ...,
    #       scaffold_n50: ...,
    #       chromosome_count: ...,
    #       organelles_count: ...,
    #     },
    #     annoinfo:
    #     {
    #       name: ...,
    #       genes_total_count: ...,
    #       genes_protein_coding: ...,
    #       genes_non_coding: ...,
    #       genes_pseudogene: ...,
    #       genes_other: ...,
    #       busco_lineage: ...,
    #       busco_version: ...,
    #       busco_result:
    #       {
    #         total: ...,
    #         complete: ...,
    #         complete_ratio: ...,
    #         singlecopy: ...,
    #         singlecopy_ratio: ...,
    #         duplicated: ...,
    #         duplicated_ratio: ...,
    #         fragmented: ...,
    #         fragmented_ratio: ...,
    #         missing: ...,
    #         missing_ratio: ...,
    #     },
    #   },
    #   ...
    # ]
    res = []
    with open(jsonl_path) as jsonl:
        for json_obj in jsonl:
            json_obj = json_obj.strip()
            if json_obj:
                tgt_d = {}
                src_d = json.loads(json_obj)
                tgt_d["accession"] = src_d["accession"]
                tgt_d["source"] = src_d["sourceDatabase"]

                #!todo get lineage by entrez

                # asm info
                tgt_d["asminfo"] = {}
                tgt_d["asminfo"]["name"] = src_d["assemblyInfo"]["assemblyName"]
                tgt_d["asminfo"]["level"] = src_d["assemblyInfo"]["assemblyLevel"]
                tgt_d["asminfo"]["type"] = src_d["assemblyInfo"]["assemblyType"]
                tgt_d["asminfo"]["atgccount"] = src_d["assemblyStats"]["atgcCount"]
                tgt_d["asminfo"]["gccount"] = src_d["assemblyStats"]["gcCount"]
                tgt_d["asminfo"]["size"] = src_d["assemblyStats"]["totalSequenceLength"]
                tgt_d["asminfo"]["ungappedsize"] = src_d["assemblyStats"][
                    "totalUngappedLength"
                ]
                tgt_d["asminfo"]["contig_count"] = src_d["assemblyStats"][
                    "numberOfContigs"
                ]
                tgt_d["asminfo"]["contig_n50"] = src_d["assemblyStats"]["contigN50"]
                tgt_d["asminfo"]["contig_l50"] = src_d["assemblyStats"]["contigL50"]
                tgt_d["asminfo"]["scaffold_count"] = src_d["assemblyStats"][
                    "numberOfScaffolds"
                ]
                tgt_d["asminfo"]["scaffold_n50"] = src_d["assemblyStats"]["scaffoldN50"]
                tgt_d["asminfo"]["scaffold_l50"] = src_d["assemblyStats"]["scaffoldL50"]
                tgt_d["asminfo"]["chromosome_count"] = src_d["assemblyStats"][
                    "totalNumberOfChromosomes"
                ]
                tgt_d["asminfo"]["organelles_count"] = src_d["assemblyStats"][
                    "numberOfOrganelles"
                ]
                tgt_d["annoinfo"]["name"] = src_d["annotationInfo"]["name"]
                tgt_d["annoinfo"]["genes_total_count"] = src_d["annotationInfo"][
                    "geneCounts"
                ]["total"]
                tgt_d["annoinfo"]["genes_protein_coding"] = src_d["annotationInfo"][
                    "geneCounts"
                ]["proteinCoding"]
                tgt_d["annoinfo"]["genes_non_coding"] = src_d["annotationInfo"][
                    "geneCounts"
                ]["nonCoding"]
                tgt_d["annoinfo"]["genes_pseudogene"] = src_d["annotationInfo"][
                    "geneCounts"
                ]["pseudogene"]
                tgt_d["annoinfo"]["genes_other"] = src_d["annotationInfo"][
                    "geneCounts"
                ]["other"]
            res.append(tgt_d)
