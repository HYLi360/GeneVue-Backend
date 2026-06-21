import re
from pathlib import Path
from typing import Literal, List, Optional, overload
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from genevue.Sequences.FASTA import FASTA
from genevue.Sequences.combine_and_split import simple_combine

import typer


# Universal functions.
def is_seq_type(
    sequence: str | Seq | SeqRecord, seqtype: Literal["DNA", "RNA", "protein"]
):
    char_dict = defaultdict(int)
    for char in sequence:
        char_dict[char.lower()] += 1

    if char_dict["a"] + char_dict["t"] + char_dict["g"] + char_dict["c"] + char_dict[
        "n"
    ] == len(sequence):
        if seqtype == "DNA":
            return True
        return False
    elif char_dict["a"] + char_dict["u"] + char_dict["g"] + char_dict["c"] + char_dict[
        "n"
    ] == len(sequence):
        if seqtype == "RNA":
            return True
        return False
    else:
        if seqtype == "protein":
            return True
        return False


def remove_any_spacechar(sequence: str):
    return re.sub(r"\s+", "", sequence)


def complement(sequence: str) -> str:
    return str(Seq(sequence).complement())


def reverse(sequence: str) -> str:
    return sequence[::-1]


def reverse_complement(sequence: str) -> str:
    return complement(reverse(sequence))


def frame6trans(sequence: str | Seq | SeqRecord):
    pass


app_sequence = typer.Typer(name="sequence", help="sequence tool")


@app_sequence.command(name="extract")
def cmd_extract(
    tgt_path: str,
    seq_path: str,
    res_path: str,
    tgt_type: Literal[*FASTA.support_target_format()],
):
    extractor = FASTA(Path(seq_path).resolve())
    extractor.filter(Path(tgt_path).resolve(), tgt_type, Path(res_path).resolve())


@app_sequence.command(name="combine")
def cmd_combine(
    in_path: List[str],
    out_path: str,
):
    in_path_ls = [Path(path_str).resolve() for path_str in in_path]
    out_path = Path(out_path).resolve()
    simple_combine(in_path_ls, out_path)


@app_sequence.command(name="getid")
def cmd_getid(
    seq_path: str,
    out_path: str,
    idx_path: Optional[str] = None,
):
    fasta = FASTA(Path(seq_path).resolve())

    if idx_path:
        fasta.sequence_idx_path = Path(idx_path).resolve()

    fasta.export_seq_ids(Path(out_path).resolve())
