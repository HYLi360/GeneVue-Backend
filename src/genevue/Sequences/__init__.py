from pathlib import Path
from typing import Literal, List, Optional
from genevue.Sequences.FASTA import FASTA
from genevue.Sequences.combine_and_split import simple_combine

import typer

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
