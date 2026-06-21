import re
from pathlib import Path
from typing import Literal, Iterable, Optional, Annotated
import typer
from genevue.Utils.FileSystem import gzip2bgzip, MAGIC_NUMBER_HEAD

pairwise_re = re.compile(
    r"^\s*\d+\s+\d+\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"
)
"""N S dN dS dN/dS(omega) t"""

mcscan_info_line = re.compile(
    r"^#\sAlignment\s(\d+):\sscore=([\d.]+)\spvalue=([\d.]+)\sN=(\d+)\s(\d+)&(\d+)\s(\w+)"
)
"""bkidx, score, pvalue, N, chr1, chr2, direction"""

yn00_res_re = re.compile(
    r"^2\s+1\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\+\-\s+\-?([\d\.]+)\s+\-?([\d\.]+)\s+\+\-\s+\-?([\d\.]+)"
)
"""S N t kappa omega dN dNSE dS dSSE"""

app_file = typer.Typer(name="file")


def header_detector(file_path: Path, format_type: Literal["gz", "bgz"]) -> bool:
    magic_number = {"gz": (4, b"\x1f\x8b\x08\x00"), "bgz": (4, b"\x1f\x8b\x08\x04")}
    leng, byte = magic_number[format_type]

    with open(file_path, "rb") as f:
        if byte == f.read(leng):
            return True
        return False


def deduplicate(iterable_obj: Iterable):
    return list(dict.fromkeys(iterable_obj))


@app_file.command()
def gz2bgz(
    gzip_path: str,
    bgzip_path: Annotated[Optional[str], typer.Argument()] = None,
    chunk_size: int = typer.Option(4096, "-c", "--chunk_size"),
):
    if bgzip_path is None:
        bgzip_path = Path(gzip_path).parent / f"{Path(gzip_path).stem}.bgz"
    gzip2bgzip(Path(gzip_path), Path(bgzip_path), chunk_size)


@app_file.command()
def echo_filetype(file_path: str):
    with open(Path(file_path), "rb") as f:
        for ftype, (rlen, byte) in MAGIC_NUMBER_HEAD.items():
            if byte == f.read(rlen):
                print(ftype)
                return
            f.seek(0)
    print("Not reg")
