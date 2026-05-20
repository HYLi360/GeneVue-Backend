import re
from pathlib import Path
from typing import Literal, Iterable

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


def header_detector(file_path: Path, format_type: Literal["gz", "bgz"]) -> bool:
    magic_number = {"gz": (4, b"\x1f\x8b\x08\x00"), "bgz": (4, b"\x1f\x8b\x08\x04")}
    leng, byte = magic_number[format_type]

    with open(file_path, "rb") as f:
        header = f.read(leng)
        if byte == header:
            return True
        return False


def deduplicate(iterable_obj: Iterable):
    return list(dict.fromkeys(iterable_obj))


def chunk_iter(data, chunk_size):
    for i in range(0, len(data), chunk_size):
        yield data[i : i + chunk_size]
