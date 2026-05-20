import sys
import urllib.parse
from pathlib import Path
from typing import Optional

from genevue import console, setup_rich_logger

logger = setup_rich_logger(__name__, console)


def _url_style_decode(string: str) -> str:
    """
    Decode string which in URL-style coding.
    """
    return urllib.parse.unquote(string)


def _parse_key_equal_value(string: str) -> dict[str, str]:
    if "=" not in string:
        return {}
    else:
        key_str, value_str = string.strip().split("=")
        return {_url_style_decode(key_str): _url_style_decode(value_str)}


def pre_check(
    gff_path: Path,
    genome_seq: Optional[Path] = None,
    cds_seq: Optional[Path] = None,
    cds_id: Optional[Path] = None,
    pep_seq: Optional[Path] = None,
    pep_id: Optional[Path] = None,
):
    logger.info("GFF3 format pre-check start.")
    logger.info(f"GFF3 path: {gff_path}")
    sym_of_pass = True

    step_counter = 1
    logger.info(f"Pre-check step {step_counter}: File head.")
    with open(gff_path) as f:
        header = f.read(15)
        if header != "##gff-version 3":
            logger.error(
                f"This GFF3 file has wrong header! except '##gff-version 3' but get {header}"
            )
            logger.error("Pre-check process halt.")
            sys.exit(0)
    logger.info(f"Step {step_counter} passed.")

    step_counter += 1
    logger.info(f"Pre-check step {step_counter}: Separator.")
    with open(gff_path) as f:
        for line in f:
            if not line.startswith("#"):
                if "\t" not in line:
                    logger.warning(
                        "This GFF3 file seems using 'space' (' ') to separate, which is not comply with the specification."
                    )
                    sym_of_pass = False

    if sym_of_pass:
        logger.info(f"Step {step_counter} passed.")
    sym_of_pass = True

    step_counter += 1
    logger.info(f"Pre-check step {step_counter}: c")


def parse(gff_path: Path):
    pass
