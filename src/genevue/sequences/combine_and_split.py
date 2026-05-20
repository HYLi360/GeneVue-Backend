from pathlib import Path
from typing import List

from genevue import console, setup_rich_logger

logger = setup_rich_logger(__name__, console)


def simple_combine(
    in_path: List[Path],
    out_path: Path,
):
    out_file_handler = open(out_path, "w")

    for in_file in in_path:
        if not in_file.exists():
            logger.warning(f"File {in_file} does not exist. Skipping.")
            continue
        else:
            with open(in_file) as f:
                for line in f:
                    out_file_handler.write(line)

    out_file_handler.close()
