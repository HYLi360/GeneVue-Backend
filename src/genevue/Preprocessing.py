#  Copyright (c) 2026 HYLi360. All rights reserved.
#
#  see license in LICENSE
#  see side-package licenses in LICENSE_OF_SIDE_PACKAGES
from pathlib import Path

import typer

from genevue.Formats.GFF3 import GFF3

app_preprocessing = typer.Typer()


@app_preprocessing.command(name="gff3_to_parquet")
def gff3_to_parquet(gff3_path: str, out_path: str = ""):
    gff3_path = Path(gff3_path)

    if str(out_path) == "":
        out_path = Path(gff3_path.parent / f"{gff3_path.stem}.parquet")

    gff3 = GFF3()
    gff3.load_from_file(gff3_path)
    gff3.to_parquet(out_path)
