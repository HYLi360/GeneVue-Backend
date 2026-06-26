#  Copyright (c) 2026 HYLi360. All rights reserved.
#
#  see license in LICENSE
#  see side-package licenses in LICENSE_OF_SIDE_PACKAGES
from typing import Literal

import typer

from genevue.Evo.Coll.powercoll import Collinearity
from genevue.Formats.GFF3 import GFF3

app_evo = typer.Typer()


@app_evo.command(name="powercoll")
def cmd_powercoll(
    gff3_left_path: str,
    gff3_right_path: str,
    blastres_path: str,
    output_basename: str,
    blastid_from: Literal["gene", "mRNA", "CDS", "protein_id"],
):
    gff3_left = GFF3()
    gff3_left.load_from_file(gff3_left_path)

    gff3_right = GFF3()
    gff3_right.load_from_file(gff3_right_path)

    c = Collinearity(gff3_left)
