#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

import typer
from genevue.Remote.NCBI_DATASETS_API import Includes, Datasets4Genome

app_remote = typer.Typer()


@app_remote.command()
def download_datasets(accessions_path):
    accessions = []
    with open(accessions_path) as f:
        for line in f:
            accessions.append(line.strip())

    Datasets4Genome(
        accessions, include=[Includes.gff3, Includes.cds, Includes.pep]
    ).download()
