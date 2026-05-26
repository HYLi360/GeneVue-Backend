#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

import typer
from genevue.configure import Configure
from genevue.Remote.NCBI_DATASETS_API import Includes, Datasets4Genome

app_remote = typer.Typer()


@app_remote.command()
def download_datasets(
    accessions_path: str,
    chunk_size: int = typer.Option(
        20, "-c", "--chunk-size", help="Accessions per download batch (max 100)."
    ),
    max_concurrent: int = typer.Option(
        3, "-p", "--parallel", help="Max parallel batch downloads (1 = sequential)."
    ),
    makesymlink: bool = typer.Option(False, "-l", "--link", is_flag=True),
):
    accessions = []
    with open(accessions_path) as f:
        for line in f:
            accessions.append(line.strip())

    configure = Configure()
    apikey = configure.get_apikey("NCBI")

    Datasets4Genome(
        accessions,
        include=[Includes.gff3, Includes.cds, Includes.pep],
        chunk_size=chunk_size,
        max_concurrent=max_concurrent,
        apikey=apikey,
        generate_symlinks=makesymlink,
    ).download()
