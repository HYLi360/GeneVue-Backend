import tomllib

import typer

app_batch = typer.Typer()


@app_batch.command()
def makeblastdb():
    pass


@app_batch.command()
def blast(
    conf_toml: str,
):
    pass
