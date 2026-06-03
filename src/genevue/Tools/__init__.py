from pathlib import Path

import typer
from genevue import setup_rich_logger, console
from genevue.Tools.filesystem import iter_path
from typing import Annotated, Optional

logger = setup_rich_logger(__name__, console)

app_tools = typer.Typer()


@app_tools.command(name="list-path")
def cmd_listpath(
    output: Annotated[str, typer.Argument],
    target_path: Optional[str] = typer.Option(None, "-t", "--target"),
    is_include_dir: bool = typer.Option(False, "-d", "--dir", is_flag=True),
    is_recursive: bool = typer.Option(False, "-r", "--recursive", is_flag=True),
):
    target_path_path = Path.cwd() if not target_path else Path(target_path)

    iter_path(Path(output), target_path_path, is_include_dir, not is_recursive)
