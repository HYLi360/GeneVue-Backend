from typing import Literal

import pandas as pd
import typer

from genevue.Taxonomy.SimpleTree import SimpleTree

app_taxonomy = typer.Typer()


@app_taxonomy.command(name="simple-tree")
def cmd_simple_tree(
    table_path: str,
    sep: Literal["tab", "comma"],
    newick_path: str,
    print_ascii: bool = typer.Option(False, "-p", "--print", is_flag=True),
    first_line_as_columns_name: bool = typer.Option(
        True, "-f", "--first", is_flag=True
    ),
):
    sep_str = {"tab": "\t", "comma": ","}.get(sep)

    if first_line_as_columns_name:
        df = pd.read_csv(table_path, sep=sep_str)
    else:
        df = pd.read_csv(table_path, sep=sep_str, header=None)

    tree = SimpleTree()
    tree.add_paths(df)
    tree.build()

    if print_ascii:
        print(tree.draw_as_ascii())

    with open(newick_path, "w") as f:
        f.write(tree.to_newick())
