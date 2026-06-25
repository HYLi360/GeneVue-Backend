#  Copyright (c) 2026 HYLi360. All rights reserved.
#
#  see LICENSE in /LICENSE
#  see side-package LICENSEs (if used) in /LICENSE_OF_SIDE_PACKAGES


"""
a very simple tool to build newick tree.
"""

from __future__ import annotations

from typing import Dict, cast

import newick
import pandas as pd
import re


class TreeNode:
    def __init__(self, name: str, level: str):
        # replace any space to single "_" to follow the newick syntax
        self.name = re.sub(r"\s+", "_", name)
        self.level = level

        self.children: Dict[str, TreeNode] = {}
        self.length: int = 0

    def add_child(self, child: TreeNode):
        self.children[child.name] = child

    def get_or_create_child(self, name: str, level: str):
        if name not in self.children:
            self.children[name] = TreeNode(name, level)
        return self.children[name]

    def to_newick(self) -> str:
        if not self.children:
            return self.name

        child_newick = ",".join([child.to_newick() for child in self.children.values()])

        if self.level == "root":
            return f"({child_newick})"
        else:
            return f"({child_newick}){self.name}"


class SimpleTree:
    def __init__(self):
        self.data = pd.DataFrame()
        self.tree = TreeNode("root", "root")

    def add_path(self, path: dict):
        self.data = pd.concat([self.data, pd.DataFrame(path)])

    def add_paths(self, paths: pd.DataFrame):
        self.data = pd.concat([self.data, paths])

    def build(self):
        self.data = self.data.fillna("?")
        for _, series in self.data.iterrows():
            current = self.tree

            for level, name in series.items():
                current = current.get_or_create_child(cast(str, name), cast(str, level))

    def to_newick(self):
        return self.tree.to_newick()

    def draw_as_ascii(self):
        return newick.loads(self.tree.to_newick())[0].ascii_art()
