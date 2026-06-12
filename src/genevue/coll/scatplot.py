"""
Plot BLAST dotplot, based on paired genes, bed, and others (if you have):
- Genome sequences, for confirm the length of chromosomes.
- Genes masker, select the genes in each genome.
- An option, how the points locale themselves (by order, or middle-point of gene. Order by default)
"""

from pathlib import Path
from typing import Literal, Optional, List

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from numpy.typing import NDArray
from sklearn.preprocessing import MinMaxScaler

from genevue.utils.parse import bed_reader, blast6reader
from genevue.GXF.GFF3tools import BlazingGFF3


class CollScatPlot:
    def __init__(
        self,
        blast_res_path: Path,
        out_plot_path: Path,
        bedx_path: Path,
        bedy_path: Path,
        genome_namex: str = "",
        genome_namey: str = "",
        locale_method: Literal["order", "mid-point"] = "order",
        max_repeat_times: int = 5,
    ):
        self.blastres = blast6reader(blast_res_path)
        self.bedx = bed_reader(bedx_path)
        self.bedy = bed_reader(bedy_path)

        self.out_plot_path = out_plot_path
        self.genome_namex = genome_namex
        self.genome_namey = genome_namey
        self.locale_method = locale_method
        self.max_repeat_times = max_repeat_times

        self.left = 0.07
        self.right = 0.97
        self.top = 0.93
        self.bottom = 0.03

        self.fig, self.ax = plt.subplots()

    def plot(self):
        # pre-process
        # get length of chromosome, and position for every genes
        self._coord_transform()

    def _coord_transform(self) -> None:
        """Transform bed to coord on canvas."""
        to_one = MinMaxScaler(feature_range=[0, 1000])

        if self.locale_method == "order":
            rangex = (
                self.bedx.groupby("chrom")
                .agg({"chrom": "count"})["count"]
                .cumsum()
                .to_numpy()
            )
            rangey = (
                self.bedy.groupby("chrom")
                .agg({"chrom": "count"})["count"]
                .cumsum()
                .to_numpy()
            )
        else:
            rangex = (
                self.bedx.groupby("chrom")
                .agg({"chrEnd": "max"})["max"]
                .cumsum()
                .to_numpy()
            )

            rangey = (
                self.bedy.groupby("chrom")
                .agg({"chrEnd": "max"})["max"]
                .cumsum()
                .to_numpy()
            )

        self.scalex: float | int = 1000 / rangex[-1]
        self.scaley: float | int = 1000 / rangey[-1]
        self.tickx = to_one.fit_transform(rangex)
        self.ticky = to_one.fit_transform(rangey)

        # coords of chromosome labels
        self.tickmidx = [
            (self.tickx[i] + self.tickx[i + 1]) / 2 for i in range(len(self.tickx) - 1)
        ]
        self.tickmidy = [
            (self.ticky[i] + self.ticky[i + 1]) / 2 for i in range(len(self.ticky) - 1)
        ]
        return

    def _initialise_canvas(self):
        # grid (as chromosomes)
        self.ax.set_xticks(self.tickx)
        self.ax.set_yticks(self.ticky)
        self.ax.set_xlim(left=0, right=1000)
        self.ax.set_ylim(bottom=0, top=1000)
        self.ax.grid(visible=True, linewidth=0.4, alpha=0)

        # chromosome labels
        for x, chrname in zip(self.tickmidx, self.bedx["chrom"].unique()):
            self.ax.text(x, 0.90, chrname)
        for y, chrname in zip(self.tickmidy, self.bedy["chrom"].unique()):
            self.ax.text(0.10, y, chrname)

        # species labels
        self.ax.set_xlabel(self.genome_namex)
        self.ax.set_ylabel(self.genome_namey)
        return

    def _plot_points(self):
        self.ax.plot(self.blastres["locx"], self.blastres["locy"], marker=".", ls="")
        return

    def _fill_ancestor_color_bar(self):
        pass

    def show(self):
        self.fig.show()

    def export(self, outpath: Path):
        """
        "print" (or export) scatter plot.
        """
        self.fig.savefig(outpath)
