"""
src/genevue/plot/CollScatPlot.py

(C) 2026 HYLi360. All rights reserved.

see LICENSE in LICENSE
see side-package LICENSEs (if used) in LICENSE_OF_SIDE_PACKAGES

--------------------
Plot BLASTp/MCScans dotplot, based on paired genes, bed, and others (if you have):
- Genome sequences, for confirm the length of chromosomes.
- Genes masker, select the genes in each genome.
- An option, how the points locale themselves (by order, or middle-point of gene. Order by default)
"""

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.preprocessing import MinMaxScaler

from genevue.Utils.Parse import blast6reader


class CollScatPlot:
    def __init__(
        self,
        blast_res_path: Path,
        bedx: pd.DataFrame,
        bedy: pd.DataFrame,
        genome_namex: str = "",
        genome_namey: str = "",
        pep_label_from: Literal["gene", "mRNA", "CDS", "protein_id"] = "gene",
        locale_method: Literal["order", "mid-point"] = "order",
        max_repeat_times: int = 5,
    ):
        self.blastres = blast6reader(blast_res_path)

        self.bedx = bedx
        self.bedy = bedy

        self.genome_namex = genome_namex
        self.genome_namey = genome_namey
        self.pep_label_from = "mRNA"
        self.locale_method = locale_method
        self.max_repeat_times = max_repeat_times

        self.fig, self.ax = plt.subplots(figsize=(20, 20))
        self.fig.subplots_adjust(left=0.07, right=0.97, top=0.93, bottom=0.03)

    def draw(self):
        self._coord_transform()
        self._initialise_canvas()
        self._locale_points()
        self._plot_points()

    def _coord_transform(self) -> None:
        """Transform bed to coord on canvas."""
        to_one = MinMaxScaler(feature_range=(0, 1000))

        if self.locale_method == "order":
            rangex = (
                self.bedx.groupby("chrom").agg({"chrom": "count"}).cumsum().to_numpy()
            )
            rangey = (
                self.bedy.groupby("chrom").agg({"chrom": "count"}).cumsum().to_numpy()
            )
        else:
            rangex = (
                self.bedx.groupby("chrom").agg({"chrEnd": "max"}).cumsum().to_numpy()
            )

            rangey = (
                self.bedy.groupby("chrom").agg({"chrEnd": "max"}).cumsum().to_numpy()
            )

        self.scalex: float | int = 1000 / rangex[-1]
        self.scaley: float | int = 1000 / rangey[-1]
        self.tickx = to_one.fit_transform(rangex).transpose()[0]
        self.ticky = to_one.fit_transform(rangey).transpose()[0]

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
        self.ax.set_xticks(self.tickx, [])
        self.ax.set_yticks(self.ticky, [])
        self.ax.set_xlim(left=0, right=1000)
        self.ax.set_ylim(bottom=0, top=1000)
        self.ax.xaxis.tick_top()
        self.ax.grid(visible=True, linewidth=0.5, alpha=1)

        # chromosome labels
        for x, chrname in zip(self.tickmidx, self.bedx["chrom"].unique()):
            self.ax.text(x, 1000, chrname, fontsize=15, ha="center", va="bottom")
        for y, chrname in zip(self.tickmidy, self.bedy["chrom"].unique()):
            self.ax.text(0, y, chrname, fontsize=15, ha="right", va="center")

        # species labels
        secax = self.ax.secondary_xaxis("top")
        secax.set_xticks([])
        secax.set_xlabel(
            self.genome_namex,
            fontsize=20,
            family="DeJaVu Sans",
            style="italic",
            ha="center",
            va="bottom",
        )
        self.ax.set_ylabel(
            self.genome_namey,
            fontsize=20,
            family="DeJaVu Sans",
            style="italic",
            ha="right",
            va="center",
        )
        return

    def _locale_points(self):
        if self.locale_method == "order":
            self.bedx["locx"] = np.arange(self.bedx.shape[0]) * self.scalex
            self.blastres = pd.merge(
                left=self.blastres,
                right=self.bedx,
                left_on="gene1",
                right_on="name",
            )
            self.bedy["locy"] = np.arange(self.bedy.shape[0]) * self.scaley
            self.blastres = pd.merge(
                left=self.blastres,
                right=self.bedy,
                left_on="gene2",
                right_on="name",
            )
        else:
            pass

    def _plot_points(self):
        self.ax.plot(
            self.blastres["locx"],
            self.blastres["locy"],
            marker=".",
            ls="",
            markersize=0.2,
        )
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
