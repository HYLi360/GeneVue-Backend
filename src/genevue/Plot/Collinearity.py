"""
src/genevue/plot/Collinearity.py

(C) 2026 HYLi360. All rights reserved.

see LICENSE in /LICENSE
see side-package LICENSEs (if used) in /LICENSE_OF_SIDE_PACKAGES

--------------------
A more powerful collinearity plot tool.
"""

from typing import Optional, Literal, Iterable, Tuple, List, Dict

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.path import Path as MplPath
from matplotlib.text import Text

# Font sizes
GENE_FS = 7
SIDE_FS = 9


class MacroCollPlot:
    def __init__(
        self,
        genome_structure_top: pd.DataFrame,
        genome_structure_bottom: pd.DataFrame,
        species_name_top: str,
        species_name_bottom: str,
        anchor: zip,
        h: tuple[float, float] = (0.6, 0.4),
        x: tuple[float, float] = (0.2, 0.9),
        spacing=0.1,
        collbox_style: Literal["polygon", "bezier"] = "bezier",
        chr_namels_top: Optional[list] = None,
        chr_namels_bottom: Optional[list] = None,
        savefile: Optional[str] = None,
    ):
        self.genome_structure_top = genome_structure_top
        self.genome_structure_bot = genome_structure_bottom
        self.species_name_top = species_name_top
        self.species_name_bot = species_name_bottom
        self.anchor = anchor
        self.h = h
        self.x = x
        self.spacing = spacing
        self.collbox_style = collbox_style
        self.chr_namels_top = chr_namels_top
        self.chr_namels_bot = chr_namels_bottom
        self.savefile = savefile

    def draw(self):
        self._chro_init()
        self._collbox_init()

        palette = plt.get_cmap("tab10")

        ax = plt.gca()

        # Chromosomes
        for start, length, chrname in self.chrposls_top:
            ax.add_patch(
                patches.Rectangle(
                    xy=(start, self.h[0] - 0.003),
                    width=length,
                    height=0.006,
                    fill=True,
                    fc=palette(0),
                    ec="none",
                )
            )
            ax.add_patch(
                patches.Rectangle(
                    xy=(start, self.h[0] - 0.01),
                    width=length,
                    height=0.02,
                    fill=False,
                    ec=(0.2, 0.2, 0.2, 0.2),
                )
            )
            ax.text(
                x=start + length / 2,
                y=self.h[0] + 0.025,
                s=chrname,
                ha="center",
                va="center",
                backgroundcolor="white",
            )

        for start, length, chrname in self.chrposls_bot:
            ax.add_patch(
                patches.Rectangle(
                    xy=(start, self.h[1] - 0.003),
                    width=length,
                    height=0.006,
                    fill=True,
                    fc=palette(1),
                    ec="none",
                )
            )
            ax.add_patch(
                patches.Rectangle(
                    xy=(start, self.h[1] - 0.01),
                    width=length,
                    height=0.02,
                    fill=False,
                    ec=(0.2, 0.2, 0.2, 0.2),
                )
            )
            ax.text(
                x=start + length / 2,
                y=self.h[1] - 0.025,
                s=chrname,
                ha="center",
                va="center",
                backgroundcolor="white",
            )

        # Add collboxes
        if self.collbox_style == "bezier":
            for collbox in self.collbox:
                collbox_patch = _get_bezier_ribbon(
                    collbox[0],
                    collbox[1],
                    collbox[2],
                    collbox[3],
                    self.h[0],
                    self.h[1],
                    "grey",
                )
                ax.add_patch(collbox_patch)
        if self.collbox_style == "polygon":
            for collbox in self.collbox:
                collbox_patch = _get_polygon_ribbon(
                    collbox[0],
                    collbox[1],
                    collbox[2],
                    collbox[3],
                    self.h[0],
                    self.h[1],
                    "grey",
                )
                ax.add_patch(collbox_patch)

        # Add species name
        text1 = Text(
            x=0.12,
            y=self.h[0],
            s=self.species_name_top,
            size=15,
            ha="center",
            va="center",
            fontstyle="italic",
        )
        text2 = Text(
            x=0.12,
            y=self.h[1],
            s=self.species_name_bot,
            size=15,
            ha="center",
            va="center",
            fontstyle="italic",
        )
        ax.add_artist(text1)
        ax.add_artist(text2)

        # Close the axis
        plt.axis("off")

        plt.show()

    def _chro_init(self):
        """Preprocessing steps prior of plotting."""
        # set chromesome list.
        chrls1 = self.genome_structure_top["chr"].drop_duplicates().tolist()
        if self.chr_namels_top is not None:
            self.chr_namels_top = [
                chro for chro in self.chr_namels_top if chro in chrls1
            ]
        else:
            self.chr_namels_top = chrls1

        chrls2 = self.genome_structure_bot["chr"].drop_duplicates().tolist()
        if self.chr_namels_bot is not None:
            self.chr_namels_bot = [
                chro for chro in self.chr_namels_bot if chro in chrls2
            ]
        else:
            self.chr_namels_bot = chrls2

        genes_num1 = len(
            self.genome_structure_top[
                self.genome_structure_top["chr"].isin(self.chr_namels_top)
            ]
        )
        genes_num2 = len(
            self.genome_structure_bot[
                self.genome_structure_bot["chr"].isin(self.chr_namels_bot)
            ]
        )

        sf1 = (self.x[1] - self.x[0]) / (genes_num1 * (1 + self.spacing))
        sf2 = (self.x[1] - self.x[0]) / (genes_num2 * (1 + self.spacing))

        if (len(self.chr_namels_top) - 1) != 0:
            space_len1 = self.spacing * genes_num1 / (len(self.chr_namels_top) - 1)
        else:
            space_len1 = 0

        if (len(self.chr_namels_top) - 1) != 0:
            space_len2 = self.spacing * genes_num2 / (len(self.chr_namels_bot) - 1)
        else:
            space_len2 = 0

        start1, length1, space_end1 = [], [], [self.x[0]]
        start2, length2, space_end2 = [], [], [self.x[0]]

        for chro in self.chr_namels_top:
            start1.append(space_end1[-1])
            length1.append(
                len(self.genome_structure_top[self.genome_structure_top["chr"] == chro])
                * sf1
            )
            space_end1.append(start1[-1] + length1[-1] + space_len1 * sf1)

        for chro in self.chr_namels_bot:
            start2.append(space_end2[-1])
            length2.append(
                len(self.genome_structure_bot[self.genome_structure_bot["chr"] == chro])
                * sf2
            )
            space_end2.append(start2[-1] + length2[-1] + space_len2 * sf2)

        self.sf_top, self.sf_bot = sf1, sf2
        self.chrposls_top, self.chrposls_bot = tuple(
            zip(start1, length1, self.chr_namels_top)
        ), tuple(zip(start2, length2, self.chr_namels_bot))
        self.startd_top, self.startd_bot = dict(zip(self.chr_namels_top, start1)), dict(
            zip(self.chr_namels_bot, start2)
        )

    def _collbox_init(self):
        res = []
        for chr1, chr2, loc1start, loc1end, loc2start, loc2end in self.anchor:
            if (str(chr1) in self.startd_top) and (str(chr2) in self.startd_bot):
                res.append(
                    [
                        self.startd_top[str(chr1)] + loc1start * self.sf_top,
                        self.startd_top[str(chr1)] + loc1end * self.sf_top,
                        self.startd_bot[str(chr2)] + loc2start * self.sf_bot,
                        self.startd_bot[str(chr2)] + loc2end * self.sf_bot,
                    ]
                )
        self.collbox = res


class MacroCollPlotMultiple:
    def __init__(
        self,
        chr_spacing_factor: float | int = 0.1,
        chr_x: Tuple[float | int, float | int] = (0.2, 0.9),
        collbox_style: Literal["polygon", "bezier"] = "bezier",
    ):
        self.chr_spacing_factor = chr_spacing_factor
        self.chr_x = chr_x
        self.collbox_style = collbox_style

        # internal attributes
        self._scale_factor = {}
        self._chr_attributes = pd.DataFrame()
        self._coll_attributes = pd.DataFrame()
        self._line_attributes = pd.DataFrame()

        self.fig, self.ax = plt.subplots(figsize=(20, 10))
        self._drown = False

    def add_genome(
        self,
        genome_name: str,
        genome_features: pd.DataFrame,
        chromosome_y: float | int,
        color,
        allow_chroms: Optional[List[str]] = None,
        rename_chroms: Optional[Dict[str, str]] = None,
    ):
        if rename_chroms is not None:
            rename_chroms_dict = rename_chroms
        else:
            rename_chroms_dict = {}

        # get chromosome name and length from any type of genome features table.
        buffer = []
        chr_length_dict = genome_features.groupby("seqid")["end"].max().to_dict()

        if allow_chroms is not None:
            chridls = [
                chro
                for chro in allow_chroms
                if chro in genome_features["seqid"].unique()
            ]

        else:
            chridls = genome_features["seqid"].unique().tolist()

        scale_factor = (self.chr_x[1] - self.chr_x[0]) / (
            sum([chr_length_dict[chro] for chro in chridls])
            * (1 + self.chr_spacing_factor)
        )

        self._scale_factor[genome_name] = scale_factor

        space_length = ((self.chr_x[1] - self.chr_x[0]) * self.chr_spacing_factor) / (
            len(chridls) - 1
        )

        s = self.chr_x[0]

        for chrid in chridls:
            buffer.append(
                {
                    "genome_name": genome_name,
                    "chr_name": rename_chroms_dict.get(chrid, chrid),
                    "y": (1 - chromosome_y),
                    "factor": scale_factor,
                    "start": s,
                    "end": s + chr_length_dict[chrid] * scale_factor,
                    "color": color,
                }
            )
            s += chr_length_dict[chrid] * scale_factor + space_length

        self._chr_attributes = pd.concat([self._chr_attributes, pd.DataFrame(buffer)])

    def add_collbox(
        self,
        genome_name1: str,
        genome_name2: str,
        anchor: pd.DataFrame,
        color,
    ):
        # read anchor file after running "coll".
        buffer = []

        factor1 = self._scale_factor[genome_name1]
        factor2 = self._scale_factor[genome_name2]
        # we only need the first 6 columns
        anchor = anchor.loc[
            :,
            [
                "chrom1",
                "chrom2",
                "loc1startbase",
                "loc1endbase",
                "loc2startbase",
                "loc2endbase",
            ],
        ]
        for _, (
            chr1,
            chr2,
            loc1startbase,
            loc1endbase,
            loc2startbase,
            loc2endbase,
        ) in anchor.iterrows():

            mask1 = (self._chr_attributes["genome_name"] == genome_name1) & (
                self._chr_attributes["chr_name"] == chr1
            )
            mask2 = (self._chr_attributes["genome_name"] == genome_name2) & (
                self._chr_attributes["chr_name"] == chr2
            )

            df1 = self._chr_attributes[mask1].reset_index()
            df2 = self._chr_attributes[mask2].reset_index()

            if (len(df1) == 0) or (len(df2) == 0):
                continue

            y1 = df1.at[0, "y"]
            chr1_start = df1.at[0, "start"]
            y2 = df2.at[0, "y"]
            chr2_start = df2.at[0, "start"]

            buffer.append(
                {
                    "genome_name1": genome_name1,
                    "chr1": chr1,
                    "y1": y1,
                    "chr1_coll_start": chr1_start + float(loc1startbase) * factor1,
                    "chr1_coll_end": chr1_start + float(loc1endbase) * factor1,
                    "genome_name2": genome_name2,
                    "chr2": chr2,
                    "y2": y2,
                    "chr2_coll_start": chr2_start + float(loc2startbase) * factor2,
                    "chr2_coll_end": chr2_start + float(loc2endbase) * factor2,
                    "color": color,
                }
            )
        self._coll_attributes = pd.concat([self._coll_attributes, pd.DataFrame(buffer)])

    def add_line(self, genome_name1, chr1, point1, genome_name2, chr2, point2, color):
        factor1 = self._scale_factor[genome_name1]
        start1 = (
            self._chr_attributes[
                (self._chr_attributes["genome_name"] == genome_name1)
                & (self._chr_attributes["chr_name"] == chr1)
            ]["start"]
            .unique()
            .item()
        )
        factor2 = self._scale_factor[genome_name2]
        start2 = (
            self._chr_attributes[
                (self._chr_attributes["genome_name"] == genome_name2)
                & (self._chr_attributes["chr_name"] == chr2)
            ]["start"]
            .unique()
            .item()
        )

        self._line_attributes = pd.concat(
            [
                self._line_attributes,
                pd.DataFrame(
                    {
                        "genome_name1": genome_name1,
                        "chr1": chr1,
                        "point1": start1 + point1 * factor1,
                        "genome_name2": genome_name2,
                        "chr2": chr2,
                        "point2": start2 + point2 * factor2,
                        "color": color,
                    }
                ),
            ],
        )

    def reorder_chromosome(self):
        pass

    def draw(self):
        self._chro_init()
        self._coll_init()
        self._line_init()

    def show(self):
        if not self._drown:
            self.draw()
        self.fig.show(warn=True)

    def save(self, path):
        self.fig.savefig(path)

    def _chro_init(self):
        for _, (
            gname,
            cname,
            y,
            factor,
            start,
            end,
            color,
        ) in self._chr_attributes.iterrows():
            self.ax = _draw_a_chromosome(self.ax, start, end, y, cname, color)

    def _coll_init(self):
        for _, (
            gname1,
            chr1,
            y1,
            collstart1,
            collend1,
            gname2,
            chr2,
            y2,
            collstart2,
            collend2,
            color,
        ) in self._coll_attributes.iterrows():
            if self.collbox_style == "bezier":
                self.ax = _draw_collbox_bezier(
                    self.ax, collstart1, collend1, collstart2, collend2, y1, y2, color
                )
            else:
                self.ax = _draw_collbox_polygon(
                    self.ax, None, None, None, None, None, y2, color
                )

    def _line_init(self):
        pass


class MicroCollPlot:
    def __init__(
        self,
        draw_method: Literal["gene_name", "gene_index", "base_range"],
        bed1: pd.DataFrame,
        bed2: pd.DataFrame,
        anchor: pd.DataFrame,
        species_name_top: str,
        species_name_bottom: str,
        list_param1: list,
        list_param2: list,
        gene_rename=None,
        y: tuple[float, float] = (0.6, 0.4),
        x: tuple[float, float] = (0.2, 0.9),
        chromosome_label: tuple[str, str] = ("", ""),
        reverse: tuple[bool, bool] = (False, False),
        colors=(
            (0.12156862745098039, 0.4666666666666667, 0.7058823529411765, 1.0),
            (1.0, 0.4980392156862745, 0.054901960784313725, 1.0),
        ),
        collbox_style: Literal["polygon", "bezier"] = "polygon",
        savefile: Optional[str] = None,
    ):
        if gene_rename is None:
            gene_rename = dict()
        bed1 = bed1.loc[:, ["gene_id", "start", "end", "strand"]].rename(
            columns={
                "gene_id": "gene_top",
                "start": "x_top_start",
                "end": "x_top_end",
                "strand": "strand_top",
            }
        )
        bed2 = bed2.loc[:, ["gene_id", "start", "end", "strand"]].rename(
            columns={
                "gene_id": "gene_bot",
                "start": "x_bot_start",
                "end": "x_bot_end",
                "strand": "strand_bot",
            }
        )
        anchor = anchor.loc[:, ["gene1", "chr1", "gene2", "chr2"]].rename(
            columns={
                "gene1": "gene_top",
                "chr1": "chr_top",
                "gene2": "gene_bot",
                "chr2": "chr_bot",
            }
        )

        anchor = pd.merge(left=anchor, right=bed1, on="gene_top", how="inner")
        anchor = pd.merge(left=anchor, right=bed2, on="gene_bot", how="inner")
        # anchor: gene_top, chr_top, gene_bot, chr_bot, x_top_start, x_top_end, strand_top, x_bot_start, x_bot_end, strand_bot
        self.anchor = anchor

        self.draw_method = draw_method
        self.list_param1 = list_param1
        self.list_param2 = list_param2
        self.gene_rename = gene_rename

        self.species_name_top = species_name_top
        self.species_name_bot = species_name_bottom

        self.ytop, self.ybot = y
        self.xleft, self.xright = x
        self.chr_label_top, chr_label_bot = chromosome_label
        self.reverse_top, self.reverse_bot = reverse
        self.color_top, self.color_bot = colors[0], colors[1]
        self.collbox_style = collbox_style
        self.savefile = savefile

    def show(self):
        self._range_init()
        self._arrows_collboxes_init()

        fig, ax = plt.subplots(figsize=(7.0, 4.0), dpi=120, constrained_layout=True)

        # Add DNA strand
        ax.add_patch(
            patches.Rectangle(
                xy=(self.xleft, self.ytop - 0.002),
                width=self.xright - self.xleft,
                height=0.004,
                fc="grey",
                ec="none",
                alpha=0.4,
            )
        )
        ax.add_patch(
            patches.Rectangle(
                xy=(self.xleft, self.ybot - 0.002),
                width=self.xright - self.xleft,
                height=0.004,
                fc="grey",
                ec="none",
                alpha=0.4,
            )
        )

        # Add collboxes
        for collbox in self.collboxes:
            ax.add_patch(collbox)

        # Add gene arrows
        for genearrow, genename, start, end in self.genearrows_top:
            try:
                genename = self.gene_rename[genename]
            except KeyError:
                genename = genename
            ax.add_patch(genearrow)
            ax.annotate(
                genename,
                ((start + end) / 2, self.ytop - 0.02),
                xycoords="data",
                ha="center",
                va="center",
                fontsize=GENE_FS,
                bbox=dict(boxstyle="round", pad=0.2, fc="white", ec="none", alpha=0.6),
            )

        for genearrow, genename, start, end in self.genearrows_bot:
            try:
                genename = self.gene_rename[genename]
            except KeyError:
                genename = genename
            ax.add_patch(genearrow)
            ax.annotate(
                genename,
                ((start + end) / 2, self.ybot + 0.02),
                xycoords="data",
                ha="center",
                va="center",
                fontsize=GENE_FS,
                bbox=dict(boxstyle="round", pad=0.2, fc="white", ec="none", alpha=0.6),
            )

        # Add texts
        # Species name
        ax.text(
            x=0.09,
            y=self.ytop,
            s=self.species_name_top,
            size=SIDE_FS,
            ha="center",
            va="bottom",
            fontstyle="italic",
        )
        ax.text(
            x=0.09,
            y=self.ybot,
            s=self.species_name_bot,
            size=SIDE_FS,
            ha="center",
            va="bottom",
            fontstyle="italic",
        )

        # Sequence info
        if abs(self.x_top_endbase - self.x_top_startbase) < 300000:
            range1 = f"{self.x_top_startbase:,}~{self.x_top_endbase:,}bp"
        else:
            range1 = f"{self.x_top_startbase / 1000000:.2f}~{self.x_top_endbase / 1000000:.2f}Mbp"

        if abs(self.x_bot_endbase - self.x_bot_startbase) < 300000:
            range2 = f"{self.x_bot_startbase:,}~{self.x_bot_endbase:,}bp"
        else:
            range2 = f"{self.x_bot_startbase / 1000000:.2f}~{self.x_bot_endbase / 1000000:.2f}Mbp"

        ax.text(x=0.09, y=self.ytop, s=range1, size=SIDE_FS, ha="center", va="top")
        ax.text(x=0.09, y=self.ybot, s=range2, size=SIDE_FS, ha="center", va="top")

        # Close the axis
        ax.axis("off")

        fig.show()

        return fig

    def _range_init(self):
        # Define top and botton ranges and scale factors
        # Method 1: genome name list
        if self.draw_method == "gene_name":
            # Get gene table with those genes
            self.anchor = self.anchor[self.anchor["gene_top"].isin(self.list_param1)]
            self.anchor = self.anchor[self.anchor["gene_bot"].isin(self.list_param2)]

        # Method 2: genome index list
        if self.draw_method == "gene_index":
            # Get gene table in that range
            self.anchor = self.anchor[self.list_param1]
            self.anchor = self.anchor[self.list_param2]

        # Method 3: base range
        if self.draw_method == "base_range":
            # Get gene table in that range
            # list param structure: [chr_id, base_start, base_end]
            self.anchor = self.anchor[self.anchor["chr_top"] == self.list_param1[0]]
            self.anchor = self.anchor[self.anchor["chr_bot"] == self.list_param2[0]]

            self.x_top_startbase, self.x_top_endbase = self.list_param1[1:]
            self.x_bot_start, self.x_bot_endbase = self.list_param2[1:]

        else:
            # Set base range
            x1, x2 = self.anchor["x_top_start"].min(), self.anchor["x_top_end"].max()
            x3, x4 = self.anchor["x_bot_start"].min(), self.anchor["x_bot_end"].max()

            # Add some offset (2%)
            self.x_top_startbase = int(x1 - (x2 - x1 + 1) * 0.02)
            self.x_top_endbase = int(x2 + (x2 - x1 + 1) * 0.02)
            self.x_bot_startbase = int(x3 - (x4 - x3 + 1) * 0.02)
            self.x_bot_endbase = int(x4 + (x4 - x3 + 1) * 0.02)

        if self.reverse_top:
            self.x_top_startbase, self.x_top_endbase = (
                self.x_top_endbase,
                self.x_top_startbase,
            )
        if self.reverse_bot:
            self.x_bot_startbase, self.x_bot_endbase = (
                self.x_bot_endbase,
                self.x_bot_startbase,
            )

        # Set scale factors
        self.sf_top = (self.xright - self.xleft) / (
            abs(self.x_top_endbase - self.x_top_startbase) + 1
        )
        self.sf_bot = (self.xright - self.xleft) / (
            abs(self.x_bot_endbase - self.x_bot_startbase) + 1
        )

    def _arrows_collboxes_init(self):
        # Prepare gene arrows and collboxes
        self.genearrows_top, self.genearrows_bot, self.collboxes = [], [], []
        # gene_top, chr_top, gene_bot, chr_bot, x_top_start, x_top_end, strand_top, x_bot_start, x_bot_end, strand_bot
        for idx, info in self.anchor.iterrows():
            (
                gene_top,
                chr_top,
                gene_bot,
                chr_bot,
                x_top_start,
                x_top_end,
                strand_top,
                x_bot_start,
                x_bot_end,
                strand_bot,
            ) = info
            x_top_start, x_top_end = (
                (x_top_start, x_top_end)
                if strand_top == "+"
                else (x_top_end, x_top_start)
            )
            x_bot_start, x_bot_end = (
                (x_bot_start, x_bot_end)
                if strand_bot == "+"
                else (x_bot_end, x_bot_start)
            )

            if self.reverse_top:
                x1 = self.xright - (x_top_start - self.x_top_endbase) * self.sf_top
                x2 = self.xright - (x_top_end - self.x_top_endbase) * self.sf_top
            else:
                x1 = self.xleft + (x_top_start - self.x_top_startbase) * self.sf_top
                x2 = self.xleft + (x_top_end - self.x_top_startbase) * self.sf_top

            if self.reverse_bot:
                x3 = self.xright - (x_bot_start - self.x_bot_endbase) * self.sf_bot
                x4 = self.xright - (x_bot_end - self.x_bot_endbase) * self.sf_bot
            else:
                x3 = self.xleft + (x_bot_start - self.x_bot_startbase) * self.sf_bot
                x4 = self.xleft + (x_bot_end - self.x_bot_startbase) * self.sf_bot

            # For top
            self.genearrows_top.append(
                (
                    _gene_arrow(x1, x2, self.ytop, self.color_top),
                    gene_top,
                    x1,
                    x2,
                )
            )

            # For bottom
            self.genearrows_bot.append(
                (
                    _gene_arrow(x3, x4, self.ybot, self.color_bot),
                    gene_bot,
                    x3,
                    x4,
                )
            )

            # For collbox
            if self.collbox_style == "polygon":
                self.collboxes.append(
                    _get_polygon_ribbon(x1, x2, x3, x4, self.ytop, self.ybot, "grey")
                )

            if self.collbox_style == "bezier":
                self.collboxes.append(
                    _get_bezier_ribbon(x1, x2, x3, x4, self.ytop, self.ybot, "grey")
                )


def _draw_a_chromosome(
    ax: plt.Axes,
    chr_start: float | int,
    chr_end: float | int,
    chr_h: float | int,
    chr_name: str,
    color,
) -> plt.Axes:
    ax.add_patch(
        patches.Rectangle(
            xy=(chr_start, chr_h - 0.003),
            width=abs(chr_end - chr_start),
            height=0.006,
            fill=True,
            fc=color,
            ec="none",
        )
    )
    ax.add_patch(
        patches.Rectangle(
            xy=(chr_start, chr_h - 0.007),
            width=abs(chr_end - chr_start),
            height=0.014,
            fill=False,
            ec=(0.2, 0.2, 0.2, 0.2),
        )
    )
    ax.text(
        x=chr_start + (chr_end - chr_start) / 2,
        y=chr_h + 0.025,
        s=chr_name,
        ha="center",
        va="center",
        backgroundcolor="white",
    )
    return ax


def _get_polygon_ribbon(x1_start, x1_end, x2_start, x2_end, y1, y2, color_code):
    return patches.Polygon(
        [
            [x1_start, y1],
            [x1_end, y1],
            [x2_start, y2],
            [x2_end, y2],
        ],
        fc=color_code,
        ec="none",
        lw=0,
        alpha=0.5,
    )


def _get_bezier_ribbon(
    x_top_start, x_top_end, x_bot_start, x_bot_end, ytop, ybot, color_code
):
    mid_y = (ytop + ybot) / 2

    codes = [
        MplPath.MOVETO,
        MplPath.LINETO,
        MplPath.CURVE4,  # verts[2]
        MplPath.CURVE4,  # verts[3]
        MplPath.CURVE4,  # verts[4]
        MplPath.LINETO,
        MplPath.CURVE4,  # verts[6]
        MplPath.CURVE4,  # verts[7]
        MplPath.CURVE4,  # verts[8]
        MplPath.CLOSEPOLY,
    ]

    verts = [
        (x_top_start, ytop),  # MOVETO
        (x_top_end, ytop),  # LINETO
        (x_top_end, mid_y),  # CURVE4
        (x_bot_end, mid_y),  # CURVE4
        (x_bot_end, ybot),  # CURVE4
        (x_bot_start, ybot),  # LINETO
        (x_bot_start, mid_y),  # CURVE4
        (x_top_start, mid_y),  # CURVE4
        (x_top_start, ytop),  # CURVE4
        (x_top_start, ytop),  # CLOSEPOLY
    ]

    path = MplPath(verts, codes)
    patch = patches.PathPatch(path, fc=color_code, ec="none", lw=0, alpha=0.5)
    return patch


def _draw_collbox_polygon(
    ax: plt.Axes, x1_start, x1_end, x2_start, x2_end, y1, y2, color_code
) -> plt.Axes:
    ax.add_patch(
        patches.Polygon(
            [
                [x1_start, y1],
                [x1_end, y1],
                [x2_start, y2],
                [x2_end, y2],
            ],
            fc=color_code,
            ec="none",
            lw=0,
            alpha=0.5,
        )
    )
    return ax


def _draw_collbox_bezier(
    ax: plt.Axes,
    x1_start,
    x1_end,
    x2_start,
    x2_end,
    y1,
    y2,
    color_code,
) -> plt.Axes:
    mid_y = (y1 + y2) / 2

    codes = [
        MplPath.MOVETO,
        MplPath.LINETO,
        MplPath.CURVE4,  # verts[2]
        MplPath.CURVE4,  # verts[3]
        MplPath.CURVE4,  # verts[4]
        MplPath.LINETO,
        MplPath.CURVE4,  # verts[6]
        MplPath.CURVE4,  # verts[7]
        MplPath.CURVE4,  # verts[8]
        MplPath.CLOSEPOLY,
    ]

    verts = [
        (x1_start, y1),  # MOVETO
        (x1_end, y1),  # LINETO
        (x1_end, mid_y),  # CURVE4
        (x2_end, mid_y),  # CURVE4
        (x2_end, y2),  # CURVE4
        (x2_start, y2),  # LINETO
        (x2_start, mid_y),  # CURVE4
        (x1_start, mid_y),  # CURVE4
        (x1_start, y1),  # CURVE4
        (x1_start, y1),  # CLOSEPOLY
    ]

    path = MplPath(verts, codes)
    ax.add_patch(patches.PathPatch(path, fc=color_code, ec="none", lw=0, alpha=0.5))
    return ax


def _gene_arrow(start: float, end: float, y: float, color):
    if abs(start - end) >= 0.02:
        if start < end:
            return patches.Polygon(
                xy=(
                    [start, y - 0.01],
                    [start, y + 0.01],
                    [end - 0.02, y + 0.01],
                    [end, y],
                    [end - 0.02, y - 0.01],
                    [start, y - 0.01],
                ),
                color=color,
            )
        else:
            return patches.Polygon(
                xy=(
                    [start, y - 0.01],
                    [start, y + 0.01],
                    [end + 0.02, y + 0.01],
                    [end, y],
                    [end + 0.02, y - 0.01],
                    [start, y - 0.01],
                ),
                color=color,
            )
    else:
        return patches.Polygon(
            xy=([start, y - 0.01], [start, y + 0.01], [end, y], [start, y - 0.01]),
            color=color,
        )
