#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.text import Text

import pandas as pd

from typing import Optional, Literal


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

        # todo: Multiple chromosome support?

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


def _get_polygon_ribbon(
    x_top_start, x_top_end, x_bot_start, x_bot_end, ytop, ybot, color_code
):
    return patches.Polygon(
        [
            [x_top_start, ytop],
            [x_top_end, ytop],
            [x_bot_start, ybot],
            [x_bot_end, ybot],
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
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE4,  # verts[2]
        Path.CURVE4,  # verts[3]
        Path.CURVE4,  # verts[4]
        Path.LINETO,
        Path.CURVE4,  # verts[6]
        Path.CURVE4,  # verts[7]
        Path.CURVE4,  # verts[8]
        Path.CLOSEPOLY,
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

    path = Path(verts, codes)
    patch = patches.PathPatch(path, fc=color_code, ec="none", lw=0, alpha=0.5)
    return patch


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
