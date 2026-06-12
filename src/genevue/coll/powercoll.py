#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Literal, Tuple, Hashable, List, cast, Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from rich.table import Table

from genevue import console, setup_rich_logger

logger = setup_rich_logger(__name__, console)


class CollRes:
    def __init__(self, chr1: Any, chr2: Any):
        self.chr1, self.chr2 = str(chr1), str(chr2)
        self.loc1ls: list[NDArray[np.int32]] = []
        self.loc2ls: list[NDArray[np.int32]] = []
        self.pvaluels: list[(float | int)] = []
        self.scorels: list[(float | int)] = []
        self.directls: list[str] = []
        self.dencels1: list[(float | int)] = []
        self.dencels2: list[(float | int)] = []

    def add_coll(
        self,
        loc1: NDArray[np.int32],
        loc2: NDArray[np.int32],
        pvalue: float | int,
        score: float | int,
        direct: Literal["plus", "minus"],
        dence1: float | int,
        dence2: float | int,
    ):
        self.loc1ls.append(loc1)
        self.loc2ls.append(loc2)
        self.pvaluels.append(pvalue)
        self.scorels.append(score)
        self.directls.append(direct)
        self.dencels1.append(dence1)
        self.dencels2.append(dence2)

    def get_coll(
        self, idx: int
    ) -> tuple[NDArray[np.int32], NDArray[np.int32], float, float, str, float, float]:
        return (
            self.loc1ls[idx],
            self.loc2ls[idx],
            self.pvaluels[idx],
            self.scorels[idx],
            self.directls[idx],
            self.dencels1[idx],
            self.dencels2[idx],
        )

    def __len__(self):
        return len(self.pvaluels)


@dataclass(
    kw_only=True,
)
class Collinearity:
    # Parameters needed
    # bed:
    # chr, gene_id, start, end, strand, order
    bed1: pd.DataFrame
    bed2: pd.DataFrame
    blast: pd.DataFrame

    # Initialize parameters with default values
    multiple: int = 1
    repeat_number: int = 10
    over_gap: int = 15
    comparison: Literal["genomes", "chromosomes"] = "genomes"
    position: Literal["order", "end"] = "order"
    grading: Tuple[int, int, int] = (50, 40, 25)
    mg: Tuple[int, int] = (40, 40)
    gap_penalty: float | int = -1
    pvalue_min: float | int = 1

    loc1: NDArray = np.array([])
    loc2: NDArray = np.array([])
    gradings: NDArray = np.array([])
    resls: List[CollRes] = field(default_factory=list)

    def show_all(self):
        table = Table(
            caption="Checkout those parameters!",
            caption_justify="left",
            highlight=True,
        )
        table.add_column("Param. Name", justify="left", style="cyan", no_wrap=True)
        table.add_column("Value", justify="right", style="green", no_wrap=True)

        for f in fields(self):
            if f.name not in [
                "bed1",
                "bed2",
                "blast",
                "savefile_name",
                "loc1",
                "loc2",
                "gradings",
                "resls",
            ]:
                table.add_row(f.name, str(getattr(self, f.name)))
        console.print(table)

    def _deal_blast_for_chromosomes(
        self,
        blast: pd.DataFrame,
        rednum: int,
        repeat_number: int,
    ):
        bluenum = rednum
        blast = blast.sort_values(by=["gene1", "bitscore"], ascending=[True, False])

        def assign_grading(group):
            group["cumcount"] = group.groupby(1).cumcount()
            group = group[group["cumcount"] <= repeat_number]
            group["grading"] = pd.cut(
                group["cumcount"],
                bins=[-1, 0, bluenum, repeat_number],
                labels=self.grading,
                right=True,
            )
            return group

        newblast: pd.DataFrame = cast(
            pd.DataFrame,
            cast(
                object,
                blast.groupby(["chr1", "chr2"])
                .apply(assign_grading)
                .reset_index(drop=True),
            ),
        )
        newblast["grading"] = newblast["grading"].astype(int)
        return newblast[newblast["grading"] > 0]

    def _deal_blast_for_genomes(
        self,
        blast: pd.DataFrame,
        rednum: int,
        repeat_number: int,
    ):
        # Define the blue number as the sum of rednum and the predefined constant
        bluenum = 4 + rednum

        # We do these steps to give grades by every pairs.
        # 1. Rip gene1 and bitscore columns and sort them.
        #    ! We keep the index (Not using reset_index()) to ensure we can restore the order
        rank = blast.loc[:, ["gene1", "bitscore"]].sort_values(
            ["gene1", "bitscore"], ascending=[True, False]
        )

        # 2. Using groupby and cumcount.
        rank = rank.groupby("gene1").cumcount()

        # 3. Mapping new values.
        grading = pd.Series(0, index=rank.index)
        grading[rank < rednum] = self.grading[0]
        grading[(rank >= rednum) & (rank < bluenum)] = self.grading[1]
        grading[(rank >= bluenum) & (rank < repeat_number)] = self.grading[2]

        # 4. And finally, combine the original blast-df.
        blast["grading"] = grading

        # Return only the rows with non-zero grading
        return blast[blast["grading"] > 0]

    def run(self):
        # Read simplified gff data
        self.bed1["strand"] = self.bed1["strand"].map({"+": 1, "-": -1})
        self.bed2["strand"] = self.bed2["strand"].map({"+": 1, "-": -1})

        # Map positions and chromosome information
        # (Using pd.merge instead of pd.map to improve perfermance)
        self.bed1 = self.bed1.rename(
            columns={
                "chr": "chr1",
                "gene_id": "gene1",
                self.position: "loc1",
                "strand": "strand1",
            }
        )
        self.bed2 = self.bed2.rename(
            columns={
                "chr": "chr2",
                "gene_id": "gene2",
                self.position: "loc2",
                "strand": "strand2",
            }
        )
        self.blast = self.blast.merge(self.bed1, on="gene1", how="left")
        self.blast = self.blast.merge(self.bed2, on="gene2", how="left")

        # Apply blast filtering and grading
        if self.comparison.lower() == "genomes":
            self.blast = self._deal_blast_for_genomes(
                self.blast, self.multiple, self.repeat_number
            )
        if self.comparison.lower() == "chromosomes":
            self.blast = self._deal_blast_for_chromosomes(
                self.blast, self.multiple, self.repeat_number
            )

        if len(self.blast) == 0:
            raise RuntimeError("SimpleGFF3 and BLAST result do not seem to match.")

        logger.info(f"The filtered homologous gene pairs are {len(self.blast)}.")

        self.blast = self.blast.sort_values(
            ["chr1", "chr2", "loc1", "loc2"]
        ).reset_index(drop=True)

        self.loc1 = self.blast["loc1"].to_numpy(dtype=np.int32, copy=False)
        self.loc2 = self.blast["loc2"].to_numpy(dtype=np.int32, copy=False)
        self.gradings = self.blast["grading"].to_numpy(dtype=np.float32, copy=False)

        self.resls = []
        for (chr1, chr2), group in self.blast.groupby(["chr1", "chr2"]):
            self.resls.append(
                self._process(
                    chr_tuple=(chr1, chr2), index=group.index.to_numpy(dtype=np.uint32)
                )
            )

    def _process(
        self,
        chr_tuple: tuple[Hashable, Hashable],
        index: NDArray[np.uint32],
    ) -> CollRes:
        # Forward and Backward scaning.
        score1, usedtimes1, parent1 = self._score_matrix(
            self.loc1[index], self.loc2[index], self.gradings[index], "forward"
        )
        score2, usedtimes2, parent2 = self._score_matrix(
            self.loc1[index], self.loc2[index], self.gradings[index], "backward"
        )

        # Collect result.
        res = self._max_path(
            self.loc1[index],
            self.loc2[index],
            score1,
            score2,
            usedtimes1,
            usedtimes2,
            parent1,
            parent2,
            chr_tuple[0],
            chr_tuple[1],
        )
        return res

    def _score_matrix(
        self,
        loc1: NDArray[np.int32],
        loc2: NDArray[np.int32],
        grading: NDArray[np.float32],
        direction: Literal["forward", "backward"],
    ) -> tuple[NDArray[np.float32], NDArray[np.int32], NDArray[np.int32]]:
        """
        Build sequence of parent.

        If you want use it directly, please note that make sure your points (DataFrame) has been sorted by

        `points.sort_values(by=['loc1', 'loc2'], kind='mergesort')`,

        whatever you select `forward` or `backward`.
        """
        scorels: NDArray[np.float32] = grading.copy()
        n: int = len(loc1)
        used: NDArray[np.int32] = cast(NDArray[np.int32], np.zeros(n, dtype=np.int32))
        parent: NDArray[np.int32] = np.full(n, -1, dtype=np.int32)

        # i = start point.
        if direction == "forward":
            for i in range(n):
                row, col = loc1[i], loc2[i]
                left = np.searchsorted(loc1, row + 1)
                right = np.searchsorted(loc1, row + self.mg[0])
                row_i_old, gap = row, self.mg[1]

                for j in range(left, right):
                    # add constraint of loc2 window range
                    if loc2[j] <= col or loc2[j] >= col + self.mg[1]:
                        continue

                    if (loc2[j] - col) > gap and loc1[j] > row_i_old:
                        break

                    score = (
                        grading[j] + (loc1[j] - row + loc2[j] - col) * self.gap_penalty
                    )
                    if score <= 0:
                        continue

                    cand = scorels[i] + score
                    if scorels[j] < cand:
                        scorels[j] = cand
                        parent[j] = i
                        used[i] += 1
                        used[j] += 1
                        dcol = loc2[j] - col
                        gap = min(gap, dcol)
                        row_i_old = loc1[j]
        else:
            for i in range(n - 1, -1, -1):
                row, col = loc1[i], loc2[i]
                left = np.searchsorted(loc1, (row - self.mg[0] + 1))
                right = np.searchsorted(loc1, row)
                row_i_old, gap = row, self.mg[1]

                for j in range(left, right):
                    if loc2[j] <= col or loc2[j] >= col + self.mg[1]:
                        continue

                    if (loc2[j] - col) > gap and loc1[j] < row_i_old:
                        break

                    score = (
                        grading[j] + (row - loc1[j] + loc2[j] - col) * self.gap_penalty
                    )
                    if score <= 0:
                        continue

                    cand = scorels[i] + score
                    if scorels[j] < cand:
                        scorels[j] = cand
                        parent[j] = i
                        used[i] += 1
                        used[j] += 1
                        dcol = loc2[j] - col
                        gap = min(gap, dcol)
                        row_i_old = loc1[j]
        return scorels, used, parent

    @staticmethod
    def _backtrack(parent: NDArray[np.int32], end: int) -> NDArray[np.int32]:
        """Return path positions from start->end ( [start, end] )."""
        path = []
        cur = end
        while cur != -1:
            path.append(cur)
            cur = int(parent[cur])
        path.reverse()
        return np.asarray(path, dtype=np.int32)

    def _max_path(
        self,
        loc1: NDArray[np.int32],
        loc2: NDArray[np.int32],
        score1: NDArray[np.float32],
        score2: NDArray[np.float32],
        usedtimes1: NDArray[np.int32],
        usedtimes2: NDArray[np.int32],
        parent1: NDArray[np.int32],
        parent2: NDArray[np.int32],
        chr1: Hashable,
        chr2: Hashable,
    ):
        # Here I still adapted a greedy algorithm, whose implementation principle is as follows:
        # - Among all unused points, locate the one with the highest score
        # - Trace back from this point to obtain a pair of collinear chains
        # - Check the chain's length against p_value; if either condition is violated, discard it and continue loop
        # - Reset the corresponding `times` flag to zero and proceed to the next search
        # - Repeat the above process until no further chains can be extracted
        n = len(loc1)
        coverage = 0
        res = CollRes(chr1, chr2)

        # use this to exclude points picked
        usedtimes1[usedtimes1 > 0], usedtimes2[usedtimes2 > 0] = 1, 1

        # use this to calculate p_value
        times = np.ones(n, dtype=np.int32)

        while True:
            cand1, cand2 = np.flatnonzero(usedtimes1 > 0), np.flatnonzero(
                usedtimes2 > 0
            )
            if (
                (len(cand1) < self.over_gap)
                and (len(cand2) < self.over_gap)
                or (len(score1[cand1]) == 0)
                or (len(score2[cand2]) == 0)
            ):
                break

            # Pick the point with the highest score
            if score1[cand1].max() > score2[cand2].max():
                # Forward has the highest score
                # forward -> left first (original argmax)
                score = score1[cand1].max()
                end = int(cand1[np.argmax(score1[cand1])])

                # Backtrack all points.
                path = self._backtrack(parent1, end)

                # Filt unused points
                path = path[usedtimes1[path] == 1]

                # If path's too short?
                if len(path) < self.over_gap:
                    usedtimes1[end] = 0
                    continue

                # bounding-box
                l1, l2, n1 = (
                    loc1[path],
                    loc2[path],
                    times[path.min() : path.max() + 1].sum(),
                )
                l1max, l1min, l2max, l2min = (
                    int(l1.max()),
                    int(l1.min()),
                    int(l2.max()),
                    int(l2.min()),
                )

                # set times
                times[path.min() : path.max() + 1] += 1

                # add coverage
                coverage += (path.max() - path.min() + 1) / n

                # calculate p-value
                pvalue = self._p_value_calc(
                    m=len(path),
                    n=path.max() - path.min() + 1,
                    n1=n1,
                    l1=l1max - l1min + 1,
                    l2=l2max - l2min + 1,
                    score=score,
                )

                # append or drop
                if pvalue < self.pvalue_min:
                    res.add_coll(
                        l1,
                        l2,
                        pvalue,
                        score,
                        "plus",
                        len(l1) / (l1max - l1min + 1),
                        len(l2) / (l2max - l2min + 1),
                    )

                # remove used points (and points among that)
                (
                    usedtimes1[path.min() : path.max() + 1],
                    usedtimes2[path.min() : path.max() + 1],
                ) = (0, 0)

            else:
                score = score2[cand2].max()
                # The issue with argmax affinity.
                # By default, argmax selects the leftmost element.
                # However, in backward branch, we want it to select the rightmost element.
                end = int(cand2[::-1][np.argmax(score2[cand2][::-1])])

                # Backtrack all points.
                path = self._backtrack(parent2, end)

                # Filt unused points
                path = path[usedtimes2[path] == 1]

                # If path's too short?
                if len(path) < self.over_gap:
                    usedtimes2[end] = 0
                    continue

                # bounding-box
                l1, l2, n1 = (
                    loc1[path],
                    loc2[path],
                    times[path.min() : path.max() + 1].sum(),
                )
                l1max, l1min, l2max, l2min = (
                    int(l1.max()),
                    int(l1.min()),
                    int(l2.max()),
                    int(l2.min()),
                )

                # set times
                times[path.min() : path.max() + 1] += 1

                # add coverage
                coverage += (path.max() - path.min() + 1) / n

                # calculate p-value
                pvalue = self._p_value_calc(
                    m=len(path),
                    n=path.max() - path.min() + 1,
                    n1=n1,
                    l1=l1max - l1min + 1,
                    l2=l2max - l2min + 1,
                    score=score,
                )

                # append or drop
                if pvalue < self.pvalue_min:
                    res.add_coll(
                        l1,
                        l2,
                        pvalue,
                        score,
                        "minus",
                        (len(l1)) / (l1max - l1min + 1),
                        (len(l2)) / (l2max - l2min + 1),
                    )

                # remove used points (and points among that)
                (
                    usedtimes1[path.min() : path.max() + 1],
                    usedtimes2[path.min() : path.max() + 1],
                ) = (0, 0)
        return res

    def _p_value_calc(self, m, n, n1, l1, l2, score) -> float:
        return (
            (1 - score / m / self.grading[0])
            * (n1 - m + 1)
            / n
            * (l1 - m + 1)
            * (l2 - m + 1)
            / l1
            / l2
        )

    def export(self, collfile_path: Path, anchor_path: Path) -> None:
        # Target format:
        # Alignment 1: score=100 pvalue=0.03 N=3 1&1 minus
        # s1g1 1 s2g1 1 1
        # s1g2 2 s2g2 2 -1
        # s1g3 3 s2g3 3 1

        # simp_gff:
        # chr, gene_id, start, end, strand, order
        self.bed1 = self.bed1.loc[:, ["loc1", "gene1", "strand1"]]
        self.bed2 = self.bed2.loc[:, ["loc2", "gene2", "strand2"]]

        # block_counter
        counter = 1

        # Start!
        fout = open(collfile_path, "w")
        anchor = open(anchor_path, "w")
        for collres in self.resls:
            for i in range(len(collres)):
                # Get block.
                loc1, loc2, pvalue, score, direct, dence1, dence2 = collres.get_coll(i)
                loc1, loc2 = loc1[::-1], loc2[::-1]

                # if sequence completely same?
                if not (loc1 != loc2).all():
                    continue

                bk = pd.DataFrame(
                    data={
                        "loc1": loc1,
                        "loc2": loc2,
                    }
                )

                # Get block length.
                n = len(loc1)

                # Give gene_id and strand based on two gff.
                chr1, chr2 = collres.chr1, collres.chr2
                bk = pd.merge(
                    left=bk,
                    right=self.bed1,
                    how="inner",
                    on="loc1",
                )
                bk = pd.merge(
                    left=bk,
                    right=self.bed2,
                    how="inner",
                    on="loc2",
                )
                # Polishing bk.
                # before: ['loc1', 'loc2', 'gene1', 'strand1', 'gene2', 'strand2']
                # after:  ['gene1', 'loc1', 'strand1', 'gene2', 'loc2', 'strand2']
                bk = bk.loc[:, ["gene1", "loc1", "strand1", "gene2", "loc2", "strand2"]]

                bk["direction"] = (bk["strand1"] + bk["strand2"]).abs() - 1

                # Write bkinfo.
                fout.write(
                    f"# Alignment {counter}: score={score} pvalue={pvalue:.4f} N={n} {chr1}&{chr2} {direct}\n"
                )

                # Write anchor info.
                anchor.write(
                    f"{chr1}\t{chr2}\t{direct}\t{score}\t{pvalue:4f}\t{dence1:.2f}\t{dence2:.2f}"
                    f"\t{n}\t{loc1[0]}\t{loc1[-1]}\t{loc2[0]}\t{loc2[-1]}\n"
                )

                # Write bk.
                fout.write(
                    bk.loc[:, ["gene1", "loc1", "gene2", "loc2", "direction"]].to_csv(
                        sep="\t",
                        index=False,
                        header=False,
                    )
                )

                counter += 1
                continue
        fout.close()
        anchor.close()
