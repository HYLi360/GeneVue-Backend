import gc
from dataclasses import dataclass, fields
from multiprocessing import Pool
from typing import Literal, Tuple

import numpy as np
import pandas as pd
from rich.progress import Progress
from rich.table import Table

from genevue import console
from genevue.Coll._corecoll import DoCollinearity, CollRes


@dataclass(
    kw_only=True,
)
class Collinearity:
    # simp_gff:
    # chr, gene_id, start, end, strand, order
    bed1: pd.DataFrame
    bed2: pd.DataFrame
    blast: pd.DataFrame
    savefile_name: str

    # Initialize parameters with default values
    repeat_number = 10
    multiple: int = 1
    over_gap: int = 15
    comparison: Literal["genomes", "chromosomes"] = "genomes"
    position: Literal["order", "end"] = "order"
    grading: Tuple[int, int, int] = (50, 40, 25)
    mg: Tuple[int, int] = (40, 40)
    coverage_ratio: float = 0.9
    processes_num: int = 4
    gap_penalty: int = -1
    over_length: int = 0
    pvalue_min: float = 1

    options = {
        "gap_penalty": gap_penalty,
        "over_length": over_length,
        "pvalue_min": pvalue_min,
        "over_gap": over_gap,
        "coverage_ratio": coverage_ratio,
        "grading": grading,
        "mg": mg,
    }

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
                "options",
                "over_length",
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

        newblast = (
            blast.groupby(["chr1", "chr2"]).apply(assign_grading).reset_index(drop=True)
        )
        newblast["grading"] = newblast["grading"].astype(int)
        return newblast[newblast["grading"] > 0]

    def _deal_blast_for_genomes(
        self,
        blast: pd.DataFrame,
        rednum: int,
        repeat_number: int,
    ):
        # Initialize the grading column
        blast["grading"] = 0

        # Define the blue number as the sum of rednum and the predefined constant
        bluenum = 4 + rednum

        # Get the indices for each group by sorting the 11th column (bitscore) in descending order
        index = [
            group.sort_values(by=["bitscore"], ascending=[False])[
                :repeat_number
            ].index.tolist()
            for name, group in blast.groupby(["gene1"])
        ]

        # Split the indices into red, blue, and gray groups
        reddata = np.array([k[:rednum] for k in index], dtype=object)
        bluedata = np.array([k[rednum:bluenum] for k in index], dtype=object)
        graydata = np.array([k[bluenum:repeat_number] for k in index], dtype=object)

        # Concatenate the results into flat lists
        redindex = np.concatenate(reddata) if reddata.size else []
        blueindex = np.concatenate(bluedata) if bluedata.size else []
        grayindex = np.concatenate(graydata) if graydata.size else []

        # Update the grading column based on the group indices
        blast.loc[redindex, "grading"] = self.grading[0]
        blast.loc[blueindex, "grading"] = self.grading[1]
        blast.loc[grayindex, "grading"] = self.grading[2]

        # Return only the rows with non-zero grading
        return blast[blast["grading"] > 0]

    def run(self):
        # Read simplified gff data
        gff1 = self.bed1
        gff1["strand"] = self.bed1["strand"].map({"+": 1, "-": -1})

        gff2 = self.bed2
        gff2["strand"] = self.bed2["strand"].map({"+": 1, "-": -1})

        # Map positions and chromosome information
        # (Using pd.merge instead of pd.map to improve perfermance)
        gff1_cols = gff1[["gene_id", self.position, "chr"]].rename(
            columns={"gene_id": "gene1", self.position: "loc1", "chr": "chr1"}
        )
        gff2_cols = gff2[["gene_id", self.position, "chr"]].rename(
            columns={"gene_id": "gene2", self.position: "loc2", "chr": "chr2"}
        )
        self.blast = self.blast.merge(gff1_cols, on="gene1", how="left")
        self.blast = self.blast.merge(gff2_cols, on="gene2", how="left")

        # Apply blast filtering and grading
        if self.comparison.lower() == "genomes":
            blast = self._deal_blast_for_genomes(
                self.blast, int(self.multiple), int(self.repeat_number)
            )
        if self.comparison.lower() == "chromosomes":
            blast = self._deal_blast_for_chromosomes(
                self.blast, int(self.multiple), int(self.repeat_number)
            )

        if len(self.blast) == 0:
            raise RuntimeError("SimpleGFF3 and BLAST result do not seem to match.")

        console.log(f"The filtered homologous gene pairs are {len(self.blast)}.")

        # Group blast data by 'chr1' and 'chr2'
        total = [
            ((chr1, chr2), group)
            for (chr1, chr2), group in self.blast.groupby(["chr1", "chr2"])
        ]
        # Determine chunk size for multiprocessing
        chunks_size = int(np.ceil(len(total) / float(self.processes_num)))

        # Running with multi-processes
        data = []

        with Progress() as progress:
            task = progress.add_task("CollinearScan running.", total=len(total))
            with Pool(processes=self.processes_num) as pool:
                for result in pool.imap(self._single, total, chunksize=chunks_size):
                    data.append(result)
                    progress.update(task, advance=1)

        self._format_output(data)

    def _single(self, args):
        chr_pair, points = args
        collinearity = DoCollinearity(
            points=points,
            chr_tuple=chr_pair,
            **self.options,
        )
        return collinearity.run()

    def _format_output(self, total_res_ls: list[CollRes]) -> None:
        # Target format:
        # Alignment 1: score=100 pvalue=0.03 N=3 1&1 minus
        # s1g1 1 s2g1 1 1
        # s1g2 2 s1g2 2 -1
        # s1g3 3 s1g3 3 1

        # simp_gff:
        # chr, gene_id, start, end, strand, order

        # block_counter
        counter = 1

        # Start!
        fout = open(self.savefile_name + ".coll", "w+")
        for chr_pair in total_res_ls:
            for i in range(chr_pair.len):
                # Get block.
                bk, p_value, score = chr_pair.get_block(i)

                # Get block length.
                n = bk.shape[0]

                # If this block not long enough?
                if n < self.over_gap:
                    continue

                # Give gene_id and strand based on two gff.
                chr1, chr2 = chr_pair.chr1, chr_pair.chr2
                bk = pd.merge(
                    left=bk,
                    right=self.bed1[self.bed1["chr"] == chr1]
                    .loc[:, ["gene_id", "strand", "order"]]
                    .rename(
                        columns={
                            "gene_id": "gene1",
                            "strand": "strand1",
                            "order": "loc1",
                        }
                    ),
                    how="inner",
                    on="loc1",
                )
                bk = pd.merge(
                    left=bk,
                    right=self.bed2[self.bed2["chr"] == chr2]
                    .loc[:, ["gene_id", "strand", "order"]]
                    .rename(
                        columns={
                            "gene_id": "gene2",
                            "strand": "strand2",
                            "order": "loc2",
                        }
                    ),
                    how="inner",
                    on="loc2",
                )

                # Polishing bk.
                # before: ['loc1', 'loc2', 'gene1', 'strand1', 'gene2', 'strand2']
                # after:  ['gene1', 'loc1', 'strand1', 'gene2', 'loc2', 'strand2']
                bk = bk.loc[:, ["gene1", "loc1", "strand1", "gene2", "loc2", "strand2"]]

                # Calculate the direction ('plus' or 'minus'.)
                bk["direction"] = (
                    (bk["strand1"] + bk["strand2"]).abs().map({0: -1, 2: 1})
                )
                mark = "plus" if bk["direction"].sum() > 0 else "minus"

                # Write bkinfo.
                fout.write(
                    f"# Alignment {counter}: score={score} pvalue={p_value:.4f} N={n} {chr1}&{chr2} {mark}\n"
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

        # Finish.
        fout.close()

        # Go Yours Garbage!
        gc.collect()
