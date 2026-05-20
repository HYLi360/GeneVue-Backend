"""
Plot BLAST dotplot, based on paired genes, bed, and others (if you have):
- Genome sequences, for confirm the length of chromosomes.
- Genes masker, select the genes in each genome.
- An option, how the points locale themselves (by order, or middle-point of gene. Order by default)
"""

from pathlib import Path
from typing import Literal, Optional

from matplotlib import pyplot as plt

from genevue.sequences import FASTA
from genevue.utils.parse import bed_reader


class CollScatPlot:
    def __init__(
        self,
        coll_path: Path,
        bed1_path: Path,
        bed2_path: Path,
        out_plot_path: Path,
        genome_seq1_path: Optional[Path] = None,
        genome_seq2_path: Optional[Path] = None,
        chromosome_order1_path: Optional[Path] = None,
        chromosome_order2_path: Optional[Path] = None,
        gene_mask1_path: Optional[Path] = None,
        gene_mask2_path: Optional[Path] = None,
        genome_name1: str = "",
        genome_name2: str = "",
        locale_method: Literal["order", "mid-point"] = "order",
    ):
        self.coll_path = coll_path
        self.bed1_path = bed1_path
        self.bed2_path = bed2_path
        self.out_plot_path = out_plot_path
        self.genome_seq1_path = genome_seq1_path
        self.genome_seq2_path = genome_seq2_path
        self.chromosome_order_path = chromosome_order1_path
        self.chromosome_order_path = chromosome_order2_path
        self.gene_mask_path = gene_mask1_path
        self.gene_mask_path = gene_mask2_path
        self.genome_name1 = genome_name1
        self.genome_name2 = genome_name2
        self.locale_method = locale_method

        self.bed1 = bed_reader(self.bed1_path)
        self.bed2 = bed_reader(self.bed2_path)

        self._pre_process()

    def _pre_process(self):
        # pre-process
        # get length of chromosome, and position for every genes
        if self.locale_method == "mid-point":
            # use the absolute position
            if self.genome_seq1_path:
                fasta = FASTA(self.genome_seq1_path)
                self.chromosome_length1 = fasta.seq_lengths
            else:
                self.chromosome_length1 = (
                    self.bed1.groupby("chrom")
                    .agg({"chromEnd": "max"})["chromEnd"]
                    .to_dict()
                )
            if self.genome_seq2_path:
                fasta = FASTA(self.genome_seq2_path)
                self.chromosome_length2 = fasta.seq_lengths
            else:
                self.chromosome_length2 = (
                    self.bed2.groupby("chrom")
                    .agg({"chromEnd": "max"})["chromEnd"]
                    .to_dict()
                )

            # position
            self.bed1[]



    def plot(self):
        pass
