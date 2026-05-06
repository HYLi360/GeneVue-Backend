# coding: utf-8
"""
GeneSetAnalysis.py

Module called for proteins mode.

Author(s): Matthew Berkeley, Mathieu Seppey, Mose Manni, Felipe Simao, Rob Waterhouse

Copyright (c) 2015-2025, Evgeny Zdobnov (ez@ezlab.org). All rights reserved.

License: Licensed under the MIT license. See LICENSE.md file.

"""

from genevue.QC.BUSCO.BuscoAnalysis import BuscoAnalysis
from genevue.QC.BUSCO.BuscoLogger import BuscoLogger
from genevue.QC.BUSCO.Analysis import ProteinAnalysis
from Bio import SeqIO

logger = BuscoLogger.get_logger(__name__)


class GeneSetAnalysis(ProteinAnalysis, BuscoAnalysis):
    """
    This class runs a BUSCO analysis on a gene set.
    """

    _mode = "proteins"

    def __init__(self):
        """
        Initialize an instance.
        """
        super().__init__()

        if self.input_file.endswith(".gz"):
            import gzip

            handle = gzip.open(self.input_file, "rt")
        else:
            handle = open(self.input_file, "rt")

        self.gene_details = {
            record.id: {"aa_seq": record}
            for record in list(SeqIO.parse(handle, "fasta"))
        }

    def cleanup(self):
        super().cleanup()

    def run_analysis(self):
        """
        This function calls all needed steps for running the analysis.
        """
        super().run_analysis()
        self.run_hmmer(self.input_file)
        self.hmmer_runner.write_buscos_to_file()
        return

    def reset(self):
        super().reset()
