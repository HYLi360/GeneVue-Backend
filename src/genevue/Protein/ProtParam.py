#  Copyright (c) 2026 HYLi360. All rights reserved.
#
#  see license in LICENSE
#  see side-package licenses in LICENSE_OF_SIDE_PACKAGES

"""
A protein phys/chem-properties predict tool, based on ExPASy ProtParam and Biopython.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from genevue import console, setup_rich_logger
from genevue.Constants.AminoAcids import molecular_formula_protein
from genevue.Sequences import remove_any_spacechar

logger = setup_rich_logger(__name__, console)


class ProtParam:
    def __init__(self, seq: str | Seq | SeqRecord):
        if isinstance(seq, SeqRecord):
            self.seq = str(seq)
            self.seqid = seq.id
            self.description = seq.description
        else:
            self.seq = remove_any_spacechar(str(seq))
            self.seqid = "input_sequence"
            self.description = ""

        # Analyse.
        self.pa_res = ProteinAnalysis(self.seq)

    # Length
    @property
    def length(self):
        return len(self.seq)

    # Theoretical pI and Charge at pH=7
    @property
    def isoelectric_point(self):
        return self.pa_res.isoelectric_point()

    @property
    def charge_at_7(self):
        return self.pa_res.charge_at_pH(7)

    # Molecular Weight (in Da)
    @property
    def molecular_weight(self):
        return molecular_weight(self.seq, "protein")

    # Amino Acid Composition
    @property
    def composition_aa(self):
        return self.pa_res.count_amino_acids()

    @property
    def composition_aa_percent(self):
        return self.pa_res.amino_acids_percent

    # Negatively/Positively Charged Amino Acid
    @property
    def negative(self):
        neg_counter = 0
        for a in self.seq:
            if a.upper() == "D" or a.upper() == "E":
                neg_counter += 1
        return neg_counter

    @property
    def positive(self):
        pos_counter = 0
        for a in self.seq:
            if a.upper() == "R" or a.upper() == "K":
                pos_counter += 1
        return pos_counter

    # Atomic Composition
    @property
    def formula(self):
        return molecular_formula_protein(self.seq)[0]

    # Atoms Count
    @property
    def atoms_count(self):
        return molecular_formula_protein(self.seq)[1]

    # Extinction coefficients
    @property
    def extinction_coefficients_mec_reduced(self):
        return self.pa_res.molar_extinction_coefficient()[0]

    @property
    def extinction_coefficients_mec_cystines(self):
        return self.pa_res.molar_extinction_coefficient()[1]

    # Instability Index (II)
    @property
    def instability_index(self):
        return self.pa_res.instability_index()

    # check is it unstable based on ii
    @property
    def is_stable_based_on_ii(self):
        return self.instability_index <= 40

    # GRAnd AVerage of hYdropathicity (GRAVY) in KyteDoolitle scale
    @property
    def gravy_kd(self):
        return self.pa_res.gravy()

    # Aliphatic Index

    # Aromaticity
    @property
    def aromaticity(self):
        return self.pa_res.aromaticity()

    # Export the report
    @property
    def report(self) -> str:
        return f"""
Sequence ID:           {self.seqid}
Sequence description:  {self.description}

Number of amino acids: {self.length}
Theoretical pI:        {self.isoelectric_point:.2f}
Molecular weight:      {self.molecular_weight:.2f}

Amino acid composition:{self.composition_aa}

Total number of negatively charged residues (Asp + Glu): {self.negative}
Total number of positively charged residues (Arg + Lys): {self.positive}

Atomic composition:     {self.formula}
Total number of atoms:  {self.atoms_count}

Extinction coefficients:
{self.extinction_coefficients_mec_cystines} (when all Cys become Cys-Cys or Cystines)
{self.extinction_coefficients_mec_reduced} (when all Cys stay reduced)

Instability Index:      {self.instability_index:.2f} ({"Stable" if self.is_stable_based_on_ii else "Unstable"})

GRand AVerage of hYdropathicity (GRAVY), based on Kyte-Doolittle scale: {self.gravy_kd:.2f}
"""

    def __str__(self):
        return f"[START]\n{self.report}\n[END]"


class ProtParamBatch:
    def __init__(self, seq_list: list[str | Seq | SeqRecord]):
        self.seq_list = seq_list

        self.protparam_res = [ProtParam(seq) for seq in seq_list]

    def report(self):
        pass
