#  Copyright (c) 2026 HYLi360. All rights reserved.
#
#  see license in LICENSE
#  see side-package licenses in LICENSE_OF_SIDE_PACKAGES

from typing import Dict, Literal

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AACmap:
    def __init__(self, cmap: Literal["Clustal2", "RasMol"] = "Clustal2"):
        self.cmap = cmap

        match self.cmap:
            case "Clustal2":
                self.palette = {
                    "A": "#80A0F0",
                    "C": "#80A0F0",
                    "I": "#80A0F0",
                    "L": "#80A0F0",
                    "M": "#80A0F0",
                    "F": "#80A0F0",
                    "W": "#80A0F0",
                    "V": "#80A0F0",
                    "N": "#00FF00",
                    "Q": "#00FF00",
                    "S": "#00FF00",
                    "T": "#00FF00",
                    "K": "#F01505",
                    "R": "#F01505",
                    "D": "#C048C0",
                    "E": "#C048C0",
                    "H": "#15A4A4",
                    "Y": "#15A4A4",
                    "G": "#F09048",
                    "P": "#FFFF00",
                }
            case "RasMol":
                self.palette = {
                    "D": "#E60A0A",
                    "E": "#E60A0A",
                    "C": "#E6E600",
                    "M": "#E6E600",
                    "K": "#145AFF",
                    "R": "#145AFF",
                    "S": "#FA9600",
                    "T": "#FA9600",
                    "F": "#3232AA",
                    "Y": "#3232AA",
                    "N": "#00DCDC",
                    "Q": "#00DCDC",
                    "G": "#EBEBEB",
                    "L": "#0F820F",
                    "V": "#0F820F",
                    "I": "#0F820F",
                    "A": "#C8C8C8",
                    "W": "#B45AB4",
                    "H": "#8282D2",
                    "P": "#DC9682",
                }
        self.aalabels = list(self.palette.keys())

    def color(self, amino_acid):
        return self.palette.get(amino_acid, "F0F0F0")


BLOSUM62 = {
    "A": 0.078,
    "C": 0.024,
    "D": 0.052,
    "E": 0.059,
    "F": 0.044,
    "G": 0.083,
    "H": 0.025,
    "I": 0.062,
    "K": 0.056,
    "L": 0.092,
    "M": 0.024,
    "N": 0.041,
    "P": 0.043,
    "Q": 0.034,
    "R": 0.051,
    "S": 0.059,
    "T": 0.055,
    "V": 0.072,
    "W": 0.014,
    "Y": 0.034,
}

MOLECULAR_FORMULA_AMINOACID: Dict[str, Dict[str, int]] = {
    "A": {"C": 3, "H": 7, "N": 1, "O": 2, "S": 0},
    "C": {"C": 3, "H": 7, "N": 1, "O": 2, "S": 1},
    "D": {"C": 4, "H": 7, "N": 1, "O": 4, "S": 0},
    "E": {"C": 5, "H": 9, "N": 1, "O": 4, "S": 0},
    "F": {"C": 9, "H": 11, "N": 1, "O": 2, "S": 0},
    "G": {"C": 2, "H": 5, "N": 1, "O": 2, "S": 0},
    "H": {"C": 6, "H": 9, "N": 3, "O": 2, "S": 0},
    "I": {"C": 6, "H": 13, "N": 1, "O": 2, "S": 0},
    "K": {"C": 6, "H": 14, "N": 2, "O": 2, "S": 0},
    "L": {"C": 6, "H": 13, "N": 1, "O": 2, "S": 0},
    "M": {"C": 5, "H": 11, "N": 1, "O": 2, "S": 1},
    "N": {"C": 4, "H": 8, "N": 2, "O": 3, "S": 0},
    "P": {"C": 5, "H": 9, "N": 1, "O": 2, "S": 0},
    "Q": {"C": 5, "H": 10, "N": 2, "O": 3, "S": 0},
    "R": {"C": 6, "H": 14, "N": 4, "O": 2, "S": 0},
    "S": {"C": 3, "H": 7, "N": 1, "O": 3, "S": 0},
    "T": {"C": 4, "H": 9, "N": 1, "O": 3, "S": 0},
    "V": {"C": 5, "H": 11, "N": 1, "O": 2, "S": 0},
    "W": {"C": 11, "H": 12, "N": 2, "O": 2, "S": 0},
    "Y": {"C": 9, "H": 11, "N": 1, "O": 3, "S": 0},
}


def molecular_formula_protein(protein_seq: str | Seq | SeqRecord) -> tuple[str, int]:
    res = {"C": 0, "H": 0, "N": 0, "O": 0, "S": 0}
    for aa in protein_seq:
        faa: Dict[str, int] = MOLECULAR_FORMULA_AMINOACID.get(
            aa.upper(), {"C": 0, "H": 0, "N": 0, "O": 0, "S": 0}
        )
        res = {
            "C": res["C"] + faa["C"],
            "H": res["H"] + faa["H"],
            "N": res["N"] + faa["N"],
            "O": res["O"] + faa["O"],
            "S": res["S"] + faa["S"],
        }
    # minus mass of H2O
    res["H"] = res["H"] - (len(protein_seq) - 1) * 2
    res["O"] = res["O"] - (len(protein_seq) - 1)

    res_ls, total_atoms = [], 0

    for element in res:
        if res[element] == 0:
            res_ls.append(f"")
        elif res[element] == 1:
            res_ls.append(f"{element}")
            total_atoms += 1
        else:
            res_ls.append(f"{element}_{res[element]}")
            total_atoms += res[element]

    return " ".join(res_ls), total_atoms
