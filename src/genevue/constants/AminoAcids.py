from typing import Literal


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


from typing import Literal


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
