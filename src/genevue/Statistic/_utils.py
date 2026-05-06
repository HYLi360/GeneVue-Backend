#   Copyright (C) 2025-2026, HYLi360.
#   Free software distributed under the terms of the GNU GPL-3.0 license,
#   and comes with ABSOLUTELY NO WARRANTY.
#   See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

import numpy as np
from rich.table import Table


def num_formatter(column_number_sequence: list, max_digit: int = 3):
    integer, decimal = list(
        zip(*[str(np.float64(num)).split(".") for num in column_number_sequence])
    )
    decimal_digits = min(max_digit, max([len(dec) for dec in decimal]))

    _, decimal2 = list(
        zip(*[f"{num:.{decimal_digits}f}".split(".") for num in column_number_sequence])
    )
    new_number_sequence = []

    for i, d, d2 in zip(integer, decimal, decimal2):
        if int(i) == 0:
            i = ""
        if int(f"{d}{"0"*(decimal_digits-len(d))}") == int(d2):
            new_number_sequence.append(f"{i}.{d}{" "*(decimal_digits-len(d))}")
        else:
            new_number_sequence.append(f"{i}.{d2}")

    return new_number_sequence


def pvalue_formatter(column_number_sequence: list):
    new_number_sequence = [
        "<0.001" if i < 0.001 else f"{i:.3f}" for i in column_number_sequence
    ]
    return new_number_sequence
