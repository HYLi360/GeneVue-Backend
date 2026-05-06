#   Copyright (C) 2025-2026, HYLi360.
#   Free software distributed under the terms of the GNU GPL-3.0 license,
#   and comes with ABSOLUTELY NO WARRANTY.
#   See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

import numpy as np
from scipy import stats
from rich.table import Table
from genevue import console
from rich import box
from typing import Optional


class SingleDiscriptive:
    def __init__(self, data: np.typing.NDArray):
        self.data = data
        self.n = len(data)
        self.mean = data.mean()
        data.sort()
        self.data_sorted = data
        self.median = (
            self.data_sorted[(self.n - 1) // 2]
            if self.n % 2 == 1
            else (
                self.data_sorted[(self.n - 1) // 2]
                + self.data_sorted[(self.n - 1) // 2 + 1]
            )
            / 2
        )
        self.max = data[-1]
        self.min = data[0]
        self.var = data.var(ddof=1)
        self.sst = data.var(ddof=1) * len(data)
        self.std = data.std(ddof=1)
        self.ser = data.std(ddof=1) / np.sqrt(len(data))
        self.ci95 = stats.norm.ppf(q=0.975) * self.ser
        self.show()

    def show(self):
        print()
        table_desc = Table(
            title="Descriptive Statistics",
            title_justify="left",
            show_header=True,
            header_style="bold magenta",
            highlight=True,
        )
        table_desc.add_column("N", justify="right", no_wrap=True)
        table_desc.add_column("Mean", justify="right", no_wrap=True)
        table_desc.add_column("Median", justify="right", no_wrap=True)
        table_desc.add_column("Maximum", justify="right", no_wrap=True)
        table_desc.add_column("Minimum", justify="right", no_wrap=True)
        table_desc.add_column("Var", justify="right", no_wrap=True)
        table_desc.add_column("Total SS", justify="right", no_wrap=True)
        table_desc.add_column("Std", justify="right", no_wrap=True)
        table_desc.add_column("SE", justify="right", no_wrap=True)
        table_desc.add_column("CI95", justify="right", no_wrap=True)
        table_desc.add_row(
            f"{self.n}",
            f"{self.mean:.2f}",
            f"{self.median:.2f}",
            f"{self.max:.2f}",
            f"{self.min:.2f}",
            f"{self.var:.2f}",
            f"{self.sst:.2f}",
            f"{self.std:.2f}",
            f"{self.ser:.2f}",
            f"[{self.mean-self.ci95:.2f}, {self.mean+self.ci95:.2f}]",
        )
        table_desc.caption = (
            "CI95 is based on the assumption of a normal distribution\n"
        )
        table_desc.caption_justify = "left"
        table_desc.box = box.SIMPLE_HEAVY
        console.print(table_desc)


class DoubleDiscriptive:
    def __init__(
        self,
        data1: np.typing.NDArray,
        data2: np.typing.NDArray,
        data1_name: Optional[str] = None,
        data2_name: Optional[str] = None,
    ):
        self.var_name1 = "X1" if not data1_name else data1_name
        self.var_name2 = "X2" if not data2_name else data2_name
        self.data1 = data1
        self.n1 = len(data1)
        self.mean1 = data1.mean()
        data1.sort()
        self.data_sorted1 = data1
        self.median1 = (
            self.data_sorted1[(len(data1) - 1) // 2]
            if len(data1) % 2 == 1
            else (
                self.data_sorted1[(len(data1) - 1) // 2]
                + self.data_sorted1[(len(data1) - 1) // 2 + 1]
            )
            / 2
        )
        self.max1 = data1.max()
        self.min1 = data1.min()
        self.var1 = data1.var(ddof=1)
        self.sst1 = data1.var(ddof=1) * len(data1)
        self.std1 = data1.std(ddof=1)
        self.ser1 = data1.std(ddof=1) / np.sqrt(len(data1))
        self.ci95_1 = stats.norm.ppf(q=0.975) * self.ser1

        self.data2 = data2
        self.n2 = len(data2)
        self.mean2 = data2.mean()
        data2.sort()
        self.data_sorted2 = data2
        self.median2 = (
            self.data_sorted2[(len(data2) - 1) // 2]
            if len(data2) % 2 == 1
            else (
                self.data_sorted2[(len(data2) - 1) // 2]
                + self.data_sorted2[(len(data2) - 1) // 2 + 1]
            )
            / 2
        )
        self.max2 = data2.max()
        self.min2 = data2.min()
        self.var2 = data2.var(ddof=1)
        self.sst2 = data2.var(ddof=1) * len(data2)
        self.std2 = data2.std(ddof=1)
        self.ser2 = data2.std(ddof=1) / np.sqrt(len(data2))
        self.ci95_2 = stats.norm.ppf(q=0.975) * self.ser2
        self.show()

    def show(self):
        print()
        table_desc = Table(
            title="Descriptive Statistics",
            show_header=True,
            header_style="bold magenta",
        )
        table_desc.add_column("Name", justify="right", no_wrap=True)
        table_desc.add_column("N", justify="right", no_wrap=True)
        table_desc.add_column("Mean", justify="right", no_wrap=True)
        table_desc.add_column("Median", justify="right", no_wrap=True)
        table_desc.add_column("Maximum", justify="right", no_wrap=True)
        table_desc.add_column("Minimum", justify="right", no_wrap=True)
        table_desc.add_column("Var", justify="right", no_wrap=True)
        table_desc.add_column("Total SS", justify="right", no_wrap=True)
        table_desc.add_column("Std", justify="right", no_wrap=True)
        table_desc.add_column("SE", justify="right", no_wrap=True)
        table_desc.add_column("CI95", justify="right", no_wrap=True)
        table_desc.add_row(
            f"{self.var_name1}",
            f"{self.n1}",
            f"{self.mean1:.2f}",
            f"{self.median1:.2f}",
            f"{self.max1:.2f}",
            f"{self.min1:.2f}",
            f"{self.var1:.2f}",
            f"{self.sst1:.2f}",
            f"{self.std1:.2f}",
            f"{self.ser1:.2f}",
            f"[{self.mean1-self.ci95_1:.2f}, {self.mean1+self.ci95_1:.2f}]",
        )
        table_desc.add_row(
            f"{self.var_name2}",
            f"{self.n2}",
            f"{self.mean2:.2f}",
            f"{self.median2:.2f}",
            f"{self.max2:.2f}",
            f"{self.min2:.2f}",
            f"{self.var2:.2f}",
            f"{self.sst2:.2f}",
            f"{self.std2:.2f}",
            f"{self.ser2:.2f}",
            f"[{self.mean2-self.ci95_2:.2f}, {self.mean2+self.ci95_2:.2f}]",
        )
        table_desc.caption = "CI95 is based on the assumption of a normal distribution"
        table_desc.caption_justify = "left"
        table_desc.box = box.SIMPLE_HEAVY
        console.print(table_desc)
        print()


def _num_formatter(column_number_sequence: list, max_digit: int = 3):
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

    new_number_sequence = (
        new_number_sequence[0] if len(new_number_sequence) == 1 else new_number_sequence
    )
    return new_number_sequence


def _pvalue_formatter(column_number_sequence: list):
    new_number_sequence = [
        "<0.001" if i < 0.001 else f"{i:.3f}" for i in column_number_sequence
    ]
    return new_number_sequence
