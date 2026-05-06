# Copyright (C) 2025-2026, HYLi360.
# Free software distributed under the terms of the GNU GPL-3.0 license,
# and comes with ABSOLUTELY NO WARRANTY.
# See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

"""
A basic framework for statistical computation and analysis. Its current functionality involves performing a series of
statistical tests and returning results to the user based on the corresponding distributions.

It makes use of the SciPy stats module to some extent; however, for educational purposes, the actual computations do not
always require it.

Results should be presented at three confidence levels (alpha value): 0.950, 0.990, and 0.999. If necessary, results can
also be visualized by calling the `show()` function.
"""

import numpy as np
from scipy import stats
from rich.table import Table
from rich.text import Text
from rich import box

from genevue.Statistic import SingleDiscriptive, DoubleDiscriptive, _num_formatter
from genevue import console


class TTest1SMP(SingleDiscriptive):
    def __init__(self, data: np.typing.NDArray, mu: float):
        super().__init__(data)
        self.mu = mu
        self.df = len(data) - 1
        self.shapiro = stats.shapiro(self.data)
        self.t = (self.mean - self.mu) / (self.std / np.sqrt(self.n))
        self.left = stats.t.cdf(x=self.t, df=self.df)
        self.right = 1 - self.left
        self.both = 2 * min(self.left, self.right)
        self.ci95 = stats.t.ppf(df=self.df, q=0.975) * (self.std / np.sqrt(self.n))
        self.show_all()

    def show_all(self):
        table_norm = Table(
            title="Normality Hypothesis Test",
            title_justify="left",
            show_header=True,
            header_style="bold magenta",
            highlight=True,
            width=64,
        )
        table_norm.add_column("Statistic", justify="right", no_wrap=True)
        table_norm.add_column("df", justify="right", no_wrap=True)
        table_norm.add_column("Sig.", justify="right", no_wrap=True)
        table_norm.add_row(
            _num_formatter([self.shapiro.statistic]),
            f"{self.df}",
            _num_formatter([self.shapiro.pvalue]),
        )
        table_norm.caption = (
            "based on Shapiro-Wilk method\n"
            "Low p-value suggests a violation of the assumption of normality\n"
        )
        table_norm.caption_justify = "left"
        table_norm.box = box.SIMPLE_HEAVY
        console.print(table_norm)

        table_ttest = Table(
            title="One-Sample t-Test",
            title_justify="left",
            show_header=True,
            header_style="bold magenta",
            highlight=True,
        )
        table_ttest.add_column("Statistic", justify="right", no_wrap=True)
        table_ttest.add_column("df", justify="right", no_wrap=True)
        table_ttest.add_column("Sig.", justify="right", no_wrap=True)
        table_ttest.add_column("Mean", justify="right", no_wrap=True)
        table_ttest.add_column("mu", justify="right", no_wrap=True)
        table_ttest.add_column("diff", justify="right", no_wrap=True)
        table_ttest.add_column("CI95", justify="right", no_wrap=True)
        table_ttest.add_row(
            f"{self.t:.3f}",
            f"{self.df}",
            f"{self.both:.3f}",
            f"{_num_formatter([self.mean])}",
            f"{_num_formatter([self.mu])}",
            f"{_num_formatter([self.mean-self.mu])}",
            f"[{self.mean-self.mu-self.ci95:.3f}, {self.mean-self.mu+self.ci95:.3f}]",
        )
        table_ttest.box = box.SIMPLE_HEAVY
        console.print(table_ttest)


class TTestIndeSMP(DoubleDiscriptive):
    def __init__(self, data1: np.typing.NDArray, data2: np.typing.NDArray):
        super().__init__(data1, data2)

        self.shapiro1 = stats.shapiro(self.data1)
        self.shapiro2 = stats.shapiro(self.data2)
        self.levene = stats.levene(self.data1, self.data2, center="mean")

        self.t = stats.ttest_ind(self.data1, self.data2, equal_var=True)
        self.var = ((self.n1 - 1) * self.var1 + (self.n2 - 1) * self.var2) / (
            self.n1 + self.n2 - 2
        )
        self.ser = np.sqrt(self.var * (1 / self.n1 + 1 / self.n2))
        self.left = stats.t.cdf(x=self.t.statistic, df=self.t.df)
        self.right = 1 - self.left
        self.both = 2 * min(self.left, self.right)
        self.ci95 = stats.t.ppf(df=self.t.df, q=0.975) * self.ser

        self.welth_t = stats.ttest_ind(self.data1, self.data2, equal_var=False)
        self.welth_left = stats.t.cdf(x=self.welth_t.statistic, df=self.welth_t.df)
        self.welth_right = 1 - self.welth_left
        self.welth_both = 2 * min(self.welth_left, self.welth_right)
        self.welth_ci95 = stats.t.ppf(df=self.welth_t.df, q=0.975) * self.ser

        print(self)

    def __str__(self):
        return (
            f"+------------------------+\n"
            "| Descriptive Statistics |\n"
            "+------------------------+\n"
            f"N      =  {self.n1},  {self.n2}\n"
            f"mean   = {self.mean1: .2f}, {self.mean2: .2f}\n"
            f"median = {self.median1: .2f}, {self.median2: .2f}\n"
            f"max    = {self.max1: .2f}, {self.max2: .2f}\n"
            f"min    = {self.min1: .2f}, {self.max2: .2f}\n"
            f"Var    = {self.var1: .2f}, {self.var2: .2f}\n"
            f"SST    = {self.sst1: .2f}, {self.sst2: .2f}\n"
            f"SD     = {self.std1: .2f}, {self.std2: .2f}\n"
            f"SE     = {self.ser1: .2f}, {self.ser2: .2f}\n"
            f"CI95   =[{self.mean1-self.ci95_1: .2f}, {self.mean1+self.ci95_1: .2f}], "
            f"[{self.mean2-self.ci95_2: .2f}, {self.mean2+self.ci95_2: .2f}]\n"
            "\n"
            "+------------------------------------------------------------+\n"
            "| Normality Hypothesis Test (based on Shapiro-Wilk's method) |\n"
            "+------------------------------------------------------------+\n"
            f"W1 = {self.shapiro1.statistic:.3} (p-value = {self.shapiro1.pvalue:.3}){" (!)" if self.shapiro1.pvalue < 0.1 else ""}\n"
            f"W2 = {self.shapiro2.statistic:.3} (p-value = {self.shapiro2.pvalue:.3}){" (!)" if self.shapiro2.pvalue < 0.1 else ""}\n"
            "* The small p-value suggests a violation of the assumption of normality.\n"
            "\n"
            "+---------------------------------------------------------+\n"
            "| Homogeneity of Variance Test (based on Levene's method) |\n"
            "+---------------------------------------------------------+\n"
            f"L = {self.levene.statistic:.3} (p-value = {self.levene.pvalue:.3}){" (!)" if self.levene.pvalue < 0.1 else ""}\n"
            "* The small p-value suggests that the populations do not have equal variances.\n"
            "\n"
            "+----------------------------+\n"
            "| Independent Samples t-Test |\n"
            "+----------------------------+\n"
            f"Student's t-Test:      df = {self.t.df:.3f}\n"
            f"              t Statistic = {self.t.statistic:.3f}\n"
            f"        p-value for  left = {self.left:.3f} (mean2<mean1)\n"
            f"                for right = {self.right:.3f} (mean2>mean1)\n"
            f"                for  both = {self.both:.3f} (mean2!=mean1)\n"
            f"            Mean-diff SER = {self.ser:.3f}\n"
            f"         95% Mean-diff CI = [{self.mean2-self.mean1-self.ci95:.3f}, {self.mean2-self.mean1:.3f}, "
            f"{self.mean2-self.mean1+self.ci95:.3f}]\n"
            f"Welth's t-Test:        df = {self.welth_t.df:.3f}\n"
            f"              t Statistic = {self.welth_t.statistic:.3f}\n"
            f"        p-value for  left = {self.welth_left:.3f} (mean2<mean1)\n"
            f"                for right = {self.welth_right:.3f} (mean2>mean1)\n"
            f"                for  both = {self.welth_both:.3f} (mean2!=mean1)\n"
            f"         95% Mean-diff CI = [{self.mean2-self.mean1-self.welth_ci95:.3f}, {self.mean2-self.mean1:.3f}, "
            f"{self.mean2-self.mean1+self.welth_ci95:.3f}]\n"
        )
