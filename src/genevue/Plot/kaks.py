"""
genevue.plot.kaks 的 Docstring
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from typing import Literal


def tribe_peak(
    ka_seq: np.ndarray,
    ks_seq: np.ndarray,
    omega_seq: np.ndarray,
    yscale: Literal["symlog", "linear"] = "linear",
) -> None:
    # two Kernel Density Estimation (KDE) methods
    kde_ka = stats.gaussian_kde(ka_seq)
    kde_ks = stats.gaussian_kde(ks_seq)
    kde_omega = stats.gaussian_kde(omega_seq)

    # plot three KDE curves
    plt.yscale(yscale)
    x_eval = np.linspace(0, 1.5, num=1000)
    plt.plot(x_eval, kde_ka(x_eval), label="Ka")
    plt.plot(x_eval, kde_ks(x_eval), label="Ks")
    plt.plot(x_eval, kde_omega(x_eval), label="Omega")

    plt.legend()
    plt.grid()
    plt.savefig("f.png")
    plt.show()
