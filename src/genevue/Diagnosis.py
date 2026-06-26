# Copyright (c) 2026 HYLi360. All rights reserved.
#
# see LICENSE in /LICENSE
# see side-package LICENSEs (if used) in /LICENSE_OF_SIDE_PACKAGES

import importlib.metadata
import os
import platform
import sys
import time

from genevue import console
from genevue.QC.BUSCO import __version__ as version_busco


class Diagnosis:
    def __init__(self):
        # OS infomation
        self._fetch_os_info()

        # architecture of CPU
        self._fetch_hardware_info()

        # python runtime infomation
        self._fetch_runtime_info()

        # version of key dependences
        self._fetch_dependences_version()

    def _fetch_os_info(self):
        if os.name == "nt":
            """Windows"""
            self.os_info = " ".join(
                [
                    " ".join(platform.win32_ver()),
                    (
                        platform.win32_edition()
                        if platform.win32_edition() is not None
                        else ""
                    ),
                    "IoT" if platform.win32_is_iot() else "",
                ]
            )
        elif os.name == "posix":
            """Linux, macOS or BSD"""
            self.os_info = " ".join(
                [platform.system(), platform.release(), *platform.libc_ver()]
            )
        else:
            self.os_info = "Unknown"

    def _fetch_hardware_info(self):
        self.arch = platform.machine()

    def _fetch_runtime_info(self):
        self.python_version = platform.python_version()
        self.python_implementation = platform.python_implementation()
        self.python_executable_path = sys.executable.replace(
            os.path.expanduser("~"), "~"
        )
        self.python_venv_prefix = sys.prefix.replace(os.path.expanduser("~"), "~")
        self.python_venv_base_prefix = sys.base_prefix.replace(
            os.path.expanduser("~"), "~"
        )
        self.python_build_time = time.strptime(
            platform.python_build()[1], "%b %d %Y %H:%M:%S"
        )
        self.python_build_time_formated = time.strftime(
            "%Y-%m-%d %H:%M:%S", self.python_build_time
        )
        self.python_compiler = platform.python_compiler()
        self.current_wd = os.getcwd().replace(os.path.expanduser("~"), "~")
        self.is_venv = os.environ.get("VIRTUAL_ENV")
        self.is_conda_venv = os.environ.get("CONDA_DEFAULT_ENV")
        self.is_64bit = sys.maxsize > 2**32
        self.is_isolated = bool(sys.flags.isolated)
        self.is_optimized = bool(sys.flags.optimize)
        self.is_repl = hasattr(sys, "ps1")

    def _fetch_dependences_version(self):
        self.dependences_version = {}
        for key_package in [
            "aiohttp",
            "bgzip",
            "biopython",
            "logomaker",
            "lxml",
            "matplotlib",
            "newick",
            "numpy",
            "pandas",
            "polars",
            "pyarrow",
            "requests",
            "rich",
            "scikit-learn",
            "scipy",
            "seaborn",
            "typer",
        ]:
            self.dependences_version[key_package] = importlib.metadata.version(
                key_package
            )

        self.dependences_version["busco"] = version_busco

    def export_simple(self):
        console.print(
            "\n"
            f"[bold]Operator System[/bold]    [cyan]{self.os_info}[/cyan]\n"
            f"[bold]Python[/bold]             [cyan]{self.python_version}[/cyan]\n"
            f"                   on {self.python_executable_path}\n"
            f"[bold]Biopython[/bold]          [cyan]{self.dependences_version['biopython']}[/cyan]\n"
            f"[bold]NumPy[/bold]              [cyan]{self.dependences_version['numpy']}[/cyan]\n"
            f"[bold]pandas[/bold]             [cyan]{self.dependences_version['pandas']}[/cyan]\n"
            f"[bold]Polars[/bold]             [cyan]{self.dependences_version['polars']}[/cyan]\n"
            f"[bold]SciPy[/bold]              [cyan]{self.dependences_version['scipy']}[/cyan]\n"
            f"[bold]BUSCO (integrated)[/bold] [cyan]{self.dependences_version['busco']}[/cyan]\n",
        )
