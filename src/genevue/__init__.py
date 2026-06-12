#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

from __future__ import annotations
from pathlib import Path
from genevue.logsystem import LogLevel, setup_rich_logger, _setup_busco_bridge
from genevue.console import console

__version__ = "0.0.11"
__nickname__ = "TESTING"
__full_version__ = f"{__version__} {__nickname__}"


# ---all custom class (inheritances from BaseException and Warning) we used---------------------------
class BaseException4GeneVue(BaseException):
    pass


class BaseWarning4GeneVue(UserWarning):
    pass


class ProgramNotFoundError(BaseException4GeneVue):
    """Called when the external executable file does not exist."""

    def __init__(self, program: str) -> None:
        self.message = f"Can't Find the Executable File of This Program: {program}"

    def __str__(self) -> str:
        return f"{self.message}"


class SeqFileNotFoundError(BaseException4GeneVue):
    """Called when the external executable file does not exist."""

    def __init__(self, seq_file_path: str | Path) -> None:
        self.message = f"Can't Find This Sequence File: {seq_file_path}"

    def __str__(self) -> str:
        return f"{self.message}"


class GFF3FileNotFoundError(BaseException4GeneVue):
    """Called when the external executable file does not exist."""

    def __init__(self, gff3_file_path: str | Path) -> None:
        self.message = f"Can't Find This SimpleGFF3 File: {gff3_file_path}"

    def __str__(self) -> str:
        return f"{self.message}"


class FileNotExistsError(BaseException4GeneVue):

    def __init__(self, file_path: Path) -> None:
        self.message = f"Can not found this file: {file_path}"

    def __str__(self) -> str:
        return f"{self.message}"


class FormatNotSuitableError(BaseException4GeneVue):

    def __init__(self, message: str) -> None:
        self.message = message

    def __str__(self) -> str:
        return self.message


class NotImplementedMethodError(BaseException4GeneVue):

    def __init__(self, method_name: str):
        self.message = f"Not implemented method: {method_name}"


class DirNotExistsError(BaseException4GeneVue):
    pass


class DirNotEmptyError(BaseException4GeneVue):
    pass


class AllFieldsEmptyError(BaseException4GeneVue):

    def __init__(self, fields: list):
        self.message = f"Fill in at least one field: {', '.join(fields)}"

    def __str__(self):
        return self.message


class NothingFoundError(BaseException4GeneVue):

    def __init__(self, message: str):
        self.message = message

    def __str__(self):
        return self.message


class ResultIsNotSpeciesError(BaseException4GeneVue):

    def __init__(self, result: str):
        self.message = f"Expect get information of species, but get {result}"

    def __str__(self):
        return self.message


__all__ = [
    "__full_version__",
    "console",
    "setup_rich_logger",
    "BaseException4GeneVue",
    "SeqFileNotFoundError",
    "GFF3FileNotFoundError",
    "FileNotExistsError",
    "NotImplementedMethodError",
    "DirNotExistsError",
    "DirNotEmptyError",
    "AllFieldsEmptyError",
    "NothingFoundError",
    "ResultIsNotSpeciesError",
]
