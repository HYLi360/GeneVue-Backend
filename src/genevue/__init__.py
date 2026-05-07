#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

from __future__ import annotations

import logging
import os
from pathlib import Path
from enum import Enum

from genevue.configure import Configure
from rich.console import Console as RichConsole
from rich.logging import RichHandler
from rich.traceback import install

__version__ = "0.0.2"
__nickname__ = "TESTING"

__full_version__ = f"{__version__} {__nickname__}"


class LogLevel(Enum):
    DEBUG = -1
    INFO = 11
    WARNING = 21
    ERROR = 31
    CRITICAL = 41


# Start Rich Engine.
# Add `warn` `error` `debug` method to make it much like JavaScript
# and give me more control of terminal output.
class GeneVueConsole(RichConsole):
    def __init__(self, name="app"):
        super().__init__()
        self.loglevel: LogLevel = LogLevel.INFO
        self.logger = logging.getLogger(name)
        self.logger.setLevel(self.loglevel.value)
        self.logger.propagate = False
        self.ppid = str(os.getppid())
        self.pid = str(os.getcwd())

        # register for logging file handler
        self._file_hdlr = None

        if not self.logger.handlers:
            handler = RichHandler(
                show_time=True,
                show_level=True,
                show_path=False,
                rich_tracebacks=True,
            )
            handler.setFormatter(logging.Formatter("%(message)s"))
            self.logger.addHandler(handler)

    def debug(self, *args):
        self.logger.debug(" ".join(map(str, args)), stacklevel=2)

    def info(self, *args):
        self.logger.info(" ".join(map(str, args)), stacklevel=2)

    def warn(self, *args):
        self.logger.warning(" ".join(map(str, args)), stacklevel=2)

    def error(self, *args):
        self.logger.error(" ".join(map(str, args)), stacklevel=2)

    def exception(self, exception: BaseException4GeneVue, *args):
        if not args:
            self.logger.error(exception, stacklevel=2)
        raise exception

    def set_new_loglevel(self, new_loglevel: LogLevel):
        self.loglevel = new_loglevel
        self.logger.setLevel(self.loglevel.value)

    # write log to one file
    def add_file_handler(self, path: str, level: int = logging.DEBUG):
        """Attach a FileHandler to this console's logger."""
        if self._file_hdlr is not None:
            self.remove_file_handler()
        self._file_hdlr = logging.FileHandler(path, mode="a")
        self._file_hdlr.setLevel(level)
        self._file_hdlr.setFormatter(
            logging.Formatter(
                "%(asctime)s %(levelname)s:%(name)s\t%(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )
        self.logger.addHandler(self._file_hdlr)

    def remove_file_handler(self):
        """Remove and close the FileHandler if one is attached."""
        if self._file_hdlr is not None:
            self.logger.removeHandler(self._file_hdlr)
            self._file_hdlr.close()
            self._file_hdlr = None

    def has_file_handler(self):
        return self._file_hdlr is not None


console = GeneVueConsole()


def _setup_busco_bridge():
    """Configure a parent logger so all BUSCO module loggers propagate
    to a RichHandler for unified console output."""
    busco_parent = logging.getLogger("genevue.QC.BUSCO")
    busco_parent.setLevel(logging.DEBUG)
    busco_parent.propagate = False
    if not busco_parent.handlers:
        handler = RichHandler(
            show_time=True,
            show_level=True,
            show_path=False,
            rich_tracebacks=True,
        )
        handler.setFormatter(logging.Formatter("%(message)s"))
        busco_parent.addHandler(handler)


_setup_busco_bridge()
install(show_locals=True)


# ---all custom class (inheritances from BaseException and Warning) we used---------------------------
class BaseException4GeneVue(BaseException):
    def __init__(self, *args: object) -> None:
        super().__init__(*args)


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


local_configure = Configure()


__all__ = [
    "__full_version__",
    "console",
    "SeqFileNotFoundError",
    "GFF3FileNotFoundError",
    "FileNotExistsError",
    "NotImplementedMethodError",
    "DirNotExistsError",
    "DirNotEmptyError",
    "AllFieldsEmptyError",
    "NothingFoundError",
    "ResultIsNotSpeciesError",
    "local_configure",
]
