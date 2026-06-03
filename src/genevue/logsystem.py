from __future__ import annotations

import logging
import os
import time
from enum import Enum
from pathlib import Path
from typing import Optional

import atexit
from rich.console import Console
from rich.logging import RichHandler

_log_file_path: Optional[Path] = None


class LogLevel(Enum):
    DEBUG = 10
    INFO = 20
    WARNING = 30
    ERROR = 40
    CRITICAL = 50


def _setup_busco_bridge(rich_console: Console):
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


def _clean_empty_log_file():
    """
    Clean empty log file. ...may useless?
    """
    global _log_file_path
    if _log_file_path is None:
        return
    else:
        try:
            if _log_file_path.exists() and _log_file_path.stat().st_size == 0:
                os.remove(_log_file_path)
        except OSError:
            # not found log file
            # expective result
            pass


def _get_log_file_handler(log_dir: Path, loglevel: int | str) -> logging.FileHandler:
    global _log_file_path
    if _log_file_path is None:
        _log_file_path = log_dir / f"genevue_{time.strftime("%Y%m%d%H%M%S")}.log"
        atexit.register(_clean_empty_log_file)

    # add delay=True make file handler activated only when needed
    file_handler = logging.FileHandler(_log_file_path, delay=True)
    file_handler.setLevel(loglevel)
    file_handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s %(name)s %(levelname)s %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    return file_handler


def setup_rich_logger(
    name: str,
    rich_console: Console,
    total_loglevel: int | str = "DEBUG",
    disp_loglevel: int | str = "INFO",
    file_loglevel: int | str = "DEBUG",
    propagate: bool = False,
    *,
    log_dir: Path = Path(".").resolve(),
) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(total_loglevel)

    if getattr(logger, "_genevue_configured", False):
        return logger

    rich_handler = RichHandler(
        console=rich_console,
        show_time=True,
        rich_tracebacks=True,
        markup=True,
        show_path=False,
        tracebacks_show_locals=True,
    )
    rich_handler.setLevel(disp_loglevel)

    file_handler = _get_log_file_handler(log_dir, file_loglevel)

    logger.addHandler(file_handler)
    logger.addHandler(rich_handler)
    logger.propagate = propagate
    logger._genevue_configured = True

    return logger
