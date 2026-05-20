"""
Localisation system.

Usage:
  _ = init_translations(lang="zh_CN", ts_path="/path/to/ts_file")

- lang: Literal
- ts_path: can set manually, but "configure.ts_path" is recommended
"""

from __future__ import annotations

import gettext
from pathlib import Path
from typing import Literal, Callable
from genevue import console

_translator = gettext.gettext


def _(message: str) -> str:
    return _translator(message)


def init_translations(lang: Literal["zh_CN", "en_US"], ts_path: Path) -> Callable:
    global _translator

    gettext.bindtextdomain("messages", str(ts_path))
    gettext.textdomain("messages")

    try:
        translation = gettext.translation(
            domain="messages",
            localedir=str(ts_path),
            languages=[lang],
        )
        translation.install("messages")
        _translator = translation.gettext
    except FileNotFoundError:
        console.print(
            f"Translation file (.mo) not found at {ts_path}. "
            f"Falling back to English."
        )

    return _
