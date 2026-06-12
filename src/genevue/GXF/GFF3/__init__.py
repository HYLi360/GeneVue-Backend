import sys
import urllib.parse
from pathlib import Path
from typing import Optional

from genevue import console, setup_rich_logger

logger = setup_rich_logger(__name__, console)


def _url_style_decode(string: str) -> str:
    """
    Decode string which in URL-style coding.
    """
    return urllib.parse.unquote(string)


def _parse_key_equal_value(string: str) -> dict[str, str]:
    if "=" not in string:
        return {}
    else:
        key_str, value_str = string.strip().split("=")
        return {_url_style_decode(key_str): _url_style_decode(value_str)}
