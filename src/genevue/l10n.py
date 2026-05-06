"""
Reading the config file and binding translation function before import.
"""

import gettext
from pathlib import Path

from genevue import local_configure
from genevue import console

# Get current language
LANGUAGE = local_configure.language
TRANSLATE_PATH = Path(local_configure.ts_path)

gettext.bindtextdomain("messages", str(TRANSLATE_PATH))
gettext.textdomain("messages")

# Initialize the gettext
# Fetch the .mo (compiled .po)
try:
    translation = gettext.translation(
        domain="messages",
        localedir=f"{TRANSLATE_PATH}",
        languages=[LANGUAGE],
    )
    translation.install("messages")
    _ = translation.gettext
except FileNotFoundError:
    console.warn(
        f"Not found the translate file (.mo). Rollback to English.\nTried path: {TRANSLATE_PATH}"
    )
    _ = gettext.gettext
