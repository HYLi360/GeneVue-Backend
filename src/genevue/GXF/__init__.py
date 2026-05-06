#  Copyright (C) 2025-2026, HYLi360.
#  Free software distributed under the terms of the GNU GPL-3.0 license,
#  and comes with ABSOLUTELY NO WARRANTY.
#  See at <https://www.gnu.org/licenses/gpl-3.0.en.html>

import re

feature_line = re.compile(
    rb"^(?P<chr>[^\t#]*)\t[^\t]*\t(?P<type>[^\t#]*)\t(?P<start>\d+)\t(?P<end>\d+)\t[^\t]*\t(?P<strand>[+-])\t("
    rb"?P<phase>[^\t]*)\t(?P<attrs>[^\n]*)$",
    re.M,
)
