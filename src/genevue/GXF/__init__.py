#  Copyright (c) 2026 HYLi360. All rights reserved.
#
#  see LICENSE in /LICENSE
#  see side-package LICENSEs (if used) in /LICENSE_OF_SIDE_PACKAGES

import re

feature_line = re.compile(
    rb"^(?P<chr>[^\t#]*)\t[^\t]*\t(?P<type>[^\t#]*)\t(?P<start>\d+)\t(?P<end>\d+)\t[^\t]*\t(?P<strand>[+-])\t("
    rb"?P<phase>[^\t]*)\t(?P<attrs>[^\n]*)$",
    re.M,
)
