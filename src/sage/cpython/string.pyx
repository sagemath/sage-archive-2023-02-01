# -*- encoding: utf-8 -*-
"""
String <-> bytes encoding/decoding
"""

#*****************************************************************************
#       Copyright (C) 2017 Erik M. Bray <erik.bray@lri.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import

import sys


# Provide this as a shortcut to calling sys.getfilesystemencoding(), which
# after interpreter initialization is constant.
FS_ENCODING = sys.getfilesystemencoding()

# Functions in this module are implemented in the .pxd file for inlining.
