# -*- coding: utf-8 -*-
"""
Environment Variables

This module defines the following subset of the Sage environment
variables:

* ``SAGE_ROOT``
* ``SAGE_SRC``
"""


#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os



try:
    SAGE_SRC = os.environ['SAGE_SRC']
except KeyError:
    SAGE_SRC = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))

try:
    SAGE_ROOT = os.environ['SAGE_ROOT']
except KeyError:
    SAGE_ROOT = os.path.dirname(SAGE_SRC)


assert os.path.isfile(os.path.join(SAGE_ROOT, 'configure.ac'))
assert os.path.isfile(os.path.join(SAGE_SRC, 'sage_bootstrap', 'setup.py'))
