"""
Precision Error.

The errors in this file indicate various styles of precision problems
that can go wrong for p-adics and power series.

AUTHORS::

    - David Roe
"""

#*****************************************************************************
#       Copyright (C) 2008 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

class PrecisionError(Exception):
    pass

class HaltingError(PrecisionError):
    pass

class PrecisionLimitError(PrecisionError):
    pass
