"""
Torsion points on modular abelan varieties.

AUTHOR:
    -- William Stein (2007-03)
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.structure.sage_object import SageObject

class TorsionPoint(SageObject):
    def __init__(self, abvar, x):
        self._abvar = abvar
        self._x = x

    def _repr_(self):
        return "Torsion point defined by %s on %s"%\
               (self._abvar, self._x)

