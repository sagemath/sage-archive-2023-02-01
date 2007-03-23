"""
Morphisms between modular abelian varieties.

AUTHOR:
    -- William Stein (2007-03)
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

import sage.categories.all

class Morphism(sage.categories.all.Morphism):
    def __init__(self, parent, x):
        self._parent = parent
        self._x = x

    def _repr_(self):
        return "Morphism defined by %s in %s"%(self._x, self._parent)
