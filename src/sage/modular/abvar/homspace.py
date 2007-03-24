"""
Spaces of homomorphisms between modular abelian varieties.

AUTHOR:
    -- William Stein (2007-03)
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.categories.homset import HomsetWithBase

class Homspace(HomsetWithBase):
    def __init__(self, domain, codomain):
        self._domain = domain
        self._codomain = codomain

    def _repr_(self):
        return "Space of homomorphisms from %s to %s"%\
               (self._domain, self._codomain)
