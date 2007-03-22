"""
L-series of modular abelian varieties.

AUTHOR:
    -- William Stein (2007-03)
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.structure.sage_object import SageObject

class LSeries(SageObject):
    def __init__(self, abvar):
        self._abvar = abvar

    def _repr_(self):
        return "L-Series attached to %s"%self

class LSeries_complex(LSeries):

    def __call__(self, x):
        raise NotImplementedError

    def rational_part(self, x):
        raise NotImplementedError


class LSeries_padic(LSeries):

    def __call__(self, x):
        raise NotImplementedError

    def rational_part(self, x):
        raise NotImplementedError
