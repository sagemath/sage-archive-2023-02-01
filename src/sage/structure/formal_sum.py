"""
Formal sums
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.misc.misc
import element
import sage.misc.latex

from sage.ext.element import AdditiveGroupElement
from sage.ext.group import AbelianGroup

class FormalSums(AbelianGroup):
    def _repr_(self):
        return "Abelian Group of all Formal Finite Sums"

formal_sums = FormalSums()

class FormalSum(AdditiveGroupElement, list):
    def __init__(self, x, parent=None, check=True, reduce=True):
        if x == 0:
            x = []
        if check:
            for t in self:
                if not (isinstance(t, tuple) and len(t) == 2):
                    raise TypeError, "Invalid formal sum"
        list.__init__(self, x)
        if parent is None:
            parent = formal_sums
        AdditiveGroupElement.__init__(self, parent)
        self.reduce()

    def _repr_(self):
        symbols = [z[1] for z in self]
        coeffs= [z[0] for z in self]
        return sage.misc.misc.repr_lincomb(symbols, coeffs)

    def _latex_(self):
        symbols = [z[1] for z in self]
        coeffs= [z[0] for z in self]
        return sage.misc.latex.repr_lincomb(symbols, coeffs)

    def _add_(self, other):
        return self.__class__(list.__add__(self,other), parent=self.parent())

    def __mul__(self, s):
        return self.__class__([(c*s, x) for c,x in self], check=False, parent=self.parent())

    def __rmul__(self, s):
        return self.__class__([(s*c, x) for c,x in self], parent=self.parent())

    def reduce(self):
        if len(self) == 0:
            return
        v = [(x,c) for c, x in self]
        v.sort()
        w = []
        last = v[0][0]
        coeff = v[0][1]
        for x, c in v[1:]:
            if x == last:
                coeff += c
            else:
                if coeff != 0:
                    w.append((coeff, last))
                last = x
                coeff = c
        w.append((coeff,last))
        list.__init__(self, w)


