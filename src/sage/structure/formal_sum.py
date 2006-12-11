"""
Formal sums

AUTHORS:
   -- David Harvey (2006-09-20): changed FormalSum not to derive from
      "list" anymore, because that breaks new Element interface
   -- Nick Alexander (2006-12-06): added test cases.

EXAMPLES:
    sage: A = FormalSum([(1, 2/3)]); A
    2/3
    sage: B = FormalSum([(3, 1/5)]); B
    3*1/5
    sage: -B
    -3*1/5
    sage: A + B
    3*1/5 + 2/3
    sage: A - B
    -3*1/5 + 2/3
    sage: B*3
    9*1/5
    sage: 2*A
    2*2/3
    sage: list(2*A + A)
    [(3, 2/3)]
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

from sage.groups.group import AbelianGroup
from sage.structure.element import AdditiveGroupElement
from sage.rings.integer_ring import ZZ
from sage.structure.parent_base import ParentWithBase

class FormalSums(AbelianGroup):
    def __init__(self):
        ParentWithBase.__init__(self, ZZ)

    def _repr_(self):
        return "Abelian Group of all Formal Finite Sums"

    def _an_element_impl(self):
        return FormalSum([], check=False, reduce=False)

formal_sums = FormalSums()

class FormalSum(AdditiveGroupElement):
    # NOTE: originally FormalSum also inherited from "list". But when
    # we made Element (base type of AdditiveGroupElement) have a cdef'ed
    # _parent attribute, it became impossible to multiply inherit.
    # So I added a "_data" attribute to FormalSum which contains the
    # list information. I also gave this class some iterator semantics.
    # So currently this class is very inefficiently implemented, and
    # also I think in the long run it would be better to require
    # users to call a FormalSum.list() method explicitly to get the list.
    # For now I'm trying to keep things from breaking with minimal effort.
    # (David H, 2006-09-20)
    def __init__(self, x, parent=None, check=True, reduce=True):
        if x == 0:
            x = []
        if check:
            for t in x:
                if not (isinstance(t, tuple) and len(t) == 2):
                    raise TypeError, "Invalid formal sum"
        self._data = x
        if parent is None:
            parent = formal_sums
        AdditiveGroupElement.__init__(self, parent)
        self.reduce()

    def __iter__(self):
        return iter(self._data)

    def __getitem__(self, n):
        return self._data[n]

    def __len__(self):
        return len(self._data)

    def _repr_(self):
        symbols = [z[1] for z in self]
        coeffs= [z[0] for z in self]
        return sage.misc.misc.repr_lincomb(symbols, coeffs)

    def _latex_(self):
        symbols = [z[1] for z in self]
        coeffs= [z[0] for z in self]
        return sage.misc.latex.repr_lincomb(symbols, coeffs)

    def _neg_(self):
        return self.__class__([(-c, s) for (c, s) in self._data], check=False, parent=self.parent())

    def _add_(self, other):
        return self.__class__(self._data + other._data, check=False, parent=self.parent())

    def _lmul_(self, s):
        return self._rmul_(s) # formal sums are abelian

    def _rmul_(self, s):
        return self.__class__([(s*c, x) for (c, x) in self], check=False, parent=self.parent())

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
        self._data = w
