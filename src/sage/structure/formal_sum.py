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

from sage.modules.module import Module
from sage.structure.element import ModuleElement
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.structure.parent_base import ParentWithBase

class FormalSums_generic(Module):
    def __init__(self, base=ZZ):
        ParentWithBase.__init__(self, base)

    def _repr_(self):
        return "Abelian Group of all Formal Finite Sums over %s"%self.base_ring()

    def _an_element_impl(self):
        return FormalSum([], check=False, reduce=False)

    def __call__(self, x, check=True, reduce=True):
        if isinstance(x, FormalSum):
            P = x.parent()
            if P is self:
                return x
            elif P == self:
                return FormalSum(x._data, check=False, reduce=False, parent=self)
            else:
                x = x._data
        if isinstance(x, list):
            return FormalSum(x, check=check,reduce=reduce,parent=self)
        if x == 0:
            return FormalSum([], check=False, reduce=False, parent=self)
        else:
            return FormalSum([(self.base_ring()(1), x)], check=False, reduce=False, parent=self)

    def _coerce_impl(self, x):
        if x == 0:
            return FormalSum([], check=False, reduce=False, parent=self)
        elif isinstance(x, FormalSum) and x.parent().has_coerce_map_from(self.base_ring()):
            return self(x)
        else:
            return FormalSum([(self.base_ring()(1), x)], check=False, reduce=False, parent=self)

    def base_extend(self, R):
        if R.has_coerce_map_from(self.base_ring()):
            return FormalSums(R)

import weakref
cache = {}
def FormalSums(R):
    try:
        F = cache[R]()
        if not F is None:
            return F
    except KeyError:
        pass
    F = FormalSums_generic(R)
    cache[R] = weakref.ref(F)
    return F


formal_sums = FormalSums_generic()

class FormalSum(ModuleElement):
    def __init__(self, x, parent=formal_sums, check=True, reduce=True):
        if x == 0:
            x = []
        if check:
            k = parent.base_ring()
            try:
                x = [(k(t[0]), t[1]) for t in x]
            except (IndexError, KeyError), msg:
                raise TypeError, "%s\nInvalid formal sum"%msg
        self._data = x
        if parent is None:
            parent = formal_sums
        ModuleElement.__init__(self, parent)
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

    def __cmp__(self, right):
        return cmp(self._data, right._data)

    def _neg_(self):
        return self.__class__([(-c, s) for (c, s) in self._data], check=False, parent=self.parent())

    def _add_(self, other):
        return self.__class__(self._data + other._data, check=False, parent=self.parent())

    def _lmul_(self, s):
        return self.__class__([(c*s, x) for (c, x) in self], check=False, parent=self.parent())

    def _rmul_(self, s):
        return self.__class__([(s*c, x) for (c, x) in self], check=False, parent=self.parent())

    def __nonzero__(self):
        if len(self._data) == 0:
            return False
        for c, _ in self._data:
            if not c.is_zero():
                return True
        return False

    def reduce(self):
        if len(self) == 0:
            return
        v = [(x,c) for c, x in self if not c.is_zero()]
        if len(v) == 0:
            self._data = v
            return
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
        if not coeff.is_zero():
            w.append((coeff,last))
        self._data = w
