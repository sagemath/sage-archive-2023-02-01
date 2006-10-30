"""
Quotient Ring Elements

AUTHOR:
    -- William Stein
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

import operator

import ring_element
import quotient_ring

from sage.rings.coerce import bin_op

from sage.structure.all import Element_cmp_

class QuotientRingElement(Element_cmp_, ring_element.RingElement):
    """
    An element of a quotient ring $R/I$.

    EXAMPLES:
        sage: R = PolynomialRing(ZZ,'x')
        sage: S = R/(4 + 3*x + x^2, 1 + x^2); S
        Quotient of Univariate Polynomial Ring in x over Integer Ring by the ideal (x^2 + 1, x^2 + 3*x + 4)
        sage: v = S.gens(); v
        (x,)

        sage: loads(v[0].dumps()) == v[0]  # todo: not implemented
        True

        sage: R = PolynomialRing(QQ, 'x,y')
        sage: S = R/(x^2 + y^2); S
        Quotient of Polynomial Ring in x, y over Rational Field by the ideal (y^2 + x^2)
        sage: S.gens()
        (x, y)

    We name each of the generators.
        sage: S = R.quotient(x^2 + y^2, 'a,b')
        sage: a
        a
        sage: b
        b
        sage: a^2 + b^2 == 0
        True
        sage: b.lift()
        y
        sage: (a^3 + b^2).lift()
        y^2 - x*y^2
    """
    def __init__(self, parent, rep, reduce=True):
        ring_element.RingElement.__init__(self, parent)
        self.__rep = rep
        if reduce:
            self._reduce_()

    def _reduce_(self):
        I = self.parent().defining_ideal()
        self.__rep = I.reduce(self.__rep)

    def copy(self):
        return QutientRingElement(self.parent(), self.__rep)

    def lift(self):
        return self.__rep

    def is_zero(self):
        return self.__rep in self.parent().defining_ideal()

    def is_unit(self):
        if self.__rep.is_unit():
            return True
        raise NotImplementedError

    def _repr_(self):
        P = self.parent()
        R = P.cover_ring()
        tmp = R._names_from_obj(P)
        s = str(self.__rep)
        R._names_from_obj(tmp)
        return s

    def _add_(self, right):
        return QuotientRingElement(self.parent(), self.__rep + right.__rep)

    def _sub_(self, right):
        return QuotientRingElement(self.parent(), self.__rep - right.__rep)

    def _mul_(self, right):
        return QuotientRingElement(self.parent(), self.__rep * right.__rep)

    def _div_(self, right):
        if not right.is_unit():
            raise ZeroDivisionError
        raise NotImplementedError

    def __int__(self):
        return int(self.lift())

    def _integer_(self):
        try:
            return self.lift()._integer_()
        except AttributeError:
            raise NotImplementedError

    def _rational_(self):
        try:
            return self.lift()._rational_()
        except AttributeError:
            raise NotImplementedError

    def __long__(self):
        return int(self.lift())

    def __rdiv__(self, left):
        return self.parent()(left)/self

    def __pow__(self, right):
        right = int(right)
        return QuotientRingElement(self.parent(), self.__rep**right)

    def __neg__(self):
        return QuotientRingElement(self.parent(), -self.__rep)

    def __pos__(self):
        return self

    def __invert__(self):
        inv = self.__rep.inverse_mod(self.parent().defining_ideal())
        return QuotientRingElement(self.parent(), inv)

    def __float__(self):
        return float(self.lift())

    def _cmp_(self, other):
        if (self.__rep - other.__rep) in self.parent().defining_ideal():
            return 0
        return -1

    def lt(self):
        return self.__rep.lt()

    def lm(self):
        return self.__rep.lm()

    def lc(self):
        return self.__rep.lc()

    def variables(self):
        return self.__rep.variables()

    def monomials(self):
        return self.__rep.monomials()
