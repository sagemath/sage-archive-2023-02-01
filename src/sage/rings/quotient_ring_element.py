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

from __future__ import with_statement

import operator

import ring_element
import quotient_ring

class QuotientRingElement(ring_element.RingElement):
    """
    An element of a quotient ring $R/I$.

    EXAMPLES:
        sage: R.<x> = PolynomialRing(ZZ)
        sage: S.<xbar> = R.quo((4 + 3*x + x^2, 1 + x^2)); S
        Quotient of Univariate Polynomial Ring in x over Integer Ring by the ideal (x^2 + 1, x^2 + 3*x + 4)
        sage: v = S.gens(); v
        (xbar,)

        sage: loads(v[0].dumps()) == v[0]
        True

        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S = R.quo(x^2 + y^2); S
        Quotient of Polynomial Ring in x, y over Rational Field by the ideal (y^2 + x^2)
        sage: S.gens()
        (xbar, ybar)

    We name each of the generators.
        sage: S.<a,b> = R.quotient(x^2 + y^2)
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
        from sage.structure.parent_gens import localvars
        P = self.parent()
        R = P.cover_ring()
        # We print by temporarily (and safely!) changing the variable
        # names of the covering structure R to those of P.
        # These names get changed back, since we're using "with".
        with localvars(R, P.variable_names(), normalize=False):
            return str(self.__rep)

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

    def __neg__(self):
        return QuotientRingElement(self.parent(), -self.__rep)

    def __pos__(self):
        return self

    def __invert__(self):
        inv = self.__rep.inverse_mod(self.parent().defining_ideal())
        return QuotientRingElement(self.parent(), inv)

    def __float__(self):
        return float(self.lift())

    def __cmp__(self, other):
        if self.__rep == other.__rep or ((self.__rep - other.__rep) in self.parent().defining_ideal()):
            return 0
        return -1

    def lt(self):
        return QuotientRingElement(self.parent(),self.__rep.lt())

    def lm(self):
        return QuotientRingElement(self.parent(),self.__rep.lm())

    def lc(self):
        return QuotientRingElement(self.parent(),self.__rep.lc())

    def variables(self):
        return [QuotientRingElement(self.parent(),v) for v in self.__rep.variables()]

    def monomials(self):
        return [QuotientRingElement(self.parent(),m) for m in self.__rep.monomials()]
