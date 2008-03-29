"""
Quotient Ring Elements

AUTHOR:
    -- William Stein

TODO: This implementation is very basic.
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

from sage.interfaces.singular import singular as singular_default

class QuotientRingElement(ring_element.RingElement):
    def __init__(self, parent, rep, reduce=True):
        """
        An element of a quotient ring $R/I$.

        EXAMPLES:
            sage: R.<x> = PolynomialRing(ZZ)
            sage: S.<xbar> = R.quo((4 + 3*x + x^2, 1 + x^2)); S
            Quotient of Univariate Polynomial Ring in x over Integer Ring by the ideal (x^2 + 3*x + 4, x^2 + 1)
            sage: v = S.gens(); v
            (xbar,)

            sage: loads(v[0].dumps()) == v[0]
            True

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S = R.quo(x^2 + y^2); S
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
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
            -x*y^2 + y^2
        """
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

    def __nonzero__(self):
        return self.__rep not in self.parent().defining_ideal()

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
        """
        Return the leading term of this quotient ring element.

        EXAMPLE:
            sage: R.<x,y,z>=PolynomialRing(GF(7),3,order='lex')
            sage: I = sage.rings.ideal.FieldIdeal(R)
            sage: Q = R.quo( I )
            sage: f = Q( z*y + 2*x )
            sage: f.lt()
            2*xbar
        """
        return QuotientRingElement(self.parent(),self.__rep.lt())

    def lm(self):
        """
        Return the leading monomial of this quotient ring element.

        EXAMPLE:
            sage: R.<x,y,z>=PolynomialRing(GF(7),3,order='lex')
            sage: I = sage.rings.ideal.FieldIdeal(R)
            sage: Q = R.quo( I )
            sage: f = Q( z*y + 2*x )
            sage: f.lm()
            xbar
        """
        return QuotientRingElement(self.parent(),self.__rep.lm())

    def lc(self):
        """
        Return the leading coefficent of this quotient ring element.

        EXAMPLE:
            sage: R.<x,y,z>=PolynomialRing(GF(7),3,order='lex')
            sage: I = sage.rings.ideal.FieldIdeal(R)
            sage: Q = R.quo( I )
            sage: f = Q( z*y + 2*x )
            sage: f.lc()
            2
        """
        return self.__rep.lc()

    def variables(self):
        return [QuotientRingElement(self.parent(),v) for v in self.__rep.variables()]

    def monomials(self):
        return [QuotientRingElement(self.parent(),m) for m in self.__rep.monomials()]

    def _singular_(self, singular=singular_default):
        """
        Return Singular representation of self.

        INPUT:
            singular -- a non-standard interpreter may be provided

        EXAMPLE:
            sage: P.<x,y>  = PolynomialRing(GF(2),2)
            sage: I = sage.rings.ideal.FieldIdeal(P)
            sage: Q = P.quo(I)
            sage: Q._singular_()
            //   characteristic : 2
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C
            // quotient ring from ideal
            _[1]=x2+x
            _[2]=y2+y
            sage: xbar = Q(x); xbar
            xbar
            sage: xbar._singular_()
            x
            sage: Q(xbar._singular_()) # a round-trip
            xbar
        """
        return self.__rep._singular_(singular)

    def reduce(self, G):
        r"""
        Reduce this quotient ring element by a set of quotient ring
        elements \var{G}.

        INPUT:
            G -- a list of quotient ring elements

        EXAMPLE:
            sage: P.<a,b,c,d,e> = PolynomialRing(GF(2), 5, order='lex')
            sage: I1 = ideal([a*b + c*d + 1, a*c*e + d*e, a*b*e + c*e, b*c + c*d*e + 1])
            sage: Q = P.quotient( sage.rings.ideal.FieldIdeal(P) )
            sage: I2 = ideal([Q(f) for f in I1.gens()])
            sage: f = Q((a*b + c*d + 1)^2  + e)
            sage: f.reduce(I2.gens())
            ebar
        """
        try:
            G = [f.lift() for f in G]
        except TypeError:
            pass
        # reduction w.r.t. the defining ideal is performed in the
        # constructor
        return QuotientRingElement(self.parent(), self.__rep.reduce(G))
