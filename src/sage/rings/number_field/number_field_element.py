"""
Number Field Elements
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@ucsd.edu>
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
import sage.rings.coerce as coerce

import sage.rings.field_element as field_element
import sage.rings.infinity as infinity
import sage.rings.polynomial_element as polynomial
import sage.rings.polynomial_ring as polynomial_ring
import sage.rings.rational_field as rational_field
import sage.rings.rational as rational
import sage.rings.integer_ring as integer_ring
import sage.rings.arith as arith
import sage.misc.misc as misc

import number_field

from sage.libs.all import pari_gen

QQ = rational_field.RationalField()

def is_NumberFieldElement(x):
    return isinstance(x, NumberFieldElement)

class NumberFieldElement(field_element.FieldElement):
    """
    An element of a number field.
    """
    def __init__(self, parent, f):
        """
        INPUT:
            parent -- a number field
            f -- defines an element of a number field.

        EXAMPLES:
        The following examples illustrate creation of elements of
        number fields, and some basic arithmetic.

        First we define a polynomial over Q.
            sage: x = PolynomialRing(QQ).0
            sage: f = x^2 + 1

        Next we use f to define the number field.
            sage: K = NumberField(f, "a"); K
            Number Field in a with defining polynomial x^2 + 1
            sage: a = K.gen()
            sage: a^2
            -1
            sage: (a+1)^2
            2*a
            sage: a^2
            -1
            sage: z = K(5); 1/z
            1/5

        We create a cube root of 2.
            sage: K = NumberField(x^3 - 2, "b")
            sage: b = K.gen()
            sage: b^3
            2
            sage: (b^2 + b + 1)^3
            12*b^2 + 15*b + 19

        This example illustrates save and load:
            sage: K, a = NumberField(x^17 - 2, 'a').objgen()
            sage: s = a^15 - 19*a + 3
            sage: loads(s.dumps()) == s
            True
        """
        ppr = parent.polynomial_ring()
        if isinstance(parent, number_field.NumberField_extension):
            ppr = parent.base_field().polynomial_ring()

        if isinstance(parent, str):
            print 'parent = ', parent
            print 'f = ', f
        if isinstance(f, pari_gen):
            f = f.lift()
            f = ppr(f)
        if not isinstance(f, polynomial.Polynomial):
            f = ppr(f)
        if f.degree() >= parent.degree():
            if isinstance(parent, number_field.NumberField_extension):
                f %= parent.absolute_polynomial()
            else:
                f %= parent.polynomial()
        field_element.FieldElement.__init__(self, parent)
        self.__element = f
        self.__field = parent

    def __repr__(self):
        x = self.__element
        return str(x).replace(x.parent().variable_name(),self.parent().variable_name())

    def _im_gens_(self, codomain, im_gens):
        # NOTE -- if you ever want to change this so relative number fields are
        # in terms of a root of a poly.
        # The issue is that elements of a relative number field are represented in terms
        # of a generator for the absolute field.  However the morphism gives the image
        # of gen, which need not be a generator for the absolute field.  The morphism
        # has to be *over* the relative element.
        return codomain(self.__element(im_gens[0]))

    def _latex_(self):
        return self.polynomial()._latex_(name=self.__field.variable_name())

    def _pari_(self, var=None):
        """
        Return PARI C-library object representation of self.
        """
        try:
            return self.__pari
        except AttributeError:
            if var == None:
                var = self.parent().variable_name()
            f = self.__element._pari_().subst("x",var)
            g = self.parent().polynomial()._pari_().subst("x",var)
            return f.Mod(g)

    def _pari_init_(self, var=None):
        """
        Return GP/PARI string representation of self.
        """
        if var == None:
            var = self.parent().variable_name()
        f = self.__element._pari_().subst("x",var)
        g = self.parent().polynomial()._pari_().subst("x",var)
        return 'Mod(%s, %s)'%(f,g)

    def __cmp__(self, other):
        if not isinstance(other, NumberFieldElement) or self.parent() != other.parent():
            return coerce.cmp(self, other)
        return misc.generic_cmp(self.__element, other.__element)

    def _add_(self, other):
        return NumberFieldElement(self.parent(), self.__element+other.__element)

    def __getitem__(self, n):
        return self.polynomial()[n]

    def _mul_(self, other):
        """
        Returns the product of self and other as elements of a number field.
        """
        return NumberFieldElement(self.parent(), self.__element * other.__element)

        #NOTES: In LiDIA, they build a multiplication table for the
        #number field, so it's not necessary to reduce modulo the
        #defining polynomial every time:
        #     src/number_fields/algebraic_num/order.cc: compute_table
        # but asymptotically fast poly multiplication means it's
        # actually faster to *not* build a table!

    def __pow__(self, right):
        right = int(right)
        if right < 0:
            x = self.__invert__()
            right *= -1
            return arith.generic_power(x, right, one=self.parent()(1))
        return arith.generic_power(self, right, one=self.parent()(1))

    def _sub_(self, other):
        return NumberFieldElement(self.parent(), self.__element - other.__element)

    def _div_(self, other):
        return self * other.__invert__()

    def __floordiv__(self, other):
        return self / other

    def __neg__(self):
        return NumberFieldElement(self.parent(), -self.__element)

    def __int__(self):
        return int(self.__element)

    def __long__(self):
        return long(self.__element)

    def __invert__(self):
        if self.__element.is_zero():
            raise ZeroDivisionError
        g, _, a = self.parent().polynomial().xgcd(self.__element)
        assert g.degree() == 0
        c = g[0]
        return NumberFieldElement(self.parent(), (1/c)*a)

    def _integer_(self):
        return integer_ring.IntegerRing()(int(self))

    def _rational_(self):
        if self.__element.degree() >= 1:
            raise TypeError, "Unable to coerce %s to a rational"%self
        return self.__element[0]

    def polynomial(self):
        return self.__element

    def denominator(self):
        """
        Return the denominator of this element, which is by definition
        the denominator of the corresponding polynomial
        representation.  I.e., elements of number fields are
        represented as a polynomial (in reduced form) modulo the
        modulus of the number field, and the denominator is the
        denominator of this polynomial.

        EXAMPLES:
            sage: K.<z> = CyclotomicField(3)
            sage: a = 1/3 + (1/5)*z
            sage: print a.denominator()
            15
        """
        return self.__element.denominator()

    def _set_multiplicative_order(self, n):
        self.__multiplicative_order = n

    def multiplicative_order(self):
        try:
            return self.__multiplicative_order
        except AttributeError:
            pass

        if self.__element == 1:
            self.__multiplicative_order = 1
            return self.__multiplicative_order
        if self.__element == -1:
            self.__multiplicative_order = 2
            return self.__multiplicative_order

        if isinstance(self.parent(), number_field.NumberField_cyclotomic):
            t = self.parent().multiplicative_order_table()
            f = self.__element
            if t.has_key(f):
                self.__multiplicative_order = t[f]
                return self.__multiplicative_order

        ####################################################################
        # VERY DUMB Algorithm to compute the multiplicative_order of
        # an element x of a number field K.
        #
        # 1. Find an integer B such that if n>=B then phi(n) > deg(K).
        #    For this use that for n>6 we have phi(n) >= log_2(n)
        #    (to see this think about the worst prime factorization
        #    in the multiplicative formula for phi.)
        # 2. Compute x, x^2, ..., x^B in order to determine the multiplicative_order.
        #
        # todo-- Alternative: Only do the above if we don't require an optional
        # argument which gives a multiple of the order, which is usually
        # something available in any actual application.
        #
        # BETTER TODO: Factor cyclotomic polynomials over K to determine
        # possible orders of elements?  Is there something even better?
        #
        ####################################################################
        d = self.parent().degree()
        B = max(7, 2**d+1)
        x = self
        i = 1
        while i < B:
            if x == 1:
                self.__multiplicative_order = i
                return self.__multiplicative_order
            x *= self
            i += 1

        # it must have infinite order
        self.__multiplicative_order = infinity.infinity
        return self.__multiplicative_order

    def trace(self):
        return QQ(self._pari_().trace())
        # return self.matrix().trace()

    def norm(self):
        return QQ(self._pari_().norm())
        #return self.matrix().determinant()

    def charpoly(self):
        r"""
        The characteristic polynomial of this element over $\Q$.

        EXAMPLES:

        We compute the charpoly of cube root of $3$.

            sage: R.<x> = QQ['x']
            sage: K.<a> = NumberField(x^3-2)
            sage: a.charpoly()
            x^3 - 2

        We construct a relative extension and find the characteristic
        polynomial over $\Q$.

            sage: S.<X> = K['X']
            sage: L.<b> = NumberField(X^3 + 17)
            sage: L
            Extension by X^3 + 17 of the Number Field in a with defining polynomial x^3 - 2
            sage: a = L.0; a
            b
            sage: a.charpoly()
            x^9 + 57*x^6 + 165*x^3 + 6859
        """
        R = self.parent().polynomial_ring()
        if not isinstance(self.parent(), number_field.NumberField_extension):
            return R(self._pari_().charpoly())
        else:
            g = self.polynomial()  # in QQ[x]
            f = self.parent().pari_polynomial()  # # field is QQ[x]/(f)
            R = g.parent()
            return R( (g._pari_().Mod(f)).charpoly() )

## This might be useful for computing relative charpoly.
## BUT -- currently I don't even know how to view elements
## as being in terms of the right thing, i.e., this code
## below as is lies.
##             nf = self.parent()._pari_base_nf()
##             prp = self.parent().pari_relative_polynomial()
##             elt = str(self.polynomial()._pari_())
##             return R(nf.rnfcharpoly(prp, elt))
##         # return self.matrix().charpoly()

    def minpoly(self):
        # The minimal polynomial is square-free and
        # divisible by same irreducible factors as
        # the characteristic polynomial.
        # TODO: factoring to find the square-free part is idiotic.
        # Instead use a GCD algorithm!
        f = polynomial_ring.PolynomialRing(QQ)(1)
        for g, _ in self.charpoly().factor():
            f *= g
        return f

    def matrix(self):
        r"""
        The matrix of right multiplication by the element on the power
        basis $1, x, x^2, \ldots, x^{d-1}$ for the number field.  Thus
        the {\em rows} of this matrix give the images of each of the $x^i$.
        """
        # Mutiply each power of field generator on
        # the left by this element; make matrix
        # whose rows are the coefficients of the result,
        # and transpose.
        try:
            return self.__matrix
        except AttributeError:
            K = self.parent()
            v = []
            x = K.gen()
            a = K(1)
            d = K.degree()
            for n in range(d):
                v += (a*self).list()
                a *= x
            Q = rational_field.RationalField()
            import sage.matrix.matrix_space
            M = sage.matrix.matrix_space.MatrixSpace(Q, d)
            self.__matrix = M(v)
            return self.__matrix


    def list(self):
        """
        EXAMPLE:
            sage: K.<z> = CyclotomicField(3)
            sage: (2+3/5*z).list()
            [2, 3/5]
            sage: (5*z).list()
            [0, 5]
            sage: K(3).list()
            [3, 0]
        """
        n = self.parent().degree()
        v = self.__element[:n]
        z = rational.Rational(0)
        return v + [z]*(n - len(v))

