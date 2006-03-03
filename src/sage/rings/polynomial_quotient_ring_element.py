"""
Elements of Quotients of Univariate Polynomial Rings
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator
import coerce
import sage.structure.element as element
import sage.rings.arith as arith
import sage.misc.misc as misc
import commutative_ring_element

import polynomial_element
import polynomial_quotient_ring



class PolynomialQuotientRingElement(commutative_ring_element.CommutativeRingElement):
    """
    Element of a quotient of a polynomial ring.
    """
    def __init__(self, parent, polynomial, check=True):
        """
        Create an element of the quotient of a polynomial ring.

        INPUT:
            parent -- a quotient of a polynomial ring
            polynomial -- a polynomial
            check -- bool (optional): whether or not to verify
                     that x is a valid element of the polynomial
                     ring and reduced (mod the modulus).
        """
        commutative_ring_element.CommutativeRingElement.__init__(self, parent)
        if check:
            if not isinstance(parent, polynomial_quotient_ring.PolynomialQuotientRing_generic):
                raise TypeError, "parent must be a polynomial quotient ring"

            if not isinstance(polynomial, polynomial_element.Polynomial):
                raise TypeError, "polynomial must be a polynomial"

            if not polynomial in parent.polynomial_ring():
                raise TypeError, "polynomial must be in the polynomial ring of the parent"

        f = parent.modulus()
        if polynomial.degree() >= f.degree():
            try:
                polynomial %= f
            except AttributeError:
                A = polynomial
                B = f
                R = A
                Q = B.parent()(0)
                X = B.parent().gen()
                while R.degree() >= B.degree():
                    S = (R.leading_coefficient()/B.leading_coefficient()) * X**(R.degree()-B.degree())
                    Q = Q + S
                    R = R - S*B
                polynomial = R
        self.__polynomial = polynomial

    def _im_gens_(self, codomain, im_gens):
        return self.__polynomial._im_gens_(codomain, im_gens)

    def __reduce__(self):
        return PolynomialQuotientRingElement, (self.parent(), self.__polynomial, False)

    def _repr_(self):
        return self.__polynomial._repr(self.parent().variable_name())

    def __add__(self, right):
        """
        Return the sum of two polynomial ring quotient elements.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3-2, 'a'); a = S.gen()
            sage: (a^2 - 4) + (a+2)
            a^2 + a - 2
            sage: int(1) + a
            a + 1
        """
        if not isinstance(right, PolynomialQuotientRingElement) \
               or self.parent() != right.parent():
            return coerce.bin_op(self, right, operator.add)
        return PolynomialQuotientRingElement(self.parent(),
                               self.__polynomial + right.__polynomial, check=False)

    def __cmp__(self, other):
        """
        Compare this element with something else, where equality
        testing coerces the object on the right, if possible (and
        necessary).

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3-2, 'a'); a = S.gen()
            sage: S(1) == 1
            True
            sage: a^3 == 2
            True

        For the purposes of comparison in SAGE the quotient element
        $a^3$ is equal to $x^3$.  This is because when the comparison
        is performed, the right element is coerced into the parent of
        the left element, and $x^3$ coerces to $a^3$.

            sage: a == x
            True
            sage: a^3 == x^3
            True
            sage: x^3
            x^3
            sage: S(x^3)
            2
        """
        if not isinstance(other, PolynomialQuotientRingElement) or other.parent() != self.parent():
            return coerce.cmp(self, other)
        if self.parent() != other.parent():
            return -1
        return misc.generic_cmp(self.__polynomial, other.__polynomial)

    def __div__(self, right):
        """
        Return the quotient of two polynomial ring quotient elements.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3-2, 'a'); a = S.gen()
            sage: (a^2 - 4) / (a+2)
            a - 2
        """
        if not isinstance(right, PolynomialQuotientRingElement) \
               or self.parent() != right.parent():
            return coerce.bin_op(self, right, operator.div)
        return self * ~right

    def __getitem__(self, n):
        return self.__polynomial[n]

    def __int__(self):
        """
        Coerce this element to an int if possible.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3-2, 'a'); a = S.gen()
            sage: int(S(10))
            10
            sage: int(a)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to int
        """
        return int(self.__polynomial)

    def __invert__(self):
        if self.__polynomial.is_zero():
            raise ZeroDivisionError, \
               "element %s of quotient polynomial ring not invertible"%self
        g, _, a = self.parent().modulus().xgcd(self.__polynomial)
        if g.degree() != 0:
            raise ZeroDivisionError, \
               "element %s of quotient polynomial ring not invertible"%self
        c = g[0]
        return PolynomialQuotientRingElement(self.parent(), (~c)*a, check=False)

    def __long__(self):
        """
        Coerce this element to a long if possible.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3-2, 'a'); a = S.gen()
            sage: long(S(10))
            10L
            sage: long(a)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to long
        """
        return long(self.__polynomial)

    def __mul__(self, right):
        """
        Return the product of two polynomial ring quotient elements.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3-2, 'a'); a = S.gen()
            sage: (a^2 - 4) * (a+2)
            2*a^2 - 4*a - 6
        """
        if not isinstance(right, PolynomialQuotientRingElement) \
               or self.parent() != right.parent():
            return coerce.bin_op(self, right, operator.mul)
        R = self.parent()
        prod = self.__polynomial * right.__polynomial
        return PolynomialQuotientRingElement(R, prod, check=False)

    def __neg__(self):
        return PolynomialQuotientRingElement(self.parent(), -self.__polynomial)

    def __pow__(self, n):
        """
        Return a power of a polynomial ring quotient element.

        EXAMPLES:
            sage: R = PolynomialRing(Integers(9), 'x'); x = R.gen()
            sage: S = R.quotient(x^4 + 2*x^3 + x + 2, 'a'); a = S.gen()
            sage: a^100
            7*a^3 + 8*a + 7
        """
        n = int(n)
        if n < 0:
            x = self.__invert__()
            n *= -1
            return arith.generic_power(x, n, one=self.parent()(1))
        return arith.generic_power(self, n, one=self.parent()(1))

    def __sub__(self, right):
        """
        Return the difference of two polynomial ring quotient elements.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3-2, 'a'); a = S.gen()
            sage: (a^2 - 4) - (a+2)
            a^2 - a - 6
            sage: int(1) - a
            -a + 1
        """
        if not isinstance(right, PolynomialQuotientRingElement) \
               or self.parent() != right.parent():
            return coerce.bin_op(self, right, operator.sub)
        return PolynomialQuotientRingElement(self.parent(),
                               self.__polynomial - right.__polynomial, check=False)

    def charpoly(self):
        """
        The characteristic polynomial of this element, which is by
        definition the characteristic polynomial of right multiplication
        by this element.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3 -389*x^2 + 2*x - 5, 'a'); a = S.gen()
            sage: a.charpoly()
            x^3 - 389*x^2 + 2*x - 5
        """
        return self.matrix().charpoly()

    def fcp(self):
        """
        Return the factorization of the characteristic polynomial
        of this element.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3 -389*x^2 + 2*x - 5, 'a'); a = S.gen()
            sage: a.fcp()
            (x^3 - 389*x^2 + 2*x - 5)
            sage: S(1).fcp()
            (x - 1)^3
        """
        return self.charpoly().factor()

    def lift(self):
        """
        Return lift of this polynomial quotient ring element to the unique
        equivalent polynomial of degree less than the modulus.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3-2, 'a'); a = S.gen()
            sage: b = a^2 - 3
            sage: b
            a^2 - 3
            sage: b.lift()
            x^2 - 3
        """
        return self.__polynomial

    def list(self):
        """
        Return list of the elements of self, of length the same as
        the degree of the quotient polynomial ring.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3 + 2*x - 5, 'a'); a = S.gen()
            sage: a^10
            -134*a^2 - 35*a + 300
            sage: (a^10).list()
            [300, -35, -134]
        """
        v = self.__polynomial.list()
        R = self.parent()
        n = R.degree()
        return v + [R.base_ring()(0)]*(n - len(v))

    def matrix(self):
        """
        The matrix of right multiplication by this element on the
        power basis for the quotient ring.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3 + 2*x - 5, 'a'); a = S.gen()
            sage: a.matrix()
            [ 0  1  0]
            [ 0  0  1]
            [ 5 -2  0]
        """
        # Mutiply each power of field generator on the right by this
        # element, then return the matrix whose rows are the
        # coefficients of the result.
        try:
            return self.__matrix
        except AttributeError:
            R = self.parent()
            v = []
            x = R.gen()
            a = R(1)
            d = R.degree()
            for n in range(d):
                v += (a*self).list()
                a *= x
            S = R.base_ring()
            import sage.matrix.matrix_space
            M = sage.matrix.matrix_space.MatrixSpace(S, d)
            self.__matrix = M(v)
            return self.__matrix

    def minpoly(self):
        """
        The minimal polynomial of this element, which is by definition
        the minimal polynomial of right multiplication by this
        element.
        """
        return self.matrix().minpoly()

    def norm(self):
        """
        The norm of this element, which is the norm of the matrix
        of right multiplication by this element.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3 -389*x^2 + 2*x - 5, 'a'); a = S.gen()
            sage: a.norm()
            5
        """
        return self.matrix().determinant()

    def trace(self):
        """
        The trace of this element, which is the trace of the matrix
        of right multiplication by this element.

        EXAMPLES:
            sage: R, x = PolynomialRing(QQ).objgen()
            sage: S = R.quotient(x^3 -389*x^2 + 2*x - 5, 'a'); a = S.gen()
            sage: a.trace()
            389
        """
        return self.matrix().trace()


