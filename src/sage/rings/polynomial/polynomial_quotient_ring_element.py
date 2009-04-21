r"""
Elements of Quotients of Univariate Polynomial Rings

EXAMPLES: We create a quotient of a univariate polynomial ring over
`\ZZ`.

::

    sage: R.<x> = ZZ[]
    sage: S.<a> = R.quotient(x^3 + 3*x -1)
    sage: 2 * a^3
    -6*a + 2

Next we make a univeriate polynomial ring over
`\ZZ[x]/(x^3+3x-1)`.

::

    sage: S.<y> = S[]

And, we quotient out that by `y^2 + a`.

::

    sage: T.<z> = S.quotient(y^2+a)

In the quotient `z^2` is `-a`.

::

    sage: z^2
    -a

And since `a^3 = -3x + 1`, we have::

    sage: z^6
    3*a - 1

::

    sage: R.<x> = PolynomialRing(Integers(9))
    sage: S.<a> = R.quotient(x^4 + 2*x^3 + x + 2)
    sage: a^100
    7*a^3 + 8*a + 7

::

    sage: R.<x> = PolynomialRing(QQ)
    sage: S.<a> = R.quotient(x^3-2)
    sage: a
    a
    sage: a^3
    2

For the purposes of comparison in Sage the quotient element
`a^3` is equal to `x^3`. This is because when the
comparison is performed, the right element is coerced into the
parent of the left element, and `x^3` coerces to
`a^3`.

::

    sage: a == x
    True
    sage: a^3 == x^3
    True
    sage: x^3
    x^3
    sage: S(x^3)
    2

AUTHORS:

- William Stein
"""

####################################################################################
#       Copyright (C) 2005, 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
####################################################################################

import sage.rings.commutative_ring_element as commutative_ring_element
import sage.rings.number_field.number_field as number_field
import sage.rings.number_field.number_field_rel as number_field_rel


class PolynomialQuotientRingElement(commutative_ring_element.CommutativeRingElement):
    """
    Element of a quotient of a polynomial ring.
    """
    def __init__(self, parent, polynomial, check=True):
        """
        Create an element of the quotient of a polynomial ring.

        INPUT:


        -  ``parent`` - a quotient of a polynomial ring

        -  ``polynomial`` - a polynomial

        -  ``check`` - bool (optional): whether or not to
           verify that x is a valid element of the polynomial ring and reduced
           (mod the modulus).
        """
        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
        from sage.rings.polynomial.polynomial_element import Polynomial

        commutative_ring_element.CommutativeRingElement.__init__(self, parent)
        if check:
            if not isinstance(parent, PolynomialQuotientRing_generic):
                raise TypeError, "parent must be a polynomial quotient ring"

            if not isinstance(polynomial, Polynomial):
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
                P = B.parent()
                Q = P(0)
                X = P.gen()
                while R.degree() >= B.degree():
                    S = P((R.leading_coefficient()/B.leading_coefficient())) * X**(R.degree()-B.degree())
                    Q = Q + S
                    R = R - S*B
                polynomial = R
        self._polynomial = polynomial

    def _im_gens_(self, codomain, im_gens):
        return self._polynomial._im_gens_(codomain, im_gens)

    def __hash__(self):
        return self._polynomial.__hash__()

    def __reduce__(self):
        """
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S.<a> = R.quotient(2*x^3 + 3/2*x -1/3)
            sage: 2 * a^3
            -3/2*a + 1/3
            sage: loads(dumps(2*a^3)) == 2*a^3
            True
        """
        return PolynomialQuotientRingElement, (self.parent(), self._polynomial, False)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S.<a> = R.quotient(3*x^3 + 3/2*x -1/3)
            sage: 3 * a^3 + S.modulus()
            -3/2*a + 1/3
        """
        # We call _repr since _repr_ does not have a name variable.
        # This is very fragile!
        return self._polynomial._repr(self.parent().variable_name())

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S.<a> = R.quotient(3*x^3 + 3/2*x -1/3)
            sage: latex(a*(3 * a^3) + S.modulus())
            -\frac{3}{2} a^{2} + \frac{1}{3} a
        """
        return self._polynomial._latex_(self.parent().variable_name())

    ##################################################
    # Arithmetic
    ##################################################

    def _mul_(self, right):
        """
        Return the product of two polynomial ring quotient elements.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3-2)
            sage: (a^2 - 4) * (a+2)
            2*a^2 - 4*a - 6
        """
        R = self.parent()
        prod = self._polynomial * right._polynomial
        return PolynomialQuotientRingElement(R, prod, check=False)

    def _sub_(self, right):
        """
        Return the difference of two polynomial ring quotient elements.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 - 2)
            sage: (a^2 - 4) - (a+2)
            a^2 - a - 6
            sage: int(1) - a
            -a + 1
        """
        return PolynomialQuotientRingElement(self.parent(),
                                             self._polynomial - right._polynomial, check=False)

    def _add_(self, right):
        """
        Return the sum of two polynomial ring quotient elements.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3-2)
            sage: (a^2 - 4) + (a+2)
            a^2 + a - 2
            sage: int(1) + a
            a + 1
        """
        return PolynomialQuotientRingElement(self.parent(),
                                             self._polynomial + right._polynomial, check=False)

    def _div_(self, right):
        """
        Return the quotient of two polynomial ring quotient elements.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3-2)
            sage: (a^2 - 4) / (a+2)
            a - 2
        """
        return self * ~right

    def __neg__(self):
        return PolynomialQuotientRingElement(self.parent(), -self._polynomial)

    def __cmp__(self, other):
        """
        Compare this element with something else, where equality testing
        coerces the object on the right, if possible (and necessary).

        EXAMPLES:
        """
        return cmp(self._polynomial, other._polynomial)



    def __getitem__(self, n):
        return self._polynomial[n]

    def __int__(self):
        """
        Coerce this element to an int if possible.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3-2)
            sage: int(S(10))
            10
            sage: int(a)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to int
        """
        return int(self._polynomial)

    def __invert__(self):
        if self._polynomial.is_zero():
            raise ZeroDivisionError, \
               "element %s of quotient polynomial ring not invertible"%self
        g, _, a = self.parent().modulus().xgcd(self._polynomial)
        if g.degree() != 0:
            raise ZeroDivisionError, \
               "element %s of quotient polynomial ring not invertible"%self
        c = g[0]
        return PolynomialQuotientRingElement(self.parent(), (~c)*a, check=False)

    def __long__(self):
        """
        Coerce this element to a long if possible.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3-2)
            sage: long(S(10))
            10L
            sage: long(a)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to long
        """
        return long(self._polynomial)

    def field_extension(self, names):
        r"""
        Given a polynomial with base ring a quotient ring, return a
        3-tuple: a number field defined by the same polynomial, a
        homomorphism from its parent to the number field sending the
        generators to one another, and the inverse isomorphism.

        INPUT:

        - ``names`` - name of generator of output field


        OUTPUT:

        -  field

        -  homomorphism from self to field

        -  homomorphism from field to self


        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<alpha> = R.quotient(x^3-2)
            sage: F.<a>, f, g = alpha.field_extension()
            sage: F
            Number Field in a with defining polynomial x^3 - 2
            sage: a = F.gen()
            sage: f(alpha)
            a
            sage: g(a)
            alpha

        Over a finite field, the corresponding field extension is not a
        number field::

            sage: R.<x> = GF(25,'b')['x']
            sage: S.<a> = R.quo(x^3 + 2*x + 1)
            sage: F.<b>, g, h = a.field_extension()
            sage: h(b^2 + 3)
            a^2 + 3
            sage: g(x^2 + 2)
            b^2 + 2

        We do an example involving a relative number field::

            sage: R.<x> = QQ['x']
            sage: K.<a> = NumberField(x^3-2)
            sage: S.<X> = K['X']
            sage: Q.<b> = S.quo(X^3 + 2*X + 1)
            sage: F, g, h = b.field_extension('c')

        Another more awkward example::

            sage: R.<x> = QQ['x']
            sage: K.<a> = NumberField(x^3-2)
            sage: S.<X> = K['X']
            sage: f = (X+a)^3 + 2*(X+a) + 1
            sage: f
            X^3 + 3*a*X^2 + (3*a^2 + 2)*X + 2*a + 3
            sage: Q.<z> = S.quo(f)
            sage: F.<w>, g, h = z.field_extension()
            sage: c = g(z)
            sage: f(c)
            0
            sage: h(g(z))
            z
            sage: g(h(w))
            w

        AUTHORS:

        - Craig Citro (2006-08-06)

        - William Stein (2006-08-06)
        """
        #TODO: is the return order backwards from the magma convention?

##         We do another example over $\ZZ$.
##             sage: R.<x> = ZZ['x']
##             sage: S.<a> = R.quo(x^3 - 2)
##             sage: F.<b>, g, h = a.field_extension()
##             sage: h(b^2 + 3)
##             a^2 + 3
##             sage: g(x^2 + 2)
##             a^2 + 2
##         Note that the homomorphism is not defined on the entire
##         ''domain''.   (Allowing creation of such functions may be
##         disallowed in a future version of SAGE.):        <----- INDEED!
##             sage: h(1/3)
##             Traceback (most recent call last):
##             ...
##             TypeError: Unable to coerce rational (=1/3) to an Integer.
##         Note that the parent ring must be an integral domain:
##             sage: R.<x> = GF(25,'b')['x']
##             sage: S.<a> = R.quo(x^3 - 2)
##             sage: F, g, h = a.field_extension()
##             Traceback (most recent call last):
##             ...
##             ValueError: polynomial must be irreducible


        R = self.parent()
        x = R.gen()

        F = R.modulus().root_field(names)
        alpha = F.gen()

        f = R.hom([alpha], F, check=False)

        if number_field_rel.is_RelativeNumberField(F):

            base_hom = F.base_field().hom([R.base_ring().gen()])
            g = F.Hom(R)(x, base_hom)

        else:
            g = F.hom([x], R, check=False)

        return F, f, g


    def charpoly(self, var):
        """
        The characteristic polynomial of this element, which is by
        definition the characteristic polynomial of right multiplication by
        this element.

        INPUT:


        -  ``var`` - string - the variable name


        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quo(x^3 -389*x^2 + 2*x - 5)
            sage: a.charpoly('X')
            X^3 - 389*X^2 + 2*X - 5
        """
        return self.matrix().charpoly(var)

    def fcp(self, var='x'):
        """
        Return the factorization of the characteristic polynomial of this
        element.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 -389*x^2 + 2*x - 5)
            sage: a.fcp('x')
            x^3 - 389*x^2 + 2*x - 5
            sage: S(1).fcp('y')
            (y - 1)^3
        """
        return self.charpoly(var).factor()

    def lift(self):
        """
        Return lift of this polynomial quotient ring element to the unique
        equivalent polynomial of degree less than the modulus.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3-2)
            sage: b = a^2 - 3
            sage: b
            a^2 - 3
            sage: b.lift()
            x^2 - 3
        """
        return self._polynomial

    def __iter__(self):
        return iter(self.list())

    def list(self):
        """
        Return list of the elements of self, of length the same as the
        degree of the quotient polynomial ring.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 + 2*x - 5)
            sage: a^10
            -134*a^2 - 35*a + 300
            sage: (a^10).list()
            [300, -35, -134]
        """
        v = self._polynomial.list()
        R = self.parent()
        n = R.degree()
        return v + [R.base_ring()(0)]*(n - len(v))

    def matrix(self):
        """
        The matrix of right multiplication by this element on the power
        basis for the quotient ring.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 + 2*x - 5)
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
            for _ in xrange(d):
                v += (a*self).list()
                a *= x
            S = R.base_ring()
            import sage.matrix.matrix_space
            M = sage.matrix.matrix_space.MatrixSpace(S, d)
            self.__matrix = M(v)
            return self.__matrix

    def minpoly(self):
        """
        The minimal polynomial of this element, which is by definition the
        minimal polynomial of right multiplication by this element.
        """
        return self.matrix().minpoly()

    def norm(self):
        """
        The norm of this element, which is the norm of the matrix of right
        multiplication by this element.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 -389*x^2 + 2*x - 5)
            sage: a.norm()
            5
        """
        return self.matrix().determinant()

    def trace(self):
        """
        The trace of this element, which is the trace of the matrix of
        right multiplication by this element.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: S.<a> = R.quotient(x^3 -389*x^2 + 2*x - 5)
            sage: a.trace()
            389
        """
        return self.matrix().trace()


