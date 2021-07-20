r"""
Lazy Laurent Series

A lazy Laurent series is a Laurent series whose coefficients are computed as
demanded or needed. Unlike the usual Laurent series in Sage, lazy Laurent
series do not have precisions because a lazy Laurent series knows (can be
computed, lazily) all its coefficients.

EXAMPLES:

Generating functions are Laurent series over the integer ring::

    sage: L.<z> = LLSRing(ZZ)

This defines the generating function of Fibonacci sequence::

    sage: f = 1 / (1 - z - z^2)
    sage: f
    1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...

The 100th element of Fibonacci sequence can be obtained from the generating
function::

    sage: f.coefficient(100)
    573147844013817084101

Coefficients are computed and cached only when necessary::

    sage: f._aux._cache[100]
    573147844013817084101
    sage: f._aux._cache[101]
    Traceback (most recent call last):
    ...
    IndexError: list index out of range

You can do arithmetic with lazy power series::

    sage: f
    1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
    sage: f^-1
    1 - z - z^2 + ...
    sage: f + f^-1
    2 + z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
    sage: g = (f + f^-1)*(f - f^-1); g
    4*z + 6*z^2 + 8*z^3 + 19*z^4 + 38*z^5 + 71*z^6 + ...

You may need to change the base ring::

    sage: h = g.change_ring(QQ)
    sage: h.parent()
    Lazy Laurent Series Ring in z over Rational Field
    sage: h
    4*z + 6*z^2 + 8*z^3 + 19*z^4 + 38*z^5 + 71*z^6 + ...
    sage: h^-1
    1/4*z^-1 - 3/8 + 1/16*z - 17/32*z^2 + 5/64*z^3 - 29/128*z^4 + 165/256*z^5 + ...
    sage: _.valuation()
    -1

AUTHORS:

- Kwankyu Lee (2019-02-24): initial version

"""

# ****************************************************************************
#       Copyright (C) 2019 Kwankyu Lee <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from .infinity import infinity
from sage.structure.element import ModuleElement
from .integer_ring import ZZ
from sage.structure.richcmp import op_EQ, op_NE
from sage.arith.power import generic_power
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.integer import Integer

class LLS(ModuleElement):
    r"""
    A Laurent series where the coefficients are computed lazily.

    INPUT:

    - ``coefficient_function`` -- Python function that computes coefficients

    - ``issparse`` -- Boolean that determines whether the implementation is sparse or dense

    - ``approximate_valuation`` -- integer; approximate valuation of the series

    - ``constant`` -- either ``None`` or pair of an element of the base ring and an integer

    Let the coefficient of index `i` mean the coefficient of the term of the
    series with exponent `i`.

    Python function ``coefficient`` returns the value of the coefficient of
    index `i` from input.

    Let ``approximate_valuation`` be `n`. All coefficients of index below `n` are zero.  If
    ``constant`` is ``None``, then the ``coefficient`` function is responsible
    to compute the values of all coefficients of index `\ge n`. If ``constant``
    is a pair `(c,m)`, then the ``coefficient`` function is responsible to
    compute the values of all coefficients of index `\ge n` and `< m` and all
    the coefficients of index `\ge m` is the constant `c`.

    EXAMPLES::

        sage: L.<z> = LLSRing(ZZ)
        sage: L(lambda i: i, valuation=-3, constant=(-1,3))
        -3*z^-3 - 2*z^-2 - z^-1 + z + 2*z^2 - z^3 - z^4 - z^5 + ...
        sage: L(lambda i: i, valuation=-3, constant=-1, degree=3)
        -3*z^-3 - 2*z^-2 - z^-1 + z + 2*z^2 - z^3 - z^4 - z^5 + ...

    ::

        sage: f = 1 / (1 - z - z^2); f
        1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
        sage: f.coefficient(100)
        573147844013817084101

    Lazy Laurent series is picklable::

        sage: g = loads(dumps(f))
        sage: g
        1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
        sage: g == f
        True
    """
    def __init__(self, parent, aux):
        """
        Initialize the series.

        TESTS::

            sage: L = LLSRing(GF(2), 'z')
            sage: z = L.gen()
            sage: TestSuite(z).run()
        """
        ModuleElement.__init__(self, parent)
        self._aux = aux

    def __getitem__(self, n):
        """
        Return the coefficient of the term with exponent ``n`` of the series.

        INPUT:

        - ``n`` -- integer

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = z/(1 - 2*z^3)
            sage: [f[n] for n in range(20)]
            [0, 1, 0, 0, 2, 0, 0, 4, 0, 0, 8, 0, 0, 16, 0, 0, 32, 0, 0, 64]
            sage: f[0:20]
            [0, 1, 0, 0, 2, 0, 0, 4, 0, 0, 8, 0, 0, 16, 0, 0, 32, 0, 0, 64]

            sage: M = L(lambda n: n)
            sage: [M[n] for n in range(20)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n)
            sage: [M[n] for n in range(20)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        """
        R = self.base_ring()
        if isinstance(n, slice):
            if n.stop is None:
                raise NotImplementedError("cannot list an infinite set")
            start = n.start if n.start is not None else self._aux.valuation()
            step = n.step if n.step is not None else 1
            return [R(self._aux[k]) for k in range(start, n.stop, step)]
        return R(self._aux[n])


    def __call__(self, other):
        """
        Return the composition of the series with ``other``.

        INPUT:

        - ``other`` -- other series

        TESTS::

            sage: L.<z> = LLSRing(QQ)
            sage: LS.<y> = LLSRing(QQ)
            sage: f = z^2 + 1 + z
            sage: g = 0
            sage: f(g)
            1
            sage: g = L(0)
            sage: f(g)
            1
            sage: g = y - y                                                                             
            sage: f(g)                                                                                    
            1
            sage: f = L(lambda n: n, 0); f                                                                
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: f(g)                                                                                    
            0
            sage: L.<z> = LLSRing(QQ)                                                                                                                                                                     
            sage: LS.<y> = LLSRing(QQ, sparse=True)                                                                                                                                                       
            sage: f = L(lambda n : 1); f                                                                                                                                                                  
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: g = LS(lambda n: 0 if n != 2 else 1, 0); g                                                                                                                                              
            y^2 + ...
            sage: f(g)                                                                                                                                                                                    
            1 + z^2 + z^4 + z^6 + ...
            sage: g = LS(lambda n: 1 if n > 1 else 0, 0); g                                                                                                                                               
            y^2 + y^3 + y^4 + y^5 + y^6 + ...
            sage: f(g)                                                                                                                                                                                    
            1 + z^2 + z^3 + 2*z^4 + 3*z^5 + 5*z^6 + ...
            sage: L.<z> = LLSRing(QQ, sparse=True)                                                                                                                                                        
            sage: f = L(lambda n : 1); f                                                                                                                                                                  
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: g = LS(lambda n: 1, 1); g                                                                                                                                                               
            y + y^2 + y^3 + y^4 + y^5 + y^6 + y^7 + ...
            sage: f(g)                                                                                                                                                                                    
            1 + z + 2*z^2 + 4*z^3 + 8*z^4 + 16*z^5 + 32*z^6 + ...
        """
        # g = 0, val(f) >= 0
        try:
            if isinstance(other._aux, LLS_zero) and self.valuation() >= 0:
                return self[0]
        except AttributeError:
            if other == 0:
                return self[0]
            raise ValueError()
        
        # g != 0, val(g) > 0 and g != 0, f has finitely many non-zero coefficients
        P = self.parent()
        return P.element_class(P, LLS_com(self._aux, other._aux, P))

    def _mul_(self, other):
        """
        Return the product of this series with ``other``.

        INPUT:

        - ``other`` -- other series

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: (1 - z)*(1 - z)
            1 - 2*z + z^2
            sage: (1 - z)*(1 - z)*(1 - z)
            1 - 3*z + 3*z^2 - z^3
            sage: M = L(lambda n: n)
            sage: M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = M * (1 - M)
            sage: N
            z + z^2 - z^3 - 6*z^4 - 15*z^5 - 29*z^6 + ...
            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L(lambda n: 1); N
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: M * N
            z + 3*z^2 + 6*z^3 + 10*z^4 + 15*z^5 + 21*z^6 + ...
        """
        P = self.parent()
        left = self._aux
        right = other._aux
        if isinstance(left, LLS_zero) or isinstance(right, LLS_zero):
            return P.zero()

        if (isinstance(left, LLS_eventually_geometric) and isinstance(right, LLS_eventually_geometric)):
            if not left._constant and not right._constant:
                p = left._laurent_polynomial * right._laurent_polynomial
                c = left._constant
                return P.element_class(P, LLS_eventually_geometric(p, P._sparse, c))
        return P.element_class(P, LLS_mul(self._aux, other._aux))

    def _add_(self, other):
        """
        Return the sum of this series with ``other``.

        INPUT:

        - ``other`` -- other series

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: (1 - z)*(1 - z)
            1 - 2*z + z^2
            sage: (1 - z)*(1 - z)*(1 - z)
            1 - 3*z + 3*z^2 - z^3
            sage: z + z
            2*z
            sage: z^2 + 3*z^2
            4*z^2
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L(lambda n: 1); N
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M + N; P
            1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + ...

            sage: A = L(1, constant=2, degree=3)
            sage: B = L(2, constant=-2, degree=5)
            sage: A + B
            3 + 2*z^3 + 2*z^4

            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L(lambda n: 1); N
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M + N; P
            1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + ...
        """
        P = self.parent()
        left = self._aux
        right = other._aux
        if (isinstance(left, LLS_eventually_geometric)
            and isinstance(right, LLS_eventually_geometric)):
            R = P._laurent_poly_ring
            c = left._constant + right._constant
            pl = left._laurent_polynomial
            pr = right._laurent_polynomial
            d = max(left._degree, right._degree)
            pl += R([left._constant]*(d-left._degree)).shift(left._degree)
            pr += R([right._constant]*(d-right._degree)).shift(right._degree)
            p = pl + pr
            if not p and not c:
                return P.zero()
            return P.element_class(P, LLS_eventually_geometric(p, P._sparse, c, d))
        return P.element_class(P, LLS_add(self._aux, other._aux))

    def _sub_(self, other):
        """
        Return the series of this series minus ``other`` series.

        INPUT:

        - ``other`` -- other series

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: z - z
            0
            sage: 3*z - 2*z
            z
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L(lambda n: 1); N
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M - N; P
            -1 + z^2 + 2*z^3 + 3*z^4 + 4*z^5 + 5*z^6 + ...

            sage: A = L(1, constant=2, degree=3)
            sage: B = L(2, constant=3, degree=5)
            sage: A - B
            -1 + 2*z^3 + 2*z^4 - z^5 - z^6 - z^7 + ...

            sage: A = L(1, constant=2, degree=3)
            sage: B = L([1,0,0,2,2], constant=2)
            sage: X = A - B; X
            0
            sage: type(X._aux)
            <class 'sage.rings.lazy_laurent_series_new.LLS_zero'>

            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L(lambda n: 1); N
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M - N; P
            -1 + z^2 + 2*z^3 + 3*z^4 + 4*z^5 + 5*z^6 + ...
        """
        P = self.parent()
        left = self._aux
        right = other._aux
        if (isinstance(left, LLS_eventually_geometric) and isinstance(right, LLS_eventually_geometric)):
            R = P._laurent_poly_ring
            c = left._constant - right._constant
            pl = left._laurent_polynomial
            pr = right._laurent_polynomial
            d = max(left._degree, right._degree)
            pl += R([left._constant]*(d-left._degree)).shift(left._degree)
            pr += R([right._constant]*(d-right._degree)).shift(right._degree)
            p = pl - pr
            if not p and not c:
                return P.zero()
            return P.element_class(P, LLS_eventually_geometric(p, P._sparse, c, d))
        if left == right:
            return P.zero()
        return P.element_class(P, LLS_sub(self._aux, other._aux))

    def _div_(self, other):
        """
        Return ``self`` divided by ``other``.

        INPUT:

        - ``other`` -- nonzero series

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: z/(1 - z)
            z + z^2 + z^3 + z^4 + z^5 + z^6 + z^7 + ...
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L(lambda n: 1); N
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M / N; P
            z + z^2 + z^3 + z^4 + z^5 + z^6 + ...

            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = L(lambda n: 1); N
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: P = M / N; P
            z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
        """
        if isinstance(other._aux, LLS_zero):
            raise ZeroDivisionError("cannot divide by 0")

        P = self.parent()
        left = self._aux
        if isinstance(left, LLS_zero):
            return P.zero()
        right = other._aux
        if (isinstance(left, LLS_eventually_geometric)
            and isinstance(right, LLS_eventually_geometric)):
            if not left._constant and not right._constant:
                ret = left._laurent_polynomial / right._laurent_polynomial
                try:
                    ret = P._laurent_poly_ring(ret)
                    return P.element_class(P, LLS_eventually_geometric(ret, P._sparse, left._constant))
                except (TypeError, ValueError):
                    # We cannot divide the polynomials, so the result must be a series
                    pass

        return P.element_class(P, LLS_mul(left, LLS_inv(right)))
        # return P.element_class(P, LLS_div(left, right))

    def _rmul_(self, scalar):
        """
        Return the scalar multiplication of this series by ``scalar``.

        INPUT:

        - ``scalar`` -- an element of the base ring

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: 2*z
            2*z
            sage: -1*z
            -z
            sage: 0*z
            0
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M * 3
            3*z + 6*z^2 + 9*z^3 + 12*z^4 + 15*z^5 + 18*z^6 + ...

            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M * 3
            3*z + 6*z^2 + 9*z^3 + 12*z^4 + 15*z^5 + 18*z^6 + ...
            sage: N = L(lambda n: 1); N
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: N * 4
            4 + 4*z + 4*z^2 + 4*z^3 + 4*z^4 + 4*z^5 + 4*z^6 + ...
        """
        P = self.parent()
        if not scalar:
            return P.zero()

        if isinstance(self._aux, LLS_eventually_geometric):
            c = scalar * self._aux._constant
            p = scalar * self._aux._laurent_polynomial
            return P.element_class(P, LLS_eventually_geometric(p, P._sparse, c, self._aux._degree))

        return P.element_class(P, LLS_scalar(self._aux, scalar))

    def _neg_(self):
        """
        Return the negative of this series.

        TESTS::

            sage: L = LLSRing(ZZ, 'z')
            sage: z = L.gen()
            sage: -(1 - z)
            -1 + z
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: P = -M; P
            -z - 2*z^2 - 3*z^3 - 4*z^4 - 5*z^5 - 6*z^6 + ...

            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: P = -M; P
            -z - 2*z^2 - 3*z^3 - 4*z^4 - 5*z^5 - 6*z^6 + ...
            sage: -(z^2 + 3*z - 4*z^3)
            -3*z - z^2 + 4*z^3
        """
        P = self.parent()
        if isinstance(self._aux, LLS_eventually_geometric):
            p = -self._aux._laurent_polynomial
            c = -self._aux._constant
            d = self._aux._degree
            return P.element_class(P, LLS_eventually_geometric(p, P._sparse, c, d))
        return P.element_class(P, LLS_neg(self._aux))

    def __invert__(self):
        """
        Return the multiplicative inverse of the element.

        TESTS::

            sage: L.<z> = LLSRing(ZZ, sparse=False)
            sage: ~(1 - z)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: P = ~M; P
            z^-1 - 2 + z + ...
            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: P = ~M; P
            z^-1 - 2 + z + ...
        """
        P = self.parent()
        return P.element_class(P, LLS_inv(self._aux))

    def coefficient(self, n):
        """
        Return the coefficient of the term with exponent ``n`` of the series.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: L.<z> = LLSRing(ZZ)
            sage: F = L(None)
            sage: F.define(1 + z*F^2)
            sage: F.coefficient(10)
            16796
            sage: F
            1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + ...

        TESTS::

            sage: L.<z> = LLSRing(ZZ, sparse=False)
            sage: e = L(None)
            sage: e.define(1 + z*e^2)
            sage: e
            1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + ...
            sage: e._aux._cache
            [1, 1, 2, 5, 14, 42, 132]
            sage: e.coefficient(10)
            16796
            sage: e._aux._cache
            [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796]
            sage: M = L(lambda n: n^2); M
            z + 4*z^2 + 9*z^3 + 16*z^4 + 25*z^5 + 36*z^6 + ...
            sage: M._aux._cache
            [0, 1, 4, 9, 16, 25, 36]
            sage: M.coefficient(9)
            81
            sage: M._aux._cache
            [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

            sage: L = LLSRing(ZZ, 'z', sparse=True)
            sage: M = L(lambda n: n^2); M
            z + 4*z^2 + 9*z^3 + 16*z^4 + 25*z^5 + 36*z^6 + ...
            sage: M._aux._cache
            {0: 0, 1: 1, 2: 4, 3: 9, 4: 16, 5: 25, 6: 36}
            sage: M.coefficient(10)
            100
            sage: M._aux._cache
            {0: 0, 1: 1, 2: 4, 3: 9, 4: 16, 5: 25, 6: 36, 10: 100}
        """
        return self.__getitem__(n)

    def map_coefficients(self, func, ring=None):
        """
        Return the series with ``func`` applied to each coefficient of this series.

        INPUT:

        - ``func`` -- Python function that takes in a coefficient and returns
          a new coefficient

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: s = z/(1 - 2*z)
            sage: t = s.map_coefficients(lambda c: c + 1)
            sage: s
            z + 2*z^2 + 4*z^3 + 8*z^4 + 16*z^5 + 32*z^6 + 64*z^7 + ...
            sage: t
            2*z + 3*z^2 + 5*z^3 + 9*z^4 + 17*z^5 + 33*z^6 + 65*z^7 + ...
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = M.map_coefficients(lambda c: c + 1); N
            1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + ...

            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = M.map_coefficients(lambda c: c + 1); N
            1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + ...
        """
        P = self.parent()
        R = P.base_ring()
        if isinstance(self._aux, LLS_eventually_geometric):
            p = p.map_coefficients(func)
            c = func(c)
            if not p and not c:
                return P.zero()
            return P.element_class(P, LLS_eventually_geometric(p, self._aux._is_sparse, c, d))
        return P.element_class(P, LLS_apply_coeff(self._aux, func, R))

    def change_ring(self, ring):
        """
        Return this series with coefficients converted to elements of ``ring``.

        INPUT:

        - ``ring`` -- a ring

        TESTS::

            sage: L.<z> = LLSRing(ZZ, sparse=False)
            sage: s = 2 + z
            sage: t = s.change_ring(QQ)
            sage: t^-1
            1/2 - 1/4*z + 1/8*z^2 - 1/16*z^3 + 1/32*z^4 - 1/64*z^5 + 1/128*z^6 + ...
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: N = M.change_ring(QQ)
            sage: N.parent()
            Lazy Laurent Series Ring in z over Rational Field
            sage: M.parent()
            Lazy Laurent Series Ring in z over Integer Ring

            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M.parent()
            Lazy Laurent Series Ring in z over Integer Ring
            sage: N = M.change_ring(QQ)
            sage: N.parent()
            Lazy Laurent Series Ring in z over Rational Field
            sage: M ^-1
            z^-1 - 2 + z + ...
        """
        from .lazy_laurent_series_ring_new import LLSRing
        Q = LLSRing(ring, names=self.parent().variable_names())
        return Q.element_class(Q, self._aux)

    def truncate(self, d):
        """
        Return this series with its terms of degree >= ``d`` truncated.

        INPUT:

        - ``d`` -- integer

        TESTS::

            sage: L.<z> = LLSRing(ZZ, sparse=False)
            sage: alpha = 1/(1-z)
            sage: alpha
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: beta = alpha.truncate(5)
            sage: beta
            1 + z + z^2 + z^3 + z^4
            sage: alpha - beta
            z^5 + z^6 + ...
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M.truncate(4)
            z + 2*z^2 + 3*z^3

            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M.truncate(4)
            z + 2*z^2 + 3*z^3
        """
        P = self.parent()
        R = P._laurent_poly_ring
        z = R.gen()
        p = R.sum(self[i] * z**i for i in range(self._aux._approximate_valuation, d))
        return P.element_class(P, LLS_eventually_geometric(p, P._sparse, ZZ.zero(), d))

    def __pow__(self, n):
        """
        Return the ``n``-th power of the series.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: (1 - z)^-1
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: (1 - z)^0
            1
            sage: (1 - z)^3
            1 - 3*z + 3*z^2 - z^3
            sage: (1 - z)^-3
            1 + 3*z + 6*z^2 + 10*z^3 + 15*z^4 + 21*z^5 + 28*z^6 + ...
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M ^ 2
            z^2 + 4*z^3 + 10*z^4 + 20*z^5 + 35*z^6 + ...

            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: M = L(lambda n: n); M
            z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + ...
            sage: M ^ 2
            z^2 + 4*z^3 + 10*z^4 + 20*z^5 + 35*z^6 + ...
        """
        if n == 0:
            return self.parent().one()

        return generic_power(self, n)

    def approximate_series(self, prec, name=None):
        """
        Return the Laurent series with absolute precision ``prec`` approximated
        from this series.

        INPUT:

        - ``prec`` -- an integer

        - ``name`` -- name of the variable; if it is ``None``, the name of the variable
          of the series is used

        OUTPUT: a Laurent series with absolute precision ``prec``

        TESTS::

            sage: L = LLSRing(ZZ, 'z')
            sage: z = L.gen()
            sage: f = (z - 2*z^3)^5/(1 - 2*z)
            sage: f
            z^5 + 2*z^6 - 6*z^7 - 12*z^8 + 16*z^9 + 32*z^10 - 16*z^11 + ...
            sage: g = f.approximate_series(10)
            sage: g
            z^5 + 2*z^6 - 6*z^7 - 12*z^8 + 16*z^9 + O(z^10)
            sage: g.parent()
            Power Series Ring in z over Integer Ring
            sage: h = (f^-1).approximate_series(3)
            sage: h
            z^-5 - 2*z^-4 + 10*z^-3 - 20*z^-2 + 60*z^-1 - 120 + 280*z - 560*z^2 + O(z^3)
            sage: h.parent()
            Laurent Series Ring in z over Integer Ring
        """
        S = self.parent()

        if name is None:
            name = S.variable_name()

        if self.valuation() < 0:
            from sage.rings.all import LaurentSeriesRing
            R = LaurentSeriesRing(S.base_ring(), name=name)
            n = self.valuation()
            return R([self[i] for i in range(n, prec)], n).add_bigoh(prec)
        else:
            from sage.rings.all import PowerSeriesRing
            R = PowerSeriesRing(S.base_ring(), name=name)
            return R([self[i] for i in range(prec)]).add_bigoh(prec)

    def prec(self):
        """
        Return the precision of the series, which is infinity.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = 1/(1 - z)
            sage: f.prec()
            +Infinity
        """
        return infinity

    def polynomial(self, degree=None, name=None):
        """
        Return the polynomial or Laurent polynomial if the series is actually so.

        INPUT:

        - ``degree`` -- ``None`` or an integer

        - ``name`` -- name of the variable; if it is ``None``, the name of the variable
          of the series is used

        OUTPUT: a Laurent polynomial if the valuation of the series is negative or
        a polynomial otherwise.

        If ``degree`` is not ``None``, the terms of the series of degree
        greater than ``degree`` are truncated first. If ``degree`` is ``None``
        and the series is not a polynomial or a Laurent polynomial, a
        ``ValueError`` is raised.

        EXAMPLES::

            sage: L = LLSRing(ZZ, 'z')
            sage: f = L([1,0,0,2,0,0,0,3], 5); f
            z^5 + 2*z^8 + 3*z^12
            sage: f.polynomial()
            3*z^12 + 2*z^8 + z^5

        TESTS::

            sage: g = L([1,0,0,2,0,0,0,3], -5); g
            z^-5 + 2*z^-2 + 3*z^2
            sage: g.polynomial()
            z^-5 + 2*z^-2 + 3*z^2
            sage: z = L.gen()
            sage: f = (1 + z)/(z^3 - z^5)
            sage: f
            z^-3 + z^-2 + z^-1 + 1 + z + z^2 + z^3 + ...
            sage: f.polynomial(5)
            z^-3 + z^-2 + z^-1 + 1 + z + z^2 + z^3 + z^4 + z^5
            sage: f.polynomial(0)
            z^-3 + z^-2 + z^-1 + 1
            sage: f.polynomial(-5)
            0
            sage: M = L(lambda n: n^2, 0)
            sage: M.polynomial(3)
            9*z^3 + 4*z^2 + z
            sage: M = L(lambda n: n^2, 0)
            sage: M.polynomial(5)
            25*z^5 + 16*z^4 + 9*z^3 + 4*z^2 + z
        """
        if degree is None:
            if isinstance(self._aux, LLS_zero):
                from sage.rings.all import PolynomialRing
                return PolynomialRing(S.base_ring(), name=name).zero()
            elif isinstance(self._aux, LLS_eventually_geometric) and self._aux._constant:
                m = self._aux._degree
            else:
                raise ValueError("not a polynomial")
        else:
            m = degree + 1

        S = self.parent()

        if name is None:
            name = S.variable_name()

        if self.valuation() < 0:
            from sage.rings.all import LaurentPolynomialRing
            R = LaurentPolynomialRing(S.base_ring(), name=name)
            n = self.valuation()
            return R([self[i] for i in range(n, m)]).shift(n)
        else:
            from sage.rings.all import PolynomialRing
            R = PolynomialRing(S.base_ring(), name=name)
            return R([self[i] for i in range(m)])

    def valuation(self):
        """
        Return the valuation of the series.

        This method determines the valuation of the series by looking for a
        nonzero coefficient. Hence if the series happens to be zero, then it
        may run forever.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: s = 1/(1 - z) - 1/(1 - 2*z)
            sage: s.valuation()
            1
            sage: t = z - z
            sage: t.valuation()
            +Infinity
            sage: M = L(lambda n: n^2, 0)
            sage: M.valuation()
            1
            sage: M = L(lambda n: n^2, 0)
            sage: M.valuation()
            1
        """
        return self._aux.valuation()

    def _repr_(self):
        """
        Return the string representation of this Laurent series.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: -1/(1 + 2*z)
            -1 + 2*z - 4*z^2 + 8*z^3 - 16*z^4 + 32*z^5 - 64*z^6 + ...
        """
        if isinstance(self._aux, LLS_zero):
            return '0'
        if isinstance(self._aux, LLS_uninitialized) and self._aux._target is None:
            return 'Uninitialized LLS'

        atomic_repr = self.base_ring()._repr_option('element_is_atomic')
        X = self.parent().variable_name()
        v = self._aux._approximate_valuation

        if not isinstance(self._aux, LLS_eventually_geometric):
            m = v + 7  # long enough
        elif not self._aux._constant:
            # Just a polynonial, so let that print itself
            return repr(self._aux._laurent_polynomial)
        else:
            m = self._aux._degree + 3

        # Use the polynomial printing
        R = self.parent()._laurent_poly_ring
        ret = repr(R([self._aux[i] for i in range(v, m)]).shift(v))
        # TODO: Better handling when ret == 0 but we have not checked up to the constant term
        return ret + ' + ...'

    def _richcmp_(self, other, op):
        """
        Compare ``self` with ``other`` with respect to the comparison operator ``op``.

        Equality is verified if the corresponding coefficients of both series
        can be checked for equality without computing coefficients
        indefinitely.  Otherwise an exception is raised to declare that
        equality is not decidable.

        Inequality is not defined for lazy Laurent series.

        TESTS::

            sage: L.<z> = LLSRing(QQ)
            sage: z + z^2 == z^2 + z
            True
            sage: z + z^2 != z^2 + z
            False
            sage: z + z^2 > z^2 + z
            False
            sage: z + z^2 < z^2 + z
            False
        """
        if op is op_EQ:
            if isinstance(self._aux, LLS_zero):  # self == 0
                return isinstance(other._aux, LLS_zero)
            if isinstance(other._aux, LLS_zero):  # self != 0 but other == 0
                return False

            if (not isinstance(self._aux, LLS_eventually_geometric)
                or not isinstance(other._aux, LLS_eventually_geometric)):
                # One of the lazy laurent series is not known to eventually be constant
                # Implement the checking of the caches here.
                n = min(self._aux._approximate_valuation, other._aux._approximate_valuation)
                m = max(self._aux._approximate_valuation, other._aux._approximate_valuation)
                for i in range(n, m):
                    if self[i] != other[i]:
                        return False
                if self._aux == other._aux:
                    return True
                raise ValueError("undecidable as lazy Laurent series")

            # Both are LLS_eventually_geometric, which implements a full check
            return self._aux == other._aux

        if op is op_NE:
            return not (self == other)

        return False

    def __hash__(self):
        """
        Return the hash of ``self``

        TESTS::

            sage: L = LLSRing(ZZ, 'z')
            sage: f = L([1,2,3,4], -5)
            sage: g = (1 + f)/(1 - f)^2
            sage: {g: 1}
            {z^5 - 2*z^6 + z^7 + 5*z^9 - 11*z^10 + z^11 + ...: 1}
        """
        return hash(self._aux)

    def __bool__(self):
        """
        Test whether ``self`` is not zero.

        TESTS::

            sage: L.<z> = LLSRing(GF(2))
            sage: (z-z).is_zero()
            True
            sage: f = 1/(1 - z)
            sage: f.is_zero()
            False
            sage: M = L(lambda n: n, 0); M
            z + z^3 + z^5 + ...
            sage: M.is_zero()
            False
        """
        if isinstance(self._aux, LLS_zero):
            return False
        if isinstance(self._aux, LLS_eventually_geometric):
            # This should always end up being True, but let's be careful about it for now...
            return self._aux._laurent_polynomial or self._aux._constant

        for a in self._aux._cache:
            if a:
                return True
        if self[self._aux._approximate_valuation]:
            return True
        raise ValueError("undecidable as lazy Laurent series")

    def define(self, s):
        r"""
        Define an equation by ``self = s``.

        EXAMPLES:

        We begin by constructing the Catalan numbers::

            sage: L.<z> = LLSRing(ZZ)
            sage: C = L(None)
            sage: C.define(1 + z*C^2)
            sage: C
            1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + ...

        The Catalan numbers but with a valuation 1::

            sage: B = L(None, 1)
            sage: B.define(z + B^2)
            sage: B
            z + z^2 + 2*z^3 + 5*z^4 + 14*z^5 + 42*z^6 + 132*z^7 + ...

        We can define multiple series that are linked::

            sage: s = L(None)
            sage: t = L(None)
            sage: s.define(1 + z*t^3)
            sage: t.define(1 + z*s^2)
            sage: s[:9]
            [1, 1, 3, 9, 34, 132, 546, 2327, 10191]
            sage: t[:9]
            [1, 1, 2, 7, 24, 95, 386, 1641, 7150]

        An bigger example::

            sage: L.<z> = LLSRing(ZZ)
            sage: A = L(None, 5)
            sage: B = L(None)
            sage: C = L(None, 2)
            sage: A.define(z^5 + B^2)
            sage: B.define(z^5 + C^2)
            sage: C.define(z^2 + C^2 + A^2)
            sage: A[0:15]
            [0, 0, 0, 0, 0, 1, 0, 0, 1, 2, 5, 4, 14, 10, 48]
            sage: B[0:15]
            [0, 0, 0, 0, 1, 1, 2, 0, 5, 0, 14, 0, 44, 0, 138]
            sage: C[0:15]
            [0, 0, 1, 0, 1, 0, 2, 0, 5, 0, 15, 0, 44, 2, 142]

        We count unlabeled ordered trees by total number of nodes
        and number of internal nodes::

            sage: R.<q> = QQ[]
            sage: Q.<z> = LLSRing(R)
            sage: leaf = z
            sage: internal_node = q * z
            sage: L = Q(constant=1, degree=1)
            sage: T = Q(None, 1)
            sage: T.define(leaf + internal_node * L(T))  # not tested - composition need
            sage: [T[i] for i in range(7)]  # not tested - composition need
            [0, 1, q, q^2 + q, q^3 + 3*q^2 + q, q^4 + 6*q^3 + 6*q^2 + q]

        TESTS::

            sage: L.<z> = LLSRing(ZZ, sparse=True)
            sage: s = L(None)
            sage: s.define(1 + z*s^3)
            sage: s[:10]
            [1, 1, 3, 12, 55, 273, 1428, 7752, 43263, 246675]

            sage: e = L(None)
            sage: e.define(1 + z*e)
            sage: e.define(1 + z*e)
            Traceback (most recent call last):
            ...
            ValueError: series already defined
            sage: z.define(1 + z^2)
            Traceback (most recent call last):
            ...
            ValueError: series already defined
        """
        if not isinstance(self._aux, LLS_uninitialized) or self._aux._target is not None:
            raise ValueError("series already defined")
        self._aux._target = s._aux

class LLS_aux():
    """
    Abstract base class for all auxillary LLS.
    """
    def __init__(self, sparse, approximate_valuation):
        self._is_sparse = sparse
        self._approximate_valuation = approximate_valuation

class LLS_inexact(LLS_aux):
    """
    LLS aux class when it is not or we do not know if it is
    eventually geometric.
    """
    def __init__(self, is_sparse, approximate_valuation):
        super().__init__(is_sparse, approximate_valuation)

        if self._is_sparse:
            self._cache = dict()  # cache of known coefficients
        else:
            self._cache = list()
            self._offset = approximate_valuation
            self._iter = self.iterate_coefficients()

    def __getstate__(self):
        d = dict(self.__dict__)
        if not self._is_sparse:
            # We cannot pickle a generator object, so we remove it and
            #   the cache from the pickle information.
            del d["_iter"]
            del d["_cache"]
        return d

    def __setstate__(self, d):
        self.__dict__ = d
        if not self._is_sparse:
            self._iter = self.iterate_coefficients()
            self._cache = []

    def __getitem__(self, n):
        if n < self._approximate_valuation:
            return ZZ.zero()

        if self._is_sparse:
            try:
                c = self._cache[n]
            except KeyError:
                c = self.get_coefficient(n)
                self._cache[n] = c
        else:
            i = n - self._offset
            if i >= len(self._cache):
                a = len(self._cache) + self._offset
                # it is important to extend by generator:
                # self._coefficient_function might recurse, and
                # thereby extend the cache itself, too
                self._cache.extend(next(self._iter) for _ in range(a, n+1))
            c = self._cache[i]

        return c

    def valuation(self):
        if self._is_sparse:
            n = self._approximate_valuation
            cache = self._cache
            while True:
                if n in cache:
                    if cache[n]:
                        self._approximate_valuation = n
                        return n
                    n += 1
                else:
                    if self[n] != 0:
                        self._approximate_valuation = n
                        return n
                    n += 1
        else:
            n = self._approximate_valuation
            cache = self._cache
            while True:
                if n - self._offset < len(cache):
                    if cache[n - self._offset]:
                        self._approximate_valuation = n
                        return n
                    n += 1
                else:
                    if self[n] != 0:
                        self._approximate_valuation = n
                        return n
                    n += 1


class LLS_unary(LLS_inexact):
    """
    Abstract base class for unary operators.

    INPUT:

    - ``series`` -- series upon which the operator operates

    """
    def __init__(self, series, *args, **kwargs):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = -1/(1 - z)
            sage: f
            -1 - z - z^2 - z^3 - z^4 - z^5 - z^6 + ...
            sage: loads(dumps(f)) == f
            True
        """
        self._series = series
        super().__init__(*args, **kwargs)

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = ~(1 - z)
            sage: {f: 1}
            {1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...: 1}
        """
        return hash((type(self), self._series))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: g = 1/(1 - z) + 1/(1 + z)
            sage: f == g
            True
            sage: f = ~(1 - z)
            sage: g = ~(1 - z)
            sage: f == g
            True
        """
        return isinstance(other, type(self)) and self._series == other._series


class LLS_binary(LLS_inexact):

    def __init__(self, left, right, *args, **kwargs):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: loads(dumps(f)) == f
            True
            sage: f = 1/(1 - z) - 1/(1 + z)
            sage: loads(dumps(f)) == f
            True
        """
        self._left = left
        self._right = right
        super().__init__(*args, **kwargs)

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: {f: 1}
            {2 + 2*z^2 + 2*z^4 + 2*z^6 + ...: 1}
        """
        return hash((type(self), self._left, self._right))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: g = 1/(1 - z) + 1/(1 + z)
            sage: f == g
            True
        """
        if not isinstance(other, type(self)):
            return False
        return self._left == other._left and self._right == other._right


class LLS_binary_commutative(LLS_binary):

    def __hash__(self):
        return hash((type(self), frozenset([self._left, self._right])))

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False
        if self._left == other._left and self._right == other._right:
            return True
        if self._left == other._right and self._right == other._left:
            return True
        return False


class LLS_eventually_geometric(LLS_aux):
    def __init__(self, laurent_polynomial, is_sparse, constant=None, degree=None):
        if constant is None:
            constant = ZZ.zero()
        if degree is None:
            if not laurent_polynomial:
                raise ValueError("you must specify the degree for the polynomial 0")
            degree = laurent_polynomial.degree() + 1
        else:
            # Consistency check
            assert not laurent_polynomial or laurent_polynomial.degree() < degree

        self._constant = constant
        self._degree = degree
        self._laurent_polynomial = laurent_polynomial
        if not laurent_polynomial:
            valuation = degree
        else:
            valuation = laurent_polynomial.valuation()
        super().__init__(is_sparse, valuation)

    def __getitem__(self, n):
        if n >= self._degree:
            return self._constant
        return self._laurent_polynomial[n]

    def valuation(self):
        return self._approximate_valuation

    def __hash__(self):
        return hash((self._laurent_polynomial, self._degree, self._constant))

    def __eq__(self, other):
        return (isinstance(other, type(self))
                and self._degree == other._degree
                and self._laurent_polynomial == other._laurent_polynomial
                and self._constant == other._constant)

class LLS_zero(LLS_aux):
    def __init__(self, sparse):
        return super().__init__(sparse, 0)

    def __getitem__(self, n):
        return ZZ.zero()

    def valuation(self):
        return infinity

    def __hash__(self):
        return 0

class LLS_coefficient_function(LLS_inexact):
    def __init__(self, coefficient_function, ring, is_sparse, approximate_valuation):
        self._coefficient_function = coefficient_function
        self._ring = ring
        super().__init__(is_sparse, approximate_valuation)

    def get_coefficient(self, n):
        return self._ring(self._coefficient_function(n))

    def iterate_coefficients(self):
        n = self._offset
        ring = self._ring
        while True:
            yield ring(self._coefficient_function(n))
            n += 1

class LLS_uninitialized(LLS_inexact):
    def __init__(self, is_sparse, approximate_valuation):
        self._target = None
        super().__init__(is_sparse, approximate_valuation)

    def get_coefficient(self, n):
        return self._target[n]

    def iterate_coefficients(self):
        n = self._approximate_valuation
        while True:
            yield self._target[n]
            n += 1

#####################################################################
## Binary operations

class LLS_add(LLS_binary):
    """
    Operator for addition.
    """
    def __init__(self, left, right):
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = min(left._approximate_valuation, right._approximate_valuation)
        super().__init__(left, right, left._is_sparse, a)

    def get_coefficient(self, n):
        return self._left[n] + self._right[n]

    def iterate_coefficients(self):
        n = self._offset
        while True:
            yield self._left[n] + self._right[n]
            n += 1

class LLS_sub(LLS_binary):
    """
    Operator for subtraction.
    """
    def __init__(self, left, right):
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = min(left._approximate_valuation, right._approximate_valuation)
        super().__init__(left, right, left._is_sparse, a)

    def get_coefficient(self, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L.<z> = LLSRing(ZZ)
            sage: f = (1 + z)*(1 - z)
            sage: f.coefficient(2)
            -1
        """
        return self._left[n] - self._right[n]

    def iterate_coefficients(self):
        n = self._offset
        while True:
            yield self._left[n] - self._right[n]
            n += 1

class LLS_mul(LLS_binary):
    """
    Operator for multiplication.

    We are assuming commutativity of the coefficient ring here.
    """
    def __init__(self, left, right):
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = left._approximate_valuation + right._approximate_valuation
        super().__init__(left, right, left._is_sparse, a)

    def get_coefficient(self, n):
        c = ZZ.zero()
        for k in range(self._left._approximate_valuation,
                       n - self._right._approximate_valuation + 1):
            l = self._left[k]
            if l:
                c += l * self._right[n-k]
        return c

    def iterate_coefficients(self):
        n = self._offset
        while True:
            c = ZZ.zero()
            for k in range(self._left._approximate_valuation,
                           n - self._right._approximate_valuation + 1):
                l = self._left[k]
                if l:
                    c += l * self._right[n-k]
            yield c
            n += 1

class LLS_div(LLS_binary):
    """
    Return ``left`` divided by ``right``.
    """
    def __init__(self, left, right):
        lv = left.valuation()
        rv = right.valuation()
        self._lv = lv
        self._rv = rv
        self._ainv = ~right[rv]
        super().__init__(left, right, left._is_sparse, lv - rv)

    def get_coefficient(self, n):
        lv = self._lv
        rv = self._rv
        if n == lv - rv:
            return self._left[lv] / self._right[rv]
        c = self._left[n + rv]
        for k in range(lv - rv, n):
            c -= self[k] * self._right[n + rv - k]
        return c * self._ainv

    def iterate_coefficients(self):
        n = self._offset
        lv = self._lv
        rv = self._rv
        while True:
            if n == lv - rv:
                yield self._left[lv] / self._right[rv]
                n += 1
                continue
            c = self._left[n + rv]
            for k in range(lv - rv, n):
                c -= self[k] * self._right[n + rv - k]
            yield c * self._ainv
            n += 1


class LLS_com(LLS_binary):
    """
    Return ``left`` composed by ``right``.
    """
    def __init__(self, left, right, ring):
        lv = left.valuation()
        rv = right.valuation()
        self._ring = ring
        self._lv = lv
        self._rv = rv
        self._ainv = ~right[rv]
        super().__init__(left, right, left._is_sparse, lv * rv)
    
    def get_coefficient(self, n):
        # g != 0, val(g) > 0
        if self._rv > 0:
            if n == 0:
                return self._left[0]
            tail = self._ring(lambda n: self._left[n + 1], 0)
            right_series = self._ring(lambda n: self._right[n], self._rv)
            z = tail(right_series) * right_series
            z.coefficient(1)
            return z[n]
    
    def iterate_coefficients(self):
        # g != 0, val(g) > 0
        if self._rv > 0:
            yield self._left[0]
            tail = self._ring(lambda n: self._left[n + 1], 0)
            right_series = self._ring(lambda n: self._right[n], self._rv)
            z = tail(right_series) * right_series
            z.coefficient(1)
            n = 1
            while True:
                yield z[n]
                n += 1


#####################################################################
## Unary operations

class LLS_scalar(LLS_unary):
    """
    Operator for multiplying with a scalar.
    """
    def __init__(self, series, scalar):
        self._scalar = scalar

        super().__init__(series, series._is_sparse, series._approximate_valuation)

    def get_coefficient(self, n):
        return self._series[n] * self._scalar

    def iterate_coefficients(self):
        n = self._offset
        while True:
            yield self._series[n] * self._scalar
            n += 1

class LLS_neg(LLS_unary):
    """
    Operator for negative of the series.
    """
    def __init__(self, series):
        super().__init__(series, series._is_sparse, series._approximate_valuation)

    def get_coefficient(self, n):
        return -self._series[n]

    def iterate_coefficients(self):
        n = self._offset
        while True:
            yield -self._series[n]
            n += 1

class LLS_inv(LLS_unary):
    """
    Operator for multiplicative inverse of the series.
    """
    def __init__(self, series):
        v = series.valuation()
        super().__init__(series, series._is_sparse, -v)

        self._ainv = ~series[v]
        self._zero = ZZ.zero()

    def get_coefficient(self, n):
        v = self._approximate_valuation
        if n == v:
            return self._ainv
        c = self._zero
        for k in range(v, n):
            c += self[k] * self._series[n - v - k]
        return -c * self._ainv

    def iterate_coefficients(self):
        n = self._offset
        while True:
            v = self._approximate_valuation
            if n == v:
                yield self._ainv
                n += 1
                continue
            c = self._zero
            for k in range(v, n):
                c += self[k] * self._series[n - v - k]
            yield -c * self._ainv
            n += 1

class LLS_apply_coeff(LLS_unary):
    """
    Return the series with ``function`` applied to each coefficient of this series.
    """
    def __init__(self, series, function, ring):
        self._function = function
        self._ring = ring
        super().__init__(series, series._is_sparse, series._approximate_valuation)

    def get_coefficient(self, n):
        c = self._ring(self._function(self._series[n]))
        return c

    def iterate_coefficients(self):
        n = self._offset
        while True:
            c = self._ring(self._function(self._series[n]))
            yield c
            n += 1
