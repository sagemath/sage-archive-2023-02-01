r"""
Lazy Laurent Series

A lazy laurent series is a laurent series whose coefficients are computed as
demanded or needed. Unlike the usual laurent series in Sage, lazy laurent
series do not have precisions because a lazy laurent series knows (can be
computed, lazily) all its coefficients.

EXAMPLES:

Generating functions are laurent series over the integer ring::

    sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
    sage: L = LazyLaurentSeriesRing(ZZ, 'z')

This defines the generating function of Fibonacci sequence::

    sage: def coeff(s, i):
    ....:     if i in [0, 1]:
    ....:         return 1
    ....:     else:
    ....:         return s.coefficient(i - 1) + s.coefficient(i - 2)
    ....:
    sage: f = L.series(coeff, valuation=0); f
    1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...

The 100th element of Fibonacci sequence can be obtained from the generating
function::

    sage: f.coefficient(100)
    573147844013817084101

Coefficients are computed and cached only when necessary::

    sage: f._cache[100]
    573147844013817084101
    sage: f._cache[101]
    Traceback (most recent call last):
    ...
    KeyError: 101

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

from sage.structure.element import Element

from sage.arith.power import generic_power

class LazyLaurentSeries(Element):
    r"""
    Return a lazy laurent series.

    INPUT:

    - ``coefficient`` -- Python function that computes coefficients

    - ``valuation`` -- integer; approximate valuation of the series

    - ``constant`` -- either ``None`` or pair of an element of the base ring and an integer

    Let the coefficient of index `i` mean the coefficient of the term of the
    series with exponent `i`.

    Python function ``coefficient`` returns the value of the coefficient of
    index `i` from input `s` and `i` where `s` is the series itself.

    Let ``valuation`` be `n`. All coefficients of index below `n` are zero.  If
    ``constant`` is ``None``, then the ``coefficient`` function is responsible
    to compute the values of all coefficients of index `\ge n`. If ``constant``
    is a pair `(c,m)`, then the ``coefficient`` function is responsible to
    compute the values of all coefficients of index `\ge n` and `< m` and all
    the coefficients of index `\ge m` is the constant `c`.

    EXAMPLES::

        sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
        sage: L = LazyLaurentSeriesRing(ZZ, 'z')
        sage: L.series(lambda s, i: i, valuation=-3, constant=(-1,3))
        -3*z^-3 - 2*z^-2 - z^-1 + z + 2*z^2 - z^3 - z^4 - z^5 + ...

        sage: def coeff(s, i):
        ....:     if i in [0, 1]:
        ....:         return 1
        ....:     else:
        ....:         return s.coefficient(i - 1) + s.coefficient(i - 2)
        ....:
        sage: f = L.series(coeff, valuation=0); f
        1 + z + 2*z^2 + 3*z^3 + 5*z^4 + 8*z^5 + 13*z^6 + ...
        sage: f.coefficient(100)
        573147844013817084101
    """
    def __init__(self, parent, coefficient=None, valuation=0, constant=None):
        """
        Initialize.

        TESTS::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: z = L.gen()
            sage: TestSuite(z).run(skip='_test_pickling')
        """
        Element.__init__(self, parent)

        self._coefficient_function = coefficient
        self._approximate_valuation = valuation
        self._constant = constant

        self._cache = dict() # cache of known coefficients

    def __eq__(self, other):
        """
        Return ``True`` if this series is known to be equal to ``other``.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: z = L.gen()
            sage: z == z^2 + z - z^2
            True
        """
        other = self.parent()(other)

        n = self._approximate_valuation

        if self._constant is None or other._constant is None:
            raise ValueError("undecidable as lazy laurent series")

        sc, sm = self._constant
        oc, om = other._constant

        if sc != oc:
            return False

        m = max(sm, om)

        for i in range(n, m):
            if self.coefficient(i) != other.coefficient(i):
                return False

        return True

    def __nonzero__(self):
        """
        Test whether ``self`` is not zero.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: z = L.gen()
            sage: (z-z).is_zero()
            True
        """
        if self._constant is None:
            raise ValueError("undecidable as lazy laurent series")

        sc, sm = self._constant

        if sc:
            return True

        for i in range(self._approximate_valuation, sm):
            if self.coefficient(i):
                return True

        return False

    def _repr_(self):
        """
        Return the string representation of this laurent series.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: -1/(1 + 2*z)
            -1 + 2*z - 4*z^2 + 8*z^3 - 16*z^4 + 32*z^5 - 64*z^6 + ...
        """
        atomic_repr = self.base_ring()._repr_option('element_is_atomic')
        X = self.parent().variable_name()

        try:
            n = self.valuation()
        except ValueError:
            n = self._approximate_valuation

        if self._constant is None:
            m = n + 7 # long enough
        elif self._constant[0] != 0:
            m = self._constant[1] + 3
        else:
            m = self._constant[1]

        s = ' '
        first = True
        while n < m:
            x = repr(self.coefficient(n))
            if x != '0':
                if not first:
                    s += ' + '
                if not atomic_repr and n > 0 and (x[1:].find('+') != -1 or x[1:].find('-') != -1):
                    x = '({})'.format(x)
                if n > 1 or n < 0:
                    var = '*%s^%s'%(X,n)
                elif n == 1:
                    var = '*%s'%X
                else:  # n == 0
                    var = ''
                s += '{}{}'.format(x,var)
                first = False
            n += 1

        s = s.replace(" + -", " - ").replace(" 1*"," ").replace(" -1*", " -")[1:]

        if len(s) == 0:
            s = '0'

        if self._constant is None or self._constant[1] > m or self._constant[0] != 0:
            s += ' + {}'.format('...')

        return s

    def coefficient(self, n):
        """
        Return the coefficient of the term of the sereis with exponent `n`.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: def g(s, i):
            ....:     if i == 0:
            ....:         return 1
            ....:     else:
            ....:         return sum(s.coefficient(j)*s.coefficient(i - 1 -j) for j in [0..i-1])
            ....:
            sage: e = L.series(g, valuation=0)
            sage: e.coefficient(10)
            16796
            sage: e
            1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + ...
        """
        R = self.base_ring()

        if self._approximate_valuation == infinity:
            return R.zero()
        elif n < self._approximate_valuation:
            return R.zero()
        elif self._constant is not None and n >= self._constant[1]:
            return self._constant[0]

        try:
            c = self._cache[n]
        except:
            c = R(self._coefficient_function(self, n))
            self._cache[n] = c

        return c

    def _mul_(self, other):
        """
        Return the product of this series with ``other``.

        INPUT:

        - ``other`` -- other series

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: (1 - z)*(1 - z)
            1 - 2*z + z^2
            sage: (1 - z)*(1 - z)*(1 - z)
            1 - 3*z + 3*z^2 - z^3
        """

        R = self.parent()

        def mul(s, n):
            c = self.base_ring().zero()
            for k in range(self._approximate_valuation, n - other._approximate_valuation + 1):
                c += self.coefficient(k) * other.coefficient(n-k)
            return c

        a = self._approximate_valuation + other._approximate_valuation

        c = None
        if self._constant is not None and other._constant is not None:
            if self._constant[0] == 0 and other._constant[0] == 0:
                c = (self._constant[0], self._constant[1] + other._constant[1] - 1)

        return R.element_class(R, coefficient=mul, valuation=a, constant=c)

    def _add_(self, other):
        """
        Return the sum of this series with ``other``.

        INPUT:

        - ``other`` -- other series

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: (1 - z)*(1 - z)
            1 - 2*z + z^2
            sage: (1 - z)*(1 - z)*(1 - z)
            1 - 3*z + 3*z^2 - z^3
        """
        R = self.parent()

        def add(s, n):
            return self.coefficient(n) + other.coefficient(n)

        a = min(self._approximate_valuation, other._approximate_valuation)

        if self._constant is not None and other._constant is not None:
            c = (self._constant[0] + other._constant[0],
                 max(self._constant[1], other._constant[1]))
        else:
            c = None

        return R.element_class(R, coefficient=add, valuation=a, constant=c)

    def _sub_(self, other):
        """
        Return the series of this series minus ``other`` series.

        INPUT:

        - ``other`` -- other series

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: z - z
            0
        """
        R = self.parent()

        def sub(s, n):
            return self.coefficient(n) - other.coefficient(n)

        a = min(self._approximate_valuation, other._approximate_valuation)

        if self._constant is not None and other._constant is not None:
            c = (self._constant[0] - other._constant[0],
                 max(self._constant[1], other._constant[1]))
        else:
            c = None

        return R.element_class(R, coefficient=sub, valuation=a, constant=c)

    def _neg_(self):
        """
        Return the negative of this series.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: -(1 - z)
            -1 + z
        """
        R = self.parent()

        def neg(s, i):
            return -self.coefficient(i)

        a = self._approximate_valuation

        if self._constant is not None:
            c = (-self._constant[0], self._constant[1])
        else:
            c = None

        return R.element_class(R, coefficient=neg, valuation=a, constant=c)

    def __invert__(self):
        """
        Return the multiplicative inverse of the element.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: ~(1 - z)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
        """
        v = self.valuation()

        if v == infinity:
            raise ZeroDivisionError('cannot invert zero')

        R = self.parent()

        a0 = self.coefficient(v)
        zero = self.base_ring().zero()

        def inv(s, n):
            if n == -v:
                return ~a0
            c = zero
            for k in range(-v, n):
                c += s.coefficient(k) * self.coefficient(n + v - k)
            return -c / a0

        return R.element_class(R, coefficient=inv, valuation=-v, constant=None)

    def __pow__(self, n):
        """
        Return the `n`-th power of the series.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: (1 - z)^-1
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: (1 - z)^0
            1
            sage: (1 - z)^3
            1 - 3*z + 3*z^2 - z^3
            sage: (1 - z)^-3
            1 + 3*z + 6*z^2 + 10*z^3 + 15*z^4 + 21*z^5 + 28*z^6 + ...
        """
        if n == 0:
            return self.parent()(1)

        return generic_power(self, n)

    def _div_(self, other):
        """
        Return the multiplicative inverse of the element.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: z/(1 - z)
            z + z^2 + z^3 + z^4 + z^5 + z^6 + z^7 + ...
        """
        return self * ~other

    def apply_to_coefficients(self, function):
        """
        Return the series with ``function`` applied to each coefficient of this series.

        INPUT:

        - ``function`` -- Python function

        Python function ``function`` returns a new coefficient for input coefficient.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: s = z/(1 - 2*z)
            sage: t = s.apply_to_coefficients(lambda c: c + 1)
            sage: s
            z + 2*z^2 + 4*z^3 + 8*z^4 + 16*z^5 + 32*z^6 + 64*z^7 + ...
            sage: t
            2*z + 3*z^2 + 5*z^3 + 9*z^4 + 17*z^5 + 33*z^6 + 65*z^7 + ...
        """
        R = self.parent()

        def apply_map(s, n):
            return function(self.coefficient(n))

        a = self._approximate_valuation

        if self._constant:
            c = (function(self._constant[0]), self._constant[1])
        else:
            c = None

        return R.element_class(R, coefficient=apply_map, valuation=a, constant=c)

    def change_ring(self, ring):
        """
        Return this series with coefficients converted to elements of ``ring``.

        INPUT:

        - ``ring`` -- a ring

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: s = 2 + z
            sage: t = s.change_ring(QQ)
            sage: t^-1
            1/2 - 1/4*z + 1/8*z^2 - 1/16*z^3 + 1/32*z^4 - 1/64*z^5 + 1/128*z^6 + ...
        """
        R = self.parent()

        def change(s, n):
            return ring(self.coefficient(n))

        a = self._approximate_valuation

        if self._constant:
            c = (ring(self._constant[0]), self._constant[1])
        else:
            c = None

        from .lazy_laurent_series_ring import LazyLaurentSeriesRing
        Q = LazyLaurentSeriesRing(ring, names=R.variable_name())
        return Q.element_class(Q, coefficient=change, valuation=a, constant=c)

    def valuation(self):
        """
        Return the valuation of the series.

        This method determines the valuation of the series by looking for a
        nonzero coefficient. Hence if the series happens to be zero, then it
        may run forever.

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: s = 1/(1 - z) - 1/(1 - 2*z)
            sage: s.valuation()
            1
            sage: t = z - z
            sage: t.valuation()
            +Infinity
        """
        if self._constant is None:
            n = self._approximate_valuation
            cache = self._cache
            while True:
                if n in cache:
                    if cache[n]:
                        self._approximate_valuation = n
                        return n
                    n += 1
                else:
                    if self.coefficient(n) != 0:
                        self._approximate_valuation = n
                        return n
                    n += 1
        else:
            n = self._approximate_valuation
            m = self._constant[1]
            while n <= m:
                if self.coefficient(n) != 0:
                    self._approximate_valuation = n
                    return n
                n += 1
            return infinity

    def truncate(self, d):
        """
        Return this series with its terms of degree > ``d`` truncated.

        INPUT:

        - ``d`` -- integer

        EXAMPLES::

            sage: from sage.rings.lazy_laurent_series_ring import LazyLaurentSeriesRing
            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: z = L.gen()
            sage: alpha = 1/(1-z)
            sage: alpha
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
            sage: beta = alpha.truncate(5)
            sage: beta
            1 + z + z^2 + z^3 + z^4 + z^5
            sage: alpha - beta
            z^6 + z^7 + z^8 + z^9 + z^10 + z^11 + z^12 + ...
        """
        R = self.parent()
        zero = self.base_ring().zero()

        def trunc(s, n):
            if n <= d:
                return self.coefficient(n)
            else:
                return zero

        return R.element_class(R, coefficient=trunc, valuation=self._approximate_valuation,
                               constant=(zero, d+1))
