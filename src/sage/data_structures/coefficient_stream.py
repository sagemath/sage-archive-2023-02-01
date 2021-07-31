r"""
Coefficient Stream

This module provides lazy class implementations of basic operators
on coefficient streams. The classes implemented in this module
can be used to build up more complex streams for different kinds of
series (Laurent, Dirichlet, etc).

EXAMPLES:

The coefficient stream can be used to build up a Lazy laurent series::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ)
    sage: f = L(lambda n: n, True)
    sage: f
    z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + 7*z^7 + ...
    sage: type(f._coeff_stream)
    <class 'sage.data_structures.coefficient_stream.CoefficientStream_coefficient_function'>

    There are basic unary and binary operators available for the coefficient streams.
    For example, we can add two streams together::

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
    sage: from sage.data_structures.coefficient_stream import CoefficientStream_add
    sage: f = CoefficientStream_coefficient_function(lambda n: n, QQ, True, 0)
    sage: [f[i] for i in range(10)]
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    sage: g = CoefficientStream_coefficient_function(lambda n: 1, QQ, True, 0)
    sage: [g[i] for i in range(10)]
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    sage: h = CoefficientStream_add(f, g)
    sage: [h[i] for i in range(10)]
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

Coefficient streams can be subtracted::

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_sub
    sage: h = CoefficientStream_sub(f, g)
    sage: [h[i] for i in range(10)]
    [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]

Coefficient streams can be multiplied::

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_product
    sage: h = CoefficientStream_cauchy_product(f, g)
    sage: [h[i] for i in range(10)]
    [0, 1, 3, 6, 10, 15, 21, 28, 36, 45]

Coefficient streams can be divided::

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_div
    sage: h = CoefficientStream_div(f, g)
    sage: [h[i] for i in range(10)]
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1]

Two coefficient streams can be composed (depending on whether it exists)::

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_composition
    sage: g = CoefficientStream_coefficient_function(lambda n: n, QQ, True, 1)
    sage: h = CoefficientStream_composition(f, g)
    sage: [h[i] for i in range(10)]
    [0, 1, 4, 14, 46, 145, 444, 1331, 3926, 11434]

We can also use the unary negation operator on a coefficient stream::

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_neg
    sage: h = CoefficientStream_neg(f)
    sage: [h[i] for i in range(10)]
    [0, -1, -2, -3, -4, -5, -6, -7, -8, -9]

Coefficient streams can be multiplied by a scalar::

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_scalar
    sage: h = CoefficientStream_scalar(f, 2)
    sage: [h[i] for i in range(10)]
    [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]

The multiplicative inverse of a series can also be obtained::

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_inverse
    sage: h = CoefficientStream_cauchy_inverse(g)
    sage: [h[i] for i in range(10)]
    [-2, 1, 0, 0, 0, 0, 0, 0, 0, 0]

Functions can also be applied to a coefficient stream::

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_apply_coeff
    sage: h = CoefficientStream_apply_coeff(f, lambda n: n^2, QQ)
    sage: [h[i] for i in range(10)]
    [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

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

from sage.rings.integer_ring import ZZ
from sage.rings.infinity import infinity
from sage.arith.misc import divisors

class CoefficientStream():
    """
    Abstract base class for all coefficient streams.

    INPUT:

    - ``sparse`` -- boolean; whether the implementation of the series is sparse
    - ``approximate_valuation`` -- integer; a lower bound for the valuation of the series
    """
    def __init__(self, sparse, approximate_valuation):
        """
        Initialize the auxillary class for any series.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream
            sage: type(CoefficientStream(True, 1))
            <class 'sage.data_structures.coefficient_stream.CoefficientStream'>
        """
        self._is_sparse = sparse
        self._approximate_valuation = approximate_valuation


class CoefficientStream_inexact(CoefficientStream):
    """
    An abstract base class for the stream when we do not know it is
    eventually geometric.

    INPUT:

    - ``sparse`` -- boolean; whether the implementation of the series is sparse
    - ``approximate_valuation`` -- integer; a lower bound for the valuation of the series
    """
    def __init__(self, is_sparse, approximate_valuation):
        """
        Initialize the stream class for a CoefficientStream when it is not
        or it cannot be determined if it is eventually geometric.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_inexact
            sage: f = CoefficientStream_inexact(True, 0)
            sage: f._is_sparse
            True
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
            sage: g = CoefficientStream_coefficient_function(lambda n: n, QQ, False, 0)
            sage: g._offset
            0
        """
        super().__init__(is_sparse, approximate_valuation)

        if self._is_sparse:
            self._cache = dict()  # cache of known coefficients
        else:
            self._cache = list()
            self._offset = approximate_valuation
            self._iter = self.iterate_coefficients()

    def __getstate__(self):
        """
        Remove the cache from the pickle information so that it can be pickled.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: h = CoefficientStream_exact([1], True)
            sage: g = CoefficientStream_exact([1, -1, -1], True)
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_div
            sage: u = CoefficientStream_div(h, g)
            sage: [u[i] for i in range(10)]
            [1, 1, 2, 3, 5, 8, 13, 21, 34, 55]
            sage: m = loads(dumps(u))
            sage: [m[i] for i in range(10)]
            [1, 1, 2, 3, 5, 8, 13, 21, 34, 55]
        """
        d = dict(self.__dict__)
        if not self._is_sparse:
            # We cannot pickle a generator object, so we remove it and
            #   the cache from the pickle information.
            del d["_iter"]
            del d["_cache"]
        return d

    def __setstate__(self, d):
        """
        Re-create the cache and the generator object when unpickling.

        INPUT:

        - ``d`` -- a dictionary that needs to be unpickled

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: h = CoefficientStream_exact([-1], True)
            sage: g = CoefficientStream_exact([1, -1], True)
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_div
            sage: u = CoefficientStream_div(h, g)
            sage: [u[i] for i in range(10)]
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
            sage: loads(dumps(u)) == u
            True
        """
        self.__dict__ = d
        if not self._is_sparse:
            self._iter = self.iterate_coefficients()
            self._cache = []

    def __getitem__(self, n):
        """
        Return the `n`-th coefficient of the series.

        INPUT:

        - ``n`` -- integer; the index of the coefficient to return

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
            sage: f = CoefficientStream_coefficient_function(lambda n: n^2, QQ, True, 0)
            sage: f[3]
            9
            sage: f._cache
            {3: 9}
            sage: [f[i] for i in range(10)]
            [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
            sage: f._cache
            {0: 0, 1: 1, 2: 4, 3: 9, 4: 16, 5: 25, 6: 36, 7: 49, 8: 64, 9: 81}

            sage: f = CoefficientStream_coefficient_function(lambda n: n^2, QQ, False, 0)
            sage: f[3]
            9
            sage: f._cache
            [0, 1, 4, 9]
            sage: [f[i] for i in range(10)]
            [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
            sage: f._cache
            [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
        """
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
        """
        Return the valuation of the series.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
            sage: f = CoefficientStream_coefficient_function(lambda n: n, QQ, True, 0)
            sage: f.valuation()
            1
        """
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


class CoefficientStream_exact(CoefficientStream):
    r"""
    A stream of eventually constant coefficients.

    INPUT:

        - ``initial_values`` -- a list of initial values

        - ``is_sparse`` -- a boolean, which specifies whether the
          series is sparse

        - ``valuation`` -- (default: 0), an integer, determining the
          degree of the first element of ``initial_values``

        - ``degree`` -- (default: None), an integer, determining the
          degree of the first element which is known to be equal to
          ``constant``

        - ``constant`` -- (default: 0), an integer, the coefficient
          of every index larger than or equal to ``degree``
    """
    def __init__(self, initial_coefficients, is_sparse, constant=None, degree=None, valuation=None):
        """
        Initialize a series that is known to be eventually geometric.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: CoefficientStream_exact([], False)
            Traceback (most recent call last):
            ...
            AssertionError: CoefficientStream_exact should only be used for non-zero streams

        """
        if constant is None:
            self._constant = ZZ.zero()
        else:
            self._constant = constant
        if valuation is None:
            valuation = 0
        if degree is None:
            self._degree = valuation + len(initial_coefficients)
        else:
            self._degree = degree

        assert valuation + len(initial_coefficients) <= self._degree

        for i, v in enumerate(initial_coefficients):
            if v:
                valuation += i
                initial_coefficients = initial_coefficients[i:]
                for j, w in enumerate(reversed(initial_coefficients)):
                    if w:
                        break
                    initial_coefficients.pop()
                self._initial_coefficients = tuple(initial_coefficients)
                break
        else:
            valuation = self._degree
            self._initial_coefficients = tuple()

        assert self._initial_coefficients or self._constant, "CoefficientStream_exact should only be used for non-zero streams"

        super().__init__(is_sparse, valuation)

    def __getitem__(self, n):
        """
        Return the coefficient of the term with exponent ``n`` of the series.

        INPUT:

        - ``n`` -- integer, the degree for which the coefficient is required

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: s = CoefficientStream_exact([1], False)
            sage: [s[i] for i in range(-2, 5)]
            [0, 0, 1, 0, 0, 0, 0]

            sage: s = CoefficientStream_exact([], False, constant=1)
            sage: [s[i] for i in range(-2, 5)]
            [0, 0, 1, 1, 1, 1, 1]

            sage: s = CoefficientStream_exact([2], False, constant=1)
            sage: [s[i] for i in range(-2, 5)]
            [0, 0, 2, 1, 1, 1, 1]

            sage: s = CoefficientStream_exact([2], False, valuation=-1, constant=1)
            sage: [s[i] for i in range(-2, 5)]
            [0, 2, 1, 1, 1, 1, 1]

            sage: s = CoefficientStream_exact([2], False, valuation=-1, degree=2, constant=1)
            sage: [s[i] for i in range(-2, 5)]
            [0, 2, 0, 0, 1, 1, 1]

            sage: t = CoefficientStream_exact([0, 2, 0], False, valuation=-2, degree=2, constant=1)
            sage: t == s
            True
        """
        if n >= self._degree:
            return self._constant
        i = n - self._approximate_valuation
        if i < 0 or i >= len(self._initial_coefficients):
            return 0
        return self._initial_coefficients[i]

    def valuation(self):
        return self._approximate_valuation

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: f = 1 + z + z^2 + z^3
            sage: {f: 1}
            {1 + z + z^2 + z^3: 1}
        """
        return hash((self._initial_coefficients, self._degree, self._constant))

    def __eq__(self, other):
        """
        Test the equality between ``self`` and ``other``.

        INPUT:

        - ``other`` -- a lazy Laurent series which is known to be eventaully geometric

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: f = 1 + z + z^2 + z^3
            sage: m = 1 + z + z^2 + z^3
            sage: f == m
            True
        """
        return (isinstance(other, type(self))
                and self._degree == other._degree
                and self._initial_coefficients == other._initial_coefficients
                and self._constant == other._constant)


class CoefficientStream_coefficient_function(CoefficientStream_inexact):
    r"""
    Class that returns the elements in the coefficient stream.

    INPUT:

    - ``coefficient_function`` -- a python function that generates the
      coefficients of the series
    - ``ring`` -- the base ring of the series
    - ``is_sparse`` -- boolean; specifies whether the series is sparse
    - ``approximate_valuation`` -- integer; a lower bound for the valuation of the series

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
        sage: f = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, False, 1)
        sage: f[3]
        9
        sage: [f[i] for i in range(10)]
        [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
    """

    def __init__(self, coefficient_function, ring, is_sparse, approximate_valuation):
        """
        Initialize.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, False, 1)
        """
        self._coefficient_function = coefficient_function
        self._ring = ring
        super().__init__(is_sparse, approximate_valuation)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
            sage: f = CoefficientStream_coefficient_function(lambda n: n, QQ, True, 0)
            sage: f.get_coefficient(4)
            4
        """
        return self._ring(self._coefficient_function(n))

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, QQ, False, 0)
            sage: n = f.iterate_coefficients()
            sage: [next(n) for _ in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        n = self._offset
        ring = self._ring
        while True:
            yield ring(self._coefficient_function(n))
            n += 1


class CoefficientStream_uninitialized(CoefficientStream_inexact):
    r"""
    Coefficient stream for an uninitialized series.

    INPUT:

    - ``is_sparse`` -- boolean; which specifies whether the series is sparse
    - ``approximate_valuation`` -- integer; a lower bound for the valuation of the series

    EXAMPLES::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ)
        sage: N = L(None)
        sage: N
        Uninitialized Lazy Laurent Series
    """
    def __init__(self, is_sparse, approximate_valuation):
        """
        Initialize an uninitialized lazy laurent series.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: N = L(None, 2)
            sage: N
            Uninitialized Lazy Laurent Series
        """
        self._target = None
        super().__init__(is_sparse, approximate_valuation)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
            sage: C = L(None)
            sage: C.define(1 + z*C^2)
            sage: C._coeff_stream.get_coefficient(4)
            14

        """
        return self._target[n]

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: C = L(None)
            sage: C.define(1 + z*C^2)
            sage: n = C._coeff_stream.iterate_coefficients()
            sage: [next(n) for _ in range(6)]
            [1, 1, 2, 5, 14, 42]

        """
        n = self._approximate_valuation
        while True:
            yield self._target[n]
            n += 1


class CoefficientStream_unary(CoefficientStream_inexact):
    r"""
    Class for unary operators for the coefficient stream.

    INPUT:

    - ``series`` -- :class:`CoefficientStream` the operator acts on

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_cauchy_inverse, CoefficientStream_scalar)
        sage: f = CoefficientStream_coefficient_function(lambda n: 2*n, ZZ, False, 1)
        sage: g = CoefficientStream_cauchy_inverse(f)
        sage: [g[i] for i in range(10)]
        [-1, 1/2, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: g = CoefficientStream_scalar(f, 2)
        sage: [g[i] for i in range(10)]
        [0, 4, 8, 12, 16, 20, 24, 28, 32, 36]
    """

    def __init__(self, series, *args, **kwargs):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_neg
            sage: CoefficientStream_neg.__base__
            <class 'sage.data_structures.coefficient_stream.CoefficientStream_unary'>
        """
        self._series = series
        super().__init__(*args, **kwargs)

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 - z)
            sage: {f: 1}
            {1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...: 1}
        """
        return hash((type(self), self._series))

    def __eq__(self, other):
        """
        Test equality.

        INPUT:

        - ``other`` -- a stream of coefficients

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_scalar)
            sage: f = CoefficientStream_coefficient_function(lambda n: 2*n, ZZ, False, 1)
            sage: g = CoefficientStream_coefficient_function(lambda n: n, ZZ, False, 1)
            sage: h = CoefficientStream_scalar(f, 2)
            sage: n = CoefficientStream_scalar(g, 2)
            sage: h == n
            False
            sage: n == n
            True
            sage: h == h
            True
        """
        return isinstance(other, type(self)) and self._series == other._series


class CoefficientStream_binary(CoefficientStream_inexact):
    """
    Class for binary operators for the coefficient stream.

    INPUT:

    - ``left`` -- :class:`CoefficientStream` for the left side of the operator
    - ``right`` -- :class:`CoefficientStream` for the right side of the operator

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_add, CoefficientStream_sub)
        sage: f = CoefficientStream_coefficient_function(lambda n: 2*n, ZZ, True, 0)
        sage: g = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
        sage: h = CoefficientStream_add(f, g)
        sage: [h[i] for i in range(10)]
        [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
        sage: h = CoefficientStream_sub(f, g)
        sage: [h[i] for i in range(10)]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    """

    def __init__(self, left, right, *args, **kwargs):
        """
        Initialize.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_product
            sage: (CoefficientStream_cauchy_product.__base__).__base__
            <class 'sage.data_structures.coefficient_stream.CoefficientStream_binary'>
        """
        self._left = left
        self._right = right
        super().__init__(*args, **kwargs)

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: {f: 1}
            {2 + 2*z^2 + 2*z^4 + 2*z^6 + ...: 1}
        """
        return hash((type(self), self._left, self._right))

    def __eq__(self, other):
        """
        Test equality.

        INPUT:

        - ``other`` -- a stream of coefficients

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_cauchy_product)
            sage: f = CoefficientStream_coefficient_function(lambda n: 2*n, ZZ, False, 1)
            sage: g = CoefficientStream_coefficient_function(lambda n: n, ZZ, False, 1)
            sage: h = CoefficientStream_coefficient_function(lambda n: 1, ZZ, False, 1)
            sage: t = CoefficientStream_cauchy_product(f, g)
            sage: u = CoefficientStream_cauchy_product(g, h)
            sage: v = CoefficientStream_cauchy_product(h, f)
            sage: t == u
            False
            sage: t == t
            True
            sage: u == v
            False
        """
        if not isinstance(other, type(self)):
            return False
        return self._left == other._left and self._right == other._right


class CoefficientStream_binary_commutative(CoefficientStream_binary):
    r"""
    Abstract base class for commutative binary operators for the
    coefficient stream.

    INPUT:

    - ``other`` -- a :class:`CoefficientStream`

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_add)
        sage: f = CoefficientStream_coefficient_function(lambda n: 2*n, ZZ, True, 0)
        sage: g = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
        sage: h = CoefficientStream_add(f, g)
        sage: [h[i] for i in range(10)]
        [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
        sage: u = CoefficientStream_add(g, f)
        sage: [u[i] for i in range(10)]
        [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
        sage: h == u
        True
    """
    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = z^2 + z
            sage: {f: 1}
            {z + z^2: 1}
        """
        return hash((type(self), frozenset([self._left, self._right])))

    def __eq__(self, other):
        """
        Test the equality between ``self`` and ``other``.

        INPUT:

        - ``other`` -- a stream of coefficients

        TESTS::
            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_cauchy_product)
            sage: f = CoefficientStream_coefficient_function(lambda n: 2*n, ZZ, True, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
            sage: h = CoefficientStream_cauchy_product(f, g)
            sage: [h[i] for i in range(10)]
            [0, 0, 2, 8, 20, 40, 70, 112, 168, 240]
            sage: u = CoefficientStream_cauchy_product(g, f)
            sage: [u[i] for i in range(10)]
            [0, 0, 2, 8, 20, 40, 70, 112, 168, 240]
            sage: h == u
            True
        """
        if not isinstance(other, type(self)):
            return False
        if self._left == other._left and self._right == other._right:
            return True
        if self._left == other._right and self._right == other._left:
            return True
        return False


class CoefficientStream_zero(CoefficientStream):
    """
    A coefficient Stream that is exactly equal to zero.

    INPUT:

    - ``sparse`` -- boolean; whether the coefficient stream is sparse or not

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import CoefficientStream_zero
        sage: s = CoefficientStream_zero(True)
        sage: s[5]
        0
    """

    def __init__(self, sparse):
        """
        Initialize.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_zero
            sage: s = CoefficientStream_zero(False)
            sage: [s[i] for i in range(10)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        """
        return super().__init__(sparse, 0)

    def __getitem__(self, n):
        """
        Return the ``n``-th coefficient of the series.

        INPUT::

        - ``n`` -- integer; the index of the coefficient to be returned

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_zero
            sage: s = CoefficientStream_zero(True)
            sage: s[1]
            0
            sage: sum([s[i] for i in range(10)])
            0
        """
        return ZZ.zero()

    def valuation(self):
        """
        Return the valuation of the series.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_zero
            sage: s = CoefficientStream_zero(True)
            sage: s.valuation()
            +Infinity
        """
        return infinity

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_zero
            sage: s = CoefficientStream_zero(False)
            sage: a = s.__hash__(); a
            0
            sage: t = CoefficientStream_zero(False)
            sage: b = t.__hash__(); b
            0
            sage: b == a
            True
        """
        return 0


#####################################################################
# Binary operations

class CoefficientStream_add(CoefficientStream_binary_commutative):
    """
    Operator for addition of two coefficient streams.

    INPUT:

    - ``left`` -- stream of coefficients on the left side of the operator
    - ``right`` -- stream of coefficients on the right side of the operator

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_add, CoefficientStream_coefficient_function)
        sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 0)
        sage: g = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 0)
        sage: h = CoefficientStream_add(f, g)
        sage: [h[i] for i in range(10)]
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        sage: u = CoefficientStream_add(g, f)
        sage: [u[i] for i in range(10)]
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    """
    def __init__(self, left, right):
        """
        Initalize.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_add)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 0)
            sage: h = CoefficientStream_add(f, g)
        """
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = min(left._approximate_valuation, right._approximate_valuation)
        super().__init__(left, right, left._is_sparse, a)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_add)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 0)
            sage: h = CoefficientStream_add(f, g)
            sage: h.get_coefficient(5)
            30
            sage: [h.get_coefficient(i) for i in range(10)]
            [0, 2, 6, 12, 20, 30, 42, 56, 72, 90]
        """
        return self._left[n] + self._right[n]

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_add)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, False, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^3, ZZ, False, 0)
            sage: h = CoefficientStream_add(f, g)
            sage: n = h.iterate_coefficients()
            sage: [next(n) for i in range(10)]
            [1, 2, 9, 28, 65, 126, 217, 344, 513, 730]
        """
        n = self._offset
        while True:
            yield self._left[n] + self._right[n]
            n += 1


class CoefficientStream_sub(CoefficientStream_binary):
    """
    Operator for subtraction of two coefficient streams.

    INPUT:

    - ``left`` -- stream of coefficients on the left side of the operator
    - ``right`` -- stream of coefficients on the right side of the operator

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_sub, CoefficientStream_coefficient_function)
        sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 0)
        sage: g = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 0)
        sage: h = CoefficientStream_sub(f, g)
        sage: [h[i] for i in range(10)]
        [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]
        sage: u = CoefficientStream_sub(g, f)
        sage: [u[i] for i in range(10)]
        [1, 0, -1, -2, -3, -4, -5, -6, -7, -8]
    """

    def __init__(self, left, right):
        """
        Initalize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_sub)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 0)
            sage: h = CoefficientStream_sub(f, g)
        """
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = min(left._approximate_valuation, right._approximate_valuation)
        super().__init__(left, right, left._is_sparse, a)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_sub)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 0)
            sage: h = CoefficientStream_sub(f, g)
            sage: h.get_coefficient(5)
            -20
            sage: [h.get_coefficient(i) for i in range(10)]
            [0, 0, -2, -6, -12, -20, -30, -42, -56, -72]
        """
        return self._left[n] - self._right[n]

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_sub)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, False, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^3, ZZ, False, 0)
            sage: h = CoefficientStream_sub(f, g)
            sage: n = h.iterate_coefficients()
            sage: [next(n) for i in range(10)]
            [1, 0, -7, -26, -63, -124, -215, -342, -511, -728]
        """
        n = self._offset
        while True:
            yield self._left[n] - self._right[n]
            n += 1


class CoefficientStream_cauchy_product(CoefficientStream_binary_commutative):
    """
    Operator for multiplication of two coefficient streams.

    We are assuming commutativity of the coefficient ring here.

    INPUT:

    - ``left`` -- stream of coefficients on the left side of the operator
    - ``right`` -- stream of coefficients on the right side of the operator

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_cauchy_product, CoefficientStream_coefficient_function)
        sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 0)
        sage: g = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 0)
        sage: h = CoefficientStream_cauchy_product(f, g)
        sage: [h[i] for i in range(10)]
        [0, 1, 3, 6, 10, 15, 21, 28, 36, 45]
        sage: u = CoefficientStream_cauchy_product(g, f)
        sage: [u[i] for i in range(10)]
        [0, 1, 3, 6, 10, 15, 21, 28, 36, 45]
    """
    def __init__(self, left, right):
        """
        Initalize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_cauchy_product)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 0)
            sage: h = CoefficientStream_cauchy_product(f, g)
        """
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = left._approximate_valuation + right._approximate_valuation
        super().__init__(left, right, left._is_sparse, a)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_cauchy_product)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 0)
            sage: h = CoefficientStream_cauchy_product(f, g)
            sage: h.get_coefficient(5)
            50
            sage: [h.get_coefficient(i) for i in range(10)]
            [0, 0, 1, 6, 20, 50, 105, 196, 336, 540]
        """
        c = ZZ.zero()
        for k in range(self._left._approximate_valuation,
                       n - self._right._approximate_valuation + 1):
            val = self._left[k]
            if val:
                c += val * self._right[n-k]
        return c

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_cauchy_product)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, False, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^3, ZZ, False, 0)
            sage: h = CoefficientStream_cauchy_product(f, g)
            sage: n = h.iterate_coefficients()
            sage: [next(n) for i in range(10)]
            [0, 1, 9, 36, 100, 225, 441, 784, 1296, 2025]
        """
        n = self._offset
        while True:
            c = ZZ.zero()
            for k in range(self._left._approximate_valuation,
                           n - self._right._approximate_valuation + 1):
                val = self._left[k]
                if val:
                    c += val * self._right[n-k]
            yield c
            n += 1

class CoefficientStream_dirichlet_convolution(CoefficientStream_binary_commutative):
    """Operator for the convolution of two coefficient streams.

    We are assuming commutativity of the coefficient ring here.
    Moreover, the valuation must be non-negative.

    INPUT:

    - ``left`` -- stream of coefficients on the left side of the operator
    - ``right`` -- stream of coefficients on the right side of the operator

    The coefficient of `n^{-s}` in the convolution of `l` and `r`
    equals `\sum_{k | n} l_k r_{n/k}`.  Note that `l[n]` yields the
    coefficient of `(n+1)^{-s}`!

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_dirichlet_convolution, CoefficientStream_coefficient_function, CoefficientStream_exact)
        sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
        sage: g = CoefficientStream_exact([0], True, constant=1)
        sage: h = CoefficientStream_dirichlet_convolution(f, g)
        sage: [h[i] for i in range(1, 10)]
        [1, 3, 4, 7, 6, 12, 8, 15, 13]
        sage: [sigma(n) for n in range(1, 10)]
        [1, 3, 4, 7, 6, 12, 8, 15, 13]

        sage: u = CoefficientStream_dirichlet_convolution(g, f)
        sage: [u[i] for i in range(1, 10)]
        [1, 3, 4, 7, 6, 12, 8, 15, 13]

    """
    def __init__(self, left, right):
        """
        Initalize ``self``.
        """
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        assert left._approximate_valuation > 0 and right._approximate_valuation > 0, "Dirichlet convolution is only defined for coefficient streams with valuation at least 1"

        vl = left._approximate_valuation
        vr = right._approximate_valuation
        a = vl * vr
        super().__init__(left, right, left._is_sparse, a)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        """
        c = ZZ.zero()
        for k in divisors(n):
            val = self._left[k]
            if val:
                c += val * self._right[n//k]
        return c

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.
        """
        n = self._offset
        while True:
            c = ZZ.zero()
            for k in divisors(n):
                val = self._left[k]
                if val:
                    c += val * self._right[n//k]
            yield c
            n += 1

class CoefficientStream_dirichlet_inv(CoefficientStream_unary):
    """
    Operator for multiplicative inverse of the stream.

    INPUT:

    - ``series`` -- a :class:`CoefficientStream`

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_dirichlet_inv, CoefficientStream_coefficient_function)
        sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
        sage: g = CoefficientStream_dirichlet_inv(f)
        sage: [g[i] for i in range(10)]
        [0, 1, -1, -1, 0, -1, 1, -1, 0, 0]
        sage: [moebius(i) for i in range(10)]
        [0, 1, -1, -1, 0, -1, 1, -1, 0, 0]
    """
    def __init__(self, series):
        """
        Initialize.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_exact, CoefficientStream_dirichlet_inv)
            sage: f = CoefficientStream_exact([0, 0], True, constant=1)
            sage: g = CoefficientStream_dirichlet_inv(f)
            Traceback (most recent call last):
            ...
            AssertionError: the Dirichlet inverse only exists if the coefficient with index 1 is non-zero
        """
        assert series[1], "the Dirichlet inverse only exists if the coefficient with index 1 is non-zero"
        super().__init__(series, series._is_sparse, 1)

        self._ainv = ~series[1]
        self._zero = ZZ.zero()

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient
        """
        if n == 1:
            return self._ainv
        c = self._zero
        for k in divisors(n):
            if k < n:
                val = self._series[n//k]
                if val:
                    c += self[k] * val
        return -c * self._ainv

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.
        """
        n = 1
        yield self._ainv
        while True:
            n += 1
            c = self._zero
            for k in divisors(n):
                if k < n:
                    val = self._series[n//k]
                    if val:
                        c += self[k] * val
            yield -c * self._ainv



class CoefficientStream_div(CoefficientStream_binary):
    """
    Operator for division of two coefficient streams.

    INPUT:

    - ``left`` -- stream of coefficients on the left side of the operator
    - ``right`` -- stream of coefficients on the right side of the operator

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_div, CoefficientStream_coefficient_function)
        sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
        sage: g = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
        sage: h = CoefficientStream_div(f, g)
        sage: [h[i] for i in range(10)]
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        sage: u = CoefficientStream_div(g, f)
        sage: [u[i] for i in range(10)]
        [1, -1, 0, 0, 0, 0, 0, 0, 0, 0]
    """

    def __init__(self, left, right):
        """
        Initalize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_div)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 1)
            sage: h = CoefficientStream_div(f, g)
        """
        lv = left.valuation()
        rv = right.valuation()
        self._lv = lv
        self._rv = rv
        self._ainv = ~right[rv]
        super().__init__(left, right, left._is_sparse, lv - rv)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_div)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 1)
            sage: h = CoefficientStream_div(f, g)
            sage: h.get_coefficient(5)
            -2
            sage: [h.get_coefficient(i) for i in range(10)]
            [1, -2, 2, -2, 2, -2, 2, -2, 2, -2]
        """
        lv = self._lv
        rv = self._rv
        if n == lv - rv:
            return self._left[lv] / self._right[rv]
        c = self._left[n + rv]
        for k in range(lv - rv, n):
            c -= self[k] * self._right[n + rv - k]
        return c * self._ainv

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_div)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, False, 1)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^3, ZZ, False, 1)
            sage: h = CoefficientStream_div(f, g)
            sage: n = h.iterate_coefficients()
            sage: [next(n) for i in range(10)]
            [1, -7, 30, -114, 426, -1590, 5934, -22146, 82650, -308454]
        """
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


class CoefficientStream_composition(CoefficientStream_binary):
    r"""
    Return ``f`` composed by ``g``.

    This is the composition `(f \circ g)(z) = f(g(z))`.

    INPUT:

    - ``f`` -- a :class:`CoefficientStream`
    - ``g`` -- a :class:`CoefficientStream` with positive valuation

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_composition, CoefficientStream_coefficient_function)
        sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
        sage: g = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
        sage: h = CoefficientStream_composition(f, g)
        sage: [h[i] for i in range(10)]
        [0, 1, 3, 8, 20, 48, 112, 256, 576, 1280]
        sage: u = CoefficientStream_composition(g, f)
        sage: [u[i] for i in range(10)]
        [0, 1, 3, 8, 21, 55, 144, 377, 987, 2584]
    """
    def __init__(self, f, g):
        """
        Initalize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_composition)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 1)
            sage: h = CoefficientStream_composition(f, g)
        """
        assert g._approximate_valuation > 0
        self._fv = f._approximate_valuation
        self._gv = g._approximate_valuation
        if self._fv < 0:
            ginv = CoefficientStream_cauchy_inverse(g)
            # the constant part makes no contribution to the negative
            # we need this for the case so self._neg_powers[0][n] => 0
            self._neg_powers = [CoefficientStream_zero(f._is_sparse), ginv]
            for i in range(1, -self._fv):
                self._neg_powers.append(CoefficientStream_cauchy_product(self._neg_powers[-1], ginv))
        # Placeholder None to make this 1-based
        self._pos_powers = [None, g]
        val = self._fv * self._gv
        super().__init__(f, g, f._is_sparse, val)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_composition)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 1)
            sage: h = CoefficientStream_composition(f, g)
            sage: h.get_coefficient(5)
            527
            sage: [h.get_coefficient(i) for i in range(10)]
            [0, 1, 6, 28, 124, 527, 2172, 8755, 34704, 135772]
        """
        if n < 0:
            return sum(self._left[i] * self._neg_powers[-i][n] for i in range(self._fv, n // self._gv + 1))
        # n > 0
        while len(self._pos_powers) <= n // self._gv:
            self._pos_powers.append(CoefficientStream_cauchy_product(self._pos_powers[-1], self._right))
        ret = sum(self._left[i] * self._neg_powers[-i][n] for i in range(self._fv, 0))
        if n == 0:
            ret += self._left[0]
        return ret + sum(self._left[i] * self._pos_powers[i][n] for i in range(1, n // self._gv+1))

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_composition)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, False, 1)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^3, ZZ, False, 1)
            sage: h = CoefficientStream_composition(f, g)
            sage: n = h.iterate_coefficients()
            sage: [next(n) for i in range(10)]
            [1, 9, 44, 207, 991, 4752, 22769, 109089, 522676, 2504295]
        """
        n = self._approximate_valuation
        while True:
            yield self.get_coefficient(n)
            n += 1


#####################################################################
# Unary operations

class CoefficientStream_scalar(CoefficientStream_unary):
    """
    Operator for multiplying a coefficient stream with a scalar.

    INPUT:

    - ``series`` -- a :class:`CoefficientStream`
    - ``scalar`` -- a scalar

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_scalar, CoefficientStream_coefficient_function)
        sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
        sage: g = CoefficientStream_scalar(f, 2)
        sage: [g[i] for i in range(10)]
        [0, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    """
    def __init__(self, series, scalar):
        """
        Initialize.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_scalar, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: -1, ZZ, True, 0)
            sage: g = CoefficientStream_scalar(f, 3)
            sage: [g[i] for i in range(10)]
            [-3, -3, -3, -3, -3, -3, -3, -3, -3, -3]
        """
        self._scalar = scalar
        super().__init__(series, series._is_sparse, series._approximate_valuation)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_scalar, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
            sage: g = CoefficientStream_scalar(f, 3)
            sage: g.get_coefficient(5)
            15
            sage: [g.get_coefficient(i) for i in range(10)]
            [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
        """
        return self._series[n] * self._scalar

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_scalar, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, False, 1)
            sage: g = CoefficientStream_scalar(f, 4)
            sage: n = g.iterate_coefficients()
            sage: [next(n) for i in range(10)]
            [4, 16, 36, 64, 100, 144, 196, 256, 324, 400]
        """
        n = self._offset
        while True:
            yield self._series[n] * self._scalar
            n += 1


class CoefficientStream_neg(CoefficientStream_unary):
    """
    Operator for negative of the stream.

    INPUT:

    - ``series`` -- a :class:`CoefficientStream`

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_neg, CoefficientStream_coefficient_function)
        sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
        sage: g = CoefficientStream_neg(f)
        sage: [g[i] for i in range(10)]
        [0, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    """

    def __init__(self, series):
        """
        Initialize.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_neg, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: -1, ZZ, True, 0)
            sage: g = CoefficientStream_neg(f)
            sage: [g[i] for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        super().__init__(series, series._is_sparse, series._approximate_valuation)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_neg, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
            sage: g = CoefficientStream_neg(f)
            sage: g.get_coefficient(5)
            -5
            sage: [g.get_coefficient(i) for i in range(10)]
            [0, -1, -2, -3, -4, -5, -6, -7, -8, -9]
        """
        return -self._series[n]

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_neg, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, False, 1)
            sage: g = CoefficientStream_neg(f)
            sage: n = g.iterate_coefficients()
            sage: [next(n) for i in range(10)]
            [-1, -4, -9, -16, -25, -36, -49, -64, -81, -100]
        """
        n = self._offset
        while True:
            yield -self._series[n]
            n += 1


class CoefficientStream_cauchy_inverse(CoefficientStream_unary):
    """
    Operator for multiplicative inverse of the stream.

    INPUT:

    - ``series`` -- a :class:`CoefficientStream`

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_cauchy_inverse, CoefficientStream_coefficient_function)
        sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
        sage: g = CoefficientStream_cauchy_inverse(f)
        sage: [g[i] for i in range(10)]
        [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    """
    def __init__(self, series):
        """
        Initialize.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_cauchy_inverse, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: -1, ZZ, True, 0)
            sage: g = CoefficientStream_cauchy_inverse(f)
            sage: [g[i] for i in range(10)]
            [-1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        """
        v = series.valuation()
        super().__init__(series, series._is_sparse, -v)

        self._ainv = ~series[v]
        self._zero = ZZ.zero()

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_cauchy_inverse, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
            sage: g = CoefficientStream_cauchy_inverse(f)
            sage: g.get_coefficient(5)
            0
            sage: [g.get_coefficient(i) for i in range(10)]
            [-2, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        """
        v = self._approximate_valuation
        if n == v:
            return self._ainv
        c = self._zero
        for k in range(v, n):
            c += self[k] * self._series[n - v - k]
        return -c * self._ainv

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_cauchy_inverse, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, False, 1)
            sage: g = CoefficientStream_cauchy_inverse(f)
            sage: n = g.iterate_coefficients()
            sage: [next(n) for i in range(10)]
            [1, -4, 7, -8, 8, -8, 8, -8, 8, -8]
        """
        v = self._approximate_valuation  # shorthand name
        n = 0  # Counts the number of places from the valuation
        yield self._ainv
        # Note that first entry of the cache will correspond to z^v
        while True:
            n += 1
            c = self._zero
            m = min(len(self._cache), n)
            for k in range(m):
                c += self._cache[k] * self._series[n - v - k]
            for k in range(v+m, v+n):
                c += self[k] * self._series[n - k]
            yield -c * self._ainv


class CoefficientStream_apply_coeff(CoefficientStream_unary):
    """
    Return the stream with ``function`` applied to each coefficient of this series.

    INPUT:

    - ``series`` -- a :class:`CoefficientStream`
    - ``function`` -- a function that modifies the elements of the stream
    - ``ring`` -- the base ring of the stream

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_apply_coeff, CoefficientStream_coefficient_function)
        sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
        sage: g = CoefficientStream_apply_coeff(f, lambda n: -n, ZZ)
        sage: [g[i] for i in range(10)]
        [0, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    """
    def __init__(self, series, function, ring):
        """
        Initialize.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_apply_coeff, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: -1, ZZ, True, 0)
            sage: g = CoefficientStream_apply_coeff(f, lambda n: n + 1, ZZ)
            sage: [g[i] for i in range(10)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        """
        self._function = function
        self._ring = ring
        super().__init__(series, series._is_sparse, series._approximate_valuation)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_apply_coeff, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
            sage: g = CoefficientStream_apply_coeff(f, lambda n: n^2, ZZ)
            sage: g.get_coefficient(5)
            25
            sage: [g.get_coefficient(i) for i in range(10)]
            [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
        """
        c = self._ring(self._function(self._series[n]))
        return c

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_apply_coeff, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, False, 1)
            sage: g = CoefficientStream_apply_coeff(f, lambda n: 2*n, ZZ)
            sage: n = g.iterate_coefficients()
            sage: [next(n) for i in range(10)]
            [2, 8, 18, 32, 50, 72, 98, 128, 162, 200]
        """
        n = self._offset
        while True:
            c = self._ring(self._function(self._series[n]))
            yield c
            n += 1
