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

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_inverse
    sage: ginv = CoefficientStream_cauchy_inverse(g)
    sage: h = CoefficientStream_cauchy_product(f, ginv)
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

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_lmul
    sage: h = CoefficientStream_lmul(f, 2)
    sage: [h[i] for i in range(10)]
    [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]

The multiplicative inverse of a series can also be obtained::

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_inverse
    sage: h = CoefficientStream_cauchy_inverse(g)
    sage: [h[i] for i in range(10)]
    [-2, 1, 0, 0, 0, 0, 0, 0, 0, 0]

Functions can also be applied to a coefficient stream::

    sage: from sage.data_structures.coefficient_stream import CoefficientStream_map_coefficients
    sage: h = CoefficientStream_map_coefficients(f, lambda n: n^2, QQ)
    sage: [h[i] for i in range(10)]
    [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

AUTHORS:

- Kwankyu Lee (2019-02-24): initial version
- Tejasvi Chebrolu, Martin Rubey, Travis Scrimshaw (2021-08):
  refactored and expanded functionality
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


class CoefficientStream():
    """
    Abstract base class for all coefficient streams.

    INPUT:

    - ``sparse`` -- boolean; whether the implementation of the stream is sparse
    - ``approximate_order`` -- integer; a lower bound for the order
      of the stream
    """
    def __init__(self, sparse, approximate_order):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream
            sage: CS = CoefficientStream(True, 1)
        """
        self._is_sparse = sparse
        self._approximate_order = approximate_order


class CoefficientStream_inexact(CoefficientStream):
    """
    An abstract base class for the stream when we do not know it is
    eventually geometric.

    INPUT:

    - ``sparse`` -- boolean; whether the implementation of the stream is sparse
    - ``approximate_order`` -- integer; a lower bound for the order
      of the stream
    """
    def __init__(self, is_sparse, approximate_order):
        """
        Initialize the stream class for a CoefficientStream when it is not
        or it cannot be determined if it is eventually geometric.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_inexact
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
            sage: g = CoefficientStream_coefficient_function(lambda n: n, QQ, False, 0)
            sage: isinstance(g, CoefficientStream_inexact)
            True
        """
        super().__init__(is_sparse, approximate_order)

        if self._is_sparse:
            self._cache = dict()  # cache of known coefficients
        else:
            self._cache = list()
            self._offset = approximate_order
            self._iter = self.iterate_coefficients()

    def __getstate__(self):
        """
        Build the dictionary for pickling ``self``.

        We remove the cache from the pickle information when it is a dense
        implementation as iterators cannot be pickled.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_product
            sage: h = CoefficientStream_exact([1], True)
            sage: g = CoefficientStream_exact([1, -1, -1], True)
            sage: u = CoefficientStream_cauchy_product(h, g)
            sage: [u[i] for i in range(10)]
            [1, -1, -1, 0, 0, 0, 0, 0, 0, 0]
            sage: u._cache
            {0: 1, 1: -1, 2: -1, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0}
            sage: m = loads(dumps(u))
            sage: m._cache
            {0: 1, 1: -1, 2: -1, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0}
            sage: [m[i] for i in range(10)]
            [1, -1, -1, 0, 0, 0, 0, 0, 0, 0]

            sage: h = CoefficientStream_exact([1], False)
            sage: g = CoefficientStream_exact([1, -1, -1], False)
            sage: u = CoefficientStream_cauchy_product(h, g)
            sage: [u[i] for i in range(10)]
            [1, -1, -1, 0, 0, 0, 0, 0, 0, 0]
            sage: u._cache
            [1, -1, -1, 0, 0, 0, 0, 0, 0, 0]
            sage: m = loads(dumps(u))
            sage: m._cache
            []
            sage: [m[i] for i in range(10)]
            [1, -1, -1, 0, 0, 0, 0, 0, 0, 0]
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
        Build an object from ``d``.

        INPUT:

        - ``d`` -- a dictionary that needs to be unpickled

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: h = CoefficientStream_exact([-1], True)
            sage: g = CoefficientStream_exact([1, -1], True)
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_product
            sage: u = CoefficientStream_cauchy_product(h, g)
            sage: [u[i] for i in range(10)]
            [-1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
            sage: loads(dumps(u)) == u
            True
        """
        self.__dict__ = d
        if not self._is_sparse:
            self._iter = self.iterate_coefficients()
            self._cache = []

    def __getitem__(self, n):
        """
        Return the `n`-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the index

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
        if n < self._approximate_order:
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
        n = self._approximate_order
        while True:
            yield self.get_coefficient(n)
            n += 1

    def order(self):
        r"""
        Return the order of ``self``, which is the minimum index ``n`` such
        that ``self[n]`` is nonzero.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
            sage: f = CoefficientStream_coefficient_function(lambda n: n, QQ, True, 0)
            sage: f.order()
            1
        """
        if self._is_sparse:
            n = self._approximate_order
            cache = self._cache
            while True:
                if n in cache:
                    if cache[n]:
                        self._approximate_order = n
                        return n
                    n += 1
                else:
                    if self[n] != 0:
                        self._approximate_order = n
                        return n
                    n += 1
        else:
            n = self._approximate_order
            cache = self._cache
            while True:
                if n - self._offset < len(cache):
                    if cache[n - self._offset]:
                        self._approximate_order = n
                        return n
                    n += 1
                else:
                    if self[n] != 0:
                        self._approximate_order = n
                        return n
                    n += 1


class CoefficientStream_exact(CoefficientStream):
    r"""
    A stream of eventually constant coefficients.

    INPUT:

    - ``initial_values`` -- a list of initial values
    - ``is_sparse`` -- boolean; specifies whether the stream is sparse
    - ``order`` -- integer (default: 0); determining the degree
      of the first element of ``initial_values``
    - ``degree`` -- integer (optional); determining the degree
      of the first element which is known to be equal to ``constant``
    - ``constant`` -- integer (default: 0); the coefficient
      of every index larger than or equal to ``degree``
    """
    def __init__(self, initial_coefficients, is_sparse, constant=None, degree=None, order=None):
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
        if order is None:
            order = 0
        if degree is None:
            self._degree = order + len(initial_coefficients)
        else:
            self._degree = degree

        assert order + len(initial_coefficients) <= self._degree

        # We do not insist that the last entry of initial_coefficients
        #   is different from constant in case comparisons can be
        #   expensive such as in the symbolic ring
        for i, v in enumerate(initial_coefficients):
            if v:
                order += i
                initial_coefficients = initial_coefficients[i:]
                for j, w in enumerate(reversed(initial_coefficients)):
                    if w:
                        break
                    initial_coefficients.pop()
                self._initial_coefficients = tuple(initial_coefficients)
                break
        else:
            order = self._degree
            self._initial_coefficients = tuple()

        assert self._initial_coefficients or self._constant, "CoefficientStream_exact should only be used for non-zero streams"

        super().__init__(is_sparse, order)

    def __getitem__(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the index

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

            sage: s = CoefficientStream_exact([2], False, order=-1, constant=1)
            sage: [s[i] for i in range(-2, 5)]
            [0, 2, 1, 1, 1, 1, 1]

            sage: s = CoefficientStream_exact([2], False, order=-1, degree=2, constant=1)
            sage: [s[i] for i in range(-2, 5)]
            [0, 2, 0, 0, 1, 1, 1]

            sage: t = CoefficientStream_exact([0, 2, 0], False, order=-2, degree=2, constant=1)
            sage: t == s
            True
        """
        if n >= self._degree:
            return self._constant
        i = n - self._approximate_order
        if i < 0 or i >= len(self._initial_coefficients):
            return ZZ.zero()
        return self._initial_coefficients[i]

    def order(self):
        r"""
        Return the order of ``self``, which is the minimum index ``n`` such
        that ``self[n]`` is nonzero.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: s = CoefficientStream_exact([1], False)
            sage: s.order()
            0
        """
        return self._approximate_order

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: s = CoefficientStream_exact([1], False)
            sage: hash(s) == hash(s)
            True
        """
        return hash((self._initial_coefficients, self._degree, self._constant))

    def __eq__(self, other):
        """
        Test the equality between ``self`` and ``other``.

        INPUT:

        - ``other`` -- a coefficient stream

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: s = CoefficientStream_exact([2], False, order=-1, degree=2, constant=1)
            sage: t = CoefficientStream_exact([0, 2, 0], False, 1, 2, -2)
            sage: [s[i] for i in range(10)]
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: [t[i] for i in range(10)]
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: s == t
            True
            sage: s = CoefficientStream_exact([2], False, constant=1)
            sage: t = CoefficientStream_exact([2], False, order=-1, constant=1)
            sage: [s[i] for i in range(10)]
            [2, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: [t[i] for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: s == t
            False
            sage: t == t
            True
        """
        return (isinstance(other, type(self))
                and self._degree == other._degree
                and self._initial_coefficients == other._initial_coefficients
                and self._constant == other._constant)

    def polynomial_part(self, R):
        """
        Return the initial part of ``self`` as a Laurent polynomial in ``R``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: s = CoefficientStream_exact([2], False, order=-1, degree=2, constant=1)
            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: s.polynomial_part(L._laurent_poly_ring)
            2*z^-1
        """
        v = self._approximate_order
        return R(self._initial_coefficients).shift(v)


class CoefficientStream_coefficient_function(CoefficientStream_inexact):
    r"""
    Class that returns the elements in the coefficient stream.

    INPUT:

    - ``coefficient_function`` -- a function that generates the
      coefficients of the stream
    - ``ring`` -- the base ring
    - ``is_sparse`` -- boolean; specifies whether the stream is sparse
    - ``approximate_order`` -- integer; a lower bound for the order
      of the stream

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
        sage: f = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, False, 1)
        sage: f[3]
        9
        sage: [f[i] for i in range(10)]
        [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
    """

    def __init__(self, coefficient_function, ring, is_sparse, approximate_order):
        """
        Initialize.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, False, 1)
            sage: TestSuite(f).run(skip="_test_pickling")
        """
        self._coefficient_function = coefficient_function
        self._ring = ring
        super().__init__(is_sparse, approximate_order)

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

    - ``is_sparse`` -- boolean; which specifies whether the stream is sparse
    - ``approximate_order`` -- integer; a lower bound for the order
      of the stream

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import CoefficientStream_uninitialized
        sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
        sage: one = CoefficientStream_exact([1], True)
        sage: C = CoefficientStream_uninitialized(True, 0)
        sage: C._target
        sage: C._target = one
        sage: C.get_coefficient(4)
        0
    """
    def __init__(self, is_sparse, approximate_order):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_uninitialized
            sage: C = CoefficientStream_uninitialized(False, 0)
            sage: TestSuite(C).run(skip="_test_pickling")
        """
        self._target = None
        super().__init__(is_sparse, approximate_order)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_uninitialized
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: one = CoefficientStream_exact([1], True)
            sage: C = CoefficientStream_uninitialized(True, 0)
            sage: C._target
            sage: C._target = one
            sage: C.get_coefficient(0)
            1
        """
        return self._target[n]

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_uninitialized
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_exact
            sage: z = CoefficientStream_exact([1], True, order=1)
            sage: C = CoefficientStream_uninitialized(True, 0)
            sage: C._target
            sage: C._target = z
            sage: n = C.iterate_coefficients()
            sage: [next(n) for _ in range(10)]
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        """
        n = self._approximate_order
        while True:
            yield self._target[n]
            n += 1


class CoefficientStream_unary(CoefficientStream_inexact):
    r"""
    Base class for unary operators on coefficient streams.

    INPUT:

    - ``series`` -- :class:`CoefficientStream` the operator acts on

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_cauchy_inverse, CoefficientStream_lmul)
        sage: f = CoefficientStream_coefficient_function(lambda n: 2*n, ZZ, False, 1)
        sage: g = CoefficientStream_cauchy_inverse(f)
        sage: [g[i] for i in range(10)]
        [-1, 1/2, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: g = CoefficientStream_lmul(f, 2)
        sage: [g[i] for i in range(10)]
        [0, 4, 8, 12, 16, 20, 24, 28, 32, 36]
    """

    def __init__(self, series, *args, **kwargs):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_unary
            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_cauchy_inverse, CoefficientStream_exact)
            sage: f = CoefficientStream_exact([1, -1], False)
            sage: g = CoefficientStream_cauchy_inverse(f)
            sage: isinstance(g, CoefficientStream_unary)
            True
            sage: TestSuite(g).run()
        """
        self._series = series
        super().__init__(*args, **kwargs)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_unary
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
            sage: M = CoefficientStream_unary(CoefficientStream_coefficient_function(lambda n: 1, ZZ, False, 1), True, 0)
            sage: hash(M) == hash(M)
            True
        """
        return hash((type(self), self._series))

    def __eq__(self, other):
        """
        Test equality.

        INPUT:

        - ``other`` -- a stream of coefficients

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_rmul)
            sage: f = CoefficientStream_coefficient_function(lambda n: 2*n, ZZ, False, 1)
            sage: g = CoefficientStream_coefficient_function(lambda n: n, ZZ, False, 1)
            sage: h = CoefficientStream_rmul(f, 2)
            sage: n = CoefficientStream_rmul(g, 2)
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
    Base class for binary operators on coefficient streams.

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
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_binary
            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_add, CoefficientStream_cauchy_inverse, CoefficientStream_exact)
            sage: f1 = CoefficientStream_exact([1, -1], False)
            sage: g1 = CoefficientStream_cauchy_inverse(f1)
            sage: f2 = CoefficientStream_exact([1, 1], False)
            sage: g2 = CoefficientStream_cauchy_inverse(f2)
            sage: O = CoefficientStream_add(g1, g2)
            sage: isinstance(O, CoefficientStream_binary)
            True
            sage: TestSuite(O).run()
        """
        self._left = left
        self._right = right
        super().__init__(*args, **kwargs)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_binary
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function
            sage: M = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 0)
            sage: N = CoefficientStream_coefficient_function(lambda n: -2*n, ZZ, True, 0)
            sage: O = CoefficientStream_binary(M, N, True, 0)
            sage: hash(O) == hash(O)
            True
        """
        return hash((type(self), self._left, self._right))

    def __eq__(self, other):
        """
        Test equality.

        INPUT:

        - ``other`` -- a stream of coefficients

        EXAMPLES::

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
    Base class for commutative binary operators on coefficient streams.

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

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_add)
            sage: f = CoefficientStream_coefficient_function(lambda n: 2*n, ZZ, True, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
            sage: h = CoefficientStream_add(f, g)
            sage: u = CoefficientStream_add(g, f)
            sage: hash(h) == hash(u)
            True
        """
        return hash((type(self), frozenset([self._left, self._right])))

    def __eq__(self, other):
        """
        Test the equality between ``self`` and ``other``.

        INPUT:

        - ``other`` -- a stream of coefficients

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
        if not isinstance(other, type(self)):
            return False
        if self._left == other._left and self._right == other._right:
            return True
        if self._left == other._right and self._right == other._left:
            return True
        return False


class CoefficientStream_zero(CoefficientStream):
    """
    A coefficient stream that is exactly equal to zero.

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
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_zero
            sage: s = CoefficientStream_zero(False)
            sage: TestSuite(s).run()
        """
        return super().__init__(sparse, 0)

    def __getitem__(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the index

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_zero
            sage: s = CoefficientStream_zero(True)
            sage: s[1]
            0
            sage: sum([s[i] for i in range(10)])
            0
        """
        return ZZ.zero()

    def order(self):
        r"""
        Return the order of ``self``, which is ``infinity``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_zero
            sage: s = CoefficientStream_zero(True)
            sage: s.order()
            +Infinity
        """
        return infinity

    def __eq__(self, other):
        """
        Check equality of ``self`` and ``other`` ignoring sparsity.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_zero
            sage: CoefficientStream_zero(True) == CoefficientStream_zero(False)
            True
        """
        return self is other or isinstance(other, CoefficientStream_zero)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_zero
            sage: s = CoefficientStream_zero(False)
            sage: a = hash(s); a
            0
            sage: t = CoefficientStream_zero(False)
            sage: b = hash(t); b
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
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_add)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 0)
            sage: h = CoefficientStream_add(f, g)
        """
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = min(left._approximate_order, right._approximate_order)
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
        initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_sub)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 0)
            sage: h = CoefficientStream_sub(f, g)
        """
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = min(left._approximate_order, right._approximate_order)
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


class CoefficientStream_cauchy_product(CoefficientStream_binary):
    """
    Operator for multiplication of two coefficient streams using the
    Cauchy product.

    We are *not* assuming commutativity of the coefficient ring here,
    only that the coefficient ring commutes with the (implicit) variable.

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
        initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_coefficient_function, CoefficientStream_cauchy_product)
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 0)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 0)
            sage: h = CoefficientStream_cauchy_product(f, g)
        """
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = left._approximate_order + right._approximate_order
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
        for k in range(self._left._approximate_order,
                       n - self._right._approximate_order + 1):
            val = self._left[k]
            if val:
                c += val * self._right[n-k]
        return c


class CoefficientStream_composition(CoefficientStream_binary):
    r"""
    Return ``f`` composed by ``g``.

    This is the composition `(f \circ g)(z) = f(g(z))`.

    INPUT:

    - ``f`` -- a :class:`CoefficientStream`
    - ``g`` -- a :class:`CoefficientStream` with positive order

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import CoefficientStream_composition, CoefficientStream_coefficient_function
        sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_product as CS_prod
        sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_inverse as CS_inv
        sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
        sage: g = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
        sage: h = CoefficientStream_composition(f, g, CS_prod, CS_inv)
        sage: [h[i] for i in range(10)]
        [0, 1, 3, 8, 20, 48, 112, 256, 576, 1280]
        sage: u = CoefficientStream_composition(g, f, CS_prod, CS_inv)
        sage: [u[i] for i in range(10)]
        [0, 1, 3, 8, 21, 55, 144, 377, 987, 2584]
    """
    def __init__(self, f, g, prod_stream, inv_stream):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function, CoefficientStream_composition
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_product as CS_prod
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_inverse as CS_inv
            sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 1)
            sage: h = CoefficientStream_composition(f, g, CS_prod, CS_inv)
        """
        assert g._approximate_order > 0
        self._fv = f._approximate_order
        self._gv = g._approximate_order
        self._prod_stream = prod_stream
        if self._fv < 0:
            ginv = inv_stream(g)
            # the constant part makes no contribution to the negative
            # we need this for the case so self._neg_powers[0][n] => 0
            self._neg_powers = [CoefficientStream_zero(f._is_sparse), ginv]
            for i in range(1, -self._fv):
                self._neg_powers.append(self._prod_stream(self._neg_powers[-1], ginv))
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

            sage: from sage.data_structures.coefficient_stream import CoefficientStream_coefficient_function, CoefficientStream_composition
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_product as CS_prod
            sage: from sage.data_structures.coefficient_stream import CoefficientStream_cauchy_inverse as CS_inv
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
            sage: g = CoefficientStream_coefficient_function(lambda n: n^2, ZZ, True, 1)
            sage: h = CoefficientStream_composition(f, g, CS_prod, CS_inv)
            sage: h.get_coefficient(5)
            527
            sage: [h.get_coefficient(i) for i in range(10)]
            [0, 1, 6, 28, 124, 527, 2172, 8755, 34704, 135772]
        """
        if n < 0:
            return sum(self._left[i] * self._neg_powers[-i][n] for i in range(self._fv, n // self._gv + 1))
        # n > 0
        while len(self._pos_powers) <= n // self._gv:
            self._pos_powers.append(self._prod_stream(self._pos_powers[-1], self._right))
        ret = sum(self._left[i] * self._neg_powers[-i][n] for i in range(self._fv, 0))
        if n == 0:
            ret += self._left[0]
        return ret + sum(self._left[i] * self._pos_powers[i][n] for i in range(1, n // self._gv+1))


#####################################################################
# Unary operations

class CoefficientStream_rmul(CoefficientStream_unary):
    """
    Operator for multiplying a coefficient stream with a scalar
    as ``scalar * self``.

    INPUT:

    - ``series`` -- a :class:`CoefficientStream`
    - ``scalar`` -- a scalar

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_rmul, CoefficientStream_coefficient_function)
        sage: W = algebras.DifferentialWeyl(QQ, names=('x',))                                     
        sage: x, dx = W.gens()
        sage: f = CoefficientStream_coefficient_function(lambda n: x^n, W, True, 1)
        sage: g = CoefficientStream_rmul(f, dx)
        sage: [g[i] for i in range(5)]
        [0, x*dx + 1, x^2*dx + 2*x, x^3*dx + 3*x^2, x^4*dx + 4*x^3]
    """
    def __init__(self, series, scalar):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_rmul, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: -1, ZZ, True, 0)
            sage: g = CoefficientStream_rmul(f, 3)
        """
        self._scalar = scalar
        super().__init__(series, series._is_sparse, series._approximate_order)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_rmul, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
            sage: g = CoefficientStream_rmul(f, 3)
            sage: g.get_coefficient(5)
            15
            sage: [g.get_coefficient(i) for i in range(10)]
            [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
        """
        return self._scalar * self._series[n]


class CoefficientStream_lmul(CoefficientStream_unary):
    """
    Operator for multiplying a coefficient stream with a scalar
    as ``self * scalar``.

    INPUT:

    - ``series`` -- a :class:`CoefficientStream`
    - ``scalar`` -- a scalar

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_lmul, CoefficientStream_coefficient_function)
        sage: W = algebras.DifferentialWeyl(QQ, names=('x',))                                     
        sage: x, dx = W.gens()
        sage: f = CoefficientStream_coefficient_function(lambda n: x^n, W, True, 1)
        sage: g = CoefficientStream_lmul(f, dx)
        sage: [g[i] for i in range(5)]
        [0, x*dx, x^2*dx, x^3*dx, x^4*dx]
    """
    def __init__(self, series, scalar):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_lmul, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: -1, ZZ, True, 0)
            sage: g = CoefficientStream_lmul(f, 3)
        """
        self._scalar = scalar
        super().__init__(series, series._is_sparse, series._approximate_order)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_lmul, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 1)
            sage: g = CoefficientStream_lmul(f, 3)
            sage: g.get_coefficient(5)
            15
            sage: [g.get_coefficient(i) for i in range(10)]
            [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
        """
        return self._series[n] * self._scalar


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
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_neg, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: -1, ZZ, True, 0)
            sage: g = CoefficientStream_neg(f)
        """
        super().__init__(series, series._is_sparse, series._approximate_order)

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
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_cauchy_inverse, CoefficientStream_exact)
            sage: f = CoefficientStream_exact([1, -1], False)
            sage: g = CoefficientStream_cauchy_inverse(f)
        """
        v = series.order()
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
        v = self._approximate_order
        if n == v:
            return self._ainv
        c = self._zero
        for k in range(v, n):
            l = self[k]
            if l:
                c += l * self._series[n - v - k]
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
        v = self._approximate_order  # shorthand name
        n = 0  # Counts the number of places from the order
        yield self._ainv
        # Note that first entry of the cache will correspond to z^v
        while True:
            n += 1
            c = self._zero
            m = min(len(self._cache), n)
            for k in range(m):
                l = self._cache[k]
                if l:
                    c += l * self._series[n - v - k]
            for k in range(v+m, v+n):
                l = self[k]
                if l:
                    c += l * self._series[n - k]
            yield -c * self._ainv


class CoefficientStream_map_coefficients(CoefficientStream_unary):
    r"""
    Return the stream with ``function`` applied to each nonzero
    coefficient of ``series``.

    INPUT:

    - ``series`` -- a :class:`CoefficientStream`
    - ``function`` -- a function that modifies the elements of the stream
    - ``ring`` -- the base ring of the stream

    EXAMPLES::

        sage: from sage.data_structures.coefficient_stream import (CoefficientStream_map_coefficients, CoefficientStream_coefficient_function)
        sage: f = CoefficientStream_coefficient_function(lambda n: 1, ZZ, True, 1)
        sage: g = CoefficientStream_map_coefficients(f, lambda n: -n, ZZ)
        sage: [g[i] for i in range(10)]
        [0, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    """
    def __init__(self, series, function, ring):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_map_coefficients, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: -1, ZZ, True, 0)
            sage: g = CoefficientStream_map_coefficients(f, lambda n: n + 1, ZZ)
            sage: TestSuite(g).run(skip="_test_pickling")
        """
        self._function = function
        self._ring = ring
        super().__init__(series, series._is_sparse, series._approximate_order)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.coefficient_stream import (CoefficientStream_map_coefficients, CoefficientStream_coefficient_function)
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, -1)
            sage: g = CoefficientStream_map_coefficients(f, lambda n: n^2 + 1, ZZ)
            sage: g.get_coefficient(5)
            26
            sage: [g.get_coefficient(i) for i in range(-1, 10)]
            [2, 0, 2, 5, 10, 17, 26, 37, 50, 65, 82]

            sage: R.<x,y> = ZZ[]
            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, -1)
            sage: g = CoefficientStream_map_coefficients(f, lambda n: n.degree() + 1, R)
            sage: [g.get_coefficient(i) for i in range(-1, 3)]
            [1, 0, 1, 1]

            sage: f = CoefficientStream_coefficient_function(lambda n: n, ZZ, True, 0)
            sage: g = CoefficientStream_map_coefficients(f, lambda n: 5, GF(3))
            sage: [g.get_coefficient(i) for i in range(10)]
            [0, 5, 5, 0, 5, 5, 0, 5, 5, 0]
        """
        c = self._ring(self._series[n])
        if c:
            return self._function(c)
        return c

