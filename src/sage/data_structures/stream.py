r"""
Streams

This module provides lazy implementations of basic operators on
streams. The classes implemented in this module can be used to build
up more complex streams for different kinds of series (Laurent,
Dirichlet, etc).

EXAMPLES:

Streams can be used as data structure for lazy Laurent series::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ)
    sage: f = L(lambda n: n, True)
    sage: f
    z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + 7*z^7 + ...
    sage: type(f._coeff_stream)
    <class 'sage.data_structures.stream.StreamFunction'>

There are basic unary and binary operators available for streams. For
example, we can add two streams::

    sage: from sage.data_structures.stream import *
    sage: f = StreamFunction(lambda n: n, QQ, True, 0)
    sage: [f[i] for i in range(10)]
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    sage: g = StreamFunction(lambda n: 1, QQ, True, 0)
    sage: [g[i] for i in range(10)]
    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    sage: h = StreamAdd(f, g)
    sage: [h[i] for i in range(10)]
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

We can subtract one stream from another::

    sage: h = StreamSub(f, g)
    sage: [h[i] for i in range(10)]
    [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]

There is a Cauchy product on streams::

    sage: h = StreamCauchyProduct(f, g)
    sage: [h[i] for i in range(10)]
    [0, 1, 3, 6, 10, 15, 21, 28, 36, 45]

We can compute the inverse corresponding to the Cauchy product::

    sage: ginv = StreamCauchyInverse(g)
    sage: h = StreamCauchyProduct(f, ginv)
    sage: [h[i] for i in range(10)]
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1]

Two streams can be composed::

    sage: g = StreamFunction(lambda n: n, QQ, True, 1)
    sage: h = StreamCauchyComposition(f, g)
    sage: [h[i] for i in range(10)]
    [0, 1, 4, 14, 46, 145, 444, 1331, 3926, 11434]

There is a unary negation operator::

    sage: h = StreamNeg(f)
    sage: [h[i] for i in range(10)]
    [0, -1, -2, -3, -4, -5, -6, -7, -8, -9]

More generally, we can multiply by a scalar::

    sage: h = StreamLmul(f, 2)
    sage: [h[i] for i in range(10)]
    [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]

Finally, we can apply an arbitrary functions to the elements of a stream::

    sage: h = StreamMapCoefficients(f, lambda n: n^2, QQ)
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


class Stream():
    """
    Abstract base class for all streams.

    INPUT:

    - ``sparse`` -- boolean; whether the implementation of the stream is sparse
    - ``approximate_order`` -- integer; a lower bound for the order
      of the stream
    """
    def __init__(self, sparse, approximate_order):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import Stream
            sage: CS = Stream(True, 1)
        """
        self._is_sparse = sparse
        self._approximate_order = approximate_order

    def __ne__(self, other):
        """
        Check inequality of ``self`` and ``other``.

        The default is to always return ``False`` as it usually
        cannot be decided whether they are equal.

        EXAMPLES::

            sage: from sage.data_structures.stream import Stream
            sage: CS = Stream(True, 1)
            sage: CS != CS
            False
            sage: CS != Stream(False, -2)
            False

        """
        return False

    def is_nonzero(self):
        r"""
        Return ``True`` if and only if this stream is known
        to be nonzero.

        The default implementation is ``False``.

        EXAMPLES::

            sage: from sage.data_structures.stream import Stream
            sage: CS = Stream(True, 1)
            sage: CS.is_nonzero()
            False
        """
        return False


class StreamInexact(Stream):
    """
    An abstract base class for the stream when we do not know it is
    eventually constant.

    INPUT:

    - ``sparse`` -- boolean; whether the implementation of the stream is sparse
    - ``approximate_order`` -- integer; a lower bound for the order
      of the stream
    """
    def __init__(self, is_sparse, approximate_order):
        """
        Initialize the stream class for a Stream when it is not
        or it cannot be determined if it is eventually geometric.

        TESTS::

            sage: from sage.data_structures.stream import StreamInexact
            sage: from sage.data_structures.stream import StreamFunction
            sage: g = StreamFunction(lambda n: n, QQ, False, 0)
            sage: isinstance(g, StreamInexact)
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

            sage: from sage.data_structures.stream import StreamExact
            sage: from sage.data_structures.stream import StreamCauchyProduct
            sage: h = StreamExact([1], True)
            sage: g = StreamExact([1, -1, -1], True)
            sage: u = StreamCauchyProduct(h, g)
            sage: [u[i] for i in range(10)]
            [1, -1, -1, 0, 0, 0, 0, 0, 0, 0]
            sage: u._cache
            {0: 1, 1: -1, 2: -1, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0}
            sage: m = loads(dumps(u))
            sage: m._cache
            {0: 1, 1: -1, 2: -1, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0}
            sage: [m[i] for i in range(10)]
            [1, -1, -1, 0, 0, 0, 0, 0, 0, 0]

            sage: h = StreamExact([1], False)
            sage: g = StreamExact([1, -1, -1], False)
            sage: u = StreamCauchyProduct(h, g)
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
            # We cannot pickle a generator object, so we remove it
            # and the cache from the pickle information.
            del d["_iter"]
            del d["_cache"]
        return d

    def __setstate__(self, d):
        """
        Build an object from ``d``.

        INPUT:

        - ``d`` -- a dictionary that needs to be unpickled

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamExact
            sage: h = StreamExact([-1], True)
            sage: g = StreamExact([1, -1], True)
            sage: from sage.data_structures.stream import StreamCauchyProduct
            sage: u = StreamCauchyProduct(h, g)
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

            sage: from sage.data_structures.stream import StreamFunction
            sage: f = StreamFunction(lambda n: n^2, QQ, True, 0)
            sage: f[3]
            9
            sage: f._cache
            {3: 9}
            sage: [f[i] for i in range(10)]
            [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
            sage: f._cache
            {0: 0, 1: 1, 2: 4, 3: 9, 4: 16, 5: 25, 6: 36, 7: 49, 8: 64, 9: 81}

            sage: f = StreamFunction(lambda n: n^2, QQ, False, 0)
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
                # It is important to extend by generator:
                # self._iter might recurse, and thereby extend the
                # cache itself, too.
                self._cache.extend(next(self._iter) for _ in range(a, n+1))
            c = self._cache[i]

        return c

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamFunction, StreamCauchyComposition
            sage: f = StreamFunction(lambda n: 1, ZZ, False, 1)
            sage: g = StreamFunction(lambda n: n^3, ZZ, False, 1)
            sage: h = StreamCauchyComposition(f, g)
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

            sage: from sage.data_structures.stream import StreamFunction
            sage: f = StreamFunction(lambda n: n, QQ, True, 0)
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
                    if self[n]:
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
                    if self[n]:
                        self._approximate_order = n
                        return n
                    n += 1

    def __ne__(self, other):
        """
        Check inequality of ``self`` and ``other``.

        Check if there are any differences in the caches to see if they
        are known to be not equal.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamFunction
            sage: f = StreamFunction(lambda n: n, QQ, True, 0)
            sage: g = StreamFunction(lambda n: n^2, QQ, True, 0)
            sage: f != g
            False
            sage: f[1], g[1]
            (1, 1)
            sage: f != g
            False
            sage: f[3], g[4]
            (3, 16)
            sage: f != g
            False
            sage: f[2], g[2]
            (2, 4)
            sage: f != g
            True

        Checking the dense implementation::

            sage: f = StreamFunction(lambda n: n if n > 0 else 0, QQ, False, -3)
            sage: g = StreamFunction(lambda n: n^2, QQ, False, 0)
            sage: f != g
            False
            sage: g != f
            False
            sage: _ = f[1], g[1]
            sage: f != g
            False
            sage: g != f
            False
            sage: _ = f[2], g[2]
            sage: f != g
            True
            sage: g != f
            True

            sage: f = StreamFunction(lambda n: n if n > 0 else 0, QQ, False, -3)
            sage: g = StreamFunction(lambda n: n^2, QQ, False, 0)
            sage: _ = f[5], g[1]
            sage: f != g
            False
            sage: g != f
            False
            sage: _ = g[2]
            sage: f != g
            True
            sage: g != f
            True

            sage: f = StreamFunction(lambda n: n if n > 0 else 0, QQ, False, -3)
            sage: g = StreamFunction(lambda n: n^2, QQ, False, 0)
            sage: _ = g[5], f[1]
            sage: f != g
            False
            sage: g != f
            False
            sage: _ = f[2]
            sage: f != g
            True
            sage: g != f
            True
        """
        if not isinstance(other, StreamInexact):
            return False

        if self._is_sparse:
            for i in self._cache:
                if i in other._cache and other._cache[i] != self._cache[i]:
                    return True
        else: # they are dense
            # Make ``self`` have the smaller approximate order.
            if self._approximate_order > other._approximate_order:
                self, other = other, self
            saorder = self._approximate_order
            soffset = self._offset
            oaorder = other._approximate_order
            ooffset = other._offset
            end = min(oaorder, soffset + len(self._cache))
            for i in range(saorder, end):
                if self._cache[i-soffset]:
                    return True
            # now check all known values
            end = min(soffset + len(self._cache), ooffset + len(other._cache))
            for i in range(oaorder, end):
                if self._cache[i-soffset] != other._cache[i-ooffset]:
                    return True

        return False

class StreamExact(Stream):
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

            sage: from sage.data_structures.stream import StreamExact
            sage: StreamExact([], False)
            Traceback (most recent call last):
            ...
            AssertionError: StreamExact should only be used for non-zero streams
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

        # We do not insist that the last entry of
        # initial_coefficients is different from constant in case
        # comparisons can be expensive such as in the symbolic ring
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

        assert self._initial_coefficients or self._constant, "StreamExact should only be used for non-zero streams"

        super().__init__(is_sparse, order)

    def __getitem__(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the index

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamExact
            sage: s = StreamExact([1], False)
            sage: [s[i] for i in range(-2, 5)]
            [0, 0, 1, 0, 0, 0, 0]

            sage: s = StreamExact([], False, constant=1)
            sage: [s[i] for i in range(-2, 5)]
            [0, 0, 1, 1, 1, 1, 1]

            sage: s = StreamExact([2], False, constant=1)
            sage: [s[i] for i in range(-2, 5)]
            [0, 0, 2, 1, 1, 1, 1]

            sage: s = StreamExact([2], False, order=-1, constant=1)
            sage: [s[i] for i in range(-2, 5)]
            [0, 2, 1, 1, 1, 1, 1]

            sage: s = StreamExact([2], False, order=-1, degree=2, constant=1)
            sage: [s[i] for i in range(-2, 5)]
            [0, 2, 0, 0, 1, 1, 1]

            sage: t = StreamExact([0, 2, 0], False, order=-2, degree=2, constant=1)
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

            sage: from sage.data_structures.stream import StreamExact
            sage: s = StreamExact([1], False)
            sage: s.order()
            0
        """
        return self._approximate_order

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamExact
            sage: s = StreamExact([1], False)
            sage: hash(s) == hash(s)
            True
        """
        return hash((self._initial_coefficients, self._degree, self._constant))

    def __eq__(self, other):
        """
        Test the equality between ``self`` and ``other``.

        INPUT:

        - ``other`` -- a stream

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamExact
            sage: s = StreamExact([2], False, order=-1, degree=2, constant=1)
            sage: t = StreamExact([0, 2, 0], False, 1, 2, -2)
            sage: [s[i] for i in range(10)]
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: [t[i] for i in range(10)]
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: s == t
            True
            sage: s = StreamExact([2], False, constant=1)
            sage: t = StreamExact([2], False, order=-1, constant=1)
            sage: [s[i] for i in range(10)]
            [2, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: [t[i] for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: s == t
            False
            sage: t == t
            True

            sage: s = StreamExact([2], False, order=0, degree=5, constant=1)
            sage: t = StreamExact([2], False, order=-1, degree=5, constant=1)
            sage: s == t
            False
        """
        return (isinstance(other, type(self))
                and self._degree == other._degree
                and self._approximate_order == other._approximate_order
                and self._initial_coefficients == other._initial_coefficients
                and self._constant == other._constant)

    def __ne__(self, other):
        """
        Test inequality between ``self`` and ``other``.

        INPUT:

        - ``other`` -- a stream

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamExact
            sage: s = StreamExact([2], False, order=-1, degree=2, constant=1)
            sage: t = StreamExact([0, 2, 0], False, 1, 2, -2)
            sage: s != t
            False
            sage: s = StreamExact([2], False, constant=1)
            sage: t = StreamExact([2], False, order=-1, constant=1)
            sage: s != t
            True

        When it is not known, then both equality and inequality
        return ``False``::

            sage: from sage.data_structures.stream import StreamFunction
            sage: f = StreamFunction(lambda n: 2 if n == 0 else 1, ZZ, False, 0)
            sage: s == f
            False
            sage: s != f
            False
            sage: [s[i] for i in range(-3, 5)]
            [0, 0, 0, 2, 1, 1, 1, 1]
            sage: [f[i] for i in range(-3, 5)]
            [0, 0, 0, 2, 1, 1, 1, 1]
        """
        if isinstance(other, type(self)):
            return (self._degree != other._degree
                    or self._approximate_order != other._approximate_order
                    or self._initial_coefficients != other._initial_coefficients
                    or self._constant != other._constant)
        return False

    def is_nonzero(self):
        r"""
        Return ``True`` if and only if this stream is known
        to be nonzero.

        An assumption of this class is that it is nonzero.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamExact
            sage: s = StreamExact([2], False, order=-1, degree=2, constant=1)
            sage: s.is_nonzero()
            True
        """
        return True

    def _polynomial_part(self, R):
        """
        Return the initial part of ``self`` as a Laurent polynomial in ``R``.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamExact
            sage: s = StreamExact([2], False, order=-1, degree=2, constant=1)
            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: s._polynomial_part(L._laurent_poly_ring)
            2*z^-1
        """
        v = self._approximate_order
        return R(self._initial_coefficients).shift(v)


class StreamFunction(StreamInexact):
    r"""
    Class that creates a stream from a function on the integers.

    INPUT:

    - ``function`` -- a function that generates the
      coefficients of the stream
    - ``ring`` -- the base ring
    - ``is_sparse`` -- boolean; specifies whether the stream is sparse
    - ``approximate_order`` -- integer; a lower bound for the order
      of the stream

    EXAMPLES::

        sage: from sage.data_structures.stream import StreamFunction
        sage: f = StreamFunction(lambda n: n^2, ZZ, False, 1)
        sage: f[3]
        9
        sage: [f[i] for i in range(10)]
        [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
    """

    def __init__(self, function, ring, is_sparse, approximate_order):
        """
        Initialize.

        TESTS::

            sage: from sage.data_structures.stream import StreamFunction
            sage: f = StreamFunction(lambda n: 1, ZZ, False, 1)
            sage: TestSuite(f).run(skip="_test_pickling")
        """
        self._function = function
        self._ring = ring
        super().__init__(is_sparse, approximate_order)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamFunction
            sage: f = StreamFunction(lambda n: n, QQ, True, 0)
            sage: f.get_coefficient(4)
            4
        """
        return self._ring(self._function(n))

    def iterate_coefficients(self):
        """
        A generator for the coefficients of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamFunction
            sage: f = StreamFunction(lambda n: 1, QQ, False, 0)
            sage: n = f.iterate_coefficients()
            sage: [next(n) for _ in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        n = self._offset
        ring = self._ring
        while True:
            yield ring(self._function(n))
            n += 1


class StreamUninitialized(StreamInexact):
    r"""
    Coefficient stream for an uninitialized series.

    INPUT:

    - ``is_sparse`` -- boolean; which specifies whether the stream is sparse
    - ``approximate_order`` -- integer; a lower bound for the order
      of the stream

    EXAMPLES::

        sage: from sage.data_structures.stream import StreamUninitialized
        sage: from sage.data_structures.stream import StreamExact
        sage: one = StreamExact([1], True)
        sage: C = StreamUninitialized(True, 0)
        sage: C._target
        sage: C._target = one
        sage: C.get_coefficient(4)
        0
    """
    def __init__(self, is_sparse, approximate_order):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import StreamUninitialized
            sage: C = StreamUninitialized(False, 0)
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

            sage: from sage.data_structures.stream import StreamUninitialized
            sage: from sage.data_structures.stream import StreamExact
            sage: one = StreamExact([1], True)
            sage: C = StreamUninitialized(True, 0)
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

            sage: from sage.data_structures.stream import StreamUninitialized
            sage: from sage.data_structures.stream import StreamExact
            sage: z = StreamExact([1], True, order=1)
            sage: C = StreamUninitialized(True, 0)
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


class StreamUnary(StreamInexact):
    r"""
    Base class for unary operators on coefficient streams.

    INPUT:

    - ``series`` -- :class:`Stream` the operator acts on

    EXAMPLES::

        sage: from sage.data_structures.stream import (StreamFunction, StreamCauchyInverse, StreamLmul)
        sage: f = StreamFunction(lambda n: 2*n, ZZ, False, 1)
        sage: g = StreamCauchyInverse(f)
        sage: [g[i] for i in range(10)]
        [-1, 1/2, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: g = StreamLmul(f, 2)
        sage: [g[i] for i in range(10)]
        [0, 4, 8, 12, 16, 20, 24, 28, 32, 36]
    """

    def __init__(self, series, *args, **kwargs):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import StreamUnary
            sage: from sage.data_structures.stream import (StreamCauchyInverse, StreamExact)
            sage: f = StreamExact([1, -1], False)
            sage: g = StreamCauchyInverse(f)
            sage: isinstance(g, StreamUnary)
            True
            sage: TestSuite(g).run()
        """
        self._series = series
        super().__init__(*args, **kwargs)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamUnary
            sage: from sage.data_structures.stream import StreamFunction
            sage: M = StreamUnary(StreamFunction(lambda n: 1, ZZ, False, 1), True, 0)
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

            sage: from sage.data_structures.stream import (StreamFunction, StreamRmul)
            sage: f = StreamFunction(lambda n: 2*n, ZZ, False, 1)
            sage: g = StreamFunction(lambda n: n, ZZ, False, 1)
            sage: h = StreamRmul(f, 2)
            sage: n = StreamRmul(g, 2)
            sage: h == n
            False
            sage: n == n
            True
            sage: h == h
            True
        """
        return isinstance(other, type(self)) and self._series == other._series


class StreamBinary(StreamInexact):
    """
    Base class for binary operators on coefficient streams.

    INPUT:

    - ``left`` -- :class:`Stream` for the left side of the operator
    - ``right`` -- :class:`Stream` for the right side of the operator

    EXAMPLES::

        sage: from sage.data_structures.stream import (StreamFunction, StreamAdd, StreamSub)
        sage: f = StreamFunction(lambda n: 2*n, ZZ, True, 0)
        sage: g = StreamFunction(lambda n: n, ZZ, True, 1)
        sage: h = StreamAdd(f, g)
        sage: [h[i] for i in range(10)]
        [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
        sage: h = StreamSub(f, g)
        sage: [h[i] for i in range(10)]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    """

    def __init__(self, left, right, *args, **kwargs):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import StreamBinary
            sage: from sage.data_structures.stream import (StreamAdd, StreamCauchyInverse, StreamExact)
            sage: f1 = StreamExact([1, -1], False)
            sage: g1 = StreamCauchyInverse(f1)
            sage: f2 = StreamExact([1, 1], False)
            sage: g2 = StreamCauchyInverse(f2)
            sage: O = StreamAdd(g1, g2)
            sage: isinstance(O, StreamBinary)
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

            sage: from sage.data_structures.stream import StreamBinary
            sage: from sage.data_structures.stream import StreamFunction
            sage: M = StreamFunction(lambda n: n, ZZ, True, 0)
            sage: N = StreamFunction(lambda n: -2*n, ZZ, True, 0)
            sage: O = StreamBinary(M, N, True, 0)
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

            sage: from sage.data_structures.stream import (StreamFunction, StreamCauchyProduct)
            sage: f = StreamFunction(lambda n: 2*n, ZZ, False, 1)
            sage: g = StreamFunction(lambda n: n, ZZ, False, 1)
            sage: h = StreamFunction(lambda n: 1, ZZ, False, 1)
            sage: t = StreamCauchyProduct(f, g)
            sage: u = StreamCauchyProduct(g, h)
            sage: v = StreamCauchyProduct(h, f)
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


class StreamBinaryCommutative(StreamBinary):
    r"""
    Base class for commutative binary operators on coefficient streams.

    EXAMPLES::

        sage: from sage.data_structures.stream import (StreamFunction, StreamAdd)
        sage: f = StreamFunction(lambda n: 2*n, ZZ, True, 0)
        sage: g = StreamFunction(lambda n: n, ZZ, True, 1)
        sage: h = StreamAdd(f, g)
        sage: [h[i] for i in range(10)]
        [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
        sage: u = StreamAdd(g, f)
        sage: [u[i] for i in range(10)]
        [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
        sage: h == u
        True
    """
    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamFunction, StreamAdd)
            sage: f = StreamFunction(lambda n: 2*n, ZZ, True, 0)
            sage: g = StreamFunction(lambda n: n, ZZ, True, 1)
            sage: h = StreamAdd(f, g)
            sage: u = StreamAdd(g, f)
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

            sage: from sage.data_structures.stream import (StreamFunction, StreamAdd)
            sage: f = StreamFunction(lambda n: 2*n, ZZ, True, 0)
            sage: g = StreamFunction(lambda n: n, ZZ, True, 1)
            sage: h = StreamAdd(f, g)
            sage: [h[i] for i in range(10)]
            [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
            sage: u = StreamAdd(g, f)
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


class StreamZero(Stream):
    """
    A coefficient stream that is exactly equal to zero.

    INPUT:

    - ``sparse`` -- boolean; whether the coefficient stream is sparse or not

    EXAMPLES::

        sage: from sage.data_structures.stream import StreamZero
        sage: s = StreamZero(True)
        sage: s[5]
        0
    """

    def __init__(self, sparse):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import StreamZero
            sage: s = StreamZero(False)
            sage: TestSuite(s).run()
        """
        return super().__init__(sparse, 0)

    def __getitem__(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the index

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamZero
            sage: s = StreamZero(True)
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

            sage: from sage.data_structures.stream import StreamZero
            sage: s = StreamZero(True)
            sage: s.order()
            +Infinity
        """
        return infinity

    def __eq__(self, other):
        """
        Check equality of ``self`` and ``other`` ignoring sparsity.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamZero
            sage: StreamZero(True) == StreamZero(False)
            True
        """
        return self is other or isinstance(other, StreamZero)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamZero
            sage: s = StreamZero(False)
            sage: a = hash(s); a
            0
            sage: t = StreamZero(False)
            sage: b = hash(t); b
            0
            sage: b == a
            True
        """
        return 0


#####################################################################
# Binary operations

class StreamAdd(StreamBinaryCommutative):
    """
    Operator for addition of two coefficient streams.

    INPUT:

    - ``left`` -- stream of coefficients on the left side of the operator
    - ``right`` -- stream of coefficients on the right side of the operator

    EXAMPLES::

        sage: from sage.data_structures.stream import (StreamAdd, StreamFunction)
        sage: f = StreamFunction(lambda n: n, ZZ, True, 0)
        sage: g = StreamFunction(lambda n: 1, ZZ, True, 0)
        sage: h = StreamAdd(f, g)
        sage: [h[i] for i in range(10)]
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        sage: u = StreamAdd(g, f)
        sage: [u[i] for i in range(10)]
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    """
    def __init__(self, left, right):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import (StreamFunction, StreamAdd)
            sage: f = StreamFunction(lambda n: 1, ZZ, True, 0)
            sage: g = StreamFunction(lambda n: n^2, ZZ, True, 0)
            sage: h = StreamAdd(f, g)
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

            sage: from sage.data_structures.stream import (StreamFunction, StreamAdd)
            sage: f = StreamFunction(lambda n: n, ZZ, True, 0)
            sage: g = StreamFunction(lambda n: n^2, ZZ, True, 0)
            sage: h = StreamAdd(f, g)
            sage: h.get_coefficient(5)
            30
            sage: [h.get_coefficient(i) for i in range(10)]
            [0, 2, 6, 12, 20, 30, 42, 56, 72, 90]
        """
        return self._left[n] + self._right[n]


class StreamSub(StreamBinary):
    """
    Operator for subtraction of two coefficient streams.

    INPUT:

    - ``left`` -- stream of coefficients on the left side of the operator
    - ``right`` -- stream of coefficients on the right side of the operator

    EXAMPLES::

        sage: from sage.data_structures.stream import (StreamSub, StreamFunction)
        sage: f = StreamFunction(lambda n: n, ZZ, True, 0)
        sage: g = StreamFunction(lambda n: 1, ZZ, True, 0)
        sage: h = StreamSub(f, g)
        sage: [h[i] for i in range(10)]
        [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8]
        sage: u = StreamSub(g, f)
        sage: [u[i] for i in range(10)]
        [1, 0, -1, -2, -3, -4, -5, -6, -7, -8]
    """

    def __init__(self, left, right):
        """
        initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import (StreamFunction, StreamSub)
            sage: f = StreamFunction(lambda n: 1, ZZ, True, 0)
            sage: g = StreamFunction(lambda n: n^2, ZZ, True, 0)
            sage: h = StreamSub(f, g)
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

            sage: from sage.data_structures.stream import (StreamFunction, StreamSub)
            sage: f = StreamFunction(lambda n: n, ZZ, True, 0)
            sage: g = StreamFunction(lambda n: n^2, ZZ, True, 0)
            sage: h = StreamSub(f, g)
            sage: h.get_coefficient(5)
            -20
            sage: [h.get_coefficient(i) for i in range(10)]
            [0, 0, -2, -6, -12, -20, -30, -42, -56, -72]
        """
        return self._left[n] - self._right[n]


class StreamCauchyProduct(StreamBinary):
    """
    Operator for multiplication of two coefficient streams using the
    Cauchy product.

    We are *not* assuming commutativity of the coefficient ring here,
    only that the coefficient ring commutes with the (implicit) variable.

    INPUT:

    - ``left`` -- stream of coefficients on the left side of the operator
    - ``right`` -- stream of coefficients on the right side of the operator

    EXAMPLES::

        sage: from sage.data_structures.stream import (StreamCauchyProduct, StreamFunction)
        sage: f = StreamFunction(lambda n: n, ZZ, True, 0)
        sage: g = StreamFunction(lambda n: 1, ZZ, True, 0)
        sage: h = StreamCauchyProduct(f, g)
        sage: [h[i] for i in range(10)]
        [0, 1, 3, 6, 10, 15, 21, 28, 36, 45]
        sage: u = StreamCauchyProduct(g, f)
        sage: [u[i] for i in range(10)]
        [0, 1, 3, 6, 10, 15, 21, 28, 36, 45]
    """
    def __init__(self, left, right):
        """
        initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import (StreamFunction, StreamCauchyProduct)
            sage: f = StreamFunction(lambda n: 1, ZZ, True, 0)
            sage: g = StreamFunction(lambda n: n^2, ZZ, True, 0)
            sage: h = StreamCauchyProduct(f, g)
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

            sage: from sage.data_structures.stream import (StreamFunction, StreamCauchyProduct)
            sage: f = StreamFunction(lambda n: n, ZZ, True, 0)
            sage: g = StreamFunction(lambda n: n^2, ZZ, True, 0)
            sage: h = StreamCauchyProduct(f, g)
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

    def is_nonzero(self):
        r"""
        Return ``True`` if and only if this stream is known
        to be nonzero.

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamFunction,
            ....:     StreamCauchyProduct, StreamCauchyInverse)
            sage: f = StreamFunction(lambda n: n, ZZ, True, 1)
            sage: g = StreamCauchyProduct(f, f)
            sage: g.is_nonzero()
            False
            sage: fi = StreamCauchyInverse(f)
            sage: h = StreamCauchyProduct(fi, fi)
            sage: h.is_nonzero()
            True
        """
        return self._left.is_nonzero() and self._right.is_nonzero()


class StreamCauchyComposition(StreamBinary):
    r"""
    Return ``f`` composed by ``g``.

    This is the composition `(f \circ g)(z) = f(g(z))`.

    INPUT:

    - ``f`` -- a :class:`Stream`
    - ``g`` -- a :class:`Stream` with positive order

    EXAMPLES::

        sage: from sage.data_structures.stream import StreamCauchyComposition, StreamFunction
        sage: f = StreamFunction(lambda n: n, ZZ, True, 1)
        sage: g = StreamFunction(lambda n: 1, ZZ, True, 1)
        sage: h = StreamCauchyComposition(f, g)
        sage: [h[i] for i in range(10)]
        [0, 1, 3, 8, 20, 48, 112, 256, 576, 1280]
        sage: u = StreamCauchyComposition(g, f)
        sage: [u[i] for i in range(10)]
        [0, 1, 3, 8, 21, 55, 144, 377, 987, 2584]
    """
    def __init__(self, f, g):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import StreamFunction, StreamCauchyComposition
            sage: f = StreamFunction(lambda n: 1, ZZ, True, 1)
            sage: g = StreamFunction(lambda n: n^2, ZZ, True, 1)
            sage: h = StreamCauchyComposition(f, g)
        """
        assert g._approximate_order > 0
        self._fv = f._approximate_order
        self._gv = g._approximate_order
        if self._fv < 0:
            ginv = StreamCauchyInverse(g)
            # The constant part makes no contribution to the negative.
            # We need this for the case so self._neg_powers[0][n] => 0.
            self._neg_powers = [StreamZero(f._is_sparse), ginv]
            for i in range(1, -self._fv):
                self._neg_powers.append(StreamCauchyProduct(self._neg_powers[-1], ginv))
        # Placeholder None to make this 1-based.
        self._pos_powers = [None, g]
        val = self._fv * self._gv
        super().__init__(f, g, f._is_sparse, val)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamFunction, StreamCauchyComposition
            sage: f = StreamFunction(lambda n: n, ZZ, True, 1)
            sage: g = StreamFunction(lambda n: n^2, ZZ, True, 1)
            sage: h = StreamCauchyComposition(f, g)
            sage: h.get_coefficient(5)
            527
            sage: [h.get_coefficient(i) for i in range(10)]
            [0, 1, 6, 28, 124, 527, 2172, 8755, 34704, 135772]
        """
        if n < 0:
            return sum(self._left[i] * self._neg_powers[-i][n] for i in range(self._fv, n // self._gv + 1))
        # n > 0
        while len(self._pos_powers) <= n // self._gv:
            self._pos_powers.append(StreamCauchyProduct(self._pos_powers[-1], self._right))
        ret = sum(self._left[i] * self._neg_powers[-i][n] for i in range(self._fv, 0))
        if n == 0:
            ret += self._left[0]
        return ret + sum(self._left[i] * self._pos_powers[i][n] for i in range(1, n // self._gv+1))


#####################################################################
# Unary operations

class StreamScalar(StreamInexact):
    """
    Base class for operators multiplying a coeffeicient stream
    by a scalar.
    """
    def __init__(self, series, scalar):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import (StreamRmul, StreamFunction)
            sage: f = StreamFunction(lambda n: -1, ZZ, True, 0)
            sage: g = StreamRmul(f, 3)
        """
        self._series = series
        self._scalar = scalar
        assert scalar, "the scalar must not be equal to 0"
        super().__init__(series._is_sparse, series._approximate_order)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamFunction
            sage: from sage.data_structures.stream import StreamRmul
            sage: a = StreamFunction(lambda n: 2*n, ZZ, False, 1)
            sage: f = StreamRmul(a, 2)
            sage: hash(f) == hash(f)
            True
        """
        return hash((type(self), self._series, self._scalar))

    def __eq__(self, other):
        """
        Test equality.

        INPUT:

        - ``other`` -- a stream of coefficients

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamFunction
            sage: from sage.data_structures.stream import StreamRmul, StreamLmul
            sage: a = StreamFunction(lambda n: 2*n, ZZ, False, 1)
            sage: b = StreamFunction(lambda n: n, ZZ, False, 1)
            sage: f = StreamRmul(a, 2)
            sage: f == StreamRmul(b, 2)
            False
            sage: f == StreamRmul(a, 2)
            True
            sage: f == StreamRmul(a, 3)
            False
            sage: f == StreamLmul(a, 3)
            False
        """
        return (isinstance(other, type(self)) and self._series == other._series
                and self._scalar == other._scalar)

    def is_nonzero(self):
        r"""
        Return ``True`` if and only if this stream is known
        to be nonzero.

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamRmul, StreamFunction)
            sage: f = StreamFunction(lambda n: n, ZZ, True, 1)
            sage: g = StreamRmul(f, 2)
            sage: g.is_nonzero()
            False

            sage: from sage.data_structures.stream import StreamCauchyInverse
            sage: fi = StreamCauchyInverse(f)
            sage: g = StreamRmul(fi, 2)
            sage: g.is_nonzero()
            True
        """
        return self._series.is_nonzero()


class StreamRmul(StreamScalar):
    """
    Operator for multiplying a coefficient stream with a scalar
    as ``scalar * self``.

    INPUT:

    - ``series`` -- a :class:`Stream`
    - ``scalar`` -- a scalar

    EXAMPLES::

        sage: from sage.data_structures.stream import (StreamRmul, StreamFunction)
        sage: W = algebras.DifferentialWeyl(QQ, names=('x',))                                     
        sage: x, dx = W.gens()
        sage: f = StreamFunction(lambda n: x^n, W, True, 1)
        sage: g = StreamRmul(f, dx)
        sage: [g[i] for i in range(5)]
        [0, x*dx + 1, x^2*dx + 2*x, x^3*dx + 3*x^2, x^4*dx + 4*x^3]
    """
    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamRmul, StreamFunction)
            sage: f = StreamFunction(lambda n: n, ZZ, True, 1)
            sage: g = StreamRmul(f, 3)
            sage: g.get_coefficient(5)
            15
            sage: [g.get_coefficient(i) for i in range(10)]
            [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
        """
        return self._scalar * self._series[n]


class StreamLmul(StreamScalar):
    """
    Operator for multiplying a coefficient stream with a scalar
    as ``self * scalar``.

    INPUT:

    - ``series`` -- a :class:`Stream`
    - ``scalar`` -- a scalar

    EXAMPLES::

        sage: from sage.data_structures.stream import (StreamLmul, StreamFunction)
        sage: W = algebras.DifferentialWeyl(QQ, names=('x',))                                     
        sage: x, dx = W.gens()
        sage: f = StreamFunction(lambda n: x^n, W, True, 1)
        sage: g = StreamLmul(f, dx)
        sage: [g[i] for i in range(5)]
        [0, x*dx, x^2*dx, x^3*dx, x^4*dx]
    """
    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamLmul, StreamFunction)
            sage: f = StreamFunction(lambda n: n, ZZ, True, 1)
            sage: g = StreamLmul(f, 3)
            sage: g.get_coefficient(5)
            15
            sage: [g.get_coefficient(i) for i in range(10)]
            [0, 3, 6, 9, 12, 15, 18, 21, 24, 27]
        """
        return self._series[n] * self._scalar


class StreamNeg(StreamUnary):
    """
    Operator for negative of the stream.

    INPUT:

    - ``series`` -- a :class:`Stream`

    EXAMPLES::

        sage: from sage.data_structures.stream import (StreamNeg, StreamFunction)
        sage: f = StreamFunction(lambda n: 1, ZZ, True, 1)
        sage: g = StreamNeg(f)
        sage: [g[i] for i in range(10)]
        [0, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    """
    def __init__(self, series):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import (StreamNeg, StreamFunction)
            sage: f = StreamFunction(lambda n: -1, ZZ, True, 0)
            sage: g = StreamNeg(f)
        """
        super().__init__(series, series._is_sparse, series._approximate_order)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamNeg, StreamFunction)
            sage: f = StreamFunction(lambda n: n, ZZ, True, 1)
            sage: g = StreamNeg(f)
            sage: g.get_coefficient(5)
            -5
            sage: [g.get_coefficient(i) for i in range(10)]
            [0, -1, -2, -3, -4, -5, -6, -7, -8, -9]
        """
        return -self._series[n]

    def is_nonzero(self):
        r"""
        Return ``True`` if and only if this stream is known
        to be nonzero.

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamNeg, StreamFunction)
            sage: f = StreamFunction(lambda n: n, ZZ, True, 1)
            sage: g = StreamNeg(f)
            sage: g.is_nonzero()
            False

            sage: from sage.data_structures.stream import StreamCauchyInverse
            sage: fi = StreamCauchyInverse(f)
            sage: g = StreamNeg(fi)
            sage: g.is_nonzero()
            True
        """
        return self._series.is_nonzero()

class StreamCauchyInverse(StreamUnary):
    """
    Operator for multiplicative inverse of the stream.

    INPUT:

    - ``series`` -- a :class:`Stream`

    EXAMPLES::

        sage: from sage.data_structures.stream import (StreamCauchyInverse, StreamFunction)
        sage: f = StreamFunction(lambda n: 1, ZZ, True, 1)
        sage: g = StreamCauchyInverse(f)
        sage: [g[i] for i in range(10)]
        [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    """
    def __init__(self, series):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import (StreamCauchyInverse, StreamExact)
            sage: f = StreamExact([1, -1], False)
            sage: g = StreamCauchyInverse(f)
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

            sage: from sage.data_structures.stream import (StreamCauchyInverse, StreamFunction)
            sage: f = StreamFunction(lambda n: n, ZZ, True, 1)
            sage: g = StreamCauchyInverse(f)
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

            sage: from sage.data_structures.stream import (StreamCauchyInverse, StreamFunction)
            sage: f = StreamFunction(lambda n: n^2, ZZ, False, 1)
            sage: g = StreamCauchyInverse(f)
            sage: n = g.iterate_coefficients()
            sage: [next(n) for i in range(10)]
            [1, -4, 7, -8, 8, -8, 8, -8, 8, -8]
        """
        v = self._approximate_order
        n = 0  # Counts the number of places from v.
        yield self._ainv
        # Note that the first entry of the cache will correspond to
        # z^v, when the stream corresponds to a Laurent series.
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

    def is_nonzero(self):
        r"""
        Return ``True`` if and only if this stream is known
        to be nonzero.

        An assumption of this class is that it is nonzero.

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamCauchyInverse, StreamFunction)
            sage: f = StreamFunction(lambda n: n^2, ZZ, False, 1)
            sage: g = StreamCauchyInverse(f)
            sage: g.is_nonzero()
            True
        """
        return True

class StreamMapCoefficients(StreamInexact):
    r"""
    The stream with ``function`` applied to each nonzero
    coefficient of ``series``.

    INPUT:

    - ``series`` -- a :class:`Stream`
    - ``function`` -- a function that modifies the elements of the stream
    - ``ring`` -- the base ring of the stream

    EXAMPLES::

        sage: from sage.data_structures.stream import (StreamMapCoefficients, StreamFunction)
        sage: f = StreamFunction(lambda n: 1, ZZ, True, 1)
        sage: g = StreamMapCoefficients(f, lambda n: -n, ZZ)
        sage: [g[i] for i in range(10)]
        [0, -1, -1, -1, -1, -1, -1, -1, -1, -1]
    """
    def __init__(self, series, function, ring):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.data_structures.stream import (StreamMapCoefficients, StreamFunction)
            sage: f = StreamFunction(lambda n: -1, ZZ, True, 0)
            sage: g = StreamMapCoefficients(f, lambda n: n + 1, ZZ)
            sage: TestSuite(g).run(skip="_test_pickling")
        """
        self._function = function
        self._ring = ring
        self._series = series
        super().__init__(series._is_sparse, series._approximate_order)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        INPUT:

        - ``n`` -- integer; the degree for the coefficient

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamMapCoefficients, StreamFunction)
            sage: f = StreamFunction(lambda n: n, ZZ, True, -1)
            sage: g = StreamMapCoefficients(f, lambda n: n^2 + 1, ZZ)
            sage: g.get_coefficient(5)
            26
            sage: [g.get_coefficient(i) for i in range(-1, 10)]
            [2, 0, 2, 5, 10, 17, 26, 37, 50, 65, 82]

            sage: R.<x,y> = ZZ[]
            sage: f = StreamFunction(lambda n: n, ZZ, True, -1)
            sage: g = StreamMapCoefficients(f, lambda n: n.degree() + 1, R)
            sage: [g.get_coefficient(i) for i in range(-1, 3)]
            [1, 0, 1, 1]

            sage: f = StreamFunction(lambda n: n, ZZ, True, 0)
            sage: g = StreamMapCoefficients(f, lambda n: 5, GF(3))
            sage: [g.get_coefficient(i) for i in range(10)]
            [0, 5, 5, 0, 5, 5, 0, 5, 5, 0]
        """
        c = self._ring(self._series[n])
        if c:
            return self._function(c)
        return c

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamMapCoefficients, StreamFunction)
            sage: f = StreamFunction(lambda n: -1, ZZ, True, 0)
            sage: g = StreamMapCoefficients(f, lambda n: n + 1, ZZ)
            sage: hash(g) == hash(g)
            True
        """
        # We don't hash the function as it might not be hashable.
        return hash((type(self), self._series, self._ring))

    def __eq__(self, other):
        """
        Test equality.

        INPUT:

        - ``other`` -- a stream of coefficients

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamMapCoefficients, StreamFunction)
            sage: f = StreamFunction(lambda n: -1, ZZ, True, 0)
            sage: def plus_one(n): return n + 1
            sage: g = StreamMapCoefficients(f, plus_one, ZZ)
            sage: g == f
            False
            sage: g == StreamMapCoefficients(f, plus_one, QQ)
            False
            sage: g == StreamMapCoefficients(f, plus_one, ZZ)
            True
        """
        return (isinstance(other, type(self)) and self._series == other._series
                and self._ring == other._ring and self._function == other._function)

class StreamShift(StreamInexact):
    """
    Operator for shifting the stream.

    INPUT:

    - ``series`` -- a :class:`Stream`
    - ``shift`` -- an integer
    """
    def __init__(self, series, shift):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamShift
            sage: from sage.data_structures.stream import StreamExact
            sage: h = StreamExact([1], False, constant=3)
            sage: M = StreamShift(h, 2)
            sage: TestSuite(M).run()
        """
        self._series = series
        self._shift = shift
        super().__init__(series._is_sparse, series._approximate_order + shift)

    def __getitem__(self, n):
        """
        Return the ``n``-th coefficient of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamShift
            sage: from sage.data_structures.stream import StreamFunction
            sage: F = StreamFunction(lambda n: n, ZZ, False, 1)
            sage: M = StreamShift(F, 2)
            sage: [F[i] for i in range(6)]
            [0, 1, 2, 3, 4, 5]
            sage: [M[i] for i in range(6)]
            [0, 0, 0, 1, 2, 3]
        """
        return self._series[n-self._shift]

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.data_structures.stream import StreamShift
            sage: from sage.data_structures.stream import StreamFunction
            sage: F = StreamFunction(lambda n: n, ZZ, False, 1)
            sage: M = StreamShift(F, 2)
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

            sage: from sage.data_structures.stream import StreamShift
            sage: from sage.data_structures.stream import StreamFunction
            sage: F = StreamFunction(lambda n: 1, ZZ, False, 1)
            sage: M2 = StreamShift(F, 2)
            sage: M3 = StreamShift(F, 3)
            sage: M2 == M3
            False
            sage: M2 == StreamShift(F, 2)
            True
        """
        return (isinstance(other, type(self)) and self._shift == other._shift
                and self._series == other._series)

    def is_nonzero(self):
        r"""
        Return ``True`` if and only if this stream is known
        to be nonzero.

        An assumption of this class is that it is nonzero.

        EXAMPLES::

            sage: from sage.data_structures.stream import (StreamCauchyInverse, StreamFunction)
            sage: f = StreamFunction(lambda n: n^2, ZZ, False, 1)
            sage: g = StreamCauchyInverse(f)
            sage: g.is_nonzero()
            True
        """
        return self._series.is_nonzero()

