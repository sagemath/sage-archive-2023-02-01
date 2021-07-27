r"""
Coefficient Stream

This module provides lazy class implementations of basic operators
on coefficient streams. The classes implemented in this module
can be used to build up more complex streams for different kinds of
series (Laurent, Dirichlet, etc).

EXAMPLES::

    The coefficient stream can be used to build up a Lazy laurent series::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ)
    sage: f = L(lambda n: n, True)
    sage: f
    z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + 7*z^7 + ...
    sage: type(f._coeff_stream)
    <class 'sage.data_structures.coefficient_stream.CoefficientStream_coefficient_function'>

    There are basic unary and binary operators available for the coefficient streams.
    For example, we can add two streams together::

    sage: g = L(lambda n: -n, True)
    sage: f + g
    0 + ...

    Coefficient streams can be subtracted::

    sage: f - g
    2*z + 4*z^2 + 6*z^3 + 8*z^4 + 10*z^5 + 12*z^6 + 14*z^7 + ...

    Coefficient streams can be multiplied::

    sage: g = L(lambda n: n^2, True)
    sage: f * g
    z^2 + 6*z^3 + 20*z^4 + 50*z^5 + 105*z^6 + 196*z^7 + 336*z^8 + ...

    Coefficient streams can be divided::

    sage: f / g
    1 - 2*z + 2*z^2 - 2*z^3 + 2*z^4 - 2*z^5 + 2*z^6 + ...

    Two coefficient streams can be composed (depending on whether it exists)::

    sage: f(g)
    z + 6*z^2 + 28*z^3 + 124*z^4 + 527*z^5 + 2172*z^6 + 8755*z^7 + ...

    We can also use the unary negation operator on a coefficient stream::

    sage: -f
    -z - 2*z^2 - 3*z^3 - 4*z^4 - 5*z^5 - 6*z^6 - 7*z^7 + ...

    Coefficient streams can be multiplied by a scalar::

    sage: f * 2
    2*z + 4*z^2 + 6*z^3 + 8*z^4 + 10*z^5 + 12*z^6 + 14*z^7 + ...

    The multiplicative inverse of a series can also be obtained::

    sage: f = L(lambda n: 1, True, 1); f
    z + z^2 + z^3 + z^4 + z^5 + z^6 + z^7 + ...
    sage: ~f
    z^-1 - 1 + ...

    Functions can also be applied to a coefficient stream::

    sage: f.map_coefficients(lambda n: n*3)
    3*z + 3*z^2 + 3*z^3 + 3*z^4 + 3*z^5 + 3*z^6 + 3*z^7 + ...

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


class CoefficientStream():
    """
    Abstract base class for all streams.

    INPUT:

    - ``sparse`` -- boolean; whether the implementation of the series is sparse

    - ``approximate_valuation`` -- the approximate valuation of the series

    EXAMPLES::

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
    CoefficientStream stream class when it is not or we do not know if it is
    eventually geometric.
    """

    def __init__(self, is_sparse, approximate_valuation):
        """
        Initialize the stream class for a CoefficientStream when it is not
        or it cannot be determined if it is eventually geometric.
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
        """
        self.__dict__ = d
        if not self._is_sparse:
            self._iter = self.iterate_coefficients()
            self._cache = []

    def __getitem__(self, n):
        """
        Return the `n`-th coefficient of the series.
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


class CoefficientStream_eventually_geometric(CoefficientStream):
    r"""
    A lazy Laurent series which is known to be eventually geometric

    INPUT:

        - ``laurent_polynomial`` -- a Laurent polynomial

        - ``is_sparse`` -- a boolean, which specifies whether the series is sparse

        - ``constant`` -- either ``None`` (default: ``None``) or pair of an element of the base ring and an integer

        - ``degree`` -- either ``None`` (default: ``None``) or an integer, the degree of the polynomial

    EXAMPLES::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ)
        sage: f = 1 + z^2
        sage: f[2]
        1

    You can do arithmetic operations with eventually geometric series::

        sage: g = z^3 + 1
        sage: s = f + g
        sage: s[2]
        1
        sage: s
        2 + z^2 + z^3
        sage: s = f - g
        sage: s[2]
        1
        sage: s
        z^2 - z^3
    """

    def __init__(self, laurent_polynomial, is_sparse, constant=None, degree=None):
        """
        Initialize a series that is known to be eventually geometric.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: M = z^2 + z
            sage: type(M._coeff_stream)
            <class 'sage.data_structures.coefficient_stream.CoefficientStream_eventually_geometric'>
        """
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
        """
        Return the coefficient of the term with exponent ``n`` of the series.

        INPUT:

        - ``n`` -- integer, the degree for which the coefficient is required

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: f = 1 + z + z^2 + z^3
            sage: f[3]
            1
            sage: f[8]
            0
        """
        if n >= self._degree:
            return self._constant
        return self._laurent_polynomial[n]

    def valuation(self):
        """
        Return the valuation of the series.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(QQ)
            sage: f = 1 + z + z^2 + z^3
            sage: f.valuation()
            0
        """
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
        return hash((self._laurent_polynomial, self._degree, self._constant))

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
                and self._laurent_polynomial == other._laurent_polynomial
                and self._constant == other._constant)


class CoefficientStream_coefficient_function(CoefficientStream_inexact):
    r"""
    The coefficient function of a lazy Laurent series.

    INPUT:

    - ``coefficient_function`` -- a python function that generates the coefficients of the series

    - ``ring`` -- the base ring of the series

    - ``is_sparse`` -- a boolean, which specifies whether the series is sparse

    - ``approximate_valuation`` -- the approximate valuation of the series

    EXAMPLES::

        sage: L.<z> = LazyLaurentSeriesRing(ZZ)
        sage: M = L(lambda n: n, True); M
        z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + 7*z^7 + ...
        sage: M[4]
        4
    """

    def __init__(self, coefficient_function, ring, is_sparse, approximate_valuation):
        """
        Initialize the coefficient function of a lazy Laurent series.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: type(L(lambda n:n, True)._coeff_stream)
            <class 'sage.data_structures.coefficient_stream.CoefficientStream_coefficient_function'>
        """
        self._coefficient_function = coefficient_function
        self._ring = ring
        super().__init__(is_sparse, approximate_valuation)

    def get_coefficient(self, n):
        """
        Return the coefficient of the term with exponent ``n`` of the series.

        INPUT:

        - ``n`` -- integer, the degree for which the coefficient is required

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: M = L(lambda n: n, True)
            sage: M._coeff_stream.get_coefficient(4)
            4
        """
        return self._ring(self._coefficient_function(n))

    def iterate_coefficients(self):
        """
        Return a generator for the coefficients of the series.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: M = L(lambda n: n, False)
            sage: n = M._coeff_stream.iterate_coefficients()
            sage: [next(n) for _ in range(10)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        n = self._offset
        ring = self._ring
        while True:
            yield ring(self._coefficient_function(n))
            n += 1


class CoefficientStream_uninitialized(CoefficientStream_inexact):
    r"""
    An uninitialized lazy Laurent series.

    INPUT:

    - ``is_sparse`` -- a boolean, which specifies whether the series is sparse

    - ``approximate_valuation`` -- the approximate valuation of the series

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
        Return the ``n``-th coefficient of the series.

        INPUT:

        - ``n`` -- integer, the degree for which the coefficient is required

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
        Return a generator for the coefficients of the series.

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
    """
    Abstract base class for unary operators.

    INPUT:

    - ``series`` -- series upon which the operator operates

    """

    def __init__(self, series, *args, **kwargs):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
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

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = ~(1 - z)
            sage: {f: 1}
            {1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...: 1}
        """
        return hash((type(self), self._series))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
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


class CoefficientStream_binary(CoefficientStream_inexact):

    def __init__(self, left, right, *args, **kwargs):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
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

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: {f: 1}
            {2 + 2*z^2 + 2*z^4 + 2*z^6 + ...: 1}
        """
        return hash((type(self), self._left, self._right))

    def __eq__(self, other):
        """
        Test equality.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z) + 1/(1 + z)
            sage: g = 1/(1 - z) + 1/(1 + z)
            sage: f == g
            True
        """
        if not isinstance(other, type(self)):
            return False
        return self._left == other._left and self._right == other._right


class CoefficientStream_binary_commutative(CoefficientStream_binary):

    def __hash__(self):
        """
        Return the hash of ``self``.
        """
        return hash((type(self), frozenset([self._left, self._right])))

    def __eq__(self, other):
        """
        Test the equality between ``self`` and ``other``.
        """
        if not isinstance(other, type(self)):
            return False
        if self._left == other._left and self._right == other._right:
            return True
        if self._left == other._right and self._right == other._left:
            return True
        return False


class CoefficientStream_zero(CoefficientStream):
    def __init__(self, sparse):
        """
        Initialise a lazy Laurent series which is known to be zero.
        """
        return super().__init__(sparse, 0)

    def __getitem__(self, n):
        """
        Return the ``n``-th coefficient of the series.
        """
        return ZZ.zero()

    def valuation(self):
        """
        Return the valuation of the series.
        """
        return infinity

    def __hash__(self):
        """
        Return the hash of ``self``.
        """
        return 0

#####################################################################
# Binary operations


class CoefficientStream_add(CoefficientStream_binary):
    """
    Operator for addition.
    """

    def __init__(self, left, right):
        """
        Initalize.
        """
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = min(left._approximate_valuation, right._approximate_valuation)
        super().__init__(left, right, left._is_sparse, a)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of the series when ``left`` is added to ``right``.
        """
        return self._left[n] + self._right[n]

    def iterate_coefficients(self):
        """
        Return a generator for the coefficients of the series when ``left`` is added to ``right``.
        """
        n = self._offset
        while True:
            yield self._left[n] + self._right[n]
            n += 1


class CoefficientStream_sub(CoefficientStream_binary):
    """
    Operator for subtraction.
    """

    def __init__(self, left, right):
        """
        Initialize.
        """
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = min(left._approximate_valuation, right._approximate_valuation)
        super().__init__(left, right, left._is_sparse, a)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of the series when ``right`` is subtracted from ``left``.
        """
        return self._left[n] - self._right[n]

    def iterate_coefficients(self):
        """
        Return the generator for the coefficients of the series when ``right`` is subtracted from ``left``.
        """
        n = self._offset
        while True:
            yield self._left[n] - self._right[n]
            n += 1


class CoefficientStream_mul(CoefficientStream_binary):
    """
    Operator for multiplication.

    We are assuming commutativity of the coefficient ring here.
    """

    def __init__(self, left, right):
        """
        Initialize.
        """
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError

        a = left._approximate_valuation + right._approximate_valuation
        super().__init__(left, right, left._is_sparse, a)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of the series when ``left`` is multiplied by ``right``.
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
        Return the generator for the coefficients of the series when ``left`` is multiplied by ``right``.
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


class CoefficientStream_div(CoefficientStream_binary):
    """
    Return ``left`` divided by ``right``.
    """

    def __init__(self, left, right):
        """
        Initialize.
        """
        lv = left.valuation()
        rv = right.valuation()
        self._lv = lv
        self._rv = rv
        self._ainv = ~right[rv]
        super().__init__(left, right, left._is_sparse, lv - rv)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of the series when ``left`` is divided by ``right``.
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
        Return the generator for the coefficients of the series when ``left`` is divided by ``right``.
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
    """

    def __init__(self, f, g):
        """
        Initialize.
        """
        assert g._approximate_valuation > 0
        self._fv = f._approximate_valuation
        self._gv = g._approximate_valuation
        if self._fv < 0:
            ginv = CoefficientStream_inv(g)
            # the constant part makes no contribution to the negative
            # we need this for the case so self._neg_powers[0][n] => 0
            self._neg_powers = [CoefficientStream_zero(f._is_sparse), ginv]
            for i in range(1, -self._fv):
                self._neg_powers.append(CoefficientStream_mul(self._neg_powers[-1], ginv))
        # Placeholder None to make this 1-based
        self._pos_powers = [None, g]
        val = self._fv * self._gv
        super().__init__(f, g, f._is_sparse, val)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of the series when ``f`` is composed by ``g``.
        """
        if n < 0:
            return sum(self._left[i] * self._neg_powers[-i][n] for i in range(self._fv, n // self._gv + 1))
        # n > 0
        while len(self._pos_powers) <= n // self._gv:
            self._pos_powers.append(CoefficientStream_mul(self._pos_powers[-1], self._right))
        ret = sum(self._left[i] * self._neg_powers[-i][n] for i in range(self._fv, 0))
        if n == 0:
            ret += self._left[0]
        return ret + sum(self._left[i] * self._pos_powers[i][n] for i in range(1, n // self._gv+1))

    def iterate_coefficients(self):
        """
        Return the generator for the coefficients of the series when ``f`` is composed by ``g``.
        """
        n = self._approximate_valuation
        while True:
            yield self.get_coefficient(n)
            n += 1

#####################################################################
# Unary operations


class CoefficientStream_scalar(CoefficientStream_unary):
    """
    Operator for multiplying with a scalar.
    """

    def __init__(self, series, scalar):
        """
        Initialize.
        """
        self._scalar = scalar

        super().__init__(series, series._is_sparse, series._approximate_valuation)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of the ``series`` when multiplied by the ``scalar``.
        """
        return self._series[n] * self._scalar

    def iterate_coefficients(self):
        """
        Return the generator for the coefficients of the ``series`` when multiplied by the ``scalar``.
        """
        n = self._offset
        while True:
            yield self._series[n] * self._scalar
            n += 1


class CoefficientStream_neg(CoefficientStream_unary):
    """
    Operator for negative of the series.
    """

    def __init__(self, series):
        """
        Initialize.
        """
        super().__init__(series, series._is_sparse, series._approximate_valuation)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of the ``series`` when negated.
        """
        return -self._series[n]

    def iterate_coefficients(self):
        """
        Return the generator for the coefficients of the ``series`` when negated.
        """
        n = self._offset
        while True:
            yield -self._series[n]
            n += 1


class CoefficientStream_inv(CoefficientStream_unary):
    """
    Operator for multiplicative inverse of the series.
    """

    def __init__(self, series):
        """
        Initialize.
        """
        v = series.valuation()
        super().__init__(series, series._is_sparse, -v)

        self._ainv = ~series[v]
        self._zero = ZZ.zero()

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of the multiplicative inverse of the ``series``.
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
        Return the generator for the coefficients of the multiplicative inverse of the ``series``.
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
    Return the series with ``function`` applied to each coefficient of this series.
    """

    def __init__(self, series, function, ring):
        """
        Initialize.
        """
        self._function = function
        self._ring = ring
        super().__init__(series, series._is_sparse, series._approximate_valuation)

    def get_coefficient(self, n):
        """
        Return the ``n``-th coefficient of the series with ``function`` applied to each coefficient.
        """
        c = self._ring(self._function(self._series[n]))
        return c

    def iterate_coefficients(self):
        """
        Return the generator for the coefficients of the series with ``function`` applied to each coefficient.
        """
        n = self._offset
        while True:
            c = self._ring(self._function(self._series[n]))
            yield c
            n += 1
