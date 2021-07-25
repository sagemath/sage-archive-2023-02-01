from .integer_ring import ZZ
from .infinity import infinity

class CoefficientStream():
    """
    Abstract base class for all auxillary LazyLaurentSeries.
    """

    def __init__(self, sparse, approximate_valuation):
        """
        Initialize the auxillary class for a LazyLaurentSeries.
        """
        self._is_sparse = sparse
        self._approximate_valuation = approximate_valuation


class LazyLaurentSeries_inexact(CoefficientStream):
    """
    LazyLaurentSeries aux class when it is not or we do not know if it is
    eventually geometric.
    """

    def __init__(self, is_sparse, approximate_valuation):
        """
        Initialize the auxillary class for a LazyLaurentSeries when it is not or it cannot be determined if it is eventually geometric.
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


class LazyLaurentSeries_unary(LazyLaurentSeries_inexact):
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


class LazyLaurentSeries_binary(LazyLaurentSeries_inexact):

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


class LazyLaurentSeries_binary_commutative(LazyLaurentSeries_binary):

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


class LazyLaurentSeries_zero(CoefficientStream):
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


class CoefficientStream_add(LazyLaurentSeries_binary):
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


class CoefficientStream_sub(LazyLaurentSeries_binary):
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


class CoefficientStream_mul(LazyLaurentSeries_binary):
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


class CoefficientStream_div(LazyLaurentSeries_binary):
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


class CoefficientStream_composition(LazyLaurentSeries_binary):
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
            self._neg_powers = [LazyLaurentSeries_zero(f._is_sparse), ginv]
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


class CoefficientStream_scalar(LazyLaurentSeries_unary):
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


class CoefficientStream_neg(LazyLaurentSeries_unary):
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


class CoefficientStream_inv(LazyLaurentSeries_unary):
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


class CoefficientStream_apply_coeff(LazyLaurentSeries_unary):
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