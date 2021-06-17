from .infinity import infinity
from sage.structure.element import ModuleElement
from .integer_ring import ZZ

class LLS(ModuleElement):
    def __init__(self, parent, aux):
        ModuleElement.__init__(self, parent)
        self._aux = aux

    def __getitem__(self, n):
        return self.base_ring()(self._aux[n])

    def _mul_(self, other):
        P = self.parent()
        return P.element_class(P, LLS_mul(self._aux, other._aux))
    
    def _add_(self, other):
        P = self.parent()
        return P.element_class(P, LLS_add(self._aux, other._aux))

class LLS_aux():
    def __init__(self, is_sparse, approximate_valuation, constant=None):
        self._approximate_valuation = approximate_valuation
        self._constant = constant
        self._is_sparse = is_sparse

        if self._is_sparse:
            self._cache = dict()  # cache of known coefficients
        else:
            self._cache = list()
            self._offset = approximate_valuation
            self._iter = self.iterate_coefficients()

    def __getitem__(self, n):
        if self._approximate_valuation is infinity:
            return ZZ.zero()
        elif n < self._approximate_valuation:
            return ZZ.zero()
        elif self._constant is not None and n >= self._constant[1]:
            return self._constant[0]

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

class LLS_coefficient_function(LLS_aux):
    """

    EXAMPLES::

        sage: s = LLS_coefficient_function(lambda n: 1, True, 0)
        sage: [s[i] for i in range(-5, 5)]
        [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
        sage: t = LLS_mul(s, s)
        sage: [t[i] for i in range(-5, 5)]
        [0, 0, 0, 0, 0, 1, 2, 3, 4, 5]

    TESTS::

        sage: s = LLS_coefficient_function(lambda n: 1, False, 0)
        sage: t = LLS_mul(s, s)
        sage: [t[i] for i in range(-5, 5)]
        [0, 0, 0, 0, 0, 1, 2, 3, 4, 5]
   
        sage: s._cache
        [1, 1, 1, 1, 1]

    """
    def __init__(self, coefficient_function, is_sparse, approximate_valuation, constant=None):
        self._coefficient_function = coefficient_function
        super().__init__(is_sparse, approximate_valuation, constant)

    def get_coefficient(self, n):
        return self._coefficient_function(n)

    def iterate_coefficients(self):
        n = self._offset
        while True:
            yield self._coefficient_function(n)
            n += 1

           
class LLS_mul(LLS_aux):
    """
    Operator for multiplication.
    """
    def __init__(self, left, right):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z) * 1/(1 + z)
            sage: loads(dumps(f)) == f
            True
        """
        self._left = left
        self._right = right
        a = left._approximate_valuation + right._approximate_valuation
        c = None
        if left._constant is not None and right._constant is not None:
            if left._constant[0] == 0 and right._constant[0] == 0:
                c = (left._constant[0], left._constant[1] + right._constant[1] - 1)

        if left._is_sparse != right._is_sparse:
            raise NotImplementedError
        super().__init__(left._is_sparse, a, c)


    def get_coefficient(self, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = (1 + z)*(1 - z)
            sage: f.coefficient(2)
            -1
        """
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


class LLS_add(LLS_aux):
    """
    Operator for addition.
    """
    def __init__(self, left, right):
        """
        Initialize.

        TESTS::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = 1/(1 - z) * 1/(1 + z)
            sage: loads(dumps(f)) == f
            True
        """
        self._left = left
        self._right = right
        a = min(left._approximate_valuation, right._approximate_valuation)

        if left._constant is not None and right._constant is not None:
            c = (left._constant[0] + right._constant[0],
                 max(left._constant[1], right._constant[1]))
        else:
            c = None
        
        if left._is_sparse != right._is_sparse:
            raise NotImplementedError
        
        super().__init__(left._is_sparse, a, c)
    
    def get_coefficient(self, n):
        """
        Return the `n`-th coefficient of the series ``s``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: f = (1 + z)*(1 - z)
            sage: f.coefficient(2)
            -1
        """
        c = ZZ.zero()
        c = self._left[n] + self._right[n]
        return c

    def iterate_coefficients(self):
        n = self._offset
        while True:
            c = ZZ.zero()
            c = self._left[n] + self._right[n]
            yield c
            n += 1