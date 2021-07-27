r"""
Lazy Laurent Series Rings

The ring of lazy Laurent series over a ring has usual arithmetic operations,
but it is actually not a ring in the usual sense since every
arithmetic operation gives a new series.

EXAMPLES:

The definition of Laurent series rings is not initially imported into the
global namespace. You need to import it explicitly to use it::

    sage: L.<z> = LazyLaurentSeriesRing(QQ)
    sage: L.category()
    Category of magmas and additive magmas
    sage: 1/(1 - z)
    1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
    sage: 1/(1 - z) == 1/(1 - z)
    True

Lazy Laurent series ring over a finite field::

    sage: L.<z> = LazyLaurentSeriesRing(GF(3)); L
    Lazy Laurent Series Ring in z over Finite Field of size 3
    sage: e = 1/(1 + z)
    sage: e.coefficient(100)
    1
    sage: e.coefficient(100).parent()
    Finite Field of size 3

Generating functions of integer sequences are Laurent series over the integer
ring::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ); L
    Lazy Laurent Series Ring in z over Integer Ring
    sage: 1/(1 - 2*z)^3
    1 + 6*z + 24*z^2 + 80*z^3 + 240*z^4 + 672*z^5 + 1792*z^6 + ...

Power series can be defined recursively::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
    sage: L._sparse
    True
    sage: s = L(None)
    sage: s._coeff_stream._is_sparse
    True
    sage: s._coeff_stream._approximate_valuation
    0
    sage: s.define(1 + z*s^2)
    sage: s
    1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + ...

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
import random

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

from sage.categories.magmas_and_additive_magmas import MagmasAndAdditiveMagmas
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets

from sage.misc.cachefunc import cached_method

from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from .integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing, LaurentPolynomialRing_generic

from .lazy_laurent_series import LazyLaurentSeries

from sage.data_structures.coefficient_stream import (
    CoefficientStream_zero,
    CoefficientStream_coefficient_function,
    CoefficientStream_eventually_geometric,
    CoefficientStream_uninitialized
)


class LazyLaurentSeriesRing(UniqueRepresentation, Parent):
    """
    Lazy Laurent series ring.

    INPUT:

    - ``base_ring`` -- base ring of this Laurent series ring

    - ``names`` -- name of the generator of this Laurent series ring

    EXAMPLES::

        sage: LazyLaurentSeriesRing(ZZ, 't')
        Lazy Laurent Series Ring in t over Integer Ring
    """
    Element = LazyLaurentSeries

    def __init__(self, base_ring, names, sparse=False, category=None):
        """
        Initialize the ring.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 't')
            sage: TestSuite(L).run(skip='_test_elements')
        """
        self._sparse = sparse
        self._laurent_poly_ring = LaurentPolynomialRing(base_ring, names, sparse=sparse)
        Parent.__init__(self, base=base_ring, names=names,
                        category=MagmasAndAdditiveMagmas().or_subcategory(category))

    def _repr_(self):
        """
        String representation of this Laurent series ring.

        EXAMPLES::

            sage: LazyLaurentSeriesRing(GF(2), 'z')
            Lazy Laurent Series Ring in z over Finite Field of size 2
        """
        return "Lazy Laurent Series Ring in {} over {}".format(self.variable_name(), self.base_ring())

    @cached_method
    def gen(self, n=0):
        """
        Return the generator of this Laurent series ring.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.gen()
            z
            sage: L.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: there is only one generator
        """
        if n != 0:
            raise IndexError("there is only one generator")
        R = self._laurent_poly_ring
        coeff_stream = CoefficientStream_eventually_geometric(R.gen(n), self._sparse, ZZ.zero(), 2)
        return self.element_class(self, coeff_stream)

    def ngens(self):
        """
        Return the number of generators of this Laurent series ring.

        This is always 1.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: L.ngens()
            1
        """
        return 1

    @cached_method
    def gens(self):
        """
        Return the tuple of the generator.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: L.gens()
            (z,)
            sage: 1/(1 - z)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
        """
        return tuple([self.gen(n) for n in range(self.ngens())])

    def _coerce_map_from_(self, S):
        """
        Return ``True`` if a coercion from ``S`` exists.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L.has_coerce_map_from(ZZ)
            True
            sage: L.has_coerce_map_from(GF(2))
            True
        """
        if self.base_ring().has_coerce_map_from(S):
            return True

        R = self._laurent_poly_ring
        if R.has_coerce_map_from(S):
            def make_series_from(poly):
                return self.element_class(self, CoefficientStream_eventually_geometric(R(poly), self._sparse))
            return SetMorphism(Hom(S, self, Sets()), make_series_from)

        return False

    def _element_constructor_(self, x=None, valuation=None, constant=None, degree=None):
        """
        Construct a Laurent series from ``x``.

        INPUT:

        - ``x`` -- a Laurent series, a Laurent polynomial, a Python function, or a list of elements in the base ring

        - ``valuation`` -- an integer or ``None`` (default: ``None``), the approximate valuation of the series

        - ``constant`` -- either ``None`` (default: ``None``) or pair of an element of the base ring and an integer

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L(2)
            0
            sage: L(3)
            1

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')

            sage: L(lambda i: i, 5, 1, 10)
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + ...
            sage: L(lambda i: i, 5, (1, 10))
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + ...

            sage: X = L(constant=5, degree=2); X
            5*z^2 + 5*z^3 + 5*z^4 + ...
            sage: X.valuation()
            2

            sage: def g(i):
            ....:     if i < 0:
            ....:         return 1
            ....:     else:
            ....:         return 1 + sum(k for k in range(i+1))
            sage: e = L(g, -5); e
            z^-5 + z^-4 + z^-3 + z^-2 + z^-1 + 1 + 2*z + ...
            sage: f = e^-1; f
            z^5 - z^6 - z^11 + ...
            sage: f.coefficient(10)
            0
            sage: f[20]
            9
            sage: f[30]
            -219

            sage: L(valuation=2, constant=1)
            z^2 + z^3 + z^4 + ...
            sage: L(constant=1)
            Traceback (most recent call last):
            ...
            ValueError: you must specify the degree for the polynomial 0

        Alternatively, ``x`` can be a list of elements of the base ring.
        Then these elements are read as coefficients of the terms of
        degrees starting from the ``valuation``. In this case, ``constant``
        may be just an element of the base ring instead of a tuple or can be
        simply omitted if it is zero::

            sage: f = L([1,2,3,4], -5)
            sage: f
            z^-5 + 2*z^-4 + 3*z^-3 + 4*z^-2
            sage: g = L([1,3,5,7,9], 5, -1)
            sage: g
            z^5 + 3*z^6 + 5*z^7 + 7*z^8 + 9*z^9 - z^10 - z^11 - z^12 + ...

        TODO::

            Add a method to make a copy of self._sparse.
        """
        if x is None:
            if valuation is None:
                valuation = 0
            return self.element_class(self, CoefficientStream_uninitialized(self._sparse, valuation))

        R = self._laurent_poly_ring
        BR = self.base_ring()
        try:
            # Try to build stuff using the polynomial ring constructor
            x = R(x)
        except (TypeError, ValueError):
            pass
        if isinstance(constant, (tuple, list)):
            constant, degree = constant
        if constant is not None:
            constant = BR(constant)
        if x in R:
            if not x and not constant:
                coeff_stream = CoefficientStream_zero(self._sparse)
            else:
                if x and valuation:
                    x = x.shift(valuation - x.valuation())
                if degree is None and not x:
                    degree = valuation
                coeff_stream = CoefficientStream_eventually_geometric(R(x), self._sparse, constant, degree)
            return self.element_class(self, coeff_stream)
        if isinstance(x, LazyLaurentSeries):
            if x._coeff_stream._is_sparse is self._sparse:
                return self.element_class(self, x._coeff_stream)
            # TODO: Implement a way to make a self._sparse copy
            raise NotImplementedError("cannot convert between sparse and dense")
        if callable(x):
            if valuation is None:
                valuation = 0
            if degree is not None:
                if constant is None:
                    constant = ZZ.zero()
                z = R.gen()
                p = R.sum(x(i) * z**i for i in range(valuation, degree))
                return self.element_class(self, CoefficientStream_eventually_geometric(p, self._sparse, constant, degree))
            return self.element_class(self, CoefficientStream_coefficient_function(x, self.base_ring(), self._sparse, valuation))
        raise ValueError(f"unable to convert {x} into a lazy Laurent series")

    def _an_element_(self):
        """
        Return a Laurent series in this ring.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.an_element()
            z^-10 + z^-9 + z^-8 + ...
        """
        c = self.base_ring().an_element()
        R = self._laurent_poly_ring
        return self.element_class(self, CoefficientStream_eventually_geometric(R.zero(), self._sparse, c, -10))

    @cached_method
    def one(self):
        """
        Return the constant series `1`.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.one()
            1
        """
        R = self._laurent_poly_ring
        return self.element_class(self, CoefficientStream_eventually_geometric(R.one(), self._sparse, ZZ.zero(), 1))

    @cached_method
    def zero(self):
        """
        Return the zero series.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.zero()
            0
        """
        return self.element_class(self, CoefficientStream_zero(self._sparse))
