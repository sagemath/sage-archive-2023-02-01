r"""
Lazy Laurent Series Rings

The ring of Laurent series over a ring with the usual arithmetic where the
coefficients are computed lazily.

EXAMPLES::

    sage: L.<z> = LazyLaurentSeriesRing(QQ)
    sage: L.category()
    Category of magmas and additive magmas
    sage: 1/(1 - z)
    1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
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

Generating functions of integer sequences are Laurent series over the
integer ring::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ); L
    Lazy Laurent Series Ring in z over Integer Ring
    sage: 1/(1 - 2*z)^3
    1 + 6*z + 24*z^2 + 80*z^3 + 240*z^4 + 672*z^5 + 1792*z^6 + O(z^7)

Power series can be defined recursively::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ)
    sage: s = L(None)
    sage: s.define(1 + z*s^2)
    sage: s
    1 + z + 2*z^2 + 5*z^3 + 14*z^4 + 42*z^5 + 132*z^6 + O(z^7)

The implementation of the ring can be either be a sparse or a dense one.
The default is a dense implementation::

    sage: L.<z> = LazyLaurentSeriesRing(ZZ)
    sage: L.is_sparse()
    False
    sage: L.<z> = LazyLaurentSeriesRing(ZZ, sparse=True)
    sage: L.is_sparse()
    True

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

from inspect import signature

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

from sage.categories.magmas_and_additive_magmas import MagmasAndAdditiveMagmas
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets

from sage.misc.cachefunc import cached_method

from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from .integer_ring import ZZ
from .infinity import infinity
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing, LaurentPolynomialRing_generic
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.global_options import GlobalOptions
from sage.symbolic.ring import SR

from .lazy_laurent_series import LazyLaurentSeries, LazyDirichletSeries, LazyTaylorSeries

from sage.data_structures.coefficient_stream import (
    CoefficientStream_zero,
    CoefficientStream_coefficient_function,
    CoefficientStream_exact,
    CoefficientStream_uninitialized
)


class LazyLaurentSeriesRing(UniqueRepresentation, Parent):
    """
    The ring of lazy Laurent series.

    INPUT:

    - ``base_ring`` -- base ring
    - ``names`` -- name of the generator
    - ``sparse`` -- (default: ``False``) whether the implementation of
      the series is sparse or not

    EXAMPLES::

        sage: LazyLaurentSeriesRing(ZZ, 't')
        Lazy Laurent Series Ring in t over Integer Ring
    """
    Element = LazyLaurentSeries

    def __init__(self, base_ring, names, sparse=False, category=None):
        """
        Initialize ``self``.

        TESTS::

            sage: L = LazyLaurentSeriesRing(ZZ, 't')
            sage: TestSuite(L).run(skip='_test_elements')
        """
        self._sparse = sparse
        self._coeff_ring = base_ring
        self._laurent_poly_ring = LaurentPolynomialRing(base_ring, names, sparse=sparse)
        Parent.__init__(self, base=base_ring, names=names,
                        category=MagmasAndAdditiveMagmas().or_subcategory(category))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LazyLaurentSeriesRing(GF(2), 'z')
            Lazy Laurent Series Ring in z over Finite Field of size 2
        """
        return "Lazy Laurent Series Ring in {} over {}".format(self.variable_name(), self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: latex(L)
            \Bold{F}_{2} [\![z]\!]
        """
        from sage.misc.latex import latex
        return latex(self.base_ring()) + r"[\![{}]\!]".format(self.variable_name())

    def is_sparse(self):
        """
        Return if ``self`` is sparse or not.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.is_sparse()
            False

            sage: L = LazyLaurentSeriesRing(ZZ, 'z', sparse=True)
            sage: L.is_sparse()
            True
        """
        return self._sparse

    @cached_method
    def monomial(self, c, n):
        r"""
        Return the interpretation of the coefficient ``c`` at index ``n``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.monomial(1, 3)
            z^3
            sage: L.monomial(2, -4)
            2*z^-4
        """
        L = self._laurent_poly_ring
        return L(c) * L.gen() ** n

    @cached_method
    def gen(self, n=0):
        r"""
        Return the ``n``-th generator of ``self``.

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
        R = self.base_ring()
        coeff_stream = CoefficientStream_exact([R.one()], self._sparse, constant=ZZ.zero(), valuation=1)
        return self.element_class(self, coeff_stream)

    def ngens(self):
        r"""
        Return the number of generators of ``self``.

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
        Return the generators of ``self``.

        EXAMPLES::

            sage: L.<z> = LazyLaurentSeriesRing(ZZ)
            sage: L.gens()
            (z,)
            sage: 1/(1 - z)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
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
        return R.has_coerce_map_from(S)

    def _element_constructor_(self, x=None, valuation=None, constant=None, degree=None):
        """
        Construct a Laurent series from ``x``.

        INPUT:

        - ``x`` -- data used to the define a Laurent series
        - ``valuation`` -- integer (optional); integer; a lower bound for
          the valuation of the series
        - ``constant`` -- (optional) the eventual constant of the series
        - ``degree`` -- (optional) the degree when the series is ``constant``

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(GF(2), 'z')
            sage: L(2)
            0
            sage: L(3)
            1

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')

            sage: L(lambda i: i, 5, 1, 10)
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)
            sage: L(lambda i: i, 5, (1, 10))
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)

            sage: X = L(constant=5, degree=2); X
            5*z^2 + 5*z^3 + 5*z^4 + O(z^5)
            sage: X.valuation()
            2

            sage: def g(i):
            ....:     if i < 0:
            ....:         return 1
            ....:     else:
            ....:         return 1 + sum(k for k in range(i+1))
            sage: e = L(g, -5); e
            z^-5 + z^-4 + z^-3 + z^-2 + z^-1 + 1 + 2*z + O(z^2)
            sage: f = e^-1; f
            z^5 - z^6 - z^11 + O(z^12)
            sage: f.coefficient(10)
            0
            sage: f[20]
            9
            sage: f[30]
            -219

            sage: L(valuation=2, constant=1)
            z^2 + z^3 + z^4 + O(z^5)
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
            z^5 + 3*z^6 + 5*z^7 + 7*z^8 + 9*z^9 - z^10 - z^11 - z^12 + O(z^13)

        .. TODO::

            Add a method to change the sparse/dense implementation.
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
                    if valuation is None:
                        raise ValueError("you must specify the degree for the polynomial 0")
                    degree = valuation
                if x == R.zero():
                    coeff_stream = CoefficientStream_exact([x], self._sparse, valuation=degree-1, constant=constant)
                    return self.element_class(self, coeff_stream)
                initial_coefficients = [x[i] for i in range(x.valuation(), x.degree() + 1)]
                coeff_stream = CoefficientStream_exact(initial_coefficients, self._sparse, valuation=x.valuation(), constant=constant, degree=degree)
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
                p = [x(i) for i in range(valuation, degree)]
                coeff_stream = CoefficientStream_exact(p, self._sparse, valuation=valuation, constant=constant, degree=degree)
                return self.element_class(self, coeff_stream)
            return self.element_class(self, CoefficientStream_coefficient_function(x, self.base_ring(), self._sparse, valuation))
        raise ValueError(f"unable to convert {x} into a lazy Laurent series")

    def _an_element_(self):
        """
        Return a Laurent series in ``self``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.an_element()
            z^-2 + 3*z^-1 + 2*z + z^2 + z^3 + z^4 + z^5 + O(z^6)
        """
        R = self.base_ring()
        coeff_stream = CoefficientStream_exact([R.an_element(), 3, 0, 2*R.an_element(), 1],
                                               self._sparse, valuation=-2, constant=1)
        return self.element_class(self, coeff_stream)

    @cached_method
    def one(self):
        r"""
        Return the constant series `1`.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.one()
            1
        """
        R = self.base_ring()
        coeff_stream = CoefficientStream_exact([R.one()], self._sparse, constant=ZZ.zero(), degree=1)
        return self.element_class(self, coeff_stream)

    @cached_method
    def zero(self):
        r"""
        Return the zero series.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.zero()
            0
        """
        return self.element_class(self, CoefficientStream_zero(self._sparse))

    # add options to class
    class options(GlobalOptions):
        r"""
        Set and display the options for Lazy Laurent series.

        If no parameters are set, then the function returns a copy of
        the options dictionary.

        The ``options`` to Lazy Laurent series can be accessed as the method
        :meth:`LazyLaurentSeriesRing.options` of :class:`LazyLaurentSeriesRing`.

        @OPTIONS@

        EXAMPLES::

            sage: LLS.<z> = LazyLaurentSeriesRing(QQ)
            sage: LLS.options.display_length
            7
            sage: f = 1/(1-z)
            sage: f
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + O(z^7)
            sage: LLS.options.display_length = 10
            sage: f
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + z^7 + z^8 + z^9 + O(z^10)
            sage: g = LLS(lambda n: n^2, valuation=-2, degree=5, constant=42)
            sage: g
            4*z^-2 + z^-1 + z + 4*z^2 + 9*z^3 + 16*z^4 + 42*z^5 + 42*z^6 + 42*z^7 + O(z^8)
            sage: LLS.options.constant_length = 1
            sage: g
            4*z^-2 + z^-1 + z + 4*z^2 + 9*z^3 + 16*z^4 + 42*z^5 + O(z^6)
            sage: LazyLaurentSeriesRing.options._reset()
            sage: LazyLaurentSeriesRing.options.display_length
            7
        """
        NAME = 'LazyLaurentSeriesRing'
        module = 'sage.rings.lazy_laurent_series_ring'
        display_length = dict(default=7,
                              description='the number of coefficients to display from the valuation',
                              checker=lambda x: x in ZZ and x > 0)
        constant_length = dict(default=3,
                               description='the number of coefficients to display for nonzero constant series',
                               checker=lambda x: x in ZZ and x > 0)

    def series(self, coefficient, valuation, constant=None, degree=None):
        r"""
        Return a lazy Laurent series.

        INPUT:

        - ``coefficient`` -- Python function that computes coefficients or a list
        - ``valuation`` -- integer; approximate valuation of the series
        - ``constant`` -- (optional) an element of the base ring or a
          pair of an element of the base ring and an integer
        - ``degree`` -- (optional) integer

        Let the coefficient of index `i` mean the coefficient of the term
        of the series with exponent `i`.

        Python function ``coefficient`` returns the value of the coefficient
        of index `i` from input `s` and `i` where `s` is the series itself.

        Let ``valuation`` be `n`. All coefficients of index below `n` are zero.
        If ``constant`` is not specified, then the ``coefficient`` function is
        responsible to compute the values of all coefficients of index `\ge n`.
        If ``constant`` is a pair `(c,m)`, then the ``coefficient`` function
        is responsible to compute the values of all coefficients of index
        `\ge n` and `< m` and all the coefficients of index `\ge m` is
        the constant `c`.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.series(lambda s, i: i, 5, (1,10))
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)

            sage: def g(s, i):
            ....:     if i < 0:
            ....:         return 1
            ....:     else:
            ....:         return s.coefficient(i - 1) + i
            sage: e = L.series(g, -5); e
            z^-5 + z^-4 + z^-3 + z^-2 + z^-1 + 1 + 2*z + ...
            sage: f = e^-1; f
            z^5 - z^6 - z^11 + ...
            sage: f.coefficient(10)
            0
            sage: f.coefficient(20)
            9
            sage: f.coefficient(30)
            -219

        Alternatively, the ``coefficient`` can be a list of elements of the
        base ring. Then these elements are read as coefficients of the terms of
        degrees starting from the ``valuation``. In this case, ``constant``
        may be just an element of the base ring instead of a tuple or can be
        simply omitted if it is zero. ::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: f = L.series([1,2,3,4], -5); f
            z^-5 + 2*z^-4 + 3*z^-3 + 4*z^-2
            sage: g = L.series([1,3,5,7,9], 5, -1); g
            z^5 + 3*z^6 + 5*z^7 + 7*z^8 + 9*z^9 - z^10 - z^11 - z^12 + O(z^13)
        """
        if isinstance(constant, (list, tuple)):
            constant, degree = constant

        if isinstance(coefficient, (tuple, list)):
            if constant is None:
                constant = self.base_ring().zero()
            if degree is None:
                degree = valuation + len(coefficient)
            coeff_stream = CoefficientStream_exact(coefficient, self._sparse, valuation=valuation,
                                                   constant=constant, degree=degree)
            return self.element_class(self, coeff_stream)

        if degree is not None and valuation > degree and constant:
            raise ValueError('inappropriate valuation')

        t = None
        t = self(lambda n: coefficient(t, n), valuation=valuation,
                 constant=constant, degree=degree)
        return t

######################################################################

class LazyTaylorSeriesRing(UniqueRepresentation, Parent):
    """
    Lazy Taylor series ring.

    INPUT:

    - ``base_ring`` -- base ring of this Taylor series ring
    - ``names`` -- name(s) of the generator of this Taylor series ring
    - ``sparse`` -- (default: ``False``) whether this series is sparse or not

    EXAMPLES::

        sage: LazyTaylorSeriesRing(ZZ, 't')
        Lazy Taylor Series Ring in t over Integer Ring

        sage: L.<x, y> = LazyTaylorSeriesRing(QQ); L
        Multivariate Lazy Taylor Series Ring in x, y over Rational Field
    """
    Element = LazyTaylorSeries

    def __init__(self, base_ring, names, sparse=False, category=None):
        """
        Initialize ``self``.

        TESTS::

            sage: L = LazyTaylorSeriesRing(ZZ, 't')
            sage: TestSuite(L).run(skip='_test_elements')
        """
        self._sparse = sparse
        if len(names) == 1:
            self._coeff_ring = base_ring
        else:
            self._coeff_ring = PolynomialRing(base_ring, names)
        self._laurent_poly_ring = PolynomialRing(base_ring, names)
        Parent.__init__(self, base=base_ring, names=names,
                        category=MagmasAndAdditiveMagmas().or_subcategory(category))

    def _repr_(self):
        """
        String representation of this Taylor series ring.

        EXAMPLES::

            sage: LazyTaylorSeriesRing(GF(2), 'z')
            Lazy Taylor Series Ring in z over Finite Field of size 2
        """
        if len(self.variable_names()) == 1:
            return "Lazy Taylor Series Ring in {} over {}".format(self.variable_name(), self.base_ring())
        generators_rep = ", ".join(self.variable_names())
        return "Multivariate Lazy Taylor Series Ring in {} over {}".format(generators_rep, self.base_ring())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: L = LazyTaylorSeriesRing(GF(2), 'z')
            sage: latex(L)
            \Bold{F}_{2} [\![z]\!]
        """
        from sage.misc.latex import latex
        generators_rep = ", ".join(self.variable_names())
        return latex(self.base_ring()) + r"[\![{}]\!]".format(generators_rep)

    @cached_method
    def monomial(self, c, n):
        r"""
        Return the interpretation of the coefficient ``c`` at index ``n``.

        EXAMPLES::

            sage: L = LazyLaurentSeriesRing(ZZ, 'z')
            sage: L.monomial(2, 3)
            2*z^3
        """
        m = len(self.variable_names())
        L = self._laurent_poly_ring
        if m == 1:
            return L(c) * L.gen() ** n
        return L(c)

    @cached_method
    def gen(self, n=0):
        """
        Return the ``n``-th generator of ``self``.

        EXAMPLES::

            sage: L = LazyTaylorSeriesRing(ZZ, 'z')
            sage: L.gen()
            z
            sage: L.gen(3)
            Traceback (most recent call last):
            ...
            IndexError: there is only one generator
        """
        m = len(self.variable_names())
        if n > m:
            if m == 1:
                raise IndexError("there is only one generator")
            raise IndexError("there are only %s generators" % m)

        R = self._laurent_poly_ring
        if len(self.variable_names()) == 1:
            coeff_stream = CoefficientStream_exact([1], self._sparse, constant=ZZ.zero(), valuation=1)
        else:
            coeff_stream = CoefficientStream_exact([R.gen(n)], self._sparse, constant=ZZ.zero(), valuation=1)
        return self.element_class(self, coeff_stream)

    def ngens(self):
        r"""
        Return the number of generators of ``self``.

        This is always 1.

        EXAMPLES::

            sage: L.<z> = LazyTaylorSeriesRing(ZZ)
            sage: L.ngens()
            1
        """
        return len(self.variable_names())

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: L.<z> = LazyTaylorSeriesRing(ZZ)
            sage: L.gens()
            (z,)
            sage: (1+z)^2
            1 + 2*z + z^2
            sage: 1/(1 - z)
            1 + z + z^2 + z^3 + z^4 + z^5 + z^6 + ...
        """
        return tuple([self.gen(n) for n in range(self.ngens())])

    def _coerce_map_from_(self, S):
        """
        Return ``True`` if a coercion from ``S`` exists.

        EXAMPLES::

            sage: L = LazyTaylorSeriesRing(GF(2), 'z')
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
                p_dict = poly.homogeneous_components()
                v = min(p_dict.keys())
                d = max(p_dict.keys())
                p_list = [p_dict.get(i, 0) for i in range(v, d + 1)]
                coeff_stream = CoefficientStream_exact(p_list, self._sparse, valuation=v)
                return self.element_class(self, coeff_stream)
            return SetMorphism(Hom(S, self, Sets()), make_series_from)

        return False

    def _element_constructor_(self, x=None, valuation=None, constant=None, degree=None):
        """
        Construct a Taylor series from ``x``.

        INPUT:

        - ``x`` -- data used to the define a Taylor series
        - ``valuation`` -- integer (optional); integer; a lower bound for the valuation of the series
        - ``constant`` -- (optional) the eventual constant of the series
        - ``degree`` -- (optional) the degree when the series is ``constant``

        EXAMPLES::

            sage: L = LazyTaylorSeriesRing(GF(2), 'z')
            sage: L(2)
            0
            sage: L(3)
            1

            sage: L = LazyTaylorSeriesRing(ZZ, 'z')

            sage: L(lambda i: i, 5, 1, 10)
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)
            sage: L(lambda i: i, 5, (1, 10))
            5*z^5 + 6*z^6 + 7*z^7 + 8*z^8 + 9*z^9 + z^10 + z^11 + z^12 + O(z^13)

            sage: X = L(constant=5, degree=2); X
            5*z^2 + 5*z^3 + 5*z^4 + O(z^5)
            sage: X.valuation()
            2

            sage: e = L(lambda n: n+1); e
            1 + 2*z + 3*z^2 + 4*z^3 + 5*z^4 + 6*z^5 + 7*z^6 + O(z^7)
            sage: f = e^-1; f
            1 - 2*z + z^2 + O(z^7)
            sage: f.coefficient(10)
            0
            sage: f[20]
            0

            sage: L(valuation=2, constant=1)
            z^2 + z^3 + z^4 + O(z^5)
            sage: L(constant=1)
            1 + z + z^2 + O(z^3)

        Alternatively, ``x`` can be a list of elements of the base ring.
        Then these elements are read as coefficients of the terms of
        degrees starting from the ``valuation``. In this case, ``constant``
        may be just an element of the base ring instead of a tuple or can be
        simply omitted if it is zero::

            sage: f = L([1,2,3,4], 1); f
            z + 2*z^2 + 3*z^3 + 4*z^4

            sage: g = L([1,3,5,7,9], 5, -1); g
            z^5 + 3*z^6 + 5*z^7 + 7*z^8 + 9*z^9 - z^10 - z^11 - z^12 + O(z^13)

        .. TODO::

            Add a method to change the sparse/dense implementation.
        """
        if valuation is None:
            valuation = 0
        assert valuation >= 0, "the valuation of a Taylor series must be positive"

        R = self._laurent_poly_ring
        if x is None:
            assert degree is None
            coeff_stream = CoefficientStream_uninitialized(self._sparse, valuation)
            return self.element_class(self, coeff_stream)

        try:
            # Try to build stuff using the polynomial ring constructor
            x = R(x)
        except (TypeError, ValueError):
            pass
        if isinstance(constant, (tuple, list)):
            constant, degree = constant
        if constant is not None:
            constant = R(constant)
        if x in R:
            if not x and not constant:
                coeff_stream = CoefficientStream_zero(self._sparse)
            else:
                if not x:
                    coeff_stream = CoefficientStream_exact([], self._sparse,
                                                           valuation=valuation,
                                                           degree=degree,
                                                           constant=constant)
                    return self.element_class(self, coeff_stream)

                if len(self.variable_names()) == 1:
                    v = x.valuation()
                    d = x.degree()
                    p_list = [x[i] for i in range(v, d + 1)]
                else:
                    p_dict = x.homogeneous_components()
                    v = min(p_dict.keys())
                    d = max(p_dict.keys())
                    p_list = [p_dict.get(i, 0) for i in range(v, d + 1)]

                coeff_stream = CoefficientStream_exact(p_list, self._sparse,
                                                       valuation=valuation,
                                                       constant=constant,
                                                       degree=degree)
            return self.element_class(self, coeff_stream)

        if isinstance(x, LazyTaylorSeries):
            if x._coeff_stream._is_sparse is self._sparse:
                return self.element_class(self, x._coeff_stream)
            # TODO: Implement a way to make a self._sparse copy
            raise NotImplementedError("cannot convert between sparse and dense")
        if callable(x):
            if degree is not None:
                if constant is None:
                    constant = ZZ.zero()
                z = R.gen()
                p = [x(i) for i in range(valuation, degree)]
                coeff_stream = CoefficientStream_exact(p, self._sparse,
                                                       valuation=valuation,
                                                       constant=constant,
                                                       degree=degree)
                return self.element_class(self, coeff_stream)
            coeff_stream = CoefficientStream_coefficient_function(x, self._coeff_ring, self._sparse, valuation)
            return self.element_class(self, coeff_stream)
        raise ValueError(f"unable to convert {x} into a lazy Taylor series")

    def _an_element_(self):
        """
        Return a Taylor series in ``self``.

        EXAMPLES::

            sage: L = LazyTaylorSeriesRing(ZZ, 'z')
            sage: L.an_element()
            z + z^2 + z^3 + z^4 + ...
        """
        c = self.base_ring()(1)
        R = self._laurent_poly_ring
        coeff_stream = CoefficientStream_exact([R.one()], self._sparse, valuation=1, constant=c)
        return self.element_class(self, coeff_stream)

    @cached_method
    def one(self):
        r"""
        Return the constant series `1`.

        EXAMPLES::

            sage: L = LazyTaylorSeriesRing(ZZ, 'z')
            sage: L.one()
            1
        """
        R = self._laurent_poly_ring
        coeff_stream = CoefficientStream_exact([R.one()], self._sparse, constant=ZZ.zero(), degree=1)
        return self.element_class(self, coeff_stream)

    @cached_method
    def zero(self):
        r"""
        Return the zero series.

        EXAMPLES::

            sage: L = LazyTaylorSeriesRing(ZZ, 'z')
            sage: L.zero()
            0
        """
        return self.element_class(self, CoefficientStream_zero(self._sparse))

    # add options to class
    class options(GlobalOptions):
        NAME = 'LazyTaylorSeriesRing'
        module = 'sage.rings.lazy_laurent_series_ring'
        display_length = dict(default=7,
                              description='the number of coefficients to display from the valuation',
                              checker=lambda x: x in ZZ and x > 0)
        constant_length = dict(default=3,
                               description='the number of coefficients to display for nonzero constant series',
                               checker=lambda x: x in ZZ and x > 0)


######################################################################

class LazyDirichletSeriesRing(UniqueRepresentation, Parent):
    """
    Lazy Dirichlet series ring.

    INPUT:

    - ``base_ring`` -- base ring of this Dirichlet series ring

    - ``names`` -- name of the generator of this Dirichlet series ring

    EXAMPLES::

        sage: LazyDirichletSeriesRing(ZZ, 't')
        Lazy Dirichlet Series Ring in t over Integer Ring
    """
    Element = LazyDirichletSeries

    def __init__(self, base_ring, names, sparse=False, category=None):
        """
        Initialize the ring.

        TESTS::

            sage: L = LazyDirichletSeriesRing(ZZ, 't')
            sage: TestSuite(L).run(skip='_test_elements')
        """
        if base_ring.characteristic() > 0:
            raise ValueError("positive characteristic not allowed for Dirichlet series")

        self._sparse = sparse
        self._coeff_ring = base_ring
        self._laurent_poly_ring = SR
        Parent.__init__(self, base=base_ring, names=names,
                        category=MagmasAndAdditiveMagmas().or_subcategory(category))

    def _repr_(self):
        """
        String representation of this Dirichlet series ring.

        EXAMPLES::

            sage: LazyDirichletSeriesRing(QQbar, 'z')
            Lazy Dirichlet Series Ring in z over Algebraic Field
        """
        return "Lazy Dirichlet Series Ring in {} over {}".format(self.variable_name(), self.base_ring())

    @cached_method
    def monomial(self, c, n):
        r"""
        Return the interpretation of the coefficient ``c`` at index ``n``.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.monomial(5, 3)
            5/3^z
        """
        L = self._laurent_poly_ring
        return L(c) * L(n) ** -L(self.variable_name())

    @cached_method
    def gen(self, n=1):
        """
        Return the `n`-th generator of this Dirichlet series ring.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.gen()
            1
            sage: L.gen(4)
            1/(4^z)

        """
        assert n >= 1
        coeff_stream = CoefficientStream_exact([1], self._sparse, valuation=n, constant=ZZ.zero())
        return self.element_class(self, coeff_stream)

    def ngens(self):
        """
        Return the number of generators of this Dirichlet series ring.

        This is always `\infty`.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, "s")
            sage: L.ngens()
            +Infinity
        """
        return infinity

    @cached_method
    def gens(self):
        """
        Return the tuple of the generator.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, "z")
            sage: L.gens()
            Lazy family (<lambda>(i))_{i in Positive integers}
            sage: L.gens()[2]
            1/(2^z)

        """
        from sage.sets.family import Family
        from sage.sets.positive_integers import PositiveIntegers

        return Family(PositiveIntegers(), lambda n: self.gen(n))

    def _coerce_map_from_(self, S):
        """
        Return ``True`` if a coercion from ``S`` exists.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.has_coerce_map_from(ZZ)
            True
            sage: L.has_coerce_map_from(QQ)
            False
        """
        if self.base_ring().has_coerce_map_from(S):
            return True

        return False

    def _element_constructor_(self, x=None, valuation=None, constant=None, degree=None):
        """
        Construct a Dirichlet series from ``x``.

        INPUT:

        - ``x`` -- a Dirichlet series, a Dirichlet polynomial, a Python function, or a list of elements in the base ring

        - ``constant`` -- either ``None`` (default: ``None``) or pair of an element of the base ring and an integer

        EXAMPLES::


            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L(3)
            3
            sage: L(lambda i: i, constant=1, degree=6)
            1 + 2/2^z + 3/3^z + 4/4^z + 5/5^z + 1/(6^z) + 1/(7^z) + 1/(8^z) + O(1/(9^z))

            sage: X = L(constant=5, degree=3); X
            5/3^z + 5/4^z + 5/5^z + O(1/(6^z))
            sage: X.valuation()
            3
            sage: e = L(moebius); e
            1 - 1/(2^z) - 1/(3^z) - 1/(5^z) + 1/(6^z) - 1/(7^z) + O(1/(8^z))

            sage: L([0], constant=1)
            1/(2^z) + 1/(3^z) + 1/(4^z) + O(1/(5^z))

            sage: L(constant=1)
            1 + 1/(2^z) + 1/(3^z) + O(1/(4^z))

        Alternatively, ``x`` can be a list of elements of the base ring.
        Then these elements are read as coefficients of the terms of
        degrees starting from the ``valuation``. In this case, ``constant``
        may be just an element of the base ring instead of a tuple or can be
        simply omitted if it is zero::

            sage: f = L([1,2,3,4], 4); f
            1/(4^z) + 2/5^z + 3/6^z + 4/7^z
            sage: g = L([1,3,5,7,9], 6, -1); g
            1/(6^z) + 3/7^z + 5/8^z + 7/9^z + 9/10^z - 1/(11^z) - 1/(12^z) - 1/(13^z) + O(1/(14^z))

        TESTS::

            sage: L = LazyDirichletSeriesRing(GF(2), 'z')
            Traceback (most recent call last):
            ...
            ValueError: positive characteristic not allowed for Dirichlet series

        TODO::

            Add a method to make a copy of self._sparse.
        """
        if valuation is None:
            valuation = 1
        assert valuation > 0, "the valuation of a Dirichlet series must be positive"

        if x is None:
            return self.element_class(self, CoefficientStream_uninitialized(self._sparse, valuation))

        BR = self.base_ring()
        if constant is None:
            constant = ZZ.zero()
        elif isinstance(constant, (tuple, list)):
            constant, degree = constant
        constant = BR(constant)

        if x in BR:
            x = BR(x)
            if not x and not constant:
                coeff_stream = CoefficientStream_zero(self._sparse)
                return self.element_class(self, coeff_stream)
            elif not x:
                x = []
            else:
                x = [x]
        if isinstance(x, (tuple, list)):
            coeff_stream = CoefficientStream_exact(x, self._sparse,
                                                   valuation=valuation,
                                                   constant=constant,
                                                   degree=degree)
            return self.element_class(self, coeff_stream)

        if isinstance(x, LazyDirichletSeries):
            if x._coeff_stream._is_sparse is self._sparse:
                return self.element_class(self, x._coeff_stream)
            # TODO: Implement a way to make a self._sparse copy
            raise NotImplementedError("cannot convert between sparse and dense")
        if callable(x):
            if degree is not None:
                if constant is None:
                    constant = ZZ.zero()
                x = [BR(x(i)) for i in range(1, degree)]
                coeff_stream = CoefficientStream_exact(x, self._sparse,
                                                       valuation=valuation,
                                                       constant=constant,
                                                       degree=degree)
                return self.element_class(self, coeff_stream)
            coeff_stream = CoefficientStream_coefficient_function(x, BR, self._sparse, valuation)
            return self.element_class(self, coeff_stream)
        raise ValueError(f"unable to convert {x} into a lazy Dirichlet series")

    def _an_element_(self):
        """
        Return a Dirichlet series in this ring.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.an_element()
            1/(4^z) + 1/(5^z) + 1/(6^z) + ...
        """
        c = self.base_ring().an_element()
        return self.element_class(self, CoefficientStream_exact([], self._sparse, constant=c, valuation=4))

    @cached_method
    def one(self):
        """
        Return the constant series `1`.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.one()
            1
        """
        return self.element_class(self, CoefficientStream_exact([1], self._sparse, valuation=1))

    @cached_method
    def zero(self):
        """
        Return the zero series.

        EXAMPLES::

            sage: L = LazyDirichletSeriesRing(ZZ, 'z')
            sage: L.zero()
            0
        """
        return self.element_class(self, CoefficientStream_zero(self._sparse))

    # add options to class
    class options(GlobalOptions):
        NAME = 'LazyDirichletSeriesRing'
        module = 'sage.rings.lazy_laurent_series_ring'
        display_length = dict(default=7,
                              description='the number of coefficients to display from the valuation',
                              checker=lambda x: x in ZZ and x > 0)
        constant_length = dict(default=3,
                               description='the number of coefficients to display for nonzero constant series',
                               checker=lambda x: x in ZZ and x > 0)
