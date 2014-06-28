"""
Cartesian products

AUTHORS:

- Nicolas Thiery (2010-03): initial version

"""
#*****************************************************************************
#       Copyright (C) 2008 Nicolas Thiery <nthiery at users.sf.net>,
#                          Mike Hansen <mhansen@gmail.com>,
#                          Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.misc import attrcall
from sage.misc.cachefunc import cached_method
from sage.misc.superseded import deprecated_function_alias
from sage.categories.sets_cat import Sets
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper

class CartesianProduct(UniqueRepresentation, Parent):
    """
    A class implementing a raw data structure for cartesian products
    of sets (and elements thereof). See :obj:`cartesian_product` for
    how to construct full fledge cartesian products.

    .. automethod:: _cartesian_product_of_elements
    """

    def __init__(self, sets, category, flatten = False):
        r"""
        INPUT:

         - ``sets`` -- a tuple of parents
         - ``category`` -- a subcategory of ``Sets().CartesianProducts()``
         - ``flatten`` -- a boolean (default: ``False``)

        ``flatten`` is current ignored, and reserved for future use.

        TESTS::

            sage: from sage.sets.cartesian_product import CartesianProduct
            sage: C = CartesianProduct((QQ, ZZ, ZZ), category = Sets().CartesianProducts())
            sage: C
            The cartesian product of (Rational Field, Integer Ring, Integer Ring)
            sage: C.an_element()
            (1/2, 1, 1)
            sage: TestSuite(C).run()
        """
        self._sets = sets
        Parent.__init__(self, category = category)

    def _element_constructor_(self,x):
        r"""
        Construct an element of a cartesian product from a list or iterable

        INPUT:

        - ``x`` -- a list (or iterable)

        Each component of `x` is converted to the corresponding
        cartesian factor.

        EXAMPLES::

            sage: C = cartesian_product([GF(5), GF(3)])
            sage: x = C((1,3)); x
            (1, 0)
            sage: x.parent()
            The cartesian product of (Finite Field of size 5, Finite Field of size 3)
            sage: x[0].parent()
            Finite Field of size 5
            sage: x[1].parent()
            Finite Field of size 3

        An iterable is also accepted as input::

            sage: C(i for i in range(2))
            (0, 1)

        TESTS::

            sage: C((1,3,4))
            Traceback (most recent call last):
            ...
            ValueError: (1, 3, 4) should be of length 2
        """
        x = tuple(x)
        if len(x) != len(self._sets):
            raise ValueError(
                "{} should be of length {}".format(x, len(self._sets)))
        x = tuple(c(xx) for c,xx in zip(self._sets,x))
        return self.element_class(self, x)

    def _repr_(self):
        """
        EXAMPLES::

            sage: cartesian_product([QQ, ZZ, ZZ]) # indirect doctest
            The cartesian product of (Rational Field, Integer Ring, Integer Ring)
        """
        return "The cartesian product of %s"%(self._sets,)

    def cartesian_factors(self):
        """
        Return the cartesian factors of ``self``.

        .. SEEALSO::

            :meth:`Sets.CartesianProducts.ParentMethods.cartesian_factors()
            <sage.categories.sets_cat.Sets.CartesianProducts.ParentMethods.cartesian_factors>`.

        EXAMPLES::

            sage: cartesian_product([QQ, ZZ, ZZ]).cartesian_factors()
            (Rational Field, Integer Ring, Integer Ring)
        """
        return self._sets

    def _sets_keys(self):
        """
        Return the indices of the cartesian factors of ``self``
        as per
        :meth:`Sets.CartesianProducts.ParentMethods._sets_keys()
        <sage.categories.sets_cat.Sets.CartesianProducts.ParentMethods._sets_keys>`.

        EXAMPLES::

            sage: cartesian_product([QQ, ZZ, ZZ])._sets_keys()
            [0, 1, 2]
        """
        return range(len(self._sets))

    @cached_method
    def cartesian_projection(self, i):
        """
        Return the natural projection onto the `i`-th cartesian
        factor of ``self`` as per
        :meth:`Sets.CartesianProducts.ParentMethods.cartesian_projection()
        <sage.categories.sets_cat.Sets.CartesianProducts.ParentMethods.cartesian_projection>`.

        INPUT:

        - ``i`` -- the index of a cartesian factor of ``self``

        EXAMPLES::

            sage: C = Sets().CartesianProducts().example(); C
            The cartesian product of (Set of prime numbers (basic implementation), An example of an infinite enumerated set: the non negative integers, An example of a finite enumerated set: {1,2,3})
            sage: x = C.an_element(); x
            (47, 42, 1)
            sage: pi = C.cartesian_projection(1)
            sage: pi(x)
            42
        """
        assert i in self._sets_keys()
        return attrcall("cartesian_projection", i)

    summand_projection = deprecated_function_alias(10963, cartesian_projection)

    def __iter__(self):
        r"""
        Iterates over the elements of self.

        EXAMPLE::

            sage: F33 = GF(2).cartesian_product(GF(2))
            sage: list(F33)
            [(0, 0), (0, 1), (1, 0), (1, 1)]
        """
        from itertools import product
        for x in product(*self._sets):
            yield self(x)

    def _cartesian_product_of_elements(self, elements):
        """
        Return the cartesian product of the given ``elements``.

        This implements :meth:`Sets.CartesianProducts.ParentMethods._cartesian_product_of_elements`.

        INPUT:

        - ``elements`` -- a tuple (or iterable) with one element of
          each cartesian factor of ``self``

        .. WARNING::

            This is meant as a fast low-level method. In particular,
            no coercion is attempted. When coercion or sanity checks
            are desirable, please use instead ``self(elements)`` or
            ``self._element_constructor(elements)``.

        EXAMPLES::

            sage: S1 = Sets().example()
            sage: S2 = InfiniteEnumeratedSets().example()
            sage: C = cartesian_product([S2, S1, S2])
            sage: C._cartesian_product_of_elements([S2.an_element(), S1.an_element(), S2.an_element()])
            (42, 47, 42)
        """
        elements = tuple(elements)
        assert len(elements) == len(self._sets)
        return self.element_class(self, elements)

    an_element = Sets.CartesianProducts.ParentMethods.an_element

    class Element(ElementWrapper):

        wrapped_class = tuple

        def cartesian_projection(self, i):
            r"""
            Return the projection of ``self`` on the `i`-th factor of
            the cartesian product, as per
            :meth:`Sets.CartesianProducts.ElementMethods.cartesian_projection()
            <sage.categories.sets_cat.Sets.CartesianProducts.ElementMethods.cartesian_projection>`.

            INPUTS:

            - ``i`` -- the index of a factor of the cartesian product

            EXAMPLES::

                sage: C = Sets().CartesianProducts().example(); C
                The cartesian product of (Set of prime numbers (basic implementation), An example of an infinite enumerated set: the non negative integers, An example of a finite enumerated set: {1,2,3})
                sage: x = C.an_element(); x
                (47, 42, 1)
                sage: x.cartesian_projection(1)
                42

                sage: x.summand_projection(1)
                doctest:...: DeprecationWarning: summand_projection is deprecated. Please use cartesian_projection instead.
                See http://trac.sagemath.org/10963 for details.
                42
            """
            return self.value[i]

        __getitem__ = cartesian_projection

        def __iter__(self):
            r"""
            Iterate over the components of an element.

            EXAMPLES::

                sage: C = Sets().CartesianProducts().example(); C
                The cartesian product of
                (Set of prime numbers (basic implementation),
                 An example of an infinite enumerated set: the non negative integers,
                 An example of a finite enumerated set: {1,2,3})
                sage: c = C.an_element(); c
                (47, 42, 1)
                sage: for i in c:
                ....:     print i
                47
                42
                1
            """
            return iter(self.value)
