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
from sage.categories.sets_cat import Sets
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

class CartesianProduct(UniqueRepresentation, Parent):
    """
    A class implementing a raw data structure for cartesian products
    of sets (and elements thereof). See
    :const:`~sage.categories.cartesian_product.cartesian_product` for
    how to construct full fledge cartesian products.
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

    def _repr_(self):
        """
        EXAMPLES::

            sage: cartesian_product([QQ, ZZ, ZZ]) # indirect doctest
            The cartesian product of (Rational Field, Integer Ring, Integer Ring)
        """
        return "The cartesian product of %s"%(self._sets,)

    def summands(self):
        """
        Return the summands of ``self``

        .. see also:: :meth:`Sets.CartesianProducts.ParentMethods.summands()
            <sage.categories.sets_cat.Sets.CartesianProducts.ParentMethods.summands>`.

        EXAMPLES::

            sage: cartesian_product([QQ, ZZ, ZZ]).summands()
            (Rational Field, Integer Ring, Integer Ring)
        """
        return self._sets

    def _sets_keys(self):
        """
        Returns the indices of the summands of ``self``
        as per
        :meth:`Sets.CartesianProducts.ParentMethods._sets_keys()
        <sage.categories.sets_cat.Sets.CartesianProducts.ParentMethods._sets_keys>`.

        EXAMPLES::

            sage: cartesian_product([QQ, ZZ, ZZ])._sets_keys()
            [0, 1, 2]
        """
        return range(len(self._sets))

    @cached_method
    def summand_projection(self, i):
        """
        Returns the natural projection onto the `i`-th summand of self
        as per
        :meth:`Sets.CartesianProducts.ParentMethods.summand_projection()
        <sage.categories.sets_cat.Sets.CartesianProducts.ParentMethods.summand_projection>`.

        INPUTS:

         - ``i`` -- the index of a summand of self

        EXAMPLES::

            sage: C = Sets().CartesianProducts().example(); C
            The cartesian product of (Set of prime numbers (basic implementation), An example of an infinite enumerated set: the non negative integers, An example of a finite enumerated set: {1,2,3})
            sage: x = C.an_element(); x
            (47, 42, 1)
            sage: pi = C.summand_projection(1)
            sage: pi(x)
            42
        """
        assert i in self._sets_keys()
        return attrcall("summand_projection", i)

    def _cartesian_product_of_elements(self, elements):
        """
        Returns the cartesian product of the given ``elements``,
        as per
        :meth:`Sets.CartesianProducts.ParentMethods._cartesian_product_of_elements()`.
        <sage.categories.sets_cat.Sets.CartesianProducts.ParentMethods._cartesian_product_of_elements>`.

        INPUT:

         - ``elements`` - a tuple with one element of each summand of self

        EXAMPLES::

            sage: S1 = Sets().example()
            sage: S2 = InfiniteEnumeratedSets().example()
            sage: C = cartesian_product([S2, S1, S2])
            sage: C._cartesian_product_of_elements([S2.an_element(), S1.an_element(), S2.an_element()])
            (42, 47, 42)
        """
        elements = tuple(elements)
        assert len(elements) == len(self._sets)
        # assert all(zip(self._sets[i], elements, operator.__contains__))
        return self(tuple(elements))

    an_element = Sets.CartesianProducts.ParentMethods.an_element

    from sage.structure.element_wrapper import ElementWrapper
    class Element(ElementWrapper):

        def summand_projection(self, i):
            """
            Returns the projection of ``self`` on the `i`-th
            summand of the cartesian product, as per
            :meth:`Sets.CartesianProducts.ElementMethods.summand_projection()
            <sage.categories.sets_cat.Sets.CartesianProducts.ElementMethods.summand_projection>`.

            INPUTS:

             - ``i`` -- the index of a summand of the cartesian product

            EXAMPLES::

                sage: C = Sets().CartesianProducts().example(); C
                The cartesian product of (Set of prime numbers (basic implementation), An example of an infinite enumerated set: the non negative integers, An example of a finite enumerated set: {1,2,3})
                sage: x = C.an_element(); x
                (47, 42, 1)
                sage: x.summand_projection(1)
                42
            """
            return self.value[i]

        wrapped_class = tuple
