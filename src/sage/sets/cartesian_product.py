"""
Cartesian products

AUTHORS:

- Nicolas Thiery (2010-03): initial version
- Daniel Krenn (2015-06): cartesian products for posets with different orders

"""
#*****************************************************************************
#       Copyright (C) 2008 Nicolas Thiery <nthiery at users.sf.net>,
#                          Mike Hansen <mhansen@gmail.com>,
#                          Florent Hivert <Florent.Hivert@univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import itertools

from sage.misc.misc import attrcall
from sage.misc.cachefunc import cached_method
from sage.misc.superseded import deprecated_function_alias
from sage.misc.misc_c import prod

from sage.categories.sets_cat import Sets

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper

from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity

class CartesianProduct(UniqueRepresentation, Parent):
    """
    A class implementing a raw data structure for cartesian products
    of sets (and elements thereof). See :obj:`cartesian_product` for
    how to construct full fledged cartesian products.

    EXAMPLES::

        sage: G = cartesian_product([GF(5), Permutations(10)])
        sage: G.cartesian_factors()
        (Finite Field of size 5, Standard permutations of 10)
        sage: G.cardinality()
        18144000
        sage: G.random_element()    # random
        (1, [4, 7, 6, 5, 10, 1, 3, 2, 8, 9])
        sage: G.category()
        Join of Category of finite monoids
            and Category of Cartesian products of monoids
            and Category of Cartesian products of finite enumerated sets

    .. automethod:: _cartesian_product_of_elements
    """
    def __init__(self, sets, category, flatten=False, **kwargs):
        r"""
        INPUT:

         - ``sets`` -- a tuple of parents
         - ``category`` -- a subcategory of ``Sets().CartesianProducts()``
         - ``flatten`` -- a boolean (default: ``False``)

        ``flatten`` is current ignored, and reserved for future use.

        No other keyword arguments (``kwargs``) are accepted.

        TESTS::

            sage: from sage.sets.cartesian_product import CartesianProduct
            sage: C = CartesianProduct((QQ, ZZ, ZZ), category = Sets().CartesianProducts())
            sage: C
            The cartesian product of (Rational Field, Integer Ring, Integer Ring)
            sage: C.an_element()
            (1/2, 1, 1)
            sage: TestSuite(C).run()
            sage: cartesian_product([ZZ, ZZ], blub=None)
            Traceback (most recent call last):
            ...
            TypeError: unknown parameters: blub
        """
        if kwargs:
            raise TypeError('unknown parameters: %s' %
                            ', '.join(str(k) for k in kwargs.iterkeys()))
        self._sets = tuple(sets)
        Parent.__init__(self, category=category)

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
        x = tuple(c(xx) for c,xx in itertools.izip(self._sets,x))
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
            {0, 1, 2}
            sage: cartesian_product([ZZ]*100)._sets_keys()
            {0, ..., 99}
        """
        from sage.sets.integer_range import IntegerRange
        return IntegerRange(len(self._sets))

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

            sage: C.cartesian_projection('hey')
            Traceback (most recent call last):
            ...
            ValueError: i (=hey) must be in {0, 1, 2}
        """
        if i not in self._sets_keys():
            raise ValueError("i (={}) must be in {}".format(i, self._sets_keys()))
        return attrcall("cartesian_projection", i)

    summand_projection = deprecated_function_alias(10963, cartesian_projection)

    def _cartesian_product_of_elements(self, elements):
        """
        Return the cartesian product of the given ``elements``.

        This implements :meth:`Sets.CartesianProducts.ParentMethods._cartesian_product_of_elements`.
        INPUT:

        - ``elements`` -- an iterable (e.g. tuple, list) with one element of
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

    def construction(self):
        r"""
        Return the construction functor and its arguments for this
        cartesian product.

        OUTPUT:

        A pair whose first entry is a cartesian product functor and
        its second a list of the cartesian factors.

        EXAMPLES::

            sage: cartesian_product([ZZ, QQ]).construction()
            (The cartesian_product functorial construction,
             (Integer Ring, Rational Field))
        """
        from sage.categories.cartesian_product import cartesian_product
        return cartesian_product, self.cartesian_factors()

    def _coerce_map_from_(self, S):
        r"""
        Return ``True`` if ``S`` coerces into this cartesian product.

        TESTS::

            sage: Z = cartesian_product([ZZ])
            sage: Q = cartesian_product([QQ])
            sage: Z.has_coerce_map_from(Q)  # indirect doctest
            False
            sage: Q.has_coerce_map_from(Z)  # indirect doctest
            True
        """
        if isinstance(S, CartesianProduct):
            S_factors = S.cartesian_factors()
            R_factors = self.cartesian_factors()
            if len(S_factors) == len(R_factors):
                if all(r.has_coerce_map_from(s) for r,s in zip(R_factors, S_factors)):
                    return True

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


class CartesianProductPosets(CartesianProduct):
    r"""
    A class implementing cartesian products of posets (and elements
    thereof). Compared to :class:`CartesianProduct` you are able to
    specify an order for comparison of the elements.

    INPUT:

    - ``sets`` -- a tuple of parents.

    - ``category`` -- a subcategory of
      ``Sets().CartesianProducts() & Posets()``.

    - ``order`` -- a string or function specifing an order less or equal.
      It can be one of the following:

      - ``'lex'`` -- elements are ordered lexicographically.

      - ``'components'`` -- an element is less or equal to another
        element, if less or equal is true for all its components
        (cartesian projections).

      - a function ``order_le(left, right)``, which performs the comparison.

    Other keyword arguments (``kwargs``) are passed to the constructor
    of :class:`CartesianProduct`.

    EXAMPLES::

        sage: P = Poset((srange(3), lambda left, right: left <= right))
        sage: P.CartesianProduct = sage.sets.cartesian_product.CartesianProductPosets
        sage: Cl = cartesian_product((P, P), order='lex')
        sage: Cl((1, 1)) <= Cl((2, 0))
        True
        sage: Cc = cartesian_product((P, P), order='components')
        sage: Cc((1, 1)) <= Cc((2, 0))
        False
        sage: def le_sum(left, right):
        ....:     return (sum(left) < sum(right) or
        ....:             sum(left) == sum(right) and left[0] <= right[0])
        sage: Cs = cartesian_product((P, P), order=le_sum)
        sage: Cs((1, 1)) <= Cs((2, 0))
        True

    TESTS::

        sage: Cc.category()
        Join of Category of finite posets and
        Category of Cartesian products of finite enumerated sets
        sage: Cs.category()
        Join of Category of finite posets and
        Category of Cartesian products of finite enumerated sets

    .. SEEALSO:

        :class:`CartesianProduct`
    """

    def __init__(self, sets, category, order, **kwargs):
        r"""
        See :class:`CartesianProductPosets` for details,

        TESTS::

            sage: P = Poset((srange(3), lambda left, right: left <= right))
            sage: P.CartesianProduct = sage.sets.cartesian_product.CartesianProductPosets
            sage: Cl = cartesian_product((P, P), order='notexisting')
            Traceback (most recent call last):
            ...
            ValueError: No order 'notexisting' known.
        """
        if isinstance(order, str):
            try:
                self._le_ = getattr(self, 'le_' + order)
            except AttributeError:
                raise ValueError("No order '%s' known." % (order,))
        else:
            self._le_ = order

        super(CartesianProductPosets, self).__init__(
            sets, category, **kwargs)
        from sage.categories.posets import Posets
        self._refine_category_(Posets())


    def le(self, left, right):
        r"""
        Tests if ``left`` is smaller or equal to ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method uses the order defined on creation of this
            cartesian product. See :class:`CartesianProductPosets`

        EXAMPLES::

            sage: QQ.CartesianProduct = sage.sets.cartesian_product.CartesianProductPosets
            sage: def le_sum(left, right):
            ....:     return (sum(left) < sum(right) or
            ....:             sum(left) == sum(right) and left[0] <= right[0])
            sage: C = cartesian_product((QQ, QQ), order=le_sum)
            sage: C.le(C((1/3, 2)), C((2, 1/3)))
            True
            sage: C.le(C((2, 1/3)), C((1/3, 2)))
            False
            sage: C.le(C((1/3, 2)), C((2, 2)))
            True
            sage: C.le(C((2, 2)), C((1/3, 2)))
            False
        """
        return self._le_(left, right)


    def le_lex(self, left, right):
        r"""
        Tests if ``left`` is lexicographically smaller or equal
        to ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: P = Poset((srange(2), lambda left, right: left <= right))
            sage: P.CartesianProduct = sage.sets.cartesian_product.CartesianProductPosets
            sage: Q = cartesian_product((P, P), order='lex')
            sage: T = [Q((0, 0)), Q((1, 1)), Q((0, 1)), Q((1, 0))]
            sage: for a in T:
            ....:     for b in T:
            ....:         assert(Q.le(a, b) == (a <= b))
            ....:         print '%s <= %s = %s' % (a, b, a <= b)
            (0, 0) <= (0, 0) = True
            (0, 0) <= (1, 1) = True
            (0, 0) <= (0, 1) = True
            (0, 0) <= (1, 0) = True
            (1, 1) <= (0, 0) = False
            (1, 1) <= (1, 1) = True
            (1, 1) <= (0, 1) = False
            (1, 1) <= (1, 0) = False
            (0, 1) <= (0, 0) = False
            (0, 1) <= (1, 1) = True
            (0, 1) <= (0, 1) = True
            (0, 1) <= (1, 0) = True
            (1, 0) <= (0, 0) = False
            (1, 0) <= (1, 1) = True
            (1, 0) <= (0, 1) = False
            (1, 0) <= (1, 0) = True
        """
        for l, r, S in \
                zip(left.value, right.value, self.cartesian_factors()):
            if l == r:
                continue
            if S.le(l, r):
                return True
            if S.le(r, l):
                return False
        return True  # equal


    def le_components(self, left, right):
        r"""
        Tests if ``left`` is component-wise smaller or equal
        to ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        A boolean.

        The comparison is ``True`` if the result of the
        comparision in each component is ``True``.

        EXAMPLES::

            sage: P = Poset((srange(2), lambda left, right: left <= right))
            sage: P.CartesianProduct = sage.sets.cartesian_product.CartesianProductPosets
            sage: Q = cartesian_product((P, P), order='components')
            sage: T = [Q((0, 0)), Q((1, 1)), Q((0, 1)), Q((1, 0))]
            sage: for a in T:
            ....:     for b in T:
            ....:         assert(Q.le(a, b) == (a <= b))
            ....:         print '%s <= %s = %s' % (a, b, a <= b)
            (0, 0) <= (0, 0) = True
            (0, 0) <= (1, 1) = True
            (0, 0) <= (0, 1) = True
            (0, 0) <= (1, 0) = True
            (1, 1) <= (0, 0) = False
            (1, 1) <= (1, 1) = True
            (1, 1) <= (0, 1) = False
            (1, 1) <= (1, 0) = False
            (0, 1) <= (0, 0) = False
            (0, 1) <= (1, 1) = True
            (0, 1) <= (0, 1) = True
            (0, 1) <= (1, 0) = False
            (1, 0) <= (0, 0) = False
            (1, 0) <= (1, 1) = True
            (1, 0) <= (0, 1) = False
            (1, 0) <= (1, 0) = True
        """
        return all(
            S.le(l, r)
            for l, r, S in
            zip(left.value, right.value, self.cartesian_factors()))


    class Element(CartesianProduct.Element):

        def _le_(self, other):
            r"""
            Return if this element is less or equal to ``other``.

            INPUT:

            - ``other`` -- an element.

            OUTPUT:

            A boolean.

            .. NOTE::

                This method calls :meth:`CartesianProductPosets.le`. Override
                it in inherited class to change this.

                It can be assumed that this element and ``other`` have
                the same parent.

            TESTS::

                sage: QQ.CartesianProduct = sage.sets.cartesian_product.CartesianProductPosets
                sage: def le_sum(left, right):
                ....:     return (sum(left) < sum(right) or
                ....:             sum(left) == sum(right) and left[0] <= right[0])
                sage: C = cartesian_product((QQ, QQ), order=le_sum)
                sage: C((1/3, 2)) <= C((2, 1/3))  # indirect doctest
                True
                sage: C((1/3, 2)) <= C((2, 2))  # indirect doctest
                True
            """
            return self.parent().le(self, other)


        def __le__(self, other):
            r"""
            Return if this element is less than or equal to ``other``.

            INPUT:

            - ``other`` -- an element.

            OUTPUT:

            A boolean.

            .. NOTE::

                This method uses the coercion framework to find a
                suitable common parent.

                This method can be deleted once :trac:`10130` is fixed and
                provides these methods automatically.

            TESTS::

                sage: QQ.CartesianProduct = sage.sets.cartesian_product.CartesianProductPosets
                sage: def le_sum(left, right):
                ....:     return (sum(left) < sum(right) or
                ....:             sum(left) == sum(right) and left[0] <= right[0])
                sage: C = cartesian_product((QQ, QQ), order=le_sum)
                sage: C((1/3, 2)) <= C((2, 1/3))
                True
                sage: C((1/3, 2)) <= C((2, 2))
                True
            """
            from sage.structure.element import have_same_parent
            if have_same_parent(self, other):
                return self._le_(other)

            from sage.structure.element import get_coercion_model
            import operator
            try:
                return get_coercion_model().bin_op(self, other, operator.le)
            except TypeError:
                return False


        def __ge__(self, other):
            r"""
            Return if this element is greater than or equal to ``other``.

            INPUT:

            - ``other`` -- an element.

            OUTPUT:

            A boolean.

            .. NOTE::

                This method uses the coercion framework to find a
                suitable common parent.

                This method can be deleted once :trac:`10130` is fixed and
                provides these methods automatically.

            TESTS::

                sage: QQ.CartesianProduct = sage.sets.cartesian_product.CartesianProductPosets
                sage: def le_sum(left, right):
                ....:     return (sum(left) < sum(right) or
                ....:             sum(left) == sum(right) and left[0] <= right[0])
                sage: C = cartesian_product((QQ, QQ), order=le_sum)
                sage: C((1/3, 2)) >= C((2, 1/3))
                False
                sage: C((1/3, 2)) >= C((2, 2))
                False
            """
            return other.__le__(self)


        def __lt__(self, other):
            r"""
            Return if this element is less than ``other``.

            INPUT:

            - ``other`` -- an element.

            OUTPUT:

            A boolean.

            .. NOTE::

                This method uses the coercion framework to find a
                suitable common parent.

                This method can be deleted once :trac:`10130` is fixed and
                provides these methods automatically.

            TESTS::

                sage: QQ.CartesianProduct = sage.sets.cartesian_product.CartesianProductPosets
                sage: def le_sum(left, right):
                ....:     return (sum(left) < sum(right) or
                ....:             sum(left) == sum(right) and left[0] <= right[0])
                sage: C = cartesian_product((QQ, QQ), order=le_sum)
                sage: C((1/3, 2)) < C((2, 1/3))
                True
                sage: C((1/3, 2)) < C((2, 2))
                True
            """
            return not self == other and self.__le__(other)


        def __gt__(self, other):
            r"""
            Return if this element is greater than ``other``.

            INPUT:

            - ``other`` -- an element.

            OUTPUT:

            A boolean.

            .. NOTE::

                This method uses the coercion framework to find a
                suitable common parent.

                This method can be deleted once :trac:`10130` is fixed and
                provides these methods automatically.

            TESTS::

                sage: QQ.CartesianProduct = sage.sets.cartesian_product.CartesianProductPosets
                sage: def le_sum(left, right):
                ....:     return (sum(left) < sum(right) or
                ....:             sum(left) == sum(right) and left[0] <= right[0])
                sage: C = cartesian_product((QQ, QQ), order=le_sum)
                sage: C((1/3, 2)) > C((2, 1/3))
                False
                sage: C((1/3, 2)) > C((2, 2))
                False
            """
            return not self == other and other.__le__(self)
