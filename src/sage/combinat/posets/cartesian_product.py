"""
Cartesian products of Posets

AUTHORS:

- Daniel Krenn (2015)

"""
#*****************************************************************************
#  Copyright (C) 2015 Daniel Krenn <dev@danielkrenn.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************

from sage.sets.cartesian_product import CartesianProduct


class CartesianProductPosets(CartesianProduct):
    r"""
    A class implementing cartesian products of posets (and elements
    thereof). Compared to :class:`CartesianProduct` you are able to
    specify an order for comparison of the elements.

    INPUT:

    - ``sets`` -- a tuple of parents.

    - ``category`` -- a subcategory of
      ``Sets().CartesianProducts() & Posets()``.

    - ``order`` -- a string or function specifying an order less or equal.
      It can be one of the following:

      - ``'lex'`` -- elements are ordered lexicographically.

      - ``'product'`` -- an element is less or equal to another
        element, if less or equal is true for all its components
        (cartesian projections).

      - a function ``order_le(left, right)``, which performs the comparison.

    Other keyword arguments (``kwargs``) are passed to the constructor
    of :class:`CartesianProduct`.

    EXAMPLES::

        sage: P = Poset((srange(3), lambda left, right: left <= right))
        sage: Cl = cartesian_product((P, P), order='lex')
        sage: Cl((1, 1)) <= Cl((2, 0))
        True
        sage: Cc = cartesian_product((P, P), order='product')
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

    def __init__(self, sets, category, order=None, **kwargs):
        r"""
        See :class:`CartesianProductPosets` for details.

        TESTS::

            sage: P = Poset((srange(3), lambda left, right: left <= right))
            sage: Cl = cartesian_product((P, P), order='notexisting')
            Traceback (most recent call last):
            ...
            ValueError: No order 'notexisting' known.
        """
        if order is None:
            self._le_ = self.le_product
        elif isinstance(order, str):
            try:
                self._le_ = getattr(self, 'le_' + order)
            except AttributeError:
                raise ValueError("No order '%s' known." % (order,))
        else:
            self._le_ = order

        from sage.categories.posets import Posets
        category = category & Posets()
        super(CartesianProductPosets, self).__init__(
            sets, category, **kwargs)


    def le(self, left, right):
        r"""
        Tests if ``left`` is less than or equal to ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method uses the order defined on creation of this
            cartesian product. See :class:`CartesianProductPosets`.

        EXAMPLES::

            sage: QQ.CartesianProduct = sage.combinat.posets.cartesian_product.CartesianProductPosets
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


    def le_product(self, left, right):
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
            sage: Q = cartesian_product((P, P), order='product')
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

                sage: QQ.CartesianProduct = sage.combinat.posets.cartesian_product.CartesianProductPosets
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

                sage: QQ.CartesianProduct = sage.combinat.posets.cartesian_product.CartesianProductPosets
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

                sage: QQ.CartesianProduct = sage.combinat.posets.cartesian_product.CartesianProductPosets
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

                sage: QQ.CartesianProduct = sage.combinat.posets.cartesian_product.CartesianProductPosets
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

                sage: QQ.CartesianProduct = sage.combinat.posets.cartesian_product.CartesianProductPosets
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
