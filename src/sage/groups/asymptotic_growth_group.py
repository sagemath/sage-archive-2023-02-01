r"""
(Asymptotic) Growth Groups

This module adds support for (asymptotic) growth groups. Such groups
are equipped with a partial order: the elements can be seen as
functions and their behavior as the argument(s) get large (tend to
`\infty`) is compared.

Besides an abstract base class :class:`GenericGrowthGroup`, this module
contains concrete realizations of growth groups. At the moment there
is

- :class:`MonomialGrowthGroup` (whose elements are powers of a fixed symbol).

More complex growth groups can be constructed via cartesian products.

These growth groups are used behind the scenes when performing
calculations in an asymptotic ring (to be implemented).

AUTHORS:

- Benjamin Hackl (2015-01): initial version
- Daniel Krenn (2015-05-29): initial version and review
- Daniel Krenn (2015-06-02): cartesian products
- Benjamin Hackl (2015-07): growth group factory

.. WARNING::

    As this code is experimental, warnings are thrown when a growth
    group is created for the first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: G = agg.GenericGrowthGroup(ZZ); G
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        Growth Group Generic(ZZ)
        sage: G = agg.MonomialGrowthGroup(ZZ, 'x'); G
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        Growth Group x^ZZ


.. NOTE::

    By using the following short notation for growth groups, their
    creation is very simple: *Monomial growth groups* (i.e. the
    group for powers of a fixed symbol;
    :class:`~sage.groups.asymptotic_growth_group.MonomialGrowthGroup`)
    are denoted as ``variable^base``, e.g. ``x^ZZ`` and ``y^QQ`` for
    the group of integer powers of `x`, and the group of rational
    powers of `y`, respectively.

    This also enables us to construct *logarithmic growth groups*,
    e.g. ``log(x)^ZZ``.

    This notation will also be extended to *Exponential growth
    groups*.

EXAMPLES::

    sage: import sage.groups.asymptotic_growth_group as agg
    sage: G_x = agg.GrowthGroup('x^ZZ'); repr(G_x)
    'Growth Group x^ZZ'
    sage: G_xy = agg.GrowthGroup('x^ZZ * y^ZZ'); G_xy
    Growth Group x^ZZ * y^ZZ
    sage: G_xy.an_element()
    x * y
    sage: x = G_xy('x'); y = G_xy('y')
    sage: elem = x^21 * y^21; elem^2
    x^42 * y^42

A monomial growth group itself is totally ordered, all elements
are comparable. However, this does **not** hold for cartesian
products::

    sage: e1 = x^2 * y; e2 = x * y^2
    sage: e1 <= e2 or e2 <= e1
    False

In terms of uniqueness, the order of factors in the cartesian
product matters::

    sage: agg.GrowthGroup('x^ZZ * y^ZZ') is agg.GrowthGroup('y^ZZ * x^ZZ')
    False
"""

#*****************************************************************************
# Copyright (C) 2014--2015 Benjamin Hackl <benjamin.hackl@aau.at>
#               2014--2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage

def repr_short_to_parent(s):
    r"""
    Helper method for the growth group factory, which converts a short
    representation string to a parent.

    INPUT:

    A string.

    OUTPUT:

    A parent.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: agg.repr_short_to_parent('ZZ')
        Integer Ring
        sage: agg.repr_short_to_parent('QQ')
        Rational Field
        sage: agg.repr_short_to_parent('SR')
        Symbolic Ring

    TESTS::

        sage: agg.repr_short_to_parent('abcdef')
        Traceback (most recent call last):
        ...
        ValueError: Cannot create a parent out of 'abcdef'.
    """
    if s == 'ZZ':
        return sage.rings.integer_ring.ZZ
    elif s == 'QQ':
        return sage.rings.rational_field.QQ
    elif s == 'SR':
        return sage.symbolic.ring.SR
    else:
        raise ValueError("Cannot create a parent out of '%s'." % s)


def parent_to_repr_short(P):
    r"""
    Helper method, which generates a short(er) representation string
    out of a parent.

    INPUT:

    A parent.

    OUTPUT:

    A string.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: agg.parent_to_repr_short(ZZ)
        'ZZ'
        sage: agg.parent_to_repr_short(QQ)
        'QQ'
        sage: agg.parent_to_repr_short(SR)
        'SR'
        sage: agg.parent_to_repr_short(ZZ[x])
        '(Univariate Polynomial Ring in x over Integer Ring)'
    """
    if P is sage.rings.integer_ring.ZZ:
        return 'ZZ'
    elif P is sage.rings.rational_field.QQ:
        return 'QQ'
    elif P is sage.symbolic.ring.SR:
        return 'SR'
    else:
        rep = repr(P)
        if ' ' in rep:
            rep = '(' + rep + ')'
        return rep


from sage.sets.cartesian_product import CartesianProductPosets
class CartesianProductGrowthGroups(CartesianProductPosets):
    r"""
    A cartesian product of growth groups.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.MonomialGrowthGroup(QQ, 'x')
        sage: L = agg.MonomialGrowthGroup(ZZ, 'log(x)')
        sage: C = cartesian_product([P, L], order='lex'); C
        Growth Group x^QQ * log(x)^ZZ
        sage: C.an_element()
        x^(1/2) * log(x)

    ::

        sage: Px = agg.MonomialGrowthGroup(QQ, 'x')
        sage: Lx = agg.MonomialGrowthGroup(ZZ, 'log(x)')
        sage: Cx = cartesian_product([Px, Lx], order='lex')
        sage: Py = agg.MonomialGrowthGroup(QQ, 'y')
        sage: C = cartesian_product([Cx, Py], order='components'); C
        Growth Group x^QQ * log(x)^ZZ * y^QQ
        sage: C.an_element()
        x^(1/2) * log(x) * y^(1/2)

    .. SEEALSO:

        :class:`~sage.sets.cartesian_product.CartesianProduct`,
        :class:`~sage.sets.cartesian_product.CartesianProductPosets`.
    """

    def _element_constructor_(self, data):
        r"""
        Converts the given object to an element of this cartesian
        product.

        EXAMPLES::

            sage: from sage.groups.asymptotic_growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ * y^ZZ')
            sage: G_log = GrowthGroup('x^ZZ * log(x)^ZZ * y^ZZ')

        Conversion from the symbolic ring works::

            sage: x,y = var('x y')
            sage: G(x^-3 * y^2)
            x^(-3) * y^2
            sage: G(x^4), G(y^2)
            (x^4, y^2)
            sage: G(1)
            1

        Even more complex expressions can be parsed::

            sage: G_log(x^42 * log(x)^-42 * y^42)
            x^42 * log(x)^(-42) * y^42

        TESTS::

            sage: G = GrowthGroup('x^ZZ * y^ZZ')
            sage: G('x'), G('y')
            (x, y)
        """
        if data == 1:
            return self.one()

        if isinstance(data, list):
            try:
                obj = super(CartesianProductGrowthGroups,
                            self)._element_constructor_(data)
                return obj
            except ValueError:
                factors = self.cartesian_factors()
                one = self.one().value
                for k in range(len(data)):
                    for l in range(len(factors)):
                        try:
                            conv_data = factors[l](data[k])
                            data[k] = self([one[j] if j != l else conv_data
                                            for j in range(len(one))])
                            break
                        except (ValueError, TypeError):
                            continue
                return self.prod(data)


        if hasattr(data, 'parent'):
            if data.parent() is self:
                return data

            elif data.parent() is sage.symbolic.ring.SR:
                import operator
                from sage.symbolic.operators import mul_vararg
                if data.operator() == operator.pow or data.is_symbol():
                    return self([data])
                elif data.operator() == mul_vararg:
                    return self(data.operands())
            # room for other parents (e.g. polynomial ring et al.)

        # final attempt: try to parse the representation string
        else:
            str_lst = str(data).replace(' ', '').split('*')
            return self(str_lst)


    def _repr_(self):
        r"""
        A representation string for this cartesian product of growth groups.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(QQ, 'x')
            sage: L = agg.MonomialGrowthGroup(ZZ, 'log(x)')
            sage: cartesian_product([P, L], order='lex')._repr_()
            'Growth Group x^QQ * log(x)^ZZ'
        """
        return 'Growth Group ' + self._repr_short_()


    def _repr_short_(self):
        r"""
        A short (shorter than :meth:`._repr_`) representation string
        for this cartesian product of growth groups.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(QQ, 'x')
            sage: L = agg.MonomialGrowthGroup(ZZ, 'log(x)')
            sage: cartesian_product([P, L], order='lex')._repr_short_()
            'x^QQ * log(x)^ZZ'
        """
        return ' * '.join(S._repr_short_() for S in self.cartesian_factors())


    class Element(CartesianProductPosets.Element):
        def _repr_(self):
            r"""
            A representation string for this cartesian product element.

            INPUT:

            Nothing.

            OUTPUT:

            A string.

            EXAMPLES::

                sage: import sage.groups.asymptotic_growth_group as agg
                sage: P = agg.MonomialGrowthGroup(QQ, 'x')
                sage: L = agg.MonomialGrowthGroup(ZZ, 'log(x)')
                sage: cartesian_product([P, L], order='lex').an_element()._repr_()
                'x^(1/2) * log(x)'
            """
            s = ' * '.join(repr(v) for v in self.value if not v.is_one())
            if s == '':
                return '1'
            return s


CartesianProductGrowthGroups.CartesianProduct = CartesianProductGrowthGroups


class GenericGrowthElement(sage.structure.element.MultiplicativeGroupElement):
    r"""
    An abstract implementation of a generic growth element.

    Growth elements form a group by multiplication and (some of) the
    elements can be compared to each other, i.e., all elements form a
    poset.

    INPUT:

    - ``parent`` -- a :class:`GenericGrowthGroup`.

    - ``raw_element`` -- an element from the base of the parent.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: G = agg.GenericGrowthGroup(ZZ)
        sage: g = agg.GenericGrowthElement(G, 42); g
        GenericGrowthElement(42)
        sage: g.parent()
        Growth Group Generic(ZZ)
        sage: G(raw_element=42) == g
        True
    """

    def __init__(self, parent, raw_element):
        r"""
        See :class:`GenericGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G(raw_element=42)
            GenericGrowthElement(42)

        TESTS::

            sage: G(raw_element=42).category()
            Category of elements of Growth Group Generic(ZZ)

        ::

            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G(raw_element=42).category()
            Category of elements of Growth Group Generic(ZZ)

        ::

            sage: agg.GenericGrowthElement(None, 0)
            Traceback (most recent call last):
            ...
            ValueError: The parent must be provided
        """
        if parent is None:
            raise ValueError('The parent must be provided')
        super(GenericGrowthElement, self).__init__(parent=parent)

        self._raw_element_ = parent.base()(raw_element)


    def _repr_(self):
        r"""
        A representation string for this abstract generic element.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G(raw_element=42)  # indirect doctest
            GenericGrowthElement(42)
        """
        return 'GenericGrowthElement(%s)' % (self._raw_element_,)


    def __hash__(self):
        r"""
        Return the hash of this element.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ);
            sage: hash(G(raw_element=42))  # random
            5656565656565656
        """
        return hash((self.parent(), self._raw_element_))


    def _mul_(self, other):
        r"""
        Abstract multiplication method for generic elements.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`.

        OUTPUT:

        A class:`GenericGrowthElement` representing the product with
        ``other``.

        .. NOTE::

            Inherited classes must override this.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: g = G.an_element()
            sage: g * g
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
        """
        raise NotImplementedError('Only implemented in concrete realizations.')


    def _div_(self, other):
        r"""
        Divide this growth element by another one.

        INPUT:

        - ``other`` -- an instance of :class:`GenericGrowthElement`.

        OUTPUT:

        An instance of :class:`GenericGrowthElement`.

        .. NOTE::

            This method is called by the coercion framework, thus, it can be
            assumed that this element, as well as ``other`` are of the same
            type. The output will have this type.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: e1 = P(raw_element=2)
            sage: e2 = e1._div_(P.gen()); e2
            x
            sage: e2 == e1 / P.gen()
            True
        """
        return self._mul_(~other)


    def __pow__(self, power):
        r"""
        Takes this growth element to the given ``power``.

        INPUT:

        - ``power`` -- a number. This can anything that is valid to be
          on the right hand side of ``*`` with an elements of the
          parent's base.

        OUTPUT:

        The result of this exponentiation a :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G.an_element()^7
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
        """
        raise NotImplementedError('Only implemented in concrete realizations.')


    def __eq__(self, other):
        r"""
        Return if this growth element is equal to ``other``.

        INPUT:

        - ``other`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

            The comparison of two elements with the same parent is done in
            :meth:`_eq_`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G.an_element() == G.an_element()
            True
            sage: G(raw_element=42) == G(raw_element=7)
            False

        ::

            sage: G_ZZ = agg.GenericGrowthGroup(ZZ)
            sage: G_QQ = agg.GenericGrowthGroup(QQ)
            sage: G_ZZ(raw_element=1) == G_QQ(raw_element=1)
            True

        ::

            sage: P_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P_ZZ.gen() == P_QQ.gen()
            True
            sage: ~P_ZZ.gen() == P_ZZ.gen()
            False
            sage: ~P_ZZ(1) == P_ZZ(1)
            True
        """
        from sage.structure.element import have_same_parent
        if have_same_parent(self, other):
            return self._eq_(other)

        from sage.structure.element import get_coercion_model
        import operator
        try:
            return get_coercion_model().bin_op(self, other, operator.eq)
        except TypeError:
            return False


    def _eq_(self, other):
        r"""
        Return if this :class:`GenericGrowthElement` is equal to ``other``.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function compares two instances of
            :class:`GenericGrowthElement`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: e1 = P(raw_element=1)
            sage: e1._eq_(P.gen())
            True
            sage: e2 = e1^4
            sage: e2 == e1^2 * e1 * e1
            True
            sage: e2 == e1
            False
        """
        return self._raw_element_ == other._raw_element_


    def __le__(self, other):
        r"""
        Return if this growth element is at most (less than or equal
        to) ``other``.

        INPUT:

        - ``other`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

            The comparison of two elements with the same parent is done in
            :meth:`_le_`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P_ZZ.gen() <= P_QQ.gen()^2
            True
            sage: ~P_ZZ.gen() <= P_ZZ.gen()
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


    def _le_(self, other):
        r"""
        Return if this :class:`GenericGrowthElement` is at most (less
        than or equal to) ``other``.

        INPUT:

        - ``other`` -- a :class:`GenericGrowthElement`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function compares two instances of
            :class:`GenericGrowthElement`.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: e1 = G(raw_element=1); e2 = G(raw_element=2)
            sage: e1 <= e2  # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented in concrete realizations.
        """
        raise NotImplementedError('Only implemented in concrete realizations.')


class GenericGrowthGroup(
        sage.structure.parent.Parent,
        sage.structure.unique_representation.UniqueRepresentation):
    r"""
    An abstract implementation for growth groups.

    INPUT:

    - ``base`` -- one of SageMath's parents, out of which the elements
      get their data (``raw_element``).

    - ``category`` -- (default: ``None``) the category of the newly
      created growth group. It has to be a subcategory of ``Join of
      Category of groups and Category of posets``. This is also the
      default category if ``None`` is specified.

    .. NOTE::

        This class should be derived to get concrete implementations.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: G = agg.GenericGrowthGroup(ZZ); G
        Growth Group Generic(ZZ)

    .. SEEALSO::

        :class:`MonomialGrowthGroup`
    """
    # TODO: implement some sort of 'assume', where basic assumptions
    # for the variables can be stored. --> within the cartesian product

    # enable the category framework for elements
    Element = GenericGrowthElement


    @sage.misc.superseded.experimental(trac_number=17601)
    def __init__(self, base, category=None):
        r"""
        See :class:`GenericGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GenericGrowthGroup(ZZ).category()
            Join of Category of groups and Category of posets

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G.is_parent_of(G(raw_element=42))
            True
            sage: G2 = agg.GenericGrowthGroup(ZZ, category=FiniteGroups() & Posets())
            sage: G2.category()
            Join of Category of finite groups and Category of finite posets
            sage: G3 = agg.GenericGrowthGroup(ZZ, category=Rings())
            Traceback (most recent call last):
            ...
            ValueError: (Category of rings,) is not a subcategory of Join of Category of groups and Category of posets

        ::

            sage: G = agg.GenericGrowthGroup('42')
            Traceback (most recent call last):
            ...
            TypeError: 42 is not a valid base
        """
        if not isinstance(base, sage.structure.parent.Parent):
            raise TypeError('%s is not a valid base' % (base,))
        from sage.categories.groups import Groups
        from sage.categories.posets import Posets

        if category is None:
            category = Groups() & Posets()
        else:
            if not isinstance(category, tuple):
                category = (category,)
            if not any(cat.is_subcategory(Groups() & Posets()) for cat in
                       category):
                raise ValueError('%s is not a subcategory of %s'
                                 % (category, Groups() & Posets()))
        super(GenericGrowthGroup, self).__init__(category=category,
                                                 base=base)


    def _repr_short_(self):
        r"""
        A short representation string of this abstract growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GenericGrowthGroup(QQ)._repr_short_()
            'Generic(QQ)'
            sage: agg.GenericGrowthGroup(QQ)
            Growth Group Generic(QQ)
        """
        return 'Generic(%s)' % (parent_to_repr_short(self.base()),)


    def _repr_(self, condense=False):
        r"""
        A representations string of this growth group.

        INPUT:

        - ``condense`` -- (default: ``False``) if set, then a shorter
          output is returned, e.g. the prefix-string ``Growth Group``
          is not show in this case.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'x')  # indirect doctest
            Growth Group x^ZZ
            sage: agg.MonomialGrowthGroup(QQ, 'log(x)')  # indirect doctest
            Growth Group log(x)^QQ

        TESTS::

            sage: agg.MonomialGrowthGroup(QQ, 'log(x)')._repr_(condense=True)
            'log(x)^QQ'
        """
        pre = 'Growth Group ' if not condense else ''
        return '%s%s' % (pre, self._repr_short_())


    def __hash__(self):
        r"""
        Return the hash of this group.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: hash(agg.GenericGrowthGroup(ZZ))  # random
            4242424242424242
        """
        return hash((self.__class__, self.base()))


    def _an_element_(self):
        r"""
        Return an element of ``self``.

        INPUT:

        Nothing.

        OUTPUT:

        An element of ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ);
            sage: G.an_element()  # indirect doctest
            GenericGrowthElement(1)
        """
        return self.element_class(self, self.base().an_element())


    def le(self, left, right):
        r"""
        Return if the growth of ``left`` is at most (less than or
        equal to) the growth of ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: x = G.gen()
            sage: G.le(x, x^2)
            True
            sage: G.le(x^2, x)
            False
            sage: G.le(x^0, 1)
            True
        """
        return self(left) <= self(right)


    def one(self):
        r"""
        Return the neutral element of this growth group.

        INPUT:

        Nothing.

        OUTPUT:

        An element of this group.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: e1 = agg.MonomialGrowthGroup(ZZ, 'x').one(); e1
            1
            sage: e1.is_idempotent()
            True
        """
        return self(1)


    def _element_constructor_(self, data, raw_element=None):
        r"""
        Converts given object to this growth group.

        INPUT:

        - ``data`` -- an object representing the element to be
          initialized.

        - ``raw_element`` -- (default: ``None``) if given, then this is
          directly passed to the element constructor (i.e., no conversion
          is performed).

        OUTPUT:

        An element of this growth group.

        .. NOTE::

            This method calls :meth:`_convert_`, which does the actual
            conversion from ``data``.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G_ZZ = agg.GenericGrowthGroup(ZZ)
            sage: z = G_ZZ(raw_element=42); z
            GenericGrowthElement(42)
            sage: z is G_ZZ(z)
            True

        ::

            sage: G_QQ = agg.GenericGrowthGroup(QQ)
            sage: q = G_QQ(raw_element=42)
            sage: q is z
            False
            sage: G_ZZ(q)
            GenericGrowthElement(42)
            sage: G_QQ(z)
            GenericGrowthElement(42)
            sage: q is G_ZZ(q)
            False

        ::

            sage: G_ZZ()
            Traceback (most recent call last):
            ...
            ValueError: No input specified. Cannot continue.
            sage: G_ZZ('blub')
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert blub.
            sage: G_ZZ('x', raw_element=42)
            Traceback (most recent call last):
            ...
            ValueError: Input is ambigous: x as well as raw_element=42 are specified
        """
        if raw_element is None:
            if type(data) == self.element_class and data.parent() == self:
                return data
            elif isinstance(data, GenericGrowthElement):
                raw_element = data._raw_element_
            elif type(data) == int and data == 0:
                raise ValueError('No input specified. Cannot continue.')
            else:
                raw_element = self._convert_(data)
            if raw_element is None:
                raise ValueError('Cannot convert %s.' % (data,))
        elif type(data) != int or data != 0:
            raise ValueError('Input is ambigous: '
                             '%s as well as raw_element=%s '
                             'are specified' % (data, raw_element))

        return self.element_class(self, raw_element)


    def _convert_(self, data):
        r"""
        Converts given ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        .. NOTE::

            This method always returns ``None`` in this abstract base
            class and should be overridden in inherited class.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G._convert_('icecream') is None
            True
        """
        pass


    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this growth group.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            Another growth group ``S`` coerces into this growth group
            if and only if the base of ``S`` coerces into the base of
            this growth group.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: G_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: bool(G_ZZ.has_coerce_map_from(G_QQ))  # indirect doctest
            False
            sage: bool(G_QQ.has_coerce_map_from(G_ZZ))  # indirect doctest
            True
        """
        if isinstance(S, GenericGrowthGroup):
            if self.base().has_coerce_map_from(S.base()):
                return True


    def gens_monomial(self):
        r"""
        Return a generator of this growth group, in case one exists.

        INPUT:

        Nothing.

        OUTPUT:

        An element of this growth group or ``None``.

        .. NOTE::

            This method is only implemented for concrete growth
            group implementations.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GenericGrowthGroup(ZZ).gens_monomial()
            Traceback (most recent call last):
            ...
            NotImplementedError: Only implemented for concrete growth group
            implementations.
        """
        raise NotImplementedError("Only implemented for concrete growth group"
                                  " implementations.")


    CartesianProduct = CartesianProductGrowthGroups


class MonomialGrowthElement(GenericGrowthElement):
    r"""
    An implementation of monomial growth elements.

    INPUT:

    - ``parent`` -- a :class:`GenericGrowthGroup`.

    - ``raw_element`` -- an element from the base ring of the parent.

      This ``raw_element`` is the exponent of the created monomial
      growth element.

    A monomial growth element represents a term of the type
    `\operatorname{variable}^{\operatorname{exponent}}`. The multiplication
    corresponds to the addition of the exponents.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
        sage: e1 = P(1); e1
        1
        sage: e2 = P(raw_element=2); e2
        x^2
        sage: e1 == e2
        False
        sage: P.le(e1, e2)
        True
        sage: P.le(e1, P.gen()) and P.le(P.gen(), e2)
        True
    """

    @property
    def exponent(self):
        r"""
        The exponent of this growth element.

        EXAMPLES:

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P(x^42).exponent
            42
        """
        return self._raw_element_


    def _repr_(self):
        r"""
        A representation string for this monomial growth element.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P(1)._repr_()
            '1'
            sage: P(x^5)  # indirect doctest
            x^5
            sage: P(x^(1/2))  # indirect doctest
            x^(1/2)

        TESTS::

            sage: P(x^-1)  # indirect doctest
            1/x
            sage: P(x^-42)  # indirect doctest
            x^(-42)
        """
        from sage.rings.integer_ring import ZZ

        if self.exponent == 0:
            return '1'
        elif self.exponent == 1:
            return self.parent()._var_
        elif self.exponent == -1:
            return '1/' + self.parent()._var_
        elif self.exponent in ZZ and self.exponent > 0:
            return self.parent()._var_ + '^' + str(self.exponent)
        else:
            return self.parent()._var_ + '^(' + str(self.exponent) + ')'


    def _mul_(self, other):
        r"""
        Multiply this monomial growth element with another.

        INPUT:

        - ``other`` -- a :class:`MonomialGrowthElement`

        OUTPUT:

        The product as a :class:`MonomialGrowthElement`.

        .. NOTE::

            Two monomial growth elements are multiplied by adding
            their exponents.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: a = P(x^2)
            sage: b = P(x^3)
            sage: c = a._mul_(b); c
            x^5
            sage: c == a * b
            True
            sage: a * b * a  # indirect doctest
            x^7
        """
        return self.parent()(raw_element=self.exponent + other.exponent)


    def __invert__(self):
        r"""
        Return the multiplicative inverse of this monomial growth element.

        INPUT:

        Nothing.

        OUTPUT:

        The multiplicative inverse as a :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: e1 = P(raw_element=2)
            sage: e2 = e1.__invert__(); e2
            x^(-2)
            sage: e2 == ~e1
            True
        """
        return self.parent()(raw_element=-self.exponent)


    def __pow__(self, power):
        r"""
        Takes this growth element to the given ``power``.

        INPUT:

        - ``power`` -- a number. This can anything that is valid to be
          on the right hand side of ``*`` with an elements of the
          parent's base.

        OUTPUT:

        The result of this exponentiation a :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: x = P.gen()
            sage: a = x^7; a
            x^7
            sage: b = a^(1/2); b
            x^(7/2)
            sage: b.parent()
            Growth Group x^QQ
            sage: b^12
            x^42
        """
        new_exponent = self.exponent * power
        try:
            return self.parent()(raw_element=new_exponent)
        except (ValueError, TypeError):
            pass

        from sage.rings.real_mpfr import RR
        if new_exponent not in RR:
            raise NotImplementedError('Only real exponents are implemented.')

        new_parent = MonomialGrowthGroup(new_exponent.parent(),
                                         self.parent()._var_)
        return new_parent(raw_element=new_exponent)


    def _le_(self, other):
        r"""
        Return if this :class:`MonomialGrowthElement` is at most
        (less than or equal to) ``other``.

        INPUT:

        - ``other`` -- a :class:`MonomialGrowthElement`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function compares two instances of
            :class:`MonomialGrowthElement`.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P_ZZ.gen() <= P_QQ.gen()^2  # indirect doctest
            True
        """
        return self.exponent <= other.exponent


class MonomialGrowthGroup(GenericGrowthGroup):
    r"""
    A growth group dealing with powers of a fixed object/symbol.

    The elements :class:`MonomialGrowthElement` of this group represent powers
    of a fixed base; the group law is the multiplication, which corresponds
    to the addition of the exponents of the monomials.

    INPUT:

    - ``base`` -- one of SageMath's parents, out of which the elements
      get their data (``raw_element``).

      As monomials are represented by this group, the elements in
      ``base`` are the exponents of these monomials.

    - ``var`` -- an object.

      The string representation of ``var`` acts as a base of the
      monomials represented by this group.

    - ``category`` -- (default: ``None``) the category of the newly
      created growth group. It has to be a subcategory of ``Join of
      Category of groups and Category of posets``. This is also the
      default category if ``None`` is specified.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.MonomialGrowthGroup(ZZ, 'x'); P
        Growth Group x^ZZ
        sage: agg.MonomialGrowthGroup(ZZ, log(SR.var('y')))
        Growth Group log(y)^ZZ

    .. SEEALSO::

        :class:`GenericGrowthGroup`
    """

    # enable the category framework for elements
    Element = MonomialGrowthElement


    @staticmethod
    def __classcall__(cls, base, var, category=None):
        r"""
        Normalizes the input in order to ensure a unique
        representation.

        For more information see :class:`MonomialGrowthGroup`.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P1 = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P2 = agg.MonomialGrowthGroup(ZZ, ZZ['x'].gen())
            sage: P3 = agg.MonomialGrowthGroup(ZZ, SR.var('x'))
            sage: P1 is P2 and P2 is P3
            True
            sage: P4 = agg.MonomialGrowthGroup(ZZ, buffer('xylophone', 0, 1))
            sage: P1 is P4
            True
            sage: P5 = agg.MonomialGrowthGroup(ZZ, 'x ')
            sage: P1 is P5
            True

        ::

            sage: L1 = agg.MonomialGrowthGroup(QQ, log(x))
            sage: L2 = agg.MonomialGrowthGroup(QQ, 'log(x)')
            sage: L1 is L2
            True
        """
        var = str(var).strip()
        return super(MonomialGrowthGroup, cls).__classcall__(
            cls, base, var, category)


    @sage.misc.superseded.experimental(trac_number=17601)
    def __init__(self, base, var, category):
        r"""
        For more information see :class:`MonomialGrowthGroup`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'x')
            Growth Group x^ZZ
            sage: agg.MonomialGrowthGroup(QQ, SR.var('n'))
            Growth Group n^QQ
            sage: agg.MonomialGrowthGroup(ZZ, ZZ['y'].gen())
            Growth Group y^ZZ
            sage: agg.MonomialGrowthGroup(QQ, 'log(x)')
            Growth Group log(x)^QQ

        TESTS::

            sage: agg.MonomialGrowthGroup('x', ZZ)
            Traceback (most recent call last):
            ...
            TypeError: x is not a valid base
        """
        if not var:
            raise ValueError('Empty var is not allowed.')
        if var[0] in '0123456789=+-*/^%':
            # This restriction is mainly for optical reasons on the
            # representation. Feel free to relax this if needed.
            raise ValueError("The variable name '%s' is inappropriate." %
                             (var,))
        self._var_ = var

        super(MonomialGrowthGroup, self).__init__(category=category, base=base)


    def _repr_short_(self):
        r"""
        A short representation string of this monomial growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'a')  # indirect doctest
            Growth Group a^ZZ


        TESTS::

            sage: agg.MonomialGrowthGroup(ZZ, 'a')._repr_short_()
            'a^ZZ'
            sage: agg.MonomialGrowthGroup(QQ, 'a')._repr_short_()
            'a^QQ'
            sage: agg.MonomialGrowthGroup(PolynomialRing(QQ, 'x'), 'a')._repr_short_()
            'a^(Univariate Polynomial Ring in x over Rational Field)'
        """
        return '%s^%s' % (self._var_, parent_to_repr_short(self.base()))


    def __hash__(self):
        r"""
        Return the hash of this group.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: hash(P)  # random
            -1234567890123456789
        """
        return hash((super(MonomialGrowthGroup, self).__hash__(), self._var_))


    def _convert_(self, data):
        r"""
        Converts given ``data`` to something the constructor of the
        element class accepts (``raw_element``).

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        An element of the base ring or ``None`` (when no such element
        can be constructed).

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P._convert_('icecream') is None
            True
            sage: P(1)  # indirect doctest
            1
            sage: P('x')  # indirect doctest
            x

        ::

            sage: P(x)  # indirect doctest
            x
            sage: P(x^-333)  # indirect doctest
            x^(-333)
            sage: P(log(x)^2)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert log(x)^2.

        ::

            sage: PR.<x> = ZZ[]; x.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: P(x^2)  # indirect doctest
            x^2

        ::

            sage: PSR.<x> = ZZ[[]]
            sage: P(x^42)  # indirect doctest
            x^42
            sage: P(x^12 + O(x^17))
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert x^12 + O(x^17).

        ::

            sage: R.<w,x> = ZZ[]
            sage: P(x^4242)  # indirect doctest
            x^4242
            sage: P(w^4242)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert w^4242.

        ::

            sage: PSR.<w,x> = ZZ[[]]
            sage: P(x^7)  # indirect doctest
            x^7
            sage: P(w^7)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert w^7.

        ::

            sage: P('x^7')
            x^7
            sage: P('1/x')
            1/x
            sage: P('x^(-2)')
            x^(-2)
            sage: P('x^-2')
            x^(-2)
        """
        if data == 1:
            return self.base().zero()
        if str(data) == self._var_:
            return self.base().one()

        try:
            P = data.parent()
        except AttributeError:
            if self._var_ not in str(data):
                return  # this has to end here

            elif str(data) == '1/' + self._var_:
                return self.base()(-1)
            elif str(data).startswith(self._var_ + '^'):
                return self.base()(str(data).replace(self._var_ + '^', '')
                                   .replace('(', '').replace(')', ''))
            else:
                return  # end of parsing


        from sage.symbolic.ring import SR
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
        from sage.rings.polynomial.multi_polynomial_ring_generic import \
            MPolynomialRing_generic
        from sage.rings.power_series_ring import PowerSeriesRing_generic
        import operator
        if P is SR:
            if data.operator() == operator.pow:
                base, exponent = data.operands()
                if str(base) == self._var_:
                    return exponent
        elif isinstance(P, (PolynomialRing_general, MPolynomialRing_generic)):
            if data.is_monomial() and len(data.variables()) == 1:
                if self._var_ == str(data.variables()[0]):
                    return data.degree()
        elif isinstance(P, PowerSeriesRing_generic):
            if hasattr(data, 'variables') and len(data.variables()) == 1:
                from sage.rings.integer_ring import ZZ
                if data.is_monomial() and data.precision_absolute() not in ZZ:
                    if self._var_ == str(data.variables()[0]):
                        return data.degree()
            elif self._var_ == str(data.variable()[0]):
                from sage.rings.integer_ring import ZZ
                if data.is_monomial() and data.precision_absolute() not in ZZ:
                    return data.degree()


    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this growth group.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P_x_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_x_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: bool(P_x_ZZ.has_coerce_map_from(P_x_QQ))  # indirect doctest
            False
            sage: bool(P_x_QQ.has_coerce_map_from(P_x_ZZ))  # indirect doctest
            True
            sage: P_y_ZZ = agg.MonomialGrowthGroup(ZZ, 'y')
            sage: bool(P_y_ZZ.has_coerce_map_from(P_x_ZZ))  # indirect doctest
            False
            sage: bool(P_x_ZZ.has_coerce_map_from(P_y_ZZ))  # indirect doctest
            False
            sage: bool(P_y_ZZ.has_coerce_map_from(P_x_QQ))  # indirect doctest
            False
            sage: bool(P_x_QQ.has_coerce_map_from(P_y_ZZ))  # indirect doctest
            False
        """
        if super(MonomialGrowthGroup, self)._coerce_map_from_(S):
            if self._var_ == S._var_:
                return True


    def gens_monomial(self):
        r"""
        Return a tuple containing generators of this growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple containing elements of this growth group.

        .. NOTE::

            If a :class:`MonomialGrowthGroup` models a logarithmic
            growth group (by having a variable name of the form
            ``log(...)``), an empty tuple is returned.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'x').gens_monomial()
            (x,)
            sage: agg.MonomialGrowthGroup(QQ, 'log(x)').gens_monomial()
            ()
        """
        if self._var_.startswith('log(') and self._var_.endswith(')'):
            return ()
        return (self(raw_element=self.base().one()),)


    def gens(self):
        r"""
        Return a tuple of all generators of this monomial growth
        group, even if the growth group is logarithmic.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple whose entries are instances of
        :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.gens()
            (x,)
            sage: agg.MonomialGrowthGroup(ZZ, 'log(x)').gens()
            (log(x),)
        """
        return (self(raw_element=self.base().one()),)


    def gen(self, n=0):
        r"""
        Return the `n`-th generator of this growth group.

        INPUT:

        - ``n`` -- default: `0`.

        OUTPUT:

        A :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.gen()
            x
        """
        return self.gens()[n]

    def ngens(self):
        r"""
        Return the number of generators of this monomial growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A Python integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.ngens()
            1
            sage: agg.MonomialGrowthGroup(ZZ, 'log(x)').ngens()
            1
        """
        return len(self.gens())


class GrowthGroupFactory(sage.structure.factory.UniqueFactory):
    r"""
    A factory creating asymptotic growth groups.

    INPUT:

    - ``specification`` -- a string.

    OUTPUT:

    An asymptotic growth group.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: agg.GrowthGroup('x^ZZ')
        Growth Group x^ZZ
        sage: agg.GrowthGroup('log(x)^QQ')
        Growth Group log(x)^QQ

    This factory can also be used to construct Cartesian products
    of growth groups::

        sage: agg.GrowthGroup('x^ZZ * y^ZZ')
        Growth Group x^ZZ * y^ZZ
        sage: agg.GrowthGroup('x^ZZ * log(x)^ZZ')
        Growth Group x^ZZ * log(x)^ZZ
        sage: agg.GrowthGroup('x^ZZ * log(x)^ZZ * y^QQ')
        Growth Group x^ZZ * log(x)^ZZ * y^QQ
    """
    def create_key_and_extra_args(self, specification, **kwds):
        r"""
        Given the arguments and keyword, create a key that uniquely
        determines this object.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GrowthGroup.create_key_and_extra_args('x^ZZ')
            (('x^ZZ',), {})
            sage: agg.GrowthGroup.create_key_and_extra_args('asdf')
            Traceback (most recent call last):
            ...
            ValueError: 'asdf' is not a valid string describing a growth group.
        """
        factors = tuple(s.strip() for s in specification.split('*'))
        for f in factors:
            if '^' not in f:
                raise ValueError("'%s' is not a valid string describing "
                                 "a growth group." % (f,))

        return factors, kwds


    def create_object(self, version, factors, **kwds):
        r"""
        Create an object from the given arguments.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GrowthGroup('as^df')
            Traceback (most recent call last):
            ...
            ValueError: 'as^df' is not a valid string describing a growth group.
            sage: agg.GrowthGroup('x^y^z')
            Traceback (most recent call last):
            ...
            ValueError: Cannot decode x^y^z.
        """

        groups = []
        for factor in factors:
            b_and_e = factor.split('^')
            if len(b_and_e) != 2:
                raise ValueError('Cannot decode %s.' % (factor,))
            (b, e) = b_and_e

            try:
                # monomial growth group: 'var^base'
                groups.append(
                    MonomialGrowthGroup(repr_short_to_parent(e), b, **kwds))
                continue
            except (TypeError, ValueError):
                pass

            raise ValueError("'%s' is not a valid string describing "
                             "a growth group." % (factor,))
            # todo: once exponential growth groups are implemented,
            #       move line above to the bottom of this loop

            try:
                # exponential growth group: 'base^var'
                groups.append(
                    ExponentialGrowthGroup(repr_short_to_parent(b), e, **kwds))
                continue
            except (TypeError, ValueError):
                pass

        if len(groups) == 1:
            return groups[0]

        # otherwise, a cartesian product is created. the growth elements
        # in this products are
        # - ordered lexicographically (over the same variable (or a
        #   function thereof)
        # - ordered component-wise (for different variables)

        if len(set(groups)) != len(groups):
            raise ValueError('The cartesian product of equal growth '
                             'groups is not supported')

        equal_var_groups = []
        vars = []
        for group in groups:
            var = group._var_
            if not vars or '(' + vars[-1] + ')' not in var:
                vars.append(var)
                equal_var_groups.append([group])
            else:
                equal_var_groups[-1].append(group)

        from sage.categories.cartesian_product import cartesian_product
        for k in range(len(equal_var_groups)):
            if len(equal_var_groups[k]) > 1:
                equal_var_groups[k] = cartesian_product(equal_var_groups[k],
                                                        order='lex')
            else:
                equal_var_groups[k] = equal_var_groups[k][0]

        if len(equal_var_groups) == 1:
            return equal_var_groups[0]

        return cartesian_product(equal_var_groups, order='components')

GrowthGroup = GrowthGroupFactory("GrowthGroup")
