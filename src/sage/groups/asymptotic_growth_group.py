r"""
(Asymptotic) Growth Groups

AUTHORS:

- Benjamin Hackl (2015-01): initial version
- Daniel Krenn (2015-05-29): initial version and review
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

class GenericGrowthElement(sage.structure.element.MultiplicativeGroupElement):
    r"""
    Class for a generic asymptotic growth group element. These
    elements hold exactly one asymptotic term and can be compared to
    each other, multiplied and divided, but possess no explicit
    coefficient. At this stage, just the order of magnitude shall be
    managed. In this class, only base structure is handled.
    For a concrete realization see :class:`MonomialGrowthElement`.

    INPUT:

    - ``parent`` -- a :class:`GenericGrowthGroup`.

    - ``raw_element`` -- an element from the base ring of the parent.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: G = agg.GenericGrowthGroup(ZZ)
        sage: g = agg.GenericGrowthElement(G, 42); g
        GenericGrowthElement(42)
        sage: g.parent()
        Generic Growth Group over Integer Ring
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
            Category of elements of Generic Growth Group over Integer Ring

        ::

            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G(raw_element=42).category()
            Category of elements of Generic Growth Group over Integer Ring

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
            NotImplementedError: Only implemented in concrete realizations
        """
        raise NotImplementedError('Only implemented in concrete realizations')


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

        TESTS::

            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: G.an_element() <= G.an_element()
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert 1.
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
            sage: P_ZZ = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P_QQ = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P_ZZ.gen() <= P_QQ.gen()^2  # indirect doctest
            True
        """
        return (self / other).is_le_one()


    def is_le_one(self):
        r"""
        Abstract method for comparison with one.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: G.gen().is_le_one()
            False
            sage: (~G.gen()).is_le_one()
            True
        """
        raise NotImplementedError('Only implemented in concrete realizations')


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
        Generic Growth Group over Integer Ring

    .. SEEALSO::

        :class:`MonomialGrowthGroup`
    """
    # TODO: implement some sort of 'assume', where basic assumptions
    # for the variables can be stored. --> within the cartesian product

    # enable the category framework for elements
    Element = GenericGrowthElement


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
        """
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


    def _repr_(self):
        r"""
        A representation string for this generic growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GenericGrowthGroup(ZZ)  # indirect doctest
            Generic Growth Group over Integer Ring
        """
        return 'Generic Growth Group over %s' % (self.base(),)


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
        Return if the growth element ``left`` is at most (less than or
        equal to) the growth element ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

        Return whether the asymptotic order of magnitude of `x` is less
        than or equal to the asymptotic order of magnitude of `y`.

        INPUT:

        - ``x``, ``y`` -- elements of ``self``.

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
            sage: G1 = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: G2 = agg.MonomialGrowthGroup(QQ, 'x')
            sage: bool(G1._coerce_map_from_(G2))
            False
            sage: bool(G2._coerce_map_from_(G1))
            True
        """
        if isinstance(S, GenericGrowthGroup):
            if self.base().has_coerce_map_from(S.base()):
                return True


class MonomialGrowthElement(GenericGrowthElement):
    r"""
    Class for the concrete realization of asymptotic growth group
    elements in the case of polynomial growth. These elements hold
    exactly one asymptotic term.

    A power growth element represents a polynomial term
    `\operatorname{variable}^{\operatorname{exponent}}`.
    More complex constructions including logarithmic or exponential
    terms can be constructed via a cartesian product of the related
    growth groups. Asymptotic growth elements can be multiplied,
    divided, inverted, and compared to each other. However, they
    possess no explicit coefficient.

    The elements can be specified by either an expression ``x`` being
    a string, an element from the symbolic or a polynomial ring or the
    integer `1`. On the other hand, elements can also be specified
    by their exponent.

    INPUT:

    - ``parent`` -- a :class:`MonomialGrowthGroup`, the
      parent of the element.
    - ``x`` -- an expression (string, polynomial ring element,
      symbolic ring element, or the integer `1`) representing
      the element to be initialized.
    - ``exponent`` -- the exponent of the power element.

    OUTPUT:

    An asymptotic power growth element with the specified
    parent and magnitude of growth, i.e. exponent.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
        sage: e1 = P(x=1); e1
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
        return self._raw_element_


    def _repr_(self):
        r"""
        Represent the asymptotic power growth element as
        ``variable^{exponent}``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(QQ, 'x')
            sage: P(x=1)._repr_()
            '1'
            sage: P(raw_element=5)._repr_()
            'x^5'
            sage: P(raw_element=1/2)._repr_()
            'x^(1/2)'
        """
        from sage.rings.integer_ring import ZZ

        if self.exponent == 0:
            return '1'
        elif self.exponent == 1:
            return self.parent()._var_
        elif self.exponent in ZZ:
            return self.parent()._var_ + '^' + str(self.exponent)
        else:
            return self.parent()._var_ + '^(' + str(self.exponent) + ')'


    def _mul_(self, other):
        r"""
        Multiply two asymptotic power growth elements from the
        same parent by adding their exponents.

        INPUT:

        - ``other`` -- the asymptotic growth element to be
          multiplied with ``self``.

        OUTPUT:

        An asymptotic power growth element representing the product
        of ``self`` and ``other``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: a = P(raw_element=2)
            sage: b = P(raw_element=3)
            sage: c = a._mul_(b); c
            x^5
            sage: c == a * b
            True
            sage: a * b * a
            x^7
        """
        cls = self.__class__
        return cls(self.parent(), self.exponent + other.exponent)


    def __invert__(self):
        r"""
        Return the multiplicative inverse from a given asymptotic power
        growth element.

        INPUT:

        Nothing.

        OUTPUT:

        The multiplicative inverse asymptotic power growth element
        of ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: e1 = P(raw_element=2)
            sage: e2 = e1.__invert__(); e2
            x^-2
            sage: e2 == ~e1
            True
        """
        cls = self.__class__
        return cls(self.parent(), -self.exponent)


    def _div_(self, other):
        r"""
        Divide two asymptotic power growth elements from the same
        parent by subtracting their exponents.

        INPUT:

        - ``other`` -- the asymptotic growth element which ``self``
          is divided by.

        OUTPUT:

        The result of the division of ``self`` by ``other``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: e1 = P(raw_element=2)
            sage: e1._div_(P.gen())
            x
            sage: e1._div_(P.gen()) == e1 / P.gen()
            True
        """
        cls = self.__class__
        return cls(self.parent(), self.exponent - other.exponent)


    def __pow__(self, power):
        r"""
        Return a asymptotic power element to the power of
        ``power``.

        INPUT:

        - ``power`` -- a rational number.

        OUTPUT:

        The asymptotic power growth element ``self`` to the power of
        ``power``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.gen().__pow__(5)
            x^5
            sage: P.gen().__pow__(1/2)
            x^(1/2)
            sage: P.gen()^7
            x^7
        """
        new_exponent = self.exponent * power
        if new_exponent in self.parent().base():
            return self.parent()(raw_element=self.exponent * power)

        from sage.rings.real_mpfr import RR
        if new_exponent in RR:
            pnt = MonomialGrowthGroup(new_exponent.parent(),
                                   self.parent()._var_)
            return pnt(raw_element=new_exponent)
        else:
            raise NotImplementedError('Only real exponents are implemented.')


    def is_le_one(self):
        r"""
        Return whether or not the growth of the asymptotic power
        growth element ``self`` is less than or equal to the
        (constant) growth of `1`.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: e1 = P.gen()
            sage: e1.is_le_one()
            False
            sage: (P.one() / P.gen()).is_le_one()
            True
        """
        return self.exponent <= 0


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
        Monomial Growth Group in x over Integer Ring
        sage: agg.MonomialGrowthGroup(ZZ, log(SR.var('y')))
        Monomial Growth Group in log(y) over Integer Ring

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
            sage: P2 = agg.MonomialGrowthGroup(ZZ, PolynomialRing(ZZ, 'x').gen())
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


    def __init__(self, base, var, category):
        r"""
        For more information see :class:`MonomialGrowthGroup`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'x')
            Monomial Growth Group in x over Integer Ring
            sage: agg.MonomialGrowthGroup(QQ, SR.var('n'))
            Monomial Growth Group in n over Rational Field
            sage: agg.MonomialGrowthGroup(ZZ, PolynomialRing(ZZ, 'y').gen())
            Monomial Growth Group in y over Integer Ring
            sage: agg.MonomialGrowthGroup(QQ, 'log(x)')
            Monomial Growth Group in log(x) over Rational Field
        """
        if not var:
            raise ValueError('Empty var is not allowed.')
        if var[0] in '0123456789=+-*/^%':
            # This restriction is mainly for optical reasons on the
            # representation. Feel free to relax this if needed.
            raise ValueError('The inapproproate variable name %s.' % (var,))
        self._var_ = var

        super(MonomialGrowthGroup, self).__init__(category=category, base=base)


    def _repr_(self):
        r"""
        A representation string for this monomial growth group.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.MonomialGrowthGroup(ZZ, 'x')  # indirect doctest
            Monomial Growth Group in x over Integer Ring

        TESTS::

            sage: agg.MonomialGrowthGroup(QQ, 'v_107')._repr_()
            'Monomial Growth Group in v_107 over Rational Field'
        """
        return 'Monomial Growth Group in %s over %s' % (self._var_, self.base())


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
            x^-333
            sage: P(log(x)^2)  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Cannot convert log(x)^2.
        """
        if data == 1:
            return self.base().zero()
        if str(data) == self._var_:
            return self.base().one()

        try:
            P = data.parent()
        except AttributeError:
            return  # this has to end here

        from sage.symbolic.ring import SR
        import operator
        if P is SR:
            if data.operator() == operator.pow:
                base, exponent = data.operands()
                if str(base) == self._var_:
                    return exponent
        #elif ...
        #TODO: PolynomialRing, PowerSeriesRing


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
            sage: bool(P_x_ZZ._coerce_map_from_(P_x_QQ))
            False
            sage: bool(P_x_QQ._coerce_map_from_(P_x_ZZ))
            True
            sage: P_y_ZZ = agg.MonomialGrowthGroup(ZZ, 'y')
            sage: bool(P_y_ZZ._coerce_map_from_(P_x_ZZ))
            False
            sage: bool(P_x_ZZ._coerce_map_from_(P_y_ZZ))
            False
            sage: bool(P_y_ZZ._coerce_map_from_(P_x_QQ))
            False
            sage: bool(P_x_QQ._coerce_map_from_(P_y_ZZ))
            False
        """
        if super(MonomialGrowthGroup, self)._coerce_map_from_(S):
            if self._var_ == S._var_:
                return True


    def gen(self):
        r"""
        Return the monomial growth element with exponent `1`.

        INPUT:

        Nothing.

        OUTPUT:

        A :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: e1 = agg.MonomialGrowthGroup(ZZ, 'x').gen(); e1
            x
            sage: e1.exponent == 1
            True
        """
        return self(raw_element=self.base().one())


    def gens(self):
        r"""
        Return a tuple of all generators of this monomial growth
        group, which is exactly consisting of one element, namely the
        monomial growth element with exponent `1`.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple whose entries are :class:`MonomialGrowthElement`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.gens()
            (x,)
        """
        return (self.gen(),)


    def ngens(self):
        r"""
        Return the number of generators of this monomial growth group,
        which is exactly `1`.

        INPUT:

        Nothing.

        OUTPUT:

        A Python integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: P.ngens()
            1
        """
        return 1
