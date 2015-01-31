r"""
Asymptotic Term Monoid

AUTHORS:

- Benjamin Hackl (2015-01): initial version

"""

# *****************************************************************************
# Copyright (C) 2014--2015 Benjamin Hackl <benjamin.hackl@aau.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************
from sage.structure.element import MonoidElement

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class GenericTerm(MonoidElement):
    r"""
    Generic Asymptotic Term class.

    TODO: Description.
    """

    def __init__(self, parent, growth):
        r"""
        See :class:`GenericTerm` for more information.

        INPUT:

        - ``parent`` -- the parent of the asymptotic term.

        - ``growth`` -- an element of the parent's ``base``,
          specifying the growth of the asymptotic term.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: P = atm.GenericTermMonoid(PG)
            sage: P(x^2)
            Generic Term with growth x^2
        """
        from sage.groups.asymptotic_growth_group import GenericGrowthElement

        if parent is None:
            raise ValueError("The parent must be provided")
        if growth is None or not isinstance(growth, GenericGrowthElement):
            raise ValueError("The growth must be provided and must inherit "
                             "from GenericGrowthElement")
        else:
            if growth not in parent.base():
                raise ValueError("%s is not in the parent's "
                                 "specified base" % growth)
            else:
                self.growth = growth
        super(GenericTerm, self).__init__(parent=parent)


    def _mul_(self, other):
        r"""
        Abstract multiplication method for generic terms.

        INPUT:

        - ``other`` -- a asymptotic term from the same parent as
          ``self``.

        OUTPUT:

        A :class:`GenericTerm` representing the product of ``self``
        and ``other``.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: P = atm.GenericTermMonoid(PG)
            sage: t1, t2 = P(x), P(x^2)
            sage: t1, t2
            (Generic Term with growth x, Generic Term with growth x^2)
            sage: t1 * t2  # indirect doctest
            Generic Term with growth x^3
        """
        cls = self.__class__
        return cls(self.parent(), self.growth * other.growth)


    def can_absorb(self, other):
        r"""
        Check, whether the asymptotic term ``self`` is able to absorb
        the asymptotic term ``other``. For example, an :class:`OTerm`
        is able to *absorb* another :class:`OTerm` or an
        :class:`ExactTerm` with lower (or equal) growth. For more
        information see :class:`OTerm`, :class:`ExactTerm`, and
        :class:`LTermGeneric`.

        INPUT:

        - ``other`` -- some asymptotic term.

        OUTPUT:

        A boolean.

        EXAMPLES:

        We want to show step by step which terms can be absorbed
        by which other terms. We start by defining the necessary
        term monoids and some terms::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: OT = atm.OTermMonoid(base=PG)
            sage: ET = atm.ExactTermMonoid(base=PG, coefficient_ring=QQ)
            sage: LT = atm.LTermGenericMonoid(base=PG, coefficient_ring=QQ)
            sage: ot1, ot2 = OT(x), OT(x^2)
            sage: et1 = ET(x^2, 2)
            sage: lt1 = LT(x^2, 2, start=0)

        :class:`OTerm` is able to absorb other :class:`OTerm`,
        :class:`LTermGeneric` (and descendants thereof) as well
        as :class:`ExactTerm`, as long as the growth of the other
        term is less than or equal to the growth of ``self``::

            sage: ot1, ot2
            (O(x), O(x^2))
            sage: ot1.can_absorb(ot2), ot2.can_absorb(ot1)
            (False, True)
            sage: et1
            2 * x^2
            sage: ot1.can_absorb(et1)
            False
            sage: ot2.can_absorb(et1)
            True
            sage: lt1
            2 * L(x^2, 0)
            sage: ot1.can_absorb(lt1)
            False
            sage: ot2.can_absorb(lt1)
            True

        :class:`ExactTerm` can only absorb another
        :class:`ExactTerm` if the growth coincides with the
        growth of ``self``::

            sage: et1.can_absorb(ET(x^2, 5))
            True
            sage: any(et1.can_absorb(t) for t in [ot1, ot2, lt1])
            False

        :class:`LTermGeneric` can absorb arbitrary other
        :class:`LTermGeneric`, and :class:`ExactTerm` whose
        growth is less than or equal to the growth of ``self``::

            sage: any(lt1.can_absorb(t) for t in [ot1, ot2])
            False
            sage: lt1.can_absorb(LT(x^5, 1, start=0))
            True
            sage: lt1.can_absorb(et1)
            True
        """
        if not isinstance(other, GenericTerm):
            raise TypeError("%s is not an asymptotic term" % other)
        if isinstance(self, OTerm):
            if isinstance(other, (OTerm, LTermGeneric, ExactTerm)):
                return other <= self
            else:
                return False
        elif isinstance(self, LTermGeneric):
            if isinstance(other, ExactTerm):
                return other <= self
            elif isinstance(other, LTermGeneric):
                return True
            else:
                return False
        elif isinstance(self, ExactTerm):
            if isinstance(other, ExactTerm):
                return self.growth == other.growth
            else:
                return False


    def absorb(self, other):
        r"""
        Absorb the asymptotic term ``other``, yielding a new
        asymptotic term (or ``None``).

        INPUT:

        - ``other`` -- an asymptotic term that can be absorbed
          by ``self``.

        OUTPUT:

        An asymptotic term.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: # todo: write examples!
        """
        # TODO: implement this.


    def __le__(self, other):
        r"""
        Return whether the growth of the term ``self`` is less than
        or equal to the growth of the term ``other``.

        INPUT:

        - ``other`` -- an asymptotic term (inherited from
          :class:`GenericTerm`) to be compared to ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: # todo: examples
        """
        from sage.structure.sage_object import have_same_parent

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
        Return whether the generic term ``self`` is of lesser or
        equal growth as the generic term ``other`` by calling
        :meth:`GenericTermMonoid.le`.

        INPUT:

        - ``other`` -- a growth power element from the same parent
          to be compared to ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: # todo: examples
        """
        return self.parent().le(self, other)


    def _repr_(self):
        r"""
        Represent the generic term as ``Generic Term with growth
        ...``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: P = atm.GenericTermMonoid(base=PG)
            sage: P(x)._repr_()
            'Generic Term with growth x'
            sage: P(x^7)._repr_()
            'Generic Term with growth x^7'
        """
        return "Generic Term with growth " + repr(self.growth)


class GenericTermMonoid(Parent, UniqueRepresentation):
    r"""
    Parent for generic asymptotic terms. In this class the base
    structure for asymptotic term monoids will be handled. These
    monoids are the parents for asymptotic terms (for example, see
    :class:`GenericTerm` or :class:`OTerm`). Basically, asymptotic
    terms consist of a ``growth`` (which is an asymptotic growth
    group element, for example
    :class:`~sage.groups.asymptotic_growth_group.GrowthElementPower`).

    INPUT:

    - ``base`` -- an asymptotic growth group. This is the parent of
      the elements specifying the growth of the asymptotic terms.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a subcategory
      of ``Join of Category of Monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    OUTPUT:

    A generic asymptotic term monoid.

    EXAMPLES::

        sage:  # todo: write example.
    """

    # enable the category framework for elements
    Element = GenericTerm


    def __init__(self, base=None, category=None):
        r"""
        See :class:`GenericTermMonoid` for more information.
        
        EXAMPLES::
        
            sage: # not tested
        """
        # TODO: examples.

        from sage.categories.monoids import Monoids
        from sage.categories.posets import Posets
        from sage.groups.asymptotic_growth_group import GenericGrowthGroup

        if category is None:
            category = Monoids() & Posets()
        else:
            if not isinstance(category, tuple):
                category = (category, )
            if not any(cat.is_subcategory(Monoids() & Posets()) for cat in
                       category):
                raise ValueError("%s is not a subcategory of %s"
                                 % (category, Monoids() & Posets()))
        if base is None:
            base = GenericGrowthGroup()
        else:
            if not isinstance(base, GenericGrowthGroup):
                raise ValueError("%s does not inherit from %s"
                                 % (base, GenericGrowthGroup()))
        super(GenericTermMonoid, self).__init__(category=category, base=base)


    def _repr_(self):
        r"""
        Represent the generic term monoid as "Generic Term Monoid
        over `base`".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: atm.GenericTermMonoid()._repr_()
            'Generic Term Monoid over Generic Asymptotic Growth Group'
            sage: PG = agg.GrowthGroupPower("x")
            sage: atm.GenericTermMonoid(base=PG)._repr_()
            'Generic Term Monoid over Asymptotic Power Growth Group in x over Integer Ring'
        """
        return "Generic Term Monoid over %s" % repr(self.base())


    def _coerce_map_from_(self, S):
        r"""
        Another GenericTermMonoid ``S`` coerces into ``self`` if and
        only if the base of ``S`` coerces into the base of ``self``.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG1 = agg.GrowthGroupPower("x", base=ZZ)
            sage: T1 = atm.GenericTermMonoid(base=PG1); T1
            Generic Term Monoid over Asymptotic Power Growth Group in x over Integer Ring
            sage: PG2 = agg.GrowthGroupPower("x", base=QQ)
            sage: T2 = atm.GenericTermMonoid(base=PG2); T2
            Generic Term Monoid over Asymptotic Power Growth Group in x over Rational Field
            sage: T2._coerce_map_from_(T1)
            True
        """
        if isinstance(S, self.__class__):
            if self.base().coerce_map_from(S.base()) is not None:
                return True

    def _element_constructor_(self, x):
        r"""
        Find the coercion of object ``x`` in ``self``.

        INPUT:

        - ``x`` -- either an asymptotic term to be coerced into
          ``self``, or an asymptotic growth element used for the
          construction of an asymptotic term in ``self`` with
          growth ``x``.

        OUTPUT:

        An element of the term monoid ``self``.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG1.<x> = agg.GrowthGroupPower(base=ZZ)
            sage: PG2.<x> = agg.GrowthGroupPower(base=QQ)
            sage: P1 = atm.GenericTermMonoid(base=PG1)
            sage: P2 = atm.GenericTermMonoid(base=PG2)
            sage: term1 = P1(PG1.gen())
            sage: term2 = P2(PG2.gen()^2)
            sage: term1 <= term2  # in order for two terms to be compared,
            ....:                 # a coercion into one of the parents
            ....:                 # has to be found.
            True
            sage: P2.coerce(term1)  # coercion does not fail
            Generic Term with growth x
        """
        from sage.groups.asymptotic_growth_group import GenericGrowthElement

        if isinstance(x, GenericGrowthElement):
            return self.element_class(self, x)

        if x is None:
            raise ValueError("The growth of the term has to be specified!")
        elif x.parent() is self:
            return x
        elif isinstance(x, GenericTerm):
            return self.element_class(self, self.base().coerce(x.growth))
        else:
            raise ValueError("Input is not an asymptotic growth element.")


    def le(self, x, y):
        r"""
        Return whether the growth of term `x` is less than or equal
        to the growth of term `y`.

        INPUT:

        - ``x``, ``y`` -- elements of ``self``.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG = agg.GrowthGroupPower("x")
            sage: P = atm.GenericTermMonoid(base=PG)
            sage: x, y = P(PG.gen()), P(PG.gen()^2)
            sage: P.le(x,y)
            True
        """
        return x.growth <= y.growth


class OTerm(GenericTerm):
    r"""
    Class for an asymptotic term representing a *Big O* term with
    specified growth.

    TODO: Description.
    """

    def _repr_(self):
        r"""
        Represent ``self`` as ``O(growth)``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: OT = atm.OTermMonoid(PG)
            sage: t1, t2, t3 = OT(x), OT(x^2), OT(x^3)
            sage: t1._repr_(), t2._repr_()
            ('O(x)', 'O(x^2)')
            sage: t3  # indirect doctest
            O(x^3)
        """
        return "O(%s)" % self.growth


class OTermMonoid(GenericTermMonoid):
    r"""
    Parent for asymptotic big `O` terms.

    INPUT:

    - ``base`` -- the underlying asymptotic growth group, compare
      :class:`GenericTermMonoid`.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    OUTPUT:

    A monoid for asymptotic `O` terms.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage:  # todo: examples!
    """
    # enable the category framework for elements
    Element = OTerm


    def _coerce_map_from_(self, S):
        r"""
        In order for ``S`` to coerce into ``self``, ``S`` may be
        an instance of the following classes:

        - :class:`OTermMonoid`

        - :class:`LTermGenericMonoid` or a descendant thereof

        - :class:`ExactTermMonoid`

        In all those cases, ``S`` coerces into ``self`` if the
        base of ``S`` coerces into the base of ``self``.

        INPUT:

        - ``S`` -- a parent which is tested for coercion into
          ``self``.

        OUTPUT:

        ``True`` or ``None``.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage:  # todo: examples
        """
        if isinstance(S, (ExactTermMonoid, LTermGenericMonoid)):
            if self.base().coerce_map_from(S.base()) is not None:
                return True
        else:
            return super(OTermMonoid, self)._coerce_map_from_(S)


    def _repr_(self):
        r"""
        Represent ``self`` as "Big O term monoid over ``base``".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: atm.OTermMonoid(PG)._repr_()
            'Big O term monoid over Asymptotic Power Growth Group in x over Integer Ring'
        """
        return "Big O term monoid over %s" % self.base()


class TermWithCoefficient(GenericTerm):
    r"""
    Base class for asymptotic terms possessing a coefficient. For
    example, the classes :class:`ExactTerm` and :class:`LTermGeneric`
    both directly inherit from this class.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element from
      ``parent.base()``.

    - ``coefficient`` -- an element from ``parent.coefficient_ring``.

    OUTPUT:

    An asymptotic term with coefficient.

    EXAMPLES::

        sage:  # todo: examples
    """

    def __init__(self, parent, growth, coefficient=1):
        r"""
        See :class:`TermWithCoefficient` for more information.

        EXAMPLES::

            sage:  # todo: examples
        """
        if coefficient not in parent.coefficient_ring:
            raise ValueError("%s is not in %s" % (coefficient,
                                                  parent.coefficient_ring))
        else:
            self.coefficient = coefficient
        super(TermWithCoefficient, self).__init__(parent=parent, growth=growth)


    def _repr_(self):
        r"""
        Represent the :class:`TermWithCoefficient` as "Asymptotic
        Term with coefficient ``coefficient`` and growth ``growth``".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: P = atm.TermWithCoefficientMonoid(PG, ZZ)
            sage: P(x^2, 5)._repr_()
            'Asymptotic Term with coefficient 5 and growth x^2'
        """
        return "Asymptotic Term with coefficient %s and growth %s" % \
               (self.coefficient, self.growth)

    def _mul_(self, other):
        r"""
        Multiplication method for asymptotic terms with coefficients.

        INPUT:

        - ``other`` -- an asymptotic term from the same parent
          as ``self``.

        OUTPUT:

        An asymptotic term representing the product of ``self`` and
        ``other``.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: CT = atm.TermWithCoefficientMonoid(PG, ZZ)
            sage: ET = atm.ExactTermMonoid(PG, ZZ)

            This method handels the multiplication of abstract terms
            with coefficient (i.e. :class:`TermWithCoefficient`) and
            exact terms (i.e. :class:`ExactTerm`). First, an example
            for abstract terms::

            sage: t1, t2 = CT(x^2, 2), CT(x^3, 3)
            sage: t1 * t2
            Asymptotic Term with coefficient 6 and growth x^5

            And now, an example for exact terms::

            sage: t1, t2 = ET(x^2, 2), ET(x^3, 3)
            sage: t1 * t2
            6 * x^5
        """
        result = super(TermWithCoefficient, self)._mul_(other)
        result.coefficient = self.coefficient * other.coefficient
        return result


class TermWithCoefficientMonoid(GenericTermMonoid):
    r"""
    This class implements the base structure for parents of
    asymptotic terms possessing a coefficient from some coefficient
    ring. In particular, this is also the parent for
    :class:`TermWithCoefficient`.

    INPUT:

    - ``base`` -- The underlying asymptotic growth group, compare
      :class:`GenericTermMonoid`.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    - ``coefficient_ring`` -- the ring which contains the
      coefficients of the elements.

    OUTPUT:

    A monoid for asymptotic terms with coefficients.

    EXAMPLES::

        sage: # todo: examples.
    """

    # enable the category framework for elements
    Element = TermWithCoefficient

    def __init__(self, base=None, coefficient_ring=None, category=None):
        r"""
        For more information see :class:`TermWithCoefficientMonoid`.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: P1 = atm.TermWithCoefficientMonoid(PG, ZZ); P1
            Monoid for asymptotic terms with coefficients from Integer Ring over Asymptotic Power Growth Group in x over Integer Ring
            sage: P2 = atm.TermWithCoefficientMonoid(PG, QQ); P2
            Monoid for asymptotic terms with coefficients from Rational Field over Asymptotic Power Growth Group in x over Integer Ring
            sage: P2.category()
            Join of Category of monoids and Category of posets
        """
        if coefficient_ring is None:
            raise ValueError("Base ring is not specified")
        self.coefficient_ring = coefficient_ring
        super(TermWithCoefficientMonoid, self).__init__(base=base,
                                                        category=category)


    def _coerce_map_from_(self, S):
        r"""
        :class:`TermWithCoefficientMonoid` is still a generic base
        class handling structure. Thus, ``S`` coerces into ``self``
        if and only if the base and the coefficient ring of ``S``
        coerce into the base and coefficient ring of ``self``,
        respectively.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG1 = agg.GrowthGroupPower("x")
            sage: PG2 = agg.GrowthGroupPower("x", base=QQ)
            sage: M1 = atm.TermWithCoefficientMonoid(PG1, ZZ); M1
            Monoid for asymptotic terms with coefficients from Integer Ring over Asymptotic Power Growth Group in x over Integer Ring
            sage: M2 = atm.TermWithCoefficientMonoid(PG2, QQ); M2
            Monoid for asymptotic terms with coefficients from Rational Field over Asymptotic Power Growth Group in x over Rational Field
            sage: M2._coerce_map_from_(M1)
            True
            sage: M1._coerce_map_from_(M2) is None
            True
        """
        if isinstance(S, TermWithCoefficientMonoid):
            return (super(TermWithCoefficientMonoid, self).
                    _coerce_map_from_(S) and self.coefficient_ring.
                    coerce_map_from(S.coefficient_ring) is not None)


    def _element_constructor_(self, x, coefficient=1):
        r"""
        Construct a asymptotic term with coefficient.

        INPUT:

        - ``x`` -- an asymptotic term with coefficient to be coerced
          into ``self``, or an asymptotic growth group element
          to be used to construct an asymptotic term of ``self``.

        - ``coefficient`` -- an element of the ``coefficient_ring``
          of ``self``.

        OUTPUT:

        An asymptotic term with parent ``self``.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: P = atm.TermWithCoefficientMonoid(PG, ZZ)
            sage: t1 = P(x^2, 5); t1  # indirect doctest
            Asymptotic Term with coefficient 5 and growth x^2
        """
        if isinstance(x, TermWithCoefficient):
            return self.element_class(self, x.growth, x.coefficient)
        else:
            return self.element_class(self, x, coefficient)


    def _repr_(self):
        r"""
        Represent ``self`` as ``Monoid for asymptotic terms with
        coefficients from coefficient_ring over base``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG = agg.GrowthGroupPower("x")
            sage: atm.TermWithCoefficientMonoid(PG, ZZ)._repr_()
            'Monoid for asymptotic terms with coefficients from Integer Ring over Asymptotic Power Growth Group in x over Integer Ring'
        """
        return "Monoid for asymptotic terms with coefficients from %s " \
               "over %s" % (self.coefficient_ring, self.base())


class LTermGeneric(TermWithCoefficient):
    r"""
    Base class for asymptotic `L` terms, i.e. big `O` terms with
    explicitly specified constant.

    .. NOTE::

        When adding ("absorbing") various different `L` terms, some
        sort of "lifting" occurs in general. This lifting affects the
        coefficient of the resulting `L` term, and works
        differently for various growth classes. This makes different
        implementations for different growth classes necessary.

    """

    def __init__(self, parent, growth, coefficient=1, start=0):
        from sage.rings.real_mpfr import RR

        if start not in RR:
            raise ValueError("%s is not a real number" % start)
        else:
            self.start = start
        super(LTermGeneric, self).__init__(parent=parent, growth=growth,
                                           coefficient=coefficient)


    def _repr_(self):
        r"""
        Represent the generic `L` term as "``coefficient * L(growth,
        start)``".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: LT = atm.LTermGenericMonoid(PG, ZZ)
            sage: LT(x^2, 5, start=0)._repr_()
            '5 * L(x^2, 0)'
        """
        return "%s * L(%s, %s)" % (self.coefficient, self.growth, self.start)


class LTermGenericMonoid(TermWithCoefficientMonoid):
    r"""
    Base class for parents of asymptotic `L` terms. Also, the parent
    for :class:`LTermGeneric`.

    TODO: Description
    """
    # enable the category framework for elements
    Element = LTermGeneric


    def _element_constructor_(self, x, coefficient=1, start=0):
        r"""
        Construct a generic `L` term.

        INPUT:

        - ``x`` -- an asymptotic `L` term to be coerced into
          ``self``, or an asymptotic growth group element to be used
           to construct a generic `L` term of ``self``.

        - ``coefficient`` -- an element of the ``coefficient_ring``
          of ``self``.

        - ``start`` -- a real number indicating the point where
          the `L` term is valid.

        OUTPUT:

        A generic `L` term.

        EXAMPLES::

            sage:  # todo: examples.
        """
        if isinstance(x, LTermGeneric):
            return self.element_class(self, x.growth, x.coefficient, x.start)
        elif isinstance(x, ExactTerm):
            return self.element_class(self, x.growth, x.coefficient, 0)
        else:
            return self.element_class(self, x, coefficient, start)

    def _coerce_map_from_(self, S):
        r"""
        Another :class:`LTermGenericMonoid` or
        :class:`ExactTermMonoid` ``S`` coerces into ``self`` if the
        base of ``S`` as well as the coefficient ring of ``S`` coerce
        into the base and the coefficient ring of ``self``,
        respectively.

        INPUT:

        - ``S`` -- a parent which is tested for coercion into
          ``self``.

        OUTPUT:

        ``True`` or ``None``.

        EXAMPLES::

            sage: # todo: examples.
        """
        if isinstance(S, ExactTermMonoid):
            if self.base().coerce_map_from(S.base()) is not None and \
                self.coefficient_ring.coerce_map_from(
                    S.coefficient_ring) is not None:
                return True
        return super(LTermGenericMonoid, self)._coerce_map_from_(S)


class ExactTerm(TermWithCoefficient):
    r"""
    Class for asymptotic exact terms. These terms primarily consist of
    an asymptotic growth element as well as a coefficient specifying
    the growth of the asymptotic term.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element from
      ``parent.base()``.

    - ``coefficient`` -- an element from ``parent.coefficient_ring``.

    OUTPUT:

    An asymptotic exact term.

    EXAMPLES::

        sage:  # todo: examples
    """

    def _repr_(self):
        r"""
        Represent the exact term ``self`` as ``coefficient *
        growth``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: ET = atm.ExactTermMonoid(base=PG, coefficient_ring=ZZ)
            sage: et1 = ET(x^2, 2); et1
            2 * x^2
        """
        return "%s * %s" % (self.coefficient, self.growth)


class ExactTermMonoid(TermWithCoefficientMonoid):
    r"""
    Parent for asymptotic exact terms, implemented in
    :class:`ExactTerm`.

    INPUT:

    - ``base`` -- The underlying asymptotic growth group, compare
      :class:`GenericTermMonoid`.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    - ``coefficient_ring`` -- the ring which contains the
      coefficients of the elements.

    OUTPUT:

    A monoid for asymptotic exact terms.

    EXAMPLES::

        sage:  # todo: examples.
    """
    # enable the category framework for elements
    Element = ExactTerm