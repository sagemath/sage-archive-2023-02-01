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
    Base class for asymptotic terms. Mainly the structure and
    several properties asymptotic terms have are handled here.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an element of the parent's ``base``, specifying
      the growth of the asymptotic term.

    OUTPUT:

    A generic asymptotic term.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: PG.<x> = agg.GrowthGroupPower()
        sage: P = atm.GenericTermMonoid(PG)
        sage: t1, t2 = P(x), P(x^2); (t1, t2)
        (Generic Term with growth x, Generic Term with growth x^2)
        sage: t1 * t2
        Generic Term with growth x^3
        sage: t1.can_absorb(t2)
        False
        sage: t1.absorb(t2)
        Traceback (most recent call last):
        ...
        ArithmeticError: Generic Term with growth x cannot absorb Generic Term with growth x^2
        sage: t1.can_absorb(t1)
        False
    """

    def __init__(self, parent, growth):
        r"""
        See :class:`GenericTerm` for more information.

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
            sage: t1 * t2
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
        else:
            return False


    def absorb(self, other):
        r"""
        Absorb the asymptotic term ``other``, yielding a new
        asymptotic term (or ``None``). For a more detailed
        explanation of the *absorption* of asymptotic terms see
        the introduction of this module, or the following examples.

        INPUT:

        - ``other`` -- an asymptotic term that can be absorbed
          by ``self``.

        OUTPUT:

        An asymptotic term or ``None``.

        EXAMPLES:

        We want to demonstrate in which cases an asymptotic term
        is able to absorb another term, as well as explain the output
        of this operation. We start by defining some parents and
        elements::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower(base=QQ)
            sage: OT = atm.OTermMonoid(base=PG)
            sage: ET = atm.ExactTermMonoid(base=PG, coefficient_ring=QQ)
            sage: LT = atm.LTermGenericMonoid(base=PG, coefficient_ring=QQ)
            sage: ot1, ot2 = OT(x), OT(x^2)
            sage: et1, et2 = ET(x, 100), ET(x^2, 2)
            sage: et3, et4 = ET(x^2, 1), ET(x^2, -2)
            sage: lt1 = LT(x, 5, start=0)

        Because of the definition of `O` terms (see
        :wikipedia:`Big_O_notation`), they are able to absorb all
        other asymptotic terms with weaker or equal growth. The
        result of the absorption is the original `O` Term::

            sage: ot1.absorb(ot1)
            O(x)
            sage: ot1.absorb(et1)
            O(x)
            sage: ot1.absorb(lt1)
            O(x)

        This corresponds to `O(x) + O(x) = O(x)`,
        `O(x) + 100x = O(x)`, and `O(x) + 5\cdot L(x,0) = O(x)`.
        If absorb is called on a term that cannot be absorbed, an
        ``ArithmeticError`` is raised::

            sage: ot1.absorb(ot2)
            Traceback (most recent call last):
            ...
            ArithmeticError: O(x) cannot absorb O(x^2)

        This would only work the other way around::

            sage: ot2.absorb(ot1)
            O(x^2)

        :class:`ExactTerm` is able to absorb another
        :class:`ExactTerm` if the terms have the same growth. In this
        case, *absorption* is nothing else than an addition of the
        respective coefficients::

            sage: et2.absorb(et3)
            3 * x^2
            sage: et3.absorb(et2)
            3 * x^2
            sage: et3.absorb(et4)
            -1 * x^2

        Note that, for technical reasons, the coefficient `0` is not
        allowed, and thus ``None`` is returned if two exact terms
        cancel each other out::

            sage: et2.absorb(et4)
            sage: repr(et4.absorb(et2))
            'None'

        .. TODO:

            The absorption of `L` terms is implemented at a later
            point.
        """
        from sage.structure.sage_object import have_same_parent
        if not self.can_absorb(other):
            raise ArithmeticError("%s cannot absorb %s" % (self, other))

        if have_same_parent(self, other):
            return self._absorb_(other)

        from sage.structure.element import get_coercion_model
        try:
            return get_coercion_model().bin_op(self, other,
                                               lambda self, other:
                                               self._absorb_(other))
        except TypeError:
            return False

    def _absorb_(self, other):
        r"""
        Absorption of ``other`` by ``self``, where ``self`` and
        ``other`` have the same parent. This is not implemented
        for abstract base classes, for concrete realizations
        see :meth:`OTerm._absorb_` or :meth:`ExactTerm._absorb_`.

        INPUT:

        - ``other`` -- an asymptotic term from the same parent as
          ``self``.

        OUTPUT:

        An asymptotic term representing the result of the absorption
        of ``other`` by ``self``, or ``None``.

        EXAMPLES:

        First, we define some asymptotic terms::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: P = atm.GenericTermMonoid(PG)
            sage: t1, t2 = P(x), P(x^2)

        When it comes to absorption, note that the method
        :meth:`can_absorb` (which is called before absorption takes
        place) does not allow the absorption of generic terms. Thus,
        an ``ArithmeticError`` is raised::

            sage: t2.absorb(t1)
            Traceback (most recent call last):
            ...
            ArithmeticError: Generic Term with growth x^2 cannot absorb Generic Term with growth x
        """
        raise NotImplementedError("Not implemented in abstract base classes")


    def __le__(self, other):
        r"""
        Return whether the growth of the term ``self`` is less than
        or equal to the growth of the term ``other``.

        INPUT:

        - ``other`` -- an asymptotic term (inherited from
          :class:`GenericTerm`) to be compared to ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES:

        First, we define some asymptotic terms (and their parents)::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: GT = atm.GenericTermMonoid(PG)
            sage: OT = atm.OTermMonoid(PG)
            sage: ET1 = atm.ExactTermMonoid(PG, ZZ)
            sage: ET2 = atm.ExactTermMonoid(PG, QQ)
            sage: g1, g2 = GT(x), GT(x^2); g1, g2
            (Generic Term with growth x, Generic Term with growth x^2)
            sage: o1, o2 = OT(x^-1), OT(x^3); o1, o2
            (O(x^-1), O(x^3))
            sage: t1, t2 = ET1(x^2, 5), ET2(x^3, 2/7); t1, t2
            (5 * x^2, 2/7 * x^3)

        In order for the comparison to work, the terms have come from
        or coerce into the same parent. Concretely, comparing
        :class:`GenericTerm` to, for example, an :class:`OTerm`
        always yields ``False``::

            sage: g1 <= g2
            True
            sage: o1, g1
            (O(x^-1), Generic Term with growth x)
            sage: o1 <= g1
            False

        If the elements of the common parent do not possess
        coefficients, then only the growth is compared::

            sage: o1 <= o1
            True
            sage: o1 <= o2
            True
            sage: o1 <= t1 and t1 <= o2
            True

        Otherwise, for asymptotic exact terms and `L` terms, in case
        of equal growth, also the coefficient is compared::

            sage: t1 <= t2
            True
            sage: t1 <= ET1(x^2, 3)
            False
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

        - ``other`` -- an asymptotic term to be compared to ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: GT = atm.GenericTermMonoid(PG)
            sage: t1, t2 = GT(x^-2), GT(x^5); t1, t2
            (Generic Term with growth x^-2, Generic Term with growth x^5)
            sage: t1._le_(t2)
            True
            sage: t2._le_(t1)
            False
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

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: PG1.<x> = agg.GrowthGroupPower()
        sage: PG2.<y> = agg.GrowthGroupPower(base=QQ)
        sage: GT1, GT2 = atm.GenericTermMonoid(PG1), atm.GenericTermMonoid(PG2)
        sage: GT1
        Generic Term Monoid over Asymptotic Power Growth Group in x over Integer Ring
        sage: GT2
        Generic Term Monoid over Asymptotic Power Growth Group in y over Rational Field
    """

    # enable the category framework for elements
    Element = GenericTerm


    def __init__(self, base=None, category=None):
        r"""
        See :class:`GenericTermMonoid` for more information.
        
        EXAMPLES::
        
            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG1.<x> = agg.GrowthGroupPower()
            sage: GT1 = atm.GenericTermMonoid(PG1); GT1
            Generic Term Monoid over Asymptotic Power Growth Group in x over Integer Ring
            sage: GT1.base()
            Asymptotic Power Growth Group in x over Integer Ring
            sage: PG2.<y> = agg.GrowthGroupPower()
            sage: GT2 = atm.GenericTermMonoid(PG2); GT2
            Generic Term Monoid over Asymptotic Power Growth Group in y over Integer Ring
            sage: GT1 is GT2
            False
        """

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

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: P = atm.GenericTermMonoid(base=PG)
            sage: t1, t2 = P(x), P(x^2)
            sage: P.le(t1,t2)
            True
        """
        return x.growth <= y.growth


class OTerm(GenericTerm):
    r"""
    Class for an asymptotic term representing an `O` term with
    specified growth. For the mathematical properties of `O` terms
    see :wikipedia:`Big_O_Notation`.

    `O` terms may *absorb* terms of weaker or equal growth.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an element of the parent's ``base``, specifying
      the growth of the asymptotic term.

    OUTPUT:

    An asymptotic `O` term.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: PG.<x> = agg.GrowthGroupPower()
        sage: OT = atm.OTermMonoid(PG)
        sage: t1, t2, t3 = OT(x^-7), OT(x^5), OT(x^42)
        sage: t1, t2, t3
        (O(x^-7), O(x^5), O(x^42))
        sage: t1.can_absorb(t2)
        False
        sage: t2.can_absorb(t1)
        True
        sage: t2.absorb(t1)
        O(x^5)
        sage: t1 <= t2 and t2 <= t3
        True
        sage: t3 <= t1
        False
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
            sage: t3
            O(x^3)
        """
        return "O(%s)" % self.growth


    def _absorb_(self, other):
        r"""
        Absorption of :class:`OTerm` ``other`` by :class:`OTerm`
        ``self``. Provided that the absorption is possible, the `O`
        term ``self`` is the dominant `O` term, and thus the result
        is ``self``.

        INPUT:

        - ``other`` -- an asymptotic `O` term from the same parent
          as ``self``.

        OUTPUT:

        The asymptotic `O` term ``self``.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: OT = atm.OTermMonoid(base=PG)
            sage: ot1, ot2 = OT(x), OT(x^2)
            sage: ot1.absorb(ot1)
            O(x)
            sage: ot2.absorb(ot1)
            O(x^2)
            sage: ot1.absorb(ot2)
            Traceback (most recent call last):
            ...
            ArithmeticError: O(x) cannot absorb O(x^2)
        """
        return self


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
        sage: PG1.<x> = agg.GrowthGroupPower(base=ZZ)
        sage: PG2.<y> = agg.GrowthGroupPower(base=QQ)
        sage: OT1 = atm.OTermMonoid(PG1); OT1
        Asymptotic O Term Monoid over Asymptotic Power Growth Group in x over Integer Ring
        sage: OT2 = atm.OTermMonoid(PG2); OT2
        Asymptotic O Term Monoid over Asymptotic Power Growth Group in y over Rational Field
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
            sage: PG1.<x> = agg.GrowthGroupPower(base=ZZ); x1 = x
            sage: PG2.<x> = agg.GrowthGroupPower(base=QQ); x2 = x
            sage: OT1 = atm.OTermMonoid(PG1)
            sage: OT2 = atm.OTermMonoid(PG2)
            sage: ET = atm.ExactTermMonoid(PG1, ZZ)

        Now, the ``OTermMonoid`` whose base is over the Integer Ring
        has to coerce into the ``OTermMonoid`` with the base over the
        rational field, and the ``ExactTermMonoid`` also has to
        coerce in both the ``OTermMonoid``'s::

            sage: OT2._coerce_map_from_(OT1)
            True
            sage: OT2._coerce_map_from_(ET)
            True
            sage: ET._coerce_map_from_(OT1) is None
            True
        """
        if isinstance(S, (ExactTermMonoid, LTermGenericMonoid)):
            if self.base().coerce_map_from(S.base()) is not None:
                return True
        else:
            return super(OTermMonoid, self)._coerce_map_from_(S)


    def _repr_(self):
        r"""
        Represent ``self`` as "Asymptotic O term monoid over ``base``".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: atm.OTermMonoid(PG)._repr_()
            'Asymptotic O Term Monoid over Asymptotic Power Growth Group in x over Integer Ring'
        """
        return "Asymptotic O Term Monoid over %s" % self.base()


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

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: PG.<x> = agg.GrowthGroupPower()
        sage: CT1 = atm.TermWithCoefficientMonoid(PG, ZZ)
        sage: CT2 = atm.TermWithCoefficientMonoid(PG, QQ)
        sage: CT1(x^2, 5)
        Asymptotic Term with coefficient 5 and growth x^2
        sage: CT2(x^3, 3/8)
        Asymptotic Term with coefficient 3/8 and growth x^3
    """

    def __init__(self, parent, growth, coefficient=1):
        r"""
        See :class:`TermWithCoefficient` for more information.

        EXAMPLES:

        First, we define some parent monoids::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: CT1 = atm.TermWithCoefficientMonoid(PG, ZZ)
            sage: CT2 = atm.TermWithCoefficientMonoid(PG, QQ)

        The coefficients have to be from the given coefficient ring::

            sage: t = CT1(x, 1/2)
            Traceback (most recent call last):
            ...
            ValueError: 1/2 is not in Integer Ring
            sage: t = CT2(x, 1/2); t
            Asymptotic Term with coefficient 1/2 and growth x

        For technical reasons, the coefficient 0 is not allowed::

            sage: t = CT1(x^42, 0)
            Traceback (most recent call last):
            ...
            ValueError: 0 is not a valid coefficient
        """
        if coefficient not in parent.coefficient_ring:
            raise ValueError("%s is not in %s" % (coefficient,
                                                  parent.coefficient_ring))
        elif coefficient == 0:
            raise ValueError("0 is not a valid coefficient")
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


    def _le_(self, other):
        r"""
        Return whether the asymptotic term with coefficient ``self``
        is of lesser or equal growth as the asymptotic term with
        coefficient ``other``.

        INPUT:

        - ``other`` -- an asymptotic term with coefficient to be
          compared to ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: ET = atm.ExactTermMonoid(PG, QQ)
            sage: t1, t2, t3 = ET(x, 5), ET(x^2, 3), ET(x^2, 42)
            sage: t1 <= t2
            True
            sage: t2 <= t1
            False
            sage: t2 <= t3
            True
            sage: t3 <= t2
            False
        """
        if self.growth == other.growth:
            return self.coefficient <= other.coefficient
        else:
            return super(TermWithCoefficient, self)._le_(other)


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

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: PG1 = agg.GrowthGroupPower("x"); x1 = PG1.gen()
        sage: PG2 = agg.GrowthGroupPower("x", base=QQ); x2 = PG2.gen()
        sage: TC1 = atm.TermWithCoefficientMonoid(PG1, QQ); TC1
        Monoid for asymptotic terms with coefficients from Rational Field over Asymptotic Power Growth Group in x over Integer Ring
        sage: TC2 = atm.TermWithCoefficientMonoid(PG2, QQ); TC2
        Monoid for asymptotic terms with coefficients from Rational Field over Asymptotic Power Growth Group in x over Rational Field
        sage: TC1 == TC2 or TC1 is TC2
        False
        sage: TC2.coerce_map_from(TC1)
        Conversion map:
          From: Monoid for asymptotic terms with coefficients from Rational Field over Asymptotic Power Growth Group in x over Integer Ring
          To:   Monoid for asymptotic terms with coefficients from Rational Field over Asymptotic Power Growth Group in x over Rational Field
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
    an explicitly specified constant and starting point.

    .. NOTE::

        When adding ("absorbing") various different `L` terms, some
        sort of "lifting" occurs in general. This lifting affects the
        coefficient of the resulting `L` term, and works
        differently for various growth classes. This makes different
        implementations for different growth classes necessary.

    .. TODO::

        Implement *absorption* (:meth:`_absorb_`) for `L` terms.

    INPUT:

        - ``parent`` -- the parent of the asymptotic term.

        - ``growth`` -- an asymptotic growth element from
          ``parent.base()``.

        - ``coefficient`` -- an element from
          ``parent.coefficient_ring``.

        - ``start`` -- a real number representing the starting point
          of the estimations in the definition of asymptotic `O`
          terms (see :wikipedia:`Big_O_Notation`).

    OUTPUT:

    An asymptotic `L` term with given coefficient and starting point.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: PG.<x> = agg.GrowthGroupPower()
        sage: LT1 = atm.LTermGenericMonoid(PG, ZZ)
        sage: LT2 = atm.LTermGenericMonoid(PG, QQ)
        sage: lt1, lt2 = LT1(x^2, 42, 42), LT2(x, 12/7, 0); lt1, lt2
        (42 * L(x^2, 42), 12/7 * L(x, 0))
        sage: lt1 <= lt2
        False
        sage: lt2 <= lt1
        True
        sage: lt1 * lt2
        72 * L(x^3, 42)
        sage: lt2.can_absorb(lt1)
        True
        sage: lt2.absorb(lt1)
        Traceback (most recent call last):
        ...
        NotImplementedError: Not yet implemented
    """


    def __init__(self, parent, growth, coefficient=1, start=0):
        r"""
        See :class:`LTermGeneric` for more information and examples.
        """
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


    def _absorb_(self, other):
        r"""
        Absorption of one `L` term ``other`` by another `L` term
        ``self``.

        INPUT:

        - ``other`` -- another asymptotic `L` term from the same
          parent as ``self``.

        OUTPUT:

        An asymptotic `L` term, representing the result of the
        absorption of ``other`` by ``self``.


        ..TODO::

            Implement this method for specific growth parents.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: LT = atm.LTermGenericMonoid(PG, QQ)
            sage: t1, t2 = LT(x^2, 5, 0), LT(x, 42, 5); t1, t2
            (5 * L(x^2, 0), 42 * L(x, 5))
            sage: t1.can_absorb(t2)
            True
            sage: t1.absorb(t2)
            Traceback (most recent call last):
            ...
            NotImplementedError: Not yet implemented
        """
        raise NotImplementedError("Not yet implemented")


    def _mul_(self, other):
        r"""
        Multiplication method for asymptotic `L` terms.

        INPUT:

        - ``other`` -- an asymptotic `L` term from the same parent
          as ``self``.

        OUTPUT:

        An asymptotic `L` term representing the product of ``self``
        and ``other``.

        .. NOTE::

            When taking the product of two asymptotic `L` terms,
            the growth of the product is the product of the growth
            elements of the factors, and analogue for the respective
            coefficient. The starting point of the resulting `L`
            term is the maximum starting point of the factors.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: LT = atm.LTermGenericMonoid(PG, QQ)
            sage: lt1, lt2 = LT(x^2, 2, 0), LT(x^-1, 3, 5)
            sage: lt1._mul_(lt1)
            4 * L(x^4, 0)
            sage: lt1._mul_(lt2)
            6 * L(x, 5)
        """
        result = super(LTermGeneric, self)._mul_(other)
        result.start = max(self.start, other.start)
        return result


class LTermGenericMonoid(TermWithCoefficientMonoid):
    r"""
    Base class for parents of asymptotic `L` terms. Also, the parent
    for :class:`LTermGeneric`. `L` terms are asymptotic terms which
    behave like `O` terms, with the difference that for `L` terms,
    the constant and the starting point for the inequality occurring
    in the definition of `O` terms (see :wikipedia:`Big_O_notation`)
    are explicitly given.

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

    A generic asymptotic `L` term monoid.

    EXAMPLES::

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: PG1 = agg.GrowthGroupPower("x"); x1 = PG1.gen()
        sage: PG2 = agg.GrowthGroupPower("x", base=QQ); x2 = PG2.gen()
        sage: LT1 = atm.LTermGenericMonoid(PG1, QQ); LT1
        Generic L Term Monoid with coefficients from Rational Field over Asymptotic Power Growth Group in x over Integer Ring
        sage: LT2 = atm.LTermGenericMonoid(PG2, QQ); LT2
        Generic L Term Monoid with coefficients from Rational Field over Asymptotic Power Growth Group in x over Rational Field
        sage: LT1 == LT2 or LT1 is LT2
        False
        sage: LT2.coerce_map_from(LT1)
        Conversion map:
          From: Generic L Term Monoid with coefficients from Rational Field over Asymptotic Power Growth Group in x over Integer Ring
          To:   Generic L Term Monoid with coefficients from Rational Field over Asymptotic Power Growth Group in x over Rational Field
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

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: LT = atm.LTermGenericMonoid(PG, QQ)
            sage: ET = atm.ExactTermMonoid(PG, QQ)
            sage: lt, et = LT(x^3, 42/5, 3), ET(x^7, 5/9)
            sage: LT(lt)
            42/5 * L(x^3, 3)
            sage: LT(lt) == lt
            True
            sage: LT(et)
            5/9 * L(x^7, 0)
        """
        if isinstance(x, LTermGeneric):
            return self.element_class(self, x.growth, x.coefficient, x.start)
        elif isinstance(x, ExactTerm):
            return self.element_class(self, x.growth, x.coefficient, 0)
        else:
            return self.element_class(self, x, coefficient, start)


    def _repr_(self):
        r"""
        Represent the generic `L` term monoid as ``Generic L Term
        Monoid with coefficients from coefficient_ring over base``.

        INPUT:

        Noting.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower(base=QQ)
            sage: atm.LTermGenericMonoid(PG, ZZ)._repr_()
            'Generic L Term Monoid with coefficients from Integer Ring over Asymptotic Power Growth Group in x over Rational Field'
        """
        return "Generic L Term Monoid with coefficients from %s over %s" % \
               (self.coefficient_ring, self.base())


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

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: LT1 = atm.LTermGenericMonoid(PG, QQ)
            sage: LT2 = atm.LTermGenericMonoid(PG, ZZ)
            sage: ET = atm.ExactTermMonoid(PG, ZZ)
            sage: LT1.coerce_map_from(LT2)
            Conversion map:
              From: Generic L Term Monoid with coefficients from Integer Ring over Asymptotic Power Growth Group in x over Integer Ring
              To:   Generic L Term Monoid with coefficients from Rational Field over Asymptotic Power Growth Group in x over Integer Ring
            sage: LT1.coerce_map_from(ET)
            Conversion map:
              From: Exact Term Monoid with coefficients from Integer Ring over Asymptotic Power Growth Group in x over Integer Ring
              To:   Generic L Term Monoid with coefficients from Rational Field over Asymptotic Power Growth Group in x over Integer Ring
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

        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: PG.<x> = agg.GrowthGroupPower()
        sage: ET = atm.ExactTermMonoid(PG, QQ)

    Asymptotic exact terms may be multiplied (with the usual rules
    applying)::

        sage: ET(x^2, 3) * ET(x, -1)
        -3 * x^3
        sage: ET(x^0, 4) * ET(x^5, 2)
        8 * x^5

    They may also be multiplied with `L` or `O` terms::

        sage: OT = atm.OTermMonoid(PG)
        sage: LT = atm.LTermGenericMonoid(PG, QQ)
        sage: ET(x^2, 42) * OT(x)
        O(x^3)
        sage: ET(x^2, 42) * LT(x, 1, 5)
        42 * L(x^3, 5)

    Absorption for asymptotic exact terms relates to addition::

        sage: ET(x^2, 5).can_absorb(ET(x^5, 12))
        False
        sage: ET(x^2, 5).can_absorb(ET(x^2, 1))
        True
        sage: ET(x^2, 5).absorb(ET(x^2, 1))
        6 * x^2

    Note that, as for technical reasons, `0` is not allowed as a
    coefficient for an asymptotic term with coefficient, ``None``
    is returned if two asymptotic exact terms cancel out each other
    during absorption::

        sage: ET(x^2, 42).can_absorb(ET(x^2, -42))
        True
        sage: ET(x^2, 42).absorb(ET(x^2, -42)) is None
        True
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


    def _absorb_(self, other):
        r"""
        Absorption of ``other`` by ``self``, where ``self`` and
        ``other`` are asymptotic exact terms over the same parent.
        For exact terms, absorption translates to adding the
        respective coefficients. For technical reasons, we return
        ``None`` if the resulting coefficient is `0`.

        INPUT:

        - ``other`` -- an asymptotic exact term from the same parent
          as ``self``.

        OUTPUT:

        An asymptotic exact term or ``None``.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: ET = atm.ExactTermMonoid(PG, QQ)

        Asymptotic exact terms can absorb other asymptotic exact
        terms with the same growth::

            sage: et1, et2 = ET(x^2, 2), ET(x^2, -2)
            sage: et1.absorb(et1)
            4 * x^2
            sage: repr(et1.absorb(et2))
            'None'

        If the growth differs, an ``ArithmeticException`` is raised::

            sage: ET(x^5, 1).absorb(et1)
            Traceback (most recent call last):
            ...
            ArithmeticError: 1 * x^5 cannot absorb 2 * x^2
        """
        cls = self.__class__
        coef_new = self.coefficient + other.coefficient
        if coef_new == 0:
            return None
        else:
            return cls(self.parent(), self.growth, coef_new)


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


        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: PG1 = agg.GrowthGroupPower("x"); x1 = PG1.gen()
        sage: PG2 = agg.GrowthGroupPower("x", base=QQ); x2 = PG2.gen()
        sage: ET1 = atm.ExactTermMonoid(PG1, ZZ); ET1
        Exact Term Monoid with coefficients from Integer Ring over Asymptotic Power Growth Group in x over Integer Ring
        sage: ET2 = atm.ExactTermMonoid(PG2, QQ); ET2
        Exact Term Monoid with coefficients from Rational Field over Asymptotic Power Growth Group in x over Rational Field
        sage: ET2.coerce_map_from(ET1)
        Conversion map:
          From: Exact Term Monoid with coefficients from Integer Ring over Asymptotic Power Growth Group in x over Integer Ring
          To:   Exact Term Monoid with coefficients from Rational Field over Asymptotic Power Growth Group in x over Rational Field
    """
    # enable the category framework for elements
    Element = ExactTerm


    def _repr_(self):
        r"""
        Represent the asymptotic exact term monoid as ``Exact Term
        Monoid with coefficients from coefficient_ring over base``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: import sage.groups.asymptotic_growth_group as agg
            sage: PG.<x> = agg.GrowthGroupPower()
            sage: atm.ExactTermMonoid(PG, QQ)._repr_()
            'Exact Term Monoid with coefficients from Rational Field over Asymptotic Power Growth Group in x over Integer Ring'
        """
        return "Exact Term Monoid with coefficients from %s over %s" % \
               (self.coefficient_ring, self.base())