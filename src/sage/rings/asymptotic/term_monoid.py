r"""
Asymptotic Term Monoid

This module implements asymptotic term monoids. The elements of these
monoids are used behind the scenes when performing calculations in an
asymptotic ring (to be implemented).

The monoids build upon the (asymptotic) growth groups. While growth
elements only model the growth of a function as it tends towards
infinity (or tends towards another fixed point; see
:mod:`~sage.rings.asymptotic.growth_group` for more details), an
asymptotic term additionally specifies its "type" and performs the
actual arithmetic operations (multiplication and partial
addition/absorption of terms).

Besides an abstract base term :class:`GenericTerm`, this module
implements the following types of terms:

- :class:`OTerm` -- `O`-terms at infinity, see
  :wikipedia:`Big_O_notation`.
- :class:`TermWithCoefficient` -- abstract base class for
  asymptotic terms with coefficients.
- :class:`ExactTerm` -- this class represents a growth element
  multiplied with some non-zero coefficient from a base ring.

.. WARNING::

    As this code is experimental, a warning is thrown when a term
    monoid is created for the first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: G = agg.MonomialGrowthGroup(ZZ, 'x')
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        sage: T = atm.TermWithCoefficientMonoid(G, ZZ)
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.

.. _term_absorption:

Absorption of Asymptotic Terms
==============================

A characteristic property of asymptotic terms is that some terms are
able to "absorb" other terms. This is realized with the method
:meth:`~sage.rings.asymptotic.term_monoid.GenericTerm.absorb`.

For instance, `O(x^2)` is able to absorb `O(x)` (with result
`O(x^2)`). This is because the functions bounded by linear growth
are bounded by quadratic growth as well. Another example would be
that `3x^5` is able to absorb `-2x^5` (with result `x^5`), which
simply corresponds to addition.

Essentially, absorption can be interpreted
as the addition of "compatible" terms (partial addition).

We want to show step by step which terms can be absorbed
by which other terms. We start by defining the necessary
term monoids and some terms::

    sage: import sage.rings.asymptotic.term_monoid as atm
    sage: import sage.rings.asymptotic.growth_group as agg
    sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
    sage: OT = atm.OTermMonoid(growth_group=G)
    sage: ET = atm.ExactTermMonoid(growth_group=G, base_ring=QQ)
    sage: ot1 = OT(x); ot2 = OT(x^2)
    sage: et1 = ET(x^2, 2)

- Because of the definition of `O`-terms (see
  :wikipedia:`Big_O_notation`), :class:`OTerm` are able to absorb all
  other asymptotic terms with weaker or equal growth. In our
  implementation, this means that :class:`OTerm` is able to absorb
  other :class:`OTerm`, as well as :class:`ExactTerm`, as long as the
  growth of the other term is less than or equal to the growth of this
  element::

      sage: ot1, ot2
      (O(x), O(x^2))
      sage: ot1.can_absorb(ot2), ot2.can_absorb(ot1)
      (False, True)
      sage: et1
      2*x^2
      sage: ot1.can_absorb(et1)
      False
      sage: ot2.can_absorb(et1)
      True

  The result of this absorption always is the dominant
  (absorbing) :class:`OTerm`::

      sage: ot1.absorb(ot1)
      O(x)
      sage: ot2.absorb(ot1)
      O(x^2)
      sage: ot2.absorb(et1)
      O(x^2)

  These examples correspond to `O(x) + O(x) = O(x)`,
  `O(x^2) + O(x) = O(x^2)`, and `O(x^2) + 2x^2 = O(x^2)`.

- :class:`ExactTerm` can only absorb another
  :class:`ExactTerm` if the growth coincides with the
  growth of this element::

      sage: et1.can_absorb(ET(x^2, 5))
      True
      sage: any(et1.can_absorb(t) for t in [ot1, ot2])
      False

  As mentioned above, absorption directly corresponds
  to addition in this case::

      sage: et1.absorb(ET(x^2, 5))
      7*x^2

  When adding two exact terms, they might cancel out.
  For technical reasons, ``None`` is returned in this
  case::

      sage: ET(x^2, 5).can_absorb(ET(x^2, -5))
      True
      sage: ET(x^2, 5).absorb(ET(x^2, -5)) is None
      True

- The abstract base terms :class:`GenericTerm` and
  :class:`TermWithCoefficient` can neither absorb any
  other term, nor be absorbed by any other term.

If ``absorb`` is called on a term that cannot be absorbed, an
:python:`ArithmeticError<library/exceptions.html#exceptions.ArithmeticError>`
is raised::

    sage: ot1.absorb(ot2)
    Traceback (most recent call last):
    ...
    ArithmeticError: O(x) cannot absorb O(x^2)

This would only work the other way around::

    sage: ot2.absorb(ot1)
    O(x^2)

Comparison
==========

The comparison of asymptotic terms with `\leq` is implemented as follows:

- When comparing ``t1 <= t2``, the coercion framework first tries to
  find a common parent for both terms. If this fails, ``False`` is
  returned.

- In case the coerced terms do not have a coefficient in their
  common parent (e.g. :class:`OTerm`), the growth of the two terms
  is compared.

- Otherwise, if the coerced terms have a coefficient (e.g.
  :class:`ExactTerm`), we compare whether ``t1`` has a growth that is
  strictly weaker than the growth of ``t2``. If so, we return
  ``True``. If the terms have equal growth, then we return ``True``
  if and only if the coefficients coincide as well.

  In all other cases, we return ``False``.

Long story short: we consider terms with different coefficients that
have equal growth to be incomparable.

Various
=======

.. TODO::

    - Implementation of more term types (e.g. `L` terms,
      `\Omega` terms, `o` terms, `\Theta` terms).

AUTHORS:

- Benjamin Hackl (2015-01): initial version
- Benjamin Hackl, Daniel Krenn (2015-05): conception of the asymptotic ring
- Benjamin Hackl (2015-06): refactoring caused by refactoring growth groups
- Daniel Krenn (2015-07): extensive review and patches
- Benjamin Hackl (2015-07): cross-review; short notation

Classes and Methods
===================
"""

# *****************************************************************************
# Copyright (C) 2014--2015 Benjamin Hackl <benjamin.hackl@aau.at>
#               2014--2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************

import sage

def absorption(left, right):
    r"""
    Let one of the two passed terms absorb the other.

    Helper function used by
    :class:`~sage.rings.asymptotic_ring.AsymptoticExpression`.

    .. NOTE::

        If neither of the terms can absorb the other, an
        :python:`ArithmeticError<library/exceptions.html#exceptions.ArithmeticError>`
        is raised.

        See the :ref:`module description <term_absorption>` for a
        detailed explanation of absorption.

    INPUT:

    - ``left`` -- an asymptotic term.

    - ``right`` -- an asymptotic term.

    OUTPUT:

    An asymptotic term or ``None``.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: G = agg.GrowthGroup('x^ZZ')
        sage: T = atm.TermMonoid('O', G)
        sage: atm.absorption(T(x^2), T(x^3))
        O(x^3)
        sage: atm.absorption(T(x^3), T(x^2))
        O(x^3)

    ::

        sage: T = atm.TermMonoid('exact', G, ZZ)
        sage: atm.absorption(T(x^2), T(x^3))
        Traceback (most recent call last):
        ...
        ArithmeticError: Absorption between x^2 and x^3 is not possible.
    """
    try:
        return left.absorb(right)
    except ArithmeticError:
        try:
            return right.absorb(left)
        except ArithmeticError:
            raise ArithmeticError('Absorption between %s and %s is not possible.' % (left, right))


def can_absorb(left, right):
    r"""
    Return whether one of the two input terms is able to absorb the
    other.

    Helper function used by
    :class:`~sage.rings.asymptotic_ring.AsymptoticExpression`.

    INPUT:

    - ``left`` -- an asymptotic term.

    - ``right`` -- an asymptotic term.

    OUTPUT:

    A boolean.

    .. NOTE::

        See the :ref:`module description <term_absorption>` for a
        detailed explanation of absorption.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: G = agg.GrowthGroup('x^ZZ')
        sage: T = atm.TermMonoid('O', G)
        sage: atm.can_absorb(T(x^2), T(x^3))
        True
        sage: atm.can_absorb(T(x^3), T(x^2))
        True
    """
    return left.can_absorb(right) or right.can_absorb(left)


class GenericTerm(sage.structure.element.MonoidElement):
    r"""
    Base class for asymptotic terms. Mainly the structure and
    several properties of asymptotic terms are handled here.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element.

    EXAMPLES::

        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
        sage: T = atm.GenericTermMonoid(G)
        sage: t1 = T(x); t2 = T(x^2); (t1, t2)
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

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.GenericTermMonoid(G)
            sage: T(x^2)
            Generic Term with growth x^2

        ::

            sage: atm.GenericTerm(parent=None, growth=x)
            Traceback (most recent call last):
            ...
            ValueError: The parent must be provided
            sage: atm.GenericTerm(T, agg.GrowthGroup('y^ZZ').gen())
            Traceback (most recent call last):
            ...
            ValueError: y is not in the parent's specified growth group
        """
        from sage.rings.asymptotic.growth_group import GenericGrowthElement

        if parent is None:
            raise ValueError('The parent must be provided')
        if growth is None or not isinstance(growth, GenericGrowthElement):
            raise ValueError('The growth must be provided and must inherit '
                             'from GenericGrowthElement')
        elif growth not in parent.growth_group:
            raise ValueError("%s is not in the parent's "
                             "specified growth group" % growth)

        self.growth = growth
        super(GenericTerm, self).__init__(parent=parent)


    def _mul_(self, other):
        r"""
        Multiplication of this term by another.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A :class:`GenericTerm`.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element, as well as ``other``
            are from a common parent.

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.GenericTermMonoid(G)
            sage: t1 = T(x); t2 = T(x^2)
            sage: t1, t2
            (Generic Term with growth x, Generic Term with growth x^2)
            sage: t1 * t2
            Generic Term with growth x^3
        """
        return self.parent()(self.growth * other.growth)


    def __div__(self, other):
        r"""
        Division of this term by another.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A :class:`GenericTerm`.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

            The comparison of two elements with the same parent is done in
            :meth:`_div_`.

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.GenericTermMonoid(G)
            sage: t1 = T(x); t2 = T(x^2)
            sage: t1 / t2  # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Inversion of Generic Term with growth x^2
            not implemented (in this abstract method).
        """
        from sage.structure.element import have_same_parent
        if have_same_parent(self, other):
            return self._div_(other)

        from sage.structure.element import get_coercion_model
        import operator
        return get_coercion_model().bin_op(self, other, operator.div)


    def _div_(self, other):
        r"""
        Division of this term by another.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A :class:`GenericTerm`.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element, as well as ``other``
            are from a common parent.

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.GenericTermMonoid(G)
            sage: t1 = T(x); t2 = T(x^2)
            sage: t1 / t2  # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Inversion of Generic Term with growth x^2
            not implemented (in this abstract method).
        """
        return self * ~other


    def __invert__(self):
        r"""
        Invert this term.

        OUTPUT:

        A :class:`GenericTerm`.

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.GenericTermMonoid(G)
            sage: ~T(x) # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Inversion of Generic Term with growth x
            not implemented (in this abstract method).
        """
        raise NotImplementedError('Inversion of %s not implemented '
                                  '(in this abstract method).' % (self,))


    def can_absorb(self, other):
        r"""
        Check whether this asymptotic term is able to absorb
        the asymptotic term ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            A :class:`GenericTerm` cannot absorb any other term.

            See the :ref:`module description <term_absorption>` for a
            detailed explanation of absorption.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: T = atm.GenericTermMonoid(G)
            sage: g1 = G(raw_element=21); g2 = G(raw_element=42)
            sage: t1 = T(g1); t2 = T(g2)
            sage: t1.can_absorb(t2)  # indirect doctest
            False
            sage: t2.can_absorb(t1)  # indirect doctest
            False
        """
        return False


    def absorb(self, other, check=True):
        r"""
        Absorb the asymptotic term ``other`` and return the resulting
        asymptotic term.

        INPUT:

        - ``other`` -- an asymptotic term.

        - ``check`` -- a boolean. If ``check`` is ``True`` (default),
          then ``can_absorb`` is called before absorption.

        OUTPUT:

        An asymptotic term or ``None`` if a cancellation occurs. If no
        absorption can be performed, an :python:`ArithmeticError<library/exceptions.html#exceptions.ArithmeticError>`
        is raised.

        .. NOTE::

            Setting ``check`` to ``False`` is meant to be used in
            cases where the respective comparison is done externally
            (in order to avoid duplicate checking).

            For a more detailed explanation of the *absorption* of
            asymptotic terms see
            the :ref:`module description <term_absorption>`.

        EXAMPLES:

        We want to demonstrate in which cases an asymptotic term
        is able to absorb another term, as well as explain the output
        of this operation. We start by defining some parents and
        elements::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_QQ = agg.GrowthGroup('x^QQ'); x = G_QQ.gen()
            sage: OT = atm.OTermMonoid(G_QQ)
            sage: ET = atm.ExactTermMonoid(growth_group=G_QQ, base_ring=QQ)
            sage: ot1 = OT(x); ot2 = OT(x^2)
            sage: et1 = ET(x, 100); et2 = ET(x^2, 2)
            sage: et3 = ET(x^2, 1); et4 = ET(x^2, -2)

        `O`-Terms are able to absorb other `O`-terms and exact terms
        with weaker or equal growth. ::

            sage: ot1.absorb(ot1)
            O(x)
            sage: ot1.absorb(et1)
            O(x)
            sage: ot1.absorb(et1) is ot1
            True

        :class:`ExactTerm` is able to absorb another
        :class:`ExactTerm` if the terms have the same growth. In this
        case, *absorption* is nothing else than an addition of the
        respective coefficients::

            sage: et2.absorb(et3)
            3*x^2
            sage: et3.absorb(et2)
            3*x^2
            sage: et3.absorb(et4)
            -x^2

        Note that, for technical reasons, the coefficient `0` is not
        allowed, and thus ``None`` is returned if two exact terms
        cancel each other out::

            sage: et2.absorb(et4)
            sage: et4.absorb(et2) is None
            True

        TESTS:

        When disabling the ``check`` flag, absorb might produce
        wrong results::

            sage: et1.absorb(ot2, check=False)
            O(x)
        """
        from sage.structure.element import have_same_parent

        if check:
            if not self.can_absorb(other):
                raise ArithmeticError('%s cannot absorb %s' % (self, other))

        if have_same_parent(self, other):
            return self._absorb_(other)

        from sage.structure.element import get_coercion_model

        return get_coercion_model().bin_op(self, other,
                                           lambda left, right:
                                           left._absorb_(right))


    def _absorb_(self, other):
        r"""
        Let this element absorb ``other``.

        INPUT:

        - ``other`` -- an asymptotic term from the same parent as
          this element.

        OUTPUT:

        An asymptotic term or ``None``.

        .. NOTE::

            This is not implemented for abstract base classes. For
            concrete realizations see, for example, :meth:`OTerm._absorb_`
            or :meth:`ExactTerm._absorb_`.
            Override this in derived class.

        EXAMPLES:

        First, we define some asymptotic terms::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.GenericTermMonoid(G)
            sage: t1 = T(x); t2 = T(x^2)

        When it comes to absorption, note that the method
        :meth:`can_absorb` (which is called before absorption takes
        place) does not allow the absorption of generic terms. Thus,
        an :python:`ArithmeticError<library/exceptions.html#exceptions.ArithmeticError>`
        is raised::

            sage: t2.absorb(t1)
            Traceback (most recent call last):
            ...
            ArithmeticError: Generic Term with growth x^2 cannot absorb Generic Term with growth x

        TESTS::

            sage: t2._absorb_(t1)
            Traceback (most recent call last):
            ...
            NotImplementedError: Not implemented in abstract base classes
        """
        raise NotImplementedError('Not implemented in abstract base classes')


    def __le__(self, other):
        r"""
        Return whether the growth of this term is less than
        or equal to the growth of ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method **only** compares the growth of the input
            terms!

        EXAMPLES:

        First, we define some asymptotic terms (and their parents)::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: GT = atm.GenericTermMonoid(G)
            sage: OT = atm.OTermMonoid(G)
            sage: ET_ZZ = atm.ExactTermMonoid(G, ZZ)
            sage: ET_QQ = atm.ExactTermMonoid(G, QQ)
            sage: g1 = GT(x); g2 = GT(x^2); g1, g2
            (Generic Term with growth x, Generic Term with growth x^2)
            sage: o1 = OT(x^-1); o2 = OT(x^3); o1, o2
            (O(1/x), O(x^3))
            sage: t1 = ET_ZZ(x^2, 5); t2 = ET_QQ(x^3, 2/7); t1, t2
            (5*x^2, 2/7*x^3)

        In order for the comparison to work, the terms have come from
        or coerce into the same parent. In particular, comparing
        :class:`GenericTerm` to, for example, an :class:`OTerm`
        always yields ``False``::

            sage: g1 <= g2
            True
            sage: o1, g1
            (O(1/x), Generic Term with growth x)
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

        For terms with coefficient (like exact terms), comparison
        works similarly, with the sole exception that terms with
        equal growth are considered incomparable. Thus, `\leq`
        only holds if the coefficients are equal as well::

            sage: t1 <= t2
            True
            sage: ET_ZZ(x, -5) <= ET_ZZ(x, 42)
            False
            sage: ET_ZZ(x, 5) <= ET_ZZ(x, 5)
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
        Return whether this generic term grows at most (i.e. less than
        or equal) like ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element, as well as ``other``
            are from the same parent.

            Also, this method **only** compares the growth of the
            input terms!

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.GenericTermMonoid(G)
            sage: t1 = T(x^-2); t2 = T(x^5); t1, t2
            (Generic Term with growth x^(-2), Generic Term with growth x^5)
            sage: t1._le_(t2)
            True
            sage: t2._le_(t1)
            False
        """
        return self.growth <= other.growth


    def __eq__(self, other):
        r"""
        Return whether this asymptotic term is equal to ``other``.

        INPUT:

        - ``other`` -- an object.

        OUTPUT:

        A boolean.

        .. NOTE::

            This function uses the coercion model to find a common
            parent for the two operands.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import (GenericTermMonoid,
            ....:      ExactTermMonoid, OTermMonoid)
            sage: GT = GenericTermMonoid(GrowthGroup('x^ZZ'))
            sage: ET = ExactTermMonoid(GrowthGroup('x^ZZ'), ZZ)
            sage: OT = OTermMonoid(GrowthGroup('x^ZZ'))
            sage: g = GT.an_element(); e = ET.an_element(); o = OT.an_element()
            sage: g, e, o
            (Generic Term with growth x, x, O(x))
            sage: e == e^2  # indirect doctest
            False
            sage: e == ET(x,1)  # indirect doctest
            True
            sage: o == OT(x^2)  # indirect doctest
            False
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
        Return whether this asymptotic term is the same as ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method gets called by the coercion framework, so it
            can be assumed that this asymptotic term is from the
            same parent as ``other``.

            Only implemented in concrete realizations.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: T = GenericTermMonoid(GrowthGroup('x^ZZ'))
            sage: t = T.an_element()
            sage: t == t
            True

        ::

            sage: from sage.rings.asymptotic.term_monoid import OTermMonoid
            sage: OT = OTermMonoid(GrowthGroup('x^ZZ'))
            sage: t = OT.an_element(); t
            O(x)
            sage: t == OT(x)  # indirect doctest
            True
            sage: t == OT(x^2)  # indirect doctest
            False
        """
        return self.growth == other.growth


    def _repr_(self):
        r"""
        A representation string for this generic term.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.GenericTermMonoid(growth_group=G)
            sage: T(x)._repr_()
            'Generic Term with growth x'
            sage: T(x^7)._repr_()
            'Generic Term with growth x^7'
        """
        return 'Generic Term with growth ' + repr(self.growth)


class GenericTermMonoid(sage.structure.parent.Parent,
                        sage.structure.unique_representation.UniqueRepresentation):
    r"""
    Parent for generic asymptotic terms.

    INPUT:

    - ``growth_group`` -- a growth group (i.e. an instance of
      :class:`~sage.rings.asymptotic.growth_group.GenericGrowthGroup`).

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. It has to be a subcategory
      of ``Join of Category of Monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    In this class the base
    structure for asymptotic term monoids will be handled. These
    monoids are the parents of asymptotic terms (for example, see
    :class:`GenericTerm` or :class:`OTerm`). Basically, asymptotic
    terms consist of a ``growth`` (which is an asymptotic growth
    group element, for example
    :class:`~sage.rings.asymptotic.growth_group.MonomialGrowthElement`);
    additional structure and properties are added by the classes inherited
    from :class:`GenericTermMonoid`.

    EXAMPLES::

        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G_x = agg.GrowthGroup('x^ZZ'); x = G_x.gen()
        sage: G_y = agg.GrowthGroup('y^QQ'); y = G_y.gen()
        sage: T_x_ZZ = atm.GenericTermMonoid(G_x); T_y_QQ = atm.GenericTermMonoid(G_y)
        sage: T_x_ZZ
        Generic Term Monoid x^ZZ
        sage: T_y_QQ
        Generic Term Monoid y^QQ
    """

    # enable the category framework for elements
    Element = GenericTerm

    @sage.misc.superseded.experimental(trac_number=17601)
    def __init__(self, growth_group, category=None):
        r"""
        See :class:`GenericTermMonoid` for more information.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_x = agg.GrowthGroup('x^ZZ')
            sage: T_x = atm.GenericTermMonoid(G_x); T_x
            Generic Term Monoid x^ZZ
            sage: T_x.growth_group
            Growth Group x^ZZ
            sage: G_y = agg.GrowthGroup('y^QQ')
            sage: T_y = atm.GenericTermMonoid(G_y); T_y
            Generic Term Monoid y^QQ
            sage: T_x is T_y
            False

        ::

            sage: atm.GenericTermMonoid(None)
            Traceback (most recent call last):
            ...
            ValueError: Growth Group has to be specified
        """
        from sage.categories.monoids import Monoids
        from sage.categories.posets import Posets
        from sage.rings.asymptotic.growth_group import GenericGrowthGroup

        if category is None:
            category = Monoids() & Posets()
        else:
            if not isinstance(category, tuple):
                category = (category,)
            if not any(cat.is_subcategory(Monoids() & Posets()) for cat in
                       category):
                raise ValueError('%s is not a subcategory of %s'
                                 % (category, Monoids() & Posets()))
        if growth_group is None:
            raise ValueError('Growth Group has to be specified')
        else:
            if not isinstance(growth_group, GenericGrowthGroup):
                raise ValueError('%s does not inherit from %s'
                                 % (growth_group, GenericGrowthGroup()))
        self._growth_group_ = growth_group
        super(GenericTermMonoid, self).__init__(category=category)

    @property
    def growth_group(self):
        r"""
        The growth group underlying this term monoid.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.ExactTermMonoid(G, ZZ).growth_group
            Growth Group x^ZZ
        """
        return self._growth_group_


    def _repr_(self):
        r"""
        A representation string for this generic term monoid.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: atm.GenericTermMonoid(G)._repr_()
            'Generic Term Monoid Generic(ZZ)'
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.GenericTermMonoid(growth_group=G)._repr_()
            'Generic Term Monoid x^ZZ'
        """
        return 'Generic Term Monoid %s' % (self.growth_group._repr_short_(), )


    def _coerce_map_from_(self, S):
        r"""
        Return whether ``S`` coerces into this term monoid.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            Another generic term monoid ``S`` coerces into this term
            monoid if and only if the growth group of ``S`` coerces
            into the growth group of this term monoid.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_ZZ = agg.GrowthGroup('x^ZZ')
            sage: T_ZZ = atm.GenericTermMonoid(growth_group=G_ZZ); T_ZZ
            Generic Term Monoid x^ZZ
            sage: G_QQ = agg.GrowthGroup('x^QQ')
            sage: T_QQ = atm.GenericTermMonoid(growth_group=G_QQ); T_QQ
            Generic Term Monoid x^QQ
            sage: T_QQ.has_coerce_map_from(T_ZZ)  # indirect doctest
            True
        """
        if isinstance(S, self.__class__):
            if self.growth_group.has_coerce_map_from(S.growth_group):
                return True


    def _element_constructor_(self, data):
        r"""
        Convert the given object to this term monoid.

        INPUT:

        - ``data`` -- an object representing the element to be
          initialized.

        OUTPUT:

        An element of this term monoid.

        .. NOTE::

            The object ``data`` is either an asymptotic term that is
            to be coerced into this term monoid, or an asymptotic
            growth element that is used for creating an element
            of this term monoid.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_ZZ = agg.GrowthGroup('x^ZZ')
            sage: G_QQ = agg.GrowthGroup('x^QQ')
            sage: T_ZZ = atm.GenericTermMonoid(growth_group=G_ZZ)
            sage: T_QQ = atm.GenericTermMonoid(growth_group=G_QQ)
            sage: term1 = T_ZZ(G_ZZ.gen())
            sage: term2 = T_QQ(G_QQ.gen()^2)

        In order for two terms to be compared, a coercion into
        a common parent has to be found::

            sage: term1.parent()
            Generic Term Monoid x^ZZ
            sage: term2.parent()
            Generic Term Monoid x^QQ
            sage: term1 <= term2
            True

        In this case, this works because ``T_ZZ``, the parent of
        ``term1``, coerces into ``T_QQ``::

            sage: T_QQ.coerce(term1)
            Generic Term with growth x

        The conversion of growth elements also works for the creation
        of terms::

            sage: x = SR('x'); x.parent()
            Symbolic Ring
            sage: T_ZZ(x^42)
            Generic Term with growth x^42
            sage: x = PolynomialRing(ZZ, 'x').gen(); x.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: T_ZZ(x^10)
            Generic Term with growth x^10
            sage: T_ZZ(10 * x^2)
            Traceback (most recent call last):
            ...
            ValueError: Input is ambiguous: cannot convert 10*x^2 to a generic term.
        """
        if isinstance(data, self.element_class) and data.parent() == self:
            return data
        elif isinstance(data, GenericTerm):
            return self.element_class(self, data.growth)
        elif isinstance(data, int) and data == 0:
            raise ValueError('No input specified. Cannot continue.')
        else:
            try:
                data = self.growth_group(data)
                return self.element_class(self, data)
            except:
                raise ValueError('Input is ambiguous: cannot convert %s to a '
                                 'generic term.' % (data,))


    def _an_element_(self):
        r"""
        Return an element of this term monoid.

        INPUT:

        Nothing.

        OUTPUT:

        An element of this term monoid.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.OTermMonoid(G).an_element()  # indirect doctest
            O(x)
            sage: atm.GenericTermMonoid(G).an_element()  # indirect doctest
            Generic Term with growth x
        """
        return self(self.growth_group.an_element())


    def some_elements(self):
        r"""
        Return some elements of this term monoid.

        See :class:`TestSuite` for a typical use case.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: tuple(atm.OTermMonoid(G).some_elements())
            (O(1), O(x), O(1/x), O(x^2), O(x^(-2)), O(x^3), ...)
        """
        return iter(self(g) for g in self.growth_group.some_elements())


    def le(self, left, right):
        r"""
        Return whether the term ``left`` is at most (less than or equal
        to) the term ``right``.

        INPUT:

        - ``left`` -- an element.

        - ``right`` -- an element.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.GenericTermMonoid(growth_group=G)
            sage: t1 = T(x); t2 = T(x^2)
            sage: T.le(t1, t2)
            True
        """
        return self(left) <= self(right)


class OTerm(GenericTerm):
    r"""
    Class for an asymptotic term representing an `O`-term with
    specified growth. For the mathematical properties of `O`-terms
    see :wikipedia:`Big_O_Notation`.

    `O`-terms can *absorb* terms of weaker or equal growth.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- a growth element.

    EXAMPLES::

        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
        sage: OT = atm.OTermMonoid(G)
        sage: t1 = OT(x^-7); t2 = OT(x^5); t3 = OT(x^42)
        sage: t1, t2, t3
        (O(x^(-7)), O(x^5), O(x^42))
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

    The conversion of growth elements also works for the creation
    of `O`-terms::

        sage: x = SR('x'); x.parent()
        Symbolic Ring
        sage: OT(x^17)
        O(x^17)
    """

    def _repr_(self):
        r"""
        A representation string for this `O`-term.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: OT = atm.OTermMonoid(G)
            sage: t1 = OT(x); t2 = OT(x^2); t3 = OT(x^3)
            sage: t1._repr_(), t2._repr_()
            ('O(x)', 'O(x^2)')
            sage: t3
            O(x^3)
        """
        return 'O(%s)' % self.growth


    def __invert__(self):
        r"""
        Invert this term.

        OUTPUT:

        A :class:`ZeroDivisionError` since `O`-terms cannot be inverted.

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.OTermMonoid(G)
            sage: ~T(x) # indirect doctest
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot invert O(x).
        """
        raise ZeroDivisionError('Cannot invert %s.' % (self,))


    def can_absorb(self, other):
        r"""
        Check whether this `O`-term can absorb ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            An :class:`OTerm` can absorb any other asymptotic term
            with weaker or equal growth.

            See the :ref:`module description <term_absorption>` for a
            detailed explanation of absorption.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: OT = atm.TermMonoid('O', agg.GrowthGroup('x^ZZ'))
            sage: t1 = OT(x^21); t2 = OT(x^42)
            sage: t1.can_absorb(t2)
            False
            sage: t2.can_absorb(t1)
            True
        """
        return other <= self


    def _absorb_(self, other):
        r"""
        Let this `O`-term absorb another `O`-term ``other``.

        INPUT:

        - ``other`` -- an asymptotic `O`-term.

        OUTPUT:

        An asymptotic `O`-term.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element and ``other``
            have the same parent.

            Also, observe that the result of a "dominant" `O`-term
            absorbing another `O`-term always is the "dominant"
            `O`-term again.

            See the :ref:`module description <term_absorption>` for a
            detailed explanation on absorption.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: OT = atm.OTermMonoid(growth_group=G)
            sage: ot1 = OT(x); ot2 = OT(x^2)
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
    Parent for asymptotic big `O`-terms.

    INPUT:

    - ``growth_group`` -- a growth group.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. It has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    EXAMPLES::

        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G_x_ZZ = agg.GrowthGroup('x^ZZ')
        sage: G_y_QQ = agg.GrowthGroup('y^QQ')
        sage: OT_x_ZZ = atm.OTermMonoid(G_x_ZZ); OT_x_ZZ
        Asymptotic O-Term Monoid x^ZZ
        sage: OT_y_QQ = atm.OTermMonoid(G_y_QQ); OT_y_QQ
        Asymptotic O-Term Monoid y^QQ

    `O`-term monoids can also be created by using the
    :class:`term factory <TermMonoidFactory>`::

        sage: atm.TermMonoid('O', G_x_ZZ) is OT_x_ZZ
        True
        sage: atm.TermMonoid('O', agg.GrowthGroup('x^QQ'))
        Asymptotic O-Term Monoid x^QQ
    """
    # enable the category framework for elements
    Element = OTerm


    def _coerce_map_from_(self, S):
        r"""
        Return whether ``S`` coerces into this term monoid.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        ``True`` or ``None``.

        .. NOTE::

            Another term monoid ``S`` coerces into this term monoid
            if ``S`` is an instance of one of the following classes:

            - :class:`OTermMonoid`

            - :class:`ExactTermMonoid`

            Additionally, the growth group underlying ``S`` has to
            coerce into the growth group of this term monoid.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_ZZ = agg.GrowthGroup('x^ZZ'); x_ZZ = G_ZZ.gen()
            sage: G_QQ = agg.GrowthGroup('x^QQ'); x_QQ = G_QQ.gen()
            sage: OT_ZZ = atm.OTermMonoid(G_ZZ)
            sage: OT_QQ = atm.OTermMonoid(G_QQ)
            sage: ET = atm.ExactTermMonoid(G_ZZ, ZZ)

        Now, the :class:`OTermMonoid` whose growth group is over the
        integer ring has to coerce into the :class:`OTermMonoid` with
        the growth group over the rational field, and the
        :class:`ExactTermMonoid` also has to coerce in each of the
        given :class:`OTermMonoid`::

            sage: OT_QQ.has_coerce_map_from(OT_ZZ)  # indirect doctest
            True
            sage: OT_QQ.has_coerce_map_from(ET)  # indirect doctest
            True
            sage: ET.has_coerce_map_from(OT_ZZ)  # indirect doctest
            False
        """
        if isinstance(S, (ExactTermMonoid,)):
            if self.growth_group.has_coerce_map_from(S.growth_group):
                return True
        else:
            return super(OTermMonoid, self)._coerce_map_from_(S)


    def _repr_(self):
        r"""
        A representation string for this `O`-term monoid.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: atm.OTermMonoid(G)._repr_()
            'Asymptotic O-Term Monoid x^ZZ'
        """
        return 'Asymptotic O-Term Monoid %s' % (self.growth_group._repr_short_(),)


class TermWithCoefficient(GenericTerm):
    r"""
    Base class for asymptotic terms possessing a coefficient. For
    example, :class:`ExactTerm` directly inherits from this class.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element of
      the parent's growth group.

    - ``coefficient`` -- an element of the parent's base ring.

    EXAMPLES::

        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
        sage: CT_ZZ = atm.TermWithCoefficientMonoid(G, ZZ)
        sage: CT_QQ = atm.TermWithCoefficientMonoid(G, QQ)
        sage: CT_ZZ(x^2, 5)
        Asymptotic Term with coefficient 5 and growth x^2
        sage: CT_QQ(x^3, 3/8)
        Asymptotic Term with coefficient 3/8 and growth x^3
    """

    def __init__(self, parent, growth, coefficient):
        r"""
        See :class:`TermWithCoefficient` for more information.

        EXAMPLES:

        First, we define some monoids::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: CT_ZZ = atm.TermWithCoefficientMonoid(G, ZZ)
            sage: CT_QQ = atm.TermWithCoefficientMonoid(G, QQ)

        The coefficients have to be from the given base ring::

            sage: t = CT_ZZ(x, 1/2)
            Traceback (most recent call last):
            ...
            ValueError: 1/2 is not in Integer Ring
            sage: t = CT_QQ(x, 1/2); t
            Asymptotic Term with coefficient 1/2 and growth x

        For technical reasons, the coefficient 0 is not allowed::

            sage: t = CT_ZZ(x^42, 0)
            Traceback (most recent call last):
            ...
            ValueError: 0 is not a valid coefficient.

        The conversion of growth elements also works for the creation
        of terms with coefficient::

            sage: x = SR('x'); x.parent()
            Symbolic Ring
            sage: CT_ZZ(x^42, 42)
            Asymptotic Term with coefficient 42 and growth x^42
        """
        if coefficient not in parent.base_ring():
            raise ValueError('%s is not in %s' % (coefficient,
                                                  parent.base_ring()))
        elif coefficient == 0:
            raise ValueError('0 is not a valid coefficient')

        self.coefficient = parent.base_ring()(coefficient)
        super(TermWithCoefficient, self).__init__(parent=parent, growth=growth)


    def _repr_(self):
        r"""
        A representation string for this term with coefficient.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.TermWithCoefficientMonoid(G, ZZ)
            sage: T(x^2, 5)._repr_()
            'Asymptotic Term with coefficient 5 and growth x^2'
        """
        return 'Asymptotic Term with coefficient %s and growth %s' % \
               (self.coefficient, self.growth)


    def _mul_(self, other):
        r"""
        Multiplication method for asymptotic terms with coefficients.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        An asymptotic term representing the product of this element
        and ``other``.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element and ``other`` have
            the same parent.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: CT = atm.TermWithCoefficientMonoid(G, ZZ)
            sage: ET = atm.ExactTermMonoid(G, ZZ)

        This method handles the multiplication of abstract terms
        with coefficient (i.e. :class:`TermWithCoefficient`) and
        exact terms (i.e. :class:`ExactTerm`). First, an example
        for abstract terms::

            sage: t1 = CT(x^2, 2); t2 = CT(x^3, 3)
            sage: t1 * t2
            Asymptotic Term with coefficient 6 and growth x^5

        And now, an example for exact terms::

            sage: t1 = ET(x^2, 2); t2 = ET(x^3, 3)
            sage: t1 * t2
            6*x^5
        """
        return self.parent()(self.growth * other.growth,
                             self.coefficient * other.coefficient)


    def _le_(self, other):
        r"""
        Return whether this asymptotic term with coefficient grows
        at most (less than or equal) like ``other``.

        INPUT:

        - ``other`` -- an asymptotic term with coefficient.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method is called by the coercion framework, thus,
            it can be assumed that this element and ``other`` are
            from the same parent.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import TermMonoid
            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: ET = TermMonoid('exact', G, QQ)
            sage: t1 = ET(x, 5); t2 = ET(x^2, 3); t3 = ET(x^2, 42)
            sage: t1 <= t2
            True
            sage: t2 <= t1
            False
            sage: t2 <= t3
            False
            sage: t3 <= t2
            False
            sage: t2 <= t2
            True

        TESTS::

            sage: ET(x, -2) <= ET(x, 1)
            False
        """
        if self.growth == other.growth:
            return self.coefficient == other.coefficient
        else:
            return super(TermWithCoefficient, self)._le_(other)


    def _eq_(self, other):
        r"""
        Return whether this :class:`TermWithCoefficient` is the same as
        ``other``.

        INPUT:

        - ``other`` -- an :class:`TermWithCoefficient`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method gets called by the coercion model, so it can
            be assumed that this term and ``other`` come from the
            same parent.

        EXAMPLES::

            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: from sage.rings.asymptotic.term_monoid import TermWithCoefficientMonoid
            sage: T = TermWithCoefficientMonoid(GrowthGroup('x^ZZ'), ZZ)
            sage: t = T.an_element(); t
            Asymptotic Term with coefficient 1 and growth x
            sage: t == T(x, 1)
            True
            sage: t == T(x, 2)
            False
            sage: t == T(x^2, 1)
            False
        """
        return super(TermWithCoefficient, self)._eq_(other) and \
            self.coefficient == other.coefficient


class TermWithCoefficientMonoid(GenericTermMonoid):
    r"""
    This class implements the base structure for parents of
    asymptotic terms possessing a coefficient from some coefficient
    ring. In particular, this is also the parent for
    :class:`TermWithCoefficient`.

    INPUT:

    - ``growth_group`` -- a growth group.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. It has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    - ``base_ring`` -- the ring which contains the
      coefficients of the elements.

    EXAMPLES::

        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G_ZZ = agg.GrowthGroup('x^ZZ'); x_ZZ = G_ZZ.gen()
        sage: G_QQ = agg.GrowthGroup('x^QQ'); x_QQ = G_QQ.gen()
        sage: TC_ZZ = atm.TermWithCoefficientMonoid(G_ZZ, QQ); TC_ZZ
        Term Monoid x^ZZ with coefficients from Rational Field
        sage: TC_QQ = atm.TermWithCoefficientMonoid(G_QQ, QQ); TC_QQ
        Term Monoid x^QQ with coefficients from Rational Field
        sage: TC_ZZ == TC_QQ or TC_ZZ is TC_QQ
        False
        sage: TC_QQ.coerce_map_from(TC_ZZ)
        Conversion map:
          From: Term Monoid x^ZZ with coefficients from Rational Field
          To:   Term Monoid x^QQ with coefficients from Rational Field
    """

    # enable the category framework for elements
    Element = TermWithCoefficient

    @sage.misc.superseded.experimental(trac_number=17601)
    def __init__(self, growth_group, base_ring, category=None):
        r"""
        For more information see :class:`TermWithCoefficientMonoid`.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T_ZZ = atm.TermWithCoefficientMonoid(G, ZZ); T_ZZ
            Term Monoid x^ZZ with coefficients from Integer Ring
            sage: T_QQ = atm.TermWithCoefficientMonoid(G, QQ); T_QQ
            Term Monoid x^ZZ with coefficients from Rational Field
            sage: T_QQ.category()
            Join of Category of monoids and Category of posets

        TESTS::

            sage: T = atm.TermWithCoefficientMonoid(G, None)
            Traceback (most recent call last):
            ...
            ValueError: None is not a valid base ring.
        """
        from sage.categories.rings import Rings
        if base_ring not in Rings():
            raise ValueError('%s is not a valid base ring.' % (base_ring,))
        self._base_ring_ = base_ring
        super(TermWithCoefficientMonoid,
              self).__init__(growth_group=growth_group, category=category)

    def base_ring(self):
        r"""
        The base ring of this term monoid, i.e. the ring where
        the coefficients are from.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.ExactTermMonoid(G, ZZ).base_ring()
            Integer Ring
        """
        return self._base_ring_


    def _coerce_map_from_(self, S):
        r"""
        Return whether ``S`` coerces into this term monoid.

        INPUT:

        - ``S`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            Another term monoid ``S`` coerces into this
            :class:`TermWithCoefficientMonoid`
            if both, the base ring as well as the growth
            group underlying ``S`` coerce into the base ring and the
            growth group underlying this term monoid, respectively.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_ZZ = agg.GrowthGroup('x^ZZ')
            sage: G_QQ = agg.GrowthGroup('x^QQ')
            sage: TC_ZZ = atm.TermWithCoefficientMonoid(G_ZZ, ZZ); TC_ZZ
            Term Monoid x^ZZ with coefficients from Integer Ring
            sage: TC_QQ = atm.TermWithCoefficientMonoid(G_QQ, QQ); TC_QQ
            Term Monoid x^QQ with coefficients from Rational Field
            sage: TC_QQ.has_coerce_map_from(TC_ZZ)  # indirect doctest
            True
            sage: TC_ZZ.has_coerce_map_from(TC_QQ)  # indirect doctest
            False
        """
        if isinstance(S, TermWithCoefficientMonoid):
            return (super(TermWithCoefficientMonoid, self)._coerce_map_from_(S) and
                    self.base_ring().has_coerce_map_from(S.base_ring()))


    def _element_constructor_(self, data, coefficient=None):
        r"""
        Construct an asymptotic term with coefficient or convert
        the given object ``data`` to this term monoid.

        INPUT:

        - ``data`` -- a growth element or an object representing the
          element to be initialized.

        - ``coefficient`` -- an element of the base ring.

        OUTPUT:

        An asymptotic term.

        .. NOTE::

            The object ``data`` is either an asymptotic term with
            coefficient that is to be coerced into this term monoid,
            or an asymptotic growth element that is used together
            with ``coefficient`` in order to create an element of
            this term monoid.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: T = atm.TermWithCoefficientMonoid(G, ZZ)
            sage: x.parent() == SR
            True
            sage: t1 = T(x^2, 5); t1  # indirect doctest
            Asymptotic Term with coefficient 5 and growth x^2

        TESTS::

            sage: T(5 * x^5)
            Asymptotic Term with coefficient 5 and growth x^5
            sage: T(G.gen()^10)
            Traceback (most recent call last):
            ...
            ValueError: Coefficient is not specified. Cannot continue.
            sage: T(G.gen()^10, coefficient=10)
            Asymptotic Term with coefficient 10 and growth x^10
            sage: T(x^123)
            Asymptotic Term with coefficient 1 and growth x^123

        ::

            sage: T(x)
            Asymptotic Term with coefficient 1 and growth x

        ::

            sage: G_log = agg.GrowthGroup('log(x)^ZZ')
            sage: T_log = atm.TermWithCoefficientMonoid(G_log, ZZ)
            sage: T_log(log(x))
            Asymptotic Term with coefficient 1 and growth log(x)
        """
        if isinstance(data, self.element_class) and data.parent() == self:
            return data
        elif isinstance(data, TermWithCoefficient):
            return self.element_class(self, data.growth, data.coefficient)
        elif isinstance(data, int) and data == 0:
            raise ValueError('No input specified. Cannot continue.')

        try:
            if coefficient is not None:
                data = self.growth_group(data)
                return self.element_class(self, data, coefficient)
            else:
                P = data.parent()
                from sage.symbolic.ring import SR
                import operator
                from sage.symbolic.operators import mul_vararg
                if P is SR:
                    op = data.operator()
                    if op == mul_vararg:
                        data, coef_tmp = data.operands()
                        data = self.growth_group(data)
                    elif op in (operator.pow, None) or\
                            isinstance(op, sage.functions.log.Function_log):
                        coef_tmp = 1
                        data = self.growth_group(data)
                else:
                    coeffs = data.coefficients()
                    if type(coeffs) == list:
                        # (multivariate) polynomial ring
                        coef_tmp = coeffs[0]
                        data = self.growth_group(data / coef_tmp)
                    elif type(coeffs) == dict:
                        # power series ring
                        coef_tmp = coeffs.values()[0]
                        data = self.growth_group(data / coef_tmp)

                return self.element_class(self, data, coef_tmp)
        except (ValueError, AttributeError):
            if coefficient is None:
                raise ValueError('Coefficient is not specified. '
                                 'Cannot continue.')
            elif coefficient not in self.base_ring():
                raise ValueError('%s is not in %s'
                                 % (coefficient, self.base_ring()))
            elif coefficient == 0:
                raise ValueError('0 is not a valid coefficient.')
            raise ValueError('Input is ambiguous: cannot convert %s with '
                             'coefficient %s to a term with coefficient.'
                             % (data, coefficient))


    def _repr_(self):
        r"""
        A representation string for this monoid of terms with
        coefficient.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.TermWithCoefficientMonoid(G, ZZ)._repr_()
            'Term Monoid x^ZZ with coefficients from Integer Ring'
        """
        return 'Term Monoid %s with coefficients from ' \
               '%s' % (self.growth_group._repr_short_(), self.base_ring())


    def _an_element_(self):
        r"""
        Return an element of this term with coefficient monoid.

        INPUT:

        Nothing.

        OUTPUT:

        An element of this term monoid.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.TermWithCoefficientMonoid(G, ZZ).an_element()  # indirect doctest
            Asymptotic Term with coefficient 1 and growth x
            sage: atm.ExactTermMonoid(G, ZZ).an_element()  # indirect doctest
            x
            sage: atm.ExactTermMonoid(G, QQ).an_element()  # indirect doctest
            1/2*x
        """
        return self(self.growth_group.an_element(),
                    self.base_ring().an_element())


    def some_elements(self):
        r"""
        Return some elements of this term with coefficient monoid.

        See :class:`TestSuite` for a typical use case.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: from itertools import islice
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GrowthGroup('z^QQ')
            sage: T = atm.ExactTermMonoid(G, ZZ)
            sage: tuple(islice(T.some_elements(), 10))
            (z^(1/2), z^(-1/2), -z^(1/2), z^2, -z^(-1/2), 2*z^(1/2),
             z^(-2), -z^2, 2*z^(-1/2), -2*z^(1/2))
        """
        from sage.misc.mrange import cantor_product
        return iter(self(g, c) for g, c in cantor_product(
            self.growth_group.some_elements(),
            iter(c for c in self.base_ring().some_elements() if c != 0)))


class ExactTerm(TermWithCoefficient):
    r"""
    Class for asymptotic exact terms. These terms primarily consist of
    an asymptotic growth element as well as a coefficient specifying
    the growth of the asymptotic term.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element from
      ``parent.growth_group``.

    - ``coefficient`` -- an element from ``parent.base_ring()``.

    EXAMPLES::

        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
        sage: ET = atm.ExactTermMonoid(G, QQ)

    Asymptotic exact terms may be multiplied (with the usual rules
    applying)::

        sage: ET(x^2, 3) * ET(x, -1)
        -3*x^3
        sage: ET(x^0, 4) * ET(x^5, 2)
        8*x^5

    They may also be multiplied with `O`-terms::

        sage: OT = atm.OTermMonoid(G)
        sage: ET(x^2, 42) * OT(x)
        O(x^3)

    Absorption for asymptotic exact terms relates to addition::

        sage: ET(x^2, 5).can_absorb(ET(x^5, 12))
        False
        sage: ET(x^2, 5).can_absorb(ET(x^2, 1))
        True
        sage: ET(x^2, 5).absorb(ET(x^2, 1))
        6*x^2

    Note that, as for technical reasons, `0` is not allowed as a
    coefficient for an asymptotic term with coefficient. Instead ``None``
    is returned if two asymptotic exact terms cancel out each other
    during absorption::

        sage: ET(x^2, 42).can_absorb(ET(x^2, -42))
        True
        sage: ET(x^2, 42).absorb(ET(x^2, -42)) is None
        True

    Exact terms can also be created by converting monomials with
    coefficient from the symbolic ring, or a suitable polynomial
    or power series ring::

        sage: x = var('x'); x.parent()
        Symbolic Ring
        sage: ET(5*x^2)
        5*x^2
        sage: x = ZZ['x'].gen(); x.parent()
        Univariate Polynomial Ring in x over Integer Ring
        sage: ET(5*x^2)
        5*x^2
        sage: x = ZZ[['x']].gen(); x.parent()
        Power Series Ring in x over Integer Ring
        sage: ET(5*x^2)
        5*x^2
    """

    def _repr_(self):
        r"""
        A representation string for this exact term.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: ET = atm.ExactTermMonoid(G, ZZ)
            sage: et1 = ET(x^2, 2); et1
            2*x^2

        TESTS::

            sage: ET(x^2, 1)
            x^2
            sage: ET(x^2, -1)
            -x^2
            sage: ET(x^0, 42)
            42
        """
        if self.growth._raw_element_ == 0:
            return '%s' % self.coefficient
        else:
            if self.coefficient == 1:
                return '%s' % self.growth
            elif self.coefficient == -1:
                return '-%s' % self.growth
            else:
                return '%s*%s' % (self.coefficient, self.growth)


    def __invert__(self):
        r"""
        Invert this term.

        OUTPUT:

        A term.

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.ExactTermMonoid(G, QQ)
            sage: ~T(x, 1/2)  # indirect doctest
            2*1/x
        """
        try:
            c = ~self.coefficient
        except ZeroDivisionError:
            raise ZeroDivisionError('Cannot invert %s since its coefficient %s '
                                    'cannot be inverted.' % (self, self.coefficient))
        return self.parent()(~self.growth, c)


    def can_absorb(self, other):
        r"""
        Check whether this exact term can absorb ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            For :class:`ExactTerm`, absorption corresponds to
            addition. This means that an exact term can absorb
            only other exact terms with the same growth.

            See the :ref:`module description <term_absorption>` for a
            detailed explanation of absorption.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: ET = atm.TermMonoid('exact', agg.GrowthGroup('x^ZZ'), ZZ)
            sage: t1 = ET(x^21, 1); t2 = ET(x^21, 2); t3 = ET(x^42, 1)
            sage: t1.can_absorb(t2)
            True
            sage: t2.can_absorb(t1)
            True
            sage: t1.can_absorb(t3) or t3.can_absorb(t1)
            False
        """
        return isinstance(other, ExactTerm) and self.growth == other.growth


    def _absorb_(self, other):
        r"""
        Let this exact term absorb another exact term ``other``.

        INPUT:

        - ``other`` -- an exact term.

        OUTPUT:

        An exact term or ``None``.

        .. NOTE::

            In the context of exact terms, absorption translates
            to addition. As the coefficient `0` is not allowed,
            ``None`` is returned instead if the terms cancel out.

            See the :ref:`module description <term_absorption>` for a
            detailed explanation on absorption.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: ET = atm.ExactTermMonoid(G, QQ)

        Asymptotic exact terms can absorb other asymptotic exact
        terms with the same growth::

            sage: et1 = ET(x^2, 2); et2 = ET(x^2, -2)
            sage: et1.absorb(et1)
            4*x^2
            sage: et1.absorb(et2) is None
            True

        If the growth differs, an ``ArithmeticException`` is raised::

            sage: ET(x^5, 1).absorb(et1)
            Traceback (most recent call last):
            ...
            ArithmeticError: x^5 cannot absorb 2*x^2
        """
        coeff_new = self.coefficient + other.coefficient
        if coeff_new == 0:
            return None
        else:
            return self.parent()(self.growth, coeff_new)


class ExactTermMonoid(TermWithCoefficientMonoid):
    r"""
    Parent for asymptotic exact terms, implemented in
    :class:`ExactTerm`.

    INPUT:

    - ``growth_group`` -- a growth group.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. It has to be a subcategory
      of ``Join of Category of monoids and Category of posets``. This
      is also the default category if ``None`` is specified.

    - ``base_ring`` -- the ring which contains the coefficients of
      the elements.

    EXAMPLES::

        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G_ZZ = agg.GrowthGroup('x^ZZ'); x_ZZ = G_ZZ.gen()
        sage: G_QQ = agg.GrowthGroup('x^QQ'); x_QQ = G_QQ.gen()
        sage: ET_ZZ = atm.ExactTermMonoid(G_ZZ, ZZ); ET_ZZ
        Exact Term Monoid x^ZZ with coefficients from Integer Ring
        sage: ET_QQ = atm.ExactTermMonoid(G_QQ, QQ); ET_QQ
        Exact Term Monoid x^QQ with coefficients from Rational Field
        sage: ET_QQ.coerce_map_from(ET_ZZ)
        Conversion map:
          From: Exact Term Monoid x^ZZ with coefficients from Integer Ring
          To:   Exact Term Monoid x^QQ with coefficients from Rational Field

    Exact term monoids can also be created using the
    :class:`term factory <TermMonoidFactory>`::

        sage: atm.TermMonoid('exact', G_ZZ, ZZ) is ET_ZZ
        True
        sage: atm.TermMonoid('exact', agg.GrowthGroup('x^ZZ'), QQ)
        Exact Term Monoid x^ZZ with coefficients from Rational Field
    """
    # enable the category framework for elements
    Element = ExactTerm

    def _repr_(self):
        r"""
        A representation string for this exact term monoid.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: atm.ExactTermMonoid(G, QQ)._repr_()
            'Exact Term Monoid x^ZZ with coefficients from Rational Field'
        """
        return 'Exact Term Monoid %s with coefficients from %s' % \
               (self.growth_group._repr_short_(), self.base_ring())


class TermMonoidFactory(sage.structure.factory.UniqueFactory):
    r"""
    Factory for asymptotic term monoids. It can generate the following
    term monoids:

    - :class:`OTermMonoid`,

    - :class:`ExactTermMonoid`.

    INPUT:

    - ``term`` -- the kind of term that shall be created. Either
      ``'exact'`` or ``'O'`` (capital letter ``O``).

    - ``growth_group`` -- a growth group.

    - ``base_ring`` -- the base ring for coefficients.

    OUTPUT:

    An asymptotic term monoid.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: G = agg.GrowthGroup('x^ZZ')
        sage: OT = atm.TermMonoid('O', G); OT
        Asymptotic O-Term Monoid x^ZZ
        sage: ET = atm.TermMonoid('exact', G, ZZ); ET
        Exact Term Monoid x^ZZ with coefficients from Integer Ring
    """
    def create_key_and_extra_args(self, term, growth_group, base_ring=None,
                                  **kwds):
        r"""
        Given the arguments and keyword, create a key that uniquely
        determines this object.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.TermMonoid.create_key_and_extra_args('O', G)
            (('O', Growth Group x^ZZ, None), {})
            sage: atm.TermMonoid.create_key_and_extra_args('exact', G, ZZ)
            (('exact', Growth Group x^ZZ, Integer Ring), {})
            sage: atm.TermMonoid.create_key_and_extra_args('exact', G)
            Traceback (most recent call last):
            ...
            ValueError: A base ring has to be specified

        TESTS::

            sage: atm.TermMonoid.create_key_and_extra_args('icecream', G)
            Traceback (most recent call last):
            ...
            ValueError: icecream has to be either 'exact' or 'O'
            sage: atm.TermMonoid.create_key_and_extra_args('O', ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Integer Ring has to be an asymptotic growth group
        """
        if term not in ['O', 'exact']:
            raise ValueError("%s has to be either 'exact' or 'O'" % term)

        from sage.rings.asymptotic.growth_group import GenericGrowthGroup
        if not isinstance(growth_group, GenericGrowthGroup):
            raise ValueError("%s has to be an asymptotic growth group"
                             % growth_group)

        if term == 'exact' and base_ring is None:
            raise ValueError("A base ring has to be specified")
        elif term == 'O':
            base_ring = None

        return (term, growth_group, base_ring), kwds


    def create_object(self, version, key, **kwds):
        r"""
        Create a object from the given arguments.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.TermMonoid('O', G)  # indirect doctest
            Asymptotic O-Term Monoid x^ZZ
            sage: atm.TermMonoid('exact', G, ZZ)  # indirect doctest
            Exact Term Monoid x^ZZ with coefficients from Integer Ring
        """

        term, growth_group, base_ring = key
        if term == 'O':
            return OTermMonoid(growth_group, **kwds)
        else:
            return ExactTermMonoid(growth_group, base_ring, **kwds)


TermMonoid = TermMonoidFactory("TermMonoid")
