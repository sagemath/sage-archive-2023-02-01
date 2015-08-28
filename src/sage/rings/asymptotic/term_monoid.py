r"""
Asymptotic Term Monoid

This module implements asymptotic term monoids. The elements of these
monoids are used behind the scenes when performing calculations in an
:mod:`asymptotic ring <sage.rings.asymptotic.asymptotic_ring>`.

The monoids build upon the (asymptotic) growth groups. While growth
elements only model the growth of a function as it tends towards
infinity (or tends towards another fixed point; see
:mod:`~sage.rings.asymptotic.growth_group` for more details), an
asymptotic term additionally specifies its "type" and performs the
actual arithmetic operations (multiplication and partial
addition/absorption of terms).

Besides an abstract base term :class:`GenericTerm`, this module
implements the following types of terms:

- :class:`OTerm` -- `O`-terms in infinity, see
  :wikipedia:`Big_O_notation`.
- :class:`TermWithCoefficient` -- abstract base class for
  asymptotic terms with coefficients.
- :class:`ExactTerm` -- this class represents a growth element
  multiplied with some non-zero coefficient from a coefficient ring.

A characteristic property of asymptotic terms is that some terms are
able to "absorb" other terms (see
:meth:`~sage.rings.asymptotic.term_monoid.GenericTerm.absorb`). For
instance, `O(x^2)` is able to absorb `O(x)` (with result
`O(x^2)`), and `3*x^5` is able to absorb `-2*x^5` (with result
`x^5`). Essentially, absorption can be interpreted as the
addition of "compatible" terms (partial addition).

.. TODO::

    - Implementation of more term types (e.g. `L` terms,
      `\Omega` terms, `o` terms, `\Theta` terms).

AUTHORS:

- Benjamin Hackl (2015-01): initial version
- Benjamin Hackl, Daniel Krenn (2015-05): conception of the asymptotic ring
- Benjamin Hackl (2015-06): refactoring caused by refactoring growth groups
- Daniel Krenn (2015-07): extensive review and patches
- Benjamin Hackl (2015-07): cross-review; short notation

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
        sage: T = atm.GenericTermMonoid(G, ZZ)
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


def product_diagonal(A, B):
    r"""
    Return an iterator over the product of `A` and `B` which iterates
    along the diagonal.

    INPUT:

    - ``A`` and ``B`` -- iterables (over a finite number of elements)

    OUTPUT:

    An iterator over `(a,b)` for `a \in A` and `b \in B`.

    EXAMPLES::

        sage: from sage.rings.asymptotic.term_monoid import product_diagonal
        sage: tuple(product_diagonal(srange(2), srange(2)))
        ((0, 0), (0, 1), (1, 0), (1, 1))
        sage: tuple(product_diagonal(srange(4), srange(2)))
        ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1))
        sage: tuple(product_diagonal(srange(2), srange(3)))
        ((0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (1, 2))
        sage: tuple(''.join(p) for p in product_diagonal('abc', 'xyz'))
        ('ax', 'ay', 'bx', 'az', 'by', 'cx', 'bz', 'cy', 'cz')

    TESTS:

    Check that all pairs are returned::

        sage: all(len(tuple(product_diagonal(srange(m), srange(n)))) == m*n
        ....:     for m in srange(5) for n in srange(5))
        True

    Check that everthing is loaded in the correct order::

        sage: def it(s, n):
        ....:     for i in srange(n):
        ....:         print '%s loads item number %s' % (s, i)
        ....:         yield i
        sage: for p in product_diagonal(it('A', 2), it('B', 2)):
        ....:     print p
        A loads item number 0
        B loads item number 0
        (0, 0)
        B loads item number 1
        (0, 1)
        A loads item number 1
        (1, 0)
        (1, 1)
        sage: for p in product_diagonal(it('A', 3), it('B', 2)):
        ....:     print p
        A loads item number 0
        B loads item number 0
        (0, 0)
        B loads item number 1
        (0, 1)
        A loads item number 1
        (1, 0)
        (1, 1)
        A loads item number 2
        (2, 0)
        (2, 1)
        sage: for p in product_diagonal(it('A', 2), it('B', 4)):
        ....:     print p
        A loads item number 0
        B loads item number 0
        (0, 0)
        B loads item number 1
        (0, 1)
        A loads item number 1
        (1, 0)
        B loads item number 2
        (0, 2)
        (1, 1)
        B loads item number 3
        (0, 3)
        (1, 2)
        (1, 3)
    """
    # when writing this code I thought the solution would be shorter...

    class iter_as_list(list):
        def __init__(self, iterable):
            self.it = iter(iterable)
            self.newdata = True
        def __getitem__(self, i):
            self.newdata = False
            try:
                while len(self) <= i:
                    self.append(next(self.it))
                    self.newdata = True
            except StopIteration:
                raise
            return list.__getitem__(self, i)

    from itertools import count
    A = iter_as_list(A)
    B = iter_as_list(B)
    for s in count():
        for i in range(s+1):
            stopped = False
            try:
                a = A[i]
            except StopIteration:
                stopped = True
            try:
                b = B[s-i]
            except StopIteration:
                stopped = True
            if stopped:
                continue
            yield a, b
        if not A.newdata and not B.newdata and s >= len(A) + len(B):
            return


def absorption(left, right):
    r"""
    Helper method used by
    :class:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticExpression`.

    INPUT:

    - ``left`` -- an asymptotic term.

    - ``right`` -- an asymptotic term.

    OUTPUT:

    An asymptotic term or ``None``.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: G = agg.GrowthGroup('x^ZZ')
        sage: T = atm.TermMonoid('O', G, ZZ)
        sage: atm.absorption(T(x^2), T(x^3))
        O(x^3)
        sage: atm.absorption(T(x^3), T(x^2))
        O(x^3)
    """
    try:
        return left.absorb(right)
    except ArithmeticError:
        return right.absorb(left)


def can_absorb(left, right):
    r"""
    Helper method used by
    :class:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticExpression`.

    INPUT:

    - ``left`` -- an asymptotic term.

    - ``right`` -- an asymptotic term.

    OUTPUT:

    A boolean.

    .. NOTE::

        This method returns whether one of the two input terms is
        able to absorb the other.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: G = agg.GrowthGroup('x^ZZ')
        sage: T = atm.TermMonoid('O', G, ZZ)
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
        sage: T = atm.GenericTermMonoid(G, QQ)
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
            sage: T = atm.GenericTermMonoid(G, ZZ)
            sage: T(x^2)
            Generic Term with growth x^2
        """
        if parent is None:
            raise ValueError('The parent must be provided')
        try:
            self.growth = parent.growth_group(growth)
        except ValueError, TypeError:
            raise ValueError("%s is not in %s" % (growth, parent.growth_group))

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
            sage: T = atm.GenericTermMonoid(G, ZZ)
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
            sage: T = atm.GenericTermMonoid(G, QQ)
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
            sage: T = atm.GenericTermMonoid(G, QQ)
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
            sage: T = atm.GenericTermMonoid(G, QQ)
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
        Check, whether this asymptotic term is able to absorb
        the asymptotic term ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        EXAMPLES:

        We want to show step by step which terms can be absorbed
        by which other terms. We start by defining the necessary
        term monoids and some terms::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: OT = atm.OTermMonoid(G, ZZ)
            sage: ET = atm.ExactTermMonoid(G, coefficient_ring=QQ)
            sage: ot1 = OT(x); ot2 = OT(x^2)
            sage: et1 = ET(x^2, 2)

        :class:`OTerm` is able to absorb other :class:`OTerm`,
        as well as :class:`ExactTerm`, as long as the growth of
        the other term is less than or equal to the growth of this
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

        :class:`ExactTerm` can only absorb another
        :class:`ExactTerm` if the growth coincides with the
        growth of this element::

            sage: et1.can_absorb(ET(x^2, 5))
            True
            sage: any(et1.can_absorb(t) for t in [ot1, ot2])
            False
        """
        return self._can_absorb_(other)


    def _can_absorb_(self, other):
        r"""
        Check, whether this generic term can absorb ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            Generic terms are not able to absorb any other term.
            Therefore, this method just returns ``False``.
            Override this in derived class.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GenericGrowthGroup(ZZ)
            sage: T = atm.GenericTermMonoid(G, QQ)
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
        Absorb the asymptotic term ``other`` and returns the resulting
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

            For a more detailed explanation of the *absorption* of
            asymptotic terms see the introduction of :mod:`this module
            <sage.rings.asymptotic.term_monoid>`, or the examples
            below.

        EXAMPLES:

        We want to demonstrate in which cases an asymptotic term
        is able to absorb another term, as well as explain the output
        of this operation. We start by defining some parents and
        elements::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_QQ = agg.GrowthGroup('x^QQ'); x = G_QQ.gen()
            sage: OT = atm.OTermMonoid(G_QQ, ZZ)
            sage: ET = atm.ExactTermMonoid(G_QQ, coefficient_ring=QQ)
            sage: ot1 = OT(x); ot2 = OT(x^2)
            sage: et1 = ET(x, 100); et2 = ET(x^2, 2)
            sage: et3 = ET(x^2, 1); et4 = ET(x^2, -2)

        Because of the definition of `O`-terms (see
        :wikipedia:`Big_O_notation`), they are able to absorb all
        other asymptotic terms with weaker or equal growth. The
        result of the absorption is the original `O`-term::

            sage: ot1.absorb(ot1)
            O(x)
            sage: ot1.absorb(et1)
            O(x)
            sage: ot1.absorb(et1) is ot1
            True

        The first example above corresponds to `O(x) + O(x) = O(x)`, and
        the second to `O(x) + 100x = O(x)`.
        If absorb is called on a term that cannot be absorbed, an
        :python:`ArithmeticError<library/exceptions.html#exceptions.ArithmeticError>`
        is raised::

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
                                           lambda self, other:
                                           self._absorb_(other))


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
            sage: T = atm.GenericTermMonoid(G, QQ)
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

        EXAMPLES:

        First, we define some asymptotic terms (and their parents)::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: GT = atm.GenericTermMonoid(G, QQ)
            sage: OT = atm.OTermMonoid(G, QQ)
            sage: ET_ZZ = atm.ExactTermMonoid(G, ZZ)
            sage: ET_QQ = atm.ExactTermMonoid(G, QQ)
            sage: g1 = GT(x); g2 = GT(x^2); g1, g2
            (Generic Term with growth x, Generic Term with growth x^2)
            sage: o1 = OT(x^-1); o2 = OT(x^3); o1, o2
            (O(x^(-1)), O(x^3))
            sage: t1 = ET_ZZ(x^2, 5); t2 = ET_QQ(x^3, 2/7); t1, t2
            (5*x^2, 2/7*x^3)

        In order for the comparison to work, the terms have come from
        or coerce into the same parent. Concretely, comparing
        :class:`GenericTerm` to, for example, an :class:`OTerm`
        always yields ``False``::

            sage: g1 <= g2
            True
            sage: o1, g1
            (O(x^(-1)), Generic Term with growth x)
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

        For terms with coefficient (like exact terms), the
        coefficient is compared as well if the growth is equal::

            sage: t1 <= t2
            True
            sage: t1 <= ET_ZZ(x^2, 3)
            False
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

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.GenericTermMonoid(G, QQ)
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
        Return if this asymptotic term is equal to ``other``.

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
            sage: GT = GenericTermMonoid(GrowthGroup('x^ZZ'), QQ)
            sage: ET = ExactTermMonoid(GrowthGroup('x^ZZ'), ZZ)
            sage: OT = OTermMonoid(GrowthGroup('x^ZZ'), QQ)
            sage: g = GT.an_element(); e = ET.an_element(); o = OT.an_element()
            sage: g, e, o
            (Generic Term with growth x, x, O(x))
            sage: g == g^2  # indirect doctest
            False
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
        Return if this asymptotic term is the same as ``other``.

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
            sage: T = GenericTermMonoid(GrowthGroup('x^ZZ'), QQ)
            sage: t = T.an_element()
            sage: t == t
            True

        ::

            sage: from sage.rings.asymptotic.term_monoid import OTermMonoid
            sage: OT = OTermMonoid(GrowthGroup('x^ZZ'), QQ)
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
            sage: T = atm.GenericTermMonoid(G, QQ)
            sage: T(x)._repr_()
            'Generic Term with growth x'
            sage: T(x^7)._repr_()
            'Generic Term with growth x^7'
        """
        return 'Generic Term with growth ' + repr(self.growth)


class GenericTermMonoid(sage.structure.unique_representation.UniqueRepresentation,
                        sage.structure.parent.Parent):
    r"""
    Parent for generic asymptotic terms.

    INPUT:

    - ``growth_group`` -- a partially ordered group (e.g. an instance of
      :class:`~sage.rings.asymptotic.growth_group.GenericGrowthGroup`).

    - ``coefficient_ring`` -- a ring which contains the (maybe implicit)
      coefficients of the elements.

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
        sage: T_x_ZZ = atm.GenericTermMonoid(G_x, QQ)
        sage: T_y_QQ = atm.GenericTermMonoid(G_y, QQ)
        sage: T_x_ZZ
        Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
        sage: T_y_QQ
        Generic Term Monoid y^QQ with (implicit) coefficients in Rational Field
    """

    # enable the category framework for elements
    Element = GenericTerm


    @staticmethod
    def __classcall__(cls, growth_group, coefficient_ring, category=None):
        r"""
        Normalizes the input in order to ensure a unique
        representation of the parent.

        TESTS::

            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G = GrowthGroup('x^ZZ')
            sage: T = GenericTermMonoid(G, QQ)
            sage: T.__class__(G, QQ) is T
            True
        """
        from sage.categories.sets_cat import Sets
        Sets_parent_class = Sets().parent_class
        while issubclass(cls, Sets_parent_class):
            cls = cls.__base__

        from sage.rings.asymptotic.growth_group import GenericGrowthGroup
        if growth_group is None:
            raise ValueError('No growth group specified.')
        if not isinstance(growth_group, sage.structure.parent.Parent):
            raise TypeError('%s is not a valid growth group.' % (growth_group,))

        if coefficient_ring is None:
            raise ValueError('No coefficient ring specified.')
        if not isinstance(coefficient_ring, sage.structure.parent.Parent):
            raise TypeError('%s is not a valid coefficient ring.' % (coefficient_ring,))

        if category is None:
            from sage.categories.monoids import Monoids
            from sage.categories.posets import Posets
            category = Monoids() & Posets()

        return super(GenericTermMonoid, cls).__classcall__(
            cls, growth_group, coefficient_ring, category)


    @sage.misc.superseded.experimental(trac_number=17601)
    def __init__(self, growth_group, coefficient_ring, category):
        r"""
        See :class:`GenericTermMonoid` for more information.

        EXAMPLES::

            sage: from sage.rings.asymptotic.term_monoid import GenericTermMonoid, TermWithCoefficientMonoid
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G_x = GrowthGroup('x^ZZ')
            sage: T_x = GenericTermMonoid(G_x, QQ); T_x
            Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
            sage: T_x.growth_group
            Growth Group x^ZZ
            sage: G_y = GrowthGroup('y^QQ')
            sage: T_y = GenericTermMonoid(G_y, QQ); T_y
            Generic Term Monoid y^QQ with (implicit) coefficients in Rational Field
            sage: T_x is T_y
            False

        ::

            sage: GenericTermMonoid(None, None)
            Traceback (most recent call last):
            ...
            ValueError: No growth group specified.

        ::

            sage: G = GrowthGroup('x^ZZ'); x = G.gen()
            sage: T_ZZ = TermWithCoefficientMonoid(G, ZZ); T_ZZ
            Generic Term Monoid x^ZZ with (implicit) coefficients in Integer Ring
            sage: T_QQ = TermWithCoefficientMonoid(G, QQ); T_QQ
            Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
            sage: T_QQ.category()
            Join of Category of monoids and Category of posets
        """
        self._growth_group_ = growth_group
        self._coefficient_ring_ = coefficient_ring
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


    @property
    def coefficient_ring(self):
        r"""
        The coefficient ring of this term monoid, i.e. the ring where
        the coefficients are from.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.GenericTermMonoid(G, ZZ).coefficient_ring
            Integer Ring
        """
        return self._coefficient_ring_


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
            sage: atm.GenericTermMonoid(G, QQ)._repr_()
            'Generic Term Monoid Generic(ZZ) with (implicit) coefficients in Rational Field'
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.GenericTermMonoid(G, QQ)._repr_()
            'Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field'
        """
        return 'Generic Term Monoid %s with (implicit) coefficients in %s' % \
            (self.growth_group._repr_short_(), self.coefficient_ring)


    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this term monoid.

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
            sage: G_ZZ = agg.GrowthGroup('x^ZZ')
            sage: T_ZZ = atm.GenericTermMonoid(G_ZZ, QQ); T_ZZ
            Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
            sage: G_QQ = agg.GrowthGroup('x^QQ')
            sage: T_QQ = atm.GenericTermMonoid(G_QQ, QQ); T_QQ
            Generic Term Monoid x^QQ with (implicit) coefficients in Rational Field
            sage: T_QQ.has_coerce_map_from(T_ZZ)  # indirect doctest
            True

        ::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G_ZZ = agg.GrowthGroup('x^ZZ')
            sage: G_QQ = agg.GrowthGroup('x^QQ')
            sage: TC_ZZ = atm.TermWithCoefficientMonoid(G_ZZ, ZZ); TC_ZZ
            Generic Term Monoid x^ZZ with (implicit) coefficients in Integer Ring
            sage: TC_QQ = atm.TermWithCoefficientMonoid(G_QQ, QQ); TC_QQ
            Generic Term Monoid x^QQ with (implicit) coefficients in Rational Field
            sage: TC_QQ.has_coerce_map_from(TC_ZZ)  # indirect doctest
            True
            sage: TC_ZZ.has_coerce_map_from(TC_QQ)  # indirect doctest
            False
        """
        if isinstance(S, self.__class__):
            if self.growth_group.has_coerce_map_from(S.growth_group) and \
                    self.coefficient_ring.has_coerce_map_from(S.coefficient_ring):
                return True


    def _element_constructor_(self, data, coefficient=None):
        r"""
        Convert the given object to this term monoid.

        INPUT:

        - ``data`` -- a growth element or an object representing the
          element to be initialized.

        - ``coefficient`` -- (default: ``None``)
          an element of the coefficient ring.

        OUTPUT:

        An element of this term monoid.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: G_QQ = GrowthGroup('x^QQ')
            sage: T_ZZ = atm.GenericTermMonoid(G_ZZ, QQ)
            sage: T_QQ = atm.GenericTermMonoid(G_QQ, QQ)
            sage: term1 = T_ZZ(G_ZZ.gen())
            sage: term2 = T_QQ(G_QQ.gen()^2)
            sage: term1 <= term2  # in order for two terms to be compared,
            ....:                 # a coercion into one of the parents
            ....:                 # has to be found.
            True
            sage: T_QQ.coerce(term1)  # coercion does not fail
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
            ValueError: 10*x^2 is not in Generic Term Monoid x^ZZ
            with (implicit) coefficients in Rational Field.
            > *previous* ValueError: Factor 10*x^2 of 10*x^2 is neither a
            coefficient (in Rational Field) nor growth (in Growth Group x^ZZ).

        ::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: T = atm.TermWithCoefficientMonoid(G, ZZ)
            sage: t1 = T(x^2, 5); t1  # indirect doctest
            Asymptotic Term with coefficient 5 and growth x^2

        TESTS::

            sage: O_ZZ = atm.OTermMonoid(G_ZZ, QQ)
            sage: O_ZZ(x^11)
            O(x^11)

        ::

            sage: T(G.gen()^10)
            Asymptotic Term with coefficient 1 and growth x^10
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
        if type(data) == self.element_class and data.parent() == self:
            return data
        elif isinstance(data, TermWithCoefficient):
            return self._create_element_(data.growth, data.coefficient)
        elif isinstance(data, GenericTerm):
            return self._create_element_(data.growth, None)
        elif type(data) == int and data == 0:
            raise ValueError('No input specified. Cannot continue '
                             'creating an element of %s.' % (self,))

        from growth_group import combine_exceptions
        if coefficient is not None:
            try:
                data = self.growth_group(data)
            except (ValueError, TypeError) as e:
                raise combine_exceptions(
                    ValueError('Growth %s is not in %s.' % (data, self)), e)
            return self._create_element_(data, coefficient)

        try:
            growth, coefficient = self._split_growth_and_coefficient_(data)
        except ValueError as e:
            raise combine_exceptions(
                ValueError('%s is not in %s.' % (data, self)), e)

        return self._create_element_(growth, coefficient)


    def _create_element_(self, growth, coefficient):
        r"""
        Helper method which creates an element by using the ``element_class``.

        INPUT:

        - ``growth`` -- a growth element.

        - ``coefficient`` -- an element of the coefficient ring.

        OUTPUT:

        A term.

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: T_ZZ = atm.GenericTermMonoid(G_ZZ, QQ)
            sage: T_ZZ(G_ZZ.gen())  # indirect doctest
            Generic Term with growth x
        """
        if coefficient is not None and coefficient != self.coefficient_ring.one():
            raise ValueError('Coefficient %s is not 1, but %s does not '
                             'support coefficients.' % (coefficient, self))
        return self.element_class(self, growth)


    def _split_growth_and_coefficient_(self, data):
        r"""
        Split given ``data`` into a growth element and a coefficient.

        INPUT:

        ``data`` -- an element.

        OUTPUT:

        A pair ``(growth, coefficient``).

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: T_ZZ = atm.ExactTermMonoid(G_ZZ, QQ)
            sage: T_ZZ._split_growth_and_coefficient_('2*x^3')
            (x^3, 2)
        """
        factors = self._get_factors_(data)

        growth_group = self.growth_group
        coefficient_ring = self.coefficient_ring

        coefficients = []
        growths = []
        for f in factors:
            try:
                growths.append(growth_group(f))
                continue
            except (ValueError, TypeError):
                pass

            try:
                coefficients.append(coefficient_ring(f))
                continue
            except (ValueError, TypeError):
                pass

            raise ValueError('Factor %s of %s is neither a coefficient (in %s) '
                             'nor growth (in %s).' %
                             (f, data, coefficient_ring, growth_group))

        from sage.misc.misc_c import prod
        growth = prod(growths) if growths else growth_group.one()
        coefficient = prod(coefficients) if coefficients else coefficient_ring.one()
        return (growth, coefficient)


    def _get_factors_(self, data):
        r"""
        Split given ``data`` into separate factors.

        INPUT:

        - ``data`` -- an object.

        OUTPUT:

        A tuple.

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: T_ZZ = atm.ExactTermMonoid(G_ZZ, QQ)
            sage: T_ZZ._get_factors_(x^2 * log(x))
            (x^2, log(x))
        """
        if isinstance(data, str):
            from growth_group import split_str_by_mul
            return split_str_by_mul(data)

        try:
            P = data.parent()
        except AttributeError:
            return (data,)

        from sage.symbolic.ring import SR
        if P is SR:
            from sage.symbolic.operators import mul_vararg
            if data.operator() == mul_vararg:
                return tuple(data.operands())

        return (data,)


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
            sage: atm.OTermMonoid(G, QQ).an_element()  # indirect doctest
            O(x)
            sage: atm.GenericTermMonoid(G, QQ).an_element()  # indirect doctest
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
            sage: tuple(atm.OTermMonoid(G, QQ).some_elements())
            (O(1), O(x), O(x^(-1)), O(x^2), O(x^(-2)), O(x^3), ...)
        """
        return iter(self(g) for g in self.growth_group.some_elements())


    def le(self, left, right):
        r"""
        Return if the term ``left`` is at most (less than or equal
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
            sage: T = atm.GenericTermMonoid(G, QQ)
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
        sage: OT = atm.OTermMonoid(G, QQ)
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
            sage: OT = atm.OTermMonoid(G, QQ)
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
            sage: T = atm.OTermMonoid(G, QQ)
            sage: ~T(x) # indirect doctest
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Cannot invert O(x).
        """
        raise ZeroDivisionError('Cannot invert %s.' % (self,))


    def _can_absorb_(self, other):
        r"""
        Check, whether this `O`-term can absorb ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            An :class:`OTerm` can absorb any other asymptotic term
            with weaker or equal growth.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: OT = atm.TermMonoid('O', agg.GrowthGroup('x^ZZ'), QQ)
            sage: t1 = OT(x^21); t2 = OT(x^42)
            sage: t1.can_absorb(t2)  # indirect doctest
            False
            sage: t2.can_absorb(t1)  # indirect doctest
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
            absorbing another `O`-term, always is the "dominant"
            `O`-term again.

        EXAMPLES::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: OT = atm.OTermMonoid(G, QQ)
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
        sage: OT_x_ZZ = atm.OTermMonoid(G_x_ZZ, QQ); OT_x_ZZ
        O-Term Monoid x^ZZ with implicit coefficients in Rational Field
        sage: OT_y_QQ = atm.OTermMonoid(G_y_QQ, QQ); OT_y_QQ
        O-Term Monoid y^QQ with implicit coefficients in Rational Field

    `O`-term monoids can also be created by using the
    :class:`term factory <TermMonoidFactory>`::

        sage: atm.TermMonoid('O', G_x_ZZ, QQ) is OT_x_ZZ
        True
        sage: atm.TermMonoid('O', agg.GrowthGroup('x^QQ'), QQ)
        O-Term Monoid x^QQ with implicit coefficients in Rational Field
    """

    # enable the category framework for elements
    Element = OTerm


    def _create_element_(self, growth, coefficient):
        r"""
        Helper method which creates an element by using the ``element_class``.

        INPUT:

        - ``growth`` -- a growth element.

        - ``coefficient`` -- an element of the coefficient ring (will be
          ignored since we create an O-Term).

        OUTPUT:

        An O-term.

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: T_ZZ = atm.OTermMonoid(G_ZZ, QQ)
            sage: T_ZZ(G_ZZ.gen())  # indirect doctest
            O(x)
        """
        return self.element_class(self, growth)


    def _coerce_map_from_(self, S):
        r"""
        Return if ``S`` coerces into this term monoid.

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
            sage: OT_ZZ = atm.OTermMonoid(G_ZZ, QQ)
            sage: OT_QQ = atm.OTermMonoid(G_QQ, QQ)
            sage: ET = atm.ExactTermMonoid(G_ZZ, ZZ)

        Now, the :class:`OTermMonoid` whose growth group is over the
        interger ring has to coerce into the :class:`OTermMonoid` with
        the growth group over the rational field, and the
        :class:`ExactTermMonoid` also has to coerce in both the
        :class:`OTermMonoid`s::

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
            sage: atm.OTermMonoid(G, QQ)._repr_()
            'O-Term Monoid x^ZZ with implicit coefficients in Rational Field'
        """
        return 'O-Term Monoid %s with implicit coefficients in %s' % \
            (self.growth_group._repr_short_(), self.coefficient_ring)


class TermWithCoefficient(GenericTerm):
    r"""
    Base class for asymptotic terms possessing a coefficient. For
    example, :class:`ExactTerm` directly inherits from this class.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element of
      the parent's growth group.

    - ``coefficient`` -- an element of the parent's coefficient ring.

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

        The coefficients have to be from the given coefficient ring::

            sage: t = CT_ZZ(x, 1/2)
            Traceback (most recent call last):
            ...
            ValueError: 1/2 is not a coefficient in
            Generic Term Monoid x^ZZ with (implicit) coefficients in Integer Ring.
            sage: t = CT_QQ(x, 1/2); t
            Asymptotic Term with coefficient 1/2 and growth x

        For technical reasons, the coefficient 0 is not allowed::

            sage: t = CT_ZZ(x^42, 0)
            Traceback (most recent call last):
            ...
            ValueError:  Zero coefficient 0 is not allowed in
            Generic Term Monoid x^ZZ with (implicit) coefficients in Integer Ring.

        The conversion of growth elements also works for the creation
        of terms with coefficient::

            sage: x = SR('x'); x.parent()
            Symbolic Ring
            sage: CT_ZZ(x^42, 42)
            Asymptotic Term with coefficient 42 and growth x^42
        """
        try:
            coefficient = parent.coefficient_ring(coefficient)
        except (ValueError, TypeError):
            raise ValueError('%s is not a coefficient in %s.' %
                             (coefficient, parent))
        if coefficient == 0:
            raise ValueError('Zero coefficient %s is not allowed in %s.' %
                             (coefficient, parent))

        self.coefficient = coefficient
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

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: ET = atm.ExactTermMonoid(G, QQ)
            sage: t1 = ET(x, 5); t2 = ET(x^2, 3); t3 = ET(x^2, 42)
            sage: t1 <= t2
            True
            sage: t2 <= t1
            False
            sage: t2 <= t3
            True
            sage: t3 <= t2
            False

        TESTS::

            sage: ET(x, -2) <= ET(x, 1)
            False
        """


        if self.growth == other.growth:
            return abs(self.coefficient) < abs(other.coefficient) or \
                self.coefficient == other.coefficient
        else:
            return super(TermWithCoefficient, self)._le_(other)


    def _eq_(self, other):
        r"""
        Return if this :class:`TermWithCoefficient` is the same as
        ``other``.

        INPUT:

        - ``other`` -- an :class:`TermWithCoefficient`.

        OUTPUT:

        A boolean.

        .. NOTE::

            This method gets called by the coercion model, so it can
            be assumed that this :class:`TermWithCoefficient` and
            ``other`` come from the same parent.

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

    - ``coefficient_ring`` -- the ring which contains the
      coefficients of the elements.

    EXAMPLES::

        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G_ZZ = agg.GrowthGroup('x^ZZ'); x_ZZ = G_ZZ.gen()
        sage: G_QQ = agg.GrowthGroup('x^QQ'); x_QQ = G_QQ.gen()
        sage: TC_ZZ = atm.TermWithCoefficientMonoid(G_ZZ, QQ); TC_ZZ
        Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
        sage: TC_QQ = atm.TermWithCoefficientMonoid(G_QQ, QQ); TC_QQ
        Generic Term Monoid x^QQ with (implicit) coefficients in Rational Field
        sage: TC_ZZ == TC_QQ or TC_ZZ is TC_QQ
        False
        sage: TC_QQ.coerce_map_from(TC_ZZ)
        Conversion map:
          From: Generic Term Monoid x^ZZ with (implicit) coefficients in Rational Field
          To:   Generic Term Monoid x^QQ with (implicit) coefficients in Rational Field
    """

    # enable the category framework for elements
    Element = TermWithCoefficient


    def _create_element_(self, growth, coefficient):
        r"""
        Helper method which creates an element by using the ``element_class``.

        INPUT:

        - ``growth`` -- a growth element.

        - ``coefficient`` -- an element of the coefficient ring.

        OUTPUT:

        A term.

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: from sage.rings.asymptotic.growth_group import GrowthGroup
            sage: G_ZZ = GrowthGroup('x^ZZ')
            sage: T_ZZ = atm.ExactTermMonoid(G_ZZ, QQ)
            sage: T_ZZ(G_ZZ.gen(), 4/3)  # indirect doctest
            4/3*x
        """
        return self.element_class(self, growth, coefficient)


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
                    self.coefficient_ring.an_element())


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
            (z^(1/2), -z^(1/2), z^(-1/2), 2*z^(1/2), -z^(-1/2),
             z^2, -2*z^(1/2), 2*z^(-1/2), -z^2, z^(-2))

        """
        return iter(self(g, c) for g, c in product_diagonal(
            self.growth_group.some_elements(),
            iter(c for c in self.coefficient_ring.some_elements() if c != 0)))


class ExactTerm(TermWithCoefficient):
    r"""
    Class for asymptotic exact terms. These terms primarily consist of
    an asymptotic growth element as well as a coefficient specifying
    the growth of the asymptotic term.

    INPUT:

    - ``parent`` -- the parent of the asymptotic term.

    - ``growth`` -- an asymptotic growth element from
      ``parent.growth_group``.

    - ``coefficient`` -- an element from ``parent.coefficient_ring``.

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

        sage: OT = atm.OTermMonoid(G, QQ)
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
        g = repr(self.growth)
        c = repr(self.coefficient)
        if g == '1':
            return c
        elif c == '1':
            return '%s' % (g,)
        elif c == '-1':
            return '-%s' % (g,)
        else:
            return '%s*%s' % (c, g)


    def __invert__(self):
        r"""
        Invert this term.

        OUTPUT:

        A term.

        TESTS::

            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: import sage.rings.asymptotic.growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ'); x = G.gen()
            sage: T = atm.ExactTermMonoid(G, ZZ)
            sage: ~T(x, 2)  # indirect doctest
            1/2*x^(-1)
            sage: (~T(x, 2)).parent()
            Exact Term Monoid x^ZZ with coefficients in Rational Field
        """
        try:
            c = ~self.coefficient
        except ZeroDivisionError:
            raise ZeroDivisionError('Cannot invert %s since its coefficient %s '
                                    'cannot be inverted.' % (self, self.coefficient))
        g = ~self.growth
        if c.parent() is self.coefficient.parent() and g.parent() is self.growth.parent():
            return self.parent()(g, c)
        else:
            new_parent = self.parent().__class__(g.parent(), c.parent())
            return new_parent(g, c)


    def _can_absorb_(self, other):
        r"""
        Check, whether this exact term can absorb ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        A boolean.

        .. NOTE::

            For :class:`ExactTerm`, absorption corresponds to
            addition. This means that an exact term can absorb
            only other exact terms with the same growth.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: ET = atm.TermMonoid('exact', agg.GrowthGroup('x^ZZ'), ZZ)
            sage: t1 = ET(x^21, 1); t2 = ET(x^21, 2); t3 = ET(x^42, 1)
            sage: t1.can_absorb(t2)  # indirect doctest
            True
            sage: t2.can_absorb(t1)  # indirect doctest
            True
            sage: t1.can_absorb(t3) or t3.can_absorb(t1) # indirect doctest
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
            to addition. As the coefficient `0` is not allowed. Instead
            ``None`` is returned if the terms cancel out.

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

    - ``coefficient_ring`` -- the ring which contains the coefficients of
      the elements.

    EXAMPLES::

        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: import sage.rings.asymptotic.growth_group as agg
        sage: G_ZZ = agg.GrowthGroup('x^ZZ'); x_ZZ = G_ZZ.gen()
        sage: G_QQ = agg.GrowthGroup('x^QQ'); x_QQ = G_QQ.gen()
        sage: ET_ZZ = atm.ExactTermMonoid(G_ZZ, ZZ); ET_ZZ
        Exact Term Monoid x^ZZ with coefficients in Integer Ring
        sage: ET_QQ = atm.ExactTermMonoid(G_QQ, QQ); ET_QQ
        Exact Term Monoid x^QQ with coefficients in Rational Field
        sage: ET_QQ.coerce_map_from(ET_ZZ)
        Conversion map:
          From: Exact Term Monoid x^ZZ with coefficients in Integer Ring
          To:   Exact Term Monoid x^QQ with coefficients in Rational Field

    Exact term monoids can also be created using the term factory::

        sage: atm.TermMonoid('exact', G_ZZ, ZZ) is ET_ZZ
        True
        sage: atm.TermMonoid('exact', agg.GrowthGroup('x^ZZ'), QQ)
        Exact Term Monoid x^ZZ with coefficients in Rational Field
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
            'Exact Term Monoid x^ZZ with coefficients in Rational Field'
        """
        return 'Exact Term Monoid %s with coefficients in %s' % \
               (self.growth_group._repr_short_(), self.coefficient_ring)


class TermMonoidFactory(sage.structure.factory.UniqueFactory):
    r"""
    Factory for asymptotic term monoids. It can generate the following
    term monoids:

    - :class:`OTermMonoid`,

    - :class:`ExactTermMonoid`.

    INPUT:

    - ``term`` -- the kind of term that shall be created. Either a string
      ``'exact'`` or ``'O'``, or an existing instance of a term.

    - ``growth_group`` -- a growth group.

    - ``coefficient_ring`` -- a ring.

    - ``asymptotic_ring`` -- if specified, then ``growth_group`` and
      ``coefficient_ring`` are taken from this asymptotic ring.

    OUTPUT:

    An asymptotic term monoid.

    EXAMPLES::

        sage: import sage.rings.asymptotic.growth_group as agg
        sage: import sage.rings.asymptotic.term_monoid as atm
        sage: G = agg.GrowthGroup('x^ZZ')
        sage: OT = atm.TermMonoid('O', G, QQ); OT
        O-Term Monoid x^ZZ with implicit coefficients in Rational Field
        sage: ET = atm.TermMonoid('exact', G, ZZ); ET
        Exact Term Monoid x^ZZ with coefficients in Integer Ring
    """
    def create_key_and_extra_args(self, term,
                                  growth_group=None, coefficient_ring=None,
                                  asymptotic_ring=None,
                                  **kwds):
        r"""
        Given the arguments and keyword, create a key that uniquely
        determines this object.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.TermMonoid.create_key_and_extra_args('O', G, QQ)
            ((<class 'sage.rings.asymptotic.term_monoid.OTermMonoid'>,
              Growth Group x^ZZ, Rational Field), {})
            sage: atm.TermMonoid.create_key_and_extra_args('exact', G, ZZ)
            ((<class 'sage.rings.asymptotic.term_monoid.ExactTermMonoid'>,
              Growth Group x^ZZ, Integer Ring), {})
            sage: atm.TermMonoid.create_key_and_extra_args('exact', G)
            Traceback (most recent call last):
            ...
            ValueError: A coefficient ring has to be specified
            to create a term monoid of type 'exact'
        """
        if isinstance(term, GenericTermMonoid):
            term_class = term.__class__
        elif term == 'O':
            term_class = OTermMonoid
        elif term == 'exact':
            term_class = ExactTermMonoid
        else:
            raise ValueError("Term specification '%s' has to be either 'exact' or 'O' "
                             "or an instance of an existing term." % term)

        if asymptotic_ring is not None and \
                (growth_group is not None or coefficient_ring is not None):
            raise ValueError("Input ambiguous: asymptotic ring %s as well as "
                             "growth group %s or coefficient ring %s are given." %
                             (asymptotic_ring, growth_group, coefficient_ring))

        if asymptotic_ring is not None:
            growth_group = asymptotic_ring.growth_group
            coefficient_ring = asymptotic_ring.coefficient_ring

        from sage.rings.asymptotic.growth_group import GenericGrowthGroup
        if not isinstance(growth_group, GenericGrowthGroup):
            raise ValueError("%s has to be an asymptotic growth group"
                             % growth_group)

        if coefficient_ring is None:
            raise ValueError("A coefficient ring has to be specified to "
                             "create a term monoid of type '%s'" % (term,))

        return (term_class, growth_group, coefficient_ring), kwds


    def create_object(self, version, key, **kwds):
        r"""
        Create a object from the given arguments.

        EXAMPLES::

            sage: import sage.rings.asymptotic.growth_group as agg
            sage: import sage.rings.asymptotic.term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: atm.TermMonoid('O', G, QQ)  # indirect doctest
            O-Term Monoid x^ZZ with implicit coefficients in Rational Field
            sage: atm.TermMonoid('exact', G, ZZ)  # indirect doctest
            Exact Term Monoid x^ZZ with coefficients in Integer Ring
        """
        term_class, growth_group, coefficient_ring = key
        return term_class(growth_group, coefficient_ring, **kwds)


TermMonoid = TermMonoidFactory("TermMonoid")
