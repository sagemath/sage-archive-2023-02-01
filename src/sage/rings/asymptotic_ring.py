r"""
Asymptotic Ring

This module implements a ring (called :class:`AsymptoticRing`) for
computations with :class:`asymptotic expressions
<AsymptoticExpression>`.

Definition
==========

An asymptotic expression is a sum; its summands are the following:

- Exact terms `c\cdot g` with a coefficient `c` and an element `g` of
  an :ref:`growth group <asymptotic_ring_growth>`.

- `O`-terms `O(g)` (see :wikipedia:`Big_O_notation`) for some
  :mod:`growth group element <sage.groups.asymptotic_growth_group>`
  `g` (:ref:`see below <asymptotic_ring_growth>`).

Examples of such elements can found :ref:`below <asymptotic_ring_intro>`.

.. _asymptotic_ring_growth:

Growth Elements
---------------

The elements of a :mod:`growth group
<sage.groups.asymptotic_growth_group>` are equipped with a partial
ordering and usually contains a variable. Examples are (among many
other possibilities)

- elements of the form `z^q` for some integer or rational `q` (growth
  groups ``z^ZZ`` or ``z^QQ``),

- elements of the form `log(z)^q` for some integer or rational `q` (growth
  groups ``log(z)^ZZ`` or ``log(z)^QQ``),

- elements of the form `a^z` for some
  rational `a` (growth group ``QQ^z``), or

- more sophisticated constructions like products `x^r log(x)^s \cdot
  a^y \cdot y^q` (this corresponds to an element of the growth group
  ``x^QQ * \log(x)^ZZ * QQ^y * y^QQ``).

The ordering in all these examples is the growth as `x`, `y`, or `z`
(independently) tend to `\infty`. For elements only using the
variable `z` this means, `g_1 \leq g_2` if

.. MATH::

    \lim_{z\to\infty} \frac{g_2}{g_1} \leq 1.

.. WARNING::

    As this code is experimental, a warning is thrown when an
    asymptotic ring (or an associated structure) is created for the
    first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: G = agg.MonomialGrowthGroup(ZZ, 'x')
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        sage: T = atm.ExactTermMonoid(G, ZZ)
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)

.. _asymptotic_ring_intro:

Introductory Examples
=====================

First, we construct the following (very simple) asymptotic ring in the variable `z`::

    sage: A.<z> = AsymptoticRing(growth_group='z^QQ', coefficient_ring=ZZ); A
    Asymptotic Ring <z^QQ> over Integer Ring

A typical element of this ring is

::

    sage: A.an_element()  # not tested
    -z^(3/2) + O(z^(1/2))

This element consists of two summands: the exact term with coefficient
`-1` and growth `x^{3/2}` and the `O`-term `O(x^{1/2})`. Note that the
growth of `x^{3/2}` is larger than the growth of `x^{1/2}` as
`x\to\infty`, thus this expression cannot be simplified (which would
be done automatically, see below).

Next, we construct a more sophisticated asymptotic ring in the
variables `x` and `y` by

::

    sage: B.<x, y> = AsymptoticRing(growth_group='x^QQ * \log(x)^ZZ * QQ^y * y^QQ', coefficient_ring=QQ); B  # not tested

Again, we can look at a typical element::

    sage: B.an_element()  # not tested

Arithemtical Operations
-----------------------

With the asymptotic rings constructed above (or more precisely with
their elements) we can do a lot of different arithmetical
calculations.

We start our calculations in the ring

::

    sage: A
    Asymptotic Ring <z^QQ> over Integer Ring

Of course, we can perform the usual ring operations `+` and `*`::

    sage: z^2 + 3*z*(1 - z)
    -2*z^2 + 3*z
    sage: (3*z + 2)^3
    27*z^3 + 54*z^2 + 36*z + 8

In addition to that, special powers---our growth group ``z^QQ`` allows
the exponents to be out of `\mathbb{Q}`---can also be computed::

    sage: (z^(5/2) + z^(1/7)) * z^(-1/5)
    z^(23/10) + z^(-2/35)

The central concepts of computations with asymptotic expressions is
that the `O`-notation can be used. For example, we have

::

    sage: z^3 + z^2 + z + O(z^2)
    z^3 + O(z^2)

and more advanced

::

    sage: (z + 2*z^2 + 3*z^3 + 4*z^4) * (O(z) + z^2)
    4*z^6 + O(z^5)

.. TODO::

   inversions

.. TODO::

    arithmetic in the ring

    ::

        sage: B  # not tested

More Examples
=============

.. TODO::

    write more examples

AUTHORS:

- Benjamin Hackl (2015-06): initial version
- Benjamin Hackl (2015-07): improvement user interface (short notation)
- Daniel Krenn (2015-08): various improvents, review; documentation
"""

# *****************************************************************************
# Copyright (C) 2015 Benjamin Hackl <benjamin.hackl@aau.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************

import sage

class AsymptoticExpression(sage.rings.ring_element.RingElement):
    r"""
    Class for asymptotic expressions, i.e., the elements of an
    :class:`AsymptoticRing`.

    INPUT:

    - ``parent`` -- the parent of the asymptotic expression.

    - ``summands`` -- the summands as a
      :class:`~sage.data_structures.mutable_poset.MutablePoset`, which
      represents the underlying structure.

    - ``simplify`` -- a boolean (default: ``True``). It controls
      automatic simplification (absorption) of the asymptotic expression.

    EXAMPLES:

    There are several ways to create asymptotic expressions; usually
    this is done by using the corresponding rings/parents::

        sage: R_x.<x> = AsymptoticRing('x^QQ', QQ); R_x
        Asymptotic Ring <x^QQ> over Rational Field
        sage: R_y.<y> = AsymptoticRing('y^ZZ', ZZ); R_y
        Asymptotic Ring <y^ZZ> over Integer Ring

    At this point, `x` and `y` are already asymptotic expressions::

        sage: type(x)
        <class 'sage.rings.asymptotic_ring.AsymptoticRing_with_category.element_class'>

    The usual ring operations, but allowing rational exponents (growth
    group ``x^QQ``) can be performed::

        sage: x^2 + 3*(x - x^(2/5))
        x^2 + 3*x - 3*x^(2/5)
        sage: (3*x^(1/3) + 2)^3
        27*x + 54*x^(2/3) + 36*x^(1/3) + 8

    One of the central ideas behind computing with asymptotic
    expressions is that the `O`-notation (see
    :wikipedia:`Big_O_notation`) can be used. For example, we have::

        sage: (x + 2*x^2 + 3*x^3 + 4*x^4) * (O(x) + x^2)
        4*x^6 + O(x^5)

    In particular, :meth:`~sage.rings.big_oh.O` can be used to
    construct the asymptotic expressions. With the help of the
    :meth:`summands`, we can also have a look at the inner structure
    of an asymptotic expression::

        sage: expr1 = x + 2*x^2 + 3*x^3 + 4*x^4; expr2 = O(x) + x^2
        sage: print(expr1.summands.repr_full())
        poset(x, 2*x^2, 3*x^3, 4*x^4)
        +-- null
        |   +-- no predecessors
        |   +-- successors:   x
        +-- x
        |   +-- predecessors:   null
        |   +-- successors:   2*x^2
        +-- 2*x^2
        |   +-- predecessors:   x
        |   +-- successors:   3*x^3
        +-- 3*x^3
        |   +-- predecessors:   2*x^2
        |   +-- successors:   4*x^4
        +-- 4*x^4
        |   +-- predecessors:   3*x^3
        |   +-- successors:   oo
        +-- oo
        |   +-- predecessors:   4*x^4
        |   +-- no successors
        sage: print(expr2.summands.repr_full())
        poset(O(x), x^2)
        +-- null
        |   +-- no predecessors
        |   +-- successors:   O(x)
        +-- O(x)
        |   +-- predecessors:   null
        |   +-- successors:   x^2
        +-- x^2
        |   +-- predecessors:   O(x)
        |   +-- successors:   oo
        +-- oo
        |   +-- predecessors:   x^2
        |   +-- no successors
        sage: print((expr1 * expr2).summands.repr_full())
        poset(O(x^5), 4*x^6)
        +-- null
        |   +-- no predecessors
        |   +-- successors:   O(x^5)
        +-- O(x^5)
        |   +-- predecessors:   null
        |   +-- successors:   4*x^6
        +-- 4*x^6
        |   +-- predecessors:   O(x^5)
        |   +-- successors:   oo
        +-- oo
        |   +-- predecessors:   4*x^6
        |   +-- no successors

    In addition to the monomial growth elements from above, we can
    also compute with logarithmic terms (simply by constructing the
    appropriate growth group)::

        sage: R_log = AsymptoticRing('log(x)^QQ', QQ)
        sage: lx = R_log(log(SR.var('x')))
        sage: (O(lx) + lx^3)^4
        log(x)^12 + O(log(x)^10)
    """
    def __init__(self, parent, summands, simplify=True):
        r"""
        See :class:`AsymptoticExpression` for more information.

        TESTS::

            sage: R_x.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: R_y.<y> = AsymptoticRing('y^ZZ', ZZ)
            sage: R_x is R_y
            False
            sage: ex1 = x + 2*x^2 + 3*x^3 + 4*x^4 + 5*x^5
            sage: ex2 = x + O(R_x(1))
            sage: ex1 * ex2
            5*x^6 + O(x^5)
        """
        super(AsymptoticExpression, self).__init__(parent=parent)

        self._summands_ = summands
        if simplify:
            self._simplify_()


    @property
    def summands(self):
        r"""
        The summands of this asymptotic expression stored in the
        underlying data structure (a
        :class:`~sage.data_structures.mutable_poset.MutablePoset`).

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: expr = 7 * x^12 + x^5 + O(x^3)
            sage: expr.summands
            poset(O(x^3), x^5, 7*x^12)
        """
        return self._summands_


    def __nonzero__(self):
        r"""
        Return if this asymptotic expression is not identially zero.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        TESTS::

            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: bool(R(0))  # indirect doctest
            False
            sage: bool(x)  # indirect doctest
            True
            sage: bool(7 * x^12 + x^5 + O(x^3))  # indirect doctest
            True
        """
        return bool(self._summands_)


    def _simplify_(self):
        r"""
        Simplify this asymptotic expression.

        INPUT:

        Nothing.

        OUTPUT:

        Nothing, but modifies this asymptotic expression.

        .. NOTE::

            This method is usually called during initialization of
            this asymptotic expression.

        .. NOTE::

            This asymptotic expression is simplified by letting
            `O`-terms that are included in this expression absorb all
            terms with smaller growth.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: OT = atm.TermMonoid('O', G); ET = atm.TermMonoid('exact', G, ZZ)
            sage: R = AsymptoticRing(G, ZZ)
            sage: lst = [ET(x,1), ET(x^2, 2), OT(x^3), ET(x^4, 4)]
            sage: expr = R(lst, simplify=False); expr  # indirect doctest
            4*x^4 + O(x^3) + 2*x^2 + x
            sage: expr._simplify_(); expr
            4*x^4 + O(x^3)
            sage: R(lst)  # indirect doctest
            4*x^4 + O(x^3)
        """
        self.summands.merge(reverse=True)


    def _repr_(self):
        r"""
        A representation string for this asymptotic expression.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: (5*x^2 + 12*x) * (x^3 + O(x))  # indirect doctest
            5*x^5 + 12*x^4 + O(x^3)
            sage: (5*x^2 - 12*x) * (x^3 + O(x))  # indirect doctest
            5*x^5 - 12*x^4 + O(x^3)
        """
        s = ' + '.join(repr(elem) for elem in
                       self.summands.elements_topological(reverse=True))
        s = s.replace('+ -', '- ')
        if not s:
            return '0'
        return s


    def _add_(self, other):
        r"""
        Add ``other`` to this asymptotic expression.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpression`.

        OUTPUT:

        The sum as an :class:`AsymptoticExpression`.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: expr1 = x^123; expr2 = x^321
            sage: expr1._add_(expr2)
            x^321 + x^123
            sage: expr1 + expr2  # indirect doctest
            x^321 + x^123

        If an `O`-term is added to an asymptotic expression, then
        the `O`-term absorbs everything it can::

            sage: x^123 + x^321 + O(x^555)  # indirect doctest
            O(x^555)

        TESTS::

            sage: x + O(x)
            O(x)
            sage: O(x) + x
            O(x)
        """
        smds = self.summands.copy().union(other.summands)
        return self.parent()(summands=smds)


    def _sub_(self, other):
        r"""
        Subtract ``other`` from this asymptotic expression.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpression`.

        OUTPUT:

        The difference as an :class:`AsymptoticExpression`.

        .. NOTE::

            Subtraction of two asymptotic expressions is implemented
            by means of addition: `e_1 - e_2 = e_1 + (-1)\cdot e_2`.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: expr1 = x^123; expr2 = x^321
            sage: expr1 - expr2  # indirect doctest
            -x^321 + x^123
            sage: O(x) - O(x)
            O(x)
        """
        return self + (-1) * other


    def _mul_term_(self, term):
        r"""
        Helper method: multiply this asymptotic expression with the
        asymptotic term ``term``.

        INPUT:

        - ``term`` -- an asymptotic term (see
          :mod:`~sage.monoids.asymptotic_term_monoid`).

        OUTPUT:

        The product as an :class:`AsymptoticExpression`.

        TESTS::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: T = atm.OTermMonoid(R.growth_group)
            sage: expr = 10*x^2 + O(x)
            sage: t = T(R.growth_group.gen())
            sage: expr._mul_term_(t)
            O(x^3)
        """
        return self.parent()([term * elem for elem in self.summands.elements()])


    def _mul_(self, other):
        r"""
        Multiply ``other`` to this asymptotic expression.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpression`.

        OUTPUT:

        The product as an :class:`AsymptoticExpression`.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: ex1 = 5*x^12
            sage: ex2 = x^3 + O(x)
            sage: ex1 * ex2  # indirect doctest
            5*x^15 + O(x^13)

        .. TODO::

            The current implementation is the school book
            multiplication. More efficient variants like Karatsuba
            multiplication, or methods that exploit the structure
            of the underlying poset shall be implemented at a later
            point.
        """
        return self.parent()(sum(self._mul_term_(term_other) for
                                 term_other in other.summands.elements()))


    def __pow__(self, power):
        r"""
        Takes this element to the given ``power``.

        INPUT:

        - ``power`` -- an element.

        OUTPUT:

        An asymptotic expression.

        TESTS::

            sage: R_QQ.<x> = AsymptoticRing('x^QQ', QQ)
            sage: x^(1/7)
            x^(1/7)
            sage: R_ZZ.<y> = AsymptoticRing('y^ZZ', ZZ)
            sage: y^(1/7)
            Traceback (most recent call last):
            ...
            ValueError: Growth Group y^ZZ disallows taking y to the power of 1/7.
            sage: (x^(1/2) + O(x^0))^15
            x^(15/2) + O(x^7)
        """
        if len(self.summands) > 1:
            from sage.rings.integer_ring import ZZ
            if power not in ZZ:
                raise NotImplementedError('Taking the sum %s to the '
                                          'non-integer power %s not '
                                          'implemented.' % (self, power))
            return super(AsymptoticExpression, self).__pow__(power)

        P = self.parent()
        if power not in P.growth_group.base():
            raise ValueError('%s disallows taking %s '
                             'to the power of %s.' %
                             (P.growth_group, self, power))

        from sage.monoids.asymptotic_term_monoid import TermWithCoefficient
        expr = self.summands.elements().next()
        if isinstance(expr, TermWithCoefficient):
            new_growth = expr.growth**power
            new_coeff = expr.coefficient**power
            return P(expr.parent()(new_growth, new_coeff))
        else:
            new_growth = expr.growth**power
            return P(expr.parent()(new_growth))



    def O(self):
        r"""
        Convert all terms in this asymptotic expression to `O`-terms.

        INPUT:

        Nothing.

        OUTPUT:

        An asymptotic expression.

        EXAMPLES::

            sage: AR.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: O(x)
            O(x)
            sage: expr = 42 * x^42 + x^10 + O(x^2); expr
            42*x^42 + x^10 + O(x^2)
            sage: expr.O()
            O(x^42)

        TESTS::

            sage: O(AR(0))
            Traceback (most recent call last):
            ...
            ValueError: Cannot build O(0).
        """
        if not self:
            raise ValueError('Cannot build O(%s).' % (self,))
        return sum(self.parent().create_summand('O', growth=element)
                   for element in self.summands.maximal_elements())



class AsymptoticRing(sage.rings.ring.Ring,
                     sage.structure.unique_representation.UniqueRepresentation):
    r"""
    A ring consisting of :class:`asymptotic expressions <AsymptoticExpression>`.

    INPUT:

    - ``growth_group`` -- either a partially ordered group (see
      :mod:`~sage.groups.asymptotic_growth_group`) or a string
      describing such a growth group (see
      :class:`~sage.groups.asymptotic_growth_group.GrowthGroupFactory`).

    - ``coefficient_ring`` -- the ring which contains the
      coefficients of the expressions.

    - ``category`` -- the category of the parent can be specified
      in order to broaden the base structure. It has to be a
      subcategory of ``Category of rings``. This is also the default
      category if ``None`` is specified.

    EXAMPLES:

    We begin with the construction of an asymptotic ring in various
    ways. First, we simply pass a string specifying the underlying
    growth group::

        sage: R1_x.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); R1_x
        Asymptotic Ring <x^QQ> over Rational Field
        sage: x
        x

    This is equivalent to the following code, which explicitly
    specifies the underlying growth group::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: G_QQ = agg.GrowthGroup('x^QQ')
        sage: R2_x.<x> = AsymptoticRing(growth_group=G_QQ, coefficient_ring=QQ); R2_x
        Asymptotic Ring <x^QQ> over Rational Field

    Of course, the coefficient ring of the asymptotic ring and the
    base ring of the underlying growth group do not need to
    coincide::

        sage: R_ZZ_x.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=ZZ); R_ZZ_x
        Asymptotic Ring <x^QQ> over Integer Ring

    Note, we can also create and use logarithmic growth groups::

        sage: R_log = AsymptoticRing('log(x)^ZZ', QQ); R_log
        Asymptotic Ring <log(x)^ZZ> over Rational Field

    Other growth groups are available. See :mod:`~sage.rings.asymptotic_ring` for
    a lot more examples.

    Below there are some technical details.

    According to the conventions for parents, uniqueness is ensured::

        sage: R1_x is R2_x
        True

    Furthermore, the coercion framework is also involved. Coercion
    between two asymptotic rings is possible (given that the
    underlying growth groups and coefficient rings are chosen
    appropriately)::

        sage: R1_x.has_coerce_map_from(R_ZZ_x)
        True

    Additionally, for the sake of convenience, the coefficient ring
    also coerces into the asymptotic ring (representing constant
    quantities)::

        sage: R1_x.has_coerce_map_from(QQ)
        True

    TESTS::

        sage: R3_x = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); R3_x
        Asymptotic Ring <x^QQ> over Rational Field
        sage: R1_x is R2_x is R3_x
        True
    """
    # enable the category framework for elements
    Element = AsymptoticExpression


    @staticmethod
    def __classcall__(cls, growth_group, coefficient_ring, names=None,
                      category=None):
        r"""
        Normalizes the input in order to ensure a unique
        representation of the parent.

        For more information see :class:`AsymptoticRing`.

        EXAMPLES:

        ``__classcall__`` unifies the input to the constructor of
        :class:`AsymptoticRing` such that the instances generated
        are unique. Also, this enables the use of the generation
        framework::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.GrowthGroup('x^ZZ')
            sage: AR1 = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: AR2.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR1 is AR2
            True

        The bracket notation can only be used if the growth group
        has a generator::

            sage: AR.<lx> = AsymptoticRing('log(x)^ZZ', ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Growth Group log(x)^ZZ does not have a generator.
        """
        if isinstance(growth_group, str):
            from sage.groups.asymptotic_growth_group import GrowthGroup
            growth_group = GrowthGroup(growth_group)

        if names is not None and not growth_group.gens_monomial():
            raise ValueError("%s does not have a generator." % (growth_group,))

        return super(AsymptoticRing, cls).__classcall__(cls, growth_group,
                                                        coefficient_ring,
                                                        category)


    @sage.misc.superseded.experimental(trac_number=17601)
    def __init__(self, growth_group, coefficient_ring, category=None):
        r"""
        See :class:`AsymptoticRing` for more information.

        TESTS::

            sage: R1 = AsymptoticRing('x^ZZ', ZZ); R1
            Asymptotic Ring <x^ZZ> over Integer Ring
            sage: R2.<x> = AsymptoticRing('x^QQ', QQ); R2
            Asymptotic Ring <x^QQ> over Rational Field
            sage: R1 is R2
            False

        ::

            sage: R3 = AsymptoticRing('x^ZZ')
            Traceback (most recent call last):
            ...
            TypeError: __classcall__() takes at least 3 arguments (2 given)
        """
        from sage.categories.rings import Rings

        if growth_group is None:
            raise ValueError('Growth group not specified. Cannot continue.')
        elif coefficient_ring is None:
            raise ValueError('Coefficient ring not specified. Cannot continue.')
        elif coefficient_ring not in Rings():
            raise ValueError('%s is not a ring. Cannot continue.' % (coefficient_ring,))

        if category is None:
            category = Rings()
        else:
            if not isinstance(category, tuple):
                category = (category,)
            if not any(cat.is_subcategory(Rings()) for cat in category):
                raise ValueError('%s is not a subcategory of %s' % (category,
                                 Rings()))

        self._coefficient_ring_ = coefficient_ring
        self._growth_group_ = growth_group
        super(AsymptoticRing, self).__init__(base=coefficient_ring,
                                             category=category)


    @property
    def growth_group(self):
        r"""
        The growth group of this asymptotic ring.

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.growth_group
            Growth Group x^ZZ
        """
        return self._growth_group_


    @property
    def coefficient_ring(self):
        r"""
        The coefficient ring of this asymptotic ring.

        EXAMPLES::

            sage: AR = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR.coefficient_ring
            Integer Ring
        """
        return self._coefficient_ring_


    @staticmethod
    def _create_empty_summands_():
        r"""
        Create an empty data structure suitable for storing and working
        with summands.

        INPUT:

        Nothing.

        OUTPUT:

        A :class:`~sage.data_structures.mutable_poset.MutablePoset`.

        TESTS::

            sage: AsymptoticRing._create_empty_summands_()
            poset()
        """
        from sage.data_structures.mutable_poset import MutablePoset
        from sage.monoids.asymptotic_term_monoid import \
            can_absorb, absorption
        return MutablePoset(key=lambda element: element.growth,
                            can_merge=can_absorb,
                            merge=absorption)


    def _element_constructor_(self, data, summands=None, simplify=True):
        r"""
        Convert a given object to this asymptotic ring.

        INPUT:

        - ``data`` -- an object representing the element to be
          initialized.

        - ``summands`` -- (default: ``None``) if given, then this is
          directly passed to the element constructor (i.e., no
          conversion is performed).

        - ``simplify`` -- (default: ``True``) if set, then the constructed
          element is simplified (terms are absorbed) automatically.

        OUTPUT:

        An element of this asymptotic ring.

        TESTS::

            sage: AR = AsymptoticRing('x^ZZ', ZZ)
            sage: AR(5)
            5
            sage: AR(3*x^2)
            3*x^2
            sage: x = ZZ['x'].gen(); x.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: AR(x)
            x
            sage: y = ZZ['y'].gen(); AR(y)
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert y to an asymptotic expression.
        """
        if summands is not None:
            if type(data) != int or data != 0:
                raise ValueError('Input is ambigous: '
                                 '%s as well as summands=%s '
                                 'are specified.' % (data, summands))
            return self.element_class(self, summands, simplify=simplify)

        if type(data) == self.element_class and data.parent() == self:
            return data

        if isinstance(data, AsymptoticExpression):
            return self.element_class(self, data.summands, simplify=simplify)

        from sage.monoids.asymptotic_term_monoid import GenericTerm
        if isinstance(data, GenericTerm):
            data = (data,)

        if isinstance(data, (list, tuple)):
            if not all(isinstance(elem, GenericTerm) for elem in data):
                raise TypeError('Not all list entries of %s '
                                'are asymptotic terms.' % (data,))
            summands = AsymptoticRing._create_empty_summands_()
            summands.union_update(data)
            return self.element_class(self, summands, simplify=simplify)

        if data == 0:
            summands = AsymptoticRing._create_empty_summands_()
            return self.element_class(self, summands, simplify=simplify)

        try:
            summand = self.create_summand('exact', growth=data)
        except (TypeError, ValueError):
            pass
        else:
            return summand

        try:
            coefficient = self.coefficient_ring(data)
        except (TypeError, ValueError):
            pass
        else:
            return self.create_summand('exact', growth=1, coefficient=coefficient)

        raise TypeError('Cannot convert %s to an asymptotic '
                        'expression.' % (data,))


    def _coerce_map_from_(self, R):
        r"""
        Return if ``R`` coerces into this asymptotic ring.

        INPUT:

        - ``R`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            There are two possible cases: either ``R`` coerces in the
            :meth:`coefficient_ring` of this asymptotic ring, or ``R``
            itself is an asymptotic ring, where both the
            meth:`growth_group` and the :meth:`coefficient_ring` coerce into
            the :meth:`growth_group` and the :meth:`coefficient_ring` of this
            asymptotic ring, respectively.

        TESTS::

            sage: AR_ZZ = AsymptoticRing('x^ZZ', coefficient_ring=ZZ); AR_ZZ
            Asymptotic Ring <x^ZZ> over Integer Ring
            sage: x_ZZ = AR_ZZ.gen()
            sage: AR_QQ = AsymptoticRing('x^QQ', coefficient_ring=QQ); AR_QQ
            Asymptotic Ring <x^QQ> over Rational Field
            sage: x_QQ = AR_QQ.gen()
            sage: AR_QQ.has_coerce_map_from(AR_ZZ)  # indirect doctest
            True
            sage: x_ZZ * x_QQ
            x^2

        ::

            sage: AR_QQ.has_coerce_map_from(QQ)
            True
            sage: AR_QQ.has_coerce_map_from(ZZ)
            True
        """
        if self.coefficient_ring.has_coerce_map_from(R):
            return True
        elif isinstance(R, AsymptoticRing):
            if self.growth_group.has_coerce_map_from(R.growth_group) and \
                    self.coefficient_ring.has_coerce_map_from(R.coefficient_ring):
                return True


    def _repr_(self):
        r"""
        A representation string of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.GrowthGroup('x^ZZ')
            sage: AR = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: repr(AR)  # indirect doctest
            'Asymptotic Ring <x^ZZ> over Integer Ring'
        """
        try:
            G = '<' + self.growth_group._repr_(condense=True) + '>'
        except TypeError:
            G = repr(self.growth_group)
        return 'Asymptotic Ring %s over %s' % (G, self.coefficient_ring)


    def gens(self):
        r"""
        Return a tuple with generators of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple of asymptotic expressions.

        .. NOTE::

            Generators do not necessarily exist. This depends on the
            underlying growth group. For example, monomial growth
            groups have a generator, and exponential growth groups
            don't.

        EXAMPLES::

            sage: AR.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: AR.gens()
            (x,)
        """
        from sage.groups.asymptotic_growth_group import MonomialGrowthGroup
        if isinstance(self.growth_group, MonomialGrowthGroup):
            return self.create_summand('exact', growth=self.growth_group.gen(), coefficient=1),


    def gen(self, n=0):
        r"""
        Return the ``n``-th generator of this asymptotic ring.

        INPUT:

        - ``n`` -- (default: `0`) a non-negative integer.

        OUTPUT:

        An asymptotic expression.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: R.gen()
            x
        """
        return self.gens()[n]


    def ngens(self):
        r"""
        Return the number of generators of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: AR.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: AR.ngens()
            1
        """
        from sage.groups.asymptotic_growth_group import MonomialGrowthGroup
        if isinstance(self.growth_group, MonomialGrowthGroup):
            return 1
        else:
            return 0


    def create_summand(self, type, growth, **kwds):
        r"""
        Create a simple asymptotic expression consisting of a single
        summand.

        INPUT:

        - ``type`` -- 'O' or 'exact'.

        - ``growth`` -- an element of the :meth:`growth_group`.

        - ``coefficient`` -- an element of the :meth:`coefficient_ring`.

        OUTPUT:

        An asymptotic expression.

        .. NOTE::

            This method calls the factory :class:`TermMonoid
            <sage.monoids.asymptotic_term_monoid.TermMonoidFactory>`
            with the appropriate arguments.

        EXAMPLES::

            sage: R = AsymptoticRing('x^ZZ', ZZ)
            sage: R.create_summand('O', growth=x^2)
            O(x^2)
            sage: R.create_summand('exact', growth=x^456, coefficient=123)
            123*x^456
        """
        from sage.monoids.asymptotic_term_monoid import TermMonoid
        TM = TermMonoid(type, self.growth_group, self.coefficient_ring)

        if type == 'exact' and kwds.get('coefficient') == 0:
            return self(kwds['coefficient'])

        return self(TM(growth, **kwds))
