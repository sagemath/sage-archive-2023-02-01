r"""
Asymptotic Ring

This module implements the central classes for computing with
asymptotic expressions. It provides the following two classes:

- :class:`AsymptoticExpression` -- this class essentially represents
  a sum of asymptotic terms (see
  :mod:`Asymptotic Terms <~sage.monoids.asymptotic_term_monoid>`).

- :class:`AsymptoticRing` -- parent structure for
  :class:`AsymptoticExpression`.

AUTHORS:

- Benjamin Hackl (2015-06): initial version
- Benjamin Hackl (2015-07): improvement user interface (short notation)


.. WARNING::

    As this code is experimental, a warning is thrown when an
    asymptotic ring (or an associated structure) is created for the
    first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: import sage.monoids.asymptotic_term_monoid as atm
        sage: G = agg.MonomialGrowthGroup(ZZ, 'x')
        doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
        See http://trac.sagemath.org/17600 for details.
        sage: T = atm.ExactTermMonoid(G, ZZ)
        doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
        See http://trac.sagemath.org/17715 for details.
        sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
        doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
        See http://trac.sagemath.org/17716 for details.

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
    Class for asymptotic expressions. ...

    INPUT:

    - ``parent`` -- the parent of the asymptotic expression.

    - ``poset`` -- a mutable poset representing the underlying
      growth structure.

    - ``simplify`` -- a boolean, controls automatic simplification
      (absorption) of asymptotic expressions. Default: ``True``.

    EXAMPLES:

    There are ways to create an asymptotic expression. First,
    we construct the corresponding parents::

        sage: R_x.<x> = AsymptoticRing('x^QQ', QQ)
        sage: import sage.groups.asymptotic_growth_group as agg
        sage: G = agg.GrowthGroup('y^ZZ')
        sage: R_y = AsymptoticRing(G, ZZ); y = R_y.gen()

    At this point, `x` and `y` are already asymptotic expressions::

        sage: (x, y)
        (x, y)
        sage: isinstance(x, sage.rings.asymptotic_ring.AsymptoticExpression)
        True
        sage: type(x)
        <class 'sage.rings.asymptotic_ring.AsymptoticRing_with_category.element_class'>

    The usual ring operations can be performed::

        sage: x^2 + 3*(x - x^2)
        3*x + -2*x^2
        sage: (3*x + 2)^3
        8 + 36*x + 54*x^2 + 27*x^3

    In addition to that, special powers (determined by the base ring
    of the growth group) can also be computed::

        sage: (x^(5/2) + x^(1/7)) * x^(-1/5)
        x^(-2/35) + x^(23/10)

    One of the central ideas behind computing with asymptotic
    expressions is that the `O`-notation (see
    :wikipedia:`Big_O_notation`) can be used. For example, we have::

        sage: (x + 2*x^2 + 3*x^3 + 4*x^4) * (O(x) + x^2)
        O(x^5) + 4*x^6

    In particular, :meth:`~sage.rings.big_oh.O` can be used to
    construct the asymptotic expressions. With the help of the
    ``poset``, we can also have a look at the inner structure
    of an asymptotic expression::

        sage: expr1 = x + 2*x^2 + 3*x^3 + 4*x^4; expr2 = O(x) + x^2
        sage: print(expr1.poset.repr_full())
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
        sage: print(expr2.poset.repr_full())
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
        sage: print((expr1 * expr2).poset.repr_full())
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

        sage: R_log.<xl> = AsymptoticRing('log(x)^QQ', QQ)
        sage: (O(xl) + xl^3)^4
        O(log(x)^10) + log(x)^12
    """
    def __init__(self, parent, poset, simplify=True):
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
            O(x^5) + 5*x^6
        """
        self._poset_ = poset
        if simplify:
            self._simplify_()
        super(AsymptoticExpression, self).__init__(parent=parent)


    @property
    def poset(self):
        r"""
        The poset of this asymptotic expression.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: expr = 7 * x^12 + x^5 + O(x^3)
            sage: expr.poset
            poset(O(x^3), x^5, 7*x^12)
        """
        return self._poset_


    def _simplify_(self):
        r"""
        Simplify this asymptotic expression.

        INPUT:

        Nothing.

        OUTPUT:

        An :class:`AsymptoticExpression`.

        .. NOTE::

            This asymptotic expression is simplified by letting
            `O`-terms that are included in this expression absorb all
            terms of lesser growth.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: OT = atm.TermMonoid('O', G); ET = atm.TermMonoid('exact', G, ZZ)
            sage: R = AsymptoticRing(G, ZZ)
            sage: lst = [ET(x,1), ET(x^2, 2), OT(x^3), ET(x^4, 4)]
            sage: expr = R(lst, simplify=False); expr  # indirect doctest
            x + 2*x^2 + O(x^3) + 4*x^4
            sage: expr._simplify_(); expr
            O(x^3) + 4*x^4
            sage: R(lst)  # indirect doctest
            O(x^3) + 4*x^4
        """
        from sage.monoids.asymptotic_term_monoid import OTerm
        for shell in self.poset.shells_topological(reverse=True):
            if shell.element.growth in self.poset and isinstance(shell.element,
                                                                 OTerm):
                self.poset.merge(shell.key)


    def _repr_(self, reverse=False):
        r"""
        A representation string for this asymptotic expression.

        INPUT:

        - ``reverse`` -- a boolean (default: ``False``).

        OUTPUT:

        A string.

        .. NOTE::

            By default, the elements with the weakest growth are
            printed first. If ``reverse`` is ``True``, then the
            printing order is reversed.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: expr = (5*x^2 + 12*x) * (x^3 + O(x))
            sage: repr(expr)  # indirect doctest
            'O(x^3) + 12*x^4 + 5*x^5'
            sage: expr._repr_(reverse=True)
            '5*x^5 + 12*x^4 + O(x^3)'
        """
        s = ' + '.join(repr(elem) for elem in
                       self.poset.elements_topological(include_special=False,
                                                       reverse=reverse))
        if s == '':
            return '0'
        else:
            return s


    def _add_(self, other):
        r"""
        Add ``other`` to this asymptotic expression.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpression`.

        OUTPUT:

        An :class:`AsymptoticExpression`.

        EXAMPLES::
            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: expr1 = x^123; expr2 = x^321
            sage: expr1._add_(expr2)
            x^123 + x^321
            sage: expr1 + expr2  # indirect doctest
            x^123 + x^321

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
        pst = self.poset.copy().union(other.poset)
        return self.parent()(poset=pst)


    def _sub_(self, other):
        r"""
        Subtract ``other`` from this asymptotic expression.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpression`.

        OUTPUT:

        An :class:`AsymptoticExpression`.

        .. NOTE::

            Subtraction of two asymptotic expressions is implemented
            by means of addition: `e_1 - e_2 = e_1 + (-1)\cdot e_2`.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: expr1 = x^123; expr2 = x^321
            sage: expr1 - expr2  # indirect doctest
            x^123 + -x^321
        """
        return self + (-1) * other

    def _mul_term_(self, other):
        r"""
        Helper method: multiply this asymptotic expression with the
        asymptotic term ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        An :class:`AsymptoticExpression`.

        TESTS::

            sage: import sage.monoids.asymptotic_term_monoid as atm
            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: T = atm.OTermMonoid(R.growth_group)
            sage: expr = 10*x^2 + O(x)
            sage: t = T(R.growth_group.gen())
            sage: expr._mul_term_(t)
            O(x^3)
        """
        return self.parent()([other * elem for elem in self.poset.elements()])


    def _mul_(self, other):
        r"""
        Multiply ``other`` to this asymptotic expression.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpression`.

        OUTPUT:

        An :class:`AsymptoticExpression`.

        .. TODO::

            The current implementation is the simple school book
            multiplication. More efficient variants like Karatsuba
            multiplication, or methods that exploit the structure
            of the underlying poset shall be implemented at a later
            point.

        EXAMPLES::

            sage: R.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: ex1 = 5*x^12
            sage: ex2 = x^3 + O(x)
            sage: ex1 * ex2  # indirect doctest
            O(x^13) + 5*x^15
        """
        return self.parent()(sum(self._mul_term_(term_other) for
                                 term_other in other.poset.elements()))


    def __pow__(self, power):
        r"""
        Return this element to the power of ``power``.

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
            ValueError: Exponent 1/7 not in base of growth group y^ZZ
            sage: (x^(1/2) + O(x^0))^15
            O(x^7) + x^(15/2)
        """
        P = self.parent()
        if len(self.poset._shells_) == 1:
            expr = self.poset.elements().next()
            if power in P.growth_group.base():
                from sage.monoids.asymptotic_term_monoid import TermWithCoefficient
                if isinstance(expr, TermWithCoefficient):
                    new_growth = expr.growth ** power
                    new_coef = expr.coefficient ** power
                    return P(expr.parent()(new_growth, new_coef))
                else:
                    new_growth = expr.growth ** power
                    return P(expr.parent()(new_growth))
            else:
                raise ValueError('Exponent %s not in base of growth group %s' %
                                 (power, P.growth_group))

        return super(AsymptoticExpression, self).__pow__(power)


    def O(self):
        r"""
        Convert all terms in this asymptotic expression to `O`-terms.

        INPUT:

        Nothing.

        OUTPUT:

        An asymptotic expression.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.GrowthGroup('x^ZZ')
            sage: AR = AsymptoticRing(MG, ZZ)
            sage: x = AR.create_term('exact', x)
            sage: expr = 42 * x^42 + x^10 + O(x^2); expr
            O(x^2) + x^10 + 42*x^42
            sage: expr.O()
            O(x^42)
        """
        if self.poset.null in self.poset.oo.predecessors():
            raise ValueError('O(%s) not defined' % self)
        else:
            return sum(self.parent().create_term('O', shell.element) for shell
                       in self.poset.oo.predecessors())



class AsymptoticRing(sage.rings.ring.Ring,
                     sage.structure.unique_representation.UniqueRepresentation):
    r"""
    Parent for asymptotic expressions.

    INPUT:

    - ``growth_group`` -- a partially ordered group (e.g. an instance
      of :class:`~sage.groups.asymptotic_growth_group.MonomialGrowthGroup`),
      or a short representation string of a growth group.

    - ``coefficient_ring`` -- the ring which contains the
      coefficients of the expressions.

    - ``category`` -- the category of the parent can be specified
      in order to broaden the base structure. It has to be a
      subcategory of ``Category of rings``. This is also the default
      category if ``None`` is specified.

    OUTPUT:

    An asymptotic ring.

    .. NOTE::

        See the following examples for more information on how to
        create an asymptotic ring.

    EXAMPLES:

    We begin with the explicit construction of an asymptotic ring,
    i.e. by explicitly specifying the underlying growth group::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: G_QQ = agg.GrowthGroup('x^QQ')
        sage: R_x = AsymptoticRing(growth_group=G_QQ, coefficient_ring=QQ); R_x
        Asymptotic Ring over x^QQ with coefficients from Rational Field

    Note that the coefficient ring of the asymptotic ring and the
    base ring of the underlying growth group do not need to
    coincide::

        sage: R_ZZ_x = AsymptoticRing(growth_group=G_QQ, coefficient_ring=ZZ); R_ZZ_x
        Asymptotic Ring over x^QQ with coefficients from Integer Ring

    As mentioned above, the short notation for growth groups can also
    be used to specify the underlying growth group. For now,
    representation strings of the form ``"variable^base"`` are
    allowed, where ``variable`` is some string, and ``base`` is
    either ``ZZ`` (for `\mathbb{Z}`), ``QQ`` (for `\mathbb{Q}`),
    or ``SR`` (for the symbolic ring). These strings correspond to
    monomial growth groups (see
    :class:`~sage.groups.asymptotic_growth_group.MonomialGrowthGroup`)::

        sage: R2_x = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); R2_x
        Asymptotic Ring over x^QQ with coefficients from Rational Field

    Alternatively, the preparser allows us to write::

        sage: R3_x.<x> = AsymptoticRing(growth_group='x^QQ', coefficient_ring=QQ); R3_x
        Asymptotic Ring over x^QQ with coefficients from Rational Field

    Note that this allows us to create logarithmic and polynomial
    growth groups::

        sage: R.<x> = AsymptoticRing('x^ZZ', QQ); R
        Asymptotic Ring over x^ZZ with coefficients from Rational Field
        sage: R_log.<lx> = AsymptoticRing('log(x)^ZZ', QQ); R_log
        Asymptotic Ring over log(x)^ZZ with coefficients from Rational Field

    According to the conventions for parents, uniqueness is ensured::

        sage: R_x is R2_x is R3_x
        True

    Furthermore, the coercion framework is also involved. Coercion
    between two asymptotic rings is possible (given that the
    underlying growth groups and coefficient rings are chosen
    appropriately)::

        sage: R_x.has_coerce_map_from(R_ZZ_x)
        True

    Additionally, for the sake of convenience, the coefficient ring
    also coerces into the asymptotic ring (representing constant
    quantities)::

        sage: R_x.has_coerce_map_from(QQ)
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

        __classcall__ unifies the input to the constructor of
        :class:`AsymptoticRing` such that the instances generated
        are unique. Also, this enables the use of the generation
        framework::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.GrowthGroup('x^ZZ')
            sage: AR1 = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: AR2.<x> = AsymptoticRing(growth_group='x^ZZ', coefficient_ring=ZZ)
            sage: AR1 is AR2
            True
        """
        if isinstance(growth_group, str):
            from sage.groups.asymptotic_growth_group import GrowthGroup
            growth_group = GrowthGroup(growth_group)

        return super(AsymptoticRing, cls).__classcall__(cls, growth_group,
                                                        coefficient_ring, category)


    @sage.misc.superseded.experimental(trac_number=17716)
    def __init__(self, growth_group=None, coefficient_ring=None, category=None):
        r"""
        See :class:`AsymptoticRing` for more information.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: R1 = AsymptoticRing(G, ZZ); R1
            Asymptotic Ring over x^ZZ with coefficients from Integer Ring
            sage: R2.<x> = AsymptoticRing('x^QQ', QQ); R2
            Asymptotic Ring over x^QQ with coefficients from Rational Field
            sage: R1 is R2
            False

        ::

            sage: R3 = AsymptoticRing(G)
            Traceback (most recent call last):
            ...
            TypeError: __classcall__() takes at least 3 arguments (2 given)
        """
        if growth_group is None:
            raise ValueError('Growth group not specified. Cannot continue.')
        elif coefficient_ring is None:
            raise ValueError('Coefficient ring not specified. Cannot continue.')
        elif not hasattr(coefficient_ring, 'is_ring') or\
                not coefficient_ring.is_ring():
            raise ValueError('%s has to be a ring.' % (coefficient_ring,))

        from sage.categories.rings import Rings
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

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.GrowthGroup('x^ZZ')
            sage: AR = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: AR.growth_group
            x^ZZ
        """
        return self._growth_group_


    @property
    def coefficient_ring(self):
        r"""
        The coefficient ring of this asymptotic ring.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.GrowthGroup('x^ZZ')
            sage: AR = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: AR.coefficient_ring
            Integer Ring
        """
        return self._coefficient_ring_


    def _element_constructor_(self, data, poset=None, simplify=True):
        r"""
        Convert a given object to this asymptotic ring.

        INPUT:

        - ``data`` -- an object representing the element to be
          initialized.

        - ``poset`` -- (default: None) if given, then this is
          directly passed to the element constructor (i.e., no
          conversion is performed).

        OUTPUT:

        An element of this growth group.

        TESTS::

            sage: AR.<x> = AsymptoticRing('x^ZZ', ZZ)
            sage: 3 * x^3
            3*x^3
            sage: 3 * x^3 + x
            x + 3*x^3
            sage: (3 * x^3) * (5 * x)
            15*x^4
            sage: 3 * x^3 + O(x^5)
            O(x^5)
        """
        if poset is None:
            if type(data) == self.element_class and data.parent() == self:
                return data
            elif isinstance(data, AsymptoticExpression):
                return self.element_class(self, data.poset, simplify)

            if data in self.coefficient_ring:
                if data != 0:
                    return self.create_term('exact', 1, data)
                else:
                    from sage.data_structures.mutable_poset import MutablePoset
                    from sage.monoids.asymptotic_term_monoid import \
                        can_absorb, absorption
                    poset = MutablePoset(key=lambda elem: elem.growth,
                                     can_merge=can_absorb,
                                     merge=absorption)
                    return self.element_class(self, poset, simplify)

            from sage.monoids.asymptotic_term_monoid import OTerm, ExactTerm
            from sage.data_structures.mutable_poset import MutablePoset
            from sage.monoids.asymptotic_term_monoid import can_absorb, \
                absorption
            if isinstance(data, (OTerm, ExactTerm)):
                data = (data,)

            try:
                data_iter = iter(data)
                if all(isinstance(elem, (OTerm, ExactTerm)) for elem in data_iter):
                    poset = MutablePoset(key=lambda elem: elem.growth,
                                         can_merge=can_absorb,
                                         merge=absorption)
                    poset = poset.union(data)
            except:
                raise TypeError('Input is ambiguous: cannot convert '
                                '%s to an asymptotic expression' % (data,))

        return self.element_class(self, poset, simplify)


    def _coerce_map_from_(self, R):
        r"""
        Return if ``R`` coerces into this asymptotic ring.

        INPUT:

        - ``R`` -- a parent.

        OUTPUT:

        A boolean.

        .. NOTE::

            There are two possible cases: either ``R`` is the
            ``coefficient_ring`` of this asymptotic ring, or ``R``
            itself is an asymptotic ring (where both the
            ``growth_group`` and the ``coefficient_ring`` coerce into
            the ``growth_group`` and the ``coefficient_ring`` of this
            asymptotic ring, respectively).

        TESTS::

            sage: AR_ZZ = AsymptoticRing('x^ZZ', coefficient_ring=ZZ); AR_ZZ
            Asymptotic Ring over x^ZZ with coefficients from Integer Ring
            sage: x_ZZ = AR_ZZ.gen()
            sage: AR_QQ = AsymptoticRing('x^QQ', coefficient_ring=QQ); AR_QQ
            Asymptotic Ring over x^QQ with coefficients from Rational Field
            sage: x_QQ = AR_QQ.gen()
            sage: AR_QQ.has_coerce_map_from(AR_ZZ)  # indirect doctest
            True
            sage: x_ZZ * x_QQ
            x^2

        ::

            sage: AR_QQ.has_coerce_map_from(QQ)
            True
            sage: 1/2 * x_QQ^2 + 7/8 * x_QQ^3
            1/2*x^2 + 7/8*x^3
        """
        if self.coefficient_ring.has_coerce_map_from(R):
            return True
        elif isinstance(R, AsymptoticRing):
            if self.growth_group.has_coerce_map_from(R.growth_group) and \
                    self.coefficient_ring.has_coerce_map_from(R.coefficient_ring):
                return True


    def _repr_(self):
        r"""
        A representation string for this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.GrowthGroup('x^ZZ')
            sage: AR = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: repr(AR)  # indirect doctest
            'Asymptotic Ring over x^ZZ with coefficients from Integer Ring'
        """
        return 'Asymptotic Ring over %s with coefficients ' \
               'from %s' % (self.growth_group, self.coefficient_ring)


    def gens(self):
        r"""
        Return a tuple with generators of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        A tuple.

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
            return self.create_term('exact', self.growth_group.gen(), 1),


    def gen(self, n=0):
        r"""
        Return the ``n``-th generator of this asymptotic ring.

        INPUT:

        - ``n`` -- a positive integer or `0`. Default: `0`.

        OUTPUT:

        An asymptotic expression.

        .. NOTE::

            Generators do not necessarily exist. This depends on the
            underlying growth group. For example, monomial growth
            groups have a generator, and exponential growth groups
            don't.

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

    def create_term(self, type, growth=None, coefficient=None):
        r"""
        Create a simple asymptotic expression consisting of a single
        term.

        INPUT:

        - ``type`` -- 'O' or 'exact',

        - ``growth`` -- a growth element,

        - ``coefficient`` -- a ring element.

        OUTPUT:

        An asymptotic expression.

        .. NOTE::

            This method calls the factory
            :class:`~sage.monoids.asymptotic_term_monoid.TermMonoid`
            with the appropriate arguments.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.GrowthGroup('x^ZZ')
            sage: R = AsymptoticRing(G, ZZ)
            sage: R.create_term('O', x^2)
            O(x^2)
            sage: R.create_term('exact', x^456, 123)
            123*x^456
        """
        from sage.monoids.asymptotic_term_monoid import TermMonoid
        if type == 'O':
            TM = TermMonoid(type, self.growth_group)
            return self(TM(growth))
        else:
            TM = TermMonoid(type, self.growth_group, self.coefficient_ring)
            return self(TM(growth, coefficient))
