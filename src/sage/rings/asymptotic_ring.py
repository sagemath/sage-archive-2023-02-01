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
        sage: R.<x> = AsymptoticRing('monomial', ZZ)
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
    ...
    """
    def __init__(self, parent, poset, simplify=True):
        r"""
        See :class:`AsymptoticExpression` for more information.

        TESTS::

            sage: R_x.<x> = AsymptoticRing('monomial', ZZ)
            sage: R_y.<y> = AsymptoticRing('monomial', ZZ)
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

            sage: R.<x> = AsymptoticRing('monomial', ZZ)
            sage: expr = 7 * x^12 + x^5 + O(x^3)
            sage: expr.poset
            poset(O(x^3), 1*x^5, 7*x^12)
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

            sage: R.<x> = AsymptoticRing('monomial', ZZ)
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
            sage: R.<x> = AsymptoticRing('monomial', ZZ)
            sage: expr1 = x^123; expr2 = x^321
            sage: expr1._add_(expr2)
            1*x^123 + 1*x^321
            sage: expr1 + expr2  # indirect doctest
            1*x^123 + 1*x^321

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

            sage: R.<x> = AsymptoticRing('monomial', ZZ)
            sage: expr1 = x^123; expr2 = x^321
            sage: expr1 - expr2  # indirect doctest
            1*x^123 + -1*x^321
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
            sage: R.<x> = AsymptoticRing('monomial', ZZ)
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

            sage: R.<x> = AsymptoticRing('monomial', ZZ)
            sage: ex1 = 5*x^12
            sage: ex2 = x^3 + O(x)
            sage: ex1 * ex2  # indirect doctest
            O(x^13) + 5*x^15
        """
        return self.parent()(sum(self._mul_term_(term_other) for
                                 term_other in other.poset.elements()))


    def O(self):
        r"""
        Convert all terms in this asymptotic expression to `O`-terms.

        INPUT:

        Nothing.

        OUTPUT:

        An asymptotic expression.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: AR = AsymptoticRing(MG, ZZ)
            sage: x = AR.create_term('exact', x)
            sage: expr = 42 * x^42 + x^10 + O(x^2); expr
            O(x^2) + 1*x^10 + 42*x^42
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
    ...
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
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: AR1 = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: AR2.<x> = AsymptoticRing(growth_group='monomial', coefficient_ring=ZZ)
            sage: AR1 is AR2
            True
        """
        if not isinstance(names, tuple):
            names = (names,)

        if len(names) > 1:
            raise NotImplementedError('Currently only one variable is supported')
        if isinstance(growth_group, str) and names is None:
            raise ValueError('names has to be specified if the growth group '
                             'is to be constructed')
        elif growth_group == 'monomial':
            from sage.groups.asymptotic_growth_group import MonomialGrowthGroup
            # coefficient ring acts as the base ring of the growth group!
            growth_group = MonomialGrowthGroup(coefficient_ring, names[0])
        elif hasattr(growth_group, '_var_') and names is not None:
            if names[0] == growth_group._var_:
                raise ValueError('Specified variable %s and growth group variable '
                                 '%s do not match' % (names[0], growth_group._var_))

        return super(AsymptoticRing, cls).__classcall__(cls, growth_group,
                                                        coefficient_ring, category)


    @sage.misc.superseded.experimental(trac_number=17716)
    def __init__(self, growth_group=None, coefficient_ring=None, category=None):
        r"""
        See :class:`AsymptoticRing` for more information.

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: G = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: R1 = AsymptoticRing(G, ZZ); R1
            Asymptotic Ring over Monomial Growth Group in x over Integer Ring with coefficients from Integer Ring
            sage: R2.<x> = AsymptoticRing('monomial', QQ); R2
            Asymptotic Ring over Monomial Growth Group in x over Rational Field with coefficients from Rational Field
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
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: AR = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: AR.growth_group
            Monomial Growth Group in x over Integer Ring
        """
        return self._growth_group_


    @property
    def coefficient_ring(self):
        r"""
        The coefficient ring of this asymptotic ring.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
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

            sage: AR.<x> = AsymptoticRing('monomial', ZZ)
            sage: 3 * x^3
            3*x^3
            sage: 3 * x^3 + x
            1*x + 3*x^3
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

            sage: AR_ZZ = AsymptoticRing('monomial', coefficient_ring=ZZ, names='x'); AR_ZZ
            Asymptotic Ring over Monomial Growth Group in x over Integer Ring with coefficients from Integer Ring
            sage: (x_ZZ,) = AR_ZZ.gens()
            sage: AR_QQ = AsymptoticRing('monomial', coefficient_ring=QQ, names='x'); AR_QQ
            Asymptotic Ring over Monomial Growth Group in x over Rational Field with coefficients from Rational Field
            sage: (x_QQ,) = AR_QQ.gens()
            sage: AR_QQ.has_coerce_map_from(AR_ZZ)  # indirect doctest
            True
            sage: x_ZZ * x_QQ
            1*x^2

        ::

            sage: AR_QQ.has_coerce_map_from(QQ)
            True
            sage: 1/2 * x_QQ^2 + 7/8 * x_QQ^3
            1/2*x^2 + 7/8*x^3
        """
        if R is self.coefficient_ring:
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
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: AR = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: repr(AR)  # indirect doctest
            'Asymptotic Ring over Monomial Growth Group in x over Integer Ring with coefficients from Integer Ring'
        """
        return 'Asymptotic Ring over %s with coefficients ' \
               'from %s' % (self.growth_group, self.coefficient_ring)


    def gens(self):
        r"""
        Return a tuple with generators of this asymptotic ring.

        INPUT:

        Noting.

        OUTPUT:

        A tuple.

        .. NOTE::

            Generators do not necessarily exist. This depends on the
            underlying growth group. For example, monomial growth
            groups have a generator, and exponential growth groups
            don't.

        EXAMPLES::

            sage: AR.<x> = AsymptoticRing('monomial', ZZ)
            sage: AR.gens()
            (1*x,)
        """
        from sage.groups.asymptotic_growth_group import MonomialGrowthGroup
        if isinstance(self.growth_group, MonomialGrowthGroup):
            return self.create_term('exact', self.growth_group.gen(), 1),

    def ngens(self):
        r"""
        Return the number of generators of this asymptotic ring.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: AR.<x> = AsymptoticRing('monomial', ZZ)
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
            sage: G = agg.MonomialGrowthGroup(ZZ, 'x')
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
