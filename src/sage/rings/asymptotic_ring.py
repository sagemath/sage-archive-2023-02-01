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
    def __init__(self, parent, poset):
        r"""
        ...
        """
        for shell in poset.shells_topological(reverse=True):
            from sage.monoids.asymptotic_term_monoid import OTerm
            if isinstance(shell.element, OTerm) and shell.element.growth in poset:
                poset.merge(shell.key)
        self._poset_ = poset
        super(AsymptoticExpression, self).__init__(parent=parent)


    @property
    def poset(self):
        r"""

        """
        return self._poset_

    def _repr_(self, reverse=False):
        r"""
        A representation string for this asymptotic expression.
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
        """
        return self + other._mul_term_(
            self.parent().create_term('exact', growth=1, coefficient=-1))

    def _mul_term_(self, other):
        r"""
        Helper method: multiply this asymptotic expression with the
        asymptotic term ``other``.

        INPUT:

        - ``other`` -- an asymptotic term.

        OUTPUT:

        An :class:`AsymptoticExpression`.
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
        ...
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


    def _element_constructor_(self, data, poset=None):
        r"""
        Converts a given object to this asymptotic ring.

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
                return self.element_class(self, data.poset)

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
                    return self.element_class(self, poset)

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

        return self.element_class(self, poset)


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
        ...
        """
        from sage.monoids.asymptotic_term_monoid import TermMonoid
        if type == 'O':
            TM = TermMonoid(type, self.growth_group)
            return self(TM(growth))
        else:
            TM = TermMonoid(type, self.growth_group, self.coefficient_ring)
            return self(TM(growth, coefficient))
