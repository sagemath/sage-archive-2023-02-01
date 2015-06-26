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

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: MG = agg.MonomialGrowthGroup(ZZ, 'x')
            sage: AR = AsymptoticRing(growth_group=MG, coefficient_ring=ZZ)
            sage: AR.create_term('exact', x^3, 3)
            3*x^3
            sage: AR.create_term('exact', x^3, 3) + AR.create_term('exact', x, 1)
            1*x + 3*x^3
            sage: AR.create_term('exact', x^3, 3) * AR.create_term('exact', x, 5)
            15*x^4
            sage: AR.create_term('exact', x^3, 3) + AR.create_term('O', x^5)
            O(x^5)
        """
        if poset is None:
            if type(data) == self.element_class and data.parent() == self:
                return data
            if data in self.coefficient_ring:
                if data != 0:
                    return self.create_term('exact', 1, data)
                else:
                    from sage.data_structures.mutable_poset import MutablePoset
                    from sage.monoids.asymptotic_term_monoid import \
                        _can_absorb_, _absorption_
                    poset = MutablePoset(key=lambda elem: elem.growth,
                                     can_merge=_can_absorb_,
                                     merge=_absorption_)
                    return self.element_class(self, poset)

            from sage.monoids.asymptotic_term_monoid import OTerm, ExactTerm
            from sage.data_structures.mutable_poset import MutablePoset
            from sage.monoids.asymptotic_term_monoid import _can_absorb_, \
                _absorption_
            if isinstance(data, (OTerm, ExactTerm)):
                data = (data,)

            try:
                data_iter = iter(data)
                if all(isinstance(elem, (OTerm, ExactTerm)) for elem in data_iter):
                    poset = MutablePoset(key=lambda elem: elem.growth,
                                         can_merge=_can_absorb_,
                                         merge=_absorption_)
                    poset = poset.union(data)
            except:
                raise TypeError('Input is ambiguous: cannot convert '
                                '%s to an asymptotic expression' % (data,))

        return self.element_class(self, poset)


    def _coerce_map_from_(self, R):
        r"""
        ...
        """
        if R is self.coefficient_ring:
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
