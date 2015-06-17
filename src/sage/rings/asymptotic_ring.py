r"""
Asymptotic Ring

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

def _absorption_(left, right):
    r"""
    Helper method for the mutable poset associated to our
    """
    try:
        return left.absorb(right)
    except ArithmeticError:
        return right.absorb(left)


class AsymptoticExpression(sage.rings.ring_element.RingElement):
    r"""
    ...
    """
    def __init__(self, parent, poset):
        r"""
        ...
        """
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
                       self._poset_.shells_topological(include_special=False,
                                                       reverse=reverse))
        return s

    # todo: implement _le_?

    def _update_poset_(self):
        r"""
        Helper method. This is called by ring operations like
        :meth:`_add_` and :meth:`_mul_`. This triggers the absorption
        of elements that are already in the poset.

        INPUT:

        Nothing.

        OUTPUT:

        Nothing.

        ...
        """
        for shell in self.poset.shells_topological(reverse=True):
            from sage.monoids.asymptotic_term_monoid import OTerm
            if isinstance(shell.element, OTerm) and shell.element in self.poset:
                while self.poset.null not in shell.predecessors():
                    for elem in shell.predecessors():
                        self.poset.remove(elem)


    def _add_(self, other):
        r"""
        Add ``other`` to this asymptotic expression.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpression`.

        OUTPUT:

        An :class:`AsymptoticExpression`.
        """
        poset = self.poset.union(other.poset)
        self._update_poset_()
        return self.parent()(poset=poset)


    def _sub_(self, other):
        r"""
        Subtract ``other`` from this asymptotic expression.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpression`.

        OUTPUT:

        An :class:`AsymptoticExpression`.
        """
        pass


    def _mul_(self, other):
        r"""
        Multiply ``other`` to this asymptotic expression.

        INPUT:

        - ``other`` -- an :class:`AsymptoticExpression`.

        OUTPUT:

        An :class:`AsymptoticExpression`.
        """
        pass



class AsymptoticRing(sage.rings.ring.Ring):
    pass