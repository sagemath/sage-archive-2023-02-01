r"""
Groups of imaginary elements

.. NOTE::

    One main purpose of such groups is in an
    :doc:`asymptotic ring's <asymptotic_ring>`
    :doc:`growth group <growth_group>` when an element like `n^z`
    (for some constant `z`) is split into
    `n^{\Re z + I\Im z}`.
    (Note that the first summand in the exponent determines the growth,
    the second does not influence the growth.)

AUTHORS:

- Daniel Krenn (2018)

Classes and Methods
===================
"""
#*****************************************************************************
# Copyright (C) 2018 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.structure.element import AdditiveGroupElement
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp_by_eq_and_lt
from sage.structure.unique_representation import UniqueRepresentation


class ImaginaryElement(AdditiveGroupElement):
    r"""
    An element of :class:`ImaginaryGroup`.

    INPUT:

    - ``parent`` -- a SageMath parent

    - ``imag`` -- an element of parent's base
    """

    def __init__(self, parent, imag):
        r"""
        See :class:`ImaginaryGroup` for more information.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=-42)  # indirect doctest
            -42*I
        """
        if parent is None:
            raise ValueError('parent must be provided')
        super(ImaginaryElement, self).__init__(parent=parent)

        try:
            self._imag_ = parent.base()(imag)
        except (TypeError, ValueError) as e:
            from sage.rings.asymptotic.misc import combine_exceptions
            from sage.structure.element import parent as parent_function
            raise combine_exceptions(
                ValueError(
                    '{} ({}) is not in {}'.format(imag,
                                                  parent_function(imag),
                                                  parent.base())),
                e)

    def imag(self):
        r"""
        Return the imaginary part of this imaginary element.

        EXAMPLES::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(I).imag()
            1
        """
        return self._imag_

    def real(self):
        r"""
        Return the real part (`=0`) of this imaginary element.

        EXAMPLES::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(I).real()
            0
        """
        return self.parent().base().zero()

    def __hash__(self):
        r"""
        Return a hash value of this imaginary element.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: hash(J(I))  # indirect doctest, random
            42
        """
        return hash((self.parent(), self._imag_))

    _richcmp_ = richcmp_by_eq_and_lt("_eq_", "_lt_")

    def _eq_(self, other):
        r"""
        Return whether this imaginary part equals ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=3) == J(imag=4)
            False
            sage: J(imag=42) == J(imag=42)
            True
        """
        return self._imag_ == other._imag_

    def _lt_(self, other):
        r"""
        Raise an error since imaginary elements are not comparable.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=1) < J(imag=1)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for imaginary elements I and I
            sage: J(imag=1) < J(imag=2)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for imaginary elements I and 2*I
            sage: J(imag=1) > J(imag=2)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for imaginary elements 2*I and I
        """
        raise RuntimeError("cannot decide '<' "
                           "for imaginary elements "
                           "{} and {}".format(self, other))

    def _repr_(self):
        r"""
        Return a representation string of this imaginary element.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=0)  # indirect doctest
            0
            sage: J(imag=1)  # indirect doctest
            I
            sage: J(imag=-1)  # indirect doctest
            -I
            sage: J(imag=42)  # indirect doctest
            42*I
            sage: J(imag=-42)  # indirect doctest
            -42*I
        """
        from sage.rings.asymptotic.misc import repr_op
        if self._imag_ == 0:
            return '0'
        if self._imag_ == 1:
            return 'I'
        if self._imag_ == -1:
            return '-I'
        return repr_op(self._imag_, '*', 'I')

    def _add_(self, other):
        r"""
        Return the sum of this imaginary element and ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=1) + J(imag=-2)  # indirect doctest
            -I
        """
        P = self.parent()
        return P.element_class(P, self._imag_ + other._imag_)

    def _sub_(self, other):
        r"""
        Return the difference of this imaginary element and ``other``.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: J(imag=1) - J(imag=-2)  # indirect doctest
            3*I
        """
        P = self.parent()
        return P.element_class(P, self._imag_ - other._imag_)

    def __neg__(self):
        r"""
        Return the negative of this imaginary element.

        TESTS::

            sage: from sage.groups.misc_gps.imaginary_groups import ImaginaryGroup
            sage: J = ImaginaryGroup(ZZ)
            sage: -J(imag=-2)  # indirect doctest
            2*I
        """
        P = self.parent()
        return P.element_class(P, -self._imag_)

