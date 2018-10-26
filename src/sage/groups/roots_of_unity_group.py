r"""
Groups of elements on the unit circle (e.g. roots of unity).

AUTHORS:

- Daniel Krenn (2018)

Classes and Methods
===================
"""
#*****************************************************************************
# Copyright (C) 2014--2015 Benjamin Hackl <benjamin.hackl@aau.at>
#               2014--2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.structure.element import MultiplicativeGroupElement
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.richcmp import richcmp_by_eq_and_lt


class UnitCirclePoint(MultiplicativeGroupElement):
    r"""
    An element of :class:`UnitCircleGroup`
    which is `e^{2\pi\mathit{exponent}}`.

    INPUT:

    - ``parent`` -- a SageMath parent

    - ``exponent`` -- a number (of a subset of the reals)

    - ``normalize`` -- a boolean
    """
    def __init__(self, parent, exponent, normalize=True):
        r"""
        See :class:`UnitCirclePoint` for more information.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=1/2)
            e^(2*pi*0.500000000000000)
            sage: C(exponent=3/2)
            e^(2*pi*0.500000000000000)

            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0)
            1
            sage: U(exponent=1)
            1
            sage: U(exponent=2/3) == U(exponent=5/3)
            True
        """
        if parent is None:
            raise ValueError('parent must be provided')
        super(UnitCirclePoint, self).__init__(parent=parent)

        exponent = parent.base()(exponent)
        if normalize:
            exponent = exponent - exponent.floor()
        self._exponent_ = exponent

    def _repr_(self):
        r"""
        Return a representation string of this point on the unit circle.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=1/3)
            e^(2*pi*0.333333333333333)
        """
        return 'e^(2*pi*{})'.format(self._exponent_)

    def __hash__(self):
        r"""
        Return a hash value of this point on the unit circle.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: C = UnitCircleGroup(RR)
            sage: hash(C(exponent=1/3))  # indirect doctest, random
            42
        """
        return hash((self.parent(), self._exponent_))

    def _mul_(self, other):
        r"""
        Return the product of this point on the unit circle and ``other``.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=0.3) * C(exponent=0.4)
            e^(2*pi*0.700000000000000)
        """
        P = self.parent()
        return P.element_class(P, self._exponent_ + other._exponent_)

    def __invert__(self):
        r"""
        Return the inverse of this point on the unit circle.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: C = UnitCircleGroup(RR)
            sage: ~C(exponent=0.4)
            e^(2*pi*0.600000000000000)
            sage: C(1) / C(exponent=0.4)
            e^(2*pi*0.600000000000000)
        """
        P = self.parent()
        return P.element_class(P, -self._exponent_)

    _richcmp_ = richcmp_by_eq_and_lt("_eq_", "_lt_")

    def _eq_(self, other):
        r"""
        Return whether this point on the unit circle equals ``other``.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0) == U(exponent=1)
            True
            sage: U(exponent=2/3) == U(exponent=5/3)
            True
            sage: U(exponent=2/3) == U(exponent=-2/3)
            False
        """
        return self._exponent_ == other._exponent_

    def _lt_(self, other):
        r"""
        Raise an error since points on the unit circle are not comparable.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0) < U(exponent=0)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for the roots of unity 1 and 1
            sage: U(exponent=0) < U(exponent=1/2)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for the roots of unity 1 and -1
            sage: U(exponent=0) > U(exponent=1/2)  # indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: cannot decide '<' for the roots of unity -1 and 1
        """
        raise RuntimeError("cannot decide '<' "
                           "for the roots of unity "
                           "{} and {}".format(self, other))

    def is_one(self):
        r"""
        Return whether this point on the unit circle is `1`.

        EXAMPLES:

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: C = UnitCircleGroup(QQ)
            sage: C(exponent=0).is_one()
            True
            sage: C(exponent=1/2).is_one()
            False
            sage: C(exponent=2/3).is_one()
            False
        """
        return self._exponent_ == 0

    def is_minus_one(self):
        r"""
        Return whether this point on the unit circle is `-1`.

        EXAMPLES:

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: C = UnitCircleGroup(QQ)
            sage: C(exponent=0).is_minus_one()
            False
            sage: C(exponent=1/2).is_minus_one()
            True
            sage: C(exponent=2/3).is_minus_one()
            False
        """
        from sage.rings.rational_field import QQ
        return self._exponent_ == QQ(1)/QQ(2)


class UnitCircleGroup(UniqueRepresentation, Parent):
    r"""
    A group for points on the unit circle.

    INPUT:

    - ``base`` -- a SageMath parent representing a subset of the reals

    - ``category`` -- a category

    EXAMPLES::

        sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
        sage: UnitCircleGroup(RR)
        Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
        sage: UnitCircleGroup(QQ)
        Unit Circle Group with Exponents in Rational Field modulo ZZ
    """

    Element = UnitCirclePoint

    @staticmethod
    def __classcall__(cls, base, category=None):
        r"""
        See :class:`UnitCircleGroup` for more information.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup(RR).category()  # indirect doctest
            Category of commutative groups
        """
        category = cls._determine_category_(category)
        return super(UnitCircleGroup, cls).__classcall__(
            cls, base, category)

    @staticmethod
    def _determine_category_(category):
        r"""
        Return the category of this unit circle group.

        INPUT:

        - ``category`` -- a category or ``None`` (in which case the output
          equals ``category``)

        OUTPUT:

        A category.

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup._determine_category_(None)
            Category of commutative groups
            sage: UnitCircleGroup._determine_category_(Groups())
            Category of groups
        """
        if category is None:
            from sage.categories.groups import Groups
            category = Groups().Commutative()
        return category

    def __init__(self, base, category):
        r"""
        See :class:`UnitCircleGroup` for more information.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup(RR).base()  # indirect doctest
            Real Field with 53 bits of precision
        """
        super(UnitCircleGroup, self).__init__(category=category,
                                              base=base)

    def _repr_(self):
        r"""
        Return a representation string of this unit circle group.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup(RR)  # indirect doctest
            Unit Circle Group with Exponents in Real Field with 53 bits of precision modulo ZZ
        """
        return 'Unit Circle Group with Exponents in {} modulo ZZ'.format(self.base())

    def _repr_short_(self):
        r"""
        Return a short representation string of this unit circle group.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup(RR)._repr_short_()
            'U_RR'
        """
        from sage.rings.asymptotic.misc import parent_to_repr_short
        s = parent_to_repr_short(self.base())
        if ' ' in s:
            s = '({})'.format(s)
        return 'U_{}'.format(s)

    def __hash__(self):
        r"""
        Return a hash value of this unit circle group.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: hash(UnitCircleGroup(RR))  # indirect doctest, random
            42
        """
        return hash((self.__class__, self.base()))

    def _an_element_(self):
        r"""
        Return an element of this unit circle group.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup
            sage: UnitCircleGroup(RR).an_element()  # indirect doctest
            e^(2*pi*0.000000000000000)
        """
        return self.element_class(self, self.base().an_element())

    def _element_constructor_(self, data, exponent=None):
        r"""
        Construct an element out of the given data.

        INPUT:

        - ``data`` -- an object

        - ``exponent`` -- a number (of a subset of the reals) or ``None``

        OUTPUT:

        A :class:`UnitCirclePoint`.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
            sage: C = UnitCircleGroup(RR)
            sage: C(exponent=1/2)
            e^(2*pi*0.500000000000000)

            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0)
            1
            sage: U(exponent=1)
            1
            sage: U(exponent=1/2)
            -1
            sage: U(exponent=1/4)
            I
            sage: U(exponent=1/3)
            zeta3

            sage: C.<z> = CyclotomicField(6)
            sage: z, U(z)
            (z, zeta6)
            sage: z^2, U(z^2)
            (z - 1, zeta3)

            sage: U(ZZ(-1))
            -1
            sage: U(QQ(-1))
            -1
            sage: U(int(-1))
            -1
        """
        from sage.groups.generic import discrete_log
        from sage.rings.asymptotic.misc import combine_exceptions
        from sage.rings.rational_field import QQ
        from sage.rings.number_field.number_field import NumberField_cyclotomic

        if exponent is None:
            if isinstance(data, int) and data == 0:
                raise ValueError('no input specified')

            elif isinstance(data, self.element_class):
                if data.parent() == self:
                    return data
                exponent = data._exponent_

            elif data == 1 or data == '1':
                exponent = 0

            elif data == -1 or data == '-1':
                exponent = QQ(1)/QQ(2)

            else:
                try:
                    P = data.parent()
                except AttributeError:
                    raise TypeError('{} is not in {}'.format(data, self))

                if isinstance(P, NumberField_cyclotomic):
                    zeta = P.gen()
                    n = zeta.multiplicative_order()
                    try:
                        exponent = QQ(discrete_log(data, zeta)) / QQ(n)
                    except ValueError as e:
                        raise combine_exceptions(
                            ValueError('{} is not in {}'.format(data, self)), e)

            if exponent is None:
                raise ValueError('{} is not in {}'.format(data, self))

        elif not isinstance(data, int) or data != 0:
            raise ValueError('input is ambigous: '
                             '{} as well as exponent={} '
                             'specified'.format(data, exponent))

        return self.element_class(self, exponent)


class RootOfUnity(UnitCirclePoint):
    r"""
    A root of unity (i.e. an element of :class:`RootsOfUnityGroup`)
    which is `e^{2\pi\mathit{exponent}}` for a rational ``exponent``.
    """

    def exponent_numerator(self):
        r"""
        Return the numerator of the rational quotient in `[0,1)`
        representing the exponent of this root of unity.

        EXAMPLES::

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: a = U(exponent=2/3); a
            zeta3^2
            sage: a.exponent_numerator()
            2
        """
        return self._exponent_.numerator()

    def exponent_denominator(self):
        r"""
        Return the denominator of the rational quotient in `[0,1)`
        representing the exponent of this root of unity.

        EXAMPLES::

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: a = U(exponent=2/3); a
            zeta3^2
            sage: a.exponent_denominator()
            3
        """
        return self._exponent_.denominator()

    def _repr_(self):
        r"""
        Return a representation string of this root of unity.

        TESTS::

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: U = RootsOfUnityGroup()
            sage: U(exponent=0)
            1
            sage: U(exponent=1/2)
            -1
            sage: U(exponent=1/4)
            I
            sage: U(exponent=3/4)
            -I
            sage: U(exponent=1/3)
            zeta3
            sage: U(exponent=2/3)
            zeta3^2
        """
        from sage.rings.rational_field import QQ
        if self._exponent_ == 0:
            return '1'
        if self._exponent_ == QQ(1)/QQ(2):
            return '-1'
        if self._exponent_ == QQ(1)/QQ(4):
            return 'I'
        if self._exponent_ == QQ(3)/QQ(4):
            return '-I'
        num = self.exponent_numerator()
        den = self.exponent_denominator()
        zeta = 'zeta{}'.format(den)
        if num == 1:
            return zeta
        return '{}^{}'.format(zeta, num)


class RootsOfUnityGroup(UnitCircleGroup):
    r"""
    The group of all roots of unity.

    INPUT:

    - ``category`` -- a category

    This is a specialized :class:`UnitCircleGroup` with base `\QQ`.
    """

    Element = RootOfUnity

    @staticmethod
    def __classcall__(cls, category=None):
        r"""
        See :class:`RootsOfUnityGroup` for more information.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: RootsOfUnityGroup().category()  # indirect doctest
            Category of commutative groups
        """
        category = cls._determine_category_(category)
        return super(UnitCircleGroup, cls).__classcall__(
            cls, category)

    def __init__(self, category):
        r"""
        See :class:`RootsOfUnityGroup` for more information.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: RootsOfUnityGroup().base()  # indirect doctest
            Rational Field
        """
        from sage.rings.rational_field import QQ
        return super(RootsOfUnityGroup, self).__init__(base=QQ,
                                                         category=category)
    def _repr_(self):
        r"""
        Return a representation string of this roots of unity group.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: RootsOfUnityGroup()  # indirect doctest
            Group of Roots of Unity
        """
        return 'Group of Roots of Unity'

    def _repr_short_(self):
        r"""
        Return a short representation string of this roots of unity group.

        TESTS:

            sage: from sage.groups.roots_of_unity_group import RootsOfUnityGroup
            sage: RootsOfUnityGroup()._repr_short_()
            'U'
        """
        return 'U'


