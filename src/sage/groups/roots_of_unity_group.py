from __future__ import absolute_import

from sage.structure.element import MultiplicativeGroupElement
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.richcmp import richcmp_by_eq_and_lt


class UnitCirclePoint(MultiplicativeGroupElement):
    def __init__(self, parent, exponent, normalize=True):
        r"""
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
        return 'e^(2*pi*{})'.format(self._exponent_)

    def __hash__(self):
        return hash((self.parent(), self._exponent_))

    def _mul_(self, other):
        P = self.parent()
        return P.element_class(P, self._exponent_ + other._exponent_)

    def __invert__(self):
        P = self.parent()
        return P.element_class(P, -self._exponent_)

    _richcmp_ = richcmp_by_eq_and_lt("_eq_", "_lt_")

    def _eq_(self, other):
        r"""
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
        return RuntimeError("cannot decide '<' "
                            "for the roots of unity "
                            "{} and {}".format(self, other))


class UnitCircleGroup(UniqueRepresentation, Parent):

    Element = UnitCirclePoint

    @staticmethod
    def __classcall__(cls, base, category=None):
        category = cls._determine_category_(category)
        return super(UnitCircleGroup, cls).__classcall__(
            cls, base, category)

    @staticmethod
    def _determine_category_(category):
        if category is None:
            from sage.categories.groups import Groups
            category = Groups().Commutative()
        return category

    def __init__(self, base, category):
        super(UnitCircleGroup, self).__init__(category=category,
                                              base=base)

    def _repr_(self):
        return 'Unit Circle Group with Exponents in {} modulo ZZ'.format(self.base())

    def __hash__(self):
        return hash((self.__class__, self.base()))

    def _an_element_(self):
        return self.element_class(self, self.base().an_element())



class RootOfUnity(UnitCirclePoint):

    def exponent_numerator(self):
        return self._exponent_.numerator()

    def exponent_denominator(self):
        return self._exponent_.denominator()

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.groups.roots_of_unity_group import UnitCircleGroup, RootsOfUnityGroup
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

    Element = RootOfUnity

    @staticmethod
    def __classcall__(cls, category=None):
        category = cls._determine_category_(category)
        return super(UnitCircleGroup, cls).__classcall__(
            cls, category)

    def __init__(self, category):
        from sage.rings.rational_field import QQ
        return super(RootsOfUnityGroup, self).__init__(base=QQ,
                                                         category=category)
    def _repr_(self):
        return 'Group of Roots of Unity'
