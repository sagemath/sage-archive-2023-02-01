"""
Image Sets
"""

# ****************************************************************************
#       Copyright (C) 2008      Mike Hansen <mhansen@gmail.com>
#                     2012      Christian Stump
#                     2020-2021 Frédéric Chapoton
#                     2021      Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from typing import Iterator

from sage.structure.parent import Parent, is_Parent
from sage.categories.map import is_Map
from sage.categories.poor_man_map import PoorManMap
from sage.categories.sets_cat import Sets
from sage.categories.enumerated_sets import EnumeratedSets
from sage.rings.integer import Integer
from sage.modules.free_module import FreeModule
from sage.symbolic.callable import is_CallableSymbolicExpression

from .set import Set_base, Set_add_sub_operators, Set_boolean_operators

class ImageSubobject(Parent):

    def __init__(self, map, domain_subset, *, category=None, is_injective=None):

        if not is_Parent(domain_subset):
            from sage.sets.set import Set
            domain_subset = Set(domain_subset)

        if not is_Map(map) and not isinstance(map, PoorManMap):
            map_name = f"The map {map}"
            if is_CallableSymbolicExpression(map):
                domain = map.parent().base()
                if len(map.arguments()) != 1:
                    domain = FreeModule(domain, len(map.arguments()))
                map = self._star(map)
            else:
                domain = domain_subset
            map = PoorManMap(map, domain, name=map_name)

        if is_Map(map):
            map_category = map.category_for()
            if is_injective is None:
                is_injective = map.is_injective()
        else:
            map_category = Sets()
            if is_injective is None:
                is_injective = False

        if category is None:
            category = map_category._meet_(domain_subset.category())

        category = category.Subobjects()
        if domain_subset in Sets().Finite():
            category = category.Finite()
        elif is_injective and domain_subset in Sets.Infinite():
            category = category.Infinite()

        if domain_subset in EnumeratedSets():
            category = category & EnumeratedSets()

        Parent.__init__(self, category=category)

        self._map = map
        self._domain_subset = domain_subset
        self._is_injective = is_injective

    @staticmethod
    def _star(function):
        def f(arg):
            return function(*arg)
        return f

    def ambient(self):
        return self._map.codomain()

    def _repr_(self) -> str:
        r"""
        TESTS::

            sage: Partitions(3).map(attrcall('conjugate'))
            Image of Partitions of the integer 3 by The map *.conjugate() from Partitions of the integer 3

        """
        return f"Image of {self._domain_subset} by {self._map}"

    def cardinality(self) -> Integer:
        r"""
        Return the cardinality of this combinatorial class

        EXAMPLES::

            sage: R = Permutations(10).map(attrcall('reduced_word'))
            sage: R.cardinality()
            3628800

        """
        if self._is_injective:
            return self._domain_subset.cardinality()
        return super().cardinality()

    def __iter__(self) -> Iterator:
        r"""
        Return an iterator over the elements of this combinatorial class

        EXAMPLES::

            sage: R = Permutations(10).map(attrcall('reduced_word'))
            sage: R.cardinality()
            3628800
        """
        if self._is_injective:
            for x in self._domain_subset:
                yield self._map(x)
        else:
            yield from super().__iter__()

    def _an_element_(self):
        r"""
        Return an element of this combinatorial class

        EXAMPLES::

            sage: R = SymmetricGroup(10).map(attrcall('reduced_word'))
            sage: R.an_element()
            [9, 8, 7, 6, 5, 4, 3, 2]
        """
        domain_element = self._domain_subset.an_element()
        return self._map(domain_element)

    def _sympy_(self):
        r"""
        Return an instance of a subclass of SymPy ``Set`` corresponding to ``self``.

        EXAMPLES::

            sage: from sage.sets.image_set import ImageSet
            sage: S = ImageSet(sin, RealSet.open(0, pi/4)); S
            Image of (0, 1/4*pi) by The map sin from (0, 1/4*pi)
            sage: S._sympy_()
            ImageSet(Lambda(x, sin(x)), Interval.open(0, pi/4))
        """
        from sympy import imageset
        try:
            sympy_map = self._map._sympy_()
        except AttributeError:
            sympy_map = self._map
        return imageset(sympy_map,
                        self._domain_subset._sympy_())


class ImageSet(ImageSubobject, Set_base, Set_add_sub_operators, Set_boolean_operators):

    r"""
    Image of a set by a map

    EXAMPLES:

        sage: from sage.sets.image_set import ImageSet

    Symbolics::

        sage: ImageSet(sin, RealSet.open(0, pi/4))
        Image of (0, 1/4*pi) by The map sin from (0, 1/4*pi)
        sage: _.an_element()
        1/2*sqrt(-sqrt(2) + 2)

        sage: sos(x,y) = x^2 + y^2; sos
        (x, y) |--> x^2 + y^2
        sage: ImageSet(sos, ZZ^2)
        Image of
         Ambient free module of rank 2 over the principal ideal domain Integer Ring by
         The map (x, y) |--> x^2 + y^2 from Vector space of dimension 2 over Symbolic Ring
        sage: _.an_element()
        1
        sage: ImageSet(sos, Set([(3, 4), (3, -4)]))
        Image of
         {(3, -4), (3, 4)} by
         The map (x, y) |--> x^2 + y^2 from Vector space of dimension 2 over Symbolic Ring
        sage: _.an_element()
        25

    """

    pass
