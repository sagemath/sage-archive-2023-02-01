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
from sage.categories.sets_cat import Sets
from sage.rings.integer import Integer

from .set import Set_base, Set_add_sub_operators, Set_boolean_operators

class ImageSubobject(Parent):

    def __init__(self, map, domain_subset, *, category=None, is_injective=None):

        if not is_Parent(domain_subset):
            from sage.sets.set import Set
            domain_subset = Set(domain_subset)

        if is_Map(map):
            map_category = map.category_for()
        else:
            map_category = Sets()

        if category is None:
            category = map_category._meet_(domain_subset.category())

        category = category.Subobjects()
        if domain_subset.is_finite():
            category = category.Finite()
        elif map.is_injective() and domain_subset.is_infinite():
            category = category.Infinite()

        Parent.__init__(self, category=category)

        self._map = map
        self._domain_subset = domain_subset

    def ambient(self):
        return self._map.codomain()

    def __repr__(self) -> str:
        r"""
        TESTS::

            sage: Partitions(3).map(attrcall('conjugate'))
            Image of Partitions of the integer 3 by *.conjugate()

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
        if self._map.is_injective():
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
        if self._map.is_injective():
            for x in self._domain_subset:
                yield self.f(x)
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
        return self._map(self._domain_subset.an_element())


class ImageSet(ImageSubobject, Set_base, Set_add_sub_operators, Set_boolean_operators):

    pass
