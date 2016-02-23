# -*- coding: utf-8 -*-
r"""
Crystal Of Mirković-Vilonen (MV) Polytopes
"""

#*****************************************************************************
#       Copyright (C) 2015 Dinakar Muthiah <your email>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType


class MVPolytope(Element):
    def __init__(self, parent, long_word, lusztig_datum):
        Element.__init__(self, parent)
        self._initial_long_word = tuple(long_word)
        self._lusztig_datum = tuple(lusztig_datum)
        self._lusztig_data_dict = {self._initial_long_word: self._lusztig_datum}

class MVPolytopes(UniqueRepresentation, Parent):
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        """
        Normalize input to ensure a unique representation.
        """
        return super(MVPolytopes, cls).__classcall__(cls, CartanType(cartan_type))

    def __init__(self, cartan_type):
        """
        Initialize ``self``.
        """
        self._cartan_type = cartan_type
        Parent.__init__(self, category=(HighestWeightCrystals(), InfiniteEnumeratedSets()))
        if not cartan_type.is_finite():
            raise NotImplementedError("only implemented for finite Cartan types")

        # There must be a better way to do the following
        i = self._cartan_type.index_set()[0]
        self._default_word = self._pbw_datum_parent._long_word_begin_with(i)
        zero_lusztig_datum = [0]*len(self._default_word)
        self.module_generators = (self.element_class(self, 
                                                     self._default_word,
                                                     zero_lusztig_datum),)

    Element = MVPolytope

