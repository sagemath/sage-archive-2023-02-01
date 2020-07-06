# ****************************************************************************
#       Copyright (C) 2020 Stefan Grosser <stefan.grosser1@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.posets.posets import Poset, FinitePoset
from sage.misc.lazy_attribute import lazy_attribute
from .linear_extensions import LinearExtensionsOfForest

class ForestPoset(FinitePoset):
    r"""
    A forest poset is a poset where the underlying Hasse diagram and is 
    directed acyclic graph. 

    """
                
    def __init__(self, hasse_diagram, elements, category, facade, key):
        FinitePoset.__init__(self, hasse_diagram=hasse_diagram, elements=elements, category=category, facade=facade, key=key)
        self._lin_ext_type = LinearExtensionsOfForest
        self._desc = 'Finite forest poset'
