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

class ForestPoset(FinitePoset):
    r"""
    A forest poset is a poset where the underlying Hasse diagram and is 
    directed acyclic graph. 

    """

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.combinat.posets.forest import ForestPoset
            sage: P = ForestPoset(DiGraph({0: [2], 1: [2], 2: [3, 4], 3: [], 4: []}))
            sage: P._repr_()
            'Finite forest poset containing 5 elements'
        """
        s = "Finite forest poset containing %s elements" % self._hasse_diagram.order()
        if self._with_linear_extension:
            s += " with distinguished linear extension"
        return s
                
    def __init__(self, hasse_diagram, elements, category, facade, key):
        FinitePoset.__init__(self, hasse_diagram=hasse_diagram, elements=elements, category=category, facade=facade, key=key)


    def linear_extensions(self, facade=False):
        r"""
        Returns the enumerated set of all the linear extensions of this forest poset

        INPUT:

        - ``facade`` -- a boolean (default: ``False``);
          whether to return the linear extensions as plain lists

        EXAMPLES::

        TESTS::
        """
        from .linear_extensions import LinearExtensionsOfForest
        return LinearExtensionsOfForest(self, facade=facade)