r"""
Relative Interiors of Polyhedra and Cones
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject

class RelativeInterior(SageObject):

    """
    The relative interior of a polyhedron or cone
    """

    def __init__(self, polyhedron):
        self._polyhedron = polyhedron

    def __contains__(self, point):
        return self._polyhedron.relative_interior_contains(point)

    def closure(self):
        return self._polyhedron

    def _repr_(self):
        """
        Return a description of ``self``.

        EXAMPLES::

            sage: P = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1]])
            sage: P.relative_interior()._repr_()
            'Relative interior of a 2-dimensional polyhedron in ZZ^4 defined as the convex hull of 3 vertices'
            sage: P.rename('A')
            sage: P.relative_interior()._repr_()
            'Relative interior of A'
        """
        repr_P = repr(self._polyhedron)
        if repr_P.startswith('A '):
            repr_P = 'a ' + repr_P[2:]
        return 'Relative interior of ' + repr_P
