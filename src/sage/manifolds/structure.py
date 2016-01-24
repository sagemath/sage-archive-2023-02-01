r"""
Manifold Structures

These classes encode the structure of a manifold.

AUTHORS:

- Travis Scrimshaw (2015-11-25): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.fast_methods import Singleton
from sage.manifolds.chart import Chart, RealChart

# This is a slight abuse by making this a Singleton, but there is no
#    need to have different copies of this object.
class TopologicalStructure(Singleton):
    """
    The structure of a topological manifold over a general topological field.
    """
    chart = Chart
    name = "topological"

    def subcategory(self, cat):
        """
        Return the subcategory of ``cat`` corresponding to the structure
        of ``self``.

        EXAMPLES::

            sage: from sage.manifolds.structure import TopologicalStructure
            sage: from sage.categories.manifolds import Manifolds
            sage: TopologicalStructure().subcategory(Manifolds(RR))
            Category of manifolds over Real Field with 53 bits of precision

        """
        return cat

class RealTopologicalStructure(Singleton):
    """
    The structure of a topological manifold over `\RR`.
    """
    chart = RealChart
    name = "topological"

    def subcategory(self, cat):
        """
        Return the subcategory of ``cat`` corresponding to the structure
        of ``self``.

        EXAMPLES::

            sage: from sage.manifolds.structure import RealTopologicalStructure
            sage: from sage.categories.manifolds import Manifolds
            sage: RealTopologicalStructure().subcategory(Manifolds(RR))
            Category of manifolds over Real Field with 53 bits of precision

        """
        return cat

