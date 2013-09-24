"""
Points of a Toric Variety
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************



from copy import copy

from sage.misc.misc import powerset
from sage.combinat.cartesian_product import CartesianProduct


class NaivePointIterator(object):
    
    def __init__(self, fan, ring):
        """
        The naive point iterator.

        This is very slow.
        """
        self.ring = ring
        self.fan = fan
        self.ker = fan.rays().matrix().integer_kernel().matrix()

    def cone_iter(self):
        """
        Iterate over all cones of the fan

        OUTPUT:

        Iterator over the cones, starting with the high-dimensional
        ones.

        EXAMPLES::

            sage: ni = toric_varieties.dP6().point_set().naive_iterator()
            sage: for cone in ni.cone_iter(): print cone.ambient_ray_indices()
            (0, 1)
            (1, 2)
            (2, 3)
            (3, 4)
            (4, 5)
            (0, 5)
            (0,)
            (1,)
            (2,)
            (3,)
            (4,)
            (5,)
            ()
        """
        fan = self.fan
        for d in range(fan.dim(), -1, -1):
            for cone in fan.cones(d):
                yield cone
            
    def coordinate_iter(self):
        """
        Iterate over all coordinates of points.

        EXAMPLES::

            sage: F2 = GF(2)
            sage: ni = toric_varieties.P2(base_ring=F2).point_set().naive_iterator()
            sage: list(ni.coordinate_iter())
            [[0, 0, 1], [1, 0, 0], [0, 1, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0], [1, 1, 1]]

            sage: ni = toric_varieties.P2(base_ring=F2).point_set().naive_iterator()
            sage: list(ni.coordinate_iter())
        """
        units = [x for x in self.ring if x != 0]
        zero = self.ring.zero()
        big_torus = [units] * self.fan.nrays()
        for cone in self.cone_iter():
            patch = copy(big_torus)
            for i in cone.ambient_ray_indices():
                patch[i] = [zero]
            for p in CartesianProduct(*patch):
                yield p

    def is_finite(self):
        return self.ring.is_finite()


    

