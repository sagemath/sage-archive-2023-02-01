"""
Points of a Toric Variety

The classes here are not meant to be instatiated manually. Instead,
you should always use the methods of the
:class:`sage.schemes.toric.homset.SchemeHomset_points_toric_field <point set>` of the variety.
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

from sage.combinat.cartesian_product import CartesianProduct
from sage.misc.misc import powerset, prod
from sage.misc.cachefunc import cached_method


class NaivePointIterator(object):
    
    def __init__(self, point_homset, fan, ring):
        """
        The naive point iterator.

        This is very slow.
        """
        self.ring = ring
        self.fan = fan
        self.ker = fan.rays().matrix().integer_kernel().matrix()
        self.point_homset = point_homset

    @cached_method
    def rescalings(self):
        """
        EXAMPLES::

            sage: ni = toric_varieties.P2_123(base_ring=GF(5)).point_set().naive_iterator()
            sage: ni.rescalings()
            ((1, 1, 1), (4, 3, 2), (4, 2, 3), (1, 4, 4))

            sage: ni = toric_varieties.dP8(base_ring=GF(3)).point_set().naive_iterator()
            sage: ni.rescalings()
            ((1, 1, 1, 1), (1, 2, 2, 2), (2, 1, 2, 1), (2, 2, 1, 2))

            sage: ni = toric_varieties.P1xP1(base_ring=GF(3)).point_set().naive_iterator()
            sage: ni.rescalings()
            ((1, 1, 1, 1), (1, 1, 2, 2), (2, 2, 1, 1), (2, 2, 2, 2))
        """
        units = [x for x in self.ring if x != 0]
        result = []
        ker = self.ker
        for phases in CartesianProduct(*([units] * ker.nrows())):
            phases = tuple(prod(mu**exponent for mu, exponent in zip(phases, column))
                           for column in ker.columns())
            result.append(phases)
        return tuple(result)                                   

    def orbit(self, point):
        """
            sage: ni = toric_varieties.P2_123(base_ring=GF(7)).point_set().naive_iterator()
            sage: ni.orbit([1, 0, 0])
            set([(1, 0, 0), (2, 0, 0), (4, 0, 0)])
            sage: ni.orbit([0, 1, 0])
            set([(0, 1, 0), (0, 6, 0)])
            sage: ni.orbit([0, 0, 1])
            set([(0, 0, 2), (0, 0, 6), (0, 0, 1), (0, 0, 5), (0, 0, 4), (0, 0, 3)])
            sage: ni.orbit([1, 1, 0])
            set([(1, 1, 0), (2, 1, 0), (4, 1, 0), (2, 6, 0), (1, 6, 0), (4, 6, 0)])
        """
        result = set()
        for phases in self.rescalings():
            p = tuple(mu*z for mu, z in zip(point, phases))
            result.add(p)
        return result

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
        Iterate over all distinct coordinates ().

        EXAMPLES::

            sage: F2 = GF(2)
            sage: ni = toric_varieties.P2(base_ring=F2).point_set().naive_iterator()
            sage: list(ni.coordinate_iter())
            [(0, 0, 1), (1, 0, 0), (0, 1, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, 1)]

            sage: ni = toric_varieties.P1xP1(base_ring=F2).point_set().naive_iterator()
            sage: list(ni.coordinate_iter())
            [(0, 1, 0, 1), (1, 0, 0, 1), (1, 0, 1, 0),
             (0, 1, 1, 0), (0, 1, 1, 1), (1, 0, 1, 1),
             (1, 1, 0, 1), (1, 1, 1, 0), (1, 1, 1, 1)]
        """
        units = [x for x in self.ring if x != 0]
        zero = self.ring.zero()
        big_torus = [units] * self.fan.nrays()
        for cone in self.cone_iter():
            patch = copy(big_torus)
            for i in cone.ambient_ray_indices():
                patch[i] = [zero]
            for p in CartesianProduct(*patch):
                yield tuple(p)

    def point_iter(self):
        """
        EXAMPLES:

            sage: ni = toric_varieties.P2(base_ring=GF(2)).point_set().naive_iterator()
            sage: list(ni.point_iter())
            [[0 : 0 : 1], [1 : 0 : 0], [0 : 1 : 0], [0 : 1 : 1], [1 : 0 : 1], [1 : 1 : 0], [1 : 1 : 1]]

            sage: ni = toric_varieties.P1xP1(base_ring=GF(3)).point_set().naive_iterator()
            sage: list(ni.point_iter()) 
            [[0 : 1 : 0 : 1], [1 : 0 : 0 : 1], [1 : 0 : 1 : 0], [0 : 1 : 1 : 0],
             [0 : 1 : 1 : 1], [0 : 1 : 1 : 2], [1 : 0 : 1 : 1], [1 : 0 : 1 : 2],
             [1 : 1 : 0 : 1], [1 : 2 : 0 : 1], [1 : 1 : 1 : 0], [1 : 2 : 1 : 0],
             [1 : 1 : 1 : 1], [1 : 1 : 1 : 2], [1 : 2 : 1 : 1], [1 : 2 : 1 : 2]]
        """
        seen = set()
        for p in self.coordinate_iter():
            if p in seen:
                continue
            seen.update(self.orbit(p))
            yield self.point_homset(p)

    def is_finite(self):
        """
        Return whether there are finitely many points.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: P2.point_set().naive_iterator().is_finite()
            False
            sage: P2.change_ring(GF(7)).point_set().naive_iterator().is_finite()
            True
        """
        return self.ring.is_finite()


    

