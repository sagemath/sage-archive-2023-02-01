# -*- coding: utf-8 -*-
"""
Enumerate Points of a Toric Variety

The classes here are not meant to be instatiated manually. Instead,
you should always use the methods of the :class:`point set
<sage.schemes.toric.homset.SchemeHomset_points_toric_field>` of the
variety.

In this module, points are always represented by tuples instead of
Sage's class for points of the toric variety. All Sage library code
must then convert it to proper point objects before returning it to
the user.

EXAMPLES::

    sage: P2 = toric_varieties.P2(base_ring=GF(3))
    sage: point_set = P2.point_set()
    sage: point_set.cardinality()
    13
    sage: iter(point_set).next()
    [0 : 0 : 1]
    sage: list(point_set)[0:5]
    [[0 : 0 : 1], [1 : 0 : 0], [0 : 1 : 0], [0 : 1 : 1], [0 : 1 : 2]]
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************



from copy import copy

from sage.combinat.cartesian_product import CartesianProduct
from sage.misc.misc import powerset, prod
from sage.misc.cachefunc import cached_method
from sage.rings.arith import gcd


class InfinitePointEnumerator(object):
    
    def __init__(self, fan, ring):
        """
        Point enumerator for infinite fields.

        INPUT:
        
        - ``fan`` -- fan of the toric variety.

        - ``ring`` -- infinite base ring over which to enumerate
          points.

        TESTS::

            sage: from sage.schemes.toric.points import InfinitePointEnumerator
            sage: fan = toric_varieties.P2().fan()
            sage: n = InfinitePointEnumerator(fan, QQ)
            sage: ni = iter(n)
            sage: [ni.next() for k in range(10)]
            [(0, 1, 1), (1, 1, 1), (-1, 1, 1), (1/2, 1, 1), (-1/2, 1, 1), 
             (2, 1, 1), (-2, 1, 1), (1/3, 1, 1), (-1/3, 1, 1), (3, 1, 1)]
          
            sage: X = ToricVariety(Fan([], lattice=ZZ^0))
            sage: X.point_set().cardinality()
            1
            sage: X.base_ring().is_finite()
            False
            sage: X.point_set().list()
            ([],)
        """
        self.ring = ring
        self.fan = fan

    def __iter__(self):
        """
        Iterate over the points.

        OUTPUT:

        Iterator over points.

        EXAMPLES::

            sage: from sage.schemes.toric.points import InfinitePointEnumerator
            sage: fan = toric_varieties.P2().fan()
            sage: n = InfinitePointEnumerator(fan, QQ)
            sage: ni = iter(n)
            sage: [ni.next() for k in range(5)]
            [(0, 1, 1), (1, 1, 1), (-1, 1, 1), (1/2, 1, 1), (-1/2, 1, 1)]
        """
        rays = self.fan().rays() + self.fan().virtual_rays()
        n = len(rays)
        if n == 0:
            yield tuple()
        else:
            R = self.ring
            p = [R.one() for k in range(n)]
            for r in R:
                p[0] = r
                yield tuple(p)


class NaiveFinitePointEnumerator(object):
    
    def __init__(self, fan, ring):
        """
        The naive point enumerator.

        This is very slow.
        
        INPUT:
        
        - ``fan`` -- fan of the toric variety.

        - ``ring`` -- finite base ring over which to enumerate points.

        EXAMPLES::

            sage: from sage.schemes.toric.points import NaiveFinitePointEnumerator
            sage: fan = toric_varieties.P2().fan()
            sage: n = NaiveFinitePointEnumerator(fan, GF(3))
            sage: iter(n).next()
            (0, 0, 1)
        """
        assert ring.is_finite()
        self.ring = ring
        self.fan = fan

    @cached_method
    def rays(self):
        """
        Return all rays (real and virtual).

        OUTPUT:

        Tuple of rays of the fan.

        EXAMPLES::

            sage: from sage.schemes.toric.points import NaiveFinitePointEnumerator
            sage: fan = toric_varieties.torus(2).fan()
            sage: fan.rays()
            Empty collection
            in 2-d lattice N
            sage: n = NaiveFinitePointEnumerator(fan, GF(3))
            sage: n.rays()
            N(1, 0),
            N(0, 1)
            in 2-d lattice N
        """
        return self.fan.rays() + self.fan.virtual_rays()

    @cached_method
    def units(self):
        """
        Return the units in the base field.

        EXAMPLES::

            sage: ne = toric_varieties.P2(base_ring=GF(5)).point_set()._naive_enumerator()
            sage: ne.units()
            (1, 2, 3, 4)
        """
        return tuple(x for x in self.ring if x != 0)
        
    @cached_method
    def roots(self, n):
        """
        Return the n-th roots in the base field

        INPUT:

        - ``n`` integer.

        OUTPUT:

        Tuple containing all n-th roots (not only the primitive
        ones). In particular, 1 is included.

        EXAMPLES::

            sage: ne = toric_varieties.P2(base_ring=GF(5)).point_set()._naive_enumerator()
            sage: ne.roots(2)
            (1, 4)
            sage: ne.roots(3)
            (1,)
            sage: ne.roots(4)
            (1, 2, 3, 4)
        """
        return tuple(x for x in self.ring if x**n == self.ring.one())
        
    def _Chow_group_free(self):
        r"""
        Return the relations coming from the free part of the Chow group

        OUTPUT:

        A tuple containing the elements of $Hom(A_{d-1,\text{free}},
        F^\times)$, including the identity.

        EXAMPLES::

            sage: fan = NormalFan(ReflexivePolytope(2, 0))
            sage: X = ToricVariety(fan, base_field=GF(7))
            sage: X.Chow_group().degree(1)
            C3 x Z
            sage: enum = X.point_set()._naive_enumerator()
            sage: enum._Chow_group_free()
            ((1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4), (5, 5, 5), (6, 6, 6))
        """
        units = self.units()
        result = []
        ker = self.rays().matrix().integer_kernel().matrix()
        for phases in CartesianProduct(*([units] * ker.nrows())):
            phases = tuple(prod(mu**exponent for mu, exponent in zip(phases, column))
                           for column in ker.columns())
            result.append(phases)
        return tuple(sorted(result))

    def _Chow_group_torsion(self):
        r"""
        Return the relations coming from the torison part of the Chow group

        OUTPUT:

        A tuple containing the non-identity elements of
        $Hom(A_{d-1,\text{tors}}, F^\times)$

        EXAMPLES::

            sage: fan = NormalFan(ReflexivePolytope(2, 0))
            sage: X = ToricVariety(fan, base_field=GF(7))
            sage: X.Chow_group().degree(1)
            C3 x Z
            sage: enum = X.point_set()._naive_enumerator()
            sage: enum._Chow_group_torsion()
            ((1, 2, 4), (1, 4, 2))
        """
        if self.fan.is_smooth():
            return tuple()
        image = self.rays().column_matrix().image()
        torsion = image.saturation().quotient(image)
        result = set()
        for t in torsion:
            t_lift = t.lift()
            for root in self.roots(t.order()):
                phases = tuple(root**exponent for exponent in t_lift)
                result.add(phases)
        result.remove(tuple(self.ring.one() for r in self.rays()))
        return tuple(sorted(result))

    @cached_method
    def rescalings(self):
        """
        Return the rescalings of homogeneous coordinates.

        OUTPUT:

        A tuple containing all points that are equivalent to
        `[1:1:\dots:1]`, the distinguished point of the big torus
        orbit.
        
        EXAMPLES::

            sage: ni = toric_varieties.P2_123(base_ring=GF(5)).point_set()._naive_enumerator()
            sage: ni.rescalings()
            ((1, 1, 1), (1, 4, 4), (4, 2, 3), (4, 3, 2))

            sage: ni = toric_varieties.dP8(base_ring=GF(3)).point_set()._naive_enumerator()
            sage: ni.rescalings()
            ((1, 1, 1, 1), (1, 2, 2, 2), (2, 1, 2, 1), (2, 2, 1, 2))

            sage: ni = toric_varieties.P1xP1(base_ring=GF(3)).point_set()._naive_enumerator()
            sage: ni.rescalings()
            ((1, 1, 1, 1), (1, 1, 2, 2), (2, 2, 1, 1), (2, 2, 2, 2))
        """
        free = self._Chow_group_free()
        tors = self._Chow_group_torsion()
        if len(tors) == 0:  # optimization for smooth fans
            return free
        result = set(free)
        for f, t in CartesianProduct(free, tors):
            phases = tuple(x*y for x, y in zip(f, t))
            result.add(phases)
        return tuple(sorted(result))

    def orbit(self, point):
        """
        Return the orbit of homogeneous coordinates under rescalings.

        OUTPUT:

        The set of all homogeneous coordinates that are equivalent to ``point``.

        EXAMPLES::

            sage: ne = toric_varieties.P2_123(base_ring=GF(7)).point_set()._naive_enumerator()
            sage: sorted(ne.orbit([1, 0, 0]))
            [(1, 0, 0), (2, 0, 0), (4, 0, 0)]
            sage: sorted(ne.orbit([0, 1, 0]))
            [(0, 1, 0), (0, 6, 0)]
            sage: sorted(ne.orbit([0, 0, 1]))
            [(0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 0, 4), (0, 0, 5), (0, 0, 6)]
            sage: sorted(ne.orbit([1, 1, 0]))
            [(1, 1, 0), (1, 6, 0), (2, 1, 0), (2, 6, 0), (4, 1, 0), (4, 6, 0)]
        """
        result = set()
        for phases in self.rescalings():
            p = tuple(mu*z for mu, z in zip(point, phases))
            result.add(p)
        return frozenset(result)

    def cone_iter(self):
        """
        Iterate over all cones of the fan

        OUTPUT:

        Iterator over the cones, starting with the high-dimensional
        ones.

        EXAMPLES::

            sage: ne = toric_varieties.dP6(base_ring=GF(11)).point_set()._naive_enumerator()
            sage: for cone in ne.cone_iter(): 
            ....:     print cone.ambient_ray_indices()
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
        Iterate over all distinct homogeneous coordinates.

        This method does NOT identify homogeneous coordinates that are
        equivalent by a homogeneous rescaling.
        
        OUTPUT:

        An iterator over the points.

        EXAMPLES::

            sage: F2 = GF(2)
            sage: ni = toric_varieties.P2(base_ring=F2).point_set()._naive_enumerator()
            sage: list(ni.coordinate_iter())
            [(0, 0, 1), (1, 0, 0), (0, 1, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, 1)]

            sage: ni = toric_varieties.P1xP1(base_ring=F2).point_set()._naive_enumerator()
            sage: list(ni.coordinate_iter())
            [(0, 1, 0, 1), (1, 0, 0, 1), (1, 0, 1, 0),
             (0, 1, 1, 0), (0, 1, 1, 1), (1, 0, 1, 1),
             (1, 1, 0, 1), (1, 1, 1, 0), (1, 1, 1, 1)]

        TESTS::

            sage: V = ToricVariety(Fan([Cone([(1,1)])]), base_ring=GF(3))
            sage: ni = V.point_set()._naive_enumerator()
            sage: list(ni.coordinate_iter())
            [(0, 1), (0, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
        """
        units = [x for x in self.ring if x != 0]
        zero = self.ring.zero()
        big_torus = [units] * len(self.rays())
        for cone in self.cone_iter():
            patch = copy(big_torus)
            for i in cone.ambient_ray_indices():
                patch[i] = [zero]
            for p in CartesianProduct(*patch):
                yield tuple(p)

    def __iter__(self):
        """
        Iterate over the distinct points of the toric variety.

        This function does identify orbits under the homogeneous
        rescalings, and returns precisely one representative per
        orbit.

        OUTPUT:

        Iterator over points.

        EXAMPLES:

            sage: ni = toric_varieties.P2(base_ring=GF(2)).point_set()._naive_enumerator()
            sage: list(ni)
            [(0, 0, 1), (1, 0, 0), (0, 1, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, 1)]

            sage: ni = toric_varieties.P1xP1(base_ring=GF(3)).point_set()._naive_enumerator()
            sage: list(ni) 
            [(0, 1, 0, 1), (1, 0, 0, 1), (1, 0, 1, 0), (0, 1, 1, 0), 
             (0, 1, 1, 1), (0, 1, 1, 2), (1, 0, 1, 1), (1, 0, 1, 2), 
             (1, 1, 0, 1), (1, 2, 0, 1), (1, 1, 1, 0), (1, 2, 1, 0), 
             (1, 1, 1, 1), (1, 1, 1, 2), (1, 2, 1, 1), (1, 2, 1, 2)]
        """
        seen = set()
        for p in self.coordinate_iter():
            if p in seen:
                continue
            seen.update(self.orbit(p))
            yield p


class FiniteFieldPointEnumerator(NaiveFinitePointEnumerator):

    @cached_method
    def multiplicative_generator(self):
        """
        Return the multiplicative generator of the finite field.

        OUTPUT:

        A finite field element.

        EXAMPLES:

            sage: point_set = toric_varieties.P2(base_ring=GF(5^2, 'a')).point_set()
            sage: ffe = point_set._finite_field_enumerator()
            sage: ffe.multiplicative_generator()
            a
        """
        return self.ring.multiplicative_generator()

    @cached_method
    def multiplicative_group_order(self):
        return self.ring.multiplicative_generator().multiplicative_order()
        
    @cached_method
    def root_generator(self, n):
        """
        Return a generator for :meth:`roots`.

        INPUT:

        - ``n`` integer.

        OUTPUT:

        A multiplicative generator for :meth:`roots`.

        EXAMPLES::

            sage: point_set = toric_varieties.P2(base_ring=GF(5)).point_set()
            sage: ffe = point_set._finite_field_enumerator()
            sage: ffe.root_generator(2)
            4
            sage: ffe.root_generator(3)
            1
            sage: ffe.root_generator(4)
            2

        TESTS::

            sage: for p in primes(10):
            ....:     for k in range(1,5):
            ....:         F = GF(p^k, 'a')
            ....:         N = F.cardinality() - 1
            ....:         ffe = point_set._finite_field_enumerator(F)
            ....:         assert N == ffe.multiplicative_group_order()
            ....:         for n in N.divisors():
            ....:             x = ffe.root_generator(n)
            ....:             assert set(x**i for i in range(N)) == set(ffe.roots(n))
        """
        N = self.multiplicative_group_order()
        k = N // gcd(n, N)
        return self.multiplicative_generator() ** k

    def _Chow_group_free_generators(self):
        r"""
        Return generators for :meth:`_Chow_group_free_generators`

        OUTPUT:

        A tuple containing generators for $Hom(A_{d-1,\text{free}},
        F^\times)$.

        EXAMPLES::

            sage: fan = NormalFan(ReflexivePolytope(2, 0))
            sage: X = ToricVariety(fan, base_field=GF(7))
            sage: X.Chow_group().degree(1)
            C3 x Z
            sage: enum = X.point_set()._finite_field_enumerator()
            sage: enum._Chow_group_free()
            ((1, 1, 1), (2, 2, 2), (3, 3, 3), (4, 4, 4), (5, 5, 5), (6, 6, 6))
            sage: enum._Chow_group_free_generators()
            ((3, 3, 3),)
        """
        result = []
        null_space = self.rays().matrix().integer_kernel()
        for ker in null_space.basis():
            phases = tuple(self.multiplicative_generator()**exponent
                           for exponent in ker)
            result.append(phases)
        return tuple(sorted(result))

    def _Chow_group_torsion_generators(self):
        r"""
        Return generators for :meth:`Chow_group_torsion`

        OUTPUT:

        A tuple containing generators for
        $Hom(A_{d-1,\text{tors}}, F^\times)$.

        EXAMPLES::

            sage: fan = NormalFan(ReflexivePolytope(2, 0))
            sage: X = ToricVariety(fan, base_field=GF(7))
            sage: X.Chow_group().degree(1)
            C3 x Z
            sage: enum = X.point_set()._finite_field_enumerator()
            sage: enum._Chow_group_torsion()
            ((1, 2, 4), (1, 4, 2))
            sage: enum._Chow_group_torsion_generators()
            ((1, 2, 4),)
        """
        if self.fan.is_smooth():
            return tuple()
        image = self.rays().column_matrix().image()
        torsion = image.saturation().quotient(image)
        result = set()
        for t in torsion.gens():
            t_lift = t.lift()
            root = self.root_generator(t.order())
            if root == 1:
                continue
            phases = tuple(root**exponent for exponent in t_lift)
            result.add(phases)
        assert tuple(self.ring.one() for r in self.rays()) not in result  # because we excluded 1 as root
        return tuple(sorted(result))
    
    def log(self, z):
        """
        Return the component-wise log of ``z``

        INPUT:

        - ``z`` -- a list/tuple/iterable of non-zero finite field
          elements.

        OUTPUT:

        Tuple of integers. The logarithm with base the
        :meth:`multiplicative_generator`.

        EXAMPLES::

            sage: F.<a> = GF(5^2)
            sage: point_set = toric_varieties.P2_123(base_ring=F).point_set()
            sage: ffe = point_set._finite_field_enumerator()
            sage: z = tuple(a^i for i in range(25));  z
            (1, a, a + 3, 4*a + 3, 2*a + 2, 4*a + 1, 2, 2*a, 2*a + 1, 3*a + 1, 
             4*a + 4, 3*a + 2, 4, 4*a, 4*a + 2, a + 2, 3*a + 3, a + 4, 3, 3*a, 
             3*a + 4, 2*a + 4, a + 1, 2*a + 3, 1)
            sage: ffe.log(z)
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 
             17, 18, 19, 20, 21, 22, 23, 0)
            sage: ffe.exp(ffe.log(z)) == z
            True
            sage: ffe.log(ffe.exp(range(24))) == tuple(range(24))
            True
        """
        base = self.multiplicative_generator()
        return tuple(zi.log(base) for zi in z)

    def exp(self, powers):
        """
        Return the component-wise exp of ``z``

        INPUT:

        - ``powers`` -- a list/tuple/iterable of integers.

        OUTPUT:

        Tuple of finite field elements. The powers of the
        :meth:`multiplicative_generator`.

        EXAMPLES::

            sage: F.<a> = GF(5^2)
            sage: point_set = toric_varieties.P2_123(base_ring=F).point_set()
            sage: ffe = point_set._finite_field_enumerator()
            sage: powers = range(24)
            sage: ffe.exp(powers)
            (1, a, a + 3, 4*a + 3, 2*a + 2, 4*a + 1, 2, 2*a, 2*a + 1, 3*a + 1, 
             4*a + 4, 3*a + 2, 4, 4*a, 4*a + 2, a + 2, 3*a + 3, a + 4, 3, 3*a, 
             3*a + 4, 2*a + 4, a + 1, 2*a + 3)
            sage: ffe.log(ffe.exp(powers)) == tuple(powers)
            True
        """
        base = self.multiplicative_generator()
        return tuple(base ** i for i in powers)

    @cached_method
    def rescaling_log_generators(self):
        """
        Return the log generators of :meth:`rescalings`.

        OUTPUT:

        A tuple containing the logarithms (see :meth:`log`) of the
        generators of the multiplicative group of :meth:`rescalings`.
        
        EXAMPLES::

            sage: point_set = toric_varieties.P2_123(base_ring=GF(5)).point_set()
            sage: ffe = point_set._finite_field_enumerator()
            sage: ffe.rescalings()
            ((1, 1, 1), (1, 4, 4), (4, 2, 3), (4, 3, 2))
            sage: map(ffe.log, ffe.rescalings())
            [(0, 0, 0), (0, 2, 2), (2, 1, 3), (2, 3, 1)]
            sage: ffe.rescaling_log_generators()
            ((2, 3, 1),)
        """
        free = self._Chow_group_free_generators()
        tors = self._Chow_group_torsion_generators()
        result = map(self.log, free + tors)
        return tuple(sorted(result))

    def cone_points_iter(self):
        """
        Iterate over the open torus orbits and yield distinct points.
        
        OUTPUT:

        For each open torus orbit (cone): A triple consisting of the
        cone, the nonzero homogeneous coordinates in that orbit (list
        of integers), and the nonzero log coordinates of distinct
        points as a cokernel.

        EXAMPLES::

            sage: fan = NormalFan(ReflexivePolytope(2, 0))
            sage: X = ToricVariety(fan, base_field=GF(7))
            sage: point_set = X.point_set()
            sage: ffe = point_set._finite_field_enumerator()
            sage: cpi = ffe.cone_points_iter()
            sage: cone, nonzero_points, cokernel = list(cpi)[5]
            sage: cone
            1-d cone of Rational polyhedral fan in 2-d lattice N
            sage: cone.ambient_ray_indices()
            (2,)
            sage: nonzero_points
            [0, 1]
            sage: cokernel
            Finitely generated module V/W over Integer Ring with invariants (2)
            sage: list(cokernel)
            [(0), (1)]
            sage: [p.lift() for p in cokernel]
            [(0, 0), (0, 1)]
        """
        from sage.matrix.constructor import matrix, block_matrix, identity_matrix
        from sage.rings.all import ZZ
        nrays = len(self.rays())
        N = self.multiplicative_group_order()
        # Want cokernel of the log rescalings in (ZZ/N)^(#rays). But
        # ZZ/N is not a integral domain. Instead: work over ZZ
        log_generators = self.rescaling_log_generators()
        log_relations = block_matrix(2, 1, [
            matrix(ZZ, len(log_generators), nrays, log_generators), 
            N * identity_matrix(ZZ, nrays)])
        for cone in self.cone_iter():
            nrays = self.fan().nrays() + len(self.fan().virtual_rays())
            nonzero_coordinates = [i for i in range(nrays)
                                   if i not in cone.ambient_ray_indices()]
            log_relations_nonzero = log_relations.matrix_from_columns(nonzero_coordinates)
            image = log_relations_nonzero.image()
            cokernel = image.ambient_module().quotient(image)
            yield cone, nonzero_coordinates, cokernel

    def __iter__(self):
        """
        Iterate over the distinct points of the toric variety.

        This function does identify orbits under the homogeneous
        rescalings, and returns precisely one representative per
        orbit.

        OUTPUT:

        Iterator over points.

        EXAMPLES:

            sage: point_set = toric_varieties.P2(base_ring=GF(2)).point_set()
            sage: ffe = point_set._finite_field_enumerator()
            sage: list(ffe)
            [(0, 0, 1), (1, 0, 0), (0, 1, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0), (1, 1, 1)]

            sage: fan = NormalFan(ReflexivePolytope(2, 0))
            sage: X = ToricVariety(fan, base_field=GF(7))
            sage: point_set = X.point_set()
            sage: ffe = point_set._finite_field_enumerator()
            sage: list(ffe)
            [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1), (0, 1, 3), (1, 0, 1), 
             (1, 0, 3), (1, 1, 0), (1, 3, 0), (1, 1, 1), (1, 1, 3), (1, 1, 2), 
             (1, 1, 6), (1, 1, 4), (1, 1, 5), (1, 3, 2), (1, 3, 6), (1, 3, 4), 
             (1, 3, 5), (1, 3, 1), (1, 3, 3)]
            sage: set(point_set._naive_enumerator()) == set(ffe)
            True
        """
        nrays = len(self.rays())
        for cone, nonzero_coordinates, cokernel in self.cone_points_iter():
            zero = [self.ring.zero()] * nrays
            for v in cokernel:
                z_nonzero = self.exp(v.lift())
                z = copy(zero)
                for i, value in zip(nonzero_coordinates, z_nonzero):
                    z[i] = value
                yield tuple(z)

    def cardinality(self):
        """
        Return the cardinality of the point set.

        OUTPUT:

        Integer. The number of points.

        EXAMPLES::

            sage: fan = NormalFan(ReflexivePolytope(2, 0))
            sage: X = ToricVariety(fan, base_field=GF(7))
            sage: point_set = X.point_set()
            sage: ffe = point_set._finite_field_enumerator()
            sage: ffe.cardinality()
            21
        """
        n = 0
        for cone, nonzero_coordinates, cokernel in self.cone_points_iter():
            n += cokernel.cardinality()
        return n
