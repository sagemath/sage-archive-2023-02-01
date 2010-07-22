r"""
Library of toric varieties

This module provides a simple way to construct often-used toric
varieties. Please see the help for the individual methods of
``toric_varieties`` for a more detailed description of which varieties
can be constructed.

AUTHORS:

- Volker Braun (2010-07-02): initial version

EXAMPLES::

    sage: toric_varieties.dP6()
    2-d CPR-Fano toric variety covered by 6 affine patches
"""

#*****************************************************************************
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject

from sage.matrix.all import matrix, identity_matrix
from sage.geometry.fan import Fan
from sage.geometry.lattice_polytope import LatticePolytope
from sage.schemes.generic.toric_variety import ToricVariety
from sage.schemes.generic.fano_toric_variety import CPRFanoToricVariety




# The combinatorial data of the toric varieties is stored separately here
# since we might want to use it later on to do the reverse lookup.
toric_varieties_rays_cones = {
    'dP6':[
        [(0, 1), (-1, 0), (-1, -1), (0, -1), (1, 0), (1, 1)],
        [[0,1],[1,2],[2,3],[3,4],[4,5],[5,0]] ],
    'dP7':[
        [(0, 1), (-1, 0), (-1, -1), (0, -1), (1, 0)],
        [[0,1],[1,2],[2,3],[3,4],[4,0]] ],
    'dP8':[
        [(1,1), (0, 1), (-1, -1), (1, 0)],
        [[0,1],[1,2],[2,3],[3,0]]
        ],
    'P1xP1':[
        [(1, 0), (-1, 0), (0, 1), (0, -1)],
        [[0,2],[2,1],[1,3],[3,0]] ],
    'P1':[
        [(1,), (-1,)],
        [[0],[1]] ],
    'P2':[
        [(1,0), (0, 1), (-1, -1)],
        [[0,1],[1,2],[2,0]] ],
    'A1':[
        [(1,)],
        [[0]] ],
    'A2':[
        [(1, 0), (0, 1)],
        [[0,1]] ],
    'A2_Z2':[
        [(1, 0), (1, 2)],
        [[0,1]] ],
    'P1xA1':[
        [(1, 0), (-1, 0), (0, 1)],
        [[0,2],[2,1]] ],
    'Conifold':[
        [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)],
        [[0,1,2,3]] ],
    'dP6xdP6':[
        [(0, 1, 0, 0), (-1, 0, 0, 0), (-1, -1, 0, 0),
         (0, -1, 0, 0), (1, 0, 0, 0), (1, 1, 0, 0),
         (0, 0, 0, 1), (0, 0, -1, 0), (0, 0, -1, -1),
         (0, 0, 0, -1), (0, 0, 1, 0), (0, 0, 1, 1)],
        [[0, 1, 6, 7], [0, 1, 7, 8], [0, 1, 8, 9], [0, 1, 9, 10],
         [0, 1, 10, 11], [0, 1, 6, 11], [1, 2, 6, 7], [1, 2, 7, 8],
         [1, 2, 8, 9], [1, 2, 9, 10], [1, 2, 10, 11], [1, 2, 6, 11],
         [2, 3, 6, 7], [2, 3, 7, 8], [2, 3, 8, 9], [2, 3, 9, 10],
         [2, 3, 10, 11], [2, 3, 6, 11], [3, 4, 6, 7], [3, 4, 7, 8],
         [3, 4, 8, 9], [3, 4, 9, 10], [3, 4, 10, 11], [3, 4, 6, 11],
         [4, 5, 6, 7], [4, 5, 7, 8], [4, 5, 8, 9], [4, 5, 9, 10],
         [4, 5, 10, 11], [4, 5, 6, 11], [0, 5, 6, 7], [0, 5, 7, 8],
         [0, 5, 8, 9], [0, 5, 9, 10], [0, 5, 10, 11], [0, 5, 6, 11]] ],
    'Cube_face_fan':[
        [(1, 1, 1), (1, -1, 1), (-1, 1, 1), (-1, -1, 1),
         (-1, -1, -1), (-1, 1, -1), (1, -1, -1), (1, 1, -1)],
        [[0,1,2,3], [4,5,6,7], [0,1,7,6], [4,5,3,2], [0,2,5,7], [4,6,1,3]] ],
    'Cube_sublattice':[
        [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 1, 1),
         (-1, 0, 0), (0, -1, 0), (0, 0, -1), (1, -1, -1)],
        [[0,1,2,3],[4,5,6,7],[0,1,7,6],[4,5,3,2],[0,2,5,7],[4,6,1,3]] ],
    'Cube_nonpolyhedral':[
        [(1, 2, 3), (1, -1, 1), (-1, 1, 1), (-1, -1, 1),
         (-1, -1, -1), (-1, 1, -1), (1, -1, -1), (1, 1, -1)],
        [[0,1,2,3],[4,5,6,7],[0,1,7,6],[4,5,3,2],[0,2,5,7],[4,6,1,3]] ],
    'BCdlOG':[
        [(-1, 0, 0, 2, 3),  #  0
         ( 0,-1, 0, 2, 3),  #  1
         ( 0, 0,-1, 2, 3),  #  2
         ( 0, 0,-1, 1, 2),  #  3
         ( 0, 0, 0,-1, 0),  #  4
         ( 0, 0, 0, 0,-1),  #  5
         ( 0, 0, 0, 2, 3),  #  6
         ( 0, 0, 1, 2, 3),  #  7
         ( 0, 0, 2, 2, 3),  #  8
         ( 0, 0, 1, 1, 1),  #  9
         ( 0, 1, 2, 2, 3),  # 10
         ( 0, 1, 3, 2, 3),  # 11
         ( 1, 0, 4, 2, 3)], # 12
        [ [0,6,7,1,4],   [0,6,10,2,4],  [0,6,1,2,4],   [0,9,7,1,5],  [0,6,7,1,5],
          [0,6,10,2,5],  [0,6,1,2,5],   [0,9,1,4,5],   [0,6,10,4,11],[0,6,7,4,11],
          [0,6,10,5,11], [0,9,7,5,11],  [0,6,7,5,11],  [0,9,4,5,11], [0,10,4,5,11],
          [0,9,7,1,8],   [0,9,1,4,8],   [0,7,1,4,8],   [0,9,7,11,8], [0,9,4,11,8],
          [0,7,4,11,8],  [0,10,2,4,3],  [0,1,2,4,3],   [0,10,2,5,3], [0,1,2,5,3],
          [0,10,4,5,3],  [0,1,4,5,3],   [12,6,7,1,4],  [12,6,10,2,4],[12,6,1,2,4],
          [12,9,7,1,5],  [12,6,7,1,5],  [12,6,10,2,5], [12,6,1,2,5], [12,9,1,4,5],
          [12,6,10,4,11],[12,6,7,4,11], [12,6,10,5,11],[12,9,7,5,11],[12,6,7,5,11],
          [12,9,4,5,11], [12,10,4,5,11],[12,9,7,1,8],  [12,9,1,4,8], [12,7,1,4,8],
          [12,9,7,11,8], [12,9,4,11,8], [12,7,4,11,8], [12,10,2,4,3],[12,1,2,4,3],
          [12,10,2,5,3], [12,1,2,5,3],  [12,10,4,5,3], [12,1,4,5,3] ]  ],
    'BCdlOG_base':[
        [(-1, 0, 0),
         ( 0,-1, 0),
         ( 0, 0,-1),
         ( 0, 0, 1),
         ( 0, 1, 2),
         ( 0, 1, 3),
         ( 1, 0, 4)],
        [[0,4,2],[0,4,5],[0,5,3],[0,1,3],[0,1,2],
         [6,4,2],[6,4,5],[6,5,3],[6,1,3],[6,1,2]] ],
    'P4_11169':[
        [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (-9, -6, -1, -1)],
        [[0,1,2,3],[0,1,2,4],[0,1,3,4],[0,2,3,4],[1,2,3,4]] ],
    'P4_11169_resolved':[
        [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (-9, -6, -1, -1), (-3, -2, 0, 0)],
        [[0, 1, 2, 3], [0, 1, 3, 4], [0, 1, 2, 4], [1, 3, 4, 5], [0, 3, 4, 5],
         [1, 2, 4, 5], [0, 2, 4, 5], [1, 2, 3, 5], [0, 2, 3, 5]] ]
}



class ToricVarietyFactory(SageObject):
    r"""
    The methods of this class construct toric varieties.

    .. WARNING::

        You need not create instances of this class. Use the
        already-provided object ``toric_varieties`` instead.
    """

    _check = True

    def _make_ToricVariety(self, name, coordinate_names):
        """
        Construct a toric variety and cache the result.

        INPUT:

        - ``name`` -- string. One of the pre-defined names in the
          ``toric_varieties_rays_cones`` data structure.

        - ``coordinate_names`` -- A string describing the names of the
          homogeneous coordinates of the toric variety.

        OUTPUT:

        A :class:`toric variety
        <sage.schemes.generic.toric_variety.ToricVariety_field>`.

        EXAMPLES::

            sage: toric_varieties.A1()           # indirect doctest
            1-d affine toric variety
        """
        dict_key = '_'+name
        if dict_key not in self.__dict__:
            rays, cones = toric_varieties_rays_cones[name]
            fan = Fan(cones, rays, check=self._check)
            self.__dict__[dict_key] = \
                ToricVariety(fan,
                             coordinate_names=coordinate_names)
        return self.__dict__[dict_key]

    def _make_CPRFanoToricVariety(self, name, coordinate_names):
        """
        Construct a (crepant partially resolved) Fano toric variety
        and cache the result.

        INPUT:

        - ``name`` -- string. One of the pre-defined names in the
          ``toric_varieties_rays_cones`` data structure.

        - ``coordinate_names`` -- A string describing the names of the
          homogeneous coordinates of the toric variety.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: toric_varieties.P2()           # indirect doctest
            2-d CPR-Fano toric variety covered by 3 affine patches
        """
        dict_key = '_'+name
        if dict_key not in self.__dict__:
            rays, cones = toric_varieties_rays_cones[name]
            polytope = LatticePolytope( matrix(rays).transpose() )
            points = map(tuple, polytope.points().columns())
            ray2point = [points.index(r) for r in rays]
            charts = [ [ray2point[i] for i in c] for c in cones ]
            self.__dict__[dict_key] = \
                CPRFanoToricVariety(Delta_polar=polytope,
                                    coordinate_points=ray2point,
                                    charts=charts,
                                    coordinate_names=coordinate_names,
                                    check=self._check)
        return self.__dict__[dict_key]

    def dP6(self):
        r"""
        Construct the del Pezzo surface of degree 6 (`\mathbb{P}^2`
        blown up at 3 points) as a toric variety.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: dP6
            2-d CPR-Fano toric variety covered by 6 affine patches
            sage: dP6.fan().ray_matrix()
            [ 0 -1 -1  0  1  1]
            [ 1  0 -1 -1  0  1]
            sage: dP6.gens()
            (x, u, y, v, z, w)
        """
        return self._make_CPRFanoToricVariety('dP6', 'x u y v z w')

    def dP7(self):
        r"""
        Construct the del Pezzo surface of degree 7 (`\mathbb{P}^2`
        blown up at 2 points) as a toric variety.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: dP7 = toric_varieties.dP7()
            sage: dP7
            2-d CPR-Fano toric variety covered by 5 affine patches
            sage: dP7.fan().ray_matrix()
            [ 0 -1 -1  0  1]
            [ 1  0 -1 -1  0]
            sage: dP7.gens()
            (x, u, y, v, z)
        """
        return self._make_CPRFanoToricVariety('dP7', 'x u y v z')

    def dP8(self):
        r"""
        Construct the del Pezzo surface of degree 8 (`\mathbb{P}^2`
        blown up at 1 point) as a toric variety.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: dP8 = toric_varieties.dP8()
            sage: dP8
            2-d CPR-Fano toric variety covered by 4 affine patches
            sage: dP8.fan().ray_matrix()
            [ 1  0 -1  1]
            [ 1  1 -1  0]
            sage: dP8.gens()
            (t, x, y, z)
        """
        return self._make_CPRFanoToricVariety('dP8', 't x y z')

    def P1xP1(self):
        r"""
        Construct the del Pezzo surface `\mathbb{P}^1 \times
        \mathbb{P}^1` as a toric variety.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1
            2-d CPR-Fano toric variety covered by 4 affine patches
            sage: P1xP1.fan().ray_matrix()
            [ 1 -1  0  0]
            [ 0  0  1 -1]
            sage: P1xP1.gens()
            (s, t, x, y)
        """
        return self._make_CPRFanoToricVariety('P1xP1', 's t x y')

    def P1(self):
        r"""
        Construct the projective line `\mathbb{P}^1` as a toric
        variety.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: P1 = toric_varieties.P1()
            sage: P1
            1-d CPR-Fano toric variety covered by 2 affine patches
            sage: P1.fan().ray_matrix()
            [ 1 -1]
            sage: P1.gens()
            (s, t)
        """
        return self._make_CPRFanoToricVariety('P1', 's t')

    def P2(self):
        r"""
        Construct the projective plane `\mathbb{P}^2` as a toric
        variety.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: P2
            2-d CPR-Fano toric variety covered by 3 affine patches
            sage: P2.fan().ray_matrix()
            [ 1  0 -1]
            [ 0  1 -1]
            sage: P2.gens()
            (x, y, z)
        """
        return self._make_CPRFanoToricVariety('P2', 'x y z')

    def P(self, n):
        r"""
        Construct the ``n``-dimensional projective space `\mathbb{P}^n`.

        INPUT:

        - ``n`` -- integer. The dimension of the projective space.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: P3 = toric_varieties.P(3)
            sage: P3
            3-d toric variety covered by 4 affine patches
            sage: P3.fan().ray_matrix()
            [-1  1  0  0]
            [-1  0  1  0]
            [-1  0  0  1]
            sage: P3.gens()
            (z0, z1, z2, z3)
        """
        rays = [ tuple([-1]*n) ] + identity_matrix(n).columns()
        cones = [ range(0,i)+range(i+1,n+1) for i in range(0,n+1) ]
        fan = Fan(cones, rays)
        return ToricVariety(fan)

    def A1(self):
        r"""
        Construct the affine line `\mathbb{A}^1` as a toric variety.

        OUTPUT:

        A :class:`toric variety
        <sage.schemes.generic.toric_variety.ToricVariety_field>`.

        EXAMPLES::

            sage: A1 = toric_varieties.A1()
            sage: A1
            1-d affine toric variety
            sage: A1.fan().ray_matrix()
            [1]
            sage: A1.gens()
            (z,)
        """
        return self._make_ToricVariety('A1', 'z')

    def A2(self):
        r"""
        Construct the affine plane `\mathbb{A}^2` as a toric variety.

        OUTPUT:

        A :class:`toric variety
        <sage.schemes.generic.toric_variety.ToricVariety_field>`.

        EXAMPLES::

            sage: A2 = toric_varieties.A2()
            sage: A2
            2-d affine toric variety
            sage: A2.fan().ray_matrix()
            [1 0]
            [0 1]
            sage: A2.gens()
            (x, y)
        """
        return self._make_ToricVariety('A2', 'x y')

    def A(self, n):
        r"""
        Construct the ``n``-dimensional affine space.

        INPUT:

        - ``n`` -- integer. The dimension of the affine space.

        OUTPUT:

        A :class:`toric variety
        <sage.schemes.generic.toric_variety.ToricVariety_field>`.

        EXAMPLES::

            sage: A3 = toric_varieties.A(3)
            sage: A3
            3-d affine toric variety
            sage: A3.fan().ray_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: A3.gens()
            (z0, z1, z2)
        """
        rays = identity_matrix(n).columns()
        cones = [ range(0,n) ]
        fan = Fan(cones, rays)
        return ToricVariety(fan)

    def A2_Z2(self):
        r"""
        Construct the orbifold `\mathbb{A}^2 / \ZZ_2` as a toric
        variety.

        OUTPUT:

        A :class:`toric variety
        <sage.schemes.generic.toric_variety.ToricVariety_field>`.

        EXAMPLES::

            sage: A2_Z2 = toric_varieties.A2_Z2()
            sage: A2_Z2
            2-d affine toric variety
            sage: A2_Z2.fan().ray_matrix()
            [1 1]
            [0 2]
            sage: A2_Z2.gens()
            (x, y)
        """
        return self._make_ToricVariety('A2_Z2', 'x y')

    def P1xA1(self):
        r"""
        Construct the cartesian product `\mathbb{P}^1 \times \mathbb{A}^1` as
        a toric variety.

        OUTPUT:

        A :class:`toric variety
        <sage.schemes.generic.toric_variety.ToricVariety_field>`.

        EXAMPLES::

            sage: P1xA1 = toric_varieties.P1xA1()
            sage: P1xA1
            2-d toric variety covered by 2 affine patches
            sage: P1xA1.fan().ray_matrix()
            [ 1 -1  0]
            [ 0  0  1]
            sage: P1xA1.gens()
            (s, t, z)
        """
        return self._make_ToricVariety('P1xA1', 's t z')

    def Conifold(self):
        r"""
        Construct the conifold as a toric variety.

        OUTPUT:

        A :class:`toric variety
        <sage.schemes.generic.toric_variety.ToricVariety_field>`.

        EXAMPLES::

            sage: Conifold = toric_varieties.Conifold()
            sage: Conifold
            3-d affine toric variety
            sage: Conifold.fan().ray_matrix()
            [0 0 1 1]
            [0 1 0 1]
            [1 1 1 1]
            sage: Conifold.gens()
            (u, x, y, v)
        """
        return self._make_ToricVariety('Conifold', 'u x y v')

    def dP6xdP6(self):
        r"""
        Construct the product of two del Pezzo surfaces of degree 6
        (`\mathbb{P}^2` blown up at 3 points) as a toric variety.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: dP6xdP6 = toric_varieties.dP6xdP6()   # long time
            sage: dP6xdP6                               # long time
            4-d CPR-Fano toric variety covered by 36 affine patches
            sage: dP6xdP6.fan().ray_matrix()            # long time
            [ 0 -1 -1  0  1  1  0  0  0  0  0  0]
            [ 1  0 -1 -1  0  1  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0 -1 -1  0  1  1]
            [ 0  0  0  0  0  0  1  0 -1 -1  0  1]
            sage: dP6xdP6.gens()                        # long time
            (x0, x1, x2, x3, x4, x5, y0, y1, y2, y3, y4, y5)
        """
        return self._make_CPRFanoToricVariety('dP6xdP6', 'x0 x1 x2 x3 x4 x5 y0 y1 y2 y3 y4 y5')

    def Cube_face_fan(self):
        r"""
        Construct the toric variety given by the face fan of the
        3-dimensional unit lattice cube.

        This variety has 6 conifold singularities but the fan is still
        polyhedral.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: Cube_face_fan = toric_varieties.Cube_face_fan()
            sage: Cube_face_fan
            3-d CPR-Fano toric variety covered by 6 affine patches
            sage: Cube_face_fan.fan().ray_matrix()
            [ 1  1 -1 -1 -1 -1  1  1]
            [ 1 -1  1 -1 -1  1 -1  1]
            [ 1  1  1  1 -1 -1 -1 -1]
            sage: Cube_face_fan.gens()
            (z0, z1, z2, z3, z4, z5, z6, z7)
        """
        return self._make_CPRFanoToricVariety('Cube_face_fan', None)

    def Cube_sublattice(self):
        r"""
        Construct the toric variety defined by a face fan over a
        3-dimensional cube, but not the unit cube in the
        N-lattice. See [FultonP65]_.

        Its Chow group is `A_2(X)=\mathbb{Z}^5`, which distinguishes
        it from the face fan of the unit cube.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: Cube_sublattice = toric_varieties.Cube_sublattice()
            sage: Cube_sublattice
            3-d CPR-Fano toric variety covered by 6 affine patches
            sage: Cube_sublattice.fan().ray_matrix()
            [ 1  0  0 -1 -1  0  0  1]
            [ 0  1  0  1  0 -1  0 -1]
            [ 0  0  1  1  0  0 -1 -1]
            sage: Cube_sublattice.gens()
            (z0, z1, z2, z3, z4, z5, z6, z7)

        REFERENCES:

        ..  [FultonP65]
            Page 65, 3rd exercise (Section 3.4) of Wiliam Fulton,
            "Introduction to Toric Varieties", Princeton University
            Press
        """
        return self._make_CPRFanoToricVariety('Cube_sublattice', None)

    def Cube_nonpolyhedral(self):
        r"""
        Construct the toric variety defined by a fan that is not the
        face fan of a polyhedron.

        This toric variety is defined by a fan that is topologically
        like the face fan of a 3-dimensional cube, but with a
        different N-lattice structure.

        OUTPUT:

        A :class:`toric variety
        <sage.schemes.generic.toric_variety.ToricVariety_field>`.

        NOTES:

        * This is an example of an non-polyhedral fan.

        * Its Chow group has torsion: `A_2(X)=\ZZ^5 \oplus \ZZ_2`

        EXAMPLES::

            sage: Cube_nonpolyhedral = toric_varieties.Cube_nonpolyhedral()
            sage: Cube_nonpolyhedral
            3-d toric variety covered by 6 affine patches
            sage: Cube_nonpolyhedral.fan().ray_matrix()
            [ 1  1 -1 -1 -1 -1  1  1]
            [ 2 -1  1 -1 -1  1 -1  1]
            [ 3  1  1  1 -1 -1 -1 -1]
            sage: Cube_nonpolyhedral.gens()
            (z0, z1, z2, z3, z4, z5, z6, z7)
        """
        return self._make_ToricVariety('Cube_nonpolyhedral', None)

    def Cube_deformation(self,k):
        r"""
        Construct, for each `k\in\ZZ`, a toric variety with
        `\ZZ_k`-torsion in the Chow group.

        The fans of this sequence of toric varieties all equal the
        face fan of a unit cube topologically, but the
        ``(1,1,1)``-vertex is moved to ``(1,1,2k+1)``. This example
        was studied in [FS]_.

        INPUT:

        - ``k`` -- integer. The case ``k=0`` is the same as
          :meth:`Cube_face_fan`.

        OUTPUT:

        A :class:`toric variety
        <sage.schemes.generic.toric_variety.ToricVariety_field>`
        `X_k`. Its Chow group is `A_1(X_k)=\ZZ_k`.

        EXAMPLES::

            sage: X_2 = toric_varieties.Cube_deformation(2)
            sage: X_2
            3-d toric variety covered by 6 affine patches
            sage: X_2.fan().ray_matrix()
            [ 1  1 -1 -1 -1 -1  1  1]
            [ 1 -1  1 -1 -1  1 -1  1]
            [ 5  1  1  1 -1 -1 -1 -1]
            sage: X_2.gens()
            (z0, z1, z2, z3, z4, z5, z6, z7)

        REFERENCES:

        ..  [FS]
            William Fulton, Bernd Sturmfels, "Intersection Theory on
            Toric Varieties", http://arxiv.org/abs/alg-geom/9403002
        """
        rays = lambda kappa: matrix([[ 1, 1, 2*kappa+1],[ 1,-1, 1],[-1, 1, 1],[-1,-1, 1],
                                       [-1,-1,-1],[-1, 1,-1],[ 1,-1,-1],[ 1, 1,-1]])
        cones = [[0,1,2,3],[4,5,6,7],[0,1,7,6],[4,5,3,2],[0,2,5,7],[4,6,1,3]]
        fan = Fan(cones, rays(k))
        return ToricVariety(fan)

    def BCdlOG(self):
        r"""
        Construct the 5-dimensional toric variety studied in
        [BCdlOG]_, [HLY]_

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: X = toric_varieties.BCdlOG()      # long time
            sage: X                                 # long time
            5-d CPR-Fano toric variety covered by 54 affine patches
            sage: X.fan().ray_matrix()              # long time
            [-1  0  0  0  0  0  0  0  0  0  0  0  1]
            [ 0 -1  0  0  0  0  0  0  0  0  1  1  0]
            [ 0  0 -1 -1  0  0  0  1  2  1  2  3  4]
            [ 2  2  2  1 -1  0  2  2  2  1  2  2  2]
            [ 3  3  3  2  0 -1  3  3  3  1  3  3  3]
            sage: X.gens()                          # long time
            (v1, v2, c1, c2, v4, v5, b, e1, e2, e3, f, g, v6)

        REFERENCES:

        ..  [BCdlOG]
            Volker Braun, Philip Candelas, Xendia de la Ossa,
            Antonella Grassi, "Toric Calabi-Yau Fourfolds, Duality
            Between N=1 Theories and Divisors that Contribute to the
            Superpotential", http://arxiv.org/abs/hep-th/0001208

        ..  [HLY]
            Yi Hu, Chien-Hao Liu, Shing-Tung Yau, "Toric morphisms and
            fibrations of toric Calabi-Yau hypersurfaces",
            http://arxiv.org/abs/math/0010082
            """
        return self._make_CPRFanoToricVariety('BCdlOG', 'v1 v2 c1 c2 v4 v5 b e1 e2 e3 f g v6')

    def BCdlOG_base(self):
        r"""
        Construct the base of the `\mathbb{P}^2(1,2,3)` fibration
        :meth:`BCdlOG`.

        OUTPUT:

        A :class:`toric variety
        <sage.schemes.generic.toric_variety.ToricVariety_field>`.

        EXAMPLES::

            sage: base = toric_varieties.BCdlOG_base()
            sage: base
            3-d toric variety covered by 10 affine patches
            sage: base.fan().ray_matrix()
            [-1  0  0  0  0  0  1]
            [ 0 -1  0  0  1  1  0]
            [ 0  0 -1  1  2  3  4]
            sage: base.gens()
            (d4, d3, r2, r1, d2, u, d1)
        """
        return self._make_ToricVariety('BCdlOG_base', 'd4 d3 r2 r1 d2 u d1')

    def P4_11169(self):
        r"""
        Construct the weighted projective space
        `\mathbb{P}^4(1,1,1,6,9)`.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: P4_11169 = toric_varieties.P4_11169()
            sage: P4_11169
            4-d CPR-Fano toric variety covered by 5 affine patches
            sage: P4_11169.fan().ray_matrix()
            [ 1  0  0  0 -9]
            [ 0  1  0  0 -6]
            [ 0  0  1  0 -1]
            [ 0  0  0  1 -1]
            sage: P4_11169.gens()
            (z0, z1, z2, z3, z4)
        """
        return self._make_CPRFanoToricVariety('P4_11169', None)

    def P4_11169_resolved(self):
        r"""
        Construct the blow-up of the weighted projective space
        `\mathbb{P}^4(1,1,1,6,9)` at its curve of `\ZZ_3` orbifold
        fixed points.

        OUTPUT:

        A :class:`CPR-Fano toric variety
        <sage.schemes.generic.fano_toric_variety.CPRFanoToricVariety_field>`.

        EXAMPLES::

            sage: P4_11169_resolved = toric_varieties.P4_11169_resolved()  # long time
            sage: P4_11169_resolved                      # long time
            4-d CPR-Fano toric variety covered by 9 affine patches
            sage: P4_11169_resolved.fan().ray_matrix()   # long time
            [ 1  0  0  0 -9 -3]
            [ 0  1  0  0 -6 -2]
            [ 0  0  1  0 -1  0]
            [ 0  0  0  1 -1  0]
            sage: P4_11169_resolved.gens()               # long time
            (z0, z1, z2, z3, z4, z5)
        """
        return self._make_CPRFanoToricVariety('P4_11169_resolved', None)


toric_varieties = ToricVarietyFactory()
