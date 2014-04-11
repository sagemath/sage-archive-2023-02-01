"""
Python Backend for Polyhedra

This module implements the python backend for polyhedra. While slower
than specialized C/C++ implementations, it is general enough to work
with any field in Sage that allows you to define polyhedra.
"""

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.matrix.constructor import matrix
from sage.modules.all import vector
from sage.geometry.polyhedron.double_description import StandardAlgorithm as Algorithm

DEBUG = True


def homogenize_ine(A, b):
    """
    Return homogeneous inequalities in matrix form

    EXAMPLES::

        sage: A = matrix([[1,2,3], [4,5,6]])
        sage: b = vector([8,9])
        sage: from sage.geometry.polyhedron.double_description_inhomogeneus import homogenize_ine
        sage: homogenize_ine(A, b)
        [-8| 1  2  3]
        [-9| 4  5  6]
    """    
    from sage.matrix.constructor import block_matrix
    return block_matrix(1, 2, [-b.column(), A])


def dehomogenize_ine(A):
    """
    Return inhomogeneous inequalities

    INPUT:

    - ``A`` -- matrix. Homogeneous inequalities `A x\geq 0`.

    OUTPUT:

    A pair consisting of a matrix and a a vector.

    EXAMPLES::

        sage: A = matrix(ZZ, [[-8, 1, 2, 3], [-9, 4, 5, 6]])
        sage: from sage.geometry.polyhedron.double_description_inhomogeneous import dehomogenize_ine
        sage: dehomogenize_ine(A)
        (
        [1 2 3]        
        [4 5 6], (8, 9)
        )
    """
    c = A.ncols()
    r = A.nrows()
    return A.submatrix(0, 1, r, c-1), -A.column(0)


def homogenize_ext(vertices, rays):
    ring = vertices.base_ring()


def dehomogenize_ext(R):
    vertices = []
    rays = []
    for r in R:
        if r[0] == 0:
            rays.append(r[1:])
        else:
            r = r / r[0]
            vertices.append(r[1:])
    return vertices, rays



def H_to_V_representation(base_ring, dim, inequalities):
    """
    Convert hyperplane to vertex representation

    EXAMPLES::
    
        sage: from sage.geometry.polyhedron.double_description_inhomogeneous import H_to_V_representation
        sage: ieqs = [(1, 1, -1, -1), (1, 1, 1, -1),   (1, 1, 1, 1),   (1, 1, -1, 1),
        ....:         (1, -1, -1, 1), (1, -1, -1, -1), (1, -1, 1, -1), (1, -1, 1, 1)]
        sage: H_to_V_representation(QQ, 3, ieqs)
        ([(-1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, -1), (0, -1, 0), (1, 0, 0)], [])
    """
    A = matrix(base_ring, inequalities)
    DD = Algorithm(A).run()
    R = DD.R
    return dehomogenize_ext(R)





class PivotedInequalities(SageObject):

    def _pivot_inequalities(self, A):
        self._linear_subspace = A.right_kernel()
        self._pivots = A.pivots()
        self._nonpivots = A.nonpivots()
        return A.matrix_from_columns(self._pivots)

    def _unpivot_ray(self, ray):
        """
        Undo the pivoting to go back to the original inequalities containing a linear subspace
        """
        result = [self.base_ring.zero()] * (self.dim + 1)
        for r, i in zip(ray, self._pivots):
            result[i] = r
        return vector(self.base_ring, result)



class Hrep2Vrep(PivotedInequalities):
    """

    EXAMPLES::
    
        sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep
        sage: Hrep2Vrep(QQ, 2, [(1,2,3), (2,4,3)])
        [-1/2|-1/2  1/2|]
        [   0| 2/3 -1/3|]
        sage: Hrep2Vrep(QQ, 2, [(1,2,3), (2,-2,-3)])
        [   1 -1/2||   1]
        [   0    0||-2/3]
        sage: Hrep2Vrep(QQ, 2, [(1,2,3), (2,2,3)])
        [-1/2| 1/2|   1]
        [   0|   0|-2/3]
        sage: Hrep2Vrep(QQ, 2, [(8,7,-2), (1,-4,3), (4,-3,-1)])
        [ 1  0 -2||]
        [ 1  4 -3||]
        sage: Hrep2Vrep(QQ, 2, [(1,2,3), (2,4,3), (5,-1,-2)])
        [-19/5  -1/2| 2/33  1/11|]
        [ 22/5     0|-1/33 -2/33|]
        sage: Hrep2Vrep(QQ, 2, [(0,2,3), (0,4,3), (0,-1,-2)])
        [   0| 1/2  1/3|]
        [   0|-1/3 -1/6|]
    """

    def __init__(self, base_ring, dim, inequalities, equations=[]):
        assert equations == []
        self.base_ring = base_ring
        self.dim = dim
        A = self._init_Vrep(inequalities)
        DD = Algorithm(A).run()
        self._extract_Vrep(DD)
        if DEBUG: 
            self.verify(inequalities)

    def _init_Vrep(self, inequalities):
        """
        Split off the linear subspace from the inequalities and select pivots
        """
        A = matrix(self.base_ring, inequalities)
        return self._pivot_inequalities(A)

    def _normalize_rays(self, rays, lines):
        """
        Shift the rays to be of the form `(0,...)` if possible.

        This is only possible if and only if there is a line
        (i.e. line generator) whose first coordinate is non-zero. It
        is a consequence of the rays not being unique if there is a
        linear subspace.

        OUTPUT:

        A quadruple consisting of
        
        * the rays with their first coordinate non-negative, if
          possible by shifting with lines.

        * the rays whose first coordinate is scaled to one
        
        * the lines whose first coordinate is zero
        
        * the lines whose first coordinate equals one

        Rays with first coordinate negative that cannot be shifted by
        lines are discarderd.

        EXAMPLES::
        
            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep
            sage: ieqs = [(1,2,3), (2,4,3)]
            sage: H = Hrep2Vrep(QQ, 2, [(1,2,3)])
            sage: r1 = vector([1,1,0])
            sage: r2 = vector([1,0,1])
            sage: l = vector([2,1,1])
            sage: H._normalize_rays([r1, r2], [l])
            ([(0, 1/2, -1/2), (0, -1/2, 1/2)], [], [], [(1, 1/2, 1/2)])

        If the line has first coordinate equal zero then it is not possible::

            sage: l = vector([0,1,1])
            sage: H._normalize_rays([r1, r2], [l])
            ([], [(1, 1, 0), (1, 0, 1)], [(0, 1, 1)], [])
        """
        L0 = []
        L1 = []
        for l in lines:
            if l[0] == 0:
                L0.append(l)
            else:
                l = l / l[0]
                L1.append(l)
        if len(L1) == 0:
            R0 = [r        for r in rays if r[0] <= 0]
            R1 = [r / r[0] for r in rays if r[0] >  0]
            return R0, R1, L0, L1
        else:
            # We have a line with first coordinate 1, use it no bring the rays into the special form
            l = L1[0]
            rays = [r - r[0]*l for r in rays]
            return rays, [], L0, L1


    def _split_linear_subspace(self):
        """
        Split the linear subspace in a generator with `x_0\not=0` and the
        remaining generators with `x_0=0`.

        EXAMPLES::
        
            sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep
            sage: H = Hrep2Vrep(QQ, 2, [(1,2,3)])
            sage: H._split_linear_subspace()
            ((1, 0, -1/3), [(0, 1, -2/3)])
            sage: H = Hrep2Vrep(QQ, 2, [(1,0,0)])
            sage: H._split_linear_subspace()
            (None, [(0, 1, 0), (0, 0, 1)])
        """
        lines = self._linear_subspace.matrix().rows()
        L0 = []
        L1 = []
        for l in lines:
            if l[0] == 0:
                L0.append(l)
            else:
                l = l / l[0]
                L1.append(l)
        if len(L1) == 0:
            return None, L0
        else:
            l1 = L1.pop()
            return l1, L0 + [l-l[0]*l1 for l in L1]

    def _extract_Vrep(self, DD):
        """
        Extract the V-representation from the extremal rays of the homogeneous cone

        The V-representation is the intersection of the cone generated
        by the rays `R` and ``self._linear_subspace`` with the
        hyperplane `(1,*,*,...,*)`.
        """
        R = map(self._unpivot_ray, DD.R)

        line1, L0 = self._split_linear_subspace()
        if line1:
            # we can shift all rays to have x_0 = 0, stay extremal
            L1 = [line1]
            R1 = []
            R0 = [r - r[0]*line1 for r in R]
        else:
            # have to really intersect with x_0 = 0
            L1 = []
            R1 = [r / r[0] for r in R if r[0] > 0]
            DD0 = DD.first_coordinate_plane()
            R0 = map(self._unpivot_ray, DD0.R)

        vertices = []
        for v in R1 + L1:
            assert v[0] == 1
            vertices.append(v[1:])
        self.vertices = vertices
        self.rays = [r[1:] for r in R0]
        self.lines = [l[1:] for l in L0]
        
        
    def _repr_(self):
        from sage.matrix.constructor import block_matrix
        def make_matrix(rows):
            return matrix(self.base_ring, len(rows), self.dim, rows).transpose()
        V = make_matrix(self.vertices)
        R = make_matrix(self.rays)
        L = make_matrix(self.lines)
        return str(block_matrix([[V, R, L]]))
    
    def verify(self, inequalities):
        """
        Debug: Compare result to PPL if the base ring is QQ
        """
        from sage.rings.all import QQ
        from sage.geometry.polyhedron.constructor import Polyhedron
        if self.base_ring != QQ:
            return
        P = Polyhedron(vertices=self.vertices, rays=self.rays, lines=self.lines, 
                       base_ring=QQ, ambient_dim=self.dim)
        Q = Polyhedron(ieqs=inequalities,
                       base_ring=QQ, ambient_dim=self.dim)
        if (P != Q) or \
           (len(self.vertices) != P.n_vertices()) or \
           (len(self.rays) != P.n_rays()) or \
           (len(self.lines) != P.n_lines()):
            print 'incorrect!', 
            print Q.Vrepresentation()
            print P.Hrepresentation()
        



class Vrep2Hrep(PivotedInequalities):
    """

    EXAMPLES::
    
        sage: from sage.geometry.polyhedron.double_description_inhomogeneous import Vrep2Hrep
        sage: Vrep2Hrep(QQ, 2, [(-1/2,0)], [(-1/2,2/3), (1/2,-1/3)], [])
        [1 2 3]
        [2 4 3]
        [-----]

        #sage: Vrep2Hrep(QQ, 2, [(1,0), (-1/2,0)], [], [(1,-2/3)])
        #[   0    0 -3/2]
        #[ 1/3  2/3  1/2]
        #[ 2/3 -2/3   -2]
        #[--------------]

        #sage: Vrep2Hrep(QQ, 2, [(-1/2,0)], [(1/2,0)], [(1,-2/3)])
        #[   0    0 -3/2]
        #[   1    0 -3/2]
        #[   1    2  3/2]
        #[--------------]

        sage: Vrep2Hrep(QQ, 2, [(1,1), (0,4), (-2,-3)], [], [])
        [ 8/13  7/13 -2/13]
        [ 1/13 -4/13  3/13]
        [ 4/13 -3/13 -1/13]
        [-----------------]

        sage: Vrep2Hrep(QQ, 2, [(-19/5,22/5), (-1/2,0)], [(2/33,-1/33), (1/11,-2/33)], [])
        [10/11 -2/11 -4/11]
        [ 66/5 132/5  99/5]
        [ 2/11  4/11  6/11]
        [-----------------]

        sage: Vrep2Hrep(QQ, 2, [(0,0)], [(1/2,-1/3), (1/3,-1/6)], [])
        [  0  -6 -12]
        [  0  12  18]
        [-----------]
    """

    def __init__(self, base_ring, dim, vertices, rays, lines):
        assert len(lines) == 0    #todo
        self.base_ring = base_ring
        self.dim = dim
        A = self._init_Vrep(vertices, rays, lines)
        DD = Algorithm(A).run()
        self._extract_Hrep(DD)
        if DEBUG: 
            self.verify(vertices, rays, lines)

    def _init_Vrep(self, vertices, rays, lines):
        """
        Split off the linear subspace from the inequalities and select pivots
        """
        homogeneous = \
            [ [1] + list(v) for v in vertices ] + \
            [ [0] + list(r) for r in rays ] + \
            [ [0] + list(l) for l in lines] + \
            [ [0] + [-x for x in l] for l in lines]
        A = matrix(self.base_ring, homogeneous)
        return self._pivot_inequalities(A)

    def _extract_Hrep(self, DD):
        def is_trivial(ray):
            # trivial Hrep output 1 >= 0
            return ray[0] > 0 and all(r == 0 for r in ray[1:])
        ieqs = map(self._unpivot_ray, DD.R)
        self.inequalities = [r for r in ieqs if not is_trivial(r)]
        self.equations = self._linear_subspace.matrix().rows()

    def _repr_(self):
        from sage.matrix.constructor import block_matrix
        def make_matrix(cols):
            return matrix(self.base_ring, len(cols), self.dim+1, cols)
        I = make_matrix(self.inequalities)
        E = make_matrix(self.equations)
        return str(block_matrix([[I], [E]]))

    def verify(self, vertices, rays, lines):
        """
        Debug: Compare result to PPL if the base ring is QQ
        """
        from sage.rings.all import QQ
        from sage.geometry.polyhedron.constructor import Polyhedron
        if self.base_ring != QQ:
            return
        P = Polyhedron(vertices=vertices, rays=rays, lines=lines, 
                       base_ring=QQ, ambient_dim=self.dim)
        Q = Polyhedron(ieqs=self.inequalities, eqns=self.equations,
                       base_ring=QQ, ambient_dim=self.dim)
        if not P == Q:
            print 'incorrect!', 
            print Q.Vrepresentation()
            print P.Hrepresentation()
