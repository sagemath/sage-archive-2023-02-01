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



def homogenize_ine(A, b):
    """
    Return homogeneous inequalities in matrix form

    EXAMPLES::

        sage: A = matrix([[1,2,3], [4,5,6]])
        sage: b = vector([8,9])
        sage: from sage.geometry.polyhedron.double_description_backend import homogenize_ine
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
        sage: from sage.geometry.polyhedron.double_description_backend import dehomogenize_ine
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
    
        sage: from sage.geometry.polyhedron.double_description_backend import H_to_V_representation
        sage: ieqs = [(1, 1, -1, -1), (1, 1, 1, -1),   (1, 1, 1, 1),   (1, 1, -1, 1),
        ....:         (1, -1, -1, 1), (1, -1, -1, -1), (1, -1, 1, -1), (1, -1, 1, 1)]
        sage: H_to_V_representation(QQ, 3, ieqs)
        ([(-1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, -1), (0, -1, 0), (1, 0, 0)], [])
    """
    A = matrix(base_ring, inequalities)
    DD = Algorithm(A).run()
    R = DD.R
    return dehomogenize_ext(R)





class Hrep2Vrep(SageObject):
    """

    EXAMPLES::
    
        sage: from sage.geometry.polyhedron.double_description_backend import Hrep2Vrep
        sage: ieqs = [(1,2,3), (2,4,3)]
        sage: H = Hrep2Vrep(QQ, 3, ieqs)
        sage: H._pivots
        sage: H.vertices
        sage: H.lines
        sage: H.rays
        sage: H
    """

    def __init__(self, base_ring, dim, inequalities):
        self.base_ring = base_ring
        self.dim = dim
        A = self._init_inequalities(inequalities)
        DD = Algorithm(A).run()
        self._init_Vrep(DD.R)

    def _init_inequalities(self, inequalities):
        """
        Split off the linear subspace from the inequalities and select pivots
        """
        A = matrix(self.base_ring, inequalities)
        self._linear_subspace = A.right_kernel()
        self._pivots = A.pivots()
        self._nonpivots = A.nonpivots()
        return A.matrix_from_columns(self._pivots)

    def _unpivot_ray(self, ray):
        """
        Undo the pivoting to go back to the original inequalities containing a linear subspace
        """
        result = [self.base_ring.zero()] * self.dim
        for r, i in zip(ray, self._pivots):
            result[i] = r
        return vector(self.base_ring, result)

    def _normalize_rays(self, rays, lines):
        """
        Shift the rays to be of the form `(0,...)` if possible.

        This is only possible if and only if there is a line
        (i.e. line generator) whose first coordinate is non-zero. It
        is a consequence of the rays not being unique if there is a
        linear subspace.

        OUTPUT:

        A quadruple consisting of
        
        * the rays with their first coordinate shifted to zero if possible

        * the rays whose first coordinate is scaled to one
        
        * the lines whose first coordinate is zero
        
        * the lines whose first coordinate equals one

        Rays with first coordinate negative that cannot be shifted by
        lines are discarderd.

        EXAMPLES::
        
            sage: from sage.geometry.polyhedron.double_description_backend import Hrep2Vrep
            sage: ieqs = [(1,2,3), (2,4,3)]
            sage: H = Hrep2Vrep(QQ, 3, [(1,2,3)])
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
            R0 = [r        for r in rays if r[0] == 0]
            R1 = [r / r[0] for r in rays if r[0] >  0]
            return R0, R1, L0, L1
        else:
            # We have a line with first coordinate 1, use it no bring the rays into the special form
            l = L1[0]
            rays = [r - r[0]*l for r in rays]
            return rays, [], L0, L1

    def _init_Vrep(self, R):
        """
        Extract the V-representation from the extremal rays of the homogeneous cone

        The V-representation is the intersection of the cone generated
        by the rays `R` and ``self._linear_subspace`` with the
        hyperplane `(1,*,*,...,*)`.
        """
        R = map(self._unpivot_ray, R)
        L = self._linear_subspace.matrix().rows()
        R0, R1, L0, L1 = self._normalize_rays(R, L)
        vertices = []
        for v in R1 + L1:
            assert v[0] == 1
            vertices.append(v[1:])
        self.vertices = vertices
        self.rays = [r[1:] for r in R0]
        self.lines = [l[1:] for l in L0]
        
        
    def _repr_(self):
        from sage.matrix.constructor import block_matrix
        V = matrix(self.vertices).transpose()
        R = matrix(self.rays).transpose()
        L = matrix(self.lines).transpose()
        return str(block_matrix([[V, R, ]]))
    
        

