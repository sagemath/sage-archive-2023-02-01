"""
Double Description Algorithm

This module implements the double description algorithm for extremal
vertex enumeration in a pointed cone. With a little bit of
preprocessing (see :mod:`double_description_backend`) this defines a
backend for polyhedral computations.

EXAMPLES::

    sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
    sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
    sage: alg = StandardAlgorithm(A);  alg
    (1, 0, 1)
    (0, 1, 1)
    (-1, -1, 1)
    sage: DD, _ = alg.initial_pair();  DD
    Double description pair (A, R) defined by
        [ 1  0  1]        [ 2/3 -1/3 -1/3]
    A = [ 0  1  1],   R = [-1/3  2/3 -1/3]
        [-1 -1  1]        [ 1/3  1/3  1/3]

The implementation works over any field that is embedded in `\RR` ::

    sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
    sage: A = matrix(AA, [(1,0,1), (0,1,1), (-AA(2).sqrt(),-AA(3).sqrt(),1), 
    ....:                 (-AA(3).sqrt(),-AA(2).sqrt(),1)])
    sage: alg = StandardAlgorithm(A)
    sage: alg.run().R
    ((-0.4177376677004119?, 0.5822623322995881?, 0.4177376677004119?), 
     (-0.2411809548974793?, -0.2411809548974793?, 0.2411809548974793?), 
     (0.07665629029830300?, 0.07665629029830300?, 0.2411809548974793?), 
     (0.5822623322995881?, -0.4177376677004119?, 0.4177376677004119?))

REFERENCES:

..  Komei Fukuda , Alain Prodon:
    Double Description Method Revisited,    
    Combinatorics and Computer Science, volume 1120 of Lecture Notes
    in Computer Science, page 91-111. Springer (1996)
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
from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import matrix, identity_matrix



def test_random(d, n):
    """
    Test random collections of inequalities.

    INPUT:

    - ``d`` -- integer. The dimension.

    - ``n``  -- integer. The number of random inequalities to generate.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.double_description import test_random
        sage: P = test_random(5, 10)
        sage: P.run().verify()
    """
    from sage.matrix.constructor import random_matrix
    from sage.rings.all import QQ
    while True:
        A = random_matrix(QQ, n, d)
        if A.rank() == min(n, d) and not any(a == 0 for a in A.rows()):
            break
    return StandardAlgorithm(A)
    

class DoubleDescriptionPair(SageObject):

    def __init__(self, problem, A_rows, R_cols):
        self.problem = problem
        self.A = tuple(A_rows)
        self.R = tuple(R_cols)

    def _make_new(self, A_rows, R_cols):
        return self.__class__(self.problem, A_rows, R_cols)

    def _repr_(self):
        from sage.misc.ascii_art import ascii_art
        s = ascii_art('Double description pair (A, R) defined by')
        A = ascii_art(matrix(self.A))
        A._baseline = (len(self.A) / 2)
        A = ascii_art('A = ') + A
        R = ascii_art(matrix(self.R).transpose())
        if len(self.R) > 0:
            R._baseline = (len(self.R[0]) / 2)
        else:
            R._baseline = 0
        R = ascii_art('R = ') + R
        return str(s * (A + ascii_art(',   ') + R))

    def inner_product_matrix(self):
        """
        Return the inner product matrix between the rows of A and the columns of R

        OUTPUT:

        A matrix. There is one row for each row of A and one column
        for each column of R.
        
        EXAMPLES::
        
            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: alg = StandardAlgorithm(A)
            sage: DD, _ = alg.initial_pair()
            sage: DD.inner_product_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        result = []
        for a in self.A:
            line = []
            for r in self.R:
                line.append(a.inner_product(r))
            result.append(line)
        return matrix(self.problem.base_ring(), result)

    def cone(self):
        """
        Return the cone defined by `A`

        This method is for debugging only.

        OUTPUT:

        A cone.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: DD, _ = StandardAlgorithm(A).initial_pair()
            sage: DD.cone().Hrepresentation()
            (An inequality (-1, -1, 1) x + 0 >= 0, 
             An inequality (0, 1, 1) x + 0 >= 0, 
             An inequality (1, 0, 1) x + 0 >= 0)
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        if len(self.A) == 0:
            return Polyhedron(vertices=[[0]*self.problem.dim()])
        else:
            ieqs = [[0] + list(a) for a in self.A]
            return Polyhedron(ieqs=ieqs, base_ring=self.problem.base_ring())

    def verify(self):
        from sage.geometry.polyhedron.constructor import Polyhedron
        from sage.rings.all import QQ
        if self.problem.base_ring() != QQ:
            return
        A_cone = self.cone()
        R_cone = Polyhedron(vertices=[[0]*self.problem.dim()], 
                            rays=self.R, base_ring=self.problem.base_ring())
        assert A_cone == R_cone
        assert A_cone.n_inequalities() <= len(self.A)
        assert R_cone.n_rays() == len(self.R)

    def R_by_sign(self, a):
        """
        Classify the rays into those that are positive, zero, and negative on `a`.
        
        INPUT:

        - ``a`` -- vector. Coefficient vector of a homogeneous inequality.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: DD, _ = StandardAlgorithm(A).initial_pair()
            sage: DD.R_by_sign(vector([1,-1,0]))
            ([(2/3, -1/3, 1/3)], [(-1/3, -1/3, 1/3)], [(-1/3, 2/3, 1/3)])
            sage: DD.R_by_sign(vector([1,1,1]))
            ([(2/3, -1/3, 1/3), (-1/3, 2/3, 1/3)], [], [(-1/3, -1/3, 1/3)])
        """
        pos = []
        nul = []
        neg = []
        for r in self.R:
            sgn = a * r
            if sgn > 0:
                pos.append(r)
            elif sgn < 0:
                neg.append(r)
            else:
                nul.append(r)
        return pos, nul, neg

    @cached_method
    def zero_set(self, ray):
        """
        Return the zero set (active set) `Z(r)`
        
        INPUT:

        - ``ray`` -- a ray vector.

        OUTPUT:

        A tuple containing the inequality vectors that are zero on `r`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
            sage: DD, _ = Problem(A).initial_pair()
            sage: DD.zero_set(DD.R[0])
            ((0, 1, 1), (-1, -1, 1))
        """
        return tuple(a for a in self.A if a.inner_product(ray) == 0)

    def is_extremal(self, ray):
        """
        Test whether the ray is extremal

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(0,1,0), (1,0,0), (0,-1,1), (-1,0,1)])
            sage: DD = StandardAlgorithm(A).run()
            sage: DD.is_extremal(DD.R[0])
            True
        """
        A_Zray = matrix(self.problem.base_ring(), self.zero_set(ray))
        return A_Zray.rank() == self.problem.dim() - 1

    def are_adjacent(self, r1, r2):
        """
        Return whether the two rays are adjacent


        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(0,1,0), (1,0,0), (0,-1,1), (-1,0,1)])
            sage: DD = StandardAlgorithm(A).run()
            sage: DD.are_adjacent(DD.R[0], DD.R[1])
            True
            sage: DD.are_adjacent(DD.R[0], DD.R[2])
            True
            sage: DD.are_adjacent(DD.R[0], DD.R[3])
            False
        """
        Z_r1 = self.zero_set(r1)
        Z_r2 = self.zero_set(r2)
        Z_12 = set(Z_r1).intersection(Z_r2)
        A_Z12 = matrix(self.problem.base_ring(), list(Z_12))
        return A_Z12.rank() == self.problem.dim() - 2

    def dual(self):
        """
        Return the dual pair `(R^T, A^T)`

        EXAMPLES::
        
            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: A = matrix(QQ, [(0,1,0), (1,0,0), (0,-1,1), (-1,0,1)])
            sage: DD, _ = Problem(A).initial_pair()
            sage: DD
            Double description pair (A, R) defined by
                [ 0  1  0]        [0 1 0]
            A = [ 1  0  0],   R = [1 0 0]
                [ 0 -1  1]        [1 0 1]
            sage: DD.dual()
             Double description pair (A, R) defined by
                [0 1 1]        [ 0  1  0]
            A = [1 0 0],   R = [ 1  0 -1]
                [0 0 1]        [ 0  0  1]
        """
        return self._make_new(self.R, self.A)

    # def remove_unnecessary_inequalities(self):
    #     """
    #     Return a new double description pair with unnecessary inequalities
    #     removed.
    #
    #     EXAMPLES::
    #     """
    #     supporting_hyperplanes = []
    #     d = self.problem.dim()
    #     ring = self.problem.base_ring()
    #     for a in self.A:
    #         R_perp = [r for r in self.R if r.inner_product(a) == 0]
    #         if matrix(ring, R_perp).rank() == d - 1:
    #             supporting_hyperplanes.append(a)
    #     if len(supporting_hyperplanes) == len(self.A):
    #         return self
    #     else:
    #         return self._make_new(supporting_hyperplanes, self.R)




class Problem(SageObject):

    pair_class = DoubleDescriptionPair

    def __init__(self, A):
        """
        Base class for implementations of the double description algorithm

        INPUT:

        - ``A`` -- a matrix. The rows of the matrix are interpreted as
          homogeneous inequalities `Ax \geq 0`. Must have maximal rank.
        """
        assert A.rank() == A.ncols()    # implementation assumes maximal rank
        self._A = A
        self._field = A.base_ring().fraction_field()

    @cached_method
    def A(self):
        """
        Return the rows of the defining matrix `A`.
        """
        rows = [a.change_ring(self._field) for a in self._A.rows()]
        map(lambda a:a.set_immutable(), rows)
        return tuple(rows)

    @cached_method
    def A_matrix(self):
        """
        Return the defining matrix `A`.
        """
        return matrix(self.base_ring(), self.A())

    def base_ring(self):
        """
        Return the base field.
        """
        return self._field

    @cached_method
    def dim(self):
        return self._A.ncols()
        
    def _repr_(self):
        return '\n'.join(map(str, self.A()))

    def initial_pair(self):
        """
        Return an initial double description pair.

        INPUT:

        - ``pair_class`` -- subclass of
          :class:`DoubleDescriptionPair`.

        OUTPUT:

        A pair consisting of a :class:`DoubleDescriptionPair` instance
        and a tuple of unused inequalities.

        EXAMPLES::

            sage: A = matrix([(-1, 1), (-1, 2), (1/2, -1/2), (1/2, 2)])
            sage: from sage.geometry.polyhedron.double_description import Problem
            sage: DD, _ = Problem(A).initial_pair()
            sage: DD.verify()
        """
        pivot_rows = self.A_matrix().pivot_rows()
        A0 = [self.A()[pivot] for pivot in pivot_rows]
        Ac = [self.A()[i] for i in range(len(self.A())) if i not in pivot_rows]
        I = identity_matrix(self.base_ring(), self.dim())
        R = matrix(A0).solve_right(I)
        return self.pair_class(self, A0, R.columns()), list(Ac)



class StandardDoubleDescriptionPair(DoubleDescriptionPair):

    def add_inequality(self, a):
        """
        EXAMPLES::

            sage: A = matrix([(-1, 1, 0), (-1, 2, 1), (1/2, -1/2, -1)])
            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: DD, _ = StandardAlgorithm(A).initial_pair()
            sage: newDD = DD.add_inequality(vector([1,0,0]));  newDD
            Double description pair (A, R) defined by
                [  -1    1    0]        [   1    1    0    0]
            A = [  -1    2    1],   R = [   1    1    1    1]
                [ 1/2 -1/2   -1]        [   0   -1 -1/2   -2]
                [   1    0    0]    
        """
        from sage.combinat.cartesian_product import CartesianProduct
        R_pos, R_nul, R_neg = self.R_by_sign(a)
        if len(R_neg) == 0:
            return self
        R_new = []
        for rp, rn in CartesianProduct(R_pos, R_neg):
            if not self.are_adjacent(rp, rn):
                continue
            r = a.inner_product(rp) * rn - a.inner_product(rn) * rp
            r.set_immutable()
            R_new.append(r)
        R_new = tuple(R_pos + R_nul + R_new)
        A_new = self.A + (a,)
        return self._make_new(A_new, R_new)
    
    
class StandardAlgorithm(Problem):
    """
    Standard implementation of the double description algorithm
    """
    pair_class = StandardDoubleDescriptionPair

    def run(self, check=True):
        """
        OUTPUT:

        A double description pair `(A, R)` as a
        :class:`DoubleDescriptionPair`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
            sage: A = matrix(QQ, [(0,1,0), (1,0,0), (0,-1,1), (-1,0,1)])
            sage: StandardAlgorithm(A).run()
            Double description pair (A, R) defined by
                [ 0  1  0]        [0 0 1 1]
            A = [ 1  0  0],   R = [1 0 1 0]
                [ 0 -1  1]        [1 1 1 1]
                [-1  0  1]    
        """
        DD, remaining = self.initial_pair()
        for a in remaining:
            DD = DD.add_inequality(a)
            if check:
                DD.verify()
        return DD
