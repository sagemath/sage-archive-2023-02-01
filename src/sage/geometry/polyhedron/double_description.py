"""
Double Description Algorithm

EXAMPLES::

    sage: from sage.geometry.polyhedron.double_description import StandardAlgorithm
    sage: A = matrix(QQ, [(1,0,1), (0,1,1), (-1,-1,1)])
    sage: alg = StandardAlgorithm(A);  alg
    (1, 0, 1)
    (0, 1, 1)
    (-1, -1, 1)
    sage: DD = alg.initial_pair();  DD
"""

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method


class DoubleDescriptionPair(object):

    def __init__(self, A_rows, R_cols):
        self.A = A_rows
        self.R = R_cols

    def cone(self):
        from sage.geometry.polyhedron.constructor import Polyhedron
        return Polyhedron(ieqs=self.A)

    def _check(self):
        assert self.cone() == Polyhedron(rays=self.R)

    def R_by_sign(self, a):
        """
        INPUT:

        - ``a`` -- vector. Coefficient vector of a homogeneous inequality.
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

    def add_inequality(self, a):
        R_pos, R_nul, R_neg = self.R_by_sign(a)

    


class StandardAlgorithm(SageObject):

    def __init__(self, A):
        """
        Standard implementation of the double description algorithm

        INPUT:

        - ``A`` -- a matrix. The rows of the matrix are interpreted as
          homogeneous inequalities `Ax \geq 0`.
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

    def base_ring(self):
        """
        Return the base field.
        """
        return self._field
        
    def _repr_(self):
        return '\n'.join(map(str, self.A()))

    def initial_pair(self):
        """
        Return an initial double description pair.

        OUTPUT:

        :class:`DoubleDescriptionPair`.
        """
        a = self.A()[0]
        

    def run(self):
        """        
        OUTPUT:

        A double description pair `(A, R)` as a
        :class:`DoubleDescriptionPair`.
        """
    
