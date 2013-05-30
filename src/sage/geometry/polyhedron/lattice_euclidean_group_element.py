"""
Lattice Euclidean Group Elements

The clases here are used to return particular isomorphisms of
:class:`PPL lattice
polytopes<sage.geometry.polyhedron.ppl_lattice_polytope.LatticePolytope_PPL_class>`.
"""
########################################################################
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.modules.all import vector
from sage.matrix.constructor import matrix



########################################################################
class LatticePolytopeError(Exception):
    """
    Base class for errors from lattice polytopes
    """
    pass


########################################################################
class LatticePolytopesNotIsomorphicError(LatticePolytopeError):
    """
    Raised when two lattice polytopes are not isomorphic.
    """
    pass


########################################################################
class LatticePolytopeNoEmbeddingError(LatticePolytopeError):
    """
    Raised when no embedding of the desired kind can be found.
    """
    pass




########################################################################
class LatticeEuclideanGroupElement(SageObject):

    def __init__(self, A, b):
        """
        An element of the lattice Euclidean group.

        Note that this is just intended as a container for results from
        LatticePolytope_PPL. There is no group-theoretic functionality to
        speak of.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope \
            ....:     import LatticePolytope_PPL, C_Polyhedron
            sage: from sage.geometry.polyhedron.lattice_euclidean_group_element \
            ....:     import LatticeEuclideanGroupElement
            sage: M = LatticeEuclideanGroupElement([[1,2],[2,3],[-1,2]], [1,2,3])
            sage: M
            The map A*x+b with A=
            [ 1  2]
            [ 2  3]
            [-1  2]
            b =
            (1, 2, 3)
            sage: M._A
            [ 1  2]
            [ 2  3]
            [-1  2]
            sage: M._b
            (1, 2, 3)
            sage: M(vector([0,0]))
            (1, 2, 3)
            sage: M( LatticePolytope_PPL((0,0),(1,0),(0,1)) )
            A 2-dimensional lattice polytope in ZZ^3 with 3 vertices
            sage: _.vertices()
            ((1, 2, 3), (2, 4, 2), (3, 5, 5))
        """
        self._A = matrix(ZZ, A)
        self._b = vector(ZZ, b)
        assert self._A.nrows() == self._b.degree()

    def __call__(self, x):
        """
        Return the image of ``x``

        INPUT:

        - ``x`` -- a vector or lattice polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope \
            ....:     import LatticePolytope_PPL, C_Polyhedron
            sage: from sage.geometry.polyhedron.lattice_euclidean_group_element \
            ....:     import LatticeEuclideanGroupElement
            sage: M = LatticeEuclideanGroupElement([[1,2],[2,3],[-1,2]], [1,2,3])
            sage: M(vector(ZZ, [11,13]))
            (38, 63, 18)
            sage: M( LatticePolytope_PPL((0,0),(1,0),(0,1)) )
            A 2-dimensional lattice polytope in ZZ^3 with 3 vertices
        """
        from sage.geometry.polyhedron.ppl_lattice_polytope import (
            LatticePolytope_PPL, LatticePolytope_PPL_class )
        if isinstance(x, LatticePolytope_PPL_class):
            if x.is_empty():
                from sage.libs.ppl import C_Polyhedron
                return LatticePolytope_PPL(C_Polyhedron(self._b.degree(), 'empty'))
            return LatticePolytope_PPL(*[self(v) for v in x.vertices()])
            pass
        v = self._A*x+self._b
        v.set_immutable()

        return v

    def _repr_(self):
        """
        Return a string representation

        EXAMPLES::

            sage: from sage.geometry.polyhedron.lattice_euclidean_group_element \
            ....:     import LatticeEuclideanGroupElement
            sage: M = LatticeEuclideanGroupElement([[1,2],[2,3],[-1,2]], [1,2,3])
            sage: M._repr_()
            'The map A*x+b with A=\n[ 1  2]\n[ 2  3]\n[-1  2]\nb = \n(1, 2, 3)'
        """
        s = 'The map A*x+b with A=\n'+str(self._A)
        s += '\nb = \n'+str(self._b)
        return s


