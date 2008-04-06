"""
Algebra orders

AUTHOR: David Kohel, 2005-09
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.ring import Ring
from sage.algebras.algebra import Algebra
from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import FreeModule, VectorSpace
from sage.algebras.algebra_element import AlgebraElement
from sage.algebras.algebra_order_element import AlgebraOrderElement

class AlgebraOrder_generic(Algebra):
    """
    An order in an algebra.
    """
    def __init__(self, A, R, gens = [], basis=None, rank = 0):
        Algebra.__init__(self, R)
        K = A.base_ring()
        self.__ambient_algebra = A
        self.__base_ring = R
        self.__gens = gens
        if rank != 0:
            self.__rank = rank
        else:
            self.__rank = A.dimension()
        n = self.__rank
        try:
            V = A.__vector_space
        except AttributeError:
            V = VectorSpace(K,n)
            A._Algebra__vector_space = V
        if basis != None:
            M = FreeModule(R,n)
            self.__algebra_basis_elements = basis
            vecs = [ x.vector() for x in self.__algebra_basis_elements ]
            # This syntax should change to V.submodule(R,vecs)
            self.__module = M.submodule(vecs)

    def __repr__(self):
        return "Ideal on generators %s over %s"%(
            self.__gens, self.__algebra)

    def __call__(self, x, check=True):
        if isinstance(x, AlgebraOrderElement) and x.parent() is self: return x
        return AlgebraOrderElement(self, x, check=check)

    def __contains__(self, x):
        if isinstance(x, AlgebraOrderElement) and x.parent() is self:
	    return True
        elif isinstance(x, AlgebraElement) \
               and x.parent() is self.ambient_algebra():
            M = self.inverse_embedding_matrix()
            R = self.base_ring()
            try:
                v = [ R(c) for c in list(x.vector() * M) ]
                return True
            except TypeError:
                return False
            # Stronger definition of 'in'...
            # return not (False in (c in R for c in list(x.vector()*M)))
        elif isinstance(x, RingElement) and x in self.base_ring():
            return True
        else:
            return False

    def ambient_algebra(self):
        """
        Return the ambient algebra of this algebra order
        """
        return self.__ambient_algebra

    def base_ring(self):
        """
        Return the base ring of this algebra order
        """
        return self.__base_ring

    def __compute_order_basis(self):
        """
        Compute a basis as an order.
        """
        A = self.__ambient_algebra
        R = self.__base_ring
        n = self.__rank
        V = A.vector_space()
        M = FreeModule(R,n)
        # syntax should change to:
        # N = V.submodule(R,[ V.gen(0) ])
        N = M.submodule([ V.gen(0) ])
        basis = [ A(1) ] + list(self.__gens)
        saturated = False
        while not saturated:
            saturated = True
            for x in self.__gens:
                for y in basis:
                    v = (x*y).vector()
                    if not v in N:
                        N += M.submodule([ v ])
                        saturated = False
                    v = (y*x).vector()
                    if not v in N:
                        N += M.submodule([ v ])
                        saturated = False
            basis = [ A(v.list()) for v in N.basis() ]
        self.__algebra_basis_elements = basis

    def basis(self):
        """
        Return the basis of this algebra order
        """
        try:
            return tuple([ self(x, check=False) for x in self.__algebra_basis_elements ])
        except AttributeError:
            self.__compute_order_basis()
        return tuple([ self(x, check=False) for x in self.__algebra_basis_elements ])

    def algebra_basis_elements(self):
        """
        Return the basis of this algebra order as a tuple
        """
        return tuple(self.__algebra_basis_elements)

    def embedding_matrix(self):
        r"""
        Return the $n\times n$ matrix whose rows are the vectors
        expressing the basis elements of this algebra order in terms
        of the ambient algebra's basis.
        """
        try:
            return self.__embedding_martix
        except AttributeError:
            n = self.__rank
            K = self.ambient_algebra().base_ring()
            M = MatrixSpace(K,n)
            self.__embedding_matrix = M([ x.ambient_algebra_element().vector() for x in self.basis() ])
        return self.__embedding_matrix

    def inverse_embedding_matrix(self):
        r"""
        Return the inverse of the $n\times n$ matrix whose rows are
        the vectors expressing the basis elements of this algebra
        order in terms of the ambient algebra's basis.
        """
        try:
            return self.__inverse_embedding_matrix
        except AttributeError:
            pass
        embedding_matrix = self.embedding_matrix()
        try:
            self.__inverse_embedding_matrix = embedding_matrix**-1
            return self.__inverse_embedding_matrix
        except ZeroDivisionError:
            raise AttributeError, "Basis for quaternion order is not of full rank."

    def module(self):
        r"""
        Convert this algebra order into a module over the base ring.
        """
        try:
            return self.__module
        except AttributeError:
            M = FreeModule(self.__base_ring,self.__rank)
            vecs = [ x.ambient_algebra_element().vector() for x in self.basis() ]
            # This syntax should change to V.submodule(R,vecs)
            self.__module = M.submodule(vecs)
        return self.__module

    def gen(self,i):
        r"""
        Return the i'th generator of this algebra order
        """
        n = len(self.__gens)
        if i < 0 or not i < n:
            raise IndexError, \
                  "Argument i (= %s) must be between 0 and %s."%(i, n-1)
        return self(self.__gens[int(i)])

    def gens(self):
        r"""
        Return the generators of this algebra order as a tuple
        """
        return tuple([ self(x) for x in self.__gens ])

    def ngens(self):
        r"""
        Return the number of generators of this algebra order
        """
        return len(self.__gens)

    def rank(self):
        r"""
        Return the rank of this algebra order
        """
        return self.__rank

