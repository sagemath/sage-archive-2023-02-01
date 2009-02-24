r"""
Homspaces between free modules

EXAMPLES: We create `\mathrm{End}(\mathbb{Z}^2)` and compute a
basis.

::

    sage: M = FreeModule(IntegerRing(),2)
    sage: E = End(M)
    sage: B = E.basis()
    sage: len(B)
    4
    sage: B[0]
    Free module morphism defined by the matrix
    [1 0]
    [0 0]
    Domain: Ambient free module of rank 2 over the principal ideal domain ...
    Codomain: Ambient free module of rank 2 over the principal ideal domain ...

We create `\mathrm{Hom}(\mathbb{Q}^3, \mathbb{Q}^2)` and
compute a basis.

::

    sage: V3 = VectorSpace(RationalField(),3)
    sage: V2 = VectorSpace(RationalField(),2)
    sage: H = Hom(V3,V2)
    sage: H
    Set of Morphisms from Vector space of dimension 3 over Rational Field
    to Vector space of dimension 2 over Rational Field in Category of
    vector spaces over Rational Field
    sage: B = H.basis()
    sage: len(B)
    6
    sage: B[0]
    Free module morphism defined by the matrix
    [1 0]
    [0 0]
    [0 0]
    Domain: Vector space of dimension 3 over Rational Field
    Codomain: Vector space of dimension 2 over Rational Field
"""

#*****************************************************************************
#  Copyright (C) 2005 William Stein <wstein@gmail.com>
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

import sage.categories.homset
import sage.modules.free_module_morphism
import sage.matrix.all as matrix
import free_module_morphism

from matrix_morphism import MatrixMorphism

def is_FreeModuleHomspace(x):
    return isinstance(x, FreeModuleHomspace)

class FreeModuleHomspace(sage.categories.homset.HomsetWithBase):
    def __call__(self, A, check=True):
        return free_module_morphism.FreeModuleMorphism(self, A)

    def _coerce_impl(self, x):
        M = self._matrix_space()
        if isinstance(x,MatrixMorphism):
            x = x.matrix()
        return self._coerce_try(x, [M])

    def _matrix_space(self):
        try:
            return self.__matrix_space
        except AttributeError:
            R = self.domain().base_ring()
            M = matrix.MatrixSpace(R, self.domain().rank(), self.codomain().rank())
            self.__matrix_space = M
            return M

    def basis(self):
        """
        Return a basis for this space of free module homomorphisms.
        """
        try:
            return self.__basis
        except AttributeError:
            M = self._matrix_space()
            B = M.basis()
            self.__basis = [self(x) for x in B]
            return self.__basis

    def identity(self):
       r"""
       Return identity morphism in an endomorphism ring.

       EXAMPLE::

           sage: V=VectorSpace(QQ,5)
           sage: H=V.Hom(V)
           sage: H.identity()
           Free module morphism defined by the matrix
           [1 0 0 0 0]
           [0 1 0 0 0]
           [0 0 1 0 0]
           [0 0 0 1 0]
           [0 0 0 0 1]
           Domain: Vector space of dimension 5 over Rational Field
           Codomain: Vector space of dimension 5 over Rational Field
       """
       if self.is_endomorphism_set():
           return self(matrix.identity_matrix(self.base_ring(),self.domain().rank()))
       else:
           raise TypeError, "Identity map only defined for endomorphisms. Try natural_map() instead."



