r"""
Homspaces between free modules

EXAMPLES: We create `\mathrm{End}(\ZZ^2)` and compute a
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

We create `\mathrm{Hom}(\QQ^3, \QQ^2)` and
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
    [0 0]...

TESTS::

    sage: H = Hom(QQ^2, QQ^1)
    sage: loads(dumps(H)) == H
    True

See trac 5886::

    sage: V = (QQ^2).span_of_basis([[1,2],[3,4]])
    sage: V.hom([V.0, V.1])
    Free module morphism defined by the matrix
    [1 0]
    [0 1]...

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
    """
    Return True if x is a Free module homspace.

    EXAMPLES::

        sage: H = Hom(QQ^3, QQ^2)
        sage: sage.modules.free_module_homspace.is_FreeModuleHomspace(H)
        True
        sage: sage.modules.free_module_homspace.is_FreeModuleHomspace(2)
        False
    """
    return isinstance(x, FreeModuleHomspace)

class FreeModuleHomspace(sage.categories.homset.HomsetWithBase):
    def __call__(self, A, check=True):
        """
        INPUT:
            A -- either a matrix or a list/tuple of images of generators
            check -- bool (default: True)

        If A is a matrix, then it is the matrix of this linear
        transformation, with respect to the basis for the domain and
        codomain.  Thus the identity matrix always defines the
        identity morphism.

        EXAMPLES::
            sage: V = (QQ^3).span_of_basis([[1,1,0],[1,0,2]])
            sage: H = V.Hom(V); H
            Set of Morphisms from ...
            sage: H([V.0,V.1])                    # indirect doctest
            Free module morphism defined by the matrix
            [1 0]
            [0 1]...
            sage: phi = H([V.1,V.0]); phi
            Free module morphism defined by the matrix
            [0 1]
            [1 0]...
            sage: phi(V.1) == V.0
            True
            sage: phi(V.0) == V.1
            True
        """
        if not matrix.is_Matrix(A):
            # Compute the matrix of the morphism that sends the
            # generators of the domain to the elements of A.
            C = self.codomain()
            try:
                v = [C(a) for a in A]
                A = matrix.matrix([C.coordinates(a) for a in v])
            except TypeError:
                pass
        return free_module_morphism.FreeModuleMorphism(self, A)

    def _matrix_space(self):
        """
        Return underlying matrix space that contains the matrices that define
        the homomorphisms in this free module homspace.

        OUTPUT:
            - matrix space

        EXAMPLES::

            sage: H = Hom(QQ^3, QQ^2)
            sage: H._matrix_space()
            Full MatrixSpace of 3 by 2 dense matrices over Rational Field
        """
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

        OUTPUT:
            - tuple

        EXAMPLES::

            sage: H = Hom(QQ^2, QQ^1)
            sage: H.basis()
            (Free module morphism defined by the matrix
            [1]
            [0]
            Domain: Vector space of dimension 2 over Rational Field
            Codomain: Vector space of dimension 1 over Rational Field,
             Free module morphism defined by the matrix
            [0]
            [1]
            Domain: Vector space of dimension 2 over Rational Field
            Codomain: Vector space of dimension 1 over Rational Field)
        """
        try:
            return self.__basis
        except AttributeError:
            M = self._matrix_space()
            B = M.basis()
            self.__basis = tuple([self(x) for x in B])
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



