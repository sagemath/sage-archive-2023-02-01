r"""
Morphisms defined by a matrix.

A matrix morphism is a morphism that is defined by multiplication by a
matrix.  Elements of domain must either have a method \code{vector()}
that returns a vector that the defining matrix can hit from the left,
or be coercible into vector space of appropriate dimension.

EXAMPLES:
    sage: from sage.modules.matrix_morphism import MatrixMorphism, is_MatrixMorphism
    sage: V = Q^3
    sage: T = End(V)
    sage: M = MatrixSpace(QQ,3)
    sage: I = M.identity_matrix()
    sage: m = MatrixMorphism(T, I); m
    Morphism defined by the matrix
    [1 0 0]
    [0 1 0]
    [0 0 1]
    sage: is_MatrixMorphism(m)
    True
    sage: m.charpoly()
    x^3 - 3*x^2 + 3*x - 1
    sage: m.base_ring()
    Rational Field
    sage: m.det()
    1
    sage: m.fcp()
    (x - 1)^3
    sage: m.matrix()
    [1 0 0]
    [0 1 0]
    [0 0 1]
    sage: m.rank()
    3
    sage: m.trace()
    3

AUTHOR:
    - William Stein: initial versions
    - David Joyner (2005-12-17): added examples
    - William Stein (2005-01-07): added __reduce__
"""


import sage.categories.all
import sage.categories.homset
import sage.matrix.all as matrix
import sage.misc.misc as misc
import sage.modules.free_module as free_module
import sage.rings.coerce
from   sage.structure.all import Sequence

def is_MatrixMorphism(x):
    return isinstance(x, MatrixMorphism)

class MatrixMorphism(sage.categories.all.Morphism):
    def __init__(self, parent, A):
        """
        INPUT:
            parent -- a homspace
            A -- matrix

        EXAMPLES:
            sage: from sage.modules.matrix_morphism import MatrixMorphism
            sage: T = End(Q^3)
            sage: M = MatrixSpace(Q,3,3)
            sage: I = M.identity_matrix()
            sage: A = MatrixMorphism(T, I)
            sage: loads(A.dumps()) == A
            True
        """
        if not sage.categories.homset.is_Homset(parent):
            raise TypeError, "parent must be a Hom space"
        if not matrix.is_Matrix(A):
            A = matrix.MatrixSpace(parent.category().base_ring(), parent.domain().rank(), parent.codomain().rank())(A)
        R = A.base_ring()
        if A.nrows() != parent.domain().rank():
            raise ArithmeticError, "number of rows of matrix (=%s) must equal rank of domain (=%s)"%(
                A.nrows(), parent.domain().rank())
        if A.ncols() != parent.codomain().rank():
                raise ArithmeticError, "number of columns of matrix (=%s) must equal rank of codomain (=%s)"%(
                    A.ncols(), parent.codomain().rank())
        self.__matrix = A
        sage.categories.all.Morphism.__init__(self, parent)

    def _repr_(self):
        if max(self.__matrix.nrows(),self.__matrix.ncols()) > 5:
            mat = "(not printing %s x %s matrix)"%(self.__matrix.nrows(), self.__matrix.ncols())
        else:
            mat = str(self.__matrix)
        return "Morphism defined by the matrix\n%s"%mat

    def __cmp__(self, other):
        if not isinstance(other, MatrixMorphism) or self.parent() != other.parent():
            return sage.rings.coerce.cmp(self, other)
        return cmp(self.__matrix, other.__matrix)

    def __call__(self, x):
        try:
            if x.parent() != self.domain():
                x = self.domain()(x)
        except AttributeError, TypeError:
            raise TypeError, "%s must be coercible into %s"%(x,self.domain())
        if self.domain().is_ambient():
            x = x.element()
        else:
            x = self.domain().coordinate_vector(x)
        v = x*self.matrix()
        C = self.codomain()
        if C.is_ambient():
            return C(v)
        return C.linear_combination_of_basis(v)

    def __invert__(self):
        if self.nrows() != self.ncols():
            raise ZeroDivisionError, "Inverse of morphism not defined."%self
        try:
            B = ~self.__matrix
        except ZeroDivisionError:
            raise ZeroDivisionError, "Inverse does not exist."
        return self.parent().reversed()(B)

    def _mul_function(self, other):
        return other.__matrix * self.__matrix

    def _add_function(self, other):
        return self.__matrix + other.__matrix

    def _sub_function(self, other):
        return self.__matrix - other.__matrix

    def __rmul__(self, left):
        R = self.base_ring()
        return self.parent()(R(left) * self.__matrix)

    def __mul__(self, right):
        if not isinstance(right, MatrixMorphism):
            R = self.base_ring()
            return self.parent()(self.__matrix * R(right))
        return sage.categories.all.Morphism.__mul__(self, right)

    def base_ring(self):
        return self.domain().base_ring()

    def charpoly(self):
        if not self.is_endomorphism():
            raise ArithmeticError, "charpoly only defined for endomorphisms " +\
                    "(i.e., domain = range)"
        return self.__matrix.charpoly()

    def decomposition(self, is_diagonalizable=False):
        if not self.is_endomorphism():
            raise ArithmeticError, "Matrix morphism must be an endomorphism."
        D = self.domain()
        E = self.__matrix.decomposition(is_diagonalizable=is_diagonalizable)
        if D.is_ambient():
            return Sequence([D.submodule(V) for V, _ in E], cr=True, check=False)
        else:
            B = D.basis_matrix()
            return Sequence([D.submodule((V.basis_matrix() * B).row_space()) for V, _ in E],
                            cr=True, check=False)

    def det(self):
        """
        Return the determinant of this endomorphism.
        """
        if not self.is_endomorphism():
            raise ArithmeticError, "Matrix morphism must be an endomorphism."
        return self.matrix().determinant()

    def fcp(self):
        """
        Return the factorization of the characteristic polynomial.
        """
        return self.charpoly().factor()

    def kernel(self):
        V = self.matrix().kernel()
        return self.domain().submodule(V)

    def image(self):
        V = self.matrix().image()
        return self.codomain().submodule(V)

    def matrix(self):
        return self.__matrix

    def rank(self):
        return self.__matrix.rank()

    def restrict_domain(self, sub):
        """
        Restrict this matrix morphism to a subspace sub of the domain.
        The subspace sub should have a basis() method and elements of the basis
        should be coercible into domain.

        The resulting morphism has the same codomain as before, but
        a new domain.
        """
        D  = self.domain()
        B  = sub.basis()
        ims= sum([(self(D(b)).coordinate_vector()).list() for b in B],[])

        MS = matrix.MatrixSpace(self.base_ring(), len(B), self.codomain().rank())
        H = sage.categories.homset.Hom(sub, self.codomain(), self.category())
        return H(MS(ims))

    def restrict(self, sub):
        """
        Restrict this matrix morphism to a subspace sub of the domain.

        The codomain and domain of the resulting matrix are both sub.
        """
        if not self.is_endomorphism():
            raise ArithmeticError, "matrix morphism must be an endomorphism"
        A = self.matrix().restrict(sub.free_module())
        H = sage.categories.homset.End(sub, self.domain().category())
        return H(A)

    def trace(self):
        return self.__matrix.trace()
