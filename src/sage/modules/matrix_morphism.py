r"""
Morphisms defined by a matrix.

A matrix morphism is a morphism that is defined by multiplication by a
matrix.  Elements of domain must either have a method \code{vector()}
that returns a vector that the defining matrix can hit from the left,
or be coercible into vector space of appropriate dimension.

EXAMPLES:
    sage: from sage.modules.matrix_morphism import MatrixMorphism, is_MatrixMorphism
    sage: V = QQ^3
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
    sage: m.charpoly('x')
    x^3 - 3*x^2 + 3*x - 1
    sage: m.base_ring()
    Rational Field
    sage: m.det()
    1
    sage: m.fcp('x')
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
            sage: T = End(QQ^3)
            sage: M = MatrixSpace(QQ,3)
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
            raise ArithmeticError, "number of rows of matrix (=%s) must equal rank of domain (=%s)"%(A.nrows(), parent.domain().rank())
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
        return cmp(self.__matrix, other.__matrix)

    def __call__(self, x):
        """
        Evaluate this matrix morphism at an element that can be
        coerced into the domain.

        EXAMPLES:
            sage: V = QQ^3; W = QQ^2
            sage: H = Hom(V, W); H
            Set of Morphisms from Vector space of dimension 3 over Rational Field to Vector space of dimension 2 over Rational Field in Category of vector spaces over Rational Field
            sage: phi = H(range(6)); phi
            Free module morphism defined by the matrix
            [0 1]
            [2 3]
            [4 5]
            Domain: Vector space of dimension 3 over Rational Field
            Codomain: Vector space of dimension 2 over Rational Field
            sage: phi(V.0)
            (0, 1)
            sage: phi([1,2,3])
            (16, 22)
            sage: phi(5)
            Traceback (most recent call last):
            ...
            TypeError: 5 must be coercible into Vector space of dimension 3 over Rational Field
            sage: phi([1,1])
            Traceback (most recent call last):
            ...
            TypeError: [1, 1] must be coercible into Vector space of dimension 3 over Rational Field
        """
        try:
            if not hasattr(x, 'parent') or x.parent() != self.domain():
                x = self.domain()(x)
        except TypeError:
            raise TypeError, "%s must be coercible into %s"%(x,self.domain())
        if self.domain().is_ambient():
            x = x.element()
        else:
            x = self.domain().coordinate_vector(x)
        v = x*self.matrix()
        C = self.codomain()
        if C.is_ambient():
            return C(v)
        return C(C.linear_combination_of_basis(v), check=False)

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

    def charpoly(self, var):
        if not self.is_endomorphism():
            raise ArithmeticError, "charpoly only defined for endomorphisms " +\
                    "(i.e., domain = range)"
        return self.__matrix.charpoly(var)

    def decomposition(self, is_diagonalizable=False):
        if not self.is_endomorphism():
            raise ArithmeticError, "Matrix morphism must be an endomorphism."
        D = self.domain()
        E = self.__matrix.decomposition(is_diagonalizable=is_diagonalizable)
        if D.is_ambient():
            return Sequence([D.submodule(V, check=False) for V, _ in E], cr=True, check=False)
        else:
            B = D.basis_matrix()
            return Sequence([D.submodule((V.basis_matrix() * B).row_space(), check=False) for V, _ in E],
                            cr=True, check=False)

    def det(self):
        """
        Return the determinant of this endomorphism.
        """
        if not self.is_endomorphism():
            raise ArithmeticError, "Matrix morphism must be an endomorphism."
        return self.matrix().determinant()

    def fcp(self, var='x'):
        """
        Return the factorization of the characteristic polynomial.
        """
        return self.charpoly(var).factor()

    def kernel(self):
        """
        Compute the kernel of this matrix.

        EXAMPLES:
            sage: V = VectorSpace(QQ,3)
            sage: id = V.Hom(V)(identity_matrix(QQ,3))
            sage: null = V.Hom(V)(0*identity_matrix(QQ,3))
            sage: id.kernel()
            Vector space of degree 3 and dimension 0 over Rational Field
            Basis matrix:
            []
            sage: phi = V.Hom(V)(matrix(QQ,3,range(9)))
            sage: phi.kernel()
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [ 1 -2  1]
        """
        V = self.matrix().kernel()
        D = self.domain()
        if not D.is_ambient():
            # Transform V to ambient space
            # This is a matrix multiply:  we take the linear combinations of the basis for
            # D given by the elements of the basis for V.
            B = V.basis_matrix() * D.basis_matrix()
            V = B.row_space()
        return self.domain().submodule(V, check=False)

    def image(self):
        V = self.matrix().image()
        D = self.codomain()
        if not D.is_ambient():
            # Transform V to ambient space
            # This is a matrix multiply:  we take the linear combinations of the basis for
            # D given by the elements of the basis for V.
            B = V.basis_matrix() * D.basis_matrix()
            V = B.row_space()
        return self.codomain().submodule(V, check=False)

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
        H = sub.Hom(self.codomain(), sub.category())
        return H(MS(ims))

        #D  = self.domain()
	#C  = self.codomain()
	#M  = self.matrix()
	#Mr = M.restrict_domain(sub)
	#return sub.Hom(C)(Mr)

    def restrict_codomain(self, sub):
        """
        Restrict this matrix morphism to a subspace sub of the codomain.

        The resulting morphism has the same domain as before, but
        a new codomain.
	"""
	A = self.matrix().restrict_codomain(sub.free_module())
        H = sage.categories.homset.Hom(self.domain(), sub, self.domain().category())
        return H(A)

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
