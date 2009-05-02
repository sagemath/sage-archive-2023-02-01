r"""
Morphisms defined by a matrix.

A matrix morphism is a morphism that is defined by multiplication
by a matrix. Elements of domain must either have a method
``vector()`` that returns a vector that the defining
matrix can hit from the left, or be coercible into vector space of
appropriate dimension.

EXAMPLES::

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

- Craig Citro (2008-03-18): refactored MatrixMorphism class
"""


import sage.categories.all
import sage.categories.homset
import sage.matrix.all as matrix
import sage.misc.misc as misc
import sage.modules.free_module as free_module
from   sage.structure.all import Sequence

def is_MatrixMorphism(x):
    """
    Return True if x is a Matrix morphism of free modules.

    EXAMPLES::

        sage: V = ZZ^2; phi = V.hom([3*V.0, 2*V.1])
        sage: sage.modules.matrix_morphism.is_MatrixMorphism(phi)
        True
        sage: sage.modules.matrix_morphism.is_MatrixMorphism(3)
        False
    """
    return isinstance(x, MatrixMorphism_abstract)

class MatrixMorphism_abstract(sage.categories.all.Morphism):
    def __init__(self, parent):
        """
        INPUT:


        -  ``parent`` - a homspace

        -  ``A`` - matrix


        EXAMPLES::

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
        sage.categories.all.Morphism.__init__(self, parent)

    def __cmp__(self, other):
        """
        Compare two matrix morphisms.

        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([3*V.0, 2*V.1])
            sage: phi == 3
            False
            sage: phi == phi
            True
        """
        return cmp(self.matrix(), other.matrix())

    def __call__(self, x):
        """
        Evaluate this matrix morphism at an element that can be coerced
        into the domain.

        EXAMPLES::

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
        """
        Invert this matrix morphism.

        EXAMPLES::

            sage: V = QQ^2; phi = V.hom([3*V.0, 2*V.1])
            sage: phi^(-1)
            Free module morphism defined by the matrix
            [1/3   0]
            [  0 1/2]
            Domain: Vector space of dimension 2 over Rational Field
            Codomain: Vector space of dimension 2 over Rational Field

        Check that a certain non-invertible morphism isn't invertible::

            sage: V = ZZ^2; phi = V.hom([3*V.0, 2*V.1])
            sage: phi^(-1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: matrix morphism not invertible
        """
        try:
            B = ~(self.matrix())
        except ZeroDivisionError:
            raise ZeroDivisionError, "matrix morphism not invertible"
        try:
            return self.parent().reversed()(B)
        except TypeError:
            raise ZeroDivisionError, "matrix morphism not invertible"

    def __rmul__(self, left):
        """
        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([V.0+V.1, 2*V.1])
            sage: 2*phi
            Free module morphism defined by the matrix
            [2 2]
            [0 4]...
            sage: phi*2
            Free module morphism defined by the matrix
            [2 2]
            [0 4]...
        """
        R = self.base_ring()
        return self.parent()(R(left) * self.matrix())

    def __mul__(self, right):
        """
        Composition of morphisms, denoted by \*.

        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([V.0+V.1, 2*V.1])
            sage: phi*phi
            Free module morphism defined by the matrix
            [1 3]
            [0 4]
            Domain: Ambient free module of rank 2 over the principal ideal domain ...
            Codomain: Ambient free module of rank 2 over the principal ideal domain ...

            sage: V = QQ^3
            sage: E = V.endomorphism_ring()
            sage: phi = E(Matrix(QQ,3,range(9))) ; phi
            Free module morphism defined by the matrix
            [0 1 2]
            [3 4 5]
            [6 7 8]
            Domain: Vector space of dimension 3 over Rational Field
            Codomain: Vector space of dimension 3 over Rational Field
            sage: phi*phi
            Free module morphism defined by the matrix
            [ 15  18  21]
            [ 42  54  66]
            [ 69  90 111]
            Domain: Vector space of dimension 3 over Rational Field
            Codomain: Vector space of dimension 3 over Rational Field
            sage: phi.matrix()**2
            [ 15  18  21]
            [ 42  54  66]
            [ 69  90 111]

        ::

            sage: W = QQ**4
            sage: E_VW = V.Hom(W)
            sage: psi = E_VW(Matrix(QQ,3,4,range(12))) ; psi
            Free module morphism defined by the matrix
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            Domain: Vector space of dimension 3 over Rational Field
            Codomain: Vector space of dimension 4 over Rational Field
            sage: psi*phi
            Free module morphism defined by the matrix
            [ 20  23  26  29]
            [ 56  68  80  92]
            [ 92 113 134 155]
            Domain: Vector space of dimension 3 over Rational Field
            Codomain: Vector space of dimension 4 over Rational Field
            sage: phi*psi
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 3 by 4 dense matrices over Rational Field' and 'Full MatrixSpace of 3 by 3 dense matrices over Rational Field'
            sage: phi.matrix()*psi.matrix()
            [ 20  23  26  29]
            [ 56  68  80  92]
            [ 92 113 134 155]
        """
        if not isinstance(right, MatrixMorphism):
            R = self.base_ring()
            return self.parent()(self.matrix() * R(right))
        M = right.matrix() * self.matrix()
        return right.domain().Hom(self.codomain())(M)

    def __add__(self, right):
        """
        Sum of morphisms, denoted by +.

        EXAMPLES::

            sage: phi = (ZZ**2).endomorphism_ring()(Matrix(ZZ,2,[2..5])) ; phi
            Free module morphism defined by the matrix
            [2 3]
            [4 5]
            Domain: Ambient free module of rank 2 over the principal ideal domain ...
            Codomain: Ambient free module of rank 2 over the principal ideal domain ...
            sage: phi + 3
            Free module morphism defined by the matrix
            [5 3]
            [4 8]
            Domain: Ambient free module of rank 2 over the principal ideal domain ...
            Codomain: Ambient free module of rank 2 over the principal ideal domain ...
            sage: phi + phi
            Free module morphism defined by the matrix
            [ 4  6]
            [ 8 10]
            Domain: Ambient free module of rank 2 over the principal ideal domain ...
            Codomain: Ambient free module of rank 2 over the principal ideal domain ...
            sage: psi = (ZZ**3).endomorphism_ring()(Matrix(ZZ,3,[22..30])) ; psi
            Free module morphism defined by the matrix
            [22 23 24]
            [25 26 27]
            [28 29 30]
            Domain: Ambient free module of rank 3 over the principal ideal domain ...
            Codomain: Ambient free module of rank 3 over the principal ideal domain ...
            sage: phi + psi
            Traceback (most recent call last):
            ...
            TypeError: entries has the wrong length
        """
        # TODO: move over to any coercion model!
        if not isinstance(right, MatrixMorphism):
            R = self.base_ring()
            return self.parent()(self.matrix() + R(right))
        if not right.parent() == self.parent():
            right = self.parent()(right)
        M = self.matrix() + right.matrix()
        return self.domain().Hom(right.codomain())(M)

    def __neg__(self):
        """
        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([V.0+V.1, 2*V.1])
            sage: -phi
            Free module morphism defined by the matrix
            [-1 -1]
            [ 0 -2]...
        """
        return self.parent()(-self.matrix())

    def __sub__(self, other):
        """
        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([V.0+V.1, 2*V.1])
            sage: phi - phi
            Free module morphism defined by the matrix
            [0 0]
            [0 0]...
        """
        # TODO: move over to any coercion model!
        if not isinstance(other, MatrixMorphism):
            R = self.base_ring()
            return self.parent()(self.matrix() - R(other))
        if not other.parent() == self.parent():
            other = self.parent()(other)
        return self.parent()(self.matrix() - other.matrix())

    def base_ring(self):
        """
        Return the base ring of self, that is, the ring over which self is
        given by a matrix.

        EXAMPLES::

            sage: sage.modules.matrix_morphism.MatrixMorphism((ZZ**2).endomorphism_ring(), Matrix(ZZ,2,[3..6])).base_ring()
            Integer Ring
        """
        return self.domain().base_ring()

    def charpoly(self, var='x'):
        """
        Return the characteristic polynomial of this endomorphism.

        INPUT:
            - var -- variable

        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([V.0+V.1, 2*V.1])
            sage: phi.charpoly()
            x^2 - 3*x + 2
            sage: phi.matrix().charpoly()
            x^2 - 3*x + 2
            sage: phi.charpoly('T')
            T^2 - 3*T + 2
        """
        if not self.is_endomorphism():
            raise ArithmeticError, "charpoly only defined for endomorphisms " +\
                    "(i.e., domain = range)"
        return self.matrix().charpoly(var)

    def decomposition(self, *args, **kwds):
        """
        Return decomposition of this endomorphism, i.e., sequence of
        subspaces obtained by finding invariant subspaces of self.

        See the documentation for self.matrix().decomposition for more
        details.  All inputs to this function are passed onto the
        matrix one.

        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([V.0+V.1, 2*V.1])
            sage: phi.decomposition()
            [
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1],
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 1 -1]
            ]
        """
        if not self.is_endomorphism():
            raise ArithmeticError, "Matrix morphism must be an endomorphism."
        D = self.domain()
        E = self.matrix().decomposition(*args,**kwds)
        if D.is_ambient():
            return Sequence([D.submodule(V, check=False) for V, _ in E],
                            cr=True, check=False)
        else:
            B = D.basis_matrix()
            return Sequence([D.submodule((V.basis_matrix() * B).row_space(),
                                         check=False) for V, _ in E],
                            cr=True, check=False)

    def det(self):
        """
        Return the determinant of this endomorphism.

        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([V.0+V.1, 2*V.1])
            sage: phi.det()
            2
        """
        if not self.is_endomorphism():
            raise ArithmeticError, "Matrix morphism must be an endomorphism."
        return self.matrix().determinant()

    def fcp(self, var='x'):
        """
        Return the factorization of the characteristic polynomial.

        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([V.0+V.1, 2*V.1])
            sage: phi.fcp()
            (x - 2) * (x - 1)
            sage: phi.fcp('T')
            (T - 2) * (T - 1)
        """
        return self.charpoly(var).factor()

    def kernel(self):
        """
        Compute the kernel of this morphism.

        EXAMPLES::

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
            sage: hom(CC^2, CC^2, 1).kernel()
            Vector space of degree 2 and dimension 0 over Complex Field with 53 bits of precision
            Basis matrix:
            []
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
        """
        Compute the image of this morphism.

        EXAMPLES::

            sage: V = VectorSpace(QQ,3)
            sage: phi = V.Hom(V)(range(9))
            sage: phi.image()
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: hom(GF(7)^3, GF(7)^2, 0).image()
            Vector space of degree 2 and dimension 0 over Finite Field of size 7
            Basis matrix:
            []


        Compute the image of the identity map on a ZZ-submodule::

            sage: V = (ZZ^2).span([[1,2],[3,4]])
            sage: phi = V.Hom(V)(identity_matrix(ZZ,2))
            sage: phi(V.0) == V.0
            True
            sage: phi(V.1) == V.1
            True
            sage: phi.image()
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 2]
            sage: phi.image() == V
            True
        """
        V = self.matrix().image()
        D = self.codomain()
        if not D.is_ambient():
            # Transform V to ambient space
            # This is a matrix multiply:  we take the linear combinations of the basis for
            # D given by the elements of the basis for V.
            B = V.basis_matrix() * D.basis_matrix()
            V = B.row_space(self.domain().base_ring())
        return self.codomain().submodule(V, check=False)

    def matrix(self):
        """
        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom(V.basis())
            sage: phi.matrix()
            [1 0]
            [0 1]
            sage: sage.modules.matrix_morphism.MatrixMorphism_abstract.matrix(phi)
            Traceback (most recent call last):
            ...
            NotImplementedError: this method must be overridden in the extension class
        """
        raise NotImplementedError, "this method must be overridden in the extension class"

    def rank(self):
        """
        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom(V.basis())
            sage: phi.rank()
            2
            sage: V = ZZ^2; phi = V.hom([V.0, V.0])
            sage: phi.rank()
            1
        """
        return self.matrix().rank()

    def restrict_domain(self, sub):
        """
        Restrict this matrix morphism to a subspace sub of the domain. The
        subspace sub should have a basis() method and elements of the basis
        should be coercible into domain.

        The resulting morphism has the same codomain as before, but a new
        domain.

        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([3*V.0, 2*V.1])
            sage: phi.restrict_domain(V.span([V.0]))
            Free module morphism defined by the matrix
            [3 0]
            Domain: Free module of degree 2 and rank 1 over Integer Ring
            Echelon ...
            Codomain: Ambient free module of rank 2 over the principal ideal domain ...
            sage: phi.restrict_domain(V.span([V.1]))
            Free module morphism defined by the matrix
            [0 2]...
        """
        H = sub.Hom(self.codomain())
        return H(self.matrix().restrict_domain(sub))

    def restrict_codomain(self, sub):
        """
        Restrict this matrix morphism to a subspace sub of the codomain.

        The resulting morphism has the same domain as before, but a new
        codomain.

        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([4*(V.0+V.1),0])
            sage: W = V.span([2*(V.0+V.1)])
            sage: phi
            Free module morphism defined by the matrix
            [4 4]
            [0 0]
            Domain: Ambient free module of rank 2 over the principal ideal domain ...
            Codomain: Ambient free module of rank 2 over the principal ideal domain ...
            sage: psi = phi.restrict_codomain(W); psi
            Free module morphism defined by the matrix
            [2]
            [0]
            Domain: Ambient free module of rank 2 over the principal ideal domain ...
            Codomain: Free module of degree 2 and rank 1 over Integer Ring
            Echelon ...
        """
        H = self.domain().Hom(sub)
        return H(self.matrix().restrict_codomain(sub.free_module()))

    def restrict(self, sub):
        """
        Restrict this matrix morphism to a subspace sub of the domain.

        The codomain and domain of the resulting matrix are both sub.

        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([3*V.0, 2*V.1])
            sage: phi.restrict(V.span([V.0]))
            Free module morphism defined by the matrix
            [3]
            Domain: Free module of degree 2 and rank 1 over Integer Ring
            Echelon ...
            Codomain: Free module of degree 2 and rank 1 over Integer Ring
            Echelon ...
        """
        if not self.is_endomorphism():
            raise ArithmeticError, "matrix morphism must be an endomorphism"
        A = self.matrix().restrict(sub.free_module())
        H = sage.categories.homset.End(sub, self.domain().category())
        return H(A)

    def trace(self):
        """
        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([V.0+V.1, 2*V.1])
            sage: phi.trace()
            3
        """
        return self.matrix().trace()

class MatrixMorphism(MatrixMorphism_abstract):
    """
    A morphism defined by a matrix.
    """
    def __init__(self, parent, A):
        """
        INPUT:


        -  ``parent`` - a homspace

        -  ``A`` - matrix


        EXAMPLES::

            sage: from sage.modules.matrix_morphism import MatrixMorphism
            sage: T = End(QQ^3)
            sage: M = MatrixSpace(QQ,3)
            sage: I = M.identity_matrix()
            sage: A = MatrixMorphism(T, I)
            sage: loads(A.dumps()) == A
            True
        """
        if not matrix.is_Matrix(A):
            A = matrix.MatrixSpace(parent.category().base_ring(),
                                   parent.domain().rank(),
                                   parent.codomain().rank())(A)
        R = A.base_ring()
        if A.nrows() != parent.domain().rank():
            raise ArithmeticError, "number of rows of matrix (=%s) must equal rank of domain (=%s)"%(A.nrows(), parent.domain().rank())
        if A.ncols() != parent.codomain().rank():
                raise ArithmeticError, "number of columns of matrix (=%s) must equal rank of codomain (=%s)"%(A.ncols(), parent.codomain().rank())
        self._matrix = A
        MatrixMorphism_abstract.__init__(self, parent)

    def matrix(self):
        """
        Return matrix that defines this morphism.

        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([3*V.0, 2*V.1])
            sage: phi.matrix()
            [3 0]
            [0 2]
        """
        return self._matrix

    def _repr_(self):
        """
        Return string representation of this morphism (this is for
        some reason currently not used at all).

        EXAMPLES::

            sage: V = ZZ^2; phi = V.hom([3*V.0, 2*V.1])
            sage: phi._repr_()
            'Morphism defined by the matrix\n[3 0]\n[0 2]'
        """
        if max(self.matrix().nrows(),self.matrix().ncols()) > 5:
            mat = "(not printing %s x %s matrix)"%(self.matrix().nrows(),
                                                   self.matrix().ncols())
        else:
            mat = str(self.matrix())
        return "Morphism defined by the matrix\n%s"%mat

