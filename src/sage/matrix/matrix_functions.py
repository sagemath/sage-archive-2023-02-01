import matrix_generic

import sage.rings.integer
from sage.rings.all import is_Finite, is_FiniteField, is_IntegerModRing, FiniteField


#############################################
## Generic matrices over an integral domain
#############################################

def eigenspaces(self):
    """
    Return a list of pairs
         (e, V)
    where e runs through all eigenvalues (up to Galois conjugation)
    of this matrix, and V is the corresponding eigenspace.

    WARNING: Uses a somewhat naive algorithm (simply factors the
    characteristic polynomial and computes kernels directly over
    the extension field).  TODO: Implement the better algorithm
    that is in dual_eigenvector in sage/hecke/module.py.

    EXAMPLES:
    We compute the eigenspaces of the matrix of the Hecke operator
    $T_2$ on a space:

        sage: A = ModularSymbols(43).T(2).matrix()
        sage: A.eigenspaces()
        [
        (3, [
        (1, 0, 1/7, 0, -1/7, 0, -2/7)
        ]),
        (-2, [
        (0, 1, 0, 1, -1, 1, -1),
        (0, 0, 1, 0, -1, 2, -1)
        ]),
        (a, [
        (0, 1, 0, -1, -a - 1, 1, -1),
        (0, 0, 1, 0, -1, 0, -a + 1)
        ])
        ]

    Next we compute the eigenspaces over the finite field
    of order 11:

        sage: A = ModularSymbols(43, base_ring=GF(11), sign=1).T(2).matrix()
        sage: A.eigenspaces()
        [
        (9, [
        (0, 0, 1, 5)
        ]),
        (3, [
        (1, 6, 0, 6)
        ]),
        (x, [
        (0, 1, 0, 5*x + 10)
        ])
        ]

    We compute the eigenspaces of a $3\times 3$ matrix.

        sage: A = Matrix(QQ,3,3,range(9))
        sage: eA = A.eigenspaces(); eA
        [
        (0, [
        (1, -2, 1)
        ]),
        (a, [
        (1, 1/15*a + 2/5, 2/15*a - 1/5)
        ])
        ]
        sage: eA[1][0].parent()
        Number Field in a with defining polynomial x^2 - 12*x - 18

    This function also allows one to compute eigenvectors and eigenvalues of a
    square matrix over a finite field:

        sage: MS = MatrixSpace(GF(7),2,2)
        sage: g = MS([[5, 1], [4, 1]])
        sage: eig = g.eigenspaces()
        sage: eigvecs = [eig[i][1][0] for i in range(2)]
        sage: eigvals = [eig[i][0] for i in range(2)]
        sage: eigvecs[0]*g == eigvals[0]*eigvecs[0]
        True
    """
    try:
        return self.__eigenvectors
    except AttributeError:
        pass
    if not self.is_square():
        raise ValueError, "matrix must be square"
    f = self.charpoly()
    G = f.factor()
    V = []
    for h, e in G:
        F = h.root_field()
        if F.degree() > 1:
            x = F.gen(0)
        else:
            x = -h[0]
        A = self.change_ring(F) - x
        W = A.kernel()
        if W.dimension() == 0:
            raise RuntimeError, "bug in eigenspaces (dimension of kernel must be nonzero!)"
        V.append((x, W.basis()))
    self.__eigenvectors = V
    return sage.structure.sequence.Sequence(V, cr=True)


def eigenvectors(self):
    """
    EXAMPLES:
        sage: MS = MatrixSpace(GF(7),2,2)
        sage: g = MS([[5, 1], [4, 1]])
        sage: g.eigenvectors()
        [(1, 5), (1, 1)]
    """
    if not self.is_square():
        raise ValueError, "matrix must be square"
    eig = self.eigenspaces()
    n = len(eig)
    return [eig[i][1][0] for i in range(n)]

def eigenvalues(self):
    """
    EXAMPLES:
        sage: MS = MatrixSpace(GF(7),2,2)
        sage: g = MS([[5, 1], [4, 1]])
        sage: g.eigenvalues()
        [4, 2]
    """
    if not self.is_square():
        raise ValueError, "matrix must be square"
    eig = self.eigenspaces()
    n = len(eig)
    return [eig[i][0] for i in range(n)]

def charpoly(self, *args, **kwds):
    r"""
    Return the characteristic polynomial of self, as a polynomial
    over the base ring.

    ALGORITHM: Compute the Hessenberg form of the matrix and read
    off the characteristic polynomial from that.  The result is
    cached.

    If the basering is not a field, use \code{self.matrix_over_field()}
    to obtain the corresponding matrix over a field, then call the
    \code{charpoly} function on it.

    EXAMPLES:
    First a matrix over $\Z$:
        sage: A = MatrixSpace(IntegerRing(),2)( [[1,2], [3,4]] )
        sage: f = A.charpoly()
        sage: f
        x^2 - 5*x - 2
        sage: f.parent()
        Univariate Polynomial Ring in x over Integer Ring

    An example over $\Q$, verifying against the trace and determinant:
        sage: A = MatrixSpace(RationalField(),3)(range(9))
        sage: A.charpoly()
        x^3 - 12*x^2 - 18*x
        sage: A.trace()
        12
        sage: A.determinant()
        0

    We compute the characteristic polynomial of a matrix over
    the polynomial ring $\Z[a]$:
        sage: R = PolynomialRing(IntegerRing(),'a'); a = R.gen()
        sage: M = MatrixSpace(R,2)([[a,1], [a,a+1]])
        sage: M
        [    a     1]
        [    a a + 1]
        sage: f = M.charpoly()
        sage: f
        x^2 + (-2*a - 1)*x + a^2
        sage: f.parent()
        Univariate Polynomial Ring in x over Univariate Polynomial Ring in a over Integer Ring
        sage: M.trace()
        2*a + 1
        sage: M.determinant()
        a^2

    We compute the characteristic polynomial of a matrix over the
    multi-variate polynomial ring $\Z[x,y]$:
        sage: R = MPolynomialRing(IntegerRing(),2); x,y = R.gens()
        sage: A = MatrixSpace(R,2)([x, y, x^2, y^2])
        sage: f = A.charpoly()
        sage: f
        x^2 + (-1*x1^2 - x0)*x + x0*x1^2 - x0^2*x1

    It's a little difficult to distinguish the variables.  To fix this,
    we rename the indeterminate $Z$:
        sage: f.parent().assign_names("Z")
        sage: f
        Z^2 + (-1*x1^2 - x0)*Z + x0*x1^2 - x0^2*x1

    We can pass parameters in, which are passed on to the charpoly
    function for matrices over a field.
        sage: A = 1000*MatrixSpace(ZZ,10)(range(100))
        sage: A.charpoly(bound=2)
        x^10 + 14707*x^9 - 21509*x^8
        sage: A = 1000*MatrixSpace(ZZ,10)(range(100))
        sage: A.charpoly()
        x^10 - 495000*x^9 - 8250000000*x^8


    """
    if not self.is_square():
        raise ValueError, "matrix must be square"
    if self.__charpoly is not None:
        return self.__charpoly
    if not is_Field(self.base_ring()):
        f = self.matrix_over_field().charpoly(*args, **kwds)
        return f.base_extend(self.base_ring())

    if self.nrows() != self.ncols():
        raise ArithmeticError, "charpoly of non-square matrix not defined."

    R = polynomial_ring.PolynomialRing(self.base_ring())
    zero = R(0)
    if self.nrows() == 0:
        self.__charpoly = zero
        return self.__charpoly
    time = misc.verbose(t=0)
    H = self.hessenberg_form()
    n = self.nrows()
    c = [zero for i in range(n+1)]
    c[0] = R(1)
    X = R.gen()
    for m in range(1,n+1):
        c[m] = (X - R(H[m-1,m-1]))*c[m-1]
        t = 1
        for i in range(1,m):
            t = t*H[m-i, m-i-1]
            c[m] = c[m] - R(t*H[m-i-1,m-1])*c[m-i-1]
    misc.verbose('computed characteristic polynomial of %sx%s matrix'%
                 (self.nrows(), self.ncols()), time)
    f = c[n]
    if self.is_immutable():
        self.__charpoly = f
    return f

def fcp(self):
    """
    Return the factorization of the characteristic polynomial of
    self.

    EXAMPLES:
        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,9,-7,4/5,4,3,6,4,3])
        sage: A.fcp()
        (x^3 - 8*x^2 + 209/5*x - 286)
        sage: A = M([3, 0, -2, 0, -2, 0, 0, 0, 0])
        sage: A.fcp()
        (x - 3) * x * (x + 2)
    """
    return self.charpoly().factor()


def determinant(self):
    """
    Return the determinant of this matrix.

    INPUT:
        -- a square matrix

    ALGORITHM: Find the characteristic polynomial and take its
    constant term (up to sign).

    EXAMPLES:
    We create a matrix over $\Z[x,y]$ and compute its determinant.
        sage: R = MPolynomialRing(IntegerRing(),2); x,y = R.gens()
        sage: A = MatrixSpace(R,2)([x, y, x**2, y**2])
        sage: A.determinant()
        x0*x1^2 - x0^2*x1
    """
    if not self.is_square():
        raise ValueError, "matrix must be square"
    # Use stupid slow but completely general method.
    d = (-1)**self.nrows() * self.charpoly()[0]
    return self.base_ring()(d)

def is_invertible(self):
    r"""
    Return True if this matrix is invertible.

    EXAMPLES:
    The following matrix is invertible over $\Q$ but not over $\Z$.
        sage: A = MatrixSpace(IntegerRing(), 2)(range(4))
        sage: A.is_invertible()
        False
        sage: A.matrix_over_field().is_invertible()
        True

    The inverse function is a constructor for matrices over the
    fraction field, so it can work even if A is not invertible.
        sage: ~A   # inverse of A
        [-3/2  1/2]
        [   1    0]

    The next matrix is invertible over $\Z$.
        sage: A = MatrixSpace(IntegerRing(),2)([1,10,0,-1])
        sage: A.is_invertible()
        True
        sage: ~A                # compute the inverse
        [ 1 10]
        [ 0 -1]

    The following nontrivial matrix is invertible over $\Z[x]$.
        sage: R = PolynomialRing(IntegerRing())
        sage: x = R.gen()
        sage: A = MatrixSpace(R,2)([1,x,0,-1])
        sage: A.is_invertible()
        True
        sage: ~A
        [ 1  x]
        [ 0 -1]
    """
    return self.is_square() and self.determinant().is_unit()

def __invert__(self):
    r"""
    Return this inverse of this matrix, as a matrix over the fraction field.

    Raises a \code{ZeroDivisionError} if the matrix has zero
    determinant, and raises an \code{ArithmeticError}, if the
    inverse doesn't exist because the matrix is nonsquare.

    EXAMPLES:
        sage: A = MatrixSpace(IntegerRing(), 2)([1,1,3,5])
        sage: ~A
        [ 5/2 -1/2]
        [-3/2  1/2]

    Even if the inverse lies in the base field, the result is still a matrix
    over the fraction field.
        sage: I = MatrixSpace(IntegerRing(),2)( 1 )  # identity matrix
        sage: ~I
        [1 0]
        [0 1]
        sage: (~I).parent()
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
    """

    if not is_Field(self.base_ring()):
        return ~self.matrix_over_field()
    if not self.is_square():
        raise ArithmeticError, "self must be a square matrix"
    A = self.augment(self.parent().identity_matrix())
    B = A.echelon_form()
    if B[self.nrows()-1,self.ncols()-1] != 1:
        raise ZeroDivisionError, "self is not invertible"
    return B.matrix_from_columns(range(self.ncols(), 2*self.ncols()))


def matrix_over_field(self):
    """
    Return this matrix, but with entries viewed as elements
    of the fraction field of the base ring.

    EXAMPLES:
        sage: A = MatrixSpace(IntegerRing(),2)([1,2,3,4])
        sage: B = A.matrix_over_field()
        sage: B
        [1 2]
        [3 4]
        sage: B.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
    """
    return self.change_ring(self.base_ring().fraction_field())

def numeric_array(self, typecode=None):
    """
    Return the Numeric array associated to this field, if possible, and Numeric
    is installed.
    """
    import Numeric
    if typecode is None:
        typecode = Numeric.Float64
    A = Numeric.array(self.list(), typecode=typecode)
    return Numeric.resize(A,(self.nrows(), self.ncols()))


#############################################
## Generic matrices over a PID
#############################################


def column_module(self):
    """
    Return the free module over the base ring spanned by the
    columns of this matrix.
    """
    return self.transpose().row_space()

def decomposition(self, is_diagonalizable=False, dual=False):
    """
    Returns the decomposition of the free module on which this
    matrix acts from the right, along with whether this matrix
    acts irreducibly on each factor.  The factors are guaranteed
    to be sorted in the same way as the corresponding factors of
    the characteristic polynomial.

    Let A be the matrix acting from the on the vector space V of
    column vectors.  Assume that A is square.  This function
    computes maximal subspaces W_1, ..., W_n corresponding to
    Galois conjugacy classes of eigenvalues of A.  More precisely,
    let f(X) be the characteristic polynomial of A.  This function
    computes the subspace $W_i = ker(g_(A)^n)$, where g_i(X) is an
    irreducible factor of f(X) and g_i(X) exactly divides f(X).
    If the optional parameter is_diagonalizable is True, then we
    let W_i = ker(g(A)), since then we know that ker(g(A)) =
    $ker(g(A)^n)$.

    If dual is True, also returns the corresponding decomposition
    of V under the action of the transpose of A.  The factors are
    guarenteed to correspond.

    OUTPUT:
        list -- list of pairs (V,t), where V is a vector spaces
                and t is a bool, and t is True exactly when the
                charpoly of self on V is irreducible.

        (optional) list -- list of pairs (W,t), where W is a vector
                space and t is a bool, and t is True exactly
                when the charpoly of the transpose of self on W
                is irreducible.

    EXAMPLES:
        sage: MS1 = MatrixSpace(ZZ,4)
        sage: MS2 = MatrixSpace(QQ,6)
        sage: A = MS1.matrix([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
        sage: B = MS2(range(36))
        sage: B*11   # random output
        [-11  22 -11 -11 -11 -11]
        [ 11 -22 -11 -22  11  11]
        [-11 -11 -11 -22 -22 -11]
        [-22  22 -22  22 -11  11]
        [ 22 -11  11 -22  11  22]
        [ 11  11  11 -22  22  22]
        sage: decomposition(A)
        [(Ambient free module of rank 4 over the principal ideal domain Integer Ring, True)]
        sage: decomposition(B)
        [(Vector space of degree 6 and dimension 4 over Rational Field
        Basis matrix:
        [ 1  0  0  0 -5  4]
        [ 0  1  0  0 -4  3]
        [ 0  0  1  0 -3  2]
        [ 0  0  0  1 -2  1],
          False),
         (Vector space of degree 6 and dimension 2 over Rational Field
        Basis matrix:
        [ 1  0 -1 -2 -3 -4]
        [ 0  1  2  3  4  5],
          True)]
    """
    if not self.is_square():
        raise ArithmeticError, "self must be a square matrix"

    if self.nrows() == 0:
        return []

    f = self.charpoly()
    E = []

    # Idea: For optimization, could compute powers of self
    #       up to max degree of any factor.  Then get g(self)
    #       by taking a linear combination.   ??????

    if dual:
        Edual = []
    F = f.factor()
    if len(F) == 1:
        V = sage.modules.free_module.FreeModule(
                          self.base_ring(), self.nrows())
        m = F[0][1]
        if dual:
            return [(V,m==1)], [(V,m==1)]
        else:
            return [(V,m==1)]
    F.sort()
    for g, m in f.factor():
        if is_diagonalizable:
            B = g(self)
        else:
            B = g(self) ** m
        E.append((B.kernel(), m==1))
        if dual:
            Edual.append((B.transpose().kernel(), m==1))
    if dual:
        return E, Edual
    return E


def image(self):
    """
    Return the image of the homomorphism on rows defined by this matrix.

    EXAMPLES:
        sage: MS1 = MatrixSpace(ZZ,4)
        sage: MS2 = MatrixSpace(QQ,6)
        sage: A = MS1.matrix([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
        sage: B = MS2.random_element()

        sage: image(A)
        Free module of degree 4 and rank 4 over Integer Ring
        Echelon basis matrix:
        [  1   0   0 426]
        [  0   1   0 518]
        [  0   0   1 293]
        [  0   0   0 687]

        sage: image(B) == B.row_module()
        True
    """
    return self.row_module()

def row_module(self):
    """
    Return the free module over the base ring spanned by the rows
    of self.

    EXAMPLES:
        sage: A = MatrixSpace(IntegerRing(), 2)([1,2,3,4])
        sage: A.row_module()
        Free module of degree 2 and rank 2 over Integer Ring
        Echelon basis matrix:
        [1 0]
        [0 2]
    """
    M = sage.modules.free_module.FreeModule(self.base_ring(), self.ncols())
    return M.span(self.rows())

def kernel_on(self, V, poly=None, check=False):
    """
    Return the kernel of self restricted to the invariant subspace V.
    The result is a vector subspace of V, which is also a subspace
    of the ambient space.

    INPUT:
        V -- vector subspace
        check -- (optional) default: False
        poly -- (optional) default: None; if not None, compute instead
                the kernel of poly(self) on V.

    OUTPUT:
        a subspace

    WARNING: This function does \emph{not} check that V is in fact
    invariant under self, unless check is True (not the default).
    """
    A = self.restrict(V, check=check)
    if not poly is None:
        A = poly(A)
    W = A.kernel()
    if V.is_ambient():
        return W
    else:
        A = V.basis_matrix()
        B = W.basis_matrix()
        C = B*A
        return C.row_module()

#############################################
## Generic matrices over integers (may also make sense over any PID)
#############################################

def rank(self):
    """
    Return the rank of self, which is the rank of the space
    spanned by the rows of self.
    """
    if self.__rank is not None:
        return self.__rank
    rank = self.echelon_form().rank()
    if self.is_immutable():
        self.__rank = rank
    return rank


def elementary_divisors(self):
    """
    Return the elementary divisors of self, in order.

    The elementary divisors are the invariants of the finite
    abelian group that is the cokernel of this matrix.  They are
    ordered in reverse by divisibility.

    INPUT:
        matrix
    OUTPUT:
        list of int's

    EXAMPLES:
        sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
        sage: A.elementary_divisors()
        [0, 3, 1]
        sage: C = MatrixSpace(ZZ,4)([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
        sage: C.elementary_divisors()
        [687, 1, 1, 1]

    SEE ALSO: smith_form
    """
    try:
        return self.__elementary_divisors
    except AttributeError:
        if self.nrows() == 0 or self.ncols() == 0:
            return []
        d = self._pari_().matsnf(0).python()
        if self.is_immutable():
            self.__elementary_divisors = d
        return d



def smith_form(self, transformation=False):
    """
    Returns matrices S, U, and V such that S = U*self*V, and S
    is in Smith normal form.  Thus S is diagonal with diagonal
    entries the ordered elementary divisors of S.

    The elementary divisors are the invariants of the finite
    abelian group that is the cokernel of this matrix.  They are
    ordered in reverse by divisibility.

    EXAMPLES:
        sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
        sage: D, U, V = A.smith_form()
        sage: D
        [0 0 0]
        [0 3 0]
        [0 0 1]
        sage: U
        [-1  2 -1]
        [ 0 -1  1]
        [ 0  1  0]
        sage: V
        [ 1  4 -1]
        [-2 -3  1]
        [ 1  0  0]
        sage: U*A*V
        [0 0 0]
        [0 3 0]
        [0 0 1]

    It also makes sense for nonsquare matrices:

        sage: A = Matrix(ZZ,3,2,range(6))
        sage: D, U, V = A.smith_form()
        sage: D
        [0 0]
        [2 0]
        [0 1]
        sage: U
        [-1  2 -1]
        [ 0 -1  1]
        [ 0  1  0]
        sage: V
        [ 3 -1]
        [-2  1]
        sage: U * A * V
        [0 0]
        [2 0]
        [0 1]

    SEE ALSO: elementary_divisors
    """
    v = self._pari_().matsnf(1).python()
    D = self.matrix_space()(v[2])
    U = self.matrix_space(ncols = self.nrows())(v[0])
    V = self.matrix_space(nrows = self.ncols())(v[1])
    return D, U, V

def frobenius(self,flag=0):
    """
    Return the Frobenius form (rational canonical form) of this matrix.

    If flag is 1, return only the elementary divisors.  If flag is
    2, return a two-components vector [F,B] where F is the
    Frobenius form and B is the basis change so that $M=B^{-1}FB$.

    INPUT:
       flag -- 0 (default), 1 or 2 as described above

    ALGORITHM: uses pari's matfrobenius()

    EXAMPLE:
       sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
       sage: A.frobenius(0)
       [ 0  0  0]
       [ 1  0 18]
       [ 0  1 12]
       sage: A.frobenius(1)
       [x^3 - 12*x^2 - 18*x]
       sage: A.frobenius(2)
       ([ 0  0  0]
       [ 1  0 18]
       [ 0  1 12],
       [    -1      2     -1]
       [     0  23/15 -14/15]
       [     0  -2/15   1/15])

    AUTHOR:
       -- 2006-04-02: Martin Albrecht

    TODO:
       -- move this to work for more general matrices than just over Z.
          This will require fixing how PARI polynomials are coerced
          to SAGE polynomials.
    """
    Z = integer_ring.IntegerRing()
    if self.base_ring() is not Z:
        raise NotImplementedError, "frobenius only implemented for matrices over Z"
    if self.nrows()!=self.ncols():
        raise ArithmeticError, "frobenius matrix of non-square matrix not defined."

    v = self._pari_().matfrobenius(flag)
    if flag==0:
        return self.matrix_space()(v.python())
    elif flag==1:
        r = polynomial_ring.PolynomialRing(self.base_ring())
        #TODO: this should be handled in PolynomialRing not here
        return [eval(str(x).replace("^","**"), {}, r.gens_dict()) for x in v]
    elif flag==2:
        F = matrix_space.MatrixSpace(rational_field.RationalField(),
                                     self.nrows())(v[0].python())
        B = matrix_space.MatrixSpace(rational_field.RationalField(),
                                     self.nrows())(v[1].python())
        return F,B


def _adjoint(self):
    """assumes self is a square matrix (checked in adjoint)"""
    return self.parent()(self._pari_().matadjoint().python())

def _lllgram(self):
    """assumes self is a square matrix (checked in lllgram)"""
    Z = integer_ring.IntegerRing()
    if self.base_ring() is not Z:
        raise NotImplementedError, "_lllgram only implemented for matrices over Z"
    n = self.nrows()
    # pari does not like negative definite forms
    if n > 0 and self[0,0] < 0:
        self = -self
    # maybe should be /unimodular/ matrices ?
    MS = matrix_space.MatrixSpace(Z,n,n)
    try:
        U = MS(self._pari_().lllgramint().python())
    except (RuntimeError, ArithmeticError):
        raise ValueError, "not a definite matrix"
    # Fix last column so that det = +1
    if U.det() == -1:
        for i in range(n):
            U[i,n-1] = - U[i,n-1]
    return U

def _ntl_(self):
    """
    ntl.mat_ZZ representation of self.

    \note{NTL only knows dense matrices, so if you provide a
    sparse matrix NTL will allocate memory for every zero entry.}
    """
    return mat_ZZ(self.nrows(),self.ncols(),self.list())



#############################################
## Generic matrices over a field
#############################################


def column_space(self):
    """
    Return the vector space over the base ring spanned by the
    columns of this matrix.

    EXAMPLES:
        sage: M = MatrixSpace(QQ,3,3)
        sage: A = M([1,9,-7,4/5,4,3,6,4,3])
        sage: A.column_space()
        Vector space of degree 3 and dimension 3 over Rational Field
        Basis matrix:
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: W = MatrixSpace(CC,2,2)
        sage: B = W([1, 2+3*I,4+5*I,9]); B
        [                       1.0000000000000000 2.0000000000000000 + 3.0000000000000000*I]
        [4.0000000000000000 + 5.0000000000000000*I                        9.0000000000000000]
        sage: B.column_space()
        Vector space of degree 2 and dimension 2 over Complex Field with 53 bits of precision
        Basis matrix:
        [                                                       1.0000000000000000 0.00000000000000044408920985006262 + 0.00000000000000088817841970012523*I]
        [                                                                        0               0.99999999999999978 - 0.000000000000000055511151231257827*I]
    """
    return self.column_module()


def decomposition_of_subspace(self, M, is_diagonalizable=False):
    """
    Suppose the right action of self on M leaves M
    invariant. Return the decomposition of M as a list of pairs
    (W, is_irred) where is_irred is True if the charpoly of self
    acting on the factor W is irreducible.
    """
    if not sage.modules.free_module.is_FreeModule(M):
        raise TypeError, "M must be a free module."
    if not self.is_square():
        raise ArithmeticError, "matrix must be square"
    if M.base_ring() != self.base_ring():
        raise ArithmeticError, "base rings are incompatible"
    if M.degree() != self.ncols():
        raise ArithmeticError, \
           "M must be a subspace of an %s-dimensional space"%self.ncols()

    time = misc.verbose(t=0)

    # 1. Restrict
    B = self.restrict(M)
    time0 = misc.verbose("restrict -- ", time)

    # 2. Decompose restriction
    D = B.decomposition(is_diagonalizable=is_diagonalizable, dual=False)

    assert sum([A.dimension() for A,_ in D]) == M.dimension(), "bug in decomposition; " + \
           "the sum of the dimensions of the factors must equal the dimension of the space acted on"

    # 3. Lift decomposition to subspaces of ambient vector space.
    # Each basis vector for an element of D defines a linear combination
    # of the basis of W, and these linear combinations define the
    # corresponding subspaces of the ambient space M.

    misc.verbose("decomposition -- ", time0)
    C = M.basis_matrix()
    Z = M.ambient_vector_space()

    D = [(Z.subspace([x*C for x in W.basis()]), is_irred) for W, is_irred in D]

    misc.verbose(t=time)
    return D

def denominator(self):
    r"""
    Return the least common multiple of the denominators of the
    elements of self.

    If there is no denominator function for the base field, or no
    LCM function for the denominators, raise a TypeError.

    EXAMPLES:
        sage: A = MatrixSpace(RationalField(),2)(['1/2', '1/3', '1/5', '1/7'])
        sage: A.denominator()
        210

    Denominators are note defined for real numbers:
        sage: A = MatrixSpace(RealField(),2)([1,2,3,4])
        sage: A.denominator()
        Traceback (most recent call last):
        ...
        TypeError: denominator not defined for elements of the base ring

    We can even compute the denominator of matrix over the fraction field
    of $\Z[x]$.
        sage: K.<x> = FractionField(PolynomialRing(IntegerRing()))
        sage: A = MatrixSpace(K,2)([1/x, 2/(x+1), 1, 5/(x^3)])
        sage: A.denominator()
        x^4 + x^3

    Here's an example involving a cyclotomic field:
        sage: K.<z> = CyclotomicField(3)
        sage: M = MatrixSpace(K,3,sparse=True)
        sage: A = M([(1+z)/3,(2+z)/3,z/3,1,1+z,-2,1,5,-1+z])
        sage: print A
        [1/3*z + 1/3 1/3*z + 2/3       1/3*z]
        [          1       z + 1          -2]
        [          1           5       z - 1]
        sage: print A.denominator()
        3
    """
    if self.nrows() == 0 or self.ncols() == 0:
        return integer.Integer(1)
    R = self.base_ring()
    x = self.list()
    try:
        d = x[0].denominator()
    except AttributeError:
        raise TypeError, "denominator not defined for elements of the base ring"
    try:
        for y in x:
            d = d.lcm(y.denominator())
    except AttributeError:
        raise TypeError, "lcm function not defined for elements of the base ring"
    return d


def echelon_form(self, include_zero_rows=True):
    """
    Returns the reduced row echelon form of self.

    INPUT:
        matrix -- an element A of a MatrixSpace

    OUTPUT:
        matrix -- The reduced row echelon form of A.
        Note that self is *not* changed by this command.

    EXAMPLES:
       sage: MS = MatrixSpace(RationalField(),2,3)
       sage: C = MS.matrix([1,2,3,4,5,6])
       sage: C.rank()
       2
       sage: C.nullity()
       1
       sage: C.echelon_form()
       [ 1  0 -1]
       [ 0  1  2]
    """
    if self.__echelon_form is not None:
        return self.__echelon_form

    if not is_Field(self.base_ring()):
        raise NotImplementedError, ""

    R = self.base_ring()
    # Fix to work with finite fields and Z/nZ, which was
    # suggested by Dan Christensen <jdc@uwo.ca>.
    if (   (is_FiniteField(R) and R.is_prime_field()) or \
           is_IntegerModRing(R)  ) and R.characteristic() < 46340:
        p = R.characteristic()
        S = sage.matrix.dense_matrix_pyx.Matrix_modint(p, self.nrows(), self.ncols(), self.list())
        S.echelon()
        A = self.parent()(S.list())    # most of time is spent here!?
        pivot_positions = S.pivots()

    else:

        t = misc.verbose("Generic echelon...")
        pivot_positions = []
        start_row = 0
        A = self.copy()
        nrows = A.nrows()
        ncols = A.ncols()
        cleared_a_column = False
        for c in range(ncols):
            misc.verbose("column %s of %s"%(c, ncols),t, level=2)
            for r in range(start_row, nrows):
                if A.get((r,c)) != 0:
                    pivot_positions.append(c)
                    # Divide row r through by 1/A[r,c], so leading coefficient
                    # is 1.
                    z = A.get((r,c))
                    if z != 1:
                        A.rescale_row(r, ~z)
                    # Swap
                    if r != start_row:
                        A.swap_rows(r,start_row)
                    # Clear column
                    cleared_a_column = True
                    for i in range(nrows):
                        if i != start_row:
                            x = A.get((i,c))
                            if x != 0:
                                # Add -x times start row to i
                                A.add_multiple_of_row(i, start_row, -x)
                if cleared_a_column:
                    start_row += 1
                    cleared_a_column = False
                    break
        # end for
        misc.verbose("Finished generic echelon.",t)
    #end if

    if not include_zero_rows:
        A = A.matrix_from_rows(range(len(pivot_positions)))
    A.__pivots = pivot_positions
    A.__rank = len(pivot_positions)
    A.set_immutable()
    return A



def hessenberg_form(self):
    """
    Return the Hessenberg form of self.

    The hessenberg form of a matrix $A$ is a matrix that is
    similar to $A$, so has the same characteristic polynomial as
    $A$, and is upper triangular except possible for entries right
    below the diagonal.

    ALGORITHM: See Henri Cohen's first book.

    EXAMPLES:
        sage: A = MatrixSpace(RationalField(),3)([2, 1, 1, -2, 2, 2, -1, -1, -1])
        sage: A.hessenberg_form()
        [  2 3/2   1]
        [ -2   3   2]
        [  0  -3  -2]

        sage: A = MatrixSpace(RationalField(),4)([2, 1, 1, -2, 2, 2, -1, -1, -1,1,2,3,4,5,6,7])
        sage: A.hessenberg_form()
        [    2  -7/2 -19/5    -2]
        [    2   1/2 -17/5    -1]
        [    0  25/4  15/2   5/2]
        [    0     0  58/5     3]
    """
    if not self.is_square():
        raise ArithmeticError, "self must be square"
    n = self.nrows()
    tm = misc.verbose("Computing Hessenberg Normal Form of %sx%s matrix"%(n,n))
    h = self.copy()
    for m in range(1,n-1):
        # Search for a non-zero entry in column m-1
        i = False
        for r in range(m+1,n):
            if h[r,m-1] != 0:
                i = r
                break
        if i:
            # Found a nonzero entry in column m-1 that is strictly below row m
            # Now set i to be the first nonzero position >= m in column m-1
            if h[m,m-1] != 0:
                i = m
            t = h[i,m-1]
            if i>m:
                h.swap_rows(i,m)
                # We must do the corresponding column swap to
                # maintain the characteristic polynomial (which is
                # an invariant of Hessenberg form)
                h.swap_columns(i,m)
            # Now the nonzero entry in position (m,m-1) is t.
            # Use t to clear the entries in column m-1 below m.
            for j in range(m+1,n):
                if h[j,m-1] != 0:
                    u = h[j,m-1]/t
                    h.add_multiple_of_row(j, m, -u)
                    # To maintain charpoly, do the corresponding column operation,
                    # which doesn't mess up the matrix, since it only changes
                    # column m, and we're only worried about column m-1 right now.
                    # Add u*column_j to column_m.
                    h.add_multiple_of_column(m, j, u)
    misc.verbose("Finished Hessenberg Normal Form of %sx%s matrix"%(n,n),tm)
    return h



def kernel(self, *args, **kwds):
    r"""
    Return the kernel of this matrix, as a vector space.

    INPUT:
        -- all additional arguments to the kernel function
           are passed directly onto the echelon call.

    \algorithm{Elementary row operations don't change the kernel,
    since they are just right multiplication by an invertible
    matrix, so we instead compute kernel of the column echelon
    form.  More precisely, there is a basis vector of the kernel
    that corresponds to each non-pivot row.  That vector has a 1
    at the non-pivot row, 0's at all other non-pivot rows, and for
    each pivot row, the negative of the entry at the non-pivot row
    in the column with that pivot element.}

    \note{Since we view matrices as acting on the right, but have
    functions for reduced \emph{row} echelon forms, we instead
    compute the reduced row echelon form of the transpose of this
    matrix, which is the reduced column echelon form.}

    EXAMPLES:

    A kernel of dimension one over $\Q$:x
        sage: A = MatrixSpace(QQ, 3)(range(9))
        sage: A.kernel()
        Vector space of degree 3 and dimension 1 over Rational Field
        Basis matrix:
        [ 1 -2  1]

    A trivial kernel:
        sage: A = MatrixSpace(QQ, 2)([1,2,3,4])
        sage: A.kernel()
        Vector space of degree 2 and dimension 0 over Rational Field
        Basis matrix:
        []

    Kernel of a zero matrix:
        sage: A = MatrixSpace(QQ, 2)(0)
        sage: A.kernel()
        Vector space of degree 2 and dimension 2 over Rational Field
        Basis matrix:
        [1 0]
        [0 1]

    Kernel of a non-square matrix:
        sage: A = MatrixSpace(QQ,3,2)(range(6))
        sage: A.kernel()
        Vector space of degree 3 and dimension 1 over Rational Field
        Basis matrix:
        [ 1 -2  1]

    The 2-dimensional kernel of a matrix over a cyclotomic field:
        sage: K = CyclotomicField(12); a=K.0
        sage: M = MatrixSpace(K,4,2)([1,-1, 0,-2, 0,-a**2-1, 0,a**2-1])
        sage: M
        [             1             -1]
        [             0             -2]
        [             0 -zeta12^2 - 1]
        [             0  zeta12^2 - 1]
        sage: M.kernel()
        Vector space of degree 4 and dimension 2 over Cyclotomic Field of order 12 and degree 4
        Basis matrix:
        [               0                1                0     -2*zeta12^2]
        [               0                0                1 -2*zeta12^2 + 1]

    A nontrivial kernel over a complicated base field.
        sage: K = FractionField(MPolynomialRing(QQ, 2))
        sage: M = MatrixSpace(K, 2)([[K.1, K.0], [K.1, K.0]])
        sage: M
        [x1 x0]
        [x1 x0]
        sage: M.kernel()
        Vector space of degree 2 and dimension 1 over Fraction Field of Polynomial Ring in x0, x1 over Rational Field
        Basis matrix:
        [ 1 -1]
    """

    R = self.base_ring()

    if self.nrows() == 0:    # from a 0 space
        V = sage.modules.free_module.VectorSpace(R, self.nrows())
        return V.zero_subspace()

    elif self.ncols() == 0:  # to a 0 space
        return sage.modules.free_module.VectorSpace(R, self.nrows())

    if isinstance(R, number_field.NumberField_generic):
        A = self._pari_().mattranspose()
        B = A.matker()
        n = self.nrows()
        V = sage.modules.free_module.VectorSpace(R, n)
        basis = [V([R(x) for x in b]) for b in B]
        return V.subspace(basis)

    E = self.transpose().echelon_form(*args, **kwds)
    pivots = E.pivots()
    pivots_set = set(pivots)
    basis = []
    VS = sage.modules.free_module.VectorSpace
    V = VS(R, self.nrows())
    ONE = R(1)
    for i in xrange(self.nrows()):
        if not (i in pivots_set):
            v = V(0)
            v[i] = ONE
            for r in range(len(pivots)):
                v[pivots[r]] = -E[r,i]
            basis.append(v)
    return V.subspace(basis)



def maxspin(self, v):
    """
    Computes the largest integer n such that the list of vectors
    $S=[v, A(v), ..., A^n(v)]$ are linearly independent, and returns
    that list.

    INPUT:
        self -- Matrix
        v -- Vector
    OUTPUT:
        list -- list of Vectors

    ALGORITHM:
        The current implementation just adds vectors to a vector
        space until the dimension doesn't grow.  This could be
        optimized by directly using matrices and doing an
        efficient Echelon form.  Also, when the base is Q, maybe
        we could simultaneously keep track of what is going on in
        the reduction modulo p, which might make things much
        faster.
    """
    if v == 0: return []
    VS = sage.modules.free_module.VectorSpace
    V = VS([v])
    w = v
    S = [v]
    while True:
        w = w*self
        W = V + VS([w])
        if W.dimension() == V.dimension():
            return S
        V = W
        S.append(w)

def nullity(self):
    # Use that rank + nullity = number of columns
    return self.ncols() - self.rank()

def pivots(self):
    """
    Return the i such that the i-th column of self is a
    pivot column of the reduced row echelon form of self.

    OUTPUT:
        list -- sorted list of integers
    """
    try:
        return self.__pivots
    except AttributeError:
        P = self.echelon_form().pivots()
        if self.is_immutable():
            self.__pivots = P
        return P



def restrict(self, V, check=True):
    """
    Returns the matrix that defines the action of self on the
    chosen basis for the invariant subspace V.  If V is an
    ambient, returns self (not a copy of self).

    INPUT:
        V -- vector subspace
        check -- (optional) default: True; if False may not check
                 that V is invariant (hence can be faster).
    OUTPUT:
        a matrix

    WARNING:
    This function returns an nxn matrix, where V has dimension n.
    It does \emph{not} check that V is in fact invariant under
    self, unless check is True (not the default).

    EXAMPLES:
        sage: V = VectorSpace(QQ, 3)
        sage: M = MatrixSpace(QQ, 3)
        sage: A = M([1,2,0, 3,4,0, 0,0,0])
        sage: W = V.subspace([[1,0,0], [0,1,0]])
        sage: A.restrict(W)
        [1 2]
        [3 4]
        sage: A.restrict(W, check=True)
        [1 2]
        [3 4]

    We illustrate the warning about invariance not being checked
    by default, by giving a non-invariant subspace.  With the default
    check=False this function returns the 'restriction' matrix, which
    is meaningless as check=True reveals.
        sage: W2 = V.subspace([[1,0,0], [0,1,1]])
        sage: A.restrict(W2, check=False)
        [1 2]
        [3 4]
        sage: A.restrict(W2, check=True)
        Traceback (most recent call last):
        ...
        ArithmeticError: subspace is not invariant under matrix
    """
    if not isinstance(V, sage.modules.free_module.FreeModule_generic):
        raise TypeError, "V must be a Vector Space"
    if V.base_field() != self.base_ring():
        raise TypeError, "base rings must be the same"
    if V.degree() != self.nrows():
        raise IndexError, "degree of V (=%s) must equal number of rows of self (=%s)"%(V.degree(), self.nrows())
    if V.rank() == 0:
        return self.new_matrix(nrows=0, ncols=0)

    if not check and V.base_ring().is_field() and not V.has_user_basis():
        B = V.echelonized_basis_matrix()
        P = B.pivots()
        return B*self.matrix_from_columns(P)
    else:
        n = V.rank()
        try:
            # todo optimize so only involves matrix multiplies ?
            C = [V.coordinate_vector(b*self) for b in V.basis()]
        except ArithmeticError:
            raise ArithmeticError, "subspace is not invariant under matrix"
        return self.new_matrix(n, n, C, sparse=False)

def restrict_domain(self, V):
    """
    Compute the matrix relative to the basis for V on the domain
    obtained by restricting self to V, but not changing the
    codomain of the matrix.  This is the matrix whose rows are the
    images of the basis for V.

    INPUT:
        V -- vector space (subspace of ambient space on which self acts)

    SEE ALSO: restrict()

    EXAMPLES:
        sage: V = VectorSpace(QQ, 3)
        sage: M = MatrixSpace(QQ, 3)
        sage: A = M([1,2,0, 3,4,0, 0,0,0])
        sage: W = V.subspace([[1,0,0], [1,2,3]])
        sage: A.restrict_domain(W)
        [1 2 0]
        [3 4 0]
        sage: W2 = V.subspace_with_basis([[1,0,0], [1,2,3]])
        sage: A.restrict_domain(W2)
        [ 1  2  0]
        [ 7 10  0]
    """
    return self.new_matrix(V.dimension(), self.ncols(), [b*self for b in V.basis()])



def wiedemann(self, i, t=0):
    """
    Application of Wiedemann's algorithm to the i-th standard
    basis vector.

    If the optimal parameter t is nonzero, use only the first t
    linear recurrence relations.
    """
    i = int(i); t=int(t)
    if self.nrows() != self.ncols():
        raise ArithmeticError, "matrix must be square."
    n = self.nrows()
    v = sage.modules.free_module.VectorSpace(self.base_ring(), n).gen(i)
    tm = misc.verbose('computing iterates...')
    cols = self.iterates(v, 2*n).columns()
    tm = misc.verbose('computed iterates', tm)
    f = None
    # Compute the minimal polynomial of the linear recurrence
    # sequence corresponding to the 0-th entries of the iterates,
    # then the 1-th entries, etc.
    if t == 0:
        R = range(n)
    else:
        R = [t]
    for i in R:
        tm = misc.verbose('applying berlekamp-massey')
        g = berlekamp_massey.berlekamp_massey(cols[i].list())
        misc.verbose('berlekamp-massey done', tm)
        if f is None:
            f = g
        else:
            f = f.lcm(g)
        if f.degree() == n:
            break
    return f
