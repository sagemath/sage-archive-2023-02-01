r"""
Matrices over specific rings

This section describes classes that implement matrices over special
rings, such as the integers or rationals.  They derive from the
generic \class{Matrix} class, but implement improved and additional
functionality.
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import copy
import operator

import sage.rings.arith
import sage.misc.misc as misc
import sage.misc.latex as latex
import sage.matrix.dense_matrix_pyx as dense_matrix_pyx
import sage.matrix.sparse_matrix_pyx as sparse_matrix_pyx

import sage.matrix.sparse_matrix
import sage.modules.free_module_element
import sage.modules.free_module
import matrix_space
import sage.libs.pari.all as pari

import berlekamp_massey
import sage.rings.polynomial_ring as polynomial_ring
import sage.rings.integer_ring as integer_ring
import sage.rings.integer as integer
import sage.rings.rational_field as rational_field
import sage.rings.rational as rational
import sage.rings.number_field.number_field as number_field
import sage.rings.coerce as coerce
from sage.rings.all import is_FiniteField, is_IntegerModRing, FiniteField

from sage.structure.mutability import Mutability

import sage.ext.dense_matrix_pyx

from sage.interfaces.all import singular as singular_default

from sage.libs.ntl.all import mat_ZZ
cputime = misc.cputime

from sage.matrix.matrix_pyx import Matrix, __reduce__Matrix_pyx

def is_Matrix(x):
    return isinstance(x, Matrix)

#############################################
## Generic matrices over an integral domain
#############################################

class Matrix_domain(Matrix):
    def __init__(self, parent):
        Matrix.__init__(self, parent)

    def charpoly(self, *args, **kwds):
        r"""
        Return the characteristic polynomial of self, as a polynomial
        over the base ring.

        ALGORITHM: Use \code{self.matrix_over_field()} to obtain the
        corresponding matrix over a field, then call the
        \code{charpoly} function on it.

        EXAMPLES:
        First a matrix over $\Z$:
            sage: A = MatrixSpace(IntegerRing(),2)( [[1,2], [3,4]] )
            sage: f = A.charpoly()
            sage: f
            x^2 - 5*x - 2
            sage: f.parent()
            Univariate Polynomial Ring in x over Integer Ring

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
        f = self.matrix_over_field().charpoly(*args, **kwds)
        return f.base_extend(self.base_ring())

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
        if self.nrows() != self.ncols():
            raise ArithmeticError, "Matrix must be square, but is %sx%s"%(
                self.nrows(), self.ncols())
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
        return ~self.matrix_over_field()


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
class Matrix_pid(Matrix_domain):
    def __init__(self, parent):
        Matrix_domain.__init__(self, parent)

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

    def echelon_form(self, include_zero_rows=True):
        """
        Return the echelon form of this matrix over the integers.

        This is a matrix over the base ring (a PID) which is, \emph{by
        definition}, what is also called the Hermite normal form.
        """
        raise NotImplementedError


    def image(self):
        """
        Return the image of the homomorphism on rows defined by this matrix.
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
## Generic matrices over the integers
#############################################

def _parimatrix_to_strlist(A):
    s = str(A)
    s = s.replace('Mat(','').replace(')','')
    s = s.replace(';',',').replace(' ','')
    s = s.replace(",", "','")
    s = s.replace("[", "['")
    s = s.replace("]", "']")
    return eval(s)

def _parimatrix_to_reversed_strlist(A):
    s = str(A)
    if s.find('Mat') != -1:
        return _parimatrix_to_strlist(A)
    s = s.replace('[','').replace(']','').replace(' ','')
    v = s.split(';')
    v.reverse()
    s = "['" + (','.join(v)) + "']"
    s = s.replace(",", "','")
    return eval(s)


class Matrix_integer(Matrix_pid):
    def __init__(self, parent):
        Matrix_pid.__init__(self, parent)

    def echelon_form(self, include_zero_rows=True):
        r"""
        Return the echelon form of this matrix over the integers.

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(),2)([1,2,3,4])
            sage: A.echelon_form()
            [1 0]
            [0 2]

            sage: A = MatrixSpace(IntegerRing(),5)(range(25))
            sage: A.echelon_form()
            [  5   0  -5 -10 -15]
            [  0   1   2   3   4]
            [  0   0   0   0   0]
            [  0   0   0   0   0]
            [  0   0   0   0   0]
        """
        if self.nrows() == 0 or self.ncols() == 0:
            self.__rank = 0
            return self
        # The following complicated sequence of column reversals
        # and transposes is needed since PARI's Hermite Normal Form
        # does column operations instead of row operations.
        n = self.ncols()
        r = [n-i for i in range(n)]
        v = self._pari_()
        v = v.vecextract(r) # this reverses the order of columns
        v = v.mattranspose()
        w = v.mathnf(1)
        def convert_parimatrix(z):
            n = z.ncols(); r = [n-i for i in range(n)]
            z = z.vecextract(r)
            z = z.mattranspose()
            n = z.ncols(); r = [n-i for i in range(n)]
            z = z.vecextract(r)
            return _parimatrix_to_strlist(z)

        H = convert_parimatrix(w[0])
        if self.ncols() == 1:
            H = [H]

        # We can do a 'fast' change of the above into a list of ints,
        # since we know all entries are ints:
        (nr,nc) = (self.nrows(), self.ncols())
        num_missing_rows = (nr*nc - len(H)) / nc
        rank = nr - num_missing_rows
        if include_zero_rows:
            H += ['0']*(num_missing_rows*nc)
            H = self.new_matrix(nrows=nr, ncols=nc, entries=H, coerce_entries=True)
        else:
            H = self.new_matrix(nrows=rank, ncols=nc, entries=H, coerce_entries=True)
        H.__rank = rank
        H.set_immutable()
        return H

    def rank(self):
        """
        Return the rank of self, which is the rank of the space
        spanned by the rows of self.
        """
        try:
            return self.__rank
        except AttributeError:
            r = self.change_ring(rational_field.RationalField()).rank()
            if self.is_immutable():
                self.__rank = r
            return r

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

    def kernel(self, LLL=False):
        r"""
        Return the kernel of this matrix, as a module over the integers.

        INPUT:
           LLL -- bool (default: False); if True the basis is an LLL
                  reduced basis; otherwise, it is an echelon basis.

        EXAMPLES:
            sage: M = MatrixSpace(IntegerRing(),4,2)(range(8))
            sage: M.kernel()
            Free module of degree 4 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -3  2]
            [ 0  1 -2  1]
        """

        Z = self.base_ring()

        if self.nrows() == 0:    # from a 0 space
            M = sage.modules.free_module.FreeModule(Z, self.nrows())
            return M.zero_subspace()

        elif self.ncols() == 0:  # to a 0 space
            return sage.modules.free_module.FreeModule(Z, self.nrows())

        A = self._pari_().mattranspose()
        B = A.matkerint()


        n = self.nrows()
        Z = integer_ring.IntegerRing()
        M = sage.modules.free_module.FreeModule(Z, n)

        if B.ncols() == 0:
            return M.zero_submodule()

        # Now B is a basis for the LLL-reduced integer kernel as a
        # PARI object.  The basis vectors or B[0], ..., B[n-1],
        # where n is the dimension of the kernel.
        B = [M([Z(x) for x in b]) for b in B]
        if LLL:
            return M.span_of_basis(B)
        else:
            return M.span(B)

    def _adjoint(self):
        """assumes self is a square matrix (checked in adjoint)"""
        return self.parent()(self._pari_().matadjoint().python())

    def _lllgram(self):
        """assumes self is a square matrix (checked in lllgram)"""
        Z = integer_ring.IntegerRing()
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
class Matrix_field(Matrix_pid):
    def __init__(self, parent):
        Matrix_pid.__init__(self, parent)

    def __invert__(self):
        """
        Return this inverse of this matrix.

        Raises a ZeroDivisionError if the matrix has zero determinant, and
        raises an ArithmeticError, if the inverse doesn't exist
        because the matrix is nonsquare.

        EXAMPLES:

        """
        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"
        A = self.augment(self.parent().identity_matrix())
        B = A.echelon_form()
        if B[self.nrows()-1,self.ncols()-1] != 1:
            raise ZeroDivisionError, "self is not invertible"
        return B.matrix_from_columns(range(self.ncols(), 2*self.ncols()))

    def charpoly(self):
        """
        Return the characteristic polynomial of this matrix.

        ALGORITHM: Compute the Hessenberg form of the matrix and read
        off the characteristic polynomial from that.  The result is
        cached.

        EXAMPLES:
            sage: A = MatrixSpace(RationalField(),3)(range(9))
            sage: A.charpoly()
            x^3 - 12*x^2 - 18*x
            sage: A.trace()
            12
            sage: A.determinant()
            0
        """
        try:
            return self.__charpoly
        except AttributeError:
            pass
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
        try:
            return self.__echelon_form
        except AttributeError:
            pass

        R = self.base_ring()
        # Fix to work with finite fields and Z/nZ, which was
        # suggested by Dan Christensen <jdc@uwo.ca>.
        if (   (is_FiniteField(R) and R.is_prime_field()) or \
               is_IntegerModRing(R)  ) and R.characteristic() < 46340:
            p = R.characteristic()
            S = sage.ext.dense_matrix_pyx.Matrix_modint(p, self.nrows(), self.ncols(), self.list())
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

    def integer_kernel(self):
        """
        Return the integral kernel of this matrix.

        Assume that the base field of this matrix has a numerator and
        denominator functions for its elements, e.g., it is the
        rational numbers or a fraction field.  This function computes
        a basis over the integers for the kernel of self.

        When kernels are implemented for matrices over general PID's,
        this function will compute kernels over PID's of matrices over
        the fraction field of the PID.  (todo)

        EXAMPLES:
            sage: A = MatrixSpace(QQ, 4)(range(16))
            sage: A.integer_kernel()
            Free module of degree 4 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -3  2]
            [ 0  1 -2  1]

        The integer kernel even makes sense for matrices with
        fractional entries:
            sage: A = MatrixSpace(QQ, 2)(['1/2',0,  0, 0])
            sage: A.integer_kernel()
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1]
        """
        d = self.denominator()
        A = self*d
        R = d.parent()
        M = matrix_space.MatrixSpace(R, self.nrows(), self.ncols())(A)
        return M.kernel()


    def is_invertible(self):
        """
        Return True if this matrix is invertible.

        EXAMPLES:

            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A.is_invertible()
            True
            sage: W = MatrixSpace(CC,2,2)
            sage: i = CC.0
            sage: B = W([1, 2+3*i, 4+5*i, 9])
            sage: B.is_invertible()
            True
            sage: N = MatrixSpace(QQ,2,2)
            sage: C = N([2,1,2,1])
            sage: C.is_invertible()
            False
        """
        return self.is_square() and self.rank() == self.nrows()

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

    def rank(self):
        try:
            return self.__rank
        except AttributeError:
            rank = self.echelon_form().rank()
            if self.is_immutable():
                self.__rank = rank
            return rank

    def _set_rank(self, r):
        self.__rank = r

    def _set_pivots(self, pivots):
        self.__pivots = pivots

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


    def row_space(self):
        """
        Return the vector space over the base field spanned by the
        rows of self.

        EXAMPLES:
            sage: A = MatrixSpace(QQ, 2,2)([1,2,3,4])
            sage: A.row_space()
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
            sage: A.row_span(IntegerRing())
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 2]
        """
        return self.row_module()

    def row_span(self, R=None):
        r"""
        Return the R-module spanned by the rows of self.

        INPUT:
            R -- (optional) principal ideal ring, such that the
                 entries of self coerce into R.  The default
                 is the base ring.
        OUTPUT:
            a free R-module.

        EXAMPLES:

        We define a $2\times 3$ matrix over $\Q$, then consider its row span
        over both $\Q$ and $\Z$.

            sage: A = MatrixSpace(QQ, 2,3)([1,2,3, '1/3', 4, 2])
            sage: A
            [  1   2   3]
            [1/3   4   2]
            sage: M1 = A.row_span(); M1
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0 12/5]
            [   0    1 3/10]
            sage: M2 = A.row_span(IntegerRing()); M2
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1/3   4   2]
            [  0  10   3]

        Note that the determinants of the inner product matrices are
        different, though their discriminants differ by a square:

            sage: M1.inner_product_matrix()
            [ 169/25   18/25]
            [  18/25 109/100]
            sage: d1 = M1.discriminant(); d1
            137/20
            sage: d2 = M2.discriminant(); d2
            685/9
            sage: d2/d1
            100/9
        """
        if R is None:
            return self.row_module()
        M = sage.modules.free_module.FreeModule(R, self.ncols())
        return M.span(self.rows())

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


#############################################
## Generic *DENSE* matrices over any field
#############################################

def _convert_dense_entries_to_list(entries):
    # Create matrix from a list of vectors
    e = []
    for v in entries:
        e += v.list()
    copy = False
    return e

class Matrix_generic_dense(Matrix):
    """
    The \\class{Matrix_generic_dense} class derives from
    \\class{Matrix}, and defines functionality for dense matrices over
    any base ring.  Matrices are represented by a list of elements in
    the base ring, and element access operations are implemented in
    this class.
    """
    def __init__(self, parent, entries=0,
                       coerce_entries=True,
                       copy=True):
        Matrix.__init__(self, parent)
        self.__nrows = parent.nrows()
        self.__ncols = parent.ncols()


        R = self.base_ring()
        if not isinstance(entries, list):
            x = R(entries)
            zero = R(0)
            entries = [zero for _ in xrange(self.nrows()*self.ncols())]
            if x != zero:
                if self.nrows() != self.ncols():
                    raise TypeError, "scalar matrix must be square"

                nc = self.ncols()
                for i in xrange(self.nrows()):
                    entries[i*nc + i] = x

        if len(entries) > 0 and hasattr(entries[0], "is_vector"):
            entries = _convert_dense_entries_to_list(entries)

        if len(entries) != self.nrows() * self.ncols():
            if len(entries) != self.nrows() * self.ncols():
                raise ArithmeticError, "entries must be a list of length %s"%\
                       (self.nrows()*self.ncols())
        if coerce_entries:
            try:
                entries = [R(x) for x in entries]
            except TypeError:
                raise TypeError, "Unable to coerce elements of entries (=%s) to %s"%(entries, R)
        elif copy:
            # Make a copy
            entries = list(entries)

        self.__entries = entries

    def __getitem__(self, ij):
        """
        INPUT:
            A[i, j] -- the i,j of A, and
            A[i]    -- the i-th row of A.
        """
        if isinstance(ij, int):
            return self.row(ij)
        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError, "index must be a 2-tuple (i,j)"
        i,j = ij
        if i < 0 or i >= self.__nrows:
            raise IndexError, "row index (i=%s) is out of range"%i
        if j < 0 or j >= self.__ncols:
            raise IndexError, "columns index (j=%s) is out of range"%j
        return self.__entries[int(i*self.ncols() + j)]

    def __setitem__(self, ij, value):
        """
        INPUT:
            A[i, j] = value -- set the (i,j) entry of A
            A[i] = value    -- set the ith row of A
        """
        self._require_mutable()
        if isinstance(ij, int):
            # Set the whole row.
            if ij < 0 or ij >= self.__nrows:
                raise IndexError, "row index (=%s) is out of range"%ij
            for j in xrange(self.ncols()):
                self[ij,j] = value[j]
            return
        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError, "index must be a 2-tuple (i,j)"
        i,j=ij
        if i < 0 or i >= self.__nrows:
            raise IndexError, "row index (i=%s) is out of range"%i
        if j < 0 or j >= self.__ncols:
            raise IndexError, "columns index (j=%s) is out of range"
        self.__entries[int(i*self.ncols() + j)] = self.base_ring()(value)

    def add_multiple_of_row(self, i, j, s):
        """
        Replace row i by s times row j.
        """
        self._require_mutable()
        s = self.base_ring()(s)
        nc = self.ncols()
        inc = int(i*nc)
        jnc = int(j*nc)
        for c in xrange(nc):
            self.__entries[inc + c] += s * self.__entries[jnc + c]

    def add_multiple_of_column(self, i, j, s):
        """
        Replace column i by s times column j.
        """
        self._require_mutable()
        s = self.base_ring()(s)
        nc = int(self.ncols())
        nr = int(self.nrows())
        i = int(i)
        j = int(j)
        for r in xrange(nr):
            self.__entries[r*nc + i] += s * self.__entries[r*nc + j]

    def _entries(self):
        return self.__entries

    def list(self):
        return list(self.__entries)

    def antitranspose(self):
        f = []
        e = self.list()
        (nc, nr) = (self.ncols(), self.nrows())
        for j in reversed(xrange(nc)):
            for i in reversed(xrange(nr)):
                f.append(e[i*nc + j])
        return self.new_matrix(nrows = nc, ncols = nr,
                               entries = f, copy=False, coerce_entries=False)

    def transpose(self):
        """
        Returns the transpose of self, without changing self.

        EXAMPLES:
        We create a matrix, compute its transpose, and note that the
        original matrix is not changed.
            sage: M = MatrixSpace(QQ,  2)
            sage: A = M([1,2,3,4])
            sage: B = A.transpose()
            sage: print B
            [1 3]
            [2 4]
            sage: print A
            [1 2]
            [3 4]
        """
        f = []
        e = self.list()
        (nc, nr) = (self.ncols(), self.nrows())
        for j in xrange(nc):
            for i in xrange(nr):
                f.append(e[i*nc + j])
        return self.new_matrix(nrows = nc, ncols = nr,
                               entries = f, copy=False, coerce_entries=False)


#############################################
## Generic sparse matrices over any field
#############################################
def _sparse_dot_product(v, w):
    """
    v and w are dictionaries with integer keys.
    """
    x = set(v.keys()).intersection(set(w.keys()))
    return sum([v[k]*w[k] for k in x])

def _convert_sparse_entries_to_dict(entries):
    e = {}
    for i in xrange(len(entries)):
        for j, x in (entries[i].dict()).iteritems():
            e[(i,j)] = x
    return e

class Matrix_generic_sparse(Matrix):
    """
    The \\class{Matrix_generic_sparse} class derives from
    \\class{Matrix}, and defines functionality for dense matrices over
    any base ring.  A generic sparse matrix is represented using a
    dictionary with keys pairs $(i,j)$ and values in the base ring.
    """
    ##WARNING: We do not check that the i,j pairs satisfy
    ##    0 <= i < nrows,  0 <= j < ncols.
    def __init__(self, parent,
                 entries=0,
                 coerce_entries=True,
                 copy=True):
        Matrix.__init__(self, parent)
        R = self.base_ring()
        self.__nrows = parent.nrows()
        self.__ncols = parent.ncols()

        if not isinstance(entries, (list, dict)):
            x = R(entries)
            zero = R(0)
            entries = {}
            if x != zero:
                if self.nrows() != self.ncols():
                    raise TypeError, "scalar matrix must be square"
                for i in xrange(self.nrows()):
                    entries[(i,i)] = x

        if isinstance(entries, list) and len(entries) > 0 and \
                hasattr(entries[0], "is_vector"):
            entries = _convert_sparse_entries_to_dict(entries)

        if isinstance(entries, list):
            if len(entries) != self.nrows() * self.ncols():
                raise TypeError, "entries has the wrong length"
            x = entries
            entries = {}
            k = 0
            for i in xrange(self.nrows()):
                for j in xrange(self.ncols()):
                    if x[k] != 0:
                        entries[(i,j)] = x[k]
                    k += 1
            copy = False

        if not isinstance(entries, dict):
            raise TypeError, "entries must be a dict"
        if coerce_entries:
            try:
                for k, x in entries.iteritems():
                    entries[k] = R(x)
            except TypeError:
                raise TypeError, "Unable to coerce entries to %s"%R
        elif copy:
            # Make a copy
            entries = dict(entries)
        self.__entries = entries
        self.__zero = R(0)

    def __getitem__(self, ij):
        """
        INPUT:
            A[i, j] -- the i,j of A, and
            A[i]    -- the i-th row of A.
        """
        if isinstance(ij, int):
            return self.row(ij)
        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError, "index must be a 2-tuple (i,j)"
        # the ij might not be int's, which would cause a problem:
        ij = (int(ij[0]), int(ij[1]))
        (i,j) = ij[0], ij[1]
        if i < 0 or i >= self.__nrows:
            raise IndexError, "row index (i=%s) is out of range"%i
        if j < 0 or j >= self.__ncols:
            raise IndexError, "columns index (j=%s) is out of range"%j
        if self.__entries.has_key(ij):
            return self.__entries[ij]
        return self.__zero

    def add_multiple_of_row(self, i, j, s):
        """
        Replace row i by s times row j.
        """
        self._require_mutable()
        s = self.base_ring()(s)

        if s == self.__zero:
            return

        i = int(i)
        j = int(j)

        for c in xrange(self.ncols()):
            jc = (j,c)
            ic = (i,c)
            if self.__entries.has_key(jc):
                iv = s*self.__entries[jc]
            else:
                continue

            if self.__entries.has_key(ic):
                self.__entries[ic] += iv
            else:
                self.__entries[ic]  = iv

    def add_multiple_of_column(self, i, j, s):
        """
        Replace column i by s times column j.
        """
        self._require_mutable()
        s = self.base_ring()(s)

        if s == self.__zero:
            return

        i = int(i)
        j = int(j)

        for r in xrange(self.nrows()):
            ri = (r,i)
            rj = (r,j)

            if self.__entries.has_key(rj):
                iv = s*self.__entries[rj]
            else:
                continue

            if self.__entries.has_key(ri):
                self.__entries[ri] += iv
            else:
                self.__entries[ri]  = iv

    def _dict(self):
        return self._entries()

    def get(self, ij):
        """
        Like __getitem__ but with no type or bounds checking.
        For (i,j) access, returns 0 if access is out of bounds.
        """
        if isinstance(ij, int):
            return self.row(ij)
        if self.__entries.has_key(ij):
            return self.__entries[ij]
        return self.__zero

    def set(self, ij, x):
        """
        Like __setitem__ but with no type or bounds checking.  Only
        works for single entries, not whole rows.
        """
        self._require_mutable()
        if x == 0:
            if self.__entries.has_key(ij):
                del self.__entries[ij]
            return
        self.__entries[ij] = x

    def __setitem__(self, ij, value):
        """
        INPUT:
            A[i, j] = value -- set the (i,j) entry of A
            A[i] = value    -- set the ith row of A
        """
        self._require_mutable()
        if isinstance(ij, int):
            # Set the whole row.
            if ij < 0 or ij >= self.__nrows:
                raise IndexError, "row index (=%s) is out of range"%ij
            try:
                for j in value.nonzero_positions():
                    self[ij,j] = value[j]
            except AttributeError:
                for j in xrange(self.ncols()):
                    self[ij,j] = value[j]
            return

        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError, "index must be a 2-tuple (i,j)"
        # the ij might not be int's, which would cause a problem, so we coerce
        ij = (int(ij[0]), int(ij[1]))
        i,j=ij
        if i < 0 or i >= self.__nrows:
            raise IndexError, "row index (i=%s) is out of range"%i
        if j < 0 or j >= self.__ncols:
            raise IndexError, "columns index (j=%s) is out of range"
        x = self.base_ring()(value)
        if x == 0:
            if self.__entries.has_key(ij):
                del self.__entries[ij]
            return
        self.__entries[ij] = x

    def __mul__(self, right):
        if isinstance(right, sage.modules.free_module_element.FreeModuleElement):
            return self.transpose().vector_matrix_multiply(right)
        if not isinstance(right, Matrix):
            return self._right_scalar_multiply(right)
        if not isinstance(right, Matrix_generic_sparse):
            right = self.matrix_space(right.nrows(), right.ncols())(right)
        if self.base_ring() != right.base_ring():
            try:
                self, right = coerce.canonical_base_coercion(self, right)
            except TypeError:
                raise ArithmeticError, "base rings must be compatible"
        if self.ncols() != right.nrows():
            raise ArithmeticError, "number of columns of self (=%s) must equal number of rows of right (=%s)."%(
                self.ncols(), right.nrows())
        nr = self.nrows()
        snc = self.ncols()
        nc = right.ncols()
        zero = self.base_ring()(0)
        E = {}
        rows = self.sparse_rows()
        cols = right.sparse_columns()
        for i in range(nr):
            row = rows[i]
            if len(row) == 0: continue
            for j in range(nc):
                x = _sparse_dot_product(row, cols[j])
                if x != 0: E[(i,j)] = x
        return self.new_matrix(nr, nc, entries = E, coerce_entries=False, copy=False)

    def hessenberg_form(self):
        """
        Return the Hessenberg form of this matrix.
        """
        # Current implementation is way too slow on sparse matrix.
        # On the dense it is reasonable.
        return self.dense_matrix().hessenberg_form()

    def sparse_columns(self):
        try:
            return self.__sparse_columns
        except AttributeError:
            C = [{} for _ in xrange(self.ncols())]
            for ij, x in self.__entries.iteritems():
                C[ij[1]][ij[0]] = x
            if self.is_immutable():
                self.__sparse_columns = C
            return C

    def sparse_rows(self):
        try:
            return self.__sparse_rows
        except AttributeError:
            R = [{} for _ in xrange(self.nrows())]
            for ij, x in self.__entries.iteritems():
                R[ij[0]][ij[1]] = x
            if self.is_immutable():
                self.__sparse_rows = R
            return R

    def __rmul__(self, left):
        left = self.base_ring()(left)
        X = {}
        for ij, x in self._entries().iteritems():
            X[ij] = left*x
        return self.new_matrix(entries=X, copy=False, coerce_entries=False)

    def denominator(self):
        R = self.base_ring()
        x = self._entries()
        if len(x) == 0:
            return integer.Integer(1)
        Z = x.iteritems()
        d = Z.next()[1].denominator()
        for _, y in Z:
            d = d.lcm(y.denominator())
        return d

    def _entries(self):
        return self.__entries

    def nonzero_positions(self):
        """
        Returns the set of pairs (i,j) such that self[i,j] != 0.
        """
        return set(self.__entries.keys())

    def matrix_from_columns(self, columns):
        """
        Return the matrix obtained from self with columns col[i] for i in
        the list of columns.

        EXAMPLES:
            sage: A = Matrix(GF(17),2, 3, range(6)); A
            [0 1 2]
            [3 4 5]
            sage: A.matrix_from_columns([0,2])
            [0 2]
            [3 5]
            sage: A.matrix_from_columns([])
            []

        The columns must all be valid:
            sage: A.matrix_from_columns([3,5])
            Traceback (most recent call last):
            ...
            IndexError: column 3 out of range
        """
        if not isinstance(columns, list):
            raise TypeError, "columns must be a list of integers"

        # Let A be this matrix and let B denote this matrix with
        # the columns deleted.
        # ALGORITHM:
        # 1. Make a table that encodes the function
        #          f : Z --> Z,
        #    f(j) = column of B where the j-th column of A goes
        # 2. Build new list of entries and return resulting matrix.
        C = set(columns)
        X = []
        j = 0
        for i in xrange(self.ncols()):
            if i in C:
                X.append(j)
                j += 1
            else:
                X.append(-1)    # column to be deleted.
        entries = {}
        E = self._entries()
        for ij in E.iterkeys():
            if ij[1] in C:
                i,j=ij
                entries[(i,X[j])] = E[ij]

        return self.new_matrix(ncols = len(columns), entries = entries,
                    copy=False, coerce_entries=False)

    def matrix_from_rows(self, rows):
        """
        Return the matrix obtained from self of rows row[i] for i in
        the list of rows.

        EXAMPLES:
            sage: A = Matrix(GF(17),3,2, range(6)); A
            [0 1]
            [2 3]
            [4 5]
            sage: A.matrix_from_rows([0,2])
            [0 1]
            [4 5]

        The matrix need not be a submatix!
            sage: A = Matrix(GF(17),3,2, range(6)); A
            [0 1]
            [2 3]
            [4 5]
            sage: A.matrix_from_rows([1,1,2,2])
            [2 3]
            [2 3]
            [4 5]
            [4 5]
        """
        if not isinstance(rows, list):
            raise TypeError, "rows must be a list of integers"
        R = set(rows)
        if not R.issubset(set(xrange(self.nrows()))):
            raise ArithmeticError, "Invalid rows."
        X = []
        i = 0
        for j in xrange(self.nrows()):
            if j in R:
                X.append(i)
                i += 1
            else:
                X.append(-1)    # row to be deleted.
        entries = {}
        E = self._entries()
        for ij in E.iterkeys():
            if ij[0] in R:
                i,j=ij
                entries[(X[i],j)] = E[ij]

        return self.new_matrix(
                    nrows = len(rows),
                    entries = entries,
                    copy=False, coerce_entries=False)

    def swap_rows(self, r1, r2):
        """
        Swap rows r1 and r2 of self.
        """
        self._require_mutable()
        nr = self.nrows()
        if r1 < 0 or r1 >= nr:
            raise IndexError, "r1 must satisfy 0<=r1<%s"%nr
        if r2 < 0 or r2 >= nr:
            raise IndexError, "r2 must satisfy 0<=r2<%s"%nr
        if r1 == r2: return
        for c in xrange(self.ncols()):
            self[r1,c], self[r2,c]  =  self[r2,c], self[r1,c]

    def transpose(self):
        """
        Returns the transpose of self, without changing self.

        EXAMPLES:
        We create a matrix, compute its transpose, and note that the
        original matrix is not changed.
            sage: M = MatrixSpace(QQ,  2, sparse=True)
            sage: A = M([1,2,3,4])
            sage: B = A.transpose()
            sage: print B
            [1 3]
            [2 4]
            sage: print A
            [1 2]
            [3 4]
        """
        X = {}
        for ij, x in self.__entries.iteritems():
            X[(ij[1],ij[0])] = x
        return self.new_matrix(nrows = self.ncols(), ncols = self.nrows(),
                           entries = X, copy=False, coerce_entries=False)


############################################################
# Generic matrices over special types of commutative rings
############################################################

class Matrix_generic_dense_domain(Matrix_domain, Matrix_generic_dense):
    def __init__(self, parent, entries=0,
                       coerce_entries=True,
                       copy=True):
        Matrix_generic_dense.__init__(self, parent, entries, coerce_entries, copy)


    def _singular_(self, singular=singular_default):
        """
        Tries to coerce this matrix to a singular matrix.
        """
        try:
            self.base_ring()._singular_()
        except (NotImplementedError, AttributeError):
            raise TypeError, "Cannot coerce to Singular"

        return singular.matrix(self.nrows(),self.ncols(),singular(self.list()))


class Matrix_generic_sparse_domain(Matrix_domain, Matrix_generic_sparse):
    def __init__(self, parent, entries=0,
                       coerce_entries=True,
                       copy=True):
        Matrix_generic_sparse.__init__(self, parent, entries, coerce_entries, copy)

    def _singular_(self, singular=singular_default):
        """
        Tries to coerce this matrix to a singular matrix.
        """
        try:
            self.base_ring()._singular_()
        except (NotImplementedError, AttributeError):
            raise TypeError, "Cannot coerce to Singular"

        As = singular.matrix(self.nrows(),self.ncols())

        sstr = "".join( "%s[%s,%s]=%s;"%(As.name(),i+1,j+1,e)
                        for (i,j),e in self._Matrix_generic_sparse__entries.iteritems() )

        singular.eval(sstr, allow_semicolon=True)

        return As


class Matrix_generic_dense_pid(Matrix_pid, Matrix_generic_dense):
    def __init__(self, parent, entries=0,
                       coerce_entries=True,
                       copy=True):
        Matrix_generic_dense.__init__(self, parent, entries, coerce_entries, copy)

class Matrix_generic_sparse_pid(Matrix_pid, Matrix_generic_sparse):
    def __init__(self, parent, entries=0,
                       coerce_entries=True,
                       copy=True):
        Matrix_generic_sparse.__init__(self, parent, entries, coerce_entries, copy)


class Matrix_generic_dense_field(Matrix_field, Matrix_generic_dense):
    def __init__(self, parent, entries=0,
                       coerce_entries=True,
                       copy=True):
        Matrix_generic_dense.__init__(self, parent, entries, coerce_entries, copy)




class Matrix_generic_sparse_field(Matrix_field, Matrix_generic_sparse):
    def __init__(self, parent, entries=0,
                       coerce_entries=True,
                       copy=True):
        Matrix_generic_sparse.__init__(self, parent, entries, coerce_entries, copy)



#############################################
## Dense matrices over the integers
#############################################
class Matrix_dense_integer(Matrix_integer, Matrix_generic_dense):
    """
    Dense matrix over the integers.

    This type is implemented mostly in Python using generic machinery,
    hence not very optimized.
    """
    def __init__(self, parent, entries=0,
                       coerce_entries=True,
                       copy=True):
        Matrix_generic_dense.__init__(self, parent, entries, coerce_entries, copy)


#############################################
## Sparse matrices over the integers
#############################################
class Matrix_sparse_integer(Matrix_integer, Matrix_generic_sparse):
    def __init__(self, parent, entries=0,
                       coerce_entries=True,
                       copy=True):
        Matrix_generic_sparse.__init__(self, parent, entries, coerce_entries, copy)


#############################################
## Generic matrices over cyclotomic fields
#############################################
class Matrix_sparse_cyclotomic(Matrix_generic_sparse_field):
    def __init__(self, parent, entries=0,
                       coerce_entries=True,
                       copy=True):
        Matrix_generic_sparse.__init__(self, parent, entries, coerce_entries, copy)

    def height(self, prec=53):
        """
        Return the height of this matrix.

        This is the maximum of the absolute values of any entries of
        self with respect to all archimedean absolute values.
        """
        K = self.base_ring()
        e = K.complex_embeddings(prec=prec)
        v = self._entries()
        return max([max([abs(f(z)) for f in e]) for z in v.itervalues()])


    def multimodular_echelon_form(self, height_guess=None, include_zero_rows=True, start_prime=3, proof=True):
        print "WARNING -- work in progress -- not finished!!!"
        d = self.denominator()
        if d != 1:
            A = d*self
        else:
            A = self

        hA = self.height()
        if height_guess is None:
            height_guess = 100000*hA**4

        if proof:
            M = self.ncols() * height_guess * hA  +  1
        else:
            M = height_guess + 1

        K = self.base_ring()
        n = K.degree()
        p = K.next_split_prime(start_prime-1)
        X = []
        best_pivots = []
        prod = 1
        f = K.defining_polynomial()

        entries = list(self._entries().iteritems())
        w = sum([z.list() for _, z in entries], [])
        B = matrix_space.MatrixSpace(integer_ring.IntegerRing(), len(entries), n)(w)
        B = B.transpose()

        while True:
            Fp = FiniteField(p)
            print 'p = ',p
            f_mod_p = f.base_extend(Fp)
            roots = f_mod_p.roots(multiplicities=False)
            print 'roots = ', roots
            M = matrix_space.MatrixSpace(Fp, n)
            v = []
            for r in roots:
                z = Fp(1)
                for i in range(n):
                    v.append(z)
                    if i < n-1:
                        z *= r    # z = z*r
            F = M(v)
            print 'F = ', F
            t = misc.cputime()
            Finv = F**(-1)
            print "time to invert", misc.cputime(t)
            # This will be *vastly* faster soon.
            print 'F^(-1) = \n', Finv

            FB = F*B
            print 'F*B = \n', FB
            MS = matrix_space.MatrixSpace(Fp, self.nrows(), self.ncols(), sparse=True)
            for i in range(n):
                row = FB.row(i)
                d = dict([(entries[j][0],row[j]) for j in range(len(entries))])
                A_mod_p = MS(d)
                print 'A modulo %sth prime is\n'%i, A_mod_p
                E_mod_p = A_mod_p.echelon_form()
                if self.nrows() == self.ncols() and E_mod_p.rank() == self.nrows():
                    # the echelon form must be the identity matrix
                    return self.parent()(1)
                print 'Echelon form is\n', E_mod_p
            break


#############################################
## Generic matrices over the rational numbers
#############################################
class Matrix_rational(Matrix_field):
    def __init__(self, parent):
        Matrix_field.__init__(self, parent)

    def _adjoint(self):
        """assumes self is a square matrix (checked in adjoint)"""
        return self.parent()(self._pari_().matadjoint().python())

    def _lllgram(self):
        """assumes self is a square matrix (checked in lllgram)"""
        Z = integer_ring.IntegerRing()
        n = self.nrows()
        # pari does not like negative definite forms
        if n > 0 and self[0,0] < 0:
            self = -self
        # maybe should be /unimodular/ matrices ?
        MS = matrix_space.MatrixSpace(Z,n,n)
        try:
            U = MS(self._pari_().lllgram().python())
        except (RuntimeError, ArithmeticError):
            raise ValueError, "not a definite matrix"
        # Fix last column so that det = +1
        if U.det() == -1:
            for i in range(n):
                U[i,n-1] = - U[i,n-1]
        return U

#############################################
## Dense matrices over the rational numbers
#############################################
class Matrix_dense_rational(Matrix_rational):
    """
    The \class{Matrix_dense_rational} class derives from
    \class{Matrix_field}, and defines functionality for dense
    matrices over the field $\Q$ of rational numbers.
    """
    def __init__(self,
                    parent,
                    entries=0,
                    coerce_entries=True,
                    copy=True):
        Matrix.__init__(self, parent)

        if isinstance(entries, dense_matrix_pyx.Matrix_rational):
            if copy:
                entries = entries.copy()
            self.__matrix = entries
            return

        if not isinstance(entries, list):
            entries = self.base_ring()(entries)

        if isinstance(entries, list) and len(entries) > 0 and \
               hasattr(entries[0],"is_vector"):
            entries = _convert_dense_entries_to_list(entries)

        if isinstance(entries, list):
            if len(entries) != parent.nrows() * parent.ncols():
                raise TypeError, "entries has wrong length"

        self.__matrix = dense_matrix_pyx.Matrix_rational(
                    parent.nrows(), parent.ncols(), entries)
        self.__entries = self.__matrix

    def __getitem__(self, ij):
        r"""
        Use \code{A[i,j]} to obtain the the $(i,j)$th entry of $A$, and
        \code{A[i]} to obtain the $i$-th row of $A$.
        """
        if isinstance(ij, int):
            return self.row(ij)
        return self.__matrix[ij]

    def __setitem__(self, ij, x):
        r"""
        Use \code{A[i,j]=x} to set the $(i,j)$th entry of $A$ to $x$,
        and \code{A[i]=v} to set the $i$th row of $A$ to the entries
        of $v$.
        """
        self._require_mutable()
        if isinstance(ij, int):
            # Set the whole row.
            if ij < 0 or ij >= self.nrows():
                raise IndexError, "row index (=%s) is out of range"%ij
            for j in range(self.ncols()):
                self.__matrix[ij,j] = x[j]
            return
        self.__matrix[ij] = x

    def __mul__(self, right):
        if isinstance(right, sage.modules.free_module_element.FreeModuleElement):
            return self.transpose().vector_matrix_multiply(right)
        if not isinstance(right, Matrix_dense_rational):
            # FIXME: this coerces too much, e.g
            #   sage: Matrix(QQ,1,1,[3]) * Matrix(GF(3),1,1,[1])
            #   [3]
            # as opposed to
            #   sage: Matrix(ZZ,1,1,[3]) * Matrix(GF(3),1,1,[1])
            #   [0]
            # IMHO, QQ should coerce down to GF(3), with a ValueError
            # if the denominator is divisible by 3, but that is not
            # how GF(3)._coerce_ is implemented...
            if isinstance(right, Matrix):
                right = self.matrix_space(right.nrows(), right.ncols())(right)
            else:
                return self._right_scalar_multiply(right)
        if self.ncols() != right.nrows():
            raise ArithmeticError, "number of columns of self (=%s) must equal number of rows of right (=%s)."%(
                self.ncols(), right.nrows())
        M = self.matrix_space(self.nrows(), right.ncols())
        A = self.__matrix * right.__matrix
        return Matrix_dense_rational(M, A, copy=False)

    def iterates(self, v, n):
        """
        Let A be this matrix.   Return a matrix with rows
        $$
           v, v*A, v*A^2, ..., v*A^(n-1).
        $$
        """
        if self.nrows() != self.ncols():
            raise ArithmeticError, "matrix must be square."
        I =  self.__matrix.iterates(list(v), n)
        M = self.matrix_space(n, self.ncols())
        return Matrix_dense_rational(M, I, copy=False)


    def hessenberg_form(self):
        """
        Return the Hessenberg form of this matrix.

        EXAMPLES:
            sage: A = Matrix(QQ, 3, 3, range(9))
            sage: A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.hessenberg_form()
            [ 0  5  2]
            [ 3 14  5]
            [ 0 -5 -2]
        """
        time = misc.verbose(t=0)
        A = self.copy()
        A.__matrix.hessenberg_form()
        misc.verbose("finished computing Hessenberg form",time)
        return A

    def charpoly(self, bound=None):
        """
        Return the characteristic polynomial of this matrix.

        INPUT:
            bound -- integer

        ALGORITHM: Use a multi-modular Hessenberg form algorithm.
        This multimodular algorithm works by first computing a bound
        B, then computing the characteristic polynomial (using
        Hessenberg form mod p) modulo enough primes so that their
        product is bigger than B.  We then uses the Chinese Remainder
        Theorem to recover the characteristic polynomial.  If the
        optional bound is specified, that bound is used for B instead
        of a potentially much worse general bound.

        EXAMPLES:
            sage: A = Matrix(QQ, 4, 4, [0, 1, -1, 0, 0, 1, -1, 1, -1, 2, -2, 1, -1, 1, 0, -1])
            sage: f = A.charpoly(); f
            x^4 + 2*x^3 - x^2 - 2*x + 1
            sage: f.factor()
            (x^2 + x - 1)^2

        Next we compute a charpoly using too low of a bound, and get an
        incorrect answer.
            sage: A = 100*Matrix(QQ, 3, 3, range(9))
            sage: A.charpoly(10)
            x^3 - 1200*x^2 + 5348*x

        Note that the incorrect answer is cached, but only with that bound:
            sage: A.charpoly()         # gives correct answer
            x^3 - 1200*x^2 - 180000*x
            sage: A.charpoly(10)       # cached incorrect answer
            x^3 - 1200*x^2 + 5348*x
        """
        if self.nrows() != self.ncols():
            raise ArithmeticError, "the matrix must be square."

        try:
            return self.__charpoly[bound]
        except AttributeError:
            self.__charpoly = {}
        except KeyError:
            pass

        if self.nrows() == 0:
            R = polynomial_ring.PolynomialRing(self.base_ring())
            f = R(0)
            if self.is_immutable():
                self.__charpoly[bound] = f
            return f

        d = self.denominator()
        if d != 1:
            A = self*d
        else:
            A = self
        if bound is None:
            # The following code (inspired by NTL's mat_poly_ZZ.c) computes
            # a bound that is vastly better than the Hadamard bound.  Victor
            # Shoup's comments 'This bound is computed via interpolation
            # through complex roots of unity'.  Reference: Mathieu and
            # Ford (1990, Section 6), 'On p-adic computation of the rational
            # form of a matrix.'.
            time = misc.verbose("computing bound")
            bound = 2
            zero = integer_ring.IntegerRing()(0)
            for i in range(A.nrows()):
                t1 = zero
                for j in range(A.ncols()):
                    t1 += (A[i,j]*A[i,j]).numerator()
                t1 += abs(A[i,i].numerator()) + 1
                if t1 > 1:
                    t1 = t1.isqrt() + 1
                bound *= t1
            bound += 1
            misc.verbose("bound = %s"%bound, time)
        else:
            bound = int(bound)

        p = 46337   # todo -- don't hardcode, so on 64-bit machine will be twice as fast!

        prod = 1
        X = []
        primes = []
        while prod < bound:
            time = misc.verbose("using p = %s"%p)
            B = A.__matrix.matrix_modint_nodenom(p)
            primes.append(p)
            fmodp = B.charpoly()   # vector of entries of the characteristic polynomial
            misc.verbose("got charpoly mod %s"%p, time)
            X.append(fmodp)
            p = sage.rings.arith.previous_prime(p)
            prod *= p

        time = misc.verbose("final CRT")
        w = sage.rings.arith.CRT_vectors(X, primes)
        misc.verbose("done with CRT",time )
        modulus = misc.prod(primes)
        m = modulus // 2
        for i in range(len(w)):
            if w[i] >= m:
                w[i] = w[i] - modulus
        g = polynomial_ring.PolynomialRing(rational_field.RationalField())(w)
        f = ((1/d)**self.nrows()) * g.rescale(d)
        if self.is_immutable():
            self.__charpoly[bound] = f
        return f

    def list(self):
        """
        Return a list of all the elements of this matrix.

        The elements of the list are the rows concatenated together.

        EXAMPLES:
            sage: A = MatrixSpace(QQ,3)(range(9))
            sage: v = A.list(); v
            [0, 1, 2, 3, 4, 5, 6, 7, 8]

        The returned list is a copy, and can be safely modified
        without changing self.
            sage: v[0] = 9999
            sage: A
            [0 1 2]
            [3 4 5]
            [6 7 8]
        """
        return list(self.__matrix.list())

    def _entries(self):
        return self.__matrix

    def _dense_matrix_mpq_(self):
        """
        Return underlying GMP-based matrix used to implement some functionality
        for this object. (Mainly for internal use.)
        """
        return self.__matrix

    def transpose(self):
        """
        Return the transpose of this matrix.

        EXAMPLES:
            sage: A = MatrixSpace(QQ,3)(range(9))
            sage: A.transpose()
            [0 3 6]
            [1 4 7]
            [2 5 8]
        """
        T = self.__matrix.transpose()
        return Matrix_dense_rational(self.matrix_space(self.ncols(), self.nrows()), T, copy=False)

    def set_row_to_multiple_of_row(self, i, j, s):
        """
        Set row i equal to s times row j.

        EXAMPLES:
            sage: A = MatrixSpace(QQ,3)(range(9))
            sage: A.set_row_to_multiple_of_row(0,1, 10)
            sage: A
            [30 40 50]
            [ 3  4  5]
            [ 6  7  8]
        """
        R = self.parent().base_ring()
        self.__matrix.set_row_to_multiple_of_row(i, j, R(s))

    def prod_of_row_sums(self, cols):
        r"""
        Calculate the product of all row sums of a submatrix of $A$ for a
        list of selected columns \code{cols}.

        This is done in C so very fast.

        AUTHOR:
            -- William Stein 2006-03-05
        """
        return self.__matrix.prod_of_row_sums(cols)

    def echelon_form(self, height_guess=None, include_zero_rows=True):
        """
        Return the echelon form of this matrix over the rational
        numbers, computed using a multi-modular algorithm.

        EXAMPLES:
            sage: A = MatrixSpace(QQ, 3)(range(9))
            sage: A.echelon_form()
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]
        """
        try:
            return self.__echelon_form[include_zero_rows]
        except AttributeError:
            self.__echelon_form = {}
        except KeyError:
            pass
        A = self.__matrix.echelon_modular(height_guess=height_guess)
        pivots = A.pivots()
        r = len(pivots)
        if not include_zero_rows:
            A = A.matrix_from_rows(range(len(pivots)))
            nr = r
        else:
            nr = self.nrows()
        E = Matrix_dense_rational(self.matrix_space(nrows=nr), A, copy=False)
        E._set_pivots(pivots)
        E._set_rank(len(pivots))
        E.set_immutable()
        if self.is_immutable():
            self.__echelon_form[include_zero_rows] = E
        return E



#############################################
## Sparse matrices over the rational numbers
#############################################
class Matrix_sparse_rational(Matrix_rational):
    r"""
    The \class{Matrix_sparse_rational} class derives from
    \class{Matrix}, and defines functionality for sparse matrices
    over the field $\Q$ of rational numbers.
    """
    def __init__(self,
                 parent,
                 entries = 0,
                 coerce_entries=True,
                 copy = True):

        Matrix.__init__(self, parent)

        if isinstance(entries, sparse_matrix_pyx.Matrix_mpq):
            if copy:
                entries = entries.copy()
            self.__matrix = entries
            return

        if isinstance(entries, list) and len(entries) > 0 and \
               hasattr(entries[0],"is_vector"):
            entries = _convert_dense_entries_to_list(entries)

        elif entries == 0:
            entries = []

        self.__matrix = sparse_matrix_pyx.Matrix_mpq(
                            parent.nrows(),
                            parent.ncols(),
                            entries,
                            coerce=coerce_entries)

    def _sparse_matrix_mpq_(self):
        return self.__matrix

    def __cmp__(self, right):
        if not isinstance(right, Matrix_sparse_rational):
            return Matrix.__cmp__(self, right)
        return self.__matrix.__cmp__(right.__matrix)

    def __getitem__(self, ij):
        if not isinstance(ij, tuple):
            return self.row(ij)
        else:
            return self.__matrix[ij]

    def __setitem__(self, ij, x):
        self._require_mutable()
        if not isinstance(ij, tuple):
            i = int(ij)
            for j in xrange(self.ncols()):
                self[i,j] = x[j]
            return
        self.__matrix[ij] = x

    def __mul__(self, B):
        if isinstance(B, sage.modules.free_module_element.FreeModuleElement):
            return self.transpose().vector_matrix_multiply(B)
        if isinstance(B, Matrix_sparse_rational):
            P = self.matrix_space(self.nrows(), B.ncols())
            return Matrix_sparse_rational(P, self.__matrix.matrix_multiply(B.__matrix),
                                          coerce_entries = False, copy=False)
        else:
            return Matrix.__mul__(self, B)

    def dense_matrix(self):
        """
        Return the dense matrix with the same entries as this sparse
        matrix.
        """
        try:
            return self.__dense_matrix
        except AttributeError:
            pass
        P = self.matrix_space(sparse=False)
        A = Matrix_dense_rational(P,
                                  entries=self.__matrix.dense_matrix(),
                                  copy = False)
        if self.is_immutable():
            self.__dense_matrix = A
            A.set_immutable()
        return A

    def echelon_form(self, height_guess=None, include_zero_rows=True, proof=True):
        """
        Return the echelon form of this sparse matrix over the
        rational numbers, computed using a sparse multi-modular
        algorithm.

        INPUT:
            height_guess --
            proof -- bool (default: True)

        The height guess is a guess for a bound on the height of the
        entries of the echelon form.  If you know for some reason that
        the entries of the echelon form are bounded, giving a good
        height bound can speed up this function.  At the end of the
        computation the result is checked for correctness and the
        height bound increased if necessary, so giving too small of a
        guess will not lead to incorrect results (though it may slow
        down the algorithm).

        EXAMPLES:
            sage: A = MatrixSpace(QQ, 3, sparse=True)(range(9))
            sage: A.echelon_form()
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]
            sage: A = 9999999*MatrixSpace(QQ, 3, sparse=True)(range(9))
            sage: A
            [       0  9999999 19999998]
            [29999997 39999996 49999995]
            [59999994 69999993 79999992]
            sage: A.echelon_form(height_guess=1)
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]
        """
        try:
            return self.__echelon_form[include_zero_rows]
        except AttributeError:
            self.__echelon_form = {}
        except KeyError:
            pass

        # TODO: If we know the echelon form with/without-out zero rows,
        # we can trivially compute the other one!  This would be an
        # easy optimization to include.

        if self.nrows() == 0:
            E = self.copy()
            E._set_pivots ([])
            E._set_rank(0)
        else:
            X = self.__matrix.echelon_multimodular(
                          height_guess = height_guess, proof = proof)
            pivots = X.pivots()
            r  = len(pivots)
            if not include_zero_rows:
                X  = X.matrix_from_rows(range(r))
                nr = r
            else:
                nr = self.nrows()
            E = Matrix_sparse_rational(self.matrix_space(nrows=nr), X,
                                       coerce_entries=False, copy=False)
            E._set_pivots(pivots)
            E._set_rank(r)

        if self.is_immutable():
            self.__echelon_form[include_zero_rows] = E
        E.set_immutable()
        return E

    def hessenberg_form(self):
        r"""
        Return the Hessenberg form of this sparse matrix.

        ALGORITHM: Compute the Hessenberg form of the corresponding
        dense matrix, obtained using \code{self.dense_matrix()}

        EXAMPLES:
            sage: A = MatrixSpace(QQ, 3, sparse=True)(range(9))
            sage: H = A.hessenberg_form(); H
            [ 0  5  2]
            [ 3 14  5]
            [ 0 -5 -2]
            sage: H.is_sparse()
            True
        """
        return self.dense_matrix().hessenberg_form().sparse_matrix()

    def charpoly(self, bound=None):
        """
        Return the characteristic polynomial of this matrix.

        See the documentation for self.dense_matrix().charpoly
        for more details.

        ALGORITHM: Compute the charpoly of the corresponding
        dense matrix, obtained using \code{self.dense_matrix()}.

        EXAMPLES:
            sage: A = MatrixSpace(QQ,4, sparse=True)(range(16))
            sage: A.charpoly()
            x^4 - 30*x^3 - 80*x^2
        """
        return self.dense_matrix().charpoly(bound = bound)

    def transpose(self):
        return self.dense_matrix().transpose().sparse_matrix()
        #P = self.matrix_space(self.ncols(), self.nrows())
        #return Matrix_sparse_rational(P, self.__matrix.transpose(),
        #                              coerce_entries = False, copy=False)

    def set_row_to_multiple_of_row(self, i, j, s):
        self._require_mutable()
        self.__matrix.set_row_to_multiple_of_row(i, j, rational.Rational(s))

    def _set_row_to_negative_of_row_of_A_using_subset_of_columns(self, i, A, r, cols):
        self.__matrix.set_row_to_negative_of_row_of_A_using_subset_of_columns(i, A.__matrix, r, cols)



