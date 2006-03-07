r"""
Matrices

Elements of matrix spaces are of class \code{Matrix}.  They can be
either sparse or dense, and can be defined over any base ring.

EXAMPLES:

We create the $2\times 3$ matrix
$$
   \left(\begin{matrix} 1&2&3\\4&5&6 \end{matrix}\right)
$$
as an element of a matrix space over $\Q$:

    sage: M = MatrixSpace(RationalField(),2,3)
    sage: A = M([1,2,3, 4,5,6])
    sage: A
    [1 2 3]
    [4 5 6]
    sage: A.parent()
    Full MatrixSpace of 2 by 3 dense matrices over Rational Field

We next change the top-right entry of $A$.  Note that matrix indexing
is $0$-based in SAGE, so the top right entry is $(0,2)$, which should
be thought of as ``row number $0$, column number 2''.

    sage: A[0,2] = 389
    sage: A
    [  1   2 389]
    [  4   5   6]

Also notice how matrices print.  All columns have the same width and
entries in a given column are right justified.   Next we compute the
reduced row echelon form of $A$.

    sage: A.echelon_form()
    [      1       0 -1933/3]
    [      0       1  1550/3]


MUTABILITY: Matrices are either immutable or not.  When initially
created, matrices are typically mutable, so one can change their
entries.  However, nothing about mutable matrices is cached.  Once a
matrix $A$ is made immutable using \code{A.set_immutable()} the
entries of $A$ cannot be changed.  However, properies of $A$ such as
its rank, characteristic polynomial, etc., are all cached so
computations involving $A$ may be more efficient.  Once $A$ is made
immutable it cannot be changed back.  However, one can obtain a
mutable copy of $A$ using \code{A.copy()}.

The echelon form method always returns immutable matrices with known
rank.
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
import sage.matrix.matrix_space
import matrix_space
import sage.libs.pari.all as pari
import sage.modules.module_element  as module_element
import berlekamp_massey
import sage.rings.polynomial_ring as polynomial_ring
import sage.rings.integer_ring as integer_ring
import sage.rings.rational_field as rational_field
import sage.rings.rational as rational
import sage.rings.number_field.number_field as number_field
import sage.rings.coerce as coerce
from sage.rings.all import is_FiniteField

from sage.structure.mutability import Mutability

import sage.ext.dense_matrix_pyx

cputime = misc.cputime

def is_Matrix(x):
    return isinstance(x, Matrix)

class Matrix(module_element.ModuleElement, Mutability):
    r"""
    The \class{Matrix} class is the base class for all matrix
    classes.  To create a \class{Matrix}, first create a
    \class{MatrixSpace}, then coerce a list of elements into the
    \class{MatrixSpace}.  See the documentation of
    \class{MatrixSpace} for more details.
    """
    def __init__(self, parent):
        """
        EXAMPLES:
            sage: import sage.matrix.matrix
            sage: A = sage.matrix.matrix.Matrix(MatrixSpace(QQ,2))
            sage: type(A)
            <class 'sage.matrix.matrix.Matrix'>
        """
        module_element.ModuleElement.__init__(self, parent)
        self.__nrows = parent.nrows()
        self.__ncols = parent.ncols()
        Mutability.__init__(self, False)

    def _matrix_(self, R):
        return self.change_ring(R)

    def __repr__(self):
        nr = self.nrows()
        nc = self.ncols()
        if nr == 0 or nc == 0:
            return "[]"
        # compute column widths
        S = [str(x) for x in self.list()]
        width = max([len(x) for x in S])
        rows = []
        m = 0

        # compute rows
        for r in xrange(nr):
            s = ""
            for c in xrange(nc):
                if c == nc-1:
                    sep=""
                else:
                    sep=" "
                entry = S[r*nc+c]
                if c == 0:
                    m = max(m, len(entry))
                entry = " "*(width-len(entry)) + entry
                s += entry + sep
            rows.append(s)

        # Remove leading spaces
##        n = width-m
##        if n > 0:
##            for r in xrange(nr):
##                rows[r] = rows[r][n:]

        # Put brackets around in a single string
        s = "\n".join(["[%s]"%row for row in rows])
        return s

    def _latex_sparse(self, variable="x"):
        """
        Return a latex string that represents this matrix as a sparse
        matrix.  The rows are printed as sums sum a_i x_i, where
        x is the variable.
        """
        nr = self.nrows()
        nc = self.ncols()
        s = "\\left(\\begin{align*}\n"
        for i in xrange(nr):
            v = [(j, self[i,j]) for j in xrange(nc) if self[i,j] != 0]
            for j in xrange(len(v)):
                s += "%s*%s_{%s}"%(v[j][1], variable, v[j][0])
                if j == 0:
                    s += "& + "
                elif j < len(v) - 1:
                    s += " + "
                else:
                    s += "\\\\\n"
        s += "\n\\end{align*}"
        return s

    def _latex_(self):
        nr = self.nrows()
        nc = self.ncols()
        if nr == 0 or nc == 0:
            return "()"
        # compute column widths
        S = [str(x) for x in self.list()]
        width = max([len(x) for x in S])
        S = self.list()
        rows = []
        m = 0

        # compute rows
        for r in xrange(nr):
            s = ""
            for c in xrange(nc):
                if c == nc-1:
                    sep=""
                else:
                    sep="&"
                entry = latex.latex(S[r*nc+c])
                if c == 0:
                    m = max(m, len(entry))
                s += entry + sep
            rows.append(s)

        # Put brackets around in a single string
        s = "\\\\\n".join(["%s"%row for row in rows])
        return "\\left(\\begin{array}{%s}\n"%('r'*nc) + s + "\n\\end{array}\\right)"

    def __str__(self):
        return self.__repr__()

    def __hash__(self):
        return hash(self.__str__())


    ###################################################
    ## Coercion to PARI
    ###################################################
    #def _pari_(self):
    #    return pari.pari.matrix(
    #        self.nrows(), self.ncols(), self.list())

    def _pari_init_(self):
        w = self.list()
        nr = self.nrows(); nc = self.ncols()
        v = [','.join([w[i*nc + j]._pari_init_() for j in range(nc)])
                      for i in range(nr)]
        return '[%s]'%(';'.join(v))

    ###################################################
    ## Coercion to GAP
    ###################################################
    def _gap_init_(self):
        """
        EXAMPLES:
            sage: A = MatrixSpace(QQ,3)([1,2,3,4/3,5/3,6/4,7,8,9])
            sage: g = gap(A); g
            [ [ 1, 2, 3 ], [ 4/3, 5/3, 3/2 ], [ 7, 8, 9 ] ]
            sage: g.IsMatrix()
            true
        """
        v = ['[%s]'%(','.join([self.get((i,j))._gap_init_() for j in xrange(self.ncols())])) for
             i in xrange(self.nrows())]
        return '[%s]'%(','.join(v))

    ###################################################
    ## Coercion to MAGMA
    ###################################################
    def _magma_init_(self):
        r"""
        EXAMPLES:
        We first coerce a square matrix.
            sage: A = MatrixSpace(QQ,3)([1,2,3,4/3,5/3,6/4,7,8,9])
            sage: B = magma(A); B                       # optional
            [  1   2   3]
            [4/3 5/3 3/2]
            [  7   8   9]
            sage: B.Type()                              # optional
            AlgMatElt
            sage: B.Parent()                            # optional
            Full Matrix Algebra of degree 3 over Rational Field

        We coerce a non-square matrix over $\Z/8\Z$.
            sage: A = MatrixSpace(Integers(8),2,3)([-1,2,3,4,4,-2])
            sage: B = magma(A); B                       # optional
            [7 2 3]
            [4 4 6]
            sage: B.Type()                              # optional
            ModMatRngElt
            sage: B.Parent()                            # optional
            Full RMatrixSpace of 2 by 3 matrices over IntegerRing(8)
        """
        K = self.base_ring()._magma_init_()
        if self.is_square():
            s = 'MatrixAlgebra(%s, %s)'%(K, self.nrows())
        else:
            s = 'RMatrixSpace(%s, %s, %s)'%(K, self.nrows(), self.ncols())
        v = [x._magma_init_() for x in self.list()]
        return s + '![%s]'%(','.join(v))

    ###################################################
    ## Access:
    ##   These *must* be overridden by each subclass
    ###################################################
    def __getitem__(self, ij):
        """
        Return element or row of self.

        INPUT:
            ij -- tuple (i,j) with i, j ints
        or
            ij -- int

        USAGE:
            A[i, j] -- the i,j of A, and
            A[i]    -- the i-th row of A.

        EXAMPLES:
            sage: A = Matrix(Integers(2006),2,2,[-1,2,3,4])
            sage: A[0,0]
            2005
            sage: A[0]
            (2005, 2)
        """
        raise NotImplementedError

    def get(self, ij):
        """
        Entry access, but potentially faster since it might be without
        bounds checking.

        INPUT:
            ij -- tuple (i,j) with i, j ints

        EXAMPLES:
            sage: A = Matrix(Integers(2006),2,2,[-1,2,3,4])
            sage: A.get(0,1)
            Traceback (most recent call last):
            ...
            TypeError: get() takes exactly 2 arguments (3 given)
            sage: A.get((0,1))
            2
        """
        return self[ij]

    def __setitem__(self, ij, x):
        """
        INPUT:
            ij -- tuple
            x -- object

        USAGE:
            A[i, j] = x -- set the (i,j) entry of A
            A[i]    = x -- set the ith row of A

        EXAMPLES:
            sage: A = Matrix(Integers(2006),2,2,[-1,2,3,4]); A
            [2005    2]
            [   3    4]
            sage: A[0,0] = 5; A
            [5 2]
            [3 4]
            sage: A[0] = [2,3]
            sage: A
            [2 3]
            [3 4]
            sage: A.set_immutable()
            sage: A[0,0] = 7
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
        """
        raise NotImplementedError

    def set(self, ij, x):
        """
        Set entry, but potentially faster because it might be without
        bounds checking.

        INPUT:
            ij -- tuple (i,j) with i, j ints
            x -- object

        EXAMPLES:
            sage: A = Matrix(Integers(2006),2,2,[-1,2,3,4])
            sage: A.set(0, 1, 5)
            Traceback (most recent call last):
            ...
            TypeError: set() takes exactly 3 arguments (4 given)
            sage: A.set((0,1), 5)
            sage: A
            [2005    5]
            [   3    4]
        """
        self[ij] = x

    ###################################################
    ## Convenience
    ###################################################
    def matrix_space(self, nrows=None, ncols=None, sparse=None):
        if nrows == None: nrows = self.nrows()
        if ncols == None: ncols = self.ncols()
        if sparse is None:
            sparse = sparse=self.is_sparse()
        return self.parent().matrix_space(nrows, ncols, sparse=sparse)

    def new_matrix(self, nrows=None, ncols=None, entries=0,
                coerce_entries=True, copy=True, sparse=None):
        """
        Create a matrix in the parent of this space with the given
        number of rows, columns, etc.

        WARNING: This function called with no arguments returns the 0
        matrix by default, not the matrix self.
        """
        return self.matrix_space(nrows, ncols, sparse=sparse).matrix(
                entries, coerce_entries, copy)

    ###################################################
    ## Properties
    ###################################################
    def ncols(self):
        """
        Return the number of rows of this matrix.

        EXAMPLES:
            sage: M = MatrixSpace(QQ, 2, 3)
            sage: A = M([1,2,3, 4,5,6])
            sage: A
            [1 2 3]
            [4 5 6]
            sage: A.ncols()
            3
            sage: A.nrows()
            2

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        return self.__ncols

    def nrows(self):
        r"""
        Return the number of rows of this matrix.

        EXAMPLES:
            sage: M = MatrixSpace(RationalField(),6,7)
            sage: A = M([1,2,3,4,5,6,7, 22,3/4,34,11,7,5,3, 99,65,1/2,2/3,3/5,4/5,5/6, 9,8/9, 9/8,7/6,6/7,76,4, 0,9,8,7,6,5,4, 123,99,91,28,6,1024,1])
            sage: A
            [   1    2    3    4    5    6    7]
            [  22  3/4   34   11    7    5    3]
            [  99   65  1/2  2/3  3/5  4/5  5/6]
            [   9  8/9  9/8  7/6  6/7   76    4]
            [   0    9    8    7    6    5    4]
            [ 123   99   91   28    6 1024    1]
            sage: A.ncols()
            7
            sage: A.nrows()
            6

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        return self.__nrows


    ###################################################
    ## Functions
    ###################################################
    def act_on_polynomial(self, f):
        """
        Returns the polynomial f(self*x).

        INPUT:
            self -- an nxn matrix
            f -- a polynomial in n variables x=(x1,...,xn)

        OUTPUT:
            The polynomial f(self*x).

        EXAMPLES:
            sage: R = MPolynomialRing(RationalField(), 2, ['x','y'])
            sage: x, y = R.gens()
            sage: f = x**2 - y**2
            sage: M = MatrixSpace(RationalField(),2)
            sage: A = M([1,2,3,4])
            sage: A.act_on_polynomial(f)
            -12*y^2 - 20*x*y - 8*x^2
        """
        if not self.is_square():
            raise ArithmeticError, "self must be square but is %s x %s"%(self.nrows(), self.ncols())

        F = f.base_ring()
        vars = f.parent().gens()
        n = len(self.rows())
        ans = tuple(sum([self.get( (i,j) )*vars[j] for j in range(n)]) for i in range(n))
        return f(ans)

    def add_multiple_of_row(self, i, j, s):
        """
        Replace row i by s times row j.
        """
        self._require_mutable()
        s = self.base_ring()(s)
        for c in xrange(self.ncols()):
            self[i,c] += s*self[j,c]

    def add_multiple_of_column(self, i, j, s):
        """
        Replace column i by s times column j.
        """
        self._require_mutable()
        s = self.base_ring()(s)
        for r in xrange(self.nrows()):
            self[r,i] += s*self[r,j]

    def augment(self, other):
        """
        Return the augmented matrix of the form [self | other].

        EXAMPLES:
            sage: M = MatrixSpace(RationalField(),2,2)
            sage: A = M([1,2, 3,4])
            sage: A
            [1 2]
            [3 4]
            sage: N = MatrixSpace(RationalField(),2,1)
            sage: B = N([9,8])
            sage: B
            [9]
            [8]
            sage: A.augment(B)
            [1 2 9]
            [3 4 8]
            sage: B.augment(A)
            [9 1 2]
            [8 3 4]
            sage: M = MatrixSpace(RationalField(),3,4)
            sage: A = M([1,2,3,4, 0,9,8,7, 2/3,3/4,4/5,9/8])
            sage: A
            [  1   2   3   4]
            [  0   9   8   7]
            [2/3 3/4 4/5 9/8]
            sage: N = MatrixSpace(RationalField(),3,2)
            sage: B = N([1,2, 3,4, 4,5])
            sage: B
            [1 2]
            [3 4]
            [4 5]
            sage: A.augment(B)
            [  1   2   3   4   1   2]
            [  0   9   8   7   3   4]
            [2/3 3/4 4/5 9/8   4   5]
            sage: B.augment(A)
            [  1   2   1   2   3   4]
            [  3   4   0   9   8   7]
            [  4   5 2/3 3/4 4/5 9/8]

        AUTHORS:
            -- Naqi Jaffery (2006-01-24): examples
        """
        if not isinstance(other, Matrix):
            raise TypeError, "other must be a Matrix"
        if self.base_ring() != other.base_ring():
            raise TypeError, "base rings must be the same"
        if self.nrows() != other.nrows():
            raise TypeError, "number of rows must be the same"
        Z = self.new_matrix(ncols = self.ncols() + other.ncols())
        for r in xrange(self.nrows()):
            for c in xrange(self.ncols()):
                Z[r,c] = self[r,c]
        nc = self.ncols()
        for r in xrange(other.nrows()):
            for c in xrange(other.ncols()):
                Z[r,c + nc] = other[r,c]
        return Z

    def block_sum(self, other):
        """
        Return the block matrix that has self and other on the diagonal:
        [self |    0  ]
        [  0  | other ]
        """
        if not isinstance(other, Matrix):
            raise TypeError, "other must be a Matrix"
        top = self.augment(self.new_matrix(ncols=other.ncols()))
        bottom = other.new_matrix(ncols=self.ncols()).augment(other)
        return top.stack(bottom)

    def change_ring(self, ring):
        """
        Return the matrix obtained by coercing the entries of this
        matrix into the given ring.
        """
        if ring == self.base_ring():
            return self
        M = sage.matrix.matrix_space.MatrixSpace(ring, self.__nrows, self.__ncols)
        return M(self.list())

    def column(self, i):
        """
        Return the i-th column of this matrix.
        """
        V = sage.modules.free_module.FreeModule(self.base_ring(), self.ncols())
        return V([self[j,i] for j in xrange(self.nrows())])

    def columns(self):
        """
        Returns the list of columns of self, as vectors.

        EXAMPLES:
            sage: A = MatrixSpace(RationalField(), 2,3) (xrange(6))
            sage: A
            [0 1 2]
            [3 4 5]
            sage: A.columns()
            [(0, 3), (1, 4), (2, 5)]
        """
        return self.transpose().rows()

    def commutator(self, other):
        """
        Return the commutator self*other - other*self.
        """
        return self*other - other*self

    def copy(self):
        """
        Return a copy of this matrix.  Changing the entries of the
        copy will not change the entries of this matrix.
        """
        return self.new_matrix(entries=self._entries(), coerce_entries=False)

    def dense_matrix(self):
        """
        If this matrix is sparse, return a dense matrix with the same
        entries.  If this matrix is dense, return this matrix (not a
        copy).

        NOTE: The definition of"dense" and "sparse" in SAGE have
        nothing to do with the number of nonzero entries.  Sparse
        and dense are properties of the underlying representation
        of the matrix.

        EXAMPLES:
            sage: A = MatrixSpace(RationalField(),2, sparse=True)([1,2,0,1])
            sage: A.is_sparse()
            True
            sage: B = A.dense_matrix()
            sage: B.is_sparse()
            False
            sage: A*B
            [1 4]
            [0 1]
            sage: A.parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
            sage: B.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: (A*B).parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
            sage: (B*A).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        """
        if self.is_dense():
            return self
        M = self.matrix_space(self.nrows(), self.ncols(), sparse=False)
        return M(self.list())

    def sparse_matrix(self):
        """
        If this matrix is dense, return a sparse matrix with
        the same entries.  If this matrix is sparse, return this
        matrix (not a copy).

        NOTE: The definition of "dense" and "sparse" in SAGE have
        nothing to do with the number of nonzero entries.  Sparse
        and dense are properties of the underlying representation
        of the matrix.

        EXAMPLES:
            sage: A = MatrixSpace(RationalField(),2, sparse=False)([1,2,0,1])
            sage: A.is_sparse()
            False
            sage: B = A.sparse_matrix()
            sage: B.is_sparse()
            True
            sage: A
            [1 2]
            [0 1]
            sage: B
            [1 2]
            [0 1]
            sage: A*B
            [1 4]
            [0 1]
            sage: A.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: B.parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
            sage: (A*B).parent()
            Full MatrixSpace of 2 by 2 dense matrices over Rational Field
            sage: (B*A).parent()
            Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        """
        if self.is_sparse():
            return self
        M = self.matrix_space(self.nrows(), self.ncols(), sparse=True)
        return M(self._dict())

    def dense_row(self, n):
        if n < 0 or n >= self.nrows():
            raise IndexError, "n must be between 0 and the number of rows"
        V = sage.modules.free_module.VectorSpace(\
                 self.base_ring(), self.ncols(), sparse=False)
        return V([self[n, i] for i in xrange(self.ncols())])

    def _dict(self):
        try:
            return self.__dict
        except AttributeError:
            v = {}
            for ij in self.nonzero_positions():
                v[ij] = self.get(ij)
            if self.is_immutable():
                self.__dict = v
            return v

    def _entries(self):
        """
        This must be defined in the derived class.  It is the
        underlying representation of elements, which may be a list or
        dict or something else.
        """
        raise NotImplementedError

    def is_dense(self):
        return self.parent().is_dense()

    def is_sparse(self):
        return self.parent().is_sparse()

    def is_square(self):
        return self.nrows() == self.ncols()

    def list(self):
        return [self[i,j] for i in xrange(self.nrows()) for j in xrange(self.ncols())]

    def nonpivots(self):
        """
        Return the list of i such that the i-th column of self
        is NOT a pivot column of the reduced row echelon form
        of self.

        OUTPUT:
            list -- sorted list of integers
        """
        X = set(self.pivots())
        return [j for j in xrange(self.ncols()) if not (j in X)]

    def nonzero_positions(self):
        """
        Returns the set of pairs (i,j) such that self[i,j] != 0.
        """
        z = self.base_ring()(0)
        return set([(i,j) for i in xrange(self.nrows()) \
                          for j in xrange(self.ncols()) \
                          if self[i,j] != z])

    def nonzero_positions_in_column(self, i):
        """
        Return the integers j such that self[j,i] is nonzero, i.e.,
        such that the j-th position of the i-th column is nonzero.
        """
        z = self.base_ring()(0)
        return set([j for j in xrange(self.nrows()) if self[j,i] != z])

    def nonzero_positions_in_row(self, i):
        """
        Return the integers j such that self[i,j] is nonzero, i.e.,
        such that the j-th position of the i-th row is nonzero.
        """
        z = self.base_ring()(0)
        return set([j for j in xrange(self.ncols()) if self[i,j] != z])

    def prod_of_row_sums(self, cols):
        r"""
        Calculate the product of all row sums of a submatrix of $A$ for a
        list of selected columns \code{cols}.

        AUTHOR:
            -- Jaap Spies (2006-02-18)
        """
        pr = 1
        for row in xrange(self.nrows()):
            pr *= sum([self[row, c] for c in cols])
        return pr

    def permanent(self):
        r"""
        Calculate and return the permanent of this $m \times n$ matrix using
        Ryser's algorithm.

        Let $A = (a_{i,j})$ be an $m \times n$ matrix over any
        commutative ring, with $m \le n$.   The permanent of $A$ is
        \[
        \text{per}(A) = \sum_\pi a_{1,\pi(1)}a_{2,\pi(2)} \cdots a_{m\,pi(m)}
        \]
        where the summation extends over all one-to-one functions $\pi$ from
        $\{1, \ldots, m\}$ to $\{1, \ldots, n\}$.

        The product $ a_{1,\pi(1)}a_{2,\pi(2)} \cdots a_{m,\pi(m)}$ is called
        diagonal product. So the permanent of an $m \times n$ matrix $A$ is the
        sum of all the diagonal products of $A$.

        Modification of theorem 7.1.1. from Brualdi and Ryser:
        Combinatorial Matrix Theory.
        Instead of deleting columns from $A$, we choose columns from $A$ and
        calculate the product of the row sums of the selected submatrix.

        INPUT:
            A -- matrix of size m x n with m <= n

        OUTPUT:
            permanent of matrix A

        EXAMPLES:
            sage: M = MatrixSpace(ZZ,4,4)
            sage: A = M([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
            sage: A.permanent()
            24

            sage: M = MatrixSpace(QQ,3,6)
            sage: A = M([1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
            sage: A.permanent()
            36

            sage: M = MatrixSpace(RR,3,6)
            sage: A = M([1.0,1.0,1.0,1.0,0,0,0,1.0,1.0,1.0,1.0,0,0,0,1.0,1.0,1.0,1.0])
            sage: A.permanent()
            36.000000000000000

        See Sloane's sequence OEIS A079908(3) = 36, "The Dancing School Problems"

            sage: print sloane_sequence(79908)                # optional (internet connection)
            Looking up in Sloane's online database...
            [79908, 'Solution to the Dancing School Problem with 3 girls: f(3,n).', [1, 4, 14, 36, 76, 140, 234, 364, 536, 756, 1030, 1364, 1764, 2236, 2786, 3420, 4144, 4964, 5886, 6916, 8060, 9324, 10714, 12236, 13896, 15700, 17654, 19764, 22036, 24476, 27090, 29884, 32864, 36036, 39406, 42980, 46764, 50764, 54986, 59436]]

            sage: M = MatrixSpace(ZZ,4,5)
            sage: A = M([1,1,0,1,1,0,1,1,1,1,1,0,1,0,1,1,1,0,1,0])
            sage: A.permanent()
            32

            See Minc: Permanents, Example 2.1, p. 5.

            sage: M = MatrixSpace(QQ,2,2)
            sage: A = M([1/5,2/7,3/2,4/5])
            sage: A.permanent()
            103/175

            sage: R = PolynomialRing(IntegerRing(),'a'); a = R.gen()
            sage: A = MatrixSpace(R,2)([[a,1], [a,a+1]])
            sage: A.permanent()
            a^2 + 2*a

            sage: R = MPolynomialRing(IntegerRing(),2); x,y = R.gens()
            sage: A = MatrixSpace(R,2)([x, y, x^2, y^2])
            sage: A.permanent()
            x_0*x_1^2 + x_0^2*x_1


        AUTHOR:
            -- Jaap Spies (2006-02-16)
                Copyright (C) 2006 Jaap Spies <j.spies@hccnet.nl>
                Copyright (C) 2006 William Stein <wstein@ucsd.edu>
            -- Jaap Spies (2006-02-21): added definition of permanent

        NOTES:
            -- Currently optimized for dense matrices over QQ.
        """
        perm = 0
        m, n = self.nrows(), self.ncols()
        if not m <= n:
            raise ValueError, "must have m <= n, but m (=%s) and n (=%s)"%(m,n)
        for r in range(1, m+1):
            lst = _combinations(range(n), r)
            s = sum([self.prod_of_row_sums(cols) for cols in lst])
            perm += (-1)**(m-r) * sage.rings.arith.binomial(n-r, m-r) * s
        return perm

    def permanental_minor(self, k):
        r"""
        Calculates the permanental $k$-minor of a $m \times n$ matrix.

        This is the sum of the permanents of all possible $k$ by $k$
        submatices of $A$.

        See Brualdi and Ryser: Combinatorial Matrix Theory, p. 203.
        Note the typo $p_0(A) = 0$ in that reference!  For
        applications see Theorem 7.2.1 and Theorem 7.2.4.

        Note that the permanental $m$-minor equals $per(A)$.

        For a (0,1)-matrix $A$ the permanental $k$-minor counts the
        number of different selections of $k$ 1's of $A$ with no two
        of the 1's on the same line.

        INPUT:
            self -- matrix of size m x n with m <= n

        OUTPUT:
            permanental k-minor of matrix A

        EXAMPLES:
            sage: M = MatrixSpace(ZZ,4,4)
            sage: A = M([1,0,1,0,1,0,1,0,1,0,10,10,1,0,1,1])
            sage: A.permanental_minor(2)
            114

            sage: M = MatrixSpace(ZZ,3,6)
            sage: A = M([1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
            sage: A.permanental_minor(0)
            1
            sage: A.permanental_minor(1)
            12
            sage: A.permanental_minor(2)
            40
            sage: A.permanental_minor(3)
            36

        Note that if k == m the permanental k-minor equals per(A)

            sage: A.permanent()
            36

            sage: A.permanental_minor(5)
            0

        For C the "complement" of A:

            sage: M = MatrixSpace(ZZ,3,6)
            sage: C = M([0,0,0,0,1,1,1,0,0,0,0,1,1,1,0,0,0,0])
            sage: m, n = 3, 6
            sage: sum([(-1)^k * C.permanental_minor(k)*factorial(n-k)/factorial(n-m) for k in range(m+1)])
            36

            See Theorem 7.2.1 of Brualdi: and Ryser: Combinatorial Matrix Theory: per(A)

        AUTHOR:
            - Jaap Spies (2006-02-19)
        """
        m, n = self.nrows(), self.ncols()
        if not m <= n:
            raise ValueError, "must have m <= n, but m (=%s) and n (=%s)"%(m,n)

        if k == 0:
            return 1
        if k > m:
            return 0
        k = int(k)
        pm = 0
        for cols in _combinations(range(n),k):
            for rows in _combinations(range(m),k):
                pm += self.matrix_from_rows_and_columns(rows, cols).permanent()
        return pm

    def rook_vector(self, check = False):
        r"""
        Returns rook vector of this matrix.

        Let $A$ be a general $m$ by $n$ (0,1)-matrix with $m \le n$.
        We identify $A$ with a chessboard where rooks can be placed on
        the fields corresponding with $a_{ij} = 1$. The number $r_k =
        p_k(A)$ (the permanental $k$-minor) counts the number of ways
        to place $k$ rooks on this board so that no two rooks can
        attack another.

        The rook vector is the list consisting of $r_0, r_1, \ldots, r_m$.

        The rook polynomial is defined by $r(x) = \sum_{k=0}^m r_k x^k$.

        INPUT:
            self -- m by n matrix with m <= n
            check -- True or False (default), optional

        OUTPUT:
            rook vector

        EXAMPLES:
            sage: M = MatrixSpace(ZZ,3,6)
            sage: A = M([1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
            sage: A.rook_vector()
            [1, 12, 40, 36]

            sage: x = PolynomialRing(IntegerRing(),'x').gen()
            sage: rv = A.rook_vector()
            sage: rook_polynomial = sum([rv[k] * x^k for k in range(len(rv))])
            sage: rook_polynomial
            36*x^3 + 40*x^2 + 12*x + 1

        AUTHOR:
            - Jaap Spies (2006-02-24)
        """
        m, n = self.nrows(), self.ncols()
        if not m <= n:
            raise ValueError, "must have m <= n, but m (=%s) and n (=%s)"%(m,n)

        if check == True:
            # verify that self[i, j] in {0, 1}
            for i in range(m):
                for j in range(n):
                    if not (self[i,j] == 0 or self[i,j] == 1):
                        raise ValueError, "must have zero or one, but we have (=%s)"%(self[i,j])

        return [self.permanental_minor(k) for k in range(m+1)]


    def determinant(self):
        r"""
        Return the determinant of self.

        ALGORITHM: This is computed using a \emph{naive generic
        algorithm}.  For matrices over more most rings more
        sophisticated algorithms can be used.  (Type
        \code{A.determinant?} to see what is done for a specific
        matrix A.)  If you're actually using the algorithm below, you
        should probably clamor for something better to be implemented
        for your base ring.

        EXAMPLES:
            sage: A = MatrixSpace(Integers(8),3)([1,7,3, 1,1,1, 3,4,5])
            sage: A.determinant()
            6
        """
        try:
            return self.__determinant
        except AttributeError:
            pass
        if not self.is_square():
            raise ValueError, "self must be square but is %s x %s"%(
                               self.nrows(), self.ncols())
        n = self.nrows()
        R = self.parent().base_ring()
        if n == 0:
            return R(1)

        d = R(0)
        s = R(1)
        A = self.matrix_from_rows(range(1, n))
        for i in range(n):
            v = range(n)
            del v[i]
            B = A.matrix_from_columns(v)
            d += s*self.get((0,i))*B.determinant()
        self.__determinant = d
        return d

    # shortcuts
    def det(self, *args, **kwds):
        """
        Synonym for self.determinant(...).

        EXAMPLES:
            sage: A = MatrixSpace(Integers(8),3)([1,7,3, 1,1,1, 3,4,5])
            sage: A.det()
            6
        """
        return self.determinant(*args, **kwds)

    def per(self, *args, **kwds):
        """
        Synonym for self.permanent(...).

        EXAMPLES:
            sage: A = MatrixSpace(Integers(8),3)([1,7,3, 1,1,1, 3,4,5])
            sage: A.per()
            6
        """
        return self.permanent(*args, **kwds)

    def characteristic_polynomial(self, *args, **kwds):
        """
        Synonym for self.charpoly(...).
        """
        try:
            return self.charpoly(*args, **kwds)
        except AttributeError:
            raise NotImplementedError

    def minimal_polynomial(self, *args, **kwds):
        """
        Synonym for self.charpoly(...).
        """
        try:
            return self.minpoly(*args, **kwds)
        except AttributeError:
            raise NotImplementedError

    ###########################

    def rescale_row(self, i, s):
        """
        Replace i-th row of self by s times i-th row of self.
        """
        if s == 1: return
        self._require_mutable()
        for j in xrange(self.ncols()):
            self[i,j] *= s

    def row(self, i):
        """
        Return the i-th row of this matrix.
        """
        V = sage.modules.free_module.FreeModule(self.base_ring(), self.ncols())
        return V([self[i,j] for j in xrange(self.ncols())])

    def rows(self):
        """
        Return a list of the rows of this matrix.
        """
        return [self.row(i) for i in xrange(self.nrows())]

    def linear_combination_of_rows(self, v):
        R = self.rows()
        return sum([v[i]*R[i] for i in xrange(len(v))], 0)

    def linear_combination_of_columns(self, v):
        C = self.columns()
        return sum([v[i]*C[i] for i in xrange(len(v))], 0)

    def set_row_to_multiple_of_row(self, i, j, s):
        """
        Set row i equal to s times row j.
        """
        self[i] = s*self[j]

    def _set_row_to_negative_of_row_of_A_using_subset_of_columns(self, i, A, r, cols):
        # this will be insanely slow for a generic matrix, and should be
        # overloaded for specific matrix classes!
        # the ints cols are assumed sorted.
        # this function exists just because it is useful for modular symbols presentations.
        l = 0
        rows = A.sparse_rows()
        v = rows[r]
        for k in cols:
            if k in v.keys():
                self[i,l] = -v[k]
            l += 1

    def sparse_columns(self):
        try:
            return self.__sparse_columns
        except AttributeError:
            C = [{} for _ in xrange(len(self.ncols()))]
            for i, j in self.nonzero_positions():
                C[i][j] = self.get(i,j)
            if self.is_immutable():
                self.__sparse_columns = C
            return C

    def sparse_rows(self):
        try:
            return self.__sparse_rows
        except AttributeError:
            R = [{} for _ in xrange(self.ncols())]
            for i, j in self.nonzero_positions():
                R[i][j] = self.get((i,j))
            if self.is_immutable():
                self.__sparse_rows = R
            return R

    def stack(self, other):
        """
        Return the augmented matrix self on top of other:
           [ self  ]
           [ other ]

        EXAMPLES:
            sage: M = Matrix(RationalField(), 2, 3, range(6))
            sage: N = Matrix(RationalField(), 1, 3, [10,11,12])
            sage: M.stack(N)
            [ 0  1  2]
            [ 3  4  5]
            [10 11 12]
        """
        if not isinstance(other, Matrix):
            raise TypeError, "other must be a matrix"
        if self.base_ring() != other.base_ring():
            raise TypeError, "base rings must be the same"
        if self.ncols() != other.ncols():
            raise TypeError, "number of columns must be the same"
        Z = self.new_matrix(nrows = self.nrows() + other.nrows())
        for r in xrange(self.nrows()):
            for c in xrange(self.ncols()):
                Z[r,c] = self[r,c]
        nr = self.nrows()
        for r in xrange(other.nrows()):
            for c in xrange(other.ncols()):
                Z[r+nr,c] = other[r,c]
        return Z

    def matrix_from_columns(self, columns):
        """
        Return the matrix constructed from self using columns with
        indices in the columns list.


        EXAMPLES:
            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_columns([2,1])
            [2 1]
            [5 4]
            [0 7]
        """
        if not isinstance(columns, (list, tuple)):
            raise TypeError, "columns (=%s) must be a list of integers"%columns
        A = self.new_matrix(ncols = len(columns))
        k = 0
        for i in columns:
            i = int(i)
            if i < 0 or i >= self.ncols():
                raise IndexError, "column %s out of range"%i
            for r in xrange(self.nrows()):
                A.set((r,k), self.get((r,i)))
            k += 1
        return A

    def matrix_from_rows(self, rows):
        """
        Return the matrix constructed from self using rows with indices
        in the rows list.

        EXAMPLES:
            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_rows([2,1])
            [6 7 0]
            [3 4 5]
        """
        if not isinstance(rows, (list, tuple)):
            raise TypeError, "rows must be a list of integers"
        A = self.new_matrix(nrows = len(rows))
        k = 0
        for i in rows:
            i = int(i)
            if i < 0 or i >= self.nrows():
                raise IndexError, "row %s out of range"%i
            for c in xrange(self.ncols()):
                A.set((k,c), self.get((i,c)))
            k += 1
        return A

    def matrix_from_rows_and_columns(self, rows, columns):
        """
        Return the matrix constructed from self from the given
        rows and columns.

        EXAMPLES:
            sage: M = MatrixSpace(Integers(8),3,3)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_rows_and_columns([1], [0,2])
            [3 5]
            sage: A.matrix_from_rows_and_columns([1,2], [1,2])
            [4 5]
            [7 0]

        Note that row and column indices can be reordered or repeated:
            sage: A.matrix_from_rows_and_columns([2,1], [2,1])
            [0 7]
            [5 4]

        For example here we take from row 1 columns 2 then 0 twice,
        and do this 3 times.
            sage: A.matrix_from_rows_and_columns([1,1,1],[2,0,0])
            [5 3 3]
            [5 3 3]
            [5 3 3]

        AUTHOR:
            -- Jaap Spies (2006-02-18)
        """
        if not isinstance(rows, list):
            raise TypeError, "rows must be a list of integers"
        if not isinstance(columns, list):
            raise TypeError, "columns must be a list of integers"
        A = self.new_matrix(nrows = len(rows), ncols = len(columns))
        r = 0
        c = len(columns)
        columns = [int(j) for j in columns if j >= 0 and j < self.ncols()]
        if c != len(columns):
            raise IndexError, "column index out of range"
        for i in rows:
            i = int(i)
            if i < 0 or i >= self.nrows():
                raise IndexError, "row %s out of range"%i
            k = 0
            for j in columns:
                A.set((r,k), self.get((i,j)))
                k += 1
            r += 1
        return A

    def swap_columns(self, c1, c2):
        """
        Swap columns c1 and c2 of self.
        """
        self._require_mutable()
        nc = self.ncols()
        if c1 < 0 or c1 >= nc:
            raise IndexError, "c1 must satisfy 0<=c1<%s"%nc
        if c2 < 0 or c2 >= nc:
            raise IndexError, "c2 must satisfy 0<=c2<%s"%nc
        if c1 == c2: return
        for r in xrange(self.nrows()):
            a, b = self.get((r,c2)), self.get((r,c1))
            self.set((r,c1), a)
            self.set((r,c2), b)

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
            a, b = self.get((r2,c)), self.get((r1,c))
            self.set((r1,c), a)
            self.set((r2,c), b)

    def trace(self):
        """
        Return the trace of self, which is the sum of the
        diagonal entries of self.
        INPUT:
            self -- a square matrix
        OUTPUT:
            element of the base ring of self
        """
        if self.nrows() != self.ncols():
            raise ArithmeticError, "matrix must be square"
        R = self.base_ring()
        return misc.add([self[i,i] for i in xrange(self.nrows())],R(0))

    def vector_matrix_multiply(self, v):
        """
        Returns the vector times matrix product.

        INPUT:
             v -- a free module element.

        OUTPUT:
            The the vector times matrix product v*A.

        EXAMPLES:
            sage: MS = MatrixSpace(RationalField(), 2,2)
            sage: B = MS.matrix([1,2,1,2])
            sage: V = VectorSpace(RationalField(), 2)
            sage: v = V([1,2])
            sage: B.vector_matrix_multiply(v)     # computes v*B
            (3, 6)
            sage: Bt = B.transpose()
            sage: Bt.vector_matrix_multiply(v)    # computes B*v
            (5, 5)
        """
        M = sage.modules.free_module.FreeModule(self.base_ring(), self.ncols())
        if not isinstance(v, sage.modules.free_module_element.FreeModuleElement):
            v = M(v)
        if self.nrows() != v.degree():
            raise ArithmeticError, "number of rows of matrix must equal degree of vector"
        s = M(0)
        for i in xrange(self.nrows()):
            if v[i] != 0:
                s += v[i]*self.row(i)
        return s

    def iterates(self, v, n):
        """
        Let $A$ be this matrix and $v$ be a free module element.
        Return a vector whose rows are the entries of the following
        vectors:
        $$
           v, v A, v A^2, \ldots, v A^{n-1}.
        $$

        INPUT:
            v -- free module element
            n -- nonnegative integer

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(), 2)([1,1,3,5]); A
            [1 1]
            [3 5]
            sage: v = FreeModule(IntegerRing(), 2)([1,0])
            sage: A.iterates(v,0)
            []
            sage: A.iterates(v,5)
            [  1   0]
            [  1   1]
            [  4   6]
            [ 22  34]
            [124 192]
        """
        n = int(n)
        if n >= 2 and self.nrows() != self.ncols():
            raise ArithmeticError, "matrix must be square if n >= 2."
        if n == 0:
            return self.matrix_space(n, self.ncols())(0)
        M = sage.modules.free_module.FreeModule(self.base_ring(), self.ncols())
        if not isinstance(v, sage.modules.free_module_element.FreeModuleElement):
            v = M(v)
        X = [v]
        for _ in range(n-1):
            X.append(X[len(X)-1]*self)
        MS = self.matrix_space(n, self.ncols())
        return MS(X)



    ###################################################
    ## Arithmetic
    ###################################################
    def __abs__(self):
        return self.determinant()

    def __is_compatible(self, other):
        return isinstance(other, Matrix) and self.parent() == other.parent()

    def _add_(self, right):
        A = self.new_matrix()
        for i in xrange(self.nrows()):
            for j in xrange(self.ncols()):
                A[i,j] = self[i,j] + right[i,j]
        return A

    def __cmp__(self, right):
        if not isinstance(right, Matrix) or right.parent() != self.parent():
            return coerce.cmp(self, right)
        return cmp(self._entries(), right._entries())

    def _div_(self, right):
        return self*(~right)

    def __mod__(self, p):
        raise NotImplementedError

    def _scalar_multiply(self, x):
        if not x in self.base_ring():
            x = self.base_ring()(x)
        M = self.new_matrix()
        nc = self.ncols()
        for i in xrange(self.nrows()):
            for j in xrange(self.ncols()):
                M[i,j] = x * self[i,j]
        return M

    def _right_scalar_multiply(self, x):
        if not x in self.base_ring():
            x = self.base_ring()(x)
        M = self.new_matrix()
        nc = self.ncols()
        for i in xrange(self.nrows()):
            for j in xrange(self.ncols()):
                M[i,j] = self[i,j] * x
        return M

    def __mul__(self, right):
        if not isinstance(right, Matrix):
            return self._right_scalar_multiply(right)
        if self.base_ring() != right.base_ring():
            raise ArithmeticError, "base rings must be equal"
        if self.ncols() != right.nrows():
            raise ArithmeticError, "number of columns of self (=%s) must equal number of rows of right (=%s)."%(
                self.ncols(), right.nrows())
        nr = self.nrows()
        snc = self.ncols()
        nc = right.ncols()
        M = self.new_matrix(nr, nc)
        if snc == 0: return M
        zero = self.base_ring()(0)
        for i in xrange(nr):
            for j in xrange(nc):
                M[i,j] = misc.add([self[i,k]*right[k,j] for k in range(snc)],zero)
        return M

    def __neg__(self):
        return self.__rmul__(self.base_ring()(-1))

    def __pos__(self):
        return self

    def __pow__(self, n):
        """
        EXAMPLES:
            sage: MS = MatrixSpace(RationalField(), 3, 3)
            sage: A = MS([0, 0, 1, 1, 0, '-2/11', 0, 1, '-3/11'])
            sage: A * A**(-1) == 1
            True
            sage: A**4
            [      -3/11     -13/121   1436/1331]
            [    127/121   -337/1331 -4445/14641]
            [    -13/121   1436/1331 -8015/14641]
        """
        if self.nrows() != self.ncols():
            raise ArithmeticError, "self must be square"
        #if not isinstance(n, int):
        #    raise TypeError, "n must be an int"
        n = int(n)
        if n == 0:
            return self.parent()(1)
        if n < 0:
            return (~self)**(-n)
        ans = self.parent()(1)
        apow = self
        while n != 0:
            if n%2 != 0:
                ans = ans * apow
            n /= 2
            if n == 0:  # to not waste time doing an extra multiplication/increment
                break
            apow = apow * apow
        return ans

    def __rmul__(self, left):
        A = self.new_matrix()
        for i, j in self.nonzero_positions():
            A[i,j] = left*self[i,j]
        return A

    def _sub_(self, right):
        A = self.new_matrix()
        for i in xrange(self.nrows()):
            for j in xrange(self.ncols()):
                A[i,j] = self[i,j] - right[i,j]
        return A


    def adjoint(self):
        """
        Returns the adjoint matrix of self (matrix of cofactors).

        ALGORITHM: Use PARI
        """
        return self.parent()(self._pari_().matadjoint().python())


    def lllgram(self):
        """
        LLL reduction of the lattice whose gram matrix is self.

        Returns the unimodular transformation matrix!

        ALGORITHM: Use PARI
        """
        if not self.is_square():
            raise ArithmeticError, "matrix must be square"
        Z = integer_ring.IntegerRing()
        n = self.nrows()
        MS = matrix_space.MatrixSpace(Z,n,n)
        if self.base_ring() == Z:
            U = MS(self._pari_().lllgramint().python())
        else:
            U = MS(self._pari_().lllgram().python())
        if U.det() == -1:
            for i in range(n):
                U[i,n-1] = - U[i,n-1]
        return U


#############################################
## Generic matrices over an integral domain
#############################################

class Matrix_domain(Matrix):
    def __init__(self, parent):
        Matrix.__init__(self, parent)

    def charpoly(self, *args, **kwds):
        r"""
        Return the characteristic polynomial of self, as a polynomial
        over the fraction field of the base ring.

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
            Univariate Polynomial Ring in x over Rational Field

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
            Univariate Polynomial Ring in x over Fraction Field of Univariate Polynomial Ring in a over Integer Ring
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
            x^2 + (-1*x_1^2 - x_0)*x + x_0*x_1^2 - x_0^2*x_1

        It's a little difficult to distinguish the variables.  To fix this,
        we rename the indeterminate $Z$:
            sage: f.parent().assign_names("Z")
            sage: f
            Z^2 + (-1*x_1^2 - x_0)*Z + x_0*x_1^2 - x_0^2*x_1

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
        return f

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
            x_0*x_1^2 - x_0^2*x_1
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

        SEE ALSO: elementary_divisors
        """
        v = self._pari_().matsnf(1).python()
        M = self.parent()
        return M(v[2]), M(v[0]), M(v[1]),

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
            sage: K = FractionField(PolynomialRing(IntegerRing()))
            sage: x = K.gen()
            sage: A = MatrixSpace(K,2)([1/x, 2/(x+1), 1, 5/(x**3)])
            sage: A.denominator()
            x^4 + x^3
        """
        if self.nrows() == 0 or self.ncols() == 0:
            return 1
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
        if is_FiniteField(R) and R.is_prime_field() and R.characteristic() < 46340:
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
                        A.rescale_row(r, ~A.get((r,c)))
                        # Swap
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
            [             0 -zeta_12^2 - 1]
            [             0  zeta_12^2 - 1]
            sage: M.kernel()
            Vector space of degree 4 and dimension 2 over Cyclotomic Field of order 12 and degree 4
            Basis matrix:
            [               0                1                0     -2*zeta_12^2]
            [               0                0                1 -2*zeta_12^2 + 1]

        A nontrivial kernel over a complicated base field.
            sage: K = FractionField(MPolynomialRing(QQ, 2))
            sage: M = MatrixSpace(K, 2)([[K.1, K.0], [K.1, K.0]])
            sage: M
            [x_1 x_0]
            [x_1 x_0]
            sage: M.kernel()
            Vector space of degree 2 and dimension 1 over Fraction Field of Polynomial Ring in x_0, x_1 over Rational Field
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
        M = sage.matrix.matrix_space.MatrixSpace(R, self.nrows(), self.ncols())(A)
        return M.kernel()


    def is_invertible(self):
        """
        Return True if this matrix is invertible.
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
        if not isinstance(right, Matrix):
            return self._right_scalar_multiply(right)
        if not isinstance(right, Matrix_generic_sparse):
            right = self.matrix_space(right.nrows(), right.ncols())(right)
        if self.base_ring() != right.base_ring():
            raise ArithmeticError, "base rings must be equal"
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
            return 1
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
        Return the submatrix of self of columns col[i] for i in
        the list of columns.
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
        Return the submatrix of self of rows row[i] for i in
        the list of rows.
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
        entries2 = []
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

class Matrix_generic_sparse_domain(Matrix_domain, Matrix_generic_sparse):
    def __init__(self, parent, entries=0,
                       coerce_entries=True,
                       copy=True):
        Matrix_generic_sparse.__init__(self, parent, entries, coerce_entries, copy)


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
## Dense matrices over the rational numbers
#############################################
class Matrix_dense_rational(Matrix_field):
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
        if not isinstance(right, Matrix_dense_rational):
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
class Matrix_sparse_rational(Matrix_field):
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
                X  = X.submatrix_from_rows(range(r))
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



def _combinations(sequence, number):
    """
    Generate all combinations of \code{number} elements from list
    \code{sequence}.

    Based on code from \code{test/test_generators.py}.

    AUTHOR:
        -- Jaap Spies (2006-02-18)
    """
    if number > len(sequence):
	return
    if number == 0:
	yield []
    else:
	first, rest = sequence[0], sequence[1:]
        # first in combination
	for result in _combinations(rest, number-1):
	    result.insert(0, first)
	    yield result
        # first not in combination
	for result in _combinations(rest, number):
	    yield result
