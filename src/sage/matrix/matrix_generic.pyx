r"""
Matrices over an arbitrary ring

AUTHORS:
    -- William Stein
    -- Martin Albrecht: conversion to Pyrex
    -- Jaap Spies: various functions
    -- Gary Zablackis: fixed a sign bug in generic determinant.

Elements of matrix spaces are of class \code{Matrix} (or a class
derived from Matrix).  They can be either sparse or dense, and can be
defined over any base ring.

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

We save and load a matrix:
    sage: A = Matrix(Integers(8),3,range(9))
    sage: loads(dumps(A)) == A
    True


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

EXAMPLES:
    sage: A = Matrix(Integers(8),3,range(9))
    sage: A.determinant()
    0
    sage: A[0,0] = 5
    sage: A.determinant()
    1
    sage: A.set_immutable()
    sage: A[0,0] = 5
    Traceback (most recent call last):
    ...
    ValueError: object is immutable; please change a copy instead.
"""


#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

doc = """
Matrices
"""

import sage.modules.free_module_element
import sage.modules.free_module
import sage.misc.latex
import sage.misc.functional
import sage.structure.coerce
from   sage.structure.sequence import _combinations

from sage.structure.mutability_pyx cimport Mutability

def is_Matrix(x):
    return isinstance(x, Matrix)

cdef class Matrix(sage.structure.element.ModuleElement):
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
            <type 'matrix_generic.Matrix'>
        """
        sage.structure.element.ModuleElement.__init__(self, parent)
        self.__nrows = parent.nrows()
        self.__ncols = parent.ncols()
        self._mutability = Mutability(False)

    def _require_mutable(self):
        self._mutability._require_mutable()

    def set_immutable(self):
        self._mutability.set_immutable()

    def is_immutable(self):
        return self._mutability.is_immutable()

    def is_mutable(self):
        return self._mutability.is_mutable()

    def _matrix_(self, R):
        return self.change_ring(R)

    def __reduce__(self):
        return sage.matrix.matrix.__reduce__Matrix_generic, \
                 (self.__class__, self.__dict__,
                  self.parent(), self.__dict, self.__determinant,
                  self.__sparse_columns, self.__sparse_rows,
                  self._mutability)

    def _init(self, attrs, parent, dict, determinant,
              sparse_columns, sparse_rows,
              mutability):
        sage.structure.element.ModuleElement.__init__(self, parent)
        self.__dict__ = attrs
        self.__dict = dict
        self.__determinant = determinant
        self.__sparse_columns = sparse_columns
        self.__sparse_rows = sparse_rows
        self.__nrows = parent.nrows()
        self.__ncols = parent.ncols()
        self._mutability = mutability

    def __repr__(self):
        nr = int(self.nrows())
        nc = int(self.ncols())
        if nr == 0 or nc == 0:
            return "[]"
        # compute column widths
        S = []
        for x in self.list():
            S.append(str(x))

        tmp = []
        for x in S:
            tmp.append(len(x))

        width = max(tmp)
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
                s = s + entry + sep
            rows.append(s)

        # Remove leading spaces
##        n = width-m
##        if n > 0:
##            for r in xrange(nr):
##                rows[r] = rows[r][n:]

        # Put brackets around in a single string
        tmp = []
        for row in rows:
            tmp.append("[%s]"%row)
        s = "\n".join(tmp)
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
            v = []
            for j in xrange(nc):
                if self[i,j] != 0:
                    v.append((j, self[i,j]))
            for j in xrange(len(v)):
                s  = s + "%s*%s_{%s}"%(v[j][1], variable, v[j][0])
                if j == 0:
                    s = s + "& + "
                elif j < len(v) - 1:
                    s =  s + " + "
                else:
                    s =  s + "\\\\\n"
        s = s + "\n\\end{align*}"
        return s

    def _latex_(self):
        nr = self.nrows()
        nc = self.ncols()
        if nr == 0 or nc == 0:
            return "()"
        # compute column widths
        S = []
        for x in self.list():
            S.append(str(x))

        tmp = []
        for x in S:
            tmp.append(len(x))

        width = max(tmp)

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
                entry = sage.misc.latex.latex(S[r*nc+c])
                if c == 0:
                    m = max(m, len(entry))
                s = s + entry + sep
            rows.append(s)

        # Put brackets around in a single string
        tmp = []
        for row in rows:
            tmp.append("%s"%row)
        s = "\\\\\n".join(tmp)
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
        v = []
        for i in range(nr):
            tmp = []
            for j in range(nc):
                tmp.append(w[i*nc + j]._pari_init_())
            v.append( ','.join(tmp))
        return 'Mat([%s])'%(';'.join(v))

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
        v = []
        for i in xrange(self.nrows()):
            tmp = []
            for j in xrange(self.ncols()):
                tmp.append(self.get((i,j))._gap_init_())
            v.append( '[%s]'%(','.join(tmp)) )
        return '[%s]'%(','.join(v))

    ###################################################
    ## Coercion to Mathematica
    ###################################################
    #def _mathematica_init_(self):
    #    """
    #    EXAMPLES:
    #        sage: A = MatrixSpace(QQ,3)([1,2,3,4/3,5/3,6/4,7,8,9])
    #        sage: g = mathematica(A); g
    #    """
    #    v = ['{%s}'%(','.join([self.get((i,j))._mathematica_init_() for j in xrange(self.ncols())])) for
    #         i in xrange(self.nrows())]
    #    return '{%s}'%(','.join(v))


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
        v = []
        for x in self.list():
            v.append(x._magma_init_())
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
                   coerce_entries=True, copy=True, sparse=None,
                   clear = True, zero=True):
        """
        Create a matrix in the parent of this space with the given
        number of rows, columns, etc.

        WARNING: This function called with no arguments returns the 0
        matrix by default, not the matrix self.
        """
        # TODO: propigate the zero flag
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
        ans = []
        for i in range(n):
            tmp = []
            for j in range(n):
                tmp.append(self.get( (i,j) )*vars[j])
            ans.append( sum(tmp) )
        return f(tuple(ans))

    def add_multiple_of_row(self, i, j, s):
        """
        Add s times row j to row i.
        """
        self._require_mutable()
        s = self.base_ring()(s)
        for c in xrange(self.ncols()):
            self[i,c] = self[i,c] + s*self[j,c]

    def add_multiple_of_column(self, i, j, s):
        """
        Add s times column j to column i.
        """
        self._require_mutable()
        s = self.base_ring()(s)
        for r in xrange(self.nrows()):
            self[r,i] = self[r,i] + s*self[r,j]

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

    def change_ring(self, ring, copy=False):
        """
        Return the matrix obtained by coercing the entries of this
        matrix into the given ring.
        """
        if ring == self.base_ring():
            if copy:
                return self.copy()
            else:
                return self
        M = sage.matrix.matrix_space.MatrixSpace(ring, self.__nrows, self.__ncols)
        return M(self.list())

    def column(self, i):
        """
        Return the i-th column of this matrix.
        """
        V = sage.modules.free_module.FreeModule(self.base_ring(), self.ncols())
        tmp = []
        for j in xrange(self.nrows()):
            tmp.append(self[j,i])
        return V(tmp)

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
        tmp = []
        for i in xrange(self.ncols()):
            tmp.append(self[n, i])

        return V(tmp)

    def _dict(self):
        if self.__dict != None:
            return self.__dict
        else:
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
        tmp = []
        for i in xrange(self.nrows()):
            for j in xrange(self.ncols()):
                tmp.append(self[i,j])
        return tmp


    def nonpivots(self):
        """
        Return the list of i such that the i-th column of self
        is NOT a pivot column of the reduced row echelon form
        of self.

        OUTPUT:
            list -- sorted list of integers
        """
        X = set(self.pivots())
        tmp = []
        for j in xrange(self.ncols()):
            if not (j in X):
                tmp.append(j)
        return tmp

    def nonzero_positions(self):
        """
        Returns the set of pairs (i,j) such that self[i,j] != 0.
        """
        z = self.base_ring()(0)
        tmp = set()
        for i in xrange(self.nrows()):
            for j in xrange(self.ncols()):
                if self[i,j] != z:
                    tmp.add((i,j))
        return tmp

    def nonzero_positions_in_column(self, i):
        """
        Return the integers j such that self[j,i] is nonzero, i.e.,
        such that the j-th position of the i-th column is nonzero.
        """
        z = self.base_ring()(0)
        tmp = set()
        for j in xrange(self.nrows()):
            if self[j,i] != z:
                tmp.add(j)
        return tmp

    def nonzero_positions_in_row(self, i):
        """
        Return the integers j such that self[i,j] is nonzero, i.e.,
        such that the j-th position of the i-th row is nonzero.
        """
        z = self.base_ring()(0)
        tmp = set()
        for j in xrange(self.ncols()):
            if self[i,j] != z:
                tmp.add(j)
        return tmp

    def prod_of_row_sums(self, cols):
        r"""
        Calculate the product of all row sums of a submatrix of $A$ for a
        list of selected columns \code{cols}.

        AUTHOR:
            -- Jaap Spies (2006-02-18)
        """
        pr = 1
        for row in xrange(self.nrows()):
            tmp = []
            for c in cols:
                tmp.append(self[row, c])
            pr = pr * sum(tmp)
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
            Searching Sloane's online database...
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
            x0*x1^2 + x0^2*x1


        AUTHOR:
            -- Jaap Spies (2006-02-16)
                Copyright (C) 2006 Jaap Spies <j.spies@hccnet.nl>
                Copyright (C) 2006 William Stein <wstein@gmail.com>
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
            tmp = []
            for cols in lst:
                tmp.append(self.prod_of_row_sums(cols))
            s = sum(tmp)
            perm = perm + (-1)**(m-r) * sage.rings.arith.binomial(n-r, m-r) * s
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
                pm = pm + self.matrix_from_rows_and_columns(rows, cols).permanent()
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
        tmp = []
        for k in range(m+1):
            tmp.append(self.permanental_minor(k))
        return tmp


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
            sage: A.determinant() is A.determinant()
            False
            sage: A.set_immutable()
            sage: A.determinant() is A.determinant()
            True

        We compute the determinant of the arbitrary 3x3 matrix:
            sage: R = PolynomialRing(QQ,9,'x')
            sage: A = matrix(R,3,R.gens())
            sage: A
            [x0 x1 x2]
            [x3 x4 x5]
            [x6 x7 x8]
            sage: A.determinant()
            -1*x2*x4*x6 + x2*x3*x7 + x1*x5*x6 - x1*x3*x8 - x0*x5*x7 + x0*x4*x8
        """
        if self.__determinant is not None:
            return self.__determinant
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
        sgn = R(-1)
        for i in range(n):
            v = range(n)
            del v[i]
            B = A.matrix_from_columns(v)
            d = d + s*self.get((0,i))*B.determinant()
            s = s*sgn
        if self.is_immutable():
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
        return self.charpoly(*args, **kwds)

    def charpoly(self, *args, **kwds):
        raise NotImplementedError

    def minimal_polynomial(self, *args, **kwds):
        """
        Synonym for self.charpoly(...).
        """
        return self.minpoly(*args, **kwds)

    def minpoly(self, *args, **kwds):
        raise NotImplementedError

    ###########################

    def rescale_row(self, i, s):
        """
        Replace i-th row of self by s times i-th row of self.
        """
        if s == 1: return
        self._require_mutable()
        for j in xrange(self.ncols()):
            self[i,j] = self[i,j] * s

    def row(self, i):
        """
        Return the i-th row of this matrix.
        """
        V = sage.modules.free_module.FreeModule(self.base_ring(), self.ncols())
        tmp = []
        for j in xrange(self.ncols()):
            tmp.append(self[i,j])
        return V(tmp)

    def rows(self):
        """
        Return a list of the rows of this matrix.
        """
        tmp = []
        for i in xrange(self.nrows()):
            tmp.append(self.row(i))
        return tmp

    def linear_combination_of_rows(self, v):
        R = self.rows()
        tmp = []
        for i in xrange(len(v)):
            tmp.append(v[i]*R[i])
        return sum(tmp, 0)

    def linear_combination_of_columns(self, v):
        C = self.columns()
        tmp = []
        for i in xrange(len(v)):
            tmp.append(v[i]*C[i])
        return sum(tmp, 0)

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
            l = l + 1

    def sparse_columns(self):
        if self.__sparse_columns!=None:
            return self.__sparse_columns
        else:
            C =[]
            for _ in xrange(self.ncols()):
                C.append({})
            for i, j in self.nonzero_positions():
                C[i][j] = self.get(i,j)
            if self.is_immutable():
                self.__sparse_columns = C
            return C

    def sparse_rows(self):
        if self.__sparse_rows != None:
            return self.__sparse_rows
        else:
            R = []
            for _ in xrange(self.ncols()):
                R.append({})
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
            k = k + 1
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
            k = k + 1
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
        tmp = []
        for j in columns:
            if j >= 0 and j < self.ncols():
                tmp.append(int(j))
        columns = tmp
        if c != len(columns):
            raise IndexError, "column index out of range"
        for i in rows:
            i = int(i)
            if i < 0 or i >= self.nrows():
                raise IndexError, "row %s out of range"%i
            k = 0
            for j in columns:
                A.set((r,k), self.get((i,j)))
                k = k + 1
            r = r + 1
        return A

    def swap_columns(self, c1, c2):
        """
        Swap columns c1 and c2 of self.

        EXAMPLES:
            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A
            [  1   9  -7]
            [4/5   4   3]
            [  6   4   3]
            sage: A.swap_columns(1,2); A
            [  1  -7   9]
            [4/5   3   4]
            [  6   3   4]
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

        EXAMPLES:
            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A
            [  1   9  -7]
            [4/5   4   3]
            [  6   4   3]
            sage: A.swap_rows(0,2); A
            [  6   4   3]
            [4/5   4   3]
            [  1   9  -7]
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
        tmp =  []
        for i in xrange(self.nrows()):
            tmp.append(self[i,i])

        return sum(tmp,R(0))

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
                s = s + v[i]*self.row(i)
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
        if not isinstance(self,Matrix):
            self,right = right,self

        A = self.new_matrix()
        for i in xrange(self.nrows()):
            for j in xrange(self.ncols()):
                A[i,j] = self[i,j] + right[i,j]
        return A

    def __cmp__(self, right):
        if not isinstance(right, Matrix) or right.parent() != self.parent():
            return sage.structure.coerce.cmp(self, right)
        return cmp(self._entries(), right._entries())

##     def __richcmp__(self,right,op):
##         res = 0
##         if not isinstance(right, Matrix) or right.parent() != self.parent():
##             res = sage.ext.coerce.cmp(self, right)
##         else:
##             res = cmp(self._entries(), right._entries())

##         if op == 0:  #<
##             return bool(res  < 0)
##         if op == 2: #==
##             return bool(res == 0)
##         if op == 4: #>
##             return bool(res  > 0)
##         if op == 1: #<=
##             return bool(res <= 0)
##         if op == 3: #!=
##             return bool(res != 0)
##         if op == 5: #>=
##             return bool(res >= 0)

    def __nonzero__(self):
        """
        EXAMPLES:
            sage: M = Matrix(ZZ, 2, 2, [0, 0, 0, 0])
            sage: bool(M)
            False
            sage: M = Matrix(ZZ, 2, 2, [1, 2, 3, 5])
            sage: bool(M)
            True
        """
        return sage.structure.coerce.cmp(self, 0)

    def _div_(self, right):
        if not isinstance(self,Matrix):
            self,right = right,self

        return self*(~right)

    def __mod__(self, p):
        r"""
        Return matrix mod $p$, returning again a matrix over the same
        base ring.

        \note{Use \code{A.Mod(p)} to obtain a matrix over the residue
        class ring modulo $p$.}

        EXAMPLES:
            sage: M = Matrix(ZZ, 2, 2, [5, 9, 13, 15])
            sage: M % 7
            [5 2]
            [6 1]
            sage: parent(M % 7)
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        """
        M = self.new_matrix()
        nc = self.ncols()
        for i in xrange(self.nrows()):
            for j in xrange(nc):
                M.set((i,j), self[i,j] % p)
        return M

    def Mod(self, p):
        """
        EXAMPLES:
            sage: M = Matrix(ZZ, 2, 2, [5, 9, 13, 15])
            sage: M.Mod(7)
            [5 2]
            [6 1]
            sage: parent(M.Mod(7))
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7
        """
        return self.change_ring(self.base_ring().quotient_ring(p))

    def lift(self):
        """
        EXAMPLES:
            sage: M = Matrix(ZZ/7, 2, 2, [5, 9, 13, 15]) ; M
            [5 2]
            [6 1]
            sage: M.lift()
            [5 2]
            [6 1]
            sage: parent(M.lift())
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        """
        return self.change_ring(self.base_ring().cover_ring())

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
        if not isinstance(self,Matrix):
            self,right = right,self

        if isinstance(right, sage.modules.free_module_element.FreeModuleElement):
            return self.transpose().vector_matrix_multiply(right)
        if not isinstance(right, Matrix):
            return self._right_scalar_multiply(right)
        if self.base_ring() != right.base_ring():
            try:
                self, right = sage.structure.coerce.canonical_base_coercion(self, right)
            except TypeError:
                raise ArithmeticError, "base rings must be compatible"
        if self.ncols() != right.nrows():
            raise ArithmeticError, "number of columns of self (=%s) must equal number of rows of right (=%s)."%(
                self.ncols(), right.nrows())
        try:
            return self._mul_(right)
        except (TypeError, AttributeError):
            pass
        nr = self.nrows()
        snc = self.ncols()
        nc = right.ncols()
        M = self.new_matrix(nr, nc)
        if snc == 0: return M
        zero = self.base_ring()(0)
        for i in xrange(nr):
            for j in xrange(nc):
                tmp = []
                for k in range(snc):
                    tmp.append(self[i,k]*right[k,j])
                M[i,j] = sum(tmp,zero)
        return M

    def __neg__(self):
        return self.__rmul__(self.base_ring()(-1))

    def __pos__(self):
        return self

    def __pow__(self, n, ignored):
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
            n = n/2
            if n == 0:  # to not waste time doing an extra multiplication/increment
                break
            apow = apow * apow
        return ans

    def __rmul__(self, left):
        if not isinstance(self,Matrix):
            self,left = left,self

        A = self.new_matrix()
        for i, j in self.nonzero_positions():
            A[i,j] = left*self[i,j]
        return A

    def _sub_(self, right):
        if not isinstance(self,Matrix):
            self,right = right,self

        A = self.new_matrix()
        for i in xrange(self.nrows()):
            for j in xrange(self.ncols()):
                A[i,j] = self[i,j] - right[i,j]
        return A


    def adjoint(self):
        """
        Returns the adjoint matrix of self (matrix of cofactors).

        INPUT:
            M -- a square matrix

        OUTPUT:
            N -- the adjoint matrix, such that

                N * M = M * N = M.parent(M.det())

        ALGORITHM:
            Use PARI

        EXAMPLES:
            sage: M = Matrix(ZZ,2,2,[5,2,3,4]) ; M
            [5 2]
            [3 4]
            sage: N = M.adjoint() ; N
            [ 4 -2]
            [-3  5]
            sage: M * N
            [14  0]
            [ 0 14]
            sage: N * M
            [14  0]
            [ 0 14]
            sage: M = Matrix(QQ,2,2,[5/3,2/56,33/13,41/10]) ; M
            [  5/3  1/28]
            [33/13 41/10]
            sage: N = M.adjoint() ; N
            [ 41/10  -1/28]
            [-33/13    5/3]
            sage: M * N
            [7363/1092         0]
            [        0 7363/1092]

        BUGS:
            only implemented for matrices over ZZ or QQ
            PARI can deal with more general base rings
        """
        if not self.is_square():
            raise ArithmeticError, "matrix must be square"
        return self._adjoint()

    def _adjoint(self):
        raise NotImplementedError

    def lllgram(self):
        """
        LLL reduction of the lattice whose gram matrix is self.

        INPUT:
            M -- gram matrix of a definite quadratic form

        OUTPUT:
            U -- unimodular transformation matrix such that

                U.transpose() * M * U

            is LLL-reduced

        ALGORITHM:
            Use PARI

        EXAMPLES:
            sage: M = Matrix(ZZ, 2, 2, [-5,3,3,-2]) ; M
            [-5  3]
            [ 3 -2]
            sage: U = M.lllgram() ; U
            [1 1]
            [1 2]
            sage: U.transpose() * M * U
            [-1  0]
            [ 0 -1]
            sage: M = Matrix(QQ,2,2,[269468, -199019/2, -199019/2, 36747]) ; M
            [   269468 -199019/2]
            [-199019/2     36747]
            sage: U = M.lllgram() ; U
            [-113  -24]
            [-306  -65]
            sage: U.transpose() * M * U
            [  2 1/2]
            [1/2   3]

        Semidefinite and indefinite forms raise a ValueError:

            sage: Matrix(ZZ,2,2,[2,6,6,3]).lllgram()
            Traceback (most recent call last):
            ...
            ValueError: not a definite matrix
            sage: Matrix(ZZ,2,2,[1,0,0,-1]).lllgram()
            Traceback (most recent call last):
            ...
            ValueError: not a definite matrix

        BUGS:
            should work for semidefinite forms (PARI is ok)
        """
        if not self.is_square():
            raise ArithmeticError, "matrix must be square"
        return self._lllgram()

    def _lllgram(self):
        raise NotImplementedError

    def _mul_submatrices(self,
                         Matrix other,
                         int self_r, int self_c, int self_nrows, int self_ncols,
                         int other_r, int other_c, int other_nrows, int other_ncols):
        """
        For internal use, correct types and dimensions are assumed in this method
        Adds two submatrices

        NOTE: self_ncols = other_nrows

        EXAMPLE:
            sage: A = Matrix(ZZ, 10, 10, range(100))
            sage: B = Matrix(ZZ, 20, 20, range(400))
            sage: M = A._submatrix(1,2,3,4) * B._submatrix(11,12,4,5)
            sage: N = A._mul_submatrices(B, 1, 2, 3, 4, 11, 12, 4, 5)
            sage: M == N
            True

        """
        M = self.new_matrix(self_nrows, other_ncols)
        zero = self.base_ring()(0)
        cdef int i, j, k
        for i from 0 <= i < self_nrows:
            for j from 0 <= j < other_ncols:
                sum = zero
                for k from 0 <= k < self_ncols:
                    sum = sum + self.get((self_r + i,self_c + k))*other.get((other_r + k,other_c + j))
                M.set((i,j), sum)
        return M



    def _add_submatrices(Matrix self, Matrix other,
                         int self_r, int self_c, int nrows, int ncols,
                         int other_r, int other_c):
        """
        For internal use, correct types and dimensions are assumed in this method
        Adds two submatrices

        EXAMPLE:
            sage: A = Matrix(ZZ, 10, 10, range(100))
            sage: B = Matrix(ZZ, 20, 20, range(400))
            sage: M = A._submatrix(1,2,3,4) + B._submatrix(11,12,3,4)
            sage: N = A._add_submatrices(B, 1, 2, 3, 4, 11, 12)
            sage: M == N
            True

        """
        M = self.new_matrix(nrows, ncols)
        cdef int i, j
        for i from 0 <= i < nrows:
            for j from 0 <= j < ncols:
                M.set((i,j), self.get((self_r + i,self_c + j)) + other.get((other_r + i, other_c + j)))
        return M


    def _sub_submatrices(Matrix self, Matrix other,
                         int self_r, int self_c, int nrows, int ncols,
                         int other_r, int other_c):
        """
        For internal use, correct types and dimensions are assumed in this method
        Subtracts two submatrices

        EXAMPLE:
            sage: A = Matrix(ZZ, 10, 10, range(100))
            sage: B = Matrix(ZZ, 20, 20, range(400))
            sage: M = A._submatrix(1,2,3,4) - B._submatrix(11,12,3,4)
            sage: N = A._sub_submatrices(B, 1, 2, 3, 4, 11, 12)
            sage: M == N
            True

        """
        M = self.new_matrix(nrows, ncols)
        cdef int i, j
        for i from 0 <= i < nrows:
            for j from 0 <= j < ncols:
                M.set((i,j), self.get((self_r + i,self_c + j)) - other.get((other_r + i, other_c + j)))
        return M


    def _submatrix(Matrix self,
                  int row, int col, int nrows, int ncols):
        """
        Creates a submatrix from this matrix, no bounds checking

        INPUT:
            M -- matrix
            row -- row of upper left entry of submatrix
            col -- column of upper left entry of submatrix
            nrows -- number of rows in the submatrix
            ncols -- number of columns in the submatrix

        OUTPUT:
            The nrows $\times$ ncols matrix with upper left entry M[row, col]

        EXAMPLE:
            sage: M = Matrix(ZZ, 3, 3, range(9))
            sage: M
            [0 1 2]
            [3 4 5]
            [6 7 8]

            sage: M._submatrix(1,1,2,2)
            [4 5]
            [7 8]

            sage: Matrix(ZZ, 6, 6, range(36))._submatrix(1,2,3,4)
            [ 8  9 10 11]
            [14 15 16 17]
            [20 21 22 23]

        """
        M = self.new_matrix(nrows, ncols)
        cdef int i, j
        for i from 0 <= i < nrows:
            for j from 0 <= j < ncols:
                M.set((i,j), self.get((row + i , col + j)))
        return M


    def _invert_submatrix(Matrix self,
                          int self_r, int self_c, int n):
        M = self._submatrix(self_r, self_c, n, n)
        return ~M.matrix_over_field()



    #####################################################################################
    # Strassen Matrix Multiplication
    #####################################################################################
    def strassen(self, other, int cutoff=64):
        # todo -- add something for nonsquare case

        if self.parent() != other.parent():
            raise IndexError, "Number of columns of self must equal number of rows of other."

        return self._mul_submatrices_strassen(other, 0, 0, self.nrows(), 0, 0, cutoff)


    def _strassen_subdivide(self, int r, int c, int n):
        cdef int m
        m = n/2
        return ((r,   c),
                (r,   c + m)
                (r+m, c)
                (r+m, c+m))


    def _mul_submatrices_strassen(self, other,
                                  int a_r, int a_c, int n,
                                  int b_r, int b_c,
                                  int cutoff):

        if n <= cutoff:
            return self._mul_submatrices(other, a_r, a_c, n, n, b_r, b_c, n, n)

        # 8 Pre-Additions:
        #
        # S0 = A10 + A11,  T0 = B01 - B00
        # S1 = S0  - A00,  T1 = B11 - T0
        # S2 = A00 - A10,  T2 = B11 - B01
        # S3 = A01 - S1,   T3 = T1 - B10

        cdef int m
        m = n/2

        cdef int a00_r, a00_c, a01_r, a01_c, a10_r, a10_c, a11_r, a11_c

        a00_r = a_r; a00_c = a_c
        a01_r = a_r; a01_c = a_c + m
        a10_r = a_r + m; a10_c = a_c
        a11_r = a_r + m; a11_c = a_c + m

        cdef int b00_r, b00_c, b01_r, b01_c, b10_r, b10_c, b11_r, b11_c

        b00_r = b_r; b00_c = b_c
        b01_r = b_r; b01_c = b_c + m
        b10_r = b_r + m; b10_c = b_c
        b11_r = b_r + m; b11_c = b_c + m

        # S0 = A10 + A11,  T0 = B01 - B00
        S0 = self._add_submatrices(self, a10_r, a10_c, m, m, a11_r, a11_c, m, m)
        T0 = other._sub_submatrices(other, b01_r, b01_c, m, m, b00_r, b00_c, m, m)

        # S1 = S0 - A00,   T1 = B11 - T0
        S1 = S0._sub_submatrices(self, 0, 0, m, m, a00_r, a00_c, m, m)
        T1 = other._sub_submatrices(T0, b11_r, b11_c, m, m, 0, 0, m, m)

        # S2 = A00 - A10,  T2 = B11 - B01
        S2 = self._sub_submatrices(self, a00_r, a00_c, m, m, a10_r, a10_c, m, m)
        T2 = other._sub_submatrices(other, b11_r, b11_c, m, m, b01_r, b01_c, m, m)

        # S3 = A01 - S1,   T3 = B10 - T1
        S3 = self._sub_submatrices(S1, a01_r, a01_c, m, m, 0, 0, m, m)
        T3 = other._sub_submatrices(T1, b10_r, b10_c, m, m, 0, 0, m, m)

        # 7. (Potentially) Recursive Multiplications
        # P0 =  A00*B00
        # P1 =  A01*B10
        # P2 =  S0*T0
        # P3 =  S1*T1
        # P4 =  S2*T2
        # P5 =  S3*B11
        # P6 =  A11*T3

        if m <= cutoff:
            # Don't do the multiplications recursively
            # P0 = A00*B00
            P0 = self._mul_submatrices(other, a00_r, a00_c, m, m, b00_r, b00_c, m, m)
            # P1 = A01*B10
            P1 = self._mul_submatrices(other, a01_r, a01_c, m, m, b10_r, b10_c, m, m)
            # P2 =  S0*T0
            P2 = S0 * T0
            # P3 = S1*T1
            P3 = S1 * T1
            # P4 =  S2*T2
            P4 = S2 * T2
            # P5 =  S3*B11
            P5 = S3._mul_submatrices(other, 0, 0, m, m, b11_r, b11_c, m, m)
            # P6 =  A11*T3
            P6 = self._mul_submatrices(T3, a11_r, a11_c, m, m, 0, 0, m, m)

        else:
            # Do the multiplications recursively
            # Don't do the multiplications recursively
            # P0 = A00*B00
            P0 = self._mul_submatrices_strassen(other, a00_r, a00_c, m, b00_r, b00_c, cutoff)
            # P1 = A01*B10
            P1 = self._mul_submatrices_strassen(other, a01_r, a01_c, m, b10_r, b10_c, cutoff)
            # P2 =  S0*T0
            P2 = S0.strassen(T0, cutoff)
            # P3 = S1*T1
            P3 = S1.strassen(T1, cutoff)
            # P4 =  S2*T2
            P4 = S2.strassen(T2, cutoff)
            # P5 =  S3*B11
            P5 = S3._mul_submatrices_strassen(other, 0, 0, m, b11_r, b11_c, cutoff)
            # P6 =  A11*T3
            P6 = self._mul_submatrices_strassen(T3, a11_r, a11_c, m, 0, 0, cutoff)

        # 7 Post Additions:

        U0 = P0 + P1
        U1 = P0 + P3
        U2 = U1 + P4
        U3 = U2 + P6
        U4 = U2 + P2
        U5 = U1 + P2
        U6 = U5 + P5

        return U0.block2_sum(U6, U3, U4)

    # strassen with windows
    # todo: matrix windows only implemented for dense matrices, but no multiple inheritance for cdef class (?)

    def strassen_multiply(self, other, cutoff):
        if self.ncols() != other.nrows():
            raise IndexError, "Number of columns of self must equal number of rows of other."

        output = self.new_matrix(self.nrows(), other.ncols())

        self_window = self.matrix_window()
        other_window = other.matrix_window()
        output_window = output.matrix_window()

        strassen_window_multiply(output_window, self_window, other_window, cutoff)
        return output

    def echelon_strassen(self, cutoff):
        pivots = strassen_echelon(self.matrix_window(), cutoff)
        self._set_pivots(pivots)
        return self


############################# end class ####################


# TODO: need MatrixWindow hierarchy
# def strassen_window_multiply(MatrixWindow C, MatrixWindow A, MatrixWindow B, int cutoff):
def strassen_window_multiply(C, A, B, int cutoff):
    """
    Multiplies the submatrices specified by A and B, places result in
    C.  Assumes that A and B have compatible dimensions to be
    multiplied, and that C is the correct size to receive the product,
    and that they are all defined over the same ring.

    Uses strassen multiplication at high levels and then uses MatrixWindow
    methods at low levels.

    todo: doc cutoff parameter as soon as I work out what it really means

    EXAMPLES:
        The following matrix dimensions are chosen especially to exercise the
        eight possible parity combinations that ocould ccur while subdividing
        the matrix in the strassen recursion. The base case in both cases will
        be a (4x5) matrix times a (5x6) matrix.

        TODO -- the doctests below are currently not
        tested/enabled/working -- enable them when linear algebra
        restructing gets going.

        sage.:P dim1 = 64; dim2 = 83; dim3 = 101
        sage.:P R = MatrixSpace(QQ, dim1, dim2)
        sage.:P S = MatrixSpace(QQ, dim2, dim3)
        sage.:P T = MatrixSpace(QQ, dim1, dim3)


        sage.:P A = R.random_element(range(-30, 30))
        sage.:P B = S.random_element(range(-30, 30))
        sage.:P C = T(0)
        sage.:P D = T(0)

        sage.: A_window = A.matrix_window(0, 0, dim1, dim2)
        sage.: B_window = B.matrix_window(0, 0, dim2, dim3)
        sage.: C_window = C.matrix_window(0, 0, dim1, dim3)
        sage.: D_window = D.matrix_window(0, 0, dim1, dim3)

        sage.: strassen_window_multiply(C, A, B, 2)   # use strassen method
        sage.: D_window.set_to_prod(A, B)             # use naive method
        sage.:P C == D
        True

        sage.:P dim1 = 79; dim2 = 83; dim3 = 101
        sage.:P R = MatrixSpace(QQ, dim1, dim2)
        sage.:P S = MatrixSpace(QQ, dim2, dim3)
        sage.:P T = MatrixSpace(QQ, dim1, dim3)

        sage.:P A = R.random_element(range(30))
        sage.:P B = S.random_element(range(30))
        sage.:P C = T(0)
        sage.:P D = T(0)

        sage.: A_window = A.matrix_window(0, 0, dim1, dim2)
        sage.: B_window = B.matrix_window(0, 0, dim2, dim3)
        sage.: C_window = C.matrix_window(0, 0, dim1, dim3)

        sage.:P strassen_window_multiply(C, A, B, 2)   # use strassen method
        sage.:P D.set_to_prod(A, B)                    # use naive method

        sage.:P C == D
        True

    AUTHOR: David Harvey
    """

    # todo -- I'm not sure how to interpret "cutoff". Should it be...
    # (a) the minimum side length of the matrices (currently implemented below)
    # (b) the maximum side length of the matrices
    # (c) the total number of entries being multiplied
    # (d) something else entirely?

    cdef int A_nrows, A_ncols, B_ncols
    A_nrows = A.nrows()
    A_ncols = A.ncols()        # this should also be the number of rows of B
    B_ncols = B.ncols()

    if (A_nrows <= cutoff) or (A_ncols <= cutoff) or (B_ncols <= cutoff):
        # note: this code is only reached if the TOP level is already beneath
        # the cutoff. In a typical large multiplication, the base case is
        # handled directly (see below).
        C.set_to_prod(A, B)
        return

    # Construct windows for the four quadrants of each matrix.
    # Note that if the side lengths are odd we're ignoring the
    # final row/column for the moment.

    cdef int A_sub_nrows, A_sub_ncols, B_sub_ncols
    A_sub_nrows = A_nrows >> 1
    A_sub_ncols = A_ncols >> 1     # this is also like B_sub_nrows
    B_sub_ncols = B_ncols >> 1

    A00 = A.matrix_window(0,           0,           A_sub_nrows, A_sub_ncols)
    A01 = A.matrix_window(0,           A_sub_ncols, A_sub_nrows, A_sub_ncols)
    A10 = A.matrix_window(A_sub_nrows, 0,           A_sub_nrows, A_sub_ncols)
    A11 = A.matrix_window(A_sub_nrows, A_sub_ncols, A_sub_nrows, A_sub_ncols)
    B00 = B.matrix_window(0,           0,           A_sub_ncols, B_sub_ncols)
    B01 = B.matrix_window(0,           B_sub_ncols, A_sub_ncols, B_sub_ncols)
    B10 = B.matrix_window(A_sub_ncols, 0,           A_sub_ncols, B_sub_ncols)
    B11 = B.matrix_window(A_sub_ncols, B_sub_ncols, A_sub_ncols, B_sub_ncols)

    # Allocate temp space.

    # todo: can I cache the bound A.new_empty_window method?

    S0 = A.new_empty_window(A_sub_nrows, A_sub_ncols, zero=False)
    S1 = A.new_empty_window(A_sub_nrows, A_sub_ncols, zero=False)
    S2 = A.new_empty_window(A_sub_nrows, A_sub_ncols, zero=False)
    S3 = A.new_empty_window(A_sub_nrows, A_sub_ncols, zero=False)

    T0 = A.new_empty_window(A_sub_ncols, B_sub_ncols, zero=False)
    T1 = A.new_empty_window(A_sub_ncols, B_sub_ncols, zero=False)
    T2 = A.new_empty_window(A_sub_ncols, B_sub_ncols, zero=False)
    T3 = A.new_empty_window(A_sub_ncols, B_sub_ncols, zero=False)

    Q0 = A.new_empty_window(A_sub_nrows, B_sub_ncols, zero=False)
    Q1 = A.new_empty_window(A_sub_nrows, B_sub_ncols, zero=False)
    Q2 = A.new_empty_window(A_sub_nrows, B_sub_ncols, zero=False)


    # Preparatory matrix additions/subtractions.

    # todo: we can probably save some memory in these
    # operations by reusing some of the buffers, if we interleave
    # these additions with the multiplications (below)

    # (I believe we can save on one S buffer and one T buffer)

    # S0 = A10 + A11,  T0 = B01 - B00
    # S1 = S0 - A00,   T1 = B11 - T0
    # S2 = A00 - A10,  T2 = B11 - B01
    # S3 = A01 - S1,   T3 = B10 - T1

    S0.set_to_sum(A10, A11)
    S1.set_to_diff(S0, A00)
    S2.set_to_diff(A00, A10)
    S3.set_to_diff(A01, S1)

    T0.set_to_diff(B01, B00)
    T1.set_to_diff(B11, T0)
    T2.set_to_diff(B11, B01)
    T3.set_to_diff(B10, T1)


    # The relations we need now are:

    # P0 = A00*B00
    # P1 = A01*B10
    # P2 =  S0*T0
    # P3 =  S1*T1
    # P4 =  S2*T2
    # P5 =  S3*B11
    # P6 = A11*T3

    # U0 = P0 + P1
    # U1 = P0 + P3
    # U2 = U1 + P4
    # U3 = U2 + P6
    # U4 = U2 + P2
    # U5 = U1 + P2
    # U6 = U5 + P5

    # We place the final answer into the matrix:
    # U0 U6
    # U3 U4

    U0 = C.matrix_window(0,           0,           A_sub_nrows, B_sub_ncols)
    U6 = C.matrix_window(0,           B_sub_ncols, A_sub_nrows, B_sub_ncols)
    U3 = C.matrix_window(A_sub_nrows, 0,           A_sub_nrows, B_sub_ncols)
    U4 = C.matrix_window(A_sub_nrows, B_sub_ncols, A_sub_nrows, B_sub_ncols)

    if (A_sub_nrows <= cutoff) or (A_sub_ncols <= cutoff) or (B_sub_ncols <= cutoff):
        # This is the base case, so we use MatrixWindow methods directly.

        # This next chunk is arranged so that each output cell gets written
        # to exactly once. This is important because the output blocks might
        # be quite fragmented in memory, whereas our temporary buffers
        # (Q0, Q1, Q2) will be quite localised, so we can afford to do a bit
        # of arithmetic in them.

        Q0.set_to_prod(A00, B00)         # now Q0 holds P0
        Q1.set_to_prod(A01, B10)         # now Q1 holds P1
        U0.set_to_sum(Q0, Q1)            # now U0 is correct
        Q0.add_prod(S1, T1)              # now Q0 holds U1
        Q1.set_to_prod(S2, T2)           # now Q1 holds P4
        Q1.add(Q0)                       # now Q1 holds U2
        Q2.set_to_prod(A11, T3)          # now Q2 holds P6
        U3.set_to_sum(Q1, Q2)            # now U3 is correct
        Q2.set_to_prod(S0, T0)           # now Q2 holds P2
        U4.set_to_sum(Q2, Q1)            # now U4 is correct
        Q0.add(Q2)                       # now Q0 holds U5
        Q2.set_to_prod(S3, B11)          # now Q2 holds P5
        U6.set_to_sum(Q0, Q2)            # now U6 is correct

    else:
        # Recurse into sub-products.

        strassen_window_multiply(Q0, A00, B00, cutoff)   # now Q0 holds P0
        strassen_window_multiply(Q1, A01, B10, cutoff)   # now Q1 holds P1
        U0.set_to_sum(Q0, Q1)                            # now U0 is correct
        strassen_window_multiply(Q1, S1, T1, cutoff)     # now Q1 holds P3
        Q0.add(Q1)                                       # now Q0 holds U1
        strassen_window_multiply(Q1, S2, T2, cutoff)     # now Q1 holds P4
        Q1.add(Q0)                                       # now Q1 holds U2
        strassen_window_multiply(Q2, A11, T3, cutoff)    # now Q2 holds P6
        U3.set_to_sum(Q1, Q2)                            # now U3 is correct
        strassen_window_multiply(Q2, S0, T0, cutoff)     # now Q2 holds P2
        U4.set_to_sum(Q2, Q1)                            # now U4 is correct
        Q0.add(Q2)                                       # now Q0 holds U5
        strassen_window_multiply(Q2, S3, B11, cutoff)    # now Q2 holds P5
        U6.set_to_sum(Q0, Q2)                            # now U6 is correct


    # Now deal with the leftover row and/or column (if they exist).

    if B_ncols & 1:
        B_last_col = B.matrix_window(0, B_ncols-1, A_ncols, 1)
        C_last_col = C.matrix_window(0, B_ncols-1, A_nrows, 1)
        C_last_col.set_to_prod(A, B_last_col)

    if A_nrows & 1:
        A_last_row = A.matrix_window(A_nrows-1, 0, 1, A_ncols)
        if B_ncols & 1:
            B_bulk = B.matrix_window(0, 0, A_ncols, B_ncols-1)
            C_last_row = C.matrix_window(A_nrows-1, 0, 1, B_ncols-1)
        else:
            B_bulk = B
            C_last_row = C.matrix_window(A_nrows-1, 0, 1, B_ncols)
        C_last_row.set_to_prod(A_last_row, B_bulk)

    if A_ncols & 1:
        A_last_col = A.matrix_window(0, A_ncols-1, A_sub_nrows << 1, 1)
        B_last_row = B.matrix_window(A_ncols-1, 0, 1, B_sub_ncols << 1)
        C_bulk = C.matrix_window(0, 0, A_sub_nrows << 1, B_sub_ncols << 1)
        C_bulk.add_prod(A_last_col, B_last_row)

# TODO: make cdef
def subtract_strassen_product(result, A, B, int cutoff):
    cutoff = 1000000 # for testing
    if (result.ncols() < cutoff or result.nrows() < cutoff):
        result.subtract_prod(A, B)
    else:
        to_sub = strassen_window_multiply(A, B, cutoff)
        result.subtract(to_sub)


def strassen_echelon(A, cutoff):
    """
    Compute echelon form, in place.
    Internal function, call with M.echelon(alg="block")
    Based on work of Robert Bradshaw and David Harvey at MSRI workshop in 2006.

    INPUT:
        A -- matrix window
        cutoff -- size at which algorithm reverts to naive gaussian elemination and multiplication

    OUTPUT:
        The list of pivot columns

    EXAMPLE:
        sage: import sage.matrix.matrix_rational_dense as m
        sage: parent = MatrixSpace(QQ, 5, 30)
        sage: data = parent.random_element(range(18), prob=.2).list() # test lots of non-pivots
        sage: A = m.Matrix_rational_dense(parent, data)
        sage: T = A.echelon(alg="gauss")
        sage: E = A.echelon_strassen(4)
        sage: A.copy() == T.copy()  # fix when parents are changed
        True

    AUTHORS:
        -- Robert Bradshaw
    """
    # The following notation will be used in the comments below, which should be understood to give
    # the general idea of what's going on, as if there were no inconvenient non-pivot columns.
    # The original matrix is given by [ A B ]
    #                                 [ C D ]
    # For compactness, let A' denote the inverse of A
    # top_left, top_right, bottom_left, and bottom_right loosely correspond to A, B, C, and D respectively,
    # however, the "cut" between the top and bottom rows need not be the same.

    cdef int nrows, ncols
    nrows = A.nrows()
    ncols = A.ncols()

    if (nrows < cutoff or ncols < cutoff):
        return A.echelon_in_place()

    cdef int top_h, bottom_cut, bottom_h, bottom_start, top_cut
    cdef int prev_pivot_count
    cdef int split
    split = nrows / 2

    top = A.matrix_window(0, 0, split, ncols)
    bottom = A.matrix_window(split, 0, nrows-split, ncols)

    top_pivots = strassen_echelon(top, cutoff)
    # effectively "multiplied" top row by A^{-1}
    #     [  I  A'B ]
    #     [  C   D  ]

    top_pivot_intervals = int_range(top_pivots)
    top_h = len(top_pivots)

    if top_h == 0:
        #                                 [ 0 0 ]
        # the whole top is a zero matrix, [ C D ]. Run echelon on the bottom
        bottom_pivots = strassen_echelon(bottom, cutoff)
        #             [  0   0  ]
        # we now have [  I  C'D ], proceed to sorting

    else:
        bottom_cut = max(top_pivots) + 1
        bottom_left = bottom.matrix_window(0, 0, nrows-split, bottom_cut)

        if top_h == ncols:
            bottom.set_to_zero()
            # [ I ]
            # [ 0 ]
            # proceed to sorting

        else:
            if bottom_cut == top_h:
                clear = bottom_left
            else:
                clear = bottom_left.to_matrix().matrix_from_cols(top_pivots).matrix_window() # TODO: read only, can I do this faster? Also below
            # Subtract off C time top from the bottom_right
            if bottom_cut < ncols:
                bottom_right = bottom.matrix_window(0, bottom_cut, nrows-split, ncols-bottom_cut)
                subtract_strassen_product(bottom_right, clear, top.matrix_window(0, bottom_cut, top_h, ncols-bottom_cut), cutoff);
            # [  I      A'B   ]
            # [  *   D - CA'B ]

            # Now subtract off C times the top from the bottom_left (pivots -> 0)
            if bottom_cut == top_h:
                bottom_left.set_to_zero()
                bottom_start = bottom_cut

            else:
                for cols in top_pivot_intervals:
                    bottom_left.matrix_window(0, cols[0], nrows-split, cols[1]).set_to_zero()
                non_pivots = int_range(0, bottom_cut) - top_pivot_intervals
                for cols in non_pivots:
                    if cols[0] == 0: continue
                    prev_pivot_count = len(top_pivot_intervals - int_range(cols[0]+cols[1], bottom_cut - cols[0]+cols[1]))
                    subtract_strassen_product(bottom_left.matrix_window(0, cols[0], nrows-split, cols[1]),
                                              clear.matrix_window(0, 0, nrows-split, prev_pivot_count),
                                              top.matrix_window(0, cols[0], prev_pivot_count, cols[1]),
                                              cutoff)
                bottom_start = non_pivots._intervals[0][0]
            # [  I      A'B   ]
            # [  0   D - CA'B ]

            # Now recursively do echelon form on the bottom
            bottom_pivots_rel = strassen_echelon(bottom.matrix_window(0, bottom_start, nrows-split, ncols-bottom_start), cutoff)
            # [  I  A'B ]
            # [  0  I F ]
            bottom_pivots = []
            for pivot in bottom_pivots_rel:
                bottom_pivots.append(pivot + bottom_start)
            bottom_h = len(bottom_pivots)

            if bottom_h + top_h == ncols:
                top.matrix_window(0, bottom_cut, split, ncols-bottom_cut).set_to_zero()
                # [ I 0 ]
                # [ 0 I ]
                # proceed to sorting

            elif bottom_h == 0:
                pass
                # [  I  A'B ]
                # [  0   0  ]
                # proceed to sorting

            else:
                #     [  I  A'B ]  =  [  I  E  G  ]
                # let [  0  I F ]  =  [  0  I  F  ]
                top_cut = max(max(bottom_pivots) + 1, bottom_cut)

                # Note: left with respect to leftmost non-zero column of bottom
                top_left = top.matrix_window(0, bottom_start, top_h, top_cut - bottom_start)

                if top_cut - bottom_start == bottom_h:
                    clear = top_left
                else:
                    clear = top_left.to_matrix().matrix_from_cols(bottom_pivots_rel).matrix_window()

                # subtract off E times bottom from top right
                if top_cut < ncols:
                    top_right = top.matrix_window(0, top_cut, top_h, ncols - top_cut)
                    subtract_strassen_product(top_right, clear, bottom.matrix_window(0, top_cut, top_h, ncols - top_cut), cutoff);
                # [  I  *  G - EF ]
                # [  0  I     F   ]

                # Now subtract of E times bottom from top left
                if top_cut - bottom_start == bottom_h:
                    top_left.set_to_zero()

                else:
                    bottom_pivot_intervals = int_range(bottom_pivots)
                    for cols in bottom_pivot_intervals:
                        top.matrix_window(0, cols[0], split, cols[1]).set_to_zero()
                    non_pivots = int_range(bottom_start, top_cut - bottom_start) - bottom_pivot_intervals - top_pivot_intervals
                    for cols in non_pivots:
                        if cols[0] == 0: continue
                        prev_pivot_count = len(bottom_pivot_intervals - int_range(cols[0]+cols[1], top_cut - cols[0]+cols[1]))
                        subtract_strassen_product(top.matrix_window(0, cols[0], split, cols[1]),
                                                  clear.matrix_window(0, 0, split, prev_pivot_count),
                                                  bottom.matrix_window(0, cols[0], prev_pivot_count, cols[1]),
                                                  cutoff)
                # [  I  0  G - EF ]
                # [  0  I     F   ]
                # proceed to sorting

    # subrows already sorted...maybe I could do this more efficiently in cases with few pivot columns (e.g. merge sort)

    pivots = top_pivots
    pivots.extend(bottom_pivots)
    pivots.sort()

    cdef int i, cur_row
    for cur_row from 0 <= cur_row < len(pivots):
        pivot = pivots[cur_row]
        for i from cur_row <= i < nrows:
            if not A.element_is_zero(i, pivot):
                break
        if i > cur_row:
            A.swap_rows(i, cur_row)

    return pivots

################################

# lots of room for optimization....
# eventually, should I just pass these around rather than lists of ints for pivots?
# would need new from_cols
class int_range:
    r"""
    Useful class for dealing with pivots in the strassen echelon, could have much more general application
    AUTHORS:
      -- Robert Bradshaw
    """
    def __init__(self, indices=None, range=None):
        if indices is None:
            self._intervals = []
            return
        elif not range is None:
            self._intervals = [(int(indices), int(range))]
        else:
            self._intervals = []
            if len(indices) == 0:
                return
            indices.sort()
            start = None
            last = None
            for ix in indices:
                if last is None:
                    start = ix
                elif ix-last > 1:
                    self._intervals.append((start, last-start+1))
                    start = ix
                last = ix
            self._intervals.append((start, last-start+1))

    def __repr__(self):
        return str(self._intervals)

    def intervals(self):
        return self._intervals

    def to_list(self):
        all = []
        for iv in self._intervals:
            for i in range(iv[0], iv[0]+iv[1]):
                all.append(i)
        return all

    def __iter__(self):
        return self._intervals.__iter__()

    def __len__(self):
        len = 0
        for iv in self._intervals:
            len = len + iv[1]
        return len

    # Yes, these two could be a lot faster...
    # Basically, this class is for abstracting away what I was trying to do by hand in several places
    def __add__(self, other):
        all = self.to_list()
        for i in other.to_list():
            all.append(i)
        return int_range(all)

    def __sub__(self, other):
        all = self.to_list()
        for i in other.to_list():
            if i in all:
                all.remove(i)
        return int_range(all)

    def __mul__(self, other):
        """
        In the boolean sense, i.e. intersection
        """
        intersection = []
        all = self.to_list()
        for i in other.to_list():
            if i in all:
                intersection.append(i)
        return int_range(intersection)

#######################

def __reduce__Matrix_generic(cls, attrs, parent, dict, determinant,
                         sparse_columns, sparse_rows,
                         mutability):
    M = cls.__new__(cls)
    M._init(attrs, parent, dict, determinant,
            sparse_columns, sparse_rows,
            mutability)
    return M


