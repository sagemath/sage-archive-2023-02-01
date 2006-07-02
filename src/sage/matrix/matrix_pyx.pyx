r"""
Matrices over an arbitrary ring

AUTHORS:
    * William Stein
    * Martin Albrecht (conversion to Pyrex)

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

doc = """
Matrices
"""

import sage.modules.free_module_element
import sage.modules.free_module
import sage.misc.latex
import sage.ext.coerce
from   sage.structure.sequence import _combinations
from   sage.rings.integer_ring import IntegerRing


cdef class Matrix(ModuleElement):
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
            <type 'matrix_pyx.Matrix'>
        """
        ModuleElement.__init__(self, parent)
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

    def __repr__(self):
        nr = self.nrows()
        nc = self.ncols()
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
        ans = []
        for i in range(n):
            tmp = []
            for j in range(n):
                tmp.append(self.get( (i,j) )*vars[j])
            ans.append( sum(tmp) )
        return f(tuple(ans))

    def add_multiple_of_row(self, i, j, s):
        """
        Replace row i by s times row j.
        """
        self._require_mutable()
        s = self.base_ring()(s)
        for c in xrange(self.ncols()):
            self[i,c] = self[i,c] + s*self[j,c]

    def add_multiple_of_column(self, i, j, s):
        """
        Replace column i by s times column j.
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
        """
        if self.__determinant != None:
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
        for i in range(n):
            v = range(n)
            del v[i]
            B = A.matrix_from_columns(v)
            d = d + s*self.get((0,i))*B.determinant()
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
            return sage.ext.coerce.cmp(self, right)
        return cmp(self._entries(), right._entries())

    def __richcmp__(self,right,op):
        res = 0
        if not isinstance(right, Matrix) or right.parent() != self.parent():
            res = sage.ext.coerce.cmp(self, right)
        else:
            res = cmp(self._entries(), right._entries())

        if op == 0:  #<
            return bool(res  < 0)
        if op == 2: #==
            return bool(res == 0)
        if op == 4: #>
            return bool(res  > 0)
        if op == 1: #<=
            return bool(res <= 0)
        if op == 3: #!=
            return bool(res != 0)
        if op == 5: #>=
            return bool(res >= 0)

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
        return sage.ext.coerce.cmp(self, 0)

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
            Full MatrixSpace of 2 by 2 dense matrices over Ring of integers modulo 7
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
                self, right = sage.ext.coerce.canonical_base_coercion(self, right)
            except TypeError:
                raise ArithmeticError, "base rings must be compatible"
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
