"""
Sparse matrices

EXAMPLES:
    sage: M = MatrixSpace(QQ, 2, 3, sparse=True)
    sage: A = M([1,2,3, 1,1,1])
    sage: A
    [1 2 3]
    [1 1 1]
    sage: A.echelon_form()
    [ 1  0 -1]
    [ 0  1  2]

    sage: M = MatrixSpace(QQ, 1000,1000, sparse=True)
    sage: A = M(0)
    sage: A[1,1] = 5


    sage: from sage.matrix.sparse_matrix import SparseMatrix
    sage: x = SparseMatrix(QQ, 5,10)
    sage: x.randomize(5)
    sage: x.echelon_form()       # random output
    [
    1, 0, 0, 0, 0, 0, 0, 0, 0, 10/29,
    0, 1, 0, 0, 0, 0, 0, 0, 0, -4/29,
    0, 0, 1, 0, 0, 0, 0, 0, 0, -12/29,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 24/29,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 4/29
    ]
"""
import random, weakref

import sage.rings.arith as arith
import sage.misc.misc as misc
import sage.rings.integer_ring as integer_ring
import sage.rings.rational_field as rational_field
import sage.rings.rational as rational
import sage.rings.integer as integer
import sage.rings.rational_field as rational_field
import sage.rings.integer_ring as integer_ring
import sage.rings.ring as ring

from sage.misc.proof import proof

START_PRIME = 20011  # used for modular algorithms

#################################################################
#
# Generic sparse matrices
#
#################################################################

def Sparse_matrix(cls, base_ring, nrows, ncols, entries=[],
                coerce=True, sort=True, copy=True):
      nrows = int(nrows)
      #if not isinstance(nrows, int):
      #    raise TypeError, "nrows must be an int"
      ncols = int(ncols)
      #if not isinstance(ncols, int):
      #    raise TypeError, "ncols must be an int"
      if not isinstance(entries, list):
          raise TypeError, "entries must be a list"
      if base_ring is rational_field.RationalField():
          return object.__new__(Sparse_matrix_rational)
      return object.__new__(Sparse_matrix_generic)
      # TODO: fit the above matrices into this tree

def Matrix_sparse_from_rows(X):
    """
    INPUT:
        X -- nonempty list of SparseVector rows

    OUTPUT:
        Sparse_matrix with those rows.

    EXAMPLES:
        sage: V = VectorSpace(QQ,20,sparse=True)
        sage: v = V(0)
        sage: v[9] = 4
        sage: from sage.matrix.sparse_matrix import *
        sage: SparseMatrix_from_rows([v])
        [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ]
        sage: SparseMatrix_from_rows([v, v, v, V(0)])
        [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ]
    """
    if not isinstance(X, (list, tuple)):
        raise TypeError, "X (=%s) must be a list or tuple"%X
    if len(X) == 0:
        raise ArithmeticError, "X must be nonempty"
    entries = []
    R = X[0].base_ring()
    ncols = X[0].degree()
    for i in range(len(X)):
        for j, x in X[i].entries().iteritems():
            entries.append((i,j,x))
    return Sparse_matrix(R, len(X), ncols, entries, coerce=False, sort=False)


#################################################################
#
# Generic sparse matrix over any ring
#
#################################################################

# TODO: should I call this Matrix_generic_sparse? (Possible name conflict with old matrix.py)
cdef class Matrix_sparse(matrix_generic.Matrix):
    """
    A generic sparse matrix.
    """
    def __init__(self, base_ring, nrows, ncols, entries=[],
                 coerce=True, sort=True, copy=True):
        self.__nrows = nrows
        self.__ncols = ncols
        self.__base_ring = base_ring
        self.set_entries(entries, coerce, sort, copy)

    def __repr__(self):
        return str(self.dense_matrix())

    def copy(self):
        E = list(self.entries())
        return SparseMatrix(self.__base_ring, self.__nrows, \
                self.__ncols, E, coerce=False, sort=False, copy=False)

    def __add__(self, other):
        if not isinstance(other, Sparse_matrix_generic):
            raise TypeError
        if self.nrows() != other.nrows() or self.ncols() != other.ncols():
            raise ArithmeticError, "Matrices must have the same number of rows and columns."
        if self.base_ring() != other.base_ring():
            raise ArithmeticError, "Matrices must have the same base ring."

        # Compute the sum of two sparse matrices.
        # This is complicated because of how we represent sparse matrices as tuples (i,j,x).
        # Algorithm:
        #   1. Sort both entry lists.
        #   2. March through building a new list, adding when the two i,j are the same.
        v = self.__entries
        v.sort()
        w = other.__entries
        w.sort()
        s = []
        i = 0  # pointer into self
        j = 0  # pointer into other
        while i < len(v) and j < len(w):
            vi = (v[i][0], v[i][1])
            wj = (w[j][0], w[j][1])
            if vi < wj:
                s.append(tuple(v[i]))
                i = i + 1
            elif vi > wj:
                s.append(tuple(w[j]))
                j = j + 1
            else:  # equal
                s.append((vi[0],vi[1],v[i][2] + w[j][2]))
                i = i + 1
                j = j + 1
        if i < len(v):
            while i < len(v):
                s.append(tuple(v[i]))
                i = i + 1
        if j < len(w):
            while j < len(w):
                s.append(tuple(w[j]))
                j = j + 1
        A = SparseMatrix(self.base_ring(), self.nrows(),
                         self.ncols(), s, coerce=False)
        return A

    def __sub__(self, right):
        return self + (-1)*right

    def __neg__(self):
        return (-1)*self

    def __getitem__(self, i):
        return self.row(i)

    def row(self, i):
        """
        Return the ith row of this matrix as a sparse vector.

        WARNING: Sparse matrices are stored as a single list of triples (i,j,x),
        so extracting the i-th row is expensive.  This command builds a redundant
        representation of the matrix as a list of sparse vectors, thus doubling
        the memory requirement.
        """
        i = int(i)
        #if not isinstance(i, int):
        #    raise TypeError, "i must be an int"
        if i < 0: i = self.nrows() + i
        return self.rows()[i]

    def rows(self):
        """
        Return a list of the sparse rows of this matrix.  The ith
        element of this list is the i-th row of this matrix.
        """
        if self.__rows is None:
            X = []
            entries = []
            k = 0
            R = self.base_ring()
            n = self.ncols()
            for i, j, x in self.entries():
                if i > k:
                    while len(X) <= k-1:
                        X.append(SparseVector(R, n))
                    X.append(SparseVector(R, n, entries, coerce=False, sort=False))
                    entries = []
                    k = i
                entries.append((j,x))
            while len(X) <= k-1:
                X.append(SparseVector(R, n))
            X.append(SparseVector(R, n, entries, coerce=False, sort=False))
            while len(X) < self.nrows():
                X.append(SparseVector(R, n))
            self.__rows = X
        return self.__rows

    def __getslice__(self, i, j):
        i = int(i)
        #if not  isinstance(i, int):
        #    raise TypeError, "i must be an int"
        j = int(j)
        #if not  isinstance(j, int):
        #    raise TypeError, "j must be an int"
        if i < 0: i = self.nrows() + i
        if j < 0: j = self.nrows() + j
        if i >= self.nrows():
            i = self.nrows() - 1
        if j >= self.nrows():
            j = self.nrows() - 1
        return self.matrix_from_rows(xrange(i,j))

    def nrows(self):
        return self.__nrows

    def ncols(self):
        return self.__ncols

    def base_ring(self):
        return self.__base_ring

    def entries(self):
        return self.__entries

    def list(self):
        return self.__entries

    def dict(self):
        X = {}
        for i, j, x in self.entries():
            X[(i,j)] = x
        return X

    def dense_matrix(self):
        """
        todo
        """
        import sage.matrix.matrix
        M = sage.matrix.matrix.MatrixSpace(self.base_ring(),
                                     self.nrows(),
                                     self.ncols(),
                                     sparse = False)(0)
        for i, j, x in self.list():
            M[i,j] = x
        return M

    def nonpivots(self):
        # We convert the pivots to a set so we have VERY FAST
        # inclusion testing
        X = set(self.pivots())
        # [j for j in xrange(self.ncols()) if not (j in X)]
        np = []
        for j in xrange(self.ncols()):
            if not (j in X):
                np.append(j)
        return np

    def pivots(self):
        if self.__pivots is None:
            self.echelon_form()
        return self.__pivots

    def _set_pivots(self, X):
        self.__pivots = X

    def matrix_from_nonpivot_columns(self):
        """
        The sparse matrix got by deleted all pivot columns.
        """
        return self.matrix_from_columns(self.nonpivots())

    def matrix_from_columns(self, columns):
        """
        Returns the sparse submatrix of self composed of the given
        list of columns.

        INPUT:
            columns -- list of int's
        OUTPUT:
            a sparse matrix.
        """
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
                j = j + 1
            else:
                X.append(-1)    # column to be deleted.
        entries2 = []
        for i, j, x in self.entries():
            if j in C:
                entries2.append((i,X[j],x))
        return SparseMatrix(self.base_ring(), self.nrows(),
                            len(C), entries2, coerce=False, sort=False)

    def matrix_from_rows(self, rows):
        """
        Returns the sparse submatrix of self composed of the given
        list of rows.

        INPUT:
            rows -- list of int's
        OUTPUT:
            a sparse matrix.
        """
        R = set(rows)
        if not R.issubset(set(xrange(self.nrows()))):
            raise ArithmeticError, "Invalid rows."
        X = []
        i = 0
        for j in xrange(self.nrows()):
            if j in R:
                X.append(i)
                i = i + 1
            else:
                X.append(-1)    # row to be deleted.
        entries2 = []
        for i, j, x in self.entries():
            if i in R:
                entries2.append((X[i],j,x))
        return SparseMatrix(self.base_ring(), len(R),
                            self.ncols(), entries2, coerce=False, sort=False)


    def set_entries(self, entries, coerce=True, sort=True, copy=True):
        """
        entries is a list of triples (i,j,x) and the
        x must be nonzero.

        This function does *not* check that i and j are in bounds
        or that the x are all nonzero.
        """
        self.__rows = None
        self.__echelon_form = None
        if not coerce:
            if copy:
                self.__entries = list(entries)   # copy
            else:
                self.__entries = entries
            if sort:
                self.__entries.sort()
            return
        self.__entries = []
        R = self.base_ring()
        # self.__entries = [(i,j,R(z)) for i,j,z in entries]
        for t in entries:
            t[2] = R(t[2])
        self.__entries = entries
        if sort:
            self.__entries.sort()

    def randomize(self, sparcity=4, exact=True):
        entries = []
        R = self.base_ring()
        for i in range(self.nrows()):
            if exact:
                r = sparcity
            else:
                r = random.randrange(sparcity)
            X = []
            for j in range(0,r+1):
                x = R.random()
                if x != R(0):
                    k = random.randrange(0,self.ncols())
                    if not (k in X):
                        entries.append((i,k, x))
                        X.append(k)
        self.set_entries(entries, coerce=False, sort=True)

    def transpose(self):
        """
        Returns the transpose of self.
        """
        entries2 = [] # [(j,i,x) for i,j,x in self.entries()]
        for i,j,x in self.entries():
            entries2.append((j,i,x))
        return SparseMatrix(self.base_ring(), self.ncols(),
                            self.nrows(), entries2, coerce=False, sort=True)

    def __rmul__(self, left):
        return self.scalar_multiple(left)

    def scalar_multiple(self, left):
        R = self.base_ring()
        left = R(left)
        if left == R(1):
            return self
        if left == R(0):
            return SparseMatrix(R, self.nrows(), self.ncols(), coerce=False, sort=False)

        X = []
        for i, j, x in self.list():
            X.append((i,j,left*x))
        return SparseMatrix(self.base_ring(), self.nrows(),
                            self.ncols(), X, coerce=False, sort=False)

    def echelon_form(self, params=None):
        """
        Returns the echelon form of this matrix.

        INPUT:
           params -- ignored.
        """
        # ALGORITHM:
        # Since we know nothing about the base field, we use a generic
        # algorithm.  Since sparse matrices are stored as triples
        # (i,j,x), which is not a suitable format for row operations,
        # we first convert to a list of sparse rows, then directly
        # perform a generic echelon algorithm on that list of rows.
        if self.__echelon_form is not None:
            return self.__echelon_form
        t = misc.verbose("Started generic sparse echelon.")
        K = self.base_ring()
        ONE = K(1)
        if not K.is_field():
            raise ArithmeticError, "The base ring must be a field."
        X = self.rows()
        nrows = self.nrows()
        ncols = self.ncols()
        pivot_positions = []
        start_row = 0
        nrows = self.nrows()
        ncols = self.ncols()
        for c in range(ncols):
#            N = [(X[r].num_nonzero(),r) for r in xrange(start_row, nrows) \
#                 if X[r].first_nonzero_position() == c]
            N = []
            for r in xrange(start_row, nrows):
                if X[r].first_nonzero_position() == c:
                    N.append((X[r].num_nonzero(),r))
            if len(N) == 0:
                continue
            N.sort()
            r = N[0][1]
            leading = X[r].first_nonzero_entry()
            if leading != 0:
                pivot_positions.append(c)
                # 1. Rescale
                X[r].rescale(ONE/leading)
                # 2. Swap
                X[r], X[start_row] = X[start_row], X[r]
                # 3. Clear column
                for i in range(nrows):
                    if i != start_row:
                        s = X[i][c]
                        if s != 0:
                            X[i] = X[i].add(X[start_row], -s)
            # endif
            start_row = start_row + 1
        #endfor
        if self.is_immutable():
            self.__pivots = pivot_positions
            E = SparseMatrix_from_rows(X)
            E.__pivots = pivot_positions
            self.__echelon_form = E
        misc.verbose("Finished generic echelon.",t)
        return E








cdef class Matrix_domain_sparse(Matrix_sparse):
    pass

cdef class Matrix_pid_sparse(Matrix_domain_sparse):
    pass

cdef class Matrix_field_sparse(Matrix_pid_sparse):
    pass
