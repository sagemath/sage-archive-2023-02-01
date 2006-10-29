"""
Sparse Matrices Base Class

EXAMPLES:
    sage: ???
"""

import sage.misc.misc as misc
import matrix_space

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
        sage: from sage.matrix.sparse_matrix import Matrix_sparse_from_rows
        sage: Matrix_sparse_from_rows([v])
        [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ]
        sage: Matrix_sparse_from_rows([v, v, v, V(0)])
        [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ]
    """
    cdef size_t i, j

    if not isinstance(X, (list, tuple)):
        raise TypeError, "X (=%s) must be a list or tuple"%X
    if len(X) == 0:
        raise ArithmeticError, "X must be nonempty"
    entries = []
    R = X[0].base_ring()
    ncols = X[0].degree()
    for i from 0 <= i < len(X):
        for j, x in X[i].entries().iteritems():
            entries.append((i,j,x))
    M = matrix_space.MatrixSpace(R, len(X), ncols)
    return M(entries, coerce=False, copy=False)


def _sparse_dot_product(v, w):
    """
    INPUT:
        v and w are dictionaries with integer keys.
    """
    x = set(v.keys()).intersection(set(w.keys()))
    return sum([v[k]*w[k] for k in x])

cdef _convert_sparse_entries_to_dict(entries):
    e = {}
    for i in xrange(len(entries)):
        for j, x in (entries[i].dict()).iteritems():
            e[(i,j)] = x
    return e

cdef class Matrix_sparse(matrix.Matrix):
    r"""
    The \class{Matrix_sparse} class derives from \class{Matrix}, and
    defines functionality for sparse matrices over any base ring.  A
    generic sparse matrix is represented using a dictionary with keys
    pairs $(i,j)$ and values in the base ring.
    """
    def __init__(self, parent, entries=0, coerce=True, copy=True):
        cdef size_t i, j

        matrix.Matrix.__init__(self, parent)
        R = parent.base_ring()
        self._base_ring = R
        self._zero = R(0)

        if not isinstance(entries, (list, dict)):
            x = R(entries)
            entries = {}
            if x != self._zero:
                if self._nrows != self._ncols:
                    raise TypeError, "scalar matrix must be square"
                for i from 0 <= i < self._nrows:
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
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    if x[k] != 0:
                        entries[(i,j)] = x[k]
                    k = k + 1
            copy = False

        if not isinstance(entries, dict):
            raise TypeError, "entries must be a dict"

        if coerce:
            try:
                for k, x in entries.iteritems():
                    entries[k] = R(x)
            except TypeError:
                raise TypeError, "Unable to coerce entries to %s"%R
        elif copy:
            # Make a copy
            entries = dict(entries)

        self._entries = entries

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

    def __copy__(self):
        cdef Matrix_sparse M
        M = Matrix_sparse.__new__(self._parent,0,0,0)
        M._entries = dict(self._entries)
        M._base_ring = self._base_ring
        M._zero = self._zero
        return M

    def __repr__(self):
        return str(self.dense_matrix())

    def __getitem__(self, t):
        cdef size_t r
        if not isinstance(t, tuple) or len(t) != 2:
            try:
                r = t
            except TypeError:
                raise IndexError, "Index of matrix item must be a row and a column."
            return self.row(r)

        cdef size_t i, j
        i, j = t
        self.check_bounds(i, j)
        try:
            return self._entries[(i,j)]
        except KeyError:
            return self._zero

    def __setitem__(self, t, x):
        """
        INPUT:
            A[i, j] -- the i,j of A, and
            A[i]    -- the i-th row of A.
        """
        cdef size_t r

        if not isinstance(t, tuple) or len(t) != 2:
            try:
                r = t
            except TypeError:
                raise IndexError, "Index of matrix item must be a row and a column."
            return self.set_row(r, x)

        cdef size_t i, j
        i, j = t
        self.check_bounds_and_mutability(i, j)

        x = self._base_ring(value)
        if x == 0:
            if self._entries.has_key(ij):
                del self._entries[ij]
            return
        self._entries[ij] = x

    def base_ring(self):
        return self._base_ring

    def row(self, size_t i):
        """
        Return the ith row of this matrix as a sparse vector.

        WARNING: Sparse matrices are stored as a single list of triples (i,j,x),
        so extracting the i-th row is expensive.  This command builds a redundant
        representation of the matrix as a list of sparse vectors, thus doubling
        the memory requirement.
        """
        if i < 0:
            i = self._nrows + i
        return self.rows()[i]

    def rows(self):
        """
        Return a list of the sparse rows of this matrix.  The i-th
        element of this list is the i-th row of this matrix.
        """
        x = self.fetch('rows')
        if not x is None: return x
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
        self.cache('rows', X)
        return X

    def list(self):
        """
        Return all entries of self as a linear list of numbers of rows
        times number of columns entries.

        This is not cached, so it is safe to change the returned list.
        """
        v = [self._base_ring(0)]*(self._nrows * self._ncols)
        for i, j, x in self._entries:
            v[i*self._ncols + j] = x
        return v

    def dict(self):
        """
        Return a copy of the underlying dictionary of self.

        It is safe to change the returned list.
        """
        return dict(self._entries)

    def dense_matrix(self):
        import sage.matrix.matrix
        M = sage.matrix.matrix.MatrixSpace(self.base_ring(),
                                           self._nrows, self._ncols,
                                           sparse = False)(0)
        for i, j, x in self._entries
            M[i,j] = x
        return M

    def nonpivots(self):
        # We convert the pivots to a set so we have fast
        # inclusion testing
        X = set(self.pivots())
        # [j for j in xrange(self.ncols()) if not (j in X)]
        np = []
        for j in xrange(self.ncols()):
            if not (j in X):
                np.append(j)
        return np

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
            E = Matrix_sparse_from_rows(X)
            E.__pivots = pivot_positions
            self.__echelon_form = E
        misc.verbose("Finished generic echelon.",t)
        return E

