r"""
Base class for sparse matrices
"""

cimport matrix
cimport matrix0
from   sage.structure.element    cimport Element

include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'
include '../ext/python.pxi'
include '../ext/interrupt.pxi'

cdef extern from "Python.h":
    PyObject* PyTuple_GET_ITEM0 "PyTuple_GET_ITEM" (PyObject*  p, Py_ssize_t pos)
    PyObject* PyList_GET_ITEM0 "PyList_GET_ITEM" (PyObject*  list, Py_ssize_t i)
    Py_ssize_t PyNumber_AsSsize_t(PyObject* o, PyObject* exc)


import sage.matrix.matrix_space

cdef class Matrix_sparse(matrix.Matrix):

    cdef int is_sparse_c(self):
        return 1

    cdef int is_dense_c(self):
        return 0

    def change_ring(self, ring):
        if ring is self._base_ring:
            if self._mutability._is_immutable:
                return self
            return self.copy()

        M = sage.matrix.matrix_space.MatrixSpace(ring, self._nrows, self._ncols, sparse=self.is_sparse())
        return M(self.dict(), coerce=True, copy=False)

    def __copy__(self):
        """
        Return a copy of this matrix.  Changing the entries of the
        copy will not change the entries of this matrix.

        EXAMPLES:
            sage: A = matrix(QQ['x,y'], 2, [0,-1,2,-2], sparse=True); A
            [ 0 -1]
            [ 2 -2]
            sage: B = copy(A); B
            [ 0 -1]
            [ 2 -2]
            sage: B is A
            False
            sage: B[0,0]=10; B
            [10 -1]
            [ 2 -2]
            sage: A
            [ 0 -1]
            [ 2 -2]
        """
        return self.new_matrix(entries=self.dict(), coerce=False, copy=False)

    def __hash__(self):
        """
        Return the hash of this matrix.

        Equal matrices should have equal hashes, even if one is sparse and
        the other is dense.

        EXAMPLES:
            sage: m = matrix(2, range(6), sparse=True)
            sage: m.set_immutable()
            sage: hash(m)
            5

        The sparse and dense hashes should agree:
            sage: d = m.dense_matrix()
            sage: d.set_immutable()
            sage: hash(d)
            5

            sage: A = Matrix(ZZ[['t']], 2, 2, range(4), sparse=True)
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()
            sage: B = A.copy(); B.set_immutable()
            sage: hash(A) == hash(B)
            True
        """
        return self._hash()

    cdef long _hash(self) except -1:
        x = self.fetch('hash')
        if not x is None: return x

        if not self._mutability._is_immutable:
            raise TypeError, "mutable matrices are unhashable"

        v = self._dict()
        cdef long i, h
        h = 0
        for ij, x in v.iteritems():
            # The following complicated line is the Python/C API optimized version
            # of the following:
            #           i = ij[0]*self._ncols + ij[1]

            i = PyInt_AS_LONG(<object>PyTuple_GET_ITEM(ij,0)) * self._ncols + \
                PyInt_AS_LONG(<object>PyTuple_GET_ITEM(ij,1))

            h = h ^ (i*PyObject_Hash(x))
        if h == -1: h = -2

        self.cache('hash', h)
        return h

    def _multiply_classical(Matrix_sparse left, Matrix_sparse right):
        """
        EXAMPLES:
            sage: A = matrix(QQ['x,y'], 2, [0,-1,2,-2], sparse=True)
            sage: type(A)
            <type 'sage.matrix.matrix_generic_sparse.Matrix_generic_sparse'>
            sage: B = matrix(QQ['x,y'], 2, [-1,-1,-2,-2], sparse=True)
            sage: A*B
            [2 2]
            [2 2]
        """
        cdef Py_ssize_t row, col, row_start, k1, k2, len_left, len_right, a, b
        left_nonzero = left.nonzero_positions(copy=False, column_order=False)
        right_nonzero = right.nonzero_positions(copy=False, column_order=True)
        len_left = len(left_nonzero)
        len_right = len(right_nonzero)

        e = {}
        k1 = 0
        _sig_on
        while k1 < len_left:
            row_start = k1
            row = get_ij(left_nonzero, row_start, 0)
            k2 = 0
            while k2 < len_right:
                col = get_ij(right_nonzero, k2, 1)
                sum = None
                k1 = row_start
                while k1 < len_left and get_ij(left_nonzero,k1,0) == row and \
                          k2 < len_right and get_ij(right_nonzero,k2,1) == col:
                    a = get_ij(left_nonzero, k1,1)
                    b = get_ij(right_nonzero,k2,0)
                    if a == b:
                        if sum is None:
                            sum = left.get_unsafe(row,a)*right.get_unsafe(a,col)
                        else:
                            sum = sum + left.get_unsafe(row,a)*right.get_unsafe(a,col)
                        k1 = k1 + 1
                        k2 = k2 + 1
                    elif a < b:
                        k1 = k1 + 1
                    else:
                        k2 = k2 + 1
                if not sum is None:
                    e[row, col] = sum
                while k2 < len_right and get_ij(right_nonzero,k2,1) == col:
                    k2 = k2 + 1
            while k1 < len_left and get_ij(left_nonzero,k1,0) == row:
                k1 = k1 + 1
        _sig_off
        return left.new_matrix(left._nrows, right._ncols, entries=e, coerce=False, copy=False)

    def _multiply_classical_with_cache(Matrix_sparse left, Matrix_sparse right):
        """
        This function computes the locations of the end of the rows/columns
        in the non-zero entries list once O(rows+cols) time and space, then
        uses these values in the inner loops. For large matrices this can be
        a 2x or more speedup, but the matrices can no longer be arbitrarily
        large as the runtime and space requirements are no longer functions
        of the total number of entries only.

        EXAMPLES:
            sage: A = matrix(QQ['x,y'], 2, [0,-1,2,-2], sparse=True)
            sage: type(A)
            <type 'sage.matrix.matrix_generic_sparse.Matrix_generic_sparse'>
            sage: B = matrix(QQ['x,y'], 2, [-1,-1,-2,-2], sparse=True)
            sage: A._multiply_classical_with_cache(B)
            [2 2]
            [2 2]
        """
        cdef Py_ssize_t row, col, row_start, k1, k2, len_left, len_right, a, b, i
        cdef Py_ssize_t* next_row
        cdef Py_ssize_t* next_col
        left_nonzero = left.nonzero_positions(copy=False, column_order=False)
        right_nonzero = right.nonzero_positions(copy=False, column_order=True)
        len_left = len(left_nonzero)
        len_right = len(right_nonzero)
        next_row = <Py_ssize_t *> sage_malloc(sizeof(Py_ssize_t) * left._nrows)
        next_col = <Py_ssize_t *> sage_malloc(sizeof(Py_ssize_t) * right._ncols)
        if next_row == NULL or next_col == NULL:
            if next_row != NULL: sage_free(next_row)
            _sig_off
            raise MemoryError, "out of memory multiplying a matrix"

        _sig_on
        i = len_left - 1
        for row from left._nrows > row >= 0:
            next_row[row] = i + 1
            while i >= 0 and get_ij(left_nonzero,i,0) == row:
                i = i - 1
        i = len_right - 1
        for col from right._ncols > col >= 0:
            next_col[col] = i + 1
            while i >= 0 and get_ij(right_nonzero,i,1) == col:
                i = i - 1

        e = {}
        k1 = 0
        while k1 < len_left:
            row_start = k1
            row = get_ij(left_nonzero, row_start, 0)
            k2 = 0
            while k2 < len_right:
                col = get_ij(right_nonzero, k2, 1)
                sum = None
                k1 = row_start
                while k1 < next_row[row] and k2 < next_col[col]:
                    a = get_ij(left_nonzero, k1,1)
                    b = get_ij(right_nonzero,k2,0)
                    if a == b:
                        if sum is None:
                            sum = left.get_unsafe(row,a)*right.get_unsafe(a,col)
                        else:
                            sum = sum + left.get_unsafe(row,a)*right.get_unsafe(a,col)
                        k1 = k1 + 1
                        k2 = k2 + 1
                    elif a < b:
                        k1 = k1 + 1
                    else:
                        k2 = k2 + 1
                if not sum is None:
                    e[row, col] = sum
                k2 = next_col[col]
            k1 = next_row[row]

        sage_free(next_row)
        sage_free(next_col)
        _sig_off

        return left.new_matrix(left._nrows, right._ncols, entries=e, coerce=False, copy=False)

    cdef int _will_use_strassen(self, matrix0.Matrix right) except -2:
        # never use Strassen for sparse matrix multiply
        return 0

    def _pickle(self):
        version = -1
        data = self._dict()  # dict of all elements
        return data, version

    def _unpickle_generic(self, data, int version):
        cdef Py_ssize_t i, j, k
        if version == -1:
            for ij, x in data.iteritems():
                i = PyNumber_AsSsize_t(PyTuple_GET_ITEM0(<PyObject*> ij, 0), NULL)
                j = PyNumber_AsSsize_t(PyTuple_GET_ITEM0(<PyObject*> ij, 1), NULL)
                self.set_unsafe(i, j, x)
        else:
            raise RuntimeError, "unknown matrix version (=%s)"%version

    cdef int _cmp_c_impl(self, Element right) except -2:
        return cmp(self._dict(), right._dict())

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
        cdef Matrix_sparse A
        A = self.new_matrix(self._ncols, self._nrows)

        nz = self.nonzero_positions(copy=False)
        cdef Py_ssize_t i, j, k
        for k from 0 <= k < len(nz):
            i = get_ij(nz, k, 0)
            j = get_ij(nz, k, 1)
            A.set_unsafe(j,i,self.get_unsafe(i,j))
        return A

    def antitranspose(self):
        cdef Matrix_sparse A
        A = self.new_matrix(self._ncols, self._nrows)

        nz = self.nonzero_positions(copy=False)
        cdef Py_ssize_t i, j, k
        for k from 0 <= k < len(nz):
            i = get_ij(nz, k, 0)
            j = get_ij(nz, k, 1)
            A.set_unsafe(self._ncols-j-1, self._nrows-i-1,self.get_unsafe(i,j))
        return A



##     def _echelon_in_place_classical(self):
##         """
##         Replace this matrix by its echelon form.

##         INPUT:
##            params -- ignored.
##         """
##         # ALGORITHM:
##         # Since we know nothing about the base field, we use a generic
##         # algorithm.  We first convert to a list of sparse rows, then
##         # directly perform a generic echelon algorithm on that list of
##         # rows.
##         if self.fetch('in_echelon_form'):
##             return
##         K = self.base_ring()
##         ONE = K(1)
##         if not K.is_field():
##             raise ArithmeticError, "The base ring must be a field."
##         X = self.rows()
##         nrows = self.nrows()
##         ncols = self.ncols()
##         pivot_positions = []
##         start_row = 0
##         nrows = self.nrows()
##         ncols = self.ncols()
##         for c in range(ncols):
##             N = []
##             for r in xrange(start_row, nrows):
##                 if X[r].first_nonzero_position() == c:
##                     N.append((X[r].num_nonzero(),r))
##             if len(N) == 0:
##                 continue
##             N.sort()
##             r = N[0][1]
##             leading = X[r].first_nonzero_entry()
##             if leading != 0:
##                 pivot_positions.append(c)
##                 # 1. Rescale
##                 X[r].rescale(ONE/leading)
##                 # 2. Swap
##                 X[r], X[start_row] = X[start_row], X[r]
##                 # 3. Clear column
##                 for i in range(nrows):
##                     if i != start_row:
##                         s = X[i][c]
##                         if s != 0:
##                             X[i] = X[i].add(X[start_row], -s)
##             # endif
##             start_row = start_row + 1
##         #endfor
##         if self.is_immutable():
##             self.__pivots = pivot_positions
##             E = Matrix_generic_sparse_from_rows(X)
##             E.__pivots = pivot_positions
##             self.__echelon_form = E
##         misc.verbose("Finished generic echelon.",t)
##         return E

##################################
# Helper code

cdef Py_ssize_t get_ij(v, Py_ssize_t i, int j):
    return PyNumber_AsSsize_t(PyTuple_GET_ITEM0(PyList_GET_ITEM0(<PyObject*>v, i), j), NULL)

#cdef Py_ssize_t get_ij(v, Py_ssize_t i, int j):
#    return PyNumber_AsSsize_t(<object>PyTuple_GET_ITEM(
#              <object>PyList_GET_ITEM(v, i), j), <object>NULL)
