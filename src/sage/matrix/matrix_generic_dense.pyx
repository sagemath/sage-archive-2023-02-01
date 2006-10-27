"""
Dense Matrices over a General Commutative Ring
"""

def _convert_dense_entries_to_list(entries):
    # Create matrix from a list of vectors
    e = []
    for v in entries:
        e = e+ v.list()
    copy = False
    return e

include "../ext/interrupt.pxi"

cimport matrix_generic
import matrix_generic

cdef class Matrix_generic_dense(matrix_generic.Matrix):
    r"""
    The \class{Matrix_dense} class derives from \class{Matrix}, and
    defines functionality for dense matrices over any base ring.
    Matrices are represented by a list of elements in the base ring,
    and element access operations are implemented in this class.
    """
    def __init__(self,
                 parent,
                 entries,
                 copy = True,
                 coerce = True):
        matrix_generic.Matrix.__init__(self, parent)

        cdef size_t i, n

        if entries is None:
            entries = 0

        if not isinstance(entries, list):
            try:
                x = parent.base_ring()(entries)
                is_list = 0
            except TypeError:
                try:
                    entries = list(entries)
                    is_list = 1
                except TypeError:
                    raise TypeError, "entries must be coercible to a list or the basering"

        else:
            is_list = 1

        if is_list:

            if len(entries) != self._nrows * self._ncols:
                raise TypeError, "entries has the wrong length"

            if not (coerce or copy):
                self.__entries = entries
            else:
                self.__entries = [None]*(self._nrows*self._ncols)
                n = len(entries)
                if coerce:
                    R = parent.base_ring()
                    for i from 0 <= i < n:
                        self.__entries[i] = R(entries[i])
                else:
                    for i from 0 <= i < n:
                        self.__entries[i] = entries[i]

        else:

            self.__entries = [None]*(self._nrows*self._ncols)
            zero = parent.base_ring()(0)
            for i from 0 <= i < self._nrows * self._ncols:
                self.__entries[i] = zero

            if x != zero:
                if self._nrows != self._ncols:
                    raise TypeError, "nonzero scalar matrix must be square"
                for i from 0 <= i < self._nrows:
                    self.__entries[i*i+i] = x

        self._init_row_indices()

    def __copy__(self):
        return self.__class__(self._parent, self.__entries, copy = True, coerce=False)

    def __deepcopy__(self):
        import copy
        return self.__class__(self._parent, copy.deepcopy(self.__entries), copy = False, coerce=False)

    def copy(self):
        """
        Make a copy of self.

        WARNING: The individual elements aren't themselves copied
        (though the list is copied).

        EXAMPLES:
            sage: R.<x> = QQ['x']
            sage: a = matrix(R,2,[x+1,2/3,  x^2/2, 1+x^3]); a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]

            sage: b = copy(a)

            sage: b[0,0] = 5

            sage: b
            [      5     2/3]
            [1/2*x^2 x^3 + 1]

            sage: a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]

            sage: b = copy(a)

            sage: f = b[0,0]; f[0] = 10

            sage: b
            [ x + 10     2/3]
            [1/2*x^2 x^3 + 1]

            sage: a
            [ x + 10     2/3]
            [1/2*x^2 x^3 + 1]
        """
        return self.__copy__()

    def _init_row_indices(self):
        self._row_indices = <int*> PyMem_Malloc(sizeof(int*) * self._nrows)
        n = 0
        for i from 0 <= i < self._nrows:
            self._row_indices[i] = n
            n = n + self._ncols

    def  __dealloc__(self):
        PyMem_Free(self._row_indices)

    def __getitem__(self, ij):
        """
        INPUT:
            A[i, j] -- the i,j of A, and
            A[i]    -- the i-th row of A.
        """
        cdef size_t i, j
        i, j = ij
        self.check_bounds(i, j)
        return self.__entries[self._row_indices[i] + j]

    def __setitem__(self, ij, value):
        """
        INPUT:
            A[i, j] = value -- set the (i,j) entry of A
            A[i] = value    -- set the ith row of A
        """
        cdef size_t i, j
        i, j = ij
        self.check_bounds_and_is_mutability(i, j)
        self.__entries[self._row_indices[i] + j] = value

    def matrix_window(self, int row=0, int col=0, int nrows=-1, int ncols=-1):
        if nrows == -1:
            nrows = self._nrows
            ncols = self._ncols
        return MatrixWindow(self, row, col, nrows, ncols)

    cdef classical_multiply_cdef(self, matrix_generic.Matrix _right):
        """
        Multiply the matrices self and right using the classical $O(n^3)$
        algorithm.

        EXAMPLES:
        We multiply two matrices over a fairly general ring:

            sage: R.<x,y> = Integers(8)['x,y']
            sage: a = matrix(R,2,[x,y,x^2,y^2]); a
            [  x   y]
            [x^2 y^2]
            sage: a*a
            [  x^2 + x^2*y     y^3 + x*y]
            [x^2*y^2 + x^3   y^4 + x^2*y]
            sage: a.det()^2 == (a*a).det()
            True

        SAGE fully supports degenerate matrices with 0 rows or 0 columns:


        """
        cdef int i, j, k, m, n, r, nr, nc, snc
        cdef object v
        cdef Matrix_dense A, right

        right = _right

        if self._ncols != right._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."

        nr = self._nrows
        nc = right._ncols
        snc = self._ncols

        R = self.base_ring()
        P = self.matrix_space(nr, nc)
        A = Matrix_dense(P)
        v = A.__entries
        zero = R(0)
        r = 0
        for i from 0 <= i < nr:
            m = self._row_indices[i]
            for j from 0 <= j < nc:
                z = zero
                for k from 0 <= k < snc:
                    z = z + self.__entries[m + k] * right.__entries[right._row_indices[k]+j]
                v[r] = z
                r = r + 1
        return A

    def _entries(self):
        return self.__entries

    def list(self):
        return self.__entries

    def antitranspose(self):
        f = []
        e = self.list()
        (nc, nr) = (self.ncols(), self.nrows())
        for j in reversed(xrange(nc)):
            for i in reversed(xrange(nr)):
                f.append(e[i*nc + j])
        return self.new_matrix(nrows = nc, ncols = nr,
                               entries = f, copy=False, coerce=False)

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
                               entries = f, copy=False,
                               coerce=False)



cdef class MatrixWindow:

    def __init__(MatrixWindow self, matrix, int row, int col, int nrows, int ncols):
        self._matrix = matrix
        self._row = row
        self._col = col
        self._nrows = nrows
        self._ncols = ncols

    def __repr__(self):
        return "Matrix window of size %s x %s at (%s,%s):\n%s"%(
            self._nrows, self._ncols, self._row, self._col, self._matrix)

    def matrix(MatrixWindow self):
        """
        Returns the underlying matrix that this window is a view of.

        EXAMPLES:

        """
        return self._matrix


    def to_matrix(MatrixWindow self):
        """
        Returns an actual matrix object representing this view. (Copy)
        """
        entries = self.list()
        return self._matrix.new_matrix(self._nrows, self._ncols, entries=entries,
                                       coerce=False, copy=False)

    def list(MatrixWindow self):
        v = self._matrix.__entries
        w = [None]*(self._nrows*self._ncols)
        cdef int i, j, k, l
        k = 0
        for i from 0 <= i < self._nrows:
            l = self._matrix._row_indices[i]
            for j from 0 <= j < self._ncols:
                w[k] = v[l + j]
                k = k + 1
        return w

    def matrix_window(MatrixWindow self, int row, int col, int n_rows, int n_cols):
        """
        Returns a matrix window relative to this window of the underlying matrix.
        """
        return self._matrix.matrix_window(self._matrix, _row + row, _col + col, n_rows, n_cols)

    def nrows(MatrixWindow self):
        return self._nrows

    def ncols(MatrixWindow self):
        return self._ncols


    def set_to(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef int start, self_ix
        cdef int A_ix
        for i from 0 <= i < self._nrows:
            A_ix = self._matrix._row_indices[i+A._row] + a._col
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix.__entries[self_ix] = A._matrix.__entries[A_ix]
                A_ix = A_ix + 1


    def set_to_zero(MatrixWindow self):
        cdef int i, j
        cdef int start, self_ix
        zero = self._matrix.base_ring(0)
        for i from 0 <= i < self._nrows:
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix.__entries[self_ix] = zero


    def add(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef int start, self_ix
        cdef int A_ix
        for i from 0 <= i < self._nrows:
            A_ix = A._matrix._row_indices[i+A._row] + A._col
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix.__entries[self_ix] = slef._matrix.__entries[self_ix] + A._matrix.__entries[A_ix]
                A_ix = A_ix + 1


    def subtract(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef int start, self_ix
        cdef int A_ix
        for i from 0 <= i < self._nrows:
            A_ix = A._matrix._row_indices[i+A._row] + a._col
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix.__entries[self_ix] = self._matrix.__entries[self_ix] - A._matrix.__entries[A_ix]
                A_ix = A_ix + 1


    def set_to_sum(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j
        cdef int start, self_ix
        cdef int A_ix
        for i from 0 <= i < self._nrows:
            A_ix = A._matrix._row_indices[i+A._row] + A._col
            B_ix = B._matrix._row_indices[i+B._row] + B._col
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix.__entries[self_ix] = A._matrix.__entries[A_ix] + B._matrix.__entries[B_ix]
                A_ix = A_ix + 1
                B_ix = B_ix + 1


    def set_to_diff(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j
        cdef int start, self_ix
        cdef int A_ix
        for i from 0 <= i < self._nrows:
            A_ix = A._matrix._row_indices[i+A._row] + A._col
            B_ix = B._matrix._row_indices[i+B._row] + B._col
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix.__entries[self_ix] = A._matrix.__entries[A_ix] - B._matrix.__entries[B_ix]
                A_ix = A_ix + 1
                B_ix = B_ix + 1


    def set_to_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j, k
        cdef int start, self_ix
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                sum = A._matrix.__entries[ A._row_indices[A._row+i]+A._col ] *B._matrix.__entries[ B._row_indices[B._row]+B._col+j ]
                for k from 1 <= k < A._ncols:
                    sum = sum + A._matrix.__entries[ A._row_indices[A._row+i]+A._col + k ] * B._matrix.__entries[ B._row_indices[B._row+k]+B._col+j ]
                self._matrix.__entries[ self._row_indices[self_.row+i]+self._col+j ] = sum


    def add_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j, k
        cdef int start, self_ix
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                sum = A._matrix.__entries[ A._row_indices[A._row+i]+A._col ] *B._matrix.__entries[ B._row_indices[B._row]+B._col+j ]
                for k from 1 <= k < A._ncols:
                    sum = sum + A._matrix.__entries[ A._row_indices[A._row+i]+A._col + k ] * B._matrix.__entries[ B._row_indices[B._row+k]+B._col+j ]
                self._matrix.__entries[ self._row_indices[self_.row+i]+self._col+j ] = self._matrix.__entries[ self._row_indices[self_.row+i]+self._col+j ] + sum


    def new_empty_window(MatrixWindow self, int nrows, int ncols, zero=True):
        return self._matrix.new_matrix(nrows,ncols, zero=zero).matrix_window()
