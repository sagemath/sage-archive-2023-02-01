
cdef class MatrixWindow:

    def __init__(MatrixWindow self, matrix, int row, int col, int nrows, int ncols):
        self._matrix = matrix
        self._row = row
        self._col = col
        self._nrows = nrows
        self._ncols = ncols


    def matrix(MatrixWindow self):
        """
        Returns the underlying matrix that this window is a view of.
        """
        return self._matrix

    def to_matrix(MatrixWindow self):
        """
        Returns a actual matrix object representing this view. (Copy)
        """
        m = self._matrix.new_matrix(self._nrows, self._ncols)
        w = m.window()
        w.set_to(self)
        return m

    def window(MatrixWindow self, int row, int col, int n_rows, int n_cols):
        """
        Returns a matrix window relative to this window of the underlying matrix.
        """
        return self._matrix.window(self._matrix, _row + row, _col + col, n_rows, n_cols)

    def nrows(MatrixWindow self):
        return _nrows

    def ncols(MatrixWindow self):
        return _ncols


    def set_to(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef int start, self_ix
        cdef int A_ix
        for i from 0 <= i < self._nrows:
            A_ix = self._matrix._row_indices[i+A._row] + a._col
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix._entries[self_ix] = A._matrix._entries[A_ix]
                A_ix = A_ix + 1


    def set_to_zero(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef int start, self_ix
        cdef int A_ix
        zero = self._matrix.base_ring(0)
        for i from 0 <= i < self._nrows:
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix._entries[self_ix] = zero


    def add(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef int start, self_ix
        cdef int A_ix
        for i from 0 <= i < self._nrows:
            A_ix = self._matrix._row_indices[i+A._row] + a._col
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix._entries[self_ix] = slef._matrix._entries[self_ix] + A._matrix._entries[A_ix]
                A_ix = A_ix + 1


    def subtract(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef int start, self_ix
        cdef int A_ix
        for i from 0 <= i < self._nrows:
            A_ix = self._matrix._row_indices[i+A._row] + a._col
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix._entries[self_ix] = self._matrix._entries[self_ix] - A._matrix._entries[A_ix]
                A_ix = A_ix + 1


    def set_to_sum(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j
        cdef int start, self_ix
        cdef int A_ix
        for i from 0 <= i < self._nrows:
            A_ix = self._matrix._row_indices[i+A._row] + a._col
            B_ix = self._matrix._row_indices[i+B._row] + b._col
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix._entries[self_ix] = A._matrix._entries[A_ix] + B._matrix._entries[B_ix]
                A_ix = A_ix + 1
                B_ix = B_ix + 1


    def set_to_diff(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j
        cdef int start, self_ix
        cdef int A_ix
        for i from 0 <= i < self._nrows:
            A_ix = self._matrix._row_indices[i+A._row] + a._col
            B_ix = self._matrix._row_indices[i+B._row] + b._col
            for self_ix from self._row_indices[i+self._row] + self._col <= self_ix < self_start + self._ncols:
                self._matrix._entries[self_ix] = A._matrix._entries[A_ix] - B._matrix._entries[B_ix]
                A_ix = A_ix + 1
                B_ix = B_ix + 1


    def set_to_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j, k
        cdef int start, self_ix
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                sum = A._matrix._entries[ A._row_indices[A._row+i]+A._col ] *B._matrix._entries[ B._row_indices[B._row]+B._col+j ]
                for k from 1 <= k < A._ncols:
                    sum = sum + A._matrix._entries[ A._row_indices[A._row+i]+A._col + k ] * B._matrix._entries[ B._row_indices[B._row+k]+B._col+j ]
                self._matrix._entries[ self._row_indices[self_.row+i]+self._col+j ] = sum


    def add_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j, k
        cdef int start, self_ix
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                sum = A._matrix._entries[ A._row_indices[A._row+i]+A._col ] *B._matrix._entries[ B._row_indices[B._row]+B._col+j ]
                for k from 1 <= k < A._ncols:
                    sum = sum + A._matrix._entries[ A._row_indices[A._row+i]+A._col + k ] * B._matrix._entries[ B._row_indices[B._row+k]+B._col+j ]
                self._matrix._entries[ self._row_indices[self_.row+i]+self._col+j ] = self._matrix._entries[ self._row_indices[self_.row+i]+self._col+j ] + sum




