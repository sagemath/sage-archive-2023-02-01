"""
Matrix windows
"""

#########################################################################
# Generic matrix windows, which are used for block echelon and strassen #
# algorithms.                                                           #
#########################################################################
cdef class MatrixWindow:

    def __init__(MatrixWindow self, Matrix matrix, Py_ssize_t row, Py_ssize_t col, Py_ssize_t nrows, Py_ssize_t ncols):

        self._matrix = matrix
        self._row = row
        self._col = col
        self._nrows = nrows
        self._ncols = ncols

    def __repr__(self):
        return "Matrix window of size %s x %s at (%s,%s):\n%s"%(
            self._nrows, self._ncols, self._row, self._col, self._matrix)

    def set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        self._matrix.set_unsafe(i + self._row, j + self._col, x)

    def get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        return self._matrix.get_unsafe(i + self._row, j + self._col)

    def __setitem__(self, ij, x):
        cdef Py_ssize_t i, j
        if PyTuple_Check(ij):
            # ij is a tuple, so we get i and j efficiently, construct corresponding integer entry.
            if PyTuple_Size(ij) != 2:
                raise IndexError, "index must be an integer or pair of integers"
            i = <object> PyTuple_GET_ITEM(ij, 0)
            j = <object> PyTuple_GET_ITEM(ij, 1)
            if i<0 or i >= self._nrows or j<0 or j >= self._ncols:
                raise IndexError, "matrix index out of range"
            self.set_unsafe(i, j, x)
        else:
            # If ij is not a tuple, coerce to an integer and set the row.
            i = ij
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, x)

    def __getitem__(self, ij):
        cdef Py_ssize_t i, j
        cdef object x

        if PyTuple_Check(ij):
            # ij is a tuple, so we get i and j efficiently, construct corresponding integer entry.
            if PyTuple_Size(ij) != 2:
                raise IndexError, "index must be an integer or pair of integers"
            i = <object> PyTuple_GET_ITEM(ij, 0)
            j = <object> PyTuple_GET_ITEM(ij, 1)
            if i<0 or i >= self._nrows or j<0 or j >= self._ncols:
                raise IndexError, "matrix index out of range"
            return self.get_unsafe(i, j)
        else:
            # If ij is not a tuple, coerce to an integer and get the row.
            i = ij
            return self.row(i)

    def matrix(MatrixWindow self):
        """
        Returns the underlying matrix that this window is a view of.
        """
        return self._matrix


    def to_matrix(MatrixWindow self):
        """
        Returns an actual matrix object representing this view.
        """
        a = self._matrix.new_matrix(self._nrows, self._ncols)
        a.matrix_window().set_to(self)
        return a

    def matrix_window(MatrixWindow self, Py_ssize_t row=0, Py_ssize_t col=0, Py_ssize_t n_rows=-1, Py_ssize_t n_cols=-1):
        """
        Returns a matrix window relative to this window of the underlying matrix.
        """
        if n_rows == -1:
            n_rows = self._nrows - row
        if n_cols == -1:
            n_cols = self._ncols - col
        if row == 0 and col == 0 and n_rows == self._nrows and n_cols == self._ncols:
            return self
        return self._matrix.matrix_window(self._row + row, self._col + col, n_rows, n_cols)

    def nrows(MatrixWindow self):
        return self._nrows

    def ncols(MatrixWindow self):
        return self._ncols

    def set_to(MatrixWindow self, MatrixWindow A):
        """
        Change self, making it equal A.
        """
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j))

    def set_to_zero(MatrixWindow self):
        cdef Py_ssize_t i, j
        z = self._matrix._base_ring(0)
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, z)

    def add(MatrixWindow self, MatrixWindow A):
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, self.get_unsafe(i, j) + A.get_unsafe(i, j))

    def subtract(MatrixWindow self, MatrixWindow A):
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, self.get_unsafe(i, j) - A.get_unsafe(i, j))

    def set_to_sum(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        if self._nrows != B._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j) + B.get_unsafe(i, j))

    def set_to_diff(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self.set_unsafe(i, j, A.get_unsafe(i, j) - B.get_unsafe(i, j))

    def set_to_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = 0
                for k from 0 <= k < A._ncols:
                    s = s + A.get_unsafe(i, k) * B.get_unsafe(k, j)
                self.set_unsafe(i, j, s)

    def add_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = self.get_unsafe(i, j)
                for k from 0 <= k < A._ncols:
                    s = s + A.get_unsafe(i, k) * B.get_unsafe(k, j)
                self.set_unsafe(i, j, s)

    def subtract_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = self.get_unsafe(i, j)
                for k from 0 <= k < A._ncols:
                    s = s - A.get_unsafe(i, k) * B.get_unsafe(k, j)
                self.set_unsafe(i, j, s)

    def swap_rows(MatrixWindow self, Py_ssize_t a, Py_ssize_t b):
        self._matrix.swap_rows(self._row + a, self._row + b)

    def echelon_in_place(MatrixWindow self):
        """
        Calculate the echelon form of this matrix, returning the list of pivot columns
        """
        echelon = self.to_matrix()
        echelon.echelonize() # TODO: read only, only need to copy pointers
        self.set_to(echelon.matrix_window())
        return echelon.pivots()

    def element_is_zero(MatrixWindow self, Py_ssize_t i, Py_ssize_t j):
        return self._matrix[i+self._row, j+self._col] == 0

    def new_empty_window(MatrixWindow self, Py_ssize_t nrows, Py_ssize_t ncols):
        return self._matrix.new_matrix(nrows, ncols).matrix_window()


