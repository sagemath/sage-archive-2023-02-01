include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"

import matrix_window
cimport matrix_window

from sage.matrix.matrix_modn_dense import Matrix_modn_dense
from sage.matrix.matrix_modn_dense cimport Matrix_modn_dense


cdef class MatrixWindow_modn_dense(matrix_window.MatrixWindow):

    cdef new_empty_window(self, Py_ssize_t nrows, Py_ssize_t ncols):
# the current code is all python, goes through inits, and possibly creates a parent...
# can we get away with something faster?
#       a = MatrixWindow_modn_dense.__new__(self._parent, )?
       a = self._matrix.new_matrix(nrows, ncols)
       return self.new_matrix_window(a, 0, 0, nrows, ncols)

    cdef set_to(MatrixWindow_modn_dense self, MatrixWindow A):
        """
        Change self, making it equal A.
        """
        cdef Py_ssize_t i, j
        cdef mod_int** self_rows
        cdef mod_int** A_rows
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        self_rows = ( <Matrix_modn_dense> self._matrix )._matrix
        A_rows = ( <Matrix_modn_dense> A._matrix )._matrix
        for i from 0 <= i < self._nrows:
            memcpy(self_rows[i+self._row] + self._col, A_rows[i+A._row] + A._col, self._ncols * sizeof(mod_int*))

    cdef set_to_zero(MatrixWindow_modn_dense self):
        cdef Py_ssize_t i, j
        cdef mod_int** rows
        rows = ( <Matrix_modn_dense> self._matrix )._matrix
        for i from self._row <= i < self._row + self._nrows:
            memset(rows[i] + self._col, 0, self._ncols * sizeof(mod_int*))

    cdef add(self, MatrixWindow A):
        cdef Py_ssize_t i, j
        cdef mod_int p
        cdef mod_int** self_matrix
        cdef mod_int** A_matrix
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        self_matrix = ( <Matrix_modn_dense> self._matrix )._matrix
        A_matrix = ( <Matrix_modn_dense> A._matrix )._matrix
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                # I really want to do in place operations here...
                self_matrix[i+self._row][j+self._col] = self_matrix[i+self._row][j+self._col] + A_matrix[i+A._row][j+A._col]
                if self_matrix[i+self._row][j+self._col] >= p:
                    self_matrix[i+self._row][j+self._col] = self_matrix[i+self._row][j+self._col] - p

    cdef subtract(self, MatrixWindow A):
        cdef Py_ssize_t i, j
        cdef mod_int p
        cdef mod_int** self_matrix
        cdef mod_int** A_matrix
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        self_matrix = ( <Matrix_modn_dense> self._matrix )._matrix
        A_matrix = ( <Matrix_modn_dense> A._matrix )._matrix
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                # I really want to do in place operations here...
                self_matrix[i+self._row][j+self._col] = p + self_matrix[i+self._row][j+self._col] - A_matrix[i+A._row][j+A._col]
                if self_matrix[i+self._row][j+self._col] >= p:
                    self_matrix[i+self._row][j+self._col] = self_matrix[i+self._row][j+self._col] - p

    cdef set_to_sum(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        cdef mod_int p
        cdef mod_int** self_matrix
        cdef mod_int** A_matrix
        cdef mod_int** B_matrix
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        if self._nrows != B._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        self_matrix = ( <Matrix_modn_dense> self._matrix )._matrix
        A_matrix = ( <Matrix_modn_dense> A._matrix )._matrix
        B_matrix = ( <Matrix_modn_dense> B._matrix )._matrix
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self_matrix[i+self._row][j+self._col] = A_matrix[i+A._row][j+A._col] + B_matrix[i+B._row][j+B._col]
                if self_matrix[i+self._row][j+self._col] >= p:
                    # I really want to do in place operations here...
                    self_matrix[i+self._row][j+self._col] = self_matrix[i+self._row][j+self._col] - p

    cdef set_to_diff(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        cdef mod_int p
        cdef mod_int** self_matrix
        cdef mod_int** A_matrix
        cdef mod_int** B_matrix
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        if self._nrows != B._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        self_matrix = ( <Matrix_modn_dense> self._matrix )._matrix
        A_matrix = ( <Matrix_modn_dense> A._matrix )._matrix
        B_matrix = ( <Matrix_modn_dense> B._matrix )._matrix
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self_matrix[i+self._row][j+self._col] =  p + A_matrix[i+A._row][j+A._col] - B_matrix[i+B._row][j+B._col]
                if self_matrix[i+self._row][j+self._col] >= p:
                    # I really want to do in place operations here...
                    self_matrix[i+self._row][j+self._col] = self_matrix[i+self._row][j+self._col] - p

    cdef set_to_prod(self, MatrixWindow A, MatrixWindow B):
        # TODO: gather
        cdef Py_ssize_t i, j, k
        cdef mod_int p, s
        cdef mod_int** self_matrix
        cdef mod_int** A_matrix
        cdef mod_int** B_matrix
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        self_matrix = ( <Matrix_modn_dense> self._matrix )._matrix
        A_matrix = ( <Matrix_modn_dense> A._matrix )._matrix
        B_matrix = ( <Matrix_modn_dense> B._matrix )._matrix
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = 0
                for k from 0 <= k < A._ncols:
                    s = (s + A_matrix[i+A._row][k+A._col] * B_matrix[k+B._row][j+B._col]) % p
                self_matrix[i+self._row][j+self._col] = s % p

    cdef add_prod(self, MatrixWindow A, MatrixWindow B):
        # TODO: gather
        cdef Py_ssize_t i, j, k
        cdef mod_int p, s
        cdef mod_int** self_matrix
        cdef mod_int** A_matrix
        cdef mod_int** B_matrix
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        self_matrix = ( <Matrix_modn_dense> self._matrix )._matrix
        A_matrix = ( <Matrix_modn_dense> A._matrix )._matrix
        B_matrix = ( <Matrix_modn_dense> B._matrix )._matrix
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = self_matrix[i+self._row][j+self._col]
                for k from 0 <= k < A._ncols:
                    s = ( s + A_matrix[i+A._row][k+A._col] * B_matrix[k+B._row][j+B._col] ) % p
                self_matrix[i+self._row][j+self._col] = s % p

    cdef subtract_prod(self, MatrixWindow A, MatrixWindow B):
        # TODO: gather
        cdef Py_ssize_t i, j, k
        cdef mod_int p, s
        cdef mod_int** self_matrix
        cdef mod_int** A_matrix
        cdef mod_int** B_matrix
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        self_matrix = ( <Matrix_modn_dense> self._matrix )._matrix
        A_matrix = ( <Matrix_modn_dense> A._matrix )._matrix
        B_matrix = ( <Matrix_modn_dense> B._matrix )._matrix
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                s = self_matrix[i+self._row][j+self._col]
                for k from 0 <= k < A._ncols:
                    s = s + p - ( A_matrix[i+A._row][k+A._col] * B_matrix[k+B._row][j+B._col] ) % p
                self_matrix[i+self._row][j+self._col] = s % p

    cdef int element_is_zero(self, Py_ssize_t i, Py_ssize_t j):
        return (<Matrix_modn_dense>self._matrix)._matrix[i+self._row][j+self._col] == 0


