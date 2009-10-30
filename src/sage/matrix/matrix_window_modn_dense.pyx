# cython: cdivision=True

"""
TESTS:
    sage: a = random_matrix(GF(11), 30, 40)
    sage: b = random_matrix(GF(11), 40, 53)
    sage: a._multiply_strassen(b, 7) == a.change_ring(ZZ)*b.change_ring(ZZ)
    True
"""


include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"

import matrix_window
cimport matrix_window

from sage.matrix.matrix_modn_dense import Matrix_modn_dense
from sage.matrix.matrix_modn_dense cimport Matrix_modn_dense


cdef class MatrixWindow_modn_dense(matrix_window.MatrixWindow):

    cpdef new_empty_window(self, Py_ssize_t nrows, Py_ssize_t ncols):
# the current code is all python, goes through inits, and possibly creates a parent...
# can we get away with something faster?
#       a = MatrixWindow_modn_dense.__new__(self._parent, )?
       a = self._matrix.new_matrix(nrows, ncols)
       return self.new_matrix_window(a, 0, 0, nrows, ncols)

    cpdef set_to(MatrixWindow_modn_dense self, MatrixWindow A):
        """
        Change self, making it equal A.
        """
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        cdef mod_int** self_rows
        cdef mod_int** A_rows
        self_rows = ( <Matrix_modn_dense> self._matrix )._matrix
        A_rows    = ( <Matrix_modn_dense> A._matrix )._matrix
        for i from 0 <= i < self._nrows:
            memcpy(self_rows[i+self._row] + self._col, A_rows[i+A._row] + A._col, self._ncols * sizeof(mod_int))

    cpdef set_to_zero(MatrixWindow_modn_dense self):
        cdef Py_ssize_t i, j
        cdef mod_int** rows
        rows = ( <Matrix_modn_dense> self._matrix )._matrix
        for i from self._row <= i < self._row + self._nrows:
            memset(rows[i] + self._col, 0, self._ncols * sizeof(mod_int))

    cpdef add(self, MatrixWindow A):
        cdef Py_ssize_t i, j
        cdef mod_int k, p
        cdef mod_int* self_row
        cdef mod_int* A_row
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            self_row = ( <Matrix_modn_dense> self._matrix )._matrix[i + self._row] + self._col
            A_row    = ( <Matrix_modn_dense>    A._matrix )._matrix[i +    A._row] + A._col
            for j from 0 <= j < self._ncols:
                k = self_row[j] + A_row[j]
                self_row[j] = k - (k >= p) * p

    cpdef subtract(self, MatrixWindow A):
        cdef Py_ssize_t i, j
        cdef mod_int k, p
        cdef mod_int* self_row
        cdef mod_int* A_row
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            self_row = ( <Matrix_modn_dense> self._matrix )._matrix[i + self._row] + self._col
            A_row    = ( <Matrix_modn_dense>    A._matrix )._matrix[i +    A._row] + A._col
            for j from 0 <= j < self._ncols:
                k = p + self_row[j] - A_row[j]
                self_row[j] = k - (k >= p) * p

    cpdef set_to_sum(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        cdef mod_int k, p
        cdef mod_int* self_row
        cdef mod_int* A_row
        cdef mod_int* B_row
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        if self._nrows != B._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            self_row = ( <Matrix_modn_dense> self._matrix )._matrix[i + self._row] + self._col
            A_row    = ( <Matrix_modn_dense>    A._matrix )._matrix[i +    A._row] + A._col
            B_row    = ( <Matrix_modn_dense>    B._matrix )._matrix[i +    B._row] + B._col
            for j from 0 <= j < self._ncols:
                k = A_row[j] + B_row[j]
                self_row[j] = k - (k >= p) * p

    cpdef set_to_diff(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        cdef mod_int k, p
        cdef mod_int* self_row
        cdef mod_int* A_row
        cdef mod_int* B_row
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        if self._nrows != B._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            self_row = ( <Matrix_modn_dense> self._matrix )._matrix[i + self._row] + self._col
            A_row    = ( <Matrix_modn_dense>    A._matrix )._matrix[i +    A._row] + A._col
            B_row    = ( <Matrix_modn_dense>    B._matrix )._matrix[i +    B._row] + B._col
            for j from 0 <= j < self._ncols:
                k = p + A_row[j] - B_row[j]
                self_row[j] = k - (k >= p) * p

    cpdef set_to_prod(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k, gather, top, A_ncols
        cdef mod_int p, s
        cdef mod_int* self_row
        cdef mod_int* A_row
        cdef mod_int** B_matrix_off
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        B_matrix_off = ( <Matrix_modn_dense> B._matrix )._matrix + B._row
        p = ( <Matrix_modn_dense> self._matrix ).p
        gather = ( <Matrix_modn_dense> self._matrix ).gather
        A_ncols = A._ncols

        cdef mod_int A_i_k
        cdef mod_int* B_row_k

        if gather <= 1:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix )._matrix[i + self._row] + self._col
                A_row    = ( <Matrix_modn_dense>    A._matrix )._matrix[i +    A._row] + A._col
                k = 0
                for j from 0 <= j < B._ncols:
                    self_row[j] = (A_row[k] * B_matrix_off[k][B._col+j]) % p
                for k from 1 <= k < A._ncols:
                    A_i_k = A_row[k]
                    B_row_k = B_matrix_off[k] + B._col
                    for j from 0 <= j < B._ncols:
                        self_row[j] = (self_row[j] + A_i_k * B_row_k[j]) % p

        else:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix )._matrix[i + self._row] + self._col
                A_row    = ( <Matrix_modn_dense>    A._matrix )._matrix[i +    A._row] + A._col
                for j from 0 <= j < B._ncols:
                    self_row[j] = 0
                k = 0
                while k < A_ncols:
                    top = k + gather
                    if top > A_ncols:
                        top = A_ncols
                    for k from k <= k < top: # = min(k+gather, A._ncols)
                        A_i_k = A_row[k]
                        B_row_k = B_matrix_off[k] + B._col
                        for j from 0 <= j < B._ncols:
                            self_row[j] += A_i_k * B_row_k[j]
                    for j from 0 <= j < B._ncols:
                        self_row[j] %= p

    cpdef add_prod(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k, gather, top, A_ncols
        cdef mod_int p, s
        cdef mod_int* self_row
        cdef mod_int* A_row
        cdef mod_int** B_matrix_off
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        B_matrix_off = ( <Matrix_modn_dense> B._matrix )._matrix + B._row
        p = ( <Matrix_modn_dense> self._matrix ).p
        gather = ( <Matrix_modn_dense> self._matrix ).gather
        A_ncols = A._ncols

        cdef mod_int A_i_k
        cdef mod_int* B_row_k

        if gather <= 1:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix )._matrix[i+self._row] + self._col
                A_row = ( <Matrix_modn_dense> A._matrix )._matrix[i+A._row] + A._col
                for k from 0 <= k < A._ncols:
                    A_i_k = A_row[k]
                    B_row_k = B_matrix_off[k] + B._col
                    for j from 0 <= j < B._ncols:
                        self_row[j] = (self_row[j] + A_i_k * B_row_k[j]) % p

        else:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix )._matrix[i+self._row] + self._col
                A_row = ( <Matrix_modn_dense> A._matrix )._matrix[i+A._row] + A._col
                k = 0
                while k < A_ncols:
                    top = k + gather
                    if top > A_ncols:
                        top = A_ncols
                    for k from k <= k < top: # = min(k+gather, A._ncols)
                        A_i_k = A_row[k]
                        B_row_k = B_matrix_off[k] + B._col
                        for j from 0 <= j < B._ncols:
                            self_row[j] += A_i_k * B_row_k[j]
                    for j from 0 <= j < B._ncols:
                        self_row[j] %= p


    cpdef subtract_prod(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k, gather, top, A_ncols
        cdef mod_int p, s, p2
        cdef mod_int* self_row
        cdef mod_int* A_row
        cdef mod_int** B_matrix_off
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        B_matrix_off = ( <Matrix_modn_dense> B._matrix )._matrix + B._row
        p = ( <Matrix_modn_dense> self._matrix ).p
        gather = ( <Matrix_modn_dense> self._matrix ).gather
        A_ncols = A._ncols
        p2 = p*(p-1)

        cdef mod_int A_i_k
        cdef mod_int* B_row_k

        if gather <= 1:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix )._matrix[i+self._row] + self._col
                A_row = ( <Matrix_modn_dense> A._matrix )._matrix[i+A._row] + A._col
                for k from 0 <= k < A._ncols:
                    A_i_k = A_row[k]
                    B_row_k = B_matrix_off[k] + B._col
                    for j from 0 <= j < B._ncols:
                        self_row[j] = ( self_row[j] + p2 - A_i_k * B_row_k[j] ) % p

        else:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix )._matrix[i+self._row] + self._col
                A_row = ( <Matrix_modn_dense> A._matrix )._matrix[i+A._row] + A._col
                k = 0
                while k < A_ncols:
                    top = k + gather
                    if top > A_ncols:
                        top = A_ncols
                    for k from k <= k < top: # = min(k+gather, A._ncols)
                        A_i_k = A_row[k]
                        B_row_k = B_matrix_off[k] + B._col
                        for j from 0 <= j < B._ncols:
                            self_row[j] += p2 - A_i_k * B_row_k[j]
                    for j from 0 <= j < B._ncols:
                        self_row[j] %= p


    cpdef bint element_is_zero(self, Py_ssize_t i, Py_ssize_t j):
        return (<Matrix_modn_dense>self._matrix)._matrix[i+self._row][j+self._col] == 0


