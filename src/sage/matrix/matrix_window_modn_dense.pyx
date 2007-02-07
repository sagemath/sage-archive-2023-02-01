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
        a = self._matrix.new_matrix(nrows, ncols)
#        a = Matrix_modn_dense.__new__(Matrix_modn_dense, False, ( <Matrix_modn_dense> self._matrix ).p, nrows, ncols)
        return self.new_matrix_window(a, 0, 0, nrows, ncols)

#        print "creating pseudo_matrix", nrows, ncols
#        cdef pseudo_matrix pm
#        pm = PY_NEW(pseudo_matrix)
#        pm.init(nrows, ncols, self._matrix.p, self._matrix.gather)

    cdef set_to(MatrixWindow_modn_dense self, MatrixWindow A):
        """
        Change self, making it equal A.
        """
        cdef Py_ssize_t i, j
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        cdef mod_int* self_row
        cdef mod_int* A_row
        for i from 0 <= i < self._nrows:
            self_row = ( <Matrix_modn_dense> self._matrix ).matrix[i + self._row] + self._col
            A_row    = ( <Matrix_modn_dense>    A._matrix ).matrix[i +    A._row] + A._col
            for j from 0 <= j < self._ncols:
                self_row[j] = A_row[j]
#            memcpy(self_rows[i+self._row] + self._col, A_rows[i+A._row] + A._col, self._ncols * sizeof(mod_int*))

    cdef set_to_zero(MatrixWindow_modn_dense self):
        cdef Py_ssize_t i, j
        cdef mod_int** rows
        rows = ( <Matrix_modn_dense> self._matrix ).matrix
        for i from self._row <= i < self._row + self._nrows:
            memset(rows[i] + self._col, 0, self._ncols * sizeof(mod_int*))

    cdef add(self, MatrixWindow A):
        cdef Py_ssize_t i, j
        cdef mod_int p
        cdef mod_int* self_row
        cdef mod_int* A_row
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            self_row = ( <Matrix_modn_dense> self._matrix ).matrix[i + self._row] + self._col
            A_row    = ( <Matrix_modn_dense>    A._matrix ).matrix[i +    A._row] + A._col
            for j from 0 <= j < self._ncols:
                self_row[j] += A_row[j]
                if self_row[j] >= p:
                    self_row[j] -= p

    cdef subtract(self, MatrixWindow A):
        cdef Py_ssize_t i, j
        cdef mod_int p
        cdef mod_int* self_row
        cdef mod_int* A_row
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            self_row = ( <Matrix_modn_dense> self._matrix ).matrix[i + self._row] + self._col
            A_row    = ( <Matrix_modn_dense>    A._matrix ).matrix[i +    A._row] + A._col
            for j from 0 <= j < self._ncols:
                if self_row[j] >= A_row[j]:
                    self_row[j] -= A_row[j]
                else:
                    self_row[j] += p - A_row[j]

    cdef set_to_sum(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        cdef mod_int p
        cdef mod_int* self_row
        cdef mod_int* A_row
        cdef mod_int* B_row
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        if self._nrows != B._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            self_row = ( <Matrix_modn_dense> self._matrix ).matrix[i + self._row] + self._col
            A_row    = ( <Matrix_modn_dense>    A._matrix ).matrix[i +    A._row] + A._col
            B_row    = ( <Matrix_modn_dense>    B._matrix ).matrix[i +    B._row] + B._col
            for j from 0 <= j < self._ncols:
                self_row[j] = A_row[j] + B_row[j]
                if self_row[j] >= p:
                    self_row[j] -= p

    cdef set_to_diff(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        cdef mod_int p
        cdef mod_int* self_row
        cdef mod_int* A_row
        cdef mod_int* B_row
        if self._nrows != A._nrows or self._ncols != A._ncols:
            raise ArithmeticError, "incompatible dimensions"
        if self._nrows != B._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        p = ( <Matrix_modn_dense> self._matrix ).p
        for i from 0 <= i < self._nrows:
            self_row = ( <Matrix_modn_dense> self._matrix ).matrix[i + self._row] + self._col
            A_row    = ( <Matrix_modn_dense>    A._matrix ).matrix[i +    A._row] + A._col
            B_row    = ( <Matrix_modn_dense>    B._matrix ).matrix[i +    B._row] + B._col
            for j from 0 <= j < self._ncols:
                if A_row[j] > B_row[j]:
                    self_row[j] = A_row[j] - B_row[j]
                else:
                    self_row[j] = p + A_row[j] - B_row[j]

    cdef set_to_prod(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k, gather, top, A_ncols
        cdef mod_int p, s
        cdef mod_int* self_row
        cdef mod_int* A_row
        cdef mod_int** B_matrix_off
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        B_matrix_off = ( <Matrix_modn_dense> B._matrix ).matrix + B._row
        p = ( <Matrix_modn_dense> self._matrix ).p
        gather = ( <Matrix_modn_dense> self._matrix ).gather
        A_ncols = A._ncols

        if gather <= 1:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix ).matrix[i + self._row] + self._col
                A_row    = ( <Matrix_modn_dense>    A._matrix ).matrix[i +    A._row] + A._col
                for j from 0 <= j < B._ncols:
                    s = 0
                    for k from 0 <= k < A._ncols:
                        s = (s + A_row[k] * B_matrix_off[k][j+B._col]) % p
                    self_row[j] = s
        else:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix ).matrix[i + self._row] + self._col
                A_row    = ( <Matrix_modn_dense>    A._matrix ).matrix[i +    A._row] + A._col
                for j from 0 <= j < B._ncols:
                    s = 0
                    k = 0
                    while k < A_ncols:
                        top = k + gather
                        if top > A_ncols:
                            top = A_ncols
                        for k from k <= k < top: # = min(k+gather, A._ncols)
                            pass
                            s += A_row[k] * B_matrix_off[k][j+B._col]
                        s %= p
                    self_row[j] = s

    cdef add_prod(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k, gather, top, A_ncols
        cdef mod_int p, s
        cdef mod_int* self_row
        cdef mod_int* A_row
        cdef mod_int** B_matrix_off
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        B_matrix_off = ( <Matrix_modn_dense> B._matrix ).matrix + B._row
        p = ( <Matrix_modn_dense> self._matrix ).p
        gather = ( <Matrix_modn_dense> self._matrix ).gather
        A_ncols = A._ncols

        if gather <= 1:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix ).matrix[i+self._row] + self._col
                A_row = ( <Matrix_modn_dense> A._matrix ).matrix[i+A._row] + A._col
                for j from 0 <= j < B._ncols:
                    s = self_row[j]
                    for k from 0 <= k < A._ncols:
                        s = ( s + A_row[k] * B_matrix_off[k][j+B._col] ) % p
                    self_row[j] = s

        else:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix ).matrix[i+self._row] + self._col
                A_row = ( <Matrix_modn_dense> A._matrix ).matrix[i+A._row] + A._col
                for j from 0 <= j < B._ncols:
                    s = self_row[j]
                    k = 0
                    while k < A_ncols:
                        top = k + gather
                        if top > A_ncols:
                            top = A_ncols
                        for k from k <= k < top: # = min(k+gather, A._ncols)
                            s += A_row[k] * B_matrix_off[k][j+B._col]
                        s %= p
                    self_row[j] = s


    cdef subtract_prod(self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k, gather, top, A_ncols
        cdef mod_int p, s, p2
        cdef mod_int* self_row
        cdef mod_int* A_row
        cdef mod_int** B_matrix_off
        if A._ncols != B._nrows or self._nrows != A._nrows or self._ncols != B._ncols:
            raise ArithmeticError, "incompatible dimensions"
        B_matrix_off = ( <Matrix_modn_dense> B._matrix ).matrix + B._row
        p = ( <Matrix_modn_dense> self._matrix ).p
        gather = ( <Matrix_modn_dense> self._matrix ).gather
        A_ncols = A._ncols
        p2 = p*p

        if gather <= 1:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix ).matrix[i+self._row] + self._col
                A_row = ( <Matrix_modn_dense> A._matrix ).matrix[i+A._row] + A._col
                for j from 0 <= j < B._ncols:
                    s = self_row[j]
                    for k from 0 <= k < A._ncols:
                        s = ( s + p2 - A_row[k] * B_matrix_off[k][j+B._col] ) % p
                    self_row[j] = s
        else:
            for i from 0 <= i < A._nrows:
                self_row = ( <Matrix_modn_dense> self._matrix ).matrix[i+self._row] + self._col
                A_row = ( <Matrix_modn_dense> A._matrix ).matrix[i+A._row] + A._col
                for j from 0 <= j < B._ncols:
                    s = self_row[j] + gather * p2 # using unsigned ints
                    k = 0
                    while k < A_ncols:
                        top = k + gather
                        if top > A_ncols:
                            top = A_ncols
                        for k from k <= k < top: # = min(k+gather, A._ncols)
                            s -= A_row[k] * B_matrix_off[k][j+B._col]
                        s = s % p + gather * p2 # using unsigned ints
                    self_row[j] = s % p


    cdef int element_is_zero(self, Py_ssize_t i, Py_ssize_t j):
        return (<Matrix_modn_dense>self._matrix).matrix[i+self._row][j+self._col] == 0
