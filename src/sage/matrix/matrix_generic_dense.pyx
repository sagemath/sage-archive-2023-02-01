"""
Dense Matrices over a general ring
"""

def _convert_dense_entries_to_list(entries):
    # Create matrix from a list of vectors
    e = []
    for v in entries:
        e = e+ v.list()
    copy = False
    return e

include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/python_list.pxi"
include "../ext/python_number.pxi"
include "../ext/python_ref.pxi"

cimport matrix_dense
import matrix_dense

cimport matrix
cimport matrix_window

cdef class Matrix_generic_dense(matrix_dense.Matrix_dense):
    r"""
    The \class{Matrix_generic_dense} class derives from \class{Matrix}, and
    defines functionality for dense matrices over any base ring.
    Matrices are represented by a list of elements in the base ring,
    and element access operations are implemented in this class.
    """
    ########################################################################
    # LEVEL 1 functionality
    # 0 * __new__   (not needed)
    # x * __init__
    # 0 * __dealloc__   (not needed)
    # x * set_unsafe
    # x * get_unsafe
    # x * def _pickle
    # x * def _unpickle
    ########################################################################
    def __init__(self, parent, entries, copy, coerce):
        matrix.Matrix.__init__(self, parent)

        cdef Py_ssize_t i, n

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
                self._entries = entries
            else:
                self._entries = [None]*(self._nrows*self._ncols)
                n = len(entries)
                if coerce:
                    R = parent.base_ring()
                    for i from 0 <= i < n:
                        self._entries[i] = R(entries[i])
                else:
                    for i from 0 <= i < n:
                        self._entries[i] = entries[i]

        else:

            self._entries = [None]*(self._nrows*self._ncols)
            zero = parent.base_ring()(0)
            for i from 0 <= i < self._nrows * self._ncols:
                self._entries[i] = zero

            if x != zero:
                if self._nrows != self._ncols:
                    raise TypeError, "nonzero scalar matrix must be square"
                for i from 0 <= i < self._nrows:
                    self._entries[i*self._ncols + i] = x

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        # TODO: make faster with Python/C API
        self._entries[i*self._ncols + j] = value

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        # TODO: make faster with Python/C API
        return self._entries[i*self._ncols + j]

    def _pickle(self):
        return self._entries, 0

    def _unpickle(self, data, int version):
        if version == 0:
            self._entries = data
        else:
            raise RuntimeError, "unknown matrix version"

    def __richcmp__(matrix.Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)
    def __hash__(self):
        return self._hash()

    ########################################################################
    # LEVEL 2 functionality
    #    * cdef _add_c_impl
    #    * cdef _mul_c_impl
    #    * cdef _cmp_c_impl
    #    * __neg__
    #    * __invert__
    # x  * __copy__
    # x  * _multiply_classical
    # x  * _list -- copy of the list of underlying elements
    #    * _dict -- copy of the sparse dictionary of underlying elements
    ########################################################################

    def __copy__(self):
        return self.__class__(self._parent, self._entries, copy = True, coerce=False)

    def _multiply_classical(left, matrix.Matrix _right):
        """
        Multiply the matrices left and right using the classical $O(n^3)$
        algorithm.

        EXAMPLES:
        We multiply two matrices over a fairly general ring:

            sage: R.<x,y> = Integers(8)['x,y']
            sage: a = matrix(R,2,[x,y,x^2,y^2]); a
            [  x   y]
            [x^2 y^2]
            sage: a*a
            [  x^2*y + x^2     y^3 + x*y]
            [x^2*y^2 + x^3   y^4 + x^2*y]
            sage: a.det()^2 == (a*a).det()
            True

            sage: A = matrix(QQ['x,y'], 2, [0,-1,2,-2])
            sage: B = matrix(QQ['x,y'], 2, [-1,-1,-2,-2])
            sage: A*B
            [2 2]
            [2 2]

        SAGE fully supports degenerate matrices with 0 rows or 0 columns:
            sage: A = matrix(QQ['x,y'], 0, 4, []); A
            []
            sage: B = matrix(QQ['x,y'], 4,0, []); B
            []
            sage: A*B
            []
            sage: B*A
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
        """
        cdef Py_ssize_t i, j, k, m, nr, nc, snc, p
        cdef object v
        cdef Matrix_generic_dense A, right
        right = _right

        if left._ncols != right._nrows:
            raise IndexError, "Number of columns of left must equal number of rows of other."

        nr = left._nrows
        nc = right._ncols
        snc = left._ncols

        R = left.base_ring()
        P = left.matrix_space(nr, nc)
        v = PyList_New(left._nrows * right._ncols)
        zero = R(0)
        p = 0
        cdef PyObject *l, *r
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                z = zero
                m = i*snc
                for k from 0 <= k < snc:
                    # The following is really:
                    #     z = z + left._entries[m + k] * right._entries[k*right._ncols + j]
                    l = PyList_GET_ITEM(left._entries, m+k)
                    r = PyList_GET_ITEM(right._entries, k*nc + j)
                    z = z + PyNumber_Multiply(<object>l, <object>r)
                Py_INCREF(z); PyList_SET_ITEM(v, p, z)   #   Basically this is "v.append(z)"
                p = p + 1

        A = Matrix_generic_dense.__new__(Matrix_generic_dense, 0, 0 ,0)
        matrix.Matrix.__init__(A, P)
        A._entries = v
        return A

    def _list(self):
        return self._entries

    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_c_impl
    # x  * __deepcopy__
    #    * __invert__
    #    * _multiply_classical
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    ########################################################################

    def __deepcopy__(self):
        import copy
        return self.__class__(self._parent, copy.deepcopy(self._entries), copy = False, coerce=False)

    def matrix_window_c(self, int row=0, int col=0, int nrows=-1, int ncols=-1):
        if nrows == -1:
            nrows = self._nrows
            ncols = self._ncols
        return MatrixWindow(self, row, col, nrows, ncols)



#######################################################################

cdef class MatrixWindow(matrix_window.MatrixWindow):

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
        v = self._matrix._entries
        w = [None]*(self._nrows*self._ncols)
        cdef Py_ssize_t i, j, k, l
        k = 0
        for i from 0 <= i < self._nrows:
            l = i * self._matrix._ncols
            for j from 0 <= j < self._ncols:
                w[k] = v[l + j]
                k = k + 1
        return w

    def matrix_window(MatrixWindow self, Py_ssize_t row, Py_ssize_t col, Py_ssize_t n_rows, Py_ssize_t n_cols):
        """
        Returns a matrix window relative to this window of the underlying matrix.
        """
        return self._matrix.matrix_window(self._row + row, self._col + col, n_rows, n_cols)

    def nrows(MatrixWindow self):
        return self._nrows

    def ncols(MatrixWindow self):
        return self._ncols

    def set_to(MatrixWindow self, MatrixWindow A):
        cdef Py_ssize_t i, j
        cdef Py_ssize_t start, self_ix
        cdef Py_ssize_t A_ix
        for i from 0 <= i < self._nrows:
            A_ix = self._matrix._ncols * (i+A._row) + A._col
            for self_ix from (i+self._row)*self._ncols + self._col <= self_ix < start + self._ncols:
                self._matrix._entries[self_ix] = A._matrix._entries[A_ix]
                A_ix = A_ix + 1


    def set_to_zero(MatrixWindow self):
        cdef Py_ssize_t i, j
        cdef Py_ssize_t start, self_ix
        zero = self._matrix.base_ring(0)
        for i from 0 <= i < self._nrows:
            for self_ix from (i+self._row)*self._ncols + self._col <= self_ix < start + self._ncols:
                self._matrix._entries[self_ix] = zero


    def add(MatrixWindow self, MatrixWindow A):
        cdef Py_ssize_t i, j
        cdef Py_ssize_t start, self_ix
        cdef Py_ssize_t A_ix
        for i from 0 <= i < self._nrows:
            A_ix = (i+A._row)*self._matrix._ncols + A._col
            for self_ix from self._ncols*(i+self._row) + self._col <= self_ix < start + self._ncols:
                self._matrix._entries[self_ix] = slef._matrix._entries[self_ix] + A._matrix._entries[A_ix]
                A_ix = A_ix + 1


    def subtract(MatrixWindow self, MatrixWindow A):
        cdef Py_ssize_t i, j
        cdef Py_ssize_t start, self_ix
        cdef Py_ssize_t A_ix
        for i from 0 <= i < self._nrows:
            A_ix = A._matrix._ncols*(i+A._row) + a._col
            for self_ix from self._ncols*(i+self._row) + self._col <= self_ix < start + self._ncols:
                self._matrix._entries[self_ix] = self._matrix._entries[self_ix] - A._matrix._entries[A_ix]
                A_ix = A_ix + 1


    def set_to_sum(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        cdef Py_ssize_t start, self_ix
        cdef Py_ssize_t A_ix
        for i from 0 <= i < self._nrows:
            A_ix = A._matrix._ncols*(i+A._row) + A._col
            B_ix = B._matrix._ncols*(i+B._row) + B._col
            for self_ix from self._ncols*(i+self._row) + self._col <= self_ix < start + self._ncols:
                self._matrix._entries[self_ix] = A._matrix._entries[A_ix] + B._matrix._entries[B_ix]
                A_ix = A_ix + 1
                B_ix = B_ix + 1


    def set_to_diff(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j
        cdef Py_ssize_t start, self_ix
        cdef Py_ssize_t A_ix
        for i from 0 <= i < self._nrows:
            A_ix = A._matrix._ncols*(i+A._row) + A._col
            B_ix = B._matrix._ncols*(i+B._row) + B._col
            for self_ix from self._ncols*(i+self._row) + self._col <= self_ix < start + self._ncols:
                self._matrix._entries[self_ix] = A._matrix._entries[A_ix] - B._matrix._entries[B_ix]
                A_ix = A_ix + 1
                B_ix = B_ix + 1


    def set_to_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k
        cdef Py_ssize_t start, self_ix
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                sum = A._matrix._entries[ A._ncols*(A._row+i) + A._col ] *B._matrix._entries[ B._ncols*(B._row)+B._col+j ]
                for k from 1 <= k < A._ncols:
                    sum = sum + A._matrix._entries[ A._ncols*(A._row+i) + A._col + k ] * B._matrix._entries[ B._ncols*(B._row+k)+B._col+j ]
                self._matrix._entries[ self._ncols*(self._row+i) + self._col+j ] = sum


    def add_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef Py_ssize_t i, j, k
        cdef Py_ssize_t start, self_ix
        for i from 0 <= i < A._nrows:
            for j from 0 <= j < B._ncols:
                sum = A._matrix._entries[ A._ncols*(A._row+i)+A._col ] *B._matrix._entries[ B._ncols* B._row  + B._col+j ]
                for k from 1 <= k < A._ncols:
                    sum = sum + A._matrix._entries[ A._ncols*(A._row+i) + A._col + k ] * B._matrix._entries[ B._ncols*(B._row+k)+B._col+j ]
                self._matrix._entries[ self._ncols*(self._row+i)+self._col+j ] = self._matrix._entries[ self._ncols*(self._row+i)+self._col+j ] + sum


    def new_empty_window(MatrixWindow self, Py_ssize_t nrows, Py_ssize_t ncols):
        return self._matrix.new_matrix(nrows,ncols).matrix_window()
