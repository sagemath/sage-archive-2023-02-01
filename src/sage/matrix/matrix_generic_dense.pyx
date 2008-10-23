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

            zero = parent.base_ring()(0)
            self._entries = [zero]*(self._nrows*self._ncols)

            if x != zero:
                if self._nrows != self._ncols:
                    raise TypeError, "nonzero scalar matrix must be square"
                for i from 0 <= i < self._nrows:
                    self._entries[i*self._ncols + i] = x

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        Py_DECREF(<object>PyList_GET_ITEM(self._entries, i*self._ncols + j))
        Py_INCREF(value)
        PyList_SET_ITEM(self._entries, i*self._ncols + j, value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        return <object>PyList_GET_ITEM(self._entries, i*self._ncols + j)

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
    #    * cdef _add_
    #    * cdef _mul_
    #    * cdef _cmp_c_impl
    #    * __neg__
    #    * __invert__
    # x  * __copy__
    # x  * _multiply_classical
    # x  * _list -- copy of the list of underlying elements
    #    * _dict -- copy of the sparse dictionary of underlying elements
    ########################################################################

    def __copy__(self):
        """
        Creates a copy of self, which may be changed without altering self.

        EXAMPLES:
            sage: A = matrix(ZZ[['t']], 2,3,range(6)); A
            [0 1 2]
            [3 4 5]
            sage: A.subdivide(1,1); A
            [0|1 2]
            [-+---]
            [3|4 5]
            sage: B = A.copy(); B
            [0|1 2]
            [-+---]
            [3|4 5]
            sage: B == A
            True
            sage: B[0,0] = 100
            sage: B
            [100|  1   2]
            [---+-------]
            [  3|  4   5]
            sage: A
            [0|1 2]
            [-+---]
            [3|4 5]
        """
        A = self.__class__(self._parent, self._entries, copy = True, coerce=False)
        if self.subdivisions is not None:
            A.subdivide(*self.get_subdivisions())
        return A

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
    #    * cdef _sub_
    # x  * __deepcopy__
    #    * __invert__
    #    * _multiply_classical
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    ########################################################################

    def __deepcopy__(self):
        import copy
        return self.__class__(self._parent, copy.deepcopy(self._entries), copy = False, coerce=False)

