cimport matrix

include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'
include '../ext/python.pxi'
include '../ext/interrupt.pxi'

cdef extern from "Python.h":
    PyObject* PyTuple_GET_ITEM0 "PyTuple_GET_ITEM" (PyObject*  p, Py_ssize_t pos)
    PyObject* PyList_GET_ITEM0 "PyList_GET_ITEM" (PyObject*  list, Py_ssize_t i)
    Py_ssize_t PyNumber_AsSsize_t(PyObject* o, PyObject* exc)



cdef class Matrix_sparse(matrix.Matrix):
    def __copy__(self):
        """
        Return a copy of this matrix.  Changing the entries of the
        copy will not change the entries of this matrix.

        EXAMPLES:
            sage: ???
        """
        return self.new_matrix(entries=self.dict(), coerce=False, copy=False)

    def __hash__(self):
        """
        Return the hash of this matrix.

        Equal matrices should have equal hashes, even if one is sparse and
        the other is dense.

        EXAMPLES:
            sage: m = matrix(2, range(6))
            sage: m.set_immutable()
            sage: hash(m)
            5
            sage: d = M.dense_matrix()
            sage: d.set_immutable()
            sage: hash(d)
            ?

            sage: hash(m) == hash(d)


            sage: A = Matrix(ZZ[['t']], 2, 2, range(4))
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()
            sage: B = A.copy(); B.set_immutable()
            sage: hash(A) == hash(B)
            True
        """
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

        self.cache('hash', h)
        return h

    def _multiply_classical(Matrix_sparse left, Matrix_sparse right):
        """
        sage: do 0 x 0 case
        """
        cdef Py_ssize_t row, col, row_start, k1, k2, len_left, len_right
        left_nonzero = left.nonzero_positions(copy=False, column_order=False)
        right_nonzero = left.nonzero_positions(copy=False, column_order=True)
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
                    if get_ij(left_nonzero,k1,1) == get_ij(right_nonzero,k2,0):
                        if sum is None:
                            sum = left[left_nonzero[k1]] * right[right_nonzero[k2]]
                        else:
                            sum = sum + left[left_nonzero[k1]] * right[right_nonzero[k2]]
                        k1 = k1 + 1
                        k2 = k2 + 1
                    elif get_ij(left_nonzero,k1,1) < get_ij(right_nonzero,k2,0):
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

    cdef int _will_use_strassen(self, matrix.Matrix right) except -1:
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

    cdef int _cmp_c_impl(self, matrix.Matrix right) except -2:
        return cmp(self._dict(), right._dict())

##################################
# Helper code

cdef Py_ssize_t get_ij(v, Py_ssize_t i, int j):
    return PyNumber_AsSsize_t(PyTuple_GET_ITEM0(PyList_GET_ITEM0(<PyObject*>v, i), j), NULL)

#cdef Py_ssize_t get_ij(v, Py_ssize_t i, int j):
#    return PyNumber_AsSsize_t(<object>PyTuple_GET_ITEM(
#              <object>PyList_GET_ITEM(v, i), j), <object>NULL)
