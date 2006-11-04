cimport matrix

include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'

cdef class Matrix_dense(matrix.Matrix):
    def __copy__(self):
        """
        Return a copy of this matrix.  Changing the entries of the
        copy will not change the entries of this matrix.

        EXAMPLES:
            sage: ???
        """
        return self.new_matrix(entries=self.list(), coerce=False, copy=False)

    def __hash__(self):
        """
        Return the hash of this matrix.

        Equal matrices should have equal hashes, even if one is sparse and
        the other is dense.

        EXAMPLES:
            sage: m = matrix(2, range(24))
            sage: m.set_immutable()
            sage: hash(m)

            sage: d = M.dense_matrix()
            sage: d.set_immutable()
            sage: hash(d)

            sage: hash(m) == hash(d)
            True
        """
        x = self.fetch('hash')
        if not x is None: return x

        if not self._mutability._is_immutable:
            raise TypeError, "mutable matrices are unhashable"

        v = self._list()
        cdef Py_ssize_t i
        cdef long h
        h = 0
        n = 1
        cdef PyObject** w
        w = FAST_SEQ_UNSAFE(v)
        for i from 0 <= i < len(v):
            h = h ^ (i * PyObject_Hash( <object> w[i] ))

        self.cache('hash', h)
        return h

    def _multiply_classical(left, Matrix_dense right):
        """
        Multiply the matrices left and right using the classical $O(n^3)$
        algorithm.

        This method assumes that left and right have the same parent and
        compatable dimensions.

        EXAMPLES
            sage: include the 0 rows and 0 columns cases
        """
        cdef Py_ssize_t i, j, k, l
        if left._ncols != right._nrows:
            raise IndexError, "Number of columns of left must equal number of rows of other."


        v = PyList_New(left._nrows * right._ncols)     # this is really sort of v = []..."
        zero = left.base_ring()(0)
        l = 0
        for i from 0 <= i < left._nrows:
            for j from 0 <= j < right._ncols:
                s = zero
                for k from 0 <= k < left._ncols:
                    s = s + left.get_unsafe(i,k) * right.get_unsafe(k,j)
                # This is really v.append(s)
                Py_INCREF(s); PyList_SET_ITEM(v, l, s)
                l = l + 1
        return left.new_matrix(left._nrows, right._ncols, entries = v, coerce=False, copy=False)


    cdef _pickle(self):
        version = -1
        data = self._list()  # linear list of all elements
        return data, version

    cdef _unpickle(self, data, int version):
        cdef Py_ssize_t i, j, k
        if version == -1:
            # data is a *list* of the entries of the matrix.
            # TODO: Change the data[k] below to use the fast list access macros from the Python/C API
            k = 0
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    print "setting %s, %s entry to %s"%(i,j,data[k])
                    self.set_unsafe(i, j, data[k])
                    k = k + 1
        else:
            raise RuntimeError, "unknown matrix version (=%s)"%version

    cdef int _cmp_c_impl(self, matrix.Matrix right) except -2:
        return cmp(self._list(), right._list())
