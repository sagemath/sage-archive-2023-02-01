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
        EXAMPLES:
            sage: A = Matrix(ZZ[['t']], 2, 2, range(4))
            sage: B = A.copy()
            sage: A.__hash__() == B.__hash__()
            True
            sage: A[0,0] = -1
            sage: A.__hash__() == B.__hash__()
            False
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

