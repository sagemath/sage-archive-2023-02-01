cimport matrix

include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'

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

        v = self._dict()
        cdef long i, h
        h = 0
        for ij, x in v.iteritems():
            # todo -- optimize with api
            #i = PyInt_AS_LONG(ij[0]) + PyInt_AS_LONG(ij[1])*self._ncols
            #i = PyInt_AS_LONG(w[0]) + PyInt_AS_LONG(w[1])*self._ncols
            #print ij
            i = ij[0] + ij[1]*self._ncols
            h = h ^ (i*PyObject_Hash(x))

        self.cache('hash', h)
        return h

    def x__hash__(self):
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

        v = self._dict()
        cdef long i, h
        h = 0
        for ij, x in v.iteritems():
            i = PyInt_AS_LONG(<object>PyTuple_GET_ITEM(ij,0)) + \
                PyInt_AS_LONG(<object>PyTuple_GET_ITEM(ij,1)) * self._ncols
            h = h ^ (i*PyObject_Hash(x))

        self.cache('hash', h)
        return h
