cimport matrix

include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'

include '../ext/interrupt.pxi'

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
        cdef Py_ssize_t row, col, row_start, k1, k2
        left_nonzero = left.nonzero_positions(copy=False, column_order=False)
        right_nonzero = left.nonzero_positions(copy=False, column_order=True)

        e = {}
        k1 = 0
        _sig_on
        while k1 < len(left_nonzero):
            row_start = k1
            row = left_nonzero[row_start][0]
            k2 = 0
            while k2 < len(right_nonzero):
                col = right_nonzero[k2][1]
                sum = 0
                k1 = row_start
                while k1 < len(left_nonzero) and left_nonzero[k1][0] == row and k2 < len(right_nonzero) and right_nonzero[k2][1] == col:
                    if left_nonzero[k1][1] == right_nonzero[k2][0]:
                        sum = sum + left[left_nonzero[k1]] * right[right_nonzero[k2]]
                        k1 = k1 + 1
                        k2 = k2 + 1
                    elif left_nonzero[k1][1] < right_nonzero[k2][0]:
                        k1 = k1 + 1
                    else:
                        k2 = k2 + 1

                if sum != 0:
                    e[row, col] = sum

                while k2 < len(right_nonzero) and right_nonzero[k2][1] == col:
                    k2 = k2 + 1

            while k1 < len(left_nonzero) and left_nonzero[k1][0] == row:
                k1 = k1 + 1

        _sig_off
        return left.new_matrix(left._nrows, right._ncols, entries=e, coerce=False, copy=False)


