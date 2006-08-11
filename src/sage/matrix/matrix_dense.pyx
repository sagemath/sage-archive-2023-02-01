#############################################
## Generic *DENSE* matrices over any field
#############################################

def _convert_dense_entries_to_list(entries):
    # Create matrix from a list of vectors
    e = []
    for v in entries:
        e = e+ v.list()
    copy = False
    return e

include "../ext/interrupt.pxi"

cimport matrix_pyx
import matrix_pyx



cdef class Matrix_dense(matrix_pyx.Matrix):
    """
    The \\class{Matrix_dense} class derives from
    \\class{Matrix}, and defines functionality for dense matrices over
    any base ring.  Matrices are represented by a list of elements in
    the base ring, and element access operations are implemented in
    this class.
    """
    def __new__(self, parent, int nrows, int ncols,
                entries = None,
                coerce = True):
        self._nrows = nrows
        self._ncols = ncols

    def __init__(self, parent, int nrows, int ncols,
                 entries = None,
                 coerce = True):
        matrix_pyx.Matrix.__init__(self, parent)
        cdef int i, n
        #self._entries = PyList_New(nrows*ncols)
        self._entries = [None]*(nrows*ncols)
        if entries:
            n = len(entries)
            if coerce:
                R = parent.base_ring()
                for i from 0 <= i < n:
                    self._entries[i] = R(entries[i])
            else:
                for i from 0 <= i < n:
                    self._entries[i] = entries[i]

        self._row_indices = <int*> PyMem_Malloc(sizeof(int*) * nrows)

        n = 0
        for i from 0 <= i < nrows:
            self._row_indices[i] = n
            n = n + ncols

    def  __dealloc__(self):
        if self._row_indices != <int*> 0:
            PyMem_Free(self._row_indices)

    def __getitem__(self, ij):
        """
        INPUT:
            A[i, j] -- the i,j of A, and
            A[i]    -- the i-th row of A.
        """
        #self._require_mutable()   # todo
        cdef int i, j
        i, j = ij
        if i < 0 or i >= self._nrows:
            raise IndexError
        return self._entries[self._row_indices[i] + j]

    def __setitem__(self, ij, value):
        """
        INPUT:
            A[i, j] = value -- set the (i,j) entry of A
            A[i] = value    -- set the ith row of A
        """
        #self._require_mutable()   # todo
        cdef int i, j
        i, j = ij
        if i < 0 or i >= self._nrows:
            raise IndexError
        self._entries[self._row_indices[i] + j] = value

    def _multiply(self, Matrix_dense right):
        cdef int i, j, k, m, n, r, nr, nc, snc
        cdef object v

        if self._ncols != right._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."

        cdef Matrix_dense A
        nr = self._nrows
        nc = right._ncols
        snc = self._ncols

        R = self.base_ring()
        P = self.matrix_space(nr, nc)
        A = Matrix_dense(P, nr, nc)
        v = A._entries
        zero = R(0)
        r = 0
        for i from 0 <= i < nr:
            m = self._row_indices[i]
            for j from 0 <= j < nc:
                z = zero
                for k from 0 <= k < snc:
                    #z = z + self._entries[m + k] * right._entries[right._row_indices[k]+j]

                    #z = z + PyList_GET_ITEM(self._entries, m + k)._mul_(
                    #    PyList_GET_ITEM(right._entries, right._row_indices[k]+j))

                    z = z + multiply_items(self._entries, m+k, right._entries, right._row_indices[k]+j)

                    #z = z + PyNumber_Multiply( PyList_GET_ITEM(self._entries, m + k),
                    #                  PyList_GET_ITEM(right._entries, right._row_indices[k]+j))

                v[r] = z
                r = r + 1
        return A

    def list(self):
        return self._entries

    def antitranspose(self):
        raise NotImplementedError

    def transpose(self):
        """
        Returns the transpose of self, without changing self.

        EXAMPLES:
        We create a matrix, compute its transpose, and note that the
        original matrix is not changed.
            sage: M = MatrixSpace(QQ,  2)
            sage: A = M([1,2,3,4])
            sage: B = A.transpose()
            sage: print B
            [1 3]
            [2 4]
            sage: print A
            [1 2]
            [3 4]
        """
        raise NotImplementedError



cdef object multiply_items(object v, int i, object w, int j):
    return PyNumber_Multiply( PyList_GET_ITEM(v, i), PyList_GET_ITEM(w, j) )
