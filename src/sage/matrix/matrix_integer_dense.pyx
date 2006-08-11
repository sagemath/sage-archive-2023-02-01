"""nodoctest
Dense matrices over the integers.
"""

######################################################################
#       Copyright (C) 2006 William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
######################################################################

from sage.misc.misc import verbose, get_verbose

include "../ext/gmp.pxi"
include "../ext/interrupt.pxi"

cimport sage.ext.integer
import  sage.ext.integer

cimport matrix_pid
import matrix_pid

cdef class Matrix_integer_dense(matrix_pid.Matrix_pid):
    """
    Matrix over the integers.
    """
    def __new__(self, parent, object entries=None, construct=False):
        cdef int i, nrows, ncols
        self.initialized = 0
        if isinstance(entries, str) and entries == LEAVE_UNINITIALIZED:
            self.matrix = <mpz_t **>0
            return
        nrows = parent.nrows()
        ncols = parent.ncols()
        self.matrix = <mpz_t **> PyMem_Malloc(sizeof(mpz_t*)*nrows)
        if self.matrix == <mpz_t**> 0:
            raise MemoryError, "Error allocating matrix."
        for i from 0 <= i < nrows:
            self.matrix[i] = <mpz_t *> PyMem_Malloc(sizeof(mpz_t)*ncols)
            if self.matrix[i] == <mpz_t *> 0:
                raise MemoryError, "Error allocating matrix."

    cdef int allocate(self) except -1:
        cdef int nrows, ncols, i
        nrows = parent.nrows()
        ncols = parent.ncols()

    def __init__(self, parent, object entries=None, construct=False):
        cdef int n, i, j, k, r, base, nrows, ncols
        cdef mpz_t *v

        matrix_field.Matrix_field.__init__(self, parent)

        nrows = parent.nrows()
        ncols = parent.ncols()
        self._nrows = nrows
        self._ncols = ncols
        self.__pivots = None
        base = 10
        if isinstance(entries, str):
            if entries == LEAVE_UNINITIALIZED:
                return
            elif construct:
                base = 32
                entries = entries.split(' ')

        if isinstance(entries, sage.ext.rational.Rational):
            if entries != 0 and nrows != ncols:
                raise TypeError, "scalar matrix must be square"
            s = str(entries)
            mpz_init(self.tmp)
            r = mpz_set_str(self.tmp, s, 0)
            if r == -1:
                raise TypeError, "Invalid rational number %s"%entries[k]
            _sig_on
            for i from 0 <= i < nrows:
                v = self.matrix[i]
                for j from 0 <= j < ncols:
                    mpz_init(v[j])
                    if i == j:
                        mpz_set(v[j], self.tmp)
                        k = k + 1
                    else:
                        mpz_set_si(v[j], 0)
            _sig_off
            self.initialized = 1
            return

        if nrows*ncols != 0:
            if entries != None and len(entries) != nrows*ncols:
                raise IndexError, "The vector of entries has length %s but should have length %s"%(len(entries), nrows*ncols)

        _sig_on
        k = 0
        for i from 0 <= i < nrows:
            v = self.matrix[i]
            for j from 0 <= j < ncols:
                mpz_init(v[j])
                if entries != None:
                    # TODO: If entries[k] is a rational,
                    # this should be WAY faster.  (Also see above)
                    s = str(entries[k])
                    r = mpz_set_str(v[j], s, base)
                    if r == -1:
                        _sig_off
                        raise TypeError, "Invalid rational number %s"%entries[k]
                    k = k + 1
                else:
                    mpz_set_si(v[j],0)
        _sig_off
        self.initialized = 1

    def nrows(self):
        return self._nrows

    def ncols(self):
        return self._ncols

    def __reduce__(self):
        import sage.matrix.reduce

        cdef int i, j, len_so_far, m, n
        cdef char *a
        cdef char *s, *t, *tmp

        if self._nrows == 0 or self._ncols == 0:
            entries = ''
        else:
            n = self._nrows*self._ncols*10
            s = <char*> PyMem_Malloc(n * sizeof(char))
            t = s
            len_so_far = 0

            _sig_on
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    m = mpz_sizeinbase (self.matrix[i][j], 32)
                    if len_so_far + m + 1 >= n:
                        # copy to new string with double the size
                        n = 2*n + m + 1
                        tmp = <char*> PyMem_Malloc(n * sizeof(char))
                        strcpy(tmp, s)
                        PyMem_Free(s)
                        s = tmp
                        t = s + len_so_far
                    #endif
                    mpz_get_str(t, 32, self.matrix[i][j])
                    m = strlen(t)
                    len_so_far = len_so_far + m + 1
                    t = t + m
                    t[0] = <char>32
                    t[1] = <char>0
                    t = t + 1
            _sig_off
            entries = str(s)[:-1]
            free(s)

        return sage.matrix.reduce.make_Matrix_integer_dense, \
               (self.parent(), entries)


    def __cmp__(self, other):
        if not isinstance(other, Matrix_integer_dense):
            return -1
        raise NotImplementedError

    def __setitem__(self, ij, x):
        i, j = ij
        if self.matrix == <mpz_t **>0:
            raise RuntimeError, "Matrix has not yet been initialized!"
        if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
            raise IndexError, "Invalid index."
        s = str(x)
        mpz_set_str(self.matrix[i][j], s, 0)

    def __getitem__(self, ij):
        raise NotImplementedError

    def  __dealloc__(self):
        cdef int i, j
        if self.matrix == <mpz_t **> 0:
            return
        for i from 0 <= i < self._nrows:
            if self.matrix[i] != <mpz_t *> 0:
                for j from 0 <= j < self._ncols:
                    if self.initialized:
                        mpz_clear(self.matrix[i][j])
                PyMem_Free(self.matrix[i])
        PyMem_Free(self.matrix)

    def _mul_(Matrix_integer_dense self, Matrix_integer_dense other):
        if self._ncols != other._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."

        cdef int i, j, k, nr, nc, snc
        cdef mpz_t *v
        cdef mpz_t s, z
        nr = self._nrows
        nc = other._ncols
        snc = self._ncols

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense(self.parent(), LEAVE_UNINITIALIZED)

        cdef mpz_t **m
        m = <mpz_t **> PyMem_Malloc(sizeof(mpz_t*)*nr)
        if m == <mpz_t**> 0:
            raise MemoryError, "Error allocating matrix"

        mpz_init(s); mpz_init(z)

        _sig_on
        for i from 0 <= i < nr:
            m[i] = <mpz_t *> PyMem_Malloc(sizeof(mpz_t)*nc)
            if m[i] == <mpz_t*> 0:
                mpz_clear(s); mpz_clear(z)
                _sig_off
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                mpz_set_si(s,0)   # set s = 0
                v = self.matrix[i]
                for k from 0 <= k < snc:
                    mpz_mul(z, v[k], other.matrix[k][j])
                    mpz_add(s, s, z)
                mpz_init(m[i][j])
                mpz_set(m[i][j], s)
        _sig_off
        M.set_matrix(m)
        mpz_clear(s); mpz_clear(z)
        return M

    def _add_(Matrix_integer_dense self, Matrix_integer_dense other):
        if self._ncols != other._ncols:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self._nrows != other._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."

        cdef int i, j, nr, nc
        nr = self._nrows
        nc = other._ncols

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense(nr, nc, LEAVE_UNINITIALIZED)

        cdef mpz_t **m
        m = <mpz_t **> PyMem_Malloc(sizeof(mpz_t*)*nr)
        if m == <mpz_t**> 0:
            raise MemoryError, "Error allocating matrix"

        _sig_on
        for i from 0 <= i < nr:
            m[i] = <mpz_t *> PyMem_Malloc(sizeof(mpz_t)*nc)
            if m[i] == <mpz_t*> 0:
                _sig_off
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                mpz_init(m[i][j])
                mpz_add(m[i][j], self.matrix[i][j], other.matrix[i][j])
        _sig_off
        M.set_matrix(m)
        return M

    cdef set_matrix(Matrix_integer_dense self, mpz_t **m):
        if self.matrix != <mpz_t **> 0:
            raise RuntimeError, "Only set matrix of uninitialized matrix."
        self.matrix = m
        self.initialized = 1

    def transpose(self):
        """
        Returns the transpose of self.
        """
        cdef int i, j
        cdef Matrix_integer_dense M

        M = Matrix_integer_dense(self._ncols, self._nrows, entries=LEAVE_UNINITIALIZED)
        cdef mpz_t **m
        m = <mpz_t **> PyMem_Malloc(sizeof(mpz_t*)*self._ncols)
        if m == <mpz_t**> 0:
            raise MemoryError, "Error allocating matrix"

        _sig_on
        for i from 0 <= i < self._ncols:
            m[i] = <mpz_t *> PyMem_Malloc(sizeof(mpz_t)*self._nrows)
            if m[i] == <mpz_t*> 0:
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < self._nrows:
                mpz_init(m[i][j])
                mpz_set(m[i][j], self.matrix[j][i])
        _sig_off
        M.set_matrix(m)
        return M

    def matrix_from_rows(self, rows):
        """
        Return the submatrix formed from the given rows.

        INPUT:
            rows -- list of int's

        OUTPUT:
            matrix created from the rows with given indexes
        """
        cdef int i, j, k, nc, nr
        cdef Matrix_integer_dense M

        if not isinstance(rows, list):
            raise TypeError, "rows (=%s) must be a list"%rows
        nr = len(rows)
        if nr == 0:
            return Matrix_integer_dense(0, self._ncols)
        nc = self._ncols
        v = []
        for i in rows:
            v.append(int(i))
        rows = v
        if min(rows) < 0 or max(rows) >= self._nrows:
            raise IndexError, "invalid row indexes; rows don't exist"

        M = Matrix_integer_dense(self.parent(), entries=LEAVE_UNINITIALIZED)
        cdef mpz_t **m
        m = <mpz_t **> PyMem_Malloc(sizeof(mpz_t*)*nr)
        if m == <mpz_t**> 0:
            raise MemoryError, "Error allocating matrix"

        for i from 0 <= i < nr:
            m[i] = <mpz_t *> PyMem_Malloc(sizeof(mpz_t)*nc)
            if m[i] == <mpz_t*> 0:
                raise MemoryError, "Error allocating matrix"
            k = rows[i]
            for j from 0 <= j < nc:
                mpz_init(m[i][j])
                mpz_set(m[i][j], self.matrix[k][j])

        M.set_matrix(m)
        return M



    def iterates(self, v, int n):
        """
        Let A be this matrix.   Return a matrix with rows
        $$
           v, Av, A^2v, ..., A^(n-1)v.
        $$
        """
        cdef int i, j, k, nr, nc
        cdef mpz_t s, z
        nr = n
        nc = self._ncols

        if self._nrows != self._ncols:
            raise ArithmeticError, "matrix must be square"
        if not isinstance(v, list):
            raise TypeError, "v must be a list"
        if len(v) != self._nrows:
            raise ArithmeticError, "incompatible matrix vector multiple"

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense(self.parent(), LEAVE_UNINITIALIZED)

        cdef mpz_t **m
        m = <mpz_t **> PyMem_Malloc(sizeof(mpz_t*)*nr)
        if m == <mpz_t**> 0:
            raise MemoryError, "Error allocating matrix"
        m[0] = <mpz_t *> PyMem_Malloc(sizeof(mpz_t)*nc)
        if m[0] == <mpz_t*> 0:
            mpz_clear(s); mpz_clear(z)
            raise MemoryError, "Error allocating matrix"
        mpz_init(self.tmp)
        for j from 0 <= j < nc:
            string = str(v[j])
            r = mpz_set_str(self.tmp, string, 0)
            if r == -1:
                raise TypeError, "Invalid integer"
            mpz_init(m[0][j])
            mpz_set(m[0][j], self.tmp)

        mpz_init(s)
        mpz_init(z)
        for i from 1 <= i < nr:
            m[i] = <mpz_t *> PyMem_Malloc(sizeof(mpz_t)*nc)
            if m[i] == <mpz_t*> 0:
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                mpz_set_si(s,0)  # set s = 0
                for k from 0 <= k < self._nrows:
                    mpz_mul(z, m[i-1][k], self.matrix[k][j])
                    mpz_add(s, s, z)
                mpz_init(m[i][j])
                mpz_set(m[i][j], s)
        M.set_matrix(m)
        mpz_clear(s); mpz_clear(z)
        return M


    def scalar_multiple(self, d):
        """
        Return the product self*d, as a new matrix.
        """
        cdef int i, j, nr, nc
        nr = self._nrows
        nc = self._ncols

        cdef mpz_t x
        mpz_init(x)
        s = str(d)
        r = mpz_set_str(x, s, 0)
        if r == -1:
            raise TypeError, "Invalid rational number %s"%entries[k]
        cdef Matrix_integer_dense M
        M = Matrix_integer_dense(self.parent(), LEAVE_UNINITIALIZED)

        cdef mpz_t **m
        m = <mpz_t **> PyMem_Malloc(sizeof(mpz_t*)*nr)
        if m == <mpz_t**> 0:
            raise MemoryError, "Error allocating matrix"

        _sig_on
        for i from 0 <= i < nr:
            m[i] = <mpz_t *> PyMem_Malloc(sizeof(mpz_t)*nc)
            if m[i] == <mpz_t*> 0:
                mpz_clear(x)
                _sig_off
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                mpz_init(m[i][j])
                mpz_mul(m[i][j], self.matrix[i][j], x)
        _sig_off
        M.set_matrix(m)
        mpz_clear(x)
        return M

    def copy(self):
        cdef int i, j, nr, nc
        nr = self._nrows; nc = self._ncols

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense(self.parent(), LEAVE_UNINITIALIZED)
        cdef mpz_t **m
        m = <mpz_t **> PyMem_Malloc(sizeof(mpz_t*)*nr)
        if m == <mpz_t**> 0:
            raise MemoryError, "Error allocating matrix"

        for i from 0 <= i < nr:
            m[i] = <mpz_t *> PyMem_Malloc(sizeof(mpz_t)*nc)
            if m[i] == <mpz_t*> 0:
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                mpz_init(m[i][j])
                mpz_set(m[i][j], self.matrix[i][j])

        M.set_matrix(m)
        return M

    def number_nonzero(self):
        cdef int i, j, n
        cdef mpz_t *v
        n = 0
        _sig_on
        for i from 0 <= i < self._nrows:
            v = self.matrix[i]
            for j from 0 <= j < self._ncols:
                if mpz_sgn(v[j]):   # if nonzero
                    n = n + 1
        _sig_off
        return n

    def list(self, int base=0):
        cdef int i, j
        cdef mpz_t *r
        cdef object v
        cdef sage.ext.integer.Integer x

        v = []
        _sig_on
        for i from 0 <= i < self._nrows:
            r = self.matrix[i]
            for j from 0 <= j < self._ncols:
                x = sage.ext.integer.Integer()
                x.set_from_mpz(r[j])
                v.append(x)
        _sig_off
        return v

    def rank(self):
        """
        Return the rank found during the last echelon operation on self.
        Of course if self is changed, and the echelon form of self is not
        recomputed, then the rank could be incorrect.
        """
        if self.__pivots == None:
            raise ArithmeticError, "Echelon form has not yet been computed."
        return len(self.__pivots)

    def pivots(self):
        """
        Return the pivots found during the last echelon operation on self.
        Of course if self is changed, and the echelon form of self is not
        recomputed, then the pivots could be incorrect.
        """
        if self.__pivots == None:
            raise ArithmeticError, "Echelon form has not yet been computed."
        return self.__pivots

    cdef int mpz_height(self, mpz_t height) except -1:
        cdef mpz_t x, h
        mpz_init(x)
        mpz_init_set_si(h, 0)
        cdef int i, j
        _sig_on
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                mpq_get_num(x,self.matrix[i][j])
                mpz_abs(x, x)
                if mpz_cmp(h,x) < 0:
                    mpz_set(h,x)
                mpq_get_den(x,self.matrix[i][j])
                mpz_abs(x, x)
                if mpz_cmp(h,x) < 0:
                    mpz_set(h,x)
        _sig_off
        mpz_set(height, h)
        mpz_clear(h)
        mpz_clear(x)
        return 0

    def height(self):
        cdef mpz_t h
        mpz_init(h)
        self.mpz_height(h)
        a = mpz_to_long(h)
        mpz_clear(h)
        return a


    def multiply_multi_modular(self, Matrix_integer_dense right):
        """
        Multiply this matrix by right using a multimodular algorithm
        and return the result.
        """
        raise NotImplementedError

cdef object mpz_to_long(mpz_t x):
    return long(mpz_to_str(x))

