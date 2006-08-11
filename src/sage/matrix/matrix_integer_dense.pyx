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

    def __init__(self, parent, object entries=None, construct=False, zero=True,
                 coerce=True):
        cdef int n, i, j, k, r, base, nrows, ncols
        cdef mpz_t *v

        matrix_pid.Matrix_pid.__init__(self, parent)

        nrows = parent.nrows()
        ncols = parent.ncols()
        self._nrows = nrows
        self._ncols = ncols

        self._entries = <mpz_t *>PyMem_Malloc(sizeof(mpz_t*) * (nrows * ncols))

        if self._entries == <mpz_t*> 0:
            raise MemoryError, "Error allocating matrix."

        for i from 0 <= i < (nrows * ncols):
            mpz_init(self._entries[i])

        self._matrix = <mpz_t **> PyMem_Malloc(sizeof(mpz_t*)*nrows)
        if self._matrix == <mpz_t**> 0:
            raise MemoryError, "Error allocating matrix."

        k = 0
        for i from 0 <= i < nrows:
            self._matrix[i] = self._entries + k
            k = k + ncols

        self.__pivots = None
        base = 10
        if isinstance(entries, str):
            if construct:
                base = 32
                entries = entries.split(' ')
                coerce = True
                raise NotImplementedError, "need to deal with base below"

        cdef sage.ext.integer.Integer z

        if isinstance(entries, sage.ext.integer.Integer):
            if entries != 0 and nrows != ncols:
                raise TypeError, "scalar matrix must be square"
            mpz_init(self.tmp)
            z = entries
            mpz_set(self.tmp, z.value)
            _sig_on
            for i from 0 <= i < nrows:
                v = self._matrix[i]
                for j from 0 <= j < ncols:
                    if i == j:
                        mpz_set(v[j], self.tmp)
                    else:
                        mpz_set_si(v[j], 0)
            _sig_off
            return

        if nrows*ncols != 0:
            if entries != None and len(entries) != nrows*ncols:
                raise IndexError, "The vector of entries has length %s but should have length %s"%(len(entries), nrows*ncols)

        if entries is None:
            if zero:
                for i from 0 <= i < nrows * ncols:
                    mpz_set_si(self._entries[i], 0)
            return

        k = 0

        if coerce:
            for i from 0 <= i < nrows*ncols:
                z = sage.ext.integer.Integer(entries[i])
                mpz_set(self._entries[i], z.value)
        else:
            for i from 0 <= i < nrows*ncols:
                z = entries[i]
                mpz_set(self._entries[i], z.value)

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
                    m = mpz_sizeinbase (self._matrix[i][j], 32)
                    if len_so_far + m + 1 >= n:
                        # copy to new string with double the size
                        n = 2*n + m + 1
                        tmp = <char*> PyMem_Malloc(n * sizeof(char))
                        strcpy(tmp, s)
                        PyMem_Free(s)
                        s = tmp
                        t = s + len_so_far
                    #endif
                    mpz_get_str(t, 32, self._matrix[i][j])
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
        cdef int i
        cdef Matrix_integer_dense c_other
        if not isinstance(other, Matrix_integer_dense):
            return -1
        c_other = other
        if self._nrows != c_other._rows or self._ncols != c_other._ncols:
            return -1
        for i from 0 <= i < self._nrows * self._ncols:
            if mpz_cmp(self._entries[i], c_other._entries[i]):
                return -1
        return 0

    def __setitem__(self, ij, x):
        i, j = ij
        cdef sage.ext.integer.Integer z
        if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
            raise IndexError, "Invalid index."
        try:
            z = x
        except (TypeError):
            z = sage.ext.integer.Integer(x)
        mpz_set(self._matrix[i][j], z.value)

    def __getitem__(self, ij):
        cdef sage.ext.integer.Integer z
        cdef int i, j

        i, j = ij

        if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
            raise IndexError, "Invalid index."

        z = sage.ext.integer.Integer()
        mpz_set(z.value, self._matrix[i][j])
        return z

    def  __dealloc__(self):
        PyMem_Free(self._entries)
        PyMem_Free(self._matrix)

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
        M = Matrix_integer_dense(self.parent(), zero=False)

        cdef mpz_t **m
        m = M._matrix

        mpz_init(s); mpz_init(z)

        _sig_on
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                mpz_set_si(s,0)   # set s = 0
                v = self._matrix[i]
                for k from 0 <= k < snc:
                    mpz_mul(z, v[k], other._matrix[k][j])
                    mpz_add(s, s, z)
                mpz_init(m[i][j])
                mpz_set(m[i][j], s)
        _sig_off
        mpz_clear(s); mpz_clear(z)
        return M

    def _add_(Matrix_integer_dense self, Matrix_integer_dense other):
        if self._ncols != other._ncols:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self._nrows != other._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."

        cdef int i, j

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense(nr, nc, zero=False)

        cdef mpz_t **m
        m = M._matrix

        _sig_on
        for i from 0 <= i < self._ncols * self._nrows:
            mpz_init(m[i])
            mpz_add(m[i], self._entries[i], other._entries[i])
        _sig_off
        return M

    def transpose(self):
        """
        Returns the transpose of self.
        """
        cdef int i, j
        cdef Matrix_integer_dense M

        M = Matrix_integer_dense(self._ncols, self._nrows, zero=False)
        cdef mpz_t **m
        m = M._matrix

        _sig_on
        for i from 0 <= i < self._ncols:
            for j from 0 <= j < self._nrows:
                mpz_init(m[i][j])
                mpz_set(m[i][j], self._matrix[j][i])
        _sig_off
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

        M = Matrix_integer_dense(self.parent(), zero=False)
        cdef mpz_t **m
        m = M._matrix

        for i from 0 <= i < nr:
            k = rows[i]
            for j from 0 <= j < nc:
                mpz_init(m[i][j])
                mpz_set(m[i][j], self._matrix[k][j])

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
        M = Matrix_integer_dense(self.parent(), zero=False)

        cdef mpz_t **m
        m = M._matrix
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
                    mpz_mul(z, m[i-1][k], self._matrix[k][j])
                    mpz_add(s, s, z)
                mpz_init(m[i][j])
                mpz_set(m[i][j], s)
        mpz_clear(s); mpz_clear(z)
        return M


    def scalar_multiple(self, d):
        """
        Return the product self*d, as a new matrix.
        """
        cdef int i, j, nr, nc
        cdef sage.ext.integer.Integer z
        nr = self._nrows
        nc = self._ncols

        cdef mpz_t x
        try:
            z = d
        except (TypeError):
            z = sage.ext.integer.Integer(d)
        mpz_init(x)
        mpz_set(x, z.value)

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense(self.parent(), zero=False)

        cdef mpz_t *e
        e = M._entries

        _sig_on
        for i from 0 <= i < nr * nc:
            mpz_init(e[i])
            mpz_mul(e[i], self._entries[i], x)
        _sig_off
        mpz_clear(x)
        return M

    def copy(self):
        cdef int i, j, nr, nc
        nr = self._nrows; nc = self._ncols

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense(self.parent(), zero=False)
v v v v v v v
*************
        cdef mpz_t *e
        e = M._entries
^ ^ ^ ^ ^ ^ ^

        for i from 0 <= i < nr * nc:
            mpz_set(M._entries[i], self._entries[i])

        for i from 0 <= i < nc:
            M._matrix[i] = self._matrix[i]

        return M

    def number_nonzero(self):
        cdef int i, j, n
        cdef mpz_t *v
        n = 0
        _sig_on
        for i from 0 <= i < self._nrows * self._ncols:
            if mpz_sgn(self._entries[i]):         # if nonzero
                n = n + 1
        _sig_off
        return n

    def list(self, int base=0):
        cdef int i, j
        cdef mpz_t* r
        cdef object v
        cdef sage.ext.integer.Integer x

        v = []
        _sig_on
        for i from 0 <= i < self._nrows:
            r = self._matrix[i]
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
                mpq_get_num(x,self._matrix[i][j])
                mpz_abs(x, x)
                if mpz_cmp(h,x) < 0:
                    mpz_set(h,x)
                mpq_get_den(x,self._matrix[i][j])
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

