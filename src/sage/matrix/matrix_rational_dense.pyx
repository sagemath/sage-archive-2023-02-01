"""
Dense matrices over the rational field.

This is a compiled implementation of dense matrix algebra over small
prime finite fields and the rational numbers, which is used mainly
internally by other classes.

TODO:
    -- do one big allocation instead of lots of small ones.
"""

#*****************************************************************************
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
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
#*****************************************************************************

from sage.misc.misc import verbose, get_verbose

include "../ext/gmp.pxi"
include "../ext/interrupt.pxi"

cimport matrix_field
import matrix_field

cimport matrix_dense
import matrix_dense

cimport sage.ext.rational
import  sage.ext.rational

cdef class Matrix_rational_dense(matrix_field.Matrix_field):
    """
    Matrix over the rational numbers.
    """
    def __init__(self, parent, object entries=None, construct=False, zero=True):
        cdef int n, i, j, k, r, base, nrows, ncols
        cdef mpq_t *v

        matrix_field.Matrix_field.__init__(self, parent)

        nrows = parent.nrows()
        ncols = parent.ncols()
        self._nrows = nrows
        self._ncols = ncols

        self._entries = <mpq_t *> PyMem_Malloc(sizeof(mpq_t)*(nrows*ncols))
        if self._entries == <mpq_t *> 0:
            raise MemoryError, "Error allocating matrix."

        for i from 0 <= i < (nrows*ncols):
            mpq_init(self._entries[i])

        self._matrix =  <mpq_t **> PyMem_Malloc(sizeof(mpq_t*)*ncols)
        if self._matrix == <mpq_t**> 0:
            raise MemoryError, "Error allocating matrix."

        k = 0
        for i from 0 <= i < nrows:
            self._matrix[i] = self._entries + k
            k = k + ncols

        if zero:
            for i from 0 <= i < (nrows*ncols):
                mpq_set_si(self._entries[i], 0, 1)

        self.__pivots = None
        base = 10
        if isinstance(entries, str):
            if construct:
                base = 32
                entries = entries.split(' ')
                raise NotImplementedError, "need to deal with base below"

        if isinstance(entries, sage.ext.rational.Rational):
            if entries != 0 and nrows != ncols:
                raise TypeError, "scalar matrix must be square"
            s = str(entries)
            mpq_init(self.tmp)
            r = mpq_set_str(self.tmp, s, 0)
            if r == -1:
                raise TypeError, "Invalid rational number"
            mpq_canonicalize(self.tmp)
            _sig_on
            for i from 0 <= i < nrows:
                v = self._matrix[i]
                for j from 0 <= j < ncols:
                    if i == j:
                        mpq_set(v[j], self.tmp)
                        k = k + 1
                    else:
                        mpq_set_si(v[j], 0, 1)
            _sig_off
            return

        if nrows*ncols != 0:
            if not (entries is None) and len(entries) != nrows*ncols:
                raise IndexError, "The vector of entries has length %s but should have length %s"%(len(entries), nrows*ncols)

        if entries is None:
            if zero:
                for i from 0 <= i < nrows * ncols:
                    mpq_set_si(self._entries[i], 0, 1)
            return

        cdef sage.ext.rational.Rational z

        if coerce:
            for i from 0 <= i < nrows*ncols:
                z = sage.ext.rational.Rational(entries[i])
                mpq_set(self._entries[i], z.value)
        else:
            for i from 0 <= i < nrows*ncols:
                z = entries[i]
                mpq_set(self._entries[i], z.value)


    def nrows(self):
        return self._nrows

    def ncols(self):
        return self._ncols

    def __reduce__(self):
        import sage._matrix.reduce

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
                    m = mpz_sizeinbase (mpq_numref(self._matrix[i][j]), 32) + \
                        mpz_sizeinbase (mpq_denref(self._matrix[i][j]), 32) + 3
                    if len_so_far + m + 1 >= n:
                        # copy to new string with double the size
                        n = 2*n + m + 1
                        tmp = <char*> PyMem_Malloc(n * sizeof(char))
                        strcpy(tmp, s)
                        PyMem_Free(s)
                        s = tmp
                        t = s + len_so_far
                    #endif
                    mpq_get_str(t, 32, self._matrix[i][j])
                    m = strlen(t)
                    len_so_far = len_so_far + m + 1
                    t = t + m
                    t[0] = <char>32
                    t[1] = <char>0
                    t = t + 1
            _sig_off
            entries = str(s)[:-1]
            free(s)

        return sage._matrix.reduce.make_Matrix_rational_dense, \
               (self.parent(), entries)


    def __cmp__(self, other):
        cdef int i, c
        cdef Matrix_rational_dense x
        if isinstance(other, Matrix_rational_dense):
            x = other
            if self._nrows != x._nrows or self._ncols != x._ncols:
                return 1
            for i from 0 <= i < self._nrows * self._ncols:
                if not mpq_equal(self._entries[i], x._entries[i]):
                    return 1
            return 0
        else:
            return -1

    def __setitem__(self, ij, x):
        i, j = ij
        if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
            raise IndexError, "Invalid index."
        cdef sage.ext.rational.Rational y
        try:
            y = x
        except TypeError:
            y = sage.ext.rational.Rational(x)
        mpq_set(self._matrix[i][j], y.value)

    def __getitem__(self, ij):
        i, j = ij
        if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
            raise IndexError, "Invalid index."
        cdef sage.ext.rational.Rational x
        x = sage.ext.rational.Rational()
        x.set_from_mpq(self._matrix[i][j])
        return x

    def  __dealloc__(self):
        PyMem_Free(self._entries)
        PyMem_Free(self._matrix)


    def _mul_(Matrix_rational_dense self, Matrix_rational_dense other):
        if self._ncols != other._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."

        cdef int i, j, k, nr, nc, snc
        cdef mpq_t *v
        cdef mpq_t s, z
        nr = self._nrows
        nc = other._ncols
        snc = self._ncols

        cdef Matrix_rational_dense M
        M = Matrix_rational_dense(self.new_matrix(self._nrows, other._ncols).parent(), zero=False) # clean up when MatrixSpace creates this class

        cdef mpq_t **m
        m = M._matrix

        mpq_init(s); mpq_init(z)

        _sig_on
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                mpq_set_si(s,0,1)   # set s = 0
                v = self._matrix[i]
                for k from 0 <= k < snc:
                    mpq_mul(z, v[k], other._matrix[k][j])
                    mpq_add(s, s, z)
                mpq_set(m[i][j], s)
        _sig_off
        mpq_clear(s); mpq_clear(z)
        return M

    def _add_(Matrix_rational_dense self, Matrix_rational_dense other):
        if self._ncols != other._ncols:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self._nrows != other._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."

        cdef int i, j, nr, nc
        nr = self._nrows
        nc = other._ncols

        cdef Matrix_rational_dense M
        M = Matrix_rational_dense(self.parent(), zero=False)

        cdef mpq_t **m
        m = M._matrix

        _sig_on
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                mpq_add(m[i][j], self._matrix[i][j], other._matrix[i][j])
        _sig_off
        return M

    def transpose(self):
        """
        Returns the transpose of self.
        """
        cdef int i, j
        cdef Matrix_rational_dense M

        M = Matrix_rational_dense(self.parent(), zero=False)
        cdef mpq_t **m
        m = M._matrix

        _sig_on
        for i from 0 <= i < self._ncols:
            for j from 0 <= j < self._nrows:
                mpq_set(m[i][j], self._matrix[j][i])
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
        cdef Matrix_rational_dense M

        if not isinstance(rows, list):
            raise TypeError, "rows (=%s) must be a list"%rows
        nr = len(rows)
        if nr == 0:
            return Matrix_rational_dense(0, self._ncols)
        nc = self._ncols
        v = []
        for i in rows:
            v.append(int(i))
        rows = v
        if min(rows) < 0 or max(rows) >= self._nrows:
            raise IndexError, "invalid row indexes; rows don't exist"

        M = Matrix_rational_dense(self.parent(), zero=False)
        cdef mpq_t **m
        m = M._matrix

        for i from 0 <= i < nr:
            k = rows[i]
            for j from 0 <= j < nc:
                mpq_init(m[i][j])
                mpq_set(m[i][j], self._matrix[k][j])

        return M


    def matrix_from_cols(self, cols):
        """
        Return the submatrix formed from the given columns.

        INPUT:
            cols -- list of int's

        OUTPUT:
            matrix created from the columns with given indexes
        """
        cdef int i, j, k, nc, nr
        cdef Matrix_rational_dense M

        if not isinstance(cols, list):
            raise TypeError, "cols (=%s) must be a list"%cols
        nc = len(cols)
        if nc == 0:
            return Matrix_rational_dense(self._nrows, 0)
        nr = self._nrows
#        v = []
#        for i in rows:
#            v.append(int(i))
#        rows = v
        if min(cols) < 0 or max(cols) >= self._ncols:
            raise IndexError, "invalid cols indexes; cols don't exist"

        M = self.new_matrix(nrows=nr, ncols=nc, zero=False) #Matrix_rational_dense(self.parent(), zero=False)
        cdef mpq_t **m
        m = M._matrix

        for j from 0 <= j < nc:
            k = int(cols[j])
            for i from 0 <= i < nr:
                mpq_init(m[i][j])
                mpq_set(m[i][j], self._matrix[i][k])

        return M




    def iterates(self, v, int n):
        """
        Let A be this matrix.   Return a matrix with rows
        $$
           v, Av, A^2v, ..., A^(n-1)v.
        $$
        """
        cdef int i, j, k, nr, nc
        cdef mpq_t s, z
        nr = n
        nc = self._ncols

        if self._nrows != self._ncols:
            raise ArithmeticError, "matrix must be square"
        if not isinstance(v, list):
            raise TypeError, "v must be a list"
        if len(v) != self._nrows:
            raise ArithmeticError, "incompatible matrix vector multiple"

        cdef Matrix_rational_dense M
        M = Matrix_rational_dense(self.parent(), zero=False)

        cdef mpq_t **m
        m = M._matrix

        mpq_init(self.tmp)
        for j from 0 <= j < nc:
            string = str(v[j])
            r = mpq_set_str(self.tmp, string, 0)
            if r == -1:
                raise TypeError, "Invalid rational number"
            mpq_set(m[0][j], self.tmp)

        mpq_init(s)
        mpq_init(z)
        for i from 1 <= i < nr:
            for j from 0 <= j < nc:
                mpq_set_si(s,0,1)  # set s = 0
                for k from 0 <= k < self._nrows:
                    mpq_mul(z, m[i-1][k], self._matrix[k][j])
                    mpq_add(s, s, z)
                mpq_set(m[i][j], s)
        mpq_clear(s); mpq_clear(z)
        return M


    def scalar_multiple(self, d):
        """
        Return the product self*d, as a new matrix.
        """
        cdef int i, j, nr, nc
        nr = self._nrows
        nc = self._ncols

        cdef mpq_t x
        mpq_init(x)
        s = str(d)
        r = mpq_set_str(x, s, 0)
        if r == -1:
            raise TypeError, "Invalid rational number"
        cdef Matrix_rational_dense M
        M = Matrix_rational_dense(self.parent(), zero=False)

        cdef mpq_t **m
        m = M._matrix

        _sig_on
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                mpq_mul(m[i][j], self._matrix[i][j], x)
        _sig_off
        mpq_clear(x)
        return M

    def copy(self):
        cdef int i, j, nr, nc
        nr = self._nrows; nc = self._ncols

        cdef Matrix_rational_dense M
        M = Matrix_rational_dense(self.parent(), zero=False)
        cdef mpq_t **m
        m = M._matrix

        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                mpq_set(m[i][j], self._matrix[i][j])

        return M

    # TODO: this function should be removed when self.parent().matrix() returns this pyx class
    def new_matrix(self, nrows=None, ncols=None, entries=0,
                   coerce_entries=True, copy=True, sparse=None,
                   clear = True, zero=True):
      return Matrix_rational_dense(self.matrix_space(nrows, ncols))


    def number_nonzero(self):
        cdef int i, j, n
        cdef mpq_t *v
        n = 0
        _sig_on
        for i from 0 <= i < self._nrows * self._ncols:
            if mpq_sgn(self._entries[i]):         # if nonzero
                n = n + 1
        _sig_off
        return n

    def list(self, int base=0):
        cdef int i, j
        cdef mpq_t *r
        cdef object v
        cdef sage.ext.rational.Rational x

        v = []
        _sig_on
        for i from 0 <= i < self._nrows:
            r = self._matrix[i]
            for j from 0 <= j < self._ncols:
                x = sage.ext.rational.Rational()
                x.set_from_mpq(r[j])
                v.append(x)
        _sig_off
        return v

    def echelon_gauss_in_place(self):
        """
        Changes self into echelon form.
        """
        cdef int start_row, c, r, nr, nc, i
        cdef mpq_t **m
        cdef mpq_t a_inverse, minus_b

        mpq_init(a_inverse)
        mpq_init(minus_b)
        start_row = 0
        m = self._matrix
        nr = self._nrows
        nc = self._ncols
        self.__pivots = []

        for c from 0 <= c < nc:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            for r from start_row <= r < nr:
                if mpq_sgn(m[r][c]):
                    self.__pivots.append(c)
                    mpq_inv(a_inverse,m[r][c])
                    self.scale_row(r, a_inverse, c)
                    self.swap_rows(r, start_row)
                    for i from 0 <= i < nr:
                        if i != start_row:
                            if mpq_sgn(m[i][c]):
                                mpq_neg(minus_b, m[i][c])
                                self.add_multiple_of_row(start_row, minus_b, i, c)
                    start_row = start_row + 1
                    break
        mpq_clear(a_inverse)
        mpq_clear(minus_b)

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

    def _set_pivots(self, v):
        self.__pivots = v

    def matrix_window(self, int row=0, int col=0, int nrows=-1, int ncols=-1):
        if nrows == -1:
            nrows = self._nrows - row
            ncols = self._ncols - col
        return MatrixWindow(self, row, col, nrows, ncols)

    def hessenberg_form(self):
        if not self.is_immutable():
            raise ValueError, "matrix must be mutable, since hessenberg form changes it"
        if self._nrows != self._ncols:
            raise ArithmeticError, "Matrix must be square to compute Hessenberg form."

        cdef int n
        n = self._nrows

        cdef mpq_t **h
        h = self._matrix

        cdef int i, j, m, p, r
        cdef mpq_t t, t_inv, u, neg_u
        mpq_init(t)
        mpq_init(t_inv)
        mpq_init(u)
        mpq_init(neg_u)

        for m from 1 <= m < n-1:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            # Search for a nonzero entry in column m-1
            i = -1
            for r from m+1 <= r < n:
                if mpq_sgn(h[r][m-1]):
                     i = r
                     break

            if i != -1:
                 # Found a nonzero entry in column m-1 that is strictly
                 # below row m.  Now set i to be the first nonzero position >=
                 # m in column m-1.
                 if mpq_sgn(h[m][m-1]): i = m
                 mpq_set(t,h[i][m-1])
                 mpq_inv(t_inv, t)
                 if i > m:
                     self.swap_rows(i,m)
                     self.swap_columns(i,m)

                 # Now the nonzero entry in position (m,m-1) is t.
                 # Use t to clear the entries in column m-1 below m.
                 for j from m+1 <= j < n:
                     if mpq_sgn(h[j][m-1]):
                         mpq_mul(u,h[j][m-1], t_inv)
                         mpq_neg(neg_u, u)
                         self.add_multiple_of_row(m, neg_u, j, 0)  # h[j] -= u*h[m]
                         # To maintain charpoly, do the corresponding
                         # column operation, which doesn't mess up the
                         # matrix, since it only changes column m, and
                         # we're only worried about column m-1 right
                         # now.  Add u*column_j to column_m.
                         self.add_multiple_of_column(j, u, m, 0)
                 # end for
            # end if
        # end for
        mpq_clear(t)
        mpq_clear(t_inv)
        mpq_clear(u)
        mpq_clear(neg_u)

    cdef scale_row(self, int row, mpq_t multiple, int start_col):
        cdef int r
        cdef mpq_t* v

        r = row*self._ncols
        v = self._matrix[row]
        for i from start_col <= i < self._ncols:
            mpq_mul(v[i], v[i], multiple)

    cdef add_multiple_of_row(self, int row_from, mpq_t multiple,
                            int row_to, int start_col):
        cdef int i
        cdef mpq_t *v_from, *v_to
        cdef mpq_t prod, x

        mpq_init(prod); mpq_init(x)
        v_from = self._matrix[row_from]
        v_to = self._matrix[row_to]
        for i from start_col <= i < self._ncols:
            mpq_mul(prod, multiple, v_from[i])
            mpq_add(x, prod, v_to[i])
            mpq_set(v_to[i], x)   # v_to[i] <-- multipe*v_from[i] + v_to[i]

        mpq_clear(prod); mpq_clear(x)

    def set_row_to_multiple_of_row(self, int row_to, int row_from, sage.ext.rational.Rational multiple):
        """
        Set row row_to equal to multiple times row row_from.
        """
        cdef int i
        cdef mpq_t *v_from, *v_to

        if row_from < 0 or row_from >= self._nrows:
            raise IndexError, "row_from is %s but must be >= 0 and < %s"%(row_from, self._nrows)
        if row_to < 0 or row_to >= self._nrows:
            raise IndexError, "row_to is %s but must be >= 0 and < %s"%(row_to, self._nrows)

        v_from = self._matrix[row_from]
        v_to = self._matrix[row_to]
        for i from 0 <= i < self._ncols:
            mpq_mul(v_to[i], multiple.value, v_from[i])


    cdef add_multiple_of_column(self, int col_from, mpq_t multiple,
                               int col_to, int start_row):
        cdef int i, p, nr
        cdef mpq_t **m
        cdef mpq_t prod, x

        mpq_init(prod); mpq_init(x)
        m = self._matrix
        nr = self._nrows
        for i from start_row <= i < self._nrows:
            mpq_mul(prod, multiple, m[i][col_from])
            mpq_add(x, m[i][col_to], prod)
            mpq_set(m[i][col_to], x)
        mpq_clear(prod); mpq_clear(x)

    cdef swap_rows(Matrix_rational_dense self, int row1, int row2):
        cdef mpq_t* temp
        temp = self._matrix[row1]
        self._matrix[row1] = self._matrix[row2]
        self._matrix[row2] = temp

    cdef swap_columns(self, int col1, int col2):
        cdef int i, nr
        cdef mpq_t **m
        cdef mpq_t t

        mpq_init(t)
        m = self._matrix
        nr = self._nrows
        for i from 0 <= i < self._nrows:
            mpq_set(t, m[i][col1])
            mpq_set(m[i][col1], m[i][col2])
            mpq_set(m[i][col2], t)
        mpq_clear(t)

    cdef int mpz_denom(self, mpz_t d) except -1:
        cdef mpz_t y
        mpz_set_si(d,1)
        mpz_init(y)
        cdef int i, j
        _sig_on
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                mpq_get_den(y,self._matrix[i][j])
                mpz_lcm(d, d, y)
        _sig_off
        mpz_clear(y)
        return 0

    def denom(self):
        cdef mpz_t d
        mpz_init(d)
        self.mpz_denom(d)
        dl = mpz_to_long(d)
        mpz_clear(d)
        return dl

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

    def prod_of_row_sums(self, cols):
        r"""
        Calculate the product of all row sums of a submatrix of $A$ for a
        list of selected columns \code{cols}.

        This is used for the computation of matrix permanents.
        """
        cdef int row, c, n, t

        n = len(cols)
        cdef int* v
        v = <int*> PyMem_Malloc(n * sizeof(int))
        for c from 0 <= c < n:
            t = cols[c]
            if t < 0 or t >= self._ncols:
                PyMem_Free(v)
                raise IndexError, "invalid column index (= %s)"%t
            v[c] = t

        cdef mpq_t pr, z
        mpq_init(pr)
        mpq_init(z)

        mpq_set_si(pr, 1, 1)
        for row from 0 <= row < self._nrows:
            mpq_set_si(z, 0, 1)
            for c from 0 <= c < n:
                mpq_add(z, z, self._matrix[row][v[c]])
            mpq_mul(pr, pr, z)

        cdef sage.ext.rational.Rational x
        x = sage.ext.rational.Rational()
        x.set_from_mpq(pr)
        mpq_clear(pr)
        mpq_clear(z)
        PyMem_Free(v)
        return x

    def _clear_denom(self):
        """
        INPUT:
            self -- a matrix
        OUTPU:
            self, D, if D=denominator is 1
            D*self, D if D > 1.

        Thus returns a copy of self only if D > 1.
        """
        cdef mpz_t d
        mpz_init(d)
        self.mpz_denom(d)
        if mpz_cmp_si(d,1) == 0:
            mpz_clear(d)
            return self, sage.ext.rational.Rational(1)
        cdef Matrix_rational_dense A
        A = self.copy()
        cdef mpq_t denom
        mpq_init(denom)
        mpq_set_z(denom, d)
        A._rescale(denom)
        mpz_clear(d)
        mpq_clear(denom)

        cdef sage.ext.rational.Rational x
        x = sage.ext.rational.Rational()
        x.set_from_mpq(denom)
        return A, x

    cdef int _rescale(self, mpq_t a) except -1:
        cdef int i, j
        _sig_on
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                mpq_mul(self._matrix[i][j], self._matrix[i][j], a)
        _sig_off

    def echelon(self, alg="modular", height_guess=None):
        """
        echelon(self, alg="modular", height_guess=None):

        Returns echelon form of self, without modifying self.
        """
        if alg=="modular":
            return self.echelon_modular(height_guess=height_guess)
        elif alg=="gauss":
            A = self.copy()
            A.echelon_gauss()
            return A
        elif alg=="block":
            A = self.copy()
            A.echelon_strassen()
            return A
        else:
            raise ValueError, "%s is not one of the allowed algorithms (modular, gauss)"%alg

    def echelon_gauss(self):
        """
        Returns echelon form of self using gaussian elimination, modifying self.
        """
        pivots = []
        cdef int row, col, i, j
        cdef mpq_t recip, front, tmp
        mpq_init(recip)
        mpq_init(front)
        mpq_init(tmp)
        row = 0
        for col from 0 <= col < self._ncols:
            if mpq_sgn(self._matrix[row][col]) == 0:
                for i from row < i < self._nrows:
                    if mpq_sgn(self._matrix[i][col]) != 0:
                        self.swap_rows(i, row)
                        break
                if i == self._nrows:
                    continue
            pivots.append(col)
            mpq_inv(recip, self._matrix[row][col])
            mpq_set_si(self._matrix[row][col], 1, 1)
            for j from col < j < self._ncols:
                mpq_mul(self._matrix[row][j], self._matrix[row][j], recip)
            for i from 0 <= i < self._nrows:
                if i == row: continue
                if mpq_sgn(self._matrix[i][col]) == 0: continue
                mpq_set(front, self._matrix[i][col])
                mpq_set_si(self._matrix[i][col], 0, 1)
                for j from col < j < self._ncols:
                    mpq_mul(tmp, front, self._matrix[row][j])
                    mpq_sub(self._matrix[i][j], self._matrix[i][j], tmp)
            row = row+1
            if row == self._nrows:
                break
        self._set_pivots(pivots)



    def echelon_modular(self, height_guess=None):
        """
        echelon_modular(self, height_guess=None):

        Returns echelon form of self, without modifying self.  Uses a
        multi-modular method.

        ALGORITHM:
        The following is a modular algorithm for computing the echelon
        form.  Define the height of a matrix to be the max of the
        absolute values of the entries.

        Input: Matrix A with n columns (this).

        0. Rescale input matrix A to have integer entries.  This does
           not change echelon form and makes reduction modulo many
           primes significantly easier if there were denominators.
           Henceforth we assume A has integer entries.

        1. Let c be a guess for the height of the echelon form.  E.g.,
           c=1000, since matrix is sparse and application is modular
           symbols.

        2. Let M = n * c * H(A) + 1,
           where n is the number of columns of A.

        3. List primes p_1, p_2, ..., such that the product of
           the p_i is at least M.

        4. Try to compute the rational reconstruction CRT echelon form
           of A mod the product of the p_i.  Throw away those A mod p_i
           whose pivot sequence is not >= all other pivot sequences of
           A mod p_j.
           If rational reconstruction fails, compute 1 more echelon
           forms mod the next prime, and attempt again.  Let E be this
           matrix.

        5. Compute the denominator d of E.
           Try to prove that result is correct by checking that

                 H(d*E) < (prod of reduction primes)/(ncols*H(A)),

           where H denotes the height.   If this fails, do step 4 with
           a few more primes.

           (TODO: Possible idea for optimization: When doing the rational_recon lift,
            keep track of the lcm d of denominators found so far, and given
                             a (mod m)
            first check to see if a*d lifts to an integer with abs <= m/2.
            If so, no nded to do rational recon.  This should be the case
            for most a after a while, and should save substantial time!!!!)
        """
        B, _ = self._clear_denom()
        hA = B.height()
        if height_guess is None:
            height_guess = (2*hA)**(self._ncols/2+1)
        verbose("height_guess=%s"%height_guess)
        M = self._ncols * height_guess * hA  +  1
        p = START_PRIME
        X = []
        best_pivots = []
        prod = 1
        while True:
            while prod < M:
                verbose("p=%s"%p)
                A = B.matrix_modint(p)
                A.echelon()
                if self._nrows == self._ncols and len(A.pivots()) == self._ncols:
                    # special case -- the echelon form must be the identity matrix.
                    return Matrix_rational_identity(self._nrows)

                c = cmp_pivots(best_pivots, A.pivots())
                if c <= 0:
                    best_pivots = A.pivots()
                    X.append(A)
                    prod = prod * p
                else:
                    if get_verbose() > 1:
                        verbose("Excluding this prime (bad pivots).", level=2)
                    pass   # do not save A since it is bad.
                #p = previous_probab_prime_int(p)
                p = next_probab_prime_int(p)

            Y = []
            prod = 1
            # We recompute product, since may drop bad matrices
            for i from 0 <= i < len(X):
                # Here best_pivots is the best collection
                # of pivots found during any echelon form computation.
                # Here cmp_pivots returns a number <= 0 if
                # X[i].pivots() is at least as good.
                if cmp_pivots(best_pivots, X[i].pivots()) <= 0:
                    # append a good matrix to the list Y.
                    Y.append(X[i])
                    # multiply the product of the good primes by this good prime
                    prod = prod * X[i].prime()
            try:
                t = verbose("start rr")
                E = Matrix_rational_using_crt_and_rr(Y)
                verbose("done",t)
            except ValueError:
                for i from 0 <= i < 10:
                    M = M * START_PRIME
                verbose("(Failed to compute rational reconstruction -- redoing with several more primes", level=2)
                continue
            d = E.denom()
            Es = E.scalar_multiple(d)
            hdE = (Es).height()
            if hdE * hA * self._ncols < prod:
                self.__pivots = best_pivots
                E._set_pivots(list(best_pivots))
                return E
            for i from 0 <= i < 3:
                M = M * START_PRIME

    def multiply_multi_modular(self, Matrix_rational_dense right):
        """
        Multiply this matrix by right using a multimodular algorithm
        and return the result.
        """
        if self._ncols != right._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of right."

        A, A_denom = self._clear_denom()
        B, B_denom = right._clear_denom()
        bound = 2 * A.height() * B.height() * A.ncols()
        p = 0
        X = []
        prod = 1
        while prod < bound:
            verbose('prod = %s, bound = %s'%(prod, bound))
            if p == 0:
                p = START_PRIME
            else:
                #p = previous_probab_prime_int(p)
                p = next_probab_prime_int(p)
            t = verbose("p=%s"%p)
            A_modp = A.matrix_modint(p)
            B_modp = B.matrix_modint(p)
            t = verbose("done reducing", t)
            C_modp = A_modp.strassen(B_modp)
            t = verbose("done multiplying", t)
            X.append(C_modp)
            prod = prod * p
        t = verbose("now doing CRT")
        C = Matrix_rational_CRT(X)
        verbose("finished CRT", t)
        return C

cdef object mpz_to_long(mpz_t x):
    return long(mpz_to_str(x))



cdef class MatrixWindow:

    def __init__(MatrixWindow self, Matrix_rational_dense matrix, int row, int col, int nrows, int ncols):
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
        """
        return self._matrix


    def to_matrix(MatrixWindow self):
        """
        Returns an actual matrix object representing this view.
        """
        a = self._matrix.new_matrix(self._nrows, self._ncols, zero=False)
        a.matrix_window().set_to(self)
        return a


    def matrix_window(MatrixWindow self, int row=0, int col=0, int n_rows=-1, int n_cols=-1):
        """
        Returns a matrix window relative to this window of the underlying matrix.
        """
        if row == 0 and col == 0 and n_rows == self._nrows and n_cols == self._ncols:
            return self
        return self._matrix.matrix_window(self._row + row, self._col + col, n_rows, n_cols)

    def nrows(MatrixWindow self):
        return self._nrows

    def ncols(MatrixWindow self):
        return self._ncols

    def set_to(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef mpq_t* s_row
        cdef mpq_t* A_row
        for i from 0 <= i < self._nrows:
            s_row = self._matrix._matrix[self._row + i] + self._col
            A_row = A._matrix._matrix[A._row + i] + A._col
            for j from 0 <= j < self._ncols:
                mpq_set(s_row[j], A_row[j])

    def set_to_zero(MatrixWindow self):
        cdef int i, j
        cdef mpq_t* s_row
        for i from 0 <= i < self._nrows:
            s_row = self._matrix._matrix[self._row + i] + self._col
            for j from 0 <= j < self._ncols:
                mpq_set_ui(s_row[j], 0, 1)

    def add(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef mpq_t* s_row
        cdef mpq_t* A_row
        for i from 0 <= i < self._nrows:
            s_row = self._matrix._matrix[self._row + i] + self._col
            A_row = A._matrix._matrix[A._row + i] + A._col
            for j from 0 <= j < self._ncols:
                mpq_add(s_row[j], s_row[j], A_row[j])

    def subtract(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef mpq_t* s_row
        cdef mpq_t* A_row
        for i from 0 <= i < self._nrows:
            s_row = self._matrix._matrix[self._row + i] + self._col
            A_row = A._matrix._matrix[A._row + i] + A._col
            for j from 0 <= j < self._ncols:
                mpq_sub(s_row[j], s_row[j], A_row[j])

    def set_to_sum(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j
        cdef mpq_t* s_row
        cdef mpq_t* A_row
        cdef mpq_t* B_row
        for i from 0 <= i < self._nrows:
            s_row = self._matrix._matrix[self._row + i] + self._col
            A_row = A._matrix._matrix[A._row + i] + A._col
            B_row = B._matrix._matrix[B._row + i] + B._col
            for j from 0 <= j < self._ncols:
                mpq_add(s_row[j], A_row[j], B_row[j])

    def set_to_diff(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j
        cdef mpq_t* s_row
        cdef mpq_t* A_row
        cdef mpq_t* B_row
        for i from 0 <= i < self._nrows:
            s_row = self._matrix._matrix[self._row + i] + self._col
            A_row = A._matrix._matrix[A._row + i] + A._col
            B_row = B._matrix._matrix[B._row + i] + B._col
            for j from 0 <= j < self._ncols:
                mpq_sub(s_row[j], A_row[j], B_row[j])

    def set_to_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j, k
        cdef mpq_t* s_row
        cdef mpq_t* A_row
        cdef mpq_t sum, prod
        mpq_init(sum)
        mpq_init(prod)
        for i from 0 <= i < A._nrows:
            A_row = A._matrix._matrix[A._row + i] + A._col
            s_row = self._matrix._matrix[self._row + i] + self._col
            for j from 0 <= j < B._ncols:
                mpq_mul(sum, A_row[0], B._matrix._matrix[B._row]+B._col+j)
                for k from 1 <= k < A._ncols:
                    mpq_mul(prod, A_row[k], B._matrix._matrix[B._row+k]+B._col+j)
                    mpq_add(sum, sum, prod)
                mpq_set(s_row[j], sum)

    def add_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j, k
        cdef mpq_t* s_row
        cdef mpq_t* A_row
        cdef mpq_t sum, prod
        mpq_init(sum)
        mpq_init(prod)
        for i from 0 <= i < A._nrows:
            A_row = A._matrix._matrix[A._row + i] + A._col
            s_row = self._matrix._matrix[self._row + i] + self._col
            for j from 0 <= j < B._ncols:
                mpq_mul(sum, A_row[0], B._matrix._matrix[B._row]+B._col+j)
                for k from 1 <= k < A._ncols:
                    mpq_mul(prod, A_row[k], B._matrix._matrix[B._row+k]+B._col+j)
                    mpq_add(sum, sum, prod)
                mpq_add(s_row[j], s_row[j], sum)

    def subtract_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j, k
        cdef mpq_t* s_row
        cdef mpq_t* A_row
        cdef mpq_t sum, prod
        mpq_init(sum)
        mpq_init(prod)
        for i from 0 <= i < A._nrows:
            A_row = A._matrix._matrix[A._row + i] + A._col
            s_row = self._matrix._matrix[self._row + i] + self._col
            for j from 0 <= j < B._ncols:
                mpq_mul(sum, A_row[0], B._matrix._matrix[B._row]+B._col+j)
                for k from 1 <= k < A._ncols:
                    mpq_mul(prod, A_row[k], B._matrix._matrix[B._row+k]+B._col+j)
                    mpq_add(sum, sum, prod)
                mpq_sub(s_row[j], s_row[j], sum)


    def swap_rows(MatrixWindow self, int a, int b):
        self._matrix.swap_rows(self._row + a, self._row + b)


    def echelon_in_place(MatrixWindow self):
        """
        calculate the echelon form of this matrix, returning the list of pivot columns
        """
        echelon = self.to_matrix().echelon(alg="gauss") # TODO: read only, only need to copy pointers
        self.set_to(echelon.matrix_window())
        return echelon.pivots()

    def element_is_zero(MatrixWindow self, int i, int j):
        return mpq_sgn(self._matrix._matrix[i+self._row][j+self._col]) == 0


    def new_empty_window(MatrixWindow self, int nrows, int ncols, zero=True):
        return self._matrix.new_matrix(nrows, ncols, zero=zero).matrix_window()
