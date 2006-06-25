r"""
Dense matrices over $\FF_p$ and $\QQ$.
"""

#*****************************************************************************
#       Copyright (C) 2004- William Stein <wstein@ucsd.edu>
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

include "gmp.pxi"
include "interrupt.pxi"

from sage.misc.misc import verbose, get_verbose


cimport rational
import rational

cimport arith
import arith
cdef arith.arith_int ai
ai = arith.arith_int()



##################################################################
# Miscellaneous useful functions
##################################################################
cdef int mpq_to_modint(mpq_t x, int n):
    """
    Convert an mpq to an integer modulo n.  There is substantial
    overhead because three mpq_t variables must be initialized and
    cleared each time this function is called.
    """
    cdef mpz_t num, den
    mpz_init(num); mpz_init(den)
    mpq_get_num(num, x)
    mpq_get_den(den, x)

    cdef int numer, denom, t
    cdef mpz_t r

    _sig_on
    mpz_init(r)
    numer = mpz_mod_ui(r, num, n)
    denom = mpz_mod_ui(r, den, n)
    mpz_clear(num)
    mpz_clear(den)
    mpz_clear(r)
    t = (numer * ai.c_inverse_mod_int(denom, n)) % n
    _sig_off

cdef int gmp_random_int(int min, int max):
    pass


##################################################################
# Matrix over integers modulo p, for p <= 46340
##################################################################

LEAVE_UNINITIALIZED = "LEAVE UNINITIALIZED"

ctypedef unsigned int uint

cdef class Matrix_modint:
    cdef uint **matrix
    cdef uint nrows, ncols, p
    cdef uint gather
    cdef object __pivots

    def __new__(self, int p, int nrows, int ncols, object entries=None):
        cdef int i
        if entries == "LEAVE UNINITIALIZED":
            self.matrix = <uint **>0
            return
        self.matrix = <uint **> PyMem_Malloc(sizeof(uint*)*nrows)
        if self.matrix == <uint**> 0:
            raise MemoryError, "Error allocating memory"
        for i from 0 <= i < nrows:
            self.matrix[i] = <uint *> PyMem_Malloc(sizeof(uint)*ncols)
            if self.matrix[i] == <uint*> 0:
               raise MemoryError, "Error allocating matrix"

    def __init__(self, uint p, uint nrows, uint ncols, object entries=None):
        cdef uint n, i, j, k
        cdef uint *v
        if p >= 46340:
            raise OverflowError, "p (=%s) must be < 46340"%p

        self.p = p
        self.nrows = nrows
        self.ncols = ncols
        self.gather = 2**32/(p*p)
        self.__pivots = None
        if entries == LEAVE_UNINITIALIZED:
            return
        if entries is None:
            for i from 0 <= i < nrows:
                v = self.matrix[i]
                for j from 0 <= j < ncols:
                    v[j] = 0
            return
        if len(entries) != nrows*ncols:
            raise IndexError, "The vector of entries has the wrong length."
        k = 0
        for i from 0 <= i < nrows:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            v = self.matrix[i]
            for j from 0 <= j < ncols:
                n = entries[k]
                v[j] = n % p
                if v[j] < 0:
                    v[j] = v[j] + p
                k = k + 1

    cdef set_matrix(Matrix_modint self, uint **m):
        if self.matrix != <uint **>0:
            raise RuntimeError, "Only set matrix of uninitialized matrix."
        self.matrix = m

    def get_entry(Matrix_modint self, int i, int j):
        if i < 0 or j < 0 or i >= self.nrows or j >= self.ncols:
            raise IndexError, "Array index (%s,%s) out of bounds."%(i,j)
        return self.matrix[i][j]

    def __call__(self, i, j):
        return self.matrix[i][j]

    def _cmp(self, Matrix_modint other):
        cdef uint i, j
        if self.p != other.p or self.nrows != other.nrows \
            or self.ncols != other.ncols:
            return -1
        _sig_on
        for i from 0 <= i < self.nrows:
            for j from 0 <= j < self.ncols:
                if self.matrix[i][j] != other.matrix[i][j]:
                    _sig_off
                    return -1
        _sig_off
        return 0

    def __cmp__(self, other):
        if not isinstance(other, Matrix_modint):
            return -1
        return self._cmp(other)

    cdef uint **get_matrix(Matrix_modint self):
        return self.matrix

    def  __dealloc__(self):
        if self.matrix == <uint **> 0:
            return
        cdef int i
        for i from 0 <= i < self.nrows:
            if self.matrix[i] != <uint *> 0:
                PyMem_Free(self.matrix[i])
        PyMem_Free(self.matrix)

    def __mul__(Matrix_modint self, Matrix_modint other):
        if self.gather >= 2:
            return self._multiply_with_delayed_mod(other)
        if self.ncols != other.nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self.p != other.p:
            raise ArithmeticError, "The base matrices must have the same modulus."

        cdef Matrix_modint M
        M = Matrix_modint(self.p, self.nrows, other.ncols, LEAVE_UNINITIALIZED)
        cdef uint **m
        m = <uint **> PyMem_Malloc(sizeof(uint*)*self.nrows)
        if m == <uint**> 0:
            raise MemoryError, "Error allocating matrix"

        cdef uint i, j, k, nr, nc, s, p
        cdef uint *v
        nr = self.nrows
        nc = other.ncols
        p = self.p

        _sig_on
        for i from 0 <= i < nr:
            m[i] = <uint *> PyMem_Malloc(sizeof(uint)*nc)
            if m[i] == <uint*> 0:
                _sig_off
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                s = 0
                v = self.matrix[i]
                for k from 0 <= k < self.ncols:
                    s = (s + (v[k] * other.matrix[k][j]))%p
                m[i][j] = s
        _sig_off
        M.set_matrix(m)
        return M

    def _multiply_with_delayed_mod(Matrix_modint self, Matrix_modint other):
        if self.ncols != other.nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self.p != other.p:
            raise ArithmeticError, "The base matrices must have the same modulus."

        cdef Matrix_modint M
        M = Matrix_modint(self.p, self.nrows, other.ncols, LEAVE_UNINITIALIZED)
        cdef uint **m
        m = <uint **> PyMem_Malloc(sizeof(int*)*self.nrows)
        if m == <uint**> 0:
            raise MemoryError, "Error allocating matrix"

        cdef uint i, j, k, nr, nc, snc, s, p, gather, w, groups, a, b
        cdef uint *v
        nr = self.nrows
        nc = other.ncols
        snc = self.ncols
        gather = self.gather
        p = self.p
        groups = (self.ncols / gather) + 1
        _sig_on
        for i from 0 <= i < nr:
            m[i] = <uint *> PyMem_Malloc(sizeof(uint)*nc)
            if m[i] == <uint*> 0:
                _sig_off
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                s = 0
                v = self.matrix[i]
                for w from 0 <= w < groups:
                    a = w*gather
                    b = (w+1)*gather
                    if b > snc or b < a: b = snc
                    for k from a <= k < b:
                        s = s + v[k]*other.matrix[k][j]
                    s = s % p
                ##for k from 0 <= k < self.ncols:
                ##    s = s  +  v[k] * other.matrix[k][j]
                m[i][j] = s
        _sig_off
        M.set_matrix(m)
        return M

    def __repr__(self):
        cdef uint i, j
        s = "[\n"
        for i from 0 <= i < self.nrows:
            for j from 0 <= j < self.ncols:
                s = s + str(self.matrix[i][j]) + ", "
            s = s + "\n"
        s = s[:-3] + "\n]"
        return s

    def numrows(self):
        return self.nrows

    def numcols(self):
        return self.ncols

    def prime(self):
        return self.p

    def number_nonzero(self):
        cdef uint i, j, n
        cdef uint *v
        n = 0
        _sig_on
        for i from 0 <= i < self.nrows:
            v = self.matrix[i]
            for j from 0 <= j < self.ncols:
                if v[j] != 0:
                    n = n + 1
        _sig_off
        return n

    def list(self):
        cdef uint i, j
        cdef uint *r
        v = []
        _sig_on
        for i from 0 <= i < self.nrows:
            r = self.matrix[i]
            for j from 0 <= j < self.ncols:
                v.append(r[j])
        _sig_off
        return v

    def echelon(self):
        cdef uint p, start_row, c, r, nr, nc, a, a_inverse, b, i
        cdef uint **m

        start_row = 0
        p = self.p
        m = self.matrix
        nr = self.nrows
        nc = self.ncols
        self.__pivots = []
        for c from 0 <= c < nc:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            for r from start_row <= r < nr:
                a = m[r][c]
                if a:
                    self.__pivots.append(c)
                    a_inverse = ai.c_inverse_mod_int(a, p)
                    self.scale_row(r, a_inverse, c)
                    self.swap_rows(r, start_row)
                    for i from 0 <= i < nr:
                        if i != start_row:
                            b = m[i][c]
                            if b != 0:
                                self.add_multiple_of_row(start_row, p-b, i, c)
                    start_row = start_row + 1
                    break

    def rank(self):
        """
        Return the rank found during the last echelon operation on self.
        Of course if self is changed, then the rank could be incorrect.
        """
        if self.__pivots == None:
            raise ArithmeticError, "Echelon form has not yet been computed."
        return len(self.__pivots)

    def pivots(self):
        """
        Return the pivots found during the last echelon operation on self.
        Of course if self is changed, then the pivots could be incorrect.
        """
        if self.__pivots == None:
            raise ArithmeticError, "Echelon form has not yet been computed."
        return self.__pivots

    def hessenberg_form(self):
        """
        Transforms self in place to its Hessenberg form.
        """
        if self.nrows != self.ncols:
            raise ArithmeticError, "Matrix must be square to compute Hessenberg form."

        cdef uint n
        n = self.nrows

        cdef uint **h
        h = self.matrix

        cdef uint p, r, t, t_inv, u
        cdef int i, j, m
        p = self.p

        _sig_on
        for m from 1 <= m < n-1:
            # Search for a nonzero entry in column m-1
            i = -1
            for r from m+1 <= r < n:
                if h[r][m-1]:
                     i = r
                     break

            if i != -1:
                 # Found a nonzero entry in column m-1 that is strictly
                 # below row m.  Now set i to be the first nonzero position >=
                 # m in column m-1.
                 if h[m][m-1]:
                     i = m
                 t = h[i][m-1]
                 t_inv = ai.c_inverse_mod_int(t,p)
                 if i > m:
                     self.swap_rows(i,m)
                     self.swap_columns(i,m)

                 # Now the nonzero entry in position (m,m-1) is t.
                 # Use t to clear the entries in column m-1 below m.
                 for j from m+1 <= j < n:
                     if h[j][m-1]:
                         u = (h[j][m-1] * t_inv) % p
                         self.add_multiple_of_row(m, p - u, j, 0)  # h[j] -= u*h[m]
                         # To maintain charpoly, do the corresponding
                         # column operation, which doesn't mess up the
                         # matrix, since it only changes column m, and
                         # we're only worried about column m-1 right
                         # now.  Add u*column_j to column_m.
                         self.add_multiple_of_column(j, u, m, 0)
                 # end for
            # end if
        # end for
        _sig_off

    def charpoly(self):
        """
        Transforms self in place to its Hessenberg form then computes
        and returns the coefficients of the characteristic polynomial of
        this matrix.

        The characteristic polynomial is represented as a vector of
        ints, where the constant term of the characteristic polynomial
        is the 0th coefficient of the vector.
        """
        if self.nrows != self.ncols:
            raise ArithmeticError, "charpoly not defined for non-square matrix."

        cdef uint i, m, n, p, t
        n = self.nrows
        p = self.p

        # Replace self by its Hessenberg form, and set H to this form
        # for notation below.
        #time = verbose('start hessenberg')
        self.hessenberg_form()
        #verbose('done with hessenberg', time)
        cdef Matrix_modint H
        H = self

        # We represent the intermediate polynomials that come up in
        # the calculations as rows of an (n+1)x(n+1) matrix, since
        # we've implemented basic arithmetic with such a matrix.
        # Please see the generic implementation of charpoly in
        # matrix.py to see more clearly how the following algorithm
        # actually works.  (The implementation is clearer if one uses
        # polynomials to represent polynomials instead of using the
        # rows of a matrix.)  Also see Cohen's first GTM, Algorithm
        # 2.2.9.

        cdef Matrix_modint c
        c = Matrix_modint(self.p, n+1, n+1)  # the 0 matrix
        c.matrix[0][0] = 1
        for m from 1 <= m <= n:
            # Set the m-th row of c to (x - H[m-1,m-1])*c[m-1] = x*c[m-1] - H[m-1,m-1]*c[m-1]
            # We do this by hand by setting the m-th row to c[m-1]
            # shifted to the right by one.  We then add
            # -H[m-1,m-1]*c[m-1] to the resulting m-th row.
            for i from 1 <= i <= n:
                c.matrix[m][i] = c.matrix[m-1][i-1]
            # the p-.. below is to keep scalar normalized between 0 and p.
            c.add_multiple_of_row(m-1, p - H.matrix[m-1][m-1], m, 0)
            t = 1
            for i from 1 <= i < m:
                t = (t*H.matrix[m-i][m-i-1]) % p
                # Set the m-th row of c to c[m] - t*H[m-i-1,m-1]*c[m-i-1]
                c.add_multiple_of_row(m-i-1, p - (t*H.matrix[m-i-1][m-1])%p, m, 0)

        # The answer is now the n-th row of c.
        v = []
        for i from 0 <= i <= n:
            v.append(int(c.matrix[n][i]))
        return v

    def decompose(self):
        raise NotImplementedError

    def rational_reconstruction(self):
        raise NotImplementedError

    cdef uint entry(self, uint i, uint j):
        return self.matrix[i][j]

    cdef scale_row(self, uint row, uint multiple, uint start_col):
        cdef uint r, p
        cdef uint* v
        r = row*self.ncols
        p = self.p
        v = self.matrix[row]
        for i from start_col <= i < self.ncols:
            v[i] = (v[i]*multiple) % p

    cdef add_multiple_of_row(self, uint row_from, uint multiple,
                            uint row_to, uint start_col):
        cdef uint i, p, nc
        cdef uint *v_from, *v_to
        p = self.p
        v_from = self.matrix[row_from]
        v_to = self.matrix[row_to]
        nc = self.ncols
        for i from start_col <= i < nc:
            v_to[i] = (multiple * v_from[i] +  v_to[i]) % p

    cdef add_multiple_of_column(self, uint col_from, uint multiple,
                               uint col_to, uint start_row):
        cdef uint i, p, nr
        cdef uint **m
        m = self.matrix
        p = self.p
        nr = self.nrows
        for i from start_row <= i < self.nrows:
            m[i][col_to] = (m[i][col_to] + multiple * m[i][col_from]) %p


    cdef swap_rows(self, uint row1, uint row2):
        cdef uint* temp
        temp = self.matrix[row1]
        self.matrix[row1] = self.matrix[row2]
        self.matrix[row2] = temp

    cdef swap_columns(self, uint col1, uint col2):
        cdef uint i, t, nr
        cdef uint **m
        m = self.matrix
        nr = self.nrows
        for i from 0 <= i < self.nrows:
            t = m[i][col1]
            m[i][col1] = m[i][col2]
            m[i][col2] = t

cdef uint **Matrix_modint_matrix(Matrix_modint A):
    return A.matrix



##################################################################
# Matrix over the rational numbers
# The rational numbers are implemented here using GMP MPQ.
##################################################################

START_PRIME = 20011  # used for modular algorithms

def cmp_pivots(x,y):
    """
    Compare two sequences of pivot columns.
    If x is short than y, return -1, i.e., x < y, "not as good".
    If x is longer than y, x > y, "better"
    If the length is the same then x is better, i.e., x > y
        if the entries of x are correspondingly >= those of y with
        one being greater.
    I
    """
    if len(x) < len(y):
        return -1
    if len(x) > len(y):
        return 1
    if x < y:
        return 1
    elif x == y:
        return 0
    else:
        return -1



#cdef class Matrix_rational

cdef Matrix_rational_cmp(Matrix_rational self, Matrix_rational other):
    if self.matrix == <mpq_t **> 0 or other.matrix == <mpq_t **> 0:
        raise RuntimeError, "Matrix has not yet been initialized!"
    if self.nrows != other.nrows or self.ncols != other.ncols:
        return False
    cdef int i, j, k
    _sig_on
    for i from 0 <= i < self.nrows:
        for j from 0 <= j < self.ncols:
            k = mpq_cmp(self.matrix[i][j], other.matrix[i][j])
            if k < 0:
                _sig_off
                return -1
            elif k > 0:
                _sig_off
                return 1
    _sig_off
    return 0

cdef class Matrix_rational:
    """
    Matrix over the rational numbers.
    """

    def __new__(self, int nrows, int ncols, object entries=None, construct=False):
        cdef int i
        self.initialized = 0
        if isinstance(entries, str) and entries == LEAVE_UNINITIALIZED:
            self.matrix = <mpq_t **>0
            return
        self.matrix = <mpq_t **> PyMem_Malloc(sizeof(mpq_t*)*nrows)
        if self.matrix == <mpq_t**> 0:
            raise MemoryError, "Error allocating matrix."
        for i from 0 <= i < nrows:
            self.matrix[i] = <mpq_t *> PyMem_Malloc(sizeof(mpq_t)*ncols)
            if self.matrix[i] == <mpq_t *> 0:
                raise MemoryError, "Error allocating matrix."

    def __init__(self, int nrows, int ncols, object entries=None, construct=False):
        cdef int n, i, j, k, r, base
        cdef mpq_t *v
        self.nrows = nrows
        self.ncols = ncols
        self.__pivots = None
        base = 10
        if isinstance(entries, str):
            if entries == LEAVE_UNINITIALIZED:
                return
            elif construct:
                base = 32
                entries = entries.split(' ')

        if isinstance(entries, rational.Rational):
            if entries != 0 and nrows != ncols:
                raise TypeError, "scalar matrix must be square"
            s = str(entries)
            mpq_init(self.tmp)
            r = mpq_set_str(self.tmp, s, 0)
            if r == -1:
                raise TypeError, "Invalid rational number %s"%entries[k]
            mpq_canonicalize(self.tmp)
            _sig_on
            for i from 0 <= i < nrows:
                v = self.matrix[i]
                for j from 0 <= j < ncols:
                    mpq_init(v[j])
                    if i == j:
                        mpq_set(v[j], self.tmp)
                        k = k + 1
                    else:
                        mpq_set_si(v[j], 0, 1)
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
                mpq_init(v[j])
                if entries != None:
                    # TODO: If entries[k] is a rational,
                    # this should be WAY faster.  (Also see above)
                    s = str(entries[k])
                    r = mpq_set_str(v[j], s, base)
                    if r == -1:
                        _sig_off
                        raise TypeError, "Invalid rational number %s"%entries[k]
                    mpq_canonicalize(v[j])
                    k = k + 1
                else:
                    mpq_set_si(v[j],0, 1)
        _sig_off
        self.initialized = 1

    def __reduce__(self):
        import sage.matrix.dense_matrix_pyx

        cdef int i, j, len_so_far, m, n
        cdef char *a
        cdef char *s, *t, *tmp

        if self.nrows == 0 or self.ncols == 0:
            entries = ''
        else:
            n = self.nrows*self.ncols*10
            s = <char*> PyMem_Malloc(n * sizeof(char))
            t = s
            len_so_far = 0

            _sig_on
            for i from 0 <= i < self.nrows:
                for j from 0 <= j < self.ncols:
                    m = mpz_sizeinbase (mpq_numref(self.matrix[i][j]), 32) + \
                        mpz_sizeinbase (mpq_denref(self.matrix[i][j]), 32) + 3
                    if len_so_far + m + 1 >= n:
                        # copy to new string with double the size
                        n = 2*n + m + 1
                        tmp = <char*> PyMem_Malloc(n * sizeof(char))
                        strcpy(tmp, s)
                        PyMem_Free(s)
                        s = tmp
                        t = s + len_so_far
                    #endif
                    mpq_get_str(t, 32, self.matrix[i][j])
                    m = strlen(t)
                    len_so_far = len_so_far + m + 1
                    t = t + m
                    t[0] = <char>32
                    t[1] = <char>0
                    t = t + 1
            _sig_off
            entries = str(s)[:-1]
            free(s)

        return sage.matrix.dense_matrix_pyx.make_rational_matrix, \
               (self.nrows, self.ncols, entries)


    def __cmp__(self, other):
        if not isinstance(other, Matrix_rational):
            return False
        return Matrix_rational_cmp(self, other)

    def __setitem__(self, ij, x):
        i, j = ij
        if self.matrix == <mpq_t **>0:
            raise RuntimeError, "Matrix has not yet been initialized!"
        if i < 0 or i >= self.nrows or j < 0 or j >= self.ncols:
            raise IndexError, "Invalid index."
        s = str(x)
        mpq_set_str(self.matrix[i][j], s, 0)

    def __getitem__(self, ij):
        i, j = ij
        if i < 0 or i >= self.nrows or j < 0 or j >= self.ncols:
            raise IndexError, "Invalid index."
        cdef rational.Rational x
        x = rational.Rational()
        x.set_from_mpq(self.matrix[i][j])
        return x

    cdef set_matrix(Matrix_rational self, mpq_t **m):
        if self.matrix != <mpq_t **> 0:
            raise RuntimeError, "Only set matrix of uninitialized matrix."
        self.matrix = m
        self.initialized = 1

    def  __dealloc__(self):
        cdef int i, j
        if self.matrix == <mpq_t **> 0:
            return
        for i from 0 <= i < self.nrows:
            if self.matrix[i] != <mpq_t *> 0:
                for j from 0 <= j < self.ncols:
                    if self.initialized:
                        mpq_clear(self.matrix[i][j])
                PyMem_Free(self.matrix[i])
        PyMem_Free(self.matrix)

    def str(self, base=10):
        cdef int i, j
        cdef char *a
        s = "[\n"
        _sig_on
        for i from 0 <= i < self.nrows:
            for j from 0 <= j < self.ncols:
                a = mpq_get_str(NULL, base, self.matrix[i][j])
                s = s + str(a) + ", "
                free(a)
            s = s + "\n"
        s = s[:-3] + "\n]"
        _sig_off
        return s


    def __repr__(self):
        return self.str(10)


    def __mul__(Matrix_rational self, Matrix_rational other):
        if self.ncols != other.nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."

        cdef int i, j, k, nr, nc, snc
        cdef mpq_t *v
        cdef mpq_t s, z
        nr = self.nrows
        nc = other.ncols
        snc = self.ncols

        cdef Matrix_rational M
        M = Matrix_rational(nr, nc, LEAVE_UNINITIALIZED)

        cdef mpq_t **m
        m = <mpq_t **> PyMem_Malloc(sizeof(mpq_t*)*nr)
        if m == <mpq_t**> 0:
            raise MemoryError, "Error allocating matrix"

        mpq_init(s); mpq_init(z)

        _sig_on
        for i from 0 <= i < nr:
            m[i] = <mpq_t *> PyMem_Malloc(sizeof(mpq_t)*nc)
            if m[i] == <mpq_t*> 0:
                mpq_clear(s); mpq_clear(z)
                _sig_off
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                mpq_set_si(s,0,1)   # set s = 0
                v = self.matrix[i]
                for k from 0 <= k < snc:
                    mpq_mul(z, v[k], other.matrix[k][j])
                    mpq_add(s, s, z)
                mpq_init(m[i][j])
                mpq_set(m[i][j], s)
        _sig_off
        M.set_matrix(m)
        mpq_clear(s); mpq_clear(z)
        return M

    def transpose(self):
        """
        Returns the transpose of self.
        """
        cdef int i, j
        cdef Matrix_rational M

        M = Matrix_rational(self.ncols, self.nrows, entries=LEAVE_UNINITIALIZED)
        cdef mpq_t **m
        m = <mpq_t **> PyMem_Malloc(sizeof(mpq_t*)*self.ncols)
        if m == <mpq_t**> 0:
            raise MemoryError, "Error allocating matrix"

        _sig_on
        for i from 0 <= i < self.ncols:
            m[i] = <mpq_t *> PyMem_Malloc(sizeof(mpq_t)*self.nrows)
            if m[i] == <mpq_t*> 0:
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < self.nrows:
                mpq_init(m[i][j])
                mpq_set(m[i][j], self.matrix[j][i])
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
        cdef Matrix_rational M

        if not isinstance(rows, list):
            raise TypeError, "rows (=%s) must be a list"%rows
        nr = len(rows)
        if nr == 0:
            return Matrix_rational(0, self.ncols)
        nc = self.ncols
        v = []
        for i in rows:
            v.append(int(i))
        rows = v
        if min(rows) < 0 or max(rows) >= self.nrows:
            raise IndexError, "invalid row indexes; rows don't exist"

        M = Matrix_rational(nr, nc, entries=LEAVE_UNINITIALIZED)
        cdef mpq_t **m
        m = <mpq_t **> PyMem_Malloc(sizeof(mpq_t*)*nr)
        if m == <mpq_t**> 0:
            raise MemoryError, "Error allocating matrix"

        for i from 0 <= i < nr:
            m[i] = <mpq_t *> PyMem_Malloc(sizeof(mpq_t)*nc)
            if m[i] == <mpq_t*> 0:
                raise MemoryError, "Error allocating matrix"
            k = rows[i]
            for j from 0 <= j < nc:
                mpq_init(m[i][j])
                mpq_set(m[i][j], self.matrix[k][j])

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
        cdef mpq_t s, z
        nr = n
        nc = self.ncols

        if self.nrows != self.ncols:
            raise ArithmeticError, "matrix must be square"
        if not isinstance(v, list):
            raise TypeError, "v must be a list"
        if len(v) != self.nrows:
            raise ArithmeticError, "incompatible matrix vector multiple"

        cdef Matrix_rational M
        M = Matrix_rational(nr, nc, LEAVE_UNINITIALIZED)

        cdef mpq_t **m
        m = <mpq_t **> PyMem_Malloc(sizeof(mpq_t*)*nr)
        if m == <mpq_t**> 0:
            raise MemoryError, "Error allocating matrix"
        m[0] = <mpq_t *> PyMem_Malloc(sizeof(mpq_t)*nc)
        if m[0] == <mpq_t*> 0:
            mpq_clear(s); mpq_clear(z)
            raise MemoryError, "Error allocating matrix"
        mpq_init(self.tmp)
        for j from 0 <= j < nc:
            string = str(v[j])
            r = mpq_set_str(self.tmp, string, 0)
            if r == -1:
                raise TypeError, "Invalid rational number %s"%v[i]
            mpq_init(m[0][j])
            mpq_set(m[0][j], self.tmp)

        mpq_init(s)
        mpq_init(z)
        for i from 1 <= i < nr:
            m[i] = <mpq_t *> PyMem_Malloc(sizeof(mpq_t)*nc)
            if m[i] == <mpq_t*> 0:
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                mpq_set_si(s,0,1)  # set s = 0
                for k from 0 <= k < self.nrows:
                    mpq_mul(z, m[i-1][k], self.matrix[k][j])
                    mpq_add(s, s, z)
                mpq_init(m[i][j])
                mpq_set(m[i][j], s)
        M.set_matrix(m)
        mpq_clear(s); mpq_clear(z)
        return M


    def scalar_multiple(self, d):
        cdef int i, j, nr, nc
        nr = self.nrows
        nc = self.ncols

        cdef mpq_t x
        mpq_init(x)
        s = str(d)
        r = mpq_set_str(x, s, 0)
        if r == -1:
            raise TypeError, "Invalid rational number %s"%entries[k]
        cdef Matrix_rational M
        M = Matrix_rational(nr, nc, LEAVE_UNINITIALIZED)

        cdef mpq_t **m
        m = <mpq_t **> PyMem_Malloc(sizeof(mpq_t*)*nr)
        if m == <mpq_t**> 0:
            raise MemoryError, "Error allocating matrix"

        _sig_on
        for i from 0 <= i < nr:
            m[i] = <mpq_t *> PyMem_Malloc(sizeof(mpq_t)*nc)
            if m[i] == <mpq_t*> 0:
                mpq_clear(x)
                _sig_off
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                mpq_init(m[i][j])
                mpq_mul(m[i][j], self.matrix[i][j], x)
        _sig_off
        M.set_matrix(m)
        mpq_clear(x)
        return M

    def copy(self):
        cdef int i, j, nr, nc
        nr = self.nrows; nc = self.ncols

        cdef Matrix_rational M
        M = Matrix_rational(nr, nc, LEAVE_UNINITIALIZED)
        cdef mpq_t **m
        m = <mpq_t **> PyMem_Malloc(sizeof(mpq_t*)*nr)
        if m == <mpq_t**> 0:
            raise MemoryError, "Error allocating matrix"

        for i from 0 <= i < nr:
            m[i] = <mpq_t *> PyMem_Malloc(sizeof(mpq_t)*nc)
            if m[i] == <mpq_t*> 0:
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                mpq_init(m[i][j])
                mpq_set(m[i][j], self.matrix[i][j])

        M.set_matrix(m)
        return M

    def matrix_modint(self, uint p):
        cdef uint nr, nc
        nr = self.nrows; nc = self.ncols

        cdef Matrix_modint M
        M = Matrix_modint(p, nr, nc, LEAVE_UNINITIALIZED)
        cdef uint **m
        cdef uint *v
        m = <uint **> PyMem_Malloc(sizeof(uint*)*nr)
        if m == <uint **> 0:
            raise MemoryError, "Error allocating matrix"

        cdef mpq_t *w

        cdef mpz_t num, den, r
        mpz_init(num); mpz_init(den); mpz_init(r)
        cdef int numer, denom

        cdef uint i, j, inv
        _sig_on
        for i from 0 <= i < nr:
            m[i] = <uint *> PyMem_Malloc(sizeof(uint)*nc)
            w = self.matrix[i]
            v = m[i]
            for j from 0 <= j < nc:
                mpq_get_num(num, w[j])
                mpq_get_den(den, w[j])
                numer = mpz_mod_ui(r, num, p)
                denom = mpz_mod_ui(r, den, p)
                if denom % p == 0:
                    _sig_off
                    raise ZeroDivisionError, \
                          "The reduction mod %s of this matrix is not defined."%p
                v[j] = (numer * ai.c_inverse_mod_int(denom, p)) % p
        _sig_off
        M.set_matrix(m)
        mpz_clear(num); mpz_clear(den); mpz_clear(r)
        return M

    def matrix_modint_nodenom(self, uint p):
        """
        matrix_modint_nodenom(self, uint p): Return the reduction of
        this matrix modulo p, assuming that all entries of this matrix
        are integers.  For the reduction when the entries might not be
        integers, use matrix_modint(p).
        """
        cdef uint nr, nc
        nr = self.nrows; nc = self.ncols

        cdef Matrix_modint M
        M = Matrix_modint(p, nr, nc, LEAVE_UNINITIALIZED)
        cdef uint **m
        cdef uint *v
        m = <uint **> PyMem_Malloc(sizeof(uint*)*nr)
        if m == <uint **> 0:
            raise MemoryError, "Error allocating matrix"

        cdef mpq_t *w

        cdef mpz_t num, rem
        mpz_init(num); mpz_init(rem)

        cdef uint i, j, inv
        _sig_on
        for i from 0 <= i < nr:
            m[i] = <uint *> PyMem_Malloc(sizeof(uint)*nc)
            w = self.matrix[i]
            v = m[i]
            for j from 0 <= j < nc:
                mpq_get_num(num, w[j])
                v[j] = mpz_mod_ui(rem, num, p)
        _sig_off
        M.set_matrix(m)
        mpz_clear(num); mpz_clear(rem)
        return M

    def number_nonzero(self):
        cdef int i, j, n
        cdef mpq_t *v
        n = 0
        _sig_on
        for i from 0 <= i < self.nrows:
            v = self.matrix[i]
            for j from 0 <= j < self.ncols:
                if mpq_sgn(v[j]):   # if nonzero
                    n = n + 1
        _sig_off
        return n

    def list(self, int base=0):
        cdef int i, j
        cdef mpq_t *r
        cdef object v
        cdef rational.Rational x

        v = []
        _sig_on
        for i from 0 <= i < self.nrows:
            r = self.matrix[i]
            for j from 0 <= j < self.ncols:
                x = rational.Rational()
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
        m = self.matrix
        nr = self.nrows
        nc = self.ncols
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

    def hessenberg_form(self):
        if self.nrows != self.ncols:
            raise ArithmeticError, "Matrix must be square to compute Hessenberg form."

        cdef int n
        n = self.nrows

        cdef mpq_t **h
        h = self.matrix

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

    def charpoly(self):
        raise NotImplementedError

    def decompose(self):
        raise NotImplementedError

    cdef scale_row(self, int row, mpq_t multiple, int start_col):
        cdef int r
        cdef mpq_t* v

        r = row*self.ncols
        v = self.matrix[row]
        for i from start_col <= i < self.ncols:
            mpq_mul(v[i], v[i], multiple)

    cdef add_multiple_of_row(self, int row_from, mpq_t multiple,
                            int row_to, int start_col):
        cdef int i
        cdef mpq_t *v_from, *v_to
        cdef mpq_t prod, x

        mpq_init(prod); mpq_init(x)
        v_from = self.matrix[row_from]
        v_to = self.matrix[row_to]
        for i from start_col <= i < self.ncols:
            mpq_mul(prod, multiple, v_from[i])
            mpq_add(x, prod, v_to[i])
            mpq_set(v_to[i], x)   # v_to[i] <-- multipe*v_from[i] + v_to[i]

        mpq_clear(prod); mpq_clear(x)

    def set_row_to_multiple_of_row(self, int row_to, int row_from, rational.Rational multiple):
        """
        Set row row_to equal to multiple times row row_from.
        """
        cdef int i
        cdef mpq_t *v_from, *v_to

        if row_from < 0 or row_from >= self.nrows:
            raise IndexError, "row_from is %s but must be >= 0 and < %s"%(row_from, self.nrows)
        if row_to < 0 or row_to >= self.nrows:
            raise IndexError, "row_to is %s but must be >= 0 and < %s"%(row_to, self.nrows)

        v_from = self.matrix[row_from]
        v_to = self.matrix[row_to]
        for i from 0 <= i < self.ncols:
            mpq_mul(v_to[i], multiple.value, v_from[i])


    cdef add_multiple_of_column(self, int col_from, mpq_t multiple,
                               int col_to, int start_row):
        cdef int i, p, nr
        cdef mpq_t **m
        cdef mpq_t prod, x

        mpq_init(prod); mpq_init(x)
        m = self.matrix
        nr = self.nrows
        for i from start_row <= i < self.nrows:
            mpq_mul(prod, multiple, m[i][col_from])
            mpq_add(x, m[i][col_to], prod)
            mpq_set(m[i][col_to], x)
        mpq_clear(prod); mpq_clear(x)

    cdef swap_rows(self, int row1, int row2):
        cdef mpq_t* temp
        temp = self.matrix[row1]
        self.matrix[row1] = self.matrix[row2]
        self.matrix[row2] = temp

    cdef swap_columns(self, int col1, int col2):
        cdef int i, nr
        cdef mpq_t **m
        cdef mpq_t t

        mpq_init(t)
        m = self.matrix
        nr = self.nrows
        for i from 0 <= i < self.nrows:
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
        for i from 0 <= i < self.nrows:
            for j from 0 <= j < self.ncols:
                mpq_get_den(y,self.matrix[i][j])
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
        for i from 0 <= i < self.nrows:
            for j from 0 <= j < self.ncols:
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

##     def prod_of_row_sums(self, cols):
##         r"""
##         Calculate the product of all row sums of a submatrix of $A$ for a
##         list of selected columns \code{cols}.

##         This is used for the computation of matrix permanents.
##         """
##         pr = 1
##         for row in xrange(self.nrows):
##             z = 0
##             for c in cols:
##                 z = z + self[row, c]
##             pr = pr * z
##         return pr

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
            if t < 0 or t >= self.ncols:
                PyMem_Free(v)
                raise IndexError, "invalid column index (= %s)"%t
            v[c] = t

        cdef mpq_t pr, z
        mpq_init(pr)
        mpq_init(z)


        mpq_set_si(pr, 1, 1)
        for row from 0 <= row < self.nrows:
            mpq_set_si(z, 0, 1)
            for c from 0 <= c < n:
                mpq_add(z, z, self.matrix[row][v[c]])
            mpq_mul(pr, pr, z)

        cdef rational.Rational x
        x = rational.Rational()
        x.set_from_mpq(pr)
        mpq_clear(pr)
        mpq_clear(z)
        PyMem_Free(v)
        return x

    def clear_denom_copy(self):
        """
        Returns self if the denominator is 1, or n*self, where n is
        the denominator of self.
        """
        cdef mpz_t d
        mpz_init(d)
        self.mpz_denom(d)
        if mpz_cmp_si(d,1) == 0:
            return self
        A = self.copy()
        A.clear_denom()
        return A

    def clear_denom(self):
        """
        Replace self by n*self, where n is the least common multiple
        of the denominators of all entries of self.
        """
        cdef mpz_t d
        mpz_init(d)
        self.mpz_denom(d)
        if mpz_cmp_si(d,1) == 0:
            return
        cdef mpq_t denom
        mpq_init(denom)
        mpq_set_z(denom, d)
        cdef int i, j
        _sig_on
        for i from 0 <= i < self.nrows:
            for j from 0 <= j < self.ncols:
                mpq_mul(self.matrix[i][j], self.matrix[i][j], denom)
        _sig_off
        mpq_clear(denom)
        mpz_clear(d)

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
        else:
            raise ValueError, "%s is not one of the allowed algorithms (modular, gauss)"%alg

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
        B = self.clear_denom_copy()
        hA = B.height()
        if height_guess is None:
            height_guess = (2*hA)**(self.ncols/2+1)
        verbose("height_guess=%s"%height_guess)
        M = self.ncols * height_guess * hA  +  1
        p = START_PRIME
        X = []
        best_pivots = []
        prod = 1
        while True:
            while prod < M:
                verbose("p=%s"%p)
                A = B.matrix_modint(p)
                A.echelon()
                #print A   # debug
                if self.nrows == self.ncols and len(A.pivots()) == self.ncols:
                    # special case -- the echelon form must be the identity matrix.
                    return Matrix_rational_identity(self.nrows)

                c = cmp_pivots(best_pivots, A.pivots())
                if c <= 0:
                    best_pivots = A.pivots()
                    X.append(A)
                    prod = prod * p
                else:
                    if get_verbose() > 1:
                        verbose("Excluding this prime (bad pivots).", level=2)
                    pass   # do not save A since it is bad.
                p = next_prime_int(p)

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
            #print "Es = ", Es   # debug
            #print "hdE = ", hdE # debug
            #print "hA = ", hA   # debug
            #print "A = ", B     # debug
            #print "prod = ", prod  # debug
            if hdE * hA * self.ncols < prod:
                self.__pivots = best_pivots
                E._set_pivots(list(best_pivots))
                return E
            for i from 0 <= i < 3:
                M = M * START_PRIME


cdef object mpz_to_long(mpz_t x):
    return long(mpz_to_str(x))

def Matrix_rational_using_crt_and_rr(X):
    """
    Matrix_rational_using_crt_and_rr(X):

    Uses the Chinese remainder theorem to combine the matrices in the input list
    into a single matrix modulo the product of the moduli, the lifts to Q using
    rational reconstruction.  Returns a Matrix_rational on success, or raises
    a ValueError exception on failure (i.e., of the rational reconstruction
    does not exists).

    INPUT:
        X -- a nonempty list of Matrix_modint's.  The moduli must be
                    distinct primes.

    OUTPUT:
        a single Matrix_rational, which reduces to each of the
        input Matrix_modint's.
    """
    if not isinstance(X, list):
        raise TypeError, "Argument 1 must be a list."
    if len(X) < 1:
        raise ValueError, "The list of matrices must be nonempty."
    for A in X:
        if not isinstance(A, Matrix_modint):
            raise TypeError, "Each argument in the input list must be a Matrix_modint"

    cdef int i, j, k, n, nr, nc
    n = len(X)
    nr = X[0].numrows()
    nc = X[0].numcols()

    for i from 1 <= i < n:
        if X[i].numrows() != nr or X[i].numcols() != nc:
            raise TypeError, "All matrices in the input must have the same dimensions"

    cdef uint *p, *v
    p = <uint*> PyMem_Malloc(sizeof(uint)*n)
    v = <uint*> PyMem_Malloc(sizeof(uint)*n)
    for i from 0 <= i < n:
        p[i] = X[i].prime()

    cdef mpq_t **m

    m = <mpq_t **> PyMem_Malloc(sizeof(mpq_t)*nr)
    if m == <mpq_t**> 0:
        raise MemoryError, "Error allocating matrix"

    # We do the following, so that in the triply-nested for
    # loops below, no access are to Python objects -- everything
    # is in pure C.
    cdef uint ***matrices
    matrices = <uint ***> PyMem_Malloc(sizeof(uint**)*n)
    for i from 0 <= i < n:
        matrices[i] = Matrix_modint_matrix(X[i])

    for i from 0 <= i < nr:
        _sig_check
        m[i] = <mpq_t *> PyMem_Malloc(sizeof(mpq_t)*nc)
        if m[i] == <mpq_t*> 0:
            raise MemoryError, "Error allocating matrix"
        for j from 0 <= j < nc:
            for k from 0 <= k < n:
                v[k] = matrices[k][i][j]
            mpq_init(m[i][j])
            mpq_using_crt_and_rr(m[i][j], v, p, n)
    PyMem_Free(p); PyMem_Free(v);
    PyMem_Free(matrices)

    cdef Matrix_rational M
    M = Matrix_rational(nr, nc, LEAVE_UNINITIALIZED)
    M.set_matrix(m)
    #print "X = ", X  # debug
    #print "M = ", M  # debug
    return M

def Matrix_rational_random(nrows, ncols, bound):
    x = []
    for i in range(nrows*ncols):
        x.append(gmp_randrange(-bound, bound))
    return Matrix_rational(nrows, ncols, x)

def Matrix_rational_identity(n):
    x = []
    new_row = True
    for i in range(n):
        x.append(1)
        for j in range(n):
            x.append(0)
    I = Matrix_rational(n,n, x[:n*n])
    I._set_pivots(range(n))
    return I
