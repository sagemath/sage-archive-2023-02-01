r"""nodoctest
Sparse matrices over $\Q$

This is an optimzed implementation of sparse matrices over
the rational numbers.

AUTHOR:
   -- William Stein (2004): first version
   -- William Stein (2006-02-12): added set_row_to_multiple_of_row
   -- William Stein (2006-03-04): added multimodular echelon, __reduce__, etc.
   -- William Stein (2007-02-20): update to new SAGE matrix class structure.

EXAMPLES:
    sage: M = MatrixSpace(QQ, 2, 3, sparse=True)
    sage: A = M([1,2,3, 1,1,1])
    sage: A
    [1 2 3]
    [1 1 1]
    sage: A.echelon_form()
    [ 1  0 -1]
    [ 0  1  2]

    sage: M = MatrixSpace(QQ, 1000,1000, sparse=True)
    sage: A = M(0)
    sage: A[1,1] = 5


    sage: from sage.matrix.sparse_matrix import SparseMatrix
    sage: x = SparseMatrix(QQ, 5,10)
    sage: x.randomize(5)
    sage: x.echelon_form()       # random output
    [
    1, 0, 0, 0, 0, 0, 0, 0, 0, 10/29,
    0, 1, 0, 0, 0, 0, 0, 0, 0, -4/29,
    0, 0, 1, 0, 0, 0, 0, 0, 0, -12/29,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 24/29,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 4/29
    ]

"""

#############################################################################
#       Copyright (C) 2004, 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

import random

import sage.misc.all
import sage.matrix.matrix_rational_sparse
import sage.rings.integer_ring
import sage.rings.finite_field
import sage.rings.arith

cimport sage.rings.rational
import  sage.rings.rational

cimport sage.rings.integer
import  sage.rings.integer

cimport sage.ext.arith
import sage.ext.arith
cdef sage.ext.arith.arith_int ai
ai = sage.ext.arith.arith_int()

cimport matrix_field_sparse
import matrix_field_sparse

cimport matrix_modn_sparse
import matrix_modn_sparse

cimport matrix_rational_dense
import matrix_rational_dense

include "../ext/gmp.pxi"
include '../ext/interrupt.pxi'

START_PRIME = 20011  # used for multi-modular algorithms

cdef class Vector_mpq

cdef void Vector_mpq_rescale(Vector_mpq w, mpq_t x):
    scale_mpq_vector(&w.v, x)

cdef class Vector_mpq:
    """
    Vector_mpq -- a sparse vector of GMP rationals.  This is a Python
    extension type that wraps the C implementation of sparse vectors
    modulo a small prime.
    """
    cdef mpq_vector v

    def __init__(self, int degree, int num_nonzero=0, entries=[], sort=True):
        cdef int i
        init_mpq_vector(&self.v, degree, num_nonzero)
        if entries != []:
            if len(entries) != num_nonzero:
                raise ValueError, "length of entries (=%s) must equal num_nonzero (=%s)"%(len(entries), num_nonzero)
            if sort:
                entries = list(entries) # copy so as not to modify passed in list
                entries.sort()
            for i from 0 <= i < num_nonzero:
                s = str(entries[i][1])
                mpq_set_str(self.v.entries[i], s, 0)
                self.v.positions[i] = entries[i][0]

    def __dealloc__(self):
        clear_mpq_vector(&self.v)

    def __getitem__(self, int n):
        cdef mpq_t x
        cdef sage.rings.rational.Rational a
        mpq_init(x)
        mpq_vector_get_entry(&x, &self.v, n)
        a = sage.rings.rational.Rational()
        a.set_from_mpq(x)
        mpq_clear(x)
        return a

    def cmp(self, Vector_mpq other):
        return mpq_vector_cmp(&self.v, &other.v)

    def __richcmp__(Vector_mpq self, x, int op):
        if not isinstance(x, Vector_mpq):
            return -1
        cdef int n
        n = self.cmp(x)
        if op == 0:
            return bool(n < 0)
        elif op == 1:
            return bool(n <= 0)
        elif op == 2:
            return bool(n == 0)
        elif op == 3:
            return bool(n != 0)
        elif op == 4:
            return bool(n > 0)
        elif op == 5:
            return bool(n >= 0)

    def __setitem__(self, int n, x):
        cdef object s
        s = str(x)
        mpq_vector_set_entry_str(&self.v, n, s)

    def __repr__(self):
        return str(list(self))

    def degree(self):
        return self.v.degree

    def num_nonzero(self):
        return self.v.num_nonzero

    def list(self):
        return mpq_vector_to_list(&self.v)

    cdef void rescale(self, mpq_t x):
        scale_mpq_vector(&self.v, x)

    def __add__(Vector_mpq self, Vector_mpq other):
        cdef mpq_vector z1, *z2
        cdef Vector_mpq w
        cdef mpq_t ONE
        mpq_init(ONE)
        mpq_set_si(ONE,1,1)

        add_mpq_vector_init(&z1, &self.v, &other.v, ONE)
        mpq_clear(ONE)
        w = Vector_mpq(self.v.degree)
        z2 = &(w.v)
        clear_mpq_vector(z2)   # free memory wasted on allocated w
        z2.entries = z1.entries
        z2.positions = z1.positions
        z2.num_nonzero = z1.num_nonzero
        # at this point we do *not* free z1, since it is referenced by w.
        return w

    def __sub__(Vector_mpq self, Vector_mpq other):
        return self + other*(-1)

    def copy(self):
        cdef int i
        cdef Vector_mpq w
        w = Vector_mpq(self.v.degree, self.v.num_nonzero)
        for i from 0 <= i < self.v.num_nonzero:
            mpq_set(w.v.entries[i], self.v.entries[i])
            w.v.positions[i] = self.v.positions[i]
        return w

    def __mul__(x, y):
        if isinstance(x, Vector_mpq):
            self = x
            other = y
        elif isinstance(y, Vector_mpq):
            self = y
            other = x
        else:
            raise TypeError, "Invalid types."
        cdef object s, z
        cdef mpq_t t
        z = self.copy()
        mpq_init(t)
        s = str(other)
        mpq_set_str(t, s, 0)
        Vector_mpq_rescale(z, t)
        mpq_clear(t)
        return z

    def randomize(self, int sparcity, bound=3):
        """
        randomize(self, int sparcity, exact=False):

        The sparcity is a bound on the number of nonzeros per row.
        """
        cdef int i
        for i from 0 <= i < sparcity:
            self[random.randrange(self.v.degree)] = random.randrange(1,bound)


#############################################################
#
#    Sparse Matrix over mpq_t (the GMP rationals)
#
#############################################################
cdef class Matrix_rational_sparse(matrix_sparse.Matrix_sparse):

    def __new__(self, int nrows, int ncols, object entries=[], init=True, coerce=False):
        # allocate memory
        cdef int i
        self.rows = <mpq_vector*> sage_malloc(nrows*sizeof(mpq_vector))
        self.is_init = init
        if self.is_init:
            self.is_init = True
            for i from 0 <= i < nrows:
                init_mpq_vector(&self.rows[i], ncols, 0)

    def __dealloc__(self):
        if not self.is_init:
            return
        cdef int i
        for i from 0 <= i < self.nr:
            clear_mpq_vector(&self.rows[i])

    def __init__(self, int nrows, int ncols, object entries=[], init=True, coerce=False):
        """
        INPUT:
            nrows -- number of rows
            ncols -- number of columns
            entries -- list of triples (i,j,x), where 0 <= i < nrows,
                       0 <= j < ncols, and x is a rational number.
                       Then the i,j entry of the matrix is set to x.
                       It is OK for some x to be zero.
                or, list of all entries of the matrix (dense representation).
            init -- bool (default: True); if False, don't allocate anything (for
                    internal use only!)
            coerce -- bool (default: False), if True, entries might not be of type Rational.
        """
        cdef object s
        cdef int ii, jj, k, n
        cdef sage.rings.rational.Rational z

        self.nr = nrows
        self.nc = ncols
        self.__pivots = None

        if isinstance(entries, str):
            v = entries.split('\n')
            for ii from 0 <= ii < nrows:
                w = v[ii].split()
                n = int(w[0])
                init_mpq_vector(&self.rows[ii], ncols, n)
                for jj from 0 <= jj < n:
                    self.rows[ii].positions[jj] = int(w[jj+1])
                    t = w[jj+n+1]
                    mpq_set_str(self.rows[ii].entries[jj], t, 32)
            self.is_init = True
            return

        if not init:
            return

        if not coerce:
            if isinstance(entries, list):
                if len(entries) == 0:
                    return
                if not isinstance(entries[0], tuple):
                    # dense input representation
                    k = 0
                    for ii from 0 <= ii < nrows:
                        _sig_check
                        for jj from 0 <= jj < ncols:
                            z = entries[k]
                            if mpq_sgn(z.value):         # if z is nonzero
                                mpq_vector_set_entry(&self.rows[ii], jj, z.value)
                            k = k + 1
                else:
                    # sparse input rep
                    for i, j, x in entries:
                        z = x
                        if mpq_sgn(z.value):         # if z is nonzero
                            mpq_vector_set_entry(&self.rows[i], j, z.value)
                return

            if isinstance(entries, dict):
                for (i,j), x in entries.iteritems():
                    z = x
                    mpq_vector_set_entry(&self.rows[i], j, z.value)
                return

        else:
            # copying the above code is very ugly.  but this
            # function must be very fast...
            if isinstance(entries, list):
                if len(entries) == 0:
                    return
                if not isinstance(entries[0], tuple):
                    # dense input representation
                    k = 0
                    for ii from 0 <= ii < nrows:
                        _sig_check
                        for jj from 0 <= jj < ncols:
                            z = sage.rings.rational.Rational(entries[k])
                            if mpq_sgn(z.value):         # if z is nonzero
                                mpq_vector_set_entry(&self.rows[ii], jj, z.value)
                            k = k + 1
                else:
                    # one sparse rep
                    for i, j, x in entries:
                        z = sage.rings.rational.Rational(x)
                        if mpq_sgn(z.value):         # if z is nonzero
                            mpq_vector_set_entry(&self.rows[i], j, z.value)
                return

            if isinstance(entries, dict):
                for (i,j), x in entries.iteritems():
                    z = sage.rings.rational.Rational(x)
                    mpq_vector_set_entry(&self.rows[i], j, z.value)
                return

        # Now assume entries is a rationl number and matrix should be scalar
        if nrows == ncols:
            x = sage.rings.rational.Rational(entries)
            for ii from 0 <= ii < nrows:
                self[ii,ii] = x
            return

        raise TypeError, "no way to make matrix from entries (=%s)"%entries


    def __cmp__(self, Matrix_rational_sparse other):
        if self.nr != other.nr:
            return -1
        if self.nc != other.nc:
            return -1
        cdef int i, j, c
        for i from 0 <= i < self.nr:
            if self.rows[i].num_nonzero != other.rows[i].num_nonzero:
                return 1
            for j from 0 <= j < self.rows[i].num_nonzero:
                c = mpq_cmp(self.rows[i].entries[j], other.rows[i].entries[j])
                if c:
                    return c
        return 0

    def __reduce__(self):
        """
        EXAMPLES:
            sage: A = MatrixSpace(QQ, 2, sparse=True)([1,2,-1/5,19/397])
            sage: loads(dumps(A)) == A
            True
        """
        # Format of reduced serialized representation of a sparse matrix
        #     number_nonzer_entries_in_row  nonzero_positions  nonzero_entries_in_base32[newline]
        #
        #     3  1 4 7 3/4 7/83  39/45
        #     2  1 193 -17 38
        #

        cdef char *s, *t, *tmp
        cdef int m, n, ln, i, j, z, len_so_far

        n = self.nr * 200 + 30
        s = <char*> sage_malloc(n * sizeof(char))
        len_so_far = 0
        t = s
        s[0] = <char>0   # make s a null-terminated string
        for i from 0 <= i < self.nr:
            ln = self.rows[i].num_nonzero
            if len_so_far + 20*ln >= n:
                # copy to new string with double the size
                n = 2*n + 20*ln
                tmp = <char*> sage_malloc(n * sizeof(char))
                strcpy(tmp, s)
                sage_free(s)
                s = tmp
                t = s + len_so_far
            #endif
            z = sprintf(t, '%d ', ln)
            t = t + z
            len_so_far = len_so_far + z
            for j from 0 <= j < ln:
                z = sprintf(t, '%d ', self.rows[i].positions[j])
                t = t + z
                len_so_far = len_so_far + z

            for j from 0 <= j < ln:
                m = mpz_sizeinbase(mpq_numref(self.rows[i].entries[j]), 32) + \
                    mpz_sizeinbase(mpq_denref(self.rows[i].entries[j]), 32) + 3
                if len_so_far + m >= n:
                    # copy to new string with double the size
                    n = 2*n + m + 1
                    tmp = <char*> sage_malloc(n)
                    strcpy(tmp, s)
                    sage_free(s)
                    s = tmp
                    t = s + len_so_far
                mpq_get_str(t, 32, self.rows[i].entries[j])
                m = strlen(t)
                len_so_far = len_so_far + m + 1
                t = t + m
                if j <= ln-1:
                    t[0] = <char>32
                    t[1] = <char>0
                    t = t + 1

            z = sprintf(t, '\n')
            t = t + z
            len_so_far = len_so_far + z
        # end for

        entries = str(s)[:-1]
        sage_free(s)
        return make_sparse_rational_matrix, (self.nr, self.nc, entries)


    def copy(self):
        """
        Return a copy of this matrix.

        EXAMPLES:
            sage: A = MatrixSpace(QQ, 2, sparse=True)([1,2,1/3,-4/39])._sparse_matrix_mpq_()
            sage: A.copy()
            [
            1, 2,
            1/3, -4/39
            ]
        """
        cdef Matrix_rational_sparse A
        A = Matrix_rational_sparse(self.nr, self.nc, init=False)
        cdef int i, j, k
        for i from 0 <= i < self.nr:
            # _sig_check      # not worth it since whole function is so fast?
            k = self.rows[i].num_nonzero
            if init_mpq_vector(&A.rows[i], self.nc, k) == -1:
                raise MemoryError, "Error allocating memory"
            for j from 0 <= j < k:
                mpq_set(A.rows[i].entries[j], self.rows[i].entries[j])
                A.rows[i].positions[j] = self.rows[i].positions[j]
        A.is_init = True
        return A

    def matrix_from_rows(self, rows):
        """
        Return the matrix obtained from the rows of self given by the
        input list of integers.

        INPUT:
            rows -- list of integers

        EXAMPLES:
            sage: M = MatrixSpace(QQ,2,2,sparse=True)(range(1,5))
            sage: m = M._Matrix_sparse_rational__matrix
            sage: s = m.matrix_from_rows([1,1,1])
            sage: s
            [
            3, 4,
            3, 4,
            3, 4
            ]
        """
        cdef int i, j, k, nr
        nr = len(rows)

        cdef Matrix_rational_sparse A
        A = Matrix_rational_sparse(nr, self.nc, init=False)

        for i from 0 <= i < nr:
            j = rows[i]
            init_mpq_vector(&A.rows[i], self.nc, self.rows[j].num_nonzero)
            for k from 0 <= k < self.rows[j].num_nonzero:
                mpq_set(A.rows[i].entries[k], self.rows[j].entries[k])
                A.rows[i].positions[k] = self.rows[j].positions[k]

        A.is_init = True

        return A

    def dense_matrix(self):
        """
        Return corresponding dense matrix.
        """
        cdef matrix_rational_dense.Matrix_rational_dense A
        A = matrix_rational_dense.Matrix_rational_dense(self.nr, self.nc)
        # now A is the zero matrix.  We fill in the entries.
        cdef int i, j, k
        for i from 0 <= i < self.nr:
            k = self.rows[i].num_nonzero
            for j from 0 <= j < k:
                # now set A[i,self.rows[i].positions[j]] equal to self.rows[i].entries[j]
                mpq_set(A._matrix[i][self.rows[i].positions[j]], self.rows[i].entries[j])
        return A

    def linear_combination_of_rows(self, Vector_mpq v):
        if self.nr != v.degree():
            raise ArithmeticError, "Incompatible vector * matrix multiply."
        cdef mpq_vector w, sum, sum2
        cdef int i, r, nr
        cdef Vector_mpq ans
        nr = self.nr
        w = v.v
        init_mpq_vector(&sum, self.nc, 0)
        _sig_on
        for i from 0 <= i < w.num_nonzero:
            r = w.positions[i]
            add_mpq_vector_init(&sum2, &sum, &self.rows[r], w.entries[i])
            # Now sum2 is initialized and equals sum + w[i]*self.rows[i]
            # We want sum to equal this.
            clear_mpq_vector(&sum)
            sum = sum2
        _sig_off
        # Now sum is a sparse C-vector that gives the linear combination of rows.
        # Convert to a Vector_mpq and return.
        ans = Vector_mpq(nr)
        clear_mpq_vector(&ans.v)
        ans.v = sum
        return ans

    def set_row_to_multiple_of_row(self, int row_to, int row_from,
                                   sage.rings.rational.Rational multiple):
        """
        Set row row_to equal to multiple times row row_from.

        EXAMPLES:
            sage: m = Matrix_rational_sparse(3,3,entries=[(1,1,10/3)])
            sage: m
            [
            0, 0, 0,
            0, 10/3, 0,
            0, 0, 0
            ]
            sage: m.set_row_to_multiple_of_row(0, 1, 6/1)   # third argument must be a rational!
            sage: m
            [
            0, 20, 0,
            0, 10/3, 0,
            0, 0, 0
            ]
            sage: m.set_row_to_multiple_of_row(2,1,-10/1)
            sage: m
            [
            0, 20, 0,
            0, 10/3, 0,
            0, -100/3, 0
            ]
            sage: m = Matrix_rational_sparse(3,3,entries=[(1,1,10/3)])
            sage: m.set_row_to_multiple_of_row(1,1,-10/1); m
            [
            0, 0, 0,
            0, -100/3, 0,
            0, 0, 0
            ]
            sage: m.set_row_to_multiple_of_row(-1, 1, 6/1)
            Traceback (most recent call last):
            ...
            IndexError: row_to is -1 but must be >= 0 and < 3
            sage: m.set_row_to_multiple_of_row(0, 3, 6/1)
            Traceback (most recent call last):
            ...
            IndexError: row_from is 3 but must be >= 0 and < 3
        """
        # A sparse matrix is an array of pointers to mpq_vector's.
        # 1. Delete the vector in position row_to
        # 2. Initialize a new one in its place.
        # 3. Fill in the entries with appropriate multiples of the entries in row_from.
        cdef int i

        if row_from < 0 or row_from >= self.nr:
            raise IndexError, "row_from is %s but must be >= 0 and < %s"%(row_from, self.nr)
        if row_to < 0 or row_to >= self.nr:
            raise IndexError, "row_to is %s but must be >= 0 and < %s"%(row_to, self.nr)

        if row_from == row_to:
            scale_mpq_vector(&self.rows[row_from], multiple.value)
            return

        clear_mpq_vector(&self.rows[row_to])
        init_mpq_vector(&self.rows[row_to], self.nc, self.rows[row_from].num_nonzero)

        for i from 0 <= i < self.rows[row_from].num_nonzero:
            mpq_mul(self.rows[row_to].entries[i], multiple.value, self.rows[row_from].entries[i])
            self.rows[row_to].positions[i] = self.rows[row_from].positions[i]


    def set_row_to_negative_of_row_of_A_using_subset_of_columns(self, int i, Matrix_rational_sparse A, int r, cols):
        # the ints cols are assumed sorted.
        # this function exists just because it is useful for modular symbols presentations.
        cdef int l
        cdef mpq_t x
        mpq_init(x)
        l = 0
        for k in cols:
            mpq_vector_get_entry(&x, &A.rows[r], k)
            if mpz_cmp_si(mpq_numref(x), 0):      # x is nonzero
                mpz_mul_si(mpq_numref(x), mpq_numref(x), -1)
                mpq_vector_set_entry(&self.rows[i], l, x)
            l = l + 1
        mpq_clear(x)

    def parent(self):
        import sage.matrix.matrix_space
        import sage.rings.rings
        return sage.matrix.matrix_space.MatrixSpace(
                  sage.rings.rings.RationalField(),self.nr,self.nc)

    def pivots(self):
        if self.__pivots is None:
            raise NotImplementedError
        return self.__pivots

    def row_to_dict(self, int i):
        """
        Return an associative arrow of pairs
               n:x
        where the keys n run through the nonzero positions of the row,
        and the x are nonzero and of type Integer.
        """
        cdef int j, n
        cdef sage.rings.rational.Rational x
        cdef object entries
        if i < 0 or i >= self.nr: raise IndexError
        X = {}
        for j from 0 <= j < self.rows[i].num_nonzero:
            n = self.rows[i].positions[j]
            x = sage.rings.rational.Rational()
            x.set_from_mpq(self.rows[i].entries[j])
            X[n] = x
        return X

    def dict(self):
        """
        Return an associative arrow of pairs
               (i,j):x
        where the keys (i,j) run through the nonzero positions of the matrix
        and the x are nonzero and of type Integer.
        """
        cdef int i, j, n
        cdef sage.rings.rational.Rational x
        cdef mpq_t t
        cdef object entries
        X = {}
        for i from 0 <= i < self.nr:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            for j from 0 <= j < self.rows[i].num_nonzero:
                n = self.rows[i].positions[j]
                x = sage.rings.rational.Rational()
                x.set_from_mpq(self.rows[i].entries[j])
                X[(i,n)] = x
        return X

    def randomize(self, int sparcity, bound=2):
        """
        randomize(self, int sparcity):

        The sparcity is a bound on the number of nonzeros per row.
        """
        cdef int i, j, k, r
        for i from 0 <= i < self.nr:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            for j from 0 <= j <= sparcity:
                self[i, random.randrange(0,self.nc)] = random.randrange(-bound,bound)

    def __repr__(self):
        cdef int i, j
        cdef mpq_t x
        cdef char *buf

        mpq_init(x)
        s = "[\n"
        for i from 0 <= i < self.nr:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            for j from 0 <= j < self.nc:
                mpq_vector_get_entry(&x, &self.rows[i], j)
                buf = mpq_get_str(NULL, 10, x)
                s = s + str(buf) + ", "
                sage_free(buf)   # use c's malloc/free
            s = s + "\n"
        s = s[:-3] + "\n]"
        mpq_clear(x)
        return s

    def list(self):
        cdef int i
        X = []
        for i from 0 <= i < self.nr:
            for j, x in mpq_vector_to_list(&self.rows[i]):
                X.append((i,j,x))
        return X

    def __getitem__(self, t):
        if not isinstance(t, tuple) or len(t) != 2:
            raise IndexError, "Index (=%s) of matrix access must be a row and a column."%t
        cdef sage.rings.rational.Rational y
        cdef int i, j
        i, j = t
        cdef mpq_t x
        mpq_init(x)
        mpq_vector_get_entry(&x, &self.rows[i], j)
        y = sage.rings.rational.Rational()
        y.set_from_mpq(x)
        mpq_clear(x)
        return y

    def __setitem__(self, t, x):
        if not isinstance(t, tuple) or len(t) != 2:
            raise IndexError, "index (=%s) for setting matrix item must be a 2-tuple."%t
        cdef int i, j
        i, j = t
        if i<0 or i >= self.nr or j<0 or j >= self.nc:
            raise IndexError, "Array index out of bounds."
        cdef object s
        s = str(x)
        mpq_vector_set_entry_str(&self.rows[i], j, s)

    def matrix_multiply(self, Matrix_rational_sparse B):
        """
        Return the matrix product of self and B.
        """
        cdef int i, j, k

        cdef mpq_t x, t
        mpq_init(x)
        mpq_init(t)

        cdef Matrix_rational_sparse A
        A = Matrix_rational_sparse(self.nr, B.nc)

        for i from 0 <= i < self.nr:
            if sage.misc.all.get_verbose()>=3 and i%50 == 0:
                sage.misc.all.verbose('row %s of %s'%(i, self.nr), level=3)
            for j from 0 <= j < B.nc:
                # dot of i-th row with j-th column
                mpq_set_si(x, 0, 1)
                for k from 0 <= k < self.rows[i].num_nonzero:
                    mpq_vector_get_entry(&t, &B.rows[self.rows[i].positions[k]],  j)
                    if mpz_cmp_si(mpq_numref(t), 0) != 0:  # is nonzero
                        mpq_mul(t, t, self.rows[i].entries[k])
                        mpq_add(x, x, t)
                if mpz_cmp_si(mpq_numref(x), 0) != 0:
                    mpq_vector_set_entry(&A.rows[i], j, x)

        mpq_clear(x)
        mpq_clear(t)
        return A


    def nrows(self):
        return self.nr

    def ncols(self):
        return self.nc

    def matrix_modint(self, int n, denoms=True):
        """
        Return reduction of this matrix modulo the integer $n$.

        INPUT:
            n -- int
            denoms -- bool (default: True) if True reduce denominators;
                      if False assume all denoms are 1.
        """
        cdef int i, j, d
        cdef matrix_modn_sparse.Matrix_modn_sparse A
        #cdef matrix_sparse.Matrix_sparse A
        cdef unsigned int num, den
        cdef mpq_vector* v
        d = denoms

        A = matrix_modn_sparse.Matrix_modn_sparse(MatrixSpace(GF(n), self.nr, self.nc), n, self.nr, self.nc)
        for i from 0 <= i < self.nr:
            v = &self.rows[i]
            for j from 0 <= j < v.num_nonzero:
                if mpz_cmp_si(mpq_denref(v.entries[j]), 1) == 0:
#matrix_rational_sparse.pyx:1480:30: Cannot assign type 'c_vector_modint (*)' to 'c_vector_modint (*)'  ?!
                    set_entry(<c_vector_modint *>&A.rows[i], v.positions[j],
                              mpz_fdiv_ui(mpq_numref(v.entries[j]), n))
                else:
                    num = mpz_fdiv_ui(mpq_numref(v.entries[j]), n)
                    if denoms:
                        den = mpz_fdiv_ui(mpq_denref(v.entries[j]), n)
                        set_entry(<c_vector_modint *>&A.rows[i], v.positions[j],
                                  int((num * ai.inverse_mod_int(den, n)) % n))
                    else:
                        set_entry(<c_vector_modint *>&A.rows[i], v.positions[j], num)

        return A

    def swap_rows(self, int n1, int n2):
        """
        Swap the rows in positions n1 and n2
        """
        if n1 < 0 or n1 >= self.nr or n2 < 0 or n2 >= self.nr:
            raise IndexError, "Invalid row number (n1=%s, n2=%s)"%(n1,n2)
        if n1 == n2:
            return
        cdef mpq_vector tmp
        tmp = self.rows[n1]
        self.rows[n1] = self.rows[n2]
        self.rows[n2] = tmp



    def height(self, scale=1):
        """
        Returns the height of scale*self, which is the maximum of the
        absolute values of all numerators and denominators of the
        entries of scale*self.

        OUTPUT:
            -- Integer

        NOTE: Since 0 = 0/1 has denominator 1, the height is at least 1.

        EXAMPLES:
            sage: A = MatrixSpace(QQ, 3,3, sparse=True)([1,2,3,4,5/13,6,-7/17,8,9])._sparse_matrix_mpq_()
            sage: A.height()
            17
            sage: A = MatrixSpace(QQ, 2,2, sparse=True)([1,2,-197/13,4])._sparse_matrix_mpq_()
            sage: A.height()
            197
            sage: A.height(26)
            394
        """
        cdef mpz_t h, a
        mpz_init(h)
        mpz_init(a)
        mpz_set_si(h, 1)

        cdef mpq_t s, v2
        cdef sage.rings.rational.Rational tmp
        tmp = sage.rings.rational.Rational(abs(scale))
        mpq_init(s)
        mpq_set(s,tmp.value)

        cdef int i, j
        cdef mpq_vector* v

        if scale == 1:
            for i from 0 <= i < self.nr:
                v = &self.rows[i]
                for j from 0 <= j < v.num_nonzero:
                    mpz_abs(a, mpq_numref(v.entries[j]))
                    if mpz_cmp(a, h) > 0:
                        mpz_set(h, a)
                    if mpz_cmp(mpq_denref(v.entries[j]), h) > 0:
                        mpz_set(h, mpq_denref(v.entries[j]))
        else:
            mpq_init(v2)
            for i from 0 <= i < self.nr:
                v = &self.rows[i]
                for j from 0 <= j < v.num_nonzero:
                    mpq_mul(v2, v.entries[j], s)
                    mpz_abs(a, mpq_numref(v2))
                    if mpz_cmp(a, h) > 0:
                        mpz_set(h, a)
                    if mpz_cmp(mpq_denref(v2), h) > 0:
                        mpz_set(h, mpq_denref(v2))
            mpq_clear(v2)

        #endif
        mpz_clear(a)
        mpq_clear(s)

        cdef sage.rings.integer.Integer r
        r = sage.rings.integer.Integer()
        r.set_from_mpz(h)
        mpz_clear(h)
        return r

    def denom(self):
        """
        Returns the denominator of self, which is the least common
        multiple of the denominators of the entries of self.

        OUTPUT:
            -- Integer

        NOTE: Since 0 = 0/1 has denominator 1, the height is at least 1.

        EXAMPLES:
            sage: A = MatrixSpace(QQ, 3,3, sparse=True)([1/17,2,3,4,5/13,6,-7/17,8,9])._sparse_matrix_mpq_()
            sage: A.denom()
            221
            sage: A = MatrixSpace(QQ, 2,2, sparse=True)([1,2,-197/13,4])._sparse_matrix_mpq_()
            sage: A.denom()
            13
        """
        cdef mpz_t d
        mpz_init(d)
        mpz_set_si(d, 1)

        cdef int i, j
        cdef mpq_vector* v

        _sig_on
        for i from 0 <= i < self.nr:
            v = &self.rows[i]
            for j from 0 <= j < v.num_nonzero:
                mpz_lcm(d, d, mpq_denref(v.entries[j]))
        _sig_off

        cdef sage.rings.integer.Integer _d
        _d = sage.rings.integer.Integer()
        _d.set_from_mpz(d)
        mpz_clear(d)
        return _d

    def clear_denom(self):
        """
        Replace self by self times the denominator of self.

        EXAMPLES:
            sage: A = MatrixSpace(QQ,2, sparse=True)([1/3,2,3/2,4])._sparse_matrix_mpq_()
            sage: A.clear_denom ()
            6
            sage: A
            [
            2, 12,
            9, 24
            ]
        """
        cdef sage.rings.integer.Integer _d
        _d = self.denom()

        cdef mpq_t d
        mpq_init(d)
        mpz_set(mpq_numref(d), _d.value)

        cdef int i, j
        cdef mpq_vector* v

        _sig_on
        for i from 0 <= i < self.nr:
            v = &self.rows[i]
            for j from 0 <= j < v.num_nonzero:
                mpq_mul(v.entries[j], v.entries[j], d)
        _sig_off

        mpq_clear(d)
        return _d


    def divide_by(self, sage.rings.integer.Integer d):
        """
        Replace self by self/d.

        EXAMPLES:
            sage: A = MatrixSpace(QQ,2, sparse=True)([1,2,3,4])._sparse_matrix_mpq_()
            sage: A.divide_by(5)
            sage: A
            [
            1/5, 2/5,
            3/5, 4/5
            ]
        """
        cdef mpq_t dd
        mpq_init(dd)
        mpq_set_si(dd, 1, 1)
        mpz_set(mpq_denref(dd), d.value)

        cdef int i, j
        cdef mpq_vector* v

        _sig_on
        for i from 0 <= i < self.nr:
            v = &self.rows[i]
            for j from 0 <= j < v.num_nonzero:
                mpq_mul(v.entries[j], v.entries[j], dd)
        _sig_off
        mpq_clear(dd)

    def echelon_multimodular(self, height_guess=None, proof=True):
        """
        Returns reduced row-echelon form using a multi-modular
        algorithm.  Does not change self.

        INPUT:
            height_guess -- integer or None
            proof -- boolean (default: True)

        EXAMPLES:
            sage: A = MatrixSpace(QQ,2, sparse=True)([1/3,2,3/2,4])._sparse_matrix_mpq_(); A
            [
            1/3, 2,
            3/2, 4
            ]
            sage: B = A.echelon_multimodular(); B
            [
            1, 0,
            0, 1
            ]

        Note that A is unchanged.
            sage: A
            [
            1/3, 2,
            3/2, 4
            ]


        ALGORITHM:
        The following is a modular algorithm for computing the echelon
        form.  Define the height of a matrix to be the max of the
        absolute values of the entries.

        Given Matrix A with n columns (self).

         0. Rescale input matrix A to have integer entries.  This does
            not change echelon form and makes reduction modulo lots of
            primes significantly easier if there were denominators.
            Henceforth we assume A has integer entries.

         1. Let c be a guess for the height of the echelon form.  E.g.,
            c=1000, e.g., if matrix is very sparse and application is to
            computing modular symbols.

         2. Let M = n * c * H(A) + 1,
            where n is the number of columns of A.

         3. List primes p_1, p_2, ..., such that the product of
            the p_i is at least M.

         4. Try to compute the rational reconstruction CRT echelon form
            of A mod the product of the p_i.  If rational
            reconstruction fails, compute 1 more echelon forms mod the
            next prime, and attempt again.  Make sure to keep the
            result of CRT on the primes from before, so we don't have
            to do that computation again.  Let E be this matrix.

         5. Compute the denominator d of E.
            Attempt to prove that result is correct by checking that

                  H(d*E)*ncols(A)*H(A) < (prod of reduction primes)

            where H denotes the height.   If this fails, do step 4 with
            a few more primes.
        """
        cdef Matrix_rational_sparse E

        if self.nr == 0 or self.nc == 0:
            return self

        cdef sage.rings.integer.Integer dd
        dd = self.clear_denom()
        #print "dd = ", dd

        hA = long(self.height())
        if height_guess is None:
            height_guess = 100000*hA**4
        tm = sage.misc.all.verbose("height_guess = %s"%height_guess, level=2)

        if proof:
            M = self.nc * height_guess * hA  +  1
        else:
            M = height_guess + 1

        p = START_PRIME
        X = []
        best_pivots = []
        prod = 1
        problem = 0
        while True:
            while prod < M:
                problem = problem + 1
                if problem > 50:
                    sage.misc.all.verbose("sparse_matrix multi-modular reduce not converging?")
                t = sage.misc.all.verbose("echelon modulo p=%s (%.2f%% done)"%(
                           p, 100*float(len(str(prod))) / len(str(M))), level=2)

                # We use denoms=False, since we made self integral by calling clear_denom above.
                A = self.matrix_modint(p, denoms=False)
                t = sage.misc.all.verbose("time to reduce matrix mod p:",t, level=2)
                A.echelon()
                t = sage.misc.all.verbose("time to put reduced matrix in echelon form:",t, level=2)
                c = matrix_rational_dense.cmp_pivots(best_pivots, A.pivots())
                if c <= 0:
                    best_pivots = A.pivots()
                    X.append(A)
                    prod = prod * p
                else:
                    # do not save A since it is bad.
                    if sage.misc.all.LEVEL > 1:
                        sage.misc.all.verbose("Excluding this prime (bad pivots).")
                p = sage.rings.arith.next_prime(p)
                t = sage.misc.all.verbose("time for pivot compare", t, level=2)
            # Find set of best matrices.
            Y = []
            # recompute product, since may drop bad matrices
            prod = 1
            for i in range(len(X)):
                if matrix_rational_dense.cmp_pivots(best_pivots, X[i].pivots()) <= 0:
                    Y.append(X[i])
                    prod = prod * X[i].prime()
            try:
                t = sage.misc.all.verbose("start crt and rr", level=2)
                E = lift_matrices_modint(Y)
                #print "E = ", E    # debug
                sage.misc.all.verbose("crt and rr time is",t, level=2)
            except ValueError, msg:
                #print msg # debug
                sage.misc.all.verbose("Redoing with several more primes", level=2)
                for i in range(3):
                    M = M * START_PRIME
                continue

            if not proof:
                sage.misc.all.verbose("Not checking validity of result (since proof=False).", level=2)
                break
            d   = E.denom()
            hdE = long(E.height(d))
            if hdE * self.ncols() * hA < prod:
                break
            for i in range(3):
                M = M * START_PRIME
        #end while
        sage.misc.all.verbose("total time",tm, level=2)
        self.__pivots = best_pivots
        E.__pivots = best_pivots
        self.divide_by(dd)
        return E

    def echelon(self):
        """
        Replace self by its reduction to reduced row echelon form.

        ALGORITHM:
        We use Gauss elimination, which is slightly intelligent, in
        these sense that we clear each column using a row with the
        minimum number of nonzero entries.

        WARNING: There is no reason to use the code below, except
                 for testing this class.  It is *vastly* faster to use
                 the multi-modular method, which is implemented in
                 sparse_matrix.Sparse_matrix_rational
        """
        cdef int i, r, c, min, min_row, start_row, sgn
        cdef mpq_vector tmp
        cdef mpq_t a, a_inverse, b, minus_b
        mpq_init(a)
        mpq_init(a_inverse)
        mpq_init(b)
        mpq_init(minus_b)

        start_row = 0
        self.__pivots = []
        for c from 0 <= c < self.nc:
            _sig_check
            if c % 10 == 0:
                sage.misc.all.verbose('clearing column %s of %s'%(c,self.nc),
                                      caller_name = 'matrix_rational_sparse echelon')
            min = self.nc + 1
            min_row = -1
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            for r from start_row <= r < self.nr:
                if self.rows[r].num_nonzero > 0 and self.rows[r].num_nonzero < min:
                    # Since there is at least one nonzero entry, the first entry
                    # of the positions list is defined.  It is the first position
                    # of a nonzero entry, and it equals c precisely if row r
                    # is a row we could use to clear column c.
                    if self.rows[r].positions[0] == c:
                        min_row = r
                        min = self.rows[r].num_nonzero
                    #endif
                #endif
            #endfor
            if min_row != -1:
                r = min_row
                self.__pivots.append(c)
                # Since we can use row r to clear column c, the
                # entry in position c in row r must be the first
                # nonzero entry.
                mpq_inv(a_inverse, self.rows[r].entries[0])
                scale_mpq_vector(&self.rows[r], a_inverse)
                self.swap_rows(r, start_row)
                for i from 0 <= i < self.nr:
                    if i != start_row:
                        mpq_vector_get_entry(&b, &self.rows[i], c)
                        if mpq_sgn(b) != 0:   # if b is nonzero
                            mpq_neg(minus_b, b)
                            add_mpq_vector_init(&tmp, &self.rows[i],
                                                &self.rows[start_row], minus_b)
                            clear_mpq_vector(&self.rows[i])
                            self.rows[i] = tmp
                start_row = start_row + 1
        mpq_clear(a)
        mpq_clear(a_inverse)
        mpq_clear(b)
        mpq_clear(minus_b)



#####################################



def Matrix_rational_sparse_from_columns(columns):
    """
    Create a sparse Matrix_rational_sparse from a list of sparse Vector_mpq's.
    Each vector must have same degree.
    INPUT:
        columns -- a list of Vector_mpq's.
    OUTPUT:
        A sparse Matrix_rational_sparse whose columns are as given.
    """
    if not isinstance(columns, list):
        raise TypeError, "columns must be a list"

    cdef int j, nc, nr
    nc = len(columns)
    if nc == 0:
        return Matrix_rational_sparse(0,0)
    if not isinstance(columns[0], Vector_mpq):
        raise TypeError, "each column must be of type Vector_mpq"
    nr = columns[0].degree()
    entries = []
    for j from 0 <= j < nc:
        v = columns[j]
        if not isinstance(v, Vector_mpq):
            raise TypeError, "each column must be of type Vector_mpq"
        if v.degree() != nr:
            raise IndexError, "each column must have degree the number of rows of self."
        for i, x in v.list():
            # now the i,j entry of our matrix should be set equal to x."
            entries.append((i,j,x))
    return Matrix_rational_sparse(nr, nc, entries)




def lift_matrices_modint(X):
    """
    INPUT:
        X -- list of sparse Matrix_modint matrices.

    OUTPUT:
        Matrix_rational_sparse -- computed if possible using rational reconstruction
    raises ValueError if there is no valid lift.

    ALGORITHM:
        0. validate input -- type of elements of X and dimensions match up.
        1. compute CRT basis as array of mpz_t's
        2. allocate new matrix
        3. for each row:
           * find union of nonzero positions as Python int set/list
           * allocate relevant memory in matrix over Q that we're constructing
           * for each nonzero position:
                - use CRT basis to lift to mod prod of moduli
                - use rational reconstruction to list to Q
                  (if fail clean up memory and raise exception)
                   -- always multiply element by lcm of denominators so far by attempting
                      rr (and only do rr if there is a denominator)
                - enter value in output matrix
        4. memory cleanup:
            * clear crt basis
    """
    cdef int i, r, c, nr, nc, lenX
    cdef sage.rings.integer.Integer a, b

    #print "X = ", X       # debug

    # 0. validate input
    lenX = len(X)   # save as C int
    # Create C level access to the entries of the Python array X
    for i from 0 <= i < lenX:
        if i == 0:
            nr = X[i].nrows()
            nc = X[i].ncols()
        else:
            if X[i].nrows() != nr or X[i].ncols() != nc:
                raise ValueError, "number of rows and columns of all input matrices must be the same"
    #print "validated input"

    # 1. crt basis as mpz_t'
    #    We use Python since this doesn't have to be fast.
    import sage.rings.integer   # I don't know why this is needed, but it is. (must be a weird pyrex thing)
    P = []
    for A in X:
        a = sage.rings.integer.Integer(A.prime())
        P.append(a)     # no list comprehensions in pyrex

    import sage.rings.integer_ring
    _B = sage.rings.integer_ring.crt_basis(P)
    #print "computed crt basis = ", _B

    cdef mpz_t* B
    cdef mpz_t prod    # product of moduli

    mpz_init(prod)
    mpz_set_si(prod, 1)

    B = <mpz_t*> sage_malloc(sizeof(mpz_t) * len(_B))
    if B == <mpz_t*> 0:
        raise MemoryError, "Error allocating memory"
    for i from 0 <= i < len(_B):
        mpz_init(B[i])
        b = _B[i]
        mpz_set(B[i], b.value)
        mpz_mul_si(prod, prod, X[i].prime())

    #print "converted crt basis to C"


    # 2. create un-allocated matrix over Q
    cdef Matrix_rational_sparse M
    M = Matrix_rational_sparse(nr, nc, init=False)
    #print "created un-allocated matrix"

    cdef int num_nonzero
    cdef mpz_t crt_lift, tmp

    mpz_init(crt_lift)
    mpz_init(tmp)

    cdef mpq_t denom
    mpq_init(denom)
    mpq_set_si(denom, 1, 1)

    # 3. for each row...
    cdef matrix_modn_sparse.Matrix_modn_sparse C

    error = None

    cdef int max_row_allocated

    try:
        for r from 0 <= r < nr:
            #print "considering row r=%s"%r
            #_sig_check
            nonzero_positions = []
            for i from 0 <= i < lenX:
                C = X[i]
                for j from 0 <= j < C.rows[r].num_nonzero:
                    nonzero_positions.append(C.rows[r].positions[j])
            nonzero_positions = list(set(nonzero_positions))
            nonzero_positions.sort()
            num_nonzero = len(nonzero_positions)
            #print "num_nonzero = ", num_nonzero

            # allocate space for that many nonzero entries.
            init_mpq_vector(&M.rows[r], nc, num_nonzero)
            max_row_allocated = r
            for j from 0 <= j < num_nonzero:
                M.rows[r].positions[j] = nonzero_positions[j]

            for j from 0 <= j < num_nonzero:
                # compute each lifted element and insert into M.rows[r]
                mpz_set_si(crt_lift, 0)
                for i from 0 <= i < lenX:
                    C = X[i]       # this could slow everything down a lot... ?
                    # tmp = (ith CRT basis) * (entry of ith modint matrix)
#matrix_rational_sparse.pyx:2058:58: Cannot assign type 'c_vector_modint (*)' to 'c_vector_modint (*)' ?!
                    mpz_mul_si(tmp,     B[i],   get_entry(<c_vector_modint *>&C.rows[r], M.rows[r].positions[j]))
                    mpz_add(crt_lift, crt_lift, tmp)

                mpz_mul(crt_lift, crt_lift, denom)
                mpq_rational_reconstruction(M.rows[r].entries[j],
                                            crt_lift, prod)
                mpq_div(M.rows[r].entries[j], M.rows[r].entries[j], denom)
                mpz_lcm(mpq_numref(denom), mpq_numref(denom), mpq_denref(M.rows[r].entries[j]))


    except ValueError, msg:
        #print "msg = ", msg
        error = ValueError
        # Delete memory in M, since M.is_init isn't true, so M would never
        # get deleted.  We do this here, since otherwise we'd have to
        # finish constructing M just so the __dealloc__ method would work correctly.
        for r from 0 <= r <= max_row_allocated:
            clear_mpq_vector(&M.rows[r])

    # 4. memory cleanup
    #print "cleaning up memory"
    mpz_clear(crt_lift)
    mpz_clear(denom)
    mpz_clear(prod)
    mpz_clear(tmp)

    for i from 0 <= i < len(_B):
        mpz_clear(B[i])
    sage_free(B)

    if not error is None:
        raise error, msg

    # done.
    M.is_init = True
    return M


def make_sparse_rational_matrix(nrows, ncols, entries):
    return Matrix_rational_sparse(nrows, ncols, entries = entries, init=False)
