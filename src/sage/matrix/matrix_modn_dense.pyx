"""
Generic matrices over the integers modulo n.
"""

cimport matrix_dense
import matrix_dense

include "../ext/interrupt.pxi"
include "../ext/cdefs.pxi"


from sage.misc.misc import verbose, get_verbose

cimport sage.ext.arith
import  sage.ext.arith
cdef sage.ext.arith.arith_int ai
ai = sage.ext.arith.arith_int()


##################################################################
# Matrix over integers modulo p, for p <= 46340
##################################################################

LEAVE_UNINITIALIZED = "LEAVE UNINITIALIZED"

MAX_MODULUS = 46340

cdef class Matrix_modn_dense(matrix_dense.Matrix_dense):
    def __init__(self, parent, object entries=None, coerce_entries=True, copy=True, clear=True):
        matrix_dense.Matrix_dense.__init__(self, parent)

        cdef int i, p, nrows, ncols
        nrows = parent.nrows()
        ncols = parent.ncols()
        p = parent.base_ring().characteristic() # Should I do a sanity check here?
        self.matrix = <uint **> PyMem_Malloc(sizeof(uint*)*nrows)
        if self.matrix == <uint**> 0:
            raise MemoryError, "Error allocating memory"

        self._entries = <uint *> PyMem_Malloc(sizeof(uint)*nrows*ncols)
        if self._entries == <uint*> 0:
           raise MemoryError, "Error allocating matrix"

        cdef uint n, j, k

        k = 0
        for i from 0 <= i < nrows:
            self.matrix[i] = self._entries + k
            k = k + ncols

        cdef uint *v
        if p >= 46340:
            raise OverflowError, "p (=%s) must be < 46340"%p

        self.p = p
        self._nrows = nrows
        self._ncols = ncols
        self.gather = 2**32/(p*p)
        self.__pivots = None
        if entries == 0:
            entries = None
            clear = True
        if (entries is None or entries == 1) and clear:
            for i from 0 <= i < nrows:
                v = self.matrix[i]
                for j from 0 <= j < ncols:
                    v[j] = 0
            return
        if entries == 1:
            if nrows != ncols:
                raise TypeError, "scalar matrix must be square"
            for i from 0 <= i < nrows:
                self.matrix[i][i] = 1
            return
        if not clear or entries is LEAVE_UNINITIALIZED:
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
                if v[j] < 0:  # WARNING testing uint for negative value!
                    v[j] = v[j] + p
                k = k + 1

    def new_matrix(Matrix_modn_dense self, uint nrows=-1, uint ncols=-1, clear=True):
      if nrows == -1:
          nrows = self._nrows
      if ncols == -1:
          ncols = self._ncols
      return Matrix_modn_dense(self.matrix_space(nrows, ncols), clear=clear)

    def copy(Matrix_modn_dense self):
        cdef Matrix_modn_dense other
        other = self.new_matrix(clear=False)
        cdef int i, len
        len = self._ncols * sizeof(uint)
        for i from 0 <= i < self._nrows:
            memcpy(other.matrix[i], self.matrix[i], len)
        return other

    cdef set_matrix(Matrix_modn_dense self, uint **m):
        if self.matrix != <uint **>0:
            raise RuntimeError, "Only set matrix of uninitialized matrix."
        self.matrix = m

    def __getitem__(self, t):
        if not isinstance(t, tuple) or len(t) != 2:
            raise IndexError, "Index of matrix item must be a row and a column."
        cdef int i, j
        i, j = t
        if i<0 or i >= self._nrows or j<0 or j >= self._ncols:
            raise IndexError, "Array index out of bounds."
        return self.matrix[i][j]

    def __setitem__(self, t, x):
        if not isinstance(t, tuple) or len(t) != 2:
            raise IndexError, "Index for setting matrix item must be a row and a column."
        cdef int i, j
        i, j = t
        if i<0 or i >= self._nrows or j<0 or j >= self._ncols:
            raise IndexError, "Array index out of bounds."
        self.matrix[i][j] = int(x) % self.p

    def get_entry(Matrix_modn_dense self, int i, int j):
        if i < 0 or j < 0 or i >= self._nrows or j >= self._ncols:
            raise IndexError, "Array index (%s,%s) out of bounds."%(i,j)
        return self.matrix[i][j]

    def __call__(self, i, j):
        return self.matrix[i][j]


    def _cmp(self, Matrix_modn_dense other):
        cdef uint i, j
        if self.p != other.p or self._nrows != other._nrows \
            or self._ncols != other._ncols:
            return -1
        _sig_on
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                if self.matrix[i][j] != other.matrix[i][j]:
                    _sig_off
                    return -1
        _sig_off
        return 0

    def __cmp__(self, other):
        if not isinstance(other, Matrix_modn_dense):
            return -1
        return self._cmp(other)

    cdef uint **get_matrix(Matrix_modn_dense self):
        return self.matrix

    def  __dealloc__(self):
        if self.matrix == <uint **> 0: # TODO: should never happen now, right
            return
        PyMem_Free(self._entries)
        PyMem_Free(self.matrix)

##     def _lift_to_Q(self):
##         cdef Matrix_rational M
##         M = Matrix_rational(self._nrows, self._ncols, LEAVE_UNINITIALIZED)

##         cdef mpq_t **m
##         m = <mpq_t **> PyMem_Malloc(sizeof(mpq_t*)*self._nrows)
##         if m == <mpq_t**> 0:
##             raise MemoryError, "Error allocating matrix"
##         _sig_on
##         cdef int i, j
##         for i from 0 <= i < self._nrows:
##             m[i] = <mpq_t *> PyMem_Malloc(sizeof(mpq_t)*self._ncols)
##             if m[i] == <mpq_t*> 0:
##                 _sig_off
##                 raise MemoryError, "Error allocating matrix"
##             for j from 0 <= j < self._ncols:
##                 mpq_init(m[i][j])
##                 mpq_set_ui(m[i][j], self.matrix[i][j], 1)
##         _sig_off
##         M.set_matrix(m)
##         return M

    def __mul__(Matrix_modn_dense self, Matrix_modn_dense other):
        if self.gather >= 2:
            return self._multiply_with_delayed_mod(other)
        if self._ncols != other._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self.p != other.p:
            raise ArithmeticError, "The base matrices must have the same modulus."

        cdef Matrix_modn_dense M
        M = Matrix_modn_dense(self.parent, LEAVE_UNINITIALIZED)
        cdef uint **m
        m = M.matrix

        cdef uint i, j, k, nr, nc, s, p
        cdef uint *v
        nr = self._nrows
        nc = other._ncols
        p = self.p

        _sig_on
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                s = 0
                v = self.matrix[i]
                for k from 0 <= k < self._ncols:
                    s = (s + (v[k] * other.matrix[k][j]))%p
                m[i][j] = s
        _sig_off
        return M

    def __add__(Matrix_modn_dense self, Matrix_modn_dense other):
        if self._ncols != other._ncols:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self._nrows != other._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self.p != other.p:
            raise ArithmeticError, "The base matrices must have the same modulus."

        cdef Matrix_modn_dense M
        M = Matrix_modn_dense(self.parent().matrix_space(self._nrows, self._ncols), LEAVE_UNINITIALIZED)
        cdef uint **m
        m = M.matrix

        cdef uint i, j, nr, nc, s, p, a
        cdef uint *v, *w
        nr = self._nrows
        nc = other._ncols
        p = self.p

        _sig_on
        for i from 0 <= i < nr:
            v = self.matrix[i]
            w = other.matrix[i]
            for j from 0 <= j < nc:
                a = v[j] + w[j]
                if a >= p:
                    m[i][j] = a - p
                else:
                    m[i][j] = a
        _sig_off
        return M

    def __sub__(Matrix_modn_dense self, Matrix_modn_dense other):
        if self._ncols != other._ncols:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self._nrows != other._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self.p != other.p:
            raise ArithmeticError, "The base matrices must have the same modulus."

        cdef Matrix_modn_dense M
        M = Matrix_modn_dense(self.parent(), LEAVE_UNINITIALIZED)
        cdef uint **m
        m = M.matrix

        cdef uint i, j, k, nr, nc, s, p, a
        cdef uint *v, *w
        nr = self._nrows
        nc = self._ncols
        p = self.p

        _sig_on
        for i from 0 <= i < nr:
            v = self.matrix[i]
            w = other.matrix[i]
            for j from 0 <= j < nc:
                a = v[j] + (p - <int> w[j])
                if a >= p:
                    m[i][j] = a - p
                else:
                    m[i][j] = a
        _sig_off
        return M


##     def _invert_submatrices(Matrix_modn_dense self,
##                             int self_r, int self_c, int n):
##         """
##         INPUT:
##              self, other -- two matrices

##         OUTPUT:
##              inverse of submatrices
##         """
##         raise NotImplementedError


    def _mul_submatrices(Matrix_modn_dense self,
                         Matrix_modn_dense other,
                         int self_r, int self_c, int self_nrows, int self_ncols,
                         int other_r, int other_c, int other_nrows, int other_ncols):
        """
        INPUT:
             self, other -- two matrices
             self_, other_ -- positions

        OUTPUT:
             product of the submatrices, as a Matrix_modn_dense
        """
        if self.gather >= 2:
            return self._mul_submatrices_with_delayed_mod(other,
                                                          self_r, self_c, self_nrows, self_ncols,
                                                          other_r, other_c, other_nrows, other_ncols)

        if self_ncols != other_nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self.p != other.p:
            raise ArithmeticError, "The base matrices must have the same modulus."

        cdef Matrix_modn_dense M
        M = Matrix_modn_dense(self.parent().matrix_space(self_nrows, other_ncols), LEAVE_UNINITIALIZED)
        cdef uint **m
        m = M.matrix

        cdef uint i, j, k, nr, nc, s, p
        cdef uint *v
        nr = self_nrows
        nc = other_ncols
        p = self.p

        _sig_on
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                s = 0
                v = self.matrix[i + self_r]
                for k from 0 <= k < self_ncols:
                    s = (s + (v[k + self_c] * other.matrix[k + other_r][j + other_c]))%p
                m[i][j] = s
        _sig_off
        return M

    def _mul_submatrices_with_delayed_mod(Matrix_modn_dense self,
                                          Matrix_modn_dense other,
                                          int self_r, int self_c, int self_nrows, int self_ncols,
                                          int other_r, int other_c, int other_nrows, int other_ncols):
        if self_ncols != other_nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."

        if self.p != other.p:
            raise ArithmeticError, "The base matrices must have the same modulus."

        cdef Matrix_modn_dense M
        M = Matrix_modn_dense(self.parent().matrix_space(self_nrows, self_ncols), LEAVE_UNINITIALIZED)
        cdef uint **m
        m = M.matrix

        cdef uint i, j, k, nr, nc, snc, s, p, gather, w, groups, a, b
        cdef uint *v
        nr = self_nrows
        nc = other_ncols
        snc = self_ncols
        gather = self.gather
        p = self.p
        groups = (self_ncols / gather) + 1
        _sig_on
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                s = 0
                v = self.matrix[i+self_r]
                for w from 0 <= w < groups:
                    a = w*gather
                    b = (w+1)*gather
                    if b > snc or b < a: b = snc
                    for k from a <= k < b:
                        s = s + v[k+self_c]*other.matrix[k+other_r][j+other_c]
                    s = s % p
                m[i][j] = s
        _sig_off
        return M


    def block2_sum(self, Matrix_modn_dense B, Matrix_modn_dense C, Matrix_modn_dense D):
        cdef Matrix_modn_dense M
        cdef uint nr, nc, i, j, s, p, a

        nr = self._nrows + C._nrows
        nc = self._ncols + B._ncols

        M = Matrix_modn_dense(self.parent().matrix_space(nr,nc), LEAVE_UNINITIALIZED)

        cdef uint **m
        m = M.matrix
        p = self.p

        _sig_on
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                m[i][j] = self.matrix[i][j]
                m[i][j+self._ncols] = B.matrix[i][j]
                m[i+self._nrows][j] = C.matrix[i][j]
                m[i+self._nrows][j+self._ncols] = D.matrix[i][j]

        _sig_off
        return M

    def _add_submatrices(Matrix_modn_dense self, Matrix_modn_dense other,
                         int self_r, int self_c, int self_nrows, int self_ncols,
                         int other_r, int other_c, int other_nrows, int other_ncols):
        """
        INPUT:
             self, other -- two matrices
             self_, other_ -- positions

        OUTPUT:
             sum of the submatrices, as a Matrix_modn_dense

        EXAMPLES:
            sage: from sage.matrix.matrix_modn_dense import Matrix_modn_dense
            sage: n = 4
            sage: A = Matrix_modn_dense(389, n,n, range(n^2))
            sage: B = Matrix_modn_dense(389, n,n, list(reversed(range(n^2))))
            sage: A._add_submatrices(B, 0, 0, 2, 2,  0, 0, 2, 2)
            [
            15, 15,
            15, 15
            ]
            sage: A._add_submatrices(B, 0, 0, 2, 2,  2, 2, 2, 2)
            [
            5, 5,
            5, 5
            ]
            sage: A._add_submatrices(B, 2,0, 2, 2,  0, 2, 2, 2)
            [
            21, 21,
            21, 21
            ]
            sage: A._add_submatrices(B, 0,0, 3, 3,  1, 1, 3, 3)
            [
            10, 10, 10,
            10, 10, 10,
            10, 10, 10
            ]
        """
        cdef Matrix_modn_dense M
        cdef uint i, j, s, p, a
        cdef uint *v, *w

        M = Matrix_modn_dense(self.parent().matrix_space(self_nrows, self_ncols), LEAVE_UNINITIALIZED)
        cdef uint **m
        m = M.matrix
        p = self.p

        _sig_on
        for i from 0 <= i < self_nrows:
            v = self.matrix[i+self_r]
            w = other.matrix[i+other_r]
            for j from 0 <= j < self_ncols:
                a = v[self_c+j] + w[other_c+j]
                if a >= p:
                    m[i][j] = a - p
                else:
                    m[i][j] = a
        _sig_off
        return M

    def _sub_submatrices(Matrix_modn_dense self, Matrix_modn_dense other,
                         int self_r, int self_c, int self_nrows, int self_ncols,
                         int other_r, int other_c, int other_nrows, int other_ncols):
        """
        INPUT:
             self, other -- two matrices
             self_, other_ -- positions

        OUTPUT:
             sum of the submatrices, as a Matrix_modn_dense

        EXAMPLES:
            sage: from sage.matrix.matrix_modn_dense import Matrix_modn_dense
            sage: n = 4
            sage: A = Matrix_modn_dense(389, n,n, range(n^2))
            sage: B = Matrix_modn_dense(389, n,n, list(reversed(range(n^2))))
            sage: A._sub_submatrices(B, 0, 0, 2, 2,  0, 0, 2, 2)
            [
            374, 376,
            382, 384
            ]
        """
        cdef Matrix_modn_dense M
        cdef uint i, j, s, p, a
        cdef uint *v, *w

        M = Matrix_modn_dense(self.parent(), LEAVE_UNINITIALIZED)
        cdef uint **m
        m = M.matrix
        p = self.p

        _sig_on
        for i from 0 <= i < self_nrows:
            v = self.matrix[i+self_r]
            w = other.matrix[i+other_r]
            for j from 0 <= j < self_ncols:
                a = v[self_c+j] + (p - <int> w[other_c+j])
                if a >= p:
                    m[i][j] = a - p
                else:
                    m[i][j] = a
        _sig_off
        return M


    def _multiply_with_delayed_mod(Matrix_modn_dense self, Matrix_modn_dense other):
        if self._ncols != other._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self.p != other.p:
            raise ArithmeticError, "The base matrices must have the same modulus."

        cdef Matrix_modn_dense M
        M = Matrix_modn_dense(self.parent().matrix_space(self._nrows, other._ncols), LEAVE_UNINITIALIZED)
        cdef uint **m
        m = M.matrix

        cdef uint i, j, k, nr, nc, snc, s, p, gather, w, groups, a, b
        cdef uint *v
        nr = self._nrows
        nc = other._ncols
        snc = self._ncols
        gather = self.gather
        p = self.p
        groups = (self._ncols / gather) + 1
        _sig_on
        for i from 0 <= i < nr:
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
                ##for k from 0 <= k < self._ncols:
                ##    s = s  +  v[k] * other.matrix[k][j]
                m[i][j] = s
        _sig_off
        return M


    def nrows(self):
        return self._nrows

    def ncols(self):
        return self._ncols

    def prime(self):
        return self.p

    def number_nonzero(self):
        cdef uint i, j, n
        cdef uint *v
        n = 0
        _sig_on
        for i from 0 <= i < self._nrows:
            v = self.matrix[i]
            for j from 0 <= j < self._ncols:
                if v[j] != 0:
                    n = n + 1
        _sig_off
        return n

    def list(self):
        cdef uint i, j
        cdef uint *r
        v = []
        _sig_on
        for i from 0 <= i < self._nrows:
            r = self.matrix[i]
            for j from 0 <= j < self._ncols:
                v.append(r[j])
        _sig_off
        return v

    def echelon(self):
        cdef uint p, start_row, c, r, nr, nc, a, a_inverse, b, i
        cdef uint **m

        start_row = 0
        p = self.p
        m = self.matrix
        nr = self._nrows
        nc = self._ncols
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

    def _set_pivots(self, pivots):
        self.__pivots = pivots

    def hessenberg_form(self):
        """
        Transforms self in place to its Hessenberg form.
        """
        if self._nrows != self._ncols:
            raise ArithmeticError, "Matrix must be square to compute Hessenberg form."

        cdef uint n
        n = self._nrows

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
        if self._nrows != self._ncols:
            raise ArithmeticError, "charpoly not defined for non-square matrix."

        cdef uint i, m, n, p, t
        n = self._nrows
        p = self.p

        # Replace self by its Hessenberg form, and set H to this form
        # for notation below.
        #time = verbose('start hessenberg')
        self.hessenberg_form()
        #verbose('done with hessenberg', time)
        cdef Matrix_modn_dense H
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

        cdef Matrix_modn_dense c
        c = Matrix_modn_dense(self.parent().matrix_space(n+1,n+1))  # the 0 matrix
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
        r = row*self._ncols
        p = self.p
        v = self.matrix[row]
        for i from start_col <= i < self._ncols:
            v[i] = (v[i]*multiple) % p

    cdef add_multiple_of_row(self, uint row_from, uint multiple,
                            uint row_to, uint start_col):
        cdef uint i, p, nc
        cdef uint *v_from, *v_to
        p = self.p
        v_from = self.matrix[row_from]
        v_to = self.matrix[row_to]
        nc = self._ncols
        for i from start_col <= i < nc:
            v_to[i] = (multiple * v_from[i] +  v_to[i]) % p

    cdef add_multiple_of_column(self, uint col_from, uint multiple,
                               uint col_to, uint start_row):
        cdef uint i, p, nr
        cdef uint **m
        m = self.matrix
        p = self.p
        nr = self._nrows
        for i from start_row <= i < self._nrows:
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
        nr = self._nrows
        for i from 0 <= i < self._nrows:
            t = m[i][col1]
            m[i][col1] = m[i][col2]
            m[i][col2] = t

    def matrix_from_cols(Matrix_modn_dense self, cols):
        """
        Return the submatrix formed from the given columns.

        INPUT:
            cols -- list of int's

        OUTPUT:
            matrix created from the columns with given indexes
        """
        cdef int i, j, k, nc, nr
        cdef Matrix_modn_dense M

        if not isinstance(cols, list):
            raise TypeError, "cols (=%s) must be a list"%cols
        nc = len(cols)
        if nc == 0:
            return Matrix_rational_dense(self._nrows, 0)
        nr = self._nrows
        if min(cols) < 0 or max(cols) >= self._ncols:
            raise IndexError, "invalid cols indexes; cols don't exist"

        M = self.new_matrix(nrows=nr, ncols=nc, clear=False)
        cdef uint **m
        m = M.matrix

        for j from 0 <= j < nc:
            k = int(cols[j])
            for i from 0 <= i < nr:
                m[i][j] = self.matrix[i][k]

        return M


    def matrix_window(self, int row=0, int col=0, int nrows=-1, int ncols=-1):
        if nrows == -1:
            nrows = self._nrows - row
            ncols = self._ncols - col
        return MatrixWindow(self, row, col, nrows, ncols)


#cdef uint **Matrix_modn_dense_matrix(Matrix_modn_dense A):
#    return A.matrix



cdef class MatrixWindow:

    def __init__(MatrixWindow self, Matrix_modn_dense matrix, int row, int col, int nrows, int ncols):
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
        a = self._matrix.new_matrix(self._nrows, self._ncols, clear=False)
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
        cdef uint* s_row
        cdef uint* A_row
        for i from 0 <= i < self._nrows:
            memcpy(self._matrix.matrix[self._row + i] + self._col, A._matrix.matrix[A._row + i] + A._col, self._ncols * sizeof(uint))

    def set_to_zero(MatrixWindow self):
        cdef int i, j
        cdef uint* s_row
        for i from 0 <= i < self._nrows:
            memset(self._matrix.matrix[self._row + i] + self._col, 0, self._ncols * sizeof(uint))

    def add(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef uint p
        cdef uint* s_row
        cdef uint* A_row
        p = self._matrix.p
        for i from 0 <= i < self._nrows:
            s_row = self._matrix.matrix[self._row + i] + self._col
            A_row = A._matrix.matrix[A._row + i] + A._col
            for j from 0 <= j < self._ncols:
                s_row[j] = s_row[j] + A_row[j]
                if s_row[j] >= p:
                    s_row[j] = s_row[j] - p

    def subtract(MatrixWindow self, MatrixWindow A):
        cdef int i, j
        cdef uint p
        cdef uint* s_row
        cdef uint* A_row
        p = self._matrix.p
        for i from 0 <= i < self._nrows:
            s_row = self._matrix.matrix[self._row + i] + self._col
            A_row = A._matrix.matrix[A._row + i] + A._col
            for j from 0 <= j < self._ncols:
                s_row[j] = s_row[j] + (p - A_row[j])
                if s_row[j] >= p:
                    s_row[j] = s_row[j] - p

    def set_to_sum(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j
        cdef uint p
        cdef uint* s_row
        cdef uint* A_row
        cdef uint* B_row
        p = self._matrix.p
        for i from 0 <= i < self._nrows:
            s_row = self._matrix.matrix[self._row + i] + self._col
            A_row = A._matrix.matrix[A._row + i] + A._col
            B_row = B._matrix.matrix[B._row + i] + B._col
            for j from 0 <= j < self._ncols:
                s_row[j] = A_row[j] + B_row[j]
                if s_row[j] >= p:
                    s_row[j] = s_row[j] - p

    def set_to_diff(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j
        cdef uint p
        cdef uint* s_row
        cdef uint* A_row
        cdef uint* B_row
        p = self._matrix.p
        for i from 0 <= i < self._nrows:
            s_row = self._matrix.matrix[self._row + i] + self._col
            A_row = A._matrix.matrix[A._row + i] + A._col
            B_row = B._matrix.matrix[B._row + i] + B._col
            for j from 0 <= j < self._ncols:
                s_row[j] = A_row[j] + (p - B_row[j])
                if s_row[j] >= p:
                    s_row[j] = s_row[j] - p

    def set_to_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j, k
        cdef uint p
        cdef uint* s_row
        cdef uint* A_row
        cdef uint sum, limit
        limit = (1 << (sizeof(uint)-1)) - p*p
        p = self._matrix.p
        for i from 0 <= i < A._nrows:
            A_row = A._matrix.matrix[A._row + i] + A._col
            s_row = self._matrix.matrix[self._row + i] + self._col
            for j from 0 <= j < B._ncols:
                sum = A_row[0] * B._matrix.matrix[B._row][B._col+j]
                for k from 1 <= k < A._ncols:
                    sum = sum + A_row[k] * B._matrix.matrix[B._row+k][B._col+j]
                    if sum > limit:
                        sum = sum % p
                s_row[j] = sum % p

    def add_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j, k
        cdef uint p
        cdef uint* s_row
        cdef uint* A_row
        cdef uint sum, limit
        limit = (1 << (sizeof(uint)-1)) - p*p
        p = self._matrix.p
        for i from 0 <= i < A._nrows:
            A_row = A._matrix.matrix[A._row + i] + A._col
            s_row = self._matrix.matrix[self._row + i] + self._col
            for j from 0 <= j < B._ncols:
                sum = s_row[j]
                for k from 0 <= k < A._ncols:
                    sum = sum + A_row[k] * B._matrix.matrix[B._row+k][B._col+j]
                    if sum > limit:
                        sum = sum % p
                s_row[j] = sum % p

    def subtract_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B):
        cdef int i, j, k
        cdef uint p
        cdef uint* s_row
        cdef uint* A_row
        cdef uint sum, limit
        limit = (1 << (sizeof(uint)-1)) - p*p
        p = self._matrix.p
        for i from 0 <= i < A._nrows:
            A_row = A._matrix.matrix[A._row + i] + A._col
            s_row = self._matrix.matrix[self._row + i] + self._col
            for j from 0 <= j < B._ncols:
                sum = A_row[0] * B._matrix.matrix[B._row][B._col+j]
                for k from 1 <= k < A._ncols:
                    sum = sum + A_row[k] * B._matrix.matrix[B._row+k][B._col+j]
                    if sum > limit:
                        sum = sum % p
                s_row[j] = s_row[j] + (p - sum % p)
                if s_row[j] >= p:
                    s_row[j] = s_row[j] - p


    def swap_rows(MatrixWindow self, int a, int b):
        self._matrix.swap_rows(self._row + a, self._row + b)


    def echelon_in_place(MatrixWindow self):
        """
        calculate the echelon form of this matrix, returning the list of pivot columns
        """
        echelon = self.to_matrix()
        echelon.echelon() # TODO: read only, only need to copy pointers
        self.set_to(echelon.matrix_window())
        return echelon.pivots()

    def element_is_zero(MatrixWindow self, int i, int j):
        return self._matrix.matrix[i+self._row][j+self._col] == 0


    def new_empty_window(MatrixWindow self, int nrows, int ncols, zero=True):
        return self._matrix.new_matrix(nrows, ncols, clear=zero).matrix_window()
