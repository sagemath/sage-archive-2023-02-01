"""
Matrix over integers modulo p, for p <= 46340
"""

include "../ext/interrupt.pxi"
include "../ext/cdefs.pxi"
include '../ext/stdsage.pxi'

cimport matrix_dense
cimport matrix

from sage.rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract

from sage.structure.element import ModuleElement

# change this to use extern cdef's methods.
from sage.ext.arith cimport arith_int
cdef arith_int ai
ai = arith_int()

##############################################################################
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

cdef class Matrix_modn_dense(matrix_dense.Matrix_dense):
    ########################################################################
    # LEVEL 1 functionality
    # x * __new__
    # x * __dealloc__
    # x * __init__
    # x * set_unsafe
    # x * get_unsafe
    # x * __richcmp__    -- always the same
    ########################################################################
    def __new__(self, parent, entries, copy, coerce):
        matrix_dense.Matrix_dense.__init__(self, parent)

        self.matrix = <uint **> sage_malloc(sizeof(uint*)*self._nrows)
        if self.matrix == NULL:
            raise MemoryError, "Error allocating memory"

        self._entries = <uint *> sage_malloc(sizeof(uint)*self._nrows*self._ncols)
        if self._entries == NULL:
           raise MemoryError, "Error allocating matrix"

        cdef uint k
        k = 0
        for i from 0 <= i < self._nrows:
            self.matrix[i] = self._entries + k
            k = k + self._ncols

    def __dealloc__(self):
        if self.matrix == NULL: # TODO: should never happen now, right
            return
        sage_free(self._entries)
        sage_free(self.matrix)

    def __init__(self, parent, entries, copy, coerce):
        """
        sage: test with negative numbers in input entries list
        """
        cdef uint p
        self.p = self._base_ring.characteristic()
        p = self.p
        if p >= 46340:
            raise OverflowError, "p (=%s) must be < 46340"%p
        self.gather = 2**32/(p*p)

        cdef uint e
        cdef Py_ssize_t i, j, k
        cdef uint *v

        # scalar?
        if not isinstance(entries, list):
            if self._nrows != self._ncols:
                raise TypeError, "scalar matrix must be square"
            e = entries   # coerce to an unsigned int
            for i from 0 <= i < self._nrows:
                v = self.matrix[i]
                for j from 0 <= j < i:
                    v[j] = 0
                v[i] = e
                for j from i+1 <= j < self._ncols:
                    v[j] = 0
            return

        # all entries are given as a long list
        if len(entries) != self._nrows * self._ncols:
            raise IndexError, "The vector of entries has the wrong length."

        k = 0
        cdef uint n

        for i from 0 <= i < self._nrows:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            v = self.matrix[i]
            for j from 0 <= j < self._ncols:
                v[j] = int(entries[k]) % p
                k = k + 1


    def __richcmp__(Matrix_modn_dense self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        self.matrix[i][j] = (<IntegerMod_int> value).ivalue

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef IntegerMod_int n
        n =  IntegerMod_int.__new__(IntegerMod_int)
        IntegerMod_abstract.__init__(n, self._base_ring)
        n.ivalue = self.matrix[i][j]
        return n

    ########################################################################
    # LEVEL 2 functionality
    #   * cdef _pickle
    #   * cdef _unpickle
    #   * cdef _add_c_impl
    #   * cdef _mul_c_impl
    #   * cdef _cmp_c_impl
    #   * __neg__
    #   * __invert__
    #   * __copy__
    #   * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################
    # cdef _pickle(self):
    # cdef _unpickle(self, data, int version):   # use version >= 0
    # cdef ModuleElement _add_c_impl(self, ModuleElement right):
    # cdef _mul_c_impl(self, Matrix right):
    # cdef int _cmp_c_impl(self, Matrix right) except -2:
    # def __neg__(self):
    # def __invert__(self):
    # def __copy__(self):
    # def _multiply_classical(left, matrix.Matrix _right):
    # def _list(self):
    # def _dict(self):


    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_c_impl
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    #        - all row/column operations, but optimized
    #        - echelon form in place
    #        - Hessenberg forms of matrices
    ########################################################################

    # fix all type conversion and Py_ssize_t's below

    def _echelon_in_place_classical(self):
        x = self.fetch('in_echelon_form')
        if not x is None and x: return  # already known to be in echelon form
        self.check_mutability()

        cdef Py_ssize_t start_row, c, r, nr, nc, i
        cdef uint p, a, a_inverse, b
        cdef uint **m

        start_row = 0
        p = self.p
        m = self.matrix
        nr = self._nrows
        nc = self._ncols
        pivots = []
        for c from 0 <= c < nc:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            for r from start_row <= r < nr:
                a = m[r][c]
                if a:
                    pivots.append(c)
                    a_inverse = ai.c_inverse_mod_int(a, p)
                    self.rescale_row_c(r, a_inverse, c)
                    self.swap_rows_c(r, start_row)
                    for i from 0 <= i < nr:
                        if i != start_row:
                            b = m[i][c]
                            if b != 0:
                                self.add_multiple_of_row_c(i, start_row, p-b, c)
                    start_row = start_row + 1
                    break
        self.cache('pivots',pivots)
        self.cache('in_echelon_form',True)

    def rescale_row(self, Py_ssize_t i, s, Py_ssize_t start_col=0):
        self.check_bounds_and_mutability(i,0)
        t = self._coerce_element(s)
        self.rescale_row_c(i, (<IntegerMod_int> t).get_int_value(), start_col)

    cdef rescale_row_c(self, Py_ssize_t row, uint multiple, Py_ssize_t start_col):
        cdef uint r, p
        cdef uint* v
        cdef Py_ssize_t i
        p = self.p
        v = self.matrix[row]
        for i from start_col <= i < self._ncols:
            v[i] = (v[i]*multiple) % p

    def rescale_col(self, Py_ssize_t i, s, Py_ssize_t start_row=0):
        """
        EXAMPLES:
            sage: n=3; b = MatrixSpace(Integers(37),n,n,sparse=False)([1]*n*n)
            sage: b
            [1 1 1]
            [1 1 1]
            [1 1 1]
            sage: b.rescale_col(1,5)
            sage: b
            [1 5 1]
            [1 5 1]
            [1 5 1]

        Bounds are checked.
            sage: b.rescale_col(3,2)
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range

        Recaling need not include the entire row.
            sage: b.rescale_col(0,2,1); b
            [1 5 1]
            [2 5 1]
            [2 5 1]
        """
        self.check_bounds_and_mutability(0,i)
        self.rescale_col_c(i, s, start_row)

    cdef rescale_col_c(self, Py_ssize_t col, uint multiple, Py_ssize_t start_row):
        cdef uint r, p, t
        cdef uint* v
        cdef Py_ssize_t i

        p = self.p
        t = self._coerce_element(multiple)
        x = (<IntegerMod_int> t).get_int_value()
        for i from start_row <= i < self._nrows:
            self.matrix[i][col] = (self.matrix[i][col]*x) % p

    def add_multiple_of_row(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t col_start=0):
        self.check_bounds_and_mutability(i,0)
        self.check_bounds_and_mutability(j,0)
        s = self._coerce_element(s)
        cdef uint x
        x = (<IntegerMod_int> s).get_int_value()
        self.add_multiple_of_row_c(i,j,x,col_start)

    cdef add_multiple_of_row_c(self,  uint row_to, uint row_from, uint multiple,
                               uint start_col):
        cdef uint i, p, nc
        cdef uint *v_from, *v_to
        p = self.p
        v_from = self.matrix[row_from]
        v_to = self.matrix[row_to]
        nc = self._ncols
        for i from start_col <= i < nc:
            v_to[i] = (multiple * v_from[i] +  v_to[i]) % p

    def add_multiple_of_column(self, Py_ssize_t i, Py_ssize_t j, s, Py_ssize_t row_start=0):
        self.check_bounds_and_mutability(0,i)
        self.check_bounds_and_mutability(0,j)
        self.add_multiple_of_column_c(i, j, s, 0)

    cdef add_multiple_of_column_c(self, uint col_to, uint col_from, uint multiple,
                               uint start_row):
        cdef uint i, p, nr
        cdef uint **m
        m = self.matrix
        p = self.p
        nr = self._nrows
        for i from start_row <= i < self._nrows:
            m[i][col_to] = (m[i][col_to] + multiple * m[i][col_from]) %p

    def swap_rows(self, r1, r2):
        self.check_bounds_and_mutability(r1,0)
        self.check_bounds_and_mutability(r2,0)
        self.swap_rows_c(r1, r2)

    cdef swap_rows_c(self, uint row1, uint row2):
        cdef uint* temp
        temp = self.matrix[row1]
        self.matrix[row1] = self.matrix[row2]
        self.matrix[row2] = temp

    def swap_columns(self, c1, c):
        self.check_bounds_and_mutability(0,c1)
        self.check_bounds_and_mutability(0,c2)
        self.swap_columns_c(c1, c)

    cdef swap_columns_c(self, uint col1, uint col2):
        cdef uint i, t, nr
        cdef uint **m
        m = self.matrix
        nr = self._nrows
        for i from 0 <= i < self._nrows:
            t = m[i][col1]
            m[i][col1] = m[i][col2]
            m[i][col2] = t

    def hessenberg_form(self):
        """
        Transforms self in place to its Hessenberg form.
        """
        self.check_mutability()
        x = self.fetch('in_hessenberg_form')
        if not x is None and x: return  # already known to be in Hessenberg form

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
                     self.swap_rows_c(i,m)
                     self.swap_columns_c(i,m)

                 # Now the nonzero entry in position (m,m-1) is t.
                 # Use t to clear the entries in column m-1 below m.
                 for j from m+1 <= j < n:
                     if h[j][m-1]:
                         u = (h[j][m-1] * t_inv) % p
                         self.add_multiple_of_row_c(j, m, p - u, 0)  # h[j] -= u*h[m]
                         # To maintain charpoly, do the corresponding
                         # column operation, which doesn't mess up the
                         # matrix, since it only changes column m, and
                         # we're only worried about column m-1 right
                         # now.  Add u*column_j to column_m.
                         self.add_multiple_of_column_c(m, j, u, 0)
                 # end for
            # end if
        # end for
        _sig_off
        self.cache('in_hessenberg_form',True)

    def _charpoly_hessenberg(self):
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
        self.hessenberg_form()

        cdef Matrix_modn_dense H
        H = self  # just for notational purposes below

        # We represent the intermediate polynomials that come up in
        # the calculations as rows of an (n+1)x(n+1) matrix, since
        # we've implemented basic arithmetic with such a matrix.
        # Please see the generic implementation of charpoly in
        # matrix.py to see more clearly how the following algorithm
        # actually works.  (The implementation is clearer (but slower)
        # if one uses polynomials to represent polynomials instead of
        # using the rows of a matrix.)  Also see Cohen's first GTM,
        # Algorithm 2.2.9.

        cdef Matrix_modn_dense c
        c = self.new_matrix(nrows=n+1,ncols=n+1)    # the 0 matrix
        c.matrix[0][0] = 1
        for m from 1 <= m <= n:
            # Set the m-th row of c to (x - H[m-1,m-1])*c[m-1] = x*c[m-1] - H[m-1,m-1]*c[m-1]
            # We do this by hand by setting the m-th row to c[m-1]
            # shifted to the right by one.  We then add
            # -H[m-1,m-1]*c[m-1] to the resulting m-th row.
            for i from 1 <= i <= n:
                c.matrix[m][i] = c.matrix[m-1][i-1]
            # the p-.. below is to keep scalar normalized between 0 and p.
            c.add_multiple_of_row_c(m, m-1, p - H.matrix[m-1][m-1], 0)
            t = 1
            for i from 1 <= i < m:
                t = (t*H.matrix[m-i][m-i-1]) % p
                # Set the m-th row of c to c[m] - t*H[m-i-1,m-1]*c[m-i-1]
                c.add_multiple_of_row_c(m, m-i-1, p - (t*H.matrix[m-i-1][m-1])%p, 0)

        # The answer is now the n-th row of c.
        v = []
        for i from 0 <= i <= n:
            v.append(int(c.matrix[n][i]))
        return v
