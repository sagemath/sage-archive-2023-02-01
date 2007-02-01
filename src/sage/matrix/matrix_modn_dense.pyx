r"""
Dense matrices over $\Z/n\Z$ for $n$ small.

This is a compiled implementation of dense matrices over $\Z/n\Z$
for $n$ small.

EXAMPLES:
    sage: a = matrix(Integers(37),3,range(9),sparse=False); a
    [0 1 2]
    [3 4 5]
    [6 7 8]
    sage: a.rank()
    2
    sage: type(a)
    <type 'sage.matrix.matrix_modn_dense.Matrix_modn_dense'>
    sage: a[0,0] = 5
    sage: a.rank()
    3
    sage: parent(a)
    Full MatrixSpace of 3 by 3 dense matrices over Ring of integers modulo 37

    sage: a^2
    [ 3 23 31]
    [20 17 29]
    [25 16  0]
    sage: a+a
    [10  2  4]
    [ 6  8 10]
    [12 14 16]

    sage: b = a.new_matrix(2,3,range(6)); b
    [0 1 2]
    [3 4 5]
    sage: a*b
    Traceback (most recent call last):
    ...
    TypeError: incompatible dimensions
    sage: b*a
    [15 18 21]
    [20 17 29]

    sage: a == loads(dumps(a))
    True
    sage: b == loads(dumps(b))
    True

    sage: a.echelonize(); a
    [1 0 0]
    [0 1 0]
    [0 0 1]
    sage: b.echelonize(); b
    [ 1  0 36]
    [ 0  1  2]

We create a matrix group and coerce it to GAP:
    sage: M = MatrixSpace(GF(3),3,3)
    sage: G = MatrixGroup([M([[0,1,0],[0,0,1],[1,0,0]]), M([[0,1,0],[1,0,0],[0,0,1]])])
    sage: G
    Matrix group over Finite Field of size 3 with 2 generators:
     [[[0, 1, 0], [0, 0, 1], [1, 0, 0]], [[0, 1, 0], [1, 0, 0], [0, 0, 1]]]
    sage: gap(G)
    Group(
    [ [ [ 0*Z(3), Z(3)^0, 0*Z(3) ], [ 0*Z(3), 0*Z(3), Z(3)^0 ], [ Z(3)^0, 0*Z(3),
               0*Z(3) ] ],
      [ [ 0*Z(3), Z(3)^0, 0*Z(3) ], [ Z(3)^0, 0*Z(3), 0*Z(3) ],
          [ 0*Z(3), 0*Z(3), Z(3)^0 ] ] ])
"""

include "../ext/interrupt.pxi"
include "../ext/cdefs.pxi"
include '../ext/stdsage.pxi'

import matrix_window_modn_dense

cimport matrix_dense
cimport matrix
cimport matrix0

from sage.rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract
cdef extern from "stdint.h":
    ctypedef int int_fast32_t

from sage.structure.element import ModuleElement

from sage.misc.misc import verbose, get_verbose

################
# TODO: change this to use extern cdef's methods.
from sage.ext.arith cimport arith_int
cdef arith_int ai
ai = arith_int()
################

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

        cdef mod_int p
        p = self._base_ring.characteristic()
        self.p = p
        if p >= MOD_INT_MAX:
            raise OverflowError, "p (=%s) must be < %s"%(p, MOD_INT_MAX)
        self.gather = MOD_INT_OVERFLOW / (p*p) - 1

        self._entries = <mod_int *> sage_malloc(sizeof(mod_int)*self._nrows*self._ncols)
        if self._entries == NULL:
           raise MemoryError, "Error allocating matrix"

        self.matrix = <mod_int **> sage_malloc(sizeof(mod_int*)*self._nrows)
        if self.matrix == NULL:
            sage_free(self._entries)
            raise MemoryError, "Error allocating memory"

        cdef mod_int k
        cdef Py_ssize_t i
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

        cdef mod_int e
        cdef Py_ssize_t i, j, k
        cdef mod_int *v

        # scalar?
        if not isinstance(entries, list):
            _sig_on
            for i from 0 <= i < self._nrows*self._ncols:
                self._entries[i] = 0
            e = entries   # coerce to an unsigned int
            if e != 0:
                for i from 0 <= i < self._nrows:
                    self.matrix[i][i] = e
            _sig_off
            return

        # all entries are given as a long list
        if len(entries) != self._nrows * self._ncols:
            raise IndexError, "The vector of entries has the wrong length."

        k = 0
        cdef mod_int n
        R = self.base_ring()

        cdef PyObject** w
        w = FAST_SEQ_UNSAFE(entries)
        if coerce:
            for i from 0 <= i < self._nrows:
                if PyErr_CheckSignals(): raise KeyboardInterrupt
                v = self.matrix[i]
                for j from 0 <= j < self._ncols:
                    v[j] = R( <object> w[k])
                    k = k + 1
        else:
            for i from 0 <= i < self._nrows:
                if PyErr_CheckSignals(): raise KeyboardInterrupt
                v = self.matrix[i]
                for j from 0 <= j < self._ncols:
                    v[j] = int( <object> w[k])
                    k = k + 1


    def __richcmp__(Matrix_modn_dense self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)
    def __hash__(self):
        return self._hash()

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
    # def _pickle(self):
    # def _unpickle(self, data, int version):   # use version >= 0
    # cdef ModuleElement _add_c_impl(self, ModuleElement right):
    # cdef _mul_c_impl(self, Matrix right):
    # cdef int _cmp_c_impl(self, Matrix right) except -2:
    # def __neg__(self):
    # def __invert__(self):
    # def __copy__(self):
    # def _multiply_classical(left, matrix.Matrix _right):
    # def _list(self):
    # def _dict(self):
    def __copy__(self):
        cdef Matrix_modn_dense A
        A = Matrix_modn_dense.__new__(Matrix_modn_dense, self._parent,
                                      0, 0, 0)
        memcpy(A._entries, self._entries, sizeof(mod_int)*self._nrows*self._ncols)
        A.p = self.p
        A.gather = self.gather
        return A



    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_c_impl
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    #        - all row/column operations, but optimized
    # x      - echelon form in place
    #        - Hessenberg forms of matrices
    ########################################################################

    # TODO TODO: fix all type conversion and Py_ssize_t's below

    def _echelon_in_place_classical(self):
        x = self.fetch('in_echelon_form')
        if not x is None: return  # already known to be in echelon form

        self.check_mutability()

        cdef Py_ssize_t start_row, c, r, nr, nc, i
        cdef mod_int p, a, a_inverse, b
        cdef mod_int **m

        start_row = 0
        p = self.p
        m = self.matrix
        nr = self._nrows
        nc = self._ncols
        pivots = []
        fifth = self._ncols / 10 + 1
        do_verb = (get_verbose() >= 2)
        for c from 0 <= c < nc:
            if do_verb and (c % fifth == 0 and c>0):
                tm = verbose('on column %s of %s'%(c, self._ncols),
                             level = 2,
                             caller_name = 'matrix_modn_dense echelon')
            #end if
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            for r from start_row <= r < nr:
                a = m[r][c]
                if a:
                    pivots.append(c)
                    a_inverse = ai.c_inverse_mod_int(a, p)
                    self._rescale_row_c(r, a_inverse, c)
                    self.swap_rows_c(r, start_row)
                    for i from 0 <= i < nr:
                        if i != start_row:
                            b = m[i][c]
                            if b != 0:
                                self._add_multiple_of_row_c(i, start_row, p-b, c)
                    start_row = start_row + 1
                    break
        self.cache('pivots',pivots)
        self.cache('in_echelon_form',True)

    cdef rescale_row_c(self, Py_ssize_t row, multiple, Py_ssize_t start_col):
        self._rescale_row_c(row, multiple, start_col)

    cdef _rescale_row_c(self, Py_ssize_t row, mod_int multiple, Py_ssize_t start_col):
        cdef mod_int r, p
        cdef mod_int* v
        cdef Py_ssize_t i
        p = self.p
        v = self.matrix[row]
        for i from start_col <= i < self._ncols:
            v[i] = (v[i]*multiple) % p

    cdef rescale_col_c(self, Py_ssize_t col, multiple, Py_ssize_t start_row):
        self._rescale_col_c(col, multiple, start_row)

    cdef _rescale_col_c(self, Py_ssize_t col, mod_int multiple, Py_ssize_t start_row):
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

        Recaling need not include the entire row.
            sage: b.rescale_col(0,2,1); b
            [1 5 1]
            [2 5 1]
            [2 5 1]

        Bounds are checked.
            sage: b.rescale_col(3,2)
            Traceback (most recent call last):
            ...
            IndexError: matrix column index out of range

        Rescaling by a negative number:
            sage: b.rescale_col(2,-3); b
            [ 1  5 34]
            [ 2  5 34]
            [ 2  5 34]
        """
        cdef mod_int r, p
        cdef mod_int* v
        cdef Py_ssize_t i
        p = self.p
        for i from start_row <= i < self._nrows:
            self.matrix[i][col] = (self.matrix[i][col]*multiple) % p

    cdef add_multiple_of_row_c(self,  Py_ssize_t row_to, Py_ssize_t row_from, multiple,
                               Py_ssize_t start_col):
        self._add_multiple_of_row_c(row_to, row_from, multiple, start_col)

    cdef _add_multiple_of_row_c(self,  Py_ssize_t row_to, Py_ssize_t row_from, mod_int multiple,
                               Py_ssize_t start_col):
        cdef mod_int p
        cdef mod_int *v_from, *v_to

        p = self.p
        v_from = self.matrix[row_from]
        v_to = self.matrix[row_to]

        cdef Py_ssize_t i, nc
        nc = self._ncols
        for i from start_col <= i < nc:
            v_to[i] = (multiple * v_from[i] +  v_to[i]) % p

    cdef add_multiple_of_column_c(self, Py_ssize_t col_to, Py_ssize_t col_from, s,
                                   Py_ssize_t start_row):
        self._add_multiple_of_column_c(col_to, col_from, s, start_row)

    cdef _add_multiple_of_column_c(self, Py_ssize_t col_to, Py_ssize_t col_from, mod_int multiple,
                                   Py_ssize_t start_row):
        cdef mod_int  p
        cdef mod_int **m

        m = self.matrix
        p = self.p

        cdef Py_ssize_t i, nr
        nr = self._nrows
        for i from start_row <= i < self._nrows:
            m[i][col_to] = (m[i][col_to] + multiple * m[i][col_from]) %p

    cdef swap_rows_c(self, Py_ssize_t row1, Py_ssize_t row2):
        cdef mod_int* temp
        temp = self.matrix[row1]
        self.matrix[row1] = self.matrix[row2]
        self.matrix[row2] = temp

    cdef swap_columns_c(self, Py_ssize_t col1, Py_ssize_t col2):
        cdef Py_ssize_t i, nr
        cdef mod_int t
        cdef mod_int **m
        m = self.matrix
        nr = self._nrows
        for i from 0 <= i < self._nrows:
            t = m[i][col1]
            m[i][col1] = m[i][col2]
            m[i][col2] = t

    def hessenbergize(self):
        """
        Transforms self in place to its Hessenberg form.
        """
        self.check_mutability()
        x = self.fetch('in_hessenberg_form')
        if not x is None and x: return  # already known to be in Hessenberg form

        if self._nrows != self._ncols:
            raise ArithmeticError, "Matrix must be square to compute Hessenberg form."

        cdef Py_ssize_t n
        n = self._nrows

        cdef mod_int **h
        h = self.matrix

        cdef mod_int p, r, t, t_inv, u
        cdef Py_ssize_t i, j, m
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
                         self._add_multiple_of_row_c(j, m, p - u, 0)  # h[j] -= u*h[m]
                         # To maintain charpoly, do the corresponding
                         # column operation, which doesn't mess up the
                         # matrix, since it only changes column m, and
                         # we're only worried about column m-1 right
                         # now.  Add u*column_j to column_m.
                         self._add_multiple_of_column_c(m, j, u, 0)
                 # end for
            # end if
        # end for
        _sig_off
        self.cache('in_hessenberg_form',True)

    def _charpoly_hessenberg(self, var):
        """
        Transforms self in place to its Hessenberg form then computes
        and returns the coefficients of the characteristic polynomial of
        this matrix.

        INPUT:
            var -- name of the indeterminate of the charpoly.

        The characteristic polynomial is represented as a vector of
        ints, where the constant term of the characteristic polynomial
        is the 0th coefficient of the vector.
        """
        if self._nrows != self._ncols:
            raise ArithmeticError, "charpoly not defined for non-square matrix."

        cdef Py_ssize_t i, m, n,
        n = self._nrows

        cdef mod_int p, t
        p = self.p

        # Replace self by its Hessenberg form, and set H to this form
        # for notation below.
        cdef Matrix_modn_dense H
        H = self.copy()
        H.hessenbergize()


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
            c._add_multiple_of_row_c(m, m-1, p - H.matrix[m-1][m-1], 0)
            t = 1
            for i from 1 <= i < m:
                t = (t*H.matrix[m-i][m-i-1]) % p
                # Set the m-th row of c to c[m] - t*H[m-i-1,m-1]*c[m-i-1]
                c._add_multiple_of_row_c(m, m-i-1, p - (t*H.matrix[m-i-1][m-1])%p, 0)

        # The answer is now the n-th row of c.
        v = []
        for i from 0 <= i <= n:
            v.append(int(c.matrix[n][i]))
        R = self._base_ring[var]    # polynomial ring over the base ring
        return R(v)

    cdef int _strassen_default_cutoff(self, matrix0.Matrix right) except -2:
        # TODO: lots of testing
        return 100

    # TODO: TEMPORARILY DISABLED due to bug on 64-bit sage.math:
    #  A = matrix(Integers(389),4,range(16)); A._echelon_strassen(4)
    # *** glibc detected *** free(): invalid next size (fast): 0x0000000000fb15e0 ***
    # error in set_to memcpy on 64-bit
    def matrix_window(self, Py_ssize_t row=0, Py_ssize_t col=0,
                      Py_ssize_t nrows=-1, Py_ssize_t ncols=-1):
        """
        Return the requested matrix window.

        EXAMPLES:
            sage: a = matrix(QQ,3,range(9))
            sage: a.matrix_window()
            Matrix window of size 3 x 3 at (0,0):
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: type(a)
            <type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
            sage: a = matrix(GF(7),3,range(9)); a
            [0 1 2]
            [3 4 5]
            [6 0 1]
            sage: type(a)
            <type 'sage.matrix.matrix_modn_dense.Matrix_modn_dense'>
        """
        if nrows == -1:
            nrows = self._nrows - row
            ncols = self._ncols - col
        return matrix_window_modn_dense.MatrixWindow_modn_dense(self, row, col, nrows, ncols)


    def matrix_window_c(self, Py_ssize_t row, Py_ssize_t col,
                        Py_ssize_t nrows, Py_ssize_t ncols):
        if nrows == -1:
            nrows = self._nrows - row
            ncols = self._ncols - col
        return matrix_window_modn_dense.MatrixWindow_modn_dense(self, row, col, nrows, ncols)