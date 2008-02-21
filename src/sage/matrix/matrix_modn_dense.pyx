r"""
Dense matrices over $\Z/n\Z$ for $n$ small.

AUTHORS:
    -- William Stein
    -- Robert Bradshaw

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
    TypeError: unsupported operand parent(s) for '*': 'Full MatrixSpace of 3 by 3 dense matrices over Ring of integers modulo 37' and 'Full MatrixSpace of 2 by 3 dense matrices over Ring of integers modulo 37'
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

TESTS:
    sage: M = MatrixSpace(GF(5),2,2)
    sage: A = M([1,0,0,1])
    sage: A - int(-1)
    [2 0]
    [0 2]
    sage: B = M([4,0,0,1])
    sage: B - int(-1)
    [0 0]
    [0 2]
"""


include "../ext/interrupt.pxi"
include "../ext/cdefs.pxi"
include '../ext/stdsage.pxi'
include '../ext/random.pxi'

import sage.ext.multi_modular
cimport sage.ext.arith
import sage.ext.arith
cdef sage.ext.arith.arith_int ArithIntObj
ArithIntObj  = sage.ext.arith.arith_int()

MAX_MODULUS = sage.ext.multi_modular.MAX_MODULUS

# TODO: DO NOT change this back until all the ints, etc., below are changed
# and get_unsafe is rewritten to return the right thing.  E.g., with
# the above uncommented, on 64-bit,
# m =matrix(Integers(101^3),2,[824362, 606695, 641451, 205942])
# m.det()
#  --> gives 0, which is totally wrong.

import matrix_window_modn_dense

from sage.modules.vector_modn_dense cimport Vector_modn_dense

from sage.rings.arith import is_prime
from sage.structure.element cimport ModuleElement

cimport matrix_dense
cimport matrix
cimport matrix0
cimport sage.structure.element


from sage.libs.linbox.linbox cimport Linbox_modn_dense
cdef Linbox_modn_dense linbox
linbox = Linbox_modn_dense()

from sage.structure.element cimport Matrix

from sage.rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract

from sage.misc.misc import verbose, get_verbose, cputime

from sage.rings.integer import Integer

from sage.structure.element cimport ModuleElement, RingElement, Element, Vector

################
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
        if p >= MAX_MODULUS:
            raise OverflowError, "p (=%s) must be < %s"%(p, MAX_MODULUS)
        self.gather = MOD_INT_OVERFLOW/(p*p)

        _sig_on
        self._entries = <mod_int *> sage_malloc(sizeof(mod_int)*self._nrows*self._ncols)
        _sig_off
        if self._entries == NULL:
           raise MemoryError, "Error allocating matrix"

        self._matrix = <mod_int **> sage_malloc(sizeof(mod_int*)*self._nrows)
        if self._matrix == NULL:
            sage_free(self._entries)
            self._entries = NULL
            raise MemoryError, "Error allocating memory"

        cdef mod_int k
        cdef Py_ssize_t i
        k = 0
        for i from 0 <= i < self._nrows:
            self._matrix[i] = self._entries + k
            k = k + self._ncols

    def __dealloc__(self):
        if self._entries == NULL:
            return
        sage_free(self._entries)
        sage_free(self._matrix)

    def __init__(self, parent, entries, copy, coerce):

        cdef mod_int e
        cdef Py_ssize_t i, j, k
        cdef mod_int *v
        cdef mod_int p
        p = self._base_ring.characteristic()

        R = self.base_ring()

        # scalar?
        if not isinstance(entries, list):
            _sig_on
            for i from 0 <= i < self._nrows*self._ncols:
                self._entries[i] = 0
            _sig_off
            if entries is None:
                # zero matrix
                pass
            else:
                e = R(entries)
                if e != 0:
                    for i from 0 <= i < self._nrows:
                        self._matrix[i][i] = e
            return

        # all entries are given as a long list
        if len(entries) != self._nrows * self._ncols:
            raise IndexError, "The vector of entries has the wrong length."

        k = 0
        cdef mod_int n

        if coerce:
            for i from 0 <= i < self._nrows:
                if PyErr_CheckSignals(): raise KeyboardInterrupt
                v = self._matrix[i]
                for j from 0 <= j < self._ncols:
                    v[j] = R(entries[k])
                    k = k + 1
        else:
            for i from 0 <= i < self._nrows:
                if PyErr_CheckSignals(): raise KeyboardInterrupt
                v = self._matrix[i]
                for j from 0 <= j < self._ncols:
                    v[j] = int(entries[k])
                    k = k + 1


    def __richcmp__(Matrix_modn_dense self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)
    def __hash__(self):
        return self._hash()

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        self._matrix[i][j] = (<IntegerMod_int> value).ivalue

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef IntegerMod_int n
        n =  IntegerMod_int.__new__(IntegerMod_int)
        IntegerMod_abstract.__init__(n, self._base_ring)
        n.ivalue = self._matrix[i][j]
        return n

    ########################################################################
    # LEVEL 2 functionality
    #   * cdef _pickle
    #   * cdef _unpickle
    #   * cdef _add_c_impl
    #   * cdef _mul_c_impl
    #   * cdef _matrix_times_matrix_c_impl
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


    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        """
        Add two dense matrices over Z/nZ

        EXAMPLES:
            sage: a = MatrixSpace(GF(19),3)(range(9))
            sage: a+a
            [ 0  2  4]
            [ 6  8 10]
            [12 14 16]
            sage: b = MatrixSpace(GF(19),3)(range(9))
            sage: b.swap_rows(1,2)
            sage: a+b
            [ 0  2  4]
            [ 9 11 13]
            [ 9 11 13]
            sage: b+a
            [ 0  2  4]
            [ 9 11 13]
            [ 9 11 13]
        """

        cdef Py_ssize_t i,j
        cdef long k
        cdef Matrix_modn_dense M

        M = Matrix_modn_dense.__new__(Matrix_modn_dense, self._parent,None,None,None)
        Matrix.__init__(M, self._parent)

        _sig_on
        cdef mod_int *row_self, *row_right, *row_ans
        for i from 0 <= i < self._nrows:
            row_self = self._matrix[i]
            row_right = (<Matrix_modn_dense> right)._matrix[i]
            row_ans = M._matrix[i]
            for j from 0<= j < self._ncols:
                k = (row_self[j]+row_right[j])
                # division is really slow, so we normalize only when it is necessary
                if k>=self.p:
                    row_ans[j]=k-self.p
                else:
                    row_ans[j]=k

        _sig_off
        return M;


    cdef Matrix _matrix_times_matrix_c_impl(self, Matrix right):
        if get_verbose() >= 2:
            verbose('mod-p multiply of %s x %s matrix by %s x %s matrix modulo %s'%(
                self._nrows, self._ncols, right._nrows, right._ncols, self.p))
        if self._will_use_strassen(right):
            return self._multiply_strassen(right)
        else:
            return self._multiply_classical(right)

    def _multiply_linbox(Matrix_modn_dense self, Matrix_modn_dense right):
        """
        Multiply matrices using LinBox.

        INPUT:
            right -- Matrix

        """
        if get_verbose() >= 2:
            verbose('linbox multiply of %s x %s matrix by %s x %s matrix modulo %s'%(
                self._nrows, self._ncols, right._nrows, right._ncols, self.p))
        cdef int e
        cdef Matrix_modn_dense ans, B

        if not self.base_ring().is_field():
            raise ArithmeticError, "LinBox only supports fields"

        ans = self.new_matrix(nrows = self.nrows(), ncols = right.ncols())

        B = right
        self._init_linbox()
        _sig_on
        linbox.matrix_matrix_multiply(ans._matrix, B._matrix, B._nrows, B._ncols)
        _sig_off
        return ans

    def _multiply_classical(left, right):
        return left._multiply_strassen(right, left._ncols + left._nrows)

    cdef Vector _vector_times_matrix_c_impl(self, Vector v):
        cdef Vector_modn_dense w, ans
        cdef Py_ssize_t i, j
        cdef mod_int k
        cdef mod_int x

        M = self._row_ambient_module()
        w = v
        ans = M.zero_vector()

        for i from 0 <= i < self._ncols:
            x = 0
            k = 0
            for j from 0 <= j < self._nrows:
                x += w._entries[j] * self._matrix[j][i]
                k += 1
                if k >= self.gather:
                    x %= self.p
                    k = 0
            ans._entries[i] = x % self.p
        return ans


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


    def charpoly(self, var='x', algorithm='linbox'):
        """
        Returns the characteristic polynomial of self.

       INPUT:
            var -- a variable name
            algorithm -- 'linbox' (default if self.base_ring() is a field)
                         'generic'

        EXAMPLES:
            sage: A = Mat(GF(7),3,3)(range(3)*3)
            sage: A.charpoly()
            x^3 + 4*x^2

            sage: A = Mat(Integers(6),3,3)(range(9))
            sage: A.charpoly()
            x^3

        ALGORITHM: Uses LinBox if self.base_ring() is a field,
        otherwise use Hessenberg form algorithm.
        """
        if algorithm == 'linbox' and (self.p == 2 or not self.base_ring().is_field()):
            algorithm = 'generic' # LinBox only supports Z/pZ (p prime)

        if algorithm == 'linbox':
            g = self._charpoly_linbox(var)
        elif algorithm == 'generic':
            g = matrix_dense.Matrix_dense.charpoly(self, var)
        else:
            raise ValueError, "no algorithm '%s'"%algorithm
        self.cache('charpoly_%s_%s'%(algorithm, var), g)
        return g

    def minpoly(self, var='x', algorithm='linbox'):
        """
        Returns the minimal polynomial of self.

        INPUT:
            var -- a variable name
            algorithm -- 'linbox' (default if self.base_ring() is a field)
                         'generic'
        """
        if algorithm == 'linbox' and (self.p == 2 or not self.base_ring().is_field()):
            algorithm='generic' #LinBox only supports fields

        if algorithm == 'linbox':
            g = self._minpoly_linbox(var)
        elif algorithm == 'generic':
            #g = self._minpoly_generic(var)
            raise NotImplementedError, "minimal polynomials are not implemented for Z/nZ"
        else:
            raise ValueError, "no algorithm '%s'"%algorithm
        self.cache('minpoly_%s_%s'%(algorithm, var), g)
        return g

    def _minpoly_linbox(self, var='x'):
        """
        Computes the minimal polynomial using LinBox. No checks are
        performed.
        """
        verbose('_minpoly_linbox...')
        return self._poly_linbox(var=var, typ='minpoly')

    def _charpoly_linbox(self, var='x'):
        """
        Computes the characteristic polynomial using LinBox. No checks
        are performed.
        """
        verbose('_charpoly_linbox...')
        return self._poly_linbox(var=var, typ='charpoly')

    def _poly_linbox(self, var='x', typ='minpoly'):
        """
        Computes either the minimal or the characteristic polynomial
        using LinBox. No checks are performed.

        INPUT:
            var -- 'x'
            typ -- 'minpoly' or 'charpoly'
        """
        if self._nrows != self._ncols:
            raise ValueError, "matrix must be square"
        if self._nrows <= 1:
            return matrix_dense.Matrix_dense.charpoly(self, var)

        self._init_linbox()
        _sig_on
        if typ == 'minpoly':
            v = linbox.minpoly()
        else:
            v = linbox.charpoly()
        _sig_off
        R = self._base_ring[var]
        return R(v)

    def echelonize(self, algorithm="linbox", **kwds):
        r"""
        Puts self in row echelon form.

        INPUT:
            self -- a mutable matrix
            algorithm -- 'linbox' -- uses the C++ linbox library
                         'gauss'  -- uses a custom slower $O(n^3)$ Gauss
                                     elimination implemented in \sage.
                         'all' -- compute using both algorithms and verify
                                  that the results are the same (for the paranoid).
            **kwds -- these are all ignored

        OUTPUT:
            -- self is put in reduced row echelon form.
            -- the rank of self is computed and cached
            -- the pivot columns of self are computed and cached.
            -- the fact that self is now in echelon form is recorded
               and cached so future calls to echelonize return
               immediately.

        EXAMPLES:
            sage: a = matrix(GF(97),3,4,range(12))
            sage: a.echelonize(); a
            [ 1  0 96 95]
            [ 0  1  2  3]
            [ 0  0  0  0]
            sage: a.pivots()
            [0, 1]
        """
        x = self.fetch('in_echelon_form')
        if not x is None: return  # already known to be in echelon form
        if not self.base_ring().is_field():
            #self._echelon_mod_n ()
            raise NotImplementedError, "Echelon form not implemented over '%s'."%self.base_ring()

        self.check_mutability()
        self.clear_cache()

        if algorithm == 'linbox':
            self._echelonize_linbox()
        elif algorithm == 'gauss':
            self._echelon_in_place_classical()
        elif algorithm == 'all':
            A = self.copy()
            self._echelonize_linbox()
            A._echelon_in_place_classical()
            if A != self:
                raise RuntimeError, "bug in echelon form"
        else:
            raise ValueError, "algorithm '%s' not known"%algorithm

    cdef _init_linbox(self):
        _sig_on
        linbox.set(self.p, self._matrix, self._nrows, self._ncols)
        _sig_off

    def _echelonize_linbox(self):
        """
        Puts self in row echelon form using LinBox.
        """
        self.check_mutability()
        self.clear_cache()

        t = verbose('calling linbox echelonize mod %s'%self.p)
        self._init_linbox()
        _sig_on
        r = linbox.echelonize()
        _sig_off
        verbose('done with linbox echelonize',t)
        self.cache('in_echelon_form',True)
        self.cache('rank', r)
        self.cache('pivots', self._pivots())

    def _pivots(self):
        if not self.fetch('in_echelon_form'):
            raise RuntimeError, "self must be in reduced row echelon form first."
        pivots = []
        cdef Py_ssize_t i, j, nc
        nc = self._ncols
        cdef mod_int* row
        i = 0
        while i < self._nrows:
            row = self._matrix[i]
            for j from i <= j < nc:
                if row[j] != 0:
                    pivots.append(j)
                    i += 1
                    break
            if j == nc:
                break
        return pivots

    def _echelon_in_place_classical(self):
        self.check_mutability()
        self.clear_cache()

        cdef Py_ssize_t start_row, c, r, nr, nc, i
        cdef mod_int p, a, s, t, b
        cdef mod_int **m

        start_row = 0
        p = self.p
        m = self._matrix
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

    cdef xgcd_eliminate (self, mod_int * row1, mod_int* row2, Py_ssize_t start_col):
        """
        Reduces row1 and row2 by a unimodular transformation using the xgcd relation between their first
        coefficients a and b.

        INPUT:
            -- self: a mutable matrix
            -- row1, row2: the two rows to be transformed (within self)
            -- start_col: the column of the pivots in row1 and row2. It is assumed that all entries before
                          start_col in row1 and row2 are zero.

        OUTPUT:
            -- g: the gcd of the first elements of row1 and row2 at column start_col
            -- put row1 = s * row1 + t * row2
                   row2 = w * row1 + t * row2
               where g = sa + tb
        """
        cdef mod_int p = self.p
        cdef mod_int * row1_p, * row2_p
        cdef mod_int tmp
        cdef int g, s, t, v, w
        cdef Py_ssize_t nc, i
        cdef mod_int a =  row1[start_col]
        cdef mod_int b =  row2[start_col]
        g = ArithIntObj.c_xgcd_int (a,b,<int*>&s,<int*>&t)
        v = a/g
        w = -<int>b/g
        nc = self.ncols()

    #    print("In wgcd_eliminate")
        for i from start_col <= i < nc:
   #         print(self)
            tmp = ( s * <int>row1[i] + t * <int>row2[i]) % p
  #          print (tmp,s, <int>row1[i],t,<int>row2[i])
 #           print (row2[i],w, <int>row1[i],v,<int>row2[i])
            row2[i] = (w* <int>row1[i] + v*<int>row2[i]) % p
#            print (row2[i],w, <int>row1[i],v,<int>row2[i])
            row1[i] = tmp
        #print(self)
       # print("sortie")
        return g


## This is code by Me and/or Clement Pernet from SD7.  Unfortunately I
## think it is still buggy, since it is so painful to implement with
## unsigned ints. Code to do basically the same thing is in
## matrix_integer_dense, by Erocal Burcin.
##     def _echelon_mod_n (self):
##         """
##         Put self in Hermite normal form modulo n
##         (echelonize the matrix over a ring $Z_n$)

##         INPUT:
##             self: a mutable matrix over $Z_n$
##         OUTPUT:
##             Transform in place the working matrix into its Hermite
##             normal form over Z, using the modulo n algorithm of
##             [Hermite Normal form computation using modulo determinant arithmetic,
##             Domich Kannan & Trotter, 1987]
##         """
##         self.check_mutability()
##         self.clear_cache()

##         cdef Py_ssize_t start_row, nr, nc,
##         cdef long c, r, i
##         cdef mod_int p, a, a_inverse, b, g
##         cdef mod_int **m
##         cdef Py_ssize_t start_row = 0
##         p = self.p
##         m = self._matrix
##         nr = self._nrows
##         nc = self._ncols
##         pivots = []
##         cdef Py_ssize_t fifth = self._ncols / 10 + 1
##         do_verb = (get_verbose() >= 2)
##         for c from 0 <= c < nc:
##             if do_verb and (c % fifth == 0 and c>0):
##                 tm = verbose('on column %s of %s'%(c, self._ncols),
##                              level = 2,
##                              caller_name = 'matrix_modn_dense echelon mod n')

##             if PyErr_CheckSignals(): raise KeyboardInterrupt
##             for r from start_row <= r < nr:
##                 a = m[r][c]
##                 if a:
##                     self.swap_rows_c(r, start_row)
##                     for i from start_row +1 <= i < nr:
##                         b = m[i][c]

##                         if b != 0:
##                             self.xgcd_eliminate (self._matrix[start_row], self._matrix[i], c)
##                             verbose('eliminating rows %s and %s', (start_row,i))
##                     for i from 0 <= i <start_row:
##                         p = -m[i][c]//m[start_row][c]
##                         self._add_multiple_of_row_c(i, start_row, p, c)

##                     pivots.append(m[start_row][c])
##                     start_row = start_row + 1
##                     break
##         self.cache('pivots',pivots)
##         self.cache('in_echelon_form',True)


    cdef rescale_row_c(self, Py_ssize_t row, multiple, Py_ssize_t start_col):
        self._rescale_row_c(row, multiple, start_col)

    cdef _rescale_row_c(self, Py_ssize_t row, mod_int multiple, Py_ssize_t start_col):
        cdef mod_int r, p
        cdef mod_int* v
        cdef Py_ssize_t i
        p = self.p
        v = self._matrix[row]
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
            self._matrix[i][col] = (self._matrix[i][col]*multiple) % p

    cdef add_multiple_of_row_c(self,  Py_ssize_t row_to, Py_ssize_t row_from, multiple,
                               Py_ssize_t start_col):
        self._add_multiple_of_row_c(row_to, row_from, multiple, start_col)

    cdef _add_multiple_of_row_c(self,  Py_ssize_t row_to, Py_ssize_t row_from, mod_int multiple,
                               Py_ssize_t start_col):
        cdef mod_int p
        cdef mod_int *v_from, *v_to

        p = self.p
        v_from = self._matrix[row_from]
        v_to = self._matrix[row_to]

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

        m = self._matrix
        p = self.p

        cdef Py_ssize_t i, nr
        nr = self._nrows
        for i from start_row <= i < self._nrows:
            m[i][col_to] = (m[i][col_to] + multiple * m[i][col_from]) %p

    cdef swap_rows_c(self, Py_ssize_t row1, Py_ssize_t row2):
        cdef mod_int* temp
        temp = self._matrix[row1]
        self._matrix[row1] = self._matrix[row2]
        self._matrix[row2] = temp

    cdef swap_columns_c(self, Py_ssize_t col1, Py_ssize_t col2):
        cdef Py_ssize_t i, nr
        cdef mod_int t
        cdef mod_int **m
        m = self._matrix
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
        h = self._matrix

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
        c._matrix[0][0] = 1
        for m from 1 <= m <= n:
            # Set the m-th row of c to (x - H[m-1,m-1])*c[m-1] = x*c[m-1] - H[m-1,m-1]*c[m-1]
            # We do this by hand by setting the m-th row to c[m-1]
            # shifted to the right by one.  We then add
            # -H[m-1,m-1]*c[m-1] to the resulting m-th row.
            for i from 1 <= i <= n:
                c._matrix[m][i] = c._matrix[m-1][i-1]
            # the p-.. below is to keep scalar normalized between 0 and p.
            c._add_multiple_of_row_c(m, m-1, p - H._matrix[m-1][m-1], 0)
            t = 1
            for i from 1 <= i < m:
                t = (t*H._matrix[m-i][m-i-1]) % p
                # Set the m-th row of c to c[m] - t*H[m-i-1,m-1]*c[m-i-1]
                c._add_multiple_of_row_c(m, m-i-1, p - (t*H._matrix[m-i-1][m-1])%p, 0)

        # The answer is now the n-th row of c.
        v = []
        for i from 0 <= i <= n:
            v.append(int(c._matrix[n][i]))
        R = self._base_ring[var]    # polynomial ring over the base ring
        return R(v)

    def rank(self):
        """
        Return the rank of this matrix.

        EXAMPLES:
            sage: m = matrix(GF(7),5,range(25))
            sage: m.rank()
            2

        Rank is not implemented over the integers modulo a composite yet.
            sage: m = matrix(Integers(4), 2, [2,2,2,2])
            sage: m.rank()
            Traceback (most recent call last):
            ...
            NotImplementedError: Echelon form not implemented over 'Ring of integers modulo 4'.
        """
        cdef Matrix_modn_dense A
        if self.p > 2 and is_prime(self.p):
            x = self.fetch('rank')
            if not x is None:
                return x
            # avoid modifying self in place!
            # TODO: it is crappy/buggy that linbox would change a matrix
            # when the rank function is called on it... Sigh.
            A = self.__copy__()
            A._init_linbox()
            _sig_on
            r = Integer(linbox.rank())
            _sig_off
            self.cache('rank', r)
            return r
        else:
            # linbox is very buggy for p=2
            return matrix_dense.Matrix_dense.rank(self)

    def determinant(self):
        """
        Return the determinant of this matrix.

        EXAMPLES:
            sage: m = matrix(GF(101),5,range(25))
            sage: m.det()
            0

            sage: m = matrix(Integers(4), 2, [2,2,2,2])
            sage: m.det()
            0

        """
        if self.p > 2 and is_prime(self.p):
            x = self.fetch('det')
            if not x is None:
                return x
            self._init_linbox()
            _sig_on
            d = linbox.det()
            _sig_off
            d2 = self._coerce_element(d)
            self.cache('det', d2)
            return d2
        else:
            return matrix_dense.Matrix_dense.determinant(self)

    def randomize(self, density=1):
        """
        Randomize density proportion of the entries of this matrix,
        leaving the rest unchanged.

        NOTE: The random() function doesn't seem to be as random as
        expected. The lower-order bits seem to have a strong bias
        towards zero. Even stranger, if you create two 1000x1000
        matrices over GF(2) they always appear to have the same
        reduced row echelon form, i.e. they span the same space. The
        higher-order bits seem to be much more random and thus we
        shift first and mod p then.
        """

        density = float(density)
        if density == 0:
            return

        self.check_mutability()
        self.clear_cache()

        cdef int nc
        if density == 1:
            for i from 0 <= i < self._nrows*self._ncols:
                # 16-bit seems safe
                self._entries[i] =  (random()>>16) % self.p
        else:
            density = float(density)
            nc = self._ncols
            num_per_row = int(density * nc)
            _sig_on
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < num_per_row:
                    k = ((random()>>16)+random())%nc # is this safe?
                    # 16-bit seems safe
                    self._matrix[i][k] = (random()>>16) % self.p
            _sig_off

    cdef int _strassen_default_cutoff(self, matrix0.Matrix right) except -2:
        # TODO: lots of testing
        return 100

    cdef matrix_window_c(self, Py_ssize_t row, Py_ssize_t col,
                        Py_ssize_t nrows, Py_ssize_t ncols):
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

    def _magma_init_(self):
        """
        Returns a string of self in MAGMA form.

        NOTE: Does not return MAGMA object but string.
        """
        cdef int i,j
        K = self._base_ring._magma_init_()
        if self._nrows == self._ncols:
            s = 'MatrixAlgebra(%s, %s)'%(K, self.nrows())
        else:
            s = 'RMatrixSpace(%s, %s, %s)'%(K, self.nrows(), self.ncols())
        v = []
        for i from 0 <= i < self._nrows*self._ncols:
                v.append(str(self._entries[i]))
        return s + '![%s]'%(','.join(v))
