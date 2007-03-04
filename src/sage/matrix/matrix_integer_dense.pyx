"""
Dense matrices over the integer ring.

AUTHORS:
    -- William Stein
    -- Robert Bradshaw

EXAMPLES:
    sage: a = matrix(ZZ, 3,3, range(9)); a
    [0 1 2]
    [3 4 5]
    [6 7 8]
    sage: a.det()
    0
    sage: a[0,0] = 10; a.det()
    -30
    sage: a.charpoly()
    x^3 - 22*x^2 + 102*x + 30
    sage: b = -3*a
    sage: a == b
    False
    sage: b < a
    True
"""

########## *** IMPORTANT ***
# If you're working on this code, we *always* assume that
#   self._matrix[i] = self._entries[i*self._ncols]
# !!!!!!!! Do not break this!
# This is assumed in the _rational_kernel_iml

######################################################################
#       Copyright (C) 2006,2007 William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
######################################################################

from sage.modules.vector_integer_dense cimport Vector_integer_dense

from sage.misc.misc import verbose, get_verbose, cputime

from sage.rings.arith import previous_prime

include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/gmp.pxi"

ctypedef unsigned int uint

from sage.ext.multi_modular import MultiModularBasis
from sage.ext.multi_modular cimport MultiModularBasis

from sage.rings.integer cimport Integer
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.integer_mod_ring import IntegerModRing
from sage.rings.polynomial_ring import PolynomialRing
from sage.structure.element cimport ModuleElement, RingElement, Element, Vector
from sage.structure.sequence import Sequence

from matrix_modn_dense import Matrix_modn_dense
from matrix_modn_dense cimport Matrix_modn_dense

import sage.modules.free_module
import sage.modules.free_module_element

from matrix cimport Matrix

cimport sage.structure.element

import matrix_space


######### linbox interface ##########
from sage.libs.linbox.linbox cimport Linbox_integer_dense
cdef Linbox_integer_dense linbox
linbox = Linbox_integer_dense()
USE_LINBOX_POLY = True


########## iml -- integer matrix library ###########

cdef extern from "iml.h":

    enum SOLU_POS:
        LeftSolu = 101
        RightSolu = 102

    long nullspaceMP(long n, long m,
                     mpz_t *A, mpz_t * *mp_N_pass)

    void nonsingSolvLlhsMM (SOLU_POS solupos, long n, \
                       long m, mpz_t *mp_A, mpz_t *mp_B, mpz_t *mp_N, \
                       mpz_t mp_D)



cdef class Matrix_integer_dense(matrix_dense.Matrix_dense):   # dense or sparse
    r"""
    Matrix over the integers.

    On a 32-bit machine, they can have at most $2^{32}-1$ rows or
    columns.  On a 64-bit machine, matrices can have at most
    $2^{64}-1$ rows or columns.

    EXAMPLES:
        sage: a = MatrixSpace(ZZ,3)(2); a
        [2 0 0]
        [0 2 0]
        [0 0 2]
        sage: a = matrix(ZZ,1,3, [1,2,-3]); a
        [ 1  2 -3]
        sage: a = MatrixSpace(ZZ,2,4)(2); a
        Traceback (most recent call last):
        ...
        TypeError: nonzero scalar matrix must be square
    """
    ########################################################################
    # LEVEL 1 functionality
    # x * __new__
    # x * __dealloc__
    # x * __init__
    # x * set_unsafe
    # x * get_unsafe
    # x * def _pickle
    # x * def _unpickle
    ########################################################################

    def __new__(self, parent, entries, coerce, copy):
        """
        Create and allocate memory for the matrix.  Does not actually initialize
        any of the memory.

        INPUT:
            parent, entries, coerce, copy -- as for __init__.

        EXAMPLES:
            sage: from sage.matrix.matrix_integer_dense import Matrix_integer_dense
            sage: a = Matrix_integer_dense.__new__(Matrix_integer_dense, Mat(ZZ,3), 0,0,0)
            sage: type(a)
            <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>

        WARNING: This is for internal use only, or if you really know what you're doing.
        """
        matrix_dense.Matrix_dense.__init__(self, parent)
        self._nrows = parent.nrows()
        self._ncols = parent.ncols()
        self._pivots = None

        # Allocate an array where all the entries of the matrix are stored.
        _sig_on
        self._entries = <mpz_t *>sage_malloc(sizeof(mpz_t) * (self._nrows * self._ncols))
        _sig_off
        if self._entries == NULL:
            raise MemoryError, "out of memory allocating a matrix"

        # Allocate an array of pointers to the rows, which is useful for
        # certain algorithms.
        ##################################
        # *IMPORTANT*: FOR MATRICES OVER ZZ, WE ALWAYS ASSUME THAT
        # THIS ARRAY IS *not* PERMUTED.  This should be OK, since all
        # algorithms are multi-modular.
        ##################################
        self._matrix = <mpz_t **> sage_malloc(sizeof(mpz_t*)*self._nrows)
        if self._matrix == NULL:
            sage_free(self._entries)
            self._entries = NULL
            raise MemoryError, "out of memory allocating a matrix"

        # Set each of the pointers in the array self._matrix to point
        # at the memory for the corresponding row.
        cdef Py_ssize_t i, k
        k = 0
        for i from 0 <= i < self._nrows:
            self._matrix[i] = self._entries + k
            k = k + self._ncols

    cdef _init_linbox(self):
        _sig_on
        linbox.set(self._matrix, self._nrows, self._ncols)
        _sig_off

    def __copy__(self):
        cdef Matrix_integer_dense A
        A = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent,
                                         0, 0, 0)
        cdef Py_ssize_t i
        _sig_on
        for i from 0 <= i < self._nrows * self._ncols:
            mpz_init_set(A._entries[i], self._entries[i])
        _sig_off
        A._initialized = True
        return A

    def __dealloc__(self):
        """
        Frees all the memory allocated for this matrix.
        EXAMPLE:
            sage: a = Matrix(ZZ,2,[1,2,3,4])
            sage: del a
        """
        cdef Py_ssize_t i
        if self._initialized:
            for i from 0 <= i < (self._nrows * self._ncols):
                mpz_clear(self._entries[i])
            sage_free(self._entries)
            sage_free(self._matrix)

    def __init__(self, parent, entries, copy, coerce):
        r"""
        Initialize a dense matrix over the integers.

        INPUT:
            parent -- a matrix space
            entries -- list - create the matrix with those entries along the rows.
                       other -- a scalar; entries is coerced to an integer and the diagonal
                                entries of this matrix are set to that integer.
            coerce -- whether need to coerce entries to the integers (program may crash
                      if you get this wrong)
            copy -- ignored (since integers are immutable)

        EXAMPLES:

        The __init__ function is called implicitly in each of the
        examples below to actually fill in the values of the matrix.

        We create a $2 \times 2$ and a $1\times 4$ matrix:
            sage: matrix(ZZ,2,2,range(4))
            [0 1]
            [2 3]
            sage: Matrix(ZZ,1,4,range(4))
            [0 1 2 3]

        If the number of columns isn't given, it is determined from the number of
        elements in the list.
            sage: matrix(ZZ,2,range(4))
            [0 1]
            [2 3]
            sage: matrix(ZZ,2,range(6))
            [0 1 2]
            [3 4 5]

        Another way to make a matrix is to create the space of
        matrices and coerce lists into it.
            sage: A = Mat(ZZ,2); A
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: A(range(4))
            [0 1]
            [2 3]

        Actually it is only necessary that the input can be coerced to a list, so
        the following also works:
            sage: v = reversed(range(4)); type(v)
            <type 'listreverseiterator'>
            sage: A(v)
            [3 2]
            [1 0]

        Matrices can have many rows or columns (in fact, on a 64-bit machine they could
        have up to $2^64-1$ rows or columns):
            sage: v = matrix(ZZ,1,10^5, range(10^5))
            sage: v.parent()
            Full MatrixSpace of 1 by 100000 dense matrices over Integer Ring
        """
        cdef Py_ssize_t i, j
        cdef int is_list
        cdef Integer x

        if not isinstance(entries, list):
            try:
                entries = list(entries)
                is_list = 1
            except TypeError:
                try:
                    # Try to coerce entries to a scalar (an integer)
                    x = ZZ(entries)
                    is_list = 0
                except TypeError:
                    raise TypeError, "entries must be coercible to a list or integer"
        else:
            is_list = 1

        if is_list:

            # Create the matrix whose entries are in the given entry list.
            if len(entries) != self._nrows * self._ncols:
                sage_free(self._entries)
                sage_free(self._matrix)
                self._entries = NULL
                raise TypeError, "entries has the wrong length"
            if coerce:
                for i from 0 <= i < self._nrows * self._ncols:
                    # TODO: Should use an unsafe un-bounds-checked array access here.
                    x = ZZ(entries[i])
                    # todo -- see integer.pyx and the TODO there; perhaps this could be
                    # sped up by creating a mpz_init_set_sage function.
                    mpz_init_set(self._entries[i], x.value)
                self._initialized = True
            else:
                for i from 0 <= i < self._nrows * self._ncols:
                    # TODO: Should use an unsafe un-bounds-checked array access here.
                    mpz_init_set(self._entries[i], (<Integer> entries[i]).value)
                self._initialized = True
        else:

            # If x is zero, make the zero matrix and be done.
            if mpz_sgn(x.value) == 0:
                self._zero_out_matrix()
                return

            # the matrix must be square:
            if self._nrows != self._ncols:
                sage_free(self._entries)
                sage_free(self._matrix)
                self._entries = NULL
                raise TypeError, "nonzero scalar matrix must be square"

            # Now we set all the diagonal entries to x and all other entries to 0.
            self._zero_out_matrix()
            j = 0
            for i from 0 <= i < self._nrows:
                mpz_init_set(self._entries[j], x.value)
                j = j + self._nrows + 1
            self._initialized = True

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        """
        Set position i,j of this matrix to x.

        (VERY UNSAFE -- value *must* be of type Integer).

        INPUT:
            ij -- tuple (i,j), where i is the row and j the column
        Alternatively, ij can be an integer, and the ij-th row is set.

        EXAMPLES:
            sage: a = matrix(ZZ,2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a[0,0] = 10
            sage: a
            [10  1  2]
            [ 3  4  5]
        """
        #cdef Integer Z
        #Z = value
        #mpz_set(self._matrix[i][j], Z.value)
        mpz_set(self._matrix[i][j], (<Integer>value).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        EXAMPLES:
            sage: a = MatrixSpace(ZZ,3)(range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a[1,2]
            5
            sage: a[0]
            (0, 1, 2)
            sage: a[4,7]
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range
            sage: a[-1,0]
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range
        """
        cdef Integer z
        z = PY_NEW(Integer)
        mpz_set(z.value, self._matrix[i][j])
        return z

    def _pickle(self):
        return self._pickle_version0(), 0

    cdef _pickle_version0(self):
        # TODO: *maybe* redo this to use mpz_import and mpz_export
        # from sec 5.14 of the GMP manual. ??
        cdef int i, j, len_so_far, m, n
        cdef char *a
        cdef char *s, *t, *tmp

        if self._nrows == 0 or self._ncols == 0:
            data = ''
        else:
            n = self._nrows*self._ncols*10
            s = <char*> sage_malloc(n * sizeof(char))
            t = s
            len_so_far = 0

            _sig_on
            for i from 0 <= i < self._nrows * self._ncols:
                m = mpz_sizeinbase (self._entries[i], 32)
                if len_so_far + m + 1 >= n:
                    # copy to new string with double the size
                    n = 2*n + m + 1
                    tmp = <char*> sage_malloc(n * sizeof(char))
                    strcpy(tmp, s)
                    sage_free(s)
                    s = tmp
                    t = s + len_so_far
                #endif
                mpz_get_str(t, 32, self._entries[i])
                m = strlen(t)
                len_so_far = len_so_far + m + 1
                t = t + m
                t[0] = <char>32
                t[1] = <char>0
                t = t + 1
            _sig_off
            data = str(s)[:-1]
            free(s)
        return data

    def _unpickle(self, data, int version):
        if version == 0:
            self._unpickle_version0(data)
        else:
            raise RuntimeError, "unknown matrix version (=%s)"%version

    cdef _unpickle_version0(self, data):
        cdef Py_ssize_t i, n
        data = data.split()
        n = self._nrows * self._ncols
        if len(data) != n:
            raise RuntimeError, "invalid pickle data."
        for i from 0 <= i < n:
            s = data[i]
            if mpz_init_set_str(self._entries[i], s, 32):
                raise RuntimeError, "invalid pickle data"
        self._initialized = True


    def __richcmp__(Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)

    def __hash__(self):
        return self._hash()

    ########################################################################
    # LEVEL 1 helpers:
    #   These function support the implementation of the level 1 functionality.
    ########################################################################
    cdef _zero_out_matrix(self):
        """
        Set this matrix to be the zero matrix.
        This is only for internal use.
        """
        # TODO: This is about 6-10 slower than MAGMA doing what seems to be the same thing.
        # Moreover, NTL can also do this quickly.  Why?   I think both have specialized
        # small integer classes.
        _sig_on
        cdef Py_ssize_t i
        for i from 0 <= i < self._nrows * self._ncols:
            mpz_init(self._entries[i])
        _sig_off
        self._initialized = True

    cdef _new_unitialized_matrix(self, Py_ssize_t nrows, Py_ssize_t ncols):
        """
        Return a new matrix over the integers with the given number of rows and columns.
        All memory is allocated for this matrix, but its entries have not yet been
        filled in.
        """
        P = self._parent.matrix_space(nrows, ncols)
        return Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)


    ########################################################################
    # LEVEL 2 functionality
    # x * cdef _add_c_impl
    # x * cdef _sub_c_impl
    # x * cdef _mul_c_impl
    # x * cdef _cmp_c_impl
    #   * __neg__
    #   * __invert__
    #   * __copy__
    # x * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################

    # cdef _mul_c_impl(self, Matrix right):
    # def __neg__(self):
    # def __invert__(self):
    # def __copy__(self):
    # def _multiply_classical(left, matrix.Matrix _right):
    # def _list(self):
    # def _dict(self):

    def is_zero(self):
        cdef mpz_t *a, *b
        cdef Py_ssize_t i, j
        cdef int k
        for i from 0 <= i < self._nrows * self._ncols:
            if mpz_cmp_si(self._entries[i], 0):
                return False
        return True

    def _multiply_linbox(self, Matrix right):
        """
        Multiply matrices over ZZ using linbox.

        WARNING: This is very slow right now, i.e., linbox is very slow.

        EXAMPLES:
            sage: A = matrix(ZZ,2,3,range(6))
            sage: A*A.transpose()
            [ 5 14]
            [14 50]
            sage: A._multiply_linbox(A.transpose())
            [ 5 14]
            [14 50]
        """
        cdef int e
        cdef Matrix_integer_dense ans, B
        B = right
        ans = self.new_matrix(nrows = self.nrows(), ncols = right.ncols())
        self._init_linbox()
        _sig_on
        linbox.matrix_matrix_multiply(ans._matrix, B._matrix, B._nrows, B._ncols)
        _sig_off
        return ans

    def _multiply_classical(self, Matrix right):
        """
        EXAMPLE:
            sage: n = 3
            sage: a = MatrixSpace(ZZ,n,n)(range(n^2))
            sage: a*a
            [ 15  18  21]
            [ 42  54  66]
            [ 69  90 111]
        """
        if self._ncols != right._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of right."

        cdef Py_ssize_t i, j, k, l, nr, nc, snc
        cdef mpz_t *v

        nr = self._nrows
        nc = right._ncols
        snc = self._ncols

        parent = self.matrix_space(nr, nc)

        cdef Matrix_integer_dense M, _right
        _right = right

        M = Matrix_integer_dense.__new__(Matrix_integer_dense, parent, None, None, None)
        Matrix.__init__(M, parent)

        # M has memory allocated but entries are not initialized
        cdef mpz_t *entries
        entries = M._entries

        cdef mpz_t s, z
        mpz_init(s)
        mpz_init(z)

        _sig_on
        l = 0
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                mpz_set_si(s,0)   # set s = 0
                v = self._matrix[i]
                for k from 0 <= k < snc:
                    mpz_mul(z, v[k], _right._matrix[k][j])
                    mpz_add(s, s, z)
                mpz_init_set(entries[l], s)
                l += 1
        _sig_off
        mpz_clear(s)
        mpz_clear(z)
        M._initialized = True
        return M

    cdef sage.structure.element.Matrix _matrix_times_matrix_c_impl(self, sage.structure.element.Matrix right):

        #############
        # see the tune_multiplication function below.
        n = max(self._nrows, self._ncols, right._nrows, right._ncols)
        if n <= 20:
            return self._multiply_classical(right)
        a = self.height(); b = right.height()
        if float(max(a,b)) / float(n) >= 0.70:
            return self._multiply_classical(right)
        else:
            return self._multiply_multi_modular(right)

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        """
        EXAMPLES:
            sage: a = matrix(QQ,2,range(6))
            sage: (3/4) * a
            [   0  3/4  3/2]
            [ 9/4    3 15/4]
        """
        cdef Py_ssize_t i
        cdef Integer _x
        _x = Integer(right)
        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent, None, None, None)
        for i from 0 <= i < self._nrows * self._ncols:
            mpz_init(M._entries[i])
            mpz_mul(M._entries[i], self._entries[i], _x.value)
        M._initialized = True
        return M

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        """
        Add two dense matrices over ZZ.

        EXAMPLES:
            sage: a = MatrixSpace(ZZ,3)(range(9))
            sage: a+a
            [ 0  2  4]
            [ 6  8 10]
            [12 14 16]
            sage: b = MatrixSpace(ZZ,3)(range(9))
            sage: b.swap_rows(1,2)
            sage: a+b
            [ 0  2  4]
            [ 9 11 13]
            [ 9 11 13]
        """
        cdef Py_ssize_t i, j

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent, None, None, None)

        _sig_on
        cdef mpz_t *row_self, *row_right, *row_ans
        for i from 0 <= i < self._nrows:
            row_self = self._matrix[i]
            row_right = (<Matrix_integer_dense> right)._matrix[i]
            row_ans = M._matrix[i]
            for j from 0 <= j < self._ncols:
                mpz_init(row_ans[j])
                mpz_add(row_ans[j], row_self[j], row_right[j])
        _sig_off
        M._initialized = True
        return M

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Subtract two dense matrices over ZZ.

        EXAMPLES:
            sage: M = Mat(ZZ,3)
            sage: a = M(range(9)); b = M(reversed(range(9)))
            sage: a - b
            [-8 -6 -4]
            [-2  0  2]
            [ 4  6  8]
        """
        cdef Py_ssize_t i, j

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent, None, None, None)

        _sig_on
        cdef mpz_t *row_self, *row_right, *row_ans
        for i from 0 <= i < self._nrows:
            row_self = self._matrix[i]
            row_right = (<Matrix_integer_dense> right)._matrix[i]
            row_ans = M._matrix[i]
            for j from 0 <= j < self._ncols:
                mpz_init(row_ans[j])
                mpz_sub(row_ans[j], row_self[j], row_right[j])
        _sig_off
        M._initialized = True
        return M


    cdef int _cmp_c_impl(self, Element right) except -2:
        cdef mpz_t *a, *b
        cdef Py_ssize_t i, j
        cdef int k
        for i from 0 <= i < self._nrows:
            a = self._matrix[i]
            b = (<Matrix_integer_dense>right)._matrix[i]
            for j from 0 <= j < self._ncols:
                k = mpz_cmp(a[j], b[j])
                if k:
                    if k < 0:
                        return -1
                    else:
                        return 1
        return 0


    cdef Vector _vector_times_matrix_c_impl(self, Vector v):
        """
        Returns the vector times matrix product.

        INPUT:
             v -- a free module element.

        OUTPUT:
            The the vector times matrix product v*A.

        EXAMPLES:
            sage: B = matrix(ZZ,2, [1,2,3,4])
            sage: V = ZZ^2
            sage: w = V([-1,5])
            sage: w*B
            (14, 18)
        """
        cdef Vector_integer_dense w, ans
        cdef Py_ssize_t i, j
        cdef mpz_t x

        M = self._row_ambient_module()
        w = <Vector_integer_dense> v
        ans = M.zero_vector()

        mpz_init(x)
        for i from 0 <= i < self._ncols:
            mpz_set_si(x, 0)
            for j from 0 <= j < self._nrows:
                mpz_addmul(x, w._entries[j], self._matrix[j][i])
            mpz_set(ans._entries[i], x)
        mpz_clear(x)
        return ans


    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_c_impl
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    #    * Specialized echelon form
    ########################################################################

    def charpoly(self, var='x', algorithm='linbox'):
        """
        INPUT:
            var -- a variable name
            algorithm -- 'linbox' (default)
                         'generic'

        NOTE: Linbox charpoly disabled on 64-bit machines, since it
        hangs in many cases.

        EXAMPLES:
            sage: A = matrix(ZZ,6, range(36))
            sage: f = A.charpoly(); f
            x^6 - 105*x^5 - 630*x^4
            sage: f(A) == 0
            True
            sage: n=20; A = Mat(ZZ,n)(range(n^2))
            sage: A.charpoly()
            x^20 - 3990*x^19 - 266000*x^18
            sage: A.minpoly()                # optional -- os x only right now
            x^3 - 3990*x^2 - 266000*x
        """
        key = 'charpoly_%s_%s'%(algorithm, var)
        x = self.fetch(key)
        if x: return x

        if algorithm == 'linbox' and not USE_LINBOX_POLY:
            algorithm = 'generic'

        if algorithm == 'linbox':
            g = self._charpoly_linbox(var)
        elif algorithm == 'generic':
            g = matrix_dense.Matrix_dense.charpoly(self, var)
        else:
            raise ValueError, "no algorithm '%s'"%algorithm

        self.cache(key, g)
        return g

    def minpoly(self, var='x', algorithm='linbox'):
        """
        INPUT:
            var -- a variable name
            algorithm -- 'linbox' (default)
                         'generic'

        NOTE: Linbox charpoly disabled on 64-bit machines, since it
        hangs in many cases.

        EXAMPLES:
            sage: A = matrix(ZZ,6, range(36))
            sage: A.minpoly()           # optional -- os x only right now
            x^3 - 105*x^2 - 630*x
            sage: n=6; A = Mat(ZZ,n)([k^2 for k in range(n^2)])
            sage: A.minpoly()           # optional -- os x only right now
            x^4 - 2695*x^3 - 257964*x^2 + 1693440*x
        """
        key = 'minpoly_%s_%s'%(algorithm, var)
        x = self.fetch(key)
        if x: return x

        if algorithm == 'linbox' and not USE_LINBOX_POLY:
            algorithm = 'generic'

        if algorithm == 'linbox':
            g = self._minpoly_linbox(var)
        elif algorithm == 'generic':
            g = matrix_dense.Matrix_dense.minpoly(self, var)
        else:
            raise ValueError, "no algorithm '%s'"%algorithm

        self.cache(key, g)
        return g

    def _minpoly_linbox(self, var='x'):
        return self._poly_linbox(var=var, typ='minpoly')

    def _charpoly_linbox(self, var='x'):
        if self.is_zero():  # program around a bug in linbox on 32-bit linux
            x = self.base_ring()[var].gen()
            return x ** self._nrows
        return self._poly_linbox(var=var, typ='charpoly')

    def _poly_linbox(self, var='x', typ='minpoly'):
        """
        INPUT:
            var -- 'x'
            typ -- 'minpoly' or 'charpoly'
        """
        time = verbose('computing %s of %s x %s matrix using linbox'%(typ, self._nrows, self._ncols))
        if self._nrows != self._ncols:
            raise ValueError, "matrix must be square"
        if self._nrows <= 1:
            return matrix_dense.Matrix_dense.charpoly(self, var)
        self._init_linbox()
        if typ == 'minpoly':
            _sig_on
            v = linbox.minpoly()
            _sig_off
        else:
            _sig_on
            v = linbox.charpoly()
            _sig_off
        R = self._base_ring[var]
        verbose('finished computing %s'%typ, time)
        return R(v)


    def height(self):
        """
        Return the height of this matrix, i.e., the max absolute value
        of the entries of the matrix.

        OUTPUT:
            A nonnegative integer.

        EXAMPLE:
            sage: a = Mat(ZZ,3)(range(9))
            sage: a.height()
            8
            sage: a = Mat(ZZ,2,3)([-17,3,-389,15,-1,0]); a
            [ -17    3 -389]
            [  15   -1    0]
            sage: a.height()
            389
        """
        cdef mpz_t h
        cdef Integer x

        self.mpz_height(h)
        x = Integer()
        x.set_from_mpz(h)
        mpz_clear(h)

        return x

    cdef int mpz_height(self, mpz_t height) except -1:
        """
        Used to compute the height of this matrix.

        INPUT:
             height -- a GMP mpz_t (that has not been initialized!)
        OUTPUT:
             sets the value of height to the height of this matrix, i.e., the max absolute
             value of the entries of the matrix.
        """
        cdef mpz_t x, h
        cdef Py_ssize_t i

        mpz_init_set_si(h, 0)
        mpz_init(x)

        _sig_on

        for i from 0 <= i < self._nrows * self._ncols:
            mpz_abs(x, self._entries[i])
            if mpz_cmp(h, x) < 0:
                mpz_set(h, x)

        _sig_off

        mpz_init_set(height, h)
        mpz_clear(h)
        mpz_clear(x)

        return 0   # no error occured.

    def _multiply_multi_modular(left, Matrix_integer_dense right):

        cdef Integer h
        cdef mod_int *moduli
        cdef int i, n, k

        h = left.height() * right.height() * left.ncols()
        verbose('multiplying matrices of height %s and %s'%(left.height(),right.height()))
        mm = MultiModularBasis(h)
        res = left._reduce(mm)
        res_right = right._reduce(mm)
        k = len(mm)
        for i in range(k):  # yes, I could do this with zip, but to conserve memory...
            t = cputime()
            res[i] *= res_right[i]
            verbose('multiplied matrices modulo a prime (%s/%s)'%(i+1,k), t)
        result = left.new_matrix(left.nrows(), right.ncols())
        _lift_crt(result, res, mm)  # changes result
        return result

    def _mod_int(self, modulus):
        return self._mod_int_c(modulus)

    cdef _mod_int_c(self, mod_int p):
        cdef Py_ssize_t i, j
        cdef Matrix_modn_dense res
        cdef mpz_t* self_row
        cdef mod_int* res_row
        res = Matrix_modn_dense.__new__(Matrix_modn_dense, matrix_space.MatrixSpace(IntegerModRing(p), self._nrows, self._ncols, sparse=False), None, None, None)
        for i from 0 <= i < self._nrows:
            self_row = self._matrix[i]
            res_row = res._matrix[i]
            for j from 0 <= j < self._ncols:
                res_row[j] = mpz_fdiv_ui(self_row[j], p)
        return res

    def _reduce(self, moduli):

        if isinstance(moduli, (int, long, Integer)):
            return self._mod_int(moduli)
        elif isinstance(moduli, list):
            moduli = MultiModularBasis(moduli)

        cdef MultiModularBasis mm
        mm = moduli

        res = [Matrix_modn_dense.__new__(Matrix_modn_dense,
                                         matrix_space.MatrixSpace(IntegerModRing(p), self._nrows, self._ncols, sparse=False),
                                         None, None, None) for p in mm]

        cdef size_t i, k, n
        cdef Py_ssize_t nr, nc

        n = len(mm)
        nr = self._nrows
        nc = self._ncols

        cdef mod_int **row_list
        row_list = <mod_int**>sage_malloc(sizeof(mod_int*) * n)
        if row_list == NULL:
            raise MemoryError, "out of memory allocating multi-modular coefficent list"

        cdef PyObject** res_seq
        res_seq = FAST_SEQ_UNSAFE(res)

        _sig_on
        for i from 0 <= i < nr:
            for k from 0 <= k < n:
                row_list[k] = (<Matrix_modn_dense>res_seq[k])._matrix[i]
            mm.mpz_reduce_vec(self._matrix[i], row_list, nc)
        _sig_off

        sage_free(row_list)
        return res

    def _echelon_in_place_classical(self):
        cdef Matrix_integer_dense E
        E = self.echelon_form()

        cdef int i
        for i from 0 <= i < self._ncols * self._nrows:
            mpz_set(self._entries[i], E._entries[i])

        self.clear_cache()

    def _echelon_strassen(self):
        raise NotImplementedError

    def echelon_form(self, algorithm="default", cutoff=0, include_zero_rows=True):
        r"""
        Return the echelon form of this matrix over the integers.

        INPUT:
            algorithm, cutoff -- ignored currently
            include_zero_rows -- (default: True) if False, don't include zero rows.

        EXAMPLES:
            sage: A = MatrixSpace(ZZ,2)([1,2,3,4])
            sage: A.echelon_form()
            [1 0]
            [0 2]

            sage: A = MatrixSpace(ZZ,5)(range(25))
            sage: A.echelon_form()
            [  5   0  -5 -10 -15]
            [  0   1   2   3   4]
            [  0   0   0   0   0]
            [  0   0   0   0   0]
            [  0   0   0   0   0]
        """
        x = self.fetch('echelon_form')
        if not x is None:
            return x

        if self._nrows == 0 or self._ncols == 0:
            self.cache('echelon_form', self)
            self.cache('pivots', [])
            self.cache('rank', 0)
            return self

        cdef Py_ssize_t nr, nc, n, i
        nr = self._nrows
        nc = self._ncols

        # The following complicated sequence of column reversals
        # and transposes is needed since PARI's Hermite Normal Form
        # does column operations instead of row operations.
        n = nc
        r = []
        for i from 0 <= i < n:
            r.append(n-i)
        v = self._pari_()
        v = v.vecextract(r) # this reverses the order of columns
        v = v.mattranspose()
        w = v.mathnf(1)

        cdef Matrix_integer_dense H_m
        H = convert_parimatrix(w[0])
        if nc == 1:
            H = [H]

        # We do a 'fast' change of the above into a list of ints,
        # since we know all entries are ints:
        num_missing_rows = (nr*nc - len(H)) / nc
        rank = nr - num_missing_rows

        if include_zero_rows:
            H = H + ['0']*(num_missing_rows*nc)
            H_m = self.new_matrix(nrows=nr, ncols=nc, entries=H, coerce=True)
        else:
            H_m = self.new_matrix(nrows=rank, ncols=nc, entries=H, coerce=True)

        H_m.set_immutable()
        H_m.cache('rank', rank)
        self.cache('rank',rank)
        H_m.cache('echelon_form',H_m)
        self.cache('echelon_form',H_m)
        return H_m

    def pivots(self):
        """
        Return the pivot column positions of this matrix as a list of Python integers.

        This returns a list, of the position of the first nonzero entry in each row
        of the echelon form.

        OUTPUT:
             list -- a list of Python ints

        EXAMPLES:
            sage: n = 3; A = matrix(ZZ,n,range(n^2)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.pivots()
            [0, 1]
            sage: A.echelon_form()
            [ 3  0 -3]
            [ 0  1  2]
            [ 0  0  0]
        """
        p = self.fetch('pivots')
        if not p is None: return p

        cdef Matrix_integer_dense E
        E = self.echelon_form()

        # Now we determine the pivots from the matrix E as quickly as we can.
        # For each row, we find the first nonzero position in that row -- it is the pivot.
        cdef Py_ssize_t i, j, k
        cdef mpz_t *row
        p = []
        k = 0
        for i from 0 <= i < E._nrows:
            row = E._matrix[i]   # pointer to ith row
            for j from k <= j < E._ncols:
                if mpz_cmp_si(row[j], 0) != 0:  # nonzero position
                    p.append(j)
                    k = j+1  # so start at next position next time
                    break
        self.cache('pivots', p)
        return p

    #### Elementary divisors

    def elementary_divisors(self, algorithm='linbox'):
        """
        Return the elementary divisors of self, in order.

        IMPLEMENTATION: Uses linbox.

        WARNING: This is MUCH faster than the smith_form function.

        The elementary divisors are the invariants of the finite
        abelian group that is the cokernel of this matrix.  They are
        ordered in reverse by divisibility.

        INPUT:
            self -- matrix
            algorithm -- 'linbox' or 'pari'

        OUTPUT:
            list of int's

        EXAMPLES:
            sage: matrix(3, range(9)).elementary_divisors()
            [1, 3, 0]
            sage: matrix(3, range(9)).elementary_divisors(algorithm='pari')
            [1, 3, 0]
            sage: C = MatrixSpace(ZZ,4)([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
            sage: C.elementary_divisors()
            [1, 1, 1, 687]

        SEE ALSO: smith_form
        """
        d = self.fetch('elementary_divisors')
        if not d is None:
            return d
        if self._nrows == 0 or self._ncols == 0:
            d = []
        else:
            if algorithm == 'linbox':
                d = self._elementary_divisors_linbox()
            elif algorithm == 'pari':
                d = self._pari_().matsnf(0).python()
                i = d.count(0)
                if i > 0:
                    d = list(reversed(d[i:])) + [d[0]]*i
            else:
                raise ValueError, "algorithm (='%s') unknown"%algorithm
        self.cache('elementary_divisors', d)
        return d

    def _elementary_divisors_linbox(self):
        self._init_linbox()
        _sig_on
        d = linbox.smithform()
        _sig_off
        return d

    def smith_form(self):
        """
        Returns matrices S, U, and V such that S = U*self*V, and S
        is in Smith normal form.  Thus S is diagonal with diagonal
        entries the ordered elementary divisors of S.

        WARNING: The elementary_divisors function, which returns
        the diagonal entries of S, is VASTLY faster than this
        function.

        The elementary divisors are the invariants of the finite
        abelian group that is the cokernel of this matrix.  They are
        ordered in reverse by divisibility.

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
            sage: D, U, V = A.smith_form()
            sage: D
            [0 0 0]
            [0 3 0]
            [0 0 1]
            sage: U
            [-1  2 -1]
            [ 0 -1  1]
            [ 0  1  0]
            sage: V
            [ 1  4 -1]
            [-2 -3  1]
            [ 1  0  0]
            sage: U*A*V
            [0 0 0]
            [0 3 0]
            [0 0 1]

        It also makes sense for nonsquare matrices:

            sage: A = Matrix(ZZ,3,2,range(6))
            sage: D, U, V = A.smith_form()
            sage: D
            [0 0]
            [2 0]
            [0 1]
            sage: U
            [-1  2 -1]
            [ 0 -1  1]
            [ 0  1  0]
            sage: V
            [ 3 -1]
            [-2  1]
            sage: U * A * V
            [0 0]
            [2 0]
            [0 1]

        SEE ALSO: elementary_divisors
        """
        v = self._pari_().matsnf(1).python()
        D = self.matrix_space()(v[2])
        U = self.matrix_space(ncols = self._nrows)(v[0])
        V = self.matrix_space(nrows = self._ncols)(v[1])
        return D, U, V

    def frobenius(self,flag=0, var='x'):
        """
        Return the Frobenius form (rational canonical form) of this matrix.

        INPUT:
            flag --an integer:
                0 -- (default) return the Frobenius form of this matrix.
                1 -- return only the elementary divisor polynomials, as
                     polynomials in var.
                2 -- return a two-components vector [F,B] where F is the
                     Frobenius form and B is the basis change so that $M=B^{-1}FB$.
            var -- a string (default: 'x')

        INPUT:
           flag -- 0 (default), 1 or 2 as described above

        ALGORITHM: uses pari's matfrobenius()

        EXAMPLE:
           sage: A = MatrixSpace(ZZ, 3)(range(9))
           sage: A.frobenius(0)
           [ 0  0  0]
           [ 1  0 18]
           [ 0  1 12]
           sage: A.frobenius(1)
           [x^3 - 12*x^2 - 18*x]
           sage: A.frobenius(1, var='y')
           [y^3 - 12*y^2 - 18*y]
           sage: A.frobenius(2)
           ([ 0  0  0]
           [ 1  0 18]
           [ 0  1 12],
           [    -1      2     -1]
           [     0  23/15 -14/15]
           [     0  -2/15   1/15])

        AUTHOR:
           -- 2006-04-02: Martin Albrecht

        TODO:
           -- move this to work for more general matrices than just over Z.
              This will require fixing how PARI polynomials are coerced
              to SAGE polynomials.
        """
        if not self.is_square():
            raise ArithmeticError, "frobenius matrix of non-square matrix not defined."

        v = self._pari_().matfrobenius(flag)
        if flag==0:
            return self.matrix_space()(v.python())
        elif flag==1:
            r = PolynomialRing(self.base_ring(), names=var)
            retr = []
            for f in v:
                retr.append(eval(str(f).replace("^","**"), {'x':r.gen()}, r.gens_dict()))
            return retr
        elif flag==2:
            F = matrix_space.MatrixSpace(QQ, self.nrows())(v[0].python())
            B = matrix_space.MatrixSpace(QQ, self.nrows())(v[1].python())
            return F, B

    def kernel(self, LLL=False):
        r"""
        Return the kernel of this matrix, as a module over the integers.

        INPUT:
           LLL -- bool (default: False); if True the basis is an LLL
                  reduced basis; otherwise, it is an echelon basis.

        EXAMPLES:
            sage: M = MatrixSpace(ZZ,4,2)(range(8))
            sage: M.kernel()
            Free module of degree 4 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -3  2]
            [ 0  1 -2  1]
        """
        if self._nrows == 0:    # from a 0 space
            M = sage.modules.free_module.FreeModule(ZZ, self._nrows)
            return M.zero_subspace()

        elif self._ncols == 0:  # to a 0 space
            return sage.modules.free_module.FreeModule(ZZ, self._nrows)

        A = self._pari_().mattranspose()
        B = A.matkerint()
        n = self._nrows
        M = sage.modules.free_module.FreeModule(ZZ, n)

        if B.ncols() == 0:
            return M.zero_submodule()

        # Now B is a basis for the LLL-reduced integer kernel as a
        # PARI object.  The basis vectors or B[0], ..., B[n-1],
        # where n is the dimension of the kernel.
        X = []
        for b in B:
            tmp = []
            for x in b:
                tmp.append(ZZ(x))
            X.append(M(tmp))

        if LLL:
            return M.span_of_basis(X)
        else:
            return M.span(X)

    def _adjoint(self):
        """
        Return the adjoint of this matrix.

        Assumes self is a square matrix (checked in adjoint).
        """
        return self.parent()(self._pari_().matadjoint().python())

    def _ntl_(self):
        r"""
        ntl.mat_ZZ representation of self.

        \note{NTL only knows dense matrices, so if you provide a
        sparse matrix NTL will allocate memory for every zero entry.}
        """
        return mat_ZZ(self._nrows,self._ncols, self.list())


    ####################################################################################
    # LLL
    ####################################################################################
    def lllgram(self):
        """
        LLL reduction of the lattice whose gram matrix is self.

        INPUT:
            M -- gram matrix of a definite quadratic form

        OUTPUT:
            U -- unimodular transformation matrix such that

                U.transpose() * M * U

            is LLL-reduced

        ALGORITHM:
            Use PARI

        EXAMPLES:
            sage: M = Matrix(ZZ, 2, 2, [5,3,3,2]) ; M
            [5 3]
            [3 2]
            sage: U = M.lllgram(); U
            [-1  1]
            [ 1 -2]
            sage: U.transpose() * M * U
            [1 0]
            [0 1]

        Semidefinite and indefinite forms raise a ValueError:

            sage: Matrix(ZZ,2,2,[2,6,6,3]).lllgram()
            Traceback (most recent call last):
            ...
            ValueError: not a definite matrix
            sage: Matrix(ZZ,2,2,[1,0,0,-1]).lllgram()
            Traceback (most recent call last):
            ...
            ValueError: not a definite matrix

        BUGS:
            should work for semidefinite forms (PARI is ok)
        """
        if self._nrows != self._ncols:
            raise ArithmeticError, "matrix must be square"

        n = self.nrows()
        # maybe should be /unimodular/ matrices ?
        P = self._pari_()
        try:
            U = P.lllgramint()
        except (RuntimeError, ArithmeticError), msg:
            raise ValueError, "not a definite matrix"
        MS = matrix_space.MatrixSpace(ZZ,n)
        U = MS(U.python())
        # Fix last column so that det = +1
        if U.det() == -1:
            for i in range(n):
                U[i,n-1] = - U[i,n-1]
        return U

    def prod_of_row_sums(self, cols):
        cdef Py_ssize_t c, row
        cdef mpz_t s, pr
        mpz_init(s)
        mpz_init(pr)

        mpz_set_si(pr, 1)
        for row from 0 <= row < self._nrows:
            tmp = []
            mpz_set_si(s, 0)
            for c in cols:
                if c<0 or c >= self._ncols:
                    raise IndexError, "matrix column index out of range"
                mpz_add(s, s, self._matrix[row][c])
            mpz_mul(pr, pr, s)
        cdef Integer z
        z = PY_NEW(Integer)
        mpz_set(z.value, pr)
        mpz_clear(s)
        mpz_clear(pr)
        return z

    def _linbox_sparse(self):
        cdef Py_ssize_t i, j
        v = ['%s %s +'%(self._nrows, self._ncols)]
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                if mpz_cmp_si(self._matrix[i][j], 0):
                    v.append('%s %s %s'%(i+1,j+1,self.get_unsafe(i,j)))
        v.append('0 0 0\n')
        return '\n'.join(v)

    def _linbox_dense(self):
        cdef Py_ssize_t i, j
        v = ['%s %s x'%(self._nrows, self._ncols)]
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                v.append(str(self.get_unsafe(i,j)))
        return ' '.join(v)

    def rational_reconstruction(self, N):
        """
        Use rational reconstruction to lift self to a matrix over the
        rational numbers (if possible), where we view self as a matrix
        modulo N.
        """
        import misc
        return misc.matrix_integer_dense_rational_reconstruction(self, N)

    def randomize(self, density=1, x=None, y=None):
        """
        Randomize density proportion of the entries of this matrix,
        leaving the rest unchanged.

        The randomized entries of this matrix to be between x and y
        and have density 1.
        """
        self.check_mutability()
        self.clear_cache()

        cdef int _min, _max
        if y is None:
            if x is None:
                min = -2
                max = 3
            else:
                min = 0
                max = x
        else:
            min = x
            max = y

        density = float(density)

        cdef int min_is_zero
        min_is_nonzero = (min != 0)

        cdef Integer n_max, n_min, n_width
        n_max = Integer(max)
        n_min = Integer(min)
        n_width = n_max - n_min

        cdef Py_ssize_t i, j, k, nc, num_per_row
        global state

        _sig_on
        if density == 1:
            for i from 0 <= i < self._nrows*self._ncols:
                mpz_urandomm(self._entries[i], state, n_width.value)
                if min_is_nonzero:
                    mpz_add(self._entries[i], self._entries[i], n_min.value)
        else:
            nc = self._ncols
            num_per_row = int(density * nc)
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < num_per_row:
                    k = random()%nc
                    mpz_urandomm(self._matrix[i][k], state, n_width.value)
                    if min_is_nonzero:
                        mpz_add(self._matrix[i][k], self._matrix[i][j], n_min.value)
        _sig_off

    #### Rank

    def rank(self):
        r = self.fetch('rank')
        if not r is None: return r

        # Can very quickly detect full rank by working modulo p.
        r = self._rank_modp()
        if r == self._nrows or r == self._ncols:
            self.cache('rank', r)
            return r
        # Detecting full rank didn't work -- use Linbox's general algorithm.
        r = self._rank_linbox()
        self.cache('rank', r)
        return r

    def _rank_linbox(self):
        """
        Compute the rank of this matrix using Linbox.
        """
        self._init_linbox()
        cdef unsigned long r
        _sig_on
        r = linbox.rank()
        _sig_off
        return Integer(r)

    def _rank_modp(self, p=46337):
        A = self._mod_int_c(p)
        return A.rank()

    #### Determinante

    def determinant(self):
        """
        Return the determinant of this matrix.

        ALGORITHM: Uses linbox.

        EXAMPLES:

        """
        d = self.fetch('det')
        if not d is None:
            return d
        d = self._det_linbox()
        self.cache('det', d)
        return d

    def _det_linbox(self):
        """
        Compute the determinant of this matrix using Linbox.
        """
        self._init_linbox()
        _sig_on
        d = linbox.det()
        _sig_off
        return Integer(d)

    #### Rational kernel, via IML
    def _rational_kernel_iml(self):
        """
        IML: Return the rational kernel of this matrix (acting from
        the left), considered as a matrix over QQ.  I.e., returns a
        matrix K such that self*K = 0, and the number of columns of K
        equals the nullity of self.

        AUTHOR:
            -- William Stein
        """
        cdef long dim
        cdef mpz_t *mp_N
        time = verbose('computing nullspace of %s x %s matrix using IML'%(self._nrows, self._ncols))
        _sig_on
        dim = nullspaceMP (self._nrows, self._ncols, self._entries, &mp_N)
        _sig_off

        # Now read the answer as a matrix.
        cdef Matrix_integer_dense M
        P = self.matrix_space(self._ncols, dim)
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)
        for i from 0 <= i < dim*self._ncols:
            mpz_init_set(M._entries[i], mp_N[i])
            mpz_clear(mp_N[i])
        free(mp_N)

        verbose("finished computing nullspace", time)
        M._initialized = True
        return M

    def _invert_iml(self, use_nullspace=False):
        """
        Invert this matrix using IML.  The output matrix is an integer
        matrix and a denominator.

        INPUT:
           self -- an invertible matrix
           use_nullspace -- (default: False): whether to use nullspace
                     algorithm, which is slower, but doesn't require
                     checking that the matrix is invertible as a
                     precondition.

        OUTPUT: A, d such that A*self = d
           A -- a matrix over ZZ
           d -- an integer

        ALGORITHM: Uses IML's p-adic nullspace function.

        EXAMPLES:
            sage: a = matrix(ZZ,3,[1,2,5, 3,7,8, 2,2,1])
            sage: b, d = a._invert_iml(); b,d
            ([  9  -8  19]
            [-13   9  -7]
            [  8  -2  -1], 23)
            sage: a*b
            [23  0  0]
            [ 0 23  0]
            [ 0  0 23]

        AUTHOR:
            -- William Stein
        """
        if self._nrows != self._ncols:
            raise TypeError, "self must be a square matrix."

        P = self.parent()
        time = verbose('computing inverse of %s x %s matrix using IML'%(self._nrows, self._ncols))
        if use_nullspace:
            A = self.augment(P.identity_matrix())
            K = A._rational_kernel_iml()
            d = -K[self._nrows,0]
            if K.ncols() != self._ncols or d == 0:
                raise ZeroDivisionError, "input matrix must be nonsingular"
            B = K[:self._nrows]
            verbose("finished computing inverse using IML", time)
            return B, d
        else:
            if self.rank() != self._nrows:
                raise ZeroDivisionError, "input matrix must be nonsingular"
            return self._solve_iml(P.identity_matrix(), right=True)


    def _solve_iml(self, Matrix_integer_dense B, right=True):
        """
        Let A equal self. Given B return an integer matrix C and an
        integer d such that self C*A == d*B if right is False or
        A*C == d*B if right is True.

        OUTPUT:
            C -- integer matrix
            d -- integer denominator

        EXAMPLES:
            sage: A = matrix(ZZ,4,4,[0, 1, -2, -1, -1, 1, 0, 2, 2, 2, 2, -1, 0, 2, 2, 1])
            sage: B = matrix(ZZ,3,4, [-1, 1, 1, 0, 2, 0, -2, -1, 0, -2, -2, -2])
            sage: C,d = A._solve_iml(B,right=False); C
            [  6 -18 -15  27]
            [  0  24  24 -36]
            [  4 -12  -6  -2]

            sage: d
            12

            sage: C*A == d*B
            True

            sage: A = matrix(ZZ,4,4,[0, 1, -2, -1, -1, 1, 0, 2, 2, 2, 2, -1, 0, 2, 2, 1])
            sage: B = matrix(ZZ,4,3, [-1, 1, 1, 0, 2, 0, -2, -1, 0, -2, -2, -2])
            sage: C,d = A._solve_iml(B)
            sage: C
            [ 12  40  28]
            [-12  -4  -4]
            [ -6 -25 -16]
            [ 12  34  16]

            sage: d
            12

            sage: A*C == d*B
            True

        ALGORITHM: Uses IML.

        AUTHOR:
            -- Martin Albrecht
        """
        cdef int i
        cdef mpz_t *mp_N, mp_D
        cdef Matrix_integer_dense M
        cdef Integer D

        if self._nrows != self._ncols:
            raise ArithmeticError, "self must be square"

        if self.nrows() == 1:
            return B, self[0,0]

        mpz_init(mp_D)


        if right:
            if self._ncols != B._nrows:
                raise ArithmeticError, "B's number of rows must match self's number of columns"

            n = self._ncols
            m = B._ncols

            if m == 0 or n == 0:
                return self.new_matrix(nrows = n, ncols = m), Integer(1)

            mp_N = <mpz_t *> sage_malloc( n * m * sizeof(mpz_t) )
            for i from 0 <= i < n * m:
                mpz_init( mp_N[i] )

            nonsingSolvLlhsMM(RightSolu, n, m, self._entries, B._entries, mp_N, mp_D)

            P = self.matrix_space(n, m)

        else: # left
            if self._nrows != B._ncols:
                raise ArithmeticError, "B's number of columns must match self's number of rows"

            n = self._ncols
            m = B._nrows

            if m == 0 or n == 0:
                return self.new_matrix(nrows = m, ncols = n), Integer(1)

            mp_N = <mpz_t *> sage_malloc( n * m * sizeof(mpz_t) )
            for i from 0 <= i < n * m:
                mpz_init( mp_N[i] )

            nonsingSolvLlhsMM(LeftSolu, n, m, self._entries, B._entries, mp_N, mp_D)

            P = self.matrix_space(m, n)

        M = Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)
        for i from 0 <= i < n*m:
            mpz_init_set(M._entries[i], mp_N[i])
            mpz_clear(mp_N[i])
        M._initialized = True

        D = PY_NEW(Integer)
        mpz_set(D.value, mp_D)
        mpz_clear(mp_D)

        return M,D

    def _rational_echelon_via_solve(self):
        r"""
        Computes information that gives the reduced row echelon form
        (over QQ!) of a matrix with integer entries.

        INPUT:
            self -- a matrix over the integers.

        OUTPUT:
            pivots -- ordered list of integers that give the pivot column positions
            nonpivots -- ordered list of the nonpivot column positions
            X -- matrix with integer entries
            d -- integer

        If you put standard basis vectors in order at the pivot
        columns, and put the matrix (1/d)*X everywhere else, then
        you get the reduced row echelon form of self, without zero
        rows at the bottom.

        NOTE: IML is the actual underlying $p$-adic solver that we use.

        EXAMPLES:



        AUTHOR:
           -- William Stein

        ALGORITHM:  I came up with this algorithm from scratch.
        As far as I know it is new.  It's pretty simple, but it
        is ... (fast?!).

        Let A be the input matrix.

        \begin{enumerate}
            \item[1] Compute r = rank(A).

            \item[2] Compute the pivot columns of the transpose $A^t$ of
                $A$.  One way to do this is to choose a random prime
                $p$ and compute the row echelon form of $A^t$ modulo
                $p$ (if the number of pivot columns is not equal to
                $r$, pick another prime).

            \item[3] Let $B$ be the submatrix of $A$ consisting of the
                rows corresponding to the pivot columns found in
                the previous step.  Note that, aside from zero rows
                at the bottom, $B$ and $A$ have the same reduced
                row echelon form.

            \item[4] Compute the pivot columns of $B$, again possibly by
                choosing a random prime $p$ as in [2] and computing
                the Echelon form modulo $p$.  If the number of pivot
                columns is not $r$, repeat with a different prime.
                Note -- in this step it is possible that we mistakenly
                choose a bad prime $p$ such that there are the right
                number of pivot columns modulo $p$, but they are
                at the wrong positions -- e.g., imagine the
                augmented matrix $[pI|I]$ -- modulo $p$ we would
                miss all the pivot columns.   This is OK, since in
                the next step we would detect this, as the matrix
                we would obtain would not be in echelon form.

            \item[5] Let $C$ be the submatrix of $B$ of pivot columns.
                Let $D$ be the complementary submatrix of $B$ of
                all all non-pivot columns.  Use a $p$-adic solver
                to find the matrix $X$ and integer $d$ such
                that $C (1/d) X=D$.  I.e., solve a bunch of linear
                systems of the form $Cx = v$, where the columns of
                $X$ are the solutions.

            \item[6] Verify that we had chosen the correct pivot
                columns.  Inspect the matrix $X$ obtained in step 5.
                If when used to construct the echelon form of $B$, $X$
                indeed gives a matrix in reduced row echelon form,
                then we are done -- output the pivot columns, $X$, and
                $d$. To tell if $X$ is correct, do the following: For
                each column of $X$ (corresponding to non-pivot column
                $i$ of $B$), read up the column of $X$ until finding
                the first nonzero position $j$; then verify that $j$
                is strictly less than the number of pivot columns of
                $B$ that are strictly less than $i$.  Otherwise, we
                got the pivots of $B$ wrong -- try again starting at
                step 4, but with a different random prime.

        \end{enumerate}
        """
        from matrix_modn_dense import MAX_MODULUS
        A = self
        # Step 1: Compute the rank

        t = verbose('computing rank', level=2, caller_name='p-adic echelon')
        r = A.rank()
        verbose('done computing rank', level=2, t=t, caller_name='p-adic echelon')

        if r == self._nrows:
            # The input matrix already has full rank.
            B = A
        else:
            # Steps 2 and 3: Extract out a submatrix of full rank.
            i = 0
            while True:
                p = previous_prime(random() % (MAX_MODULUS-15000) + 10000)
                P = A._mod_int(p).transpose().pivots()
                if len(P) == r:
                    B = A.matrix_from_rows(P)
                    break
                else:
                    i += 1
                    if i == 50:
                        raise RuntimeError, "Bug in _rational_echelon_via_solve in finding linearly independent rows."

        _i = 0
        while True:
            _i += 1
            if _i == 50:
                raise RuntimeError, "Bug in _rational_echelon_via_solve -- pivot columns keep being wrong."

            # Step 4: Now we instead worry about computing the reduced row echelon form of B.
            i = 0
            while True:
                p = previous_prime(random() % (MAX_MODULUS-15000) + 10000)
                pivots = B._mod_int(p).pivots()
                if len(pivots) == r:
                    break
                else:
                    i += 1
                    if i == 50:
                        raise RuntimeError, "Bug in _rational_echelon_via_solve in finding pivot columns."

            # Step 5: Apply p-adic solver
            C = B.matrix_from_columns(pivots)
            pivots_ = set(pivots)
            non_pivots = [i for i in range(B.ncols()) if not i in pivots_]
            D = B.matrix_from_columns(non_pivots)
            t = verbose('calling IML solver', level=2, caller_name='p-adic echelon')
            X, d = C._solve_iml(D, right=True)
            t = verbose('finished IML solver', level=2, caller_name='p-adic echelon', t=t)

            # Step 6: Verify that we had chosen the correct pivot columns.
            pivots_are_right = True
            for z in range(X.ncols()):
                if not pivots_are_right:
                    break
                i = non_pivots[z]
                np = len([k for k in pivots if k < i])
                for j in reversed(range(X.nrows())):
                    if X[j,z] != 0:
                        if j < np:
                            break # we're good -- go on to next column of X
                        else:
                            pivots_are_right = False
                            break
            if pivots_are_right:
                break

        #end while


        # Finally, return the answer.
        return pivots, non_pivots, X, d

    def decomposition(self, **kwds):
        """
        Returns the decomposition of the free module on which this
        matrix A acts from the right (i.e., the action is x goes to x
        A), along with whether this matrix acts irreducibly on each
        factor.  The factors are guaranteed to be sorted in the same
        way as the corresponding factors of the characteristic
        polynomial, and are saturated as ZZ modules.

        INPUT:
            self -- a matrix over the integers
            **kwds -- these are passed onto to the decomposition over QQ command.

        EXAMPLES:
            sage: t = ModularSymbols(11,sign=1).hecke_matrix(2)
            sage: w = t.change_ring(ZZ)
            sage: w.list()
            [3, -1, 0, -2]
        """
        F = self.charpoly().factor()
        if len(F) == 1:
            V = self.base_ring()**self.nrows()
            return decomp_seq([(V, bool(F[0][1]==1))])

        A = self.change_ring(QQ)
        X = A.decomposition(**kwds)
        V = ZZ**self.nrows()
        if isinstance(X, tuple):
            D, E = X
            D = [(W.intersection(V), t) for W, t in D]
            E = [(W.intersection(V), t) for W, t in E]
            return D, E
        else:
            return [(W.intersection(V), t) for W, t in X]


###############################################################

###########################################
# Helper code for Echelon form algorithm.
###########################################
def _parimatrix_to_strlist(A):
    s = str(A)
    s = s.replace('Mat(','').replace(')','')
    s = s.replace(';',',').replace(' ','')
    s = s.replace(",", "','")
    s = s.replace("[", "['")
    s = s.replace("]", "']")
    return eval(s)

def _parimatrix_to_reversed_strlist(A):
    s = str(A)
    if s.find('Mat') != -1:
        return _parimatrix_to_strlist(A)
    s = s.replace('[','').replace(']','').replace(' ','')
    v = s.split(';')
    v.reverse()
    s = "['" + (','.join(v)) + "']"
    s = s.replace(",", "','")
    return eval(s)

def convert_parimatrix(z):
    n = z.ncols();
    r = []
    for i from 0 <= i < n:
        r.append(n-i)
    z = z.vecextract(r)
    z = z.mattranspose()
    n = z.ncols();
    r = []
    for i from 0 <= i < n:
        r.append(n-i)
    z = z.vecextract(r)
    return _parimatrix_to_strlist(z)



def _lift_crt(Matrix_integer_dense M, residues, moduli=None):

    cdef size_t n, i, j, k
    cdef Py_ssize_t nr, nc

    n = len(residues)
    nr = residues[0].nrows()
    nc = residues[0].ncols()

    if moduli is None:
        moduli = MultiModularBasis([m.base_ring().order() for m in residues])
    else:
        if len(residues) != len(moduli):
            raise IndexError, "Number of residues (%s) does not match number of moduli (%s)"%(len(residues), len(moduli))

    cdef MultiModularBasis mm
    mm = moduli

    for b in residues:
        if not PY_TYPE_CHECK(b, Matrix_modn_dense):
            raise TypeError, "Can only perform CRT on list of type Matrix_modn_dense."
    cdef PyObject** res
    res = FAST_SEQ_UNSAFE(residues)

    cdef mod_int **row_list
    row_list = <mod_int**>sage_malloc(sizeof(mod_int*) * n)
    if row_list == NULL:
        raise MemoryError, "out of memory allocating multi-modular coefficent list"

    _sig_on
    for i from 0 <= i < nr:
        for k from 0 <= k < n:
            row_list[k] = (<Matrix_modn_dense>res[k])._matrix[i]
        mm.mpz_crt_vec(M._matrix[i], row_list, nc)
    _sig_off

    sage_free(row_list)
    return M

##########################################################
# Setup the c-library and GMP random number generators.
# seed it when module is loaded.
from random import randrange
cdef extern from "stdlib.h":
    long random()
    void srandom(unsigned int seed)
k = randrange(0,Integer(2)**(32))
srandom(k)

cdef gmp_randstate_t state
gmp_randinit_mt(state)
gmp_randseed_ui(state,k)

#######################################################

# Conclusions:
#  On OS X Intel, at least:
#    - if log_2(height) >= 0.70 * nrows, use classical

def tune_multiplication(k, nmin=10, nmax=200, bitmin=2,bitmax=64):
    from constructor import random_matrix
    from sage.rings.integer_ring import ZZ
    for n in range(nmin,nmax,10):
        for i in range(bitmin, bitmax, 4):
            A = random_matrix(ZZ, n, n, x = 2**i)
            B = random_matrix(ZZ, n, n, x = 2**i)
            t0 = cputime()
            for j in range(k//n + 1):
                C = A._multiply_classical(B)
            t0 = cputime(t0)
            t1 = cputime()
            for j in range(k//n+1):
                C = A._multiply_multi_modular(B)
            t1 = cputime(t1)
            print n, i, float(i)/float(n)
            if t0 < t1:
                print 'classical'
            else:
                print 'multimod'


#########################################################
# Interface to the integer matrix library.
#########################################################
## cdef class Matrix_IML:
##     cdef mpz_t* matrix
##     cdef long nrows
##     cdef long ncols

##     def __init__(self):
##         self.matrix = NULL

##     cdef set(self, mpz_t* matrix, long nrows, long ncols):
##         self.nrows = nrows
##         self.ncols = ncols
##         self.matrix = matrix

##     def kernel(self, mpz_t** answer):
##         """
##         INPUT:
##             self  -- assumes nrows, ncols and matrix have been set
##         OUTPUT:
##         """
##         if self.matrix == NULL:
##             raise RuntimeError, "you must set self.matrix first."
##         raise NotImplementedError

##     def solve(self, v):
##         raise NotImplementedError




cdef decomp_seq(v):
    return Sequence(v, universe=tuple, check=False, cr=True)
