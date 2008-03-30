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

TESTS:
    sage: a = matrix(ZZ,2,range(4), sparse=False)
    sage: loads(dumps(a)) == a
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
from sage.structure.element import is_Element
from sage.structure.proof.proof import get_flag as get_proof_flag

include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/gmp.pxi"
include "../ext/random.pxi"

ctypedef unsigned int uint

from sage.ext.multi_modular import MultiModularBasis
from sage.ext.multi_modular cimport MultiModularBasis

from sage.rings.integer cimport Integer
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.rings.integer_ring import ZZ, IntegerRing_class
from sage.rings.integer_ring cimport IntegerRing_class
from sage.rings.integer_mod_ring import IntegerModRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing
from sage.structure.element cimport ModuleElement, RingElement, Element, Vector
from sage.structure.element import is_Vector
from sage.structure.sequence import Sequence

from matrix_modn_dense import Matrix_modn_dense
from matrix_modn_dense cimport Matrix_modn_dense

from matrix2 import decomp_seq

import sage.modules.free_module
import sage.modules.free_module_element

from matrix cimport Matrix

cimport sage.structure.element

import matrix_space

################
# Used for modular HNF
from sage.ext.arith cimport arith_int
cdef arith_int ai
ai = arith_int()
################

######### linbox interface ##########
from sage.libs.linbox.linbox cimport Linbox_integer_dense, Linbox_modn_dense
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
        r"""
        Returns a new copy of this matrix.

        EXAMPLES:
            sage: a = matrix(ZZ,1,3, [1,2,-3]); a
            [ 1  2 -3]
            sage: b = a.__copy__(); b
            [ 1  2 -3]
            sage: b is a
            False
            sage: b == a
            True
        """
        cdef Matrix_integer_dense A
        A = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent,
                                         0, 0, 0)
        cdef Py_ssize_t i
        _sig_on
        for i from 0 <= i < self._nrows * self._ncols:
            mpz_init_set(A._entries[i], self._entries[i])
        _sig_off
        A._initialized = True
        if self.subdivisions is not None:
            A.subdivide(*self.get_subdivisions())
        return A

    def __hash__(self):
        r"""
        Returns hash of self.

        self must be immutable.

        EXAMPLES:
            sage: a = Matrix(ZZ,2,[1,2,3,4])
            sage: hash(a)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

            sage: a.set_immutable()
            sage: hash(a)
            8
        """
        return self._hash()

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
        cdef bint is_list
        cdef Integer x

        if entries is None:
            x = ZZ(0)
            is_list = 0
        elif isinstance(entries, (int,long)) or is_Element(entries):
            try:
                x = ZZ(entries)
            except TypeError:
                self._initialized = False
                raise TypeError, "unable to coerce entry to an integer"
            is_list = 0
        else:
            entries = list(entries)
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
                    x = ZZ(entries[i])
                    # todo -- see integer.pyx and the TODO there; perhaps this could be
                    # sped up by creating a mpz_init_set_sage function.
                    mpz_init_set(self._entries[i], x.value)
                self._initialized = True
            else:
                for i from 0 <= i < self._nrows * self._ncols:
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
                mpz_set(self._entries[j], x.value)
                j = j + self._nrows + 1
            self._initialized = True

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        """
        Set position i,j of this matrix to x.

        (VERY UNSAFE -- value *must* be of type Integer).

        INPUT:
        i -- row
        j -- column
        value -- The value to set self[i,j] to.  value MUST be of type Integer

        EXAMPLES:
            sage: a = matrix(ZZ,2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a[0,0] = 10
            sage: a
            [10  1  2]
            [ 3  4  5]
        """
        mpz_set(self._matrix[i][j], (<Integer>value).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Returns (j, i) entry of self as a new Integer.

        WARNING: this is very unsafe; it assumes i and j are in the right range.

        EXAMPLES:
            sage: a = MatrixSpace(ZZ,3)(range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a[1,2]
            5
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

    ########################################################################
    # LEVEL 1 helpers:
    #   These function support the implementation of the level 1 functionality.
    ########################################################################
    cdef _zero_out_matrix(self):
        """
        Set this matrix to be the zero matrix.
        This is only for internal use.
        (Note: this matrix must NOT already have initialised entries.)
        """
        # TODO: This is about 6-10 slower than MAGMA doing what seems to be the same thing.
        # Moreover, NTL can also do this quickly.  Why?   I think both have specialized
        # small integer classes. (dmharvey: yes, NTL does not allocate any memory when
        # intialising a ZZ to zero.)
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

    def __nonzero__(self):
        r"""
        Tests whether self is the zero matrix.

        EXAMPLES:
            sage: a = MatrixSpace(ZZ, 2, 3)(range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.__nonzero__()
            True
            sage: (a - a).__nonzero__()
            False

            sage: a = MatrixSpace(ZZ, 0, 3)()
            sage: a.__nonzero__()
            False
            sage: a = MatrixSpace(ZZ, 3, 0)()
            sage: a.__nonzero__()
            False
            sage: a = MatrixSpace(ZZ, 0, 0)()
            sage: a.__nonzero__()
            False

        """
        cdef mpz_t *a, *b
        cdef Py_ssize_t i, j
        cdef int k
        for i from 0 <= i < self._nrows * self._ncols:
            if mpz_sgn(self._entries[i]):
                return True
        return False

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
            sage: b = MatrixSpace(ZZ,n,n)(range(1, n^2 + 1))
            sage: a._multiply_classical(b)
            [ 18  21  24]
            [ 54  66  78]
            [ 90 111 132]
        """
        if self._ncols != right._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of right."

        cdef Py_ssize_t i, j, k, l, nr, nc, snc
        cdef mpz_t *v
        cdef object parent

        nr = self._nrows
        nc = right._ncols
        snc = self._ncols

        if self._nrows == right._nrows:
            # self acts on the space of right
            parent = right.parent()
        if self._ncols == right._ncols:
            # right acts on the space of self
            parent = self.parent()
        else:
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
        r"""
        Compares self with right, examining entries in lexicographic
        (row major) ordering.

        EXAMPLES:
            sage: Matrix(ZZ, [[0, 10], [20, 30]]).__cmp__(Matrix(ZZ, [[0, 10], [20, 30]]))
            0
            sage: Matrix(ZZ, [[0, 10], [20, 30]]).__cmp__(Matrix(ZZ, [[0, 15], [20, 30]]))
            -1
            sage: Matrix(ZZ, [[5, 10], [20, 30]]).__cmp__(Matrix(ZZ, [[0, 15], [20, 30]]))
            1
            sage: Matrix(ZZ, [[5, 10], [20, 30]]).__cmp__(Matrix(ZZ, [[0, 10], [25, 30]]))
            1
        """
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

    cdef void reduce_entry_unsafe(self, Py_ssize_t i, Py_ssize_t j, Integer modulus):
        # Used for p-adic matrices.
        if mpz_cmp(self._matrix[i][j], modulus.value) >= 0 or mpz_cmp_ui(self._matrix[i][j], 0) < 0:
            mpz_mod(self._matrix[i][j], self._matrix[i][j], modulus.value)

    def _mod_int(self, modulus):
        cdef mod_int c = modulus
        if int(c) != modulus:
            raise OverflowError
        return self._mod_int_c(modulus)

    cdef _mod_int_c(self, mod_int p):
        cdef Py_ssize_t i, j
        cdef Matrix_modn_dense res
        cdef mpz_t* self_row
        cdef mod_int* res_row
        res = Matrix_modn_dense.__new__(Matrix_modn_dense,
                        matrix_space.MatrixSpace(IntegerModRing(p), self._nrows, self._ncols, sparse=False), None, None, None)
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

        _sig_on
        for i from 0 <= i < nr:
            for k from 0 <= k < n:
                row_list[k] = (<Matrix_modn_dense>res[k])._matrix[i]
            mm.mpz_reduce_vec(self._matrix[i], row_list, nc)
        _sig_off

        sage_free(row_list)
        return res

    def _linbox_modn_det(self, mod_int n):
        """
        INPUT:
            n -- a prime (at most 67108879)

        EXAMPLES:
            sage: a = matrix(ZZ, 3, [1,2,5,-7,8,10,192,5,18])
            sage: a.det()
            -3669
            sage: a._linbox_modn_det(5077)
            1408
            sage: a._linbox_modn_det(3)
            0
            sage: a._linbox_modn_det(2)
            1
            sage: a.det()%5077
            1408
            sage: a.det()%2
            1
            sage: a.det()%3
            0
        """
        d = self._linbox_modn(n).det()
        return IntegerModRing(n)(d)

    def _linbox_modn(self, mod_int n):
        """
        Return modn linbox object associated to this integer matrix.

        EXAMPLES:
            sage: a = matrix(ZZ, 3, [1,2,5,-7,8,10,192,5,18])
            sage: b = a._linbox_modn(19); b
            <sage.libs.linbox.linbox.Linbox_modn_dense object at ...>
            sage: b.charpoly()
            [2L, 10L, 11L, 1L]
        """
        if n > 67108879:   # doesn't work for bigger primes -- experimental observation
            raise NotImplementedError, "modulus to big"
        cdef mod_int** matrix = <mod_int**>sage_malloc(sizeof(mod_int*) * self._nrows)
        if matrix == NULL:
            raise MemoryError, "out of memory allocating multi-modular coefficent list"

        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            matrix[i] = <mod_int *>sage_malloc(sizeof(mod_int) * self._ncols)
            if matrix[i] == NULL:
                raise MemoryError, "out of memory allocating multi-modular coefficent list"
            for j from 0 <= j < self._ncols:
                matrix[i][j] = mpz_fdiv_ui(self._matrix[i][j], n)

        cdef Linbox_modn_dense L = Linbox_modn_dense()
        L.set(n, matrix, self._nrows, self._ncols)
        return L

    def _echelon_in_place_classical(self):
        cdef Matrix_integer_dense E
        E = self.echelon_form()

        cdef int i
        for i from 0 <= i < self._ncols * self._nrows:
            mpz_set(self._entries[i], E._entries[i])

        self.clear_cache()

    def _echelon_strassen(self):
        raise NotImplementedError

    def hermite_form(self, *args, **kwds):
        r"""
        Return the Hermite normal form of self.

        This is a synonym for \code{self.echelon_form(...)}.  See the documentation
        for \code{self.echelon_form} for more details.

        EXAMPLES:
            sage: A = matrix(ZZ, 3, 5, [-1, -1, -2, 2, -2, -4, -19, -17, 1, 2, -3, 1, 1, -4, 1])
            sage: E, U = A.hermite_form(transformation=True)
            sage: E
            [   1    0   52 -133  109]
            [   0    1   19  -47   38]
            [   0    0   69 -178  145]
            sage: U
            [-46   3  11]
            [-16   1   4]
            [-61   4  15]
            sage: U*A
            [   1    0   52 -133  109]
            [   0    1   19  -47   38]
            [   0    0   69 -178  145]
            sage: A.hermite_form()
            [   1    0   52 -133  109]
            [   0    1   19  -47   38]
            [   0    0   69 -178  145]

        TESTS:
        This example illustrated trac 2398.
            sage: a = matrix([(0, 0, 3), (0, -2, 2), (0, 1, 2), (0, -2, 5)])
            sage: a.hermite_form()
            [0 1 2]
            [0 0 3]
            [0 0 0]
            [0 0 0]
        """
        return self.echelon_form(*args, **kwds)

    def echelon_form(self, algorithm="padic", proof=None, include_zero_rows=True,
                     transformation=False, D=None):
        r"""
        Return the echelon form of this matrix over the integers also
        known as the hermit normal form (HNF).

        INPUT:
            algorithm --
                  'padic' -- (default) a new fast p-adic modular algorithm,
                  'pari' -- use PARI with flag 1
                  'ntl' -- use NTL (only works for square matrices of full rank!)
            proof -- (default: True); if proof=False certain
                   determinants are computed using a randomized hybrid
                   p-adic multimodular strategy until it stabilizes
                   twice (instead of up to the Hadamard bound).  It is
                   *incredibly* unlikely that one would ever get an
                   incorrect result with proof=False.
            include_zero_rows -- (default: True)
                                 if False, don't include zero rows
            transformation -- if given, also compute transformation matrix; only
                              valid for padic algorithm
            D -- (default: None)  if given and the algorithm is 'ntl', then D
                   must be a multiple of the determinant and this function
                   will use that fact.

        OUTPUT:
            matrix -- the Hermite normal form (=echelon form over ZZ) of self.

        NOTE: The result is \emph{not} cached.

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

        TESTS:
        Make sure the zero matrices are handled correctly:
            sage: m = matrix(ZZ,3,3,[0]*9)
            sage: m.echelon_form()
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: m = matrix(ZZ,3,1,[0]*3)
            sage: m.echelon_form()
            [0]
            [0]
            [0]
            sage: m = matrix(ZZ,1,3,[0]*3)
            sage: m.echelon_form()
            [0 0 0]

        The ultimate border case!
            sage: m = matrix(ZZ,0,0,[])
            sage: m.echelon_form()
            []

        NOTE: If 'ntl' is chosen for a non square matrix this function raises
        a ValueError.

        Special cases: 0 or 1 rows:
            sage: a = matrix(ZZ, 1,2,[0,-1])
            sage: a.hermite_form()
            [0 1]
            sage: a.pivots()
            [1]
            sage: a = matrix(ZZ, 1,2,[0,0])
            sage: a.hermite_form()
            [0 0]
            sage: a.pivots()
            []
            sage: a = matrix(ZZ,1,3); a
            [0 0 0]
            sage: a.echelon_form(include_zero_rows=False)
            []
            sage: a.echelon_form(include_zero_rows=True)
            [0 0 0]
        """
        if self._nrows == 0 or self._ncols == 0:
            self.cache('pivots', [])
            self.cache('rank', 0)
            if transformation:
                return self, self
            return self

        cdef Py_ssize_t nr, nc, n, i
        cdef int nr = self._nrows
        cdef int nc = self._ncols
        cdef int i, j

        cdef Matrix_integer_dense H_m

        proof = get_proof_flag(proof, "linear_algebra")
        pivots = None

        if algorithm == "padic":
            import matrix_integer_dense_hnf
            if transformation:
                if not include_zero_rows:
                    raise ValueError, "if you get the transformation matrix you must include zero rows"
                H_m, U, pivots = matrix_integer_dense_hnf.hnf_with_transformation(self, proof=proof)
            else:
                H_m, pivots = matrix_integer_dense_hnf.hnf(self,
                                   include_zero_rows=include_zero_rows, proof=proof)
            self.cache('pivots', pivots)
            self.cache('rank', len(pivots))


        elif algorithm == 'pari':
            if transformation:
                raise ValueError, "transformation matrix only available with p-adic algorithm"
            # The following complicated sequence of column reversals
            # and transposes is needed since PARI's Hermite Normal Form
            # does column operations instead of row operations.
            tm = verbose("pari hermite form")
            n = nc
            r = []
            for i from 0 <= i < n:
                r.append(n-i)
            v = self._pari_()
            v = v.vecextract(r) # this reverses the order of columns
            v = v.mattranspose()
            w = v.mathnf(1)
            H = _convert_parimatrix(w[0])
            # if H=[] may occur if we start with a column of zeroes
            if nc == 1 and H!=[]:
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
            verbose("finished pari hermite form",tm)

        elif algorithm == 'ntl':
            if transformation:
                raise ValueError, "transformation matrix only available with p-adic algorithm"

            if nr != nc:
                raise ValueError, "ntl only computes HNF for square matrices of full rank."

            import sage.libs.ntl.ntl_mat_ZZ
            v =  sage.libs.ntl.ntl_mat_ZZ.ntl_mat_ZZ(self._nrows,self._ncols)
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    v[i,j] = self.get_unsafe(nr-i-1,nc-j-1)

            try:
                w = v.HNF(D=D)
            except RuntimeError: # HNF may fail if a nxm matrix has rank < m
                raise ValueError, "ntl only computes HNF for square matrices of full rank."

            rank = w.nrows()

            if include_zero_rows:
                H_m = self.new_matrix()
            else:
                H_m = self.new_matrix(nrows=w.nrows())

            nr = w.nrows()
            nc = w.ncols()

            for i from 0 <= i < w.nrows():
                for j from 0 <= j < w.ncols():
                    H_m[i,j] = w[nr-i-1,nc-j-1]
            H_m.set_immutable()

        else:
            raise TypeError, "algorithm '%s' not understood"%(algorithm)

        H_m.set_immutable()
        if pivots is None:
            from matrix_integer_dense_hnf import pivots_of_hnf_matrix
            pivots = pivots_of_hnf_matrix(H_m)
            rank = len(pivots)

        H_m.cache('pivots', pivots)
        self.cache('pivots', pivots)

        H_m.cache('rank', rank)
        self.cache('rank',rank)

        H_m.cache('in_echelon_form', True)

        if transformation:
            return H_m, U
        else:
            return H_m

    def saturation(self, p=0, proof=None, max_dets=0):
        r"""
        Return a saturation matrix of self, which is a matrix whose
        rows span the saturation of the row span of self.  This
        is not unique.

        The saturation of a $\ZZ$ module $M$ embedded in $\ZZ^n$ is
        the a module $S$ that contains $M$ with finite index such that
        $\ZZ^n/S$ is torsion free.  This function takes the row span
        $M$ of self, and finds another matrix of full rank with row
        span the saturation of $M$.

        INPUT:
            p -- (default: 0); if nonzero given, saturate only at the
                 prime $p$, i.e., return a matrix whose row span is a
                 $\ZZ$-module $S$ that contains self and such that the
                 index of $S$ in its saturation is coprime to $p$.  If
                 $p$ is None, return full saturation of self.
            proof -- (default: use proof.linear_algebra()); if False, the
                 determinant calculations are done with proof=False.
            max_dets -- (default: 10); technical parameter -- max
                 number of determinant to compute when bounding prime
                 divisor of self in its saturation.

        OUTPUT:
            matrix -- a matrix over ZZ

        NOTE: The result is \emph{not} cached.

        ALGORITHM:
            1. Replace input by a matrix of full rank got from a
               subset of the rows.
            2. Divide out any common factors from rows.
            3. Check max_dets random dets of submatrices to see if their
               gcd (with p) is 1 -- if so matrix is saturated and we're done.
            4. Finally, use that if A is a matrix of full rank, then
                   $hnf(transpose(A))^{-1}*A$
               is a saturation of A.

        EXAMPLES:
            sage: A = matrix(ZZ, 3, 5, [-51, -1509, -71, -109, -593, -19, -341, 4, 86, 98, 0, -246, -11, 65, 217])
            sage: A.echelon_form()
            [      1       5    2262   20364   56576]
            [      0       6   35653  320873  891313]
            [      0       0   42993  386937 1074825]
            sage: S = A.saturation(); S
            [  -51 -1509   -71  -109  -593]
            [  -19  -341     4    86    98]
            [   35   994    43    51   347]

        Notice that the saturation spans a different module than A.
            sage: S.echelon_form()
            [ 1  2  0  8 32]
            [ 0  3  0 -2 -6]
            [ 0  0  1  9 25]
            sage: V = A.row_space(); W = S.row_space()
            sage: V.is_submodule(W)
            True
            sage: V.index_in(W)
            85986
            sage: V.index_in_saturation()
            85986

        We illustrate each option:
            sage: S = A.saturation(p=2)
            sage: S = A.saturation(proof=False)
            sage: S = A.saturation(max_dets=2)
        """
        proof = get_proof_flag(proof, "linear_algebra")
        import matrix_integer_dense_saturation
        return matrix_integer_dense_saturation.saturation(self, p=p, proof=proof, max_dets=max_dets)

    def index_in_saturation(self, proof=None):
        """
        Return the index of self in its saturation.

        INPUT:
            proof -- (default: use proof.linear_algebra()); if False, the
                     determinant calculations are done with proof=False.

        OUTPUT:
            positive integer -- the index of the row span of this matrix
                                in its saturation

        ALGORITHM:
            Use Hermite normal form twice to find an invertible matrix whose
            inverse transforms a matrix with the same row span as self
            to its saturation, then compute the determinant of that matrix.

        EXAMPLES:
            sage: A = matrix(ZZ, 2,3, [1..6]); A
            [1 2 3]
            [4 5 6]
            sage: A.index_in_saturation()
            3
            sage: A.saturation()
            [1 2 3]
            [1 1 1]
        """
        proof = get_proof_flag(proof, "linear_algebra")
        import matrix_integer_dense_saturation
        return matrix_integer_dense_saturation.index_in_saturation(self, proof=proof)

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

    def elementary_divisors(self, algorithm='pari'):
        """
        Return the elementary divisors of self, in order.

        IMPLEMENTATION: Uses linbox, except sometimes linbox doesn't
        work (errors about pre-conditioning), in which case PARI is
        used.

        WARNING: This is MUCH faster than the smith_form function.

        The elementary divisors are the invariants of the finite
        abelian group that is the cokernel of this matrix.  They are
        ordered in reverse by divisibility.

        INPUT:
            self -- matrix
            algorithm -- (default: 'pari')
                 'pari': works robustless, but is slower.
                 'linbox' -- use linbox (currently off, broken)

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

            sage: M = matrix(ZZ, 3, [1,5,7, 3,6,9, 0,1,2])
            sage: M.elementary_divisors()
            [1, 1, 6]

        This returns a copy, which is safe to change:
            sage: edivs = M.elementary_divisors()
            sage: edivs.pop()
            6
            sage: M.elementary_divisors()
            [1, 1, 6]

        SEE ALSO: smith_form
        """
        d = self.fetch('elementary_divisors')
        if not d is None:
            return d[:]
        if self._nrows == 0 or self._ncols == 0:
            d = []
        else:
            if algorithm == 'linbox':
                raise ValueError, "linbox too broken -- currently Linbox SNF is disabled."
                # This fails in linbox: a = matrix(ZZ,2,[1, 1, -1, 0, 0, 0, 0, 1])
                try:
                    d = self._elementary_divisors_linbox()
                except RuntimeError:
                    import sys
                    sys.stderr.write("DONT PANIC -- switching to using PARI (which will work fine)\n")
                    algorithm = 'pari'
            if algorithm == 'pari':
                d = self._pari_().matsnf(0).python()
                i = d.count(0)
                d.sort()
                if i > 0:
                    d = d[i:] + [d[0]]*i
            elif not (algorithm in ['pari', 'linbox']):
                raise ValueError, "algorithm (='%s') unknown"%algorithm
        self.cache('elementary_divisors', d)
        return d[:]

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

    def kernel_matrix(self, algorithm='padic', LLL=False, proof=None):
        """
        The options are exactly like self.kernel(...), but returns a
        matrix A whose rows form a basis for the left kernel, i.e.,
        so that A*self = 0.

        This is mainly useful to avoid all overhead associated with
        creating a free module.

        EXAMPLES:
            sage: A = matrix(ZZ, 3, 3, [1..9])
            sage: A.kernel_matrix()
            [-1  2 -1]

        Note that the basis matrix returned above is not in Hermite form.
            sage: A.kernel()
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 1 -2  1]

        We compute another kernel:
            sage: A = matrix(ZZ, 4, 2, [2, -1, 1, 1, -18, -1, -1, -5])
            sage: K = A.kernel_matrix(); K
            [-17 -20  -3   0]
            [  7   3   1  -1]

        K is a basis for the left kernel:
            sage: K*A
            [0 0]
            [0 0]

        We illustrate the LLL flag:
            sage: L = A.kernel_matrix(LLL=True); L
            [  7   3   1  -1]
            [  4 -11   0  -3]
            sage: K.hermite_form()
            [ 1 64  3 12]
            [ 0 89  4 17]
            sage: L.hermite_form()
            [ 1 64  3 12]
            [ 0 89  4 17]
        """
        if self._nrows == 0:    # from a 0 space
            return self.new_matrix(0, self.nrows())
        elif self._ncols == 0:  # to a 0 space
            # n x n identity matrix with n = self.nrows()
            import constructor
            return constructor.identity_matrix(self.nrows())

        proof = get_proof_flag(proof, "linear_algebra")

        if algorithm == 'pari':
            return self._kernel_gens_using_pari()
        else:
            A = self._kernel_matrix_using_padic_algorithm(proof)
            if LLL:
                return A.LLL()
            else:
                return A

    def kernel(self, algorithm='padic', LLL=False, proof=None, echelonize=True):
        r"""
        Return the left kernel of this matrix, as a module over the
        integers.  This is the saturated ZZ-module spanned by all the
        row vectors v such that v*self = 0.

        INPUT:
            algorithm -- 'padic': a new p-adic based algorithm
                         'pari': use PARI
            LLL -- bool (default: False); if True the basis is an LLL
                   reduced basis; otherwise, it is an echelon basis.
            proof -- None (default: proof.linear_algebra()); if False,
                   impacts how determinants are computed.

        By convention if self has 0 rows, the kernel is of dimension
        0, whereas the kernel is the whole domain if self has 0 columns.

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
            return M.zero_submodule()
        elif self._ncols == 0:  # to a 0 space
            return sage.modules.free_module.FreeModule(ZZ, self._nrows)

        X = self.kernel_matrix(algorithm=algorithm, LLL=LLL, proof=proof)
        if not LLL and echelonize:
            X = X.hermite_form(proof=proof)
        X = X.rows()

        M = sage.modules.free_module.FreeModule(ZZ, self.nrows())
        if LLL:
            return M.span_of_basis(X, check=False)
        else:
            return M.span(X, check=False)

    def _kernel_matrix_using_padic_algorithm(self, proof):
        """
        Compute a list of independent generators that span the right kernel
        of self.

        ALGORITHM: Use IML to compute the kernel over QQ, clear denominators,
        then saturate.

        """
        return self.transpose()._rational_kernel_iml().transpose().saturation(proof=proof)


    def _kernel_gens_using_pari(self):
        """
        Compute an LLL reduced list of independent generators that
        span the kernel of self.

        ALGORITHM: Call pari's matkerint function.
        """
        A = self._pari_().mattranspose()
        B = A.matkerint()

        if B.ncols() == 0:
            return []

        # Now B is a basis for the LLL-reduced integer kernel as a
        # PARI object.  The basis vectors or B[0], ..., B[n-1],
        # where n is the dimension of the kernel.
        M = sage.modules.free_module.FreeModule(ZZ, self.nrows())
        X = []
        for b in B:
            tmp = []
            for x in b:
                tmp.append(ZZ(x))
            X.append(M(tmp))

        return X

    def _adjoint(self):
        """
        Return the adjoint of this matrix.

        Assumes self is a square matrix (checked in adjoint).
        """
        return self.parent()(self._pari_().matadjoint().python())

    def _ntl_(self):
        r"""
        ntl.mat_ZZ representation of self.

        EXAMPLE:
            sage: a = MatrixSpace(ZZ,200).random_element(x=-2, y=2)    # -2 to 2
            sage: A = a._ntl_()

        \note{NTL only knows dense matrices, so if you provide a
        sparse matrix NTL will allocate memory for every zero entry.}
        """
        import sage.libs.ntl.ntl_mat_ZZ
        return sage.libs.ntl.ntl_mat_ZZ.ntl_mat_ZZ(self._nrows,self._ncols, self.list())


    ####################################################################################
    # LLL
    ####################################################################################
    def LLL_gram(self):
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
            sage: U = M.LLL_gram(); U
            [-1  1]
            [ 1 -2]
            sage: U.transpose() * M * U
            [1 0]
            [0 1]

        Semidefinite and indefinite forms raise a ValueError:

            sage: Matrix(ZZ,2,2,[2,6,6,3]).LLL_gram()
            Traceback (most recent call last):
            ...
            ValueError: not a definite matrix
            sage: Matrix(ZZ,2,2,[1,0,0,-1]).LLL_gram()
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

    def LLL(self, delta=None, eta=None, algorithm="fpLLL:wrapper", fp=None, prec=0, early_red = False, use_givens = False):
        r"""
        Returns LLL reduced or approximated LLL reduced lattice R for
        this matrix interpreted as a lattice.

        A lattice $(b_1, b_2, ..., b_d)$ is $(\delta, \eta)$-LLL-reduced
        if the two following conditions hold:

        (a) For any $i>j$, we have $|mu_{i, j}| <= \eta$,
        (b) For any $i<d$, we have
        $\delta |b_i^*|^2 <= |b_{i + 1}^* + mu_{i + 1, i} b_{i + 1}^* |^2$,

        where $mu_{i,j} = <b_i, b_j^*>/<b_j^*,b_j^*>$ and $b_i^*$ is
        the $i$-th vector of the Gram-Schmidt orthogonalisation of
        $(b_1, b_2, ..., b_d)$.

        The default reduction parameters are $\delta=3/4$ and
        $eta=0.501$. The parameters $\delta$ and $\eta$ must satisfy:
        $0.25 < \delta <= 1.0$ and $0.5 <= \eta <
        sqrt(\delta)$. Polynomial time complexity is only guaranteed
        for $\delta < 1$.

        The lattice is returned as a matrix. Also the rank (and the
        determinant) of self are cached if those are computed during
        the reduction. Note that in general this only happens when
        self.rank() == self.ncols() and the exact algorithm is used.

        INPUT:
            delta -- parameter as described above (default: 3/4)
            eta -- parameter as described above (default: 0.501), ignored
                   by NTL
            algorithm -- string (default: "fpLLL:wrapper")
                         one of the algorithms mentioned below
            fp -- None -- NTL's exact reduction or fpLLL's wrapper
               'fp' -- double precision: NTL's FP or fpLLL's double
               'qd' -- quad doubles: NTL's QP
               'xd' -- extended exponent: NTL's XD or fpLLL's dpe
               'rr' -- arbitrary precision: NTL'RR or fpLLL's MPFR
            prec -- precision, ignored by NTL (default: auto choose)
            early_red -- perform early reduction, ignored by NTL (default: False)
            use_givens -- use Givens orthogonalization (default: False)
                          only applicable to approximate reductions and NTL.
                          This is more stable but slower.

        Also, if the verbose level is >= 2, some more verbose output
        is printed during the calculation if NTL is used.

        AVAILABLE ALGORITHMS:
            NTL:LLL -- NTL's LLL + fp
            fpLLL:heuristic -- fpLLL's heuristic + fp
            fpLLL:fast -- fpLLL's fast
            fpLLL:wrapper -- fpLLL's automatic choice (default)

        OUTPUT:
            a matrix over the integers

        EXAMPLE:
            sage: A = Matrix(ZZ,3,3,range(1,10))
            sage: A.LLL()
            [ 0  0  0]
            [ 2  1  0]
            [-1  1  3]

        We compute the extended GCD of a list of integers using LLL,
        this example is from the Magma handbook:

        sage: Q = [ 67015143, 248934363018, 109210, 25590011055, 74631449, \
                    10230248, 709487, 68965012139, 972065, 864972271 ]
        sage: n = len(Q)
        sage: S = 100
        sage: X = Matrix(ZZ, n, n + 1)
        sage: for i in xrange(n):
        ...       X[i,i + 1] = 1
        sage: for i in xrange(n):
        ...       X[i,0] = S*Q[i]
        sage: L = X.LLL()
        sage: M = L.row(n-1).list()[1:]
        sage: M
        [-3, -1, 13, -1, -4, 2, 3, 4, 5, -1]
        sage: add([Q[i]*M[i] for i in range(n)])
        -1

        ALGORITHM: Uses the NTL library by Victor Shoup or fpLLL
        library by Damien Stehle depending on the chosen algorithm.

        REFERENCES: \code{ntl.mat_ZZ} or \code{sage.libs.fplll.fplll}
            for details on the used algorithms.
        """
        tm = verbose("LLL of %sx%s matrix (algorithm %s)"%(self.nrows(), self.ncols(), algorithm))
        import sage.libs.ntl.all
        ntl_ZZ = sage.libs.ntl.all.ZZ

        from sage.libs.fplll.fplll import FP_LLL

        if get_verbose() >= 2: verb = True
        else: verb = False

        # auto choice

        # FP choice
        if algorithm == 'NTL:LLL':
            if fp == None:
                algorithm = 'NTL:LLL_FP'
            elif fp == 'fp':
                algorithm = 'NTL:LLL_FP'
            elif fp == 'qd':
                algorithm = 'NTL:LLL_QD'
            elif fp == 'xd':
                algorithm = 'NTL:LLL_XD'
            elif fp == 'rr':
                algorithn = 'NTL:LLL_RR'
        elif algorithm == 'fpLLL:heuristic':
            if fp == None:
                raise TypeError, "if 'fpLLL:heuristic' is chosen, a floating point number implementation must be chosen"
            elif fp == 'fp':
                fp = 'double'
            elif fp == 'qd':
                raise TypeError, "fpLLL does not support quad doubles."
            elif fp == 'xd':
                fp = 'dpe'
            elif fp == 'rr':
                fp = 'mpfr'

        if algorithm == "NTL:LLL":
            if delta is None:
                delta = ZZ(3)/ZZ(4)
            elif delta <= ZZ(1)/ZZ(4):
                raise TypeError, "delta must be > 1/4"
            elif delta > 1:
                raise TypeError, "delta must be <= 1"
            delta = QQ(delta)
            a = delta.numer()
            b = delta.denom()

        else:
            if delta is None:
                delta = 0.99
            elif delta <= 0.25:
                raise TypeError, "delta must be > 0.25"
            elif delta > 1:
                raise TypeError, "delta must be <= 1"
            delta = float(delta)

            if eta is None:
                eta = 0.501
            elif eta < 0.5:
                raise TypeError, "eta must be >= 0.5"

        if prec < 0:
            raise TypeError, "precision prec must be >= 0"
        int(prec)

        if algorithm.startswith('NTL:'):
            A = sage.libs.ntl.all.mat_ZZ(self.nrows(),self.ncols(),map(ntl_ZZ,self.list()))

            if algorithm == "NTL:LLL":
                r, det2 = A.LLL(a,b, verbose=verb)
                det2 = ZZ(det2)
                try:
                    det = ZZ(det2.sqrt_approx())
                    self.cache("det", det)
                except TypeError:
                    pass
            elif algorithm == "NTL:LLL_FP":
                if use_givens:
                    r = A.G_LLL_FP(delta, verbose=verb)
                else:
                    r = A.LLL_FP(delta, verbose=verb)
            elif algorithm == "NTL:LLL_QP":
                if use_givens:
                    r = A.G_LLL_QP(delta, verbose=verb)
                else:
                    r = A.LLL_QP(delta, verbose=verb)
            elif algorithm == "NTL:LLL_XD":
                if use_givens:
                    r = A.G_LLL_XD(delta, verbose=verb)
                else:
                    r = A.LLL_XD(delta, verbose=verb)
            elif algorithm == "NTL:LLL_RR":
                if use_givens:
                    r = A.G_LLL_RR(delta, verbose=verb)
                else:
                    r = A.LLL_XD(delta, verbose=verb)
            else:
                raise TypeError, "algorithm %s not supported"%algorithm

            r = ZZ(r)

            R = <Matrix_integer_dense>self.new_matrix(entries=map(ZZ,A.list()))
            self.cache("rank",r)

        elif algorithm.startswith('fpLLL:'):

            A = sage.libs.fplll.fplll.FP_LLL(self)
            if algorithm == 'fpLLL:wrapper':
                A.wrapper(prec, eta, delta)
            elif algorithm == 'fpLLL:heuristic':
                if early_red:
                    A.heuristic_early_red(prec,eta,delta,fp)
                else:
                    A.heuristic(prec,eta,delta,fp)
            elif algorithm == 'fpLLL:fast':
                if early_red:
                    A.fast_early_red(prec,eta,delta)
                else:
                    A.fast(prec,eta,delta)
            elif algorithm == 'fpLLL:proved':
                A.proved(prec,eta,delta)
            else:
                raise TypeError, "algorithm %s not supported"%algorithm
            R = A._sage_()
        else:
            raise TypeError, "algorithm %s not supported"%algorithm

        verbose("LLL finished", tm)
        return R

    def prod_of_row_sums(self, cols):
        """
        Return the product of the sums of the entries in the submatrix
        of self with given columns.

        INPUT:
            cols -- a list (or set) of integers representing columns of self.

        OUTPUT:
            an integer

        EXAMPLES:
            sage: a = matrix(ZZ,2,3,[1..6]); a
            [1 2 3]
            [4 5 6]
            sage: a.prod_of_row_sums([0,2])
            40
            sage: (1+3)*(4+6)
            40
            sage: a.prod_of_row_sums(set([0,2]))
            40
        """
        cdef Py_ssize_t c, row
        cdef mpz_t s, pr
        mpz_init(s)
        mpz_init(pr)

        mpz_set_si(pr, 1)
        for row from 0 <= row < self._nrows:
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
        v = ['%s %s M'%(self._nrows, self._ncols)]
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

        INPUT:
            N -- an integer

        OUTPUT:
            matrix -- over QQ or raise a ValueError

        EXAMPLES:
        We create a random 4x4 matrix over ZZ.
            sage: A = matrix(ZZ, 4, [4, -4, 7, 1, -1, 1, -1, -12, -1, -1, 1, -1, -3, 1, 5, -1])

        There isn't a unique rational reconstruction of it:
            sage: A.rational_reconstruction(11)
            Traceback (most recent call last):
            ...
            ValueError: Rational reconstruction of 4 (mod 11) does not exist.

        We throw in a denominator and reduce the matrix modulo 389 --
        it does rationally reconstruct.
            sage: B = (A/3 % 389).change_ring(ZZ)
            sage: B.rational_reconstruction(389) == A/3
            True
        """
        import misc
        return misc.matrix_integer_dense_rational_reconstruction(self, N)

    def randomize(self, density=1, x=None, y=None, distribution=None):
        """
        Randomize density proportion of the entries of this matrix,
        leaving the rest unchanged.

        The parameters are the same as the integer ring's random_element function.

        If x and y are given, randomized entries of this matrix to be between x and y
        and have density 1.

        INPUT:
            self -- a mutable matrix over ZZ
            density -- a float between 0 and 1
            x, y -- if not None are passed to ZZ.random_element function as the
                    upper and lower endpoints in the uniform distribution
            distribution -- would also be passed into ZZ.random_element if given

        OUTPUT:
            -- modifies this matrix in place

        EXAMPLES:
            sage: A = matrix(ZZ, 2,3, [1..6]); A
            [1 2 3]
            [4 5 6]
            sage: A.randomize()
            sage: A   # random output
            [ 1 -5 -1]
            [-2  4 -2]
            sage: A.randomize(x=-30,y=30)
            sage: A   # random output
            [-28 -12 -16]
            [ 24   2  21]
        """
        self.check_mutability()
        self.clear_cache()

        density = float(density)

        cdef Py_ssize_t i, j, k, nc, num_per_row
        global state, ZZ

        cdef IntegerRing_class the_integer_ring = ZZ

        _sig_on
        if density == 1:
            for i from 0 <= i < self._nrows*self._ncols:
                the_integer_ring._randomize_mpz(self._entries[i], x, y, distribution)

        else:
            nc = self._ncols
            num_per_row = int(density * nc)
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < num_per_row:
                    k = random()%nc
                    the_integer_ring._randomize_mpz(self._matrix[i][k], x, y, distribution)

        _sig_off

    #### Rank

    def rank(self):
        """
        Return the rank of this matrix.

        OUTPUT:
            nonnegative integer -- the rank

        NOTE: The rank is cached.

        ALGORITHM: First check if the matrix has maxim posible rank by
        working modulo one random prime.  If not call Linbox's rank
        function.

        EXAMPLES:
            sage: a = matrix(ZZ,2,3,[1..6]); a
            [1 2 3]
            [4 5 6]
            sage: a.rank()
            2
            sage: a = matrix(ZZ,3,3,[1..9]); a
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: a.rank()
            2

        Here's a bigger example -- the rank is of course still 2:
            sage: a = matrix(ZZ,100,[1..100^2]); a.rank()
            2
        """
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

    #### Determinant

    def determinant(self, algorithm='padic', proof=None, stabilize=2):
        r"""
        Return the determinant of this matrix.

        INPUT:
            algorithm --
                'padic' -- (the default) uses a p-adic / multimodular
                           algorithm that relies on code in IML and linbox
                'linbox' -- calls linbox det
                'ntl' -- calls NTL's det function
            proof -- bool or None; if None use proof.linear_algebra(); only
                     relevant for the padic algorithm.
                     NOTE: It would be *VERY VERY* hard for det to fail
                     even with proof=False.
            stabilize -- if proof is False, require det to be the same
                         for this many CRT primes in a row.  Ignored
                         if proof is True.

        ALGORITHM:
            The p-adic algorithm works by first finding a random
            vector v, then solving A*x = v and taking the denominator
            $d$.  This gives a divisor of the determinant.  Then we
            compute $\det(A)/d$ using a multimodular algorithm and the
            Hadamard bound, skipping primes that divide $d$.

        TIMINGS:
            This is perhaps the fastest implementation of determinants
            in the world.  E.g., for a 500x500 random matrix with
            32-bit entries on a core2 duo 2.6Ghz running OS X, Sage
            takes 4.12 seconds, whereas Magma takes 62.87 seconds
            (both with proof False).  With proof=True on the same problem
            Sage takes 5.73 seconds.  For another example, a 200x200
            random matrix with 1-digit entries takes 4.18 seconds in
            pari, 0.18 in Sage with proof True, 0.11 in Sage with
            proof False, and 0.21 seconds in Magma with proof True
            and 0.18 in Magma with proof False.

        EXAMPLES:
            sage: A = matrix(ZZ,8,8,[3..66])
            sage: A.determinant()
            0

            sage: A = random_matrix(ZZ,20,20)
            sage: D1 = A.determinant()
            sage: A._clear_cache()
            sage: D2 = A.determinant(algorithm='ntl')
            sage: D1 == D2
            True
        """
        d = self.fetch('det')
        if not d is None:
            return d
        if not self.is_square():
            raise ValueError, "self must be square"
        n = self.nrows()

        if n <= 3:
            # use generic special cased code.
            return matrix_dense.Matrix_dense.determinant(self)

        proof = get_proof_flag(proof, "linear_algebra")

        if algorithm == 'padic':
            import matrix_integer_dense_hnf
            return matrix_integer_dense_hnf.det_padic(self, proof=proof, stabilize=stabilize)
        elif algorithm == 'linbox':
            d = self._det_linbox()
        elif algorithm == 'ntl':
            d = self._det_ntl()
        else:
            raise TypeError, "algorithm '%s' not known"%(algorithm)

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

    def _det_ntl(self):
        """
        Compute the determinant of this matrix using NTL.
        """
        _sig_on
        d = self._ntl_().determinant()
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
        if self._nrows == 0 or self._ncols == 0:
            return self.matrix_space(self._ncols, 0).zero_matrix()

        cdef long dim
        cdef mpz_t *mp_N
        time = verbose('computing nullspace of %s x %s matrix using IML'%(self._nrows, self._ncols))
        _sig_on
        dim = nullspaceMP (self._nrows, self._ncols, self._entries, &mp_N)
        _sig_off
        P = self.matrix_space(self._ncols, dim)

        # Now read the answer as a matrix.
        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)
        for i from 0 <= i < dim*self._ncols:
            mpz_init_set(M._entries[i], mp_N[i])
            mpz_clear(mp_N[i])
        free(mp_N)

        verbose("finished computing nullspace", time)
        M._initialized = True
        return M

    def _invert_iml(self, use_nullspace=False, check_invertible=True):
        """
        Invert this matrix using IML.  The output matrix is an integer
        matrix and a denominator.

        INPUT:
           self -- an invertible matrix
           use_nullspace -- (default: False): whether to use nullspace
                     algorithm, which is slower, but doesn't require
                     checking that the matrix is invertible as a
                     precondition.
           check_invertible -- (default: True) whether to check that
                     the matrix is invertible.

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
            if check_invertible and self.rank() != self._nrows:
                raise ZeroDivisionError, "input matrix must be nonsingular"
            return self._solve_iml(P.identity_matrix(), right=True)

    def _solve_right_nonsingular_square(self, B, check_rank=True):
        r"""
        If self is a matrix $A$ of full rank, then this function
        returns a vector or matrix $X$ such that $A X = B$.  If $B$ is
        a vector then $X$ is a vector and if $B$ is a matrix, then $X$
        is a matrix.  The base ring of $X$ is the integers unless a
        denominator is needed in which case the base ring is the
        rational numbers.

        NOTE: In SAGE one can also write \code{A \ B} for
        \code{A.solve_right(B)}, i.e., SAGE implements the ``the
        MATLAB/Octave backslash operator''.

        NOTE: This is currently only implemented when A is square.

        INPUT:
            B -- a matrix or vector
            check_rank -- bool (default: True); if True verify that in
                          fact the rank is full.
        OUTPUT:
            a matrix or vector over $\QQ$

        EXAMPLES:
            sage: a = matrix(ZZ, 2, [0, -1, 1, 0])
            sage: v = vector(ZZ, [2, 3])
            sage: a \ v
            (3, -2)

        Note that the output vector or matrix is always over $\QQ$.
            sage: parent(a\v)
            Vector space of dimension 2 over Rational Field

        We solve a bigger system where the answer is over the rationals.
            sage: a = matrix(ZZ, 3, 3, [1,2,3,4, 5, 6, 8, -2, 3])
            sage: v = vector(ZZ, [1,2,3])
            sage: w = a \ v; w
            (2/15, -4/15, 7/15)
            sage: parent(w)
            Vector space of dimension 3 over Rational Field
            sage: a * w
            (1, 2, 3)

        We solve a system where the right hand matrix has multiple columns.
            sage: a = matrix(ZZ, 3, 3, [1,2,3,4, 5, 6, 8, -2, 3])
            sage: b = matrix(ZZ, 3, 2, [1,5, 2, -3, 3, 0])
            sage: w = a \ b; w
            [ 2/15 -19/5]
            [-4/15 -27/5]
            [ 7/15 98/15]
            sage: a * w
            [ 1  5]
            [ 2 -3]
            [ 3  0]

        TESTS:
        We create a random 100x100 matrix and solve the corresponding system,
        then verify that the result is correct.  (NOTE: This test is very
        risky without having a seeded random number generator!)
            sage: n = 100
            sage: a = random_matrix(ZZ,n)
            sage: v = vector(ZZ,n,range(n))
            sage: x = a \ v
            sage: a * x == v
            True

        """
        t = verbose('starting IML solve_right...')
        # It would probably be much better to rewrite linbox so it
        # throws an error instead of ** going into an infinite loop **
        # in the non-full rank case.  In any case, we do this for now,
        # since rank is very fast and infinite loops are evil.
        if check_rank and self.rank() < self.nrows():
            raise ValueError, "self must be of full rank."

        if not self.is_square():
            raise NotImplementedError, "the input matrix must be square."

        if is_Vector(B):
            if self.nrows() != B.degree():
                raise ValueError, "number of rows of self must equal degree of B."
        elif self.nrows() != B.nrows():
                raise ValueError, "number of rows of self must equal number of rows of B."

        if self.nrows() == 0:
            return B

        matrix = True
        C = B
        if not isinstance(B, Matrix_integer_dense):
            if is_Vector(B):
                matrix = False
                C = self.matrix_space(self.nrows(), 1)(B.list())
            else:
                raise NotImplementedError

        if C.ncols() >= 2*self.ncols():
            # likely much better to just invert then multiply
            X = self**(-1)*C
            verbose('finished solve_right (via inverse)', t)
            return X

        X, d = self._solve_iml(C, right=True)
        if d != 1:
            X = (1/d) * X
        if not matrix:
            # Convert back to a vector
            X = (X.base_ring() ** X.nrows())(X.list())
        verbose('finished solve_right via IML', t)
        return X

    def _solve_iml(self, Matrix_integer_dense B, right=True):
        """
        Let A equal self be a square matrix. Given B return an integer
        matrix C and an integer d such that self C*A == d*B if right
        is False or A*C == d*B if right is True.

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
            # This is *required* by the IML function we call below.
            raise ArithmeticError, "self must be square"

        if self.nrows() == 1:
            return B, self[0,0]

        mpz_init(mp_D)


        if right:
            if self._ncols != B._nrows:
                raise ArithmeticError, "B's number of rows must match self's number of columns"

            n = self._ncols
            m = B._ncols
            P = self.matrix_space(n, m)
            if self._nrows == 0 or self._ncols == 0:
                return P.zero_matrix(), Integer(1)

            if m == 0 or n == 0:
                return self.new_matrix(nrows = n, ncols = m), Integer(1)

            mp_N = <mpz_t *> sage_malloc( n * m * sizeof(mpz_t) )
            for i from 0 <= i < n * m:
                mpz_init( mp_N[i] )

            nonsingSolvLlhsMM(RightSolu, n, m, self._entries, B._entries, mp_N, mp_D)

        else: # left
            if self._nrows != B._ncols:
                raise ArithmeticError, "B's number of columns must match self's number of rows"

            n = self._ncols
            m = B._nrows

            P = self.matrix_space(m, n)
            if self._nrows == 0 or self._ncols == 0:
                return P.zero_matrix(), Integer(1)

            if m == 0 or n == 0:
                return self.new_matrix(nrows = m, ncols = n), Integer(1)

            mp_N = <mpz_t *> sage_malloc( n * m * sizeof(mpz_t) )
            for i from 0 <= i < n * m:
                mpz_init( mp_N[i] )

            nonsingSolvLlhsMM(LeftSolu, n, m, self._entries, B._entries, mp_N, mp_D)


        M = Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)
        for i from 0 <= i < n*m:
            mpz_init_set(M._entries[i], mp_N[i])
            mpz_clear(mp_N[i])
        sage_free(mp_N)
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
        if self._nrows == 0:
            pivots = []
            nonpivots = range(self._ncols)
            X = self.copy()
            d = Integer(1)
            return pivots, nonpivots, X, d
        elif self._ncols == 0:
            pivots = []
            nonpivots = []
            X = self.copy()
            d = Integer(1)
            return pivots, nonpivots, X, d

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
            return decomp_seq([(V, F[0][1]==1)])

        A = self.change_ring(QQ)
        X = A.decomposition(**kwds)
        V = ZZ**self.nrows()
        if isinstance(X, tuple):
            D, E = X
            D = [(W.intersection(V), t) for W, t in D]
            E = [(W.intersection(V), t) for W, t in E]
            return decomp_seq(D), decomp_seq(E)
        else:
            return decomp_seq([(W.intersection(V), t) for W, t in X])

    def _add_row_and_maintain_echelon_form(self, row, pivots):
        """
        Assuming self is a full rank n x m matrix in reduced row
        Echelon form over ZZ and row is a vector of degree m, this
        function creates a new matrix that is the echelon form of self
        with row appended to the bottom.

        WARNING: It is assumed that self is in echelon form.

        INPUT:
            row -- a vector of degree m over ZZ
            pivots -- a list of integers that are the pivot columns of self.

        OUTPUT:
            matrix -- a matrix of in reduced row echelon form over ZZ
            pivots -- list of integers

        ALGORITHM: For each pivot column of self, we use the extended
        Euclidean algorithm to clear the column.  The result is a new
        matrix B whose row span is the same as self.stack(row), and
        whose last row is 0 if and only if row is in the QQ-span of
        the rows of self.  If row is not in the QQ-span of the rows of
        self, then row is nonzero and suitable to be inserted into the
        top n rows of A to form a new matrix that is in reduced row
        echelon form.  We then clear that corresponding new pivot column.

        EXAMPLES:
            sage: a = matrix(ZZ, 3, [1, 0, 110, 0, 3, 112, 0, 0, 221]); a
            [  1   0 110]
            [  0   3 112]
            [  0   0 221]
            sage: a._add_row_and_maintain_echelon_form(vector(ZZ,[1,2,3]),[0,1,2])
            ([1 0 0]
            [0 1 0]
            [0 0 1], [0, 1, 2])
            sage: a._add_row_and_maintain_echelon_form(vector(ZZ,[0,0,0]),[0,1,2])
            ([  1   0 110]
            [  0   3 112]
            [  0   0 221], [0, 1, 2])
            sage: a = matrix(ZZ, 2, [1, 0, 110, 0, 3, 112])
            sage: a._add_row_and_maintain_echelon_form(vector(ZZ,[1,2,3]),[0,1])
            ([  1   0 110]
            [  0   1 219]
            [  0   0 545], [0, 1, 2])
        """
        from sage.all import get_memory_usage
        cdef Py_ssize_t i, j, piv, n = self._nrows, m = self._ncols

        import constructor

        # 0. Base case
        if self.nrows() == 0:
            pos = row.nonzero_positions()
            if len(pos) > 0:
                pivots = [pos[0]]
                if row[pivots[0]] < 0:
                    row *= -1
            else:
                pivots = []
            return constructor.matrix([row]), pivots


        if row == 0:
            return self, pivots
        # 1. Create a new matrix that has row as the last row.
        row_mat = constructor.matrix(row)
        A = self.stack(row_mat)
        # 2. Working from the left, clear each column to put
        #    the resulting matrix back in echelon form.
        for i, p in enumerate(pivots):
            # p is the i-th pivot
            b = A[n,p]
            if not b:
                continue

            # (a). Take xgcd of pivot positions in last row and in ith
            # row.

            # TODO (optimize) -- change to use direct call to gmp and
            # no bounds checking!
            a = A[i,p]
            if not a:
                raise ZeroDivisionError, "claimed pivot is not a pivot"
            if b % a == 0:
                # (b) Subtract a multiple of row i from row n.
                c = b // a
                if c:
                    for j in range(m):
                        A[n,j] -= c * A[i,j]
            else:
                # (b). More elaborate.
                #  Replace the ith row by s*A[i] + t*A[n], which will
                # have g in the i,p position, and replace the last row by
                # (b//g)*A[i] - (a//g)*A[n], which will have 0 in the i,p
                # position.
                g, s, t = a.xgcd(b)
                if not g:
                    raise ZeroDivisionError, "claimed pivot is not a pivot (got a 0 gcd)"

                row_i = A.row(i)
                row_n = A.row(n)

                ag = a//g; bg = b//g

                new_top = s*row_i  +  t*row_n
                new_bot = bg*row_i - ag*row_n


                # OK -- now we have to make sure the top part of the matrix
                # but with row i replaced by
                #     r = s*row_i[j]  +  t*row_n[j]
                # is put in rref.  We do this by recursively calling this
                # function with the top part of A (all but last row) and the
                # row r.

                zz = range(A.nrows()-1)
                del zz[i]
                top_mat = A.matrix_from_rows(zz)
                new_pivots = list(pivots)
                del new_pivots[i]

                top_mat, pivots = top_mat._add_row_and_maintain_echelon_form(new_top, new_pivots)
                w = top_mat._add_row_and_maintain_echelon_form(new_bot, pivots)
                return w
        # 3. If it turns out that the last row is nonzero,
        #    insert last row in A sliding other rows down.
        v = A.row(n)
        new_pivots = list(pivots)
        if v != 0:
            i = v.nonzero_positions()[0]
            assert not (i in pivots), 'WARNING: bug in add_row -- i (=%s) should not be a pivot'%i

            # If pivot entry is negative negate this row.
            if v[i] < 0:
                v = -v

            # Determine where the last row should be inserted.
            new_pivots.append(i)
            new_pivots.sort()
            import bisect
            j = bisect.bisect(pivots, i)
            # The new row should go *before* row j, so it becomes row j
            A = A.insert_row(j, v)
        try:
            _clear_columns(A, new_pivots, A.nrows())
        except RuntimeError:
            raise ZeroDivisionError, "mistake in claimed pivots"
        if A.row(A.nrows() - 1) == 0:
            A = A.matrix_from_rows(range(A.nrows()-1))
        return A, new_pivots

    #####################################################################################
    # Hermite form modulo D
    # This code below is by E. Burcin.  Thanks!
    #####################################################################################
    cdef _new_uninitialized_matrix(self, Py_ssize_t nrows, Py_ssize_t ncols):
        """
        Return a new matrix over the integers with the given number of rows and columns.
        All memory is allocated for this matrix, but its entries have not yet been
        filled in.
        """
        P = self._parent.matrix_space(nrows, ncols)
        return Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)

    def _hnf_mod(self, D):
        """
        INPUT:
            D -- a small integer that is assumed to be a multiple of 2*det(self)
        OUTPUT:
            matrix -- the Hermite normal form of self.
        """
        t = verbose('hermite mod %s'%D, caller_name='matrix_integer_dense')
        cdef Matrix_integer_dense res = self._new_uninitialized_matrix(self._nrows, self._ncols)
        self._hnf_modn(res, D)
        verbose('finished hnf mod', t, caller_name='matrix_integer_dense')
        return res

    cdef int _hnf_modn(Matrix_integer_dense self, Matrix_integer_dense res,
            mod_int det) except -1:
        """
        Puts self into HNT form modulo det.  Changes self in place.
        """
        cdef long long *res_l
        cdef Py_ssize_t i,j,k
        res_l = self._hnf_modn_impl(det, self._nrows, self._ncols)
        k = 0
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                mpz_init_set_si(res._matrix[i][j], res_l[k])
                k += 1
        res._initialized = True
        sage_free(res_l)


    cdef long long* _hnf_modn_impl(Matrix_integer_dense self, mod_int det,
            Py_ssize_t nrows, Py_ssize_t ncols):
        cdef long long *res, *T_ent, **res_rows, **T_rows, *B
        cdef Py_ssize_t i, j, k
        cdef long long R, mod, T_i_i, T_j_i, c1, c2, q, t
        cdef int u, v, d
        cdef mpz_t m

        # allocate memory for result matrix
        res = <long long*> sage_malloc(sizeof(long long)*ncols*nrows)
        if res == NULL:
            raise MemoryError, "out of memory allocating a matrix"
        res_rows = <long long**> sage_malloc(sizeof(long long*)*nrows)
        if res_rows == NULL:
            sage_free(res)
            raise MemoryError, "out of memory allocating a matrix"

        # allocate memory for temporary matrix
        T_ent = <long long*> sage_malloc(sizeof(long long)*ncols*nrows)
        if T_ent == NULL:
            sage_free(res)
            sage_free(res_rows)
            raise MemoryError, "out of memory allocating a matrix"
        T_rows = <long long**> sage_malloc(sizeof(long long*)*nrows)
        if T_rows == NULL:
            sage_free(res)
            sage_free(res_rows)
            sage_free(T_ent)
            raise MemoryError, "out of memory allocating a matrix"

        # allocate memory for temporary row vector
        B = <long long*>sage_malloc(sizeof(long long)*nrows)
        if B == NULL:
            sage_free(res)
            sage_free(res_rows)
            sage_free(T_ent)
            sage_free(T_rows)
            raise MemoryError, "out of memory allocating a matrix"

        # initialize the row pointers
        k = 0
        for i from 0 <= i < nrows:
            res_rows[i] = res + k
            T_rows[i] = T_ent + k
            k += nrows


        mpz_init(m)
        # copy entries from self to temporary matrix
        k = 0
        for i from 0 <= i < nrows:
            for j from 0 <= j < ncols:
                mpz_mod_ui(m, self._matrix[i][j], det)
                T_ent[k] = mpz_get_si(m)
                k += 1
        mpz_clear(m)


        # initialize variables
        i = 0
        j = 0
        R = det

        while 1:
            if j == nrows-1:
                T_i_i = T_rows[i][i]
                d = ai.c_xgcd_int(T_i_i, R, &u, &v)
                for k from 0 <= k < i:
                    res_rows[i][k] = 0
                for k from i <= k < ncols:
                    t = (u*T_rows[i][k])%R
                    if t < 0:
                        t += R
                    res_rows[i][k] = t
                if res_rows[i][i] == 0:
                    res_rows[i][i] = R
                d = res_rows[i][i]
                for j from 0 <= j < i:
                    q = res_rows[j][i]/d
                    for k from i <= k < ncols:
                        u = (res_rows[j][k] - q*res_rows[i][k])%R
                        if u < 0:
                            u += R
                        res_rows[j][k] = u

                R = R/d
                i += 1
                j = i
                if i == nrows :
                    break # return res
                if T_rows[i][i] == 0:
                    T_rows[i][i] = R
                continue


            j += 1
            if T_rows[j][i] == 0:
                continue

            T_i_i = T_rows[i][i]
            T_j_i = T_rows[j][i]
            d = ai.c_xgcd_int(T_i_i , T_j_i, &u, &v)
            if d != T_i_i:
                for k from i <= k < ncols:
                    B[k] = (u*T_rows[i][k] + v*T_rows[j][k])
            c1 = T_i_i/d
            c2 = -T_j_i/d
            for k from i <= k < ncols:
                T_rows[j][k] = (c1*T_rows[j][k] + c2*T_rows[i][k])%R
            if d != T_i_i:
                for k from i <= k < ncols:
                    T_rows[i][k] = B[k]%R

        sage_free(B)
        sage_free(res_rows)
        sage_free(T_ent)
        sage_free(T_rows)
        return res


    #####################################################################################
    # operations with matrices
    #####################################################################################
    def stack(self, other):
        """
        Return the matrix self on top of other:
           [ self  ]
           [ other ]

        EXAMPLES:
            sage: M = Matrix(ZZ, 2, 3, range(6))
            sage: N = Matrix(ZZ, 1, 3, [10,11,12])
            sage: M.stack(N)
            [ 0  1  2]
            [ 3  4  5]
            [10 11 12]
        """
        if self._ncols != other.ncols():
            raise TypeError, "number of columns must be the same"
        if not (self._base_ring is other.base_ring()):
            other = other.change_ring(self._base_ring)
        cdef Matrix_integer_dense A = other
        cdef Matrix_integer_dense M
        M = self.new_matrix(nrows = self._nrows + A._nrows, ncols = self.ncols())
        cdef Py_ssize_t i, k
        k = self._nrows * self._ncols
        for i from 0 <= i < k:
            mpz_set(M._entries[i], self._entries[i])
        for i from 0 <= i < A._nrows * A._ncols:
            mpz_set(M._entries[k + i], A._entries[i])
        return M

    def insert_row(self, Py_ssize_t index, row):
        """
        Create a new matrix from self with.

        INPUT:
            index -- integer
            row -- a vector

        EXAMPLES:
            sage: X = matrix(ZZ,3,range(9)); X
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: X.insert_row(1, [1,5,-10])
            [  0   1   2]
            [  1   5 -10]
            [  3   4   5]
            [  6   7   8]
            sage: X.insert_row(0, [1,5,-10])
            [  1   5 -10]
            [  0   1   2]
            [  3   4   5]
            [  6   7   8]
            sage: X.insert_row(3, [1,5,-10])
            [  0   1   2]
            [  3   4   5]
            [  6   7   8]
            [  1   5 -10]
        """
        cdef Matrix_integer_dense res = self._new_uninitialized_matrix(self._nrows + 1, self._ncols)
        cdef Py_ssize_t j, k
        cdef Integer z
        if index < 0:
            raise ValueError, "index must be nonnegative"
        if index > self._nrows:
            raise ValueError, "index must be less than number of rows"
        for j from 0 <= j < self._ncols * index:
            mpz_init_set(res._entries[j], self._entries[j])

        k = 0
        for j from self._ncols * index <= j < self._ncols * (index+1):
            z = row[k]
            mpz_init_set(res._entries[j], z.value)
            k += 1

        for j from self._ncols * (index+1) <= j < (self._nrows + 1)*self._ncols:
            mpz_init_set(res._entries[j], self._entries[j - self._ncols])

        res._initialized = True
        return res

    def _factor_out_common_factors_from_each_row(self):
        """
        Very very quickly modifies self so that the gcd of the entries
        in each row is 1 by dividing each row by the common gcd.

        EXAMPLES:
            sage: a = matrix(ZZ, 3, [-9, 3, -3, -36, 18, -5, -40, -5, -5, -20, -45, 15, 30, -15, 180])
            sage: a
            [ -9   3  -3 -36  18]
            [ -5 -40  -5  -5 -20]
            [-45  15  30 -15 180]
            sage: a._factor_out_common_factors_from_each_row()
            sage: a
            [ -3   1  -1 -12   6]
            [ -1  -8  -1  -1  -4]
            [ -3   1   2  -1  12]
        """
        self.check_mutability()

        cdef mpz_t g
        mpz_init(g)
        cdef Py_ssize_t i, j
        cdef mpz_t* row

        for i from 0 <= i < self._nrows:
            mpz_set_ui(g, 0)
            row = self._matrix[i]
            for j from 0 <= j < self._ncols:
                mpz_gcd(g, g, row[j])
                if mpz_cmp_ui(g, 1) == 0:
                    break
            if mpz_cmp_ui(g, 1) != 0:
                # divide through row
                for j from 0 <= j < self._ncols:
                    mpz_divexact(row[j], row[j], g)
        mpz_clear(g)

    def gcd(self):
        """
        Return the gcd of all entries of self; very fast.

        EXAMPLES:
            sage: a = matrix(ZZ,2, [6,15,-6,150])
            sage: a.gcd()
            3
        """
        cdef Integer g = Integer(0)
        cdef Py_ssize_t i, j
        cdef mpz_t* row

        for i from 0 <= i < self._nrows:
            row = self._matrix[i]
            for j from 0 <= j < self._ncols:
                mpz_gcd(g.value, g.value, row[j])
                if mpz_cmp_ui(g.value, 1) == 0:
                    return g
        return g

    def _change_ring(self, ring):
        """
        Return the matrix obtained by coercing the entries of this
        matrix into the given ring.

        EXAMPLES:
            sage: a = matrix(ZZ,2,[1,-7,3,5])
            sage: a._change_ring(RDF)
            [ 1.0 -7.0]
            [ 3.0  5.0]
        """
        if ring == RDF:
            import change_ring
            return change_ring.integer_to_real_double_dense(self)
        else:
            raise NotImplementedError

##     def augment(self, other):
##         """
##         """
##         if self._nrows != other.nrows():
##             raise TypeError, "number of rows must be the same"
##         if not (self._base_ring is other.base_ring()):
##             other = other.change_ring(self._base_ring)
##         cdef Matrix_integer_dense A = other
##         cdef Matrix_integer_dense M
##         M = self.new_matrix(nrows = self._nrows, ncols = self._ncols + A._ncols)
##         cdef Py_ssize_t i, k
##         k = self._nrows * self._ncols
##         for i from 0 <= i < k:
##             mpz_set(M._entries[i], self._entries[i])
##         for i from 0 <= i < A._nrows * A._ncols:
##             mpz_set(M._entries[k + i], A._entries[i])
##         return M


    #####################################################################################

cdef _clear_columns(Matrix_integer_dense A, pivots, Py_ssize_t n):
    # Clear all columns
    cdef Py_ssize_t i, k, p, l, m = A._ncols
    cdef mpz_t** matrix = A._matrix
    cdef mpz_t c, t
    _sig_on
    mpz_init(c)
    mpz_init(t)
    for i from 0 <= i < len(pivots):
        p = pivots[i]
        for k from 0 <= k < n:
            if k != i:
                if mpz_cmp_si(matrix[k][p],0):
                    mpz_fdiv_q(c, matrix[k][p], matrix[i][p])
                    # subtract off c*v from row k; resulting A[k,i] entry will be < b, hence in Echelon form.
                    for l from 0 <= l < m:
                        mpz_mul(t, c, matrix[i][l])
                        mpz_sub(matrix[k][l], matrix[k][l], t)
    mpz_clear(c)
    mpz_clear(t)
    _sig_off

###############################################################

###########################################
# Helper code for Echelon form algorithm.
###########################################
def _parimatrix_to_strlist(A):
    s = str(A)
    s = s.replace('Mat(','').replace(')','')
    # Deal correctly with an empty pari matrix [;]
    if s=='[;]':
        return []
    s = s.replace(';',',').replace(' ','')
    s = s.replace(",", "','")
    s = s.replace("[", "['")
    s = s.replace("]", "']")
    return eval(s)

def _parimatrix_to_reversed_strlist(A):
    s = str(A)
    if s.find('Mat') != -1:
        return _parimatrix_to_strlist(A)
    # Deal correctly with an empty pari matrix [;]
    if s=='[;]':
        return []
    s = s.replace('[','').replace(']','').replace(' ','')
    v = s.split(';')
    v.reverse()
    s = "['" + (','.join(v)) + "']"
    s = s.replace(",", "','")
    return eval(s)

def _convert_parimatrix(z):
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

    cdef mod_int **row_list
    row_list = <mod_int**>sage_malloc(sizeof(mod_int*) * n)
    if row_list == NULL:
        raise MemoryError, "out of memory allocating multi-modular coefficent list"

    _sig_on
    for i from 0 <= i < nr:
        for k from 0 <= k < n:
            row_list[k] = (<Matrix_modn_dense>residues[k])._matrix[i]
        mm.mpz_crt_vec(M._matrix[i], row_list, nc)
    _sig_off

    sage_free(row_list)
    return M

#######################################################

# Conclusions:
#  On OS X Intel, at least:
#    - if log_2(height) >= 0.70 * nrows, use classical

def tune_multiplication(k, nmin=10, nmax=200, bitmin=2,bitmax=64):
    """
    Compare various multiplication algorithms.

    INPUT:
        k -- integer; affects numbers of trials
        nmin -- integer; smallest matrix to use
        nmax -- integer; largest matrix to use
        bitmin -- integer; smallest bitsize
        bitmax -- integer; largest bitsize

    OUTPUT:
        prints what doing then who wins -- multimodular or classical

    EXAMPLES:
        sage: from sage.matrix.matrix_integer_dense import tune_multiplication
        sage: tune_multiplication(2, nmin=10, nmax=60, bitmin=2,bitmax=8)
        10 2 0.2
        ...
    """
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


