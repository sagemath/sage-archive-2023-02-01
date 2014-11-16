r"""
Dense matrices over `\ZZ/n\ZZ` for `n` small

AUTHORS:

- William Stein

- Robert Bradshaw

This is a compiled implementation of dense matrices over
`\ZZ/n\ZZ` for `n` small.

EXAMPLES::

    sage: a = matrix(Integers(37),3,range(9),sparse=False); a
    [0 1 2]
    [3 4 5]
    [6 7 8]
    sage: a.rank()
    2
    sage: type(a)
    <type 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>
    sage: a[0,0] = 5
    sage: a.rank()
    3
    sage: parent(a)
    Full MatrixSpace of 3 by 3 dense matrices over Ring of integers modulo 37

::

    sage: a^2
    [ 3 23 31]
    [20 17 29]
    [25 16  0]
    sage: a+a
    [10  2  4]
    [ 6  8 10]
    [12 14 16]

::

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

::

    sage: TestSuite(a).run()
    sage: TestSuite(b).run()

::

    sage: a.echelonize(); a
    [1 0 0]
    [0 1 0]
    [0 0 1]
    sage: b.echelonize(); b
    [ 1  0 36]
    [ 0  1  2]

We create a matrix group::

    sage: M = MatrixSpace(GF(3),3,3)
    sage: G = MatrixGroup([M([[0,1,0],[0,0,1],[1,0,0]]), M([[0,1,0],[1,0,0],[0,0,1]])])
    sage: G
    Matrix group over Finite Field of size 3 with 2 generators (
    [0 1 0]  [0 1 0]
    [0 0 1]  [1 0 0]
    [1 0 0], [0 0 1]
    )
    sage: G.gap()
    Group([ [ [ 0*Z(3), Z(3)^0, 0*Z(3) ], [ 0*Z(3), 0*Z(3), Z(3)^0 ], [ Z(3)^0, 0*Z(3), 0*Z(3) ] ],
            [ [ 0*Z(3), Z(3)^0, 0*Z(3) ], [ Z(3)^0, 0*Z(3), 0*Z(3) ], [ 0*Z(3), 0*Z(3), Z(3)^0 ] ] ])

TESTS::

    sage: M = MatrixSpace(GF(5),2,2)
    sage: A = M([1,0,0,1])
    sage: A - int(-1)
    [2 0]
    [0 2]
    sage: B = M([4,0,0,1])
    sage: B - int(-1)
    [0 0]
    [0 2]
    sage: Matrix(GF(5),0,0, sparse=False).inverse()
    []
"""

include "sage/ext/interrupt.pxi"
include 'sage/ext/stdsage.pxi'
include 'sage/ext/cdefs.pxi'

from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_mat cimport *
from sage.misc.randstate cimport randstate, current_randstate
from cpython.string cimport *

cimport sage.rings.fast_arith
import sage.rings.fast_arith
cdef sage.rings.fast_arith.arith_int ArithIntObj
ArithIntObj  = sage.rings.fast_arith.arith_int()


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

from sage.matrix.matrix_modn_dense_float cimport Matrix_modn_dense_float
from sage.matrix.matrix_modn_dense_double cimport Matrix_modn_dense_double

from sage.structure.element cimport Matrix

from sage.rings.finite_rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract

from sage.misc.misc import verbose, get_verbose, cputime

from sage.rings.integer cimport Integer

from sage.structure.element cimport ModuleElement, RingElement, Element, Vector

from matrix_integer_dense cimport Matrix_integer_dense
from sage.rings.integer_ring   import ZZ

################
from sage.rings.fast_arith cimport arith_int
cdef arith_int ai
ai = arith_int()
################

from sage.structure.proof.proof import get_flag as get_proof_flag

cdef long num = 1
cdef bint little_endian = (<char*>(&num))[0]

def __matrix_from_rows_of_matrices(X):
    """
    Return a matrix whose ith row is ``X[i].list()``.

    INPUT:

    - ``X`` - a nonempty list of matrices of the same size mod a
       single modulus ``p``

    OUTPUT: A single matrix mod p whose ith row is ``X[i].list()``.


    """
    # The code below is just a fast version of the following:
    ##     from constructor import matrix
    ##     K = X[0].base_ring()
    ##     v = sum([y.list() for y in X],[])
    ##     return matrix(K, len(X), X[0].nrows()*X[0].ncols(), v)

    from matrix_space import MatrixSpace
    cdef Matrix_modn_dense A, T
    cdef Py_ssize_t i, n, m
    n = len(X)

    T = X[0]
    m = T._nrows * T._ncols
    A = Matrix_modn_dense.__new__(Matrix_modn_dense, MatrixSpace(X[0].base_ring(), n, m), 0, 0, 0)
    A.p = T.p
    A.gather = T.gather

    for i from 0 <= i < n:
        T = X[i]
        memcpy(A._entries + i*m, T._entries, sizeof(mod_int)*m)
    return A

cpdef is_Matrix_modn_dense(self):
    """
    """
    return isinstance(self, Matrix_modn_dense) | isinstance(self, Matrix_modn_dense_float) | isinstance(self, Matrix_modn_dense_double)

##############################################################################
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

cdef class Matrix_modn_dense(matrix_dense.Matrix_dense):
    ########################################################################
    # LEVEL 1 functionality
    # x * __cinit__
    # x * __dealloc__
    # x * __init__
    # x * set_unsafe
    # x * get_unsafe
    # x * __richcmp__    -- always the same
    ########################################################################
    def __cinit__(self, parent, entries, copy, coerce):
        from sage.misc.superseded import deprecation
        deprecation(4260, "This class is replaced by Matrix_modn_dense_float/Matrix_modn_dense_double.")

        matrix_dense.Matrix_dense.__init__(self, parent)

        cdef long p
        p = self._base_ring.characteristic()
        self.p = p
        MAX_MODULUS = 2**23
        if p >= MAX_MODULUS:
            raise OverflowError, "p (=%s) must be < %s"%(p, MAX_MODULUS)
        self.gather = MOD_INT_OVERFLOW/<mod_int>(p*p)

        if not isinstance(entries, list):
            sig_on()
            self._entries = <mod_int *> sage_calloc(self._nrows*self._ncols,sizeof(mod_int))
            sig_off()
        else:
            sig_on()
            self._entries = <mod_int *> sage_malloc(self._nrows*self._ncols * sizeof(mod_int))
            sig_off()

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
        """
        TESTS::

            sage: matrix(GF(7), 2, 2, [-1, int(-2), GF(7)(-3), 1/4])
            [6 5]
            [4 2]
        """
        cdef mod_int e
        cdef Py_ssize_t i, j, k
        cdef mod_int *v
        cdef long p
        p = self._base_ring.characteristic()

        R = self.base_ring()

        # scalar?
        if not isinstance(entries, list):
            # sig_on()
            # for i from 0 <= i < self._nrows*self._ncols:
            #     self._entries[i] = 0
            # sig_off()
            if entries is None:
                # zero matrix
                pass
            else:
                e = R(entries)
                if e != 0:
                    for i from 0 <= i < min(self._nrows, self._ncols):
                        self._matrix[i][i] = e
            return

        # all entries are given as a long list
        if len(entries) != self._nrows * self._ncols:
            raise IndexError, "The vector of entries has the wrong length."

        k = 0
        cdef mod_int n
        cdef long tmp

        for i from 0 <= i < self._nrows:
            sig_check()
            v = self._matrix[i]
            for j from 0 <= j < self._ncols:
                x = entries[k]
                if PY_TYPE_CHECK_EXACT(x, int):
                    tmp = (<long>x) % p
                    v[j] = tmp + (tmp<0)*p
                elif PY_TYPE_CHECK_EXACT(x, IntegerMod_int) and (<IntegerMod_int>x)._parent is R:
                    v[j] = (<IntegerMod_int>x).ivalue
                elif PY_TYPE_CHECK_EXACT(x, Integer):
                    if coerce:
                        v[j] = mpz_fdiv_ui((<Integer>x).value, p)
                    else:
                        v[j] = mpz_get_ui((<Integer>x).value)
                elif coerce:
                    v[j] = R(entries[k])
                else:
                    v[j] = <long>(entries[k])
                k = k + 1

    def __richcmp__(Matrix_modn_dense self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)

    def __hash__(self):
        """
        EXAMPLE::

            sage: B = random_matrix(GF(127),3,3)
            sage: B.set_immutable()
            sage: {B:0} # indirect doctest
            {[  9  75  94]
            [  4  57 112]
            [ 59  85  45]: 0}

        ::

            sage: M = random_matrix(GF(7), 10, 10)
            sage: M.set_immutable()
            sage: hash(M)
            143
            sage: MZ = M.change_ring(ZZ)
            sage: MZ.set_immutable()
            sage: hash(MZ)
            143
            sage: MS = M.sparse_matrix()
            sage: MS.set_immutable()
            sage: hash(MS)
            143

        TEST::

            sage: A = matrix(GF(2),2,0)
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()
            sage: hash(A)
            0
        """
        if self.is_mutable():
            raise TypeError("mutable matrices are unhashable")
        x = self.fetch('hash')
        if not x is None:
            return x

        cdef long _hash = 0
        cdef mod_int *_matrix
        cdef long n = 0
        cdef Py_ssize_t i, j

        if self._nrows == 0 or self._ncols == 0:
            return 0

        sig_on()
        for i from 0 <= i < self._nrows:
            _matrix = self._matrix[i]
            for j from 0 <= j < self._ncols:
                _hash ^= (n * _matrix[j])
                n+=1
        sig_off()

        if _hash == -1:
            return -2

        self.cache('hash', _hash)

        return _hash

    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value):
        """
        Set the (i,j) entry of self to the int value.
        """
        self._matrix[i][j] = value

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        """
        Set the (i,j) entry of self to the value, which is assumed
        to be of type IntegerMod_int.
        """
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
    # x * cdef _add_
    #   * cdef _mul_
    # x * cdef _matrix_times_matrix_
    # x * cdef _cmp_c_impl
    # x * __neg__
    #   * __invert__
    # x * __copy__
    #   * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################
    # def _pickle(self):
    # def _unpickle(self, data, int version):   # use version >= 0
    # cpdef ModuleElement _add_(self, ModuleElement right):
    # cdef _mul_(self, Matrix right):
    # cdef int _cmp_c_impl(self, Matrix right) except -2:
    # def __invert__(self):
    # def _multiply_classical(left, matrix.Matrix _right):
    # def _list(self):
    # def _dict(self):

    def _pickle(self):
        """
        Utility function for pickling.

        If the prime is small enough to fit in a byte, then it is stored as
        a contiguous string of bytes (to save space). Otherwise, memcpy is
        used to copy the raw data in the platforms native format. Any
        byte-swapping or word size difference is taken care of in
        unpickling (optimizing for unpickling on the same platform things
        were pickled on).

        The upcoming buffer protocol would be useful to not have to do any
        copying.

        EXAMPLES::

            sage: m = matrix(Integers(128), 3, 3, [ord(c) for c in "Hi there!"]); m
            [ 72 105  32]
            [116 104 101]
            [114 101  33]
            sage: m._pickle()
            ((1, ..., 'Hi there!'), 10)
        """
        cdef Py_ssize_t i, j
        cdef unsigned char* ss
        cdef unsigned char* row_ss
        cdef long word_size
        cdef mod_int *row_self

        if self.p <= 0xFF:
            word_size = sizeof(char)
        else:
            word_size = sizeof(mod_int)

        s = PyString_FromStringAndSize(NULL, word_size * self._nrows * self._ncols)
        ss = <unsigned char*><char *>s

        sig_on()
        if word_size == sizeof(char):
            for i from 0 <= i < self._nrows:
                row_self = self._matrix[i]
                row_ss = &ss[i*self._ncols]
                for j from 0<= j < self._ncols:
                    row_ss[j] = row_self[j]
        else:
            for i from 0 <= i < self._nrows:
                memcpy(&ss[i*sizeof(mod_int)*self._ncols], self._matrix[i], sizeof(mod_int) * self._ncols)
        sig_off()

        return (word_size, little_endian, s), 10


    def _unpickle(self, data, int version):
        r"""
        TESTS:

        Test for char-sized modulus::

            sage: A = random_matrix(GF(7), 5, 9)
            sage: data, version = A._pickle()
            sage: B = A.parent()(0)
            sage: B._unpickle(data, version)
            sage: B == A
            True

        And for larger modulus::

            sage: A = random_matrix(GF(1009), 51, 5)
            sage: data, version = A._pickle()
            sage: B = A.parent()(0)
            sage: B._unpickle(data, version)
            sage: B == A
            True

        Now test all the bit-packing options::

            sage: A = matrix(Integers(1000), 2, 2)
            sage: A._unpickle((1, True, '\x01\x02\xFF\x00'), 10)
            sage: A
            [  1   2]
            [255   0]

        ::

            sage: A = matrix(Integers(1000), 1, 2)
            sage: A._unpickle((4, True, '\x02\x01\x00\x00\x01\x00\x00\x00'), 10)
            sage: A
            [258   1]
            sage: A._unpickle((4, False, '\x00\x00\x02\x01\x00\x00\x01\x03'), 10)
            sage: A
            [513 259]
            sage: A._unpickle((8, True, '\x03\x01\x00\x00\x00\x00\x00\x00\x05\x00\x00\x00\x00\x00\x00\x00'), 10)
            sage: A
            [259   5]
            sage: A._unpickle((8, False, '\x00\x00\x00\x00\x00\x00\x02\x08\x00\x00\x00\x00\x00\x00\x01\x04'), 10)
            sage: A
            [520 260]

        Now make sure it works in context::

            sage: A = random_matrix(Integers(33), 31, 31)
            sage: loads(dumps(A)) == A
            True
            sage: A = random_matrix(Integers(3333), 31, 31)
            sage: loads(dumps(A)) == A
            True
        """

        if version < 10:
            return matrix_dense.Matrix_dense._unpickle(self, data, version)

        cdef Py_ssize_t i, j
        cdef unsigned char* ss
        cdef unsigned char* row_ss
        cdef long word_size
        cdef mod_int *row_self
        cdef bint little_endian_data
        cdef unsigned char* udata

        if version == 10:
            word_size, little_endian_data, s = data
            ss = <unsigned char*><char *>s

            sig_on()
            if word_size == sizeof(char):
                for i from 0 <= i < self._nrows:
                    row_self = self._matrix[i]
                    row_ss = &ss[i*self._ncols]
                    for j from 0<= j < self._ncols:
                        row_self[j] = row_ss[j]

            elif word_size == sizeof(mod_int) and little_endian == little_endian_data:
                for i from 0 <= i < self._nrows:
                    memcpy(self._matrix[i], &ss[i*sizeof(mod_int)*self._ncols], sizeof(mod_int) * self._ncols)

            # Note that mod_int is at least 32 bits, and never stores more than 32 bits of info
            elif little_endian_data:
                for i from 0 <= i < self._nrows:
                    row_self = self._matrix[i]
                    for j from 0<= j < self._ncols:
                        udata = &ss[(i*self._ncols+j)*word_size]
                        row_self[j] =  ((udata[0]) +
                                        (udata[1] << 8) +
                                        (udata[2] << 16) +
                                        (udata[3] << 24))

            else:
                for i from 0 <= i < self._nrows:
                    row_self = self._matrix[i]
                    for j from 0<= j < self._ncols:
                        udata = &ss[(i*self._ncols+j)*word_size]
                        row_self[j] =  ((udata[word_size-1]) +
                                        (udata[word_size-2] << 8) +
                                        (udata[word_size-3] << 16) +
                                        (udata[word_size-4] << 24))

            sig_off()

        else:
            raise RuntimeError, "unknown matrix version"


    cdef long _hash(self) except -1:
        """
        TESTS::

            sage: a = random_matrix(GF(11), 5, 5)
            sage: a.set_immutable()
            sage: hash(a) #random
            216
            sage: b = a.change_ring(ZZ)
            sage: b.set_immutable()
            sage: hash(b) == hash(a)
            True
        """
        x = self.fetch('hash')
        if not x is None: return x

        if not self._is_immutable:
            raise TypeError, "mutable matrices are unhashable"

        cdef Py_ssize_t i
        cdef long h = 0, n = 0

        sig_on()
        for i from 0 <= i < self._nrows:
            for j from 0<= j < self._ncols:
                h ^= n * self._matrix[i][j]
                n += 1
        sig_off()

        self.cache('hash', h)
        return h


    def __neg__(self):
        """
        EXAMPLES::

            sage: m = matrix(GF(19), 3, 3, range(9)); m
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: -m
            [ 0 18 17]
            [16 15 14]
            [13 12 11]
        """
        cdef Py_ssize_t i,j
        cdef Matrix_modn_dense M
        cdef mod_int p = self.p

        M = Matrix_modn_dense.__new__(Matrix_modn_dense, self._parent,None,None,None)
        M.p = p

        sig_on()
        cdef mod_int *row_self, *row_ans
        for i from 0 <= i < self._nrows:
            row_self = self._matrix[i]
            row_ans = M._matrix[i]
            for j from 0<= j < self._ncols:
                if row_self[j]:
                    row_ans[j] = p - row_self[j]
                else:
                    row_ans[j] = 0
        sig_off()
        return M


    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        EXAMPLES::

            sage: a = random_matrix(Integers(60), 400, 500)
            sage: 3*a + 9*a == 12*a
            True
        """
        return self._rmul_(right)

    cpdef ModuleElement _rmul_(self, RingElement left):
        """
        EXAMPLES::

            sage: a = matrix(GF(101), 3, 3, range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a * 5
            [ 0  5 10]
            [15 20 25]
            [30 35 40]
            sage: a * 50
            [  0  50 100]
            [ 49  99  48]
            [ 98  47  97]
        """
        cdef Py_ssize_t i,j
        cdef Matrix_modn_dense M
        cdef mod_int p = self.p
        cdef mod_int a = left

        M = Matrix_modn_dense.__new__(Matrix_modn_dense, self._parent,None,None,None)
        M.p = p

        sig_on()
        cdef mod_int *row_self, *row_ans
        for i from 0 <= i < self._nrows:
            row_self = self._matrix[i]
            row_ans = M._matrix[i]
            for j from 0<= j < self._ncols:
                row_ans[j] = (a*row_self[j]) % p
        sig_off()
        return M

    def __copy__(self):
        cdef Matrix_modn_dense A
        A = Matrix_modn_dense.__new__(Matrix_modn_dense, self._parent,
                                      0, 0, 0)
        memcpy(A._entries, self._entries, sizeof(mod_int)*self._nrows*self._ncols)
        A.p = self.p
        A.gather = self.gather
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A


    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Add two dense matrices over Z/nZ

        EXAMPLES::

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
        cdef mod_int k, p
        cdef Matrix_modn_dense M

        M = Matrix_modn_dense.__new__(Matrix_modn_dense, self._parent,None,None,None)
        Matrix.__init__(M, self._parent)
        p = self.p

        sig_on()
        cdef mod_int *row_self, *row_right, *row_ans
        for i from 0 <= i < self._nrows:
            row_self = self._matrix[i]
            row_right = (<Matrix_modn_dense> right)._matrix[i]
            row_ans = M._matrix[i]
            for j from 0<= j < self._ncols:
                k = row_self[j] + row_right[j]
                row_ans[j] = k - (k >= p) * p
        sig_off()
        return M


    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Subtract two dense matrices over Z/nZ

        EXAMPLES::

            sage: a = matrix(GF(11), 3, 3, range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a - 4
            [7 1 2]
            [3 0 5]
            [6 7 4]
            sage: a - matrix(GF(11), 3, 3, range(1, 19, 2))
            [10  9  8]
            [ 7  6  5]
            [ 4  3  2]
        """

        cdef Py_ssize_t i,j
        cdef mod_int k, p
        cdef Matrix_modn_dense M

        M = Matrix_modn_dense.__new__(Matrix_modn_dense, self._parent,None,None,None)
        Matrix.__init__(M, self._parent)
        p = self.p

        sig_on()
        cdef mod_int *row_self, *row_right, *row_ans
        for i from 0 <= i < self._nrows:
            row_self = self._matrix[i]
            row_right = (<Matrix_modn_dense> right)._matrix[i]
            row_ans = M._matrix[i]
            for j from 0<= j < self._ncols:
                k = p + row_self[j] - row_right[j]
                row_ans[j] = k - (k >= p) * p
        sig_off()
        return M


    cdef int _cmp_c_impl(self, Element right) except -2:
        """
        Compare two dense matrices over Z/nZ

        EXAMPLES::

            sage: a = matrix(GF(17), 4, range(3, 83, 5)); a
            [ 3  8 13  1]
            [ 6 11 16  4]
            [ 9 14  2  7]
            [12  0  5 10]
            sage: a == a
            True
            sage: b = a - 3; b
            [ 0  8 13  1]
            [ 6  8 16  4]
            [ 9 14 16  7]
            [12  0  5  7]
            sage: b < a
            True
            sage: b > a
            False
            sage: b == a
            False
            sage: b + 3 == a
            True
        """

        cdef Py_ssize_t i, j
        cdef int cmp

        sig_on()
        cdef mod_int *row_self, *row_right
        for i from 0 <= i < self._nrows:
            row_self = self._matrix[i]
            row_right = (<Matrix_modn_dense> right)._matrix[i]
            for j from 0 <= j < self._ncols:
                if row_self[j] < row_right[j]:
                    sig_off()
                    return -1
                elif row_self[j] > row_right[j]:
                    sig_off()
                    return 1
            #cmp = memcmp(row_self, row_right, sizeof(mod_int)*self._ncols)
            #if cmp:
            #    return cmp
        sig_off()
        return 0


    cdef Matrix _matrix_times_matrix_(self, Matrix right):
        if get_verbose() >= 2:
            verbose('mod-p multiply of %s x %s matrix by %s x %s matrix modulo %s'%(
                self._nrows, self._ncols, right._nrows, right._ncols, self.p))
        if self._will_use_strassen(right):
            return self._multiply_strassen(right)
        else:
            return self._multiply_classical(right)

    def _multiply_classical(left, right):
        return left._multiply_strassen(right, left._ncols + left._nrows)

    cdef Vector _vector_times_matrix_(self, Vector v):
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
    #  x * cdef _sub_
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    #        - all row/column operations, but optimized
    # x      - echelon form in place
    #        - Hessenberg forms of matrices
    ########################################################################


    def charpoly(self, var='x', algorithm='generic'):
        """
        Returns the characteristic polynomial of self.

        INPUT:

        -  ``var`` - a variable name
        -  ``algorithm`` - 'generic' (default)

        EXAMPLES::

            sage: A = Mat(GF(7),3,3)(range(3)*3)
            sage: A.charpoly()
            x^3 + 4*x^2

            sage: A = Mat(Integers(6),3,3)(range(9))
            sage: A.charpoly()
            x^3
        """
        if algorithm == 'generic':
            g = matrix_dense.Matrix_dense.charpoly(self, var)
        else:
            raise ValueError, "no algorithm '%s'"%algorithm
        self.cache('charpoly_%s_%s'%(algorithm, var), g)
        return g

    def minpoly(self, var='x', algorithm='generic', proof=None):
        """
        Returns the minimal polynomial of self.

        INPUT:

           - ``var`` - a variable name
           - ``algorithm`` - 'generic' (default)

           - ``proof`` -- (default: True); whether to provably return
             the true minimal polynomial; if False, we only guarantee
             to return a divisor of the minimal polynomial.  There are
             also certainly cases where the computed results is
             frequently not exactly equal to the minimal polynomial
             (but is instead merely a divisor of it).

        WARNING: If proof=True, minpoly is insanely slow compared to
        proof=False.

        EXAMPLES::

            sage: R.<x>=GF(3)[]
            sage: A = matrix(GF(3),2,[0,0,1,2])
            sage: A.minpoly()
            x^2 + x

        Minpoly with proof=False is *dramatically* ("millions" of
        times!)  faster than minpoly with proof=True.  This matters
        since proof=True is the default, unless you first type
        ''proof.linear_algebra(False)''.::

            sage: A.minpoly(proof=False) in [x, x+1, x^2+x]
            True
        """

        proof = get_proof_flag(proof, "linear_algebra")

        if algorithm == 'generic':
            raise NotImplementedError, "minimal polynomials are not implemented for Z/nZ"
        else:
            raise ValueError, "no algorithm '%s'"%algorithm
        self.cache('minpoly_%s_%s'%(algorithm, var), g)
        return g

    def echelonize(self, algorithm="gauss", **kwds):
        """
        Puts self in row echelon form.

        INPUT:

        -  ``self`` - a mutable matrix
        - ``algorithm``- ``'gauss'`` - uses a custom slower `O(n^3)`
           Gauss elimination implemented in Sage.
        -  ``**kwds`` - these are all ignored


        OUTPUT:

        -  self is put in reduced row echelon form.
        -  the rank of self is computed and cached
        - the pivot columns of self are computed and cached.
        - the fact that self is now in echelon form is recorded and
          cached so future calls to echelonize return immediately.

        EXAMPLES::

            sage: a = matrix(GF(97),3,4,range(12))
            sage: a.echelonize(); a
            [ 1  0 96 95]
            [ 0  1  2  3]
            [ 0  0  0  0]
            sage: a.pivots()
            (0, 1)
        """
        x = self.fetch('in_echelon_form')
        if not x is None: return  # already known to be in echelon form
        if not self.base_ring().is_field():
            #self._echelon_mod_n ()
            raise NotImplementedError, "Echelon form not implemented over '%s'."%self.base_ring()

        self.check_mutability()
        self.clear_cache()

        if algorithm == 'gauss':
            self._echelon_in_place_classical()
        else:
            raise ValueError, "algorithm '%s' not known"%algorithm

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
            sig_check()
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
        self.cache('pivots', tuple(pivots))
        self.cache('in_echelon_form',True)

    cdef xgcd_eliminate (self, mod_int * row1, mod_int* row2, Py_ssize_t start_col):
        """
        Reduces row1 and row2 by a unimodular transformation using the xgcd
        relation between their first coefficients a and b.

        INPUT:


        -   self: a mutable matrix

        -````- row1, row2: the two rows to be transformed
           (within self)

        -```` - start_col: the column of the pivots in row1
           and row2. It is assumed that all entries before start_col in row1
           and row2 are zero.


        OUTPUT:


        - g: the gcd of the first elements of row1 and
          row2 at column start_col

        - put row1 = s \* row1 + t \* row2 row2 = w \*
          row1 + t \* row2 where g = sa + tb
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


## This is code by William Stein and/or Clement Pernet from SD7. Unfortunately I (W.S.)
## think it is still buggy, since it is so painful to implement with
## unsigned ints. Code to do basically the same thing is in
## matrix_integer_dense, by Burcin Erocal.
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

##             sig_check()
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
        EXAMPLES::

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

        Rescaling need not include the entire row.

        ::

            sage: b.rescale_col(0,2,1); b
            [1 5 1]
            [2 5 1]
            [2 5 1]

        Bounds are checked.

        ::

            sage: b.rescale_col(3,2)
            Traceback (most recent call last):
            ...
            IndexError: matrix column index out of range

        Rescaling by a negative number::

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
        """
        EXAMPLES::

            sage: A = matrix(Integers(8), 2,[1,2,3,4])
            sage: A.swap_rows(0,1)
            sage: A
            [3 4]
            [1 2]
        """
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

        sig_on()
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
        sig_off()
        self.cache('in_hessenberg_form',True)

    def _charpoly_hessenberg(self, var):
        """
        Transforms self in place to its Hessenberg form then computes and
        returns the coefficients of the characteristic polynomial of this
        matrix.

        INPUT:


        -  ``var`` - name of the indeterminate of the
           charpoly.


        The characteristic polynomial is represented as a vector of ints,
        where the constant term of the characteristic polynomial is the 0th
        coefficient of the vector.
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
        H = self.__copy__()
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

        EXAMPLES::

            sage: m = matrix(GF(7),5,range(25))
            sage: m.rank()
            2

        Rank is not implemented over the integers modulo a composite yet.::

            sage: m = matrix(Integers(4), 2, [2,2,2,2])
            sage: m.rank()
            Traceback (most recent call last):
            ...
            NotImplementedError: Echelon form not implemented over 'Ring of integers modulo 4'.
        """
        return matrix_dense.Matrix_dense.rank(self)

    def determinant(self):
        """
        Return the determinant of this matrix.

        EXAMPLES::

            sage: m = matrix(GF(101),5,range(25))
            sage: m.det()
            0

        ::

            sage: m = matrix(Integers(4), 2, [2,2,2,2])
            sage: m.det()
            0

        TESTS::

            sage: m = random_matrix(GF(3), 3, 4)
            sage: m.determinant()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
        """
        if self._nrows != self._ncols:
            raise ValueError, "self must be a square matrix"
        if self._nrows == 0:
            return self._coerce_element(1)

            return matrix_dense.Matrix_dense.determinant(self)

    def randomize(self, density=1, nonzero=False):
        """
        Randomize ``density`` proportion of the entries of this matrix,
        leaving the rest unchanged.

        INPUT:

        -  ``density`` - Integer; proportion (roughly) to be considered for
           changes
        -  ``nonzero`` - Bool (default: ``False``); whether the new entries
           are forced to be non-zero

        OUTPUT:

        -  None, the matrix is modified in-space

        EXAMPLES::

            sage: A = matrix(GF(5), 5, 5, 0)
            sage: A.randomize(0.5); A
            [0 0 0 2 0]
            [0 3 0 0 2]
            [4 0 0 0 0]
            [4 0 0 0 0]
            [0 1 0 0 0]
            sage: A.randomize(); A
            [3 3 2 1 2]
            [4 3 3 2 2]
            [0 3 3 3 3]
            [3 3 2 2 4]
            [2 2 2 1 4]
        """
        density = float(density)
        if density <= 0:
            return
        if density > 1:
            density = float(1)

        self.check_mutability()
        self.clear_cache()

        cdef randstate rstate = current_randstate()
        cdef int nc
        cdef long pm1

        if not nonzero:
            # Original code, before adding the ``nonzero`` option.
            if density == 1:
                for i from 0 <= i < self._nrows*self._ncols:
                    self._entries[i] = rstate.c_random() % self.p
            else:
                nc = self._ncols
                num_per_row = int(density * nc)
                sig_on()
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < num_per_row:
                        k = rstate.c_random() % nc
                        self._matrix[i][k] = rstate.c_random() % self.p
                sig_off()
        else:
            # New code, to implement the ``nonzero`` option.
            pm1 = self.p - 1
            if density == 1:
                for i from 0 <= i < self._nrows*self._ncols:
                    self._entries[i] = (rstate.c_random() % pm1) + 1
            else:
                nc = self._ncols
                num_per_row = int(density * nc)
                sig_on()
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < num_per_row:
                        k = rstate.c_random() % nc
                        self._matrix[i][k] = (rstate.c_random() % pm1) + 1
                sig_off()

    cdef int _strassen_default_cutoff(self, matrix0.Matrix right) except -2:
        # TODO: lots of testing
        return 100

    cpdef matrix_window(self, Py_ssize_t row=0, Py_ssize_t col=0,
                       Py_ssize_t nrows=-1, Py_ssize_t ncols=-1,
                       bint check=1):
        """
        Return the requested matrix window.

        EXAMPLES::

            sage: a = matrix(GF(7),3,range(9)); a
            [0 1 2]
            [3 4 5]
            [6 0 1]
            sage: type(a)
            <type 'sage.matrix.matrix_modn_dense_float.Matrix_modn_dense_float'>

        We test the optional check flag.

        ::

            sage: matrix(GF(7),[1]).matrix_window(0,1,1,1)
            Traceback (most recent call last):
            ...
            IndexError: matrix window index out of range
            sage: matrix(GF(7),[1]).matrix_window(0,1,1,1, check=False)
            Matrix window of size 1 x 1 at (0,1):
            [1]
        """
        if nrows == -1:
            nrows = self._nrows - row
            ncols = self._ncols - col
        if check and (row < 0 or col < 0 or row + nrows > self._nrows or \
           col + ncols > self._ncols):
            raise IndexError, "matrix window index out of range"
        return matrix_window_modn_dense.MatrixWindow_modn_dense(self, row, col, nrows, ncols)

    def _magma_init_(self, magma):
        """
        Returns a string representation of self in Magma form.

        INPUT:


        -  ``magma`` - a Magma session


        OUTPUT: string

        EXAMPLES::

            sage: a = matrix(GF(389),2,2,[1..4])
            sage: magma(a)                                         # optional - magma
            [  1   2]
            [  3   4]
            sage: a._magma_init_(magma)                            # optional - magma
            'Matrix(GF(389),2,2,StringToIntegerSequence("1 2 3 4"))'

        A consistency check::

            sage: a = random_matrix(GF(13),50); b = random_matrix(GF(13),50)
            sage: magma(a*b) == magma(a)*magma(b)                  # optional - magma
            True
        """
        s = self.base_ring()._magma_init_(magma)
        return 'Matrix(%s,%s,%s,StringToIntegerSequence("%s"))'%(
            s, self._nrows, self._ncols, self._export_as_string())

    cpdef _export_as_string(self):
        """
        Return space separated string of the entries in this matrix.

        EXAMPLES::

            sage: w = matrix(GF(997),2,3,[1,2,5,-3,8,2]); w
            [  1   2   5]
            [994   8   2]
            sage: w._export_as_string()
            '1 2 5 994 8 2'
        """
        cdef int ndigits = len(str(self.p))

        cdef Py_ssize_t i, n
        cdef char *s, *t

        if self._nrows == 0 or self._ncols == 0:
            data = ''
        else:
            n = self._nrows*self._ncols*(ndigits + 1) + 2  # spaces between each number plus trailing null
            s = <char*> sage_malloc(n * sizeof(char))
            t = s
            sig_on()
            for i in range(self._nrows * self._ncols):
                sprintf(t, "%d ", self._entries[i])
                t += strlen(t)
            sig_off()
            data = str(s)[:-1]
            sage_free(s)
        return data

    def _list(self):
        """
        Return list of elements of self.  This method is called by self.list().

        EXAMPLES::

            sage: w = matrix(GF(19), 2, 3, [1..6])
            sage: w.list()
            [1, 2, 3, 4, 5, 6]
            sage: w._list()
            [1, 2, 3, 4, 5, 6]
            sage: w.list()[0].parent()
            Finite Field of size 19

        TESTS::

            sage: w = random_matrix(GF(3),100)
            sage: w.parent()(w.list()) == w
            True
        """
        cdef Py_ssize_t i, j
        F = self.base_ring()
        entries = []
        for i from 0 <= i < self._nrows:
            for j from 0<= j < self._ncols:
                entries.append(F(self._matrix[i][j]))
        return entries

    def lift(self):
        """
        Return the lift of this matrix to the integers.

        EXAMPLES::

            sage: a = matrix(GF(7),2,3,[1..6])
            sage: a.lift()
            [1 2 3]
            [4 5 6]
            sage: a.lift().parent()
            Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

        Subdivisions are preserved when lifting::

            sage: a.subdivide([], [1,1]); a
            [1||2 3]
            [4||5 6]
            sage: a.lift()
            [1||2 3]
            [4||5 6]
        """
        cdef Py_ssize_t i, j

        cdef Matrix_integer_dense L
        L = Matrix_integer_dense._new_uninitialized_matrix(Matrix_integer_dense,self._nrows,self._ncols)
        cdef mod_int* A_row
        for i from 0 <= i < self._nrows:
            A_row = self._matrix[i]
            for j from 0 <= j < self._ncols:
                fmpz_set_si(fmpz_mat_entry(L._matrix,i,j),A_row[j])
        L._initialized = 1
        L.subdivide(self.subdivisions())
        return L


    def _matrices_from_rows(self, Py_ssize_t nrows, Py_ssize_t ncols):
        """
        Make a list of matrix from the rows of this matrix.  This is a
        fairly technical function which is used internally, e.g., by
        the cyclotomic field linear algebra code.

        INPUT:

        - ``nrows, ncols`` - integers

        OUTPUT:

        - ``list`` - a list of matrices
        """
        if nrows * ncols != self._ncols:
            raise ValueError, "nrows * ncols must equal self's number of columns"

        from matrix_space import MatrixSpace
        F = self.base_ring()
        MS = MatrixSpace(F, nrows, ncols)

        cdef Matrix_modn_dense M
        cdef Py_ssize_t i
        cdef Py_ssize_t n = nrows * ncols
        ans = []
        for i from 0 <= i < self._nrows:
            # Quickly construct a new mod-p matrix
            M = Matrix_modn_dense.__new__(Matrix_modn_dense, MS, 0,0,0)
            M.p = self.p
            M.gather = self.gather
            # Set the entries
            memcpy(M._entries, self._entries+i*n, sizeof(mod_int)*n)
            ans.append(M)
        return ans

    _matrix_from_rows_of_matrices = staticmethod(__matrix_from_rows_of_matrices)
