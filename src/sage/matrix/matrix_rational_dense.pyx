"""
Dense matrices over the rational field.

EXAMPLES:
We create a 3x3 matrix with rational entries and do some
operations with it.

    sage: a = matrix(QQ, 3,3, [1,2/3, -4/5, 1,1,1, 8,2, -3/19]); a
    [    1   2/3  -4/5]
    [    1     1     1]
    [    8     2 -3/19]
    sage: a.det()
    2303/285
    sage: a.charpoly()
    x^3 - 35/19*x^2 + 1259/285*x - 2303/285
    sage: b = a^(-1); b
    [ -615/2303  -426/2303   418/2303]
    [ 2325/2303  1779/2303  -513/2303]
    [-1710/2303   950/2303    95/2303]
    sage: b.det()
    285/2303
    sage: a == b
    False
    sage: a < b
    False
    sage: b < a
    True
    sage: a > b
    True
    sage: a*b
    [1 0 0]
    [0 1 0]
    [0 0 1]
"""

##############################################################################
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.modules.vector_rational_dense cimport Vector_rational_dense

include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"
include "../ext/gmp.pxi"
include "../ext/random.pxi"

cimport sage.structure.element

from sage.structure.sequence import Sequence
from sage.rings.rational cimport Rational
from matrix cimport Matrix
from matrix_integer_dense cimport Matrix_integer_dense
from matrix_integer_dense import _lift_crt
import sage.structure.coerce
from sage.structure.element cimport ModuleElement, RingElement, Element, Vector
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.finite_field import GF
from sage.rings.rational_field import QQ

import sage.ext.multi_modular
from matrix2 import cmp_pivots

from sage.misc.misc import verbose, get_verbose

cdef class Matrix_rational_dense(matrix_dense.Matrix_dense):

    ########################################################################
    # LEVEL 1 functionality
    # x * __new__
    # x * __dealloc__
    # x * __init__
    # x * set_unsafe
    # x * get_unsafe
    # x * cdef _pickle
    # x * cdef _unpickle
    ########################################################################
    def __new__(self, parent, entries, copy, coerce):
        """
        Create and allocate memory for the matrix.

        Unlike over matrix_integer_dense, mpq_init() is called (as there is no mpq_init_set function).

        INPUT:
            parent, entries, coerce, copy -- as for __init__.

        EXAMPLES:
            sage: from sage.matrix.matrix_rational_dense import Matrix_rational_dense
            sage: a = Matrix_rational_dense.__new__(Matrix_rational_dense, Mat(ZZ,3), 0,0,0)
            sage: type(a)
            <type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>

        WARNING: This is for internal use only, or if you really know what you're doing.
        """
        matrix_dense.Matrix_dense.__init__(self, parent)

        cdef Py_ssize_t i, k

        _sig_on
        self._entries = <mpq_t *> sage_malloc(sizeof(mpq_t)*(self._nrows * self._ncols))
        _sig_off
        if self._entries == NULL:
            raise MemoryError, "out of memory allocating a matrix"

        self._matrix =  <mpq_t **> sage_malloc(sizeof(mpq_t*) * self._nrows)
        if self._matrix == NULL:
            raise MemoryError, "out of memory allocating a matrix"

        # store pointers to the starts of the rows
        _sig_on
        k = 0
        for i from 0 <= i < self._nrows:
            self._matrix[i] = self._entries + k
            k = k + self._ncols

        for i from 0 <= i < self._nrows * self._ncols:
            mpq_init(self._entries[i])
        _sig_off

    def  __dealloc__(self):
        cdef Py_ssize_t i
        for i from 0 <= i < self._nrows * self._ncols:
            mpq_clear(self._entries[i])
        sage_free(self._entries)
        sage_free(self._matrix)

    def __init__(self, parent, entries=0, coerce=True, copy=True):

        cdef Py_ssize_t i
        cdef Rational z

        if isinstance(entries, list):
            if len(entries) != self._nrows * self._ncols:
                raise TypeError, "entries has the wrong length"

            _sig_on
            if coerce:
                for i from 0 <= i < self._nrows * self._ncols:
                    # TODO: Should use an unsafe un-bounds-checked array access here.
                    z = Rational(entries[i])
                    mpq_set(self._entries[i], z.value)
            else:
                for i from 0 <= i < self._nrows * self._ncols:
                    # TODO: Should use an unsafe un-bounds-checked array access here.
                    mpq_set(self._entries[i], (<Rational> entries[i]).value)
            _sig_off

        else:
            # is it a scalar?
            try:
                # Try to coerce entries to a scalar (an integer)
                z = Rational(entries)
                is_list = False
            except TypeError:
                raise TypeError, "entries must be coercible to a list or integer"

            if not z.is_zero():
                if self._nrows != self._ncols:
                    raise TypeError, "nonzero scalar matrix must be square"
                for i from 0 <= i < self._nrows:
                    mpq_set(self._entries[i*self._ncols+i], z.value)


    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        cdef Rational y
        y = value
        mpq_set(self._matrix[i][j], y.value)


    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef Rational x
        x = PY_NEW(Rational)
        mpq_set(x.value, self._matrix[i][j])
        return x

    def _pickle(self):
        return self._pickle_version0(), 0

    def _unpickle(self, data, int version):
        if version == 0:
            self._unpickle_version0(data)
        else:
            raise RuntimeError, "unknown matrix version (=%s)"%version

    cdef _pickle_version0(self):
        cdef Py_ssize_t i, j, len_so_far, m, n
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
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    m = mpz_sizeinbase (mpq_numref(self._matrix[i][j]), 32) + \
                        mpz_sizeinbase (mpq_denref(self._matrix[i][j]), 32) + 3
                    if len_so_far + m + 1 >= n:
                        # copy to new string with double the size
                        n = 2*n + m + 1
                        tmp = <char*> sage_malloc(n * sizeof(char))
                        strcpy(tmp, s)
                        sage_free(s)
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
            data = str(s)[:-1]
            free(s)
        return data

    cdef _unpickle_version0(self, data):
        cdef Py_ssize_t i, n
        data = data.split()
        n = self._nrows * self._ncols
        if len(data) != n:
            raise RuntimeError, "invalid pickle data."
        for i from 0 <= i < n:
            s = data[i]
            if mpq_set_str(self._entries[i], s, 32):
                raise RuntimeError, "invalid pickle data"

    def __richcmp__(Matrix self, right, int op):
        return self._richcmp(right, op)
    def __hash__(self):
        return self._hash()

    ########################################################################
    # LEVEL 2 functionality
    # x * cdef _add_c_impl
    # x * cdef _mul_c_impl
    # x * cdef _vector_times_matrix_c_impl
    # x * cdef _cmp_c_impl
    # x * __neg__
    #   * __invert__
    # x * __copy__
    # x * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        """
        EXAMPLES:
            sage: a = matrix(QQ,2,range(6))
            sage: (3/4) * a
            [   0  3/4  3/2]
            [ 9/4    3 15/4]
        """
        cdef Py_ssize_t i
        cdef Rational _x
        _x = Rational(right)
        cdef Matrix_rational_dense M
        M = Matrix_rational_dense.__new__(Matrix_rational_dense, self._parent, None, None, None)
        for i from 0 <= i < self._nrows * self._ncols:
            mpq_mul(M._entries[i], self._entries[i], _x.value)
        return M

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        """
        Add two dense matrices over QQ.

        EXAMPLES:
            sage: a = MatrixSpace(QQ,3)(range(9))
            sage: b = MatrixSpace(QQ,3)([1/n for n in range(1,10)])
            sage: a+b
            [   1  3/2  7/3]
            [13/4 21/5 31/6]
            [43/7 57/8 73/9]
            sage: b.swap_rows(1,2)
            sage: #a+b

        """
        cdef Py_ssize_t i, j
        cdef Matrix_rational_dense M
        M = Matrix_rational_dense.__new__(Matrix_rational_dense, self._parent, None, None, None)

        cdef mpq_t *M_row
        cdef mpq_t *self_row
        cdef mpq_t *right_row
        _sig_on
        for i from 0 <= i < self._nrows:
            M_row = M._matrix[i]
            self_row = self._matrix[i]
            right_row = (<Matrix_rational_dense>right)._matrix[i]
            for j from 0 <= j < self._ncols:
                mpq_add(M_row[0], self_row[0], right_row[0])
                M_row = M_row + 1
                self_row = self_row + 1
                right_row = right_row + 1
        _sig_off
        return M

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Subtract two dense matrices over QQ.

        EXAMPLES:
            sage: a = MatrixSpace(QQ,3)(range(9))
            sage: b = MatrixSpace(QQ,3)([1/n for n in range(1,10)])
            sage: a-b
            [  -1  1/2  5/3]
            [11/4 19/5 29/6]
            [41/7 55/8 71/9]
        """
        cdef Py_ssize_t i, j
        cdef Matrix_rational_dense M
        M = Matrix_rational_dense.__new__(Matrix_rational_dense, self._parent, None, None, None)

        cdef mpq_t *M_row
        cdef mpq_t *self_row
        cdef mpq_t *right_row
        _sig_on
        for i from 0 <= i < self._nrows:
            M_row = M._matrix[i]
            self_row = self._matrix[i]
            right_row = (<Matrix_rational_dense>right)._matrix[i]
            for j from 0 <= j < self._ncols:
                mpq_sub(M_row[0], self_row[0], right_row[0])
                M_row = M_row + 1
                self_row = self_row + 1
                right_row = right_row + 1
        _sig_off
        return M

    cdef int _cmp_c_impl(self, Element right) except -2:
        cdef mpq_t *a, *b
        cdef Py_ssize_t i, j
        cdef int k
        for i from 0 <= i < self._nrows:
            a = self._matrix[i]
            b = (<Matrix_rational_dense>right)._matrix[i]
            for j from 0 <= j < self._ncols:
                k = mpq_cmp(a[j], b[j])
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
            sage: B = matrix(QQ,2, [1,2,3,4])
            sage: V = QQ^2
            sage: w = V([-1,5/2])
            sage: w*B
            (13/2, 8)
        """
        cdef Vector_rational_dense w, ans
        cdef Py_ssize_t i, j
        cdef mpq_t x

        M = self._row_ambient_module()
        w = <Vector_rational_dense> v
        ans = M.zero_vector()

        mpq_init(x)
        mpq_init(y)
        for i from 0 <= i < self._ncols:
            mpq_set_si(x, 0,1)
            for j from 0 <= j < self._nrows:
                mpq_mul(y, w._entries[j], self._matrix[j][i])
                mpq_add(x, x, y)
            mpq_set(ans._entries[i], x)
        mpq_clear(x)
        mpq_clear(y)
        return ans


    def __neg__(self):
        """
        Negate a matrix over QQ.

        EXAMPLES:
            sage: a = MatrixSpace(QQ,3)([1/n for n in range(1,10)])
            sage: -a
            [  -1 -1/2 -1/3]
            [-1/4 -1/5 -1/6]
            [-1/7 -1/8 -1/9]
        """
        cdef Py_ssize_t i, j
        cdef Matrix_rational_dense M
        M = Matrix_rational_dense.__new__(Matrix_rational_dense, self._parent, None, None, None)

        cdef mpq_t *M_row
        cdef mpq_t *self_row
        _sig_on
        for i from 0 <= i < self._nrows:
            M_row = M._matrix[i]
            self_row = self._matrix[i]
            for j from 0 <= j < self._ncols:
                mpq_neg(M_row[0], self_row[0])
                M_row = M_row + 1
                self_row = self_row + 1
        _sig_off
        return M

    def __copy__(self):
        """
        Copy a matrix over QQ.

        EXAMPLES:
            sage: a = MatrixSpace(QQ,3)([1/n for n in range(1,10)])
            sage: -a
            [  -1 -1/2 -1/3]
            [-1/4 -1/5 -1/6]
            [-1/7 -1/8 -1/9]
        """
        cdef Py_ssize_t i, j
        cdef Matrix_rational_dense M
        M = Matrix_rational_dense.__new__(Matrix_rational_dense, self._parent, None, None, None)

        cdef mpq_t *M_row
        cdef mpq_t *self_row
        _sig_on
        for i from 0 <= i < self._nrows:
            M_row = M._matrix[i]
            self_row = self._matrix[i]
            for j from 0 <= j < self._ncols:
                mpq_set(M_row[0], self_row[0])
                M_row = M_row + 1
                self_row = self_row + 1
        _sig_off
        return M



    # cdef _mul_c_impl(self, Matrix right):
    # cdef int _cmp_c_impl(self, Matrix right) except -2:
    # def __invert__(self):
    # def _multiply_classical(left, matrix.Matrix _right):
    # def _list(self):
    # def _dict(self):


    ########################################################################
    # LEVEL 3 functionality (Optional)
    # x * cdef _sub_c_impl
    #   * __deepcopy__
    #   * __invert__
    #   * Matrix windows -- only if you need strassen for that base
    #   * Other functions (list them here):
    # x * denom(self):
    # x * mpz_denom(self, mpz_t d):
    # x * _clear_denom(self):
    # x * _multiply_multi_modular(self, Matrix_rational_dense right):
    # o * echelon_modular(self, height_guess=None):
    ########################################################################

    def __invert__(self):
        """
        OUTPUT:
           -- the inverse of self

        If self is not invertible, a ZeroDivisionError is raised.

        EXAMPLES:
            sage: a = matrix(QQ,3,range(9))
            sage: a^(-1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular
        """
        A, denom = self._clear_denom()
        B, d = A._invert_iml()
        return (denom/d)*B

    def determinant(self):
        """
        Return the determinant of this matrix.

        ALGORITHM: Clear denominators and call the integer determinant function.

        EXAMPLES:
            sage: m = matrix(QQ,3,[1,2/3,4/5, 2,2,2, 5,3,2/5])
            sage: m.determinant()
            -34/15
            sage: m.charpoly()
            x^3 - 17/5*x^2 - 122/15*x + 34/15
        """
        det = self.fetch('det')
        if not det is None: return det

        A, denom = self._clear_denom()
        det = Rational(A.determinant())
        if denom != 1:
            det = det / (denom**self.nrows())
        self.cache('det', det)
        return det


    def denom(self):
        """
        Return the denominator of this matrix.

        OUTPUT:
            -- SAGE Integer

        EXAMPLES:
            sage: b = matrix(QQ,2,range(6)); b[0,0]=-5007/293; b
            [-5007/293         1         2]
            [        3         4         5]
            sage: b.denom()
            293
        """
        cdef Integer z
        z = PY_NEW(Integer)
        self.mpz_denom(z.value)
        return z

    cdef int mpz_denom(self, mpz_t d) except -1:
        mpz_set_si(d,1)
        cdef Py_ssize_t i, j
        cdef mpq_t *self_row
        _sig_on
        for i from 0 <= i < self._nrows:
            self_row = self._matrix[i]
            for j from 0 <= j < self._ncols:
                mpz_lcm(d, d, mpq_denref(self_row[0]))
                self_row = self_row + 1
        _sig_off
        return 0

    def _clear_denom(self):
        """
        INPUT:
            self -- a matrix
        OUTPUT:
            D*self, D

        The product is a matrix over ZZ
        """
        cdef Integer D
        cdef Py_ssize_t i, j
        cdef Matrix_integer_dense A
        cdef mpq_t *self_row
        cdef mpz_t *A_row
        D = <Integer>PY_NEW(Integer)
        self.mpz_denom(D.value)
        MZ = sage.matrix.matrix_space.MatrixSpace(ZZ, self._nrows, self._ncols, sparse=self.is_sparse())
        A = Matrix_integer_dense.__new__(Matrix_integer_dense, MZ, 0, 0, 0)
        _sig_on
        for i from 0 <= i < self._nrows:
            A_row = A._matrix[i]
            self_row = self._matrix[i]
            for j from 0 <= j < self._ncols:
                mpz_init(A_row[0])
                mpz_divexact(A_row[0], D.value, mpq_denref(self_row[0]))
                mpz_mul(A_row[0], A_row[0], mpq_numref(self_row[0]))
                A_row = A_row + 1
                self_row = self_row + 1
        _sig_off
        A._initialized = 1
        return A, D

    def charpoly(self, var='x', algorithm='linbox'):
        """
        Return the characteristic polynomial of this matrix.

        INPUT:
            var -- 'x' (string)
            algorithm -- 'linbox' (default)
                         'generic'

        OUTPUT:
            a polynomial over the rational numbers.

        EXAMPLES:
            sage: a = matrix(QQ, 3, [4/3, 2/5, 1/5, 4, -3/2, 0, 0, -2/3, 3/4])
            sage: f = a.charpoly(); f
            x^3 - 7/12*x^2 - 149/40*x + 97/30
            sage: f(a)
            [0 0 0]
            [0 0 0]
            [0 0 0]
        """
        key = 'charpoly_%s_%s'%(algorithm, var)
        x = self.fetch(key)
        if x: return x

        if algorithm == 'linbox':
            A, denom = self._clear_denom()
            f = A.charpoly(var, algorithm='linbox')
            x = f.parent().gen()
            g = f(x * denom) * (1 / (denom**f.degree()))
        elif algorithm == 'generic':
            g = matrix_dense.Matrix_dense.charpoly(self, var)
        else:
            raise ValueError, "no algorithm '%s'"%algorithm

        self.cache(key, g)
        return g

    def minpoly(self, var='x', algorithm='linbox'):
        """
        Return the minimal polynomial of this matrix.

        INPUT:
            var -- 'x' (string)
            algorithm -- 'linbox' (default)
                         'generic'

        OUTPUT:
            a polynomial over the rational numbers.

        EXAMPLES:
            sage: a = matrix(QQ, 3, [4/3, 2/5, 1/5, 4, -3/2, 0, 0, -2/3, 3/4])
            sage: f = a.minpoly(); f           # optional -- os x only right now
            x^3 - 7/12*x^2 - 149/40*x + 97/30
            sage: a = Mat(ZZ,4)(range(16))
            sage: f = a.minpoly(); f.factor()  # optional -- os x only right now
            x * (x^2 - 30*x - 80)
            sage: f(a) == 0                    # optional -- os x only right now
            True
        """
        key = 'minpoly_%s_%s'%(algorithm, var)
        x = self.fetch(key)
        if x: return x

        if algorithm == 'linbox':
            A, denom = self._clear_denom()
            f = A.minpoly(var, algorithm='linbox')
            x = f.parent().gen()
            g = f(x * denom) * (1 / (denom**f.degree()))
        elif algorithm == 'generic':
            g = matrix_dense.Matrix_dense.minpoly(self, var)
        else:
            raise ValueError, "no algorithm '%s'"%algorithm

        self.cache(key, g)
        return g

    cdef sage.structure.element.Matrix _matrix_times_matrix_c_impl(self, sage.structure.element.Matrix right):
        return self._multiply_over_integers(right)

    def _multiply_over_integers(self, Matrix_rational_dense right, algorithm='default'):
        """
        Multiply this matrix by right using a multimodular algorithm
        and return the result.

        INPUT:
            self -- matrix over QQ
            right -- matrix over QQ
            algorithm -- 'default': use whatever is the defalt for A*B when A, B are over ZZ.
                         'multimodular': use a multimodular algorithm

        EXAMPLES:
            sage: a = MatrixSpace(QQ,10,5)(range(50))
            sage: b = MatrixSpace(QQ,5,12)([1/n for n in range(1,61)])
            sage: a._multiply_over_integers(b) == a._multiply_over_integers(b, algorithm='multimodular')
            True

            sage: a = MatrixSpace(QQ,3)(range(9))
            sage: b = MatrixSpace(QQ,3)([1/n for n in range(1,10)])
            sage: a._multiply_over_integers(b, algorithm = 'multimodular')
            [ 15/28   9/20   7/18]
            [  33/7 117/40   20/9]
            [249/28   27/5  73/18]

        """
        cdef Matrix_integer_dense A, B, AB
        cdef Matrix_rational_dense res
        cdef Integer D
        cdef mpz_t* AB_row,
        cdef mpq_t* res_row
        A, A_denom = self._clear_denom()
        B, B_denom = right._clear_denom()
        if algorithm == 'default':
            AB = A*B
        elif algorithm == 'multimodular':
            AB = A._multiply_multi_modular(B)
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm
        D = A_denom * B_denom
        res = Matrix_rational_dense.__new__(Matrix_rational_dense,
                                            self.matrix_space(AB._nrows, AB._ncols), 0, 0, 0)
        for i from 0 <= i < res._nrows:
            AB_row = AB._matrix[i]
            res_row = res._matrix[i]
            for j from 0 <= j < res._ncols:
                mpz_set(mpq_numref(res_row[0]), AB_row[0])
                mpz_set(mpq_denref(res_row[0]), D.value)
                mpq_canonicalize(res_row[0])
                AB_row = AB_row + 1
                res_row = res_row + 1
        _sig_off
        return res


    def height(self):
        """
        Return the height of this matrix, which is the least common
        multiple of all numerators and denominators of elements of
        this matrix.

        OUTPUT:
            -- SAGE Integer

        EXAMPLES:
            sage: b = matrix(QQ,2,range(6)); b[0,0]=-5007/293; b
            [-5007/293         1         2]
            [        3         4         5]
            sage: b.height()
            5007
        """
        cdef Integer z
        z = PY_NEW(Integer)
        self.mpz_height(z.value)
        return z

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

    cdef int _rescale(self, mpq_t a) except -1:
        cdef int i, j
        _sig_on
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                mpq_mul(self._matrix[i][j], self._matrix[i][j], a)
        _sig_off

    def _adjoint(self):
        """
        Return the adjoint of this matrix.

        Assumes self is a square matrix (checked in adjoint).
        """
        return self.parent()(self._pari_().matadjoint().python())

    def prod_of_row_sums(self, cols):
        cdef Py_ssize_t c, row
        cdef mpq_t s, pr
        mpq_init(s)
        mpq_init(pr)

        mpq_set_si(pr, 1, 1)
        for row from 0 <= row < self._nrows:
            tmp = []
            mpq_set_si(s, 0, 1)
            for c in cols:
                if c<0 or c >= self._ncols:
                    raise IndexError, "matrix column index out of range"
                mpq_add(s, s, self._matrix[row][c])
            mpq_mul(pr, pr, s)
        cdef Rational _pr
        _pr = PY_NEW(Rational)
        mpq_set(_pr.value, pr)
        mpq_clear(s)
        mpq_clear(pr)
        return _pr

    ################################################
    # Kernel
    ################################################
    def kernel(self):
        """
        Return the kernel of this matrix, as a vector space over QQ.
        """
        K = self.fetch('kernel')
        if not K is None:
            return K
        A, _ = self.transpose()._clear_denom()
        K = A._rational_kernel_iml().change_ring(QQ)
        V = K.column_space()
        self.cache('kernel', V)
        return V


    ################################################
    # Echelon form
    ################################################
    def echelonize(self, algorithm='default',
                   height_guess=None, proof=True, **kwds):
        """
        INPUT:
            algorithm -- 'default' (default): use heuristic choice
                         'padic': an algorithm based on the IML p-adic solver.
                         'multimodular': uses a multimodular algorithm the uses linbox
                                         modulo many primes.
            height_guess, proof, **kwds -- all passed to the multimodular algorithm; ignored
                                           by the p-adic algorithm.

        OUTPUT:
            matrix -- the reduced row echelon for of self.

        EXAMPLES:
            sage: a = matrix(QQ, 4, range(16)); a[0,0] = 1/19; a[0,1] = 1/5; a
            [1/19  1/5    2    3]
            [   4    5    6    7]
            [   8    9   10   11]
            [  12   13   14   15]
            sage: a.echelonize(); a
            [      1       0       0 -76/157]
            [      0       1       0  -5/157]
            [      0       0       1 238/157]
            [      0       0       0       0]

            sage: a = matrix(QQ, 4, range(16)); a[0,0] = 1/19; a[0,1] = 1/5
            sage: a.echelonize(algorithm='multimodular'); a
            [      1       0       0 -76/157]
            [      0       1       0  -5/157]
            [      0       0       1 238/157]
            [      0       0       0       0]
        """

        x = self.fetch('in_echelon_form')
        if not x is None: return  # already known to be in echelon form
        self.check_mutability()
        self.clear_cache()

        cdef Matrix_rational_dense E
        if algorithm == 'default':
            algorithm = 'padic'
        if algorithm == 'padic':
            pivots = self._echelonize_padic()
        elif algorithm == 'multimodular':
            pivots = self._echelonize_multimodular(height_guess, proof, **kwds)
        else:
            raise ValueError, "no algorithm '%s'"%algorithm
        self.cache('in_echelon_form', True)
        self.cache('pivots', pivots)


    def echelon_form(self, algorithm='default',
                     height_guess=None, proof=True, **kwds):
        """
        INPUT:
            algorithm -- 'default' (default): use heuristic choice
                         'padic': an algorithm based on the IML p-adic solver.
                         'multimodular': uses a multimodular algorithm the uses linbox
                                         modulo many primes.
            height_guess, proof, **kwds -- all passed to the multimodular algorithm; ignored
                                           by the p-adic algorithm.

        OUTPUT:
            self is no in reduced row echelon form.

        EXAMPLES:
            sage: a = matrix(QQ, 4, range(16)); a[0,0] = 1/19; a[0,1] = 1/5; a
            [1/19  1/5    2    3]
            [   4    5    6    7]
            [   8    9   10   11]
            [  12   13   14   15]
            sage: a.echelon_form()
            [      1       0       0 -76/157]
            [      0       1       0  -5/157]
            [      0       0       1 238/157]
            [      0       0       0       0]
            sage: a.echelon_form(algorithm='multimodular')
            [      1       0       0 -76/157]
            [      0       1       0  -5/157]
            [      0       0       1 238/157]
            [      0       0       0       0]
        """
        label = 'echelon_form_%s'%algorithm
        x = self.fetch(label)
        if not x is None:
            return x
        if self.fetch('in_echelon_form'): return self

        if algorithm == 'default':
            algorithm = 'padic'

        if algorithm == 'padic':
            E = self._echelon_form_padic()
        elif algorithm == 'multimodular':
            E = self._echelon_form_multimodular(height_guess, proof=proof)
        else:
            raise ValueError, "no algorithm '%s'"%algorithm
        self.cache(label, E)
        self.cache('pivots', E.pivots())
        return E


    # p-adic echelonization algorithms
    def _echelon_form_padic(self, include_zero_rows=True):
        """
        Compute and return the echelon form of self using a p-adic nullspace algorithm.
        """
        cdef Matrix_integer_dense X
        cdef Matrix_rational_dense E
        cdef Integer d
        cdef mpq_t* E_row
        cdef mpz_t* X_row

        t = verbose('Computing echelon form of %s x %s matrix over QQ using p-adic nullspace algorithm.'%(
            self.nrows(), self.ncols()))
        A, _ = self._clear_denom()
        t = verbose('  Got integral matrix', t)
        pivots, nonpivots, X, d = A._rational_echelon_via_solve()
        t = verbose('  Computed ZZ-echelon using p-adic algorithm.', t)

        nr = self.nrows() if include_zero_rows else X.nrows()
        parent = self.matrix_space(nr, self.ncols())
        E = Matrix_rational_dense.__new__(Matrix_rational_dense, parent, None, None, None)

        # Fill in the identity part of the matrix
        cdef Py_ssize_t i, j
        for i from 0 <= i < len(pivots):
            mpz_set_si(mpq_numref(E._matrix[i][pivots[i]]), 1)

        # Fill in the non-pivot part of the matrix
        for i from 0 <= i < X.nrows():
            E_row = E._matrix[i]
            X_row = X._matrix[i]
            for j from 0 <= j < X.ncols():
                mpz_set(mpq_numref(E_row[nonpivots[j]]), X_row[j])
                mpz_set(mpq_denref(E_row[nonpivots[j]]), d.value)
                mpq_canonicalize(E_row[nonpivots[j]])

        t = verbose('Reconstructed solution over QQ, thus completing the echelonize', t)
        E.cache('in_echelon_form', True)
        E.cache('pivots', pivots)
        return E

    def _echelonize_padic(self):
        """
        Echelonize self using a p-adic nullspace algorithm.
        """
        cdef Matrix_integer_dense X
        cdef Integer d
        cdef mpq_t* E_row
        cdef mpz_t* X_row

        t = verbose('Computing echelonization of %s x %s matrix over QQ using p-adic nullspace algorithm.'%
                    (self.nrows(), self.ncols()))
        A, _ = self._clear_denom()
        t = verbose('  Got integral matrix', t)
        pivots, nonpivots, X, d = A._rational_echelon_via_solve()
        t = verbose('  Computed ZZ-echelon using p-adic algorithm.', t)

        # Fill in the identity part of self.
        cdef Py_ssize_t i, j, k
        for j from 0 <= j < len(pivots):
            k = pivots[j]
            for i from 0 <= i < len(pivots):
                if i == j:
                    mpq_set_si(self._matrix[i][k], 1, 1)
                else:
                    mpq_set_si(self._matrix[i][k], 0, 1)


        # Fill in the non-pivot part of self.
        for i from 0 <= i < X.nrows():
            E_row = self._matrix[i]
            X_row = X._matrix[i]
            for j from 0 <= j < X.ncols():
                mpz_set(mpq_numref(E_row[nonpivots[j]]), X_row[j])
                mpz_set(mpq_denref(E_row[nonpivots[j]]), d.value)
                mpq_canonicalize(E_row[nonpivots[j]])

        # Fill in the 0-rows at the bottom.
        for i from len(pivots) <= i < self._nrows:
            E_row = self._matrix[i]
            for j from 0 <= j < self._ncols:
                mpq_set_si(E_row[j], 0, 1)

        t = verbose('Filled in echelonization of self, thus completing the echelonize', t)
        return pivots


    # Multimodular echelonization algorithms
    def _echelonize_multimodular(self, height_guess=None, proof=True, **kwds):
        cdef Matrix_rational_dense E
        E = self._echelon_form_multimodular(height_guess, proof=proof, **kwds)
        cdef Py_ssize_t i, j
        cdef mpq_t *row0, *row1
        for i from 0 <= i < self._nrows:
            row0 = self._matrix[i]
            row1 = E._matrix[i]
            for j from 0 <= j < self._ncols:
                mpq_set(row0[j], row1[j])
        return E.pivots()

    def _echelon_form_multimodular(self, height_guess=None, proof=True):
        """
        Returns reduced row-echelon form using a multi-modular
        algorithm.  Does not change self.

        REFERENCE: Chapter 7 of Stein's "Explicitly Computing Modular Forms".

        INPUT:
            height_guess -- integer or None
            proof -- boolean (default: True)
        """
        import misc
        return misc.matrix_rational_echelon_form_multimodular(self,
                                 height_guess=height_guess, proof=proof)

    def _decomposition_devel(self, is_diagonalizable=False, dual=False,
                            echelon_algorithm='default', height_guess=None):
        """
        Returns the decomposition of the free module on which this
        matrix A acts from the right (i.e., the action is x goes to x
        A), along with whether this matrix acts irreducibly on each
        factor.  The factors are guaranteed to be sorted in the same
        way as the corresponding factors of the characteristic
        polynomial.

        INPUT:
            self -- a matrix over a field

        OUTPUT:
            Sequence -- list of pairs (V,t), where V is a vector spaces
                    and t is a bool, and t is True exactly when the
                    charpoly of self on V is irreducible.

            (optional) list -- list of pairs (W,t), where W is a vector
                    space and t is a bool, and t is True exactly
                    when the charpoly of the transpose of self on W
                    is irreducible.
        """
        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"

        if not self.base_ring().is_field():
            raise TypeError, "self must be over a field."

        if self.nrows() == 0:
            return decomp_seq([])

        f = self.charpoly('x')
        E = decomp_seq([])

        if dual:
            Edual = decomp_seq([])

        t = verbose('decomposition -- factoring the characteristic polynomial', level=2)
        F = f.factor()
        verbose('decomposition -- done factoring', t=t, level=2)

        if len(F) == 1:
            V = QQ**self.nrows()
            m = F[0][1]
            if dual:
                return decomp_seq([(V,m==1)]), decomp_seq([(V,m==1)])
            else:
                return decomp_seq([(V,m==1)])


        V = ZZ**self.nrows()
        v = V.random_element()

        num_iterates = max([f.degree() - g.degree() for g, _ in F]) + 1

        A, _ = self._clear_denom()

        t = verbose('decomposition -- computing %s iterates of a random vector.'%num_iterates)
        S = A.iterates(v, num_iterates)
        verbose('decomposition -- done computing iterates.', t = t)

        F.sort()
        for i in range(len(F)):
            g, m = F[i]

            # Compute the complementary factor.
            h = f // (g**m)
            h = h * h.denominator()

            # Compute one element of the kernel of g(self)**m.
            t = verbose('decomposition -- compute element of kernel of g(self), for g of degree %s'%g.degree(),level=2)
            w = S.linear_combination_of_rows(h.list())
            t = verbose('decomposition -- done computing element of kernel of g(self)', t=t,level=2)

            # Get the rest of the kernel.
            t = verbose('decomposition -- fill out rest of kernel',level=2)
            W = A.iterates(w, g.degree())
            t = verbose('decomposition -- finished filling out rest of kernel',level=2, t=t)

            W = W.change_ring(QQ)

            t = verbose('decomposition -- now computing row space', level=2)
            W.echelonize(algorithm = echelon_algorithm, height_guess=height_guess)
            E.append((W.row_space(), m==1))
            verbose('decomposition -- computed row space', level=2,t=t)

        return E

    ################################################################
    # second implementation of the above, usually over twice as fast
    # even without denominator lcm trick
    # TODO: merge with the above
    def _echelon_multimodular(self, height_guess=None, proof=True):
        """
        _echelon_multimodular(self, height_guess=None):

        Returns echelon form of self, without modifying self.  Uses a
        multi-modular method.

        REFERENCE: Chapter 7 of Stein's "Explicitly Computing Modular Forms".

        INPUT:
            self -- matrix with n columns (this).
            height_guess -- integer or None
            proof -- boolean (default: True)

        ALGORITHM:
        The following is a modular algorithm for computing the echelon
        form.  Define the height of a matrix to be the max of the
        absolute values of the entries.

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

        AUTHORS:
            -- William Stein
            -- Robert Bradshaw
        """
        cdef Matrix_integer_dense A
        cdef Matrix_rational_dense E
        cdef Integer d
        cdef int problem
        A, d = self._clear_denom()
        hA = A.height()
        if height_guess is None:
            height_guess = (2*hA)**(self._ncols/2+1)
        tm = verbose("height_guess = %s"%height_guess, level=2)
        best_pivots = []

        if proof:
            M = self._ncols * height_guess * hA  +  1
        else:
            M = height_guess + 1
        mm = sage.ext.multi_modular.MutableMultiModularBasis(M)

        res = []
        # reduction via several primes can be made more efficient than reduction via each prime one at a time
        t = verbose("Reducing mod %s:"%mm, level=2)
        new_res = A._reduce(mm) # TODO: can I recognize special forms (e.g. identity) before calculating all of these?
        t = verbose("time to reduce matrix mod p:",t, level=2)
        problem = 0

        _sig_on
        while True:
            # calculate the new echelon forms
            for B in new_res:
                B.echelonize()
                t = verbose("time to put reduced matrix in echelon form:",t, level=2)
            i = len(res)
            res += new_res
            new_res = []
            # make sure they all have the same pivots
            while i < len(res):
                c = cmp_pivots(best_pivots, res[i].pivots())
                if c == 0:
                    i += 1
                elif c < 0:
                    best_pivots = res[i].pivots()
                    i = 0
                else:
                    p = mm.replace_prime(i)
                    res[i] = A._mod_int(p)
                    res[i].echelonize()
                    verbose("Excluding this prime (bad pivots).")
            t = verbose("time for pivot compare", t, level=2)
            # now try and lift
            try:
                t = verbose("start crt/rr", level=1)
                E = self._lift_crt_rr_with_lcm(res, mm)
                verbose("done crt/rr", t, level=1)
            except ValueError:
                mm._extend_moduli(3)
                new_res = A._reduce(mm[-3:])
                verbose("Failed to compute rational reconstruction -- redoing with several more primes", level=2)
                continue

            # see if we have enough clearance for the height
            if not proof:
                verbose("Not checking validity of result (since proof=False).", level=2)
                break

            dE, d = E._clear_denom()
            hE = dE.height()
            if hE * hA * self._ncols < mm.prod():
                self.cache('pivots', best_pivots)
                (<Matrix_rational_dense>E).cache('pivots', best_pivots)
                break

            # try a few more primes
            mm._extend_moduli(3)
            new_res = A._reduce(mm[-3:])
            problem += 1
            if problem > 50:
                verbose("sparse_matrix multi-modular reduce not converging?")

        #end while
        _sig_off
        verbose("total time", tm, level=2)
        self.cache('pivots', best_pivots)
        E.cache('pivots', best_pivots)
        return E

    ################################################################

    def _lift_crt_rr(self, res, mm):
        cdef Integer m
        cdef Matrix_integer_dense ZA
        cdef Matrix_rational_dense QA
        cdef Py_ssize_t i, j, nr, nc
        cdef mpz_t* Z_row
        cdef mpq_t* Q_row

        ZA = _lift_crt(res, mm)
        nr = ZA._nrows
        nc = ZA._ncols
        QA = Matrix_rational_dense.__new__(Matrix_rational_dense, self.parent(), None, None, None)
        m = mm.prod()
        for i from 0 <= i < nr:
            Z_row = ZA._matrix[i]
            Q_row = QA._matrix[i]
            for j from 0 <= j < nc:
                mpq_rational_reconstruction(Q_row[j], Z_row[j], m.value)
        return QA

    def _lift_crt_rr_with_lcm(self, res, mm):
        """
            Optimizations: When doing the rational_recon lift of a (mod m)
            first see if |a| < sqrt(m/2) in which case it lifts to
            an integer (often a=0 or 1).

            If that fails, keep track of the lcm d of denominators found so far,
            and check to see if z = a*d lifts to an integer with |z| <= sqrt(m/2).
            If so, no need to do rational recon.  This should be the case
            for most a after a while, and should saves substantial time!
        """
        cdef Integer m
        cdef Matrix_integer_dense ZA
        cdef Matrix_rational_dense QA
        cdef Py_ssize_t i, j, nr, nc
        cdef mpz_t* Z_row
        cdef mpq_t* Q_row
        cdef mpz_t lcm_denom, sqrt_m, neg_sqrt_m, z

        mpz_init(z)
        mpz_init(sqrt_m)
        mpz_init(neg_sqrt_m)
        mpz_init_set_ui(lcm_denom, 1)

        m = mm.prod()
        mpz_fdiv_q_2exp(sqrt_m, m.value, 1)
        mpz_sqrt(sqrt_m, sqrt_m)
        mpz_sub(neg_sqrt_m, m.value, sqrt_m)

        t = verbose("Starting crt", level=2)
        ZA = _lift_crt(res, mm)
        t = verbose("crt finished", t, level=2)
        nr = ZA._nrows
        nc = ZA._ncols
        QA = Matrix_rational_dense.__new__(Matrix_rational_dense, self.parent(), None, None, None)

        cdef int is_integral, lcm_trick
        is_integral = 0
        lcm_trick = 0

        t = verbose("Starting rational reconstruction", level=2)
        for i from 0 <= i < nr:
            Z_row = ZA._matrix[i]
            Q_row = QA._matrix[i]
            for j from 0 <= j < nc:
                if mpz_cmp(Z_row[j], sqrt_m) < 0:
                    mpz_set(mpq_numref(Q_row[j]), Z_row[j])
                    is_integral += 1
                elif mpz_cmp(Z_row[j], neg_sqrt_m) > 0:
                    mpz_sub(mpq_numref(Q_row[j]), Z_row[j], m.value)
                    is_integral += 1
                else:
                    mpz_mul(z, Z_row[j], lcm_denom)
                    mpz_fdiv_r(z, z, m.value)
                    if mpz_cmp(z, sqrt_m) < 0:
                        mpz_set(mpq_numref(Q_row[j]), z)
                        mpz_set(mpq_denref(Q_row[j]), lcm_denom)
                        mpq_canonicalize(Q_row[j])
                        lcm_trick += 1
                    elif mpz_cmp(z, neg_sqrt_m) > 0:
                        mpz_sub(mpq_numref(Q_row[j]), z, m.value)
                        mpz_set(mpq_denref(Q_row[j]), lcm_denom)
                        mpq_canonicalize(Q_row[j])
                        lcm_trick += 1
                    else:
                        mpq_rational_reconstruction(Q_row[j], Z_row[j], m.value)
                        mpz_lcm(lcm_denom, lcm_denom, mpq_denref(Q_row[j]))
        mpz_clear(z)
        mpz_clear(sqrt_m)
        mpz_clear(neg_sqrt_m)
        mpz_clear(lcm_denom)
        t = verbose("rr finished. integral entries: %s, lcm trick: %s, other: %s"%(is_integral, lcm_trick, nr*nc - is_integral - lcm_trick), t, level=2)
        return QA


    def randomize(self, density=1, num_bound=2, den_bound=2):
        """
        Randomize density proportion of the entries of this matrix to
        be rationals with numerators and denominators at most the
        given bounds.
        """
        density = float(density)
        if density == 0:
            return
        self.check_mutability()
        self.clear_cache()

        cdef Integer B, C
        B = Integer(num_bound+1)
        C = Integer(den_bound+1)

        cdef Py_ssize_t i, j, k, nc, num_per_row
        global state

        cdef double total
        total = self._nrows * self._ncols
        cdef int r, s
        r = self._nrows * self._ncols

        _sig_on
        if density == 1:
            if mpz_cmp_si(C.value, 2):   # denom is > 1
                for i from 0 <= i < self._nrows*self._ncols:
                    mpq_randomize_entry(self._entries[i], B.value, C.value)
            else:
                for i from 0 <= i < self._nrows*self._ncols:
                    mpq_randomize_entry_as_int(self._entries[i], B.value)
        else:
            nc = self._ncols
            num_per_row = int(density * nc)
            if mpz_cmp_si(C.value, 2):   # denom is > 1
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < num_per_row:
                        k = random()%nc
                        mpq_randomize_entry(self._matrix[i][k], B.value, C.value)
            else:
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < num_per_row:
                        k = random()%nc
                        mpq_randomize_entry_as_int(self._matrix[i][k], B.value)
        _sig_off




cdef decomp_seq(v):
    return Sequence(v, universe=tuple, check=False, cr=True)
