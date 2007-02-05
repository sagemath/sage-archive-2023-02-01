"""
Dense matrices over the rational field.
"""

##############################################################################
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

cdef extern from "matrix_rational_dense_linbox.h":
    void linbox_rational_dense_echelon_form(mpq_t** matrix, size_t nr, size_t nc)


include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"

cimport sage.structure.element
from sage.rings.rational cimport Rational
from matrix cimport Matrix
from matrix_integer_dense cimport Matrix_integer_dense
import sage.structure.coerce
from sage.structure.element cimport ModuleElement, RingElement
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.finite_field import GF

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

        self._entries = <mpq_t *> sage_malloc(sizeof(mpq_t)*(self._nrows * self._ncols))
        if self._entries == NULL:
            raise MemoryError, "out of memory allocating a matrix"

        self._matrix =  <mpq_t **> sage_malloc(sizeof(mpq_t*) * self._nrows)
        if self._matrix == NULL:
            raise MemoryError, "out of memory allocating a matrix"

        # store pointers to the starts of the rows
        k = 0
        for i from 0 <= i < self._nrows:
            self._matrix[i] = self._entries + k
            k = k + self._ncols

        for i from 0 <= i < self._nrows * self._ncols:
            mpq_init(self._entries[i])

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
    #   * cdef _mul_c_impl
    #   * cdef _cmp_c_impl
    # x * __neg__
    #   * __invert__
    # x * __copy__
    #   * _multiply_classical
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
    ########################################################################
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
        cdef int i, j
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
            sage: f = a.minpoly(); f
            x^3 - 7/12*x^2 - 149/40*x + 97/30
            sage: a = Mat(ZZ,4)(range(16))
            sage: f = a.minpoly(); f.factor()
            x * (x^2 - 30*x - 80)
            sage: f(a) == 0
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
    # Echelon form
    ################################################
    def echelonize(self, height_guess=None, proof=True, **kwds):
        x = self.fetch('in_echelon_form')
        if not x is None: return  # already known to be in echelon form
        self.check_mutability()
        cdef Matrix_rational_dense E
        E = self._echelon_form_multimodular(height_guess, proof=proof)
        cdef Py_ssize_t i, j
        cdef mpq_t *row0, *row1
        for i from 0 <= i < self._nrows:
            row0 = self._matrix[i]
            row1 = E._matrix[i]
            for j from 0 <= j < self._ncols:
                mpq_set(row0[j], row1[j])
        self.cache('in_echelon_form', True)
        self.cache('pivots', E.pivots())

    def echelon_form(self, height_guess=None, proof=True, **kwds):
        x = self.fetch('echelon_form')
        if not x is None:
            return x
        cdef Matrix_rational_dense E
        E = self._echelon_form_multimodular(height_guess, proof=proof)
        self.cache('echelon_form', E)
        self.cache('pivots', E.pivots())
        return E

    def _echelonize_linbox(self):
        raise NotImplementedError
        linbox_rational_dense_echelon_form(self._matrix, self._nrows, self._ncols)

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
