"""
Dense matrices over the integer ring.
"""

######################################################################
#       Copyright (C) 2006,2007 William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
######################################################################

from sage.misc.misc import verbose, get_verbose

include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/gmp.pxi"

cdef extern from "matrix_integer_dense_linbox.h":
    void linbox_integer_dense_minpoly(mpz_t* *minpoly, size_t* degree,
                                      size_t n, mpz_t** matrix, int do_minpoly)
    void linbox_integer_dense_charpoly(mpz_t* *charpoly, size_t* degree,
                                       size_t n, mpz_t** matrix)
    void linbox_integer_dense_delete_array(mpz_t* f)
    int linbox_integer_dense_matrix_matrix_multiply(mpz_t** ans, mpz_t **A, mpz_t **B,
                                      size_t A_nr, size_t A_nc, size_t B_nr, size_t B_nc)

ctypedef unsigned int uint

from sage.ext.multi_modular cimport MultiModularBasis
cdef MultiModularBasis mm
mm = MultiModularBasis()

from sage.rings.integer cimport Integer
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.integer_mod_ring import IntegerModRing
from sage.rings.polynomial_ring import PolynomialRing
from sage.structure.element cimport ModuleElement

from matrix_modn_dense import Matrix_modn_dense
from matrix_modn_dense cimport Matrix_modn_dense

import sage.modules.free_module

from matrix cimport Matrix

import matrix_space


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
        self._initialized = 0
        self._nrows = parent.nrows()
        self._ncols = parent.ncols()
        self._pivots = None

        # Allocate an array where all the entries of the matrix are stored.
        self._entries = <mpz_t *>sage_malloc(sizeof(mpz_t) * (self._nrows * self._ncols))
        if self._entries == NULL:
            raise MemoryError, "out of memory allocating a matrix"

        # Allocate an array of pointers to the rows, which is useful for
        # certain algorithms.
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

    def __copy__(self):
        cdef Matrix_integer_dense A
        A = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent,
                                         0, 0, 0)
        cdef Py_ssize_t i
        for i from 0 <= i < self._nrows * self._ncols:
            mpz_init_set(A._entries[i], self._entries[i])
        return A

    def  __dealloc__(self):
        """
        Frees all the memory allocated for this matrix.
        EXAMPLE:
            sage: a = Matrix(ZZ,2,[1,2,3,4])
            sage: del a
        """
        if self._entries == NULL: return
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

        if not isinstance(entries, list):  # todo -- change to PyObject_TypeCheck???
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
                raise TypeError, "entries has the wrong length"
            if coerce:
                for i from 0 <= i < self._nrows * self._ncols:
                    # TODO: Should use an unsafe un-bounds-checked array access here.
                    x = ZZ(entries[i])
                    # todo -- see integer.pyx and the TODO there; perhaps this could be
                    # sped up by creating a mpz_init_set_sage function.
                    mpz_init_set(self._entries[i], x.value)
            else:
                for i from 0 <= i < self._nrows * self._ncols:
                    # TODO: Should use an unsafe un-bounds-checked array access here.
                    mpz_init_set(self._entries[i], (<Integer> entries[i]).value)

        else:

            # If x is zero, make the zero matrix and be done.
            if mpz_cmp_si(x.value, 0) == 0:
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
            self._initialized = 1


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
        self._initialized = 1

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
    #   * cdef _mul_c_impl
    #   * cdef _cmp_c_impl
    #   * __neg__
    #   * __invert__
    #   * __copy__
    # x * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################

    # cdef _mul_c_impl(self, Matrix right):
    # cdef int _cmp_c_impl(self, Matrix right) except -2:
    # def __neg__(self):
    # def __invert__(self):
    # def __copy__(self):
    # def _multiply_classical(left, matrix.Matrix _right):
    # def _list(self):
    # def _dict(self):

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
        ans = self.new_matrix(nrows = self.nrows(), ncols = right.ncols())
        B = right
        _sig_on
        e = linbox_integer_dense_matrix_matrix_multiply(ans._matrix, self._matrix, B._matrix,
                                          self._nrows, self._ncols,
                                          right._nrows, right._ncols)
        _sig_off
        if e:
            raise RuntimeError
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
        cdef Py_ssize_t i

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent, None, None, None)
        Matrix.__init__(M, self._parent)

        _sig_on

        cdef mpz_t *entries
        entries = M._entries
        for i from 0 <= i < self._ncols * self._nrows:
            mpz_init(entries[i])
            mpz_add(entries[i], self._entries[i], (<Matrix_integer_dense> right)._entries[i])
        _sig_off
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
        cdef Py_ssize_t i

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent, None, None, None)
        Matrix.__init__(M, self._parent)

        _sig_on

        cdef mpz_t *entries
        entries = M._entries
        for i from 0 <= i < self._ncols * self._nrows:
            mpz_init(entries[i])
            mpz_sub(entries[i], self._entries[i], (<Matrix_integer_dense> right)._entries[i])

        _sig_off
        return M



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

        EXAMPLES:
            sage: A = matrix(ZZ,6, range(36))
            sage: f = A.charpoly(); f
            x^6 - 105*x^5 - 630*x^4
            sage: f(A) == 0
            True
            sage: n=20; A = Mat(ZZ,n)(range(n^2))
            sage: A.charpoly()
            x^20 - 3990*x^19 - 266000*x^18
            sage: A.minpoly()
            x^3 - 3990*x^2 - 266000*x
        """
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
        INPUT:
            var -- a variable name
            algorithm -- 'linbox' (default)
                         'generic'

        EXAMPLES:
            sage: A = matrix(ZZ,6, range(36))
            sage: A.minpoly()
            x^3 - 105*x^2 - 630*x
            sage: n=6; A = Mat(ZZ,n)([k^2 for k in range(n^2)])
            sage: A.minpoly()
            x^4 - 2695*x^3 - 257964*x^2 + 1693440*x
        """
        if algorithm == 'linbox':
            g = self._minpoly_linbox(var)
        elif algorithm == 'generic':
            g = self._minpoly_generic(var)
        else:
            raise ValueError, "no algorithm '%s'"%algorithm
        self.cache('minpoly_%s_%s'%(algorithm, var), g)
        return g

    def _minpoly_linbox(self, var='x'):
        return self._poly_linbox(var=var, typ='minpoly')

    def _charpoly_linbox(self, var='x'):
        return self._poly_linbox(var=var, typ='charpoly')

    def _poly_linbox(self, var='x', typ='minpoly'):
        """
        INPUT:
            var -- 'x'
            typ -- 'minpoly' or 'charpoly'
        """
        if self._nrows != self._ncols:
            raise ValueError, "matrix must be square"
        if self._nrows <= 1:
            return matrix_dense.Matrix_dense.charpoly(self, var)
        cdef mpz_t* poly
        cdef size_t n
        cdef size_t degree
        if typ == 'minpoly':
            _sig_on
            linbox_integer_dense_minpoly(&poly, &degree, self._nrows, self._matrix, 1)
            _sig_off
        else:
            _sig_on
            linbox_integer_dense_charpoly(&poly, &degree, self._nrows, self._matrix)
            _sig_off

        v = []
        cdef Integer k
        for n from 0 <= n <= degree:
            k = PY_NEW(Integer)
            mpz_set(k.value, poly[n])
            mpz_clear(poly[n])
            v.append(k)
        linbox_integer_dense_delete_array(poly)
        R = self._base_ring[var]
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
        cdef int i, n

        h = left.height() * right.height()
        n = mm.moduli_list_c(&moduli, h.value)
        res = []
        for i from 0 <= i < n:
            res.append(left._mod_int_c(moduli[i]) * right._mod_int_c(moduli[i]))
        sage_free(moduli)
        return left._lift_crt(res)

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
            res_row = res.matrix[i]
            for j from 0 <= j < self._ncols:
                res_row[j] = mpz_fdiv_ui(self_row[j], p)
        return res

    def _lift_crt(self, residues):

        cdef size_t n, i, j, k
        cdef Py_ssize_t nr, nc

        n = len(residues)
        nr = residues[0].nrows()
        nc = residues[0].ncols()

        for b in residues:
            if not PY_TYPE_CHECK(b, Matrix_modn_dense):
                raise TypeError, "Can only perform CRT on list of type Matrix_modn_dense."
        cdef PyObject** res
        res = FAST_SEQ_UNSAFE(residues)

        cdef mod_int **row_list
        row_list = <mod_int**>sage_malloc(sizeof(mod_int*) * n)
        if row_list == NULL:
            raise MemoryError, "out of memory allocating multi-modular coefficent list"

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, self.matrix_space(nr, nc), None, None, None)

        _sig_on
        for i from 0 <= i < nr:
            for k from 0 <= k < n:
                row_list[k] = (<Matrix_modn_dense>res[k]).matrix[i]
            mm.mpz_crt_vec(M._matrix[i], row_list, n, nc)
        _sig_off

        sage_free(row_list)
        return M


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

    def elementary_divisors(self):
        """
        Return the elementary divisors of self, in order.

        The elementary divisors are the invariants of the finite
        abelian group that is the cokernel of this matrix.  They are
        ordered in reverse by divisibility.

        INPUT:
            matrix
        OUTPUT:
            list of int's

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
            sage: A.elementary_divisors()
            [0, 3, 1]
            sage: C = MatrixSpace(ZZ,4)([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
            sage: C.elementary_divisors()
            [687, 1, 1, 1]

        SEE ALSO: smith_form
        """
        d = self.fetch('elementary_divisors')
        if not d is None:
            return d
        if self._nrows == 0 or self._ncols == 0:
            d = []
        else:
            d = self._pari_().matsnf(0).python()
        self.cache('elementary_divisors', d)
        return d

    def smith_form(self):
        """
        Returns matrices S, U, and V such that S = U*self*V, and S
        is in Smith normal form.  Thus S is diagonal with diagonal
        entries the ordered elementary divisors of S.

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

##     def _linbox(self):
##         cdef Py_ssize_t i, j
##         s = '%s %s M\n'%(self._nrows, self._ncols)
##         for i from 0 <= i < self._nrows:
##             for j from 0 <= j < self._ncols:
##                 if mpz_cmp_si(self._matrix[i][j], 0):
##                     s += '%s %s %s\n'%(i+1,j+1,self.get_unsafe(i,j))
##         return s

    def _linbox(self):
        cdef Py_ssize_t i, j
        s = '%s %s x'%(self._nrows, self._ncols)
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                if mpz_cmp_si(self._matrix[i][j], 0):
                    s += ' %s'%self.get_unsafe(i,j)
        return s


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


