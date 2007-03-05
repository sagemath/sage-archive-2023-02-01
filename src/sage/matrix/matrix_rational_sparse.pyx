"""
Sparse rational matrices.

AUTHORS:
    -- William Stein (2007-02-21)
    -- Soroosh Yazdani (2007-02-21)
"""

##############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

include '../modules/binary_search.pxi'

include '../modules/vector_integer_sparse_h.pxi'
include '../modules/vector_integer_sparse_c.pxi'
include '../modules/vector_rational_sparse_h.pxi'
include '../modules/vector_rational_sparse_c.pxi'
include '../ext/stdsage.pxi'

include '../ext/interrupt.pxi'

from sage.rings.rational cimport Rational
from sage.rings.integer  cimport Integer
from matrix cimport Matrix

from sage.rings.integer_ring import ZZ
import sage.matrix.matrix_space

from matrix_integer_sparse cimport Matrix_integer_sparse

cdef class Matrix_rational_sparse(matrix_sparse.Matrix_sparse):

    ########################################################################
    # LEVEL 1 functionality
    #   * __new__
    #   * __dealloc__
    #   * __init__
    #   * set_unsafe
    #   * get_unsafe
    #   * __richcmp__    -- always the same
    #   * __hash__       -- alway simple
    ########################################################################
    def __new__(self, parent, entries, copy, coerce):
        # set the parent, nrows, ncols, etc.
        matrix_sparse.Matrix_sparse.__init__(self, parent)

        self._matrix = <mpq_vector*> sage_malloc(parent.nrows()*sizeof(mpq_vector))
        if self._matrix == NULL:
            raise MemoryError, "error allocating sparse matrix"
        # initialize the rows
        for i from 0 <= i < parent.nrows():
            mpq_vector_init(&self._matrix[i], self._ncols, 0)

        # record that rows have been initialized
        self._initialized = True


    def __dealloc__(self):
        self._dealloc()

    cdef _dealloc(self):
        cdef Py_ssize_t i
        if self._initialized:
            for i from 0 <= i < self._nrows:
                mpq_vector_clear(&self._matrix[i])
        if self._matrix != NULL:
            sage_free(self._matrix)

    def __init__(self, parent, entries, copy, coerce):
        """
        Create a sparse matrix over the rational numbers

        INPUT:
            parent -- a matrix space
            entries -- * a Python list of triples (i,j,x), where 0 <= i < nrows,
                         0 <= j < ncols, and x is coercible to an int.  The i,j
                         entry of self is set to x.  The x's can be 0.
                       * Alternatively, entries can be a list of *all* the entries
                         of the sparse matrix (so they would be mostly 0).
            copy -- ignored
            coerce -- ignored
        """
        cdef Py_ssize_t i, j, k
        cdef Rational z
        cdef void** X

        # fill in entries in the dict case
        if isinstance(entries, dict):
            R = self.base_ring()
            for ij, x in entries.iteritems():
                z = R(x)
                if z != 0:
                    i, j = ij  # nothing better to do since this is user input, which may be bogus.
                    if i < 0 or j < 0 or i >= self._nrows or j >= self._ncols:
                        raise IndexError, "invalid entries list"
                    mpq_vector_set_entry(&self._matrix[i], j, z.value)
        elif isinstance(entries, list):
            # Dense input format -- fill in entries
            if len(entries) != self._nrows * self._ncols:
                raise TypeError, "list of entries must be a dictionary of (i,j):x or a dense list of n * m elements"
            seq = PySequence_Fast(entries,"expected a list")
            X = PySequence_Fast_ITEMS(seq)
            k = 0
            R = self.base_ring()
            # Get fast access to the entries list.
            for i from 0 <= i < self._nrows:
                for  j from 0 <= j < self._ncols:
                    z = R(<object>X[k])
                    if z != 0:
                        mpq_vector_set_entry(&self._matrix[i], j, z.value)
                    k = k + 1
        else:
            # fill in entries in the scalar case
            z = Rational(entries)
            if z == 0:
                return
            if self._nrows != self._ncols:
                raise TypeError, "matrix must be square to initialize with a scalar."
            for i from 0 <= i < self._nrows:
                mpq_vector_set_entry(&self._matrix[i], i, z.value)


    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        mpq_vector_set_entry(&self._matrix[i], j, (<Rational> x).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef Rational x
        x = Rational()
        mpq_vector_get_entry(&x.value, &self._matrix[i], j)
        return x

    def __richcmp__(Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)
    def __hash__(self):
        return self._hash()


    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    #   * cdef _add_c_impl
    #   * cdef _sub_c_impl
    #   * cdef _mul_c_impl
    #   * cdef _cmp_c_impl
    #   * __neg__
    #   * __invert__
    #   * __copy__
    #   * _multiply_classical
    #   * _matrix_times_matrix_c_impl
    #   * _list -- list of underlying elements (need not be a copy)
    #   * x _dict -- sparse dictionary of underlying elements (need not be a copy)
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

# TODO
##     cdef ModuleElement _lmul_c_impl(self, RingElement right):
##         """
##         EXAMPLES:
##             sage: a = matrix(QQ,2,range(6))
##             sage: (3/4) * a
##             [   0  3/4  3/2]
##             [ 9/4    3 15/4]
##         """

    def _dict(self):
        """
        Unsafe version of the dict method, mainly for internal use.
        This may return the dict of elements, but as an *unsafe*
        reference to the underlying dict of the object.  It is might
        be dangerous if you change entries of the returned dict.
        """
        d = self.fetch('dict')
        if not d is None:
            return d

        cdef Py_ssize_t i, j, k
        d = {}
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._matrix[i].num_nonzero:
                x = Integer()
                mpz_set((<Integer>x).value, self._matrix[i].entries[j])
                d[(int(i),int(self._matrix[i].positions[j]))] = x
        self.cache('dict', d)
        return d


    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_c_impl
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    ########################################################################

    def height(self):
        """
        Return the height of this matrix, which is the least common
        multiple of all numerators and denominators of elements of
        this matrix.

        OUTPUT:
            -- Integer

        EXAMPLES:
            sage: b = matrix(QQ,2,range(6), sparse=True); b[0,0]=-5007/293; b
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
            for j from 0 <= j < self._matrix[i].num_nonzero:
                mpq_get_num(x, self._matrix[i].entries[j])
                mpz_abs(x, x)
                if mpz_cmp(h,x) < 0:
                    mpz_set(h,x)
                mpq_get_den(x, self._matrix[i].entries[j])
                mpz_abs(x, x)
                if mpz_cmp(h,x) < 0:
                    mpz_set(h,x)
        _sig_off
        mpz_set(height, h)
        mpz_clear(h)
        mpz_clear(x)
        return 0

    cdef int mpz_denom(self, mpz_t d) except -1:
        mpz_set_si(d,1)
        cdef Py_ssize_t i, j
        cdef mpq_vector *v

        _sig_on
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._matrix[i].num_nonzero:
                mpz_lcm(d, d, mpq_denref(self._matrix[i].entries[j]))
        _sig_off
        return 0


    def denominator(self):
        """
        Return the denominator of this matrix.

        OUTPUT:
            -- SAGE Integer

        EXAMPLES:
            sage: b = matrix(QQ,2,range(6)); b[0,0]=-5007/293; b
            [-5007/293         1         2]
            [        3         4         5]
            sage: b.denominator()
            293
        """
        cdef Integer z
        z = PY_NEW(Integer)
        self.mpz_denom(z.value)
        return z

    def _clear_denom(self):
        """
        INPUT:
            self -- a matrix

        OUTPUT:
            D*self, D

        The product D*self is a matrix over ZZ

        EXAMPLES:
            sage: a = matrix(QQ,3,[-2/7, -1/4, -2, 0, 1/7, 1, 0, 1/2, 1/5],sparse=True)
            sage: a.denominator()
            140
            sage: a._clear_denom()
            ([ -40  -35 -280]
            [   0   20  140]
            [   0   70   28], 140)
        """
        cdef Integer D
        cdef Py_ssize_t i, j
        cdef Matrix_integer_sparse A
        cdef mpz_t t
        cdef mpq_vector* v

        D = Integer()
        self.mpz_denom(D.value)

        MZ = sage.matrix.matrix_space.MatrixSpace(ZZ, self._nrows, self._ncols, sparse=True)
        A = MZ.zero_matrix()

        mpz_init(t)
        _sig_on
        for i from 0 <= i < self._nrows:
            v = &(self._matrix[i])
            for j from 0 <= j < v.num_nonzero:
                mpz_divexact(t, D.value, mpq_denref(v.entries[j]))
                mpz_mul(t, t, mpq_numref(v.entries[j]))
                mpz_vector_set_entry(&(A._matrix[i]), v.positions[j], t)
        _sig_off
        mpz_clear(t)
        A._initialized = 1
        return A, D


    ################################################
    # Echelon form
    ################################################
    def echelonize(self, height_guess=None, proof=True, **kwds):
        """
        INPUT:
            height_guess, proof, **kwds -- all passed to the multimodular algorithm; ignored
                                           by the p-adic algorithm.

        OUTPUT:
            matrix -- the reduced row echelon for of self.

        ALGORITHM: a multimodular algorithm.

        EXAMPLES:
            sage: a = matrix(QQ, 4, range(16), sparse=True); a[0,0] = 1/19; a[0,1] = 1/5; a
            [1/19  1/5    2    3]
            [   4    5    6    7]
            [   8    9   10   11]
            [  12   13   14   15]
            sage: a.echelonize(); a
            [      1       0       0 -76/157]
            [      0       1       0  -5/157]
            [      0       0       1 238/157]
            [      0       0       0       0]
        """

        x = self.fetch('in_echelon_form')
        if not x is None: return  # already known to be in echelon form
        self.check_mutability()
        self.clear_cache()

        pivots = self._echelonize_multimodular(height_guess, proof, **kwds)

        self.cache('in_echelon_form', True)
        self.cache('pivots', pivots)


    def echelon_form(self, algorithm='default',
                     height_guess=None, proof=True, **kwds):
        """
        INPUT:
            height_guess, proof, **kwds -- all passed to the multimodular algorithm; ignored
                                           by the p-adic algorithm.

        OUTPUT:
            self is no in reduced row echelon form.

        EXAMPLES:
            sage: a = matrix(QQ, 4, range(16), sparse=True); a[0,0] = 1/19; a[0,1] = 1/5; a
            [1/19  1/5    2    3]
            [   4    5    6    7]
            [   8    9   10   11]
            [  12   13   14   15]
            sage: a.echelon_form()
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

        E = self._echelon_form_multimodular(height_guess, proof=proof)

        self.cache(label, E)
        self.cache('pivots', E.pivots())
        return E

    # Multimodular echelonization algorithms
    def _echelonize_multimodular(self, height_guess=None, proof=True, **kwds):
        cdef Matrix_rational_sparse E
        E = self._echelon_form_multimodular(height_guess, proof=proof, **kwds)
        # Get rid of self's data
        self._dealloc()

        # Change self's data to point to E's.
        self._matrix = E._matrix

        # Make sure that E's destructure doesn't delete self's data.
        E._matrix = NULL
        E._initialized = False


    def _echelon_form_multimodular(self, height_guess=None, proof=True):
        """
        Returns reduced row-echelon form using a multi-modular
        algorithm.  Does not change self.

        INPUT:
            height_guess -- integer or None
            proof -- boolean (default: True)
        """
        import misc
        return misc.matrix_rational_echelon_form_multimodular(self,
                                 height_guess=height_guess, proof=proof)


    def set_row_to_multiple_of_row(self, i, j, s):
        """
        Set row i equal to s times row j.

        EXAMPLES:
            sage: a = matrix(QQ,2,3,range(6), sparse=True); a
            [0 1 2]
            [3 4 5]
            sage: a.set_row_to_multiple_of_row(1,2,-3)
            sage: a
            [ 0  1  2]
            [ 0 -3 -6]
        """
        self.check_mutability()
        cdef Rational _s
        _s = Rational(s)
        mpq_vector_scalar_multiply(&self._matrix[i], &self._matrix[j], _s.value)

    def dense_matrix(self):
        import misc
        return misc.matrix_rational_sparse__dense_matrix(self)

    def _set_row_to_negative_of_row_of_A_using_subset_of_columns(self, Py_ssize_t i, Matrix A,
                                                                 Py_ssize_t r, cols):
        """
        Set row i of self to -(row r of A), but where we only take
        the given column positions in that row of A:

        INPUT:
            i -- integer, index into the rows of self
            A -- a sparse matrix
            r -- integer, index into rows of A
            cols -- a *sorted* list of integers.

        EXAMPLES:
            sage: a = matrix(QQ['x'],2,3,range(6), sparse=True); a
            [0 1 2]
            [3 4 5]
            sage: a._set_row_to_negative_of_row_of_A_using_subset_of_columns(0,a,1,[1,2])
            sage: a
            [-4 -5  2]
            [ 3  4  5]
        """
        # this function exists just because it is useful for modular symbols presentations.
        self.check_mutability()

        cdef Matrix_rational_sparse _A
        _A = A

        cdef Py_ssize_t l
        l = 0

        cdef mpq_vector *v, *w
        v = &self._matrix[i]
        w = &_A._matrix[r]
        cdef mpq_t x
        mpq_init(x)
        for k in cols:
            mpq_vector_get_entry(&x, w, k)
            mpq_vector_set_entry(v, l, x)
            l += 1
        mpq_clear(x)
