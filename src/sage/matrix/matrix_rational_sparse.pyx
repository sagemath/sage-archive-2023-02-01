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

include '../modules/vector_rational_sparse_h.pxi'
include '../modules/vector_rational_sparse_c.pxi'
include '../ext/stdsage.pxi'

from sage.rings.rational cimport Rational
from sage.rings.integer  cimport Integer
from matrix cimport Matrix

from sage.rings.integer_ring import ZZ
import sage.matrix.matrix_space

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
        self._matrix = <mpq_vector*> sage_malloc(parent.nrows()*sizeof(mpq_vector))
        if self._matrix == NULL:
            raise MemoryError, "error allocating sparse matrix"

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._initialized:
            for i from 0 <= i < self._nrows:
                clear_mpq_vector(&self._matrix[i])
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

        # set the parent, nrows, ncols, etc.
        matrix_sparse.Matrix_sparse.__init__(self, parent)

        # initialize the rows
        for i from 0 <= i < parent.nrows():
            init_mpq_vector(&self._matrix[i], self._ncols, 0)

        # record that rows have been initialized
        self._initialized = True

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


    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_c_impl
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    ########################################################################

##     cdef int mpz_denom(self, mpz_t d) except -1:
##         mpz_set_si(d,1)
##         cdef Py_ssize_t i, j
##         cdef mpq_vector *v

##         _sig_on
##         for i from 0 <= i < self._nrows:
##             v = self._matrix[i]
##             for j from 0 <= j < v.num_nonzero:
##                 mpz_lcm(d, d, mpq_denref(v.entries[j]))
##         _sig_off
##         return 0

##     def _clear_denom(self):
##         """
##         INPUT:
##             self -- a matrix

##         OUTPUT:
##             D*self, D

##         The product D*self is a matrix over ZZ
##         """
##         cdef Integer D
##         cdef Py_ssize_t i, j
##         cdef Matrix_integer_sparse A
##         cdef mpz_t t

##         D = Integer()
##         self.mpz_denom(D.value)

##         MZ = sage.matrix.matrix_space.MatrixSpace(ZZ, self._nrows, self._ncols, sparse=True)
##         A = MZ.zero_matrix()

##         mpz_init(t)
##             _sig_on
##         for i from 0 <= i < self._nrows:
##             v = self._matrix[i]
##             for j from 0 <= j < v.num_nonzero:
##                 mpz_divexact(t, D.value, mpq_denref(v.entries[j]))
##                 mpz_mul(t, t, mpq_numref(v.entries[j]))
##                 mpz_vector_set_entry(&A._matrix[i], v.positions[j], t)
##         _sig_off
##         mpz_clear(t)
##         A._initialized = 1
##         return A, D


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
