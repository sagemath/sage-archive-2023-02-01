"""
Sparse integer matrices.

AUTHORS:
    -- William Stein (2007-02-21)
    -- Soroosh Yazdani (2007-02-21)

TESTS:
    sage: a = matrix(ZZ,2,range(4), sparse=True)
    sage: loads(dumps(a)) == a
    True
    sage: Matrix(ZZ,0,0,sparse=True).inverse()
    []
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
include '../modules/vector_modn_sparse_h.pxi'
include '../modules/vector_modn_sparse_c.pxi'

include '../ext/stdsage.pxi'

from sage.rings.integer  cimport Integer
from matrix cimport Matrix

from matrix_modn_sparse cimport Matrix_modn_sparse
from sage.structure.element cimport ModuleElement, RingElement, Element, Vector

import matrix_space

from sage.rings.integer_ring import ZZ
from sage.rings.integer_mod_ring import IntegerModRing


cdef class Matrix_integer_sparse(matrix_sparse.Matrix_sparse):

    ########################################################################
    # LEVEL 1 functionality
    #   * __cinit__
    #   * __dealloc__
    #   * __init__
    #   * set_unsafe
    #   * get_unsafe
    #   * __richcmp__    -- always the same
    #   * __hash__       -- always simple
    ########################################################################
    def __cinit__(self, parent, entries, copy, coerce):
        self._initialized = False
        # set the parent, nrows, ncols, etc.
        matrix_sparse.Matrix_sparse.__init__(self, parent)

        self._matrix = <mpz_vector*> sage_malloc(parent.nrows()*sizeof(mpz_vector))
        if self._matrix == NULL:
            raise MemoryError, "error allocating sparse matrix"

        # initialize the rows
        for i from 0 <= i < parent.nrows():
            mpz_vector_init(&self._matrix[i], self._ncols, 0)
        # record that rows have been initialized
        self._initialized = True

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._initialized:
            for i from 0 <= i < self._nrows:
                mpz_vector_clear(&self._matrix[i])
        sage_free(self._matrix)

    def __init__(self, parent, entries, copy, coerce):
        """
        Create a sparse matrix over the integers.

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
        cdef Integer z
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
                    mpz_vector_set_entry(&self._matrix[i], j, z.value)

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
                        mpz_vector_set_entry(&self._matrix[i], j, z.value)
                    k = k + 1

        else:

            # fill in entries in the scalar case
            z = Integer(entries)
            if z == 0:
                return
            if self._nrows != self._ncols:
                raise TypeError, "matrix must be square to initialize with a scalar."
            for i from 0 <= i < self._nrows:
                mpz_vector_set_entry(&self._matrix[i], i, z.value)


    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        mpz_vector_set_entry(&self._matrix[i], j, (<Integer> x).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef Integer x
        x = Integer()
        mpz_vector_get_entry(&x.value, &self._matrix[i], j)
        return x

    def __richcmp__(Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)

    def __hash__(self):
        return self._hash()


    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    #   * cdef _add_
    #   * cdef _sub_
    #   * cdef _mul_
    #   * cdef _cmp_c_impl
    #   * __neg__
    #   * __invert__
    #   * __copy__
    #   * _multiply_classical
    #   * _matrix_times_matrix_
    #   * _list -- list of underlying elements (need not be a copy)
    #   * x _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################
    # def _pickle(self):
    # def _unpickle(self, data, int version):   # use version >= 0
    # cpdef ModuleElement _add_(self, ModuleElement right):
    # cdef _mul_(self, Matrix right):
    # cdef int _cmp_c_impl(self, Matrix right) except -2:
    # def __neg__(self):
    # def __invert__(self):
    # def __copy__(self):
    # def _multiply_classical(left, matrix.Matrix _right):
    # def _list(self):

    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        EXAMPLES:
            sage: a = matrix(QQ,2,range(6))
            sage: (3/4) * a
            [   0  3/4  3/2]
            [ 9/4    3 15/4]
        """
        cdef Py_ssize_t i, j
        cdef mpz_vector* self_row, *M_row
        cdef Matrix_integer_sparse M
        cdef Integer _x
        _x = Integer(right)
        M = Matrix_integer_sparse.__new__(Matrix_integer_sparse, self._parent, None, None, None)
        for i from 0 <= i < self._nrows:
            self_row = &self._matrix[i]
            M_row = &M._matrix[i]
            mpz_vector_scalar_multiply(M_row, self_row, _x.value)
            #for j from 0 <= j < self._matrix[i].num_nonzero:
            #    mpz_vector_set_entry(M_row, self_row.positions[j], self_row.entries[j])
            #    mpz_mul(M_row.entries[j], M_row.entries[j], _x.value)
        return M

    cpdef ModuleElement _add_(self, ModuleElement right):
        cdef Py_ssize_t i, j
        cdef mpz_vector* self_row, *M_row
        cdef Matrix_integer_sparse M

        M = Matrix_integer_sparse.__new__(Matrix_integer_sparse, self._parent, None, None, None)
        cdef mpz_t mul
        mpz_init_set_si(mul,1)
        for i from 0 <= i < self._nrows:
            mpz_vector_clear(&M._matrix[i])
            add_mpz_vector_init(&M._matrix[i], &self._matrix[i], &(<Matrix_integer_sparse>right)._matrix[i], mul)
        mpz_clear(mul)
        return M

    cpdef ModuleElement _sub_(self, ModuleElement right):
        cdef Py_ssize_t i, j
        cdef mpz_vector* self_row, *M_row
        cdef Matrix_integer_sparse M

        M = Matrix_integer_sparse.__new__(Matrix_integer_sparse, self._parent, None, None, None)
        cdef mpz_t mul
        mpz_init_set_si(mul,-1)
        for i from 0 <= i < self._nrows:
            mpz_vector_clear(&M._matrix[i])
            add_mpz_vector_init(&M._matrix[i], &self._matrix[i], &(<Matrix_integer_sparse>right)._matrix[i], mul)
        mpz_clear(mul)
        return M

    def _dict(self):
        """
        Unsafe version of the dict method, mainly for internal use.
        This may return the dict of elements, but as an *unsafe*
        reference to the underlying dict of the object.  It might
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
    #    * cdef _sub_
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    ########################################################################

    def _nonzero_positions_by_row(self, copy=True):
        """
        Returns the list of pairs (i,j) such that self[i,j] != 0.

        It is safe to change the resulting list (unless you give the option copy=False).

        EXAMPLE::
            sage: M = Matrix(ZZ, [[0,0,0,1,0,0,0,0],[0,1,0,0,0,0,1,0]], sparse=True); M
            [0 0 0 1 0 0 0 0]
            [0 1 0 0 0 0 1 0]
            sage: M.nonzero_positions()
            [(0, 3), (1, 1), (1, 6)]

        """
        x = self.fetch('nonzero_positions')
        if not x is None:
            if copy:
                return list(x)
            return x
        nzp = []
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._matrix[i].num_nonzero:
                nzp.append((i,self._matrix[i].positions[j]))
        self.cache('nonzero_positions', nzp)
        if copy:
            return list(nzp)
        return nzp

    def _mod_int(self, modulus):
        return self._mod_int_c(modulus)

    cdef _mod_int_c(self, mod_int p):
        cdef Py_ssize_t i, j
        cdef Matrix_modn_sparse res
        cdef mpz_vector* self_row
        cdef c_vector_modint* res_row
        res = Matrix_modn_sparse.__new__(Matrix_modn_sparse, matrix_space.MatrixSpace(
            IntegerModRing(p), self._nrows, self._ncols, sparse=False), None, None, None)
        for i from 0 <= i < self._nrows:
            self_row = &(self._matrix[i])
            res_row = &(res.rows[i])
            for j from 0 <= j < self_row.num_nonzero:
                set_entry(res_row, self_row.positions[j], mpz_fdiv_ui(self_row.entries[j], p))
        return res


    def rational_reconstruction(self, N):
        """
        Use rational reconstruction to lift self to a matrix over the
        rational numbers (if possible), where we view self as a matrix
        modulo N.

        EXAMPLES:
            sage: A = matrix(ZZ, 3, 4, [(1/3)%500, 2, 3, (-4)%500, 7, 2, 2, 3, 4, 3, 4, (5/7)%500], sparse=True)
            sage: A.rational_reconstruction(500)
            [1/3   2   3  -4]
            [  7   2   2   3]
            [  4   3   4 5/7]
        """
        import misc
        return misc.matrix_integer_sparse_rational_reconstruction(self, N)

    def _linbox_sparse(self):
        cdef Py_ssize_t i, j
        v = ['%s %s M'%(self._nrows, self._ncols)]
        d = self._dict()
        for ij, x in d.iteritems():
            v.append('%s %s %s'%(ij[0]+1,ij[1]+1,x))
        v.append('0 0 0\n')
        return '\n'.join(v)

    def right_kernel(self, algorithm='padic', LLL=False, proof=None, echelonize=True):
        r"""
        Return the right kernel of this matrix, as a module over the
        integers.  This is the saturated ZZ-module spanned by all the
        column vectors v such that self*\v = 0.

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
            sage: M = MatrixSpace(ZZ,2,4,sparse=True)(range(8))
            sage: M.right_kernel()
            Free module of degree 4 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -3  2]
            [ 0  1 -2  1]
        """
        return self.dense_matrix().right_kernel(algorithm, LLL, proof, echelonize)

    def elementary_divisors(self, algorithm='pari'):
        """
        Return the elementary divisors of self, in order.

        The elementary divisors are the invariants of the finite
        abelian group that is the cokernel of *left* multiplication by
        this matrix.  They are ordered in reverse by divisibility.

        INPUT:
            self -- matrix
            algorithm -- (default: 'pari')
                 'pari': works robustly, but is slower.
                 'linbox' -- use linbox (currently off, broken)

        OUTPUT:
            list of integers

        EXAMPLES:
            sage: matrix(3, range(9),sparse=True).elementary_divisors()
            [1, 3, 0]
            sage: M = matrix(ZZ, 3, [1,5,7, 3,6,9, 0,1,2], sparse=True)
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
        return self.dense_matrix().elementary_divisors(algorithm=algorithm)
