"""
Matrices over Cyclotomic Fields

The underlying matrix for a Matrix_cyclo_dense object is stored as
follows: given an n x m matrix over a cyclotomic field of degree d, we
store a d x (nm) matrix over QQ, each column of which corresponds to
an element of the original matrix. This can be retrieved via the
_rational_matrix method. Here is an example illustrating this:

EXAMPLES:
    sage: F.<zeta> = CyclotomicField(5)
    sage: M = Matrix(F, 2, 3, [zeta, 3, zeta**4+5, (zeta+1)**4, 0, 1])
    sage: M
    [                        zeta                            3  -zeta^3 - zeta^2 - zeta + 4]
    [3*zeta^3 + 5*zeta^2 + 3*zeta                            0                            1]

    sage: M._rational_matrix()
    [ 0  3  4  0  0  1]
    [ 1  0 -1  3  0  0]
    [ 0  0 -1  5  0  0]
    [ 0  0 -1  3  0  0]


AUTHORS:
   * William Stein
   * Craig Citro
"""

######################################################################
#       Copyright (C) 2008 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
######################################################################
from sage.structure.element cimport ModuleElement, RingElement, Element, Vector
from matrix_space import MatrixSpace
from constructor import matrix
from sage.rings.rational_field import QQ, ZZ
from sage.rings.arith import previous_prime, binomial
from matrix cimport Matrix
import matrix_dense
from matrix_integer_dense import _lift_crt
from matrix_modn_dense import _matrix_from_rows_of_matrices
from sage.misc.misc import verbose
import math

from sage.rings.integer cimport Integer

from sage.structure.element cimport Matrix as baseMatrix

from misc import matrix_integer_dense_rational_reconstruction

from sage.ext.multi_modular import MAX_MODULUS

from sage.structure.proof.proof import get_flag as get_proof_flag

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

# parameters for tuning
echelon_primes_increment = 15
echelon_verbose_level = 1

cdef class Matrix_cyclo_dense(matrix_dense.Matrix_dense):
    ########################################################################
    # LEVEL 1 functionality
    # x * __new__
    # x * __dealloc__     (not needed)
    # x * __init__
    # x * set_unsafe
    # x * get_unsafe
    # x * _pickle
    # x * _unpickle
    ########################################################################

    def __new__(self, parent, entries, coerce, copy):
        """
        Create a new dense cyclotomic matrix.

        INPUT:
            parent -- a matrix space over a cyclotomic field
            entries -- a list of entries or scalar
            coerce -- bool; if true entries are coerced to base ring
            copy -- bool; ignored due to underlying data structure

        EXAMPLES:
            sage: from sage.matrix.matrix_cyclo_dense import Matrix_cyclo_dense
            sage: A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, MatrixSpace(CyclotomicField(3),2), [0,1,2,3], True, True)
            sage: type(A)
            <type 'sage.matrix.matrix_cyclo_dense.Matrix_cyclo_dense'>

        Note that the entries of A haven't even been set yet above; that doesn't
        happen until init is called:
            sage: A[0,0]
            Traceback (most recent call last):
            ...
            AttributeError: 'NoneType' object has no attribute 'column'
        """
        Matrix.__init__(self, parent)
        self._degree = self._base_ring.degree()

    # This is not necessary, since we do not (yet) explicitly allocate
    # any memory.
    #def __dealloc__(self):
    #    pass

    def __init__(self, parent, entries, copy, coerce):
        """
        Initialize a newly created cyclotomic matrix.

        INPUT:
            parent -- a matrix space over a cyclotomic field
            entries -- a list of entries or scalar
            coerce -- bool; if true entries are coerced to base ring
            copy -- bool; ignored due to underlying data structure

        EXAMPLES:
        This function is called implicitly when you create new cyclotomic
        dense matrices.
            sage: W.<a> = CyclotomicField(100)
            sage: A = matrix(2, 3, [1, 1/a, 1-a,a, -2/3*a, a^19])
            sage: A
            [                        1 -a^39 + a^29 - a^19 + a^9                    -a + 1]
            [                        a                    -2/3*a                      a^19]

        TESTS:
        We call __init__ explicitly below.
            sage: from sage.matrix.matrix_cyclo_dense import Matrix_cyclo_dense
            sage: A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, MatrixSpace(CyclotomicField(3),2), [0,1,2,3], True, True)
            sage: A.__init__(MatrixSpace(CyclotomicField(3),2), [0,1,2,3], True, True)
            sage: A
            [0 1]
            [2 3]

        """
        z = None
        if entries == 0:
            pass
        elif isinstance(entries, list):
            # This code could be made much faster using Cython, etc.
            if coerce:
                K = parent.base_ring()
                entries = [K(a) for a in entries]
            entries = sum([a.list() for a in entries], [])
        else:
            K = self._base_ring
            z = K(entries)
            entries = 0

        self._matrix = Matrix_rational_dense(MatrixSpace(QQ, self._nrows*self._ncols, self._degree),
                                            entries, copy=False, coerce=False).transpose()
        # This could also be made much faster.
        if z is not None:
            for i in range(self._nrows):
                self.set_unsafe(i,i,z)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        """
        Set the ij-th entry of self to value.

        WARNING: As the name suggests, there is no bound or type
        checking.  Expect errors, segfaults, etc., if you give bad
        input!!  This is for internal use only.

        EXAMPLES:
            sage: W.<a> = CyclotomicField(5)
            sage: A = matrix(2, 3, [1, 1/a, 1-a,a, -2/3*a, a^19])

        This indirectly calls set_unsafe:
            sage: A[0,0] = 109308420
            sage: A
            [         109308420 -a^3 - a^2 - a - 1             -a + 1]
            [                 a             -2/3*a -a^3 - a^2 - a - 1]
        """
        # The i,j entry is the (i * self._ncols + j)'th column.
        # TODO: This could be made way faster via direct access to the
        # underlying matrix
        cdef Py_ssize_t k, c
        v = value.list()
        c = i * self._ncols + j
        for k from 0 <= k < self._degree:
            self._matrix.set_unsafe(k, c, v[k])

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Get the ij-th of self.

        WARNING: As the name suggests, expect segfaults if i,j are out
        of bounds!! This is for internal use only.

        EXAMPLES:
            sage: W.<a> = CyclotomicField(5)
            sage: A = matrix(2, 3, [9939208341, 1/a, 1-a,a, -2/3*a, a^19])

        This implicitly calls get_unsafe:
            sage: A[0,0]
            9939208341
        """
        # The i,j entry is the (i * self._ncols + j)'th column.
        # TODO: This could be made way faster via direct access to the
        # underlying matrix
        return self.base_ring()(self._matrix.column(i*self._ncols + j).list())

    def _pickle(self):
        """
        Used for pickling matrices. This function returns the
        underlying data and pickle version.

        OUTPUT:
            data -- output of pickle
            version -- int

        EXAMPLES:
            sage: K.<z> = CyclotomicField(3)
            sage: w = matrix(K, 3, 3, [0, -z, -2, -2*z + 2, 2*z, z, z, 1-z, 2+3*z])
            sage: w._pickle()
            (('0 0 -2 2 0 0 0 1 2 0 -1 0 -2 2 1 1 -1 3', 0), 0)
        """
        data = self._matrix._pickle()
        return data, 0

    def _unpickle(self, data, int version):
        """
        Called when unpickling matrices.

        INPUT:
            data -- a string
            version -- int

        This modifies self.

        EXAMPLES:
            sage: K.<z> = CyclotomicField(3)
            sage: w = matrix(K, 3, 3, [0, -z, -2, -2*z + 2, 2*z, z, z, 1-z, 2+3*z])
            sage: data, version = w._pickle()
            sage: k = w.parent()(0)
            sage: k._unpickle(data, version)
            sage: k == w
            True
        """
        self.check_mutability()
        if version == 0:
            self._matrix = Matrix_rational_dense(MatrixSpace(QQ, self._degree, self._nrows*self._ncols), None, False, False)
            self._matrix._unpickle(*data)  # data is (data, matrix_QQ_version)
        else:
            raise RuntimeError, "unknown matrix version (=%s)"%version

    ########################################################################
    # LEVEL 2 functionality
    # x * cdef _add_c_impl
    # x * cdef _sub_c_impl
    #   * cdef _mul_c_impl
    # x * cdef _lmul_c_impl    -- scalar multiplication
    # x * cdef _cmp_c_impl
    # x * __neg__
    #   * __invert__
    # x * __copy__
    #   * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################

    cdef ModuleElement _add_c_impl(self, ModuleElement right):
        """
        Return the sum of two dense cyclotomic matrices.

        INPUT:
            self, right -- dense cyclotomic matrices with the same
                           parents
        OUTPUT:
            a dense cyclotomic matrix

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(2, 3, [1,z,z^2,z^3,z^4,2/3*z]); B = matrix(2, 3, [-1,2*z,3*z^2,5*z+1,z^4,1/3*z])
            sage: A + B
            [                       0                      3*z                    4*z^2]
            [           z^3 + 5*z + 1 -2*z^3 - 2*z^2 - 2*z - 2                        z]

        Verify directly that the above output is correct:
            sage: (A+B).list() == [A.list()[i]+B.list()[i] for i in range(6)]
            True
        """
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)
        # Fix the redundancy here.
        A._matrix = self._matrix + (<Matrix_cyclo_dense>right)._matrix
        return A

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        """
        Return the difference of two dense cyclotomic matrices.

        INPUT:
            self, right -- dense cyclotomic matrices with the same
                           parent
        OUTPUT:
            a dense cyclotomic matrix

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(2, 3, [1,z,z^2,z^3,z^4,2/3*z]); B = matrix(2, 3, [-1,2*z,3*z^2,5*z+1,z^4,1/3*z])
            sage: A - B
            [            2            -z        -2*z^2]
            [z^3 - 5*z - 1             0         1/3*z]

        Verify directly that the above output is correct:
            sage: (A-B).list() == [A.list()[i]-B.list()[i] for i in range(6)]
            True
        """
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)
        A._matrix = self._matrix - (<Matrix_cyclo_dense>right)._matrix
        return A

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        """
        Multiply a dense cyclotomic matrix by a scalar.

        INPUT:
            self -- dense cyclotomic matrix
            right --- scalar in the base cyclotomic field

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(2, 3, [1,z,z^2,z^3,z^4,2/3*z])
            sage: (1 + z/3)*A
            [                      1/3*z + 1                     1/3*z^2 + z                   1/3*z^3 + z^2]
            [2/3*z^3 - 1/3*z^2 - 1/3*z - 1/3            -z^3 - z^2 - z - 2/3                 2/9*z^2 + 2/3*z]

        Verify directly that the above output is correct:
            sage: ((1+z/3)*A).list() == [(1+z/3)*x for x in A.list()]
            True
        """
        if right == 1:
            return self
        elif right == 0:
            return self.parent()(0)

        # Create a new matrix object but with the _matrix attribute not initialized:
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense,
                                               self.parent(), None, None, None)

        if right.polynomial().degree() == 0:
            # multiplication by a rational number
            A._matrix = self._matrix._lmul_c_impl(right)
        else:
            # Multiply by nontrivial element of the cyclotomic number field
            # We do this by finding the matrix of this element, then left
            # multiplying by it.  We have to tweak the matrix some to
            # get the right basis, etc.
            T = right.matrix().transpose()
            A._matrix = T * self._matrix
        return A

    cdef baseMatrix _matrix_times_matrix_c_impl(self, baseMatrix right):
        """
        Return the product of two cyclotomic dense matrices.

        INPUT:
            self, right -- cyclotomic dense matrices with compatible
                           parents (same base ring, and compatible
                           dimensions for matrix multiplication).

        OUTPUT:
            cyclotomic dense matrix

        ALGORITHM:
            Use a multimodular algorithm that involves multiplying the
            two matrices modulo split primes.

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(3, 3, [1,z,z^2,z^3,z^4,2/3*z,-3*z,z,2+z]); B = matrix(3, 3, [-1,2*z,3*z^2,5*z+1,z^4,1/3*z,2-z,3-z,5-z])
            sage: A*B
            [        -z^3 + 7*z^2 + z - 1       -z^3 + 3*z^2 + 2*z + 1              -z^3 + 25/3*z^2]
            [-2*z^3 - 5/3*z^2 + 1/3*z + 4           -z^3 - 8/3*z^2 - 2     -2/3*z^2 + 10/3*z + 10/3]
            [             4*z^2 + 4*z + 4               -7*z^2 + z + 7  -9*z^3 - 2/3*z^2 + 3*z + 10]

        Verify that the answer above is consistent with what the
        generic sparse matrix multiply gives (which is a different
        implementation).
            sage: A*B == A.sparse_matrix()*B.sparse_matrix()
            True

            sage: N1 = Matrix(CyclotomicField(6), 1, [1])
            sage: cf6 = CyclotomicField(6) ; z6 = cf6.0
            sage: N2 = Matrix(CyclotomicField(6), 1, 5, [0,1,z6,-z6,-z6+1])
            sage: N1*N2
            [         0          1      zeta6     -zeta6 -zeta6 + 1]
            sage: N1 = Matrix(CyclotomicField(6), 1, [-1])
            sage: N1*N2
            [        0        -1    -zeta6     zeta6 zeta6 - 1]
        """
        A, denom_self = self._matrix._clear_denom()
        B, denom_right = (<Matrix_cyclo_dense>right)._matrix._clear_denom()

        # conservative but correct estimate
        bound = A.height() * B.height() * self._ncols

        n = self._base_ring._n()
        p = previous_prime(MAX_MODULUS)
        prod = 1
        v = []
        while prod <= bound:
            while (n >= 2 and p % n != 1) or denom_self % p == 0 or denom_right % p == 0:
                if p == 2:
                    raise RuntimeError, "we ran out of primes in matrix multiplication."
                p = previous_prime(p)
            prod *= p
            Amodp, _ = self._reductions(p)
            Bmodp, _ = right._reductions(p)
            _,     S = self._reduction_matrix(p)
            X = _matrix_from_rows_of_matrices([Amodp[i] * Bmodp[i] for i
                                               in range(len(Amodp))])
            v.append(S*X)
            p = previous_prime(p)
        M = matrix(ZZ, self._base_ring.degree(), self._nrows*right.ncols())
        _lift_crt(M, v)
        d = denom_self * denom_right
        if d == 1:
            M = M.change_ring(QQ)
        else:
            M = (1/d)*M
        cdef Matrix_cyclo_dense C = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense,
                    MatrixSpace(self._base_ring, self._nrows, right.ncols()),
                                                               None, None, None)
        C._matrix = M
        return C

    def __richcmp__(Matrix self, right, int op):
        """
        Compare a matrix with something else. This immediately calls
        a base class _richcmp.

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,z,-z,1+z/2])

        These implicitly call richcmp:
            sage: A == 5
            False
            sage: A < 100
            True
        """
        return self._richcmp(right, op)

    cdef long _hash(self) except -1:
        """
        Return hash of this matrix.

        EXAMPLES:
        This is called implicitly by the hash function.
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,z,-z,1+z/2])

        The matrix must be immutable.
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()

        Yes, this works:
            sage: hash(A)
            -25
        """
        return self._matrix._hash()

    def __hash__(self):
        """
        Return hash of an immutable matrix. Raise a TypeError if input
        matrix is mutable.

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,2/3*z+z^2,-z,1+z/2])
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()
            sage: A.__hash__()
            -18
        """
        if self._mutability._is_immutable:
            return self._hash()
        else:
            raise TypeError, "mutable matrices are unhashable"

    cdef int _cmp_c_impl(self, Element right) except -2:
        """
        Implements comparison of two cyclotomic matrices with
        identical parents.

        INPUT:
            self, right -- matrices with same parent
        OUTPUT:
            int; either -1, 0, or 1

        EXAMPLES:
        This function is called implicitly when comparisons with matrices
        are done or the cmp function is used.
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,2/3*z+z^2,-z,1+z/2])
            sage: cmp(A,A)
            0
            sage: cmp(A,2*A)
            -1
            sage: cmp(2*A,A)
            1
        """
        return self._matrix._cmp_c_impl((<Matrix_cyclo_dense>right)._matrix)

    def __copy__(self):
        """
        Make a copy of this matrix.

        EXAMPLES:
        We create a cyclotomic matrix.
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,2/3*z+z^2,-z,1+z/2])

        We make a copy of A.
            sage: C = A.copy()

        We make another reference to A.
            sage: B = A

        Changing this reference changes A itself:
            sage: B[0,0] = 10
            sage: A[0,0]
            10

        Changing the copy does not change A.
            sage: C[0,0] = 20
            sage: C[0,0]
            20
            sage: A[0,0]
            10
        """
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)
        A._matrix = self._matrix.__copy__()
        return A

    def __neg__(self):
        """
        Return the negative of this matrix.

        OUTPUT:
            matrix

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,2/3*z+z^2,-z,1+z/2])
            sage: -A
            [          -1 -z^2 - 2/3*z]
            [           z   -1/2*z - 1]
            sage: A.__neg__()
            [          -1 -z^2 - 2/3*z]
            [           z   -1/2*z - 1]
        """
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)
        A._matrix = self._matrix.__neg__()
        return A


    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    #    * Specialized echelon form
    ########################################################################
    def set_immutable(self):
        """
        Change this matrix so that it is immutable.

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,2/3*z+z^2,-z,1+z/2])
            sage: A[0,0] = 10
            sage: A.set_immutable()
            sage: A[0,0] = 20
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead (use self.copy()).

        Note that there is no function to set a matrix to be mutable
        again, since such a function would violate the whole point.
        Instead make a copy, which is always mutable by default.
            sage: A.set_mutable()
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.matrix.matrix_cyclo_dense.Matrix_cyclo_dense' object has no attribute 'set_mutable'
            sage: B = A.copy()
            sage: B[0,0] = 20
            sage: B[0,0]
            20
        """
        self._matrix.set_immutable()
        matrix_dense.Matrix_dense.set_immutable(self)

    def _rational_matrix(self):
        """
        Return the underlying rational matrix corresponding to self.

        EXAMPLES:
            sage: Matrix(CyclotomicField(7),4,4,range(16))._rational_matrix()
            [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            sage: Matrix(CyclotomicField(7),4,4,[CyclotomicField(7).gen(0)**i for i in range(16)])._rational_matrix()
            [ 1  0  0  0  0  0 -1  1  0  0  0  0  0 -1  1  0]
            [ 0  1  0  0  0  0 -1  0  1  0  0  0  0 -1  0  1]
            [ 0  0  1  0  0  0 -1  0  0  1  0  0  0 -1  0  0]
            [ 0  0  0  1  0  0 -1  0  0  0  1  0  0 -1  0  0]
            [ 0  0  0  0  1  0 -1  0  0  0  0  1  0 -1  0  0]
            [ 0  0  0  0  0  1 -1  0  0  0  0  0  1 -1  0  0]
        """
        return self._matrix

    def denominator(self):
        """
        Return the denominator of the entries of this matrix.

        OUTPUT:
            integer -- the smallest integer d so that d * self has
                       entries in the ring of integers

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [-2/7,2/3*z+z^2,-z,1+z/19]); A
            [       -2/7 z^2 + 2/3*z]
            [         -z  1/19*z + 1]
            sage: d = A.denominator(); d
            399
        """
        return self._matrix.denominator()

    def coefficient_bound(self):
        r"""
        Return an upper bound for the (complex) absolute values of all
        entries of self with respect to all embeddings.

        Use \code{self.height()} for a sharper bound.

        This is computed using just the Cauchy-Schwarz inequality, i.e.,
        we use the fact that
        $$ \left| \sum_i a_i\zeta^i \right| \leq \sum_i |a_i|, $$
        as $|\zeta| = 1$.

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1+z, 0, 9*z+7, -3 + 4*z]); A
            [  z + 1       0]
            [9*z + 7 4*z - 3]
            sage: A.coefficient_bound()
            16

        The above bound is just $9 + 7$, coming from the lower left entry.
        A better bound would be the following:
            sage: (A[1,0]).abs()
            12.9975436637560
        """
        cdef Py_ssize_t i, j

        bound = 0
        for i from 0 <= i < self._matrix._ncols:

            n = 0
            for j from 0 <= j < self._matrix._nrows:
                n += self._matrix[j, i].abs()
            if bound < n:
                bound = n

        return bound


    def height(self):
        r"""
        Return the height of self.

        If we let $a_{ij}$ be the $i,j$ entry of self, then we define
        the height of self to be
        $$
          \max_v \max_{i,j} |a_{ij}|_v,
        $$
        where $v$ runs over all complex embeddings of \code{self.base_ring()}.

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1+z, 0, 9*z+7, -3 + 4*z]); A
            [  z + 1       0]
            [9*z + 7 4*z - 3]
            sage: A.height()
            12.9975436637560
            sage: (A[1,0]).abs()
            12.9975436637560
        """
        cdef Py_ssize_t i, j

        emb = self._base_ring.complex_embeddings()

        ht = 0
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                t = max([ x.norm().sqrt() for x in [ f(self.get_unsafe(i,j)) for f in emb ] ])
                if t > ht:
                    ht = t

        return ht

    def randomize(self, density=1, num_bound=2, den_bound=2, distribution=None):
        r"""
        Randomize the entries of self.

        Choose rational numbers according to \code{distribution},
        whose numerators are bounded by \code{num_bound} and whose
        denominators are bounded by \code{den_bound}.

        EXAMPLES:
            sage: A = Matrix(CyclotomicField(5),2,2,range(4)) ; A
            [0 1]
            [2 3]
            sage: A.randomize()
            sage: A   # random output
            [       1/2*zeta5^2 + zeta5                        1/2]
            [        -zeta5^2 + 2*zeta5 -2*zeta5^3 + 2*zeta5^2 + 2]
        """
        self._cache = {}
        self._matrix.randomize(density, num_bound, den_bound, distribution)

    def _charpoly_bound(self):
        """
        Determine a bound for the coefficients of the characteristic
        polynomial of self. We use the bound in Lemma 2.1 of:

          Dumas, J-G. "Bounds on the coefficients of characteristic
          and minimal polynomials." J. Inequal. Pure Appl. Math. 8
          (2007), no. 2.

        This bound only applies for self._nrows >= 4.

        EXAMPLES:
            sage: A = Matrix(CyclotomicField(7),3,3,range(9))
            sage: A._charpoly_bound()
            2048
            sage: A.charpoly()
            x^3 + (-12)*x^2 + (-18)*x

        An example from the above paper, where our bound is sharp:
            sage: B = Matrix(CyclotomicField(7), 5,5, [1,1,1,1,1,1,1,-1,-1,-1,1,-1,1,-1,-1,1,-1,-1,1,-1,1,-1,-1,-1,1])
            sage: B._charpoly_bound()
            81
            sage: B.charpoly()
            x^5 + (-5)*x^4 + 40*x^2 + (-80)*x + 48
        """
        cdef Py_ssize_t i, j
        cdef float alpha, delta

        # should we even bother with this check, or just say in
        # the docstring that we assume it's square?
        if self._nrows != self._ncols:
            raise ValueError, "self must be square"

        # This is an approximation to 2^(5/6*log_2(5) - 2/3*log_2(6))
        alpha = 1.15799718800731
        # This is 2*e^(1-(2(7\gamma-4))/(13(3-2\gamma))), where \gamma
        # is Euler's constant.
        delta = 5.418236

        B = self.coefficient_bound()

        # this bound is only valid for n >= 4, use naive bounds
        # in other cases.
        # TODO: should charpoly just hardcode the return value for
        # self.nrows() < 4?
        if self._nrows > 3:
            return ZZ(int(math.ceil((alpha * self._nrows * B**2)**(self._nrows/2.0))))
        elif self._nrows == 3:
            return max(6*B**2, 4*B**3)
        elif self._nrows == 2:
            return 2*B**2
        else:
            return B

    def charpoly(self, var='x', algorithm="multimodular", proof=None):
        r"""
        Return the characteristic polynomial of self, as a polynomial
        over the base ring.

        INPUT:
            algorithm -- 'multimodular' (default): reduce modulo
                                        primes, compute charpoly mod
                                        p, and lift (very fast)
                         'pari': use pari (quite slow; comparable to
                                        Magma v2.14 though)
                         'hessenberg': put matrix in Hessenberg form
                                        (double dog slow)
            proof -- bool (default: None) proof flag determined by
                                          global linalg proof.

        OUTPUT:
            polynomial

        EXAMPLES:
            sage: K.<z> = CyclotomicField(5)
            sage: a = matrix(K, 3, [1,z,1+z^2, z/3,1,2,3,z^2,1-z])
            sage: f = a.charpoly(); f
            x^3 + (z - 3)*x^2 + (-16/3*z^2 - 2*z)*x - 2/3*z^3 + 16/3*z^2 - 5*z + 5/3
            sage: f(a)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: a.charpoly(algorithm='pari')
            x^3 + (z - 3)*x^2 + (-16/3*z^2 - 2*z)*x - 2/3*z^3 + 16/3*z^2 - 5*z + 5/3
            sage: a.charpoly(algorithm='hessenberg')
            x^3 + (z - 3)*x^2 + (-16/3*z^2 - 2*z)*x - 2/3*z^3 + 16/3*z^2 - 5*z + 5/3

            sage: Matrix(K, 1, [0]).charpoly()
            x
            sage: Matrix(K, 1, [5]).charpoly(var='y')
            y - 5

            sage: Matrix(CyclotomicField(13),3).charpoly()
            x^3
            sage: Matrix(CyclotomicField(13),3).charpoly()[2].parent()
            Cyclotomic Field of order 13 and degree 12
        """
        key = 'charpoly-%s-%s'%(algorithm,proof)
        f = self.fetch(key)
        if f is not None:
            return f.change_variable_name(var)

        if self.nrows() != self.ncols():
            raise TypeError, "self must be square"

        if self.is_zero():
            R = PolynomialRing(self.base_ring(), name=var)
            f = R.gen(0)**self.nrows()
            self.cache(key, f)
            return f

        if self.nrows() == 1:
            R = PolynomialRing(self.base_ring(), name=var)
            f = R.gen(0) - self[0,0]
            self.cache(key, f)
            return f

        if algorithm == 'multimodular':
            f = self._charpoly_multimodular(var, proof=proof)
        elif algorithm == 'pari':
            f = self._charpoly_over_number_field(var)
        elif algorithm == 'hessenberg':
            f = self._charpoly_hessenberg(var)
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm
        self.cache(key, f)
        return f

    def _charpoly_mod(self, p):
        """
        Return the characteristic polynomial of self*denom modulo all
        primes over $p$.

        This is used internally by the multimodular charpoly algorithm.

        INPUT:
            p -- a prime that splits completely

        OUTPUT:
            matrix over GF(p) whose columns correspond to the entries
            of all the charpolys of the reduction of self modulo all
            the primes over p.

        EXAMPLES:
            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1+z, 0, 9*z+7, -3 + 4*z]); A
            [  z + 1       0]
            [9*z + 7 4*z - 3]
            sage: A._charpoly_mod(11)
            [8 2 1]
            [1 6 0]
            [4 0 0]
            [0 0 0]
        """
        tm = verbose("Computing characteristic polynomial of cyclomotic matrix modulo %s."%p)
        # Reduce self modulo all primes over p
        R, denom = self._reductions(p)
        # Compute the characteristic polynomial of each reduced matrix
        F = [A.charpoly('x') for A in R]
        # Put the charpolys together as the rows of a mod-p matrix
        k = R[0].base_ring()
        S = matrix(k, len(F), self.nrows()+1, [f.list() for f in F])
        # multiply by inverse of reduction matrix to lift
        _, L = self._reduction_matrix(p)
        X = L * S
        # Now the columns of the matrix X define the entries of the
        # charpoly modulo p.
        verbose("Finished computing charpoly mod %s."%p, tm)
        return X

    def _charpoly_multimodular(self, var='x', proof=None):
        """
        Compute the characteristic polynomial of self using a
        multimodular algorithm.

        INPUT:
            proof -- bool (default: global flag); if False, compute
                     using primes $p_i$ until the lift modulo all
                     primes up to $p_i$ is the same as the lift modulo
                     all primes up to $p_{i+3}$ or the bound is
                     reached.

        EXAMPLES:
            sage: K.<z> = CyclotomicField(3)
            sage: A = matrix(3, [-z, 2*z + 1, 1/2*z + 2, 1, -1/2, 2*z + 2, -2*z - 2, -2*z - 2, 2*z - 1])
            sage: A._charpoly_multimodular()
            x^3 + (-z + 3/2)*x^2 + (17/2*z + 9/2)*x - 9/2*z - 23/2
            sage: A._charpoly_multimodular('T')
            T^3 + (-z + 3/2)*T^2 + (17/2*z + 9/2)*T - 9/2*z - 23/2
            sage: A._charpoly_multimodular('T', proof=False)
            T^3 + (-z + 3/2)*T^2 + (17/2*z + 9/2)*T - 9/2*z - 23/2

        TESTS:
        We test a degenerate case:
            sage: A = matrix(CyclotomicField(1),2,[1,2,3,4]); A.charpoly()
            x^2 + (-5)*x - 2
        """
        cdef Matrix_cyclo_dense A
        A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(),
                                       None, None, None)

        proof = get_proof_flag(proof, "linear_algebra")

        n = self._base_ring._n()
        p = previous_prime(MAX_MODULUS)
        prod = 1
        v = []
        #A, denom = self._matrix._clear_denom()
        # TODO: this might be stupidly slow
        denom = self._matrix.denominator()
        A._matrix = <Matrix_rational_dense>(denom*self._matrix)
        bound = A._charpoly_bound()
        L_last = 0
        while prod <= bound:
            while (n >= 2  and p % n != 1) or denom % p == 0:
                if p == 2:
                    raise RuntimeError, "we ran out of primes in multimodular charpoly algorithm."
                p = previous_prime(p)

            X = A._charpoly_mod(p)
            v.append(X)
            prod *= p
            p = previous_prime(p)

            # if we've used enough primes as determined by bound, or
            # if we've used 3 primes, we check to see if the result is
            # the same.
            if prod >= bound or (not proof  and  (len(v) % 3 == 0)):
                M = matrix(ZZ, self._base_ring.degree(), self._nrows+1)
                L = _lift_crt(M, v)
                if not proof and L == L_last:
                    break
                L_last = L

        # Now each column of L encodes a coefficient of the output polynomial,
        # with column 0 being the constant coefficient.
        K = self.base_ring()
        R = K[var]
        coeffs = [K(w.list()) for w in L.columns()]
        f = R(coeffs)

        # Rescale to account for denominator, if necessary
        if denom != 1:
            x = R.gen()
            f = f(x * denom) * (1 / (denom**f.degree()))

        return f


    def _reductions(self, p):
        """
        Compute the reductions modulo all primes over p of denom*self,
        where denom is the denominator of self.

        INPUT:
            p -- a prime that splits completely in the base cyclotomic field.

        OUTPUT:
            list -- of r distinct matrices modulo p, where r is
                    the degree of the cyclotomic base field.
            denom -- an integer

        EXAMPLES:
            sage: K.<z> = CyclotomicField(3)
            sage: w = matrix(K, 2, 3, [0, -z/5, -2/3, -2*z + 2, 2*z, z])
            sage: R, d = w._reductions(7)
            sage: R[0]
            [0 2 4]
            [1 1 4]
            sage: R[1]
            [0 1 4]
            [5 4 2]
            sage: d
            15
        """
        # Get matrix that defines the linear reduction maps modulo
        # each prime of the base ring over p.
        T, _ = self._reduction_matrix(p)
        # Clear denominator and get matrix over the integers suitable
        # for reduction.
        A, denom = self._matrix._clear_denom()
        # Actually reduce the matrix over the integers modulo the
        # prime p.
        B = A._mod_int(p)
        # Now multiply, which computes from B all the reductions of
        # self*denom modulo each of the primes over p.
        R = T * B
        # Finally compute the actual reductions by extracting them
        # from R (note that the rows of R define the reductions).
        ans = R._matrices_from_rows(self._nrows, self._ncols)
        return ans, denom

    def _reduction_matrix(self, p):
        """
        INPUT:
            p -- a prime that splits completely in the base field.

        OUTPUT:
            -- Matrix over GF(p) whose action from the left
               gives the map from O_K to GF(p) x ... x GF(p)
               given by reducing modulo all the primes over p.
            -- inverse of this matrix

        EXAMPLES:
            sage: K.<z> = CyclotomicField(3)
            sage: w = matrix(K, 2, 3, [0, -z/5, -2/3, -2*z + 2, 2*z, z])
            sage: A, B = w._reduction_matrix(7)
            sage: A
            [1 4]
            [1 2]
            sage: B
            [6 2]
            [4 3]

        The reduction matrix is used to calculate the reductions mod primes
        above p.
            sage: K.<z> = CyclotomicField(5)
            sage: A = matrix(K, 2, 2, [1, z, z^2+1, 5*z^3]); A
            [      1       z]
            [z^2 + 1   5*z^3]
            sage: T, S = A._reduction_matrix(11)
            sage: T * A._rational_matrix().change_ring(GF(11))
            [ 1  9  5  4]
            [ 1  5  4  9]
            [ 1  4  6  1]
            [ 1  3 10  3]

        The rows of this product are the (flattened) matrices mod each prime above p:
            sage: roots = [r for r, e in K.defining_polynomial().change_ring(GF(11)).roots()]; roots
            [9, 5, 4, 3]
            sage: [r^2+1 for r in roots]
            [5, 4, 6, 10]
            sage: [5*r^3 for r in roots]
            [4, 9, 1, 3]

        The reduction matrix is cached:
            sage: w._reduction_matrix(7) is w._reduction_matrix(7)
            True
        """
        cache = self.fetch('reduction_matrices')
        if cache is None:
            cache = {}
            self.cache('reduction_matrices', cache)
        try:
            return cache[p]
        except KeyError:
            pass
        K = self.base_ring()
        phi = K.defining_polynomial()
        from sage.rings.all import GF
        from constructor import matrix
        F = GF(p)
        aa = [a for a, _ in phi.change_ring(F).roots()]
        n = K.degree()
        if len(aa) != n:
            raise ValueError, "the prime p (=%s) must split completely but doesn't"%p
        T = matrix(F, n)
        for i in range(n):
            a = aa[i]
            b = 1
            for j in range(n):
                T[i,j] = b
                b *= a
        T.set_immutable()
        ans = (T, T**(-1))
        cache[p] = ans
        return ans

    def echelon_form(self, algorithm='multimodular'):
        """
        Find the echelon form of self, using the specified algorithm.

        The result is cached for each algorithm separately.

        EXAMPLES:
            sage: W.<z> = CyclotomicField(3)
            sage: A = matrix(W, 2, 3, [1+z, 2/3, 9*z+7, -3 + 4*z, z, -7*z]); A
            [  z + 1     2/3 9*z + 7]
            [4*z - 3       z    -7*z]
            sage: A.echelon_form()
            [                  1                   0  -192/97*z - 361/97]
            [                  0                   1 1851/97*z + 1272/97]
            sage: A.echelon_form(algorithm='classical')
            [                  1                   0  -192/97*z - 361/97]
            [                  0                   1 1851/97*z + 1272/97]

        We verify that the result is cached and that the caches are separate:
            sage: A.echelon_form() is A.echelon_form()
            True
            sage: A.echelon_form() is A.echelon_form(algorithm='classical')
            False

        TESTS:
            sage: W.<z> = CyclotomicField(13)
            sage: A = Matrix(W, 2,3, [10^30*(1-z)^13, 1, 2, 3, 4, z])
            sage: B = Matrix(W, 2,3, [(1-z)^13, 1, 2, 3, 4, z])
            sage: A.echelon_form() == A.echelon_form('classical')
            True
            sage: B.echelon_form() == B.echelon_form('classical')
            True

        A degenerate case with the degree 1 cyclotomic field:
            sage: A = matrix(CyclotomicField(1),2,3,[1,2,3,4,5,6]);
            sage: A.echelon_form()
            [ 1  0 -1]
            [ 0  1  2]

        A case that checks the bug in trac #3500.
            sage: cf4 = CyclotomicField(4) ; z4 = cf4.0
            sage: A = Matrix(cf4, 1, 2, [-z4, 1])
            sage: A.echelon_form()
            [    1 zeta4]
        """
        key = 'echelon_form-%s'%algorithm
        E = self.fetch(key)
        if E is not None:
            return E

        if self._nrows == 0:
            E = self.copy()
            self.cache(key, E)
            self.cache('pivots', [])
            return E

        if algorithm == 'multimodular':
            E = self._echelon_form_multimodular()
        elif algorithm == 'classical':
            E = (self*self.denominator())._echelon_classical()
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm

        self.cache(key, E)
        return E

    def _echelon_form_multimodular(self, num_primes=10, height_guess=None):
        """
        Use a multimodular algorithm to find the echelon form of self.

        INPUT:
            num_primes -- number of primes to work modulo
            height_guess -- guess for the height of the echelon form
                            of self

        OUTPUT:
            matrix in reduced row echelon form

        EXAMPLES:
            sage: W.<z> = CyclotomicField(3)
            sage: A = matrix(W, 2, 3, [1+z, 2/3, 9*z+7, -3 + 4*z, z, -7*z]); A
            [  z + 1     2/3 9*z + 7]
            [4*z - 3       z    -7*z]
            sage: A._echelon_form_multimodular(10)
            [                  1                   0  -192/97*z - 361/97]
            [                  0                   1 1851/97*z + 1272/97]

        TESTS:
        We test a degenerate case:
            sage: A = matrix(CyclotomicField(5),0); A
            []
            sage: A._echelon_form_multimodular(10)
            []
            sage: A.pivots()
            []
        """
        cdef int i
        cdef Matrix_cyclo_dense res

        verbose("entering _echelon_form_multimodular", level=echelon_verbose_level)

        denom = self._matrix.denominator()
        A = denom * self

        # This bound is chosen somewhat arbitrarily. Changing it affects the
        # runtime, not the correctness of the result.
        if height_guess is None:
            height_guess = (A.coefficient_bound()+100)*1000000

        # This is all setup to keep track of various data
        # in the loop below.
        p = previous_prime(MAX_MODULUS)
        found = 0
        prod = 1
        n = self._base_ring._n()
        height_bound = self._ncols * height_guess * A.coefficient_bound() + 1
        mod_p_ech_ls = []
        max_pivots = []
        is_square = self._nrows == self._ncols

        verbose("using height bound %s"%height_bound, level=echelon_verbose_level)

        while True:
            # Generate primes to use, and find echelon form
            # modulo those primes.
            while found < num_primes or prod <= height_bound:
                if (n == 1) or p%n == 1:
                    try:
                        mod_p_ech, piv_ls = A._echelon_form_one_prime(p)
                    except ValueError:
                        # This means that we chose a prime which divides
                        # the denominator of the echelon form of self, so
                        # just skip it and continue
                        p = previous_prime(p)
                        continue
                    # if we have the identity, just return it, and
                    # we're done.
                    if is_square and len(piv_ls) == self._nrows:
                        self.cache('pivots', range(self._nrows))
                        return self.parent().identity_matrix()
                    if piv_ls > max_pivots:
                        mod_p_ech_ls = [mod_p_ech]
                        max_pivots = piv_ls
                    elif piv_ls == max_pivots:
                        mod_p_ech_ls.append(mod_p_ech)

                    # add this to the list of primes
                    found += 1
                    prod *= p
                p = previous_prime(p)

            if found > num_primes:
                num_primes = found

            verbose("computed echelon form mod %s primes"%num_primes, level=echelon_verbose_level)
            verbose("current product of primes used: %s"%prod, level=echelon_verbose_level)

            # Use CRT to lift back to ZZ
            mat_over_ZZ = matrix(ZZ, self._base_ring.degree(), self._nrows * self._ncols)
            _lift_crt(mat_over_ZZ, mod_p_ech_ls)

            # Attempt to use rational reconstruction to find
            # our echelon form
            try:
                verbose("attempting rational reconstruction ...", level=echelon_verbose_level)
                res = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(),
                                                 None, None, None)
                res._matrix = <Matrix_rational_dense>matrix_integer_dense_rational_reconstruction(mat_over_ZZ, prod)

            except ValueError:
                # If a ValueError is raised here, it means that the
                # rational reconstruction failed. In this case, add
                # on a few more primes, and try again.

                # TODO: can we reuse the _lift_crt computation?
                num_primes += echelon_primes_increment
                verbose("rational reconstruction failed, trying with %s primes"%num_primes, level=echelon_verbose_level)
                continue

            verbose("rational reconstruction succeeded with %s primes!"%num_primes, level=echelon_verbose_level)

            if ((res * res.denominator()).coefficient_bound() *
                self.coefficient_bound() * self.ncols()) > prod:
                # In this case, we don't know the result to sufficient
                # "precision" (here precision is just the modulus, prod)
                # to guarantee its correctness, so loop.

                # TODO: can we reuse the _lift_crt computation?
                num_primes += echelon_primes_increment
                verbose("height not sufficient to determine echelon form", level=echelon_verbose_level)
                continue

            verbose("found echelon form with %s primes, whose product is %s"%(num_primes, prod), level=echelon_verbose_level)
            self.cache('pivots', max_pivots)
            return res


    def _echelon_form_one_prime(self, p):
        """
        Find the echelon form of self mod the primes dividing p. Return
        the rational matrix representing this lift. If the pivots of the
        reductions mod the primes over p are different, then no such lift
        exists, and we raise a ValueError. If this happens, then the
        denominator of the echelon form of self is divisible by p. (Note
        that the converse need not be true.)

        INPUT:
            p -- a prime that splits completely in the cyclotomic base field.

        OUTPUT:
            matrix -- Lift via CRT of the echelon forms of self modulo
                      each of the primes over p.
            list -- the list of pivots for the echelon form of self mod the
                    primes dividing p

        EXAMPLES:
            sage: W.<z> = CyclotomicField(3)
            sage: A = matrix(W, 2, 3, [1+z, 2/3, 9*z+7, -3 + 4*z, z, -7*z]); A
            [  z + 1     2/3 9*z + 7]
            [4*z - 3       z    -7*z]
            sage: A._echelon_form_one_prime(7)
            ([1 0 4 0 1 2]
            [0 0 3 0 0 4],
            [0, 1])
            sage: Matrix(W,2,3,[2*z+3,0,1,0,1,0])._echelon_form_one_prime(7)
            Traceback (most recent call last):
            ...
            ValueError: echelon form mod 7 not defined

        """
        cdef Matrix_cyclo_dense res
        cdef int i

        # Initialize variables
        is_square = self._nrows == self._ncols
        ls, denom = self._reductions(p)

        # Find our first echelon form, and the associated list
        # of pivots
        ech_ls = [ls[0].echelon_form()]
        pivot_ls = ech_ls[0].pivots()
        # If we've found the identity matrix, we're all done.
        if self._nrows == self._ncols == len(pivot_ls):
            return (self.parent().identity_matrix(), range(self._nrows))

        # For each reduction of self (i.e. for each prime of
        # self.base_ring() over p), compute the echelon form, and
        # keep track of all reductions which have the largest
        # number of pivots seen so far.
        for i from 1 <= i < len(ls):
            ech = ls[i].echelon_form()

            # This should only occur when p divides the denominator
            # of the echelon form of self.
            if ech.pivots() != pivot_ls:
                raise ValueError, "echelon form mod %s not defined"%p

            ech_ls.append(ech)



        # Now, just lift back to ZZ and return it.

        # TODO: coercion going on here
        reduction = matrix(ZZ, len(ech_ls), self._nrows * self._ncols,
                           [ [y.lift() for y in E.list()] for E in ech_ls])

        # TODO: more coercion happening here
        _, Finv = self._reduction_matrix(p)

        lifted_matrix = Finv * reduction

        return (lifted_matrix, pivot_ls)

