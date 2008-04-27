"""
Matrices over Cyclotomic Fields
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
from sage.rings.arith import previous_prime
from matrix cimport Matrix
import matrix_dense
from matrix_integer_dense import _lift_crt
from sage.misc.misc import verbose

from sage.structure.proof.proof import get_flag as get_proof_flag


cdef class Matrix_cyclo_dense(matrix_dense.Matrix_dense):
    ########################################################################
    # LEVEL 1 functionality
    # x * __new__
    # x * __dealloc__     (not needed)
    # x * __init__
    # x * set_unsafe
    # x * get_unsafe
    # x * def _pickle
    # x * def _unpickle
    ########################################################################

    def __new__(self, parent, entries, coerce, copy):
        Matrix.__init__(self, parent)
        self._degree = self._base_ring.degree()

    def __init__(self, parent, entries, copy, coerce):
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
        # The i,j entry is the (i * self._ncols + j)'th column.
        # TODO: This could be made way faster via direct access to the underlying matrix
        cdef Py_ssize_t k, c
        v = value.list()
        c = i * self._ncols + j
        for k from 0 <= k < self._degree:
            self._matrix.set_unsafe(k, c, v[k])

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        # The i,j entry is the (i * self._ncols + j)'th column.
        # TODO: This could be made way faster via direct access to the underlying matrix
        return self.base_ring()(self._matrix.column(i*self._ncols + j).list())

    def _pickle(self):
        data = self._matrix._pickle()
        return data, 0

    def _unpickle(self, data, int version):
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
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)
        # Fix the redundancy here.
        A._matrix = self._matrix + (<Matrix_cyclo_dense>right)._matrix
        return A

    cdef ModuleElement _sub_c_impl(self, ModuleElement right):
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)
        A._matrix = self._matrix - (<Matrix_cyclo_dense>right)._matrix
        return A

    cdef ModuleElement _lmul_c_impl(self, RingElement right):
        if right == 1:
            return self
        elif right == 0:
            return self.parent()(0)

        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)

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

    def __richcmp__(Matrix self, right, int op):
        return self._richcmp(right, op)

    cdef long _hash(self) except -1:
        return self._matrix._hash()

    cdef int _cmp_c_impl(self, Element right) except -2:
        return self._matrix._cmp_c_impl((<Matrix_cyclo_dense>right)._matrix)

    def __copy__(self):
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)
        A._matrix = self._matrix.__copy__()
        return A

    def __neg__(self):
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
        self._matrix.set_immutable()
        matrix_dense.Matrix_dense.set_immutable(self)

    def _rational_matrix(self):
        return self._matrix

    def denominator(self):
        """
        Return the denominator of the intries of this matrix.
        """
        return self._matrix.denominator()

    def randomize(self, density=1, num_bound=2, den_bound=2, distribution=None):
        self._matrix.randomize(density, num_bound, den_bound, distribution)

    def charpoly(self, var='x', algorithm="multimodular"):
        r"""
        Return the characteristic polynomial of self, as a polynomial
        over the base ring.

        INPUT:
            algorithm -- 'multimodular' (default): reduce modulo primes,
                                        compute charpoly mod p, and lift (very fast)
                         'pari': use pari (quite slow; comparable to Magma v2.14 though)
                         'hessenberg': put matrix in Hessenberg form (double dog slow)

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
        """
        f = self.fetch('charpoly')
        if f is not None:
            return f.change_variable_name(var)

        if algorithm == 'multimodular':
            f = self._charpoly_multimodular(var)
        elif algorithm == 'pari':
            f = self._charpoly_over_number_field(var)
        elif algorithm == 'hessenberg':
            f = self._charpoly_hessenberg(var)
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm
        # TODO: Caching commented out only for testing.
        #self.cache('charpoly', f)
        return f

    def _charpoly_mod(self, p):
        """
        Return the characteristic polynomial of self*denom modulo p.

        INPUT:
            p -- a prime that splits completely

        OUTPUT:
            list
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
        Compute the characteristic polynomial of self using a multimodular algorithm.

        INPUT:
            proof -- bool (default: global flag); if False, computes until 3 successive
                     charpoly mod p computations stabilize; always computes modulo at
                     least 6 primes.
        """
        proof = get_proof_flag(proof, "linear_algebra")

        # TODO: only proof = False for now.
        # TODO: Need to figure out how to compute a theoretical bound on the sizes
        # of the coefficients in the charpoly
        proof = False

        n = self._base_ring.zeta_order()
        p = 45989
        prod = 1
        v = []
        L0 = 0
        denom = self.denominator()
        while True:
            while p % n != 1 or denom % p == 0:
                if p == 2:
                    raise RuntimeError, "we ran out of primes in multimodular charpoly algorithm."
                p = previous_prime(p)

            X = self._charpoly_mod(p)
            prod *= p
            v.append(X)
            if proof == False and len(v) % 3 == 0:
                M = matrix(ZZ, self._base_ring.degree(), self._nrows+1)
                L = _lift_crt(M, v)
                if len(v) > 3 and L == L0:
                    break
                else:
                    L0 = L
            p -= 2


        # Now each column of L encodes a coefficient of the output polynomial,
        # with column 0 being the constant coefficient.
        K = self.base_ring()
        coeffs = [K(w.list()) for w in L.columns()]
        R = K[var]
        f = R(coeffs)

        # Rescale to account for denominator, if necessary
        if denom != 1:
            x = R.gen()
            f = f(x * denom) * (1 / (denom**f.degree()))

        return f

    def _reductions(self, p):
        """
        Compute the reductions modulo all primes over p of denom*self, where
        denom is the denominator of self.

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
        # Clear denominator and get matrix over the integers suitable for reduction.
        A, denom = self._matrix._clear_denom()
        # Actually reduce the matrix over the integers modulo the prime p.
        B = A._mod_int(p)
        # Now multiply, which computes from B all the reductions of self*denom
        # modulo each of the primes over p.
        R = T * B
        # Finally compute the actual reductions by extract them from R (note that
        # the rows of R define the reductions).
        ans = R._matrices_from_rows(self._nrows, self._ncols)
        return ans, denom

    def _reduction_matrix(self, p):
        """
        INPUT:
            p -- a prime that splits completely in the base field.

        OUTPUT:
            -- Matrix over GF(p) whose action from the left
               gives the map from O_K to GF(p) x ... x GF(p)
               got by reducing modulo all the primes over p.
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
