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
from sage.rings.arith import previous_prime, binomial
from matrix cimport Matrix
import matrix_dense
from matrix_integer_dense import _lift_crt
from sage.misc.misc import verbose
import math

from misc import matrix_integer_dense_rational_reconstruction

from sage.ext.multi_modular import MAX_MODULUS

from sage.structure.proof.proof import get_flag as get_proof_flag

#TODO: only here for debugging
import time

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
        """
        return self._matrix.denominator()

    def coefficient_bound(self):
        r"""
        Return an upper bound for the (complex) absolute values of all
        entries of self.

        This is computed using just the Cauchy-Schwarz inequality, i.e.,
        we use the fact that
        $$ \left| \sum_i a_i\zeta^i \right| \leq \sum_i |a_i|, $$
        as $|\zeta| = 1$.
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
        Return the height of self. If we let $a_{ij}$ be the $i$,$j$
        entry of self, then this is defined to be
        $$ \operatorname{max}_v\ \operatorname{max}_{i, j} |a_{ij}|, $$
        where $v$ runs over all complex embeddings of
        \code{self.base_ring()}.
        """
        cdef Py_ssize_t i, j

        emb = self._base_ring.complex_embeddings()

        ht = 0
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                t = max([ x.norm().sqrt() for x in [ f(self[i,j]) for f in emb ] ])
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
            sage: A.randomize() ; A
            [       1/2*zeta5^2 + zeta5                        1/2]
            [        -zeta5^2 + 2*zeta5 -2*zeta5^3 + 2*zeta5^2 + 2]
        """
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
            147
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
        # This is 2*e^(1-(2(7\gamma-4))/(13(3-2\gamma))), where \gamma is
        # Euler's constant.
        delta = 5.418236

        B = self.coefficient_bound()

        # this bound is only valid for n >= 4, use naive bounds
        # in other cases
        if self._nrows > 3:
            return ZZ(int(math.ceil((alpha * self._nrows * B**2)**(self._nrows/2.0))))
        elif self._nrows == 3:
            return max(6*B**2, 4*B**3)
        elif self._nrows == 2:
            return 2*B**2
        else:
            return B

        # This is the code for computing the bound from Lemma 2.2
        # of the paper. I don't get it: this bound seems to always
        # be *much* worse than the one from Lemma 2.1 ... why would
        # you want this bound? Or did I screw something up? I checked
        # a few examples by hand, and seemed to be correctly computing
        # the bound in the paper ...

#         D = ZZ(int(math.ceil((math.sqrt(1+2*delta*self._nrows*(B**2))-1)/(delta*B**2))))

#         # TODO: we don't check anything about overflows anywhere here;
#         # should we?

#         # i = 0 case
#         i = 0
#         M = ZZ(int(math.ceil((self._nrows * B**2)**(self._nrows/2.0))))
#         for i from 1 <= i < D:
#             val = ZZ(int(math.ceil(binomial(self._nrows, i) *
#                                    ((self._nrows-i)*B**2)**((n-i)/2.0))))
#             if val > M:
#                 M = val

#         other_bound = ZZ(int(math.ceil((alpha * self._nrows * B)**(self._nrows/2.0))))
#         return min(M, other_bound)



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

    # TODO: delete this, it's only here for debugging
    def xxx_charpoly_multimodular(self, var='x', proof=None):
        """
        Compute the characteristic polynomial of self using a multimodular algorithm.

        INPUT:
            proof -- bool (default: global flag); if False, computes
                     until 3 successive charpoly mod p computations
                     stabilize; always computes modulo at least 6
                     primes.
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

    def _charpoly_multimodular(self, var='x', proof=None):
        """
        Compute the characteristic polynomial of self using a multimodular algorithm.

        INPUT:
            proof -- bool (default: global flag); if False, computes
                     until 3 successive charpoly mod p computations
                     stabilize; always computes modulo at least 6
                     primes.
        """
        proof = get_proof_flag(proof, "linear_algebra")

        # TODO: only proof = False for now.
        # TODO: Need to figure out how to compute a theoretical bound on the sizes
        # of the coefficients in the charpoly
        proof = False

        # TODO: do we want self.zeta_order() or self.n()?
        # (this only matters for n odd.)
        n = self._base_ring._n()
        p = 45989
        #p = previous_prime(MAX_MODULUS)
        prod = 1
        v = []
        denom = self.denominator()
        A = self*denom
        bound = A._charpoly_bound()
        #bound = self._charpoly_bound()
        while prod < bound:
            while p % n != 1 or denom % p == 0:
                if p == 2:
                    raise RuntimeError, "we ran out of primes in multimodular charpoly algorithm."
                p = previous_prime(p)

            X = A._charpoly_mod(p)
            #X = self._charpoly_mod(p)
            v.append(X)
            prod *= p
            p = previous_prime(p)

        M = matrix(ZZ, self._base_ring.degree(), self._nrows+1)
        L = _lift_crt(M, v)

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
        Find the echelon form of self, using the specified
        algorithm.
        """
        E = self.fetch('echelon_form')
        if E is not None:
            return E

        if algorithm == 'multimodular':
            E = self._echelon_form_multimodular()
        elif algorithm == 'classical':
            E = (self*self.denominator())._echelon_classical()
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm

        # TODO: Caching commented out only for testing.
        #self.cache('charpoly', f)
        return E

    def _echelon_form_multimodular(self, num_primes=5, max_prime=None, height_guess=None):
        """
        Use a multimodular algorithm to find the echelon
        form of self.
        """

        cdef Matrix_cyclo_dense res

        denom = self.denominator()
        A = denom * self

        # TODO: max_prime is only here for debugging purposes.
        if max_prime is None:
            max_prime = MAX_MODULUS

        # TODO: i basically copied this from matrix_rational_dense
        # ... we should think about it, at the very least.
        if height_guess is None:
            height_guess = (A.coefficient_bound()+100)*1000000

        # generate primes to use
        # TODO: I stupidly just generate 5 primes and use those.
        # this obviously needs to change.
        p = previous_prime(max_prime)
        found = 0
        prime_ls = []
        prod = 1
        n = self._base_ring._n()
        height_bound = self._ncols * height_guess * A.coefficient_bound() + 1
        while found < num_primes or prod < height_bound:
            if p%n == 1:
                prime_ls.append(p)
                found += 1
                prod *= p
            p = previous_prime(p)

        mod_p_ech_ls = []
        num_pivots = 0
        for pp in prime_ls:
            mod_p_ech, piv = A._echelon_form_one_prime(pp)
            # if we have the identity, just return it, and
            # we're done.
            if mod_p_ech.is_one():
                return mod_p_ech
            if piv > num_pivots:
                mod_p_ech_ls = [mod_p_ech]
                num_pivots = piv
            elif piv == num_pivots:
                mod_p_ech_ls.append(mod_p_ech)

        mat_over_ZZ = matrix(ZZ, self._base_ring.degree(), self._nrows * self._ncols)
        _lift_crt(mat_over_ZZ, mod_p_ech_ls)

        try:
            res = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(),
                                             None, None, None)
            res._matrix = <Matrix_rational_dense>matrix_integer_dense_rational_reconstruction(mat_over_ZZ, prod)

        except ValueError:
            if num_primes > 100:
                raise ValueError, "sorry buddy, i'm just tired."
            # TODO: this is just throwing away the info constructed so far
            return self._echelon_form_multimodular(found + 5, max_prime, height_guess)

        if ((res * res.denominator()).coefficient_bound() *
            self.coefficient_bound()) > prod:
            # TODO: this is just throwing away the info constructed so far
            return self._echelon_form_multimodular(found + 5, max_prime, height_guess)

        return res

    def _echelon_form_one_prime(self, p):

        cdef Matrix_cyclo_dense res

        is_square = self._nrows == self._ncols

        ls, denom = self._reductions(p)

        ech_ls = []
        most_pivots = 0
        longest_pivot_ls = []
        longest_pivot_index_ls = []

        for i in range(len(ls)):
            ech = ls[i].echelon_form()
            piv = ech.pivots()
            if is_square and len(piv) == self._nrows:
                return (self.parent().identity_matrix(), self._nrows)
            ech_ls.append(ech)
            if len(piv) > most_pivots:
                most_pivots = len(piv)

        #return ech_ls, longest_pivot_ls, longest_pivot_index_ls

        # TODO: coercion going on here
        reduction = matrix(ZZ, len(ech_ls), self._nrows * self._ncols,
                           [ [y.lift() for y in E.list()] for E in ech_ls])

        # TODO: more coercion happening here
        _, Finv = self._reduction_matrix(p)

        lifted_matrix = Finv * reduction

        return (lifted_matrix, most_pivots)

