# choose: dense

cimport matrix_generic_dense

cdef class Matrix_cyclo_dense(matrix_generic_dense.Matrix_generic_dense):
    def charpoly(self, var='x', algorithm="pari"):
        r"""
        Return the characteristic polynomial of self, as a polynomial
        over the base ring.

        INPUT:
            algorithm -- 'pari'
                         'hessenberg'
        """
        f = self.fetch('charpoly')
        if f is not None:
            return f.change_variable_name(var)

        if algorithm == 'pari':
            f = self._charpoly_over_number_field(var)
        elif algorithm == 'hessenberg':
            f = self._charpoly_hessenberg(var)
        else:
            raise ValueError, "unknown algorithm '%s'"%algorithm
        self.cache('charpoly', f)
        return f

    def _reductions(self, T):
        """
        INPUT:
            T -- n x m matrix over GF(p) where self is n x m.

        OUTPUT:
            list -- of r distinct matrices modulo p, where r is
                    the degree of the cyclotomic base field.
        """
        # Get the base ring GF(p) of the reduction.
        F = T.base_ring()
        # Construct a matrix with 1 row for each entry in self
        X = sum([a.list() for a in self.list()], [])

        from constructor import matrix
        S = matrix(F, self.nrows()*self.ncols(), T.nrows(), X)
        W = S * T
        # The columns of W are the n distinct matrices.
        return [matrix(F, self.nrows(), self.ncols(), w.list()) for w in W.columns()]

    def _charpoly_multimodular(self):
        """
        EXAMPLES:
        One crazy example:
            sage: eps = DirichletGroup(23*3, CyclotomicField(11)).1^2
            sage: M = ModularSymbols(eps); M
            sage: t = M.hecke_matrix(2)
            sage: f = t.charpoly()
            sage: f[0]
            88*zeta11^9 + 252*zeta11^8 + 381*zeta11^7 + 508*zeta11^6 + 594*zeta11^5 + 508*zeta11^4 + 381*zeta11^3 + 252*zeta11^2 + 88*zeta11
            sage: t._charpoly_multimodular()[0]
            88*zeta11^9 + 252*zeta11^8 + 381*zeta11^7 + 508*zeta11^6 + 594*zeta11^5 + 508*zeta11^4 + 381*zeta11^3 + 252*zeta11^2 + 88*zeta11
        """
        K = self.base_ring()
        n = K.zeta_order()

        from sage.rings.arith import previous_prime
        p = 45000
        while p % n != 1:
            p = previous_prime(p)
        # we only do 1 prime so far -- it's just a test.
        return self._charpoly_mod(p)

    def _charpoly_mod(self, p):
        K = self.base_ring()
        n = K.zeta_order()
        if p % n != 1:
            raise ValueError, "n must be 1 modulo p."
        T = reduction_matrices(K, p)
        X = self._reductions(T)
        cp = [A.charpoly() for A in X]
        # Now combine the f together using the inverse of T
        S = T**(-1)

        F = T.base_ring()
        from constructor import matrix
        X = [[f[i] for f in cp] for i in range(self.nrows())]
        W = matrix(F, self.nrows(), T.nrows(), X)
        R = W * S
        return [K(v.list()) for v in R.rows()]


def reduction_matrices(K, p):
    """
    INPUT:
        K -- a cyclotomic field
        p -- a prime that splits completely in K.

    OUTPUT:
        -- Matrix over GF(p) that gives the map from O_K to GF(p) x ... x GF(p)
           got by reducing modulo all the primes over p.
    """
    phi = K.defining_polynomial()
    from sage.rings.all import GF
    from constructor import matrix
    F = GF(p)
    aa = [a for a, _ in phi.change_ring(F).roots()]
    n = K.degree()
    T = matrix(F, n)
    for i in range(n):
        a = aa[i]
        b = 1
        for j in range(n):
            T[j,i] = b
            b *= a
    return T

