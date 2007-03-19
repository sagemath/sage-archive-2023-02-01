from finite_subgroup            import FiniteSubgroup
from sage.rings.all             import divisors, gcd, ZZ, prime_range
from sage.sets.primes           import Primes
from sage.modular.congroup      import is_Gamma0



class TorsionSubgroup(FiniteSubgroup):
    def __init__(self, abvar):
        FiniteSubgroup.__init__(self, abvar)

    def _repr_(self):
        return "Torsion subgroup of %s"%self.abelian_variety()

    def order(self):
        try:
            return self._order
        except AttributeError:
            pass
        O = self.possible_orders()
        if len(O) == 1:
            n = O[0]
            self._order = n
            return n
        raise RuntimeError, "Unable to compute order of torsion subgroup (it is in %s)"%O

    def possible_orders(self):
        try:
            return self._possible_orders
        except AttributeError:
            pass
        u = self.multiple_of_order()
        l = self.divisor_of_order()
        assert u % l == 0
        O = [l * d for d in divisors(u//l)]
        self._possible_orders = O
        return O

    def divisor_of_order(self):
        A = self.abelian_variety()
        if A.dimension() == 0:
            return ZZ(1)
        R = A.rational_cuspidal_subgroup()
        return R.order()

    def multiple_of_order(self, maxp=None):
        """
        Return a multiple of the order of this torsion group.

        The multiple is computed using characteristic polynomials of
        Hecke operators of odd index not dividing the level.

        INPUT:
            maxp -- (default: None)  If maxp is None (the default),
                    compute bound until it stabilizes for 3 successive
                    primes.  Otherwise, use all primes up to and
                    including maxp.

        EXAMPLES:
            sage: J = J0(11)
            sage: G = J.torsion_subgroup()
            sage: G.multiple_of_order(11)
            5
            sage: J = J0(389)
            sage: G = J.torsion_subgroup(); G
            Torsion subgroup of Jacobian of the modular curve associated to the congruence subgroup Gamma0(389)
            sage: G.multiple_of_order()
            97
            sage: [G.multiple_of_order(p) for p in prime_range(3,11)]
            [92645296242160800, 7275, 291]
            sage: [G.multiple_of_order(p) for p in prime_range(3,13)]
            [92645296242160800, 7275, 291, 97]
            sage: [G.multiple_of_order(p) for p in prime_range(3,19)]
            [92645296242160800, 7275, 291, 97, 97, 97]
        """
        bnd = ZZ(0)
        A = self.abelian_variety()
        if A.dimension() == 0:
            return ZZ(1)
        N = A.level()
        if not is_Gamma0(A.group()):
            # to generalize to this case, you'll need to
            # (1) define a charpoly_of_frob function:
            #       this is tricky because I don't know a simple
            #       way to do this for Gamma1 and GammaH.  Will
            #       probably have to compute explicit matrix for
            #       <p> operator (add to modular symbols code),
            #       then compute some charpoly involving
            #       that directly...
            # (2) use (1) -- see my MAGMA code.
            raise NotImplementedError, "torsion multiple only implemented for Gamma0"
        cnt = 0
        if maxp is None:
            X = Primes()
        else:
            X = prime_range(maxp+1)
        for p in X:
            if (2*N) % p == 0:
                continue

            f = A.hecke_polynomial(p)
            b = f(p+1)

            if bnd == 0:
                bnd = b
            else:
                bnd_last = bnd
                bnd = ZZ(gcd(bnd, b))
                if bnd == bnd_last:
                    cnt += 1
                else:
                    cnt = 0
                if maxp is None and cnt >= 2:
                    break

        return bnd
