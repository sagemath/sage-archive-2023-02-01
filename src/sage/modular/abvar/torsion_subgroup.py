"""
Torsion subgroups of modular abelian varieties.

\sage can compute information about the structure of the torsion
subgroup of a modular abelian variety.  \sage computes a multiple of
the order by computing the greatest common divisor of the orders of
the torsion subgroup of the reduction of the abelian variety modulo p
for various primes p.  \sage computes a divisor of the order by
computing the rational cuspidal subgroup.  When these two bounds agree
(which is often the case), we determine the exact structure of the
torsion subgroup.

AUTHOR:
    -- William Stein (2007-03)

EXAMPLES:
First we consider $J_0(50)$ where everything works out nicely:
    sage: J = J0(50)
    sage: T = J.torsion_subgroup(); T
    Torsion subgroup of Jacobian of the modular curve associated to the congruence subgroup Gamma0(50)
    sage: T.multiple_of_order()
    15
    sage: T.divisor_of_order()
    15
    sage: T.gens()
    [[(1/15, 3/5, -3/5, -1/15)]]
    sage: T.invariants()
    [15]
    sage: d = J.decomposition(); d
    [
    Modular abelian variety quotient of dimension 1 and level 50,
    Modular abelian variety quotient of dimension 1 and level 50
    ]
    sage: d[0].torsion_subgroup().order()
    5
    sage: d[1].torsion_subgroup().order()
    3

Next we make a table of the upper and lower bounds for each new factor.
    sage: for N in range(1,38):
    ...    for A in J0(N).new_quotient().decomposition():
    ...        T = A.torsion_subgroup()
    ...        print '%-5s%-5s%-5s%-5s%-5s'%(N, A.dimension(), A.factor_number(), T.divisor_of_order(), T.multiple_of_order())
    11   1    0    5    5
    14   1    0    6    6
    15   1    0    8    8
    17   1    0    4    4
    19   1    0    3    3
    20   1    0    6    6
    21   1    0    8    8
    23   2    0    11   11
    24   1    0    8    8
    26   1    0    3    3
    26   1    1    7    7
    27   1    0    3    3
    29   2    0    7    7
    30   1    0    6    12
    31   2    0    5    5
    32   1    0    4    4
    33   1    0    4    4
    34   1    0    3    6
    35   1    0    3    3
    35   2    1    16   16
    36   1    0    6    6
    37   1    0    3    3
    37   1    1    1    1

TESTS:
    sage: T = J0(54).torsion_subgroup()
    sage: loads(dumps(T)) == T
    True
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################


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
        """
        Return the order of the torsion subgroup of this modular
        abelian variety.

        This may fail if the multiple obtained by counting points
        modulo $p$ exceeds the divisor obtained from the rational
        cuspidal subgroup.

        EXAMPLES:
            sage: a = J0(11)
            sage: a.torsion_subgroup().order()
            5
            sage: a = J0(23)
            sage: a.torsion_subgroup().order()
            11
            sage: t = J0(37)[0].torsion_subgroup()
            sage: t.order()
            3
        """
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

    def _generators(self):
        A = self.abelian_variety()
        if A.dimension() == 0:
            return []
        R = A.rational_cuspidal_subgroup()
        if R.order() == self.multiple_of_order():
            return R._generators()
        else:
            raise ValueError, "no explicit presentation of this finite subgroup is known (unable to compute explicitly)"

    def possible_orders(self):
        """
        Return the possible orders of this torsion subgroup, computed from
        the divisor and multiple of the order.

        EXAMPLES:
        """
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
        """
        Return a divisor of the order of this torsion subgroup of a
        modular abelian variety.

        EXAMPLES:
           sage: t = J0(37)[0].torsion_subgroup()
           sage: t.divisor_of_order()
           3
        """
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
