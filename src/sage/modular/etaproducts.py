r"""
Eta-products on modular curves X_0(N).

This package provides a class for representing eta-products, which are
meromorphic functions on modular curves of the form
\[\prod_{d | N} \eta(q^d)^{r_d}\]
where $\eta(q)$ is Dirichlet's eta function $q^{1/24} \prod_{n =
1}^\infty(1-q^n)$. These are useful for obtaining explicit models of modular
curves.

See trac ticket #3934 for background.

AUTHOR:
    -- David Loeffler (2008-08-22): initial version
"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                     2008 David Loeffler <d.loeffler.01@cantab.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.arith import divisors, prime_divisors, is_square, euler_phi, gcd
from sage.rings.all import IntegerRing
from sage.structure.formal_sum import FormalSum
from string import join
import sage.misc.latex as latex
from sage.rings.integer_mod import Mod
from sage.rings.integer_mod_ring import IntegerModRing
from sage.matrix.constructor import matrix
from sage.modules.free_module import FreeModule

ZZ = IntegerRing()

class EtaProduct(SageObject):

    def __init__(self, N, rdict):
        r"""
        Create an EtaProduct object representing the function $\prod_{d | N}
        \eta(q^d)^{r_d}$. Checks the critera of Ligozat to ensure that this
        product really is the q-expansion of a meromorphic function on X_0(N).

        INPUT:
            -- (integer) N: a positive integer, the level.
            -- rdict: a dictionary with keys divisors of N and values the
               corresponding r_d. Divisors may be omitted (understood as zero).

        OUTPUT:
            -- an EtaProduct object

        EXAMPLES:
            sage: EtaProduct(3, {3:12, 1:-12})
            Eta product of level 3 : (eta_1)^-12 (eta_3)^12

        NOTE:
            -- The EtaProduct ``knows'' what modular curve it lives on. It is
            possible for two EtaProducts with different $N$'s to be created with
            the same dictionary, and these represent different objects
            (although they will have the same $q$-expansion at the cusp $\infty$).

        """
        # Check Ligozat criteria
        sumR = sumDR = sumNoverDr = 0
        prod = 1

        for d in rdict.keys():
            if N % d:
                raise ValueError, "%s does not divide %s" % (d, N)

        for d in divisors(N):
            if not rdict.has_key(d) or rdict[d] == 0:
                continue
            sumR += rdict[d]
            sumDR += rdict[d]*d
            sumNoverDr += rdict[d]*N/d
            prod *= (N/d)**rdict[d]

        if sumR != 0:
            raise ValueError, "sum r_d is not 0"
        if (sumDR % 24) != 0:
            raise ValueError, "sum d r_d is not 0 mod 24"
        if (sumNoverDr % 24) != 0:
            raise ValueError, "sum (N/d) r_d is not 0 mod 24"
        if not is_square(prod):
            raise ValueError, "product (N/d)^(r_d) is not a square"

        self._N = N
        self._sumDR = sumDR # this is useful to have around
        self._rdict = rdict

    def level(self):
        r""" Return the level of this eta product.
        INPUT:
            -- (none)
        OUTPUT:
            -- (integer) the level of self
        EXAMPLES:

            sage: e = EtaProduct(3, {3:12, 1:-12})
            sage: e.level()
            3
            sage: EtaProduct(12, {6:6, 2:-6}).level() # not the lcm of the d's
            12
            sage: EtaProduct(36, {6:6, 2:-6}).level() # not minimal
            36
        """
        return self._N

    def _repr_(self):
        r"""
        Return the string representation of self.

        EXAMPLES:
            sage: EtaProduct(3, {3:12, 1:-12})
            Eta product of level 3 : (eta_1)^-12 (eta_3)^12
        """
        return "Eta product of level %s : " % self.level() + join(["(eta_%s)^%s" % (d, self.r(d)) for d in divisors(self.level()) if self.r(d) != 0])

    def qexp(self, n):
        r"""
        The q-expansion of self at the cusp at infinity.

        INPUT:
            -- (integer) n: number of terms to calculate

        OUTPUT:
            -- a power series over ZZ in the variable q, with a *relative* precision of 1 + O(q^n).

        ALGORITHM:
            Calculates eta to (n/m) terms, where m is the smallest integer
            dividing self.level() such that self.r(m) != 0. Then multiplies.

        EXAMPLES:
            sage: EtaProduct(36, {6:6, 2:-6}).qexp(10)
            q + 6*q^3 + 27*q^5 + 92*q^7 + 279*q^9 + O(q^11)
            sage: R.<q> = ZZ[[]]
            sage: EtaProduct(2,{2:24,1:-24}).qexp(100) == delta_qexp(101)(q^2)/delta_qexp(101)(q)
            True
        """
        R,q = PowerSeriesRing(ZZ, 'q').objgen()
        pr = R(1)
        eta_n = max([ (n/d).floor() for d in divisors(self.level()) if self.r(d) != 0])
        eta = qexp_eta(R, eta_n)
        for d in divisors(self.level()):
            if self.r(d) != 0:
                pr *= eta(q**d)**self.r(d)
        return pr*q**(self._sumDR / ZZ(24))*( R(1).add_bigoh(n))

    def order_at_cusp(self, cusp):
        r"""
        Return the order of vanishing of self at the given cusp.

        INPUT:
            -- cusp: a CuspFamily object

        OUTPUT:
            -- an integer

        EXAMPLES:
            sage: e = EtaProduct(2, {2:24, 1:-24})
            sage: e.order_at_cusp(CuspFamily(2, 1)) # cusp at infinity
            1
            sage: e.order_at_cusp(CuspFamily(2, 2)) # cusp 0
            -1
        """
        if not isinstance(cusp, CuspFamily):
            raise TypeError, "Argument (=%s) should be a CuspFamily" % cusp
        if cusp.level() != self.level():
            raise ValueError, "Cusp not on right curve!"
        return 1/ZZ(24)/gcd(cusp.width(), self.level()/cusp.width()) * sum( [ell*self.r(ell)/cusp.width() * (gcd(cusp.width(), self.level()/ell))**2  for ell in divisors(self.level())] )

    def divisor(self):
        r"""
        Return the divisor of self, as a formal sum of CuspFamily objects.

        EXAMPLES:
            sage: e = EtaProduct(12, {1:-336, 2:576, 3:696, 4:-216, 6:-576, 12:-144})
            sage: e.divisor() # FormalSum seems to print things in a random order?
            -131*(Inf) - 50*(c_{2}) + 11*(0) + 50*(c_{6}) + 169*(c_{4}) - 49*(c_{3})
            sage: e = EtaProduct(2^8, {8:1,32:-1})
            sage: e.divisor() # random
            -(c_{2}) - (Inf) - (c_{8,2}) - (c_{8,3}) - (c_{8,4}) - (c_{4,2}) - (c_{8,1}) - (c_{4,1}) + (c_{32,4}) + (c_{32,3}) + (c_{64,1}) + (0) + (c_{32,2}) + (c_{64,2}) + (c_{128}) + (c_{32,1})
        """
        return FormalSum([ (self.order_at_cusp(c), c) for c in AllCusps(self.level())])

    def degree(self):
        r"""
        Return the degree of self as a map $X_0(N) \to \mathbb{P}^1$, which is
        equal to the sum of all the positive coefficients in the divisor of
        self.

        EXAMPLES:
             sage: e = EtaProduct(12, {1:-336, 2:576, 3:696, 4:-216, 6:-576, 12:-144})
             sage: e.degree()
             230
        """
        return sum( [self.order_at_cusp(c) for c in AllCusps(self.level()) if self.order_at_cusp(c) > 0])

    def plot(self):
        r""" Returns an error as it's not clear what plotting an eta product means. """
        raise NotImplementedError

    def r(self, d):
        r""" Return the exponent $r_d$ of $\eta(q^d)$ in self.

        EXAMPLES:
            sage: e = EtaProduct(12, {2:24, 3:-24})
            sage: e.r(3)
            -24
            sage: e.r(4)
            0
        """
        if self._rdict.has_key(d):
            return self._rdict[d]
        else:
            return 0

#    def __call__(self, cusp):
#        r""" Calculate the value of self at the given cusp. """
#        assert self.level() == cusp.level()
#        if self.order_at_cusp(cusp) < 0:
#            return Infinity
#        if self.order_at_cusp(cusp) > 0:
#            return 0
#        else:
#            s = ZZ(1)
#            for ell in divisors(self.level()):
#                s *= 1/ZZ(cusp.width())*gcd(cusp.width(), self.level() / ell)**(self.r(ell) / ZZ(2))
#            return s

def basis_eta_products(N, reduce=True):
    r"""
    Produce a basis for the free abelian group of eta-products of level N (under multiplication),
    attempting to find basis vectors of the smallest possible degree.

    INPUT:
        -- an integer $N$ (the level)
        -- a boolean (default True) indicating whether or not to apply LLL-reduction to the calculated basis

    EXAMPLE:
        sage: basis_eta_products(5)
        [Eta product of level 5 : (eta_1)^6 (eta_5)^-6]
        sage: basis_eta_products(12)
        [Eta product of level 12 : (eta_1)^2 (eta_2)^1 (eta_3)^2 (eta_4)^-1 (eta_6)^-7 (eta_12)^3,
        Eta product of level 12 : (eta_1)^3 (eta_2)^-2 (eta_3)^-1 (eta_4)^1 (eta_6)^2 (eta_12)^-3,
        Eta product of level 12 : (eta_1)^-2 (eta_2)^3 (eta_3)^6 (eta_4)^-1 (eta_6)^-9 (eta_12)^3,
        Eta product of level 12 : (eta_1)^-1 (eta_2)^1 (eta_3)^3 (eta_4)^2 (eta_6)^-7 (eta_12)^2,
        Eta product of level 12 : (eta_1)^6 (eta_2)^-9 (eta_3)^-2 (eta_4)^3 (eta_6)^3 (eta_12)^-1]
        sage: basis_eta_products(12, reduce=False) # much bigger coefficients
        [Eta product of level 12 : (eta_1)^15 (eta_2)^-24 (eta_3)^-29 (eta_4)^9 (eta_6)^24 (eta_12)^5,
        Eta product of level 12 : (eta_1)^1 (eta_2)^9 (eta_3)^13 (eta_4)^-4 (eta_6)^-15 (eta_12)^-4,
        Eta product of level 12 : (eta_1)^-8 (eta_2)^-2 (eta_6)^2 (eta_12)^8,
        Eta product of level 12 : (eta_1)^-336 (eta_2)^576 (eta_3)^696 (eta_4)^-216 (eta_6)^-576 (eta_12)^-144,
        Eta product of level 12 : (eta_2)^24 (eta_12)^-24]

    ALGORITHM:
        An eta product of level $N$ is uniquely determined by the integers
        $r_d$ for $d | N$ with $d < N$, since $\sum_{d | N} r_d = 0$. The valid
        $r_d$ are those that satisfy two congruences modulo 24, and one
        congruence modulo 2 for every prime divisor of N. We beef up the
        congruences modulo 2 to congruences modulo 24 by multiplying by 12. To
        calculate the kernel of the ensuing map $\mathbb{Z}^m \to
        (\mathbb{Z}/24\mathbb{Z})^n$ we lift it arbitrarily to an integer
        matrix and calculate its Smith normal form. This gives a basis for the
        lattice.

        This lattice typically contains ``large'' elements, so by default we
        pass it to the eta_lattice_reduce() function to give a more manageable
        basis.
    """

    divs = divisors(N)[:-1]
    s = len(divs)
    primedivs = prime_divisors(N)

    rows = []
    for i in xrange(s):
        # generate a row of relation matrix
        row = [ Mod(divs[i], 24) - Mod(N, 24), Mod(N/divs[i], 24) - Mod(1, 24)]
        for p in primedivs:
            row.append( Mod(12*(N/divs[i]).valuation(p), 24))
        rows.append(row)
    M = matrix(IntegerModRing(24), rows)
    Mlift = M.change_ring(ZZ)
    # now we compute elementary factors of Mlift
    S,U,V = Mlift.smith_form()
    good_vects = []
    for vect in U.rows():
        nf = sum(vect*Mlift*V) # has only one nonzero entry, but hard to predict
                    # which one it is!
        good_vects.append((vect * 24/gcd(nf, 24)).list())
    for v in good_vects:
        v.append(-sum([r for r in v]))
    dicts = []
    for v in good_vects:
        dicts.append({})
        for i in xrange(s):
            dicts[-1][divs[i]] = v[i]
        dicts[-1][N] = v[-1]
    if reduce:
        return eta_lattice_reduce([EtaProduct(N, d) for d in dicts])
    else:
        return [EtaProduct(N, d) for d in dicts]


def eta_lattice_reduce(EtaProducts):
    r"""
    Produce a more manageable basis via LLL-reduction.

    INPUT:
        -- a list of EtaProduct objects (which should all be of the same level)

    OUTPUT:
        -- a new list of EtaProduct objects having hopefully smaller norm

    ALGORITHM:
        We define the norm of an eta-product to be the $L^2$ norm of its
        divisor (as an element of the free $\mathbb{Z}$-module with the cusps
        as basis and the standard inner product). Applying LLL-reduction to
        this gives a basis of hopefully more tractable elements. Of course we'd
        like to use the $L^1$ norm as this is just twice the degree, which is a
        much more natural invariant, but $L^2$ norm is easier to work with!

    EXAMPLES:
        sage: eta_lattice_reduce([ EtaProduct(4, {1:8,2:24,4:-32}), EtaProduct(4, {1:8, 4:-8})])
        [Eta product of level 4 : (eta_1)^8 (eta_4)^-8,
         Eta product of level 4 : (eta_1)^-8 (eta_2)^24 (eta_4)^-16]

        """
    cusps = AllCusps(EtaProducts[0].level())
    r = matrix(ZZ, [[et.order_at_cusp(c) for c in cusps] for et in EtaProducts])
    N = EtaProducts[0].level()
    V = FreeModule(ZZ, r.ncols())
    A = V.submodule_with_basis([V(rw) for rw in r.rows()])
    rred = r.LLL()
    short_etas = []
    for shortvect in rred.rows():
        bv = A.coordinates(shortvect)
        dict = {}
        for d in divisors(N):
            dict[d] = sum( [bv[i]*EtaProducts[i].r(d) for i in xrange(r.nrows())])
        short_etas.append(EtaProduct(N, dict))
    return short_etas

def num_cusps_of_width(N, d):
    r""" Return the number of cusps on $X_0(N)$ of width d.
    INPUT:
        -- (integer) N: the level
        -- (integer) d: an integer dividing N, the cusp width

    EXAMPLES:
        sage: [num_cusps_of_width(18,d) for d in divisors(18)]
        [1, 1, 2, 2, 1, 1]
    """
    assert ((N % d) == 0)
    return euler_phi(gcd(d, N/d))

def AllCusps(N):
    r""" Return a list of CuspFamily objects corresponding to the cusps of X_0(N).

    INPUT:
        -- (integer) N: the level
    EXAMPLES:
        sage: AllCusps(18)
        [(Inf), (c_{2}), (c_{3,1}), (c_{3,2}), (c_{6,1}), (c_{6,2}), (c_{9}), (0)]
    """
    c = []
    for d in divisors(N):
        n = num_cusps_of_width(N, d)
        if n == 1:
            c.append(CuspFamily(N, d))
        elif n > 1:
            for i in xrange(n):
                c.append(CuspFamily(N, d, label=str(i+1)))
    return c

class CuspFamily(SageObject):
    r""" A family of elliptic curves parametrising a region of
    $X_0(N)$."""

    def __init__(self, N, width, label = None):
        r""" Create the cusp of width d on X_0(N) corresponding to the
        family of Tate curves $(\mathbb{C}_p/q^d, \langle \zeta q\rangle)$.
        Here $\zeta$ is a primitive root of unity of order $r$ with
        $\mathrm{lcm}(r,d) = N$. The cusp doesn't store zeta, so we store
        an arbitrary label instead.

        EXAMPLE:
            sage: CuspFamily(8, 4)
            (c_{4})
            sage: CuspFamily(16, 4, '1')
            (c_{4,1})
        """
        self._N = N
        self._width = width
        if (N % width):
            raise ValueError, "Bad width"
        if num_cusps_of_width(N, width) > 1 and label == None:
            raise ValueError, "There are %s > 1 cusps of width %s on X_0(%s): specify a label" % (num_cusps_of_width(N,width), width, N)
        if num_cusps_of_width(N, width) == 1 and label != None:
            raise ValueError, "There is only one cusp of width %s on X_0(%s): no need to specify a label" % (width, N)
        self.label = label

    def width(self):
        r"""
        The width of this cusp.

        EXAMPLES:
            sage: e = CuspFamily(10, 1)
            sage: e.width()
            1
        """
        return self._width

    def level(self):
        r"""
        The width of this cusp.

        EXAMPLES:
            sage: e = CuspFamily(10, 1)
            sage: e.level()
            10
        """
        return self._N

    def sage_cusp(self):
        """
        Return the corresponding element of $\mathbb{P}^1(\mathbb{Q})$.
        EXAMPLE:
            sage: CuspFamily(10, 1).sage_cusp() # not implemented
            Infinity
        """
        raise NotImplementedError

    def _repr_(self):
        r"""
        Return a string representation of self.
        EXAMPLE:
            sage: CuspFamily(16, 4, "1")
            (c_{4,1})
        """
        if self.width() == 1:
            return "(Inf)"
        elif self.width() == self.level():
            return "(0)"
        else:
            return "(c_{%s%s})" % (self.width(), ((self.label and (","+self.label)) or ""))

def qexp_eta(ps_ring, n):
    r"""
    Return the q-expansion of $\eta(q) / q^{1/24}$, where $\eta(q)$ is
    Dedekind's function $$\eta(q) = q^{1/24}\prod_{i=1}^\infty (1-q^i)$$, as an
    element of ps_ring, to precision n. Completely naive algorithm.

    INPUTS:
        -- (PowerSeriesRing) ps_ring: a power series ring -- we pass this as an
            argument as we frequently need to create multiple series in the same ring.
        -- (integer) n: the number of terms to compute.

    OUTPUT:
        An element of ps_ring which is the q-expansion of eta(q)/q^{1/24} truncated to n terms.

    ALGORITHM:
        Multiply out the product $\prod_{i=1}^n (1 - q^i)$. Could perhaps be speeded up by using
        the identity \[ \eta(q) = q^{1/24}( 1 + \sum_{i \ge 1} (-1)^n (q^{n(3n+1)/2} + q^{n(3n-1)/2}),\]
        but I'm lazy.

    EXAMPLES:
        sage: qexp_eta(ZZ[['q']], 100)
        1 - q - q^2 + q^5 + q^7 - q^12 - q^15 + q^22 + q^26 - q^35 - q^40 + q^51 + q^57 - q^70 - q^77 + q^92 + O(q^100)
    """
    q = ps_ring.gen()
    t = ps_ring(1).add_bigoh(n)
    for i in xrange(1,n):
        t = t*ps_ring( 1 - q**i)
    return t
