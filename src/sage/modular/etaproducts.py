r"""
Eta-products on modular curves :math:`X_0(N)`

This package provides a class for representing eta-products, which
are meromorphic functions on modular curves of the form

.. math::

    \prod_{d | N} \eta(q^d)^{r_d}

where
:math:`\eta(q)` is Dirichlet's eta function
:math:`q^{1/24} \prod_{n = 1}^\infty(1-q^n)`. These are useful
for obtaining explicit models of modular curves.

See trac ticket #3934 for background.

AUTHOR:

- David Loeffler (2008-08-22): initial version
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
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.arith import divisors, prime_divisors, is_square, euler_phi, gcd
from sage.rings.all import Integer, IntegerRing, RationalField
from sage.groups.old import AbelianGroup
from sage.structure.element import MultiplicativeGroupElement
from sage.structure.formal_sum import FormalSum
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.matrix.constructor import matrix
from sage.modules.free_module import FreeModule
from sage.misc.misc import union

from string import join
import weakref

ZZ = IntegerRing()
QQ = RationalField()

_cache = {}
def EtaGroup(level):
    r"""
    Create the group of eta products of the given level.

    EXAMPLES::

        sage: EtaGroup(12)
        Group of eta products on X_0(12)
        sage: EtaGroup(1/2)
        Traceback (most recent call last):
        ...
        TypeError: Level (=1/2) must be a positive integer
        sage: EtaGroup(0)
        Traceback (most recent call last):
        ...
        ValueError: Level (=0) must be a positive integer
    """
    if level in _cache:
        G = _cache[level]()
        if not G is None:
            return G
    G = EtaGroup_class(level)
    _cache[level] = weakref.ref(G)
    return G

class EtaGroup_class(AbelianGroup):
    r"""
    The group of eta products of a given level under multiplication.
    """

    def __init__(self, level):
        r"""
        Create the group of eta products of a given level, which must be a
        positive integer.

        EXAMPLES:
            sage: G = EtaGroup(12); G # indirect doctest
            Group of eta products on X_0(12)
            sage: G is loads(dumps(G))
            True
        """
        try:
            level = ZZ(level)
        except TypeError:
            raise TypeError("Level (=%s) must be a positive integer" % level)
        if (level < 1):
            raise ValueError("Level (=%s) must be a positive integer" % level)
        self._N = level

    def __reduce__(self):
        r"""
        Return the data used to construct self. Used for pickling.

        EXAMPLE::

            sage: EtaGroup(13).__reduce__()
            (<function EtaGroup at ...>, (13,))
        """
        return (EtaGroup, (self.level(),))

    def __cmp__(self, other):
        r"""
        Compare self to other. If other is not an EtaGroup, compare by
        type; otherwise compare by level. EtaGroups of the same level
        compare as identical.

        EXAMPLE::

            sage: EtaGroup(12) == 12
            False
            sage: EtaGroup(12) < EtaGroup(13)
            True
            sage: EtaGroup(12) == EtaGroup(12)
            True
        """
        if not isinstance(other, EtaGroup_class):
            return cmp(type(self), type(other))
        else:
            return cmp(self.level(), other.level())

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: EtaGroup(12)._repr_()
            'Group of eta products on X_0(12)'
        """
        return "Group of eta products on X_0(%s)" % self.level()

    def one(self):
        r"""
        Return the identity element of ``self``.

        EXAMPLE::

            sage: EtaGroup(12).one()
            Eta product of level 12 : 1
        """
        return self({})

    def __call__(self, dict):
        r"""
        Create an element of this group (an eta product object) with
        exponents from the given dictionary. See the docstring for the
        EtaProduct() factory function for how dict is used.

        EXAMPLE::

            sage: EtaGroup(2).__call__({1:24, 2:-24})
            Eta product of level 2 : (eta_1)^24 (eta_2)^-24
        """
        return EtaGroupElement(self, dict)

    def level(self):
        r""" Return the level of self.
        EXAMPLES::

            sage: EtaGroup(10).level()
            10
        """
        return self._N

    def basis(self, reduce=True):
        r"""
        Produce a basis for the free abelian group of eta-products of level
        N (under multiplication), attempting to find basis vectors of the
        smallest possible degree.

        INPUT:


        -  ``reduce`` - a boolean (default True) indicating
           whether or not to apply LLL-reduction to the calculated basis


        EXAMPLE::

            sage: EtaGroup(5).basis()
            [Eta product of level 5 : (eta_1)^6 (eta_5)^-6]
            sage: EtaGroup(12).basis()
            [Eta product of level 12 : (eta_1)^2 (eta_2)^1 (eta_3)^2 (eta_4)^-1 (eta_6)^-7 (eta_12)^3,
            Eta product of level 12 : (eta_1)^-4 (eta_2)^2 (eta_3)^4 (eta_6)^-2,
            Eta product of level 12 : (eta_1)^-1 (eta_2)^3 (eta_3)^3 (eta_4)^-2 (eta_6)^-9 (eta_12)^6,
            Eta product of level 12 : (eta_1)^1 (eta_2)^-1 (eta_3)^-3 (eta_4)^-2 (eta_6)^7 (eta_12)^-2,
            Eta product of level 12 : (eta_1)^-6 (eta_2)^9 (eta_3)^2 (eta_4)^-3 (eta_6)^-3 (eta_12)^1]
            sage: EtaGroup(12).basis(reduce=False) # much bigger coefficients
            [Eta product of level 12 : (eta_2)^24 (eta_12)^-24,
            Eta product of level 12 : (eta_1)^-336 (eta_2)^576 (eta_3)^696 (eta_4)^-216 (eta_6)^-576 (eta_12)^-144,
            Eta product of level 12 : (eta_1)^-8 (eta_2)^-2 (eta_6)^2 (eta_12)^8,
            Eta product of level 12 : (eta_1)^1 (eta_2)^9 (eta_3)^13 (eta_4)^-4 (eta_6)^-15 (eta_12)^-4,
            Eta product of level 12 : (eta_1)^15 (eta_2)^-24 (eta_3)^-29 (eta_4)^9 (eta_6)^24 (eta_12)^5]

        ALGORITHM: An eta product of level `N` is uniquely
        determined by the integers `r_d` for `d | N` with
        `d < N`, since `\sum_{d | N} r_d = 0`. The valid
        `r_d` are those that satisfy two congruences modulo 24,
        and one congruence modulo 2 for every prime divisor of N. We beef
        up the congruences modulo 2 to congruences modulo 24 by multiplying
        by 12. To calculate the kernel of the ensuing map
        `\ZZ^m \to (\ZZ/24\ZZ)^n`
        we lift it arbitrarily to an integer matrix and calculate its Smith
        normal form. This gives a basis for the lattice.

        This lattice typically contains "large" elements, so by default we
        pass it to the reduce_basis() function which performs
        LLL-reduction to give a more manageable basis.
        """

        N = self.level()
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
        for i in xrange(U.nrows()):
            vect = U.row(i)
            nf = (i < S.ncols() and S[i,i]) or 0
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
            return self.reduce_basis([ self(d) for d in dicts])
        else:
            return [self(d) for d in dicts]

    def reduce_basis(self, long_etas):
        r"""
        Produce a more manageable basis via LLL-reduction.

        INPUT:


        - ``long_etas`` -  a list of EtaGroupElement objects (which
          should all be of the same level)


        OUTPUT:


        - a new list of EtaGroupElement objects having
          hopefully smaller norm


        ALGORITHM: We define the norm of an eta-product to be the
        `L^2` norm of its divisor (as an element of the free
        `\ZZ`-module with the cusps as basis and the
        standard inner product). Applying LLL-reduction to this gives a
        basis of hopefully more tractable elements. Of course we'd like to
        use the `L^1` norm as this is just twice the degree, which
        is a much more natural invariant, but `L^2` norm is easier
        to work with!

        EXAMPLES::

            sage: EtaGroup(4).reduce_basis([ EtaProduct(4, {1:8,2:24,4:-32}), EtaProduct(4, {1:8, 4:-8})])
            [Eta product of level 4 : (eta_1)^8 (eta_4)^-8,
            Eta product of level 4 : (eta_1)^-8 (eta_2)^24 (eta_4)^-16]
        """
        N = self.level()
        cusps = AllCusps(N)
        r = matrix(ZZ, [[et.order_at_cusp(c) for c in cusps] for et in long_etas])
        V = FreeModule(ZZ, r.ncols())
        A = V.submodule_with_basis([V(rw) for rw in r.rows()])
        rred = r.LLL()
        short_etas = []
        for shortvect in rred.rows():
            bv = A.coordinates(shortvect)
            dict = {}
            for d in divisors(N):
                dict[d] = sum( [bv[i]*long_etas[i].r(d) for i in xrange(r.nrows())])
            short_etas.append(self(dict))
        return short_etas


def EtaProduct(level, dict):
    r"""
    Create an EtaGroupElement object representing the function
    `\prod_{d | N} \eta(q^d)^{r_d}`. Checks the criteria
    of Ligozat to ensure that this product really is the q-expansion of
    a meromorphic function on X_0(N).

    INPUT:


    -  ``level`` -  (integer): the N such that this eta
       product is a function on X_0(N).

    -  ``dict`` - (dictionary): a dictionary indexed by
       divisors of N such that the coefficient of `\eta(q^d)` is
       r[d]. Only nonzero coefficients need be specified. If Ligozat's
       criteria are not satisfied, a ValueError will be raised.


    OUTPUT:


    -  an EtaGroupElement object, whose parent is
       the EtaGroup of level N and whose coefficients are the given
       dictionary.


    .. note::

       The dictionary dict does not uniquely specify N. It is
       possible for two EtaGroupElements with different `N`'s to
       be created with the same dictionary, and these represent different
       objects (although they will have the same `q`-expansion at
       the cusp `\infty`).

    EXAMPLES::

        sage: EtaProduct(3, {3:12, 1:-12})
        Eta product of level 3 : (eta_1)^-12 (eta_3)^12
        sage: EtaProduct(3, {3:6, 1:-6})
        Traceback (most recent call last):
        ...
        ValueError: sum d r_d (=12) is not 0 mod 24
        sage: EtaProduct(3, {4:6, 1:-6})
        Traceback (most recent call last):
        ...
        ValueError: 4 does not divide 3
    """
    return EtaGroup(level)(dict)

class EtaGroupElement(MultiplicativeGroupElement):

    def __init__(self, parent, rdict):
        r"""
        Create an eta product object. Usually called implicitly via
        EtaGroup_class.__call__ or the EtaProduct factory function.

        EXAMPLE::

            sage: EtaGroupElement(EtaGroup(8), {1:24, 2:-24})
            Eta product of level 8 : (eta_1)^24 (eta_2)^-24
            sage: g = _; g == loads(dumps(g))
            True
        """
        MultiplicativeGroupElement.__init__(self, parent)

        self._N = self.parent().level()
        N = self._N

        if isinstance(rdict, EtaGroupElement):
            rdict = rdict._rdict
            # Note: This is needed because the "x in G" test tries to call G(x)
            # and see if it returns an error. So sometimes this will be getting
            # called with rdict being an eta product, not a dictionary.

        if rdict == 1:
            rdict = {}
        # Check Ligozat criteria
        sumR = sumDR = sumNoverDr = 0
        prod = 1

        for d in rdict.keys():
            if N % d:
                raise ValueError("%s does not divide %s" % (d, N))

        for d in rdict.keys():
            if rdict[d] == 0:
                rdict.pop(d)
                continue
            sumR += rdict[d]
            sumDR += rdict[d]*d
            sumNoverDr += rdict[d]*N/d
            prod *= (N/d)**rdict[d]

        if sumR != 0:
            raise ValueError("sum r_d (=%s) is not 0" % sumR)
        if (sumDR % 24) != 0:
            raise ValueError("sum d r_d (=%s) is not 0 mod 24" % sumDR)
        if (sumNoverDr % 24) != 0:
            raise ValueError("sum (N/d) r_d (=%s) is not 0 mod 24" % sumNoverDr)
        if not is_square(prod):
            raise ValueError("product (N/d)^(r_d) (=%s) is not a square" % prod)

        self._sumDR = sumDR # this is useful to have around
        self._rdict = rdict
        self._keys = rdict.keys() # avoid factoring N every time

    def _mul_(self, other):
        r"""
        Return the product of self and other.

        EXAMPLES::

            sage: eta1, eta2 = EtaGroup(4).basis() # indirect doctest
            sage: eta1 * eta2
            Eta product of level 4 : (eta_1)^8 (eta_4)^-8
        """
        newdict = {}
        for d in union(self._keys, other._keys):
            newdict[d] = self.r(d) + other.r(d)
        return EtaProduct(self.level(), newdict)

    def _div_(self, other):
        r"""
        Return `self * other^{-1}`.

        EXAMPLES::

            sage: eta1, eta2 = EtaGroup(4).basis()
            sage: eta1 / eta2 # indirect doctest
            Eta product of level 4 : (eta_1)^-24 (eta_2)^48 (eta_4)^-24
            sage: (eta1 / eta2) * eta2 == eta1
            True
        """
        newdict = {}
        for d in union(self._keys, other._keys):
            newdict[d] = self.r(d) - other.r(d)
        return EtaProduct(self.level(), newdict)

    def __cmp__(self, other):
        r"""
        Compare self to other. Eta products compare first according to
        their levels, then according to their rdicts.

        EXAMPLES::

            sage: EtaProduct(2, {2:24,1:-24}) == 1
            False
            sage: EtaProduct(2, {2:24, 1:-24}) < EtaProduct(4, {2:24, 1:-24})
            True
            sage: EtaProduct(2, {2:24, 1:-24}) == EtaProduct(4, {2:24, 1:-24})
            False
            sage: EtaProduct(2, {2:24, 1:-24}) < EtaProduct(4, {2:48, 1:-48})
            True
        """
        if not isinstance(other, EtaGroupElement):
            return cmp(type(self), type(other))
        return (cmp(self.level(), other.level()) or cmp(self._rdict, other._rdict))

    def _short_repr(self):
        r"""
        A short string representation of self, which doesn't specify the
        level.

        EXAMPLES::

            sage: EtaProduct(3, {3:12, 1:-12})._short_repr()
            '(eta_1)^-12 (eta_3)^12'
        """
        if self.degree() == 0:
            return "1"
        else:
            return join(["(eta_%s)^%s" % (d,self.r(d)) for d in self._keys])

    def _repr_(self):
        r"""
        Return the string representation of self.

        EXAMPLES::

            sage: EtaProduct(3, {3:12, 1:-12})._repr_()
            'Eta product of level 3 : (eta_1)^-12 (eta_3)^12'
        """
        return "Eta product of level %s : " % self.level() + self._short_repr()

    def level(self):
        r"""
        Return the level of this eta product.

        EXAMPLES::

            sage: e = EtaProduct(3, {3:12, 1:-12})
            sage: e.level()
            3
            sage: EtaProduct(12, {6:6, 2:-6}).level() # not the lcm of the d's
            12
            sage: EtaProduct(36, {6:6, 2:-6}).level() # not minimal
            36
        """
        return self._N

    def q_expansion(self, n):
        r"""
        The q-expansion of self at the cusp at infinity.

        INPUT:


        - ``n`` (integer): number of terms to calculate


        OUTPUT:


        -  a power series over `\ZZ` in
           the variable `q`, with a *relative* precision of
           `1 + O(q^n)`.


        ALGORITHM: Calculates eta to (n/m) terms, where m is the smallest
        integer dividing self.level() such that self.r(m) != 0. Then
        multiplies.

        EXAMPLES::

            sage: EtaProduct(36, {6:6, 2:-6}).q_expansion(10)
            q + 6*q^3 + 27*q^5 + 92*q^7 + 279*q^9 + O(q^11)
            sage: R.<q> = ZZ[[]]
            sage: EtaProduct(2,{2:24,1:-24}).q_expansion(100) == delta_qexp(101)(q^2)/delta_qexp(101)(q)
            True
        """
        R,q = PowerSeriesRing(ZZ, 'q').objgen()
        pr = R(1)
        if self == self.parent()(1):
            return pr
        eta_n = max([ (n/d).floor() for d in self._keys if self.r(d) != 0])
        eta = qexp_eta(R, eta_n)
        for d in self._keys:
            if self.r(d) != 0:
                pr *= eta(q**d)**self.r(d)
        return pr*q**(self._sumDR / ZZ(24))*( R(1).add_bigoh(n))

    def qexp(self, n):
        """
        Alias for ``self.q_expansion()``.

        EXAMPLES::

            sage: e = EtaProduct(36, {6:8, 3:-8})
            sage: e.qexp(10)
            q + 8*q^4 + 36*q^7 + O(q^10)
            sage: e.qexp(30) == e.q_expansion(30)
            True
        """
        return self.q_expansion(n)

    def order_at_cusp(self, cusp):
        r"""
        Return the order of vanishing of self at the given cusp.

        INPUT:


        -  ``cusp`` -  a CuspFamily object


        OUTPUT:


        - an integer


        EXAMPLES::

            sage: e = EtaProduct(2, {2:24, 1:-24})
            sage: e.order_at_cusp(CuspFamily(2, 1)) # cusp at infinity
            1
            sage: e.order_at_cusp(CuspFamily(2, 2)) # cusp 0
            -1
        """
        if not isinstance(cusp, CuspFamily):
            raise TypeError("Argument (=%s) should be a CuspFamily" % cusp)
        if cusp.level() != self.level():
            raise ValueError("Cusp not on right curve!")
        return 1/ZZ(24)/gcd(cusp.width(), self.level()//cusp.width()) * sum( [ell*self.r(ell)/cusp.width() * (gcd(cusp.width(), self.level()//ell))**2  for ell in self._keys] )

    def divisor(self):
        r"""
        Return the divisor of self, as a formal sum of CuspFamily objects.

        EXAMPLES::

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
        Return the degree of self as a map
        `X_0(N) \to \mathbb{P}^1`, which is equal to the sum of
        all the positive coefficients in the divisor of self.

        EXAMPLES::

            sage: e = EtaProduct(12, {1:-336, 2:576, 3:696, 4:-216, 6:-576, 12:-144})
            sage: e.degree()
            230
        """
        return sum( [self.order_at_cusp(c) for c in AllCusps(self.level()) if self.order_at_cusp(c) > 0])

#     def plot(self):
#         r""" Returns an error as it's not clear what plotting an eta product means. """
#         raise NotImplementedError

    def r(self, d):
        r"""
        Return the exponent `r_d` of `\eta(q^d)` in self.

        EXAMPLES::

            sage: e = EtaProduct(12, {2:24, 3:-24})
            sage: e.r(3)
            -24
            sage: e.r(4)
            0
        """
        return self._rdict.get(d, 0)

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
#                s *= 1/ZZ(cusp.width())*gcd(cusp.width(), self.level() // ell)**(self.r(ell) / ZZ(2))
#            return s

def num_cusps_of_width(N, d):
    r"""
    Return the number of cusps on `X_0(N)` of width d.

    INPUT:


    -  ``N`` - (integer): the level

    -  ``d`` - (integer): an integer dividing N, the cusp
       width


    EXAMPLES::

        sage: [num_cusps_of_width(18,d) for d in divisors(18)]
        [1, 1, 2, 2, 1, 1]
    """
    try:
        N = ZZ(N)
        d = ZZ(d)
        assert N>0
        assert d>0
        assert ((N % d) == 0)
    except TypeError:
        raise TypeError("N and d must be integers")
    except AssertionError:
        raise AssertionError("N and d must be positive integers with d|N")

    return euler_phi(gcd(d, N//d))

def AllCusps(N):
    r"""
    Return a list of CuspFamily objects corresponding to the cusps of
    `X_0(N)`.

    INPUT:

    -  ``N`` - (integer): the level


    EXAMPLES::

        sage: AllCusps(18)
        [(Inf), (c_{2}), (c_{3,1}), (c_{3,2}), (c_{6,1}), (c_{6,2}), (c_{9}), (0)]
    """
    try:
        N = ZZ(N)
        assert N>0
    except TypeError:
        raise TypeError("N must be an integer")
    except AssertionError:
        raise AssertionError("N must be positive")
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
    r"""
    A family of elliptic curves parametrising a region of
    `X_0(N)`.
    """
    def __init__(self, N, width, label = None):
        r"""
        Create the cusp of width d on X_0(N) corresponding to the family
        of Tate curves
        `(\CC_p/q^d, \langle \zeta q\rangle)`. Here
        `\zeta` is a primitive root of unity of order `r`
        with `\mathrm{lcm}(r,d) = N`. The cusp doesn't store zeta,
        so we store an arbitrary label instead.

        EXAMPLE::

            sage: CuspFamily(8, 4)
            (c_{4})
            sage: CuspFamily(16, 4, '1')
            (c_{4,1})
        """
        try:
            N = ZZ(N)
            assert N>0
        except TypeError:
            raise TypeError("N must be an integer")
        except AssertionError:
            raise AssertionError("N must be positive")
        self._N = N
        self._width = width
        if (N % width):
            raise ValueError("Bad width")
        if num_cusps_of_width(N, width) > 1 and label is None:
            raise ValueError("There are %s > 1 cusps of width %s on X_0(%s): specify a label" % (num_cusps_of_width(N,width), width, N))
        if num_cusps_of_width(N, width) == 1 and label is not None:
            raise ValueError("There is only one cusp of width %s on X_0(%s): no need to specify a label" % (width, N))
        self.label = label

    def width(self):
        r"""
        The width of this cusp.

        EXAMPLES::

            sage: e = CuspFamily(10, 1)
            sage: e.width()
            1
        """
        return self._width

    def level(self):
        r"""
        The level of this cusp.

        EXAMPLES::

            sage: e = CuspFamily(10, 1)
            sage: e.level()
            10
        """
        return self._N

    def sage_cusp(self):
        """
        Return the corresponding element of
        `\mathbb{P}^1(\QQ)`.

        EXAMPLE::

            sage: CuspFamily(10, 1).sage_cusp() # not implemented
            Infinity
        """
        raise NotImplementedError

    def _repr_(self):
        r"""
        Return a string representation of self.

        EXAMPLE::

            sage: CuspFamily(16, 4, "1")._repr_()
            '(c_{4,1})'
        """
        if self.width() == 1:
            return "(Inf)"
        elif self.width() == self.level():
            return "(0)"
        else:
            return "(c_{%s%s})" % (self.width(), ((self.label and (","+self.label)) or ""))

def qexp_eta(ps_ring, prec):
    r"""
    Return the q-expansion of `\eta(q) / q^{1/24}`, where
    `\eta(q)` is Dedekind's function

    .. math::

        \eta(q) = q^{1/24}\prod_{n=1}^\infty (1-q^n),


    as an element of ps_ring, to precision prec.

    INPUT:

    -  ``ps_ring`` - (PowerSeriesRing): a power series ring

    -  ``prec`` - (integer): the number of terms to compute.


    OUTPUT: An element of ps_ring which is the q-expansion of
    `\eta(q)/q^{1/24}` truncated to prec terms.

    ALGORITHM: We use the Euler identity

    .. math::

         \eta(q) = q^{1/24}( 1 + \sum_{n \ge 1} (-1)^n (q^{n(3n+1)/2} + q^{n(3n-1)/2})

    to compute the expansion.

    EXAMPLES::

        sage: qexp_eta(ZZ[['q']], 100)
        1 - q - q^2 + q^5 + q^7 - q^12 - q^15 + q^22 + q^26 - q^35 - q^40 + q^51 + q^57 - q^70 - q^77 + q^92 + O(q^100)
    """
    prec = Integer(prec)
    assert prec>0, "prec must be a positive integer"
    v = [Integer(0)] * prec
    pm = Integer(1)
    v[0] = pm
    try:
        n = 1
        while True:
            pm = -pm
            v[n*(3*n-1)//2] = pm
            v[n*(3*n+1)//2] = pm
            n += 1
    except IndexError:
        pass
    return ps_ring(v, prec=prec)

def eta_poly_relations(eta_elements, degree, labels=['x1','x2'], verbose=False):
    r"""
    Find polynomial relations between eta products.

    INPUTS:

    - ``eta_elements`` - (list): a list of EtaGroupElement objects.
      Not implemented unless this list has precisely two elements. degree

    - ``degree`` - (integer): the maximal degree of polynomial to look for.

    - ``labels`` - (list of strings): labels to use for the polynomial returned.

    - ``verbose``` - (boolean, default False): if True, prints information as
      it goes.

    OUTPUTS: a list of polynomials which is a Groebner basis for the
    part of the ideal of relations between eta_elements which is
    generated by elements up to the given degree; or None, if no
    relations were found.

    ALGORITHM: An expression of the form
    `\sum_{0 \le i,j \le d} a_{ij} x^i y^j` is zero if and
    only if it vanishes at the cusp infinity to degree at least
    `v = d(deg(x) + deg(y))`. For all terms up to `q^v`
    in the `q`-expansion of this expression to be zero is a
    system of `v + k` linear equations in `d^2`
    coefficients, where `k` is the number of nonzero negative
    coefficients that can appear.

    Solving these equations and calculating a basis for the solution
    space gives us a set of polynomial relations, but this is generally
    far from a minimal generating set for the ideal, so we calculate a
    Groebner basis.

    As a test, we calculate five extra terms of `q`-expansion
    and check that this doesn't change the answer.

    EXAMPLES::

        sage: t = EtaProduct(26, {2:2,13:2,26:-2,1:-2})
        sage: u = EtaProduct(26, {2:4,13:2,26:-4,1:-2})
        sage: eta_poly_relations([t, u], 3)
        sage: eta_poly_relations([t, u], 4)
        [x1^3*x2 - 13*x1^3 - 4*x1^2*x2 - 4*x1*x2 - x2^2 + x2]

    Use verbose=True to see the details of the computation::

        sage: eta_poly_relations([t, u], 3, verbose=True)
        Trying to find a relation of degree 3
        Lowest order of a term at infinity = -12
        Highest possible degree of a term = 15
        Trying all coefficients from q^-12 to q^15 inclusive
        No polynomial relation of order 3 valid for 28 terms
        Check: Trying all coefficients from q^-12 to q^20 inclusive
        No polynomial relation of order 3 valid for 33 terms

    ::

        sage: eta_poly_relations([t, u], 4, verbose=True)
        Trying to find a relation of degree 4
        Lowest order of a term at infinity = -16
        Highest possible degree of a term = 20
        Trying all coefficients from q^-16 to q^20 inclusive
        Check: Trying all coefficients from q^-16 to q^25 inclusive
        [x1^3*x2 - 13*x1^3 - 4*x1^2*x2 - 4*x1*x2 - x2^2 + x2]
    """
    if len(eta_elements) > 2:
        raise NotImplementedError("Don't know how to find relations between more than two elements")

    eta1, eta2 = eta_elements

    if verbose: print "Trying to find a relation of degree %s" % degree
    inf = CuspFamily(eta1.level(), 1)
    loterm = -(min([0, eta1.order_at_cusp(inf)]) + min([0,eta2.order_at_cusp(inf)]))*degree
    if verbose: print "Lowest order of a term at infinity = %s" % -loterm

    maxdeg = sum([eta1.degree(), eta2.degree()])*degree
    if verbose: print "Highest possible degree of a term = %s" % maxdeg
    m = loterm + maxdeg + 1
    oldgrob = _eta_relations_helper(eta1, eta2, degree, m, labels, verbose)
    if verbose: print "Check:",
    newgrob = _eta_relations_helper(eta1, eta2, degree, m+5, labels, verbose)
    if oldgrob != newgrob:
        if verbose:
            raise ArithmeticError("Answers different!")
        else:
            raise ArithmeticError("Check: answers different!")
    return newgrob

def _eta_relations_helper(eta1, eta2, degree, qexp_terms, labels, verbose):
    r"""
    Helper function used by eta_poly_relations. Finds a basis for the
    space of linear relations between the first qexp_terms of the
    `q`-expansions of the monomials
    `\eta_1^i * \eta_2^j` for `0 \le i,j < degree`,
    and calculates a Groebner basis for the ideal generated by these
    relations.

    Liable to return meaningless results if qexp_terms isn't at least
    `1 + d*(m_1,m_2)` where

    .. math::

       m_i = min(0, {\text degree of the pole of $\eta_i$ at $\infty$})

    as then 1 will be in the ideal.

    EXAMPLE::

        sage: from sage.modular.etaproducts import _eta_relations_helper
        sage: r,s = EtaGroup(4).basis()
        sage: _eta_relations_helper(r,s,4,100,['a','b'],False)
        [a*b - a + 16]
        sage: _eta_relations_helper(EtaProduct(26, {2:2,13:2,26:-2,1:-2}),EtaProduct(26, {2:4,13:2,26:-4,1:-2}),3,12,['a','b'],False) # not enough terms, will return rubbish
        [1]
    """

    indices = [(i,j) for j in range(degree) for i in range(degree)]
    inf = CuspFamily(eta1.level(), 1)

    pole_at_infinity = -(min([0, eta1.order_at_cusp(inf)]) + min([0,eta2.order_at_cusp(inf)]))*degree
    if verbose: print "Trying all coefficients from q^%s to q^%s inclusive" % (-pole_at_infinity, -pole_at_infinity + qexp_terms - 1)

    rows = []
    for j in xrange(qexp_terms):
        rows.append([])
    for i in indices:
        func = (eta1**i[0]*eta2**i[1]).qexp(qexp_terms)
        for j in xrange(qexp_terms):
            rows[j].append(func[j - pole_at_infinity])
    M = matrix(rows)
    V = M.right_kernel()
    if V.dimension() == 0:
        if verbose: print "No polynomial relation of order %s valid for %s terms" % (degree, qexp_terms)
        return None
    if V.dimension() >= 1:
        #print "Found relation: "
        R = PolynomialRing(QQ, 2, labels)
        x,y = R.gens()
        relations = []
        for c in V.basis():
            relations.append(sum( [ c[v] * x**indices[v][0] * y**indices[v][1] for v in xrange(len(indices))]))
            #print relations[-1], " = 0"
        id = R.ideal(relations)
        return id.groebner_basis()
