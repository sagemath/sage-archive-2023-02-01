"""
p-adic L-functions of elliptic curves

AUTHORS:
   -- William Stein (2007-01-01): first version
"""

######################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
######################################################################


from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.padics.factory import Qp
from sage.rings.infinity import infinity

from sage.rings.integer import Integer
from sage.rings.arith import valuation, binomial

from sage.structure.sage_object import SageObject

from sage.misc.all import verbose

class pAdicLseries(SageObject):
    """
    The p-adic L-series of an elliptic curve.

    EXAMPLES:
    A superingular example:
        sage: e = EllipticCurve('37a')
        sage: L = e.padic_lseries(3); L
        3-adic L-series of Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        sage: L.series(2)
        WARNING: supersingular prime -- answer is *not* to as high of precision as claimed.
        ((3^-1 + O(3^2))*alpha + (2*3^-1 + O(3^2)))*T + ((3^-1 + O(3^2))*alpha + (2*3^-1 + O(3^2)))*T^2 + O(T^3)

    An ordinary example:
        sage: e = EllipticCurve('389a')
        sage: L = e.padic_lseries(5)
        sage: L.series(0)
        Traceback (most recent call last):
        ...
        ValueError: n (=0) must be a positive integer
        sage: L.series(1)
        O(T^1)
        sage: L.series(2)
        (4 + O(5))*T^2 + (2 + O(5))*T^3 + (3 + O(5))*T^4 + O(T^5)
        sage: L.series(3)
        (4 + 4*5 + O(5^2))*T^2 + (2 + 4*5 + O(5^2))*T^3 + (3 + O(5^2))*T^4 + (1 + O(5))*T^5 + (3*5 + O(5^2))*T^6 + (4 + 5 + O(5^2))*T^7 + (2 + 5 + O(5^2))*T^8 + (5 + O(5^2))*T^11 + (2 + O(5^2))*T^12 + (1 + 3*5 + O(5^2))*T^13 + (4 + 4*5 + O(5^2))*T^14 + (3 + O(5))*T^15 + (3*5 + O(5^2))*T^16 + (3*5 + O(5^2))*T^17 + (1 + 4*5 + O(5^2))*T^19 + (3 + O(5))*T^20 + (4 + 3*5 + O(5^2))*T^21 + (2 + O(5^2))*T^22 + (2 + O(5^2))*T^23 + (2 + O(5^2))*T^24 + O(T^25)

    A prime p such that E[p] is reducible:
        sage: L = EllipticCurve('11a').padic_lseries(5)
        sage: L.series(1)
        5 + 4*5^2 + O(5^3) + O(T)
        sage: L.series(2)
        5 + 4*5^2 + 4*5^3 + O(5^4) + O(T^5)
        sage: L.series(3)
        5 + 4*5^2 + 4*5^3 + O(5^5) + (4*5 + O(5^2))*T + (5 + O(5^2))*T^2 + (4*5 + O(5^2))*T^3 + (3*5 + O(5^2))*T^4 + (3*5 + O(5^2))*T^6 + (4*5 + O(5^2))*T^8 + (5 + O(5^2))*T^9 + (5 + O(5^2))*T^11 + (3*5 + O(5^2))*T^13 + (3*5 + O(5^2))*T^14 + (3*5 + O(5^2))*T^16 + (4*5 + O(5^2))*T^18 + (2*5 + O(5^2))*T^21 + (3*5 + O(5^2))*T^22 + (2*5 + O(5^2))*T^23 + (2*5 + O(5^2))*T^24 + O(T^25)

    """
    def __init__(self, E, p, normalize):
        """
        INPUT:
            E -- an elliptic curve
            p -- a prime of good reduction
            normalize -- (bool, default: True); whether or not to correctly
                 normalize the L-series, up to a power of -1 and 2.
                 If False computations may be faster.
        """
        self._E = E
        self._p = ZZ(p)
        self._normalize = normalize
        if not self._p.is_prime():
            raise ValueError, "p (=%s) must be a prime"%p
        #if E.conductor() % self._p == 0:
        #    raise NotImplementedError, "p (=%s) must be a prime of good reduction"%p
        self._modular_symbol = E.modular_symbol(sign=1, normalize=normalize)

    def elliptic_curve(self):
        """
        Return the elliptic curve to which this p-adic L-series is associated.

        EXAMPLES:
            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.elliptic_curve()
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        """
        return self._E

    def prime(self):
        """
        EXAMPLES:
            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.prime()
            5
        """
        return self._p

    def _repr_(self):
        """
        Return print representation.

            sage: e = EllipticCurve('37a')
            sage: e.padic_lseries(3)
            3-adic L-series of Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: e.padic_lseries(3,normalize=False)
            3-adic L-series of Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field (not normalized)
            sage: L = e.padic_lseries(3,normalize=False)
            sage: L.rename('(factor)*L_3(T)')
            sage: L
            (factor)*L_3(T)
        """
        s = "%s-adic L-series of %s"%(self._p, self._E)
        if not self._normalize:
            s += ' (not normalized)'
        return s

    def modular_symbol(self, r):
        """
        Return the modular symbol used to compute this p-adic
        L-series evaluated at r.

        EXAMPLES:
            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: [L.modular_symbol(r) for r in [0,1/5,oo,1/11]]
            [1/5, 6/5, 0, 0]
        """
        return self._modular_symbol(r)

    def measure(self, a, n, prec):
        r"""
        Return the measure on $\ZZ_p^*$ defined by
           $$
             \mu_{E,\alpha}^+ ( a + p^n \ZZ_p  ) =
                   \frac{1}{\alpha^n} \modsym{a}{p^n} - \frac{1}{\alpha^{n+1}} \modsym{a}{p^{n-1}}
           $$
        that is used to define this $p$-adic $L$-function.

        INPUT:
            a -- an integer
            n -- a non-negative integer
            prec -- an integer

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: L = E.padic_lseries(5)
            sage: L.measure(1,2, prec=9)
            1 + 4*5 + 2*5^2 + 4*5^3 + 3*5^4 + 5^5 + 4*5^6 + 4*5^7 + 4*5^8 + O(5^9)
        """
        p = self._p
        alpha = self.alpha(prec)
        z = 1/(alpha**n)
        w = p**(n-1)
        f = self._modular_symbol
        return z * f(a/(p*w)) - (z/alpha) * f(a/w)

    def alpha(self, prec):
        r"""
        Return a p-adic root $\alpha$ of the polynomial $x^2 - a_p x
        + p$ with $\ord_p(\alpha) < 1$.  In the ordinary case this is
        just the unit root.

        INPUT:
            prec -- positive integer, the p-adic precision of the root.

        EXAMPLES:
        Consider the elliptic curve 37a:
            sage: E = EllipticCurve('37a')

        An ordinary prime:
            sage: L = E.padic_lseries(5)
            sage: alpha = L.alpha(10); alpha
            3 + 2*5 + 4*5^2 + 2*5^3 + 5^4 + 4*5^5 + 2*5^7 + 5^8 + 5^9 + O(5^10)
            sage: alpha^2 - E.ap(5)*alpha + 5
            O(5^10)

        A supersingular prime.
            sage: L = E.padic_lseries(3)
            sage: alpha = L.alpha(10); alpha
            (1 + O(3^10))*alpha
            sage: alpha^2 - E.ap(3)*alpha + 3
            0

        A reducible prime:
            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.alpha(5)
            1 + 4*5 + 3*5^2 + 2*5^3 + 4*5^4 + O(5^5)
        """
        try:
            return self._alpha[prec]
        except AttributeError:
            self._alpha = {}
        except KeyError:
            pass
        E = self._E
        p = self._p
        a_p = E.ap(p)
        K = Qp(p, prec, print_mode='series')

        R = ZZ['x']
        f = R([p, -a_p, 1])
        if E.is_ordinary(p):
            G = f.factor_padic(p, prec+5)
            for pr, e in G:
                a = -pr[0]
                if a.valuation() < 1:
                    self._alpha[prec] = K(a)
                    return K(a)
            raise ValueError, "bug in p-adic L-function alpha"
        else: # supersingular case
            f = f.change_ring(Qp(p, prec, print_mode='series'))
            a = f.root_field('alpha', check_irreducible=False).gen()
            self._alpha[prec] = a
            return a

    def order_of_vanishing(self, proof=True):
        """
        Return the order of vanishing of this p-adic L-series.

        The output of this function is provably correct, due to a
        theorem of Kato.  This function will terminate if and only if
        the Mazur-Tate-Teitelbaum analogue of the BSD conjecture about
        the rank of the curve is true and the subgroup of elements of
        p-power order in the Shafarevich-Tate group of this curve is
        finite.  I.e., if this function terminates (with no errors!),
        then you may conclude that the p-adic BSD rank conjecture is
        true and that the p-part of Sha is finite.

        NOTE: currently p must be a prime of good ordinary reduction.

        EXAMPLES:
            sage: L = EllipticCurve('11a').padic_lseries(3)
            sage: L.order_of_vanishing()
            0
            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.order_of_vanishing()
            0
            sage: L = EllipticCurve('37a').padic_lseries(5)
            sage: L.order_of_vanishing()
            1
            sage: L = EllipticCurve('43a').padic_lseries(3)
            sage: L.order_of_vanishing()
            1
            sage: L = EllipticCurve('37b').padic_lseries(3)
            sage: L.order_of_vanishing()
            0

        We verify that Sha(E)(p) is finite for p=3,5,7 for the
        first curve of rank 2:
            sage: e = EllipticCurve('389a')
            sage: for p in primes(3,10):
            ...    print p, e.padic_lseries(p).order_of_vanishing()
            3 2
            5 2
            7 2
        """
        try:
            return self.__ord
        except AttributeError:
            pass

        if not self.is_ordinary():
            raise NotImplementedError
        E = self.elliptic_curve()
        if not E.is_good(self.prime()):
            raise ValueError, "prime must be of good reduction"
        r = E.rank()
        n = 1
        while True:
            f = self.series(n)
            v = f.valuation()
            if v < r:
                raise RuntimeError, "while computing p-adic order of vanishing, got a contradiction: the curve is %s, the curve has rank %s, but the p-adic L-series vanishes to order <= %s"%(E, r, v)
            if v == r:
                self.__ord = v
                return v
            n += 1


    def _c_bounds(self, n):
        raise NotImplementedError

    def _prec_bounds(self, n):
        raise NotImplementedError

    def teichmuller(self, prec):
        r"""
        Return Teichmuller lifts to the given precision.

        INPUT:
            prec -- a positive integer.

        OUTPUT:
            the cached Teichmuller lifts

        EXAMPLES:
            sage: L = EllipticCurve('11a').padic_lseries(7)
            sage: L.teichmuller(1)
            [0, 1, 2, 3, 4, 5, 6]
            sage: L.teichmuller(2)
            [0, 1, 30, 31, 18, 19, 48]
        """
        p = self._p
        K = Qp(p, prec, print_mode='series')
        return [Integer(0)] + \
               [a.residue(prec).lift() for a in K.teichmuller_system()]

    def _e_bounds(self, n):
        p = self._p
        T = (ZZ['T']).gen()
        w = (1+T)**(p**n) - 1
        return [infinity] + [valuation(w[j],p) for j in range(1,w.degree()+1)]


class pAdicLseriesOrdinary(pAdicLseries):
    def series(self, n):
        """
        Return the n-th approximation to the p-adic L-series as a
        power series in T.  Each coefficient is a p-adic number whose
        precision is provably correct (as long as p is a prime of
        ordinary reduction).

        INPUT:
            n -- a positive integer

        EXAMPLES:
        We compute some $p$-adic $L$-functions associated to the elliptic
        curve 11a.

            sage: E = EllipticCurve('11a')
            sage: p = 3
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p)
            sage: L.series(3)
            2 + 3 + 3^2 + 2*3^3 + O(3^5) + (1 + 3 + O(3^2))*T + (1 + 2*3 + O(3^2))*T^2 + (2*3 + O(3^2))*T^4 + (2 + O(3^2))*T^5 + (1 + O(3))*T^6 + (2 + O(3^2))*T^7 + (2 + O(3^2))*T^8 + O(T^9)


        Another example at a prime of bad reduction, where the
        $p$-adic $L$-function has an extra 0 (compared to the non
        $p$-adic $L$-function).

            sage: E = EllipticCurve('11a')
            sage: p = 11
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p)
            sage: L.series(2)
            (10 + O(11))*T + (6 + O(11))*T^2 + (2 + O(11))*T^3 + (5 + O(11))*T^4 + (4 + O(11))*T^5 + (1 + O(11))*T^6 + (5 + O(11))*T^8 + (6 + O(11))*T^9 + (5 + O(11))*T^10 + O(T^11)

        We compute a $p$-adic $L$-function that vanishes to order $2$.
            sage: E = EllipticCurve('389a')
            sage: p = 3
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p)
            sage: L.series(1)
            O(T^1)
            sage: L.series(2)
            (2 + O(3))*T^2 + O(T^3)
            sage: L.series(3)
            (2 + 2*3 + O(3^2))*T^2 + (2 + O(3))*T^3 + (1 + 3 + O(3^2))*T^4 + (1 + O(3^2))*T^5 + (2 + O(3))*T^6 + (1 + 3 + O(3^2))*T^7 + (1 + 3 + O(3^2))*T^8 + O(T^9)
        """
        n = ZZ(n)
        if n < 1:
            raise ValueError, "n (=%s) must be a positive integer"%n

        try:
            return self.__series[n]
        except AttributeError:
            self.__series = {}
        except KeyError:
            pass

        bounds = self._prec_bounds(n)
        prec = max(bounds[1:]) + 5

        p = self._p

        K = QQ
        gamma = K(1 + p)
        R = K[['T']]
        T = R(R.gen(), p**(n-1))
        L = R(0)
        one_plus_T_factor = R(1)
        gamma_power = 1
        teich = self.teichmuller(prec)
        for j in range(p**(n-1)):
            s = K(0)
            for a in range(1,p):
                b = teich[a] * gamma_power
                s += self.measure(b, n, prec).lift()
            L += s * one_plus_T_factor
            one_plus_T_factor *= 1+T
            gamma_power *= gamma

        # Now create series but with each coefficient truncated
        # so it is proven correct:
        K = Qp(p, prec, print_mode='series')
        R = K[['T']]
        L = R(L)
        aj = L.list()
        if len(aj) > 0:
            aj = [aj[0].add_bigoh(prec-2)] + [aj[j].add_bigoh(bounds[j]) for j in range(1,len(aj))]
        L = R(aj, p**(n-1))
        self.__series[n] = L
        return L

    def is_ordinary(self):
        """
        EXAMPLES:
            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.is_ordinary()
            True
        """
        return True

    def is_supersingular(self):
        """
        EXAMPLES:
            sage: L = EllipticCurve('11a').padic_lseries(5)
            sage: L.is_supersingular()
            False
        """
        return False

    def _c_bound(self):
        try:
            return self.__c_bound
        except AttributeError:
            pass
        E = self._E
        p = self._p
        if E.is_irreducible(p):
            ans = 0
        else:
            m = E.modular_symbol_space(sign=1)
            b = m.boundary_map().codomain()
            C = b._known_cusps()  # all known, since computed the boundary map
            ans = max([valuation(self.modular_symbol(a).denominator(), p) for a in C], [0])
            if ans > 0:
                ans = 0
        self.__c_bound = ans
        return ans

    def _prec_bounds(self, n):
        p = self._p
        e = self._e_bounds(n-1)
        c = self._c_bound()
        return [e[j] - c for j in range(len(e))]


class pAdicLseriesSupersingular(pAdicLseries):
    def series(self, n):
        """
        Return the n-th approximation to the p-adic L-series as a
        power series in T.  Each coefficient is a p-adic number whose
        precision is provably correct (as long as p is a prime of
        ordinary reduction).

        INPUT:
            n -- a positive integer

        EXAMPLES:
            sage: L = EllipticCurve('37a').padic_lseries(3)
            sage: L.series(2)
            WARNING: supersingular prime -- answer is *not* to as high of precision as claimed.
            ((3^-1 + O(3^2))*alpha + (2*3^-1 + O(3^2)))*T + ((3^-1 + O(3^2))*alpha + (2*3^-1 + O(3^2)))*T^2 + O(T^3)
            sage: L.alpha(2).parent()
            Univariate Quotient Polynomial Ring in alpha over 3-adic Field with capped relative precision 2 with modulus (1 + O(3^2))*x^2 + (3 + O(3^3))*x + (3 + O(3^3))
        """
        n = ZZ(n)
        if n < 1:
            raise ValueError, "n (=%s) must be a positive integer"%n

        try:
            return self.__series[n]
        except AttributeError:
            self.__series = {}
        except KeyError:
            pass

        bounds = self._prec_bounds(n)
        prec = max(bounds[1:]) + 5

        p = self._p

        alpha = self.alpha(prec)
        K = alpha.parent()
        gamma = 1 + p
        R = K[['T']]
        T = R(R.gen(), p**(n-1))
        L = R(0)
        one_plus_T_factor = R(1)
        gamma_power = 1
        teich = self.teichmuller(prec)
        for j in range(p**(n-1)):
            s = K(0)
            for a in range(1,p):
                b = teich[a] * gamma_power
                s += self.measure(b, n, prec)
            L += s * one_plus_T_factor
            one_plus_T_factor *= 1+T
            gamma_power *= gamma

        # To fix this, just need to use p-adic extension rings, which
        # David Roe has probably not quite finished yet.
        print 'WARNING: supersingular prime -- answer is *not* to as high of precision as claimed.'
        self.__series[n] = L
        return L

    def is_ordinary(self):
        return False

    def is_supersingular(self):
        return True

    def _prec_bounds(self, n):
        p = self._p
        e = self._e_bounds(n-1)
        c = ZZ(n+2)/2
        return [infinity] + [(e[j] - c).floor() for j in range(1,len(e))]


