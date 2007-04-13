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
from sage.rings.padics.qp import Qp
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
        sage: L.approx(2)
        verbose 0 (166: padic_lseries.py, approx) warning: not explicitly truncating p-adic L-function in supersingular case.
        ((3^-1 + O(3^3))*alpha + 2*3^-1 + O(3^3))*T + ((3^-1 + O(3^3))*alpha + 2*3^-1 + O(3^3))*T^2 + O(T^3)

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
            sage:
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
            2 + 3*5 + 4*5^3 + 2*5^4 + 3*5^5 + 3*5^6 + 4*5^7 + 4*5^8 + O(5^9)
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

        EXAMPLES:
        Consider the elliptic curve 37a:
            sage: E = EllipticCurve('37a')

        An ordinary prime:
            sage: L = E.padic_lseries(5)
            sage: alpha = L.alpha(10); alpha
            3 + 2*5 + 4*5^2 + 2*5^3 + 5^4 + 4*5^5 + 2*5^7 + 5^8 + 5^9 + O(5^10)
            sage: alpha^2 - E.ap(5)*alpha + 5
            ??

        A supersingular prime.
            sage: L = E.padic_lseries(3)
            sage: alpha = L.alpha(10); alpha
            ??
            sage: alpha^2 - E.ap(3)*alpha + 3
            ??
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

    def approx(self, n):
        """
        Return the n-th approximation to the p-adic L-series as a
        power series in T.  Each coefficient is a p-adic number whose
        precision is provably correct.

        INPUT:
            n -- an integer

        EXAMPLES:
        We compute some $p$-adic $L$-functions associated to the elliptic
        curve 11a.

            sage: E = EllipticCurve('11a')
            sage: p = 3
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p, prec=10)

            sage: L.approx(4)
            1 + 2*3 + 3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 3^6 + 3^8 + 3^9 + O(3^10) + (2 + 3^3 + 3^5 + 2*3^7 + 3^9 + O(3^10))*T + (2 + 2*3 + 2*3^2 + 3^3 + 2*3^4 + 2*3^5 + 2*3^7 + 3^8 + O(3^10))*T^2 + (2*3 + 3^2 + 3^5 + 2*3^6 + 3^8 + 2*3^9 + O(3^10))*T^3 + (3 + 2*3^2 + 3^5 + 2*3^6 + 2*3^7 + 3^8 + O(3^10))*T^4 + O(T^5)
            sage: L.approx(5)
            1 + 2*3 + 3^2 + 2*3^4 + 3^5 + 2*3^6 + O(3^10) + (2 + 3^5 + 3^6 + 2*3^7 + 3^9 + O(3^10))*T + (2 + 2*3 + 2*3^3 + 2*3^4 + 2*3^6 + 3^7 + 3^9 + O(3^10))*T^2 + (2*3 + 2*3^2 + 2*3^3 + 2*3^5 + 3^6 + 2*3^8 + 2*3^9 + O(3^10))*T^3 + (3 + 2*3^2 + 2*3^5 + 2*3^6 + 3^8 + 2*3^9 + O(3^10))*T^4 + (1 + 3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^7 + 3^8 + O(3^10))*T^5 + O(T^6)

        Another example at a prime of bad reduction, where the
        $p$-adic $L$-function has an extra 0 (compared to the non
        $p$-adic $L$-function).

            sage: E = EllipticCurve('11a')
            sage: p = 11
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p, prec=10)

            sage: L.approx(2)
            10*11^2 + 7*11^3 + 9*11^5 + 6*11^8 + 10*11^9 + O(11^10) + (6 + 8*11 + 9*11^2 + 11^3 + 11^4 + 5*11^5 + 5*11^6 + 5*11^7 + 4*11^8 + 3*11^9 + O(11^10))*T + (8 + 9*11 + 3*11^2 + 7*11^3 + 9*11^4 + 8*11^6 + 4*11^7 + 6*11^8 + 10*11^9 + O(11^10))*T^2 + O(T^3)

        We compute a $p$-adic $L$-function that vanishes to order $2$.
            sage: E = EllipticCurve('389a')
            sage: p = 3
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p,prec=10)

            sage: L.approx(1)
            0
            sage: L.approx(2)
            (3 + 3^2 + 2*3^3 + 3^4 + 3^5 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10))*T + (2 + 2*3 + 3^2 + 3^5 + 3^6 + 2*3^7 + 2*3^8 + O(3^10))*T^2 + O(T^3)
            sage: L.approx(3)
            (2*3^2 + 3^5 + 3^6 + 2*3^7 + 2*3^8 + O(3^10))*T + (2 + 2*3 + 2*3^4 + 2*3^5 + 3^6 + 3^7 + 2*3^8 + O(3^10))*T^2 + (2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 3^5 + O(3^10))*T^3 + O(T^4)
            sage: L.approx(4)
            (2*3^3 + 3^6 + 3^7 + 2*3^8 + 2*3^9 + O(3^10))*T + (2 + 2*3 + 2*3^2 + 3^3 + 3^4 + 2*3^5 + 2*3^6 + 2*3^7 + O(3^10))*T^2 + (2 + 2*3^2 + 3^6 + 3^8 + O(3^10))*T^3 + (1 + 2*3 + 3^3 + 2*3^4 + 2*3^5 + 3^7 + 3^8 + O(3^10))*T^4 + O(T^5)

        """
        try:
            return self.__approx[n]
        except AttributeError:
            self.__approx = {}
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
        if self.is_ordinary():
            aj = L.list()
            aj = [aj[0]] + [aj[j].add_bigoh(bounds[j]) for j in range(1,len(aj))]
            L = R(aj, p**(n-1))
        else:
            # To fix this, just need to use p-adic extension rings, which
            # David Roe has probably not quite finished yet.
            verbose('warning: not explicitly truncating p-adic L-function in supersingular case.',level=0)
        self.__approx[n] = L
        return L

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
    def is_ordinary(self):
        return True

    def is_supersingular(self):
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
    def is_ordinary(self):
        return False

    def is_supersingular(self):
        return True

    def _prec_bounds(self, n):
        p = self._p
        e = self._e_bounds(n-1)
        c = ZZ(n+2)/2
        return [infinity] + [(e[j] - c).floor() for j in range(1,len(e))]


