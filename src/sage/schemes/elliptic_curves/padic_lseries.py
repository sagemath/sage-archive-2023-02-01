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
from sage.rings.padics.qp import Qp

from sage.structure.sage_object import SageObject

class pAdicLseries(SageObject):
    def __init__(self, E, p, prec):
        self._E = E
        self._p = ZZ(p)
        if not self._p.is_prime():
            raise ValueError, "p (=%s) must be a prime"%p
        #if E.conductor() % self._p == 0:
        #    raise NotImplementedError, "p (=%s) must be a prime of good reduction"%p
        self._prec = prec
        self._modular_symbol = E.modular_symbol()
        self._Qp = Qp(p, self._prec, print_mode='series')

    def _repr_(self):
        return "%s-adic L-series of %s"%(self._p, self._E)

    def modular_symbol(self, x):
        return self._modular_symbol(x)

    def measure(self, a, n):
        r"""
        Return the measure on $\ZZ_p^*$ defined by
           $$
             \mu_{E,\alpha}^+ ( a + p^n \ZZ_p  ) =
                   \frac{1}{\alpha^n} \modsym{a}{p^n} - \frac{1}{\alpha^{n+1}} \modsym{a}{p^{n-1}}
           $$
        that is used to define this $p$-adic $L$-function.

            sage: E = EllipticCurve('37a')
            sage: L = E.padic_lseries(5, prec=9)
            sage: L.measure(1,2)
            2 + 3*5 + 4*5^3 + 2*5^4 + 3*5^5 + 3*5^6 + 4*5^7 + 4*5^8 + O(5^9)
        """
        p = self._p
        alpha = self.alpha()
        z = 1/(alpha**n)
        w = p**(n-1)
        f = self._modular_symbol
        return z * f(a/(p*w)) - (z/alpha) * f(a/w)

    def alpha(self):
        r"""
        Return the p-adic root $\alpha$ of the polynomial $x^2 - a_p x
        + p$ with $\ord_p(\alpha) < 1$.

        EXAMPLES:
            sage: E = EllipticCurve('37a')
            sage: L = E.padic_lseries(5, prec=10)
            sage: L.alpha()
            3 + 2*5 + 4*5^2 + 2*5^3 + 5^4 + 4*5^5 + 2*5^7 + 5^8 + 5^9 + O(5^10)
        """
        try:
            return self._alpha
        except AttributeError:
            pass
        E = self._E
        p = self._p
        a_p = E.ap(p)

        R = ZZ['x']
        f = R([p, -a_p, 1])
        if E.is_ordinary(p):
            G = f.factor_padic(p, self._prec)
            for pr, e in G:
                a = -pr[0]
                if a.valuation() < 1:
                    self._alpha = self._Qp(a)
                    return self._Qp(a)
            raise ValueError, "bug in p-adic L-function alpha"
        else: # supersingular case
            f = f.change_ring(Qp(p, self._prec, print_mode='series'))
            a = f.root_field('alpha', check_irreducible=False).gen()
            self._alpha = a
            return a

    def approx(self, n):
        """
        Return the n-th approximation to the p-adic L-series as a
        power series in T.

        INPUT:
            n -- an integer

        EXAMPLES:
        We compute some $p$-adic $L$-functions associated to the elliptic curve 11a.
            sage: E = EllipticCurve('11a')
            sage: p = 3
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p, prec=10)

            #sage: L.approx(4)
            #1 + 2*3 + 3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 3^6 + 3^8 + 3^9 + O(3^10) + (2 + 3^3 + 3^5 + 2*3^7 + 3^9 + O(3^10))*T + (2 + 2*3 + 2*3^2 + 3^3 + 2*3^4 + 2*3^5 + 2*3^7 + 3^8 + O(3^10))*T^2 + (2*3 + 3^2 + 3^5 + 2*3^6 + 3^8 + 2*3^9 + O(3^10))*T^3 + (3 + 2*3^2 + 3^5 + 2*3^6 + 2*3^7 + 3^8 + O(3^10))*T^4 + O(T^5)
            #sage: L.approx(5)
            #1 + 2*3 + 3^2 + 2*3^4 + 3^5 + 2*3^6 + O(3^10) + (2 + 3^5 + 3^6 + 2*3^7 + 3^9 + O(3^10))*T + (2 + 2*3 + 2*3^3 + 2*3^4 + 2*3^6 + 3^7 + 3^9 + O(3^10))*T^2 + (2*3 + 2*3^2 + 2*3^3 + 2*3^5 + 3^6 + 2*3^8 + 2*3^9 + O(3^10))*T^3 + (3 + 2*3^2 + 2*3^5 + 2*3^6 + 3^8 + 2*3^9 + O(3^10))*T^4 + (1 + 3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^7 + 3^8 + O(3^10))*T^5 + O(T^6)

        Another example at a prime of bad reduction, where the
        $p$-adic $L$-function has an extra 0 (compared to the non
        $p$-adic $L$-function).

            sage: E = EllipticCurve('11a')
            sage: p = 11
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p, prec=10)

            #sage: L.approx(2)
            #10*11^2 + 7*11^3 + 9*11^5 + 6*11^8 + 10*11^9 + O(11^10) + (6 + 8*11 + 9*11^2 + 11^3 + 11^4 + 5*11^5 + 5*11^6 + 5*11^7 + 4*11^8 + 3*11^9 + O(11^10))*T + (8 + 9*11 + 3*11^2 + 7*11^3 + 9*11^4 + 8*11^6 + 4*11^7 + 6*11^8 + 10*11^9 + O(11^10))*T^2 + O(T^3)

        We compute a $p$-adic $L$-function that vanishes to order $2$.
            sage: E = EllipticCurve('389a')
            sage: p = 3
            sage: E.is_ordinary(p)
            True
            sage: L = E.padic_lseries(p,prec=10)

            #sage: L.approx(1)
            #0
            #sage: L.approx(2)
            #(3 + 3^2 + 2*3^3 + 3^4 + 3^5 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10))*T + (2 + 2*3 + 3^2 + 3^5 + 3^6 + 2*3^7 + 2*3^8 + O(3^10))*T^2 + O(T^3)
            #sage: L.approx(3)
            #(2*3^2 + 3^5 + 3^6 + 2*3^7 + 2*3^8 + O(3^10))*T + (2 + 2*3 + 2*3^4 + 2*3^5 + 3^6 + 3^7 + 2*3^8 + O(3^10))*T^2 + (2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 3^5 + O(3^10))*T^3 + O(T^4)
            #sage: L.approx(4)
            #(2*3^3 + 3^6 + 3^7 + 2*3^8 + 2*3^9 + O(3^10))*T + (2 + 2*3 + 2*3^2 + 3^3 + 3^4 + 2*3^5 + 2*3^6 + 2*3^7 + O(3^10))*T^2 + (2 + 2*3^2 + 3^6 + 3^8 + O(3^10))*T^3 + (1 + 2*3 + 3^3 + 2*3^4 + 2*3^5 + 3^7 + 3^8 + O(3^10))*T^4 + O(T^5)

        """
        try:
            return self.__approx[n]
        except AttributeError:
            self.__approx = {}
        except KeyError:
            pass

        p = self._p
        gamma = 1 + p
        R = self._Qp[['T']]
        T = R(R.gen(), n+1)
        L = R(0)
        one_plus_T_factor = R(1)
        gamma_power = 1
        teich = self.teichmuller(n)
        for j in range(p**(n-1)):
            s = 0
            for a in range(1,p):
                b = teich[a] * gamma_power
                s += self.measure(b, n)
            L += s * one_plus_T_factor
            one_plus_T_factor *= 1+T
            gamma_power *= gamma
        self.__approx[n] = L
        return L

    def teichmuller(self, n):
        """
        INPUT:
            n -- a positive integer.

        OUTPUT:
            the Teichmuller lifts to precision p^n as integers.
        """
        K = self._Qp
        pn = self._p ** n
        return [0] + [(K.teichmuller(a).residue(n)).lift() for a in range(1,self._p)]

class pAdicLseriesOrdinary(pAdicLseries):
    pass

class pAdicLseriesSupersingular(pAdicLseries):
    pass
