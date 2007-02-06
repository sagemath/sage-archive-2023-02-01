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
from sage.rings.padic_field import pAdicField

from sage.structure.sage_object import SageObject

class pAdicLseries(SageObject):
    def __init__(self, E, p, prec):
        self._E = E
        self._p = ZZ(p)
        if not self._p.is_prime():
            raise ValueError, "p (=%s) must be a prime"%p
        if E.conductor() % self._p == 0:
            raise NotImplementedError, "p (=%s) must be a prime of good reduction"%p
        self._prec = prec
        self._modular_symbol = E.modular_symbol()

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
                    self._alpha = a
                    return a
            raise ValueError, "bug in p-adic L-function alpha"
        else: # supersingular case
            f = f.change_ring(pAdicField(p, self._prec))
            a = f.root_field('alpha', check_irreducible=False).gen()
            self._alpha = a
            return a

    def approx(self, n, var='T'):
        """
        Return the n-th polynomial approximation to the p-adic L-series
        as a power series in the given variable.

        INPUT:
            n -- an integer
            var -- a string (default: 'T')
        """
        p = self._p
        gamma = 1 + p
        R = pAdicField(p, self._prec)[var]
        T = R.gen()
        L = R(0)
        one_plus_T_factor = R(1)
        gamma_power = 1
        for j in range(p**(n-1)-1):
            s = 0
            for a in range(1,p):
                b = self._teich(a) * gamma_power
                s += self.measure(b, n)
            L += s * one_plus_T_factor
            one_plus_T_factor *= 1+T
            gamma_power *= gamma
        return L

    def _teich(self, a):
        """
        INPUT:
            a -- an integer between 1 and p-1, inclusive
        OUTPUT:
            the Teichmuller lift of a.
        """
        try:
            return self.__teich[a]
        except AttributeError:
            pass
        v = [0]
        # compute a (p-1)st root of unity in Z_p.
        K = pAdicField(p, self._prec)
        zeta = K.zeta(p-1)

        self.__teich = v

class pAdicLseriesOrdinary(pAdicLseries):
    pass

class pAdicLseriesSupersingular(pAdicLseries):
    pass
