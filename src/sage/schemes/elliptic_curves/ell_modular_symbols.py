r"""
Modular symbols associated to elliptic curves over the rational numbers

AUTHORS:
   -- William Stein (2007): first version
"""

#*****************************************************************************
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
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.modular.modsym.all import ModularSymbols
from sage.rings.arith import next_prime
from sage.rings.infinity import unsigned_infinity as infinity
from sage.rings.integer import Integer
from sage.modular.cusps import Cusps

oo = Cusps(infinity)
zero = Integer(0)

def modular_symbol_space(E, sign, base_ring, bound=None):
    """
    INPUT:
        E -- elliptic curve
        sign -- integer, -1, 0, or 1
        base_ring -- ring
        bound -- (default: None) maximum number of Hecke operators to
                 use to cut out modular symbols factor.  If None, use
                 enough to provably get the correct answer.

    OUTPUT:
        a space of modular symbols

    EXAMPLES:
        sage: import sage.schemes.elliptic_curves.ell_modular_symbols
        sage: E=EllipticCurve('11a1')
        sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.modular_symbol_space(E,-1,GF(37))
        sage: M
        Modular Symbols space of dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Finite Field of size 37
    """
    _sign = int(sign)
    if _sign != sign:
        raise TypeError, 'sign must be an integer'
    if not (_sign in [-1,0,1]):
        raise TypeError, 'sign must -1, 0, or 1'
    N = E.conductor()
    M = ModularSymbols(N, sign=sign, base_ring=base_ring)
    if bound is None:
        bound = M.hecke_bound() + 10
    V = M
    p = 2
    while p <= bound and V.dimension() > 1:
        t = V.T(p)
        ap = E.ap(p)
        V = (t - ap).kernel()
        p = next_prime(p)

    return V

class ModularSymbol(SageObject):
    r"""
    A modular symbol attached to an elliptic curve, which is the map
    from $\QQ\to \QQ$ obtained by sending $r$ to the normalized
    symmetrized (or anti-symmetrized) integral from r to infinity.

    This is as defined in Mazur-Tate-Teitelbaum.  It's possible the
    map could be off from what you expect by -1 or +/- 2, but
    otherwise it is definitely normalized correctly.

    EXAMPLES:

    """
    def __init__(self, E, sign, normalize=True):
        """
        INPUT:
            E -- an elliptic curve
            sign -- an integer, -1 or 1
            normalize -- (default: True); if True, the modular symbol
                is correctly normalized (up to possibly a factor of
                -1 or 2).  If False, the modular symbol is almost certainly
                not correctly normalized, i.e., all values will be a
                fixed scalar multiple of what they should be.  But
                the initial computation of the modular symbol is
                much faster, though evaluation of it after computing
                it won't be any faster.

        EXAMPLES:
            sage: E=EllipticCurve('11a1')
            sage: import sage.schemes.elliptic_curves.ell_modular_symbols
            sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbol(E,+1)
            sage: M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        """
        _sign = int(sign)
        if _sign != sign:
            raise TypeError, 'sign must be an integer'
        if _sign != -1 and _sign != 1:
            raise TypeError, 'sign must -1 or 1'
        self._E = E
        self._modsym = E.modular_symbol_space(sign=_sign)
        self._ambient_modsym = self._modsym.ambient_module()
        if normalize:
            P = self._modsym.integral_period_mapping()
            e = P.matrix().transpose().row(0)
            e /= 2
        else:
            e = self._modsym.dual_eigenvector()
        self._e = e

    def sign(self):
        """
        Return the sign of this elliptic curve modular symbol.

        EXAMPLES:
            sage: E=EllipticCurve('11a1')
            sage: import sage.schemes.elliptic_curves.ell_modular_symbols
            sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbol(E,+1)
            sage: M.sign()
            1
        """

        return self._modsym.sign()

    def base_ring(self):
        """
        Return the base ring for this modular symbol.
        EXAMPLES:
            sage: E=EllipticCurve('11a1')
            sage: import sage.schemes.elliptic_curves.ell_modular_symbols
            sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbol(E,+1)
            sage: M.base_ring()
            Rational Field
       """
        return self._modsym.base_ring()

    def elliptic_curve(self):
        """
        Return the elliptic curve of this modular symbol
        EXAMPLES:
            sage: E=EllipticCurve('11a1')
            sage: import sage.schemes.elliptic_curves.ell_modular_symbols
            sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbol(E,+1)
            sage: M.elliptic_curve()
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field

        """
        return self._E

    def _call_with_caching(self, r):
        try:
            return self.__cache[r]
        except AttributeError:
            self.__cache = {}
        except KeyError:
            pass
        w = self._ambient_modsym([oo,r]).element()
        c = (self._e).dot_product(w)
        self.__cache[r] = c
        return c

    def __call__(self, r):
        # this next line takes most of the time
        w = self._ambient_modsym.modular_symbol([zero, oo, Cusps(r)], check=False)

        return (self._e).dot_product(w.element())



    def _repr_(self):
        return "Modular symbol with sign %s over %s attached to %s"%(
            self.sign(), self.base_ring(), self.elliptic_curve())

