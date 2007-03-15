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

    return V

class ModularSymbol(SageObject):
    def __init__(self, E, sign, base_ring):
        """
        INPUT:
            E -- an elliptic curve
            sign -- an integer, -1 or 1
            base_ring -- a ring
        """
        _sign = int(sign)
        if _sign != sign:
            raise TypeError, 'sign must be an integer'
        if _sign != -1 and _sign != 1:
            raise TypeError, 'sign must -1 or 1'
        self._E = E
        self._modsym = E.modular_symbol_space(sign=_sign, base_ring=base_ring)
        self._e = self._modsym.dual_eigenvector()
        # todo -- here must rescale self._e

    def sign(self):
        return self._modsym.sign()

    def base_ring(self):
        return self._modsym.base_ring()

    def elliptic_curve(self):
        return self._E

    def __call__(self, x):
        w = self._modsym([0,x]).element()
        return (self._e).dot_product(w)

    def _repr_(self):
        return "Modular symbol with sign %s over %s attached to %s"%(
            self.sign(), self.base_ring(), self.elliptic_curve())

