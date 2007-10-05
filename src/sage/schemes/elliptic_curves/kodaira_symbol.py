"""
Kodaira symbols.
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe       <roed@math.harvard.edu>
#                          William Stein   <wstein@gmail.com>
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
from sage.rings.integer import Integer
import weakref

class KodairaSymbol_class(SageObject):
    def __init__(self, symbol):
        if not isinstance(symbol, str):
            n = Integer(symbol)
            self._n = None
            if n == 1:
                self._n = 0
                self._str = 'I0'
                self._latex = '$I_0$'
            elif n == 2:
                self._str = 'II'
                self._latex = '$II$'
            elif n == 3:
                self._str = 'III'
                self._latex = '$III$'
            elif n == 4:
                self._str = 'IV'
                self._latex = '$IV$'
            elif n > 4:
                nu = n - 4
                self._n = nu
                self._str = 'I' + nu.str()
                self._latex = '$I_{' + nu.str() + '}$'
            elif n == -1:
                self._str = 'I0*'
                self._latex = '$I_0^{*}$'
            elif n == -2:
                self._str = 'II*'
                self._latex = '$II^{*}$'
            elif n == -3:
                self._str = 'III*'
                self._latex = '$III^{*}$'
            elif n == -4:
                self._str = 'IV*'
                self._latex = '$IV^{*}$'
            elif n < -4:
                nu = -n - 4
                self._n = nu
                self._str = 'I' + nu.str() +'*'
                self._latex = '$I_' + nu.str() + '^{*}$'
            self._starred = (n < 0)
            self._pari = n
            return
        elif len(symbol) == 0:
            raise TypeError, "symbol must be a nonemptystring"
        if symbol[0] == "I":
            symbol = symbol[1:]
        starred = False
        if symbol[-1] == "*":
            starred = True
            symbol = symbol[:-1]
        self._starred = starred
        if symbol in ["I", "II", "V"]:
            self._n = None
            if starred:
                sign = -1
                self._str = "I" + symbol + "*"
                self._latex = "$I" + symbol + "^*$"
            else:
                sign = 1
                self._str = "I" + symbol
                self._latex = "$" + self._str + "$"
            if symbol == "I":
                self._pari = 2 * sign
            elif symbol == "II":
                self._pari = 3 * sign
            elif symbol == "V":
                self._pari = 4 * sign
        elif symbol == "n":
            self._pari = None
            self._n = "generic"
            if starred:
                self._str = "In*"
                self._latex = "$I_n^*$"
            else:
                self._str = "In"
                self._str = "$I_n$"
        elif symbol.isdigit():
            self._n = Integer(symbol)
            if starred:
                if self._n == 0:
                    self._pari = -1
                else:
                    self._pari = -self._n - 4
                self._str = "I" + symbol + "*"
                self._latex = "$I_{%s}^*$"%(symbol)
            else:
                if self._n == 0:
                    self._pari = 1
                else:
                    self._pari = self._n + 4
                self._str = "I" + symbol
                self._latex = "$I_{%s}$"%(symbol)
        else:
            raise ValueError, "input is not a Kodaira symbol"

    def __repr__(self):
        return self._str

    def _latex_(self):
        return self._latex

    def __cmp__(self, other):
        if isinstance(other, KodairaSymbol_class):
            if (self._n == "generic" and not other._n is None) or (other._n == "generic" and not self._n is None):
                return cmp(self._starred, other._starred)
            return cmp(self._str, other._str)
        else:
            return cmp(type(self), type(other))

_ks_cache = {}
def KodairaSymbol(symbol):
    """
    Returns the specified Kodaira symbol.

    INPUT:
    symbol -- A string of the form "I0", "I1", \ldots, "In", "II", "III", "IV", "I0*", "I1*", \ldots, "In*", "II*", "III*", or "IV*",
              or an integer encoding a Kodaia symbol using Pari's conventions:
              1 = "I0" (good reduction)
              2 = "II"
              3 = "III"
              4 = "IV"
              4+n = "I_n"
              -1 = "I0*"
              -2 = "II*"
              -3 = "III*"
              -4 = "IV*"
              -4-n = "I_n^*"
    OUTPUT:
    KodairaSymbol -- the corresponding Kodaira symbol.
    """
    if _ks_cache.has_key(symbol):
        ks = _ks_cache[symbol]()
        if not ks is None:
            return ks
    ks = KodairaSymbol_class(symbol)
    _ks_cache[symbol] = weakref.ref(ks)
    return ks
