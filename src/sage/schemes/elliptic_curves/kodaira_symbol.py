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
    r"""
    Class to hold a Kodaira symbol of an elliptic curve over a $p$-adic local field.

    The standard notation for Kodaira Symbols is as a string which is
    one of Im, II, III, IV, I*m, II*, III*, IV* where m denotes a
    non-negative integer.  These have been encoded by single integers
    by different people.  For convenience we give here the conversion
    table between strings, the pari coding and the eclib encoding.

    Kodaira Symbol        Eclib coding   Pari Coding

    I0                    0              1
    I*0                   1             -1
    Im  (m>0)             10*m           m+4
    I*m (m>0)             10*m+1        -(m+4)
    II, III, IV           2, 3, 4        2,  3,  4
    II*. III*, IV*        7, 6, 5       -2, -3, -4

    """
    def __init__(self, symbol):
        r"""
        Constructor for Kodaira Symbol class.

        INPUT: symbol -- string or integer.  The string should be a
           standard string representation (e.g. 'III*') of a Kodaira
           symbol, which will be parsed.  Alternatively, use the Pari
           encoding of Kodaira symbols as integers.

        EXAMPLES:

        """
        if not isinstance(symbol, str):
            n = Integer(symbol)
            self._n = None
            if n == 0:
                raise ValueError, "Kodaira Symbol code number must be nonzero."
            if n == 1:
                self._n = 0
                self._roman = 1
                self._str = 'I0'
                self._latex = '$I_0$'
            elif n == 2:
                self._roman = 2
                self._str = 'II'
                self._latex = '$II$'
            elif n == 3:
                self._roman = 3
                self._str = 'III'
                self._latex = '$III$'
            elif n == 4:
                self._roman = 4
                self._str = 'IV'
                self._latex = '$IV$'
            elif n > 4:
                nu = n - 4
                self._n = nu
                self._roman = 1
                self._str = 'I' + nu.str()
                self._latex = '$I_{' + nu.str() + '}$'
            elif n == -1:
                self._roman = 1
                self._str = 'I0*'
                self._latex = '$I_0^{*}$'
            elif n == -2:
                self._roman = 2
                self._str = 'II*'
                self._latex = '$II^{*}$'
            elif n == -3:
                self._roman = 3
                self._str = 'III*'
                self._latex = '$III^{*}$'
            elif n == -4:
                self._roman = 4
                self._str = 'IV*'
                self._latex = '$IV^{*}$'
            elif n < -4:
                nu = -n - 4
                self._roman = 1
                self._n = nu
                self._str = 'I' + nu.str() +'*'
                self._latex = '$I_' + nu.str() + '^{*}$'
            self._starred = (n < 0)
            self._pari = n
            return
        elif len(symbol) == 0:
            raise TypeError, "symbol must be a nonempty string"
        if symbol[0] == "I":
            symbol = symbol[1:]
        starred = False
        if symbol[-1] == "*":
            starred = True
            symbol = symbol[:-1]
        self._starred = starred
        if symbol in ["I", "II", "V"]:    # NB we have already stripped off the leading 'I'
            self._roman = ["I", "II", "V"].index(symbol) + 2   # =2, 3 or 4
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
            self._roman = 1
            self._pari = None
            self._n = "generic"
            if starred:
                self._str = "In*"
                self._latex = "$I_n^*$"
            else:
                self._str = "In"
                self._str = "$I_n$"
        elif symbol.isdigit():
            self._roman = 1
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

    def pari_code(self):
        """
        Return the Pari encoding of this Kodaira Symbol.

        EXAMPLES:
            sage: KodairaSymbol('I0').pari_code()
            1
            sage: KodairaSymbol('I10').pari_code()
            14
            sage: KodairaSymbol('I10*').pari_code()
            -14
            sage: [KodairaSymbol(s).pari_code() for s in ['II','III','IV']]
            [2, 3, 4]
            sage: [KodairaSymbol(s).pari_code() for s in ['II*','III*','IV*']]
            [-2, -3, -4]
        """
        return self._pari

_ks_cache = {}
def KodairaSymbol(symbol):
    """
    Returns the specified Kodaira symbol.

    INPUT:
    symbol -- A string of the form "I0", "I1", \ldots, "In", "II", "III", "IV", "I0*", "I1*", \ldots, "In*", "II*", "III*", or "IV*",
              or an integer encoding a Kodaira symbol using Pari's conventions:
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

    EXAMPLES:
        sage: KS = KodairaSymbol
        sage: [KS(n) for n in range(1,10)]
        [I0, II, III, IV, I1, I2, I3, I4, I5]
        sage: [KS(-n) for n in range(1,10)]
        [I0*, II*, III*, IV*, I1*, I2*, I3*, I4*, I5*]
        sage: all([KS(str(KS(n)))==KS(n) for n in range(-10,10) if n!=0])
        True
    """
    if _ks_cache.has_key(symbol):
        ks = _ks_cache[symbol]()
        if not ks is None:
            return ks
    ks = KodairaSymbol_class(symbol)
    _ks_cache[symbol] = weakref.ref(ks)
    return ks
