"""
String Monoid Elements

AUTHORS:

- David Kohel <kohel@maths.usyd.edu.au>, 2007-01

Elements of free string monoids, internal representation subject to change.

These are special classes of free monoid elements with distinct printing.

The internal representation of elements does not use the exponential
compression of FreeMonoid elements (a feature), and could be packed into words.
"""

#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# import operator
from sage.rings.integer import Integer
from sage.rings.all import RealField
# from sage.structure.element import MonoidElement
from sage.probability.random_variable import DiscreteProbabilitySpace
from free_monoid_element import FreeMonoidElement
import string_monoid

def is_StringMonoidElement(x):
    r"""
    """
    return isinstance(x, StringMonoidElement)

def is_AlphabeticStringMonoidElement(x):
    r"""
    """
    return isinstance(x, StringMonoidElement) and \
           isinstance(x.parent(), string_monoid.AlphabeticStringMonoid)

def is_BinaryStringMonoidElement(x):
    r"""
    """
    return isinstance(x, StringMonoidElement) and \
           isinstance(x.parent(), string_monoid.BinaryStringMonoid)

def is_OctalStringMonoidElement(x):
    r"""
    """
    return isinstance(x, StringMonoidElement) and \
           isinstance(x.parent(), string_monoid.OctalStringMonoid)

def is_HexadecimalStringMonoidElement(x):
    r"""
    """
    return isinstance(x, StringMonoidElement) and \
           isinstance(x.parent(), string_monoid.HexadecimalStringMonoid)

def is_Radix64StringMonoidElement(x):
    r"""
    """
    return isinstance(x, StringMonoidElement) and \
           isinstance(x.parent(), string_monoid.Radix64StringMonoid)


class StringMonoidElement(FreeMonoidElement):
    """
    Element of a free string monoid.
    """

    def __init__(self, S, x, check=True):
        """
        Create the element ``x`` of the StringMonoid ``S``.

        This should typically be called by a StringMonoid.
        """
        FreeMonoidElement.__init__(self, S, [])
        if isinstance(x, list):
            if check:
                for b in x:
                    if not isinstance(b, (int, long, Integer)):
                        raise TypeError(
                            "x (= %s) must be a list of integers." % x)
            self._element_list = list(x) # make copy
        elif isinstance(x, str):
            alphabet = list(self.parent().alphabet())
            self._element_list = []
            for i in range(len(x)):
                try:
                    b = alphabet.index(x[i])
                except ValueError:
                    raise TypeError(
                        "Argument x (= %s) is not a valid string." % x)
                self._element_list += [b]
        else:
            raise TypeError("Argument x (= %s) is of the wrong type." % x)

    def __cmp__(left, right):
        """
        Compare two free monoid elements with the same parents.

        The ordering is the one on the underlying sorted list of
        (monomial,coefficients) pairs.

        EXAMPLES::

            sage: S = BinaryStrings()
            sage: (x,y) = S.gens()
            sage: x * y < y * x
            True
            sage: S("01") < S("10")
            True
        """
        return cmp(left._element_list, right._element_list)

    def _repr_(self):
        """
        The self-representation of a string monoid element. Unlike Python
        strings we print without the enclosing quotes.
        """
        S = self.parent()
        alphabet = S.alphabet()
        s = ''
        for b in self._element_list:
            c = alphabet[b]
            s += c
        return s

    def _latex_(self):
        """
        Return latex representation of self.

        EXAMPLES::

            sage: S = BinaryStrings()
            sage: s = S('101111000')
            sage: latex(s)
            101111000
        """
        return self._repr_()

    def __mul__(self, y):
        """
        Multiply 2 free string monoid elements.

        EXAMPLES::

            sage: S = BinaryStrings()
            sage: (x,y) = S.gens()
            sage: x*y
            01
        """
        if not isinstance(y, StringMonoidElement):
            raise TypeError("Argument y (= %s) is of wrong type." % y)
        S = self.parent()
        x_elt = self._element_list
        y_elt = y._element_list
        z = S('')
        z._element_list = x_elt + y_elt
        return z

    def __pow__(self, n):
        """
        Return the `n`-th power of the string element.

        EXAMPLES::

            sage: (x,y) = BinaryStrings().gens()
            sage: x**3 * y**5 * x**7
            000111110000000
            sage: x**0


        Note that raising to a negative power is *not* a constructor
        for an element of the corresponding free group (yet).

        ::

            sage: x**(-1)
            Traceback (most recent call last):
            ...
            IndexError: Argument n (= -1) must be non-negative.
        """
        if not isinstance(n, (int, long, Integer)):
            raise TypeError("Argument n (= %s) must be an integer." % n)
        if n < 0:
            raise IndexError("Argument n (= %s) must be non-negative." % n)
        elif n == 0:
            return self.parent()('')
        elif n == 1:
            return self
        z = self.parent()('')
        z._element_list = self._element_list * n
        return z

    def __len__(self):
        """
        Return the number of products that occur in this monoid element.
        For example, the length of the identity is 0, and the length
        of the monoid `x_0^2x_1` is three.

        EXAMPLES::

            sage: S = BinaryStrings()
            sage: z = S('')
            sage: len(z)
            0
            sage: (x,y) = S.gens()
            sage: len(x**2 * y**3)
            5
        """
        return len(self._element_list)

    def __getitem__(self, n):
        """
        Return the n-th string character.
        """
        try:
            c = self._element_list[n]
        except:
            raise IndexError("Argument n (= %s) is not a valid index." % n)
        if not isinstance(c, list):
            c = [c]
        return self.parent()(c)

    def decoding(self, padic=False):
        r"""
        The byte string associated to a binary or hexadecimal string
        monoid element.

        EXAMPLES::

            sage: S = HexadecimalStrings()
            sage: s = S.encoding("A..Za..z"); s
            412e2e5a612e2e7a
            sage: s.decoding()
            'A..Za..z'
            sage: s = S.encoding("A..Za..z",padic=True); s
            14e2e2a516e2e2a7
            sage: s.decoding()
            '\x14\xe2\xe2\xa5\x16\xe2\xe2\xa7'
            sage: s.decoding(padic=True)
            'A..Za..z'
            sage: S = BinaryStrings()
            sage: s = S.encoding("A..Za..z"); s
            0100000100101110001011100101101001100001001011100010111001111010
            sage: s.decoding()
            'A..Za..z'
            sage: s = S.encoding("A..Za..z",padic=True); s
            1000001001110100011101000101101010000110011101000111010001011110
            sage: s.decoding()
            '\x82ttZ\x86tt^'
            sage: s.decoding(padic=True)
            'A..Za..z'
        """
        S = self.parent()
        from Crypto.Util.number import long_to_bytes
        if isinstance(S, string_monoid.AlphabeticStringMonoid):
            return ''.join([ long_to_bytes(65+i) for i in self._element_list ])
        n = self.__len__()
        if isinstance(S, string_monoid.HexadecimalStringMonoid):
            if not n % 2 == 0:
                "String %s must have even length to determine a byte character string." % str(self)
            s = []
            x = self._element_list
            for k in range(n//2):
                m = 2*k
                if padic:
                    c = long_to_bytes(x[m]+16*x[m+1])
                else:
                    c = long_to_bytes(16*x[m]+x[m+1])
                s.append(c)
            return ''.join(s)
        if isinstance(S, string_monoid.BinaryStringMonoid):
            if not n % 8 == 0:
                "String %s must have even length 0 mod 8 to determine a byte character string." % str(self)
            pows = [ 2**i for i in range(8) ]
            s = []
            x = self._element_list
            for k in range(n//8):
                m = 8*k
                if padic:
                    c = long_to_bytes(sum([ x[m+i]*pows[i] for i in range(8) ]))
                else:
                    c = long_to_bytes(sum([ x[m+7-i]*pows[i] for i in range(8) ]))
                s.append(c)
            return ''.join(s)
        raise TypeError(
            "Argument %s must be an alphabetic, binary, or hexadecimal string." % str(self))

    def coincidence_index(self, prec=0):
        """
        Returns the probability of two randomly chosen characters being equal.
        """
        if prec == 0:
            RR = RealField()
        else:
            RR = RealField(prec)
        char_dict = {}
        for i in self._element_list:
            # if char_dict.has_key(i):
            # The method .has_key() has been deprecated since Python 2.2. Use
            # "k in Dict" instead of "Dict.has_key(k)".
            if i in char_dict:
                char_dict[i] += 1
            else:
                char_dict[i] = 1
        nn = 0
        ci_num = 0
        for i in char_dict.keys():
            ni = char_dict[i]
            nn += ni
            ci_num += ni*(ni-1)
        ci_den = nn*(nn-1)
        return RR(ci_num)/ci_den

    def frequency_distribution(self, length=1, prec=0):
        """
        Returns the probability space of character frequencies.
        """
        if not length in (1, 2):
            raise NotImplementedError("Not implemented")
        if prec == 0:
            RR = RealField()
        else:
            RR = RealField(prec)
        S = self.parent()
        n = S.ngens()
        if length == 1:
            Alph = S.gens()
        else:
            Alph = tuple([ x*y for x in S.gens() for y in S.gens() ])
        X = {}
        N = len(self)-length+1
        eps = RR(Integer(1)/N)
        for i in range(N):
            c = self[i:i+length]
            # if X.has_key(c):
            # The method .has_key() has been deprecated since Python 2.2. Use
            # "k in Dict" instead of "Dict.has_key(k)".
            if c in X:
                X[c] += eps
            else:
                X[c] = eps
        return DiscreteProbabilitySpace(Alph, X, RR)
