r"""
Free String Monoids

AUTHOR: David Kohel <kohel@maths.usyd.edu.au>, 2007-01

SAGE supports a wide range of specific free string monoids.
"""
#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer import Integer
from sage.structure.parent_gens import ParentWithGens, normalize_names
from free_monoid import FreeMonoid_class
from string_monoid_element import StringMonoidElement
from string_ops import strip_encoding

import weakref

_cache = {}

def BinaryStrings():
    r"""
    Returns the free binary string monoid on generators $\{0,1\}$.

    INPUT: None

    OUTPUT:
        Free binary string monoid

    EXAMPLES:
        sage: S = BinaryStrings(); S
        Free binary string monoid
        sage: u = S('')
        sage: u

        sage: x = S('0')
        sage: x
        0
        sage: y = S('1')
        sage: y
        1
        sage: z = S('01110')
        sage: z
        01110
        sage: x*y^3*x == z
        True
        sage: u*x == x*u
        True
    """
    # Here we cache the binary strings to make them unique
    if _cache.has_key(2):
        S = _cache[2]()
        if not S is None:
            return S
    S = BinaryStringMonoid()
    _cache[2] = weakref.ref(S)
    return S

def OctalStrings():
    r"""
    Returns the free octal string monoid on generators $\{0,1,..,7\}$.

    INPUT: None

    OUTPUT:
        Free octal string monoid

    EXAMPLES:
        sage: S = OctalStrings(); S
        Free octal string monoid
        sage: x = S.gens()
        sage: x[0]
        0
        sage: x[7]
        7
        sage: x[0] * x[3]^3 * x[5]^4 * x[6]
        033355556
    """
    # Here we cache the octal strings to make them unique
    if _cache.has_key(8):
        S = _cache[8]()
        if not S is None:
            return S
    S = OctalStringMonoid()
    _cache[8] = weakref.ref(S)
    return S

def HexadecimalStrings():
    r"""
    Returns the free hexadecimal string monoid on generators
    $\{0,1,..,9,a,b,c,d,e,f\}$.

    INPUT: None

    OUTPUT:
        Free hexadecimal string monoid

    EXAMPLES:
        sage: S = HexadecimalStrings(); S
        Free hexadecimal string monoid
        sage: x = S.gen(0)
        sage: y = S.gen(10)
        sage: z = S.gen(15)
        sage: z
        f
        sage: x*y^3*z
        0aaaf
    """
    # Here we cache the hexadecimal strings to make them unique
    if _cache.has_key(16):
        S = _cache[16]()
        if not S is None:
            return S
    S = HexadecimalStringMonoid()
    _cache[16] = weakref.ref(S)
    return S

def Radix64Strings():
    r"""
    Returns the free radix 64 string monoid on 64 generators
    $$
    \{
       A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,
       a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,
       0,1,2,3,4,5,6,7,8,9,+,/
    \}.
    $$

    INPUT: None

    OUTPUT:
        Free radix 64 string monoid

    EXAMPLES:
        sage: S = Radix64Strings(); S
        Free radix 64 string monoid
        sage: x = S.gens()
        sage: x[0]
        A
        sage: x[62]
        +
        sage: x[63]
        /
    """
    # Here we cache the radix-64 strings to make them unique
    if _cache.has_key(64):
        S = _cache[64]()
        if not S is None:
            return S
    S = Radix64StringMonoid()
    _cache[64] = weakref.ref(S)
    return S

def AlphabeticStrings():
    r"""
    Returns the string monoid on generators A-Z:
    $$
    \{ A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z \}
    $$

    INPUT: None

    OUTPUT:
        Free alphabetic string monoid on A-Z.

    EXAMPLES:
        sage: S = AlphabeticStrings(); S
        Free alphabetic string monoid on A-Z
        sage: x = S.gens()
        sage: x[0]
        A
        sage: x[25]
        Z
    """
    # Here we cache the alphabetic strings to make them unique
    if _cache.has_key(26):
        S = _cache[26]()
        if not S is None:
            return S
    S = AlphabeticStringMonoid()
    _cache[26] = weakref.ref(S)
    return S

#*****************************************************************************

class StringMonoid_class(FreeMonoid_class):
    r"""
    A free string monoid on $n$ generators.
    """
    def __init__(self,n,alphabet=()):
        r"""
        Create free binary string monoid on $n$ generators$.

        INPUT:
            n: Integer
            alphabet: String or tuple whose characters or elements denote the generators.

        EXAMPLES:
            sage: S = BinaryStrings(); S
            Free binary string monoid
            sage: x = S.gens()
            sage: x[0]*x[1]**5 * (x[0]*x[1])
            01111101
        """
        # Names must be alphabetical -- omitted since printing is defined locally
        # FreeMonoid_class.__init__(self, n, names = alphabet)
        FreeMonoid_class.__init__(self, n)
        self._alphabet = alphabet

    def __contains__(self, x):
        return isinstance(x, StringMonoidElement) and x.parent() == self

    def alphabet(self):
        return tuple(self._alphabet)

    def gen(self,i=0):
        r"""
        The $i$-th generator of the monoid.

        INPUT:
            i -- integer (default: 0)

        EXAMPLES:
            sage: S = BinaryStrings()
            sage: S.gen(0)
            0
            sage: S.gen(1)
            1
            sage: S.gen(2)
            Traceback (most recent call last):
            ...
            IndexError: Argument i (= 2) must be between 0 and 1.
            sage: S = HexadecimalStrings()
            sage: S.gen(0)
            0
            sage: S.gen(12)
            c
            sage: S.gen(16)
            Traceback (most recent call last):
            ...
            IndexError: Argument i (= 16) must be between 0 and 15.
        """
        n = self.ngens()
        if i < 0 or not i < n:
            raise IndexError, "Argument i (= %s) must be between 0 and %s."%(i,n-1)
        return StringMonoidElement(self,[int(i)])

#*****************************************************************************
# Specific global string monoids
#*****************************************************************************

class BinaryStringMonoid(StringMonoid_class):
    r"""
    The free binary string monoid on generators $\{0,1\}$.
    """
    def __init__(self):
        r"""
        Create free binary string monoid on generators $\{0,1\}$.

        INPUT: None

        EXAMPLES:
            sage: S = BinaryStrings(); S
            Free binary string monoid
            sage: x = S.gens()
            sage: x[0]*x[1]**5 * (x[0]*x[1])
            01111101
        """
        StringMonoid_class.__init__(self, 2, ['0','1'])

    def __cmp__(self, other):
        if not isinstance(other, BinaryStringMonoid):
            return -1
        return 0

    def __repr__(self):
        return "Free binary string monoid"

    def __call__(self, x, check=True):
        r"""
        Return $x$ coerced into this free monoid.

        One can create a free binary string monoid element from a
        Python string of 0's and 1's or list integers

        NOTE: Due to the ambiguity of the second generator '1' with the that the
        identity element '' of the monoid, the syntax S(1) is not permissible.

        EXAMPLES:
            sage: S = BinaryStrings()
            sage: S('101')
            101
            sage: S.gen(0)
            0
            sage: S.gen(1)
            1
        """
        ## There should really some careful type checking here...
        if isinstance(x, StringMonoidElement) and x.parent() == self:
            return x
        elif isinstance(x, list):
            return StringMonoidElement(self, x, check)
        elif isinstance(x, str):
            return StringMonoidElement(self, x, check)
        else:
            raise TypeError, "Argument x (= %s) is of the wrong type."%x

    def encoding(self,S,padic=False):
        r"""
        The binary encoding of the string S, as a binary string element.

        The default is to keep the standard ascii byte encoding, e.g.

            A = 65 -> 01000001
            B = 66 -> 01000010
            .
            .
            .
            Z = 90 -> 01001110

        rather than a 2-adic representation 65 -> 10000010.

        Set padic = True to reverse the bit string.

        EXAMPLES:
            sage: S = BinaryStrings()
            sage: S.encoding('A')
            01000001
            sage: S.encoding('A',padic=True)
            10000010
            sage: S.encoding(' ',padic=True)
            00000100
        """
        from Crypto.Util.number import bytes_to_long
        bit_string = [ ]
        for i in range(len(S)):
            n = int(bytes_to_long(S[i]))
            bits = []
            for i in range(8):
                bits.append(n%2)
                n = n >> 1
            if not padic:
                bits.reverse()
            bit_string.extend(bits)
        return self(bit_string)

class OctalStringMonoid(StringMonoid_class):
    r"""
    The free octal string monoid on generators $\{0,1,..,7\}$.
    """
    def __init__(self):
        r"""
        Create free octal string monoid on generators $\{0,1,..,7\}$.

        INPUT: None

        EXAMPLES:
            sage: S = OctalStrings(); S
            Free octal string monoid
            sage: x = S.gens()
            sage: (x[0]*x[7])**3 * (x[0]*x[1]*x[6]*x[5])**2
            07070701650165
            sage: S([ i for i in range(8) ])
            01234567
        """
        StringMonoid_class.__init__(self, 8, [ str(i) for i in range(8) ])

    def __cmp__(self, other):
        if not isinstance(other, OctalStringMonoid):
            return -1
        return 0

    def __repr__(self):
        return "Free octal string monoid"

    def __call__(self, x, check=True):
        r"""
        Return $x$ coerced into this free monoid.

        One can create a free octal string monoid element from a
        Python string of 0's to 7's or list of integers.

        EXAMPLES:
            sage: S = OctalStrings()
            sage: S('07070701650165')
            07070701650165
            sage: S.gen(0)
            0
            sage: S.gen(1)
            1
            sage: S([ i for i in range(8) ])
            01234567
        """
        ## There should really some careful type checking here...
        if isinstance(x, StringMonoidElement) and x.parent() == self:
            return x
        elif isinstance(x, list):
            return StringMonoidElement(self, x, check)
        elif isinstance(x, str):
            return StringMonoidElement(self, x, check)
        else:
            raise TypeError, "Argument x (= %s) is of the wrong type."%x

class HexadecimalStringMonoid(StringMonoid_class):
    r"""
    The free hexadecimal string monoid on generators $\{0,1,..,9,a,b,c,d,e,f\}$.
    """
    def __init__(self):
        r"""
        Create free hexadecimal string monoid on generators $\{0,1,..,9,a,b,c,d,e,f\}$.

        INPUT: None

        EXAMPLES:
            sage: S = HexadecimalStrings(); S
            Free hexadecimal string monoid
            sage: x = S.gens()
            sage: (x[0]*x[10])**3 * (x[0]*x[1]*x[9]*x[15])**2
            0a0a0a019f019f
            sage: S([ i for i in range(16) ])
            0123456789abcdef
        """
        alph = '0123456789abcdef'
        StringMonoid_class.__init__(self, 16, [ alph[i] for i in range(16) ])

    def __cmp__(self, other):
        if not isinstance(other, HexadecimalStringMonoid):
            return -1
        return 0

    def __repr__(self):
        return "Free hexadecimal string monoid"

    def __call__(self, x, check=True):
        r"""
        Return $x$ coerced into this free monoid.

        One can create a free hexadecimal string monoid element from a
        Python string of a list of integers in $\{0,..,15\}$.

        EXAMPLES:
            sage: S = HexadecimalStrings()
            sage: S('0a0a0a019f019f')
            0a0a0a019f019f
            sage: S.gen(0)
            0
            sage: S.gen(1)
            1
            sage: S([ i for i in range(16) ])
            0123456789abcdef
        """
        ## There should really some careful type checking here...
        if isinstance(x, StringMonoidElement) and x.parent() == self:
            return x
        elif isinstance(x, list):
            return StringMonoidElement(self, x, check)
        elif isinstance(x, str):
            return StringMonoidElement(self, x, check)
        else:
            raise TypeError, "Argument x (= %s) is of the wrong type."%x

    def encoding(self,S,padic=False):
        r"""
        The encoding of the string S, as a hexidecimal string element.

        The default is to keep the standard right-to-left byte encoding, e.g.

            A = '\x41' -> 41
            B = '\x42' -> 42
            .
            .
            .
            Z = '\x5a' -> 5a

        rather than a left-to-right representation A = 65 -> 14.
        Although standard (e.g., in the Python constructor '\xhh'),
        this can be confusing when the string reads left-to-right.

        Set padic = True to reverse the character encoding.

        EXAMPLES:
            sage: S = HexadecimalStrings()
            sage: S.encoding('A')
            41
            sage: S.encoding('A',padic=True)
            14
            sage: S.encoding(' ',padic=False)
            20
            sage: S.encoding(' ',padic=True)
            02
        """
        from Crypto.Util.number import bytes_to_long
        hex_string = [ ]
        for i in range(len(S)):
            n = int(bytes_to_long(S[i]))
            n0 = n % 16; n1 = n // 16
            if not padic:
                hex_chars = [n1,n0]
            else:
                hex_chars = [n0,n1]
            hex_string.extend(hex_chars)
        return self(hex_string)

class Radix64StringMonoid(StringMonoid_class):
    r"""
    The free radix 64 string monoid on 64 generators.
    """
    def __init__(self):
        r"""
        Create free radix 64 string monoid on 64 generators.

        INPUT: None

        EXAMPLES:
            sage: S = Radix64Strings(); S
            Free radix 64 string monoid
            sage: x = S.gens()
            sage: (x[50]*x[10])**3 * (x[60]*x[1]*x[19]*x[35])**2
            yKyKyK8BTj8BTj
            sage: S([ i for i in range(64) ])
            ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/
        """
        alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'
        StringMonoid_class.__init__(self, 64, [ alph[i] for i in range(64) ])


    def __cmp__(self, other):
        if not isinstance(other, Radix64StringMonoid):
            return -1
        return 0

    def __repr__(self):
        return "Free radix 64 string monoid"

    def __call__(self, x, check=True):
        r"""
        Return $x$ coerced into this free monoid.

        One can create a free radix 64 string monoid element from a
        Python string or a list of integers in $0,..,63$, as for
        generic FreeMonoids.

        EXAMPLES:
            sage: S = Radix64Strings()
            sage: S.gen(0)
            A
            sage: S.gen(1)
            B
            sage: S.gen(62)
            +
            sage: S.gen(63)
            /
            sage: S([ i for i in range(64) ])
            ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/
        """
        ## There should really some careful type checking here...
        if isinstance(x, StringMonoidElement) and x.parent() == self:
            return x
        elif isinstance(x, list):
            return StringMonoidElement(self, x, check)
        elif isinstance(x, str):
            return StringMonoidElement(self, x, check)
        else:
            raise TypeError, "Argument x (= %s) is of the wrong type."%x

class AlphabeticStringMonoid(StringMonoid_class):
    """
    The free alphabetic string monoid on generators A-Z.
    """
    def __init__(self):
        r"""
        Create free alphabetic string monoid on generators A-Z.

        INPUT: None

        EXAMPLES:
            sage: S = AlphabeticStrings(); S
            Free alphabetic string monoid on A-Z
            sage: S.gen(0)
            A
            sage: S.gen(25)
            Z
            sage: S([ i for i in range(26) ])
            ABCDEFGHIJKLMNOPQRSTUVWXYZ
        """
        alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        StringMonoid_class.__init__(self, 26, [ alph[i] for i in range(26) ])

    def __cmp__(self, other):
        if not isinstance(other, AlphabeticStringMonoid):
            return -1
        return 0

    def __repr__(self):
        return "Free alphabetic string monoid on A-Z"

    def __call__(self, x, check=True):
        r"""
        Return $x$ coerced into this free monoid.

        One can create a free alphabetic string monoid element from a
        Python string, or a list of integers in $0,..,25$.

        EXAMPLES:
            sage: S = AlphabeticStrings()
            sage: S.gen(0)
            A
            sage: S.gen(1)
            B
            sage: S.gen(25)
            Z
            sage: S([ i for i in range(26) ])
            ABCDEFGHIJKLMNOPQRSTUVWXYZ
        """
        ## There should really some careful type checking here...
        if isinstance(x, StringMonoidElement) and x.parent() == self:
            return x
        elif isinstance(x, list):
            return StringMonoidElement(self, x, check)
        elif isinstance(x, str):
            return StringMonoidElement(self, x, check)
        else:
            raise TypeError, "Argument x (= %s) is of the wrong type."%x

    def encoding(self,S):
        r"""
        The encoding of the string S in the alphabetic string monoid,
        obtained by the monoid homomorphism
        \begin{verbatim}
            A -> A, ..., Z -> Z, a -> A, ..., z -> Z
        \end{verbatim}
        and stripping away all other characters.

        It should be noted that this is a non-injective monoid homomorphism.

        EXAMPLES:
            sage: S = AlphabeticStrings()
            sage: s = S.encoding("The cat in the hat."); s
            THECATINTHEHAT
            sage: s.decoding()
            'THECATINTHEHAT'
        """
        return self(strip_encoding(S))



