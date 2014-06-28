#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
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

include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"
include 'misc.pxi'
include 'decl.pxi'

from sage.rings.integer cimport Integer

from ntl_ZZ import unpickle_class_value
from ntl_GF2 cimport ntl_GF2


##############################################################################
#
# ntl_GF2X: Polynomials over GF(2) via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2006-01: initial version (based on code by William Stein)
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2007-10: adapted to new conventions
#
##############################################################################

def GF2XHexOutput(have_hex=None):
    """
    Represent GF2X and GF2E elements in the more compact
    hexadecimal form to the user.

    If no parameter is provided the currently set value will be
    returned.

    INPUT:
        have_hex if True hex representation will be used

    EXAMPLES:
        sage: m = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
        sage: x = ntl.GF2E([1,0,1,0,1], m) ; x
        [1 0 1 0 1]

        sage: ntl.GF2XHexOutput() ## indirect doctest
        False
        sage: ntl.GF2XHexOutput(True)
        sage: ntl.GF2XHexOutput()
        True

        sage: x
        0x51

        sage: ntl.GF2XHexOutput(False)
        sage: x
        [1 0 1 0 1]
    """
    if have_hex is None:
        return bool(GF2XHexOutput_c[0])

    if have_hex:
        GF2XHexOutput_c[0] = 1
    else:
        GF2XHexOutput_c[0] = 0

cdef class ntl_GF2X:
    """
    Univariate Polynomials over GF(2) via NTL.
    """
    def __init__(self, x=[]):
        """
        Constructs a new polynomial over GF(2).

        A value may be passed to this constructor. If you pass a string
        to the constructor please note that byte sequences and the hexadecimal
        notation are little endian.  So e.g. '[0 1]' == '0x2' == x.

        Input types are ntl.ZZ_px, strings, lists of digits, FiniteFieldElements
        from extension fields over GF(2), Polynomials over GF(2), Integers, and finite
        extension fields over GF(2) (uses modulus).

        INPUT:
            x -- value to be assigned to this element. See examples.

        OUTPUT:
            a new ntl.GF2X element

        EXAMPLES:
            sage: ntl.GF2X(ntl.ZZ_pX([1,1,3],2))
            [1 1 1]
            sage: ntl.GF2X('0x1c')
            [1 0 0 0 0 0 1 1]
            sage: ntl.GF2X('[1 0 1 0]')
            [1 0 1]
            sage: ntl.GF2X([1,0,1,0])
            [1 0 1]
            sage: ntl.GF2X(GF(2**8,'a').gen()**20)
            [0 0 1 0 1 1 0 1]
            sage: ntl.GF2X(GF(2**8,'a'))
            [1 0 1 1 1 0 0 0 1]
            sage: ntl.GF2X(2)
            [0 1]
            sage: ntl.GF2X(ntl.GF2(1))
            [1]

            sage: R.<x> = GF(2)[]
            sage: f = x^5+x^2+1
            sage: ntl.GF2X(f)
            [1 0 1 0 0 1]
        """

        from sage.rings.finite_rings.element_ext_pari import FiniteField_ext_pariElement
        from sage.rings.finite_rings.element_givaro import FiniteField_givaroElement
        from sage.rings.finite_rings.element_ntl_gf2e import FiniteField_ntl_gf2eElement
        from sage.rings.finite_rings.finite_field_base import FiniteField
        from sage.rings.polynomial.polynomial_gf2x import Polynomial_GF2X

        cdef long _x

        if PY_TYPE_CHECK(x,ntl_GF2):
            GF2X_conv_GF2(self.x,(<ntl_GF2>x).x)
            return
        elif PY_TYPE_CHECK(x,ntl_GF2X):
            self.x = (<ntl_GF2X>x).x
            return
        elif PY_TYPE_CHECK(x,int):
            _x = x
            GF2XFromBytes(self.x, <unsigned char *>(&_x),sizeof(long))
            return

        if PY_TYPE_CHECK(x, Integer):
            #binary repr, reversed, and "["..."]" added
            x="["+x.binary()[::-1].replace(""," ")+"]"
        elif PY_TYPE_CHECK(x, Polynomial_GF2X):
            x = x.list() # this is slow but cimport leads to circular imports
        elif PY_TYPE_CHECK(x, FiniteField):
            if x.characteristic() == 2:
                x= list(x.modulus())
        elif PY_TYPE_CHECK(x, FiniteField_ext_pariElement):
            if x.parent().characteristic() == 2:
                x=x._pari_().centerlift().centerlift().subst('a',2).int_unsafe()
                x="0x"+hex(x)[2:][::-1]
        elif PY_TYPE_CHECK(x, FiniteField_givaroElement):
            x = "0x"+hex(x.integer_representation())[2:][::-1]
        elif PY_TYPE_CHECK(x, FiniteField_ntl_gf2eElement):
            x = x.polynomial().list()
        s = str(x).replace(","," ")
        sig_on()
        # TODO: this is very slow, but we wait until somebody complains
        GF2X_from_str(&self.x, s)
        sig_off()

    def __cinit__(self):
        GF2X_construct(&self.x)

    def __dealloc__(self):
        GF2X_destruct(&self.x)

    def __reduce__(self):
        """
        EXAMPLES:
            sage: f = ntl.GF2X(ntl.ZZ_pX([1,1,3],2))
            sage: loads(dumps(f)) == f
            True
            sage: f = ntl.GF2X('0x1c')
            sage: loads(dumps(f)) == f
            True
        """
        return unpickle_class_value, (ntl_GF2X, hex(self))

    def __repr__(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: ntl.GF2X(ntl.ZZ_pX([1,1,3],2)).__repr__()
            '[1 1 1]'
        """
        return GF2X_to_PyString(&self.x)

    def __mul__(ntl_GF2X self, other):
        """
        EXAMPLES:
            sage: f = ntl.GF2X([1,0,1,1]) ; g = ntl.GF2X([0,1])
            sage: f*g ## indirect doctest
            [0 1 0 1 1]
        """
        cdef ntl_GF2X r = PY_NEW(ntl_GF2X)
        if not PY_TYPE_CHECK(other, ntl_GF2X):
            other = ntl_GF2X(other)
        GF2X_mul(r.x, self.x, (<ntl_GF2X>other).x)
        return r

    def __div__(ntl_GF2X self, b):
        """
        EXAMPLES:
            sage: a = ntl.GF2X(4)
            sage: a / ntl.GF2X(2)
            [0 1]
            sage: a / ntl.GF2X(3)
            Traceback (most recent call last):
            ...
            ArithmeticError: self (=[0 0 1]) is not divisible by b (=[1 1])
        """
        cdef ntl_GF2X q = PY_NEW(ntl_GF2X)
        cdef int divisible

        if not PY_TYPE_CHECK(b, ntl_GF2X):
            b = ntl_GF2X(b)

        divisible = GF2X_divide(q.x, self.x, (<ntl_GF2X>b).x)
        if not divisible:
            raise ArithmeticError, "self (=%s) is not divisible by b (=%s)"%(self, b)
        return q

    def DivRem(ntl_GF2X self, b):
        """
        EXAMPLES:
            sage: a = ntl.GF2X(4)
            sage: a.DivRem( ntl.GF2X(2) )
            ([0 1], [])
            sage: a.DivRem( ntl.GF2X(3) )
            ([1 1], [1])
        """
        cdef ntl_GF2X q = PY_NEW(ntl_GF2X)
        cdef ntl_GF2X r = PY_NEW(ntl_GF2X)

        if not PY_TYPE_CHECK(b, ntl_GF2X):
            b = ntl_GF2X(b)

        GF2X_DivRem(q.x, r.x, self.x, (<ntl_GF2X>b).x)
        return q,r

    def __floordiv__(ntl_GF2X self, b):
        """
        EXAMPLES:
            sage: a = ntl.GF2X(4)
            sage: a // ntl.GF2X(2)
            [0 1]
            sage: a // ntl.GF2X(3)
            [1 1]
        """
        cdef ntl_GF2X q = PY_NEW(ntl_GF2X)

        if not PY_TYPE_CHECK(b, ntl_GF2X):
            b = ntl_GF2X(b)

        GF2X_div(q.x, self.x, (<ntl_GF2X>b).x)
        return q

    def __mod__(ntl_GF2X self, b):
        """
        EXAMPLES:
            sage: a = ntl.GF2X(4)
            sage: a % ntl.GF2X(2)
            []
            sage: a % ntl.GF2X(3)
            [1]
        """
        cdef ntl_GF2X r = PY_NEW(ntl_GF2X)

        if not PY_TYPE_CHECK(b, ntl_GF2X):
            b = ntl_GF2X(b)

        GF2X_rem(r.x, self.x, (<ntl_GF2X>b).x)
        return r

    def __sub__(ntl_GF2X self, other):
        """
        EXAMPLES:
            sage: f = ntl.GF2X([1,0,1,1]) ; g = ntl.GF2X([0,1])
            sage: f - g ## indirect doctest
            [1 1 1 1]
            sage: g - f
            [1 1 1 1]
        """
        cdef ntl_GF2X r = PY_NEW(ntl_GF2X)
        if not PY_TYPE_CHECK(other, ntl_GF2X):
            other = ntl_GF2X(other)
        GF2X_sub(r.x, self.x, (<ntl_GF2X>other).x)
        return r

    def __add__(ntl_GF2X self, other):
        """
        EXAMPLES:
            sage: f = ntl.GF2X([1,0,1,1]) ; g = ntl.GF2X([0,1,0])
            sage: f + g ## indirect doctest
            [1 1 1 1]
        """
        cdef ntl_GF2X r = PY_NEW(ntl_GF2X)
        if not PY_TYPE_CHECK(other, ntl_GF2X):
            other = ntl_GF2X(other)
        GF2X_add(r.x, self.x, (<ntl_GF2X>other).x)
        return r

    def __neg__(ntl_GF2X self):
        """
        EXAMPLES:
            sage: f = ntl.GF2X([1,0,1,1])
            sage: -f ## indirect doctest
            [1 0 1 1]
            sage: f == -f
            True
        """
        cdef ntl_GF2X r = PY_NEW(ntl_GF2X)
        GF2X_negate(r.x, self.x)
        return r

    def __pow__(ntl_GF2X self, long e, ignored):
        """
        EXAMPLES:
            sage: f = ntl.GF2X([1,0,1,1]) ; g = ntl.GF2X([0,1,0])
            sage: f**3 ## indirect doctest
            [1 0 1 1 1 0 0 1 1 1]
        """
        cdef ntl_GF2X r = PY_NEW(ntl_GF2X)
        GF2X_power(r.x, self.x, e)
        return r


    def __richcmp__(self, other, op):
        """
        EXAMPLES:
            sage: f = ntl.GF2X([1,0,1,1]) ; g = ntl.GF2X([0,1,0])
            sage: f == g ## indirect doctest
            False
            sage: f == f
            True
        """
        if op != 2 and op != 3:
            raise TypeError, "elements in GF(2)[X] are not ordered."

        if not PY_TYPE_CHECK(other, ntl_GF2X):
            other = ntl_GF2X(other)

        if not PY_TYPE_CHECK(self, ntl_GF2X):
            self = ntl_GF2X(self)

        cdef int t
        t = GF2X_equal((<ntl_GF2X>self).x, (<ntl_GF2X>other).x)
        if op == 2:
            return t == 1
        elif op == 3:
            return t == 0

    def __lshift__(ntl_GF2X self, int i):
        """
        Return left shift of self by i bits ( == multiplication by
        $X^i$).

        INPUT:
            i -- offset/power of X

        EXAMPLES:
            sage: a = ntl.GF2X(4); a
            [0 0 1]
            sage: a << 2
            [0 0 0 0 1]
        """
        cdef ntl_GF2X r = PY_NEW(ntl_GF2X)
        GF2X_LeftShift(r.x, self.x, <long>i)
        return r

    def __rshift__(ntl_GF2X self, int offset):
        """
        Return right shift of self by i bits ( == floor division by
        $X^i$).

        INPUT:
            i -- offset/power of X

        EXAMPLES:
            sage: a = ntl.GF2X(4); a
            [0 0 1]
            sage: a >> 1
            [0 1]
        """
        cdef ntl_GF2X r = PY_NEW(ntl_GF2X)
        GF2X_RightShift(r.x, self.x, <long>offset)
        return r

    def GCD(ntl_GF2X self, other):
        """
        Return GCD of self and other.

        INPUT:
            other -- ntl.GF2X

        EXAMPLE:
            sage: a = ntl.GF2X(10)
            sage: b = ntl.GF2X(4)
            sage: a.GCD(b)
            [0 1]
        """
        cdef ntl_GF2X gcd = PY_NEW(ntl_GF2X)

        if not PY_TYPE_CHECK(other, ntl_GF2X):
            other = ntl_GF2X(other)

        gcd.x = GF2X_GCD(self.x, (<ntl_GF2X>other).x)
        return gcd

    def XGCD(ntl_GF2X self, other):
        """
        Return the extended gcd of self and other, i.e., elements r, s, t such that

            r = s  * self + t  * other.

        INPUT:
            other -- ntl.GF2X

        EXAMPLE:
            sage: a = ntl.GF2X(10)
            sage: b = ntl.GF2X(4)
            sage: r,s,t = a.XGCD(b)
            sage: r == a*s + t*b
            True

        """
        cdef ntl_GF2X r = PY_NEW(ntl_GF2X)
        cdef ntl_GF2X s = PY_NEW(ntl_GF2X)
        cdef ntl_GF2X t = PY_NEW(ntl_GF2X)

        if not PY_TYPE_CHECK(other, ntl_GF2X):
            other = ntl_GF2X(other)

        GF2X_XGCD(r.x, s.x, t.x, self.x, (<ntl_GF2X>other).x)
        return r,s,t

    def deg(ntl_GF2X self):
        """
        Returns the degree of this polynomial

        EXAMPLES:
            sage: ntl.GF2X([1,0,1,1]).deg()
            3
        """
        return GF2X_deg(self.x)

    def list(ntl_GF2X self):
        """
        Represents this element as a list of binary digits.

        EXAMPLES:
             sage: e=ntl.GF2X([0,1,1])
             sage: e.list()
             [0, 1, 1]
             sage: e=ntl.GF2X('0xff')
             sage: e.list()
             [1, 1, 1, 1, 1, 1, 1, 1]

        OUTPUT:
             a list of digits representing the coefficients in this element's
             polynomial representation
        """
        return [self[i] for i in range(GF2X_deg(self.x)+1)]

    def bin(ntl_GF2X self):
        """
        Returns binary representation of this element. It is
        the same as setting \code{ntl.GF2XHexOutput(False)} and
        representing this element afterwards. However it should be
        faster and preserves the HexOutput state as opposed to
        the above code.

        EXAMPLES:
             sage: e=ntl.GF2X([1,1,0,1,1,1,0,0,1])
             sage: e.bin()
             '[1 1 0 1 1 1 0 0 1]'

        OUTPUT:
            string representing this element in binary digits
        """
        cdef long _hex = GF2XHexOutput_c[0]
        GF2XHexOutput_c[0] = 0
        s = GF2X_to_PyString(&self.x)
        GF2XHexOutput_c[0] = _hex
        return s

    def __hex__(ntl_GF2X self):
        """
        Returns hexadecimal representation of this element. It is
        the same as setting \code{ntl.GF2XHexOutput(True)} and
        representing this element afterwards. However it should be
        faster and preserves the HexOutput state as opposed to
        the above code.

        EXAMPLES:
             sage: e=ntl.GF2X([1,1,0,1,1,1,0,0,1])
             sage: hex(e)
             '0xb31'

        OUTPUT:
            string representing this element in hexadecimal

        """
        cdef long _hex = GF2XHexOutput_c[0]
        GF2XHexOutput_c[0] = 1
        s = GF2X_to_PyString(&self.x)
        GF2XHexOutput_c[0] = _hex
        return s

    def __hash__(self):
        return hash(hex(self))

    def _sage_(ntl_GF2X self, R=None):
        """
        Returns a Sage polynomial over GF(2) equivalent to
        this element. If a ring R is provided it is used
        to construct the polynomial in, otherwise
        an appropriate ring is generated.

        INPUT:
            self  -- GF2X element
            R     -- PolynomialRing over GF(2)

        OUTPUT:
            polynomial in R

        EXAMPLES:
            sage: f = ntl.GF2X([1,0,1,1,0,1])
            sage: f._sage_()
            x^5 + x^3 + x^2 + 1
            sage: f._sage_(PolynomialRing(Integers(2),'y'))
            y^5 + y^3 + y^2 + 1
        """
        if R is None:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            from sage.rings.finite_rings.constructor import FiniteField
            R = PolynomialRing(FiniteField(2), 'x')

        return R(map(int,self.list()))

    def coeff(self, int i):
        """
        Return the coefficient of the monomial $X^i$ in self.

        INPUT:
            i -- degree of X

        EXAMPLE:
            sage: e = ntl.GF2X([0,1,0,1])
            sage: e.coeff(0)
            0
            sage: e.coeff(1)
            1
        """
        cdef ntl_GF2 c = PY_NEW(ntl_GF2)
        c.x = GF2X_coeff(self.x, i)
        return c

    def __getitem__(self, int i):
        """
            sage: e = ntl.GF2X([0,1,0,1])
            sage: e[0] # indirect doctest
            0
            sage: e[1]
            1
        """
        cdef ntl_GF2 c = PY_NEW(ntl_GF2)
        c.x = GF2X_coeff(self.x, i)
        return c

    def LeadCoeff(self):
        """
        Return the leading coefficient of self. This is always 1
        except when self == 0.

        EXAMPLE:
            sage: e = ntl.GF2X([0,1])
            sage: e.LeadCoeff()
            1
            sage: e = ntl.GF2X(0)
            sage: e.LeadCoeff()
            0
        """
        cdef ntl_GF2 c = PY_NEW(ntl_GF2)
        c.x = GF2X_LeadCoeff(self.x)
        return c

    def ConstTerm(self):
        """
        Return the constant term of self.

        EXAMPLE:
            sage: e = ntl.GF2X([1,0,1])
            sage: e.ConstTerm()
            1
            sage: e = ntl.GF2X(0)
            sage: e.ConstTerm()
            0
        """
        cdef ntl_GF2 c = PY_NEW(ntl_GF2)
        c.x = GF2X_ConstTerm (self.x)
        return c

    def SetCoeff(self, int i, a):
        """
        Return the constant term of self.

        EXAMPLE:
            sage: e = ntl.GF2X([1,0,1]); e
            [1 0 1]
            sage: e.SetCoeff(1,1)
            sage: e
            [1 1 1]
        """
        cdef ntl_GF2 _a = ntl_GF2(a)

        GF2X_SetCoeff(self.x, i, _a.x)

    def __setitem__(self, int i, a):
        """
            sage: e = ntl.GF2X([1,0,1]); e
            [1 0 1]
            sage: e[1] = 1 # indirect doctest
            sage: e
            [1 1 1]
        """
        cdef ntl_GF2 _a = ntl_GF2(a)
        GF2X_SetCoeff(self.x, i, _a.x)

    def diff(self):
        """
        Differentiate self.
            sage: e = ntl.GF2X([1,0,1,1,0])
            sage: e.diff()
            [0 0 1]
        """
        cdef ntl_GF2X d = PY_NEW(ntl_GF2X)
        d.x = GF2X_diff(self.x)
        return d

    def reverse(self, int hi = -2):
        """
        Return reverse of a[0]..a[hi] (hi >= -1)
        hi defaults to deg(a)

        INPUT:
            hi -- bit position until which reverse is requested

        EXAMPLE:
            sage: e = ntl.GF2X([1,0,1,1,0])
            sage: e.reverse()
            [1 1 0 1]
        """
        cdef ntl_GF2X r = PY_NEW(ntl_GF2X)
        if hi < -1:
            hi = GF2X_deg(self.x)
        r.x = GF2X_reverse(self.x, hi)
        return r

    def weight(self):
        """
        Return the number of nonzero coefficients in self.

        EXAMPLE:
            sage: e = ntl.GF2X([1,0,1,1,0])
            sage: e.weight()
            3
        """
        return int(GF2X_weight(self.x))

    def __int__(self):
        """
            sage: e = ntl.GF2X([1,0,1,1,0])
            sage: int(e)
            Traceback (most recent call last):
            ...
            ValueError: cannot convert non-constant polynomial to integer
            sage: e = ntl.GF2X([1])
            sage: int(e)
            1
        """
        cdef long l = 0
        if GF2X_deg(self.x) != 0:
            raise ValueError, "cannot convert non-constant polynomial to integer"
        else:
            return GF2_conv_to_long(GF2X_coeff(self.x,0))

    def NumBits(self):
        """
        returns number of bits of self, i.e., deg(self) + 1.

        EXAMPLE:
            sage: e = ntl.GF2X([1,0,1,1,0])
            sage: e.NumBits()
            4
        """
        return int(GF2X_NumBits(self.x))

    def __len__(self):
        """
            sage: e = ntl.GF2X([1,0,1,1,0])
            sage: len(e)
            4
        """
        return int(GF2X_NumBits(self.x))

    def NumBytes(self):
        """
        Returns number of bytes of self, i.e., floor((NumBits(self)+7)/8)

        EXAMPLE:
            sage: e = ntl.GF2X([1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0,1,1])
            sage: e.NumBytes()
            3
        """
        return int(GF2X_NumBytes(self.x))
