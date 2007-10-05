#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include 'misc.pxi'
include 'decl.pxi'

from ntl_ZZ import unpickle_class_value

##############################################################################
#
# ntl_GF2X: Polynomials over GF(2) via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2006-01: initial version (based on code by William Stein)
#
##############################################################################

cdef bint __have_GF2X_hex_repr = False # hex representation of GF2X


cdef class ntl_GF2X:
    """
    Polynomials over GF(2) via NTL
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
        """
        from sage.rings.finite_field_element import FiniteField_ext_pariElement
        from sage.rings.finite_field import FiniteField_ext_pari
        from sage.rings.finite_field_givaro import FiniteField_givaro,FiniteField_givaroElement
        from sage.rings.polynomial.polynomial_modn_dense_ntl import Polynomial_dense_mod_p
        from sage.rings.integer import Integer

        if isinstance(x, Integer):
            #binary repr, reversed, and "["..."]" added
            x="["+x.binary()[::-1].replace(""," ")+"]"
        elif type(x) == int:
            #hex repr, "0x" stripped, reversed (!)
            x="0x"+hex(x)[2:][::-1]
        elif isinstance(x, Polynomial_dense_mod_p):
            if x.base_ring().characteristic():
                x=x._Polynomial_dense_mod_n__poly
        elif isinstance(x, (FiniteField_ext_pari,FiniteField_givaro)):
            if x.characteristic() == 2:
                x= list(x.modulus())
        elif isinstance(x, FiniteField_ext_pariElement):
            if x.parent().characteristic() == 2:
                x=x._pari_().centerlift().centerlift().subst('a',2).int_unsafe()
                x="0x"+hex(x)[2:][::-1]
        elif isinstance(x, FiniteField_givaroElement):
            x = "0x"+hex(int(x))[2:][::-1]
        s = str(x).replace(","," ")
        _sig_on
        GF2X_from_str(&self.gf2x_x, s)
        _sig_off

    def __new__(self, x=[]):
        GF2X_construct(&self.gf2x_x)

    def __dealloc__(self):
        GF2X_destruct(&self.gf2x_x)

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
        return unpickle_class_value, (ntl_GF2X, self.hex())

    def __repr__(self):
        return GF2X_to_PyString(&self.gf2x_x)

    def __mul__(ntl_GF2X self, other):
        cdef ntl_GF2X y
        cdef ntl_GF2X r = ntl_GF2X()
        if not isinstance(other, ntl_GF2X):
            other = ntl_GF2X(other)
        y = other
        _sig_on
        GF2X_mul(r.gf2x_x, self.gf2x_x, y.gf2x_x)
        _sig_off
        return r

    def __sub__(ntl_GF2X self, other):
        cdef ntl_GF2X y
        cdef ntl_GF2X r = ntl_GF2X()
        if not isinstance(other, ntl_GF2X):
            other = ntl_GF2X(other)
        y = other
        _sig_on
        GF2X_sub(r.gf2x_x, self.gf2x_x, y.gf2x_x)
        _sig_off
        return r

    def __add__(ntl_GF2X self, other):
        cdef ntl_GF2X y
        cdef ntl_GF2X r = ntl_GF2X()
        if not isinstance(other, ntl_GF2X):
            other = ntl_GF2X(other)
        y = other
        _sig_on
        GF2X_add(r.gf2x_x, self.gf2x_x, y.gf2x_x)
        _sig_off
        return r

    def __neg__(ntl_GF2X self):
        cdef ntl_GF2X r = ntl_GF2X()
        _sig_on
        GF2X_negate(r.gf2x_x, self.gf2x_x)
        _sig_off
        return r

    def __pow__(ntl_GF2X self, long e, ignored):
        cdef ntl_GF2X y
        cdef ntl_GF2X r
        if not isinstance(other, ntl_GF2X):
            other = ntl_GF2X(other)
        y = other
        r = ntl_GF2X()
        _sig_on
        GF2X_power(r.gf2x_x, self.gf2x_x, e)
        _sig_off
        return r


    def __cmp__(ntl_GF2X self, ntl_GF2X other):
        cdef int t
        _sig_on
        t = GF2X_equal(self.gf2x_x, other.gf2x_x)
        _sig_off
        if t:
            return 0
        return 1

    def degree(ntl_GF2X self):
        """
        Returns the degree of this polynomial
        """
        return GF2X_deg(self.gf2x_x)

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
        #yields e.g. "[1 1 0 0 1 1 0 1]"
        #_sig_on
        s = GF2X_to_bin(&self.gf2x_x)

        #yields e.g. [1,1,0,0,1,1,0,1]
        return map(int,list(s[1:][:len(s)-2].replace(" ","")))


    def bin(ntl_GF2X self):
        """
        Returns binary representation of this element. It is
        the same as setting \code{ntl.hex_output(False)} and
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
        return GF2X_to_bin(&self.gf2x_x)

    def hex(ntl_GF2X self):
        """
        Returns hexadecimal representation of this element. It is
        the same as setting \code{ntl.hex_output(True)} and
        representing this element afterwards. However it should be
        faster and preserves the HexOutput state as opposed to
        the above code.

        EXAMPLES:
             sage: e=ntl.GF2X([1,1,0,1,1,1,0,0,1])
             sage: e.hex()
             '0xb31'

        OUTPUT:
            string representing this element in hexadecimal

        """
        return GF2X_to_hex(&self.gf2x_x)

    def _sage_(ntl_GF2X self,R=None,cache=None):
        """
        Returns a SAGE polynomial over GF(2) equivalent to
        this element. If a ring R is provided it is used
        to construct the polynomial in, otherwise
        an appropriate ring is generated.

        INPUT:
            self  -- GF2X element
            R     -- PolynomialRing over GF(2)
            cache -- optional NTL to SAGE cache (dict)

        OUTPUT:
            polynomial in R
        """
        if R==None:
            from sage.rings.polynomial.polynomial_ring import PolynomialRing
            from sage.rings.finite_field import FiniteField
            R = PolynomialRing(FiniteField(2), 'x')

        if cache != None:
            try:
                return cache[self.hex()]
            except KeyError:
                cache[self.hex()] = R(self.list())
                return cache[self.hex()]

        return R(self.list())
