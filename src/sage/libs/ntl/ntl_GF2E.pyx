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

from ntl_GF2X cimport ntl_GF2X

cdef make_GF2X(GF2X_c *x):
    """ These make_XXXX functions are deprecated and should be phased out."""
    cdef ntl_GF2X y
    _sig_off
    y = ntl_GF2X()
    y.gf2x_x = x[0]
    GF2X_delete(x)
    return y

##############################################################################
#
# ntl_GF2E: GF(2**x) via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2006-01: initial version (based on cody by William Stein)
#
##############################################################################

def GF2X_hex_repr(have_hex=None):
    """
    Represent GF2X and GF2E elements in the more compact
    hexadecimal form to the user.

    If no parameter is provided the currently set value will be
    returned.

    INPUT:
        have_hex if True hex representation will be used
    """
    global __have_GF2X_hex_repr

    if have_hex==None:
        return __have_GF2X_hex_repr

    if have_hex==True:
        GF2X_hex(1)
    else:
        GF2X_hex(0)
    __have_GF2X_hex_repr=have_hex

def ntl_GF2E_modulus(p=None):
    """
    Initializes the current modulus to P; required: deg(P) >= 1

    The input is either ntl.GF2X or is tried to be converted to a
    ntl.GF2X element.

    If no parameter p is given: Yields copy of the current GF2E
    modulus.

    INPUT:
        p -- modulus

    EXAMPLES:
        sage: ntl.GF2E_modulus([1,1,0,1,1,0,0,0,1])
        sage: ntl.GF2E_modulus().hex()
        '0xb11'
    """
    global __have_GF2E_modulus
    cdef ntl_GF2X elem

    if p==None:
        if __have_GF2E_modulus:
            return make_GF2X(<GF2X_c*>GF2E_modulus())
        else:
            raise "NoModulus"

    if not isinstance(p,ntl_GF2X):
        elem = ntl_GF2X(p)
    else:
        elem = p

    if(elem.degree()<1):
        raise "DegreeToSmall"

    ntl_GF2E_set_modulus(<GF2X_c*>&elem.gf2x_x)
    __have_GF2E_modulus = True

def ntl_GF2E_modulus_degree():
    """
    Returns deg(modulus) for GF2E elements
    """
    if __have_GF2E_modulus:
        return GF2E_degree()
    else:
        raise "NoModulus"

def ntl_GF2E_sage(names='a'):
    """
    Returns a SAGE FiniteField element matching the current modulus.

    EXAMPLES:
        sage: ntl.GF2E_modulus([1,1,0,1,1,0,0,0,1])
        sage: ntl.GF2E_sage()
        Finite Field in a of size 2^8
    """
    from sage.rings.finite_field import FiniteField
    f = ntl_GF2E_modulus()._sage_()
    return FiniteField(int(2)**GF2E_degree(),modulus=f,name=names)


def ntl_GF2E_random():
    """
    Returns a random element from GF2E modulo the current modulus.
    """
    cdef ntl_GF2E r = ntl_GF2E()
    _sig_on
    r.gf2e_x = GF2E_random()
    _sig_off
    return r

def unpickle_GF2E(hex, mod):
    ntl_GF2E_modulus(mod)
    return ntl_GF2E(hex)


# make sure not to segfault
__have_GF2E_modulus = False


cdef class ntl_GF2E(ntl_GF2X):
    r"""
    The \\class{GF2E} represents a finite extension field over GF(2) using NTL.
    Elements are represented as polynomials over GF(2) modulo \\code{ntl.GF2E_modulus()}.

    This modulus must be set using \\code{ ntl.GF2E_modulus(p) } and is unique for
    alle elements in ntl.GF2E. So it is not possible at the moment e.g. to have elements
    in GF(2**4) and GF(2**8) at the same time. You might however be lucky and get away
    with not touch the elements in GF(2**4) while having the GF(2**8) modulus set and vice
    versa.
    """

    def __init__(self,x=[]):
        """
        Constructs a new  finite field element in GF(2**x).

        A modulus \emph{must} have been set with \code{ntl.GF2E_modulus(ntl.GF2X(<value>))}
        when calling this constructor.  A value may be passed to this constructor. If you pass a string
        to the constructor please note that byte sequences and the hexadecimal notation are Little Endian in NTL.
        So e.g. '[0 1]' == '0x2' == x.

        INPUT:
            x -- value to be assigned to this element. Same types as ntl.GF2X() are accepted.

        OUTPUT:
            a new ntl.GF2E element

        EXAMPLES:
            sage: m=ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: ntl.GF2E(ntl.ZZ_pX([1,1,3],2))
            [1 1 1]
            sage: ntl.GF2E('0x1c')
            [1 0 0 0 0 0 1 1]
            sage: ntl.GF2E('[1 0 1 0]')
            [1 0 1]
            sage: ntl.GF2E([1,0,1,0])
            [1 0 1]
            sage: ntl.GF2E(GF(2**8,'a').gen()**20)
            [0 0 1 0 1 1 0 1]
        """
        if not __have_GF2E_modulus:
            raise "NoModulus"

        s = str(ntl_GF2X(x))
        _sig_on
        GF2E_from_str(&self.gf2e_x, s)
        self.gf2x_x = GF2E_rep(self.gf2e_x)
        _sig_off

    def __new__(self, x=[]):
        GF2E_construct(&self.gf2e_x)

    def __dealloc__(self):
        GF2E_destruct(&self.gf2e_x)

    def __reduce__(self):
        """
        EXAMPES:
            sage: m=ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: a = ntl.GF2E(ntl.ZZ_pX([1,1,3],2))
            sage: m=ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,1,0,1]))
            sage: loads(dumps(a)) == a
            True
        """
        return unpickle_GF2E, (self.ntl_GF2X().hex(), ntl_GF2E_modulus())

    def __repr__(self):
        _sig_on
        return string_delete(GF2E_to_str(&self.gf2e_x))

    def __mul__(ntl_GF2E self, other):
        cdef ntl_GF2E y
        cdef ntl_GF2E r = ntl_GF2E()
        if not isinstance(other, ntl_GF2E):
            other = ntl_GF2E(other)
        y = other
        _sig_on
        GF2E_mul(r.gf2e_x, self.gf2e_x, y.gf2e_x)
        _sig_off
        return r

    def __sub__(ntl_GF2E self, other):
        cdef ntl_GF2E y
        cdef ntl_GF2E r = ntl_GF2E()
        if not isinstance(other, ntl_GF2E):
            other = ntl_GF2E(other)
        y = other
        _sig_on
        GF2E_sub(r.gf2e_x, self.gf2e_x, y.gf2e_x)
        _sig_off
        return r

    def __add__(ntl_GF2E self, other):
        cdef ntl_GF2E y
        cdef ntl_GF2E r = ntl_GF2E()
        if not isinstance(other, ntl_GF2E):
            other = ntl_GF2E(other)
        y = other
        _sig_on
        GF2E_add(r.gf2e_x, self.gf2e_x, y.gf2e_x)
        _sig_off
        return r

    def __neg__(ntl_GF2E self):
        cdef ntl_GF2E r = ntl_GF2E()
        _sig_on
        GF2E_negate(r.gf2e_x, self.gf2e_x)
        _sig_off
        return r

    def __pow__(ntl_GF2E self, long e, ignored):
        cdef ntl_GF2E y
        cdef ntl_GF2E r
        if not isinstance(other, ntl_GF2E):
            other = ntl_GF2E(other)
        y = other
        r = ntl_GF2E()
        _sig_on
        GF2E_power(r.gf2e_x, self.gf2e_x, e)
        _sig_off
        return r

    def __cmp__(ntl_GF2E self, ntl_GF2E other):
        cdef int t
        _sig_on
        t = GF2E_eq(&self.gf2e_x, &other.gf2e_x)
        _sig_off
        if t:
            return 0
        return 1

    def is_zero(ntl_GF2E self):
        """
        Returns True if this element equals zero, False otherwise.
        """
        return bool(GF2E_is_zero(&self.gf2e_x))

    def is_one(ntl_GF2E self):
        """
        Returns True if this element equals one, False otherwise.
        """
        return bool(GF2E_is_one(&self.gf2e_x))

    def __copy__(self):
        cdef ntl_GF2E r = ntl_GF2E()
        _sig_on
        r.gf2e_x = self.gf2e_x
        _sig_off
        return r

    def copy(ntl_GF2E self):
        """
        Returns a copy of this element.
        """
        cdef ntl_GF2E r = ntl_GF2E()
        _sig_on
        r.gf2e_x = self.gf2e_x
        _sig_off
        return r

    def trace(ntl_GF2E self):
        """
        Returns the trace of this element.
        """
        return GF2E_trace(&self.gf2e_x)

    def ntl_GF2X(ntl_GF2E self):
        """
        Returns a ntl.GF2X copy of this element.

        EXAMPLE:
            sage: m=ntl.GF2E_modulus(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: a = ntl.GF2E('0x1c')
            sage: a.ntl_GF2X()
            [1 0 0 0 0 0 1 1]
            sage: a.copy().ntl_GF2X()
            [1 0 0 0 0 0 1 1]
        """
        cdef ntl_GF2X x = ntl_GF2X()
        x.gf2x_x = self.gf2x_x
        return x

    def _sage_(ntl_GF2E self, k=None, cache=None):
        """
        Returns a \class{FiniteFieldElement} representation
        of this element. If a \class{FiniteField} k is provided
        it is constructed in this field if possible. A \class{FiniteField}
        will be constructed if none is provided.

        INPUT:
            self  -- \class{GF2E} element
            k     -- optional GF(2**deg)
            cache -- optional NTL to SAGE conversion dictionary

        OUTPUT:
            FiniteFieldElement over k

        EXAMPLE:
            sage: ntl.GF2E_modulus([1,1,0,1,1,0,0,0,1])
            sage: e = ntl.GF2E([0,1])
            sage: a = e._sage_(); a
            a

        """
        cdef int i
        cdef int length
        deg= GF2E_degree()

        if k==None:
            from sage.rings.finite_field import FiniteField
            f = ntl_GF2E_modulus()._sage_()
            k = FiniteField(2**deg,name='a',modulus=f)

        if cache != None:
            try:
                return cache[self.hex()]
            except KeyError:
                pass

        a=k.gen()
        l = self.list()

        length = len(l)
        ret = 0

        for i from 0 <= i < length:
            if l[i]==1:
                ret = ret + a**i

        if cache != None:
            cache[self.hex()] = ret

        return ret
