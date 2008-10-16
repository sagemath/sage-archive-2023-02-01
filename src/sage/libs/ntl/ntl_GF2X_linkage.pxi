r"""
Linkage for arithmetic with NTL's GF2X elements.

This file provides the backend for \class{Polynomial_GF2X} via
templating.

AUTHOR:
    -- Martin Albrecht (2008-10): initial version
"""
#*****************************************************************************
#       Copyright (C) 2008 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.ntl.ntl_GF2_decl cimport *, GF2_c
from sage.libs.ntl.ntl_GF2X_decl cimport *, GF2X_c, GF2XModulus_c

cdef int celement_construct(GF2X_c *e, parent):
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
    """
    GF2X_construct(e)

cdef int celement_destruct(GF2X_c *e, parent):
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: del x
    """
    GF2X_destruct(e)

cdef int celement_gen(GF2X_c *e, long i, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
    """
    cdef unsigned char g = 2
    GF2XFromBytes(e[0], <unsigned char *>(&g), 1)

cdef object celement_repr(GF2X_c *e, parent):
    """
    We ignore NTL's printing.

    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: x
        x
    """
    #return GF2X_to_PyString(e)
    raise NotImplementedError

cdef inline int celement_set(GF2X_c* res, GF2X_c* a, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: y = x; y
        x
    """
    res[0] = a[0]

cdef inline int celement_set_si(GF2X_c* res, long i, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: P(0)
        0
        sage: P(2)
        0
        sage: P(1)
        1
    """
    GF2X_conv_long(res[0], i)

cdef inline long celement_get_si(GF2X_c* res, parent) except -2:
    raise NotImplementedError

cdef inline bint celement_is_zero(GF2X_c* a, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: bool(x), x.is_zero()
        (True, False)
        sage: bool(P(0)), P(0).is_zero()
        (False, True)
    """
    return GF2X_IsZero(a[0])

cdef inline bint celement_is_one(GF2X_c *a, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: x.is_one()
        False
        sage: P(1).is_one()
        True
    """
    return GF2X_IsOne(a[0])

cdef inline bint celement_equal(GF2X_c *a, GF2X_c *b, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: x == x
        True
        sage: y = x; x == y
        True
        sage: x^2 + 1 == x^2 + x
        False
    """
    return GF2X_equal(a[0], b[0])

cdef inline int celement_cmp(GF2X_c *a, GF2X_c *b, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: x != 1
        True
        sage: x < 1
        False
        sage: x > 1
        True
    """
    cdef bint t
    cdef long diff
    diff = GF2X_NumBits(a[0]) - GF2X_NumBits(b[0])
    if diff > 0:
        return 1
    elif diff == 0:
        t = GF2X_equal(a[0], b[0])
        return not t
    else:
        return -1

cdef inline long celement_hash(GF2X_c *a, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: {x:1}
        {x: 1}
    """
    cdef long _hex = GF2XHexOutput_c[0]
    GF2XHexOutput_c[0] = 1
    s = GF2X_to_PyString(a)
    GF2XHexOutput_c[0] = _hex
    return hash(s)

cdef long celement_len(GF2X_c *a, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: len(x)
        2
        sage: len(x+1)
        2
    """
    return int(GF2X_NumBits(a[0]))

cdef inline int celement_add(GF2X_c *res, GF2X_c *a, GF2X_c *b, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: x + 1
        x + 1
    """
    GF2X_add(res[0], a[0], b[0])

cdef inline int celement_sub(GF2X_c* res, GF2X_c* a, GF2X_c* b, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: x - 1
        x + 1
    """
    GF2X_sub(res[0], a[0], b[0])

cdef inline int celement_neg(GF2X_c* res, GF2X_c* a, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: -x
        x
    """
    res[0] = a[0]

cdef inline int celement_mul(GF2X_c* res, GF2X_c* a, GF2X_c* b, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: x*(x+1)
        x^2 + x
    """
    GF2X_mul(res[0], a[0], b[0])

cdef inline int celement_div(GF2X_c* res, GF2X_c* a, GF2X_c* b, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
    """
    return GF2X_divide(res[0], a[0], b[0])

cdef inline int celement_floordiv(GF2X_c* res, GF2X_c* a, GF2X_c* b, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: x//(x + 1)
        1
        sage: (x + 1)//x
        1
    """
    GF2X_div(res[0], a[0], b[0])

cdef inline int celement_mod(GF2X_c* res, GF2X_c* a, GF2X_c* b, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: (x^2 + 1) % x^2
        1
    """
    GF2X_rem(res[0], a[0], b[0])

cdef inline int celement_quorem(GF2X_c* q, GF2X_c* r, GF2X_c* a, GF2X_c* b, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: f = x^2 + x + 1
        sage: f.quo_rem(x + 1)
        (x, 1)
    """
    GF2X_DivRem(q[0], r[0], a[0], b[0])

cdef inline int celement_inv(GF2X_c* res, GF2X_c* a, parent) except -2:
    """
    We ignore NTL here and use the fraction field constructor.

    EXAMPLE:
        sage: P.<x> = GF(2)[]
    """
    raise NotImplementedError

cdef inline int celement_pow(GF2X_c* res, GF2X_c* x, long e, GF2X_c *modulus, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: x^1000
        x^1000
        sage: (x+1)^2
        x^2 + 1
        sage: (x+1)^(-2)
        1/(x^2 + 1)
        sage: f = x^9 + x^7 + x^6 + x^5 + x^4 + x^2 + x
        sage: h = x^10 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + 1
        sage: (f^2) % h
        x^9 + x^8 + x^7 + x^5 + x^3
        sage: pow(f, 2, h)
        x^9 + x^8 + x^7 + x^5 + x^3
    """
    cdef GF2XModulus_c mod

    if modulus == NULL:
        if GF2X_IsX(x[0]):
                GF2X_LeftShift(res[0], x[0], e - 1)
        else:
            do_sig = GF2X_deg(x[0]) > 1e5
            if do_sig: _sig_on
            GF2X_power(res[0], x[0], e)
            if do_sig: _sig_off
    else:
        GF2XModulus_build(mod, modulus[0])

        do_sig = GF2X_deg(x[0]) > 1e5
        if do_sig: _sig_on
        GF2X_PowerMod_long_pre(res[0], x[0], e, mod)
        if do_sig: _sig_off

cdef inline int celement_gcd(GF2X_c* res, GF2X_c* a, GF2X_c *b, parent) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(2)[]
        sage: f = x*(x+1)
        sage: f.gcd(x+1)
        x + 1
        sage: f.gcd(x^2)
        x
    """
    GF2X_GCD(res[0], a[0], b[0])
