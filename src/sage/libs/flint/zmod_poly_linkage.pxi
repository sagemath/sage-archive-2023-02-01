r"""
Linkage for arithmetic with FLINT's zmod_poly_t elements.

This file provides the backend for \class{Polynomial_zmod_flint} via
templating.

AUTHOR:
    -- Martin Albrecht (2009-01) another initial implementation
    -- Burcin Erocal (2008-11) initial implementation
"""
#*****************************************************************************
#       Copyright (C) 2008 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2009 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.flint.zmod_poly cimport *, zmod_poly_t
from sage.libs.flint.long_extras cimport *

include "../../ext/stdsage.pxi"

cdef inline celement *celement_new(unsigned long n):
    cdef celement *g = <celement *>sage_malloc(sizeof(zmod_poly_t))
    zmod_poly_init(g, n)
    return g

cdef inline int celement_delete(zmod_poly_t e, unsigned long n):
    zmod_poly_clear(e)
    sage_free(e)

cdef inline int celement_construct(zmod_poly_t e, unsigned long n):
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]

        sage: Q.<x> = GF(7)[]
    """
    zmod_poly_init(e, n)

cdef inline int celement_destruct(zmod_poly_t e, unsigned long n):
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: del x

        sage: Q.<x> = GF(7)[]
        sage: del x
    """
    zmod_poly_clear(e)

cdef inline int celement_gen(zmod_poly_t e, long i, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]

        sage: Q.<x> = GF(7)[]
    """
    zmod_poly_zero(e)
    zmod_poly_set_coeff_ui(e, 1, 1)

cdef object celement_repr(zmod_poly_t e, unsigned long n):
    raise NotImplementedError

cdef inline int celement_set(zmod_poly_t res, zmod_poly_t a, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: y = copy(x)
        sage: y is x
        False
        sage: y == x
        True

        sage: Q.<x> = GF(7)[]
        sage: y = copy(x)
        sage: y is x
        False
        sage: y == x
        True
    """
    zmod_poly_set(res, a)

cdef inline int celement_set_si(zmod_poly_t res, long i, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: P(32003)
        0
        sage: P(1)
        1
        sage: P(32004)
        1

        sage: Q.<x> = GF(7)[]
        sage: Q(7)
        0
        sage: Q(1)
        1
        sage: Q(8)
        1
    """
    while i < 0:
        i += n
    zmod_poly_zero(res)
    if i:
        zmod_poly_set_coeff_ui(res, 0, <unsigned long>i)

cdef inline long celement_get_si(zmod_poly_t res, unsigned long n) except -2:
    raise NotImplementedError

cdef inline bint celement_is_zero(zmod_poly_t a, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: P(1).is_zero()
        False
        sage: P(0).is_zero()
        True

        sage: Q.<x> = GF(7)[]
        sage: Q(1).is_zero()
        False
        sage: Q(0).is_zero()
        True
    """
    # is_zero doesn't exist
    return zmod_poly_degree(a) == -1

cdef inline bint celement_is_one(zmod_poly_t a, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: P(1).is_one()
        True
        sage: P(0).is_one()
        False

        sage: Q.<x> = GF(7)[]
        sage: Q(1).is_one()
        True
        sage: Q(0).is_one()
        False
    """

    return zmod_poly_is_one(a)

cdef inline bint celement_equal(zmod_poly_t a, zmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: (3*2)*x == 3*(2*x)
        True
        sage: (3*2)*x + 2 == 3*(2*x) + 1 + 1
        True
        sage: (3*2)*x + 7 == 3*(2*x) + 1 + 1
        False

        sage: Q.<x> = GF(7)[]
        sage: (3*2)*x == 3*(2*x)
        True
        sage: (3*2)*x + 2 == 3*(2*x) + 1 + 1
        True
        sage: (3*2)*x + 7 == 3*(2*x) + 1 + 1
        False
    """
    return zmod_poly_equal(a, b)

cdef inline int celement_cmp(zmod_poly_t l, zmod_poly_t r, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: x > x
        False
        sage: x^2 > x
        True
        sage: 3*x > x
        True

        sage: Q.<x> = GF(7)[]
        sage: x > x
        False
        sage: x^2 > x
        True
        sage: 3*x > x
        True

        sage: f = x^64 + x^20 + 1
        sage: g = x^63 + x^20 + 1
        sage: f > g
        True

        sage: f = x^64 + x^10 + 1
        sage: g = x^64 + x^20 + 1
        sage: f < g
        True

        sage: f = x^64 + x^20
        sage: g = x^64 + x^20 + 1
        sage: f < g
        True
    """
    cdef int deg_right = zmod_poly_degree(r)
    cdef int degdiff = deg_right - zmod_poly_degree(l)
    cdef int i
    cdef unsigned long rcoeff, lcoeff
    if degdiff > 0:
        return -1
    elif degdiff < 0:
        return 1
    else:
        if zmod_poly_equal(l, r):
            return 0
        i = deg_right
        rcoeff = zmod_poly_get_coeff_ui(r, i)
        lcoeff = zmod_poly_get_coeff_ui(l, i)
        while rcoeff == lcoeff and i > 0:
            i -= 1
            rcoeff = zmod_poly_get_coeff_ui(r, i)
            lcoeff = zmod_poly_get_coeff_ui(l, i)
        return cmp(lcoeff, rcoeff)

cdef long celement_len(zmod_poly_t a, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: (x + 1).degree()
        1
        sage: (x).degree()
        1
        sage: P(0).degree()
        -1

        sage: Q.<x> = GF(7)[]
        sage: (x + 1).degree()
        1
        sage: (x).degree()
        1
        sage: P(0).degree()
        -1
    """
    return <long>zmod_poly_length(a)

cdef inline int celement_add(zmod_poly_t res, zmod_poly_t a, zmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: x + 1
        x + 1

        sage: Q.<x> = GF(7)[]
        sage: x + 1
        x + 1
    """
    zmod_poly_add(res, a, b)

cdef inline int celement_sub(zmod_poly_t res, zmod_poly_t a, zmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: x - 1
        x + 32002

        sage: Q.<x> = GF(7)[]
        sage: x - 1
        x + 6
    """
    zmod_poly_sub(res, a, b)

cdef inline int celement_neg(zmod_poly_t res, zmod_poly_t a, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: -(x + 2)
        32002*x + 32001

        sage: Q.<x> = GF(7)[]
        sage: -(x + 2)
        6*x + 5
    """
    zmod_poly_neg(res, a)

cdef inline int celement_mul(zmod_poly_t res, zmod_poly_t a, zmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: (x + 1) * (x + 2)
        x^2 + 3*x + 2

        sage: Q.<x> = GF(7)[]
        sage: (x + 1) * (x + 2)
        x^2 + 3*x + 2
    """
    zmod_poly_mul(res, a, b)

cdef inline int celement_div(zmod_poly_t res, zmod_poly_t a, zmod_poly_t b, unsigned long n) except -2:
    raise NotImplementedError

cdef inline int celement_floordiv(zmod_poly_t res, zmod_poly_t a, zmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: (x + 1) // (x + 2)
        1
        sage: (3*x^2 + 1) // (x + 2)
        3*x + 31997
        sage: (x^2 + 3*x + 2)//(x + 1)
        x + 2
        sage: (x^2 + 3*x + 2)//(x + 2)
        x + 1

        sage: Q.<x> = GF(7)[]
        sage: (x + 1) // (x + 2)
        1
        sage: (3*x^2 + 1) // (x + 2)
        3*x + 1
        sage: (x^2 + 3*x + 2)//(x + 1)
        x + 2
        sage: (x^2 + 3*x + 2)//(x + 2)
        x + 1
    """
    zmod_poly_div(res, a, b)

cdef inline int celement_mod(zmod_poly_t res, zmod_poly_t a, zmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: f = 24998*x^2 + 29761*x + 2252
        sage: g = 20778*x^2 + 15346*x + 12697

        sage: f % g
        5815*x + 10280

        sage: f^5 % g
        7231*x + 17274

        sage: R.<x> = Integers(81)[]
        sage: f = x^7 + x + 1; g = x^3
        sage: r = f % g; r
        x + 1
        sage: g * x^4 + r
        x^7 + x + 1

        sage: f = x^3 + 1
        sage: f % 3
        Traceback (most recent call last):
        ...
        ValueError: Leading coefficient of a must be invertible.

        sage: f % 0
        Traceback (most recent call last):
        ...
        ZeroDivisionError
    """
    cdef zmod_poly_t q
    cdef unsigned long leadcoeff, modulus

    zmod_poly_init(q, n)
    leadcoeff = zmod_poly_get_coeff_ui(b, zmod_poly_degree(b))
    modulus = zmod_poly_modulus(b)
    if (leadcoeff > 1 and z_gcd(modulus,leadcoeff) != 1):
        raise ValueError("Leading coefficient of a must be invertible.")

    zmod_poly_divrem(q, res, a, b)
    zmod_poly_clear(q)

cdef inline int celement_quorem(zmod_poly_t q, zmod_poly_t r, zmod_poly_t a, zmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLES:
        sage: R.<x> = Integers(125)[]
        sage: f = x^5+1; g = (x+1)^2
        sage: q, r = f.quo_rem(g)
        sage: q
        x^3 + 123*x^2 + 3*x + 121
        sage: r
        5*x + 5
        sage: q*g + r
        x^5 + 1

        sage: x.quo_rem(5*x)
        Traceback (most recent call last):
        ...
        ValueError: Leading coefficient of a must be invertible.

        sage: x.quo_rem(0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError
    """
    cdef unsigned long leadcoeff, modulus

    leadcoeff = zmod_poly_get_coeff_ui(b, zmod_poly_degree(b))
    modulus = zmod_poly_modulus(b)
    if (leadcoeff > 1 and z_gcd(modulus,leadcoeff) != 1):
        raise ValueError("Leading coefficient of a must be invertible.")

    zmod_poly_divrem(q, r, a, b)

cdef inline int celement_inv(zmod_poly_t res, zmod_poly_t a, unsigned long n) except -2:
    raise NotImplementedError

cdef inline int celement_pow(zmod_poly_t res, zmod_poly_t x, long e, zmod_poly_t modulus, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: f = 24998*x^2 + 29761*x + 2252

        sage: f*f
        9426*x^4 + 15477*x^3 + 6531*x^2 + 14980*x + 15030
        sage: f^2
        9426*x^4 + 15477*x^3 + 6531*x^2 + 14980*x + 15030

        sage: f*f*f
        25062*x^6 + 30670*x^5 + 15816*x^4 + 20746*x^3 + 9142*x^2 + 5697*x + 20389
        sage: f^3
        25062*x^6 + 30670*x^5 + 15816*x^4 + 20746*x^3 + 9142*x^2 + 5697*x + 20389

        sage: f*f*f*f*f
        20269*x^10 + 20535*x^9 + 7313*x^8 + 7311*x^7 + 16853*x^6 + 142*x^5 + 23853*x^4 + 12065*x^3 + 516*x^2 + 8473*x + 17945
        sage: f^5
        20269*x^10 + 20535*x^9 + 7313*x^8 + 7311*x^7 + 16853*x^6 + 142*x^5 + 23853*x^4 + 12065*x^3 + 516*x^2 + 8473*x + 17945

        sage: f^0
        1

        sage: f^1
        24998*x^2 + 29761*x + 2252

        sage: f^-1
        1/(24998*x^2 + 29761*x + 2252)

        sage: f^-5
        1/(20269*x^10 + 20535*x^9 + 7313*x^8 + 7311*x^7 + 16853*x^6 + 142*x^5 + 23853*x^4 + 12065*x^3 + 516*x^2 + 8473*x + 17945)

     Testing the modulus:

        sage: g = 20778*x^2 + 15346*x + 12697

        sage: pow(f, 2, g)
        15328*x + 6968
        sage: f^2 % g
        15328*x + 6968

        sage: pow(f, -2, g)
        1/(15328*x + 6968)
        sage: (f^2 % g)^-1
        1/(15328*x + 6968)

        sage: pow(f, 5, g)
        7231*x + 17274
        sage: f^5 % g
        7231*x + 17274
    """
    cdef zmod_poly_t pow2
    cdef zmod_poly_t q
    cdef zmod_poly_t tmp

    zmod_poly_init(q, n)
    zmod_poly_init(tmp, n)

    if zmod_poly_degree(x) == 1 and zmod_poly_get_coeff_ui(x,0) == 0 and zmod_poly_get_coeff_ui(x,1) == 1:
        zmod_poly_zero(res)
        zmod_poly_set_coeff_ui(res,e,1)
    elif e == 0:
        zmod_poly_zero(res)
        zmod_poly_set_coeff_ui(res,0,1)
    elif e == 1:
        zmod_poly_set(res, x)
    elif e == 2:
        zmod_poly_sqr(res, x)
    else:
        if res == x:
            zmod_poly_set(tmp, x)
            x = tmp
        zmod_poly_init(pow2, n)
        zmod_poly_set(pow2, x)
        if e % 2:
            zmod_poly_set(res, x)
        else:
            zmod_poly_zero(res)
            zmod_poly_set_coeff_ui(res, 0, 1)
        e = e >> 1
        while(e != 0):
            zmod_poly_sqr(pow2, pow2)
            if e % 2:
                zmod_poly_mul(res, res, pow2)
            e = e >> 1
            if modulus != NULL:
                zmod_poly_divrem(q, res, res, modulus)
        zmod_poly_clear(pow2)

    if modulus != NULL:
        zmod_poly_divrem(q, res, res, modulus)
    zmod_poly_clear(q)
    zmod_poly_clear(tmp)

cdef inline int celement_gcd(zmod_poly_t res, zmod_poly_t a, zmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: f = P.random_element(degree=4); f
        16660*x^4 + 10640*x^3 + 1430*x^2 + 16460*x + 3566
        sage: g = P.random_element(degree=3); g
        28452*x^3 + 2561*x^2 + 22429*x + 5847
        sage: h = P.random_element(degree=2); h
        24731*x^2 + 28238*x + 18622
        sage: F = f*g; F
        13887*x^7 + 19164*x^6 + 25146*x^5 + 25986*x^4 + 21143*x^3 + 14830*x^2 + 14916*x + 16449
        sage: G = f*h; G
        11838*x^6 + 10154*x^5 + 15609*x^4 + 26164*x^3 + 11353*x^2 + 8656*x + 31830
        sage: d = (F).gcd(G); d
        x^4 + 18557*x^3 + 22917*x^2 + 30813*x + 4914
        sage: (F//d)*d == F
        True
        sage: (G//d)*d == G
        True

        sage: Q.<x> = GF(7)[]
        sage: f = Q.random_element(degree=4); f
        5*x^4 + 3*x^3 + 6*x^2 + 6*x + 1
        sage: g = Q.random_element(degree=3); g
        2*x^3 + 5*x^2 + 2*x + 3
        sage: h = Q.random_element(degree=2); h
        4*x^2 + 4*x + 6
        sage: F = f*g; F
        3*x^7 + 3*x^6 + 2*x^5 + 4*x^3 + 6*x + 3
        sage: G = f*h; G
        6*x^6 + 4*x^5 + 3*x^4 + 3*x^3 + x^2 + 5*x + 6
        sage: d = (F).gcd(G); d
        x^4 + 2*x^3 + 4*x^2 + 4*x + 3
        sage: (F//d)*d == F
        True
        sage: (G//d)*d == G
        True
    """
    if celement_is_zero(b, n):
        zmod_poly_set(res, a)
        return 0

    zmod_poly_gcd(res, a, b)
    cdef unsigned long leadcoeff = zmod_poly_get_coeff_ui(res, zmod_poly_degree(res))
    cdef unsigned long modulus = zmod_poly_modulus(res)
    if z_gcd(modulus,leadcoeff) == 1:
        zmod_poly_make_monic(res, res)

cdef inline int celement_xgcd(zmod_poly_t res, zmod_poly_t s, zmod_poly_t t, zmod_poly_t a, zmod_poly_t b, unsigned long n) except -2:
    """
    EXAMPLE:
        sage: P.<x> = GF(32003)[]
        sage: f = P.random_element(degree=4); f
        16660*x^4 + 10640*x^3 + 1430*x^2 + 16460*x + 3566
        sage: g = P.random_element(degree=3); g
        28452*x^3 + 2561*x^2 + 22429*x + 5847
        sage: h = P.random_element(degree=2); h
        24731*x^2 + 28238*x + 18622
        sage: F = f*g; F
        13887*x^7 + 19164*x^6 + 25146*x^5 + 25986*x^4 + 21143*x^3 + 14830*x^2 + 14916*x + 16449
        sage: G = f*h; G
        11838*x^6 + 10154*x^5 + 15609*x^4 + 26164*x^3 + 11353*x^2 + 8656*x + 31830
        sage: d,s,t = (F).xgcd(G); d
        x^4 + 18557*x^3 + 22917*x^2 + 30813*x + 4914
        sage: (F//d)*d == F
        True
        sage: (G//d)*d == G
        True

        sage: Q.<x> = GF(7)[]
        sage: f = Q.random_element(degree=4); f
        5*x^4 + 3*x^3 + 6*x^2 + 6*x + 1
        sage: g = Q.random_element(degree=3); g
        2*x^3 + 5*x^2 + 2*x + 3
        sage: h = Q.random_element(degree=2); h
        4*x^2 + 4*x + 6
        sage: F = f*g; F
        3*x^7 + 3*x^6 + 2*x^5 + 4*x^3 + 6*x + 3
        sage: G = f*h; G
        6*x^6 + 4*x^5 + 3*x^4 + 3*x^3 + x^2 + 5*x + 6
        sage: d,s,t = (F).xgcd(G); d
        x^4 + 2*x^3 + 4*x^2 + 4*x + 3
        sage: (F//d)*d == F
        True
        sage: (G//d)*d == G
        True
    """
    zmod_poly_xgcd(res, s, t, a, b)


cdef inline unsigned long zmod_poly_evaluate_horner(zmod_poly_t _poly, unsigned long c):
    """
    Evaluate _poly at c using Horner's rule.

    AUTHOR:
        -- Bill Hart (2009-01)

    EXAMPLE:
        sage: P.<x> = PolynomialRing(GF(7))
        sage: f= x^2 + 1
        sage: f(0)
        1
        sage: f(2)
        5
        sage: f(3)
        3

    NOTE:
        This function should be removed as soon as we update the the
        new FLINT version which has this function natively.
    """
    cdef zmod_poly_struct *poly = <zmod_poly_struct *>_poly
    cdef long i
    cdef unsigned long bits
    if poly.length == 0:
        return 0

    cdef unsigned long p = poly.p
    cdef double pinv = poly.p_inv

    cdef unsigned long value = poly.coeffs[poly.length - 1]

    if FLINT_BITS == 64:
        bits = FLINT_BIT_COUNT(p)
        if bits < FLINT_D_BITS:
            for i in xrange(poly.length - 2, -1, -1):
                value = z_addmod(z_mulmod_precomp(value, c, p, pinv), poly.coeffs[i], p)
        else:
            for i in xrange(poly.length - 2, -1, -1):
                value = z_addmod(z_mulmod2_precomp(value, c, p, pinv), poly.coeffs[i], p)
    else:
        for i in xrange(poly.length - 2, -1, -1):
            value = z_addmod(z_mulmod2_precomp(value, c, p, pinv), poly.coeffs[i], p)
    return value
