"""
This module provides the core implementation of multiprecision
floating-point arithmetic. Operations are done in-place.
"""

include '../../ext/interrupt.pxi'
include "../../ext/stdsage.pxi"
include "../../ext/python_int.pxi"
include "../../ext/python_long.pxi"
include "../../ext/python_float.pxi"
include "../../ext/python_complex.pxi"
include "../../ext/python_number.pxi"

cdef extern from "math.h":
    cdef double fpow "pow" (double, double)
    cdef double frexp "frexp" (double, int*)

from sage.libs.gmp.all cimport *
from sage.rings.integer cimport Integer

cdef extern from "mpz_pylong.h":
    cdef mpz_get_pylong(mpz_t src)
    cdef mpz_get_pyintlong(mpz_t src)
    cdef int mpz_set_pylong(mpz_t dst, src) except -1
    cdef long mpz_pythonhash(mpz_t src)

cdef mpz_set_integer(mpz_t v, x):
    if PyInt_CheckExact(x):
        mpz_set_si(v, PyInt_AS_LONG(x))
    elif PyLong_CheckExact(x):
        mpz_set_pylong(v, x)
    elif PY_TYPE_CHECK(x, Integer):
        mpz_set(v, (<Integer>x).value)
    else:
        raise TypeError("cannot convert %s to an integer" % x)

cdef inline void mpz_add_si(mpz_t a, mpz_t b, long x):
    if x >= 0:
        mpz_add_ui(a, b, x)
    else:
        # careful: overflow when negating INT_MIN
        mpz_sub_ui(a, b, <unsigned long>(-x))

cdef inline mpzi(mpz_t n):
    return mpz_get_pyintlong(n)

# This should be done better
cdef int mpz_tstbit_abs(mpz_t z, unsigned long bit_index):
    cdef int res
    if mpz_sgn(z) < 0:
        mpz_neg(z, z)
        res = mpz_tstbit(z, bit_index)
        mpz_neg(z, z)
    else:
        res = mpz_tstbit(z, bit_index)
    return res

cdef unsigned long mpz_bitcount(mpz_t z):
    if mpz_sgn(z) == 0:
        return 0
    return mpz_sizeinbase(z, 2)

# The following limits allowed exponent shifts. We could use mpz_fits_slong_p,
# but then (-LONG_MIN) wraps around; we may also not be able to add large
# shifts safely. A higher limit could be used on 64-bit systems, but
# it is unlikely that anyone will run into this (adding numbers
# that differ by 2^(2^30), at precisions of 2^30 bits).

DEF MAX_SHIFT = 1073741824      # 2^30

cdef int mpz_reasonable_shift(mpz_t z):
    if mpz_sgn(z) > 0:
        return mpz_cmp_ui(z, MAX_SHIFT) < 0
    else:
        return mpz_cmp_si(z, -MAX_SHIFT) > 0

DEF ROUND_N = 0
DEF ROUND_F = 1
DEF ROUND_C = 2
DEF ROUND_D = 3
DEF ROUND_U = 4

DEF S_NORMAL = 0
DEF S_ZERO = 1
DEF S_NZERO = 2
DEF S_INF = 3
DEF S_NINF = 4
DEF S_NAN = 5

cdef inline str rndmode_to_python(int rnd):
    if rnd == ROUND_N: return 'n'
    if rnd == ROUND_F: return 'f'
    if rnd == ROUND_C: return 'c'
    if rnd == ROUND_D: return 'd'
    if rnd == ROUND_U: return 'u'

cdef inline rndmode_from_python(str rnd):
    if rnd == 'n': return ROUND_N
    if rnd == 'f': return ROUND_F
    if rnd == 'c': return ROUND_C
    if rnd == 'd': return ROUND_D
    if rnd == 'u': return ROUND_U

cdef MPopts opts_exact
cdef MPopts opts_double_precision
cdef MPopts opts_mini_prec

opts_exact.prec = 0
opts_exact.rounding = ROUND_N
opts_double_precision.prec = 53
opts_double_precision.rounding = ROUND_N
opts_mini_prec.prec = 5
opts_mini_prec.rounding = ROUND_D

cdef double _double_inf = float("1e300") * float("1e300")
cdef double _double_ninf = -_double_inf
cdef double _double_nan = _double_inf - _double_inf

cdef inline void MPF_init(MPF *x):
    """Allocate space and set value to zero.
    Must be called exactly once when creating a new MPF."""
    x.special = S_ZERO
    mpz_init(x.man)
    mpz_init(x.exp)

cdef inline void MPF_clear(MPF *x):
    """Deallocate space. Must be called exactly once when finished with an MPF."""
    mpz_clear(x.man)
    mpz_clear(x.exp)

cdef inline void MPF_set(MPF *dest, MPF *src):
    """Clone MPF value. Assumes source value is already normalized."""
    if src is dest:
        return
    dest.special = src.special
    mpz_set(dest.man, src.man)
    mpz_set(dest.exp, src.exp)

cdef inline void MPF_set_zero(MPF *x):
    """Set value to 0."""
    x.special = S_ZERO

cdef inline void MPF_set_nan(MPF *x):
    """Set value to NaN (not a number)."""
    x.special = S_NAN

cdef inline void MPF_set_inf(MPF *x):
    """Set value to +infinity."""
    x.special = S_INF

cdef inline void MPF_set_ninf(MPF *x):
    """Set value to -infinity."""
    x.special = S_NINF

cdef MPF_set_si(MPF *x, long n):
    """Set value to that of a given C (long) integer."""
    if n:
        x.special = S_NORMAL
        mpz_set_si(x.man, n)
        mpz_set_ui(x.exp, 0)
        MPF_normalize(x, opts_exact)
    else:
        MPF_set_zero(x)

cdef MPF_set_int(MPF *x, n):
    """Set value to that of a given Python integer."""
    x.special = S_NORMAL
    mpz_set_integer(x.man, n)
    if mpz_sgn(x.man):
        mpz_set_ui(x.exp, 0)
        MPF_normalize(x, opts_exact)
    else:
        MPF_set_zero(x)

cdef MPF_set_man_exp(MPF *x, man, exp):
    """
    Set value to man*2^exp where man, exp may be of any appropriate
    Python integer types.
    """
    x.special = S_NORMAL
    mpz_set_integer(x.man, man)
    mpz_set_integer(x.exp, exp)
    MPF_normalize(x, opts_exact)


# Temporary variables. Note: not thread-safe.
# Used by MPF_add/MPF_sub/MPF_div
cdef mpz_t tmp_exponent
mpz_init(tmp_exponent)
cdef MPF tmp0
MPF_init(&tmp0)

# Used by MPF_hypot and MPF_cmp, which may call MPF_add/MPF_sub
cdef MPF tmp1
MPF_init(&tmp1)
cdef MPF tmp2
MPF_init(&tmp2)


# Constants needed in a few places
cdef MPF MPF_C_1
MPF_init(&MPF_C_1)
MPF_set_si(&MPF_C_1, 1)
cdef Integer MPZ_ZERO = Integer(0)
cdef tuple _mpf_fzero = (0, MPZ_ZERO, 0, 0)
cdef tuple _mpf_fnan = (0, MPZ_ZERO, -123, -1)
cdef tuple _mpf_finf = (0, MPZ_ZERO, -456, -2)
cdef tuple _mpf_fninf = (1, MPZ_ZERO, -789, -3)

cdef MPF_set_tuple(MPF *x, tuple value):
    """
    Set value of an MPF to that of a normalized (sign, man, exp, bc) tuple
    in the format used by mpmath.libmp.
    """
    #cdef int sign
    cdef Integer man
    sign, _man, exp, bc = value
    if PY_TYPE_CHECK(_man, Integer):
        man = <Integer>_man
    else:
        # This is actually very unlikely; it should never happen
        # in internal code that man isn't an Integer. Maybe the check
        # can be avoided by doing checks in e.g. MPF_set_any?
        man = Integer(_man)
    if mpz_sgn(man.value):
        MPF_set_man_exp(x, man, exp)
        if sign:
            mpz_neg(x.man, x.man)
        return
    if value == _mpf_fzero:
        MPF_set_zero(x)
    elif value == _mpf_finf:
        MPF_set_inf(x)
    elif value == _mpf_fninf:
        MPF_set_ninf(x)
    else:
        MPF_set_nan(x)

cdef MPF_to_tuple(MPF *x):
    """Convert MPF value to (sign, man, exp, bc) tuple."""
    cdef Integer man
    if x.special:
        if x.special == S_ZERO: return _mpf_fzero
        #if x.special == S_NZERO: return _mpf_fnzero
        if x.special == S_INF: return _mpf_finf
        if x.special == S_NINF: return _mpf_fninf
        return _mpf_fnan
    man = PY_NEW(Integer)
    if mpz_sgn(x.man) < 0:
        mpz_neg(man.value, x.man)
        sign = 1
    else:
        mpz_set(man.value, x.man)
        sign = 0
    exp = mpz_get_pyintlong(x.exp)
    bc = mpz_sizeinbase(x.man, 2)
    return (sign, man, exp, bc)

cdef MPF_set_double(MPF *r, double x):
    """
    Set value of a C double.
    """
    cdef int exp
    cdef double man
    if x != x:
        MPF_set_nan(r)
        return
    if x == _double_inf:
        MPF_set_inf(r)
        return
    if x == _double_ninf:
        MPF_set_ninf(r)
        return
    man = frexp(x, &exp)
    man *= 9007199254740992.0
    mpz_set_d(r.man, man)
    mpz_set_si(r.exp, exp-53)
    r.special = S_NORMAL
    MPF_normalize(r, opts_exact)

import math as pymath

# TODO: implement this function safely without using the Python math module
cdef double MPF_to_double(MPF *x, bint strict):
    """Convert MPF value to a Python float."""
    if x.special == S_NORMAL:
        man = mpzi(x.man)
        exp = mpzi(x.exp)
        bc = mpz_sizeinbase(x.man, 2)
        try:
            if bc < 100:
                return pymath.ldexp(man, exp)
            # Try resizing the mantissa. Overflow may still happen here.
            n = bc - 53
            m = man >> n
            return pymath.ldexp(m, exp + n)
        except OverflowError:
            if strict:
                raise
            # Overflow to infinity
            if exp + bc > 0:
                if man < 0:
                    return _double_ninf
                else:
                    return _double_inf
            # Underflow to zero
            return 0.0
    if x.special == S_ZERO:
        return 0.0
    if x.special == S_INF:
        return _double_inf
    if x.special == S_NINF:
        return _double_ninf
    return _double_nan

cdef MPF_to_fixed(mpz_t r, MPF *x, long prec, bint truncate):
    """
    Set r = x, r being in the format of a fixed-point number with prec bits.
    Floor division is used unless truncate=True in which case
    truncating division is used.
    """
    cdef long shift
    if x.special:
        if x.special == S_ZERO or x.special == S_NZERO:
            mpz_set_ui(r, 0)
            return
        raise ValueError("cannot create fixed-point number from special value")
    if mpz_reasonable_shift(x.exp):
        # XXX: signed integer overflow
        shift = mpz_get_si(x.exp) + prec
        if shift >= 0:
            mpz_mul_2exp(r, x.man, shift)
        else:
            if truncate:
                mpz_tdiv_q_2exp(r, x.man, -shift)
            else:
                mpz_fdiv_q_2exp(r, x.man, -shift)
        return
    # Underflow
    if mpz_sgn(x.exp) < 0:
        mpz_set_ui(r, 0)
        return
    raise OverflowError("cannot convert huge number to fixed-point format")

cdef int MPF_sgn(MPF *x):
    """
    Gives the sign of an MPF (-1, 0, or 1).
    """
    if x.special:
        if x.special == S_INF:
            return 1
        if x.special == S_NINF:
            return -1
        return 0
    return mpz_sgn(x.man)

cdef void MPF_neg(MPF *r, MPF *s):
    """
    Sets r = -s. MPF_neg(x, x) negates in place.
    """
    if s.special:
        if   s.special == S_ZERO: r.special = S_ZERO #r.special = S_NZERO
        elif s.special == S_NZERO: r.special = S_ZERO
        elif s.special == S_INF: r.special = S_NINF
        elif s.special == S_NINF: r.special = S_INF
        else: r.special = s.special
        return
    r.special = s.special
    mpz_neg(r.man, s.man)
    if r is not s:
        mpz_set(r.exp, s.exp)

cdef void MPF_abs(MPF *r, MPF *s):
    if s.special:
        if    s.special == S_NINF: r.special = S_INF
        else: r.special = s.special
        return
    r.special = s.special
    mpz_abs(r.man, s.man)
    if r is not s:
        mpz_set(r.exp, s.exp)

cdef MPF_normalize(MPF *x, MPopts opts):
    """
    Normalize.

    With prec = 0, trailing zero bits are stripped but no rounding
    is performed.
    """
    cdef int sign
    cdef long trail, bc, shift
    if x.special != S_NORMAL:
        return
    sign = mpz_sgn(x.man)
    if sign == 0:
        x.special = S_ZERO
        mpz_set_ui(x.exp, 0)
        return
    bc = mpz_sizeinbase(x.man, 2)
    shift = bc - opts.prec
    # Ok if mantissa small and no trailing zero bits
    if (shift <= 0 or not opts.prec) and mpz_odd_p(x.man):
        return
    # Mantissa is too large, so divide by appropriate power of 2
    # Need to be careful about rounding
    if shift > 0 and opts.prec:
        if opts.rounding == ROUND_N:
            if mpz_tstbit_abs(&x.man, shift-1):
                if mpz_tstbit_abs(&x.man, shift) or mpz_scan1(x.man, 0) < (shift-1):
                    if sign < 0:
                        mpz_fdiv_q_2exp(x.man, x.man, shift)
                    else:
                        mpz_cdiv_q_2exp(x.man, x.man, shift)
                else:
                    mpz_tdiv_q_2exp(x.man, x.man, shift)
            else:
                mpz_tdiv_q_2exp(x.man, x.man, shift)
        elif opts.rounding == ROUND_D:
            mpz_tdiv_q_2exp(x.man, x.man, shift)
        elif opts.rounding == ROUND_F:
            mpz_fdiv_q_2exp(x.man, x.man, shift)
        elif opts.rounding == ROUND_C:
            mpz_cdiv_q_2exp(x.man, x.man, shift)
        elif opts.rounding == ROUND_U:
            if sign < 0:
                mpz_fdiv_q_2exp(x.man, x.man, shift)
            else:
                mpz_cdiv_q_2exp(x.man, x.man, shift)
        else:
            raise ValueError("bad rounding mode")
    else:
        shift = 0
    # Strip trailing bits
    trail = mpz_scan1(x.man, 0)
    if 0 < trail < bc:
        mpz_tdiv_q_2exp(x.man, x.man, trail)
        shift += trail
    mpz_add_si(x.exp, x.exp, shift)

cdef void _add_special(MPF *r, MPF *s, MPF *t):
    if s.special == S_ZERO:
        # (+0) + (-0) = +0
        if t.special == S_NZERO:
            MPF_set(r, s)
        # (+0) + x = x
        else:
            MPF_set(r, t)
    elif t.special == S_ZERO:
        # (-0) + (+0) = +0
        if s.special == S_NZERO:
            MPF_set(r, t)
        # x + (+0) = x
        else:
            MPF_set(r, s)
    # (+/- 0) + x = x
    elif s.special == S_NZERO:
        MPF_set(r, t)
    elif t.special == S_NZERO:
        MPF_set(r, s)
    # (+/- inf) + (-/+ inf) = nan
    elif ((s.special == S_INF and t.special == S_NINF) or
       (s.special == S_NINF and t.special == S_INF)):
        MPF_set_nan(r)
    # nan or +/- inf trumps any finite number
    elif s.special == S_NAN or t.special == S_NAN:
        MPF_set_nan(r)
    elif s.special:
        MPF_set(r, s)
    else:
        MPF_set(r, t)
    return

cdef void _sub_special(MPF *r, MPF *s, MPF *t):
    if s.special == S_ZERO:
        # (+0) - (+/-0) = (+0)
        if t.special == S_NZERO:
            MPF_set(r, s)
        else:
            # (+0) - x = (-x)
            MPF_neg(r, t)
    elif t.special == S_ZERO:
        # x - (+0) = x; also covers (-0) - (+0) = (-0)
        MPF_set(r, s)
    # (-0) - x = x
    elif s.special == S_NZERO:
        # (-0) - (-0) = (+0)
        if t.special == S_NZERO:
            MPF_set_zero(r)
        # (-0) - x = -x
        else:
            MPF_neg(r, t)
    elif t.special == S_NZERO:
        # x - (-0) = x
        MPF_set(r, s)
    # (+/- inf) - (+/- inf) = nan
    elif ((s.special == S_INF and t.special == S_INF) or
       (s.special == S_NINF and t.special == S_NINF)):
        MPF_set_nan(r)
    elif s.special == S_NAN or t.special == S_NAN:
        MPF_set_nan(r)
    # nan - x or (+/-inf) - x = l.h.s
    elif s.special:
        MPF_set(r, s)
    # x - nan or x - (+/-inf) = (- r.h.s)
    else:
        MPF_neg(r, t)

cdef void _mul_special(MPF *r, MPF *s, MPF *t):
    if s.special == S_ZERO:
        if t.special == S_NORMAL or t.special == S_ZERO:
            MPF_set(r, s)
        elif t.special == S_NZERO:
            MPF_set(r, t)
        else:
            MPF_set_nan(r)
    elif t.special == S_ZERO:
        if s.special == S_NORMAL:
            MPF_set(r, t)
        elif s.special == S_NZERO:
            MPF_set(r, s)
        else:
            MPF_set_nan(r)
    elif s.special == S_NZERO:
        if t.special == S_NORMAL:
            if mpz_sgn(t.man) < 0:
                MPF_set_zero(r)
            else:
                MPF_set(r, s)
        else:
            MPF_set_nan(r)
    elif t.special == S_NZERO:
        if s.special == S_NORMAL:
            if mpz_sgn(s.man) < 0:
                MPF_set_zero(r)
            else:
                MPF_set(r, t)
        else:
            MPF_set_nan(r)
    elif s.special == S_NAN or t.special == S_NAN:
        MPF_set_nan(r)
    else:
        if MPF_sgn(s) == MPF_sgn(t):
            MPF_set_inf(r)
        else:
            MPF_set_ninf(r)

cdef _div_special(MPF *r, MPF *s, MPF *t):
    # TODO: handle signed zeros correctly
    if s.special == S_NAN or t.special == S_NAN:
        MPF_set_nan(r)
    elif t.special == S_ZERO or t.special == S_NZERO:
        raise ZeroDivisionError
    elif s.special == S_ZERO or s.special == S_NZERO:
        MPF_set_zero(r)
    elif s.special == S_NORMAL:
        MPF_set_zero(r)
    elif s.special == S_INF or s.special == S_NINF:
        if t.special == S_INF or t.special == S_NINF:
            MPF_set_nan(r)
        elif MPF_sgn(s) == MPF_sgn(t):
            MPF_set_inf(r)
        else:
            MPF_set_ninf(r)
    # else:
    elif t.special == S_INF or t.special == S_NINF:
        MPF_set_zero(r)

cdef _add_perturbation(MPF *r, MPF *s, int sign, MPopts opts):
    cdef long shift
    if opts.rounding == ROUND_N:
        MPF_set(r, s)
    else:
        shift = opts.prec - mpz_sizeinbase(s, 2) + 8
        if shift < 0:
            shift = 8
        mpz_mul_2exp(r.man, s.man, shift)
        mpz_add_si(r.man, r.man, sign)
        mpz_sub_ui(r.exp, s.exp, shift)
        MPF_normalize(r, opts)

cdef MPF_add(MPF *r, MPF *s, MPF *t, MPopts opts):
    """
    Set r = s + t, with exact rounding.

    With prec = 0, the addition is performed exactly. Note that this
    may cause overflow if the exponents are huge.
    """
    cdef long shift, sbc, tbc
    #assert (r is not s) and (r is not t)
    if s.special or t.special:
        _add_special(r, s, t)
        return
    r.special = S_NORMAL
    # Difference between exponents
    mpz_sub(tmp_exponent, s.exp, t.exp)
    if mpz_reasonable_shift(tmp_exponent):
        shift = mpz_get_si(tmp_exponent)
        if shift >= 0:
            # |s| >> |t|
            if shift > 2*opts.prec and opts.prec:
                sbc = mpz_sizeinbase(s,2)
                tbc = mpz_sizeinbase(t,2)
                if shift + sbc - tbc > opts.prec+8:
                    _add_perturbation(r, s, mpz_sgn(t.man), opts)
                    return
            # |s| > |t|
            mpz_mul_2exp(tmp0.man, s.man, shift)
            mpz_add(r.man, tmp0.man, t.man)
            mpz_set(r.exp, t.exp)
            MPF_normalize(r, opts)
        elif shift < 0:
            shift = -shift
            # |s| << |t|
            if shift > 2*opts.prec and opts.prec:
                sbc = mpz_sizeinbase(s,2)
                tbc = mpz_sizeinbase(t,2)
                if shift + tbc - sbc > opts.prec+8:
                    _add_perturbation(r, t, mpz_sgn(s.man), opts)
                    return
            # |s| < |t|
            mpz_mul_2exp(tmp0.man, t.man, shift)
            mpz_add(r.man, tmp0.man, s.man)
            mpz_set(r.exp, s.exp)
            MPF_normalize(r, opts)
    else:
        if not opts.prec:
            raise OverflowError("the exact result does not fit in memory")
        # |s| >>> |t|
        if mpz_sgn(tmp_exponent) > 0:
            _add_perturbation(r, s, mpz_sgn(t.man), opts)
        # |s| <<< |t|
        else:
            _add_perturbation(r, t, mpz_sgn(s.man), opts)

cdef MPF_sub(MPF *r, MPF *s, MPF *t, MPopts opts):
    """
    Set r = s - t, with exact rounding.

    With prec = 0, the addition is performed exactly. Note that this
    may cause overflow if the exponents are huge.
    """
    cdef long shift, sbc, tbc
    #assert (r is not s) and (r is not t)
    if s.special or t.special:
        _sub_special(r, s, t)
        return
    r.special = S_NORMAL
    # Difference between exponents
    mpz_sub(tmp_exponent, s.exp, t.exp)
    if mpz_reasonable_shift(tmp_exponent):
        shift = mpz_get_si(tmp_exponent)
        if shift >= 0:
            # |s| >> |t|
            if shift > 2*opts.prec and opts.prec:
                sbc = mpz_sizeinbase(s,2)
                tbc = mpz_sizeinbase(t,2)
                if shift + sbc - tbc > opts.prec+8:
                    _add_perturbation(r, s, -mpz_sgn(t.man), opts)
                    return
            # |s| > |t|
            mpz_mul_2exp(tmp0.man, s.man, shift)
            mpz_sub(r.man, tmp0.man, t.man)
            mpz_set(r.exp, t.exp)
            MPF_normalize(r, opts)
        elif shift < 0:
            shift = -shift
            # |s| << |t|
            if shift > 2*opts.prec and opts.prec:
                sbc = mpz_sizeinbase(s,2)
                tbc = mpz_sizeinbase(t,2)
                if shift + tbc - sbc > opts.prec+8:
                    _add_perturbation(r, t, -mpz_sgn(s.man), opts)
                    MPF_neg(r, r)
                    return
            # |s| < |t|
            mpz_mul_2exp(tmp0.man, t.man, shift)
            mpz_sub(r.man, s.man, tmp0.man)
            mpz_set(r.exp, s.exp)
            MPF_normalize(r, opts)
    else:
        if not opts.prec:
            raise OverflowError("the exact result does not fit in memory")
        # |s| >>> |t|
        if mpz_sgn(tmp_exponent) > 0:
            _add_perturbation(r, s, -mpz_sgn(t.man), opts)
        # |s| <<< |t|
        else:
            _add_perturbation(r, t, -mpz_sgn(s.man), opts)
            MPF_neg(r, r)

cdef bint MPF_eq(MPF *s, MPF *t):
    if s.special == S_NAN or t.special == S_NAN:
        return False
    if s.special == t.special:
        if s.special == S_NORMAL:
            return (mpz_cmp(s.man, t.man) == 0) and (mpz_cmp(s.exp, t.exp) == 0)
        else:
            return True
    return False

cdef bint MPF_ne(MPF *s, MPF *t):
    if s.special == S_NAN or t.special == S_NAN:
        return True
    if s.special == S_NORMAL and t.special == S_NORMAL:
        return (mpz_cmp(s.man, t.man) != 0) or (mpz_cmp(s.exp, t.exp) != 0)
    return s.special != t.special

cdef int MPF_cmp(MPF *s, MPF *t):
    cdef long sbc, tbc
    cdef int cm
    if MPF_eq(s, t):
        return 0
    if s.special != S_NORMAL or t.special != S_NORMAL:
        if s.special == S_ZERO: return -MPF_sgn(t)
        if t.special == S_ZERO: return MPF_sgn(s)
        if t.special == S_NAN: return 1
        if s.special == S_INF: return 1
        if t.special == S_NINF: return 1
        return -1
    if mpz_sgn(s.man) != mpz_sgn(t.man):
        if mpz_sgn(s.man) < 0:
            return -1
        else:
            return 1
    if not mpz_cmp(s.exp, t.exp):
        return mpz_cmp(s.man, t.man)
    mpz_add_ui(tmp1.exp, s.exp, mpz_sizeinbase(s.man, 2))
    mpz_add_ui(tmp2.exp, t.exp, mpz_sizeinbase(t.man, 2))
    cm = mpz_cmp(tmp1.exp, tmp2.exp)
    if mpz_sgn(s.man) < 0:
        if cm < 0: return 1
        if cm > 0: return -1
    else:
        if cm < 0: return -1
        if cm > 0: return 1
    MPF_sub(&tmp1, s, t, opts_mini_prec)
    return MPF_sgn(&tmp1)

cdef bint MPF_lt(MPF *s, MPF *t):
    if s.special == S_NAN or t.special == S_NAN:
        return False
    return MPF_cmp(s, t) < 0

cdef bint MPF_le(MPF *s, MPF *t):
    if s.special == S_NAN or t.special == S_NAN:
        return False
    return MPF_cmp(s, t) <= 0

cdef bint MPF_gt(MPF *s, MPF *t):
    if s.special == S_NAN or t.special == S_NAN:
        return False
    return MPF_cmp(s, t) > 0

cdef bint MPF_ge(MPF *s, MPF *t):
    if s.special == S_NAN or t.special == S_NAN:
        return False
    return MPF_cmp(s, t) >= 0

cdef MPF_mul(MPF *r, MPF *s, MPF *t, MPopts opts):
    """
    Set r = s * t, with correct rounding.

    With prec = 0, the multiplication is performed exactly,
    i.e. no rounding is performed.
    """
    if s.special or t.special:
        _mul_special(r, s, t)
    else:
        r.special = S_NORMAL
        mpz_mul(r.man, s.man, t.man)
        mpz_add(r.exp, s.exp, t.exp)
        if opts.prec:
            MPF_normalize(r, opts)

cdef MPF_div(MPF *r, MPF *s, MPF *t, MPopts opts):
    """
    Set r = s / t, with correct rounding.
    """
    cdef int sign
    cdef long sbc, tbc, extra
    cdef mpz_t rem
    #assert (r is not s) and (r is not t)
    if s.special or t.special:
        _div_special(r, s, t)
        return
    r.special = S_NORMAL
    # Division by a power of two <=> shift exponents
    if mpz_cmp_si(t.man, 1) == 0:
        MPF_set(&tmp0, s)
        mpz_sub(tmp0.exp, tmp0.exp, t.exp)
        MPF_normalize(&tmp0, opts)
        MPF_set(r, &tmp0)
        return
    elif mpz_cmp_si(t.man, -1) == 0:
        MPF_neg(&tmp0, s)
        mpz_sub(tmp0.exp, tmp0.exp, t.exp)
        MPF_normalize(&tmp0, opts)
        MPF_set(r, &tmp0)
        return
    sign = mpz_sgn(s.man) != mpz_sgn(t.man)
    # Same strategy as for addition: if there is a remainder, perturb
    # the result a few bits outside the precision range before rounding
    extra = opts.prec - mpz_sizeinbase(s.man,2) + mpz_sizeinbase(t.man,2) + 5
    if extra < 5:
        extra = 5
    mpz_init(rem)
    mpz_mul_2exp(tmp0.man, s.man, extra)
    mpz_tdiv_qr(r.man, rem, tmp0.man, t.man)
    if mpz_sgn(rem):
        mpz_mul_2exp(r.man, r.man, 1)
        if sign:
            mpz_sub_ui(r.man, r.man, 1)
        else:
            mpz_add_ui(r.man, r.man, 1)
        extra  += 1
    mpz_clear(rem)
    mpz_sub(r.exp, s.exp, t.exp)
    mpz_sub_ui(r.exp, r.exp, extra)
    MPF_normalize(r, opts)

cdef MPF_sqrt(MPF *r, MPF *s, MPopts opts):
    """
    Set r = sqrt(s), with correct rounding.
    """
    cdef long shift
    cdef mpz_t rem
    #assert r is not s
    if s.special:
        if s.special == S_ZERO or s.special == S_INF:
            MPF_set(r, s)
        else:
            MPF_set_nan(r)
        return
    if mpz_sgn(s.man) < 0:
        MPF_set_nan(r)
        return
    r.special = S_NORMAL
    if mpz_odd_p(s.exp):
        mpz_sub_ui(r.exp, s.exp, 1)
        mpz_mul_2exp(r.man, s.man, 1)
    elif mpz_cmp_ui(s.man, 1) == 0:
        # Square of a power of two
        mpz_set_ui(r.man, 1)
        mpz_tdiv_q_2exp(r.exp, s.exp, 1)
        MPF_normalize(r, opts)
        return
    else:
        mpz_set(r.man, s.man)
        mpz_set(r.exp, s.exp)
    shift = 2*opts.prec - mpz_sizeinbase(r.man,2) + 4
    if shift < 4:
        shift = 4
    shift += shift & 1
    mpz_mul_2exp(r.man, r.man, shift)
    if opts.rounding == ROUND_F or opts.rounding == ROUND_D:
        mpz_sqrt(r.man, r.man)
    else:
        mpz_init(rem)
        mpz_sqrtrem(r.man, rem, r.man)
        if mpz_sgn(rem):
            mpz_mul_2exp(r.man, r.man, 1)
            mpz_add_ui(r.man, r.man, 1)
            shift += 2
        mpz_clear(rem)
    mpz_add_si(r.exp, r.exp, -shift)
    mpz_tdiv_q_2exp(r.exp, r.exp, 1)
    MPF_normalize(r, opts)

cdef MPF_hypot(MPF *r, MPF *a, MPF *b, MPopts opts):
    """
    Set r = sqrt(a^2 + b^2)
    """
    cdef MPopts tmp_opts
    if a.special == S_ZERO:
        MPF_abs(r, b)
        MPF_normalize(r, opts)
        return
    if b.special == S_ZERO:
        MPF_abs(r, a)
        MPF_normalize(r, opts)
        return
    tmp_opts = opts
    tmp_opts.prec += 30
    MPF_mul(&tmp1, a, a, opts_exact)
    MPF_mul(&tmp2, b, b, opts_exact)
    MPF_add(r, &tmp1, &tmp2, tmp_opts)
    MPF_sqrt(r, r, opts)

cdef MPF_pow_int(MPF *r, MPF *x, mpz_t n, MPopts opts):
    """
    Set r = x ** n. Currently falls back to mpmath.libmp
    unless n is tiny.
    """
    cdef long m, absm
    cdef unsigned long bc
    cdef int nsign
    if x.special != S_NORMAL:
        nsign = mpz_sgn(n)
        if x.special == S_ZERO:
            if nsign < 0:
                raise ZeroDivisionError
            elif nsign == 0:
                MPF_set(r, &MPF_C_1)
            else:
                MPF_set_zero(r)
        elif x.special == S_INF:
            if nsign > 0:
                MPF_set(r, x)
            elif nsign == 0:
                MPF_set_nan(r)
            else:
                MPF_set_zero(r)
        elif x.special == S_NINF:
            if nsign > 0:
                if mpz_odd_p(n):
                    MPF_set(r, x)
                else:
                    MPF_neg(r, x)
            elif nsign == 0:
                MPF_set_nan(r)
            else:
                MPF_set_zero(r)
        else:
            MPF_set_nan(r)
        return
    bc = mpz_sizeinbase(r.man,2)
    r.special = S_NORMAL
    if mpz_reasonable_shift(n):
        m = mpz_get_si(n)
        if m == 0:
            MPF_set(r, &MPF_C_1)
            return
        if m == 1:
            MPF_set(r, x)
            MPF_normalize(r, opts)
            return
        if m == 2:
            MPF_mul(r, x, x, opts)
            return
        if m == -1:
            MPF_div(r, &MPF_C_1, x, opts)
            return
        if m == -2:
            MPF_mul(r, x, x, opts_exact)
            MPF_div(r, &MPF_C_1, r, opts)
            return
        absm = abs(m)
        if bc * absm < 10000:
            mpz_pow_ui(r.man, x.man, absm)
            mpz_mul_ui(r.exp, x.exp, absm)
            if m < 0:
                MPF_div(r, &MPF_C_1, r, opts)
            else:
                MPF_normalize(r, opts)
            return
    r.special = S_NORMAL
    # (2^p)^n
    if mpz_cmp_si(x.man, 1) == 0:
        mpz_set(r.man, x.man)
        mpz_mul(r.exp, x.exp, n)
        return
    # (-2^p)^n
    if mpz_cmp_si(x.man, -1) == 0:
        if mpz_odd_p(n):
            mpz_set(r.man, x.man)
        else:
            mpz_neg(r.man, x.man)
        mpz_mul(r.exp, x.exp, n)
        return
    # TODO: implement efficiently here
    import mpmath.libmp
    MPF_set_tuple(r,
        mpmath.libmp.mpf_pow_int(MPF_to_tuple(x), mpzi(n),
        opts.prec, rndmode_to_python(opts.rounding)))

