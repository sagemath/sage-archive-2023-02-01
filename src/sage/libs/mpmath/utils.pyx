# Utilities for Sage-mpmath interaction
# Also patches some mpmath functions for speed

include "../../ext/stdsage.pxi"

from sage.rings.integer cimport Integer
from sage.rings.real_mpfr cimport RealNumber
from sage.rings.complex_number cimport ComplexNumber
from sage.structure.element cimport Element

import sage.all

from sage.libs.mpfr cimport *
from sage.libs.gmp.all cimport *

from sage.rings.complex_field import ComplexField
from sage.rings.real_mpfr import RealField_constructor as RealField

cpdef from_man_exp(man, exp, long prec = 0, str rnd = 'd'):
    """
    Create normalized mpf value tuple from mantissa and exponent.
    With prec > 0, rounds the result in the desired direction
    if necessary.

    EXAMPLES::

        sage: from mpmath.libmpf import from_man_exp
        sage: from_man_exp(-6, -1)
        (1, 3, 0, 2)
        sage: from_man_exp(-6, -1, 1, 'd')
        (1, 1, 1, 1)
        sage: from_man_exp(-6, -1, 1, 'u')
        (1, 1, 2, 1)
    """
    cdef Integer res
    cdef long bc
    res = Integer(man)
    bc = mpz_sizeinbase(res.value, 2)
    if not prec:
        prec = bc
    if mpz_sgn(res.value) < 0:
        mpz_neg(res.value, res.value)
        return normalize(1, res, exp, bc, prec, rnd)
    else:
        return normalize(0, res, exp, bc, prec, rnd)

cpdef normalize(long sign, Integer man, exp, long bc, long prec, str rnd):
    """
    Create normalized mpf value tuple from full list of components.

    EXAMPLES::

        sage: from mpmath.libmpf import normalize
        sage: normalize(0, 4, 5, 3, 53, 'n')
        (0, 1, 7, 1)
    """
    cdef long shift
    cdef unsigned long haverem
    cdef Integer res
    cdef unsigned long trail
    if mpz_sgn(man.value) == 0:
        from mpmath.libmpf import fzero
        return fzero
    if bc <= prec and mpz_odd_p(man.value):
        return (sign, man, exp, bc)
    shift = bc - prec
    res = PY_NEW(Integer)
    if shift > 0:
        if rnd == 'n':
            if mpz_tstbit(man.value, shift-1) and (mpz_tstbit(man.value, shift)\
                or (mpz_scan1(man.value, 0) < (shift-1))):
                mpz_cdiv_q_2exp(res.value, man.value, shift)
            else:
                mpz_fdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'd':
            mpz_fdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'f':
            if sign: mpz_cdiv_q_2exp(res.value, man.value, shift)
            else:    mpz_fdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'c':
            if sign: mpz_fdiv_q_2exp(res.value, man.value, shift)
            else:    mpz_cdiv_q_2exp(res.value, man.value, shift)
        elif rnd == 'u':
            mpz_cdiv_q_2exp(res.value, man.value, shift)
        exp += shift
    else:
        mpz_set(res.value, man.value)
    # Strip trailing bits
    trail = mpz_scan1(res.value, 0)
    if 0 < trail < bc:
        mpz_tdiv_q_2exp(res.value, res.value, trail)
        exp += trail
    bc = mpz_sizeinbase(res.value, 2)
    return (sign, res, int(exp), bc)

cdef mpfr_from_mpfval(mpfr_t res, tuple x):
    """
    Set value of an MPFR number (in place) to that of a given mpmath mpf
    data tuple.
    """
    cdef int sign
    cdef Integer man
    cdef long exp
    cdef long bc
    sign, man, exp, bc = x
    if man.__nonzero__():
        mpfr_set_z(res, man.value, GMP_RNDZ)
        if sign:
            mpfr_neg(res, res, GMP_RNDZ)
        mpfr_mul_2si(res, res, exp, GMP_RNDZ)
        return
    from mpmath.libmpf import finf, fninf
    if exp == 0:
        mpfr_set_ui(res, 0, GMP_RNDZ)
    elif x == finf:
        mpfr_set_inf(res, 1)
    elif x == fninf:
        mpfr_set_inf(res, -1)
    else:
        mpfr_set_nan(res)

cdef mpfr_to_mpfval(mpfr_t value):
    """
    Given an MPFR value, return an mpmath mpf data tuple representing
    the same number.
    """
    if mpfr_nan_p(value):
        from mpmath.libmpf import fnan
        return fnan
    if mpfr_inf_p(value):
        from mpmath.libmpf import finf, fninf
        if mpfr_sgn(value) > 0:
            return finf
        else:
            return fninf
    if mpfr_sgn(value) == 0:
        from mpmath.libmpf import fzero
        return fzero
    sign = 0
    cdef Integer man = PY_NEW(Integer)
    exp = mpfr_get_z_exp(man.value, value)
    if mpz_sgn(man.value) < 0:
        mpz_neg(man.value, man.value)
        sign = 1
    cdef unsigned long trailing
    trailing = mpz_scan1(man.value, 0)
    if trailing:
        mpz_tdiv_q_2exp(man.value, man.value, trailing)
        exp += trailing
    bc = mpz_sizeinbase(man.value, 2)
    return (sign, man, int(exp), bc)

def mpmath_to_sage(x, prec):
    """
    Convert any mpmath number (mpf or mpc) to a Sage RealNumber or
    ComplexNumber of the given precision.

    EXAMPLES::

        sage: import sage.libs.mpmath.all as a
        sage: a.mpmath_to_sage(a.mpf('2.5'), 53)
        2.50000000000000
        sage: a.mpmath_to_sage(a.mpc('2.5','-3.5'), 53)
        2.50000000000000 - 3.50000000000000*I
        sage: a.mpmath_to_sage(a.mpf('inf'), 53)
        +infinity
        sage: a.mpmath_to_sage(a.mpf('-inf'), 53)
        -infinity
        sage: a.mpmath_to_sage(a.mpf('nan'), 53)
        NaN
        sage: a.mpmath_to_sage(a.mpf('0'), 53)
        0.000000000000000
    """
    cdef RealNumber y
    cdef ComplexNumber z
    if hasattr(x, "_mpf_"):
        y = RealField(prec)()
        mpfr_from_mpfval(y.value, x._mpf_)
        return y
    elif hasattr(x, "_mpc_"):
        z = ComplexField(prec)(0)
        re, im = x._mpc_
        mpfr_from_mpfval(z.__re, re)
        mpfr_from_mpfval(z.__im, im)
        return z
    else:
        raise TypeError("cannot convert %r to Sage", x)

def sage_to_mpmath(x, prec):
    """
    Convert any Sage number that can be coerced into a RealNumber
    or ComplexNumber of the given precision into an mpmath mpf or mpc.
    Integers are currently converted to int.

    EXAMPLES::

        sage: import sage.libs.mpmath.all as a
        sage: a.mp.dps = 15
        sage: print a.sage_to_mpmath(2/3, 53)
        0.666666666666667
        sage: print a.sage_to_mpmath(2./3, 53)
        0.666666666666667
        sage: print a.sage_to_mpmath(3+4*I, 53)
        (3.0 + 4.0j)
        sage: print a.sage_to_mpmath(1+pi, 53)
        4.14159265358979
        sage: a.sage_to_mpmath(infinity, 53)
        mpf('+inf')
        sage: a.sage_to_mpmath(-infinity, 53)
        mpf('-inf')
        sage: a.sage_to_mpmath(NaN, 53)
        mpf('nan')
        sage: a.sage_to_mpmath(0, 53)
        0
    """
    from mpmath.mptypes import make_mpf, make_mpc
    cdef RealNumber y
    if isinstance(x, Element):
        if PY_TYPE_CHECK(x, Integer):
            return int(<Integer>x)
        try:
            if PY_TYPE_CHECK(x, RealNumber):
                return x._mpmath_()
            else:
                x = RealField(prec)(x)
                return x._mpmath_()
        except TypeError:
            if PY_TYPE_CHECK(x, ComplexNumber):
                return x._mpmath_()
            else:
                x = ComplexField(prec)(x)
                return x._mpmath_()
    if PY_TYPE_CHECK(x, tuple) or PY_TYPE_CHECK(x, list):
        return type(x)([sage_to_mpmath(v, prec) for v in x])
    return x

def call(func, *args, **kwargs):
    """
    Call an mpmath function with Sage objects as inputs and
    convert the result back to a Sage real or complex number. Use the
    keyword argument prec=n to set the working precision in bits
    (by default mpmath's global working precision is used).

    Arguments should be Sage objects that can be coerced into RealField
    or ComplexField elements. Arguments may also be tuples/lists (which
    are converted recursively), or any type that mpmath understands
    natively (e.g. Python floats, strings for options).

    EXAMPLES::

        sage: import sage.libs.mpmath.all as a
        sage: a.mp.prec = 53
        sage: a.call(a.erf, 3+4*I)
        -120.186991395079 - 27.7503372936239*I
        sage: a.call(a.polylog, 2, 1/3+4/5*I)
        0.153548951541433 + 0.875114412499637*I
        sage: a.call(a.barnesg, 3+4*I)
        -0.000676375932234244 - 0.0000442236140124728*I
        sage: a.call(a.barnesg, -4)
        0.000000000000000
        sage: a.call(a.hyper, [2,3], [4,5], 1/3)
        1.10703578162508
        sage: a.call(a.hyper, [2,3], [4,(2,3)], 1/3)
        1.95762943509305
        sage: a.call(a.quad, a.erf, [0,1])
        0.486064958112256
        sage: a.call(a.gammainc, 3+4*I, 2/3, 1-pi*I, prec=100)
        -274.18871130777160922270612331 + 101.59521032382593402947725236*I
        sage: x = (3+4*I).n(100)
        sage: y = (2/3).n(100)
        sage: z = (1-pi*I).n(100)
        sage: a.call(a.gammainc, x, y, z, prec=100)
        -274.18871130777160922270612331 + 101.59521032382593402947725236*I
        sage: a.call(a.erf, infinity)
        1.00000000000000
        sage: a.call(a.erf, -infinity)
        -1.00000000000000
        sage: a.call(a.gamma, infinity)
        +infinity
    """
    from mpmath import mp
    orig = mp.prec
    prec = kwargs.pop('prec', orig)
    prec2 = prec + 20
    args = [sage_to_mpmath(x, prec2) for x in args]
    if kwargs:
        kwargs = dict([(key, sage_to_mpmath(value, prec2)) for (key, value) in \
            kwargs.items()])
    try:
        mp.prec = prec
        y = func(*args, **kwargs)
    finally:
        mp.prec = orig
    return mpmath_to_sage(y, prec)


