"""
Implements mpf and mpc types, with binary operations and support
for interaction with other types. Also implements the main
context class, and related utilities.
"""

include '../../ext/interrupt.pxi'
include "../../ext/stdsage.pxi"
include "../../ext/python_int.pxi"
include "../../ext/python_long.pxi"
include "../../ext/python_float.pxi"
include "../../ext/python_complex.pxi"
include "../../ext/python_number.pxi"

from sage.libs.gmp.all cimport *
from sage.rings.integer cimport Integer

cdef extern from "mpz_pylong.h":
    cdef mpz_get_pylong(mpz_t src)
    cdef mpz_get_pyintlong(mpz_t src)
    cdef int mpz_set_pylong(mpz_t dst, src) except -1
    cdef long mpz_pythonhash(mpz_t src)

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

from ext_impl cimport *

import mpmath.rational as rationallib
import mpmath.libmp as libmp
import mpmath.function_docs as function_docs
from mpmath.libmp import to_str
from mpmath.libmp import repr_dps, prec_to_dps, dps_to_prec

DEF OP_ADD = 0
DEF OP_SUB = 1
DEF OP_MUL = 2
DEF OP_DIV = 3
DEF OP_POW = 4
DEF OP_MOD = 5
DEF OP_RICHCMP = 6
DEF OP_EQ = (OP_RICHCMP + 2)
DEF OP_NE = (OP_RICHCMP + 3)
DEF OP_LT = (OP_RICHCMP + 0)
DEF OP_GT = (OP_RICHCMP + 4)
DEF OP_LE = (OP_RICHCMP + 1)
DEF OP_GE = (OP_RICHCMP + 5)

cdef MPopts opts_exact
cdef MPopts opts_double_precision
cdef MPopts opts_mini_prec

opts_exact.prec = 0
opts_exact.rounding = ROUND_N
opts_double_precision.prec = 53
opts_double_precision.rounding = ROUND_N
opts_mini_prec.prec = 5
opts_mini_prec.rounding = ROUND_D

cdef MPF MPF_C_0
cdef MPF MPF_C_1
cdef MPF MPF_C_2

MPF_init(&MPF_C_0); MPF_set_zero(&MPF_C_0);
MPF_init(&MPF_C_1); MPF_set_si(&MPF_C_1, 1);
MPF_init(&MPF_C_2); MPF_set_si(&MPF_C_2, 2);

# Temporaries used for operands in binary operations
cdef mpz_t tmp_mpz
mpz_init(tmp_mpz)

cdef MPF tmp1
cdef MPF tmp2
cdef MPF tmp_opx_re
cdef MPF tmp_opx_im
cdef MPF tmp_opy_re
cdef MPF tmp_opy_im

MPF_init(&tmp1)
MPF_init(&tmp2)
MPF_init(&tmp_opx_re)
MPF_init(&tmp_opx_im)
MPF_init(&tmp_opy_re)
MPF_init(&tmp_opy_im)

cdef class Context
cdef class mpnumber
cdef class mpf_base
cdef class mpf
cdef class mpc
cdef class constant
cdef class wrapped_libmp_function
cdef class wrapped_specfun

cdef __isint(MPF *v):
    return v.special == S_ZERO or (v.special == S_NORMAL and mpz_sgn(v.exp) >= 0)

cdef int MPF_set_any(MPF *re, MPF *im, x, MPopts opts, bint str_tuple_ok) except -1:
    """
    Sets re + im*i = x, where x is any Python number.

    Returns 0 if unable to coerce x; 1 if x is real and re was set;
    2 if x is complex and both re and im were set.

    If str_tuple_ok=True, strings and tuples are accepted and converted
    (useful for parsing arguments, but not for arithmetic operands).
    """
    if PY_TYPE_CHECK(x, mpf):
        MPF_set(re, &(<mpf>x).value)
        return 1
    if PY_TYPE_CHECK(x, mpc):
        MPF_set(re, &(<mpc>x).re)
        MPF_set(im, &(<mpc>x).im)
        return 2
    if PyInt_Check(x) or PyLong_Check(x) or PY_TYPE_CHECK(x, Integer):
        MPF_set_int(re, x)
        return 1
    if PyFloat_Check(x):
        MPF_set_double(re, x)
        return 1
    if PyComplex_Check(x):
        MPF_set_double(re, x.real)
        MPF_set_double(im, x.imag)
        return 2
    if PY_TYPE_CHECK(x, constant):
        MPF_set_tuple(re, x.func(opts.prec, rndmode_to_python(opts.rounding)))
        return 1
    if hasattr(x, "_mpf_"):
        MPF_set_tuple(re, x._mpf_)
        return 1
    if hasattr(x, "_mpc_"):
        r, i = x._mpc_
        MPF_set_tuple(re, r)
        MPF_set_tuple(im, i)
        return 2
    if hasattr(x, "_mpmath_"):
        return MPF_set_any(re, im, x._mpmath_(opts.prec,
            rndmode_to_python(opts.rounding)), opts, False)
    if PY_TYPE_CHECK(x, rationallib.mpq):
        p, q = x._mpq_
        MPF_set_int(re, p)
        MPF_set_int(im, q)
        MPF_div(re, re, im, opts)
        #MPF_set_tuple(re, libmp.from_rational(p, q, opts.prec,
        #    rndmode_to_python(opts.rounding)))
        return 1
    if hasattr(x, '_mpi_'):
        a, b = x._mpi_
        if a == b:
            MPF_set_tuple(re, a)
            return 1
        raise ValueError("can only create mpf from zero-width interval")
    if str_tuple_ok:
        if PY_TYPE_CHECK(x, tuple):
            if len(x) == 2:
                MPF_set_man_exp(re, x[0], x[1])
                return 1
            elif len(x) == 4:
                MPF_set_tuple(re, x)
                return 1
        if isinstance(x, basestring):
            try:
                st = libmp.from_str(x, opts.prec,
                    rndmode_to_python(opts.rounding))
            except ValueError:
                return 0
            MPF_set_tuple(re, st)
            return 1
    return 0

cdef binop(int op, x, y, MPopts opts):
    cdef int typx
    cdef int typy
    cdef MPF xre, xim, yre, yim
    cdef mpf rr
    cdef mpc rc
    cdef MPopts altopts

    if PY_TYPE_CHECK(x, mpf):
        xre = (<mpf>x).value
        typx = 1
    elif PY_TYPE_CHECK(x, mpc):
        xre = (<mpc>x).re
        xim = (<mpc>x).im
        typx = 2
    else:
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, opts, False)
        if typx == 0:
            return NotImplemented
        xre = tmp_opx_re
        xim = tmp_opx_im

    if PY_TYPE_CHECK(y, mpf):
        yre = (<mpf>y).value
        typy = 1
    elif PY_TYPE_CHECK(y, mpc):
        yre = (<mpc>y).re
        yim = (<mpc>y).im
        typy = 2
    else:
        typy = MPF_set_any(&tmp_opy_re, &tmp_opy_im, y, opts, False)
        if typy == 0:
            return NotImplemented
        yre = tmp_opy_re
        yim = tmp_opy_im

    if op == OP_ADD:
        if typx == 1 and typy == 1:
            # Real result
            rr = PY_NEW(mpf)
            MPF_add(&rr.value, &xre, &yre, opts)
            return rr
        else:
            # Complex result
            rc = PY_NEW(mpc)
            MPF_add(&rc.re, &xre, &yre, opts)
            if typx == 1:
                MPF_set(&rc.im, &yim)
                #MPF_normalize(&rc.im, opts)
            elif typy == 1:
                MPF_set(&rc.im, &xim)
                #MPF_normalize(&rc.im, opts)
            else:
                MPF_add(&rc.im, &xim, &yim, opts)
            return rc

    elif op == OP_SUB:
        if typx == 1 and typy == 1:
            # Real result
            rr = PY_NEW(mpf)
            MPF_sub(&rr.value, &xre, &yre, opts)
            return rr
        else:
            # Complex result
            rc = PY_NEW(mpc)
            MPF_sub(&rc.re, &xre, &yre, opts)
            if typx == 1:
                MPF_neg(&rc.im, &yim)
                MPF_normalize(&rc.im, opts)
            elif typy == 1:
                MPF_set(&rc.im, &xim)
                MPF_normalize(&rc.im, opts)
            else:
                MPF_sub(&rc.im, &xim, &yim, opts)
            return rc

    elif op == OP_MUL:
        if typx == 1 and typy == 1:
            # Real result
            rr = PY_NEW(mpf)
            MPF_mul(&rr.value, &xre, &yre, opts)
            return rr
        else:
            # Complex result
            rc = PY_NEW(mpc)
            if typx == 1:
                MPF_mul(&rc.re, &yre, &xre, opts)
                MPF_mul(&rc.im, &yim, &xre, opts)
            elif typy == 1:
                MPF_mul(&rc.re, &xre, &yre, opts)
                MPF_mul(&rc.im, &xim, &yre, opts)
            else:
                # a*c - b*d
                MPF_mul(&rc.re, &xre,   &yre,  opts_exact)
                MPF_mul(&tmp1,  &xim,   &yim,  opts_exact)
                MPF_sub(&rc.re, &rc.re, &tmp1, opts)
                # a*d + b*c
                MPF_mul(&rc.im, &xre,   &yim,  opts_exact)
                MPF_mul(&tmp1,  &xim,   &yre,  opts_exact)
                MPF_add(&rc.im, &rc.im, &tmp1, opts)
            return rc


    elif op == OP_DIV:
        if typx == 1 and typy == 1:
            # Real result
            rr = PY_NEW(mpf)
            MPF_div(&rr.value, &xre, &yre, opts)
            return rr
        else:
            rc = PY_NEW(mpc)
            if typy == 1:
                MPF_div(&rc.re, &xre, &yre, opts)
                MPF_div(&rc.im, &xim, &yre, opts)
            else:
                if typx == 1:
                    xim = MPF_C_0
                altopts = opts
                altopts.prec += 10
                # m = c*c + d*d
                MPF_mul(&tmp1, &yre,  &yre, opts_exact)
                MPF_mul(&tmp2, &yim,  &yim, opts_exact)
                MPF_add(&tmp1, &tmp1, &tmp2, altopts)
                # (a*c+b*d)/m
                MPF_mul(&rc.re, &xre,   &yre, opts_exact)
                MPF_mul(&tmp2,  &xim,   &yim, opts_exact)
                MPF_add(&rc.re, &rc.re, &tmp2, altopts)
                MPF_div(&rc.re, &rc.re, &tmp1, opts)
                # (b*c-a*d)/m
                MPF_mul(&rc.im, &xim,   &yre, opts_exact)
                MPF_mul(&tmp2,  &xre,   &yim, opts_exact)
                MPF_sub(&rc.im, &rc.im, &tmp2, altopts)
                MPF_div(&rc.im, &rc.im, &tmp1, opts)
            return rc

    elif op == OP_POW:
        if typx == 1 and typy == 1:
            rr = PY_NEW(mpf)
            if not MPF_pow(&rr.value, &xre, &yre, opts):
                return rr
        if typx == 1: xim = MPF_C_0
        if typy == 1: yim = MPF_C_0
        rc = PY_NEW(mpc)
        MPF_complex_pow(&rc.re, &rc.im, &xre, &xim, &yre, &yim, opts)
        return rc

    elif op == OP_MOD:
        if typx != 1 or typx != 1:
            raise TypeError("mod for complex numbers")
        xret = MPF_to_tuple(&xre)
        yret = MPF_to_tuple(&yre)
        v = libmp.mpf_mod(xret, yret, opts.prec, rndmode_to_python(opts.rounding))
        rr = PY_NEW(mpf)
        MPF_set_tuple(&rr.value, v)
        return rr

    elif op == OP_EQ:
        if typx == 1 and typy == 1:
            return MPF_eq(&xre, &yre)
        if typx == 1:
            return MPF_eq(&xre, &yre) and MPF_eq(&yim, &MPF_C_0)
        if typy == 1:
            return MPF_eq(&xre, &yre) and MPF_eq(&xim, &MPF_C_0)
        return MPF_eq(&xre, &yre) and MPF_eq(&xim, &yim)

    elif op == OP_NE:
        if typx == 1 and typy == 1:
            return MPF_ne(&xre, &yre)
        if typx == 1:
            return MPF_ne(&xre, &yre) or MPF_ne(&yim, &MPF_C_0)
        if typy == 1:
            return MPF_ne(&xre, &yre) or MPF_ne(&xim, &MPF_C_0)
        return MPF_ne(&xre, &yre) or MPF_ne(&xim, &yim)

    elif op == OP_LT:
        if typx != 1 or typy != 1:
            raise ValueError("cannot compare complex numbers")
        return MPF_lt(&xre, &yre)

    elif op == OP_GT:
        if typx != 1 or typy != 1:
            raise ValueError("cannot compare complex numbers")
        return MPF_gt(&xre, &yre)

    elif op == OP_LE:
        if typx != 1 or typy != 1:
            raise ValueError("cannot compare complex numbers")
        return MPF_le(&xre, &yre)

    elif op == OP_GE:
        if typx != 1 or typy != 1:
            raise ValueError("cannot compare complex numbers")
        return MPF_ge(&xre, &yre)

    return NotImplemented


cdef MPopts global_opts

global_context = None

cdef class Context:
    cdef public mpf, mpc, constant #, def_mp_function
    cdef public trap_complex
    cdef public pretty

    def __cinit__(ctx):
        global global_opts, global_context
        global_opts = opts_double_precision
        global_context = ctx
        ctx.mpf = mpf
        ctx.mpc = mpc
        ctx.constant = constant
        #ctx.def_mp_function = def_mp_function
        ctx._mpq = rationallib.mpq

    def default(ctx):
        global_opts = opts_double_precision
        ctx.trap_complex = False
        ctx.pretty = False

    def _get_prec(ctx): return global_opts.prec
    def _set_prec(ctx, prec): global_opts.prec = prec
    def _set_dps(ctx, n): global_opts.prec = dps_to_prec(int(n))
    def _get_dps(ctx): return libmp.prec_to_dps(global_opts.prec)
    dps = property(_get_dps, _set_dps)
    prec = property(_get_prec, _set_prec)
    _dps = property(_get_dps, _set_dps)
    _prec = property(_get_prec, _set_prec)

    def _get_prec_rounding(ctx):
        return global_opts.prec, rndmode_to_python(global_opts.rounding)

    _prec_rounding = property(_get_prec_rounding)

    cpdef mpf make_mpf(ctx, tuple v):
        cdef mpf x
        x = PY_NEW(mpf)
        MPF_set_tuple(&x.value, v)
        return x

    cpdef mpc make_mpc(ctx, tuple v):
        cdef mpc x
        x = PY_NEW(mpc)
        MPF_set_tuple(&x.re, v[0])
        MPF_set_tuple(&x.im, v[1])
        return x

    def convert(ctx, x, strings=True):
        """
        Converts *x* to an ``mpf``, ``mpc`` or ``mpi``. If *x* is of type ``mpf``,
        ``mpc``, ``int``, ``float``, ``complex``, the conversion
        will be performed losslessly.

        If *x* is a string, the result will be rounded to the present
        working precision. Strings representing fractions or complex
        numbers are permitted.

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> mpmathify(3.5)
            mpf('3.5')
            >>> mpmathify('2.1')
            mpf('2.1000000000000001')
            >>> mpmathify('3/4')
            mpf('0.75')
            >>> mpmathify('2+3j')
            mpc(real='2.0', imag='3.0')

        """
        cdef mpf rr
        cdef mpc rc
        if PY_TYPE_CHECK(x, mpnumber):
            return x
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, strings)
        if typx == 1:
            rr = PY_NEW(mpf)
            MPF_set(&rr.value, &tmp_opx_re)
            return rr
        if typx == 2:
            rc = PY_NEW(mpc)
            MPF_set(&rc.re, &tmp_opx_re)
            MPF_set(&rc.im, &tmp_opx_im)
            return rc
        return ctx._convert_fallback(x, strings)

    def isnan(ctx, x):
        """
        For an ``mpf`` *x*, determines whether *x* is not-a-number (nan)::

            >>> from mpmath import *
            >>> isnan(nan), isnan(3)
            (True, False)
        """
        cdef int s, t, typ
        if PY_TYPE_CHECK(x, mpf):
            return (<mpf>x).value.special == S_NAN
        if PY_TYPE_CHECK(x, mpc):
            s = (<mpc>x).re.special
            t = (<mpc>x).im.special
            return s == S_NAN or t == S_NAN
        if PyInt_CheckExact(x) or PyLong_CheckExact(x) or PY_TYPE_CHECK(x, Integer) \
            or PY_TYPE_CHECK(x, rationallib.mpq):
            return False
        typ = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 0)
        if typ == 1:
            s = tmp_opx_re.special
            return s == S_NAN
        if typ == 2:
            s = tmp_opx_re.special
            t = tmp_opx_im.special
            return s == S_NAN or t == S_NAN
        raise TypeError("isnan() needs a number as input")

    def isinf(ctx, x):
        """
        Return *True* if the absolute value of *x* is infinite;
        otherwise return *False*::

            >>> from mpmath import *
            >>> isinf(inf)
            True
            >>> isinf(-inf)
            True
            >>> isinf(3)
            False
            >>> isinf(3+4j)
            False
            >>> isinf(mpc(3,inf))
            True
            >>> isinf(mpc(inf,3))
            True

        """
        cdef int s, t, typ
        if PY_TYPE_CHECK(x, mpf):
            s = (<mpf>x).value.special
            return s == S_INF or s == S_NINF
        if PY_TYPE_CHECK(x, mpc):
            s = (<mpc>x).re.special
            t = (<mpc>x).im.special
            return s == S_INF or s == S_NINF or t == S_INF or t == S_NINF
        if PyInt_CheckExact(x) or PyLong_CheckExact(x) or PY_TYPE_CHECK(x, Integer) \
            or PY_TYPE_CHECK(x, rationallib.mpq):
            return False
        typ = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 0)
        if typ == 1:
            s = tmp_opx_re.special
            return s == S_INF or s == S_NINF
        if typ == 2:
            s = tmp_opx_re.special
            t = tmp_opx_im.special
            return s == S_INF or s == S_NINF or t == S_INF or t == S_NINF
        raise TypeError("isinf() needs a number as input")

    def isnormal(ctx, x):
        """
        Determine whether *x* is "normal" in the sense of floating-point
        representation; that is, return *False* if *x* is zero, an
        infinity or NaN; otherwise return *True*. By extension, a
        complex number *x* is considered "normal" if its magnitude is
        normal::

            >>> from mpmath import *
            >>> isnormal(3)
            True
            >>> isnormal(0)
            False
            >>> isnormal(inf); isnormal(-inf); isnormal(nan)
            False
            False
            False
            >>> isnormal(0+0j)
            False
            >>> isnormal(0+3j)
            True
            >>> isnormal(mpc(2,nan))
            False
        """
        # TODO: optimize this
        if hasattr(x, "_mpf_"):
            return bool(x._mpf_[1])
        if hasattr(x, "_mpc_"):
            re, im = x._mpc_
            re_normal = bool(re[1])
            im_normal = bool(im[1])
            if re == libmp.fzero: return im_normal
            if im == libmp.fzero: return re_normal
            return re_normal and im_normal
        if PyInt_CheckExact(x) or PyLong_CheckExact(x) or PY_TYPE_CHECK(x, Integer) \
            or PY_TYPE_CHECK(x, rationallib.mpq):
            return bool(x)
        x = ctx.convert(x)
        if hasattr(x, '_mpf_') or hasattr(x, '_mpc_'):
            return ctx.isnormal(x)
        raise TypeError("isnormal() needs a number as input")

    def isint(ctx, x, gaussian=False):
        """
        Return *True* if *x* is integer-valued; otherwise return
        *False*::

            >>> from mpmath import *
            >>> isint(3)
            True
            >>> isint(mpf(3))
            True
            >>> isint(3.2)
            False
            >>> isint(inf)
            False

        Optionally, Gaussian integers can be checked for::

            >>> isint(3+0j)
            True
            >>> isint(3+2j)
            False
            >>> isint(3+2j, gaussian=True)
            True

        """
        cdef MPF v
        cdef MPF w
        cdef int typ
        if PyInt_CheckExact(x) or PyLong_CheckExact(x) or PY_TYPE_CHECK(x, Integer):
            return True
        if PY_TYPE_CHECK(x, mpf):
            v = (<mpf>x).value
            return __isint(&v)
        if PY_TYPE_CHECK(x, mpc):
            v = (<mpc>x).re
            w = (<mpc>x).im
            if gaussian:
                return __isint(&v) and __isint(&w)
            return (w.special == S_ZERO) and __isint(&v)
        if PY_TYPE_CHECK(x, rationallib.mpq):
            p, q = x._mpq_
            return not (p % q)
        typ = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 0)
        if typ == 1:
            return __isint(&tmp_opx_re)
        if typ == 2:
            v = tmp_opx_re
            w = tmp_opx_im
            if gaussian:
                return __isint(&v) and __isint(&w)
            return (w.special == S_ZERO) and __isint(&v)
        raise TypeError("isint() needs a number as input")

    def fsum(ctx, terms, bint absolute=False, bint squared=False):
        """
        Calculates a sum containing a finite number of terms (for infinite
        series, see :func:`nsum`). The terms will be converted to
        mpmath numbers. For len(terms) > 2, this function is generally
        faster and produces more accurate results than the builtin
        Python function :func:`sum`.

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> fsum([1, 2, 0.5, 7])
            mpf('10.5')

        With squared=True each term is squared, and with absolute=True
        the absolute value of each term is used.
        """
        cdef MPF sre, sim, tre, tim, tmp
        cdef mpf rr
        cdef mpc rc
        cdef MPopts workopts
        cdef int styp, ttyp
        workopts = global_opts
        workopts.prec = workopts.prec * 2 + 50
        workopts.rounding = ROUND_D
        unknown = global_context.zero
        try:  # Way down, there is a ``finally`` with sig_off()
            sig_on()
            MPF_init(&sre)
            MPF_init(&sim)
            MPF_init(&tre)
            MPF_init(&tim)
            MPF_init(&tmp)
            styp = 1
            for term in terms:
                ttyp = MPF_set_any(&tre, &tim, term, workopts, 0)
                if ttyp == 0:
                    if absolute: term = ctx.absmax(term)
                    if squared: term = term**2
                    unknown += term
                    continue
                if absolute:
                    if squared:
                        if ttyp == 1:
                            MPF_mul(&tre, &tre, &tre, opts_exact)
                            MPF_add(&sre, &sre, &tre, workopts)
                        elif ttyp == 2:
                            # |(a+bi)^2| = a^2+b^2
                            MPF_mul(&tre, &tre, &tre, opts_exact)
                            MPF_add(&sre, &sre, &tre, workopts)
                            MPF_mul(&tim, &tim, &tim, opts_exact)
                            MPF_add(&sre, &sre, &tim, workopts)
                    else:
                        if ttyp == 1:
                            MPF_abs(&tre, &tre)
                            MPF_add(&sre, &sre, &tre, workopts)
                        elif ttyp == 2:
                            # |a+bi| = sqrt(a^2+b^2)
                            MPF_mul(&tre, &tre, &tre, opts_exact)
                            MPF_mul(&tim, &tim, &tim, opts_exact)
                            MPF_add(&tre, &tre, &tim, workopts)
                            MPF_sqrt(&tre, &tre, workopts)
                            MPF_add(&sre, &sre, &tre, workopts)
                elif squared:
                    if ttyp == 1:
                        MPF_mul(&tre, &tre, &tre, opts_exact)
                        MPF_add(&sre, &sre, &tre, workopts)
                    elif ttyp == 2:
                        # (a+bi)^2 = a^2-b^2 + 2i*ab
                        MPF_mul(&tmp, &tre, &tim, opts_exact)
                        MPF_mul(&tmp, &tmp, &MPF_C_2, opts_exact)
                        MPF_add(&sim, &sim, &tmp, workopts)
                        MPF_mul(&tre, &tre, &tre, opts_exact)
                        MPF_add(&sre, &sre, &tre, workopts)
                        MPF_mul(&tim, &tim, &tim, opts_exact)
                        MPF_sub(&sre, &sre, &tim, workopts)
                        styp = 2
                else:
                    if ttyp == 1:
                        MPF_add(&sre, &sre, &tre, workopts)
                    elif ttyp == 2:
                        MPF_add(&sre, &sre, &tre, workopts)
                        MPF_add(&sim, &sim, &tim, workopts)
                        styp = 2
            MPF_clear(&tre)
            MPF_clear(&tim)
            if styp == 1:
                rr = PY_NEW(mpf)
                MPF_set(&rr.value, &sre)
                MPF_clear(&sre)
                MPF_clear(&sim)
                MPF_normalize(&rr.value, global_opts)
                if unknown is not global_context.zero:
                    return ctx._stupid_add(rr, unknown)
                return rr
            elif styp == 2:
                rc = PY_NEW(mpc)
                MPF_set(&rc.re, &sre)
                MPF_set(&rc.im, &sim)
                MPF_clear(&sre)
                MPF_clear(&sim)
                MPF_normalize(&rc.re, global_opts)
                MPF_normalize(&rc.im, global_opts)
                if unknown is not global_context.zero:
                    return ctx._stupid_add(rc, unknown)
                return rc
            else:
                MPF_clear(&sre)
                MPF_clear(&sim)
                return +unknown
        finally:
            sig_off()

    def fdot(ctx, A, B=None, bint conjugate=False):
        r"""
        Computes the dot product of the iterables `A` and `B`,

        .. math ::

            \sum_{k=0} A_k B_k.

        Alternatively, :func:`fdot` accepts a single iterable of pairs.
        In other words, ``fdot(A,B)`` and ``fdot(zip(A,B))`` are equivalent.

        The elements are automatically converted to mpmath numbers.

        Examples::

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> A = [2, 1.5, 3]
            >>> B = [1, -1, 2]
            >>> fdot(A, B)
            mpf('6.5')
            >>> zip(A, B)
            [(2, 1), (1.5, -1), (3, 2)]
            >>> fdot(_)
            mpf('6.5')
        """
        if B:
            A = zip(A, B)
        cdef MPF sre, sim, tre, tim, ure, uim, tmp
        cdef mpf rr
        cdef mpc rc
        cdef MPopts workopts
        cdef int styp, ttyp, utyp
        MPF_init(&sre)
        MPF_init(&sim)
        MPF_init(&tre)
        MPF_init(&tim)
        MPF_init(&ure)
        MPF_init(&uim)
        MPF_init(&tmp)
        workopts = global_opts
        workopts.prec = workopts.prec * 2 + 50
        workopts.rounding = ROUND_D
        unknown = global_context.zero
        styp = 1
        for a, b in A:
            ttyp = MPF_set_any(&tre, &tim, a, workopts, 0)
            utyp = MPF_set_any(&ure, &uim, b, workopts, 0)
            if utyp == 2 and conjugate:
                MPF_neg(&uim, &uim);
            if ttyp == 0 or utyp == 0:
                if conjugate:
                    b = b.conj()
                unknown += a * b
                continue
            styp = max(styp, ttyp)
            styp = max(styp, utyp)
            if ttyp == 1:
                if utyp == 1:
                    MPF_mul(&tre, &tre, &ure, opts_exact)
                    MPF_add(&sre, &sre, &tre, workopts)
                elif utyp == 2:
                    MPF_mul(&ure, &ure, &tre, opts_exact)
                    MPF_mul(&uim, &uim, &tre, opts_exact)
                    MPF_add(&sre, &sre, &ure, workopts)
                    MPF_add(&sim, &sim, &uim, workopts)
                    styp = 2
            elif ttyp == 2:
                styp = 2
                if utyp == 1:
                    MPF_mul(&tre, &tre, &ure, opts_exact)
                    MPF_mul(&tim, &tim, &ure, opts_exact)
                    MPF_add(&sre, &sre, &tre, workopts)
                    MPF_add(&sim, &sim, &tim, workopts)
                elif utyp == 2:
                    MPF_mul(&tmp, &tre, &ure, opts_exact)
                    MPF_add(&sre, &sre, &tmp, workopts)
                    MPF_mul(&tmp, &tim, &uim, opts_exact)
                    MPF_sub(&sre, &sre, &tmp, workopts)
                    MPF_mul(&tmp, &tim, &ure, opts_exact)
                    MPF_add(&sim, &sim, &tmp, workopts)
                    MPF_mul(&tmp, &tre, &uim, opts_exact)
                    MPF_add(&sim, &sim, &tmp, workopts)
        MPF_clear(&tre)
        MPF_clear(&tim)
        MPF_clear(&ure)
        MPF_clear(&uim)
        MPF_clear(&tmp)
        if styp == 1:
            rr = PY_NEW(mpf)
            MPF_set(&rr.value, &sre)
            MPF_clear(&sre)
            MPF_clear(&sim)
            MPF_normalize(&rr.value, global_opts)
            if unknown is not global_context.zero:
                return ctx._stupid_add(rr, unknown)
            return rr
        elif styp == 2:
            rc = PY_NEW(mpc)
            MPF_set(&rc.re, &sre)
            MPF_set(&rc.im, &sim)
            MPF_clear(&sre)
            MPF_clear(&sim)
            MPF_normalize(&rc.re, global_opts)
            MPF_normalize(&rc.im, global_opts)
            if unknown is not global_context.zero:
                return ctx._stupid_add(rc, unknown)
            return rc
        else:
            MPF_clear(&sre)
            MPF_clear(&sim)
            return +unknown

    # Doing a+b directly doesn't work with mpi, presumably due to
    # Cython trying to be clever with the operation resolution
    cdef _stupid_add(ctx, a, b):
        return a + b

    def _convert_param(ctx, x):
        cdef MPF v
        cdef bint ismpf, ismpc
        if PyInt_Check(x) or PyLong_Check(x) or PY_TYPE_CHECK(x, Integer):
            return int(x), 'Z'
        if PY_TYPE_CHECK(x, tuple):
            p, q = x
            p = int(p)
            q = int(q)
            if not p % q:
                return p // q, 'Z'
            return rationallib.mpq((p,q)), 'Q'
        if PY_TYPE_CHECK(x, basestring) and '/' in x:
            p, q = x.split('/')
            p = int(p)
            q = int(q)
            if not p % q:
                return p // q, 'Z'
            return rationallib.mpq((p,q)), 'Q'
        if PY_TYPE_CHECK(x, constant):
            return x, 'R'
        ismpf = PY_TYPE_CHECK(x, mpf)
        ismpc = PY_TYPE_CHECK(x, mpc)
        if not (ismpf or ismpc):
            x = global_context.convert(x)
            ismpf = PY_TYPE_CHECK(x, mpf)
            ismpc = PY_TYPE_CHECK(x, mpc)
            if not (ismpf or ismpc):
                return x, 'U'
        if ismpf:
            v = (<mpf>x).value
        elif ismpc:
            if (<mpc>x).im.special != S_ZERO:
                return x, 'C'
            x = (<mpc>x).real
            v = (<mpc>x).re
        # A real number
        if v.special == S_ZERO:
            return 0, 'Z'
        if v.special != S_NORMAL:
            return x, 'U'
        if mpz_sgn(v.exp) >= 0:
            return (mpzi(v.man) << mpzi(v.exp)), 'Z'
        if mpz_cmp_si(v.exp, -4) > 0:
            p = mpzi(v.man)
            q = 1 << (-mpzi(v.exp))
            return rationallib.mpq((p,q)), 'Q'
        return x, 'R'

    def mag(ctx, x):
        """
        Quick logarithmic magnitude estimate of a number. Returns an
        integer or infinity `m` such that `|x| <= 2^m`. It is not
        guaranteed that `m` is an optimal bound, but it will never
        be too large by more than 2 (and probably not more than 1).

        EXAMPLES::

            sage: from mpmath import *
            sage: mp.pretty = True
            sage: mag(10), mag(10.0), mag(mpf(10)), int(ceil(log(10,2)))
            (4, 4, 4, 4)
            sage: mag(10j), mag(10+10j)
            (4, 5)
            sage: mag(0.01), int(ceil(log(0.01,2)))
            (-6, -6)
            sage: mag(0), mag(inf), mag(-inf), mag(nan)
            (-inf, +inf, +inf, nan)

        ::

            sage: class MyInt(int):
            ...       pass
            sage: class MyLong(long):
            ...       pass
            sage: class MyFloat(float):
            ...       pass
            sage: mag(MyInt(10)), mag(MyLong(10))
            (4, 4)

        """
        cdef int typ
        if PyInt_Check(x) or PyLong_Check(x) or PY_TYPE_CHECK(x, Integer):
            mpz_set_integer(tmp_opx_re.man, x)
            if mpz_sgn(tmp_opx_re.man) == 0:
                return global_context.ninf
            else:
                return mpz_sizeinbase(tmp_opx_re.man,2)
        if PY_TYPE_CHECK(x, rationallib.mpq):
            p, q = x._mpq_
            mpz_set_integer(tmp_opx_re.man, int(p))
            if mpz_sgn(tmp_opx_re.man) == 0:
                return global_context.ninf
            mpz_set_integer(tmp_opx_re.exp, int(q))
            return 1 + mpz_sizeinbase(tmp_opx_re.man,2) + mpz_sizeinbase(tmp_opx_re.exp,2)
        typ = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, False)
        if typ == 1:
            if tmp_opx_re.special == S_ZERO:
                return global_context.ninf
            if tmp_opx_re.special == S_INF or tmp_opx_re.special == S_NINF:
                return global_context.inf
            if tmp_opx_re.special != S_NORMAL:
                return global_context.nan
            mpz_add_ui(tmp_opx_re.exp, tmp_opx_re.exp, mpz_sizeinbase(tmp_opx_re.man, 2))
            return mpzi(tmp_opx_re.exp)
        if typ == 2:
            if tmp_opx_re.special == S_NAN or tmp_opx_im.special == S_NAN:
                return global_context.nan
            if tmp_opx_re.special == S_INF or tmp_opx_im.special == S_NINF or \
               tmp_opx_im.special == S_INF or tmp_opx_im.special == S_NINF:
                return global_context.inf
            if tmp_opx_re.special == S_ZERO:
                if tmp_opx_im.special == S_ZERO:
                    return global_context.ninf
                else:
                    mpz_add_ui(tmp_opx_im.exp, tmp_opx_im.exp, mpz_sizeinbase(tmp_opx_im.man, 2))
                    return mpzi(tmp_opx_im.exp)
            elif tmp_opx_im.special == S_ZERO:
                mpz_add_ui(tmp_opx_re.exp, tmp_opx_re.exp, mpz_sizeinbase(tmp_opx_re.man, 2))
                return mpzi(tmp_opx_re.exp)
            mpz_add_ui(tmp_opx_im.exp, tmp_opx_im.exp, mpz_sizeinbase(tmp_opx_im.man, 2))
            mpz_add_ui(tmp_opx_re.exp, tmp_opx_re.exp, mpz_sizeinbase(tmp_opx_re.man, 2))
            if mpz_cmp(tmp_opx_re.exp, tmp_opx_im.exp) >= 0:
                mpz_add_ui(tmp_opx_re.exp, tmp_opx_re.exp, 1)
                return mpzi(tmp_opx_re.exp)
            else:
                mpz_add_ui(tmp_opx_im.exp, tmp_opx_im.exp, 1)
                return mpzi(tmp_opx_im.exp)
        raise TypeError("requires an mpf/mpc")

    def _wrap_libmp_function(ctx, mpf_f, mpc_f=None, mpi_f=None, doc="<no doc>"):
        name = mpf_f.__name__[4:]
        doc = function_docs.__dict__.get(name, "Computes the %s of x" % doc)
        # workaround lack of closures in Cython
        f_cls = type(name, (wrapped_libmp_function,), {'__doc__':doc})
        f = f_cls(mpf_f, mpc_f, mpi_f, doc)
        return f

    @classmethod
    def _wrap_specfun(cls, name, f, wrap):
        doc = function_docs.__dict__.get(name, getattr(f, '__doc__', '<no doc>'))
        if wrap:
            # workaround lack of closures in Cython
            f_wrapped_cls = type(name, (wrapped_specfun,), {'__doc__':doc})
            f_wrapped = f_wrapped_cls(name, f)
        else:
            f_wrapped = f
        f_wrapped.__doc__ = doc
        setattr(cls, name, f_wrapped)

    cdef MPopts _fun_get_opts(ctx, kwargs):
        """
        Helper function that extracts precision and rounding information
        from kwargs, or returns the global working precision and rounding
        if no options are specified.
        """
        cdef MPopts opts
        opts.prec = global_opts.prec
        opts.rounding = global_opts.rounding
        if kwargs:
            if 'prec' in kwargs:
                opts.prec = int(kwargs['prec'])
            if 'dps'  in kwargs:
                opts.prec = libmp.dps_to_prec(int(kwargs['dps']))
            if 'rounding' in kwargs:
                opts.rounding = rndmode_from_python(kwargs['rounding'])
        return opts

    def _sage_sqrt(ctx, x, **kwargs):
        """
        Square root of an mpmath number x, using the Sage mpmath backend.

        EXAMPLES::

            sage: from mpmath import mp
            sage: mp.dps = 15
            sage: print mp.sqrt(2)   # indirect doctest
            1.4142135623731
            sage: print mp.sqrt(-2)
            (0.0 + 1.4142135623731j)
            sage: print mp.sqrt(2+2j)
            (1.55377397403004 + 0.643594252905583j)

        """
        cdef MPopts opts
        cdef int typx
        cdef mpf rr
        cdef mpc rc
        cdef tuple rev, imv
        opts = ctx._fun_get_opts(kwargs)
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, opts, 1)
        if typx == 1:
            rr = PY_NEW(mpf)
            if MPF_sqrt(&rr.value, &tmp_opx_re, opts):
                rc = PY_NEW(mpc)
                MPF_complex_sqrt(&rc.re, &rc.im, &tmp_opx_re, &MPF_C_0, opts)
                return rc
            return rr
        elif typx == 2:
            rc = PY_NEW(mpc)
            MPF_complex_sqrt(&rc.re, &rc.im, &tmp_opx_re, &tmp_opx_im, opts)
            return rc
        else:
            raise NotImplementedError("unknown argument")

    def _sage_exp(ctx, x, **kwargs):
        """
        Exponential function of an mpmath number x, using the Sage
        mpmath backend.

        EXAMPLES::

            sage: from mpmath import mp
            sage: mp.dps = 15
            sage: print mp.exp(2)   # indirect doctest
            7.38905609893065
            sage: print mp.exp(2+2j)
            (-3.07493232063936 + 6.71884969742825j)

        """
        cdef MPopts opts
        cdef int typx
        cdef mpf rr
        cdef mpc rc
        cdef tuple rev, imv
        opts = ctx._fun_get_opts(kwargs)
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, opts, 1)
        if typx == 1:
            rr = PY_NEW(mpf)
            MPF_exp(&rr.value, &tmp_opx_re, opts)
            return rr
        elif typx == 2:
            rc = PY_NEW(mpc)
            MPF_complex_exp(&rc.re, &rc.im, &tmp_opx_re, &tmp_opx_im, opts)
            return rc
        else:
            raise NotImplementedError("unknown argument")

    def _sage_cos(ctx, x, **kwargs):
        """
        Cosine of an mpmath number x, using the Sage mpmath backend.

        EXAMPLES::

            sage: from mpmath import mp
            sage: mp.dps = 15
            sage: print mp.cos(2)   # indirect doctest
            -0.416146836547142
            sage: print mp.cos(2+2j)
            (-1.56562583531574 - 3.29789483631124j)

        """
        cdef MPopts opts
        cdef int typx
        cdef mpf rr
        cdef mpc rc
        cdef tuple rev, imv
        opts = ctx._fun_get_opts(kwargs)
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 1)
        if typx == 1:
            rr = PY_NEW(mpf)
            MPF_cos(&rr.value, &tmp_opx_re, opts)
            return rr
        elif typx == 2:
            rev = MPF_to_tuple(&tmp_opx_re)
            imv = MPF_to_tuple(&tmp_opx_im)
            cxu = libmp.mpc_cos((rev, imv), opts.prec,
                rndmode_to_python(opts.rounding))
            rc = PY_NEW(mpc)
            MPF_set_tuple(&rc.re, cxu[0])
            MPF_set_tuple(&rc.im, cxu[1])
            return rc
        else:
            raise NotImplementedError("unknown argument")

    def _sage_sin(ctx, x, **kwargs):
        """
        Sine of an mpmath number x, using the Sage mpmath backend.

        EXAMPLES::

            sage: from mpmath import mp
            sage: mp.dps = 15
            sage: print mp.sin(2)   # indirect doctest
            0.909297426825682
            sage: print mp.sin(2+2j)
            (3.42095486111701 - 1.50930648532362j)

        """
        cdef MPopts opts
        cdef int typx
        cdef mpf rr
        cdef mpc rc
        cdef tuple rev, imv
        opts = ctx._fun_get_opts(kwargs)
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 1)
        if typx == 1:
            rr = PY_NEW(mpf)
            MPF_sin(&rr.value, &tmp_opx_re, opts)
            return rr
        elif typx == 2:
            rev = MPF_to_tuple(&tmp_opx_re)
            imv = MPF_to_tuple(&tmp_opx_im)
            cxu = libmp.mpc_sin((rev, imv), opts.prec,
                rndmode_to_python(opts.rounding))
            rc = PY_NEW(mpc)
            MPF_set_tuple(&rc.re, cxu[0])
            MPF_set_tuple(&rc.im, cxu[1])
            return rc
        else:
            raise NotImplementedError("unknown argument")

    def _sage_ln(ctx, x, **kwargs):
        """
        Natural logarithm of an mpmath number x, using the Sage mpmath backend.

        EXAMPLES::

            sage: from mpmath import mp
            sage: print mp.ln(2)   # indirect doctest
            0.693147180559945
            sage: print mp.ln(-2)
            (0.693147180559945 + 3.14159265358979j)
            sage: print mp.ln(2+2j)
            (1.03972077083992 + 0.785398163397448j)

        """
        cdef MPopts opts
        cdef int typx
        cdef mpf rr
        cdef mpc rc
        cdef tuple rev, imv
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 1)
        prec = global_opts.prec
        rounding = rndmode_to_python(global_opts.rounding)
        opts.prec = global_opts.prec
        opts.rounding = global_opts.rounding
        if kwargs:
            if 'prec' in kwargs: opts.prec = int(kwargs['prec'])
            if 'dps'  in kwargs: opts.prec = libmp.dps_to_prec(int(kwargs['dps']))
            if 'rounding' in kwargs: opts.rounding = rndmode_from_python(kwargs['rounding'])
        if typx == 1:
            if MPF_sgn(&tmp_opx_re) < 0:
                rc = PY_NEW(mpc)
                MPF_log(&rc.re, &tmp_opx_re, opts)
                MPF_set_pi(&rc.im, opts)
                return rc
            else:
                rr = PY_NEW(mpf)
                MPF_log(&rr.value, &tmp_opx_re, opts)
                return rr
        elif typx == 2:
            rev = MPF_to_tuple(&tmp_opx_re)
            imv = MPF_to_tuple(&tmp_opx_im)
            cxu = libmp.mpc_log((rev, imv), opts.prec, rndmode_to_python(opts.rounding))
            rc = PY_NEW(mpc)
            MPF_set_tuple(&rc.re, cxu[0])
            MPF_set_tuple(&rc.im, cxu[1])
            return rc

            #rc = PY_NEW(mpc)
            #MPF_complex_log(&rc.re, &rc.im, &tmp_opx_re, &tmp_opx_im, opts)
            #return rc
        else:
            raise NotImplementedError("unknown argument")


cdef class wrapped_libmp_function:

    cdef public mpf_f, mpc_f, mpi_f, name, __name__, __doc__

    def __init__(self, mpf_f, mpc_f=None, mpi_f=None, doc="<no doc>"):
        self.mpf_f = mpf_f
        self.mpc_f = mpc_f
        self.mpi_f = mpi_f
        self.name = self.__name__ = mpf_f.__name__[4:]
        self.__doc__ = function_docs.__dict__.get(self.name, "Computes the %s of x" % doc)

    def __call__(self, x, **kwargs):
        cdef int typx
        cdef tuple rev, imv, reu, cxu
        cdef mpf rr
        cdef mpc rc
        prec = global_opts.prec
        rounding = rndmode_to_python(global_opts.rounding)
        if kwargs:
            if 'prec' in kwargs: prec = int(kwargs['prec'])
            if 'dps'  in kwargs: prec = libmp.dps_to_prec(int(kwargs['dps']))
            if 'rounding' in kwargs: rounding = kwargs['rounding']
        typx = MPF_set_any(&tmp_opx_re, &tmp_opx_im, x, global_opts, 1)
        if typx == 1:
            rev = MPF_to_tuple(&tmp_opx_re)
            try:
                reu = self.mpf_f(rev, prec, rounding)
                rr = PY_NEW(mpf)
                MPF_set_tuple(&rr.value, reu)
                return rr
            except libmp.ComplexResult:
                if global_context.trap_complex:
                    raise
                cxu = self.mpc_f((rev, libmp.fzero), prec, rounding)
                rc = PY_NEW(mpc)
                MPF_set_tuple(&rc.re, cxu[0])
                MPF_set_tuple(&rc.im, cxu[1])
                return rc
        if typx == 2:
            rev = MPF_to_tuple(&tmp_opx_re)
            imv = MPF_to_tuple(&tmp_opx_im)
            cxu = self.mpc_f((rev, imv), prec, rounding)
            rc = PY_NEW(mpc)
            MPF_set_tuple(&rc.re, cxu[0])
            MPF_set_tuple(&rc.im, cxu[1])
            return rc
        x = global_context.convert(x)
        if hasattr(x, "_mpf_") or hasattr(x, "_mpc_"):
            return self.__call__(x, **kwargs)
        #if hasattr(x, "_mpi_"):
        #    if self.mpi_f:
        #        return global_context.make_mpi(self.mpi_f(x._mpi_, prec))
        raise NotImplementedError("%s of a %s" % (self.name, type(x)))



cdef class wrapped_specfun:
    cdef public f, name, __name__, __doc__

    def __init__(self, name, f):
        self.name = self.__name__ = name
        self.f = f
        self.__doc__ = function_docs.__dict__.get(name, "<no doc>")

    def __call__(self, *args, **kwargs):
        cdef int origprec
        args = [global_context.convert(a) for a in args]
        origprec = global_opts.prec
        global_opts.prec += 10
        try:
            retval = self.f(global_context, *args, **kwargs)
        finally:
            global_opts.prec = origprec
        return +retval


cdef class mpnumber:
    def __richcmp__(self, other, int op):
        return binop(OP_RICHCMP+op, self, other, global_opts)
    def __add__(self, other): return binop(OP_ADD, self, other, global_opts)
    def __sub__(self, other): return binop(OP_SUB, self, other, global_opts)
    def __mul__(self, other): return binop(OP_MUL, self, other, global_opts)
    def __div__(self, other): return binop(OP_DIV, self, other, global_opts)
    def __truediv__(self, other): return binop(OP_DIV, self, other, global_opts)
    def __mod__(self, other): return binop(OP_MOD, self, other, global_opts)
    def __pow__(self, other, mod):
        if mod is not None:
            raise ValueError("three-argument pow not supported")
        return binop(OP_POW, self, other, global_opts)
    def ae(s, t, rel_eps=None, abs_eps=None):
        return global_context.almosteq(s, t, rel_eps, abs_eps)


cdef class mpf_base(mpnumber):

    # Shared methods for mpf, constant. However, somehow some methods
    # (hash?, __richcmp__?) aren't inerited, so they have to
    # be defined multiple times. TODO: fix this.

    def __hash__(self):
        return libmp.mpf_hash(self._mpf_)

    def __repr__(self):
        if global_context.pretty:
            return self.__str__()
        n = repr_dps(global_opts.prec)
        return "mpf('%s')" % to_str(self._mpf_, n)

    def __str__(self):
        return to_str(self._mpf_, global_context._str_digits)

    @property
    def real(self): return self

    @property
    def imag(self): return global_context.zero

    def conjugate(self): return self

    @property
    def man(self): return self._mpf_[1]
    @property
    def exp(self): return self._mpf_[2]
    @property
    def bc(self): return self._mpf_[3]

    # XXX: optimize
    def __int__(self): return int(libmp.to_int(self._mpf_))
    def __long__(self): return long(self.__int__())
    def __float__(self): return libmp.to_float(self._mpf_)
    def __complex__(self): return complex(float(self))
    def to_fixed(self, prec): return libmp.to_fixed(self._mpf_, prec)
    def __getstate__(self): return libmp.to_pickable(self._mpf_)
    def __setstate__(self, val): self._mpf_ = libmp.from_pickable(val)


cdef class mpf(mpf_base):
    """
    An mpf instance holds a real-valued floating-point number. mpf:s
    work analogously to Python floats, but support arbitrary-precision
    arithmetic.
    """

    cdef MPF value

    def __init__(self, x=0, **kwargs):
        cdef MPopts opts
        opts = global_opts
        if kwargs:
            if 'prec' in kwargs: opts.prec = int(kwargs['prec'])
            if 'dps'  in kwargs: opts.prec = libmp.dps_to_prec(int(kwargs['dps']))
            if 'rounding' in kwargs: opts.rounding = rndmode_from_python(kwargs['rounding'])
        if MPF_set_any(&self.value, &self.value, x, opts, 1) != 1:
            raise TypeError

    def __reduce__(self): return (mpf, (), self._mpf_)
    def _get_mpf(self): return MPF_to_tuple(&self.value)
    def _set_mpf(self, v): MPF_set_tuple(&self.value, v)
    _mpf_ = property(_get_mpf, _set_mpf)

    def __nonzero__(self): return self.value.special != S_ZERO
    def __hash__(self): return libmp.mpf_hash(self._mpf_)

    @property
    def real(self): return self

    @property
    def imag(self): return global_context.zero

    def conjugate(self): return self

    @property
    def man(self): return self._mpf_[1]

    @property
    def exp(self): return self._mpf_[2]

    @property
    def bc(self): return self._mpf_[3]

    def to_fixed(self, long prec):
        # return libmp.to_fixed(self._mpf_, prec)
        MPF_to_fixed(tmp_mpz, &self.value, prec, False)
        cdef Integer r
        r = PY_NEW(Integer)
        mpz_set(r.value, tmp_mpz)
        return r

    def __int__(self):
        MPF_to_fixed(tmp_mpz, &self.value, 0, True)
        return mpzi(tmp_mpz)

    def __long__(self):
        r"""
        Convert this mpf value to a long.

        (Due to http://bugs.python.org/issue9869, to allow NZMATH to use
        this Sage-modified version of mpmath, it is vital that we
        return a long, not an int.)

        TESTS::

            sage: import mpmath
            sage: v = mpmath.mpf(2)
            sage: class MyLong(long):
            ...       pass
            sage: MyLong(v)
            2L
        """
        MPF_to_fixed(tmp_mpz, &self.value, 0, True)
        return mpzl(tmp_mpz)

    def __float__(self):
        return MPF_to_double(&self.value, False)

    def __getstate__(self):
        return libmp.to_pickable(self._mpf_)

    def __setstate__(self, val):
        self._mpf_ = libmp.from_pickable(val)

    def __cinit__(self):
        MPF_init(&self.value)

    def  __dealloc__(self):
        MPF_clear(&self.value)

    def __neg__(s):
        cdef mpf r = PY_NEW(mpf)
        MPF_neg(&r.value, &s.value)
        MPF_normalize(&r.value, global_opts)
        return r

    def __pos__(s):
        cdef mpf r = PY_NEW(mpf)
        MPF_set(&r.value, &s.value)
        MPF_normalize(&r.value, global_opts)
        return r

    def __abs__(s):
        cdef mpf r = PY_NEW(mpf)
        MPF_abs(&r.value, &s.value)
        MPF_normalize(&r.value, global_opts)
        return r

    def sqrt(s):
        cdef mpf r = PY_NEW(mpf)
        MPF_sqrt(&r.value, &s.value, global_opts)
        return r

    def __richcmp__(self, other, int op):
        return binop(OP_RICHCMP+op, self, other, global_opts)



cdef class constant(mpf_base):
    """
    Represents a mathematical constant with dynamic precision.
    When printed or used in an arithmetic operation, a constant
    is converted to a regular mpf at the working precision. A
    regular mpf can also be obtained using the operation +x.
    """

    cdef public name, func, __doc__

    def __init__(self, func, name, docname=''):
        self.name = name
        self.func = func
        self.__doc__ = getattr(function_docs, docname, '')

    def __call__(self, prec=None, dps=None, rounding=None):
        prec2 = global_opts.prec
        rounding2 = rndmode_to_python(global_opts.rounding)
        if not prec: prec = prec2
        if not rounding: rounding = rounding2
        if dps: prec = dps_to_prec(dps)
        return global_context.make_mpf(self.func(prec, rounding))

    @property
    def _mpf_(self):
        prec = global_opts.prec
        rounding = rndmode_to_python(global_opts.rounding)
        return self.func(prec, rounding)

    def __repr__(self):
        if global_context.pretty:
            return self.__str__()
        return "<%s: %s~>" % (self.name, global_context.nstr(self))

    def __nonzero__(self): return self._mpf_ != libmp.fzero
    def __neg__(self): return -mpf(self)
    def __pos__(self): return mpf(self)
    def __abs__(self): return abs(mpf(self))
    def sqrt(self): return mpf(self).sqrt()

    # XXX: optimize
    def to_fixed(self, prec):
        return libmp.to_fixed(self._mpf_, prec)
    def __getstate__(self):
        return libmp.to_pickable(self._mpf_)
    def __setstate__(self, val):
        self._mpf_ = libmp.from_pickable(val)

    # WHY?
    def __hash__(self):
        return libmp.mpf_hash(self._mpf_)
    def __richcmp__(self, other, int op):
        return binop(OP_RICHCMP+op, self, other, global_opts)


cdef class mpc(mpnumber):
    """
    An mpc represents a complex number using a pair of mpf:s (one
    for the real part and another for the imaginary part.) The mpc
    class behaves fairly similarly to Python's complex type.
    """

    cdef MPF re
    cdef MPF im

    def __init__(self, real=0, imag=0):
        cdef int typx, typy
        typx = MPF_set_any(&self.re, &self.im, real, global_opts, 1)
        if typx == 2:
            typy = 1
        else:
            typy = MPF_set_any(&self.im, &self.im, imag, global_opts, 1)
        if typx == 0 or typy != 1:
            raise TypeError

    def __cinit__(self):
        MPF_init(&self.re)
        MPF_init(&self.im)

    def  __dealloc__(self):
        MPF_clear(&self.re)
        MPF_clear(&self.im)

    def __reduce__(self):
        return (mpc, (), self._mpc_)

    def __setstate__(self, val):
        self._mpc_ = val[0], val[1]

    def __repr__(self):
        if global_context.pretty:
            return self.__str__()
        re, im = self._mpc_
        n = repr_dps(global_opts.prec)
        return "mpc(real='%s', imag='%s')" % (to_str(re, n), to_str(im, n))

    def __str__(s):
        return "(%s)" % libmp.mpc_to_str(s._mpc_, global_context._str_digits)

    def __nonzero__(self):
        return self.re.special != S_ZERO or self.im.special != S_ZERO

    #def __complex__(self):
    #    a, b = self._mpc_
    #    return complex(libmp.to_float(a), libmp.to_float(b))

    def __complex__(self):
        return complex(MPF_to_double(&self.re, False), MPF_to_double(&self.im, False))

    def _get_mpc(self):
        return MPF_to_tuple(&self.re), MPF_to_tuple(&self.im)

    def _set_mpc(self, tuple v):
        MPF_set_tuple(&self.re, v[0])
        MPF_set_tuple(&self.im, v[1])

    _mpc_ = property(_get_mpc, _set_mpc)

    @property
    def real(self):
        cdef mpf r = PY_NEW(mpf)
        MPF_set(&r.value, &self.re)
        return r

    @property
    def imag(self):
        cdef mpf r = PY_NEW(mpf)
        MPF_set(&r.value, &self.im)
        return r

    def __hash__(self):
        return libmp.mpc_hash(self._mpc_)

    def __neg__(s):
        cdef mpc r = PY_NEW(mpc)
        MPF_neg(&r.re, &s.re)
        MPF_neg(&r.im, &s.im)
        MPF_normalize(&r.re, global_opts)
        MPF_normalize(&r.im, global_opts)
        return r

    def conjugate(s):
        cdef mpc r = PY_NEW(mpc)
        MPF_set(&r.re, &s.re)
        MPF_neg(&r.im, &s.im)
        MPF_normalize(&r.re, global_opts)
        MPF_normalize(&r.im, global_opts)
        return r

    def __pos__(s):
        cdef mpc r = PY_NEW(mpc)
        MPF_set(&r.re, &s.re)
        MPF_set(&r.im, &s.im)
        MPF_normalize(&r.re, global_opts)
        MPF_normalize(&r.im, global_opts)
        return r

    def __abs__(s):
        cdef mpf r = PY_NEW(mpf)
        MPF_hypot(&r.value, &s.re, &s.im, global_opts)
        return r

    def __richcmp__(self, other, int op):
        return binop(OP_RICHCMP+op, self, other, global_opts)


def hypsum_internal(int p, int q, param_types, str ztype, coeffs, z,
    long prec, long wp, long epsshift, dict magnitude_check, kwargs):
    """
    Internal summation routine for hypergeometric series (wraps
    extension function MPF_hypsum to handle mpf/mpc types).

    EXAMPLES::

        sage: from mpmath import mp  # indirect doctest
        sage: mp.dps = 15
        sage: print mp.hyp1f1(1,2,3)
        6.36184564106256

    TODO: convert mpf/mpc parameters to fixed-point numbers here
    instead of converting to tuples within MPF_hypsum.
    """
    cdef mpf f
    cdef mpc c
    c = PY_NEW(mpc)
    have_complex, magn = MPF_hypsum(&c.re, &c.im, p, q, param_types, \
        ztype, coeffs, z, prec, wp, epsshift, magnitude_check, kwargs)
    if have_complex:
        v = c
    else:
        f = PY_NEW(mpf)
        MPF_set(&f.value, &c.re)
        v = f
    return v, have_complex, magn
