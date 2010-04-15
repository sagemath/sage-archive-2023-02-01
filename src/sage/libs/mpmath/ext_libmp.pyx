"""
Faster versions of some key functions in mpmath.libmp.
"""

include "../../ext/stdsage.pxi"
from ext_impl cimport *
from sage.libs.gmp.all cimport *
from sage.rings.integer cimport Integer

# Note: not thread-safe
cdef MPF tmp1
cdef MPF tmp2
MPF_init(&tmp1)
MPF_init(&tmp2)

def mpf_add(tuple x, tuple y, int prec=0, str rnd='d'):
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    MPF_set_tuple(&tmp2, y)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_add(&tmp1, &tmp1, &tmp2, opts)
    return MPF_to_tuple(&tmp1)

def mpf_sub(tuple x, tuple y, int prec=0, str rnd='d'):
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    MPF_set_tuple(&tmp2, y)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_sub(&tmp1, &tmp1, &tmp2, opts)
    return MPF_to_tuple(&tmp1)

def mpf_mul(tuple x, tuple y, int prec=0, str rnd='d'):
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    MPF_set_tuple(&tmp2, y)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_mul(&tmp1, &tmp1, &tmp2, opts)
    return MPF_to_tuple(&tmp1)

def mpf_div(tuple x, tuple y, int prec, str rnd='d'):
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    MPF_set_tuple(&tmp2, y)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_div(&tmp1, &tmp1, &tmp2, opts)
    return MPF_to_tuple(&tmp1)

def mpf_sqrt(tuple x, int prec, str rnd='d'):
    if x[0]:
        import mpmath.libmp as libmp
        raise libmp.ComplexResult("square root of a negative number")
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_sqrt(&tmp1, &tmp1, opts)
    return MPF_to_tuple(&tmp1)
