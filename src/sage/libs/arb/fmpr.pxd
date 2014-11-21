from sage.libs.mpfr cimport mpfr_t, mpfr_rnd_t

cdef extern from "fmpr.h":
    ctypedef struct fmpr_struct:
        pass
    ctypedef fmpr_struct[1] fmpr_t

    void fmpr_init(fmpr_t x)
    void fmpr_clear(fmpr_t x)
    int fmpr_get_mpfr(mpfr_t x, const fmpr_t y, mpfr_rnd_t rnd);
