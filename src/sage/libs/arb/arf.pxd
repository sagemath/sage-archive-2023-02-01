from sage.libs.mpfr cimport mpfr_t, mpfr_rnd_t

cdef extern from "arf.h":
    ctypedef struct arf_struct:
        pass
    ctypedef arf_struct[1] arf_t

    int arf_get_mpfr(mpfr_t x, const arf_t y, mpfr_rnd_t rnd);
