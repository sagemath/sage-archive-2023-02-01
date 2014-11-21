from sage.libs.arb.fmpr cimport fmpr_t

cdef extern from "mag.h":
    ctypedef struct mag_struct:
        pass
    ctypedef mag_struct[1] mag_t

    void mag_get_fmpr(fmpr_t x, const mag_t r)
