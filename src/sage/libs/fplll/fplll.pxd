include "../../ext/cdefs.pxi"
include "fplll.pxi"

cdef class FP_LLL:
    cdef ZZ_mat *_lattice
    cdef int _check_precision(self, int precision) except -1
    cdef int _check_eta(self, float eta) except -1
    cdef int _check_delta(self, float delta) except -1
