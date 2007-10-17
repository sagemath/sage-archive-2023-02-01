include "decl.pxi"
include "../../ext/cdefs.pxi"

cdef class ntl_mat_GF2E:
    cdef mat_GF2E_c x
    cdef long __nrows, __ncols
