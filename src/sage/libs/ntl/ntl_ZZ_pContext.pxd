include "decl.pxi"
include "../../ext/cdefs.pxi"

cdef class ntl_ZZ_pContext_class:
    cdef ZZ_pContext_c x
    cdef void restore_c(self)