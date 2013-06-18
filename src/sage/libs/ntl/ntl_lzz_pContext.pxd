include "decl.pxi"
include "sage/ext/cdefs.pxi"

cdef class ntl_zz_pContext_class:
    cdef zz_pContext_c x
    cdef void restore_c(self)
    cdef long p
