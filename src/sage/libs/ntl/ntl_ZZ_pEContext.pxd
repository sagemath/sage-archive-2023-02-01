include "decl.pxi"
include "../../ext/cdefs.pxi"

from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class

cdef class ntl_ZZ_pEContext_class:
    cdef ZZ_pEContext_c x
    cdef ntl_ZZ_pContext_class pc
    cdef void restore_c(self)
    cdef f