
include "decl.pxi"

cdef class ntl_ZZ:
    cdef ZZ_c x
    cdef public int get_as_int(ntl_ZZ self)
    cdef public void set_from_int(ntl_ZZ self, int value)
