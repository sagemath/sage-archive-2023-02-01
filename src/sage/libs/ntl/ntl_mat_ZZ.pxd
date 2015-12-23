include "decl.pxi"

cdef class ntl_mat_ZZ(object):
    cdef mat_ZZ_c x
    cdef long __nrows, __ncols
