include "decl.pxi"
include "../../ext/cdefs.pxi"

from ntl_GF2X cimport ntl_GF2X

cdef class ntl_GF2E(ntl_GF2X):
    cdef GF2E_c gf2e_x
