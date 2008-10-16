from sage.libs.ntl.ntl_GF2X_decl cimport GF2X_c

ctypedef GF2X_c celement

include "polynomial_template_header.pxi"

cdef class Polynomial_GF2X(Polynomial_template):
    pass

