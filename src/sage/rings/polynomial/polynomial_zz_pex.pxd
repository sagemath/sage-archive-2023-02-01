from sage.libs.ntl.ntl_ZZ_pEX_decl cimport ZZ_pEX_c

ctypedef ZZ_pEX_c celement

include "polynomial_template_header.pxi"

cdef class Polynomial_ZZ_pX(Polynomial_template):
    pass

