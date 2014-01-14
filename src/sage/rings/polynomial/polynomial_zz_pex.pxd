from sage.libs.ntl.ntl_ZZ_pEX_decl cimport ZZ_pEX_c
from sage.libs.ntl.ntl_ZZ_pEContext_decl cimport ZZ_pEContext_c

ctypedef ZZ_pEX_c celement
ctypedef ZZ_pEContext_c *cparent

include "polynomial_template_header.pxi"

cdef class Polynomial_ZZ_pX(Polynomial_template):
    pass

