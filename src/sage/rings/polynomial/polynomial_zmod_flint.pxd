cdef extern from "zn_poly/zn_poly.h":
     # This header needs to appear before the flint headers.
     pass

from sage.libs.flint.zmod_poly cimport zmod_poly_t, zmod_poly_struct

ctypedef zmod_poly_struct celement

include "polynomial_template_header.pxi"

cdef class Polynomial_zmod_flint(Polynomial_template):
    cdef _set_list(self, x)

