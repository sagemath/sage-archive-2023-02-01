cdef extern from "zn_poly/zn_poly.h":
     # This header needs to appear before the flint headers.
     pass

from sage.libs.flint.zmod_poly cimport zmod_poly_t, zmod_poly_struct

ctypedef zmod_poly_struct celement

include "polynomial_template_header.pxi"

cdef class Polynomial_zmod_flint(Polynomial_template):
    cdef Polynomial_template _new(self)
    cdef _set_list(self, x)
    cpdef _mul_short(self, Polynomial_zmod_flint other, length)
    cpdef _mul_short_opposite(self, Polynomial_zmod_flint other, length)
    cpdef rational_reconstruct(self, m, n_deg=?, d_deg=?)

