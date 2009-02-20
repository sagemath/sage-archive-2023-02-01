cdef extern from "stdlib.h":
    void free(void *ptr)

from sage.libs.mpfr cimport *

include '../ext/cdefs.pxi'
include '../libs/pari/decl.pxi'

cimport sage.rings.ring
import  sage.rings.ring

cimport sage.structure.element
import  sage.structure.element

cdef class RealNumber(sage.structure.element.RingElement)  # forward decl

cdef class RealField(sage.rings.ring.Field):
    cdef int __prec
    cdef bint sci_not
    cdef mpfr_rnd_t rnd
    cdef object rnd_str
    cdef RealNumber _new(self)


cdef class RealNumber(sage.structure.element.RingElement):
    cdef mpfr_t value
    cdef char init
    cdef RealNumber _new(self)
    cdef _set(self, x, int base)
    cdef _set_from_GEN_REAL(self, GEN g)
    cdef RealNumber abs(RealNumber self)
