cdef extern from "stdlib.h":
    ctypedef int size_t
    void free(void *ptr)

include 'mpfr.pxi'

cimport sage.rings.ring
import  sage.rings.ring

cimport sage.structure.element
import  sage.structure.element

cdef class RealField(sage.rings.ring.Field):
    cdef int __prec, sci_not
    cdef mp_rnd_t rnd
    cdef object rnd_str

cdef class RealNumber(sage.structure.element.RingElement):
    cdef mpfr_t value
    cdef char init

    cdef RealNumber abs(RealNumber self)
    cdef int cmp(RealNumber self, RealNumber x)
