cdef extern from "stdlib.h":
    ctypedef int size_t
    void free(void *ptr)

cimport sage.structure.element
import  sage.structure.element

cimport sage.libs.pari.gen

cimport sage.rings.real_mpfr
import  sage.rings.real_mpfr

cdef class ComplexNumber(sage.structure.element.RingElement):
    cdef sage.rings.real_mpfr.RealNumber __re
    cdef sage.rings.real_mpfr.RealNumber __im
    cdef object __pari
    cdef object _multiplicative_order

    cdef sage.rings.real_mpfr.RealNumber abs(ComplexNumber self)
    cdef sage.rings.real_mpfr.RealNumber norm(ComplexNumber self)
