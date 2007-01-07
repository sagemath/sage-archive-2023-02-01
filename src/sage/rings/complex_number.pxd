include '../ext/cdefs.pxi'
include 'mpfr.pxi'

cimport sage.structure.element
cimport sage.libs.pari.gen
cimport real_mpfr

cdef class ComplexNumber(sage.structure.element.RingElement):
    cdef mpfr_t __re
    cdef mpfr_t __im
    #cdef sage.rings.real_mpfr.RealNumber __re
    #cdef sage.rings.real_mpfr.RealNumber __im
    cdef object _multiplicative_order
    cdef int _prec

    cdef real_mpfr.RealNumber abs(ComplexNumber self)
    cdef real_mpfr.RealNumber norm(ComplexNumber self)

    cdef ComplexNumber _new(self)
