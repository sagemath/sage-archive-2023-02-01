from sage.libs.arb.arb cimport arb_t
from sage.rings.real_mpfi cimport RealIntervalFieldElement
from sage.structure.sage_object cimport SageObject

include 'mpfi.pxi'

cdef void mpfi_to_arb(arb_t target, const mpfi_t source, const unsigned long precision)
cdef void arb_to_mpfi(mpfi_t target, arb_t source, const unsigned long precision)

cdef class RealBallElement(SageObject):
     cdef arb_t value
     cdef unsigned int _precision_
     cpdef RealIntervalFieldElement RealIntervalFieldElement(self)
     cpdef RealBallElement psi(self)
