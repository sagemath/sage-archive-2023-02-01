from sage.libs.arb.arb cimport arb_t
from sage.libs.mpfi cimport mpfi_t
from sage.rings.real_mpfi cimport RealIntervalFieldElement
from sage.structure.parent cimport Parent
from sage.structure.element cimport Element

cdef void mpfi_to_arb(arb_t target, const mpfi_t source, const long precision)
cdef void arb_to_mpfi(mpfi_t target, arb_t source, const long precision)

cdef class RealBall(Element):
    cdef arb_t value
    cdef RealBall _new(self)
    cpdef RealIntervalFieldElement _interval(self)
    cpdef RealBall psi(self)
