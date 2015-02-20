from sage.libs.arb.acb cimport acb_t
from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.rings.real_arb cimport RealBall
from sage.structure.element cimport Element
from sage.structure.parent cimport Parent

cdef void ComplexIntervalFieldElement_to_acb(
    acb_t target,
    ComplexIntervalFieldElement source)

cdef ComplexIntervalFieldElement acb_to_ComplexIntervalFieldElement(
    const acb_t source,
    Parent CIF)

cdef class ComplexBall(Element):
    cdef acb_t value
    cdef ComplexBall _new(self)
    cpdef ComplexIntervalFieldElement _interval(self)
    cpdef RealBall real(self)
    cpdef RealBall imag(self)
