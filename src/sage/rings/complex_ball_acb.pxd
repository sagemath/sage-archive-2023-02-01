from sage.libs.arb.acb cimport acb_t
from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.structure.sage_object cimport SageObject

cdef void ComplexIntervalFieldElement_to_acb(
    acb_t target,
    ComplexIntervalFieldElement source)

cdef ComplexIntervalFieldElement acb_to_ComplexIntervalFieldElement(
    const acb_t source,
    const unsigned long precision)

cdef class Acb(SageObject):
    cdef acb_t value
    cdef unsigned long _precision_
    cpdef ComplexIntervalFieldElement ComplexIntervalFieldElement(self)
