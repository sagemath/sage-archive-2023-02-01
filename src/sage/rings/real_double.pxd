cimport sage.structure.element
cimport sage.rings.ring

from sage.structure.element cimport RingElement, ModuleElement, Element, FieldElement
from sage.rings.ring cimport Field

cdef extern from "limits.h":
    int INT_MAX
    double NAN

cdef class RealDoubleField_class(Field):
    cdef _new_c(self, double value)

cdef class RealDoubleElement(FieldElement):
    cdef double _value
    cdef _new_c(self, double value)
    cdef double __nth_root(RealDoubleElement self, int n)
    cdef RealDoubleElement abs(RealDoubleElement self)


