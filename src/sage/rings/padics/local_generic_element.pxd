import sage.structure.element
cimport sage.structure.element
from sage.structure.element cimport RingElement, ModuleElement, CommutativeRingElement

cdef class LocalGenericElement(CommutativeRingElement):
    cpdef RingElement _div_(self, RingElement right)
    cpdef ModuleElement _sub_(self, ModuleElement right)
