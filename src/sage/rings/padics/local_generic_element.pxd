import sage.structure.element
cimport sage.structure.element
from sage.structure.element cimport RingElement, ModuleElement, CommutativeRingElement

cdef class LocalGenericElement(CommutativeRingElement):
    cdef RingElement _div_c_impl(self, RingElement right)
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)
