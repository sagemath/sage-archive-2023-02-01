from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from map cimport Map

cdef class Morphism(Map):
    pass

cdef class SetMorphism(Morphism):
    cdef object _function
    cpdef bool _eq_c_impl(left, Element right)
