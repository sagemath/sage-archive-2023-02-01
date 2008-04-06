from sage.structure.element cimport Element
from sage.structure.parent cimport Parent

cdef class Morphism(Element):
    cdef Parent _domain
    cdef Parent _codomain

# TODO: remove this requirement when we better understand pickling
#    cdef __dict__

    cdef _update_slots(self, _dict)
    cdef _extra_slots(self, _dict)

    # these methods assume x is an element of domain, and returns an element with parent codomain
    cdef Element _call_c(self, x)
    cdef Element _call_c_impl(self, Element x)

cdef class Section(Morphism):
    cdef Morphism _morphism

cdef class FormalCoercionMorphism(Morphism):
    cdef Element _call_c(self, x)

cdef class FormalCompositeMorphism(Morphism):
    cdef Morphism __first
    cdef Morphism __second
    pass