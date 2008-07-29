from sage.structure.element cimport Element
from sage.structure.parent cimport Parent

cdef class Morphism(Element):
    cdef Parent _domain
    cdef Parent _codomain

    cdef public int _coerce_cost # a rough measure of the cost of using this morphism in the coercion system.
                          # 10 by default, 100 if a DefaultCoercionMorphism, 10000 if inexact.

# TODO: remove this requirement when we better understand pickling
#    cdef __dict__

    cdef _update_slots(self, _dict)
    cdef _extra_slots(self, _dict)

    # these methods assume x is an element of domain, and returns an element with parent codomain
    cpdef Element _call_(self, x)
    cpdef Element _call_with_args(self, x, args=*, kwds=*)

    cpdef domain(self)

    cpdef codomain(self)

cdef class Section(Morphism):
    cdef Morphism _morphism

cdef class FormalCoercionMorphism(Morphism):
    pass

cdef class FormalCompositeMorphism(Morphism):
    cdef Morphism __first
    cdef Morphism __second
    pass

cdef class SetMorphism(Morphism):
    cdef object _function
