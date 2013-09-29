from sage.structure.parent cimport Parent
from sage.structure.element cimport Element

cdef class Map(Element):
    cdef public int _coerce_cost # a rough measure of the cost of using this morphism in the coercion system.
                          # 10 by default, 100 if a DefaultCoercionMorphism, 10000 if inexact.

    cdef _update_slots(self, _dict)
    cdef _extra_slots(self, _dict)

    # these methods require x is an element of domain, and returns an element with parent codomain
    cpdef Element _call_(self, x)
    cpdef Element _call_with_args(self, x, args=*, kwds=*)

    cdef public domain   # will be either a weakref or a constant map
    cdef public codomain

    cdef public _repr_type_str


cdef class Section(Map):
    cdef Map _inverse

cdef class FormalCompositeMap(Map):
    cdef Map __first
    cdef Map __second

