from sage.structure.element cimport Element
from morphism cimport Morphism
from functor cimport Functor

cdef class Action(Functor):
    cdef _G
    cdef _S
    cdef bint _is_left
    cdef Element _call_c(self, a, b)
    cdef Element _call_c_impl(self, Element a, Element b)

cdef class ActionEndomorphism(Morphism):
    pass