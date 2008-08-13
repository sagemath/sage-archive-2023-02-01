from sage.structure.element cimport Element
from morphism cimport Morphism
from functor cimport Functor

cdef class Action(Functor):
    cdef G
    cdef S
    cdef bint _is_left
    cdef op
    cpdef Element _call_(self, a, b)


cdef class InverseAction(Action):
    cdef Action _action
    cdef Morphism S_precomposition

cdef class PrecomposedAction(Action):
    cdef Action _action
    cdef Morphism left_precomposition
    cdef Morphism right_precomposition

cdef class ActionEndomorphism(Morphism):
    cdef Action _action
    cdef Element _g