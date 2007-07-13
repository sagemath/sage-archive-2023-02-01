from sage.structure.element cimport Element
from morphism cimport Morphism
from functor cimport Functor

cdef class Action(Functor):
    cdef _G
    cdef _S
    cdef bint _is_left
    cdef _op
    cdef Element _call_c(self, a, b)
    cdef Element _call_c_impl(self, Element a, Element b)


cdef class InverseAction(Action):
    cdef Action _action

cdef class PrecomposedAction(Action):
    cdef Action _action
    cdef Morphism _left_precomposition
    cdef Morphism _right_precomposition

cdef class ActionEndomorphism(Morphism):
    cdef Action _action
    cdef Element _g