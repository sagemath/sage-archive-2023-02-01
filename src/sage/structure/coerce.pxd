
from element cimport Element, RingElement, ModuleElement, CoercionModel

from parent cimport Parent
from parent_base cimport ParentWithBase
from sage.categories.action cimport Action
from sage.categories.morphism cimport Morphism

from coerce_dict cimport TripleDict, TripleDictIter

cdef class CoercionModel_original(CoercionModel):
    pass

cdef class CoercionModel_cache_maps(CoercionModel_original):
    # This MUST be a mapping to tuples, where each
    # tuple contains at least two elements that are either
    # None or of type Morphism.
    cdef TripleDict _coercion_maps

    # This MUST be a mapping to actions.
    cdef TripleDict _action_maps

    cdef coercion_maps_c(self, R, S)
    cdef discover_coercion_c(self, R, S)
    cpdef verify_coercion_maps(self, R, S, homs, bint fix=*)

    cdef get_action_c(self, R, S, op)
    cdef discover_action_c(self, R, S, op)
    cpdef verify_action(self, action, R, S, op, bint fix=*)

cdef class CoercionModel_profile(CoercionModel_cache_maps):
    cdef object profiling_info
    cdef object timer
    cdef void _log_time(self, xp, yp, op, time, data)


cdef class LeftModuleAction(Action):
    cdef Morphism connecting
    cdef extended_base

cdef class RightModuleAction(Action):
    cdef Morphism connecting
    cdef extended_base
    cdef bint is_inplace

cdef class PyScalarAction(Action):
    cdef Action _action