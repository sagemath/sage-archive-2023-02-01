
from element cimport Element, RingElement, ModuleElement, CoercionModel

from parent cimport Parent
from parent_base cimport ParentWithBase
from sage.categories.action cimport Action
from sage.categories.morphism cimport Morphism

cdef class CoercionModel_original(CoercionModel):
    pass

cdef class CoercionModel_cache_maps(CoercionModel_original):
    # This MUST be a mapping to tuples, where each
    # tuple contains at least two elements that are either
    # None or of type Morphism.
    cdef object _coercion_maps

    # This MUST be a mapping to actions.
    cdef object _action_maps

    cdef coercion_maps_c(self, R, S)
    cdef discover_coercion_c(self, R, S)

    cdef action_maps_c(self, R, S, op)
    cdef discover_action_c(self, R, S, op)

cdef class LeftModuleAction(Action):
    pass

cdef class RightModuleAction(Action):
    pass