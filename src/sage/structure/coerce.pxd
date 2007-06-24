
from element cimport Element, RingElement, ModuleElement, CoercionModel
from parent cimport Parent

cdef class CoercionModel_simple(CoercionModel):
    pass

cdef class CoercionModel_cache_maps(CoercionModel_simple):
    cdef object _coercion_maps
    cdef coercion_maps_c(self, R, S)
    cdef discover_coercion_c(self, R, S)

    cdef object _action_maps
    cdef action_maps_c(self, R, S, op)
    cdef discover_action_c(self, R, S, op)
