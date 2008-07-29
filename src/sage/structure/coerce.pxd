
from element cimport Element, RingElement, ModuleElement, CoercionModel

from parent cimport Parent
from sage.categories.action cimport Action
from sage.categories.morphism cimport Morphism

from coerce_dict cimport TripleDict, TripleDictIter

cdef class CoercionModel_cache_maps(CoercionModel):
    # This MUST be a mapping to tuples, where each
    # tuple contains at least two elements that are either
    # None or of type Morphism.
    cdef TripleDict _coercion_maps

    # This MUST be a mapping to actions.
    cdef TripleDict _action_maps

    cpdef coercion_maps(self, R, S)
    cpdef discover_coercion(self, R, S)
    cpdef verify_coercion_maps(self, R, S, homs, bint fix=*)
    cpdef verify_action(self, action, R, S, op, bint fix=*)

    cpdef get_action(self, xp, yp, op)
    cpdef discover_action(self, R, S, op)

    cdef readonly _exception_stack
    cdef bint _exceptions_cleared
