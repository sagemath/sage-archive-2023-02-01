
from element cimport Element, RingElement, ModuleElement
from parent cimport Parent

from sage.categories.action cimport Action
from sage.categories.morphism cimport Morphism


cdef class ModuleAction(Action):
    cdef Morphism connecting
    cdef extended_base

cdef class LeftModuleAction(ModuleAction):
    pass

cdef class RightModuleAction(ModuleAction):
    cdef public bint is_inplace

cdef class PyScalarAction(Action):
    cdef Action _action

cdef class IntegerMulAction(Action):
    pass