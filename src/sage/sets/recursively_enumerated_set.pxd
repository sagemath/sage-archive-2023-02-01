
cimport sage.structure.parent

cdef class RecursivelyEnumeratedSet_generic(sage.structure.parent.Parent):
    cdef readonly _seeds
    cdef public successors
    cdef readonly str _enumeration
    cdef readonly _max_depth
    cdef readonly _graded_component
    cdef readonly _graded_component_it

    cpdef seeds(self)
    cpdef graded_component(self, depth)

cdef class RecursivelyEnumeratedSet_symmetric(RecursivelyEnumeratedSet_generic):
    cdef set _get_next_graded_component(self, set A, set B)

cdef class RecursivelyEnumeratedSet_graded(RecursivelyEnumeratedSet_generic):
    cdef set _get_next_graded_component(self, set B)

