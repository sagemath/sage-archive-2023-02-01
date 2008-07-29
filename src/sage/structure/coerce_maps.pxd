from sage.categories.action cimport Action
from sage.categories.morphism cimport Morphism

cdef class DefaultConvertMorphism(Morphism):
    cdef public bint _force_use
    cdef public bint _is_coercion

cdef class DefaultConvertMorphism_unique(DefaultConvertMorphism):
    pass

cdef class NamedConvertMorphism(Morphism):
    cdef readonly method_name
    cdef public bint _force_use

cdef class TryMorphism(Morphism):
    cdef Morphism _morphism_p
    cdef Morphism _morphism_b
    cdef _error_types
