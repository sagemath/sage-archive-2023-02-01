cimport numpy as np
import numpy as np
from sage.rings.number_field.number_field_base cimport NumberField
from sage.rings.number_field.number_field_element cimport NumberFieldElement_absolute

cdef class KSHandler:
    cdef list obj_cache
    cdef np.ndarray ks_dat
    cdef NumberField field
    cdef public shm

    cdef bint contains(self, int idx)
    cdef NumberFieldElement_absolute get(self, int idx)
    cdef setitem(self, int idx, rhs)
    cpdef update(self, list eqns)

cdef class FvarsHandler:
    cdef dict sext_to_idx, obj_cache
    cdef unsigned int ngens
    cdef fvars_t
    cdef NumberField field
    cdef np.ndarray fvars
    cdef public shm

    cdef clear_modified(self)
