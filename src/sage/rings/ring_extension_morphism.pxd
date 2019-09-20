from sage.structure.element cimport Element
from sage.rings.morphism cimport RingHomomorphism

cpdef _backend_morphism(f)

cdef class RingExtensionHomomorphism(RingHomomorphism):
    cdef _backend
    cdef _im_gens
    cdef _base_map_construction
    cpdef Element _call_(self, x)

cdef class RingExtensionBackendIsomorphism(RingExtensionHomomorphism):
    pass

cdef class RingExtensionBackendReverseIsomorphism(RingExtensionHomomorphism):
    pass
