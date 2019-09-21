from sage.structure.element cimport Element
from sage.categories.map cimport Map
from sage.rings.morphism cimport RingHomomorphism

cpdef _backend_morphism(f)

cdef class RingExtensionHomomorphism(RingHomomorphism):
    cdef _backend
    cdef _im_gens
    cdef _base_map_construction

cdef class RingExtensionBackendIsomorphism(RingExtensionHomomorphism):
    pass

cdef class RingExtensionBackendReverseIsomorphism(RingExtensionHomomorphism):
    pass

#cdef class MapVectorSpaceToRelativeField(Map):
#    pass

#cdef class MapRelativeFieldToVectorSpace(Map):
#    cdef backend_coefficients(self, x)


