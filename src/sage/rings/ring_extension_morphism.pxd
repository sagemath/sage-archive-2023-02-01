from sage.structure.element cimport Element
from sage.categories.map cimport Map
from sage.rings.morphism cimport RingMap
from sage.rings.ring_extension_element cimport RingExtensionElement


cdef _backend_morphism(f)
cpdef backend_morphism(f, forget=*)
cdef are_equal_morphisms(f, g)


cdef class RingExtensionHomomorphism(RingMap):
    cdef _backend
    cdef _im_gens
    cdef _base_map_construction

cdef class RingExtensionBackendIsomorphism(RingExtensionHomomorphism):
    pass

cdef class RingExtensionBackendReverseIsomorphism(RingExtensionHomomorphism):
    pass

cdef class MapVectorSpaceToRelativeField(Map):
    cdef _degree
    cdef list _basis
    cdef Map _f

cdef class MapRelativeFieldToVectorSpace(Map):
    cdef _degree
    cdef list _basis
    cdef _dimK
    cdef Map _iK
    cdef Map _jL
    cdef _matrix

    cdef list backend_coefficients(self, RingExtensionElement x)
