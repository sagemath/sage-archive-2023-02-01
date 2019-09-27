from sage.rings.ring_extension cimport RingExtension_class


cpdef backend_parent(R)
cpdef from_backend_parent(R, RingExtension_class E)

cpdef backend_element(x)
cpdef from_backend_element(x, RingExtension_class E)

cdef _backend_morphism(f)
cpdef backend_morphism(f, forget=*)
cpdef from_backend_morphism(f, RingExtension_class E)

