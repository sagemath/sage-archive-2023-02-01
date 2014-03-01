from sage.rings.morphism cimport RingHomomorphism, RingHomomorphism_from_base

cdef class PolynomialRingHomomorphism_from_base(RingHomomorphism_from_base):
    pass

cdef class PolynomialRingMorphism(RingHomomorphism):
    cdef object base_morphism
    cdef object im_gen
