from hom_finite_field cimport SectionFiniteFieldHomomorphism_generic
from hom_finite_field cimport FiniteFieldHomomorphism_generic
from hom_finite_field cimport FrobeniusEndomorphism_finite_field


cdef class SectionFiniteFieldHomomorphism_prime(SectionFiniteFieldHomomorphism_generic):
    pass


cdef class FiniteFieldHomomorphism_prime(FiniteFieldHomomorphism_generic):
    pass


cdef class FrobeniusEndomorphism_prime(FrobeniusEndomorphism_finite_field):
    pass
