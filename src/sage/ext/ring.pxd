import gens
cimport gens

cdef class Ring(gens.Generators):
    pass

cdef class CommutativeRing(Ring):
    pass

cdef class IntegralDomain(CommutativeRing):
#    cdef object _fraction_field
    pass

cdef class DedekindDomain(IntegralDomain):
    pass


cdef class PrincipalIdealDomain(IntegralDomain):
    pass

cdef class EuclideanDomain(PrincipalIdealDomain):
    pass

cdef class Field(PrincipalIdealDomain):
    pass

cdef class Algebra(Ring):
    pass

cdef class CommutativeAlgebra(Algebra):
    pass
