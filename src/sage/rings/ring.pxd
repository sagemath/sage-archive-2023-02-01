cimport sage.structure.parent_gens

cdef class Ring(sage.structure.parent_gens.ParentWithGens):
    pass

cdef class CommutativeRing(Ring):
    cdef public object __ideal_monoid

cdef class IntegralDomain(CommutativeRing):
    cdef public object __fraction_field
    pass

cdef class DedekindDomain(IntegralDomain):
    pass


cdef class PrincipalIdealDomain(IntegralDomain):
    pass

cdef class EuclideanDomain(PrincipalIdealDomain):
    pass

cdef class Field(PrincipalIdealDomain):
    pass

cdef class FiniteField(Field):
    cdef public object __multiplicative_generator
    cdef public object __polynomial_ring
    cdef public object __vector_space

cdef class Algebra(Ring):
    pass

cdef class CommutativeAlgebra(CommutativeRing):
    pass
