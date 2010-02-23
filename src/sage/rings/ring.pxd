from sage.structure.parent_gens cimport ParentWithGens

cdef class Ring(ParentWithGens):
    cdef public object _zero_element
    cdef public object _one_element
    cdef public object _zero_ideal
    cdef public object _unit_ideal
    cdef _an_element_c_impl(self)

cdef class CommutativeRing(Ring):
    cdef public object __fraction_field
    cdef public object __ideal_monoid

cdef class IntegralDomain(CommutativeRing):
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

cdef class CommutativeAlgebra(CommutativeRing):
    pass
