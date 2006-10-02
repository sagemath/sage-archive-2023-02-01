import  sage.ext.sage_object
cimport sage.ext.sage_object

cdef class Element(sage.ext.sage_object.SageObject):
    cdef sage.ext.sage_object.SageObject _parent

cdef class ModuleElement(Element):
    pass

cdef class MonoidElement(Element):
    pass

cdef class MultiplicativeGroupElement(MonoidElement):
    pass

cdef class AdditiveGroupElement(ModuleElement):
    pass

cdef class RingElement(Element):
    pass

cdef class CommutativeRingElement(RingElement):
    pass

cdef class IntegralDomainElement(CommutativeRingElement):
    pass

cdef class DedekindDomainElement(IntegralDomainElement):
    pass

cdef class PrincipalIdealDomainElement(DedekindDomainElement):
    pass

cdef class EuclideanDomainElement(PrincipalIdealDomainElement):
    pass

cdef class FieldElement(CommutativeRingElement):
    pass

cdef class AlgebraElement(RingElement):
    pass

cdef class CommutativeAlgebra(AlgebraElement):
    pass

cdef class InfinityElement(RingElement):
    pass
