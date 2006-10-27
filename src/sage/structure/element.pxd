
# It is important to keep this line here, basically to trick Pyrex.
# If you remove this line then other modules that cimport element
# from other directories will fail.

cimport sage.structure.sage_object


cimport sage_object
import  sage_object

cdef class Element(sage_object.SageObject):
    #cdef sage_object.SageObject _parent
    cdef object _parent

cdef class ModuleElement(Element):
    cdef ModuleElement _add_sibling_cdef(self, ModuleElement right)

cdef class MonoidElement(Element):
    pass

cdef class MultiplicativeGroupElement(MonoidElement):
    pass

cdef class AdditiveGroupElement(ModuleElement):
    pass

cdef class RingElement(Element):
    cdef RingElement _add_sibling_cdef(self, RingElement right)

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
