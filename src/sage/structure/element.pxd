
# It is important to keep this line here, basically to trick Pyrex.
# If you remove this line then other modules that cimport element
# from other directories will fail.

cimport sage.structure.sage_object
cimport sage.structure.parent

from sage.structure.structure.element import Element

cimport sage_object
import  sage_object

cdef class Element(sage_object.SageObject):
    cdef sage.structure.parent.Parent _parent
    #cdef object _parent
    cdef _rich_to_bool(self, int op, int n)
    cdef int _cmp_c_impl(left, Element right) except -2
    cdef public _richcmp(self, right, int op)

cdef class ModuleElement(Element):
    cdef ModuleElement _add_c(self, ModuleElement right)
    cdef ModuleElement _add_c_impl(self, ModuleElement right)
    cdef ModuleElement _sub_c(self, ModuleElement right)
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)
    cdef ModuleElement _neg_c(self)
    cdef ModuleElement _neg_c_impl(self)

cdef class MonoidElement(Element):
    pass

cdef class MultiplicativeGroupElement(MonoidElement):
    pass

cdef class AdditiveGroupElement(ModuleElement):
    pass

cdef class RingElement(Element):
    cdef RingElement _add_c(self, RingElement right)
    cdef RingElement _add_c_impl(self, RingElement right)
    cdef RingElement _sub_c(self, RingElement right)
    cdef RingElement _sub_c_impl(self, RingElement right)
    cdef RingElement _neg_c(self)
    cdef RingElement _neg_c_impl(self)

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

cdef class FiniteFieldElement(FieldElement):
    pass

cdef class AlgebraElement(RingElement):
    pass

cdef class CommutativeAlgebra(AlgebraElement):
    pass

cdef class InfinityElement(RingElement):
    pass
