
# It is important to keep this line here, basically to trick Pyrex.
# If you remove this line then other modules that cimport element
# from other directories will fail.

cimport sage.structure.sage_object
cimport sage.structure.parent

cimport sage_object
import  sage_object

cdef class Element(sage_object.SageObject):
    cdef sage.structure.parent.Parent _parent
    cdef int _cmp_c_impl(left, Element right) except -2
    cdef public _richcmp(self, right, int op)

cdef class ModuleElement(Element):
    cdef ModuleElement _add_c(self, ModuleElement right)             # do *NOT* override
    cdef ModuleElement _sub_c(self, ModuleElement right)             # do *NOT* override
    cdef ModuleElement _neg_c(self)                                  # do *NOT* override
    cdef ModuleElement _mul_left_scalar_c(self, RingElement left)    # do *NOT* override
    cdef ModuleElement _mul_right_scalar_c(self, RingElement right)  # do *NOT* override

    cdef ModuleElement _add_c_impl(self, ModuleElement right)        # OK to override
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)        # OK to override
    cdef ModuleElement _neg_c_impl(self)                             # OK to override
    cdef ModuleElement _mul_left_scalar_c_impl(self, RingElement left)    # OK to override
    cdef ModuleElement _mul_right_scalar_c_impl(self, RingElement right)  # OK to override

cdef class MonoidElement(Element):
    pass

cdef class MultiplicativeGroupElement(MonoidElement):
    pass

cdef class AdditiveGroupElement(ModuleElement):
    pass

cdef class RingElement(Element):
    cdef RingElement _add_c(self, RingElement right)          # do *NOT* override
    cdef RingElement _sub_c(self, RingElement right)          # do *NOT* override
    cdef RingElement _neg_c(self)                             # do *NOT* override
    cdef RingElement _mul_c(self, RingElement right)          # do *NOT* override
    cdef RingElement _div_c(self, RingElement right)          # do *NOT* override

    cdef RingElement _add_c_impl(self, RingElement right)     # do *NOT* override
    cdef RingElement _sub_c_impl(self, RingElement right)     # do *NOT* override
    cdef RingElement _neg_c_impl(self)                        # do *NOT* override
    cdef RingElement _mul_c_impl(self, RingElement right)     # do *NOT* override
    cdef RingElement _div_c_impl(self, RingElement right)     # do *NOT* override

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
    cdef AlgebraElement _mul_left_scalar_c(self, RingElement left)      # do not override
    cdef AlgebraElement _mul_right_scalar_c(self, RingElement right)    # do not override

    cdef AlgebraElement _mul_left_scalar_c_impl(self, RingElement left)     # ok to override
    cdef AlgebraElement _mul_right_scalar_c_impl(self, RingElement right)   # ok to override


cdef class CommutativeAlgebra(AlgebraElement):
    pass

cdef class InfinityElement(RingElement):
    pass
