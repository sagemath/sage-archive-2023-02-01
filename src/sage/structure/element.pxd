
# It is important to keep this line here, basically to trick Pyrex.
# If you remove this line then other modules that cimport element
# from other directories will fail.

cimport sage.structure.sage_object
from sage.structure.parent cimport Parent

cimport sage_object
import  sage_object

cdef class Element(sage_object.SageObject):
    cdef Parent _parent
    cdef int _cmp_c_impl(left, Element right) except -2
    cdef public _richcmp(self, right, int op)
    cdef _set_parent_c(self, Parent parent)

    cdef base_extend_c(self, Parent R)       # do *NOT* override, but OK to call directly
    cdef base_extend_c_impl(self, Parent R)  # OK to override, but do NOT call

cdef class ModuleElement(Element)       # forward declaration
cdef class RingElement(ModuleElement)   # forward declaration

cdef class ModuleElement(Element):
    cdef ModuleElement _add_c(self, ModuleElement right)             # do *NOT* override, but OK to call directly
    cdef ModuleElement _sub_c(self, ModuleElement right)             # do *NOT* override, but OK to call directly
    cdef ModuleElement _neg_c(self)                                  # do *NOT* override, but OK to call directly
    # self._rmul_c(x) is x * self
    cdef ModuleElement _rmul_c(self, RingElement left)               # do *NOT* override, but OK to call directly
    # self._lmul_c(x) is self * x, to abide with Python conventions.
    cdef ModuleElement _lmul_c(self, RingElement right)              # do *NOT* override, but OK to call directly

    cdef ModuleElement _add_c_impl(self, ModuleElement right)        # OK to override, but do NOT call
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)        # OK to override, but do NOT call
    cdef ModuleElement _neg_c_impl(self)                             # OK to override, but do *NOT* call directly
    cdef ModuleElement _lmul_c_impl(self, RingElement right)          # OK to override, but do *NOT* call directly
    cdef ModuleElement _rmul_c_impl(self, RingElement left)         # OK to override, but do *NOT* call directly

    # Coerce x to the base ring of self and return the result.
    cdef RingElement coerce_to_base_ring(self, x)

cdef class MonoidElement(Element):
    cdef MonoidElement _mul_c(self, MonoidElement right)             # do *NOT* override, but OK to call directly
    cdef MonoidElement _mul_c_impl(self, MonoidElement right)        # OK to override, but do *NOT* call directly

cdef class MultiplicativeGroupElement(MonoidElement):
    pass

cdef class AdditiveGroupElement(ModuleElement):
    pass

cdef class RingElement(ModuleElement):
    cdef RingElement _mul_c(self, RingElement right)          # do *NOT* override, but OK to call directly
    cdef RingElement _div_c(self, RingElement right)          # do *NOT* override, but OK to call directly

    cdef RingElement _mul_c_impl(self, RingElement right)     # OK to override, but do *NOT* call directly
    cdef RingElement _div_c_impl(self, RingElement right)     # OK to override, but do *NOT* call directly

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
