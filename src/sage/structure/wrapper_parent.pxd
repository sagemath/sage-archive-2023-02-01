from sage.structure.parent cimport Parent
from sage.categories.morphism cimport Morphism
from sage.structure.element cimport AlgebraElement, Element, ModuleElement, RingElement

cdef class StealMorphism(Morphism):
    cdef Element _call_c_impl(self, Element x)

cdef class WrapperParent_model0(Parent):
    cdef Parent R

cdef class WrapperParent_model1(Parent):
    cdef Parent R

cdef class WrapperElement(AlgebraElement):
    cdef Element val
    cdef _richcmp_c_impl(left, Element right, int op)
    cdef int _cmp_c_impl(left, Element right) except -2
    cdef base_extend_c_impl(self, Parent R)
    cdef ModuleElement _add_c_impl(self, ModuleElement right)
    cdef ModuleElement _sub_c_impl(self, ModuleElement right)
    cdef ModuleElement _neg_c_impl(self)
    cdef ModuleElement _lmul_c_impl(self, RingElement right)
    cdef ModuleElement _rmul_c_impl(self, RingElement left)
    cdef RingElement coerce_to_base_ring(self, x)
    cdef ModuleElement _lmul_nonscalar_c_impl(left, right)
    cdef ModuleElement _rmul_nonscalar_c_impl(right, left)
    cdef RingElement _mul_c_impl(self, RingElement right)
    cdef RingElement _div_c_impl(self, RingElement right)
