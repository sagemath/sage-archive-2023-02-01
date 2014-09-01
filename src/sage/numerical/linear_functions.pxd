from sage.structure.parent cimport Parent
from sage.structure.element cimport ModuleElement, RingElement, Element

cpdef is_LinearFunction(x)

cdef class LinearFunctionsParent_class(Parent):
    cpdef _element_constructor_(self, x)
    cpdef _coerce_map_from_(self, R)
    cdef object _multiplication_symbol
    cpdef object _get_multiplication_symbol(self)

cdef class LinearFunction(ModuleElement):
    cdef dict _f
    cpdef iteritems(self)
    cpdef ModuleElement _add_(self, ModuleElement b)
    cpdef ModuleElement _sub_(self, ModuleElement b)
    cpdef ModuleElement _lmul_(self, RingElement b)
    cpdef ModuleElement _rmul_(self, RingElement b)
    cpdef ModuleElement _neg_(self)
    cpdef _acted_upon_(self, x, bint self_on_left)
    cdef _richcmp(left, right, int op)
    cpdef is_zero(self)
    cpdef equals(LinearFunction left, LinearFunction right)
    cdef int _cmp_c_impl(left, Element right) except -2

cdef class LinearConstraintsParent_class(Parent):
    cdef LinearFunctionsParent_class _LF
    cpdef _element_constructor_(self, left, right=?, equality=?)
    cpdef _coerce_map_from_(self, R)

cdef class LinearConstraint(Element):
    cdef bint equality
    cdef list constraints
    cdef _richcmp(left, right, int op)
    cdef LinearConstraint _chained_comparator_hack_part1(LinearConstraint left, LinearConstraint right)
    cdef _chained_comparator_hack_part2(self)
    cpdef equals(LinearConstraint left, LinearConstraint right)

cdef LinearConstraint _chained_comparator_hack_search
cdef LinearConstraint _chained_comparator_hack_replace
