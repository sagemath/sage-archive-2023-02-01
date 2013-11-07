from sage.structure.element cimport Element, MonoidElement, MultiplicativeGroupElement
from sage.rings.finite_rings.integer_mod cimport IntegerMod_int

cdef class SemimonomialTransformation(MultiplicativeGroupElement):
    cdef tuple v
    cdef object perm, alpha

    cpdef MonoidElement _mul_(left, MonoidElement _right)
    cdef int _cmp_c_impl(left, Element _right) except -2
    cdef _new_c(self)