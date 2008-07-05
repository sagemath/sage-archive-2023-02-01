import sage.rings.ring
cimport sage.rings.ring

from sage.structure.parent cimport Parent

cdef class MPolynomialRing_generic(sage.rings.ring.CommutativeRing):
    cdef object __ngens
    cdef object __term_order
    cdef object _has_singular
    cdef object __magma
    cdef object __magma_gens

    cdef _coerce_c_impl(self, x)
    cdef int _cmp_c_impl(left, Parent right) except -2
