from sage.structure.element cimport Element, RingElement

cdef class IndexedFreeModuleElement(Element):
    cdef public dict _monomial_coefficients
    cdef long _hash
    cdef bint _hash_set

    cpdef dict monomial_coefficients(self, bint copy=*)
    cpdef _coefficient_fast(self, m)
    cpdef _lmul_(self, Element right)
    cpdef _rmul_(self, Element left)

