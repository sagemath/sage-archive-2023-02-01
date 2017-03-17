from sage.structure.element cimport Element

cdef class IndexedFreeModuleElement(Element):
    cdef public dict _monomial_coefficients
    cdef long _hash
    cdef bint _hash_set

    cpdef dict monomial_coefficients(self, bint copy=*)
    cpdef _coefficient_fast(self, m)

