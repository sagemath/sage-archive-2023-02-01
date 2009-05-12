from sage.structure.element cimport MultiplicativeGroupElement, MonoidElement, Element


cdef class PermutationGroupElement(MultiplicativeGroupElement):
    cdef int* perm
    cdef int n
    cdef int perm_buf[15] # to avoid malloc for small elements
    cdef __gap
    cdef Element _gap_element
    cdef __tuple
    cdef PermutationGroupElement _new_c(self)
    cpdef list(self)
    cdef public __custom_name
