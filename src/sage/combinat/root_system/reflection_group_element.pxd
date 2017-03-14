from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement

cdef class ComplexReflectionGroupElement(PermutationGroupElement):
    pass

cdef class RealReflectionGroupElement(ComplexReflectionGroupElement):
    cpdef bint has_left_descent(self, i)
    cpdef bint has_descent(self, i, side=*, positive=*)
