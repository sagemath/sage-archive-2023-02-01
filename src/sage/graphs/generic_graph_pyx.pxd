
from sage.structure.sage_object cimport SageObject

cdef run_spring(int, int, double*, int*, int, bint)

cdef class GenericGraph_pyx(SageObject):
    pass

cdef inline bint vectors_equal(int n, int *a, int *b)
cdef inline bint vectors_inferior(int n, int *a, int *b)
