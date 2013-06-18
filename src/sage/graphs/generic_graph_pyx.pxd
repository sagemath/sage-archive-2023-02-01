
from sage.structure.sage_object cimport SageObject
from sage.graphs.base.dense_graph cimport DenseGraph

cdef run_spring(int, int, double*, int*, int, bint)

cdef class GenericGraph_pyx(SageObject):
    pass



cdef class SubgraphSearch:
    cdef int ng
    cdef int nh
    cdef (bint) (*is_admissible) (int, int *, int *)
    cdef DenseGraph g
    cdef DenseGraph h
    cdef int *busy
    cdef int *stack
    cdef int active
    cdef int *vertices
    cdef int **line_h_out
    cdef int **line_h_in
    cdef list g_vertices
    cdef int i
    cdef bint directed


cdef inline bint vectors_equal(int n, int *a, int *b)
cdef inline bint vectors_inferior(int n, int *a, int *b)
