from c_graph cimport CGraph
from static_sparse_graph cimport short_digraph, ushort
from c_graph import CGraphBackend

include 'sage/ext/stdsage.pxi'

cdef class StaticSparseCGraph(CGraph):
    cdef short_digraph g
    cdef short_digraph g_rev
    cdef bint _directed

    cpdef bint has_vertex(self, int n)
    cdef int add_vertex_unsafe(self, int k)
    cdef int del_vertex_unsafe(self, int v)
    cpdef list verts(self)
    cdef int has_arc_unsafe(self, int u, int v)
    cpdef bint has_arc(self, int u, int v)
    cdef int out_neighbors_unsafe(self, int u, int *neighbors, int size) except? -2
    cpdef list out_neighbors(self, int u)
    cpdef int out_degree(self, int u)
    cdef int in_neighbors_unsafe(self, int u, int *neighbors, int size) except? -2
    cpdef list in_neighbors(self, int u)
    cpdef int in_degree(self, int u)
