from c_graph cimport CGraph
from static_sparse_graph cimport short_digraph, ushort
from c_graph import CGraphBackend

cdef class StaticSparseCGraph(CGraph):
    cdef short_digraph g
    cdef short_digraph g_rev
    cdef bint _directed

    cpdef int out_degree(self, int u) except -1
    cpdef int in_degree(self, int u) except -1
