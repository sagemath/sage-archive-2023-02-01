from .c_graph cimport CGraph, CGraphBackend
from .static_sparse_graph cimport short_digraph, ushort
from libc.stdint cimport uint64_t, uint32_t, INT32_MAX, UINT32_MAX
from sage.data_structures.bitset cimport *

cdef class StaticSparseCGraph(CGraph):
    cdef short_digraph g
    cdef short_digraph g_rev
    cdef bint _directed
    cdef int * number_of_loops

    cpdef int out_degree(self, int u) except -1
    cpdef int in_degree(self, int u) except -1

cdef class StaticSparseBackend(CGraphBackend):
    cdef int _order
    cdef bint _multiedges
    cdef list _vertex_to_labels
    cdef dict _vertex_to_int
    cdef bitset_t _seen

cdef uint32_t simple_BFS(short_digraph g,
                         uint32_t source,
                         uint32_t *distances,
                         uint32_t *predecessors,
                         uint32_t *waiting_list,
                         bitset_t seen)