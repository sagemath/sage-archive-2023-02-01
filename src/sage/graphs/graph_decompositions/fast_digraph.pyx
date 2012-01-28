r"""
Compact structure for fast operations on less than 32 vertices

This module implements a digraph structure meant to be used in Cython in
**highly enumerative** algorithms. It can store graphs on less than
``sizeof(int)`` vertices and perform several basic operations **quickly**
(add/remove arcs, count the out-neighborhood of a set of vertices or return its
cardinality).

**Sets and integers :**

In the following code, sets are represented as integers, where the ith bit is
set if element i belongs to the set.
"""

include '../../ext/stdsage.pxi'
include '../../ext/cdefs.pxi'
include '../../ext/interrupt.pxi'

from libc.stdint cimport uint8_t

cdef class FastDigraph:
    cdef uint8_t n
    cdef int * graph

    def __cinit__(self, D):
        if D.order() > 8*sizeof(int):
            raise Exception("Too many vertices. This structure can only encode digraphs on at most sizeof(int) vertices.")

        self.n = D.order()
        self.graph = NULL

        self.graph = <int *> sage_malloc(self.n*sizeof(int))

        memset(self.graph, 0, self.n * sizeof(int))

        cdef int i, j
        cdef int tmp

        # When the vertices are not consecutive integers
        cdef dict vertices_to_int = {}
        for i,v in enumerate(D.vertices()):
            vertices_to_int[v] = i

        for u in D:
            tmp = 0
            for v in D.neighbors_out(u):
                tmp |= 1 << vertices_to_int[v]
            self.graph[vertices_to_int[u]] = tmp

    def __dealloc__(self):
        if self.graph != NULL:
            sage_free(self.graph)

    def print_adjacency_matrix(self):
        cdef int i,j
        for 0<= i<self.n:
            for 0<= j <self.n:
                print ((self.graph[i]>>j)&1),
            print ""

cdef inline int compute_out_neighborhood_cardinality(FastDigraph g, int S):
    r"""
    Returns the cardinality of `N^+(S)\S`.

    INPUT:

    - ``g`` a FastDigraph
    - S (integer) an integer describing the set
    """
    cdef int i
    cdef int tmp = 0
    for 0<= i<g.n:
        tmp |= g.graph[i] * ((S >> i)&1)

    tmp &= (~S)
    return popcount(tmp)


cdef extern from *:
    int __builtin_popcount(unsigned int)

cdef inline int popcount(int i):
    #return __builtin_popcount(i)
    i = i - ((i >> 1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return ((i + (i >> 4) & 0x0F0F0F0F) * 0x01010101) >> 24;
