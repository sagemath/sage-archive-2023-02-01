
#*******************************************************************************
#        Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

from c_graph cimport CGraph
include '../../ext/stdsage.pxi'

cdef struct SparseGraphLLNode:
    int label
    SparseGraphLLNode *next

cdef struct SparseGraphBTNode:
    int vertex
    SparseGraphLLNode *labels
    SparseGraphBTNode *left
    SparseGraphBTNode *right

cdef class SparseGraph(CGraph):
    # Values inherited from CGraph:
    # cdef int num_verts
    # cdef int num_arcs
    # cdef int *in_degrees
    # cdef int *out_degrees
    # Values specific to SparseGraph:
    cdef int hash_length
    cdef int hash_mask
    cdef SparseGraphBTNode **vertices

    # Method declarations inherited from CGraph:
    # cdef int add_arc_unsafe(self, int, int)
    # cdef int has_arc_unsafe(self, int, int)
    # cdef int del_arc_unsafe(self, int, int)
    # cdef int out_neighbors_unsafe(self, int, int *, int)
    # cdef int in_neighbors_unsafe(self, int, int *, int)
    # Methods specific to SparseGraph only (labels!):
    cdef int add_arc_label_unsafe(self, int, int, int)
    cdef int arc_label_unsafe(self, int, int)
    cdef int all_arcs_unsafe(self, int, int, int *, int)
    cdef int del_arc_label_unsafe(self, int, int, int)
    cdef int has_arc_label_unsafe(self, int, int, int)
    cdef int del_vertex_unsafe(self, int)
    cdef int add_vertices_unsafe(self, int)





