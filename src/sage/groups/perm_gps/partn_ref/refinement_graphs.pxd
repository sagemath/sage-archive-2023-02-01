
#*****************************************************************************
#      Copyright (C) 2006 - 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include '../../../ext/cdefs.pxi'
include '../../../ext/stdsage.pxi'

from sage.graphs.base.c_graph cimport CGraph
from sage.graphs.base.sparse_graph cimport SparseGraph
from sage.graphs.base.dense_graph cimport DenseGraph
from sage.rings.integer cimport Integer
from sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label cimport PartitionStack, PS_is_discrete, PS_move_min_to_front, PS_num_cells, get_aut_gp_and_can_lab, aut_gp_and_can_lab_return

cdef class GraphStruct:
    cdef CGraph G
    cdef bint directed
    cdef bint use_indicator
    cdef int *scratch # length 3n+1

cdef int refine_by_degree(PartitionStack *, object, int *, int)
cdef int compare_graphs(int *, int *, object)
cdef bint all_children_are_equivalent(PartitionStack *, object)
cdef inline int degree(PartitionStack *, CGraph, int, int, bint)
cdef inline int sort_by_function(PartitionStack *, int, int *)



