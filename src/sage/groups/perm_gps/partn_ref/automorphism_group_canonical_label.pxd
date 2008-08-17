
#*****************************************************************************
#      Copyright (C) 2006 - 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include '../../../ext/cdefs.pxi'
include '../../../ext/stdsage.pxi'
include '../../../misc/bitset_pxd.pxi'

from sage.rings.integer cimport Integer

cdef struct OrbitPartition:
    # Disjoint set data structure for representing the orbits of the generators
    # so far found. Also keeps track of the minimum elements of the cells and
    # their sizes.
    int degree
    int *parent
    int *rank
    int *mcr # minimum cell representatives - only valid at the root of a cell
    int *size # also only valid at the root of a cell

cdef inline OrbitPartition *OP_new(int)
cdef inline int OP_dealloc(OrbitPartition *)
cdef inline int OP_find(OrbitPartition *, int)
cdef inline int OP_join(OrbitPartition *, int, int)
cdef inline int OP_merge_list_perm(OrbitPartition *, int *)

cdef struct PartitionStack:
    # Representation of a node of the search tree. A sequence of partitions of
    # length depth + 1, each of which is finer than the last. Partition k is
    # represented as PS.entries in order, broken up immediately after each
    # entry of levels which is at most k.
    int *entries
    int *levels
    int depth
    int degree

cdef inline PartitionStack *PS_new(int, bint)
cdef inline PartitionStack *PS_copy(PartitionStack *)
cdef inline int PS_copy_from_to(PartitionStack *, PartitionStack *)
cdef inline PartitionStack *PS_from_list(object)
cdef inline int PS_dealloc(PartitionStack *)
cdef inline int PS_print(PartitionStack *)
cdef inline int PS_print_partition(PartitionStack *, int)
cdef inline bint PS_is_discrete(PartitionStack *)
cdef inline int PS_num_cells(PartitionStack *)
cdef inline int PS_move_min_to_front(PartitionStack *, int, int)
cdef inline bint PS_is_mcr(PartitionStack *, int)
cdef inline bint PS_is_fixed(PartitionStack *, int)
cdef inline int PS_clear(PartitionStack *)
cdef inline int PS_split_point(PartitionStack *, int)
cdef inline int PS_first_smallest(PartitionStack *, bitset_t)
cdef inline int PS_get_perm_from(PartitionStack *, PartitionStack *, int *)

cdef inline int split_point_and_refine(PartitionStack *, int, object,
    int (*)(PartitionStack *, object, int *, int), int *)

cdef struct aut_gp_and_can_lab_return:
    int *generators
    int num_gens
    int *relabeling
    int *base
    int base_size
    mpz_t order

cdef aut_gp_and_can_lab_return *get_aut_gp_and_can_lab( object, int **, int,
    bint (*)(PartitionStack *, object),
    int (*)(PartitionStack *, object, int *, int),
    int (*)(int *, int *, object), bint, bint, bint)



