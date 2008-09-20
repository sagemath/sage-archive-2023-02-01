
#*****************************************************************************
#      Copyright (C) 2006 - 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include '../../../ext/cdefs.pxi'
include '../../../ext/stdsage.pxi'
include 'data_structures_pxd.pxi' # includes bitsets

from sage.rings.integer cimport Integer
from sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label cimport get_aut_gp_and_can_lab, aut_gp_and_can_lab_return
from double_coset cimport double_coset

cdef class BinaryCodeStruct:
    cdef bitset_s *alpha_is_wd # single bitset of length nwords + degree
    cdef int degree
    cdef int nwords
    cdef bint first_time
    cdef PartitionStack *word_ps
    cdef int *alpha # length nwords + degree
    cdef int *scratch # length 3*nwords + 3*degree + 2
    cdef aut_gp_and_can_lab_return *output
    cdef int (*ith_word)(BinaryCodeStruct self, int, bitset_s *)

cdef class LinearBinaryCodeStruct(BinaryCodeStruct):
    cdef bitset_s *basis
    cdef bitset_s *scratch_bitsets # length 2*dimension + 2
    cdef int dimension
cdef int ith_word_linear(BinaryCodeStruct, int, bitset_s *)

cdef class NonlinearBinaryCodeStruct(BinaryCodeStruct):
    cdef bitset_s *words
    cdef bitset_s *scratch_bitsets # length 4*nwords + 1
cdef int ith_word_nonlinear(BinaryCodeStruct, int, bitset_s *)

cdef int refine_by_bip_degree(PartitionStack *, object, int *, int)
cdef int compare_linear_codes(int *, int *, object, object)
cdef int compare_nonlinear_codes(int *, int *, object, object)
cdef bint all_children_are_equivalent(PartitionStack *, object)
cdef inline int word_degree(PartitionStack *, BinaryCodeStruct, int, int, PartitionStack *)
cdef inline int col_degree(PartitionStack *, BinaryCodeStruct, int, int, PartitionStack *)
cdef inline int sort_by_function(PartitionStack *, int, int *, int *, int *, int)



