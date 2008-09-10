
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
from sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label cimport PartitionStack, PS_new, PS_dealloc, PS_copy_from_to, get_aut_gp_and_can_lab, aut_gp_and_can_lab_return
from sage.groups.perm_gps.partn_ref.refinement_binary cimport NonlinearBinaryCodeStruct, refine_by_bip_degree
from sage.groups.perm_gps.partn_ref.refinement_binary cimport all_children_are_equivalent as all_binary_children_are_equivalent

cdef class MatrixStruct:
    cdef list symbol_structs
    cdef object matrix
    cdef int degree
    cdef int nwords
    cdef list symbols
    cdef int nsymbols
    cdef PartitionStack *temp_col_ps
    cdef aut_gp_and_can_lab_return *output

cdef int refine_matrix(PartitionStack *, object, int *, int)
cdef int compare_matrices(int *, int *, object)
cdef bint all_matrix_children_are_equivalent(PartitionStack *, object)

