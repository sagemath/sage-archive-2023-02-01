
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
    int (*)(int *, int *, object, object), bint, bint, bint)



