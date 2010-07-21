
#*****************************************************************************
#      Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include '../../../ext/cdefs.pxi'
include '../../../ext/stdsage.pxi'
include 'data_structures_pxd.pxi' # includes bitsets

from sage.rings.integer cimport Integer

cdef struct aut_gp_and_can_lab:
    int *generators
    int num_gens
    StabilizerChain *group
    int *relabeling

cdef aut_gp_and_can_lab *get_aut_gp_and_can_lab( void *,
    PartitionStack *, int,
    bint (*)(PartitionStack *, void *),
    int (*)(PartitionStack *, void *, int *, int),
    int (*)(int *, int *, void *, void *), bint, StabilizerChain *) except NULL



