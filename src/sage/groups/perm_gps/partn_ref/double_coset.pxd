
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

cdef int *double_coset( object, object, int **, int *, int,
    bint (*)(PartitionStack *, object),
    int (*)(PartitionStack *, object, int *, int),
    int (*)(int *, int *, object, object))



