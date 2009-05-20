
#*****************************************************************************
#      Copyright (C) 2006 - 2008 Robert L. Miller <rlmillster@gmail.com>
#      Copyright (C) 2009 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************


include '../../../ext/cdefs.pxi'
include '../../../ext/stdsage.pxi'
include 'data_structures_pxd.pxi' # includes bitsets


from sage.rings.integer cimport Integer
from double_coset cimport double_coset

# name of the three functions to customize
cdef int refine_list(PartitionStack *, object, int *, int)
cdef int compare_lists(int *, int *, object, object)
cdef bint all_list_children_are_equivalent(PartitionStack *, object)