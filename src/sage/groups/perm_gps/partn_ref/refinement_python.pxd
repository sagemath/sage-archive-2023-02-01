
#*****************************************************************************
#      Copyright (C) 2006 - 2009 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************


include '../../../ext/cdefs.pxi'
include '../../../ext/stdsage.pxi'
include 'data_structures_pxd.pxi' # includes bitsets


from sage.rings.integer cimport Integer
from automorphism_group_canonical_label cimport get_aut_gp_and_can_lab, aut_gp_and_can_lab_return
from double_coset cimport double_coset


cdef class PythonPartitionStack:
    cdef PartitionStack *c_ps

