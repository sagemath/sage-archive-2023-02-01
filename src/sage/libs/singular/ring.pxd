"""
Wrapper for Singular's Rings

AUTHOR:

- Martin Albrecht (2009-07): initial implementation
"""
#*****************************************************************************
#       Copyright (C) 2009 Martin Albrecht <malb@informatik.uni-bremen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.singular.decl cimport ring

# create a new singular ring
cdef ring *singular_ring_new(base_ring, n, names, term_order) except NULL

# carefully delete a ring
cdef void singular_ring_delete(ring *ring)
