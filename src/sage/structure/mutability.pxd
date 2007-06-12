"""
Mutability -- Pyrex Implementation
"""

##########################################################################
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
##########################################################################

cdef class Mutability:
    cdef bint _is_immutable
    cdef _require_mutable_cdef(self)
