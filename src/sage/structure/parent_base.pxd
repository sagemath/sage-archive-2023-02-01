###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

cimport parent_old

cdef class ParentWithBase(parent_old.Parent):
    # DO NOT OVERRIDE ANY OF THE FOLLOWING
    cdef base_extend_recursive_c(self, ParentWithBase X)
    cdef base_extend_canonical_c(self, ParentWithBase X)
    cdef base_extend_canonical_sym_c(self, ParentWithBase X)
    # THIS ONE IS PRIVATE, FOR RECURSION SAKE
    cdef _base_extend_canonical_rec(self, ParentWithBase X)

