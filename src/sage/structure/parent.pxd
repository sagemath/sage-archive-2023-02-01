###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

cimport sage_object
cdef class Parent(sage_object.SageObject):
    cdef public object _has_coerce_map_from

    #########################################
    # Canonical Coercion Methods
    cdef has_coerce_map_from_c(self, S)
    cdef has_coerce_map_from_c_impl(self, S)
    cdef _coerce_c(self, x)
    cdef _coerce_c_impl(self, x)
    cdef _coerce_self_c(self, x)
    cdef public object __an_element
    cdef _an_element_c_impl(self)
    cdef _an_element_c(self)

    ################################################
    # Comparison of parent objects
    cdef _richcmp(left, right, int op)
    cdef int _cmp_c_impl(left, Parent right) except -2




