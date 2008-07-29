###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

cimport parent
cdef class Parent(parent.Parent):

    # List consisting of Morphisms (from anything to self)
    # and Parents for which the __call__ method of self
    # results in natural coercion.
    # Initalized at ring creation.
    cdef _coerce_from_list
    # Hashtable of everything we've (possibliy recursively) discovered so far.
    cdef _coerce_from_hash

    # List consisting of Actions (either by or on self)
    # and Parents for which self._rmul_ and/or self._lmul_
    # do the correct thing.
    # Initalized at ring creation.
    cdef _action_list
    # Hashtable of everything we've (possibliy recursively) discovered so far.
    cdef _action_hash

    # returns a Morphism from S to self, or None
    cdef coerce_map_from_c(self, S)
    cdef coerce_map_from_c_impl(self, S)

    # returns the Action by/on self on/by S
    # corresponding to op and self_on_left
    cdef get_action_c(self, S, op, bint self_on_left)
    cdef get_action_c_impl(self, S, op, bint self_on_left)



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




