
#*******************************************************************************
#        Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

cdef class CGraph:
    cdef int num_verts
    cdef int num_arcs
    cdef int *in_degrees
    cdef int *out_degrees

    cdef int add_arc_unsafe(self, int, int)
    cdef int has_arc_unsafe(self, int, int)
    cdef int del_arc_unsafe(self, int, int)
    cdef int out_neighbors_unsafe(self, int, int *, int)
    cdef int in_neighbors_unsafe(self, int, int *, int)


# TODO: edge functions!


