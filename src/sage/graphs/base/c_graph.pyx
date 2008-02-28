
#*******************************************************************************
#        Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

cdef class CGraph:
    cdef int add_arc_unsafe(self, int u, int v) except? -1:
        raise NotImplementedError()
    cdef int has_arc_unsafe(self, int u, int v) except? -1:
        raise NotImplementedError()
    cdef int del_arc_unsafe(self, int u, int v) except? -1:
        raise NotImplementedError()
    cdef int out_neighbors_unsafe(self, int u, int *neighbors, int size) except? -2:
        raise NotImplementedError()
    cdef int in_neighbors_unsafe(self, int u, int *neighbors, int size) except? -2:
        raise NotImplementedError()

