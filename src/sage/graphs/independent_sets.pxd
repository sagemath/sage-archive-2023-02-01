from libc.stdint cimport uint32_t, uint64_t
include "sage/misc/binary_matrix.pxi"

cdef class IndependentSets:
    cdef binary_matrix_t g
    cdef list vertices
    cdef dict vertex_to_int
    cdef int n
    cdef int i
    cdef int count_only
    cdef int maximal

