cimport cython
from libc.stdint cimport uint64_t
from sage.ext.memory_allocator cimport MemoryAllocator

cdef int vertex_list_to_bit_repr(tuple vertex_list, uint64_t *output,
                                 size_t face_length) except -1

cdef int incidences_to_bit_repr(tuple incidences, uint64_t *output,
                                size_t face_length) except -1

cdef size_t bit_repr_to_vertex_list(uint64_t *face, size_t *output,
                                    size_t face_length) except -1

@cython.final
cdef class ListOfFaces:
    cdef MemoryAllocator _mem

    # ``face_length`` is the length of each face in terms of ``uint64_t``.
    cdef readonly size_t nr_faces, face_length, nr_vertices

    # ``data`` points to the raw data.
    # It will be of "type" ``uint64_t[nr_faces][face_length]``
    cdef uint64_t **data

    cpdef int calculate_dimension(self) except -2
    cdef int calculate_dimension_loop(self, uint64_t **faces, size_t nr_faces,
                                      size_t face_length) except -2
