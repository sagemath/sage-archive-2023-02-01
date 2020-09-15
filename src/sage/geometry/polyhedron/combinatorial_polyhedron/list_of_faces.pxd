cimport cython
from libc.stdint cimport uint64_t
from sage.ext.memory_allocator cimport MemoryAllocator
from .face_data_structure   cimport *

cdef struct face_list_s:
    face_t* faces
    size_t n_faces
    size_t total_n_faces
    size_t n_atoms
    size_t n_coatoms
    bint polyhedron_is_simple
    bint* is_not_new_face

ctypedef face_list_s face_list_t[1]

@cython.final
cdef class ListOfFaces:
    cdef MemoryAllocator _mem

    # ``face_length`` is the length of each face in terms of ``uint64_t``.
    cdef readonly size_t n_faces, face_length, n_atoms

    # ``data`` points to the raw data.
    # It will be of "type" ``uint64_t[n_faces][face_length]``
    cdef uint64_t **data

    cpdef int compute_dimension(self) except -2
    cdef int compute_dimension_loop(self, uint64_t **faces, size_t n_faces,
                                      size_t face_length) except -2

    cpdef ListOfFaces pyramid(self)
