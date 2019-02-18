cimport cython
from libc.stdint cimport uint64_t
from sage.ext.memory_allocator cimport MemoryAllocator
from sage.structure.sage_object cimport SageObject

cdef int char_from_tuple(tuple tup, uint64_t *output,
                         size_t face_length) except 0

cdef int char_from_incidence(tuple incidences, uint64_t *output,
                             size_t face_length) except 0

cdef size_t vertex_repr_from_bitrep(uint64_t *face, size_t *output,
                                    size_t face_length) except? 0


cdef ListOfFaces get_facets_from_incidence_matrix(tuple matrix)

cdef ListOfFaces get_vertices_from_incidence_matrix(tuple matrix)

cdef ListOfFaces get_facets_bitrep_from_facets_tuple(tuple facets_input,
                                                     size_t nr_vertices)

cdef ListOfFaces get_vertices_bitrep_from_facets_tuple(tuple facets_input,
                                                       size_t nr_vertices)

@cython.final
cdef class ListOfFaces:
    cdef uint64_t **data
    cdef MemoryAllocator _mem
    cdef size_t nr_faces, face_length, nr_vertices

cdef int calculate_dimension(ListOfFaces faces) except -2

@cython.final
cdef class FaceIterator:
    cdef uint64_t *face
    cdef int current_dimension, dimension, record_dimension, lowest_dimension
    cdef MemoryAllocator _mem
    cdef uint64_t ***newfaces2
    cdef tuple newfaces_lists
    cdef uint64_t ***newfaces
    cdef uint64_t **forbidden
    cdef size_t *nr_faces
    cdef size_t *nr_forbidden
    cdef int *first_time
    cdef size_t yet_to_yield, face_length, nr_facets, nr_vertices
    cdef size_t *output1
    cdef size_t *output2
    cdef int nr_lines

    cdef void set_record_dimension(self, int dim)
    cdef inline int next_face(self) except -1
    cdef inline int next_face_loop(self) except -1
    cdef size_t length_vertex_repr(self) except -1
    cdef size_t facet_repr(self, size_t *output) except -1
    cdef size_t vertex_repr(self, size_t *output) except -1
    cdef size_t *get_output1_array(self) except NULL
    cdef size_t *get_output2_array(self) except NULL

@cython.final
cdef class ListOfAllFaces:
    # cdef tuple lists_facet_repr #might need this for flag-vector
    cdef MemoryAllocator _mem
    cdef tuple lists_vertex_repr
    cdef size_t nr_facets
    cdef size_t nr_vertices
    cdef size_t face_length_vertex
    # cdef size_t face_length_facet #might need this for flag-vector
    cdef int dimension
    cdef size_t *face_counter
    cdef size_t *f_vector
    cdef int is_sorted
    cdef size_t *output1
    cdef size_t *output2
    cdef int incidence_dim_one
    cdef int incidence_dim_two
    cdef size_t incidence_counter_one
    cdef size_t incidence_counter_two
    cdef ListOfFaces incidence_face_mem
    cdef uint64_t *incidence_face
    cdef uint64_t **facets
    cdef int incidence_is_initialized
    cdef uint64_t ***data_vertex

    cdef int add_face(self, int face_dim, uint64_t *face) except 0
    cdef int sort(self) except 0
    cdef int _sort_one_list(self, uint64_t **faces, size_t nr_faces) except 0
    cdef int _sort_one_list_loop(
            self, uint64_t **inp, uint64_t **output1,
            uint64_t **output2, size_t nr_faces) except 0
    cdef inline size_t find_face(self, int dimension, uint64_t *face) except -1
    cdef inline int is_smaller(self, uint64_t *one, uint64_t *two)
    cdef inline int is_equal(self, int dimension, size_t index,
                             uint64_t *face) except -1
    cdef size_t facet_repr(self, int dimension, size_t index,
                           size_t *output) except -1
    cdef size_t vertex_repr(self, int dimension, size_t index,
                            size_t *output) except -1
    cdef size_t *get_output1_array(self) except NULL
    cdef size_t *get_output2_array(self) except NULL
    cdef void incidence_init(self, int dimension_one, int dimension_two)
    cdef inline int next_trivial_incidence(self, size_t *one, size_t *two)
    cdef inline int next_trivial_incidence2(self, size_t *one, size_t *two)
    cdef int next_incidence(self, size_t *one, size_t *two)

@cython.final
cdef class CombinatorialPolyhedron(SageObject):
    cdef tuple _V
    cdef tuple _H
    cdef tuple _equalities
    cdef dict _Hinv
    cdef dict _Vinv
    cdef int is_trivial  # in some instances the polyhedron might not
    # have facets or otherwise produce errors in the C function
    cdef int _dimension  # manually set, if `is_trivial`
    cdef unsigned int _length_Hrep
    cdef unsigned int _length_Vrep
    cdef int _unbounded
    cdef int _nr_lines
    cdef ListOfFaces bitrep_facets
    cdef ListOfFaces bitrep_vertices
    cdef int _polar
    cdef size_t *_f_vector
    cdef size_t _length_edges_list
    cdef size_t **_edges
    cdef size_t **_ridges
    cdef size_t **_face_lattice_incidences
    cdef size_t _nr_edges
    cdef size_t _nr_ridges
    cdef size_t _nr_face_lattice_incidences
    cdef ListOfAllFaces _all_faces
    cdef size_t _nr_facets
    cdef tuple _MemoryAllocators

    cdef int _calculate_f_vector(self) except 0
    cdef int _calculate_edges(self) except 0
    cdef int _calculate_ridges(self) except 0
    cdef int _calculate_face_lattice_incidences(self) except 0
    cdef int _record_all_faces_helper(self) except 0
