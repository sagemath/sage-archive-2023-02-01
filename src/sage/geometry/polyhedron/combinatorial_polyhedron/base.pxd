cimport cython
from libc.stdint cimport uint64_t
from sage.ext.memory_allocator cimport MemoryAllocator
from sage.structure.sage_object cimport SageObject

cdef int char_from_tuple(tuple tup, uint64_t *output,
                         size_t face_length) except -1

cdef int char_from_incidence(tuple incidences, uint64_t *output,
                             size_t face_length) except -1

cdef size_t vertex_repr_from_bitrep(uint64_t *face, size_t *output,
                                    size_t face_length) except -1

@cython.final
cdef class ListOfFaces(SageObject):
    cdef uint64_t **data
    cdef MemoryAllocator _mem
    cdef readonly size_t nr_faces, face_length, nr_vertices
    # face_length is the length of the face in terms of uint64_t

cpdef int calculate_dimension(ListOfFaces faces) except -2

@cython.final
cdef class FaceIterator(SageObject):
    cdef readonly bint dual  # if 1, then iterate over dual Polyhedron
    cdef uint64_t *face  # the current face of the iterator
    cdef size_t *atom_repr_face  # a place where repr will be stored
    cdef size_t *coatom_repr_face  # a place where repr will be stored
    cdef int current_dimension  # of current face, dual dimension if ``dual``
    cdef int dimension  # total dimension
    cdef int nr_lines  # _nr_lines of CombinatorialPolyhedron
    cdef int record_dimension  # only faces of this (dual) dimension are considered
    cdef int lowest_dimension  # don't consider faces bewow this (dual) dimension
    cdef MemoryAllocator _mem
    cdef tuple newfaces_lists  # tuple to hold the ListOfFaces corresponding to maybe_newfaces
    cdef size_t face_length  # stores length of the faces in terms of uint64_t
    cdef ListOfFaces atoms, coatoms
    # atoms and coatoms are the vertices/facets of the Polyedron
    # if dual == 0, then coatoms are facets, atoms vertices and vice versa

    # some copies from CombinatorialPolyhedron
    cdef tuple _V, _H, _equalities


    cdef uint64_t **forbidden
    cdef size_t *nr_forbidden
    # forbidden stores faces, of which we have visited all faces already
    # the length might alter:
    # consider we visit the facets A,B of some face F
    # we will first visit all faces of A and then add A to forbidden,
    # then visit all faces of B and add B to forbidden
    # then we have visited F completely, so "delete" A and B from forbidden
    # and add F to forbidden
    # this is why there might be different nr_forbidden, depending on dimension
    cdef uint64_t ***maybe_newfaces
    # maybe_newfaces is where all possible facets of a face are stored
    # e.g. when visiting all faces of ``newfaces2[dim][i]``
    # the possible facets of ``newfaces2[dim][i]`` will be stored in
    # ``newfaces[dim-1]``
    cdef uint64_t ***newfaces
    cdef size_t *nr_newfaces
    # ``newfaces`` will point to those faces in ``maybe_newfaces``, that are
    # of codimension 1 and not contained in a face in forbidden
    # (so not visited before)
    cdef bint *first_time
    # adding a face to forbidden is a bit more complicated
    # when having visited all faces of A we want to add A to to forbidden
    # but nr_newfaces[dim A] has gone down alread,
    # if first_time[dim A] == False, we know, that there is one more face,
    # we have visited already
    cdef size_t yet_to_yield
    # the number of elements in newfaces[current_dimension],
    # that have not been yielded yet

    cdef inline int next_face(self) except -1
    cdef inline int next_face_loop(self) except -1
    cdef size_t length_atom_repr(self) except -1
    cdef size_t coatom_repr(self) except -1
    cdef size_t atom_repr(self) except -1

@cython.final
cdef class ListOfAllFaces(SageObject):
    cdef MemoryAllocator _mem
    cdef int dimension  # dimension of Polyhedron
    cdef readonly bint dual  # if True, then List of all faces by dual Polyhedron
    cdef size_t face_length  # stores length of the faces in terms of uint64_t
    cdef ListOfFaces atoms, coatoms
    # atoms and coatoms are the vertices/facets of the Polyedron
    # if dual == 0, then coatoms are facets, atoms vertices and vice versa
    cdef size_t *atom_repr_face  # a place where repr will be stored
    cdef size_t *coatom_repr_face  # a place where repr will be stored
    # some copies from CombinatorialPolyhedron
    cdef tuple _V, _H, _equalities
    cdef size_t *f_vector  # a copy of the f_vector, is reversed if dual
    cdef size_t *face_counter  # how many faces of each dimension have been initialized
    cdef tuple faces_mem # tuple to hold the ListOfFaces corresponding to faces
    cdef uint64_t *** faces
    # faces is the data of all faces, faces[dim] stores all faces of dimension dim

    # it follows a number of variables to determine incidences
    cdef int incidence_is_initialized
    cdef int incidence_dim_one
    cdef int incidence_dim_two
    # inicidences are determined by intersecting each face of incidence_dim_one
    # with each coatom
    cdef size_t incidence_counter_one  # walks trough faces of incidence_dim_one
    cdef size_t incidence_counter_two
    # incidences_counter_two walks through the coatoms for each face in
    # incidence_dim_one
    cdef uint64_t *incidence_face  # this is the place, where each intersection is stored

    cdef int _add_face(self, int face_dim, uint64_t *face) except -1
    cdef int _sort(self) except -1
    cdef int _sort_one_list(self, uint64_t **faces, size_t nr_faces) except -1
    cdef int _sort_one_list_loop(
            self, uint64_t **inp, uint64_t **output1,
            uint64_t **output2, size_t nr_faces) except -1
    cdef inline size_t find_face(self, int dimension, uint64_t *face) except -1
    cdef inline bint is_smaller(self, uint64_t *one, uint64_t *two)
    cdef inline int is_equal(self, int dimension, size_t index,
                             uint64_t *face) except -1
    cdef size_t coatom_repr(self, int dimension, size_t index) except -1
    cdef size_t atom_repr(self, int dimension, size_t index) except -1
    cdef void incidence_init(self, int dimension_one, int dimension_two)
    cdef inline int next_trivial_incidence(self, size_t *one, size_t *two)
    cdef inline int next_trivial_incidence2(self, size_t *one, size_t *two)
    cdef int next_incidence(self, size_t *one, size_t *two)

@cython.final
cdef class CombinatorialPolyhedron(SageObject):
    cdef tuple _V  # the names of VRep, if they exist
    cdef dict _Vinv  # dictionary to look up enumeration of vertices
    cdef tuple _H  # the names of HRep, if they exist
    cdef tuple _equalities
    cdef int _dimension  # stores dimension, -2 on init
    cdef unsigned int _length_Hrep  # Hrep might include equalities
    cdef unsigned int _length_Vrep  # Vrep might include rays/lines
    cdef size_t _nr_facets  # length Hrep without equalities
    cdef bint _unbounded  # 1 iff Polyhedron is unbounded
    cdef int _nr_lines  # nr of affinely independet lines in the Polyhedron
    cdef ListOfFaces bitrep_facets  # facets in Bit-Representation
    cdef ListOfFaces bitrep_vertices  # vertices in Bit-Representation
    cdef tuple _f_vector
    cdef size_t _length_edges_list
    # edges, ridges and incidences are stores in a pointer of pointers,
    # this number determines how many edges are stored in ``edges[0]`` etc.
    cdef size_t **_edges
    cdef size_t **_ridges
    cdef size_t **_face_lattice_incidences
    cdef size_t _nr_edges
    cdef size_t _nr_ridges
    cdef size_t _nr_face_lattice_incidences
    cdef ListOfAllFaces _all_faces
    cdef tuple _mem_tuple  # stores MemoryAllocators
    cdef FaceIterator _face_iter(self, bint dual)
    cdef int _calculate_f_vector(self) except -1
    cdef int _calculate_edges(self, dual) except -1
    cdef int _calculate_ridges(self, dual) except -1
    cdef int _calculate_face_lattice_incidences(self) except -1
