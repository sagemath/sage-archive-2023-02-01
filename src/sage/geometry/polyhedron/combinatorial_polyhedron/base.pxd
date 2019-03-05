cimport cython
from libc.stdint cimport uint64_t
from sage.ext.memory_allocator cimport MemoryAllocator
from sage.structure.sage_object cimport SageObject

cdef int vertex_list_to_bit_repr(tuple vertex_list, uint64_t *output,
                                 size_t face_length) except -1

cdef int incidences_to_bit_repr(tuple incidences, uint64_t *output,
                                size_t face_length) except -1

cdef size_t bit_repr_to_vertex_list(uint64_t *face, size_t *output,
                                    size_t face_length) except -1

@cython.final
cdef class ListOfFaces(SageObject):
    cdef MemoryAllocator _mem

    # ``face_length`` is the length of each face in terms of ``uint64_t``.
    cdef readonly size_t nr_faces, face_length, nr_vertices

    # ``data`` points to the raw data.
    # It will be of "type" ``uint64_t[nr_faces][face_length]``
    cdef uint64_t **data

cpdef int calculate_dimension(ListOfFaces faces) except -2

@cython.final
cdef class FaceIterator(SageObject):
    cdef readonly bint dual         # if 1, then iterate over dual Polyhedron
    cdef uint64_t *face             # the current face of the iterator
    cdef size_t *atom_repr          # a place where atom-representaion of face will be stored
    cdef size_t *coatom_repr        # a place where coatom-representaion of face will be stored
    cdef int current_dimension      # dimension of current face, dual dimension if ``dual``
    cdef int dimension              # dimension of the polyhedron
    cdef int nr_lines               # ``_nr_lines`` of ``CombinatorialPolyhedron``
    cdef int request_dimension      # only faces of this (dual?) dimension are considered
    cdef int lowest_dimension       # don't consider faces bewow this (dual?) dimension
    cdef MemoryAllocator _mem
    cdef tuple newfaces_lists       # tuple to hold the ListOfFaces corresponding to maybe_newfaces
    cdef size_t face_length         # stores length of the faces in terms of uint64_t
    cdef tuple _V, _H, _equalities  # some copies from ``CombinatorialPolyhedron``

    # Atoms and coatoms are the vertices/facets of the Polyedron.
    # If ``dual == 0``, then coatoms are facets, atoms vertices and vice versa.
    cdef ListOfFaces atoms, coatoms

    # ``visited_all`` points to faces, of which we have visited all faces already.
    # The number of faces in ``visited_all` might depend on the current dimension:
    #     Consider we visit the facets A,B of some face F.
    #     We will first visit all faces of A and then add A to visited_all.
    #     Then we visit all faces of B and add B to visited_all.
    #     Then we have visited F completely.
    #     Instead of having A and B in ``visited_all`` we will point to F.
    #     In this way, we will append ``visited_all`` in lower dimension, but
    #     will ignore those changes when going up in dimension again.
    #     This is why the number of faces in ``visited_all``depends on dimension.
    cdef uint64_t **visited_all
    cdef size_t *nr_visited_all

    # ``maybe_newfaces`` is where all possible facets of a face are stored.
    # In dimension ``dim`` when visiting all faces of some face,
    # the intersections with other faces are stored in ``newfaces2[dim]``.
    cdef uint64_t ***maybe_newfaces

    # ``newfaces`` will point to those faces in ``maybe_newfaces``
    # that are of codimension 1 and not already visited.
    cdef uint64_t ***newfaces
    cdef size_t *nr_newfaces  # number of newfaces for each dimension

    # After having visited a face completely, we want to add it to ``visited_all``.
    # ``first_dim[i]`` will indicate, wether there is one more face in
    # ``newfaces[i]`` then ``nr_newfaces[i]`` suggests
    # that has to be added to ``visited_all``.
    # If ``first_time[i] == False``, we still need to
    # add ``newfaces[i][nr_newfaces[i]]`` to ``visited_all``.
    cdef bint *first_time

    # The number of elements in newfaces[current_dimension],
    # that have not been visited yet.
    cdef size_t yet_to_visit

    cdef inline int next_face(self) except -1
    cdef inline int next_face_loop(self) except -1
    cdef size_t length_atom_repr(self) except -1
    cdef size_t set_coatom_repr(self) except -1
    cdef size_t set_atom_repr(self) except -1

@cython.final
cdef class ListOfAllFaces(SageObject):
    cdef MemoryAllocator _mem
    cdef int dimension              # dimension of Polyhedron
    cdef readonly bint dual         # if True, then List of all faces by dual Polyhedron
    cdef size_t face_length         # stores length of the faces in terms of uint64_t
    cdef tuple _V, _H, _equalities  # some copies from CombinatorialPolyhedron
    cdef size_t *f_vector           # a copy of the f-vector, is reversed if dual
    cdef size_t *face_counter       # how many faces of each dimension have been initialized
    cdef size_t *atom_repr          # a place where atom-representaion of face will be stored
    cdef size_t *coatom_repr        # a place where coatom-representaion of face will be stored

    # Atoms and coatoms are the vertices/facets of the Polyedron.
    # If ``dual == 0``, then coatoms are facets, atoms vertices and vice versa.
    cdef ListOfFaces atoms, coatoms

    # All faces are stored in ``faces``. ``faces[i]`` stores all faces of
    # dimension `i` in Bit-representation (of atoms).
    cdef uint64_t *** faces
    cdef tuple faces_mem  # tuple to hold the ListOfFaces corresponding to faces

    # It follows a number of attributes to iterate over all incidences.
    # After initialization, ``ListOfAllFaces`` can iterate over all incidences
    # of faces of dimension ``incidence_dim_one`` and ``incidence_dim_two``.
    cdef int is_incidence_initialized
    cdef int incidence_dim_one
    cdef int incidence_dim_two
    cdef size_t incidence_counter_one  # walks trough faces of incidence_dim_one
    cdef size_t incidence_counter_two  # walks through all indices of coatoms (for each face in incidence_dim_one)

    # Intersection of ``faces[incidence_dim_one][incidence_counter_one]`` with
    # ``coatoms[incidence_counter_two]``.
    # If this is contained in ``faces[incidence_dim_two]``, there is an incidence.
    cdef uint64_t *incidence_face

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
    cdef size_t set_coatom_repr(self, int dimension, size_t index) except -1
    cdef size_t set_atom_repr(self, int dimension, size_t index) except -1
    cdef void incidence_init(self, int dimension_one, int dimension_two)
    cdef inline bint next_incidence(self, size_t *one, size_t *two)
    cdef inline bint next_incidence_loop(self, size_t *one, size_t *two)
    cdef inline bint next_trivial_incidence(self, size_t *one, size_t *two)
    cdef inline bint next_trivial_incidence2(self, size_t *one, size_t *two)

@cython.final
cdef class CombinatorialPolyhedron(SageObject):
    cdef tuple _V                       # the names of VRep, if they exist
    cdef dict _Vinv                     # dictionary to look up enumeration of vertices
    cdef tuple _H                       # the names of HRep, if they exist
    cdef tuple _equalities              # stores equalities, given on input (might belong to Hrep)
    cdef int _dimension                 # stores dimension, -2 on init
    cdef unsigned int _length_Hrep      # Hrep might include equalities
    cdef unsigned int _length_Vrep      # Vrep might include rays/lines
    cdef size_t _nr_facets              # length Hrep without equalities
    cdef bint _unbounded                # ``True`` iff Polyhedron is unbounded
    cdef int _nr_lines                  # number of affinely independet lines in the Polyhedron
    cdef ListOfFaces bitrep_facets      # facets in Bit-Representation
    cdef ListOfFaces bitrep_vertices    # vertices in Bit-Representation
    cdef tuple _f_vector

    # Edges, ridges and incidences are stored in a pointer of pointers,
    # This number determines how many edges are stored in ``edges[0]``,
    # how many ridges are stored in ``ridges[0]`` etc.
    cdef size_t _length_edges_list


    cdef size_t **_edges                    # stores edges
    cdef size_t _nr_edges
    cdef size_t **_ridges                   # stores ridges
    cdef size_t _nr_ridges
    cdef size_t **_face_lattice_incidences  # stores incidences in Hasse diagram
    cdef size_t _nr_face_lattice_incidences
    cdef ListOfAllFaces _all_faces          # class to generate Hasse diagram incidences
    cdef tuple _mem_tuple                   # stores MemoryAllocators for edges, ridges etc.

    cdef FaceIterator _face_iter(self, bint dual)
    cdef int _calculate_f_vector(self) except -1
    cdef int _calculate_edges(self, dual) except -1
    cdef int _calculate_ridges(self, dual) except -1
    cdef int _calculate_face_lattice_incidences(self) except -1
