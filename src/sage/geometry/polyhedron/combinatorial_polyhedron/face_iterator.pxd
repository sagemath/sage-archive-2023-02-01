cimport cython
from libc.stdint                cimport uint64_t
from sage.ext.memory_allocator  cimport MemoryAllocator
from sage.structure.sage_object cimport SageObject
from .list_of_faces             cimport ListOfFaces

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
