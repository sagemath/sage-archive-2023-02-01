cimport cython
from libc.stdint                cimport uint64_t
from sage.ext.memory_allocator  cimport MemoryAllocator
from sage.structure.sage_object cimport SageObject
from .face_iterator             cimport FaceIterator
from .list_of_faces             cimport ListOfFaces
from .list_of_all_faces         cimport ListOfAllFaces

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
