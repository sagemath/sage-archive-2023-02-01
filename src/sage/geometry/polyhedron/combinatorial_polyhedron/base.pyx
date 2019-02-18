#distutils: language = c++

r"""
Several algorithms working implicitly with the hasse_diagram
(of a polytope), including calculating the f_vector, dimension,
flags, level-sets and even the face lattice.

This computes implicitely a finite atomic and coatomic lattices,
where every interval of length two has at exactly 4 elements
(known as the diamond property).
In particular this module calculates quickly the f_vector of polytopes.
The input must be a tuple of coatoms given each by a tuple of atoms.
The atoms must be labeled 0,...,n.


AUTHOR:

- Jonathan Kliem (2019-02)
"""

#*****************************************************************************
#       Copyright (C) 2019 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, division, print_function
from sage.rings.integer import Integer
from sage.graphs.graph import Graph
from sage.graphs.digraph import DiGraph
from sage.combinat.posets.lattices import FiniteLatticePoset
from sage.geometry.polyhedron.base import is_Polyhedron
from sage.geometry.lattice_polytope import is_LatticePolytope
from sage.structure.element import is_Matrix
from sage.misc.misc import is_iterator

from libc.stdint cimport uint64_t
from libc.string cimport memcmp, memcpy
from cysignals.memory cimport sig_realloc
from cysignals.signals cimport sig_check, sig_on, sig_off, sig_block, \
                               sig_unblock

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef extern from "helper.cc":
    cdef const size_t chunksize

    cdef void intersection(uint64_t *A, uint64_t *B, uint64_t *C,
                           size_t face_length)
    # C = A & B
    # will set C to be the intersection of A and B
    # `face_length` is the length of A, B and C in terms of uint64_t

    cdef size_t get_next_level(\
        uint64_t **faces, size_t lenfaces, uint64_t **nextfaces,
        uint64_t **nextfaces2, uint64_t **forbidden,
        size_t nr_forbidden, size_t face_length)
    # intersects the first `lenfaces - 1` faces of `faces`
    # with `faces[lenfaces-1]`
    # stores the faces in `nextfaces`
    # (so `nextfaces` must be the data of `ListOfFaces`)
    # determines which ones are exactly of one dimension less
    # by considering containment
    # newfaces2 will point at those of exactly one dimension less
    # which are not contained in any of the faces in `forbidden`
    # returns the number of those faces

    cdef size_t CountFaceBits(uint64_t *A, size_t face_length)

    cdef size_t facet_repr_from_bitrep(
        uint64_t *face, uint64_t **facets, size_t *output,
        size_t nr_facets, size_t face_length)
    # Writes the facet_repr of the current face in output.
    # Returns the length of the representation.

cdef uint64_t vertex_to_bit_dictionary[64]
# this dictionary helps storing a vector of 64 incidences as uint64_t,
# where each bit represents an incidence
for i in range(64):
    vertex_to_bit_dictionary[i] = 2**(64-i-1)

cdef int char_from_tuple(tuple tup, uint64_t *output,
                         size_t face_length) except 0:
    r"""
    Converts a tuple into a bit-representation. Stores it in `output`.

    The first bit represent the entry `0` and is set to one, iff `0` is in
    `tup`. The second bit represents `1` and so on.

    INPUT:

    - `tup` -- tuple of pairwise distinct positive integers in
      `range(face_length*64)`.
    - `output` -- array of `uint64_t` of length `face_lenght`.
    - `face_length`.

    OUTPUT:

    - `1` on sucess
    - fills and initializes `output`.

    EXAMPLES::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport char_from_tuple
        ....:
        ....: cdef uint64_t output
        ....: cdef uint64_t[2] output2
        ....:
        ....: char_from_tuple((62,63), &output, 1)
        ....: print(str(output))
        ....:
        ....: char_from_tuple((61,63,125), output2, 2)
        ....: print(str(output2))
        ....:
        ....: tup = tuple(i for i in range(64))
        ....: char_from_tuple(tup, &output, 1)
        ....: print(str(output))
        ....: print(str(output+1))
        ....:
        ....: try:
        ....:     char_from_tuple((62,70), &output, 1)
        ....: except:
        ....:     print('out of range')
        ....:
        ....: try:
        ....:     char_from_tuple((-1,12), &output, 1)
        ....: except:
        ....:     print('negatives not allowed')
        ....:
        ....: try:
        ....:     char_from_tuple((0,0), &output, 1)
        ....: except:
        ....:     print('duplicates not allowed')
        ....: ''') # long time
        3
        [5L, 4L]
        18446744073709551615
        0
        out of range
        negatives not allowed
        duplicates not allowed
    """
    cdef size_t entry, position, value, i
    for i in range(face_length):
        output[i] = 0
    if unlikely(len(tup) != len(set(tup))):
        raise ValueError("entries of `tup` are not distinct")
    for entry in tup:
        value = entry % 64
        position = entry//64
        if unlikely(position >= face_length):
            raise IndexError("output too small to represent %s"%position)
        output[position] += vertex_to_bit_dictionary[value]
    return 1

cdef int char_from_incidence(tuple incidences, uint64_t *output,
                             size_t face_length) except 0:

    r"""
    Converts a tuple of incidences into a bit-representation.

    Stores it in `output`. Each entry in `incidences` represents a bit in
    `output`. It is set to `1`, iff the entry in `incidences` in non-zero.

    INPUT:

    - `tup` -- tuple of integers representing incidences of length at most
      `face_length*64`.
    - `output` -- array of `uint64_t` of length `face_lenght`.
    - `face_length`.

    OUTPUT:

    - `1` on sucess
    - fills and initializes `output`.

    EXAMPLES::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport char_from_incidence
        ....:
        ....: cdef uint64_t output
        ....: cdef uint64_t[2] output2
        ....:
        ....: tup = tuple(0 for _ in range(62)) + (1,1)
        ....: char_from_incidence(tup, &output, 1)
        ....: print(str(output))
        ....:
        ....: tup = tuple(0 for _ in range(61)) + (1,0,1)
        ....: tup += tuple(0 for _ in range(61)) + (1,)
        ....: char_from_incidence(tup, output2, 2)
        ....: print(str(output2))
        ....:
        ....: tup = tuple(1 for _ in range(64))
        ....: char_from_incidence(tup, &output, 1)
        ....: print(str(output))
        ....: print(str(output+1))
        ....:
        ....: tup = tuple(0 for _ in range(70))
        ....: try:
        ....:     raise IndexError()
        ....: except:
        ....:     print('out of range')
        ....: ''')  # long time
        3
        [5L, 4L]
        18446744073709551615
        0
        out of range
    """
    cdef size_t position, value, i, entry
    cdef size_t length = len(incidences)
    for i in range(face_length):
        output[i] = 0
    if unlikely(length > 64*face_length):
        raise IndexError("output too small to represent all incidences")
    for entry in range(length):
        if incidences[entry]:
            value = entry % 64
            position = entry//64
            output[position] += vertex_to_bit_dictionary[value]
    return 1

cdef ListOfFaces get_facets_from_incidence_matrix(tuple matrix):
    r"""
    Initializes facets in bit-representation as `ListOfFaces`.

    INPUT:

    - `matrix` -- an incidence matrix as tuple of tuples.

    .. NOTE::

        :meth:`sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`
        also has equalities appearing the incidence matrix. Those must be
        removed.

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, vertex_repr_from_bitrep, \
        ....:         get_facets_from_incidence_matrix
        ....: from sage.geometry.polyhedron.library import polytopes
        ....:
        ....: cdef ListOfFaces facets
        ....: cdef size_t[24] output
        ....: cdef size_t length
        ....: P = polytopes.permutahedron(4)
        ....: data = P.incidence_matrix()
        ....: rg = range(data.nrows())
        ....: matrix = tuple(tuple(data[i,j] for i in rg)
        ....:                for j in range(data.ncols())
        ....:                if not all(data[i,j] for i in rg))
        ....:
        ....: facets = get_facets_from_incidence_matrix(matrix)
        ....: for i in range(facets.nr_faces):
        ....:     length = vertex_repr_from_bitrep(facets.data[i], output,
        ....:                                      facets.face_length)
        ....:     print(tuple(output[j] for j in range(length)))
        ....: ''') # long time
        (18, 19, 20, 21, 22, 23)
        (9, 11, 15, 17, 21, 23)
        (16, 17, 22, 23)
        (0, 1, 2, 3, 4, 5)
        (2, 4, 8, 10)
        (0, 1, 6, 7)
        (6, 7, 12, 13, 18, 19)
        (3, 5, 9, 11)
        (1, 3, 7, 9, 13, 15)
        (0, 2, 6, 8, 12, 14)
        (12, 14, 18, 20)
        (4, 5, 10, 11, 16, 17)
        (8, 10, 14, 16, 20, 22)
        (13, 15, 19, 21)
    """
    cdef ListOfFaces facets = ListOfFaces(len(matrix), len(matrix[0]))
    cdef uint64_t **facets_data = facets.data
    cdef int i
    for i in range(len(matrix)):
        char_from_incidence(matrix[i], facets_data[i],
                            facets.face_length)
    return facets

cdef ListOfFaces get_vertices_from_incidence_matrix(tuple matrix):
    r"""
    Initializes vertices in bit-representation as `ListOfFaces`.

    Those are the facets of the polar polyhedron.

    INPUT:

    - `matrix` -- an incidence matrix as tuple of tuples.

    .. NOTE::

        :meth:`sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`
        also has equalities appearing the incidence matrix. Those must be
        removed.

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, vertex_repr_from_bitrep, \
        ....:         get_vertices_from_incidence_matrix
        ....: from sage.geometry.polyhedron.library import polytopes
        ....:
        ....: cdef ListOfFaces vertices
        ....: cdef size_t[14] output
        ....: cdef size_t length
        ....: P = polytopes.permutahedron(4)
        ....: data = P.incidence_matrix()
        ....: rg = range(data.nrows())
        ....: matrix = tuple(tuple(data[i,j] for i in rg)
        ....:                for j in range(data.ncols())
        ....:                if not all(data[i,j] for i in rg))
        ....:
        ....: vertices = get_vertices_from_incidence_matrix(matrix)
        ....: for i in range(vertices.nr_faces):
        ....:     length = vertex_repr_from_bitrep(vertices.data[i], output,
        ....:                                      vertices.face_length)
        ....:     print(tuple(output[j] for j in range(length)))
        ....: ''') # long time
        (3, 5, 9)
        (3, 5, 8)
        (3, 4, 9)
        (3, 7, 8)
        (3, 4, 11)
        (3, 7, 11)
        (5, 6, 9)
        (5, 6, 8)
        (4, 9, 12)
        (1, 7, 8)
        (4, 11, 12)
        (1, 7, 11)
        (6, 9, 10)
        (6, 8, 13)
        (9, 10, 12)
        (1, 8, 13)
        (2, 11, 12)
        (1, 2, 11)
        (0, 6, 10)
        (0, 6, 13)
        (0, 10, 12)
        (0, 1, 13)
        (0, 2, 12)
        (0, 1, 2)
    """
    # expects the incidence matrix to be given as tuple of tuples
    cdef int i
    cdef int j
    cdef tuple newmatrix
    newmatrix = tuple(tuple(matrix[i][j] for i in range(len(matrix)))
                      for j in range(len(matrix[0])))
    return get_facets_from_incidence_matrix(newmatrix)

cdef ListOfFaces get_facets_bitrep_from_facets_tuple(tuple facets_input,
                                                     size_t nr_vertices):
    r"""
    Initializes facets in bit-representation as `ListOfFaces`.

    INPUT:

    - `facets_input` -- a tuple of facets, each facet a tuple of vertices.
      The vertices must be exactly range(nr_vertices).
    - `nr_vertices`.

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, vertex_repr_from_bitrep, \
        ....:         get_facets_bitrep_from_facets_tuple
        ....:
        ....: cdef ListOfFaces facets
        ....: cdef size_t[6] output
        ....: cdef size_t length
        ....:
        ....: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        ....: facets = get_facets_bitrep_from_facets_tuple(bi_pyr, 6)
        ....: for i in range(facets.nr_faces):
        ....:     length = vertex_repr_from_bitrep(facets.data[i], output,
        ....:                                      facets.face_length)
        ....:     print(tuple(output[j] for j in range(length)))
        ....: ''')
        (0, 1, 4)
        (1, 2, 4)
        (2, 3, 4)
        (0, 3, 4)
        (0, 1, 5)
        (1, 2, 5)
        (2, 3, 5)
        (0, 3, 5)
    """
    cdef int i
    cdef ListOfFaces facets = ListOfFaces(len(facets_input),
                                          nr_vertices)
    cdef size_t face_length = facets.face_length
    cdef uint64_t **facets_data = facets.data
    for i in range(len(facets_input)):
        char_from_tuple(facets_input[i], facets_data[i], face_length)
    return facets

cdef ListOfFaces get_vertices_bitrep_from_facets_tuple(tuple facets_input,
                                                       size_t nr_vertices):
    r"""
    Initializes vertices in bit-representation as `ListOfFaces`.

    Those are the facets of the polar polyhedron.

    INPUT:

    - `facets_input` -- a tuple of facets, each facet a tuple of vertices.
      The vertices must be exactly range(nr_vertices).
    - `nr_vertices`.

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, vertex_repr_from_bitrep, \
        ....:         get_vertices_bitrep_from_facets_tuple
        ....:
        ....: cdef ListOfFaces vertices
        ....: cdef size_t[8] output
        ....: cdef size_t length
        ....:
        ....: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        ....: vertices = get_vertices_bitrep_from_facets_tuple(bi_pyr, 6)
        ....: for i in range(vertices.nr_faces):
        ....:     length = vertex_repr_from_bitrep(vertices.data[i], output,
        ....:                                      vertices.face_length)
        ....:     print(tuple(output[j] for j in range(length)))
        ....: ''')
        (0, 3, 4, 7)
        (0, 1, 4, 5)
        (1, 2, 5, 6)
        (2, 3, 6, 7)
        (0, 1, 2, 3)
        (4, 5, 6, 7)
    """
    cdef tuple new_input
    cdef size_t j
    cdef size_t i
    cdef size_t inputlength
    cdef size_t value
    cdef size_t position
    cdef int k
    cdef ListOfFaces vertices = ListOfFaces(nr_vertices,
                                            len(facets_input))
    cdef uint64_t **vertices_data = vertices.data
    cdef size_t face_length = vertices.face_length
    for i in range(nr_vertices):
        for j in range(face_length):
            vertices_data[i][j] = 0
    inputlength = len(facets_input)
    for i in range(inputlength):
        value = i % 64
        position = i//64
        for j in facets_input[i]:
            vertices_data[j][position] += \
                vertex_to_bit_dictionary[value]
    return vertices

cdef size_t vertex_repr_from_bitrep(uint64_t *face, size_t *output,
                                    size_t face_length) except? 0:
    r"""
    Takes a bitrep-representation and converts it to an array.

    Basically this is an inverse to `char_from_tuple`. Instead of a tuple,
    it stores the vertices in `output`. Returns length of representation.

    INPUT::

    - `face` -- a bit-representation of a face.
    - `output` -- an array of `size_t` long enough to contain all vertices of
      that face. `face_length*64` will suffice.
    - `face_length` -- the length of `face`.

    OUTPUT:

    - returns length of `output`
    - stores vertices in `ouput`

    EXAMPLES::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport vertex_repr_from_bitrep
        ....:
        ....: cdef uint64_t[2] face
        ....: cdef size_t[128] output
        ....: cdef size_t length
        ....:
        ....: face[0] = 17
        ....: face[1] = 31
        ....: length = vertex_repr_from_bitrep(face, output, 2)
        ....: print(tuple(output[i] for i in range(length)))
        ....:
        ....: face[0] = 0
        ....: face[1] = 61
        ....: length = vertex_repr_from_bitrep(face, output, 2)
        ....: print(tuple(output[i] for i in range(length)))
        ....: ''')
        (59, 63, 123, 124, 125, 126, 127)
        (122, 123, 124, 125, 127)

    TESTS::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport vertex_repr_from_bitrep, char_from_tuple
        ....: from sage.misc.prandom import randint
        ....:
        ....: cdef uint64_t[2] face
        ....: cdef size_t length
        ....: cdef size_t[128] output
        ....:
        ....: for _ in range(10):
        ....:     st = set(randint(0,127) for i in range(40))
        ....:     tup = tuple(sorted(tuple(st)))
        ....:     char_from_tuple(tup, face, 2)
        ....:     length = vertex_repr_from_bitrep(face, output, 2)
        ....:     tup2 = tuple(output[i] for i in range(length))
        ....:     if not tup == tup2:
        ....:         print('`vertex_repr_from_bitrep` does not behave',
        ....:               'as the inverse of `char_from_tuple`')
        ....: ''')
    """
    cdef size_t i
    cdef size_t j
    cdef size_t counter = 0
    cdef uint64_t copy
    for i in range(face_length):
        if face[i]:
            copy = face[i]
            for j in range(64):
                if copy >= vertex_to_bit_dictionary[j]:
                    output[counter] = i*64 + j
                    counter += 1
                    copy -= vertex_to_bit_dictionary[j]
    return counter

cdef void *aligned_malloc(MemoryAllocator mem, size_t size,
                          size_t align) except NULL:
    r"""
    Allocates `align`-byte aligned memory of size at least `size`.

    Allocates memory of length `size + align - 1` with `mem.malloc`.
    Returns the pointer the the first address that is `align`-byte aligned.

    `align` needs to be a power of two.

    Taken from https://github.com/xjw/cpp/blob/master/cpp/memory_alignment.cpp
    """
    # alignment could not be less than 0
    if (size<0):
        return NULL
    # allocate necessary memory for
    # alignment +
    # area to store the address of memory returned by malloc
    cdef void *p = mem.malloc(size + align-1)
    if (p == NULL):
        raise MemoryError()

    # address of the aligned memory according to the align parameter
    cdef void *ptr = <void *> ((<unsigned long>p + align-1) & ~(align-1))

    # return the address of aligned memory
    return ptr

cdef class ListOfFaces:
    r"""
    A class to store the bit-representation of faces in.

    This class will allocate the memory to store the faces in.

    INPUT:

    - `nr_faces` -- the number of faces to be stored.
    - `nr_vertices` -- the total number of vertices of the Polyhedron.

    .. SEEALSO::

        :meth:`get_facets_from_incidence_matrix`,
        :meth:`get_vertices_from_incidence_matrix`,
        :meth:`get_facets_bitrep_from_facets_tuple`,
        :meth:`get_vertices_bitrep_from_facets_tuple`,
        :class:`FaceIterator`,
        :class:`CombinatorialPolyhedron`.

    EXAMPLES::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces
        ....:
        ....: cdef ListOfFaces facets
        ....:
        ....: facets = ListOfFaces(5, 13)
        ....: print(facets.face_length in (1, 2, 4))
        ....: print(facets.nr_vertices)
        ....: print(facets.nr_faces)
        ....: ''')
        True
        13
        5
    """
    # Those are the definitions done in the header:
    # cdef uint64_t **data
    # cdef MemoryAllocator _mem
    # cdef size_t nr_faces, face_length, nr_vertices

    def __cinit__(self, size_t nr_faces, size_t nr_vertices):
        r"""
        Initializes `ListOfFaces`.

        See :class:`ListOfFaces`.

        TESTS:

        Checking for correct alignment of the data::

            sage: cython('''
            ....: from libc.stdint cimport uint64_t
            ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
            ....: cimport ListOfFaces
            ....:
            ....: cdef ListOfFaces facets
            ....: cdef size_t address
            ....: cdef size_t required_alignment
            ....:
            ....: facets = ListOfFaces(10, 13)
            ....: required_alignment = facets.face_length*8
            ....: for i in range(10):
            ....:     address = <size_t> facets.data[i]
            ....:     if not address == address & ~(required_alignment - 1):
            ....:         print('Alignment not correct')
            ....: ''')
        """
        cdef size_t i
        self.nr_faces = nr_faces
        self.face_length = ((nr_vertices - 1)//chunksize + 1)*chunksize//64
        # `face_length` is the length in terms of `uint64_t
        # WARNING: This needs to be divisible by 2, if chunksize is 128
        #          and divisible by 4, if chunksize is 256.
        self._mem = MemoryAllocator()
        self.data = <uint64_t **> \
            self._mem.malloc(nr_faces * sizeof(uint64_t *))
        self.nr_vertices = nr_vertices
        for i in range(nr_faces):
            self.data[i] = \
                <uint64_t *> aligned_malloc(self._mem, self.face_length*8,
                                            chunksize//8)
            # we must allocate the memory for ListOfFaces overaligned:
            #     must be 16-byte aligned if chunksize = 128
            #     must be 32-byte aligned if chunksize = 256

cdef int calculate_dimension(ListOfFaces faces) except -2:
    r"""
    Calculates the dimension of a polyhedron by its facets.

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, calculate_dimension, \
        ....:         get_facets_bitrep_from_facets_tuple, \
        ....:         get_vertices_bitrep_from_facets_tuple
        ....:
        ....: cdef ListOfFaces facets
        ....: cdef ListOfFaces vertices
        ....:
        ....: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        ....: facets = get_facets_bitrep_from_facets_tuple(bi_pyr, 6)
        ....: vertices = get_vertices_bitrep_from_facets_tuple(bi_pyr, 6)
        ....: print(calculate_dimension(facets))
        ....: print(calculate_dimension(vertices))
        ....: ''')
        3
        3

    ALGORITHM:

    This is done by iteration:

    Calculates the facets of one of the facets (i.e. the ridges contained in
    one of the facets). Then calculates the dimension of the facet, by
    considering its facets.

    Repeats until a face has only one facet. Usually this is a  vertex.

    However, in the unbounded case, this might be different. The face with only
    one facet might be a ray or a line. So the correct dimension of a
    Polyhedron with one facet is the number of "lines/rays/vertices"
    that the facet contains.

    Hence, we know the dimension of a face, which has only one facet and
    iteratively we know the dimension of entire Polyhedron we started with.

    TESTS::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, calculate_dimension, \
        ....:         get_facets_from_incidence_matrix, \
        ....:         get_vertices_from_incidence_matrix
        ....: from sage.geometry.polyhedron.base import Polyhedron
        ....: from sage.misc.prandom import randint
        ....:
        ....: cdef ListOfFaces facets
        ....: cdef ListOfFaces vertices
        ....: for _ in range(10):
        ....:     points = tuple(tuple(randint(-1000,1000) for _ in range(10))
        ....:                    for _ in range(randint(3,15)))
        ....:     P = Polyhedron(vertices=points)
        ....:     data = P.incidence_matrix()
        ....:     rg = range(data.nrows())
        ....:     matrix = tuple(tuple(data[i,j] for i in rg)
        ....:                    for j in range(data.ncols())
        ....:                    if not all(data[i,j] for i in rg))
        ....:
        ....:     d1 = P.dimension()
        ....:     facets = get_facets_from_incidence_matrix(matrix)
        ....:     vertices = get_vertices_from_incidence_matrix(matrix)
        ....:     d2 = calculate_dimension(facets)
        ....:     d3 = calculate_dimension(vertices)
        ....:     if not d1 == d2 == d3:
        ....:         print('calculation_dimension() seems to be incorrect')
        ....: ''') # long time
    """
    cdef size_t nr_faces
    cdef int dim
    nr_faces = faces.nr_faces
    if nr_faces == 0:
        raise TypeError("at least one face needed")

    return calculate_dimension_loop(faces.data, nr_faces, faces.face_length)

cdef int calculate_dimension_loop(uint64_t **facesdata, size_t nr_faces,
                                  size_t face_length) except -2:
    r"""
    Calculates the dimension of a polyhedron by its facets.

    INPUT:

    - `facesdata` -- facets in bit-representation.
    - `nr_faces` -- length of facesdata.
    - `face_length` -- the length of each face in terms of `uint64_t`.

    OUTPUT:

    - Dimension of the Polyhedron.

    .. SEEALSO::

        :meth:`calculate_dimension`
    """
    cdef size_t bitcount, new_nr_faces
    cdef ListOfFaces newfaces
    cdef uint64_t **newfacesdata
    cdef MemoryAllocator newfaces2
    cdef uint64_t **newfaces2data
    cdef int dim
    cdef int returnvalue

    if nr_faces == 0:
            raise TypeError("wrong usage of `calculate_dimension_loop`,\n" +
                            "at least one face needed.")

    if nr_faces == 1:
        # we expect the face to be the empty Polyhedron
        # possibly it contains more than one vertex/rays/lines
        # the dimension of the Polyhedron with this face as only facet is
        # `bitcount`
        bitcount = CountFaceBits(facesdata[0], face_length)
        return int(bitcount)

    # initializing newfaces and newfaces2
    newfaces = ListOfFaces(nr_faces, face_length*64)
    newfaces2 = MemoryAllocator()
    newfacesdata = newfaces.data
    newfaces2data = <uint64_t **> newfaces2.malloc(nr_faces*sizeof(uint64_t *))

    sig_on()
    # newfacesdata will contain all intersections of the form
    # facesdata[i] \cap facesdata[nr_faces -1] with i < nr_faces - 1
    # newfaces2data will contain all inclusion-maximal faces of
    # newfacesdata
    new_nr_faces = get_next_level(facesdata, nr_faces, newfacesdata,
                                  newfaces2data, newfacesdata, 0,
                                  face_length)
    sig_off()
    sig_check()
    return calculate_dimension_loop(newfaces2data, new_nr_faces,
                                    face_length) + 1

cdef class FaceIterator:
    r"""
    A class to iterate over all faces of a Polyhedron.

    INPUT:

    - `facets` -- the facets as `ListOfFaces`.
      If bounded, one can also use the `vertices` and iterate over the
      polar/dual faces. This will be faster if `#vertices < #facets`.
    - `dimension` -- the dimension of the `Polyhedron`.
    - `nr_lines` -- the dimension of the space of all lines the Polyhedron
      contains. If the `Polyhedron is bounded or a cone this is zero.

    .. SEEALSO::

        :meth:`get_facets_from_incidence_matrix`,
        :meth:`get_vertices_from_incidence_matrix`,
        :meth:`get_facets_bitrep_from_facets_tuple`,
        :meth:`get_vertices_bitrep_from_facets_tuple`,
        :class:`ListOfFaces`,
        :class:`CombinatorialPolyhedron`.

    EXAMPLES::

        sage: cython('''
        ....: from __future__ import print_function
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport FaceIterator, ListOfFaces, \
        ....:         get_facets_bitrep_from_facets_tuple
        ....:
        ....: cdef ListOfFaces facets
        ....: cdef size_t *output_vertex_repr
        ....: cdef size_t *output_facet_repr
        ....: cdef size_t leng
        ....: cdef int d
        ....: cdef FaceIterator face_iter
        ....:
        ....: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        ....: facets = get_facets_bitrep_from_facets_tuple(bi_pyr, 6)
        ....: face_iter = FaceIterator(facets, 3, 0)
        ....: output_vertex_repr = face_iter.get_output1_array()
        ....: output_facet_repr = face_iter.get_output2_array()
        ....:
        ....: d = face_iter.next_face()
        ....: while d < 3:
        ....:     print('Next face:')
        ....:
        ....:     print('dimension is', d)
        ....:
        ....:     leng = face_iter.length_vertex_repr()
        ....:     print('nr_vertices it contains:', leng)
        ....:
        ....:     leng = face_iter.vertex_repr(output_vertex_repr)
        ....:     tup = tuple(output_vertex_repr[i] for i in range(leng))
        ....:     print('vertex representation:', tup)
        ....:
        ....:     leng = face_iter.facet_repr(output_facet_repr)
        ....:     tup = tuple(output_facet_repr[i] for i in range(leng))
        ....:     print('facet representation:', tup)
        ....:
        ....:     d = face_iter.next_face()
        ....: ''')
        Next face:
        dimension is 2
        nr_vertices it contains: 3
        vertex representation: (0, 3, 5)
        facet representation: (7,)
        Next face:
        dimension is 2
        nr_vertices it contains: 3
        vertex representation: (2, 3, 5)
        facet representation: (6,)
        Next face:
        dimension is 2
        nr_vertices it contains: 3
        vertex representation: (1, 2, 5)
        facet representation: (5,)
        Next face:
        dimension is 2
        nr_vertices it contains: 3
        vertex representation: (0, 1, 5)
        facet representation: (4,)
        Next face:
        dimension is 2
        nr_vertices it contains: 3
        vertex representation: (0, 3, 4)
        facet representation: (3,)
        Next face:
        dimension is 2
        nr_vertices it contains: 3
        vertex representation: (2, 3, 4)
        facet representation: (2,)
        Next face:
        dimension is 2
        nr_vertices it contains: 3
        vertex representation: (1, 2, 4)
        facet representation: (1,)
        Next face:
        dimension is 2
        nr_vertices it contains: 3
        vertex representation: (0, 1, 4)
        facet representation: (0,)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (3, 5)
        facet representation: (6, 7)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (0, 5)
        facet representation: (4, 7)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (0, 3)
        facet representation: (3, 7)
        Next face:
        dimension is 0
        nr_vertices it contains: 1
        vertex representation: (5,)
        facet representation: (4, 5, 6, 7)
        Next face:
        dimension is 0
        nr_vertices it contains: 1
        vertex representation: (3,)
        facet representation: (2, 3, 6, 7)
        Next face:
        dimension is 0
        nr_vertices it contains: 1
        vertex representation: (0,)
        facet representation: (0, 3, 4, 7)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (2, 5)
        facet representation: (5, 6)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (2, 3)
        facet representation: (2, 6)
        Next face:
        dimension is 0
        nr_vertices it contains: 1
        vertex representation: (2,)
        facet representation: (1, 2, 5, 6)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (1, 5)
        facet representation: (4, 5)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (1, 2)
        facet representation: (1, 5)
        Next face:
        dimension is 0
        nr_vertices it contains: 1
        vertex representation: (1,)
        facet representation: (0, 1, 4, 5)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (0, 1)
        facet representation: (0, 4)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (3, 4)
        facet representation: (2, 3)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (0, 4)
        facet representation: (0, 3)
        Next face:
        dimension is 0
        nr_vertices it contains: 1
        vertex representation: (4,)
        facet representation: (0, 1, 2, 3)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (2, 4)
        facet representation: (1, 2)
        Next face:
        dimension is 1
        nr_vertices it contains: 2
        vertex representation: (1, 4)
        facet representation: (0, 1)
    """
    # Those are the definitions done in the header:
    # cdef uint64_t *face
    # cdef int current_dimension, dimension, record_dimension, lowest_dimension
    # cdef MemoryAllocator _mem
    # cdef uint64_t ***newfaces2
    # cdef tuple newfaces_lists
    # cdef uint64_t ***newfaces
    # cdef uint64_t **forbidden
    # cdef size_t *nr_faces
    # cdef size_t *nr_forbidden
    # cdef int *first_time
    # cdef size_t yet_to_yield, face_length, nr_facets, nr_vertices
    # cdef size_t *output1
    # cdef size_t *output2
    # cdef int nr_lines

    def __cinit__(self, ListOfFaces facets, int dimension, int nr_lines):
        r"""
        Initializes `FaceIterator`.

        See :class:`FaceIterator`.

        EXAMPLES::

            sage: cython('''
            ....: from __future__ import print_function
            ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
            ....: cimport FaceIterator, ListOfFaces, \
            ....:         get_facets_from_incidence_matrix
            ....: from sage.geometry.polyhedron.library import polytopes
            ....:
            ....: cdef ListOfFaces facets
            ....: cdef int d
            ....: cdef FaceIterator face_iter
            ....:
            ....: P = polytopes.permutahedron(4)
            ....: data = P.incidence_matrix()
            ....: rg = range(data.nrows())
            ....: matrix = tuple(tuple(data[i,j] for i in rg)
            ....:                for j in range(data.ncols())
            ....:                if not all(data[i,j] for i in rg))
            ....: facets = get_facets_from_incidence_matrix(matrix)
            ....:
            ....: face_iter = FaceIterator(facets, 3, 0)
            ....: f_vector = [1, 0, 0, 0, 1]
            ....:
            ....: d = face_iter.next_face()
            ....: while d < 3:
            ....:     f_vector[d+1] += 1
            ....:     d = face_iter.next_face()
            ....: print ('f_vector of permutahedron(4): ', f_vector)
            ....: ''') # long time
            f_vector of permutahedron(4):  [1, 24, 36, 14, 1]
        """
        cdef int i
        cdef ListOfFaces some_list
        if dimension <= 0:
            raise TypeError("FaceIterator expects positive dimensions")
        if facets.nr_faces <= 0:
            raise TypeError("FaceIterator expects non-empty `ListOfFaces`")
        self.dimension = dimension
        self.current_dimension = dimension - 1
        self.face_length = facets.face_length
        self._mem = MemoryAllocator()

        # in each dimension we need to know the number of facets in newfaces2
        self.nr_faces = \
            <size_t *> self._mem.malloc(dimension * sizeof(size_t))
        self.nr_facets = facets.nr_faces
        self.nr_faces[dimension - 1] = self.nr_facets

        # in each dimension we need to know how many faces are in forbidden
        # if we have just visited all faces of `face`, then we need only `face`
        # in forbidden to know, that we shall not visit any of them again
        self.nr_forbidden = \
            <size_t *> self._mem.malloc(dimension * sizeof(size_t))
        self.nr_forbidden[dimension -1] = 0

        # the place where the new facets will be stored
        self.newfaces2 = <uint64_t ***> \
            self._mem.malloc(dimension * sizeof(uint64_t **))
        for i in range(dimension - 1):
            self.newfaces2[i] = \
                <uint64_t **> self._mem.malloc(self.nr_facets *
                                               sizeof(uint64_t *))
        self.newfaces2[dimension - 1] = facets.data

        # the place where the possible new faces are being stored
        # note that they are actually located here, whereas newfaces2 only
        # points to the correct ones
        self.newfaces_lists = \
            tuple(ListOfFaces(facets.nr_faces, facets.nr_vertices)
                  for i in range(dimension -1))
        self.newfaces = <uint64_t ***> \
            self._mem.malloc((dimension -1) * sizeof(uint64_t **))
        for i in range(dimension -1):
            some_list = self.newfaces_lists[i]
            self.newfaces[i] = some_list.data

        # forbidden will point to the faces, we shall not visit again
        self.forbidden = <uint64_t **> \
            self._mem.malloc(facets.nr_faces*sizeof(uint64_t *))

        # at each step we first want to yield all facets, before considering
        # faces farther down the lattices
        self.yet_to_yield = facets.nr_faces

        # the user can set this to be a different number, such that only faces
        # of a fixed dimension are yield
        self.record_dimension = -2

        # we will not yield the empty face, if there are `n` lines, than there
        # are no faces below dimension `n`
        self.lowest_dimension = nr_lines
        self.nr_lines = nr_lines

        # we can only store a face `face1` in forbidden,
        # after we have visited all its faces,
        # when we consider the next face `face2`, we have to add `face1` to
        # forbidden, the only way to know, is to remember, that we have been
        # here before and there is actually a face more in `newfaces2`, than
        # what the counter seems to say
        self.first_time = \
            <int *> self._mem.malloc(dimension * sizeof(int))
        self.first_time[dimension - 1] = 1

        self.output1 = NULL
        self.output2 = NULL
        self.nr_vertices = facets.nr_vertices

    cdef void set_record_dimension(self, int dim):
        r"""
        This will have the iterator only yield faces in dimension `dim`.

        EXAMPLES::

            sage: cython('''
            ....: from __future__ import print_function
            ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
            ....: cimport FaceIterator, ListOfFaces, \
            ....:         get_facets_from_incidence_matrix
            ....: from sage.geometry.polyhedron.library import polytopes
            ....:
            ....: cdef ListOfFaces facets
            ....: cdef FaceIterator face_iter
            ....:
            ....: P = polytopes.permutahedron(5)
            ....: data = P.incidence_matrix()
            ....: rg = range(data.nrows())
            ....: matrix = tuple(tuple(data[i,j] for i in rg)
            ....:                for j in range(data.ncols())
            ....:                if not all(data[i,j] for i in rg))
            ....: facets = get_facets_from_incidence_matrix(matrix)
            ....:
            ....: face_iter = FaceIterator(facets, 4, 0)
            ....: twofaces = 0
            ....:
            ....: face_iter.set_record_dimension(2)
            ....: while face_iter.next_face() == 2:
            ....:     twofaces += 1
            ....: print ('permutahedron(5) has', twofaces,
            ....:        'faces of dimension 2')
            ....: ''') # long time
            permutahedron(5) has 150 faces of dimension 2
        """
        self.record_dimension = dim
        self.lowest_dimension = max(self.nr_lines, dim)

    cdef inline int next_face(self) except -1:
        r"""
        Sets `FaceIterator.face` to the next face and returns the dimension.

        Returns the dimension of the `Polyhedron` on failure.

        The function calls `self.next_face_loop` until it set a new face or
        until it is consumed.

        **** Messing with the face_iterator *****
        suppose face_iterator returns `face` and you do not want
        to visit any farther faces of `face` you can do the following:

        self.forbidden[self.face_iterator_nr_forbidden] = self.face
        face_iterator_nr_forbidden += 1

        This prevents any faces of `face` of appearing in the face iterator.
        """
        self.face = NULL
        cdef int dim = self.dimension
        while (not self.next_face_loop()) and (self.current_dimension < dim):
            sig_check()
        return self.current_dimension

    cdef inline int next_face_loop(self) except -1:
        r"""
        Sets `self._face` to the next face. On success returns 1. Otherwise 0.
        Need to be recalled then.

        If `self.current_dimension == self.dimension`, then the iterator is
        consumed.
        """
        cdef size_t nr_faces
        cdef size_t nr_forbidden
        cdef uint64_t **faces
        cdef size_t i, newfacescounter
        if unlikely(self.current_dimension == self.dimension):
            # the function is not supposed to be called in this case
            # just to prevent it from crashing
            raise ValueError("calling a next_face_loop, " +
                             "but iterator is consumed")

        nr_faces = self.nr_faces[self.current_dimension]
        nr_forbidden = self.nr_forbidden[self.current_dimension]
        faces = self.newfaces2[self.current_dimension]

        if (self.record_dimension > -2) and \
                (self.record_dimension != self.current_dimension):
            # if we are not in dimension `record_dimension`,
            # then we should not yield any faces
            # (in case `face_iterator_dimension == -2` we yield all faces)
            self.yet_to_yield = 0

        if self.yet_to_yield:
            # return the next face
            self.yet_to_yield -= 1
            self.face = faces[self.yet_to_yield]
            return 1

        if self.current_dimension == self.lowest_dimension:
            # we will not yield the empty face
            # we will not yield below what is wanted
            self.current_dimension += 1
            return 0

        if nr_faces <= 1:
            # there will be no more faces from intersections
            self.current_dimension += 1
            return 0

        i = nr_faces - 1
        self.nr_faces[self.current_dimension] -= 1
        if not self.first_time[self.current_dimension]:
            # if there exists faces[i+1], we have visited all its faces already
            # hence we should not visit any of them again
            self.forbidden[nr_forbidden] = faces[i+1]
            self.nr_forbidden[self.current_dimension] += 1
            nr_forbidden = self.nr_forbidden[self.current_dimension]

        else:
            # we will visit all the faces of faces[nr_faces]
            # once we have done so, we want to add this face to forbidden
            self.first_time[self.current_dimension] = 0

        # get the facets contained in faces[i] but not in any of the forbidden
        sig_check()
        sig_on()
        newfacescounter = get_next_level(
            faces, i+1, self.newfaces[self.current_dimension-1],
            self.newfaces2[self.current_dimension-1],
            self.forbidden, nr_forbidden, self.face_length)
        sig_off()
        if newfacescounter:
            # if there are new faces contained in faces[i+1],
            # we will set up the variables to correctly visit them on the next
            # call of `next_face_loop`
            self.current_dimension -= 1
            self.first_time[self.current_dimension] = 1
            self.nr_faces[self.current_dimension] = newfacescounter
            self.nr_forbidden[self.current_dimension] = nr_forbidden
            self.yet_to_yield = newfacescounter
            return 0

        else:
            # if there are no faces in lower dimension,
            # then there is no need to add the face to forbidden
            # this might become important when calculating simpliness
            # and simpliciality, where we will mess with the iterator
            # and add some faces to forbidden in order to not consider subfaces
            self.first_time[self.current_dimension] = 1
            return 0

    cdef size_t length_vertex_repr(self) except -1:
        r"""
        Calculates the number of vertices in the current face by counting the
        number of set bits.

        .. SEEALSO::

            :class:`ListOfFaces`,
        """
        if self.face:
            return CountFaceBits(self.face, self.face_length)

        # the face was not initialized properly
        raise LookupError("`FaceIterator` does not point to a face")

    cdef size_t facet_repr(self, size_t *output) except -1:
        r"""
        Writes facet_repr of the current face in output. Returns its length.

        Output must have length of that facet representation.
        Usually one should allocate output to be of length
        `nr_facets` of the `Polyhedron`.

        The method :meth:`ListOfFaces.get_output2_array` allocates an array of
        correct size.

        .. SEEALSO::

            :class:`ListOfFaces`,
            :meth:`ListOfFaces.get_output2_array`
        """
        cdef size_t nr_facets = self.nr_facets
        cdef uint64_t **facets = self.newfaces2[self.dimension - 1]
        cdef size_t face_length = self.face_length
        return facet_repr_from_bitrep(self.face, facets, output,
                                      nr_facets, face_length)

    cdef size_t vertex_repr(self, size_t *output) except -1:
        r"""
        Writes vertex_repr of the current face in output. Returns its length.

        Output must have length of that vertex representation.
        Usually one should allocate output to be of length
        `nr_vertices` of the `Polyhedron`.

        The method :meth:`ListOfFaces.get_output1_array` allocates an array of
        correct size.
        """
        cdef size_t face_length = self.face_length
        return vertex_repr_from_bitrep(self.face, output, face_length)

    cdef size_t *get_output1_array(self) except NULL:
        r"""
        Allocates array to store vertex_repr of a face in. Returns a pointer.

        A FaceIterator will have only one such array.
        """
        if self.output1 is NULL:
            self.output1 = \
                <size_t *> self._mem.malloc(self.nr_vertices *
                                            sizeof(size_t))
        return self.output1

    cdef size_t *get_output2_array(self) except NULL:
        r"""
        Allocates array to store facet_repr of a face in. Returns a pointer.

        A FaceIterator will have only one such array.
        """
        if self.output2 is NULL:
            self.output2 = \
                <size_t *> self._mem.malloc(self.nr_facets *
                                            sizeof(size_t))
        return self.output2

cdef class ListOfAllFaces:
    r"""
    A class to store all faces of :class:`CombinatorialPolyhedron`.

    Once all faces are added, they can be sorted and then the incidences
    for the `face_lattice` can be generated.

    INPUT:

    - `facets` -- the facets of the `CombinatorialPolyhedron` in
      Bit-representation. In the polar case also the vertices can be given.
    - `dimension` -- the dimension of the `CombinatorialPolyhedron`.
    - `f_vector` -- the f_vector of the `CombinatorialPolyhedron`.
      If in the polar case the vertices of the `Polyhedron` are given as
      `facets`, then the `f_vector` must be flipped.

    .. SEEALSO::

        :meth:`CombinatorialPolyhedron._record_all_faces`,
        :meth:`CombinatorialPolyhedron._record_all_faces_helper`,
        :meth:`CombinatorialPolyhedron.face_lattice`,
        :meth:`CombinatorialPolyhedron._calculate_face_lattice_incidences`.

    EXAMPLES::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfAllFaces, FaceIterator, ListOfFaces, \
        ....:         get_facets_bitrep_from_facets_tuple
        ....:
        ....: cdef ListOfFaces facets
        ....: cdef FaceIterator face_iter
        ....: cdef ListOfAllFaces all_faces
        ....:
        ....: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        ....: facets = get_facets_bitrep_from_facets_tuple(bi_pyr, 6)
        ....: face_iter = FaceIterator(facets, 3, 0)
        ....: all_faces = ListOfAllFaces(facets, 3, (1, 6, 12, 8, 1))
        ....:
        ....: d = face_iter.next_face()
        ....: while d < 3:
        ....:     if d < 2:
        ....:         all_faces.add_face(d, face_iter.face)
        ....:     d = face_iter.next_face()
        ....: print('ListOfAllFaces initialized')
        ....: all_faces.sort()
        ....: print('ListOfAllFaces sorted')
        ....: ''')
        ListOfAllFaces initialized
        ListOfAllFaces sorted
    """
    # those are the defintions from the header
    # cdef tuple lists_facet_repr #might need this for flag-vector
    # cdef MemoryAllocator _mem
    # cdef tuple lists_vertex_repr
    # cdef size_t nr_facets
    # cdef size_t nr_vertices
    # cdef size_t face_length_vertex
    # cdef size_t face_length_facet #might need this for flag-vector
    # cdef int dimension
    # cdef size_t *face_counter
    # cdef size_t *f_vector
    # cdef int is_sorted
    # cdef size_t *output1
    # cdef size_t *output2
    # cdef int incidence_dim_one
    # cdef int incidence_dim_two
    # cdef size_t incidence_counter_one
    # cdef size_t incidence_counter_two
    # cdef ListOfFaces incidence_face_mem
    # cdef uint64_t *incidence_face
    # cdef uint64_t **facets
    # cdef int incidence_is_initialized
    # cdef uint64_t ***data_vertex

    def __cinit__(self, ListOfFaces facets, int dimension, f_vector):
        r"""
        Initializes `ListOfAllFaces`.

        See :class:`ListOfAllFaces`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C._record_all_faces() # indirect doctests
            sage: C.face_lattice(vertices=True, facets=True, names=True)
            Finite lattice containing 28 elements
        """
        cdef int i
        cdef size_t j
        cdef ListOfFaces some_list
        self._mem = MemoryAllocator()
        self.nr_facets = facets.nr_faces
        self.nr_vertices = facets.nr_vertices
        self.face_length_vertex = facets.face_length
        self.dimension = dimension
        # keeping the vertex_repr of each dimension in a tuple
        self.lists_vertex_repr = \
            tuple(ListOfFaces(f_vector[i+1], self.nr_vertices)
                  for i in range(-1, dimension-1))
        self.lists_vertex_repr += (facets,)
        self.lists_vertex_repr += \
            (ListOfFaces(1, self.nr_vertices),)
        # setting up a pointer for direct access to the data
        self.data_vertex = <uint64_t ***> \
            self._mem.malloc((dimension + 2)*(sizeof(uint64_t **)))
        for i in range(dimension + 2):
            some_list = self.lists_vertex_repr[i]
            self.data_vertex[i] = some_list.data
        # initialize the empty face
        char_from_tuple((), self.data_vertex[0][0], self.face_length_vertex)
        # intialize the full polyhedron
        char_from_tuple(tuple(j for j in range(self.nr_vertices)),
                        self.data_vertex[dimension + 1][0],
                        self.face_length_vertex)
        # face_counter keeps track, if all faces have been added already
        self.face_counter = \
            <size_t *> self._mem.malloc((dimension + 2) *
                                        sizeof(size_t))
        self.face_counter[0] = 1
        self.face_counter[dimension + 1] = 1
        self.face_counter[dimension] = self.nr_facets
        for i in range(1, dimension):
            self.face_counter[i] = 0
        # copy f_vector for later use
        self.f_vector = <size_t *> self._mem.malloc((dimension + 2) *
                                                    sizeof(size_t))
        for i in range(dimension + 2):
            self.f_vector[i] = f_vector[i]
        self.is_sorted = 0
        self.output1 = NULL
        self.output2 = NULL
        self.facets = facets.data
        self.incidence_is_initialized = 0

    cdef int add_face(self, int face_dim, uint64_t *face) except 0:
        r"""
        Adds a face to `ListOfAllFaces`.

        EXAMPLES::

            sage: cython('''
            ....: from libc.stdint cimport uint64_t
            ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
            ....: cimport ListOfAllFaces, FaceIterator, ListOfFaces, \
            ....:         get_vertices_bitrep_from_facets_tuple
            ....:
            ....: cdef ListOfFaces vertices
            ....: cdef FaceIterator face_iter
            ....: cdef ListOfAllFaces all_faces
            ....:
            ....: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            ....: vertices = get_vertices_bitrep_from_facets_tuple(bi_pyr, 6)
            ....: face_iter = FaceIterator(vertices, 3, 0)
            ....: all_faces = ListOfAllFaces(vertices, 3, (1, 8, 12, 6, 1))
            ....:
            ....: d = face_iter.next_face()
            ....: try:
            ....:     all_faces.add_face(d, face_iter.face)
            ....: except:
            ....:     print('facets are initialized already')
            ....: while d < 3:
            ....:     if d < 2:
            ....:         all_faces.add_face(d, face_iter.face)
            ....:     d = face_iter.next_face()
            ....: print('ListOfAllFaces initialized')
            ....: face_iter = FaceIterator(vertices, 3, 0)
            ....: d = face_iter.next_face()
            ....: try:
            ....:     all_faces.add_face(d, face_iter.face)
            ....: except:
            ....:     print('trying to add too many faces to `ListOfAllFaces`')
            ....: ''')
            facets are initialized already
            ListOfAllFaces initialized
            trying to add too many faces to `ListOfAllFaces`
        """
        cdef size_t counter = self.face_counter[face_dim + 1]
        cdef size_t max_number = self.f_vector[face_dim + 1]
        if unlikely(counter >= max_number):
            raise IOError("trying to add too many faces to `ListOfAllFaces`")
        memcpy(self.data_vertex[face_dim + 1][counter], face,
               self.face_length_vertex*8)
        self.face_counter[face_dim + 1] += 1
        return 1

    cdef int sort(self) except 0:
        r"""
        Sorts the list faces in vertex-representation (except for facets).

        This way one can fastly find a certain face in the list later.

        EXAMPLES::

            sage: cython('''
            ....: from libc.stdint cimport uint64_t
            ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
            ....: cimport ListOfAllFaces, FaceIterator, ListOfFaces, \
            ....:         get_facets_bitrep_from_facets_tuple
            ....:
            ....: cdef ListOfFaces facets
            ....: cdef FaceIterator face_iter
            ....: cdef ListOfAllFaces all_faces
            ....:
            ....: pyr = ((0,1,5), (1,2,5), (2,3,5), (3,4,5), (4,0,5), (0,1,2,3,4))
            ....: facets = get_facets_bitrep_from_facets_tuple(pyr, 6)
            ....: face_iter = FaceIterator(facets, 3, 0)
            ....: all_faces = ListOfAllFaces(facets, 3, (1, 6, 10, 6, 1))
            ....:
            ....: try:
            ....:     all_faces.sort()
            ....: except:
            ....:     print('`ListOfAllFaces` does not contain all faces')
            ....:
            ....: d = face_iter.next_face()
            ....: while d < 3:
            ....:     if d < 2:
            ....:         all_faces.add_face(d, face_iter.face)
            ....:     d = face_iter.next_face()
            ....: print('ListOfAllFaces initialized')
            ....: all_faces.sort()
            ....: print('ListOfAllFaces sorted')
            ....: ''')
            `ListOfAllFaces` does not contain all faces
            ListOfAllFaces initialized
            ListOfAllFaces sorted
        """
        cdef int dim = self.dimension
        cdef int i
        if unlikely(self.is_sorted):
            return 1
        for i in range(dim + 2):
            if unlikely(self.f_vector[i] != self.face_counter[i]):
                print (tuple(i, self.f_vector[i], self.face_counter[i],
                             i+1, self.f_vector[i+1], self.face_counter[i+1]))
                raise ValueError("`ListOfAllFaces` does not contain all faces")
        for i in range(0, dim):
            self._sort_one_list(self.data_vertex[i], self.f_vector[i])
        self.is_sorted = 1
        return 1

    cdef int _sort_one_list(self, uint64_t **faces, size_t nr_faces) except 0:
        r"""
        Sorts `faces` of length `nr_faces`. Each face in `faces`
        is supposed to be in vertex-repr, i.e. of length `face_length_vertex`.

        See :meth:`sort`.
        """
        cdef MemoryAllocator mem = MemoryAllocator()
        cdef uint64_t **extra_mem = \
            <uint64_t **> mem.malloc(nr_faces*sizeof(uint64_t *))
        self._sort_one_list_loop(faces, faces, extra_mem, nr_faces)
        return 1

    cdef int _sort_one_list_loop(
            self, uint64_t **inp, uint64_t **output1,
            uint64_t **output2, size_t nr_faces) except 0:
        r"""
        This is mergesort.

        Sorts `inp` and returns it in `output1`.

        BEWARE: Input is the same as output1 or output2

        See :meth:`sort`.
        """
        cdef size_t middle = nr_faces//2
        cdef size_t other = nr_faces - middle
        cdef size_t i = 0
        cdef size_t j = middle
        cdef size_t counter = 0
        if nr_faces == 1:
            output1[0] = inp[0]
            return 1

        # sort the upper and lower half of inp iteratively into `output2`
        self._sort_one_list_loop(inp, output2, output1, middle)
        self._sort_one_list_loop(&(inp[middle]), &(output2[middle]),
                                 &(output1[middle]), other)
        while i < middle and j < nr_faces:
            if self.is_smaller(output2[i], output2[j]):
                output1[counter] = output2[i]
                i += 1
                counter += 1
            else:
                output1[counter] = output2[j]
                j += 1
                counter += 1
        if i < middle:
            while i < middle:
                output1[counter] = output2[i]
                i += 1
                counter += 1
        else:
            while j < nr_faces:
                output1[counter] = output2[j]
                j += 1
                counter += 1
        return 1

    cdef inline size_t find_face(self, int dimension, uint64_t *face) except -1:
        r"""
        Returns the index of `face`, if it is of dimension `dimension`.
        Assumes `face` in vertex-representation.

        .. NOTE::

            Will give an index no matter if `face` is actual of dimension
            `dimension`. Check the result with belows `is_equal`.

        EXAMPLES::

            sage: cython('''
            ....: from libc.stdint cimport uint64_t
            ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
            ....: cimport ListOfAllFaces, FaceIterator, ListOfFaces, \
            ....:         get_vertices_bitrep_from_facets_tuple
            ....:
            ....: cdef ListOfFaces vertices
            ....: cdef FaceIterator face_iter
            ....: cdef ListOfAllFaces all_faces
            ....: cdef size_t position
            ....:
            ....: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            ....: vertices = get_vertices_bitrep_from_facets_tuple(bi_pyr, 6)
            ....: face_iter = FaceIterator(vertices, 3, 0)
            ....: all_faces = ListOfAllFaces(vertices, 3, (1, 8, 12, 6, 1))
            ....:
            ....: d = face_iter.next_face()
            ....: while d < 3:
            ....:     if d < 2:
            ....:         all_faces.add_face(d, face_iter.face)
            ....:     d = face_iter.next_face()
            ....: print('ListOfAllFaces initialized')
            ....: try:
            ....:     all_faces.find_face(2, vertices.data[0])
            ....: except:
            ....:     print('`ListOfAllFaces` needs to be sorted first')
            ....: all_faces.sort()
            ....: try:
            ....:     all_faces.find_face(2, vertices.data[0])
            ....: except:
            ....:     print('cannot find a facet, as those are not sorted')
            ....: position = all_faces.find_face(1, vertices.data[0])
            ....: print(all_faces.is_equal(1, position, vertices.data[0]))
            ....: ''')
            ListOfAllFaces initialized
            `ListOfAllFaces` needs to be sorted first
            cannot find a facet, as those are not sorted
            0
        """
        if unlikely(not self.is_sorted):
            raise ValueError("`ListOfAllFaces` needs to be sorted first")
        if unlikely(dimension == self.dimension -1):
            raise ValueError("cannot find a facet, as those are not sorted")
            # of course one can easily add a function to search for a facet as
            # well, but there seems to be no need for that
        if unlikely(dimension < -1 or dimension > self.dimension):
            raise IndexError("dimension out of range")
        cdef size_t start = 0
        cdef size_t middle
        cdef nr_faces = self.f_vector[dimension + 1]
        cdef uint64_t **list_faces = self.data_vertex[dimension + 1]
        while (nr_faces > 1):
            middle = nr_faces//2
            if self.is_smaller(face, list_faces[middle + start]):
                nr_faces = middle
            else:
                nr_faces -= middle
                start += middle
        return start

    cdef inline int is_smaller(self, uint64_t *one, uint64_t *two):
        r"""
        Returns 1 if `one` is smaller than `two`, otherwise 0.
        Expects `one` and `two` to be in vertex-representation.
        """
        return memcmp(one, two, self.face_length_vertex*8) < 0

    cdef inline int is_equal(self, int dimension, size_t index,
                             uint64_t *face) except -1:
        r"""
        Checks wether `face` is in the list with dimension `dimension`
        and index `index`.
        """
        if unlikely(dimension < -1 or dimension > self.dimension
                    or index >= self.f_vector[dimension + 1]):
            raise IndexError()
        cdef uint64_t * face2 = self.data_vertex[dimension+1][index]
        cdef size_t i
        cdef size_t length = self.face_length_vertex
        return (0 == memcmp(face, face2,
                            length*8))

    cdef size_t facet_repr(self, int dimension, size_t index,
                           size_t *output) except -1:
        r"""
        Writes the facet_repr of the face of dimension `dimension`
        and index `index` in `output`. Returns its length.

        Output must have length of that facet representation.
        Usually one should allocate output to be of length
        `nr_facets` of the `Polyhedron`.

        The method :meth:`ListOfAllFaces.get_output2_array` allocates
        an array of correct size.

        .. SEEALSO::

            :class:`ListOfAllFaces`,
            :meth:`ListOfAllFaces.get_output2_array`
        """
        cdef size_t nr_facets = self.nr_facets
        cdef ListOfFaces faces = self.lists_vertex_repr[dimension + 1]
        cdef ListOfFaces facets = self.lists_vertex_repr[self.dimension]
        cdef size_t face_length = self.face_length_vertex
        cdef uint64_t **facesdata = faces.data
        cdef uint64_t *face = facesdata[index]
        return facet_repr_from_bitrep(face, facets.data, output,
                                      nr_facets, face_length)

    cdef size_t vertex_repr(self, int dimension, size_t index,
                            size_t *output) except -1:
        r"""
        Writes the vertex_repr of the face of dimension `dimension`
        and index `index` in `output`. Returns its length.

        Output must have length of that vertex representation.
        Usually one should allocate output to be of length
        `nr_vertices` of the `Polyhedron`.

        The method :meth:`ListOfAllFaces.get_output1_array` allocates
        an array of correct size.

        .. SEEALSO::

            :class:`ListOfAllFaces`,
            :meth:`ListOfAllFaces.get_output1_array`
        """
        cdef size_t face_length = self.face_length_vertex
        cdef ListOfFaces faces = self.lists_vertex_repr[dimension + 1]
        cdef uint64_t **facesdata = faces.data
        cdef uint64_t *face = facesdata[index]
        return vertex_repr_from_bitrep(face, output, face_length)

    cdef size_t *get_output1_array(self) except NULL:
        r"""
        Allocates array to store vertex_repr of a face in. Returns a pointer.

        A FaceIterator will have only one such array.
        """
        if self.output1 is NULL:
            self.output1 = \
                <size_t *> self._mem.malloc(self.nr_vertices *
                                            sizeof(size_t))
        return self.output1

    cdef size_t *get_output2_array(self) except NULL:
        r"""
        Allocates array to store facet_repr of a face in. Returns a pointer.

        A FaceIterator will have only one such array.
        """
        if self.output2 is NULL:
            self.output2 = \
                <size_t *> self._mem.malloc(self.nr_facets *
                                            sizeof(size_t))
        return self.output2

    cdef void incidence_init(self, int dimension_one, int dimension_two):
        r"""
        Initialized the ListOfAllFaces to give incidences between
        `dimension_one` and `dimension_two`.

        Currently only `dimension_one == dimension_two + 1` is properly
        implmented. This is the relevant case to get all incidences for the
        `hasse_diagram`/`face_lattice`.

        :meth:`next_incidence`.
        """
        cdef size_t nr_facets = self.nr_facets
        cdef size_t i
        if dimension_one == self.dimension:
            # the entire polyhedron is incident to every face
            if dimension_two < -1:
                raise ValueError("no faces of dimension %s"%dimension_two)
            if dimension_two > self.dimension:
                raise ValueError("no faces of dimension %s"%dimension_two)
            self.incidence_dim_one = dimension_one
            self.incidence_dim_two = dimension_two
            self.incidence_counter_one = 0
            self.incidence_counter_two = 0
            self.incidence_is_initialized = 2
            return
        if dimension_two == -1:
            # the entire polyhedron is incident to every face
            if dimension_one < -1:
                raise ValueError("no faces of dimension %s"%dimension_two)
            if dimension_one > self.dimension:
                raise ValueError("no faces of dimension %s"%dimension_two)
            self.incidence_dim_one = dimension_one
            self.incidence_dim_two = dimension_two
            self.incidence_counter_one = 0
            self.incidence_counter_two = 0
            self.incidence_is_initialized = 3
            return
        if dimension_one != dimension_two + 1:
            raise ValueError("`dimension_one` = `dimension_two` + 1 must hold")
            # we give this function in more generality,
            # so that we can later calculate more than just incidences of
            # neighbor-dimensions
        if dimension_one > self.dimension:
            raise ValueError("no faces of dimension %s"%dimension_one)
        if dimension_two < -1:
            raise ValueError("no faces of dimension %s"%dimension_two)
        if not self.is_sorted:
            raise ValueError("Allfaces need to be sorted with sort() yet")
        if not self.incidence_face_mem:
            self.incidence_face_mem = \
                ListOfFaces(1, self.nr_vertices)
        self.incidence_face = self.incidence_face_mem.data[0]
        self.incidence_dim_one = dimension_one
        self.incidence_dim_two = dimension_two
        self.incidence_counter_one = 0
        self.incidence_counter_two = 0
        self.incidence_is_initialized = 1

    cdef inline int next_trivial_incidence(self, size_t *one, size_t *two):
        r"""
        Handling the case where dimension_one is the full polyhedron.

        See :meth:`next_incidence`.
        """
        one[0] = 0
        two[0] = self.incidence_counter_two
        self.incidence_counter_two += 1
        return (two[0] < self.f_vector[self.incidence_dim_two + 1])

    cdef inline int next_trivial_incidence2(self, size_t *one, size_t *two):
        r"""
        Handling the case where dimension_two is the empty face

        See :meth:`next_incidence`.
        """
        two[0] = 0
        one[0] = self.incidence_counter_one
        self.incidence_counter_one += 1
        return (one[0] < self.f_vector[self.incidence_dim_one + 1])

    cdef int next_incidence(self, size_t *one, size_t *two):
        r"""
        Sets `one[0]` and `two[0]` to be the next incidence. Returns 1 until
        no more incidences.

        After initialization with :meth:`next_incidence`, this method will give
        all incidences of faces of `dimension_one` and `dimension_two`.
        `one[0]` will represent the index of a face in `dimension_one` and
        `two[0]` will represent the index of a face in `dimension_two`
        according to their order in `ListOfAllFaces`.

        Use :meth:`vertex_repr` and :meth:`facet_repr` to interprete the
        output.
        """
        cdef ListOfFaces dimension_one_faces
        cdef uint64_t *dimension_one_face
        cdef uint64_t **facets
        cdef uint64_t *face_one
        cdef size_t location
        cdef int is_it_equal
        if unlikely(not self.incidence_is_initialized):
            return 0
        if unlikely(self.incidence_is_initialized == 2):
            return self.next_trivial_incidence(one, two)
        if unlikely(self.incidence_is_initialized == 3):
            return self.next_trivial_incidence2(one, two)
        if unlikely(self.incidence_counter_one
                    == self.f_vector[self.incidence_dim_one + 1]):
            # in this case there are no more incidences
            self.incidence_is_initialized = 0
            return 0

        one[0] = self.incidence_counter_one
        dimension_one_face = self.data_vertex[self.incidence_dim_one + 1]\
                                             [self.incidence_counter_one]
        facets = self.data_vertex[self.dimension]
        while (self.incidence_counter_two < self.nr_facets):
            intersection(dimension_one_face,
                         facets[self.incidence_counter_two],
                         self.incidence_face,
                         self.face_length_vertex)
            location = \
                self.find_face(self.incidence_dim_two, self.incidence_face)
            is_it_equal = self.is_equal(self.incidence_dim_two,
                                        location, self.incidence_face)
            self.incidence_counter_two += 1
            if is_it_equal:
                two[0] = location
                if self.incidence_counter_two == self.nr_facets:
                    self.incidence_counter_one += 1
                    self.incidence_counter_two = 0
                return 1
        self.incidence_counter_one += 1
        self.incidence_counter_two = 0
        return self.next_incidence(one, two)

cdef class CombinatorialPolyhedron(SageObject):
    r"""
    The class of the Combinatorial Type of a Polyehdron, a Polytope.

    INPUT:

    - ``data`` -- a ``Polyhedron``, i.e. an instance of
      :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`.

    or

    - ``data`` -- a ``LatticePolytope``, i.e. an instance of
      :class:`~sage.geometry.lattice_polytope.LatticePolytopeClass`.

    or

    - ``data`` -- an ``incidence_matrix`` as in
      :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`
      of :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`.

      * ``vertices`` -- a list of ``[vertices, rays, lines]``, if
        the rows in the incidence_matrix should correspond to names.

      * ``facets`` -- a list of facets, if
        the columns in the incidence_matrix should correspond to names.

      * ``nr_lines`` -- for bounded Polyhedra, this should be the
      default: ``None``. For unbounded Polyhedra, this needs to be set
      to the correct nr of lines, i.e. the the maximum nr of lines with
      linearly independent directions in the Polyehdron.

    or

    - ``data`` -- a list of facets,
      each facet given as a list of ``[vertices, rays, lines]``.
      If the Polyhedron is unbounded, then rays and lines are required.
      If the Polyehdron contains no lines the rays can be thought of as
      the vertices of the facets deleted from a bounded Polyhedron. See
      :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`
      on how to use rays and lines.

      * ``facets`` -- a list of names of the facets, if
        the facets given should correspond to names.

      * ``unbounded`` -- for bounded Polyhedra, this should be the
      default ``None``. For unbounded Polyhedra, this needs to be set
      to the correct nr of lines, i.e. the the maximum nr of lines with
      linearly independent directions in the Polyehdron.


    EXAMPLES:

    Input is Polyhedron::

        sage: P = polytopes.cube()
        sage: CombinatorialPolyhedron(P)
        Combinatorial Type of a Polyhedron of dimension 3 with 8 vertices

    Input is a LatticePolytope::

        sage: points = [(1,0,0), (0,1,0), (0,0,1),
        ....: (-1,0,0), (0,-1,0), (0,0,-1)]
        sage: L = LatticePolytope(points)
        sage: CombinatorialPolyhedron(L)
        Combinatorial Type of a Polyhedron of dimension 3 with 6 vertices

    Input is an incidence matrix::

        sage: data = Polyhedron(rays=[[0,1]]).incidence_matrix()
        sage: CombinatorialPolyhedron(data, nr_lines=0)
        Combinatorial Type of a half-space of dimension 1
        sage: C = CombinatorialPolyhedron(data, vertices=['myvertex'],
        ....: facets=['myfacet'], nr_lines=0)
        sage: C.Vrepresentation()
        ('myvertex',)
        sage: C.Hrepresentation()
        ('myfacet',)

    You can also give the facets explicitely::

        sage: CombinatorialPolyhedron(((1,2,3),(1,2,4),(1,3,4),(2,3,4)))
        Combinatorial Type of a Polyhedron of dimension 3 with 4 vertices
        sage: facetnames = ['facet0', 'facet1', 'facet2', 'myfacet3']
        sage: facetinc = ((1,2,3),(1,2,4),(1,3,4),(2,3,4))
        sage: C = CombinatorialPolyhedron(facetinc, facets=facetnames)
        sage: C.Vrepresentation()
        (1, 2, 3, 4)
        sage: C.Hrepresentation()
        ('facet0', 'facet1', 'facet2', 'myfacet3')

    Specifying the nr of lines is important::

        sage: P = Polyhedron(ieqs=[[1,-1,0],[1,1,0]])
        sage: C = CombinatorialPolyhedron(P) #this works fine
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices
        sage: data = P.incidence_matrix()
        sage: C = CombinatorialPolyhedron(data) #wrong usage!
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 3 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C = CombinatorialPolyhedron(data, nr_lines=1) #correct
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices
        sage: C.f_vector()
        (1, 0, 2, 1)
    """
    # Those are the definitions done in the header:
    # cdef tuple _V
    # cdef tuple _H
    # cdef tuple _equalities
    # cdef dict _Hinv
    # cdef dict _Vinv
    # cdef int is_trivial  # in some instances the polyhedron might not
    # # have facets or otherwise produce errors in the C function
    # cdef int _dimension  # manually set, if `is_trivial`
    # cdef unsigned int _length_Hrep
    # cdef unsigned int _length_Vrep
    # cdef int _unbounded
    # cdef int _nr_lines
    # cdef ListOfFaces bitrep_facets
    # cdef ListOfFaces bitrep_vertices
    # cdef int polar
    # cdef size_t *_f_vector
    # cdef size_t _length_edges_list
    # cdef size_t **_edges
    # cdef size_t **_ridges
    # cdef size_t **_face_lattice_incidences
    # cdef size_t _nr_edges
    # cdef size_t _nr_ridges
    # cdef size_t _nr_face_lattice_incidences
    # cdef ListOfAllFaces _all_faces
    # cdef size_t _nr_facets
    # cdef tuple _MemoryAllocators

    def __init__(self, data, vertices=None, facets=None, nr_lines=None):
        r"""
        Initializes the combinatorial polyhedron.

        See :class:`CombinatorialPolyhedron`.

        TESTS::

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],
            ....: [0,2,3],[1,2,3]])    # indirect doctests
        """
        cdef unsigned int **incidence_matrix
        cdef unsigned int **facets_pointer
        cdef unsigned int *len_facets
        self._dimension = -2
        self._f_vector = NULL
        self._edges = NULL
        self._ridges = NULL
        self._face_lattice_incidences = NULL
        self._polar = 0
        self._nr_lines = 0
        self._length_edges_list = 16348
        self._all_faces = None
        self._MemoryAllocators = ()
        if nr_lines is None:
            self._unbounded = 0
            self._nr_lines = 0
        else:
            self._unbounded = 1 + int(nr_lines)
            self._nr_lines = int(nr_lines)
        self.is_trivial = 0
        self._equalities = ()

        if is_Polyhedron(data):
            # Input is `Polyhedron`
            if data.is_empty():
                self.is_trivial = 1
                self._dimension = -1
                return
            vertices = data.Vrepresentation()
            facets = tuple(inequality for inequality in data.Hrepresentation())

            if (len(vertices) == data.n_lines() + 1) and (data.n_lines > 0):
                # in this case the Polyhedron does not have facets
                self.is_trivial = 1
                self._dimension = data.n_lines()
                self._V = tuple(vertices)
                return

            if not data.is_compact():
                self._unbounded = 1 + data.n_lines()
                self._nr_lines = int(data.n_lines())

            data = data.incidence_matrix()
            self._length_Hrep = data.ncols()
            self._length_Vrep = data.nrows()

        if is_LatticePolytope(data):
            # Input is LatticePolytope
            if data.npoints() == 0:
                self.is_trivial = 1
                self._dimension = -1
                return
            if data.npoints() == 1:
                self.is_trivial = 1
                self._dimension = 0
                self._V = data.vertices()
                return

            self._unbounded = 0
            self._nr_lines = 0
            vertices = data.vertices()
            self._length_Vrep = len(vertices)
            facets = data.facets()
            self._length_Hrep = len(facets)
            data = tuple(tuple(vert for vert in facet.vertices())
                         for facet in facets)

        if vertices:
            self._V = tuple(vertices)
            self._Vinv = {v: i for i,v in enumerate(self._V)}
        else:
            self._V = None
            self._Vinv = None

        if facets:
            facets = tuple(facets)
            test = [0] * len(facets)
            for i in range(len(facets)):
                test[i] = 1
                if hasattr(facets[i], "is_inequality"):
                    if not facets[i].is_inequality():
                        test[i] = 0
            self._H = \
                tuple(facets[i] for i in range(len(facets)) if test[i])
            # only keeping those that are actual inequalities
            self._equalities = \
                tuple(facets[i] for i in range(len(facets)) if not test[i])
            # the inequalities are saved here
            self._Hinv = {v: i for i,v in enumerate(self._H)}
        else:
            self._H = None
            self._Hinv = None

        if is_Matrix(data):
            self._length_Hrep = data.ncols()
            self._length_Vrep = data.nrows()
            rg = range(data.nrows())
            matrix = tuple(tuple(data[i,j] for i in rg)
                           for j in range(data.ncols())
                           if not all(data[i,j] for i in rg))
            # transpose and get rid of equalities (which all vertices satisfie)
            self._nr_facets = len(matrix)

            if len(matrix) == 0:  # the case of the empty Polyhedron
                self.is_trivial = 1
                self._dimension = -1 + data.nrows()
                # the elements in Vrep are assumed to be one vertex
                # and otherwise lines
                return

            if len(matrix) == 1:  # the case of a half space
                self.is_trivial = 2
                self._dimension = -1 + data.nrows()
                # the elements in Vrep are assumed to be one vertex
                # and otherwise lines
                return

            nr_vertices = self._length_Vrep
            nr_facets = len(matrix)
            if self._unbounded or nr_facets <= nr_vertices:
                self._polar = 0
                # initializing the facets as BitVectors
                self.bitrep_facets = \
                    get_facets_from_incidence_matrix(matrix)
                # initializing the vertices as BitVectors
                self.bitrep_vertices = \
                    get_vertices_from_incidence_matrix(matrix)
            else:
                # for speed we will do all calculations with the polar/dual
                self._polar = 1
                # initializing the vertices as BitVectors
                self.bitrep_vertices = \
                    get_facets_from_incidence_matrix(matrix)
                # initializing the facets as BitVectors
                self.bitrep_facets = \
                    get_vertices_from_incidence_matrix(matrix)

        elif isinstance(data, Integer):  # intput for a trivial Polyhedron
            if data < -1:
                TypeError("any polyhedron must have dimension at least -1")
            self._nr_facets = 0
            self.is_trivial = 1
            self._dimension = data

        else:
            if is_iterator(data):
                data = tuple(data)
            # assumes the facets are given as a list of vertices/rays/lines
            if len(data) == 0:
                # the empty Polyhedron
                self.is_trivial = 1
                self._dimension = -1
                return

            if len(data) == 1:
                self.is_trivial = 2
                self._dimension = len(data[0]) - 1
                if self._dimension <= 0:
                    self.is_trivial = 1
                    # we are treating a polyhedron equal to its affine hull
                return

            if self._V is None:
                vertices = sorted(set.union(*map(set, data)))
                nr_vertices = len(vertices)
                if vertices != range(len(vertices)):
                    self._V = tuple(vertices)
                    self._Vinv = {v: i for i,v in enumerate(self._V)}
            else:
                nr_vertices = len(self._V)
            self._length_Vrep = nr_vertices

            if self._V is not None:
                def f(v): return self._Vinv[v]
            else:
                def f(v): return int(v)

            facets = tuple(tuple(f(i) for i in j) for j in data)
            self._nr_facets = len(facets)
            self._length_Hrep = len(facets)

            if self._unbounded or len(facets) <= nr_vertices:
                self._polar = 0
                # initializing the facets as BitVectors
                self.bitrep_facets = \
                    get_facets_bitrep_from_facets_tuple(
                        facets, nr_vertices)
                # initializing the vertices as BitVectors

                self.bitrep_vertices = \
                    get_vertices_bitrep_from_facets_tuple(
                        facets, nr_vertices)
            else:
                # for speed we will do all calculations with the polar/dual
                self._polar = 1
                # initializing the vertices as BitVectors
                self.bitrep_vertices = \
                    get_facets_bitrep_from_facets_tuple(
                        facets, nr_vertices)
                # initializing the facets as BitVectors
                self.bitrep_facets = \
                    get_vertices_bitrep_from_facets_tuple(
                        facets, nr_vertices)

    def _repr_(self):
        r"""
        Returns a description of the Combinatorial Polyhedron.

        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a Polyhedron of dimension 3 with 4 vertices'

            sage: P = Polyhedron(vertices=[])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of the empty Polyhedron'

            sage: P = Polyhedron(vertices=[[0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of the Polyhedron with one vertex'

            sage: P = Polyhedron(lines=[[0,0,1],[0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a trivial Polyhedron of dimension 2'

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[-1,0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a half-space of dimension 2'
        """
        if self.is_trivial == 1:
            if self._dimension == 0:
                return "Combinatorial Type of the Polyhedron with one vertex"

            if self._dimension > 0:
                return "Combinatorial Type of a trivial Polyhedron of \
                    dimension %s" % self._dimension

            return "Combinatorial Type of the empty Polyhedron"

        if self.is_trivial == 2:
            return "Combinatorial Type of a half-space of dimension %s"\
                % self._dimension

        return "Combinatorial Type of a Polyhedron of "\
                "dimension %s with %s vertices" \
                % (self.dimension(), len(self.vertices()))

    def __reduce__(self):
        r"""
        Override __reduce__ to correctly pickle/unpickle.

        TESTS::

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: tuple(i for i in C.face_iter(facet_repr=True)) \
            ....: == tuple(i for i in C1.face_iter(facet_repr=True))
            True

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: tuple(i for i in C.face_iter(facet_repr=True)) \
            ....: == tuple(i for i in C1.face_iter(facet_repr=True))
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: tuple(i for i in C.face_iter(facet_repr=True)) \
            ....: == tuple(i for i in C1.face_iter(facet_repr=True))
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0],
            ....:                      [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: C.f_vector() == C1.f_vector()
            True
        """
        nr_lines = None
        if self._unbounded:
            nr_lines = Integer(self._nr_lines)

        if self.is_trivial == 1:
            return (CombinatorialPolyhedron, (Integer(self._dimension),
                    self.Vrepresentation(), self.Hrepresentation(), nr_lines))

        if self.is_trivial == 2:
            pickletuple = tuple(0 for _ in range(self._dimension + 1))
            return (CombinatorialPolyhedron, ((pickletuple,),
                    self.Vrepresentation(), self.Hrepresentation(), nr_lines))
            # if `data` is a `tuple` containining exactly one tuple
            # of some length, than __init__ will figure it is a halfspace

        return (CombinatorialPolyhedron, (self.facets(),
                self.Vrepresentation(), self.Hrepresentation(), nr_lines))

    def Vrepresentation(self):
        r"""
        Return a list of names of ``[vertices, rays, lines]``.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0,0], [0,1,0], \
            ....:                      [0,0,1],[0,0,-1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.Vrepresentation()
            (A line in the direction (0, 0, 1),
             A ray in the direction (1, 0, 0),
             A vertex at (0, 0, 0),
             A ray in the direction (0, 1, 0))
        """
        if self._V is not None:
            return self._V
        else:
            return tuple(Integer(i) for i in range(self._length_Vrep))

    def Hrepresentation(self):
        r"""
        Returns a list of names of facets and possibly some equalities.

        EXAMPLES::

            sage: P = polytopes.permutahedron(3)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.Hrepresentation()
            (An equation (1, 1, 1) x - 6 == 0,
             An inequality (0, -1, -1) x + 5 >= 0,
             An inequality (0, 0, -1) x + 3 >= 0,
             An inequality (0, -1, 0) x + 3 >= 0,
             An inequality (0, 1, 0) x - 1 >= 0,
             An inequality (0, 1, 1) x - 3 >= 0,
             An inequality (0, 0, 1) x - 1 >= 0)
        """
        if self._H is not None:
            return self._equalities + self._H
        else:
            return tuple(Integer(i) for i in range(self._length_Hrep))

    def vertices(self, names=True):
        r"""
        Returns the elements in the ``Vrepresentation`` that are vertices.

        In the case of an unbounded Polyhedron, there might be lines and
        rays in the Vrepresentation.

        If ``names`` is set to ``False``, then the vertices are given by
        their indices in the Vrepresentation.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[0,0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.vertices()
            (A vertex at (0, 0, 0),)
            sage: C.Vrepresentation()
            (A vertex at (0, 0, 0),
             A ray in the direction (0, 0, 1),
             A ray in the direction (0, 1, 0),
             A ray in the direction (1, 0, 0))
            sage: P = polytopes.cross_polytope(3)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.vertices()
            (A vertex at (-1, 0, 0),
             A vertex at (0, -1, 0),
             A vertex at (0, 0, -1),
             A vertex at (0, 0, 1),
             A vertex at (0, 1, 0),
             A vertex at (1, 0, 0))

            sage: C.vertices(names=False)
            (0, 1, 2, 3, 4, 5)

            sage: points = [(1,0,0), (0,1,0), (0,0,1),
            ....:           (-1,0,0), (0,-1,0), (0,0,-1)]
            sage: L = LatticePolytope(points)
            sage: C = CombinatorialPolyhedron(L)
            sage: C.vertices()
            (M(1, 0, 0), M(0, 1, 0), M(0, 0, 1), M(-1, 0, 0), M(0, -1, 0), M(0, 0, -1))
            sage: C.vertices(names=False)
            (0, 1, 2, 3, 4, 5)
        """
        return tuple(i[0] for i in self.faces(0, names=names))

    def facets(self, names=True):
        r"""
        Returns the facets as lists of ``[vertices, rays, lines]``.

        If ``names`` is ``False``, then the vertices in the facets
        are given by their indices in the Vrepresentation.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.facets()
            ((A vertex at (-1, -1, 1),
              A vertex at (-1, 1, 1),
              A vertex at (1, -1, 1),
              A vertex at (1, 1, 1)),
             (A vertex at (-1, 1, -1),
              A vertex at (-1, 1, 1),
              A vertex at (1, 1, -1),
              A vertex at (1, 1, 1)),
             (A vertex at (1, -1, -1),
              A vertex at (1, -1, 1),
              A vertex at (1, 1, -1),
              A vertex at (1, 1, 1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, -1, 1),
              A vertex at (-1, 1, -1),
              A vertex at (-1, 1, 1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, 1, -1),
              A vertex at (1, -1, -1),
              A vertex at (1, 1, -1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, -1, 1),
              A vertex at (1, -1, -1),
              A vertex at (1, -1, 1)))
            sage: C.facets(names=False)
            ((1, 3, 5, 7),
             (2, 3, 6, 7),
             (4, 5, 6, 7),
             (0, 1, 2, 3),
             (0, 2, 4, 6),
             (0, 1, 4, 5))
        """
        tup = self.faces(self.dimension()-1, names=names)
        # it is important to have the facets in the exact same order as
        # on input
        # every facet knows its order by the facet representation
        order = self.faces(self.dimension()-1, facet_repr=True, names=False)
        dic = {}
        for i in range(len(tup)):
            dic[order[i][0]] = tup[i]
        return tuple(dic[i] for i in range(len(tup)))

    def edges(self, names=True):
        r"""
        Returns the edges of the CombinatorialPolyhedron,
        i.e. the rank 1 faces, which contain 2 vertices.

        If ``names`` is set to ``False``, then the vertices in the edges
        are given by their indices in the Vrepresentation.

        If you want to compute all faces of dimension 1,
        use :meth:`CombinatorialPolyhedron.faces` instead.

        .. NOTE::

            To compute edges and f_vector, first compute the edges.
            This might be faster.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(3,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ((A vertex at (3, 9, 27), A vertex at (4, 16, 64)),
             (A vertex at (2, 4, 8), A vertex at (4, 16, 64)),
             (A vertex at (1, 1, 1), A vertex at (4, 16, 64)),
             (A vertex at (0, 0, 0), A vertex at (4, 16, 64)),
             (A vertex at (2, 4, 8), A vertex at (3, 9, 27)),
             (A vertex at (0, 0, 0), A vertex at (3, 9, 27)),
             (A vertex at (1, 1, 1), A vertex at (2, 4, 8)),
             (A vertex at (0, 0, 0), A vertex at (2, 4, 8)),
             (A vertex at (0, 0, 0), A vertex at (1, 1, 1)))

            sage: C.edges(names=False)
            ((3, 4), (2, 4), (1, 4), (0, 4), (2, 3), (0, 3), (1, 2), (0, 2), (0, 1))

            sage: P = Polyhedron(rays=[[-1,0],[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ()
            sage: C.faces(1)
            ((A line in the direction (1, 0), A vertex at (0, 0)),)

            sage: P = Polyhedron(vertices=[[0,0],[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edges()
            ((A vertex at (0, 0), A vertex at (1, 0)),)

            sage: from itertools import combinations
            sage: N = combinations(['a','b','c','d','e'], 4)
            sage: C = CombinatorialPolyhedron(N)
            sage: C.edges()
            (('d', 'e'),
             ('c', 'e'),
             ('c', 'd'),
             ('b', 'e'),
             ('b', 'd'),
             ('b', 'c'),
             ('a', 'e'),
             ('a', 'd'),
             ('a', 'c'),
             ('a', 'b'))

        TESTS::

            sage: from itertools import combinations
            sage: N = combinations(range(20),19)
            sage: C = CombinatorialPolyhedron(N)
            sage: try:
            ....:     alarm(0.001)
            ....:     C.edges()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: len(C.edges())
            190
        """
        cdef size_t **edges
        cdef size_t nr_edges
        cdef size_t len_edgelist = self._length_edges_list
        cdef size_t j
        if self._polar:
            if not self._ridges:
                self._calculate_ridges()
            if self._ridges is NULL:
                raise KeyboardInterrupt()
            edges = self._ridges
            nr_edges = self._nr_ridges
        else:
            if not self._edges:
                self._calculate_edges()
            if self._edges is NULL:
                raise KeyboardInterrupt()
            edges = self._edges
            nr_edges = self._nr_edges

        # the edges are being saved in a list basically
        # with the first entry the first vertex of the first edges,
        # the second entry the second vertex of that edge

        # possibly there are many edges,
        # hence they are stored in an array of arrays,
        # with each array containing maxnumberedges of edges

        if self._V is not None and names is True:
            def f(size_t i): return self._V[i]
        else:
            def f(size_t i): return Integer(i)

        def vertex_one(size_t i):
            return f(edges[i // len_edgelist][2*(i % len_edgelist)])

        def vertex_two(size_t i):
            return f(edges[i // len_edgelist][2*(i % len_edgelist)+1])

        return tuple((vertex_one(j), vertex_two(j)) for j in range(nr_edges))

    def edge_graph(self, names=True):
        r"""
        Returns the edge graph.

        If ``names`` is set to ``False``, the vertices will carry names
        according to the indexing of the Vrepresentation.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(3,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.edge_graph()
            Graph on 5 vertices
            sage: G = C.edge_graph()
            sage: G.degree()
            [4, 3, 4, 3, 4]
        """
        return Graph(self.edges(names=names), format="list_of_edges")

    def dimension(self):
        r"""
        Returns the dimension of the ``CombinatorialPolyehdron``.

        EXAMPLES::

            sage: C = CombinatorialPolyhedron([(1,2,3), (1,2,4),
            ....:                              (1,3,4), (2,3,4)])
            sage: C.dimension()
            3

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
            sage: CombinatorialPolyhedron(P).dimension()
            3

        TESTS::

            sage: from itertools import combinations
            sage: N = combinations(range(1200), 1199)
            sage: C = CombinatorialPolyhedron(N)
            sage: try:
            ....:     alarm(0.1)
            ....:     C.dimension()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: try:
            ....:     alarm(0.1)
            ....:     C.f_vector()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
        """
        if self._dimension == -2:
            self._dimension = calculate_dimension(self.bitrep_facets)
        return Integer(self._dimension)

    def ridges(self, add_equalities=False, names=True):
        r"""
        Returns the ridges.

        The ridges of the CombinatorialPolyhedron are the faces
        contained in exactly two facets.

        If you want to compute all faces of codimension 1,
        use :meth:`CombinatorialPolyhedron.faces` instead.

        The ridges will be given by the facets, they are contained in.

        - If ``add_equalities`` is ``True``, then equalities the entire
          Polyhedron satisfies, are added.

        - If ``names`` is ``False``, then the facets in the ridges are
          given by their indices in the Hrepresentation.

        .. NOTE::

            To compute ridges and f_vector, compute the ridges first.
            This might be faster.

        EXAMPLES::

            sage: P = polytopes.permutahedron(2)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridges()
            ((An inequality (0, -1) x + 2 >= 0, An inequality (0, 1) x - 1 >= 0),)
            sage: C.ridges(add_equalities=True)
            (((An equation (1, 1) x - 3 == 0, An inequality (0, -1) x + 2 >= 0),
              (An equation (1, 1) x - 3 == 0, An inequality (0, 1) x - 1 >= 0)),)

            sage: P = polytopes.cyclic_polytope(4,5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridges()
            ((An inequality (24, -26, 9, -1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-12, 19, -8, 1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (-50, 35, -10, 1) x + 24 >= 0),
             (An inequality (-12, 19, -8, 1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (24, -26, 9, -1) x + 0 >= 0),
             (An inequality (8, -14, 7, -1) x + 0 >= 0,
              An inequality (-12, 19, -8, 1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (-12, 19, -8, 1) x + 0 >= 0),
             (An inequality (-6, 11, -6, 1) x + 0 >= 0,
              An inequality (8, -14, 7, -1) x + 0 >= 0))
            sage: C.ridges(names=False)
            ((3, 4),
             (2, 4),
             (1, 4),
             (0, 4),
             (2, 3),
             (1, 3),
             (0, 3),
             (1, 2),
             (0, 2),
             (0, 1))

            sage: P = Polyhedron(rays=[[1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C
            Combinatorial Type of a half-space of dimension 1
            sage: C.ridges()
            ()
            sage: C.faces(0, facet_repr=True)
            ((An equation (0, 1) x + 0 == 0, An inequality (1, 0) x + 0 >= 0),)

        TESTS::

            sage: from itertools import combinations
            sage: N = combinations(range(200),199)
            sage: C = CombinatorialPolyhedron(N)
            sage: try:
            ....:     alarm(0.01)
            ....:     C.ridges()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: len(C.ridges())
            19900
        """
        cdef size_t **ridges
        cdef size_t nr_ridges
        cdef size_t len_edgelist = self._length_edges_list
        cdef size_t j
        if self._polar:
            if not self._edges:
                self._calculate_edges()
            if self._edges is NULL:
                raise KeyboardInterrupt()
            ridges = self._edges
            nr_ridges = self._nr_edges
        else:
            if not self._ridges:
                self._calculate_ridges()
            if self._ridges is NULL:
                raise KeyboardInterrupt()
            ridges = self._ridges
            nr_ridges = self._nr_ridges

        # the ridges are being saved in a list basically
        # with the first entry the first facet of the first ridge,
        # the second entry the second facet of that ridge

        # possibly there are many ridges,
        # hence they are stored in an array of arrays,
        # with each array containing maxnumberedges of ridges

        if self._H is not None and names is True:
            def f(size_t i): return self._H[i]
        else:
            def f(size_t i): return Integer(i)

        def facet_one(size_t i):
            return f(ridges[i // len_edgelist][2*(i % len_edgelist)])

        def facet_two(size_t i):
            return f(ridges[i // len_edgelist][2*(i % len_edgelist)+1])

        if add_equalities:
            return tuple(
                ((self._equalities + (facet_one(i),)),
                 (self._equalities + (facet_two(i),)))
                for i in range(nr_ridges))
        else:
            return tuple((facet_one(i), facet_two(i))
                         for i in range(nr_ridges))

    def ridge_graph(self, names=True):
        r"""
        Returns the ridge graph.

        The ridge graph of the CombinatorialPolyhedron consists of
        ridges as edges and facets as vertices.

        If ``names`` is ``False``, the ``vertices`` of the graph  will
        be the incidences of the facets in the Hrepresentation.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(4,6)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.ridge_graph()
            Graph on 9 vertices
        """
        return Graph(self.ridges(names=names), format="list_of_edges")

    def f_vector(self):
        r"""
        Calculates the ``f_vector`` of the CombinatorialPolyhedron.

        The ``f_vector`` contains the number of faces of dimension ``k``
        for each ``k`` in ``range(-1, self.dimension() + 1)``.

        .. NOTE::

            If you also want to compute edges and/or ridges, do so
            first.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.f_vector()
            (1, 120, 240, 150, 30, 1)

            sage: P = polytopes.cyclic_polytope(6,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.f_vector()
            (1, 10, 45, 120, 185, 150, 50, 1)

        TESTS::

            sage: from itertools import combinations
            sage: N = combinations(range(20),19)
            sage: C = CombinatorialPolyhedron(N)
            sage: try:
            ....:     alarm(0.001)
            ....:     C.f_vector()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: C.f_vector()
            (1,
             20,
             190,
             1140,
             4845,
             15504,
             38760,
             77520,
             125970,
             167960,
             184756,
             167960,
             125970,
             77520,
             38760,
             15504,
             4845,
             1140,
             190,
             20,
             1)

        ALGORITHM:

        The number of facets is assumed to be at least two here.
        The algorithm to visit all proper faces exactly once is roughly
        equivalent to:

        facets = [set(facet) for facet in self.facets()]
        ComputeNextStep(facets, [])

        # this algorithm assumes at each step to receive all facets of
        # some face except those contained in a face of forbidden
        def ComputeNextStep(faces, forbidden):

            for face in faces:
                pass #do something here with that face

            while len(faces) > 1:
                one_face = faces.pop()
                newfaces = [one_face.intersection(face) for face in faces]
                #newfaces contains all intersection

                newfaces2 = []
                for face1 in newfaces:
                    # face1 is a facet of one_face iff
                    # it is not contained in another facet
                    if all(not face1 < face2 for face2 in newfaces):
                        newfaces2.append(face1)
                # newfaces2 contains all facets of one_face not contained
                # in any one of forbidden and maybe some that are
                # contained in one of forbidden

                newfaces3 = []
                for face1 in newfaces2:
                    if all(not face1 < face2 for face2 in forbidden):
                        newfaces3.append(face1)
                # newfaces3 contains exactly all facets of one_face but
                # those contained in one face of forbidden

                # visit all faces in one_face that are not contained in
                # one of forbidden
                ComputeNextStep(newfaces3, forbidden)

                # we have visited all faces in one_face, so we should not
                # visit one ever again

                forbidden.append(one_face)

            return
        """
        cdef int dim = self.dimension()
        if self._f_vector is NULL:
            self._calculate_f_vector()
        if self._f_vector is NULL:
            raise ValueError("could not determine f_vector")
        if self._polar:
            f = tuple(Integer(self._f_vector[dim+1-i]) for i in range(dim + 2))
        else:
            f = tuple(Integer(self._f_vector[i]) for i in range(dim + 2))
        return f

    def faces(self, dimension, facet_repr=False, names=True):
        r"""
        Gets all k-faces for specified dimenion k.

        By default faces are given as tuple of vertices.

        If ``facet_repr`` is set to ``True``, then vertices are given in
        as tuple of facets.

        If ``names`` is set to ``False``, then the vertices are given by
        their indices in the Vrepresentation.


        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.faces(2)
             ((A vertex at (0, 0, 1, 0),
              A vertex at (0, 1, 0, 0),
              A vertex at (1, 0, 0, 0)),
             (A vertex at (0, 0, 0, 1),
              A vertex at (0, 1, 0, 0),
              A vertex at (1, 0, 0, 0)),
             (A vertex at (0, 0, 0, 1),
              A vertex at (0, 0, 1, 0),
              A vertex at (1, 0, 0, 0)),
             (A vertex at (0, 0, 0, 1),
              A vertex at (0, 0, 1, 0),
              A vertex at (0, 1, 0, 0)))
            sage: C.faces(2, facet_repr=True)
            ((An equation (1, 1, 1, 1) x - 1 == 0, An inequality (0, 0, 0, 1) x + 0 >= 0),
             (An equation (1, 1, 1, 1) x - 1 == 0, An inequality (0, 0, 1, 0) x + 0 >= 0),
             (An equation (1, 1, 1, 1) x - 1 == 0, An inequality (0, 1, 0, 0) x + 0 >= 0),
             (An equation (1, 1, 1, 1) x - 1 == 0,
              An inequality (0, -1, -1, -1) x + 1 >= 0))

            sage: P = Polyhedron(rays=[[1,0],[0,1]], vertices=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.faces(1)
            ((A vertex at (0, 1), A vertex at (1, 0)),
             (A ray in the direction (1, 0), A vertex at (1, 0)),
             (A ray in the direction (0, 1), A vertex at (0, 1)))
            sage: C.faces(1, names=False)
            ((1, 3), (2, 3), (0, 1))
            sage: C.faces(1, facet_repr=True)
            ((An inequality (1, 1) x - 1 >= 0,),
             (An inequality (0, 1) x + 0 >= 0,),
             (An inequality (1, 0) x + 0 >= 0,))

        ALGORITHM:

        See :meth:`f_vector` for a description
        on how all faces are visited.
        """
        if self.is_trivial == 1:
            # the case of a Polyhedron equal to its affine hull
            if dimension == -1:
                return ((),)
            if dimension == self._dimension:
                if facet_repr is True:
                    return ((),)
                if self._V is not None and names is True:
                    return (tuple(self._V),)
                else:
                    if dimension == 0:
                        return (Integer(0),)
                    return ((),)
            return ()

        if self.is_trivial == 2:  # the case of a half-space
            if dimension == -1:
                if facet_repr is True:
                    if self._H is not None and names is True:
                        return (self._equalities + tuple(self._H),)
                    else:
                        return ((Integer(0),),)
                return ((),)
            if dimension == self._dimension - 1:
                if facet_repr is True:
                    if self._H is not None and names is True:
                        return (self._equalities + tuple(self._H),)
                    else:
                        return ((Integer(0),),)
                return ((),)
            if dimension == self._dimension:
                return ((),)
            return ()

        dim = self.dimension()
        if dimension not in range(-1, dim + 1):
            return ()

        return tuple(i for i in self.face_iter(dimension=dimension,
                                               vertex_repr=(not facet_repr),
                                               facet_repr=facet_repr,
                                               names=names))

    def face_iter(self, dimension=None, vertex_repr=True,
                  facet_repr=False, give_dimension=False, names=True):
        r"""
        Iterator over all faces of specified dimension.

        If ``dimension`` is not specified then iterate over all faces.

        If ``vertex_repr`` is ``True``, then give the faces as lists of
        elements in ``Vrepresentation``.

        If ``facet_repr`` is ``True``, then give the faces as lists of
        elements in ``Hrepresentation``.

        - Both ``vertex_repr`` and ``facet_repr`` can be set to ``True,
          this will give each face as tuple of the form
          (``vertex_repr``, ``facet_rerpr``).

        If ``names`` is ``False``, then vertices and facets are labeled
        by their indexes.

        If ``give_dimension`` is ``True``, then the dimension of the face is
        printed as well.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(dimension=2)
            sage: next(it)
            (A vertex at (4, 1, 5, 2, 3),
             A vertex at (4, 2, 5, 1, 3),
             A vertex at (5, 1, 4, 2, 3),
             A vertex at (5, 2, 4, 1, 3))
            sage: next(it)
            (A vertex at (4, 1, 5, 2, 3),
             A vertex at (4, 1, 5, 3, 2),
             A vertex at (5, 1, 4, 2, 3),
             A vertex at (5, 1, 4, 3, 2))
            sage: it = C.face_iter(dimension=2, names=False)
            sage: next(it)
            (76, 82, 100, 106)
            sage: next(it)
            (76, 77, 100, 101)
            sage: it = C.face_iter(dimension=2, facet_repr=True)
            sage: next(it)
            ((A vertex at (4, 1, 5, 2, 3),
              A vertex at (4, 2, 5, 1, 3),
              A vertex at (5, 1, 4, 2, 3),
              A vertex at (5, 2, 4, 1, 3)),
             (An equation (1, 1, 1, 1, 1) x - 15 == 0,
              An inequality (0, 1, 0, 1, 0) x - 3 >= 0,
              An inequality (0, 1, 0, 1, 1) x - 6 >= 0))
            sage: it = C.face_iter(dimension=2, vertex_repr=False,
            ....:                  facet_repr=True, names=False)
            sage: next(it)
            (28, 29)
            sage: next(it)
            (25, 29)

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],[0,2,3],[1,2,3]])
            sage: it = C.face_iter(give_dimension=True)
            sage: for i in it: i
            ((), -1)
            (((0,), (1,), (2,), (3,)), 3)
            ((1, 2, 3), 2)
            ((0, 2, 3), 2)
            ((0, 1, 3), 2)
            ((0, 1, 2), 2)
            ((2, 3), 1)
            ((1, 3), 1)
            ((1, 2), 1)
            ((3,), 0)
            ((2,), 0)
            ((1,), 0)
            ((0, 3), 1)
            ((0, 2), 1)
            ((0,), 0)
            ((0, 1), 1)

        TESTS::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: f = C.f_vector()
            sage: altf = tuple(len(tuple(C.face_iter(i)))
            ....:              for i in range(-1,C.dimension()+1))
            sage: altf == f
            True
            sage: allfaces = tuple(tuple(C.face_iter(i))
            ....:                  for i in range(-1,C.dimension()+1))
            sage: C._record_all_faces()
            sage: allfaces2 = tuple(tuple(C.face_iter(i))
            ....:                   for i in range(-1,C.dimension()+1))
            sage: all(sorted(sorted(j) for j in allfaces[i]) == \
            ....:   sorted(sorted(j) for j in allfaces2[i])
            ....:          for i in range(C.dimension()+2))
            True

            sage: P = polytopes.cyclic_polytope(5,20)
            sage: C = CombinatorialPolyhedron(P)
            sage: f = C.f_vector()
            sage: altf = tuple(len(tuple(C.face_iter(i)))
            ....:              for i in range(-1,C.dimension()+1))
            sage: altf == f
            True
            sage: allfaces = tuple(tuple(C.face_iter(i))
            ....:                  for i in range(-1,C.dimension()+1))
            sage: C._record_all_faces()
            sage: allfaces2 = tuple(tuple(C.face_iter(i))
            ....:                   for i in range(-1,C.dimension()+1))
            sage: all(sorted(sorted(j) for j in allfaces[i]) == \
            ....: sorted(sorted(j) for j in allfaces2[i])
            ....:        for i in range(C.dimension()+2))
            True

        ALGORITHM:

        See :meth:`f_vector` for a description
        on how all faces are visited.
        """
        cdef FaceIterator face_iter
        cdef ListOfAllFaces all_faces
        cdef ListOfFaces facets
        cdef ListOfFaces vertices
        cdef int dim
        cdef int next_dim
        cdef size_t *output1
        cdef size_t *output2
        cdef size_t *f_vector
        cdef size_t index

        if not vertex_repr and not facet_repr and not give_dimension:
            # the user wants no output for some reason
            return
        if dimension is not None:
            dimension = int(dimension)
            dimensionrange = (dimension,)
        else:
            dimensionrange = range(-1, self.dimension()+1)
            dimension = -2

        if self.is_trivial > 0:  # taking care of the trivial polynomial
            for dim in dimensionrange:
                if vertex_repr and facet_repr:
                    vert = self.faces(dim, names=names)
                    fac = self.faces(dim, facet_repr=True, names=names)
                    for i in range(len(vert)):
                        if give_dimension:
                            yield (vert[i], fac[i], Integer(dim))
                        else:
                            yield (vert[i], fac[i])
                elif vertex_repr:
                    vert = self.faces(dim, names=names)
                    for i in range(len(vert)):
                        if give_dimension:
                            yield (vert[i], Integer(dim))
                        else:
                            yield vert[i]
                elif facet_repr:
                    fac = self.faces(dim, facet_repr=True, names=names)
                    for i in range(len(fac)):
                        if give_dimension:
                            yield (fac[i], Integer(dim))
                        else:
                            yield fac[i]
                else:
                    fac = self.faces(dim, facet_repr=True, names=False)
                    for i in range(len(fac)):
                        yield Integer(dim)
            return

        facets = self.bitrep_facets
        addtuple = ()
        # translating the result to the desired representation
        if names:
            addtuple = self._equalities
        if self._H is not None and names is True:
            def h(i): return self._H[i]
        else:
            def h(i): return Integer(i)

        if self._V is not None and names is True:
            def v(i): return self._V[i]
        else:
            def v(i): return Integer(i)

        if -1 in dimensionrange:
            # yield the empty face
            vert = ()
            fac = tuple((h(i),) + addtuple
                        for i in range(self._nr_facets))
            if vertex_repr and facet_repr:
                if give_dimension:
                    yield (vert, fac, Integer(-1))
                else:
                    yield (vert, fac)
            elif vertex_repr:
                if give_dimension:
                    yield (vert, Integer(-1))
                else:
                    yield vert
            elif facet_repr:
                if give_dimension:
                    yield (fac, Integer(-1))
                else:
                    yield fac
            else:
                yield Integer(-1)
            if -1 == dimension:
                return

        if 0 == dimension and not facet_repr and not self._unbounded:
            # in most cases, there is an easy way to yield all vertices
            if vertex_repr and give_dimension:
                for j in range(self._length_Vrep):
                    yield (v(j),Integer(0))
            elif vertex_repr:
                for j in range(self._length_Vrep):
                    yield (v(j),)
            else:
                for j in range(self._length_Vrep):
                    yield Integer(0)
            return

        dim = self.dimension()
        if dim in dimensionrange:
            # yield the full Polyhedron
            fac = addtuple
            vert = tuple((v(i),) for i in range(self._length_Vrep))
            if vertex_repr and facet_repr:
                if give_dimension:
                    yield (vert, fac, Integer(dim))
                else:
                    yield (vert, fac)
            elif vertex_repr:
                if give_dimension:
                    yield (vert, Integer(dim))
                else:
                    yield vert
            elif facet_repr:
                if give_dimension:
                    yield (fac, Integer(dim))
                else:
                    yield fac
            else:
                yield Integer(dim)
            if dim == dimension:
                return

        specialcase = 0

        # figuring out how to determin vertex-repr, facet_repr and dimension
        # of a face
        if self._all_faces:
            # if there is a already a list of all faces,
            # we do not need to start
            # the iterator again
            all_faces = self._all_faces
            if self._polar:
                def get_vertex_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = all_faces.facet_repr(next_dim, index, output2)
                    return tuple(v(output2[t]) for t in range(leng))

                def get_facet_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = all_faces.vertex_repr(next_dim, index, output1)
                    return addtuple + tuple(h(output1[t]) for t in range(leng))

                def print_dim():
                    return Integer(dim - 1 - next_dim)
                if dimension != -2:
                    dimension = dim - 1 - dimension
            else:
                def get_facet_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = all_faces.facet_repr(next_dim, index, output2)
                    return addtuple + tuple(h(output2[t]) for t in range(leng))

                def get_vertex_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = all_faces.vertex_repr(next_dim, index, output1)
                    return tuple(v(output1[t]) for t in range(leng))

                def print_dim():
                    return Integer(next_dim)
        else:
            if dimension == 0 and not self._unbounded and \
                    not self._polar:
                # we have stored the vertices, so we should use them
                # (if the Polyhedron is unbounded, we don't know,
                # whats a vertex and whats a ray or line)
                def get_vertex_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = face_iter.facet_repr(output2)
                    return tuple(v(output2[t]) for t in range(leng))

                def get_facet_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = face_iter.vertex_repr(output1)
                    return addtuple + tuple(h(output1[t]) for t in range(leng))

                def print_dim():
                        return Integer(dim - 1 - next_dim)
                dimension = dim - 1 - dimension
                specialcase = 1
            elif dimension == dim - 1 and self._polar:
                # we have stored the actual facets in bitrep_vertices
                def get_facet_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = face_iter.facet_repr(output2)
                    return addtuple + tuple(h(output2[t]) for t in range(leng))

                def get_vertex_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = face_iter.vertex_repr(output1)
                    return tuple(v(output1[t]) for t in range(leng))

                def print_dim():
                    return Integer(next_dim)
                specialcase = 1
            elif self._polar:
                def get_vertex_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = face_iter.facet_repr(output2)
                    return tuple(v(output2[t]) for t in range(leng))

                def get_facet_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = face_iter.vertex_repr(output1)
                    return addtuple + tuple(h(output1[t]) for t in range(leng))

                def print_dim():
                        return Integer(dim - 1 - next_dim)
                if dimension > -2:
                    dimension = dim - 1 - dimension
            else:
                def get_facet_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = face_iter.facet_repr(output2)
                    return addtuple + tuple(h(output2[t]) for t in range(leng))

                def get_vertex_repr():
                    cdef size_t leng
                    cdef size_t t
                    leng = face_iter.vertex_repr(output1)
                    return tuple(v(output1[t]) for t in range(leng))

                def print_dim():
                    return Integer(next_dim)

        # settling the kind of output the user wants once and for all
        if vertex_repr and facet_repr and give_dimension:
            def generate_output():
                return (get_vertex_repr(), get_facet_repr(),
                        print_dim())
        elif vertex_repr and facet_repr:
            def generate_output():
                return (get_vertex_repr(), get_facet_repr())
        elif vertex_repr and give_dimension:
            def generate_output():
                return (get_vertex_repr(), print_dim())
        elif vertex_repr:
            def generate_output():
                return get_vertex_repr()
        elif facet_repr and give_dimension:
            def generate_output():
                return (get_facet_repr(), print_dim())
        elif facet_repr:
            def generate_output():
                return get_facet_repr()
        else:
            def generate_output():
                return print_dim()

        if specialcase:
            # in this case it makes much more sense to use the Bitrep vertices
            vertices = self.bitrep_vertices
            face_iter = FaceIterator(vertices, dim, self._nr_lines)
            face_iter.set_record_dimension(dim - 1)
            output1 = face_iter.get_output1_array()
            output2 = face_iter.get_output2_array()
            next_dim = face_iter.next_face()
            while next_dim is not dim:
                yield generate_output()
                next_dim = face_iter.next_face()
            return

        if self._all_faces:
            f_vector = self._f_vector
            output1 = all_faces.get_output1_array()
            output2 = all_faces.get_output2_array()
            if dimension == -2:
                for next_dim in range(dim):
                    for index in range(f_vector[next_dim + 1]):
                        yield generate_output()
            else:
                next_dim = dimension
                for index in range(f_vector[next_dim + 1]):
                    yield generate_output()
        else:
            face_iter = FaceIterator(facets, dim, self._nr_lines)
            face_iter.set_record_dimension(dimension)
            output1 = face_iter.get_output1_array()
            output2 = face_iter.get_output2_array()
            next_dim = face_iter.next_face()
            while next_dim is not dim:
                yield generate_output()
                next_dim = face_iter.next_face()

    def face_lattice(self, vertices=False, facets=False, names=False):
        r"""
        Generates the face-lattice.

        INPUT:

        - ``vertices`` -- if set to ``True`` the elements in the lattice
          will be named according to the vertices they contain.

        - ``facets`` -- if set to ``True`` the elements in the lattice
          will be named according to the facets they are contained in.

        - ``names`` -- if set to ``False``, facets and vertices will be
          named according to their index.

        Both vertices and facets can be set to ``True``.

        In the case of a trivial Polyhedron, which is equal to its own
        affine hull, ``facets`` will be set to ``False``, as the
        elements need distinct names.

        In the case of a half-space ``vertices`` and ``facets`` will be
        set to ``False``.

        OUTPUT:

        - a :class:'~sage.combinat.posets.lattices.FiniteLatticePoset'

        .. NOTE::

            As :class:'~sage.combinat.posets.lattices.FiniteLatticePoset'
            is awfully slow with elements having meaningful labels,
            the default of this function is to not do so.


        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice()
            Finite lattice containing 5 elements
            sage: C.face_lattice(vertices=True, names=True).atoms()
            [(A vertex at (0, 0),)]

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: P1 = Polyhedron(rays=[[1,0], [-1,0]])
            sage: C1 = CombinatorialPolyhedron(P1)
            sage: C.face_lattice().is_isomorphic(C1.face_lattice())
            True

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice()
            Finite lattice containing 542 elements

        TESTS::

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice().is_isomorphic(P.face_lattice())
            True

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice().is_isomorphic(P.face_lattice())
            True
        """
        cdef size_t **incidences
        cdef size_t nr_incidences
        cdef size_t len_edgelist = self._length_edges_list
        cdef size_t *f_vector
        cdef int k
        cdef int l
        cdef int dim
        if not self._face_lattice_incidences:
            self._calculate_face_lattice_incidences()
        if self._face_lattice_incidences is NULL:
            raise TypeError("could not determine face lattice")
        incidences = self._face_lattice_incidences
        nr_incidences = self._nr_face_lattice_incidences

        # in trivial case,
        # we must ignore part of the input to ensure an injective relabeling
        if self.is_trivial == 1:
            facets = False
        if self.is_trivial == 2:
            vertices = False
            facets = False

        dim = self._dimension
        self._calculate_f_vector()
        f_vector = self._f_vector
        polar_correct = tuple((sum(f_vector[k] for k in range(0, l)),
                               sum(f_vector[k] for k in range(l+1, dim+2)))
                              for l in range(dim+2))

        # the incidences recorded assume we are not in the polar case
        # in the polar case, we want to somewhat reverse the order of the faces
        # the j-th face of dimension i will be the j-th face of dimension
        # `self.dimension -1 -j`

        if self._polar:
            def pol(size_t i):
                for l in range(1, dim + 2):
                    if i < polar_correct[l][0]:
                        # so i is of dimension `l` in the polar version
                        # we will correct its entry to correspond to dimension
                        # `self.dimension() - 1 - l`
                        return Integer(i - polar_correct[l-1][0]
                                       + polar_correct[l-1][1])
                return Integer(i - polar_correct[dim+1][0]
                               + polar_correct[dim+1][1])
        else:
            def pol(size_t i): return Integer(i)

        def face_one(size_t i):
            return pol(incidences[i // len_edgelist][2*(i % len_edgelist)])

        def face_two(size_t i):
            return pol(incidences[i // len_edgelist][2*(i % len_edgelist)+1])

        if self._polar:
            edges = tuple((face_two(j), face_one(j))
                          for j in range(nr_incidences))
        else:
            edges = tuple((face_one(j), face_two(j))
                          for j in range(nr_incidences))

        V = tuple(range(sum(f_vector[k] for k in range(dim+2))))
        D = DiGraph([V, edges], format='vertices_and_edges')
        if vertices or facets:
            dic = {}
            real_f = self.f_vector()
            counters = list(sum(real_f[:i]) for i in range(dim+2))
            for i in self.face_iter(vertex_repr=vertices, facet_repr=facets,
                                    names=names, give_dimension=True):
                d = i[-1]
                if len(i) > 2:
                    dic[counters[d+1]] = i[:-1]
                else:
                    dic[counters[d+1]] = i[0]
                counters[d+1] += 1
            D.relabel(dic)
        return FiniteLatticePoset(D)

    cdef int _calculate_f_vector(self) except 0:
        r"""
        Calculates the f_vector of the `CombinatorialPolyhedron`.

        In the polar case, calculates the f_vector of the polar.

        See :meth:`f_vector`.
        """
        if self._f_vector:
            return 1  # there is no need to recalculate the f_vector
        cdef FaceIterator face_iter
        cdef ListOfFaces facets
        cdef int dim = self.dimension()
        cdef int d
        cdef size_t *f_vector
        cdef MemoryAllocator mem = MemoryAllocator()
        f_vector = <size_t *> mem.malloc((dim + 2) * sizeof(size_t))
        for i in range(dim + 2):
            f_vector[i] = 0
        f_vector[0] = 1
        f_vector[dim + 1] = 1

        if self.is_trivial == 1:
            # Polyhedron equal its affine hull
            # on sucess copying the f_vector to `CombinatorialPolyhedron`
            sig_block()
            self._f_vector = f_vector
            self._MemoryAllocators += (mem,)
            sig_unblock()
            return 1
        if self.is_trivial == 2:
            # Polyhedron of a half space
            f_vector[dim] = 1
            # on sucess copying the f_vector to `CombinatorialPolyhedron`
            sig_block()
            self._f_vector = f_vector
            self._MemoryAllocators += (mem,)
            sig_unblock()
            return 1

        facets = self.bitrep_facets
        face_iter = FaceIterator(facets, dim, self._nr_lines)
        f_vector[0] = 1
        f_vector[dim + 1] = 1
        d = face_iter.next_face()
        while (d < dim):
            f_vector[d+1] += 1
            d = face_iter.next_face()

        # on sucess copying the f_vector to `CombinatorialPolyhedron`
        sig_block()
        self._f_vector = f_vector
        self._MemoryAllocators += (mem,)
        sig_unblock()
        return 1

    cdef int _calculate_edges(self) except 0:
        r"""
        Calculates the edges of the `CombinatorialPolyhedron`.

        In the polar case, calculates the edges of the polar.
        This will also calculate the `f_vector`, if its not already calculated.

        See :meth:`edges` and :meth:`ridges`.
        """
        if self._edges:
            return 1  # there is no need to recalculate the edges
        cdef size_t len_edgelist = self._length_edges_list
        cdef int dim = self.dimension()
        cdef size_t counter
        cdef size_t one
        cdef size_t two
        cdef size_t *output
        cdef int d
        cdef int j
        cdef size_t i
        cdef size_t **edges
        cdef size_t *f_vector
        cdef FaceIterator face_iter
        cdef ListOfFaces facets
        cdef ListOfFaces vertices
        cdef int is_f_vector
        cdef MemoryAllocator mem = MemoryAllocator()
        cdef size_t location_of_mem
        cdef size_t current_length # dynamically enlarge **edges

        self._edges = NULL
        if self._f_vector:
            is_f_vector = 1
        else:
            # in this case we will calculate the f_vector while we're at it
            is_f_vector = 0
        counter = 0
        output = NULL


        edges = <size_t**> mem.malloc(sizeof(size_t*))
        location_of_mem = mem.n - 1
        current_length = 1
        # I will be messing with MemoryAllocator, as it lacks realloc
        # long term one should add realloc to it

        if self.is_trivial > 0:
            self._nr_edges = 0
            sig_block()
            self._MemoryAllocators += (mem,)
            self._edges = edges
            sig_unblock()
            return 1

        if dim < 1:
            self._nr_edges = 0
            sig_block()
            self._MemoryAllocators += (mem,)
            self._edges = edges
            sig_unblock()
            return 1

        if dim == 1:
            self._nr_edges = 1
            edges[0] = <size_t *> mem.malloc(2 * sizeof(size_t))
            edges[0][0] = 0
            edges[0][1] = 1
            sig_block()
            self._MemoryAllocators += (mem,)
            self._edges = edges
            sig_unblock()
            return 1

        facets = self.bitrep_facets
        vertices = self.bitrep_vertices
        face_iter = FaceIterator(facets, dim, self._nr_lines)
        output = face_iter.get_output1_array()

        if is_f_vector:
            face_iter.set_record_dimension(1)
            while (face_iter.next_face() == 1):
                face_iter.vertex_repr(output)
                one = counter // len_edgelist
                two = counter % len_edgelist
                if unlikely(two == 0):
                    if unlikely(one + 1 > current_length):
                        sig_block()
                        current_length *= 2
                        mem.pointers[location_of_mem] = \
                            sig_realloc(edges, current_length*sizeof(size_t*))
                        edges = <size_t **> mem.pointers[location_of_mem]
                        # messing with MemoryAllocator
                        sig_unblock()
                    edges[one] = <size_t *> \
                        mem.malloc(2 * len_edgelist * sizeof(size_t))
                edges[one][2*two] = output[0]
                edges[one][2*two + 1] = output[1]
                counter += 1

            self._nr_edges = counter
            sig_block()
            self._edges = edges
            self._MemoryAllocators += (mem,)
            sig_unblock()
            return 1
        else:
            # while doing the edges one might as well do the f_vector
            f_vector = <size_t *> mem.malloc((dim + 2) * sizeof(size_t))
            for j in range(dim + 2):
                f_vector[j] = 0
            f_vector[0] = 1
            f_vector[dim + 1] = 1

            counter = 0
            d = face_iter.next_face()
            while (d < dim):
                f_vector[d+1] += 1
                if d == 1:
                    face_iter.vertex_repr(output)
                    one = counter // len_edgelist
                    two = counter % len_edgelist
                    if unlikely(two == 0):
                        if unlikely(one + 1 > current_length):
                            sig_block()
                            current_length *= 2
                            mem.pointers[location_of_mem] = \
                                sig_realloc(edges, current_length*sizeof(size_t*))
                            edges = <size_t **> mem.pointers[location_of_mem]
                            # messing with MemoryAllocator
                            sig_unblock()
                        edges[one] = <size_t *> \
                            mem.malloc(2 * len_edgelist * sizeof(size_t))
                    edges[one][2*two] = output[0]
                    edges[one][2*two + 1] = output[1]
                    counter += 1

                d = face_iter.next_face()

            self._nr_edges = counter
            sig_block()
            self._f_vector = f_vector
            self._edges = edges
            self._MemoryAllocators += (mem,)
            sig_unblock()
            return 1

    cdef int _calculate_ridges(self) except 0:
        r"""
        Calculates the ridges of the `CombinatorialPolyhedron`.

        In the polar case, calculates the ridges of the polar.

        See :meth:`edges` and :meth:`ridges`.
        """
        if self._ridges:
            return 1  # there is no need to recalculate the ridges
        cdef size_t len_edgelist = self._length_edges_list
        cdef int dim = self.dimension()
        cdef size_t counter
        cdef size_t one
        cdef size_t two
        cdef size_t *output
        cdef FaceIterator face_iter
        cdef ListOfFaces facets
        cdef int is_f_vector
        cdef MemoryAllocator mem = MemoryAllocator()
        cdef size_t location_of_mem
        cdef size_t **ridges
        cdef size_t current_length # dynamically enlarge **ridges

        self._ridges = NULL
        counter = 0

        ridges = <size_t**> mem.malloc(sizeof(size_t*))
        location_of_mem = mem.n - 1
        current_length = 1
        # I will be messing with MemoryAllocator, as it lacks realloc
        # long term one should add realloc to it

        if self.is_trivial > 0:
            self._nr_ridges = 0
            sig_block()
            self._MemoryAllocators += (mem,)
            self._ridges = ridges
            sig_unblock()
            return 1

        if dim == 1:
            self._nr_ridges = 1
            ridges[0] = <size_t *> mem.malloc(2 * sizeof(size_t))
            ridges[0][0] = 0
            ridges[0][1] = 1
            sig_block()
            self._MemoryAllocators += (mem,)
            self._ridges = ridges
            sig_unblock()
            return 1

        facets = self.bitrep_facets
        face_iter = FaceIterator(facets, dim, self._nr_lines)
        face_iter.set_record_dimension(dim - 2)
        output = face_iter.get_output2_array()
        while (face_iter.next_face() == dim - 2):
            # iterating through all ridges and recording them
            face_iter.facet_repr(output)
            one = counter // len_edgelist
            two = counter % len_edgelist
            if unlikely(two == 0):
                if unlikely(one + 1 > current_length):
                    sig_block()
                    current_length *= 2
                    mem.pointers[location_of_mem] = \
                        sig_realloc(ridges, current_length*sizeof(size_t*))
                    ridges = <size_t **> mem.pointers[location_of_mem]
                    # messing with MemoryAllocator
                    sig_unblock()
                ridges[one] = \
                    <size_t *> mem.malloc(2*len_edgelist*sizeof(size_t))
            ridges[one][2*two] = output[0]
            ridges[one][2*two + 1] = output[1]
            counter += 1

        self._nr_ridges = counter
        sig_block()
        self._ridges = ridges
        self._MemoryAllocators += (mem,)
        sig_unblock()
        return 1

    cdef int _calculate_face_lattice_incidences(self) except 0:
        r"""
        Calculates all incidences for the face lattice.

        See :meth:`face_lattice`.
        """
        if self._face_lattice_incidences:
            return 1  # there is no need to recalculate the incidences
        cdef size_t len_edgelist = self._length_edges_list
        cdef int dim = self.dimension()
        cdef size_t counter
        cdef size_t one
        cdef size_t two
        cdef size_t *output
        cdef int j
        cdef size_t i
        cdef ListOfAllFaces all_faces
        cdef size_t first = 0
        cdef size_t second = 0
        cdef int dimension_two
        cdef int dimension_one
        cdef size_t *f_vector
        cdef size_t already_seen
        cdef size_t already_seen_next
        cdef MemoryAllocator mem = MemoryAllocator()
        cdef size_t location_of_mem
        cdef size_t **incidences
        cdef size_t current_length # dynamically enlarge **incidences

        self._record_all_faces()
        all_faces = self._all_faces
        f_vector = self._f_vector
        if all_faces is None and self.is_trivial == 0:
            raise ValueError("could not determine a list of all faces")

        incidences = <size_t**> mem.malloc(sizeof(size_t*))
        location_of_mem = mem.n - 1
        current_length = 1
        # I will be messing with MemoryAllocator, as it lacks realloc
        # long term one should add realloc to it

        if self.is_trivial == 1:
            # the case of a Polyhedron with only two faces
            incidences[0] = <size_t *> mem.malloc(2 * sizeof(size_t))
            incidences[0][0] = 0
            incidences[0][1] = 1
            self._nr_face_lattice_incidences = 1
            sig_block()
            self._MemoryAllocators += (mem,)
            self._face_lattice_incidences = incidences
            sig_unblock()
            return 1

        if self.is_trivial == 2:
            # the case of a Polyhedron with only three face
            # the Polyhedron is a half-space
            incidences[0] = <size_t *> mem.malloc(4 * sizeof(size_t))
            incidences[0][0] = 0
            incidences[0][1] = 1
            incidences[0][2] = 1
            incidences[0][3] = 2
            self._nr_face_lattice_incidences = 2
            sig_block()
            self._MemoryAllocators += (mem,)
            self._face_lattice_incidences = incidences
            sig_unblock()
            return 1

        counter = 0

        dimension_one = 0
        while (f_vector[dimension_one + 1] == 0):
            # taking care of cases, where there might be no faces
            # of dimension 0
            dimension_one += 1
        dimension_two = -1
        while (dimension_one < dim + 1):
            already_seen = \
                sum(f_vector[j] for j in range(dimension_two + 1))
            already_seen_next = already_seen + f_vector[dimension_two + 1]
            dimension_one = dimension_two + 1
            all_faces.incidence_init(dimension_one, dimension_two)
            while all_faces.next_incidence(&second, &first):
                second += already_seen_next
                first += already_seen
                one = counter // len_edgelist
                two = counter % len_edgelist
                if unlikely(two == 0):
                    if unlikely(one + 1 > current_length):
                        sig_block()
                        current_length *= 2
                        mem.pointers[location_of_mem] = \
                            sig_realloc(incidences,
                                        (current_length)*sizeof(size_t*))
                        incidences = <size_t **> mem.pointers[location_of_mem]
                        # messing with MemoryAllocator
                        sig_unblock()
                    incidences[one] = <size_t *> \
                        mem.malloc(2 * len_edgelist * sizeof(size_t))
                incidences[one][2*two] = first
                incidences[one][2*two + 1] = second
                counter += 1
                sig_check()
            dimension_one += 1
            dimension_two = dimension_one - 1

        self._nr_face_lattice_incidences = counter
        sig_block()
        self._MemoryAllocators += (mem,)
        self._face_lattice_incidences = incidences
        sig_unblock()
        return 1

    def _record_all_faces(self):
        r"""
        Records and sorts all faces of the Polyhedron for quicker acces later.

        EXAMPLES::

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C._record_all_faces()

        TESTS::

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: tup = tuple(i for i in C.face_iter(facet_repr=True))
            sage: C._record_all_faces()
            sage: tup1 = tuple(i for i in C.face_iter(facet_repr=True))
            sage: sorted(tup) == sorted(tup1)
            True

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: tup = tuple(i for i in C.face_iter(facet_repr=True))
            sage: C._record_all_faces()
            sage: tup1 = tuple(i for i in C.face_iter(facet_repr=True))
            sage: sorted(tup) == sorted(tup1)
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: tup = tuple(i for i in C.face_iter(facet_repr=True))
            sage: C._record_all_faces()
            sage: tup1 = tuple(i for i in C.face_iter(facet_repr=True))
            sage: sorted(tup) == sorted(tup1)
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0],
            ....:                      [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: tup = tuple(i for i in C.face_iter(facet_repr=True))
            sage: C._record_all_faces()
            sage: tup1 = tuple(i for i in C.face_iter(facet_repr=True))
            sage: sorted(tup) == sorted(tup1)
            True
        """
        if self.is_trivial >= 1:
            return

        if self._all_faces:
            return

        self._record_all_faces_helper()
        if self.is_trivial == 0 and self._all_faces is None:
            raise ValueError("could not determine a list of all faces")

    cdef int _record_all_faces_helper(self) except 0:
        r"""
        Records and sorts all faces of `CombinatorialPolyhedron`.

        See :meth:`_record_all_faces`.
        """
        cdef int dim = self.dimension()
        cdef tuple f_tuple
        cdef ListOfAllFaces all_faces
        cdef uint64_t *face
        cdef FaceIterator face_iter
        cdef ListOfFaces facets
        cdef int d
        self._calculate_f_vector()
        if not self._f_vector:
            raise TypeError("could not determine f_vector")
        if self.is_trivial:
            return 1  # in this case we will not record all faces
        f_tuple = tuple(self._f_vector[i] for i in range(dim + 2))
        all_faces = ListOfAllFaces(self.bitrep_facets, dim, f_tuple)
        facets = self.bitrep_facets
        face_iter = FaceIterator(facets, dim, self._nr_lines)
        d = face_iter.next_face()
        while (d == dim - 1):
            d = face_iter.next_face()
        while (d < dim):
            all_faces.add_face(d, face_iter.face)
            d = face_iter.next_face()
        all_faces.sort()
        self._all_faces = all_faces
        return 1
