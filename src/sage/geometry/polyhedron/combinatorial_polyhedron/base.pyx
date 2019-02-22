#distutils: language = c++

r"""
CombinatorialPolyhedron gathers several algorithms of Polyhedra depending only
on the vertex-facet incidences.

Most importanly, this module offers a quick face iterator, which can be used
to calculate f_vector, edges, ridges and even the face lattice.

The FaceIterator would work on every atomic and coatomic lattice, where every
inteverl of length two has exactly 4 elements (known as the diamond property).

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

from libc.string cimport memcmp, memcpy, memset
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

cdef inline uint64_t vertex_to_bit_dictionary(size_t i):
    r"""
    Returns an uint64_t with exactly the i-th bit set to 1.

    This "dictionary" helps storign a vector of 64 incidences as uint64_t
    """
    return (<uint64_t>1) << (64 - i - 1)

cdef int char_from_tuple(tuple tup, uint64_t *output,
                         size_t face_length) except -1:
    r"""
    Converts a tuple into a bit-representation. Stores it in `output`.

    The first bit represent the entry `0` and is set to one, iff `0` is in
    `tup`. The second bit represents `1` and so on.

    INPUT:

    - `tup` -- tuple of pairwise distinct positive integers in
      `range(face_length*64)`.
    - `output` -- array of `uint64_t` of length `face_length`.
    - `face_length`.

    OUTPUT:

    - `1` on sucess
    - fills and initializes `output`.

    EXAMPLES::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport char_from_tuple
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from sage.rings.integer import Integer
        ....:
        ....: def char_from_tuple_wrapper(tup):
        ....:     cdef size_t face_length = max(tup)//64 + 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = <uint64_t *> mem.allocarray(face_length, 8)
        ....:     char_from_tuple(tup, output, face_length)
        ....:     return tuple(Integer(output[i]) for i in range(face_length))
        ....:
        ....: def char_from_tuple_wrong_size(tup):
        ....:     cdef size_t face_length = 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = <uint64_t *> mem.allocarray(face_length, 8)
        ....:     char_from_tuple(tup, output, face_length)
        ....:     return tuple(Integer(output[i]) for i in range(face_length))
        ....: ''')  # long time

        sage: char_from_tuple_wrapper((62, 63))  # long time
        (3,)
        sage: char_from_tuple_wrapper((61, 63, 125))  # long time
        (5, 4)
        sage: char_from_tuple_wrong_size((62, 70))  # long time
        Traceback (most recent call last):
        ..
        IndexError: output too small to represent 70
        sage: char_from_tuple_wrapper((-1, 12))  # long time
        Traceback (most recent call last):
        ..
        OverflowError: can't convert negative value to size_t
        sage: char_from_tuple_wrapper((0, 0))  # long time
        Traceback (most recent call last):
        ..
        ValueError: entries of `tup` are not distinct
    """
    cdef size_t entry  # will iterate over tup
    cdef size_t position  # determines the position in output of entry
    cdef size_t value  # determines which bit will be set in output[position]

    memset(output, 0, face_length*8)
    if unlikely(len(tup) != len(set(tup))):
        raise ValueError("entries of `tup` are not distinct")
    for entry in tup:
        value = entry % 64
        position = entry//64
        if unlikely(position >= face_length):
            raise IndexError("output too small to represent %s"%entry)
        output[position] += vertex_to_bit_dictionary(value)

cdef int char_from_incidence(tuple incidences, uint64_t *output,
                             size_t face_length) except -1:

    r"""
    Converts a tuple of incidences into a bit-representation.

    Stores it in `output`. Each entry in `incidences` represents a bit in
    `output`. It is set to `1`, iff the entry in `incidences` in non-zero.

    INPUT:

    - `tup` -- tuple of integers representing incidences of length at most
      `face_length*64`.
    - `output` -- array of `uint64_t` of length `face_length`.
    - `face_length`.

    OUTPUT:

    - `1` on sucess
    - fills and initializes `output`.

    EXAMPLES::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport char_from_incidence
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from sage.rings.integer import Integer
        ....:
        ....: def char_from_incidences_wrapper(tup):
        ....:     cdef size_t face_length = (len(tup)-1)//64 + 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = \
        ....:          <uint64_t *> mem.allocarray(face_length, 8)
        ....:     char_from_incidence(tup, output, face_length)
        ....:     return tuple(Integer(output[i]) for i in range(face_length))
        ....:
        ....: def char_from_incidences_wrong_size(tup):
        ....:     cdef size_t face_length = 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = \
        ....:         <uint64_t *> mem.allocarray(face_length, 8)
        ....:     char_from_incidence(tup, output, face_length)
        ....:     return tuple(Integer(output[i]) for i in range(face_length))
        ....: ''')  # long time

        sage: char_from_incidences_wrapper((0,) * 62 + (1,1))  # long time
        (3,)
        sage: char_from_incidences_wrapper((0,) * 61 + (1,0,1) +  # long time
        ....:                              (0,) * 61 + (1,))
        (5, 4)
        sage: char_from_incidences_wrapper((1,) * 64)  # long time
        (18446744073709551615,)
        sage: char_from_incidences_wrong_size((1,) * 70)  # long time
        Traceback (most recent call last):
        ..
        IndexError: output too small to represent all incidences
    """
    cdef size_t entry  # will iterate over entries in tup
    cdef size_t position  # determines the position in output of entry
    cdef size_t value  # determines which bit will be set in output[position]
    cdef size_t length = len(incidences)

    memset(output, 0, face_length*8)
    if unlikely(length > 64*face_length):
        raise IndexError("output too small to represent all incidences")
    for entry in range(length):
        if incidences[entry]:
            value = entry % 64
            position = entry//64
            output[position] += vertex_to_bit_dictionary(value)

def get_facets_from_incidence_matrix(matrix):
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
        ....: cimport ListOfFaces, vertex_repr_from_bitrep
        ....: from sage.rings.integer import Integer
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def vertex_repr_from_bitrep_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.nr_vertices,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = vertex_repr_from_bitrep(
        ....:         data, output, faces.face_length)
        ....:     return tuple(Integer(output[i]) for i in range(length))
        ....: ''')

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import get_facets_from_incidence_matrix
        sage: P = polytopes.permutahedron(4)
        sage: facets = get_facets_from_incidence_matrix(P.incidence_matrix())
        sage: facets.nr_faces
        14
        sage: facets.nr_vertices
        24
        sage: for i in range(facets.nr_faces):
        ....:     print(vertex_repr_from_bitrep_wrapper(facets, i))
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
    if unlikely(not is_Matrix(matrix)):
        raise ValueError("input must be matrix")
    # transpose and get rid of equalities (which all vertices satisfie)
    matrix = matrix.transpose()
    rg = range(matrix.nrows())
    matrix = matrix[list(i for i in rg if not all(j for j in matrix[i]))]
    cdef ListOfFaces facets
    cdef size_t length = matrix.nrows()
    facets = ListOfFaces(length, matrix.ncols())
    cdef uint64_t **facets_data = facets.data
    cdef size_t i
    for i in range(length):
        # filling each facet with the corresponding data, which is the
        # "i-th column" of the original matrix
        char_from_incidence(tuple(matrix[i]), facets_data[i],
                            facets.face_length)
    return facets

def get_vertices_from_incidence_matrix(matrix):
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
        ....: cimport ListOfFaces, vertex_repr_from_bitrep
        ....: from sage.rings.integer import Integer
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def vertex_repr_from_bitrep_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.nr_vertices,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = vertex_repr_from_bitrep(
        ....:         data, output, faces.face_length)
        ....:     return tuple(Integer(output[i]) for i in range(length))
        ....: ''')

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import get_vertices_from_incidence_matrix
        sage: P = polytopes.permutahedron(4)
        sage: vertices = get_vertices_from_incidence_matrix(P.incidence_matrix())
        sage: vertices.nr_faces
        24
        sage: vertices.nr_vertices
        14
        sage: for i in range(vertices.nr_faces):
        ....:     print(vertex_repr_from_bitrep_wrapper(vertices, i))
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
    if unlikely(not is_Matrix(matrix)):
        raise ValueError("input must be matrix")
    # transpose and get rid of equalities (which all vertices satisfie)
    matrix = matrix.transpose()
    rg = range(matrix.nrows())
    matrix = matrix[list(i for i in rg if not all(j for j in matrix[i]))]
    cdef ListOfFaces vertices
    cdef size_t length = matrix.ncols()
    vertices = ListOfFaces(length, matrix.nrows())
    cdef uint64_t **vertices_data = vertices.data
    cdef size_t i
    for i in range(length):
        # filling each vertex with the corresponding data, which is the
        # "i-th row" of the original matrix
        char_from_incidence(tuple(matrix.column(i)), vertices_data[i],
                            vertices.face_length)
    return vertices

def get_facets_bitrep_from_facets_tuple(tuple facets_input, size_t nr_vertices):
    r"""
    Initializes facets in bit-representation as `ListOfFaces`.

    INPUT:

    - `facets_input` -- a tuple of facets, each facet a tuple of vertices.
      The vertices must be exactly range(nr_vertices).
    - `nr_vertices`.

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, vertex_repr_from_bitrep
        ....: from sage.rings.integer import Integer
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def vertex_repr_from_bitrep_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.nr_vertices,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = vertex_repr_from_bitrep(
        ....:         data, output, faces.face_length)
        ....:     return tuple(Integer(output[i]) for i in range(length))
        ....: ''')

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import get_facets_bitrep_from_facets_tuple
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: facets = get_facets_bitrep_from_facets_tuple(bi_pyr, 6)
        sage: for i in range(8):
        ....:     print(vertex_repr_from_bitrep_wrapper(facets, i))
        (0, 1, 4)
        (1, 2, 4)
        (2, 3, 4)
        (0, 3, 4)
        (0, 1, 5)
        (1, 2, 5)
        (2, 3, 5)
        (0, 3, 5)
    """
    cdef Py_ssize_t i
    cdef ListOfFaces facets = ListOfFaces(len(facets_input),
                                          nr_vertices)
    cdef size_t face_length = facets.face_length
    cdef uint64_t **facets_data = facets.data
    for i in range(len(facets_input)):
        # filling each facet with the the data from the corresponding facet
        char_from_tuple(facets_input[i], facets_data[i], face_length)
    return facets

def get_vertices_bitrep_from_facets_tuple(tuple facets_input, size_t nr_vertices):
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
        ....: cimport ListOfFaces, vertex_repr_from_bitrep
        ....: from sage.rings.integer import Integer
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def vertex_repr_from_bitrep_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.nr_vertices,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = vertex_repr_from_bitrep(
        ....:         data, output, faces.face_length)
        ....:     return tuple(Integer(output[i]) for i in range(length))
        ....: ''')

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import get_vertices_bitrep_from_facets_tuple
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: vertices = get_vertices_bitrep_from_facets_tuple(bi_pyr, 6)
        sage: for i in range(6):
        ....:     print(vertex_repr_from_bitrep_wrapper(vertices, i))
        (0, 3, 4, 7)
        (0, 1, 4, 5)
        (1, 2, 5, 6)
        (2, 3, 6, 7)
        (0, 1, 2, 3)
        (4, 5, 6, 7)
    """
    cdef size_t i
    cdef tuple new_input
    cdef size_t input_vertex
    # each input_vertex will correspond to a "facet" in the facet-representation
    # of the vertices
    cdef size_t input_facet
    # each input_facet will correspond to a "vertex" in the facet-representation
    # of the vertices
    cdef size_t inputlength
    cdef size_t position
    # determining at which position to set a bit for the input_facet
    cdef size_t value  # determing which bit to set for the input_facet
    cdef ListOfFaces vertices = ListOfFaces(nr_vertices,
                                            len(facets_input))
    cdef uint64_t **vertices_data = vertices.data
    cdef size_t face_length = vertices.face_length
    inputlength = len(facets_input)

    for i in range(nr_vertices):
        # initializing the data of ListOfFaces
        memset(vertices_data[i], 0, face_length*8)

    for input_facet in range(inputlength):
        value = input_facet % 64
        position = input_facet//64
        for input_vertex in facets_input[input_facet]:
            # if the input_vertex is in the input_facet,
            # than in facet_repr of vertices
            # input_facet is a vertex of input_vertex
            vertices_data[input_vertex][position] += \
                vertex_to_bit_dictionary(value)
    return vertices

cdef inline size_t vertex_repr_from_bitrep(uint64_t *face, size_t *output,
                                    size_t face_length) except -1:
    r"""
    Takes a bitrep-representation and converts it to an array.

    Basically this is an inverse to `char_from_tuple`. Instead of a tuple,
    it stores the vertices in `output`. Returns length of representation.

    INPUT:

    - `face` -- a bit-representation of a face.
    - `output` -- an array of `size_t` long enough to contain all vertices of
      that face. `face_length*64` will suffice.
    - `face_length` -- the length of `face`.

    OUTPUT:

    - returns length of `output`
    - stores vertices in `ouput`

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, vertex_repr_from_bitrep, char_from_tuple
        ....: from sage.rings.integer import Integer
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def vertex_repr_from_bitrep_wrapper(tup):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef length = len(tup)
        ....:     output = <size_t *> mem.allocarray(length*64,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data
        ....:     data = <uint64_t *> mem.allocarray(length, 8)
        ....:     for i in range(len(tup)):
        ....:         data[i] = tup[i]
        ....:     outputlength = vertex_repr_from_bitrep(
        ....:         data, output, length)
        ....:     return tuple(Integer(output[i]) for i in range(outputlength))
        ....: ''')  # long time

        sage: vertex_repr_from_bitrep_wrapper((17, 31))  # long time
        (59, 63, 123, 124, 125, 126, 127)
        sage: vertex_repr_from_bitrep_wrapper((13,))  # long time
        (60, 61, 63)
        sage: vertex_repr_from_bitrep_wrapper((0, 61))  # long time
        (122, 123, 124, 125, 127)

    TESTS:

    Testing that :meth`vertex_repr_from_bitrep` is the
    inverse to :meth:`char_from_tuple`::

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
        ....: ''') # long time
    """
    cdef size_t i, j
    cdef size_t output_length = 0
    cdef uint64_t copy
    for i in range(face_length):
        if face[i]:
            copy = face[i]
            for j in range(64):
                if copy >= vertex_to_bit_dictionary(j):
                    # then face[i] has the j-th bit set to 1
                    # this corresponds face containing vertex i*64 + j
                    output[output_length] = i*64 + j
                    output_length += 1
                    copy -= vertex_to_bit_dictionary(j)
    return output_length

cdef class ListOfFaces(SageObject):
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

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import ListOfFaces
        sage: facets = ListOfFaces(5, 13)
        sage: facets.face_length in (1, 2, 4)
        True
        sage: facets.nr_vertices
        13
        sage: facets.nr_faces
        5
    """
    def __init__(self, size_t nr_faces, size_t nr_vertices):
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
            self.data[i] = <uint64_t *> \
                self._mem.aligned_malloc(chunksize//8, self.face_length*8)
            # we must allocate the memory for ListOfFaces overaligned:
            #     must be 16-byte aligned if chunksize = 128
            #     must be 32-byte aligned if chunksize = 256

cpdef int calculate_dimension(ListOfFaces faces) except -2:
    r"""
    Calculates the dimension of a polyhedron by its facets.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import calculate_dimension, get_facets_bitrep_from_facets_tuple, \
        ....:            get_vertices_bitrep_from_facets_tuple
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: facets = get_facets_bitrep_from_facets_tuple(bi_pyr, 6)
        sage: vertices = get_vertices_bitrep_from_facets_tuple(bi_pyr, 6)
        sage: calculate_dimension(facets)
        3
        sage: calculate_dimension(vertices)
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

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import calculate_dimension, get_facets_from_incidence_matrix, \
        ....:            get_vertices_from_incidence_matrix
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: for _ in range(10):
        ....:     points = tuple(tuple(randint(-1000,1000) for _ in range(10))
        ....:                    for _ in range(randint(3,15)))
        ....:     P = Polyhedron(vertices=points)
        ....:     facets = get_facets_from_incidence_matrix(P.incidence_matrix())
        ....:     vertices = get_vertices_from_incidence_matrix(P.incidence_matrix())
        ....:     d1 = P.dimension()
        ....:     d2 = calculate_dimension(facets)
        ....:     d3 = calculate_dimension(vertices)
        ....:     if not d1 == d2 == d3:
        ....:         print('calculation_dimension() seems to be incorrect')
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
    # newfaces will contain the intersections of the faces with the last face
    cdef MemoryAllocator newfaces2
    cdef uint64_t **newfaces2data
    # newfaces2 will point to those faces in newfaces, that are facets

    if nr_faces == 0:
            raise TypeError("wrong usage of `calculate_dimension_loop`,\n" +
                            "at least one face needed.")

    if nr_faces == 1:
        # we expect the face to be the empty Polyhedron
        # possibly it contains more than one vertex/rays/lines
        # the dimension of the Polyhedron with this face as only facet is
        # `bitcount`
        bitcount = CountFaceBits(facesdata[0], face_length)
        return bitcount

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
    return calculate_dimension_loop(newfaces2data, new_nr_faces,
                                    face_length) + 1

cdef class FaceIterator(SageObject):
    r"""
    A class to iterate over all faces of a Polyhedron.

    INPUT:

    - `C` -- a `CombinatorialPolyhedron`.
    - `dual` -- if True, then the dual Polyhedron is used for iteration.
      Only possible if Polyhedron is bounded.
      If `dual`, then vertices will be visited first and other faces
      constructed from the vertices.
      This method is faster for less vertices then facets.

    .. SEEALSO::

        :class:`CombinatorialPolyhedron`.

    EXAMPLES:

    Construct a FaceIterator::

        sage: P = polytopes.cuboctahedron()
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_iter()

    By default it will give the dimension of each face::

        sage: next(it)
        0
        sage: next(it)
        0
        sage: next(it)
        0
        sage: next(it)
        0
        sage: next(it)
        0
        sage: next(it)
        0
        sage: next(it)
        0
        sage: next(it)
        0
        sage: next(it)
        0
        sage: next(it)
        0
        sage: next(it)
        0
        sage: next(it)
        0
        sage: next(it)
        1
        sage: next(it)
        1

    Get more knowledge about current face::

        sage: it.length_vertex_repr()
        2
        sage: it.vertex_repr()
        (A vertex at (1, 0, -1), A vertex at (1, 1, 0))
        sage: it.length_facet_repr()
        2
        sage: it.facet_repr()
        (An inequality (-1, 0, 0) x + 1 >= 0, An inequality (-1, -1, 1) x + 2 >= 0)
        sage: it.get_dimension()
        1

    Ignore faces the current face contains::

        sage: it.ignore_supfaces()
        sage: next(it)
        1
        sage: next(it)
        1
        sage: next(it)
        2
        sage: next(it)
        2
        sage: next(it)
        1
        sage: it.length_facet_repr()
        2

    Construct faces by the dual or not::

        sage: it = C.face_iter(dual=False)
        sage: it.next()
        2
        sage: it.next()
        2
        sage: it.ignore_subfaces()
        sage: it.ignore_supfaces()
        Traceback (most recent call last):
        ..
        ValueError: only possible when in dual mode
        sage: it = C.face_iter(dual=True)
        sage: it.next()
        0
        sage: it.next()
        0
        sage: it.ignore_supfaces()
        sage: it.ignore_subfaces()
        Traceback (most recent call last):
        ..
        ValueError: only possible when not in dual mode

    ALGORITHM:

    The algorithm to visit all proper faces exactly once is roughly
    equivalent to:

    faces = [set(facet) for facet in self.facets()]
    ComputeNextStep(facets, [])

    # this algorithm assumes at each step to receive all facets of
    # some face except those contained in a face of forbidden

    def ComputeNextStep(faces, forbidden):

        for face in faces:
            pass #do something here with that face

        while len(faces) > 1:
            one_face = faces.pop()
            maybe_newfaces = [one_face.intersection(face) for face in faces]
            # maybe_newfaces contains all intersection with one_face

            newfaces2 = []
            for face1 in newfaces:
                # face1 is a facet of one_face iff
                # it is not contained in another facet
                if all(not face1 < face2 for face2 in newfaces):
                    newfaces2.append(face1)
            # newfaces2 contains all facets of one_face not contained
            # in any one of forbidden and maybe some that are
            # contained in one of forbidden

            newfaces = []
            for face1 in newfaces2:
                if all(not face1 < face2 for face2 in forbidden):
                    newfaces.append(face1)
            # newfaces contains exactly all facets of one_face but
            # those contained in one face of forbidden

            # visit all faces in one_face that are not contained in
            # one of forbidden
            ComputeNextStep(newfaces, forbidden)

            # we have visited all faces in one_face, so we should not
            # visit one ever again

            forbidden.append(one_face)

            return
    """
    def __init__(self, CombinatorialPolyhedron C, int dual):
        r"""
        Initializes `FaceIterator`.

        See :class:`FaceIterator`.

        EXAMPLES::

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: f_vector = [1, 0, 0, 0, 1]
            sage: for d in it: f_vector[d+1] += 1
            sage: print ('f_vector of permutahedron(4): ', f_vector)
            f_vector of permutahedron(4):  [1, 24, 36, 14, 1]
        """
        if dual and C._unbounded:
            raise ValueError("cannot iterate over dual of unbounded Polyedron")
        cdef int i
        cdef ListOfFaces some_list  # just to make Cython aware of type

        self.dual = dual
        self.face = NULL
        self.dimension = C.dimension()
        self.current_dimension = self.dimension -1
        self.nr_lines = C._nr_lines
        self.record_dimension = -2
        self._mem = MemoryAllocator()
        self.lowest_dimension = self.nr_lines
        # we will not yield the empty face, if there are `n` lines, than there
        # are no faces below dimension `n`
        if dual:
            self.atoms = C.bitrep_facets
            self.coatoms = C.bitrep_vertices
        else:
            self.coatoms = C.bitrep_facets
            self.atoms = C.bitrep_vertices
        self.face_length = self.coatoms.face_length
        self.atom_repr_face = <size_t *> \
            self._mem.allocarray(self.coatoms.nr_vertices, sizeof(size_t))
        self.coatom_repr_face = <size_t *> \
            self._mem.allocarray(self.coatoms.nr_faces, sizeof(size_t))
        self._V = C._V
        self._H = C._H
        self._equalities = C._equalities

        if self.dimension == 0 or self.coatoms.nr_faces == 0:
            # as we will only yield prober faces, there is nothing in those cases
            # we have to discontinue with init,
            # as it assume dimension > 0 and nr_faces > 0
            self.current_dimension = self.dimension
            return

        # the place where the possible new faces are being stored
        # note that they are actually located here, whereas newfaces2 only
        # points to the correct ones
        self.newfaces_lists = \
            tuple(ListOfFaces(self.coatoms.nr_faces, self.coatoms.nr_vertices)
                  for i in range(self.dimension -1))
        self.maybe_newfaces = <uint64_t ***> \
            self._mem.allocarray((self.dimension -1), sizeof(uint64_t **))
        for i in range(self.dimension -1):
            some_list = self.newfaces_lists[i]
            self.maybe_newfaces[i] = some_list.data

        # initializing forbidden
        self.forbidden = <uint64_t **> \
            self._mem.allocarray(self.coatoms.nr_faces, sizeof(uint64_t *))
        self.nr_forbidden = \
            <size_t *> self._mem.allocarray(self.dimension, sizeof(size_t))
        self.nr_forbidden[self.dimension -1] = 0

        # initializing newfaces, the pointers to the faces we want to visit
        self.newfaces = <uint64_t ***> \
            self._mem.allocarray(self.dimension, sizeof(uint64_t **))
        for i in range(self.dimension - 1):
            self.newfaces[i] = \
                <uint64_t **> self._mem.allocarray(self.coatoms.nr_faces,
                                                   sizeof(uint64_t *))
        self.newfaces[self.dimension - 1] = self.coatoms.data  # we start with coatoms

        # initialize nr_newfaces
        self.nr_newfaces = \
            <size_t *> self._mem.allocarray(self.dimension, sizeof(size_t))
        self.nr_newfaces[self.dimension - 1] = self.coatoms.nr_faces

        # initializing first time
        self.first_time = \
            <int *> self._mem.allocarray(self.dimension, sizeof(int))
        self.first_time[self.dimension - 1] = 1

        self.yet_to_yield = self.coatoms.nr_faces

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: P = polytopes.associahedron(['A',3])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_iter()
            An Iterator over the faces of a Polyhedron of dimension 3
        """
        return "An Iterator over the faces of a Polyhedron of dimension %s"%self.dimension

    def __next__(self):
        r"""
        Visit the next face and return its dimension.

        EXAMPLES::
            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: next(it)
            2
            sage: next(it)
            2
            sage: next(it)
            2
            sage: next(it)
            2
            sage: next(it)
            2
            sage: next(it)
            2
            sage: next(it)
            1
        """
        cdef int d = self.next_face()
        if unlikely(d == self.dimension):
            raise StopIteration()
        return Integer(self.dual*(self.dimension-1-d) + (1-self.dual)*d)

    next = __next__

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: for d in it: d
            2
            2
            2
            2
            1
            1
            1
            0
            0
            0
            1
            1
            0
            1
        """
        return self

    def set_record_dimension(self, dim):
        r"""
        This will have the iterator only yield faces in dimension `dim`.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: next(it)
            3
            sage: counter = 0
            sage: it.set_record_dimension(2)
            sage: for _ in it: counter += 1
            sage: print ('permutahedron(5) has', counter,
            ....:        'faces of dimension 2')
            permutahedron(5) has 150 faces of dimension 2
            sage: C.f_vector()
            (1, 120, 240, 150, 30, 1)
        """
        if self.dual:
            self.record_dimension = self.dimension - 1 - dim
        else:
            self.record_dimension = dim
        self.lowest_dimension = max(self.nr_lines, self.record_dimension)

    def get_dimension(self):
        r"""
        Returns the dimension of the current face.

        EXAMPLES::

            sage: P = polytopes.associahedron(['A', 3])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: next(it)
            2
            sage: it.get_dimension()
            2
            sage: all(d == it.get_dimension() for d in it)
            True
        """
        if unlikely(self.face is NULL):
            raise ValueError("iterator not set to a face yet")
        if unlikely(self.current_dimension == self.dimension):
            raise ValueError("iterator consumed")
        return Integer(self.dual*(self.dimension-1-self.current_dimension) +
                       (1-self.dual)*self.current_dimension)

    def vertex_repr(self, names=True):
        r"""
        Returns the vertex-representation of the current face.

        The vertex-representation consists of
        the vertices/rays/lines that face contains.

        INPUT:

        - `names` -- if `True` returns the names of the vertices/rays/lines
          as given on initialization of the `CombinatorialPolyhedron`.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(dimension=2)
            sage: next(it)
            2
            sage: it.vertex_repr()
            (A vertex at (4, 1, 5, 2, 3),
             A vertex at (4, 2, 5, 1, 3),
             A vertex at (5, 1, 4, 2, 3),
             A vertex at (5, 2, 4, 1, 3))
            sage: next(it)
            2
            sage: it.vertex_repr()
            (A vertex at (4, 1, 5, 2, 3),
             A vertex at (4, 1, 5, 3, 2),
             A vertex at (5, 1, 4, 2, 3),
             A vertex at (5, 1, 4, 3, 2))
            sage: next(it)
            2
            sage: it.vertex_repr(False)
            (76, 77, 82, 83, 88, 89)
            sage: next(it)
            2
            sage: it.vertex_repr(False)
            (77, 83, 101, 107)

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],[0,2,3],[1,2,3]])
            sage: it = C.face_iter()
            sage: for i in it: (i, it.vertex_repr())
            (2, (1, 2, 3))
            (2, (0, 2, 3))
            (2, (0, 1, 3))
            (2, (0, 1, 2))
            (1, (2, 3))
            (1, (1, 3))
            (1, (1, 2))
            (0, (3,))
            (0, (2,))
            (0, (1,))
            (1, (0, 3))
            (1, (0, 2))
            (0, (0,))
            (1, (0, 1))
        """
        cdef size_t length
        if unlikely(self.face is NULL):
            raise ValueError("iterator not set to a face yet")
        if self.dual:
            length = self.coatom_repr()
            if names and self._V:
                return tuple(self._V[self.coatom_repr_face[i]]
                             for i in range(length))
            else:
                return tuple(Integer(self.coatom_repr_face[i])
                             for i in range(length))
        else:
            length = self.atom_repr()
            if names and self._V:
                return tuple(self._V[self.atom_repr_face[i]]
                             for i in range(length))
            else:
                return tuple(Integer(self.atom_repr_face[i])
                             for i in range(length))

    def length_vertex_repr(self):
        r"""
        Returns the length of the vertex_repr.

        Might be faster then `len(self.vertex_repr())`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: for i in it: (i, it.length_vertex_repr())
            (2, 4)
            (2, 4)
            (2, 4)
            (2, 4)
            (2, 4)
            (2, 4)
            (1, 2)
            (1, 2)
            (1, 2)
            (1, 2)
            (0, 1)
            (0, 1)
            (0, 1)
            (0, 1)
            (1, 2)
            (1, 2)
            (1, 2)
            (0, 1)
            (0, 1)
            (1, 2)
            (1, 2)
            (0, 1)
            (1, 2)
            (1, 2)
            (0, 1)
            (1, 2)
        """
        if unlikely(self.face is NULL):
            raise ValueError("iterator not set to a face yet")
        if self.dual:
            return Integer(self.coatom_repr())
        else:
            return Integer(self.length_atom_repr())

    def facet_repr(self, names=True):
        r"""
        Returns the facet-representation of the current face.

        The facet-representation consists of the facets
        that contain the face and of the equalities of the Polyhedron.

        INPUT:

        - `names` -- if `True` returns the names of the facets/equations
          as given on initialization of the `CombinatorialPolyhedron`.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(2)
            sage: next(it)
            2
            sage: it.facet_repr()
            (An inequality (0, 1, 0, 1, 0) x - 3 >= 0,
             An inequality (0, 1, 0, 1, 1) x - 6 >= 0,
             An equation (1, 1, 1, 1, 1) x - 15 == 0)
            sage: next(it)
            2
            sage: it.facet_repr()
            (An inequality (0, 1, 0, 0, 0) x - 1 >= 0,
             An inequality (0, 1, 0, 1, 1) x - 6 >= 0,
             An equation (1, 1, 1, 1, 1) x - 15 == 0)
            sage: next(it)
            2
            sage: it.facet_repr(False)
            (12, 29)
            sage: next(it)
            2
            sage: it.facet_repr(False)
            (6, 29)

            sage: P = polytopes.cyclic_polytope(4,6)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: next(it)
            0
            sage: it.facet_repr()
            (An inequality (-20, 29, -10, 1) x + 0 >= 0,
             An inequality (60, -47, 12, -1) x + 0 >= 0,
             An inequality (30, -31, 10, -1) x + 0 >= 0,
             An inequality (10, -17, 8, -1) x + 0 >= 0,
             An inequality (-154, 71, -14, 1) x + 120 >= 0,
             An inequality (-78, 49, -12, 1) x + 40 >= 0)
            sage: next(it)
            0
            sage: it.facet_repr()
            (An inequality (-50, 35, -10, 1) x + 24 >= 0,
             An inequality (-12, 19, -8, 1) x + 0 >= 0,
             An inequality (-20, 29, -10, 1) x + 0 >= 0,
             An inequality (60, -47, 12, -1) x + 0 >= 0,
             An inequality (-154, 71, -14, 1) x + 120 >= 0,
             An inequality (-78, 49, -12, 1) x + 40 >= 0)
            sage: next(it)
            0
            sage: it.facet_repr(False)
            (0, 1, 2, 4, 5, 7)
            sage: next(it)
            0
            sage: it.facet_repr(False)
            (0, 1, 5, 6, 7, 8)
            sage: next(it)
            0
            sage: it.facet_repr(False)
            (0, 1, 2, 3, 6, 8)
            sage: next(it)
            0
            sage: next(it)
            1
            sage: next(it)
            1
            sage: it.facet_repr(False)
            (4, 5, 7)
            sage: it.facet_repr()
            (An inequality (60, -47, 12, -1) x + 0 >= 0,
             An inequality (30, -31, 10, -1) x + 0 >= 0,
             An inequality (-154, 71, -14, 1) x + 120 >= 0)
        """
        cdef size_t length
        if unlikely(self.face is NULL):
            raise ValueError("iterator not set to a face yet")
        if not self.dual:
            length = self.coatom_repr()
            if names and self._H:
                return tuple(self._H[self.coatom_repr_face[i]]
                             for i in range(length)) + self._equalities
            else:
                return tuple(Integer(self.coatom_repr_face[i])
                             for i in range(length))
        else:
            length = self.atom_repr()
            if names and self._H:
                return tuple(self._H[self.atom_repr_face[i]]
                             for i in range(length)) + self._equalities
            else:
                return tuple(Integer(self.atom_repr_face[i])
                             for i in range(length))

    def length_facet_repr(self):
        r"""
        Returns the length of the facet_repr.

        Might be faster then `len(self.facet_repr())`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: for i in it: (i, it.length_facet_repr())
            (2, 1)
            (2, 1)
            (2, 1)
            (2, 1)
            (2, 1)
            (2, 1)
            (1, 2)
            (1, 2)
            (1, 2)
            (1, 2)
            (0, 3)
            (0, 3)
            (0, 3)
            (0, 3)
            (1, 2)
            (1, 2)
            (1, 2)
            (0, 3)
            (0, 3)
            (1, 2)
            (1, 2)
            (0, 3)
            (1, 2)
            (1, 2)
            (0, 3)
            (1, 2)
        """
        if unlikely(self.face is NULL):
            raise ValueError("iterator not set to a face yet")
        if not self.dual:
            return Integer(self.coatom_repr())
        else:
            return Integer(self.length_atom_repr())

    def ignore_subfaces(self):
        r"""
        Will not visit any faces of the current face.

        Only possible when not in dual mode.

        EXAMPLES::

            sage: P = polytopes.Gosset_3_21()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(dual=False)
            sage: nr_non_simplex_faces = 1
            sage: for d in it:
            ....:     if it.length_vertex_repr() > d + 1:
            ....:         nr_non_simplex_faces += 1
            ....:     else:
            ....:         it.ignore_subfaces()
            ....:
            sage: nr_non_simplex_faces
            127
        """
        if unlikely(self.dual):
            raise ValueError("only possible when not in dual mode")
        if unlikely(self.face is NULL):
            raise ValueError("iterator not set to a face yet")
        self.forbidden[self.nr_forbidden[self.current_dimension]] = self.face
        self.nr_forbidden[self.current_dimension] += 1

    def ignore_supfaces(self):
        r"""
        Will not visit any faces of the current face.

        Only possible when not in dual mode.

        EXAMPLES::

            sage: P = polytopes.Gosset_3_21()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(dual=True)
            sage: nr_faces_with_non_simplex_quotient = 1
            sage: for d in it:
            ....:     if it.length_facet_repr() > C.dimension() - d + 1:
            ....:         nr_faces_with_non_simplex_quotient += 1
            ....:     else:
            ....:         it.ignore_supfaces()
            ....:
            sage: nr_faces_with_non_simplex_quotient
            4845
        """
        if unlikely(not self.dual):
            raise ValueError("only possible when in dual mode")
        if unlikely(self.face is NULL):
            raise ValueError("iterator not set to a face yet")
        self.forbidden[self.nr_forbidden[self.current_dimension]] = self.face
        self.nr_forbidden[self.current_dimension] += 1

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
        cdef size_t nr_newfaces
        cdef size_t nr_forbidden
        cdef uint64_t **faces
        cdef size_t newfacescounter
        if unlikely(self.current_dimension == self.dimension):
            # the function is not supposed to be called in this case
            # just to prevent it from crashing
            raise StopIteration()

        nr_newfaces = self.nr_newfaces[self.current_dimension]
        nr_forbidden = self.nr_forbidden[self.current_dimension]
        faces = self.newfaces[self.current_dimension]

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

        if self.current_dimension <= self.lowest_dimension:
            # we will not yield the empty face
            # we will not yield below what is wanted
            self.current_dimension += 1
            return 0

        if nr_newfaces <= 1:
            # there will be no more faces from intersections
            self.current_dimension += 1
            return 0

        self.nr_newfaces[self.current_dimension] -= 1
        nr_newfaces -= 1
        if not self.first_time[self.current_dimension]:
            # if there exists faces[i+1], we have visited all its faces already
            # hence we should not visit any of them again
            self.forbidden[nr_forbidden] = faces[nr_newfaces + 1]
            self.nr_forbidden[self.current_dimension] += 1
            nr_forbidden = self.nr_forbidden[self.current_dimension]

        else:
            # we will visit all the faces of faces[nr_faces]
            # once we have done so, we want to add this face to forbidden
            self.first_time[self.current_dimension] = 0

        # get the facets contained in faces[i] but not in any of the forbidden
        sig_on()
        newfacescounter = get_next_level(
            faces, nr_newfaces + 1,
            self.maybe_newfaces[self.current_dimension-1],
            self.newfaces[self.current_dimension-1],
            self.forbidden, nr_forbidden, self.face_length)
        sig_off()
        if newfacescounter:
            # if there are new faces contained in faces[i+1],
            # we will set up the variables to correctly visit them on the next
            # call of `next_face_loop`
            self.current_dimension -= 1
            self.first_time[self.current_dimension] = 1
            self.nr_newfaces[self.current_dimension] = newfacescounter
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

    cdef size_t length_atom_repr(self) except -1:
        r"""
        Calculates the number of atoms in the current face by counting the
        number of set bits.
        """
        if self.face:
            return CountFaceBits(self.face, self.face_length)

        # the face was not initialized properly
        raise LookupError("`FaceIterator` does not point to a face")

    cdef size_t coatom_repr(self) except -1:
        r"""
        Writes coatom_repr of the current face in `self.coatom_repr_face`.
        Returns its length.
        """
        cdef size_t nr_coatoms = self.coatoms.nr_faces
        cdef uint64_t **coatoms = self.coatoms.data
        cdef size_t face_length = self.face_length
        return facet_repr_from_bitrep(self.face, coatoms, self.coatom_repr_face,
                                      nr_coatoms, face_length)

    cdef size_t atom_repr(self) except -1:
        r"""
        Writes atom_repr of the current face in `self.atom_repr_face`.
        Returns its length.
        """
        cdef size_t face_length = self.face_length
        return vertex_repr_from_bitrep(self.face, self.atom_repr_face, face_length)

cdef class ListOfAllFaces(SageObject):
    r"""
    A class to store all faces of :class:`CombinatorialPolyhedron`.

    Once all faces are added, they can be sorted and then the incidences
    for the `face_lattice` can be generated.

    INPUT:

    - `CombinatorialPolyhedron`.

    .. SEEALSO::

        :meth:`CombinatorialPolyhedron._record_all_faces`,
        :meth:`CombinatorialPolyhedron._record_all_faces_helper`,
        :meth:`CombinatorialPolyhedron.face_lattice`,
        :meth:`CombinatorialPolyhedron._calculate_face_lattice_incidences`.

    EXAMPLES::

        sage: P = polytopes.Birkhoff_polytope(3)
        sage: C = CombinatorialPolyhedron(P)
        sage: C.face_lattice() # indirect doctests
        Finite lattice containing 50 elements

    ALGORITHM:

    The faces are recorded with :class:`FaceIterator` in bit-representation.
    Once created, all level-sets but the coatoms are sorted with merge sort.
    Non-trivial incidences of elements whos rank differs by 1 are determinded
    by intersecting with all coatoms. Then each intersection is looked up in
    the sorted level sets.
    """
    def __init__(self, CombinatorialPolyhedron C):
        r"""
        Initializes `ListOfAllFaces`.

        See :class:`ListOfAllFaces`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C._record_all_faces() # indirect doctests
            sage: C.face_lattice()
            Finite lattice containing 28 elements
        """
        self._mem = MemoryAllocator()
        self.dimension = C.dimension()
        self.dual = 0
        if C.bitrep_facets.nr_faces > C.bitrep_vertices.nr_faces:
            self.dual = 1
        if C._unbounded:
            self.dual = 0
        cdef FaceIterator face_iter = C._face_iter(self.dual)
        self.face_length = face_iter.face_length
        if self.dimension == 0:
            # in the case of the 0-dimensional Polyhedron, we have to fix
            # atoms and coatoms, so far this didn't matter, because we only
            # cared about proper faces
            self.atoms = get_vertices_bitrep_from_facets_tuple(((),), 1)
            self.coatoms = get_facets_bitrep_from_facets_tuple(((),), 1)
            self.face_length = self.coatoms.face_length
        else:
            self.atoms = face_iter.atoms
            self.coatoms = face_iter.coatoms
        cdef size_t nr_atoms = self.atoms.nr_faces
        self.atom_repr_face = <size_t *> \
            self._mem.allocarray(self.coatoms.nr_vertices, sizeof(size_t))
        self.coatom_repr_face = <size_t *> \
            self._mem.allocarray(self.coatoms.nr_faces, sizeof(size_t))
        self._V = C._V
        self._H = C._H
        self._equalities = C._equalities

        # copy f_vector for later use
        f_vector = C.f_vector()
        self.f_vector = <size_t *> self._mem.allocarray(self.dimension + 2,
                                                        sizeof(size_t))
        if self.dual:
            for i in range(-1, self.dimension + 1):
                self.f_vector[i+1] = f_vector[-i-2]
        else:
            for i in range(-1, self.dimension + 1):
                self.f_vector[i+1] = f_vector[i+1]

        # face_counter keeps track, if all faces have been added already
        self.face_counter = \
            <size_t *> self._mem.calloc(self.dimension + 2, sizeof(size_t))
        self.face_counter[0] = 1
        self.face_counter[self.dimension + 1] = 1
        if self.dimension > -1:
            self.face_counter[self.dimension] = self.f_vector[self.dimension]

        # initializing `faces`
        cdef ListOfFaces some_list  # assuming a type
        cdef ListOfFaces coatoms_mem
        self.faces_mem = \
            tuple(ListOfFaces(self.f_vector[i+1], nr_atoms)
                  for i in range(-1, self.dimension-1))
        if self.dimension > -1:
            self.faces_mem += (self.coatoms,)
        self.faces_mem += (ListOfFaces(1, nr_atoms),)
        # setting up a pointer for direct access to the data
        self.faces = <uint64_t ***> \
            self._mem.allocarray(self.dimension + 2, sizeof(uint64_t **))
        for i in range(self.dimension + 2):
            some_list = self.faces_mem[i]
            self.faces[i] = some_list.data

        if self.dimension != 0:
            # in case of dimension == 0, we would overwrite the coatoms
            # initialize the empty face
            char_from_tuple((), self.faces[0][0], self.face_length)
        # intialize the full polyhedron
        char_from_tuple(tuple(j for j in range(nr_atoms)),
                        self.faces[self.dimension + 1][0],
                        self.face_length)

        self.incidence_is_initialized = 0
        cdef ListOfFaces incidence_face_mem = ListOfFaces(1, nr_atoms)
        self.incidence_face = incidence_face_mem.data[0]
        self.faces_mem += (incidence_face_mem,)  # needs to be stored somewhere

        cdef int d
        # adding all faces from the iterator
        if face_iter.current_dimension != self.dimension:
            # in case there are proper faces
            d = face_iter.next_face()
            while (d == self.dimension - 1):
                d = face_iter.next_face()
            while (d < self.dimension):
                self._add_face(d, face_iter.face)
                d = face_iter.next_face()
        self._sort()

    cdef int _add_face(self, int face_dim, uint64_t *face) except -1:
        r"""
        Adds a face to `ListOfAllFaces`.

        This method is used at initialization only.
        """
        cdef size_t counter = self.face_counter[face_dim + 1]
        cdef size_t max_number = self.f_vector[face_dim + 1]
        if unlikely(counter >= max_number):
            raise IOError("trying to add too many faces to `ListOfAllFaces`")
        memcpy(self.faces[face_dim + 1][counter], face,
               self.face_length*8)
        self.face_counter[face_dim + 1] += 1

    cdef int _sort(self) except -1:
        r"""
        Sorts the list faces in vertex-representation (except for facets).

        This way one can fastly find a certain face in the list later.

        This method is used an initialization only.
        """
        cdef int dim = self.dimension
        cdef int i
        for i in range(dim + 2):
            if unlikely(self.f_vector[i] != self.face_counter[i]):
                print (tuple((i, self.f_vector[i], self.face_counter[i],
                              i+1, self.f_vector[i+1], self.face_counter[i+1])))
                raise ValueError("`ListOfAllFaces` does not contain all faces")
        for i in range(0, dim):
            self._sort_one_list(self.faces[i], self.f_vector[i])

    cdef int _sort_one_list(self, uint64_t **faces, size_t nr_faces) except -1:
        r"""
        Sorts `faces` of length `nr_faces`. Each face in `faces`
        is supposed to be in vertex-repr, i.e. of length `face_length_vertex`.

        See :meth:`sort`.
        """
        cdef MemoryAllocator mem = MemoryAllocator()
        cdef uint64_t **extra_mem = \
            <uint64_t **> mem.malloc(nr_faces*sizeof(uint64_t *))
        self._sort_one_list_loop(faces, faces, extra_mem, nr_faces)

    cdef int _sort_one_list_loop(
            self, uint64_t **inp, uint64_t **output1,
            uint64_t **output2, size_t nr_faces) except -1:
        r"""
        This is mergesort.

        Sorts `inp` and returns it in `output1`.

        BEWARE: Input is the same as output1 or output2

        See :meth:`sort`.
        """
        if unlikely(nr_faces == 0):
            # prevent it from crashing
            # in this case there is nothing to do anyway
            return 0
        cdef size_t middle = nr_faces//2
        cdef size_t other = nr_faces - middle
        cdef size_t i = 0
        cdef size_t j = middle
        cdef size_t counter = 0
        if nr_faces == 1:
            output1[0] = inp[0]
            return 0

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
            ....: cimport CombinatorialPolyhedron, FaceIterator, ListOfAllFaces
            ....:
            ....: def find_face_from_iterator(it, C1):
            ....:     cdef FaceIterator face_iter = it
            ....:     cdef CombinatorialPolyhedron C = C1
            ....:     C._record_all_faces()
            ....:     cdef ListOfAllFaces all_faces = C._all_faces
            ....:     if not (all_faces.dual == it.dual):
            ....:         raise ValueError("iterator and allfaces not in same mode")
            ....:     return all_faces.find_face(face_iter.current_dimension, face_iter.face)
            ....: ''')
            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: next(it)
            2
            sage: find_face_from_iterator(it, C)
            Traceback (most recent call last):
            ..
            ValueError: cannot find a facet, as those are not sorted
            sage: it.set_record_dimension(1)
            sage: set(find_face_from_iterator(it, C) for _ in it)
            {0,
             1,
             2,
             3,
             4,
             5,
             6,
             7,
             8,
             9,
             10,
             11,
             12,
             13,
             14,
             15,
             16,
             17,
             18,
             19,
             20,
             21,
             22,
             23,
             24,
             25,
             26,
             27,
             28,
             29,
             30,
             31,
             32,
             33,
             34,
             35}
        """
        if unlikely(dimension == self.dimension -1):
            raise ValueError("cannot find a facet, as those are not sorted")
            # of course one can easily add a function to search for a facet as
            # well, but there seems to be no need for that
        if unlikely(dimension < -1 or dimension > self.dimension):
            raise IndexError("dimension out of range")
        cdef size_t start = 0
        cdef size_t middle
        cdef nr_faces = self.f_vector[dimension + 1]
        cdef uint64_t **list_faces = self.faces[dimension + 1]
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
        return memcmp(one, two, self.face_length*8) < 0

    cdef inline int is_equal(self, int dimension, size_t index,
                             uint64_t *face) except -1:
        r"""
        Checks wether `face` is in the list with dimension `dimension`
        and index `index`.
        """
        if unlikely(dimension < -1 or dimension > self.dimension
                    or index >= self.f_vector[dimension + 1]):
            raise IndexError()
        cdef uint64_t * face2 = self.faces[dimension+1][index]
        cdef size_t i
        return (0 == memcmp(face, face2, self.face_length*8))

    def vertex_repr(self, dimension, index, names=True):
        r"""
        Returns the vertex-representation of the face of dimension `dimension`
        and index `index`.

        The vertex-representation consists of
        the vertices/rays/lines that face contains.

        INPUT:

        - `dimension` -- dimension of the face.
        - `index` -- index of the face.
        - `names` -- if `True` returns the names of the vertices/rays/lines
          as given on initialization of the `CombinatorialPolyhedron`.

        EXAMPLES::

            sage: cython('''
            ....: from libc.stdint cimport uint64_t
            ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
            ....: cimport CombinatorialPolyhedron, FaceIterator, ListOfAllFaces
            ....:
            ....: def vertex_repr_via_all_faces_from_iterator(it, C1, names):
            ....:     cdef FaceIterator face_iter = it
            ....:     cdef CombinatorialPolyhedron C = C1
            ....:     cdef int dimension = face_iter.current_dimension
            ....:     C._record_all_faces()
            ....:     cdef ListOfAllFaces all_faces = C._all_faces
            ....:     if not (all_faces.dual == it.dual):
            ....:         raise ValueError("iterator and allfaces not in same mode")
            ....:     index = all_faces.find_face(dimension, face_iter.face)
            ....:     return all_faces.vertex_repr(dimension, index, names)
            ....: ''')
            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: it.set_record_dimension(1)
            sage: next(it)
            1
            sage: vertex_repr_via_all_faces_from_iterator(it, C, True)
            (A vertex at (3, 1, 4, 2), A vertex at (3, 2, 4, 1))
            sage: it.vertex_repr()
            (A vertex at (3, 1, 4, 2), A vertex at (3, 2, 4, 1))
            sage: all(vertex_repr_via_all_faces_from_iterator(it, C, True) ==
            ....:     it.vertex_repr() for _ in it)
            True

            sage: P = polytopes.twenty_four_cell()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: d = next(it)
            sage: while (d == 3): d = next(it)
            sage: vertex_repr_via_all_faces_from_iterator(it, C, True)
            (A vertex at (-1/2, 1/2, -1/2, -1/2),
             A vertex at (-1/2, 1/2, 1/2, -1/2),
             A vertex at (0, 0, 0, -1))
            sage: vertex_repr_via_all_faces_from_iterator(it, C, False)
            (5, 7, 11)
            sage: all(vertex_repr_via_all_faces_from_iterator(it, C, False) ==
            ....:     it.vertex_repr(False) for _ in it)
            True
        """
        cdef size_t length
        if self.dual:
            dimension = self.dimension - 1 - dimension
            length = self.coatom_repr(dimension, index)
            if names and self._V:
                return tuple(self._V[self.coatom_repr_face[i]]
                             for i in range(length))
            else:
                return tuple(Integer(self.coatom_repr_face[i])
                             for i in range(length))
        else:
            length = self.atom_repr(dimension, index)
            if names and self._V:
                return tuple(self._V[self.atom_repr_face[i]]
                             for i in range(length))
            else:
                return tuple(Integer(self.atom_repr_face[i])
                             for i in range(length))

    def facet_repr(self, dimension, index, names=True):
        r"""
        Returns the facet-representation of the face of dimension `dimension`
        and index `index`.

        The facet-representation consists of the facets
        that contain the face and of the equalities of the Polyhedron.

        INPUT:

        - `dimension` -- dimension of the face.
        - `index` -- index of the face.
        - `names` -- if `True` returns the names of the facets/equations
          as given on initialization of the `CombinatorialPolyhedron`.

        EXAMPLES::

            sage: cython('''
            ....: from libc.stdint cimport uint64_t
            ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
            ....: cimport CombinatorialPolyhedron, FaceIterator, ListOfAllFaces
            ....:
            ....: def facet_repr_via_all_faces_from_iterator(it, C1, names):
            ....:     cdef FaceIterator face_iter = it
            ....:     cdef CombinatorialPolyhedron C = C1
            ....:     cdef int dimension = face_iter.current_dimension
            ....:     C._record_all_faces()
            ....:     cdef ListOfAllFaces all_faces = C._all_faces
            ....:     if not (all_faces.dual == it.dual):
            ....:         raise ValueError("iterator and allfaces not in same mode")
            ....:     index = all_faces.find_face(dimension, face_iter.face)
            ....:     return all_faces.facet_repr(dimension, index, names)
            ....: ''')
            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: it.set_record_dimension(1)
            sage: next(it)
            1
            sage: facet_repr_via_all_faces_from_iterator(it, C, True)
            (An inequality (0, 0, -1, 0) x + 4 >= 0,
             An inequality (0, 1, 0, 1) x - 3 >= 0,
             An equation (1, 1, 1, 1) x - 10 == 0)
            sage: all(facet_repr_via_all_faces_from_iterator(it, C, True) ==
            ....:     it.facet_repr() for _ in it)
            True

            sage: P = polytopes.twenty_four_cell()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: d = next(it)
            sage: while (d == 3): d = next(it)
            sage: facet_repr_via_all_faces_from_iterator(it, C, True)
            (An inequality (1, 0, 0, 1) x + 1 >= 0, An inequality (0, -1, 0, 1) x + 1 >= 0)
            sage: facet_repr_via_all_faces_from_iterator(it, C, False)
            (16, 23)
            sage: all(facet_repr_via_all_faces_from_iterator(it, C, False) ==
            ....:     it.facet_repr(False) for _ in it)
            True
        """
        cdef size_t length
        if not self.dual:
            length = self.coatom_repr(dimension, index)
            if names and self._H:
                return tuple(self._H[self.coatom_repr_face[i]]
                             for i in range(length)) + self._equalities
            else:
                return tuple(Integer(self.coatom_repr_face[i])
                             for i in range(length))
        else:
            dimension = self.dimension - 1 - dimension
            length = self.atom_repr(dimension, index)
            if names and self._H:
                return tuple(self._H[self.atom_repr_face[i]]
                             for i in range(length)) + self._equalities
            else:
                return tuple(Integer(self.atom_repr_face[i])
                             for i in range(length))

    cdef size_t coatom_repr(self, int dimension, size_t index) except -1:
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
        if unlikely(dimension < -1 or dimension > self.dimension):
            raise ValueError("no face of dimension %s"%dimension)
        if unlikely(index >= self.f_vector[dimension + 1]):
            raise IndexError("no %s-th face of dimension %s"%(index, dimension))
        cdef size_t nr_coatoms = self.f_vector[self.dimension]
        cdef uint64_t **coatoms = self.faces[self.dimension]
        cdef size_t face_length = self.face_length
        cdef uint64_t *face = self.faces[dimension+1][index]
        return facet_repr_from_bitrep(face, coatoms, self.coatom_repr_face,
                                      nr_coatoms, face_length)

    cdef size_t atom_repr(self, int dimension, size_t index) except -1:
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
        if unlikely(dimension < -1 or dimension > self.dimension):
            raise ValueError("no face of dimension %s"%dimension)
        if unlikely(index >= self.f_vector[dimension + 1]):
            raise IndexError("no %s-th face of dimension %s"%(index, dimension))
        cdef size_t face_length = self.face_length
        cdef uint64_t *face = self.faces[dimension+1][index]
        return vertex_repr_from_bitrep(face, self.atom_repr_face, face_length)

    cdef void incidence_init(self, int dimension_one, int dimension_two):
        r"""
        Initialized the ListOfAllFaces to give incidences between
        `dimension_one` and `dimension_two`.

        Currently only `dimension_one == dimension_two + 1` is properly
        implmented. This is the relevant case to get all incidences for the
        `hasse_diagram`/`face_lattice`.

        :meth:`next_incidence`.
        """
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
            # the empty polyhedron is incident to every face
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
        self.incidence_dim_one = dimension_one
        self.incidence_dim_two = dimension_two
        self.incidence_counter_one = 0
        self.incidence_counter_two = 0
        self.incidence_is_initialized = 1

    cdef inline int next_trivial_incidence(self, size_t *one, size_t *two):
        r"""
        Handling the case where dimension_one is dimension of Polyhedron.

        See :meth:`next_incidence`.
        """
        one[0] = 0
        two[0] = self.incidence_counter_two
        self.incidence_counter_two += 1
        return (two[0] < self.f_vector[self.incidence_dim_two + 1])

    cdef inline int next_trivial_incidence2(self, size_t *one, size_t *two):
        r"""
        Handling the case where dimension_two is -1.

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
        cdef uint64_t *dimension_one_face
        cdef uint64_t **coatoms
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
        dimension_one_face = self.faces[self.incidence_dim_one + 1]\
                                       [self.incidence_counter_one]
        coatoms = self.faces[self.dimension]
        while (self.incidence_counter_two < self.f_vector[self.dimension]):
            intersection(dimension_one_face,
                         coatoms[self.incidence_counter_two],
                         self.incidence_face,
                         self.face_length)
            location = \
                self.find_face(self.incidence_dim_two, self.incidence_face)
            is_it_equal = self.is_equal(self.incidence_dim_two,
                                        location, self.incidence_face)
            self.incidence_counter_two += 1
            if is_it_equal:
                two[0] = location
                if self.incidence_counter_two == self.f_vector[self.dimension]:
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
        Combinatorial Type of a Polyhedron of dimension 1 with 1 vertices
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
        sage: vert = P.Vrepresentation()
        sage: C = CombinatorialPolyhedron(data, vertices=vert) # wrong usage!
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 3 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C.vertices()
        (A line in the direction (0, 1),
         A line in the direction (0, 1),
         A line in the direction (0, 1))
        sage: C = CombinatorialPolyhedron(data, vertices=vert, nr_lines=1) #correct
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices
        sage: C.f_vector()
        (1, 0, 2, 1)
        sage: C.vertices()
        ()

        sage: P = Polyhedron(rays=[[1,0],[0,1]])
        sage: C = CombinatorialPolyhedron(P) # this works fine
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 1 vertices
        sage: data = P.incidence_matrix()
        sage: vert = P.Vrepresentation()
        sage: C = CombinatorialPolyhedron(data, vertices=vert) # wrong usage!
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 3 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C.vertices()
        (A vertex at (0, 0), A vertex at (0, 0), A vertex at (0, 0))
        sage: C = CombinatorialPolyhedron(data, vertices=vert, nr_lines=0) #correct
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 1 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C.vertices()
        (A vertex at (0, 0),)
    """
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
        self._edges = NULL
        self._ridges = NULL
        self._face_lattice_incidences = NULL
        self._nr_lines = 0
        self._length_edges_list = 16348
        self._all_faces = None
        self._mem_tuple = ()
        if nr_lines is None:
            self._unbounded = 0
            self._nr_lines = 0
        else:
            self._unbounded = 1 + int(nr_lines)
            self._nr_lines = int(nr_lines)
        self._equalities = ()

        if is_Polyhedron(data):
            # Input is `Polyhedron`
            vertices = data.Vrepresentation()
            facets = tuple(inequality for inequality in data.Hrepresentation())

            if not data.is_compact():
                self._unbounded = 1 + data.n_lines()
                self._nr_lines = int(data.n_lines())

            data = data.incidence_matrix()

        if is_LatticePolytope(data):
            # Input is LatticePolytope
            self._unbounded = 0
            self._nr_lines = 0
            vertices = data.vertices()
            self._length_Vrep = len(vertices)
            facets = data.facets()
            self._length_Hrep = len(facets)
            data = tuple(tuple(vert for vert in facet.vertices())
                         for facet in facets)

        if vertices:
            # the vertices have names, which the user might want later
            self._V = tuple(vertices)
            self._Vinv = {v: i for i,v in enumerate(self._V)}
        else:
            self._V = None
            self._Vinv = None

        if facets:
            # the facets have names, which the user might want later
            facets = tuple(facets)
            test = [1] * len(facets) # 0 if that facet is an equality
            for i in range(len(facets)):
                if hasattr(facets[i], "is_inequality"):
                    if not facets[i].is_inequality():
                        test[i] = 0
            self._H = \
                tuple(facets[i] for i in range(len(facets)) if test[i])
            # only keeping those that are actual inequalities
            self._equalities = \
                tuple(facets[i] for i in range(len(facets)) if not test[i])
            # the inequalities are saved here
        else:
            self._H = None

        if is_Matrix(data):
            self._length_Hrep = data.ncols()
            self._length_Vrep = data.nrows()

            # initializing the facets as BitVectors
            self.bitrep_facets = \
                get_facets_from_incidence_matrix(data)
            # initializing the vertices as BitVectors
            self.bitrep_vertices = \
                get_vertices_from_incidence_matrix(data)

            self._nr_facets = self.bitrep_facets.nr_faces

        elif isinstance(data, Integer):  # input for a trivial Polyhedron
            if data < -1:
                TypeError("any polyhedron must have dimension at least -1")
            self._nr_facets = 0
            self._dimension = data
            # initializing the facets as BitVectors
            self.bitrep_facets = \
                get_facets_bitrep_from_facets_tuple((), 0)
            # initializing the vertices as BitVectors
            self.bitrep_vertices = \
                get_vertices_bitrep_from_facets_tuple((), 0)

        else:
            if is_iterator(data):
                data = tuple(data)
            # assumes the facets are given as a list of vertices/rays/lines

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

            # initializing the facets as BitVectors
            self.bitrep_facets = \
                get_facets_bitrep_from_facets_tuple(facets, nr_vertices)
            # initializing the vertices as BitVectors
            self.bitrep_vertices = \
                get_vertices_bitrep_from_facets_tuple(facets, nr_vertices)

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
            'Combinatorial Type of a Polyhedron of dimension -1 with 0 vertices'

            sage: P = Polyhedron(vertices=[[0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a Polyhedron of dimension 0 with 1 vertices'

            sage: P = Polyhedron(lines=[[0,0,1],[0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices'

            sage: P = Polyhedron(rays=[[1,0,0],[0,1,0],[-1,0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C._repr_()
            'Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices'
        """
        return "Combinatorial Type of a Polyhedron of "\
               "dimension %s with %s vertices" \
               % (self.dimension(), self.nr_vertices())

    def __reduce__(self):
        r"""
        Override __reduce__ to correctly pickle/unpickle.

        TESTS::

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_iter()
            sage: it1 = C1.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(),it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_iter()
            sage: it1 = C1.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(),it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_iter()
            sage: it1 = C1.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(),it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0],
            ....:                      [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_iter()
            sage: it1 = C1.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(),it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True
        """
        nr_lines = None
        if self._unbounded:
            nr_lines = Integer(self._nr_lines)
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
            if self._nr_facets == 0:
                self._dimension = self._length_Vrep - 1
            elif self._unbounded or self._nr_facets <= self._length_Vrep:
                self._dimension = calculate_dimension(self.bitrep_facets)
            else:
                # if the Polyhedron is bounded an has many facets,
                # calculating the dimenion of the dual will be faster
                self._dimension = calculate_dimension(self.bitrep_vertices)
        return Integer(self._dimension)

    def nr_vertices(self):
        r"""
        Returns the number of vertices.

        Is equivalent to `len(self.vertices())`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_vertices()
            8

            sage: P = polytopes.cyclic_polytope(4,20)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_vertices()
            20

            sage: P = Polyhedron(lines=[[0,1]], vertices=[[1,0], [-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_vertices()
            0

            sage: P = Polyhedron(rays=[[1,0,0], [0,1,0]], lines=[[0,0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_vertices()
            0

            sage: C = CombinatorialPolyhedron(4)
            sage: C.f_vector()
            (1, 0, 0, 0, 0, 1)
            sage: C.nr_vertices()
            0

            sage: C = CombinatorialPolyhedron(0)
            sage: C.f_vector()
            (1, 1)
            sage: C.nr_vertices()
            1
        """
        if self.dimension() == 0:
            # this trivial Polyhedron needs special attention
            return Integer(1)
        elif not self._unbounded:
            return Integer(self._length_Vrep)
        else:
            return len(self.vertices())

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
            (A vertex at (1, 0, 0),
             A vertex at (0, 1, 0),
             A vertex at (0, 0, 1),
             A vertex at (0, 0, -1),
             A vertex at (0, -1, 0),
             A vertex at (-1, 0, 0))

            sage: C.vertices(names=False)
            (5, 4, 3, 2, 1, 0)

            sage: points = [(1,0,0), (0,1,0), (0,0,1),
            ....:           (-1,0,0), (0,-1,0), (0,0,-1)]
            sage: L = LatticePolytope(points)
            sage: C = CombinatorialPolyhedron(L)
            sage: C.vertices()
            (M(0, 0, -1), M(0, -1, 0), M(-1, 0, 0), M(0, 0, 1), M(0, 1, 0), M(1, 0, 0))
            sage: C.vertices(names=False)
            (5, 4, 3, 2, 1, 0)

            sage: P = Polyhedron(vertices=[[0,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.vertices()
            (A vertex at (0, 0),)
        """
        if unlikely(self.dimension() == 0):
            if names and self._V:
                return (self._V[0],)
            else:
                return (Integer(0),)
        dual = False
        if not self._unbounded:
            # it is much easier to get the vertices, that the already know
            dual = True
        face_iter = self.face_iter(0, dual=dual)
        return tuple(face_iter.vertex_repr(names=names)[0] for _ in face_iter)

    def nr_facets(self):
        r"""
        Returns the number of facets.

        Is equivalent to `len(self.facets())`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_facets()
            6

            sage: P = polytopes.cyclic_polytope(4,20)
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_facets()
            170

            sage: P = Polyhedron(lines=[[0,1]], vertices=[[1,0], [-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_facets()
            2

            sage: P = Polyhedron(rays=[[1,0], [-1,0], [0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.nr_facets()
            1

            sage: C = CombinatorialPolyhedron(-1)
            sage: C.f_vector()
            (1,)
            sage: C.nr_facets()
            0

            sage: C = CombinatorialPolyhedron(0)
            sage: C.f_vector()
            (1, 1)
            sage: C.nr_facets()
            1
        """
        if unlikely(self.dimension() == 0):
            # this trivial Polyhedron needs special attention
            return Integer(1)
        return Integer(self._nr_facets)

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
        if unlikely(self.dimension() == 0):
            # special attention for this trivial case
            return ((),)
        face_iter = self.face_iter(self.dimension() - 1, dual=False)
        tup = tuple(face_iter.vertex_repr(names=names) for _ in face_iter)

        # it is important to have the facets in the exact same order as
        # on input
        # every facet knows its order by the facet representation
        face_iter = self.face_iter(self.dimension() - 1, dual=False)
        order = tuple(face_iter.facet_repr(names=False)[0] for _ in face_iter)
        dic = {}
        for i in range(len(tup)):
            dic[order[i]] = tup[i]
        return tuple(dic[i] for i in range(len(tup)))

    def edges(self, names=True):
        r"""
        Returns the edges of the CombinatorialPolyhedron,
        i.e. the rank 1 faces.

        If ``names`` is set to ``False``, then the vertices in the edges
        are given by their indices in the Vrepresentation.

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
             ('b', 'e'),
             ('a', 'e'),
             ('c', 'd'),
             ('b', 'd'),
             ('a', 'd'),
             ('b', 'c'),
             ('a', 'c'),
             ('a', 'b'))

        TESTS::

            sage: from itertools import combinations
            sage: N = combinations(range(200),199)
            sage: C = CombinatorialPolyhedron(N)
            sage: try:
            ....:     alarm(0.001)
            ....:     C.edges()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: len(C.edges())
            19900
        """
        cdef size_t **edges
        cdef size_t nr_edges
        cdef size_t len_edgelist = self._length_edges_list
        cdef size_t j
        if self._edges is NULL:
            if self._unbounded:
                self._calculate_edges(dual=False)
            elif self._length_Vrep > self._nr_facets*self._nr_facets:
                # this is a wild estimate
                # that in this case it is better not to use the dual
                self._calculate_edges(dual=False)
            else:
                self._calculate_ridges(dual=True)
        if self._edges is NULL:
            raise ValueError('could not determine edges')
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
            return f(self._edges[i // len_edgelist][2*(i % len_edgelist)])

        def vertex_two(size_t i):
            return f(self._edges[i // len_edgelist][2*(i % len_edgelist)+1])

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

    def ridges(self, add_equalities=False, names=True):
        r"""
        Returns the ridges.

        The ridges of the CombinatorialPolyhedron are the faces
        contained in exactly two facets.

        If you want to compute all faces of codimension 1,
        use :meth:`CombinatorialPolyhedron.face_iter` instead.

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
            Combinatorial Type of a Polyhedron of dimension 1 with 1 vertices
            sage: C.ridges()
            ()
            sage: it = C.face_iter(0)
            sage: for d in it: it.facet_repr()
            (An inequality (1, 0) x + 0 >= 0, An equation (0, 1) x + 0 == 0)

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
        if self._edges is NULL:
            if self._unbounded:
                self._calculate_ridges(dual=False)
            elif self._length_Vrep*self._length_Vrep < self._nr_facets:
                # this is a wild estimate
                # that in this case it is better to use the dual
                self._calculate_edges(dual=True)
            else:
                self._calculate_ridges(dual=False)
        if self._ridges is NULL:
            raise ValueError('could not determine ridges')
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
            return f(self._ridges[i // len_edgelist][2*(i % len_edgelist)])

        def facet_two(size_t i):
            return f(self._ridges[i // len_edgelist][2*(i % len_edgelist)+1])

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

        """
        if not self._f_vector:
            self._calculate_f_vector()
        if not self._f_vector:
            raise ValueError("could not determine f_vector")
        return self._f_vector

    def face_iter(self, dimension=None, dual=None):
        r"""
        Iterator over all proper faces of specified dimension.

        INPUT:

        - `dimension` -- if specified, then iterate over only this dimension.
        - `dual` -- if `True`, then iteration starting with the vertices.
          If `False`, then iteration starting with the facets.
          If `None`, this will be determined automatically to be the fastest.

        OUTPUT: A face iterator :class:`FaceIterator`.

        .. NOTE::

            :class:`FaceIterator` is more than just a plain iterator.
            By default it will iterate over the dimensions of the faces, but
            more information can be received.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(dimension=2)
            sage: next(it)
            2
            sage: it.vertex_repr()
            (A vertex at (4, 1, 5, 2, 3),
             A vertex at (4, 2, 5, 1, 3),
             A vertex at (5, 1, 4, 2, 3),
             A vertex at (5, 2, 4, 1, 3))
            sage: next(it)
            2
            sage: it.vertex_repr()
            (A vertex at (4, 1, 5, 2, 3),
             A vertex at (4, 1, 5, 3, 2),
             A vertex at (5, 1, 4, 2, 3),
             A vertex at (5, 1, 4, 3, 2))
            sage: it.facet_repr()
            (An inequality (0, 1, 0, 0, 0) x - 1 >= 0,
             An inequality (0, 1, 0, 1, 1) x - 6 >= 0,
             An equation (1, 1, 1, 1, 1) x - 15 == 0)
            sage: it.facet_repr(names=False)
            (25, 29)
            sage: next(it)
            2
            sage: it.facet_repr(names=False)
            (12, 29)
            sage: it.vertex_repr(names=False)
            (76, 77, 82, 83, 88, 89)

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],[0,2,3],[1,2,3]])
            sage: it = C.face_iter()
            sage: for _ in it: it.vertex_repr()
            (1, 2, 3)
            (0, 2, 3)
            (0, 1, 3)
            (0, 1, 2)
            (2, 3)
            (1, 3)
            (1, 2)
            (3,)
            (2,)
            (1,)
            (0, 3)
            (0, 2)
            (0,)
            (0, 1)

            sage: P = Polyhedron(rays=[[1,0],[0,1]], vertices=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter(1)
            sage: for _ in it: it.vertex_repr()
            (A vertex at (0, 1), A vertex at (1, 0))
            (A ray in the direction (1, 0), A vertex at (1, 0))
            (A ray in the direction (0, 1), A vertex at (0, 1))

        .. SEEALSO::

            :class:`FaceIterator`.
        """
        cdef FaceIterator face_iter
        if dual is None:
            if self._unbounded or self._nr_facets <= self._length_Vrep:
                dual = False
            else:
                dual = True
        face_iter = self._face_iter(int(dual))
        if dimension is not None:
            face_iter.set_record_dimension(dimension)
        return face_iter

    cdef FaceIterator _face_iter(self, int dual):
        r"""
        A method to obtain the FaceIterator as Cython object.

        See :meth:`CombinatorialPolyhedron.face_iter`
        """
        if dual == 1 and self._unbounded:
            raise ValueError("cannot iterate over dual of unbounded Polyhedron")
        return FaceIterator(self, dual)

    def face_lattice(self):
        r"""
        Generates the face-lattice.

        OUTPUT:

        - a :class:'~sage.combinat.posets.lattices.FiniteLatticePoset'

        .. NOTE::

            Use :meth:`CombinatorialPolyhedron.face_lattice_dimension` to get
            the dimension for each element in the Face Lattice.
            Use :meth:`CombinatorialPolyhedron.face_lattice_vertex_repr` to get
            the vertex representation for each element in the Face Lattice.
            Use :meth:`CombinatorialPolyhedron.face_lattice_facet_repr` to get
            the facet_repr for each element in the Face Lattice.

        EXAMPLES::

            sage: P = Polyhedron(rays=[[1,0],[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice()
            Finite lattice containing 5 elements

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
        cdef int k
        cdef int dim
        if not self._face_lattice_incidences:
            self._calculate_face_lattice_incidences()
        if self._face_lattice_incidences is NULL:
            raise TypeError("could not determine face lattice")
        incidences = self._face_lattice_incidences
        nr_incidences = self._nr_face_lattice_incidences

        dim = self._dimension

        # the incidences recorded assume we are not in the polar case
        # in the polar case, we want to somewhat reverse the order of the faces
        # the j-th face of dimension i will be the j-th face of dimension
        # `self.dimension -1 -j`

        def face_one(size_t i):
            return Integer(incidences[i // len_edgelist][2*(i % len_edgelist)])

        def face_two(size_t i):
            return Integer(incidences[i // len_edgelist][2*(i % len_edgelist)+1])

        edges = tuple((face_one(j), face_two(j))
                          for j in range(nr_incidences))

        V = tuple(range(sum(self._f_vector)))
        D = DiGraph([V, edges], format='vertices_and_edges')
        return FiniteLatticePoset(D)

    def face_lattice_dimension(self, index):
        r"""
        Returns for each element in :meth:`CombinatorialPolyhedron.face_lattice`
        its dimension.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: def f(i):
            ....:     return (i, C.face_lattice_dimension(i))
            ....:
            sage: G = F.relabel(f)
            sage: G._elements
            ((0, -1),
             (1, 0),
             (2, 0),
             (9, 1),
             (3, 0),
             (10, 1),
             (4, 0),
             (11, 1),
             (12, 1),
             (23, 2),
             (5, 0),
             (13, 1),
             (6, 0),
             (14, 1),
             (15, 1),
             (22, 2),
             (7, 0),
             (16, 1),
             (17, 1),
             (21, 2),
             (8, 0),
             (18, 1),
             (19, 1),
             (25, 2),
             (20, 1),
             (24, 2),
             (26, 2),
             (27, 3))
        """
        f_vector = self.f_vector()
        dim = self.dimension()
        return Integer(max(k for k in range(dim+2)
                       if sum(f_vector[:k]) <= index)) - 1

    def face_lattice_vertex_repr(self, index, names=True):
        r"""
        Returns for each element in :meth:`CombinatorialPolyhedron.face_lattice`
        its vertex-representation.

        If `names` is set to True, then names of the vertices are used.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: F
            Finite lattice containing 28 elements
            sage: G = F.relabel(C.face_lattice_vertex_repr)
            sage: G._elements
            ((),
             (A vertex at (-1, -1, -1),),
             (A vertex at (-1, -1, 1),),
             (A vertex at (-1, -1, -1), A vertex at (-1, -1, 1)),
             (A vertex at (-1, 1, -1),),
             (A vertex at (-1, -1, -1), A vertex at (-1, 1, -1)),
             (A vertex at (-1, 1, 1),),
             (A vertex at (-1, -1, 1), A vertex at (-1, 1, 1)),
             (A vertex at (-1, 1, -1), A vertex at (-1, 1, 1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, -1, 1),
              A vertex at (-1, 1, -1),
              A vertex at (-1, 1, 1)),
             (A vertex at (1, -1, -1),),
             (A vertex at (-1, -1, -1), A vertex at (1, -1, -1)),
             (A vertex at (1, -1, 1),),
             (A vertex at (-1, -1, 1), A vertex at (1, -1, 1)),
             (A vertex at (1, -1, -1), A vertex at (1, -1, 1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, -1, 1),
              A vertex at (1, -1, -1),
              A vertex at (1, -1, 1)),
             (A vertex at (1, 1, -1),),
             (A vertex at (-1, 1, -1), A vertex at (1, 1, -1)),
             (A vertex at (1, -1, -1), A vertex at (1, 1, -1)),
             (A vertex at (-1, -1, -1),
              A vertex at (-1, 1, -1),
              A vertex at (1, -1, -1),
              A vertex at (1, 1, -1)),
             (A vertex at (1, 1, 1),),
             (A vertex at (-1, 1, 1), A vertex at (1, 1, 1)),
             (A vertex at (1, -1, 1), A vertex at (1, 1, 1)),
             (A vertex at (-1, -1, 1),
              A vertex at (-1, 1, 1),
              A vertex at (1, -1, 1),
              A vertex at (1, 1, 1)),
             (A vertex at (1, 1, -1), A vertex at (1, 1, 1)),
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
              A vertex at (-1, 1, 1),
              A vertex at (1, -1, -1),
              A vertex at (1, -1, 1),
              A vertex at (1, 1, -1),
              A vertex at (1, 1, 1)))

            sage: P = Polyhedron(rays=[[0,1], [1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: G = F.relabel(C.face_lattice_vertex_repr)
            sage: G._elements
            ((),
             (A vertex at (0, 0),),
             (A vertex at (0, 0), A ray in the direction (0, 1)),
             (A vertex at (0, 0), A ray in the direction (1, 0)),
             (A vertex at (0, 0),
              A ray in the direction (0, 1),
              A ray in the direction (1, 0)))

            sage: def f(i): return C.face_lattice_vertex_repr(i, False)
            sage: G = F.relabel(f)
            sage: G._elements
            ((), (0,), (0, 1), (0, 2), (0, 1, 2))
        """
        self._record_all_faces()
        dim = self.face_lattice_dimension(index)
        newindex = index - sum(self._f_vector[:dim + 1])
        return self._all_faces.vertex_repr(dim, newindex, names=names)

    def face_lattice_facet_repr(self, index, names=True):
        r"""
        Returns for each element in :meth:`CombinatorialPolyhedron.face_lattice`
        its facet-representation.

        If `names` is set to True, then names of the vertices are used.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: F
            Finite lattice containing 28 elements
            sage: G = F.relabel(C.face_lattice_facet_repr)
            sage: G._elements
            ((An inequality (0, 0, -1) x + 1 >= 0,
              An inequality (0, -1, 0) x + 1 >= 0,
              An inequality (-1, 0, 0) x + 1 >= 0,
              An inequality (1, 0, 0) x + 1 >= 0,
              An inequality (0, 0, 1) x + 1 >= 0,
              An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0,
              An inequality (0, -1, 0) x + 1 >= 0,
              An inequality (1, 0, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0,
              An inequality (1, 0, 0) x + 1 >= 0,
              An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0, An inequality (1, 0, 0) x + 1 >= 0),
             (An inequality (1, 0, 0) x + 1 >= 0,
              An inequality (0, 0, 1) x + 1 >= 0,
              An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (1, 0, 0) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0,
              An inequality (0, -1, 0) x + 1 >= 0,
              An inequality (-1, 0, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0, An inequality (0, -1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0,
              An inequality (-1, 0, 0) x + 1 >= 0,
              An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0, An inequality (-1, 0, 0) x + 1 >= 0),
             (An inequality (0, 0, -1) x + 1 >= 0,),
             (An inequality (0, -1, 0) x + 1 >= 0,
              An inequality (1, 0, 0) x + 1 >= 0,
              An inequality (0, 0, 1) x + 1 >= 0),
             (An inequality (1, 0, 0) x + 1 >= 0, An inequality (0, 0, 1) x + 1 >= 0),
             (An inequality (0, -1, 0) x + 1 >= 0, An inequality (1, 0, 0) x + 1 >= 0),
             (An inequality (1, 0, 0) x + 1 >= 0,),
             (An inequality (0, -1, 0) x + 1 >= 0,
              An inequality (-1, 0, 0) x + 1 >= 0,
              An inequality (0, 0, 1) x + 1 >= 0),
             (An inequality (0, -1, 0) x + 1 >= 0, An inequality (0, 0, 1) x + 1 >= 0),
             (An inequality (0, -1, 0) x + 1 >= 0, An inequality (-1, 0, 0) x + 1 >= 0),
             (An inequality (0, -1, 0) x + 1 >= 0,),
             (An inequality (-1, 0, 0) x + 1 >= 0,
              An inequality (0, 0, 1) x + 1 >= 0,
              An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (-1, 0, 0) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 0, 1) x + 1 >= 0, An inequality (0, 1, 0) x + 1 >= 0),
             (An inequality (0, 1, 0) x + 1 >= 0,),
             (An inequality (-1, 0, 0) x + 1 >= 0, An inequality (0, 0, 1) x + 1 >= 0),
             (An inequality (0, 0, 1) x + 1 >= 0,),
             (An inequality (-1, 0, 0) x + 1 >= 0,),
             ())

            sage: P = Polyhedron(rays=[[0,1], [1,0]], vertices=[[0,1], [1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: G = F.relabel(C.face_lattice_facet_repr)
            sage: G._elements
            ((An inequality (1, 0) x + 0 >= 0,
              An inequality (0, 1) x + 0 >= 0,
              An inequality (1, 1) x - 1 >= 0),
             (An inequality (1, 0) x + 0 >= 0, An inequality (1, 1) x - 1 >= 0),
             (An inequality (1, 0) x + 0 >= 0,),
             (An inequality (0, 1) x + 0 >= 0, An inequality (1, 1) x - 1 >= 0),
             (An inequality (0, 1) x + 0 >= 0,),
             (An inequality (1, 1) x - 1 >= 0,),
             ())

            sage: def f(i): return C.face_lattice_facet_repr(i, False)
            sage: G = F.relabel(f)
            sage: G._elements
            ((0, 1, 2), (0, 2), (1, 2), (1,), (2,), (0,), ())
        """
        self._record_all_faces()
        dim = self.face_lattice_dimension(index)
        newindex = index - sum(self._f_vector[:dim + 1])
        return self._all_faces.facet_repr(dim, newindex, names=names)

    cdef int _calculate_f_vector(self) except -1:
        r"""
        Calculates the f_vector of the `CombinatorialPolyhedron`.

        In the polar case, calculates the f_vector of the polar.

        See :meth:`f_vector`.
        """
        if self._f_vector:
            return 0  # there is no need to recalculate the f_vector
        cdef FaceIterator face_iter
        cdef int dual
        if self._unbounded or self._nr_facets <= self._length_Vrep:
            dual = 0
            face_iter = self._face_iter(0)
        else:
            dual = 1
            face_iter = self._face_iter(1)
        cdef int dim = self.dimension()
        cdef int d
        cdef size_t *f_vector
        cdef MemoryAllocator mem = MemoryAllocator()
        f_vector = <size_t *> mem.calloc((dim + 2), sizeof(size_t))
        f_vector[0] = 1
        f_vector[dim + 1] = 1
        if self._nr_facets > 0 and dim > 0:
            d = face_iter.next_face()
            while (d < dim):
                f_vector[d+1] += 1
                d = face_iter.next_face()

        if dual:
            self._f_vector = \
                tuple(Integer(f_vector[dim+1-i]) for i in range(dim+2))
        else:
            self._f_vector = tuple(Integer(f_vector[i]) for i in range(dim+2))

    cdef int _calculate_edges(self, dual) except -1:
        r"""
        Calculates the edges of the `CombinatorialPolyhedron`.

        If `dual` is `True`, calculates the edges of the polar.
        This will also calculate the `f_vector`, if its not already calculated.

        See :meth:`edges` and :meth:`ridges`.
        """
        if (self._edges is not NULL and not dual) \
                or (self._ridges is not NULL and dual):
            return 0  # there is no need to recalculate
        cdef FaceIterator face_iter = self._face_iter(dual)
        cdef size_t len_edgelist = self._length_edges_list
        cdef int dim = self.dimension()
        cdef int d  # dimension of the current face of FaceIterator
        cdef MemoryAllocator mem = MemoryAllocator()
        cdef size_t *f_vector
        cdef int is_f_vector  # True if f_vector was calculated previously
        cdef size_t **edges
        cdef size_t counter  # the number of edges so far
        # for each edge we determine its location in `edges` by edges[one][two]
        cdef size_t one
        cdef size_t two
        cdef size_t current_length # dynamically enlarge **edges

        if self._f_vector:
            is_f_vector = 1
        else:
            # in this case we will calculate the f_vector while we're at it
            is_f_vector = 0
        counter = 0

        edges = <size_t**> mem.malloc(sizeof(size_t*))
        current_length = 1

        if dim == 1:
            # in this case there is an edge, but its not a proper face
            edges[0] = <size_t *> mem.malloc(2 * sizeof(size_t))
            edges[0][0] = 0
            edges[0][1] = 1
            counter = 1
            if dual:
                sig_block()
                self._nr_ridges = counter
                self._ridges = edges
                self._mem_tuple += (mem,)
                sig_unblock()
            else:
                sig_block()
                self._nr_edges = counter
                self._edges = edges
                self._mem_tuple += (mem,)
                sig_unblock()
            return 0

        if is_f_vector:
            if not dual:
                face_iter.set_record_dimension(1)
            else:
                face_iter.set_record_dimension(dim - 2)
            if self._nr_facets > 0 and dim > 0:
                # if not, there won't even be any edges
                while (face_iter.next_face() == 1):
                    one = counter // len_edgelist
                    two = counter % len_edgelist
                    if unlikely(two == 0):
                        if unlikely(one + 1 > current_length):
                            # enlarge **edges
                            current_length *= 2
                            edges = <size_t **> \
                                mem.realloc(edges, current_length*sizeof(size_t*))
                        edges[one] = <size_t *> \
                            mem.malloc(2 * len_edgelist * sizeof(size_t))

                    # setting up face_iter.atom_repr_face
                    face_iter.atom_repr()
                    # now copying the information to edges
                    edges[one][2*two] = face_iter.atom_repr_face[0]
                    edges[one][2*two + 1] = face_iter.atom_repr_face[1]
                    counter += 1

            # success, we can copy the data to the CombinatorialPolyhedron
            if dual:
                sig_block()
                self._nr_ridges = counter
                self._ridges = edges
                self._mem_tuple += (mem,)
                sig_unblock()
            else:
                sig_block()
                self._nr_edges = counter
                self._edges = edges
                self._mem_tuple += (mem,)
                sig_unblock()
        else:
            # while doing the edges one might as well do the f_vector
            f_vector = <size_t *> mem.malloc((dim + 2) * sizeof(size_t))
            for i in range(dim + 2):
                f_vector[i] = 0
            f_vector[0] = 1
            f_vector[dim + 1] = 1

            counter = 0
            if self._nr_facets > 0 and dim > 0:
                # if not, there won't even be any edges
                d = face_iter.next_face()
                while (d < dim):
                    f_vector[d+1] += 1
                    if d == 1:
                        one = counter // len_edgelist
                        two = counter % len_edgelist
                        if unlikely(two == 0):
                            if unlikely(one + 1 > current_length):
                                # enlarge **edges
                                current_length *= 2
                                edges = <size_t **> \
                                    mem.realloc(edges, current_length*sizeof(size_t*))
                            edges[one] = <size_t *> \
                                mem.malloc(2 * len_edgelist * sizeof(size_t))

                        # setting up face_iter.atom_repr_face
                        face_iter.atom_repr()
                        # now copying the information to edges
                        edges[one][2*two] = face_iter.atom_repr_face[0]
                        edges[one][2*two + 1] = face_iter.atom_repr_face[1]
                        counter += 1

                    d = face_iter.next_face()

            # success, we can copy the data to the CombinatorialPolyhedron
            if dual:
                sig_block()
                self._f_vector = \
                    tuple(Integer(f_vector[dim+1-i]) for i in range(dim+2))
                self._nr_ridges = counter
                self._ridges = edges
                self._mem_tuple += (mem,)
                sig_unblock()
            else:
                sig_block()
                self._f_vector = \
                    tuple(Integer(f_vector[i]) for i in range(dim+2))
                self._nr_edges = counter
                self._edges = edges
                self._mem_tuple += (mem,)
                sig_unblock()

    cdef int _calculate_ridges(self, dual) except -1:
        r"""
        Calculates the ridges of the `CombinatorialPolyhedron`.

        If `dual` is `True`, calculates the ridges of the polar.

        See :meth:`edges` and :meth:`ridges`.
        """
        if (self._edges is not NULL and dual) \
                or (self._ridges is not NULL and not dual):
            return 0  # there is no need to recalculate
        cdef FaceIterator face_iter = self._face_iter(dual)
        cdef size_t len_edgelist = self._length_edges_list
        cdef int dim = self.dimension()
        cdef int d  # dimension of the current face of FaceIterator
        cdef MemoryAllocator mem = MemoryAllocator()
        cdef size_t **ridges
        cdef size_t counter  # the number of riddges so far
        # for each ridge we determine its location in `ridges` by ridges[one][two]
        cdef size_t one
        cdef size_t two
        cdef size_t current_length # dynamically enlarge **ridges

        counter = 0

        ridges = <size_t**> mem.malloc(sizeof(size_t*))
        current_length = 1

        if dim == 1 and self._nr_facets > 1:
            # in this case there is a ridge, but its not a proper face
            ridges[0] = <size_t *> mem.malloc(2 * sizeof(size_t))
            ridges[0][0] = 0
            ridges[0][1] = 1
            counter = 1
            if not dual:
                sig_block()
                self._nr_ridges = counter
                self._ridges = ridges
                self._mem_tuple += (mem,)
                sig_unblock()
            else:
                sig_block()
                self._nr_edges = counter
                self._edges = ridges
                self._mem_tuple += (mem,)
                sig_unblock()
            return 0

        if dual:
            face_iter.set_record_dimension(1)
        else:
            face_iter.set_record_dimension(dim - 2)
        if self._nr_facets > 1 and dim > 0:
            # if not, there won't even be any ridges (as intersection of two distince facets)
            while (face_iter.next_face() == dim - 2):
                one = counter // len_edgelist
                two = counter % len_edgelist
                if unlikely(two == 0):
                    if unlikely(one + 1 > current_length):
                        # enlarge **ridges
                        current_length *= 2
                        ridges = <size_t **> \
                            mem.realloc(ridges, current_length*sizeof(size_t*))
                    ridges[one] = <size_t *> \
                        mem.malloc(2 * len_edgelist * sizeof(size_t))

                # setting up face_iter.coatom_repr_face
                face_iter.coatom_repr()
                # now copying the information to edges
                ridges[one][2*two] = face_iter.coatom_repr_face[0]
                ridges[one][2*two + 1] = face_iter.coatom_repr_face[1]
                counter += 1

        # success, we can copy the data to the CombinatorialPolyhedron
        if not dual:
            sig_block()
            self._nr_ridges = counter
            self._ridges = ridges
            self._mem_tuple += (mem,)
            sig_unblock()
        else:
            sig_block()
            self._nr_edges = counter
            self._edges = ridges
            self._mem_tuple += (mem,)
            sig_unblock()

    cdef int _calculate_face_lattice_incidences(self) except -1:
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
        cdef size_t already_seen
        cdef size_t already_seen_next
        cdef MemoryAllocator mem = MemoryAllocator()
        cdef size_t **incidences
        cdef size_t current_length # dynamically enlarge **incidences

        self._record_all_faces()
        all_faces = self._all_faces
        f_vector = self.f_vector()
        if all_faces is None:
            raise ValueError("could not determine a list of all faces")

        incidences = <size_t**> mem.malloc(sizeof(size_t*))
        current_length = 1

        counter = 0

        dimension_one = 0
        if dim > -1:
            while (f_vector[dimension_one + 1] == 0):
                # taking care of cases, where there might be no faces
                # of dimension 0
                dimension_one += 1
            dimension_two = -1
        while (dimension_one < dim + 1):
            already_seen = \
                sum(f_vector[j] for j in range(dimension_two + 1))
            already_seen_next = already_seen + f_vector[dimension_two + 1]
            if all_faces.dual:
                all_faces.incidence_init(dim - 1 - dimension_two, dim - 1 - dimension_one)
            else:
                all_faces.incidence_init(dimension_one, dimension_two)
            while all_faces.next_incidence(&second, &first):
                one = counter // len_edgelist
                two = counter % len_edgelist
                if unlikely(two == 0):
                    if unlikely(one + 1 > current_length):
                        # enlarge **incidences
                        current_length *= 2
                        incidences = <size_t **> \
                            mem.realloc(incidences,
                                        (current_length)*sizeof(size_t*))
                    incidences[one] = <size_t *> \
                        mem.malloc(2 * len_edgelist * sizeof(size_t))
                if all_faces.dual:
                    second += already_seen
                    first += already_seen_next
                    incidences[one][2*two] = second
                    incidences[one][2*two + 1] = first
                else:
                    second += already_seen_next
                    first += already_seen
                    incidences[one][2*two] = first
                    incidences[one][2*two + 1] = second
                counter += 1
                sig_check()
            dimension_one += 1
            dimension_two = dimension_one - 1

        self._nr_face_lattice_incidences = counter
        sig_block()
        self._mem_tuple += (mem,)
        self._face_lattice_incidences = incidences
        sig_unblock()

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
            sage: it = C.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple((C.face_lattice_vertex_repr(i),
            ....:               C.face_lattice_facet_repr(i)) for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple((C.face_lattice_vertex_repr(i),
            ....:               C.face_lattice_facet_repr(i)) for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple((C.face_lattice_vertex_repr(i),
            ....:               C.face_lattice_facet_repr(i)) for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0],
            ....:                      [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: tup = tuple((it.vertex_repr(),it.facet_repr()) for _ in it)
            sage: rg = range(1,sum(C.f_vector()) - 1)
            sage: tup2 = tuple((C.face_lattice_vertex_repr(i),
            ....:               C.face_lattice_facet_repr(i)) for i in rg)
            sage: sorted(tup) == sorted(tup2)
            True
        """
        if self._all_faces:
            return
        self._all_faces = ListOfAllFaces(self)
        if self._all_faces is None:
            raise ValueError("could not determine a list of all faces")
