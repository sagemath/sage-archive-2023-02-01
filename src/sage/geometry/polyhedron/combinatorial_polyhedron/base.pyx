# distutils: language = c++

r"""
CombinatorialPolyhedron gathers several algorithms of Polyhedra depending only
on the vertex-facet incidences.

Most importanly, this module offers a quick face iterator, which can be used
to calculate f-vector, edges, ridges and even the face lattice.

The FaceIterator would work on every atomic and coatomic lattice, where every
interval of length two has exactly 4 elements (known as the diamond property).

It also works on every lattice that is obtained by from such a lattice by
deleting all elements but the empty face contained in some of the coatoms.
Important examples are face lattices of unbounded polyhedra.

Representations in this module:

- Vertices              -- ``[vertices, rays, lines]`` of the polyhedron.
- Facets                -- facets of the polyhedron.
- Coatoms               -- the faces from which all others are constructed in
                           the algorithm. This will be facets or vertices.
                           In non-dual mode, faces are constructed as
                           intersections of the facets. In dual mode, the are
                           constructed theoretically as joins of vertices.
                           The coatoms are reprsented as incidences with the
                           atoms they contain.
- Atoms                 -- facets or vertices depending on application of algorithm.
                           Atoms are reprsented as incidences of coatoms they
                           are contained in.

- Vertex-Representation -- represents a face by a list of vertices it contains.
- Facet-Representation  -- represents a face by a list of facets it is contained in.
- Bit-Representation    -- represents incidences as ``uint64_t``-array, where
                           each Bit represents one incidences. There might
                           be trailing zeros, to fit alignment-requirements.
                           In most instances, faces are represented by the
                           Bit-representation, where each bit corresponds to
                           an atom.

AUTHOR:

- Jonathan Kliem (2019-03)
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

cdef extern from "bit_vector_operations.cc":
    # Any Bit-representation is assumed to be `chunksize`-Bit aligned.
    cdef const size_t chunksize

    cdef void intersection(uint64_t *A, uint64_t *B, uint64_t *C,
                           size_t face_length)
#    Return ``A & ~B == 0``.
#    A is not subset of B, iff there is a vertex in A, which is not in B.
#    ``face_length`` is the length of A and B in terms of uint64_t.

    cdef size_t get_next_level(
        uint64_t **faces, const size_t nr_faces, uint64_t **nextfaces,
        uint64_t **nextfaces2, uint64_t **visited_all,
        size_t nr_visited_all, size_t face_length)
#        Set ``newfaces`` to be the facets of ``faces[nr_faces -1]``
#        that are not contained in a face of ``visited_all``.

#        INPUT:

#        - ``maybe_newfaces`` -- quasi of type ``uint64_t[nr_faces -1][face_length]``,
#          needs to be ``chunksize``-Bit aligned
#        - ``newfaces`` -- quasi of type ``*uint64_t[nr_faces -1]
#        - ``visited_all`` -- quasi of type ``*uint64_t[nr_visited_all]
#        - ``face_length`` -- length of the faces

#        OUTPUT:

#        - return number of ``newfaces``
#        - set ``newfaces`` to point to the new faces

#        ALGORITHM:

#        To get all facets of ``faces[nr_faces-1]``, we would have to:
#        - Intersect the first ``nr_faces-1`` faces of ``faces`` with the last face.
#        - Add all the intersection of ``visited_all`` with the last face
#        - Out of both the inclusion-maximal ones are of codimension 1, i.e. facets.

#        As we have visited all faces of ``visited_all``, we alter the algorithm
#        to not revisit:
#        Step 1: Intersect the first ``nr_faces-1`` faces of ``faces`` with the last face.
#        Step 2: Out of thosse the inclusion-maximal ones are some of the facets.
#                At least we obtain all of those, that we have not already visited.
#                Maybe, we get some more.
#        Step 3: Only keep those that we have not already visited.
#                We obtain exactly the facets of ``faces[nr_faces-1]`` that we have
#                not visited yet.

    cdef size_t count_atoms(uint64_t *A, size_t face_length)
#        Return the number of atoms/vertices in A.
#        This is the number of set bits in A.
#        ``face_length`` is the length of A in terms of uint64_t.

    cdef size_t bit_repr_to_coatom_repr(
            uint64_t *face, uint64_t **coatoms, size_t nr_coatoms,
            size_t face_length, size_t *output)
#        Write the coatom-representation of face in output. Return length.
#        ``face_length`` is the length of ``face`` and ``coatoms[i]``
#        in terms of uint64_t.
#        ``nr_coatoms`` length of ``coatoms``.

cdef inline uint64_t vertex_to_bit_dictionary(size_t i):
    r"""
    Return an uint64_t with exactly the i-th bit set to 1.

    This "dictionary" helps storing a vector of 64 incidences as ``uint64_t``.
    """
    return (<uint64_t>1) << (64 - i - 1)

cdef int vertex_list_to_bit_repr(tuple vertex_list, uint64_t *output,
                                 size_t face_length) except -1:
    r"""
    Convert a vertex list into Bit-representation. Store it in ``output``.

    The first bit represent the entry ``0`` and is set to one, iff ``0`` is in
    ``vertex_list``. The second bit represents ``1`` and so on.

    INPUT:

    - ``vertex_list`` -- tuple of pairwise distinct positive integers in
      ``range(face_length*64)``
    - ``output`` -- array of ``uint64_t`` of length ``face_length``
    - ``face_length``

    OUTPUT:

    - ``output`` is filled

    EXAMPLES::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport vertex_list_to_bit_repr
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from sage.rings.integer import Integer
        ....:
        ....: def vertex_list_to_bit_repr_wrapper(tup):
        ....:     cdef size_t face_length = max(tup)//64 + 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = <uint64_t *> mem.allocarray(face_length, 8)
        ....:     vertex_list_to_bit_repr(tup, output, face_length)
        ....:     return tuple(Integer(output[i]) for i in range(face_length))
        ....:
        ....: def vertex_list_to_bit_repr_wrong_size(tup):
        ....:     cdef size_t face_length = 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = <uint64_t *> mem.allocarray(face_length, 8)
        ....:     vertex_list_to_bit_repr(tup, output, face_length)
        ....:     return tuple(Integer(output[i]) for i in range(face_length))
        ....: ''')  # long time

        sage: vertex_list_to_bit_repr_wrapper((62, 63))  # long time
        (3,)
        sage: vertex_list_to_bit_repr_wrapper((61, 63, 125))  # long time
        (5, 4)
        sage: vertex_list_to_bit_repr_wrong_size((62, 70))  # long time
        Traceback (most recent call last):
        ...
        IndexError: output too small to represent 70
        sage: vertex_list_to_bit_repr_wrapper((-1, 12))  # long time
        Traceback (most recent call last):
        ...
        OverflowError: can't convert negative value to size_t
        sage: vertex_list_to_bit_repr_wrapper((0, 0))  # long time
        Traceback (most recent call last):
        ...
        ValueError: entries of ``tup`` are not distinct
    """
    cdef size_t entry     # will iterate over tup
    cdef size_t position  # determines the position in output of entry
    cdef size_t value     # determines which bit will be set in output[position]

    memset(output, 0, face_length*8)  # initialize output
    if unlikely(len(vertex_list) != len(set(vertex_list))):
        raise ValueError("entries of ``tup`` are not distinct")
    for entry in vertex_list:
        value = entry % 64
        position = entry//64
        if unlikely(position >= face_length):
            raise IndexError("output too small to represent %s"%entry)
        output[position] += vertex_to_bit_dictionary(value)

cdef int incidences_to_bit_repr(tuple incidences, uint64_t *output,
                                size_t face_length) except -1:

    r"""
    Convert a tuple of incidences into Bit-representation.

    Store it in ``output``. Each entry in ``incidences`` represents a bit in
    ``output``. It is set to ``1``, iff the entry in ``incidences`` is non-zero.

    INPUT:

    - ``incidences`` -- tuple of integers representing incidences
      of length at most ``face_length*64``
    - ``output`` -- array of ``uint64_t`` of length ``face_length``
    - ``face_length``

    OUTPUT:

    - ``output`` is filled

    EXAMPLES::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport incidences_to_bit_repr
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from sage.rings.integer import Integer
        ....:
        ....: def incidences_to_bit_reprs_wrapper(tup):
        ....:     cdef size_t face_length = (len(tup)-1)//64 + 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = \
        ....:          <uint64_t *> mem.allocarray(face_length, 8)
        ....:     incidences_to_bit_repr(tup, output, face_length)
        ....:     return tuple(Integer(output[i]) for i in range(face_length))
        ....:
        ....: def incidences_to_bit_reprs_wrong_size(tup):
        ....:     cdef size_t face_length = 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = \
        ....:         <uint64_t *> mem.allocarray(face_length, 8)
        ....:     incidences_to_bit_repr(tup, output, face_length)
        ....:     return tuple(Integer(output[i]) for i in range(face_length))
        ....: ''')  # long time

        sage: incidences_to_bit_reprs_wrapper((0,) * 62 + (1,1))  # long time
        (3,)
        sage: incidences_to_bit_reprs_wrapper((0,) * 61 + (1,0,1) +  # long time
        ....:                              (0,) * 61 + (1,))
        (5, 4)
        sage: incidences_to_bit_reprs_wrapper((1,) * 64)  # long time
        (18446744073709551615,)
        sage: incidences_to_bit_reprs_wrong_size((1,) * 70)  # long time
        Traceback (most recent call last):
        ...
        IndexError: output too small to represent all incidences
    """
    cdef size_t entry       # index for the entries in tup
    cdef size_t position    # determines the position in output of entry
    cdef size_t value       # determines which bit will be set in output[position]
    cdef size_t length = len(incidences)

    memset(output, 0, face_length*8)  #initialize
    if unlikely(length > 64*face_length):
        raise IndexError("output too small to represent all incidences")
    for entry in range(length):
        if incidences[entry]:
            # vertex ``entry`` is contained in the face, so set the corresponding bit
            value = entry % 64
            position = entry//64
            output[position] += vertex_to_bit_dictionary(value)

def incidence_matrix_to_bit_repr_of_facets(matrix):
    r"""
    Initialize facets in Bit-representation as :class:`ListOfFaces`.

    INPUT:

    - ``matrix`` -- an incidence matrix as in
      :meth:`sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`

    OUTPUT:

    - :class:`ListOfFaces`

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, bit_repr_to_vertex_list
        ....: from sage.rings.integer import Integer
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def bit_repr_to_vertex_list_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.nr_vertices,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = bit_repr_to_vertex_list(
        ....:         data, output, faces.face_length)
        ....:     return tuple(Integer(output[i]) for i in range(length))
        ....: ''')

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import incidence_matrix_to_bit_repr_of_facets
        sage: P = polytopes.permutahedron(4)
        sage: facets = incidence_matrix_to_bit_repr_of_facets(P.incidence_matrix())
        sage: facets.nr_faces
        14
        sage: facets.nr_vertices
        24
        sage: for i in range(facets.nr_faces):
        ....:     print(bit_repr_to_vertex_list_wrapper(facets, i))
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

    # The incidence-matrix of a polyhedron, might contain columns with all 1's.
    # Those correspond to hyperplanes, the polyhedron lies in, not facets.
    # We delete those from the matrix. Also we transpose.
    matrix = matrix.transpose()
    rg = range(matrix.nrows())
    matrix = matrix[list(i for i in rg if not all(j for j in matrix[i]))]

    # Output will be a ``ListOfFaces`` with ``matrix.nrows()`` faces and
    # ``matrix.ncols()`` vertices.
    cdef size_t length = matrix.nrows()
    cdef ListOfFaces facets = ListOfFaces(length, matrix.ncols())
    cdef uint64_t **facets_data = facets.data

    cdef size_t i
    for i in range(length):
        # Filling each facet with its vertex-incidences, which "is" the
        # "i-th column" of the original matrix (but we have transposed).
        incidences_to_bit_repr(tuple(matrix[i]), facets_data[i], facets.face_length)
    return facets

def incidence_matrix_to_bit_repr_of_vertices(matrix):
    r"""
    Initialize vertices in Bit-representation as :class:`ListOfFaces`.

    Each vertex is represented as the facets it is contained in.
    Those are the facets of the polar polyhedron, if it exists.

    INPUT:

    - ``matrix`` -- an incidence matrix as in
      :meth:`sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`

    OUTPUT:

    - :class:`ListOfFaces`

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, bit_repr_to_vertex_list
        ....: from sage.rings.integer import Integer
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def bit_repr_to_vertex_list_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.nr_vertices,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = bit_repr_to_vertex_list(
        ....:         data, output, faces.face_length)
        ....:     return tuple(Integer(output[i]) for i in range(length))
        ....: ''')

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import incidence_matrix_to_bit_repr_of_vertices
        sage: P = polytopes.permutahedron(4)
        sage: vertices = incidence_matrix_to_bit_repr_of_vertices(P.incidence_matrix())
        sage: vertices.nr_faces
        24
        sage: vertices.nr_vertices
        14
        sage: for i in range(vertices.nr_faces):
        ....:     print(bit_repr_to_vertex_list_wrapper(vertices, i))
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

    # The incidence-matrix of a polyhedron, might contain columns with all 1's.
    # Those correspond to hyperplanes, the polyhedron lies in, not facets.
    # We delete those from the matrix. Also we transpose.
    matrix = matrix.transpose()
    rg = range(matrix.nrows())
    matrix = matrix[list(i for i in rg if not all(j for j in matrix[i]))]

    # Output will be a ``ListOfFaces`` with ``matrix.ncols()`` faces and
    # ``matrix.nrows()`` vertices.
    cdef size_t length = matrix.ncols()
    cdef ListOfFaces vertices = ListOfFaces(length, matrix.nrows())
    cdef uint64_t **vertices_data = vertices.data

    cdef size_t i
    for i in range(length):
        # Filling each facet with its vertex-incidences, which "is" the
        # "i-th row" of the original matrix (but we have transposed).
        incidences_to_bit_repr(tuple(matrix.column(i)), vertices_data[i], vertices.face_length)
    return vertices

def facets_tuple_to_bit_repr_of_facets(tuple facets_input, size_t nr_vertices):
    r"""
    Initializes facets in Bit-representation as :class:`ListOfFaces`.

    INPUT:

    - ``facets_input`` -- tuple of facets, each facet a tuple of vertices,
      vertices must be exactly ``range(nr_vertices)``
    - ``nr_vertices``

    OUTPUT:

    - :class:`ListOfFaces`

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, bit_repr_to_vertex_list
        ....: from sage.rings.integer import Integer
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def bit_repr_to_vertex_list_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.nr_vertices,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = bit_repr_to_vertex_list(
        ....:         data, output, faces.face_length)
        ....:     return tuple(Integer(output[i]) for i in range(length))
        ....: ''')

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import facets_tuple_to_bit_repr_of_facets
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: facets = facets_tuple_to_bit_repr_of_facets(bi_pyr, 6)
        sage: for i in range(8):
        ....:     print(bit_repr_to_vertex_list_wrapper(facets, i))
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
        vertex_list_to_bit_repr(facets_input[i], facets_data[i], face_length)
    return facets

def facets_tuple_to_bit_repr_of_vertices(tuple facets_input, size_t nr_vertices):
    r"""
    Initialize vertices in Bit-representation as :class:`ListOfFaces`.

    Each vertex is represented as the facets it is contained in.
    Those are the facets of the polar polyhedron, if it exists.

    INPUT:

    - ``facets_input`` -- tuple of facets, each facet a tuple of vertices,
      vertices must be exactly ``range(nr_vertices)``
    - ``nr_vertices``

    OUTPUT:

    - :class:`ListOfFaces`


    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, bit_repr_to_vertex_list
        ....: from sage.rings.integer import Integer
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def bit_repr_to_vertex_list_wrapper(list_of_faces, index):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef ListOfFaces faces = list_of_faces
        ....:     output = <size_t *> mem.allocarray(faces.nr_vertices,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data = faces.data[index]
        ....:     length = bit_repr_to_vertex_list(
        ....:         data, output, faces.face_length)
        ....:     return tuple(Integer(output[i]) for i in range(length))
        ....: ''')

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import facets_tuple_to_bit_repr_of_vertices
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: vertices = facets_tuple_to_bit_repr_of_vertices(bi_pyr, 6)
        sage: for i in range(6):
        ....:     print(bit_repr_to_vertex_list_wrapper(vertices, i))
        (0, 3, 4, 7)
        (0, 1, 4, 5)
        (1, 2, 5, 6)
        (2, 3, 6, 7)
        (0, 1, 2, 3)
        (4, 5, 6, 7)
    """
    cdef size_t nr_facets = len(facets_input)

    # Vertices in facet-representation will be a ``ListOfFaces``
    # with number of vertices faces and
    # number of facets "vertices"/atoms.
    cdef ListOfFaces vertices = ListOfFaces(nr_vertices, nr_facets)
    cdef uint64_t **vertices_data = vertices.data

    # Initializing the data of ListOfFaces.
    cdef size_t face_length = vertices.face_length
    cdef size_t i
    for i in range(nr_vertices):
        memset(vertices_data[i], 0, face_length*8)

    cdef size_t input_facet   # will iterate over indices of facets_input
    cdef size_t position      # determines the position in output of entry
    cdef size_t value         # determines which bit will be set in output[position]
    cdef size_t input_vertex  # will iterate over vertices in facet ``input_facet``

    for input_facet in range(nr_facets):
        value = input_facet % 64
        position = input_facet//64
        for input_vertex in facets_input[input_facet]:
            # Iff the input-vertex is in the input-facet,
            # then in facet-representation of the vertices
            # input-facet is a vertex of intput-vertex.
            vertices_data[input_vertex][position] += vertex_to_bit_dictionary(value)
    return vertices

cdef inline size_t bit_repr_to_vertex_list(uint64_t *face, size_t *output,
                                           size_t face_length) except -1:
    r"""
    Convert a bitrep-representation to a list of vertices. Return length of representation.

    Basically this is an inverse to :meth:`vertex_list_to_bit_repr`.
    Instead of returning a tuple, it stores the vertices in ``output``.

    INPUT:

    - ``face`` -- a Bit-representation of a face
    - ``output`` -- an array of ``size_t`` long enough to contain all vertices
      of that face (``face_length*64`` will suffice)
    - ``face_length`` -- the length of ``face``

    OUTPUT:

    - store vertices in ``output``
    - return "length" of ``output``

    EXAMPLES::

        sage: cython('''
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport ListOfFaces, bit_repr_to_vertex_list, vertex_list_to_bit_repr
        ....: from sage.rings.integer import Integer
        ....: from sage.ext.memory_allocator cimport MemoryAllocator
        ....: from libc.stdint cimport uint64_t
        ....:
        ....: def bit_repr_to_vertex_list_wrapper(tup):
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef size_t *output
        ....:     cdef length = len(tup)
        ....:     output = <size_t *> mem.allocarray(length*64,
        ....:                                        sizeof(size_t))
        ....:     cdef uint64_t * data
        ....:     data = <uint64_t *> mem.allocarray(length, 8)
        ....:     for i in range(len(tup)):
        ....:         data[i] = tup[i]
        ....:     outputlength = bit_repr_to_vertex_list(
        ....:         data, output, length)
        ....:     return tuple(Integer(output[i]) for i in range(outputlength))
        ....: ''')  # long time

        sage: bit_repr_to_vertex_list_wrapper((17, 31))  # long time
        (59, 63, 123, 124, 125, 126, 127)
        sage: bit_repr_to_vertex_list_wrapper((13,))  # long time
        (60, 61, 63)
        sage: bit_repr_to_vertex_list_wrapper((0, 61))  # long time
        (122, 123, 124, 125, 127)

    TESTS:

    Testing that :meth`bit_repr_to_vertex_list` is the
    inverse to :meth:`vertex_list_to_bit_repr`::

        sage: cython('''
        ....: from libc.stdint cimport uint64_t
        ....: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....: cimport bit_repr_to_vertex_list, vertex_list_to_bit_repr
        ....: from sage.misc.prandom import randint
        ....:
        ....: cdef uint64_t[2] face
        ....: cdef size_t length
        ....: cdef size_t[128] output
        ....:
        ....: for _ in range(10):
        ....:     st = set(randint(0,127) for i in range(40))
        ....:     tup = tuple(sorted(tuple(st)))
        ....:     vertex_list_to_bit_repr(tup, face, 2)
        ....:     length = bit_repr_to_vertex_list(face, output, 2)
        ....:     tup2 = tuple(output[i] for i in range(length))
        ....:     if not tup == tup2:
        ....:         print('``bit_repr_to_vertex_list`` does not behave',
        ....:               'as the inverse of ``vertex_list_to_bit_repr``')
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
                    # Then face[i] has the j-th bit set to 1.
                    # This corresponds to face containing vertex i*64 + j.
                    output[output_length] = i*64 + j
                    output_length += 1
                    copy -= vertex_to_bit_dictionary(j)
    return output_length

cdef class ListOfFaces(SageObject):
    r"""
    A class to store the Bit-representation of faces in.

    This class will allocate the memory for the faces.

    INPUT:

    - ``nr_faces`` -- the number of faces to be stored
    - ``nr_vertices`` -- the total number of vertices of the Polyhedron

    .. SEEALSO::

        :meth:`incidence_matrix_to_bit_repr_of_facets`,
        :meth:`incidence_matrix_to_bit_repr_of_vertices`,
        :meth:`facets_tuple_to_bit_repr_of_facets`,
        :meth:`facets_tuple_to_bit_repr_of_vertices`,
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
        Initialize :class:`ListOfFaces`.

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

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.base.ListOfFaces).run()
        """
        self.nr_faces = nr_faces
        self.nr_vertices = nr_vertices
        self._mem = MemoryAllocator()

        # ``data`` will point to the faces as ``*uint64_t``.
        self.data = <uint64_t **> self._mem.allocarray(nr_faces, sizeof(uint64_t *))

        # ``face_length`` is the length in terms of ``uint64_t``
        # NOTE: This needs to be divisible by 2, if chunksize is 128
        #       and divisible by 4, if chunksize is 256.
        self.face_length = ((nr_vertices - 1)//chunksize + 1)*chunksize//64


        cdef size_t i
        for i in range(nr_faces):
            # Allocate the memory for the i-th face.
            # We must allocate the memory for ListOfFaces overaligned:
            # - must be 16-byte aligned if chunksize = 128
            # - must be 32-byte aligned if chunksize = 256
            self.data[i] = <uint64_t *> \
                self._mem.aligned_malloc(chunksize//8, self.face_length*8)

    def __reduce__(self):
        r"""
        Override __reduce__ to indicate that pickle/unpickle will not work.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base import ListOfFaces
            sage: faces = ListOfFaces(3,5)
            sage: faces1 = loads(faces.dumps())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

cpdef int calculate_dimension(ListOfFaces faces) except -2:
    r"""
    Calculate the dimension of a Polyhedron by its facets.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import calculate_dimension, facets_tuple_to_bit_repr_of_facets, \
        ....:            facets_tuple_to_bit_repr_of_vertices
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: facets = facets_tuple_to_bit_repr_of_facets(bi_pyr, 6)
        sage: vertices = facets_tuple_to_bit_repr_of_vertices(bi_pyr, 6)
        sage: calculate_dimension(facets)
        3
        sage: calculate_dimension(vertices)
        3

    ALGORITHM:

    This is done by iteration:

    Calculates the facets of one of the facets (i.e. the ridges contained in
    one of the facets). Then calculates the dimension of the facet, by
    considering its facets.

    Repeats until a face has only one facet. Usually this is a vertex.

    However, in the unbounded case, this might be different. The face with only
    one facet might be a ray or a line. So the correct dimension of a
    Polyhedron with one facet is the number of ``[lines, rays, vertices]``
    that the facet contains.

    Hence, we know the dimension of a face, which has only one facet and
    iteratively we know the dimension of entire Polyhedron we started with.

    TESTS::

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
        ....:     import calculate_dimension, incidence_matrix_to_bit_repr_of_facets, \
        ....:            incidence_matrix_to_bit_repr_of_vertices
        sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
        ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
        sage: for _ in range(10):
        ....:     points = tuple(tuple(randint(-1000,1000) for _ in range(10))
        ....:                    for _ in range(randint(3,15)))
        ....:     P = Polyhedron(vertices=points)
        ....:     facets = incidence_matrix_to_bit_repr_of_facets(P.incidence_matrix())
        ....:     vertices = incidence_matrix_to_bit_repr_of_vertices(P.incidence_matrix())
        ....:     d1 = P.dimension()
        ....:     if d1 == 0:
        ....:         continue
        ....:     d2 = calculate_dimension(facets)
        ....:     d3 = calculate_dimension(vertices)
        ....:     if not d1 == d2 == d3:
        ....:         print('calculation_dimension() seems to be incorrect')
    """
    if faces.nr_faces == 0:
        raise TypeError("at least one face needed")
    return calculate_dimension_loop(faces.data, faces.nr_faces, faces.face_length)

cdef int calculate_dimension_loop(uint64_t **faces, size_t nr_faces,
                                  size_t face_length) except -2:
    r"""
    Calculate the dimension of a Polyhedron by its facets.

    INPUT:

    - ``facesdata`` -- facets in Bit-representation
    - ``nr_faces`` -- length of facesdata
    - ``face_length`` -- the length of each face in terms of ``uint64_t``

    OUTPUT:

    - dimension of the Polyhedron

    .. SEEALSO::

        :meth:`calculate_dimension`
    """
    if nr_faces == 0:
        raise TypeError("wrong usage of ``calculate_dimension_loop``,\n" +
                        "at least one face needed.")

    if nr_faces == 1:
        # We expect the face to be the empty Polyhedron.
        # Possibly it contains more than one vertex/rays/lines.
        # The dimension of a polyhedron with this face as only facet is
        # the number of atoms it contains.
        return count_atoms(faces[0], face_length)

    # ``maybe_newfaces`` are all intersection of ``faces[nr_faces -1]`` with previous faces.
    # It needs to be allcoated to store those faces.
    cdef ListOfFaces maybe_newfaces_mem = ListOfFaces(nr_faces, face_length*64)
    cdef uint64_t **maybe_newfaces = maybe_newfaces_mem.data

    # ``newfaces`` point to the actual facets of ``faces[nr_faces -1]``.
    cdef MemoryAllocator newfaces_mem = MemoryAllocator()
    cdef uint64_t **newfaces = <uint64_t **> newfaces_mem.allocarray(nr_faces, sizeof(uint64_t *))

    # Calculating ``maybe_newfaces`` and ``newfaces``
    # such that ``newfaces`` points to all facets of ``faces[nr_faces -1]``.
    cdef size_t new_nr_faces
    sig_on()
    new_nr_faces = get_next_level(faces, nr_faces, maybe_newfaces,
                                  newfaces, NULL, 0, face_length)
    sig_off()

    # Calculate the dimension of the polyhedron,
    # by calculating dimension of one of its faces.
    return calculate_dimension_loop(newfaces, new_nr_faces, face_length) + 1

cdef class FaceIterator(SageObject):
    r"""
    A class to iterate over all faces of a Polyhedron.

    Constructs all proper from the facets. In dual mode, constructs all proper
    faces from the vertices. Dual will be faster for less vertices than facets.

    INPUT:

    - ``C`` -- a :class:`CombinatorialPolyhedron`
    - ``dual`` -- if True, then dual Polyhedron is used for iteration
      (only possible for bounded Polyhedra)

    .. SEEALSO::

        :class:`CombinatorialPolyhedron`.

    EXAMPLES:

    Construct a FaceIterator::

        sage: P = polytopes.cuboctahedron()
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_iter()

    By default it will give the dimension of each face::

        sage: [next(it) for _ in range(14)]
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]

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
        sage: [next(it) for _ in range(5)]
        [1, 1, 2, 2, 1]
        sage: it.length_facet_repr()
        2

    Construct faces by the dual or not::

        sage: it = C.face_iter(dual=False)
        sage: next(it)
        2
        sage: next(it)
        2
        sage: it.ignore_subfaces()
        sage: it.ignore_supfaces()
        Traceback (most recent call last):
        ...
        ValueError: only possible when in dual mode
        sage: it = C.face_iter(dual=True)
        sage: next(it)
        0
        sage: next(it)
        0
        sage: it.ignore_supfaces()
        sage: it.ignore_subfaces()
        Traceback (most recent call last):
        ...
        ValueError: only possible when not in dual mode

    ALGORITHM:

    The algorithm to visit all proper faces exactly once is roughly
    equivalent to::

        faces = [set(facet) for facet in P.facets()]
        face_iterator(faces, [])

        def face_iterator(faces, visited_all):
            # Visit all faces of a Polyhedron `P`, except those contained in
            # any of the visited all.

            # Assumes ``faces`` to be excactly those facets of `P`
            # that are not contained in any of the ``visited_all``.

            # Assumes ``visited_all`` to be some list of faces of
            # a Polyhedron `P_2`, which contains `P` as one of its faces.

            while len(facets) > 0:
                one_face = faces.pop()
                maybe_newfaces = [one_face.intersection(face) for face in faces]

                # ``maybe_newfaces`` contains all facets of ``one_face``,
                # which we have not visited before.
                # Proof: Let `F` be a facet of ``one_face``.
                # We have a chain:
                # `P` ⊃ ``one_face`` ⊃ `F`.
                # By diamond property there exists ``second_face`` with:
                # `P` ⊃ ``second_face`` ⊃ `F`.

                # Either ``second_face`` is not an element of ``faces``:
                #     Hence ``second_face`` is contained in one of ``visited_all``.
                #     In particular, `F` is contained in ``visited_all``.
                # Or ``second_face`` is an element of ``faces``:
                #     Then, intersecting ``one_face`` with ``second_face`` gives
                #     ``F``. ∎

                # Let ``maybe_newfaces2`` be the inclusion maximal faces of
                # ``maybe_newfaces``.
                # If an element in ``maybe_newfaces`` is inclusion maximal
                # and not contained in any of the ``visited_all``,
                # it is a facet of ``one_face``.
                # Any facet in ``maybe_newfaces`` of ``one_face``
                # is inlcusion maximal.
                maybe_newfaces2 = []
                for face1 in maybe_newfaces:
                    # ``face1`` is a facet of ``one_face``,
                    # iff it is not contained in another facet.
                    if all(not face1 < face2 for face2 in maybe_newfaces):
                        maybe_newfaces2.append(face1)

                # ``maybe_newfaces2`` contains only facets of ``one_face``
                # and some faces contained in any of ``visited_all``.
                # It also contains all the facets not contained in any of ``visited_all``.
                # Let ``newfaces`` be the list of all facets of ``one_face``
                # not contained in any of ``visited_all``.
                newfaces = []
                for face1 in maybe_newfaces2:
                    if all(not face1 < face2 for face2 in visited_all):
                        newfaces.append(face1)

                # By induction we can apply the algorithm, to visit all
                # faces of ``one_face`` not contained in ``visited_all``:
                face_iterator(newfaces, visited_all)

                # Finally visit ``one_face`` and add it to ``visited_all``:
                visit(one_face)
                visited_all.append(one_face)

                # Note: At this point, we have visited exactly those faces,
                # contained in any of the ``visited_all``.
    """
    def __init__(self, CombinatorialPolyhedron C, bint dual):
        r"""
        Initialize :class:`FaceIterator`.

        See :class:`FaceIterator`.

        EXAMPLES::

            sage: P = polytopes.permutahedron(4)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: f_vector = [1, 0, 0, 0, 1]
            sage: for d in it: f_vector[d+1] += 1
            sage: print ('f_vector of permutahedron(4): ', f_vector)
            f_vector of permutahedron(4):  [1, 24, 36, 14, 1]

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.base.FaceIterator).run()
        """
        if dual and C._unbounded:
            raise ValueError("cannot iterate over dual of unbounded Polyedron")
        cdef int i
        cdef ListOfFaces some_list  # make Cython aware of type

        self.dual = dual
        self.face = NULL
        self.dimension = C.dimension()
        self.current_dimension = self.dimension -1
        self.nr_lines = C._nr_lines
        self.request_dimension = -2
        self._mem = MemoryAllocator()
        self.lowest_dimension = self.nr_lines
        # We will not yield the empty face.
        # If there are `n` lines, than there
        # are no faces below dimension `n`.
        # The dimension of the level-sets in the face lattice jumps from `n` to `-1`.
        if dual:
            self.atoms = C.bitrep_facets
            self.coatoms = C.bitrep_vertices
        else:
            self.coatoms = C.bitrep_facets
            self.atoms = C.bitrep_vertices
        self.face_length = self.coatoms.face_length
        self._V = C._V
        self._H = C._H
        self._equalities = C._equalities

        self.atom_repr = <size_t *> self._mem.allocarray(self.coatoms.nr_vertices, sizeof(size_t))
        self.coatom_repr = <size_t *> self._mem.allocarray(self.coatoms.nr_faces, sizeof(size_t))

        if self.dimension == 0 or self.coatoms.nr_faces == 0:
            # As we will only yield proper faces,
            # there is nothing to yield in those cases.
            # We have to discontinue initialization,
            # as it assumes ``self.dimension > 0`` and ``self.nr_faces > 0``.
            self.current_dimension = self.dimension
            return
        # We may assume ``dimension > 0`` and ``nr_faces > 0``.

        # Initialize ``maybe_newfaces``,
        # the place where the new faces are being stored.
        self.newfaces_lists = tuple(ListOfFaces(self.coatoms.nr_faces, self.coatoms.nr_vertices)
                                    for i in range(self.dimension -1))
        self.maybe_newfaces = <uint64_t ***> self._mem.allocarray((self.dimension -1), sizeof(uint64_t **))
        for i in range(self.dimension -1):
            some_list = self.newfaces_lists[i]
            self.maybe_newfaces[i] = some_list.data

        # Initialize ``visited_all``.
        self.visited_all = <uint64_t **> self._mem.allocarray(self.coatoms.nr_faces, sizeof(uint64_t *))
        self.nr_visited_all = <size_t *> self._mem.allocarray(self.dimension, sizeof(size_t))
        self.nr_visited_all[self.dimension -1] = 0

        # Initialize ``newfaces``, which will point to the new faces of codimension 1,
        # which have not been visited yet.
        self.newfaces = <uint64_t ***> self._mem.allocarray(self.dimension, sizeof(uint64_t **))
        for i in range(self.dimension - 1):
            self.newfaces[i] = <uint64_t **> self._mem.allocarray(self.coatoms.nr_faces, sizeof(uint64_t *))
        self.newfaces[self.dimension - 1] = self.coatoms.data  # we start with coatoms

        # Initialize ``nr_newfaces``.
        self.nr_newfaces = <size_t *> self._mem.allocarray(self.dimension, sizeof(size_t))
        self.nr_newfaces[self.dimension - 1] = self.coatoms.nr_faces

        # Initialize ``first_time``.
        self.first_time = <bint *> self._mem.allocarray(self.dimension, sizeof(bint))
        self.first_time[self.dimension - 1] = True

        self.yet_to_visit = self.coatoms.nr_faces

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: P = polytopes.associahedron(['A',3])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_iter()
            Iterator over the faces of a Polyhedron of dimension 3
        """
        return "Iterator over the faces of a Polyhedron of dimension %s"%self.dimension

    def __next__(self):
        r"""
        Visit the next face and return its dimension.

        EXAMPLES::
            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: [next(it) for _ in range(7)]
            [2, 2, 2, 2, 2, 2, 1]
        """
        cdef int d = self.next_face()
        if unlikely(d == self.dimension):
            raise StopIteration

        # If ``dual == 0`` return current dimension,
        # if ``dual == 1`` translate current dimension to dual and then return.
        return Integer(self.dual*(self.dimension-1-d) + (1-self.dual)*d)

    next = __next__

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: [d for d in it]
            [2, 2, 2, 2, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1]
        """
        return self

    def __reduce__(self):
        r"""
        Override __reduce__ to indicate that pickle/unpickle will not work.

        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: it1 = loads(it.dumps())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def set_request_dimension(self, dim):
        r"""
        Set the iterator to only yield faces of dimension ``dim``.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: next(it)
            3
            sage: counter = 0
            sage: it.set_request_dimension(2)
            sage: for _ in it: counter += 1
            sage: print ('permutahedron(5) has', counter,
            ....:        'faces of dimension 2')
            permutahedron(5) has 150 faces of dimension 2
            sage: C.f_vector()
            (1, 120, 240, 150, 30, 1)
        """
        if self.dual:
            # In dual mode, the dimensions are reversed.
            self.request_dimension = self.dimension - 1 - dim
        else:
            self.request_dimension = dim
        self.lowest_dimension = max(self.nr_lines, self.request_dimension)

    def get_dimension(self):
        r"""
        Return the dimension of the current face.

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
        # If ``dual == 0`` return current dimension,
        # if ``dual == 1`` translate current dimension to dual and then return.
        return Integer(self.dual*(self.dimension-1-self.current_dimension) +
                       (1-self.dual)*self.current_dimension)

    def vertex_repr(self, names=True):
        r"""
        Return the vertex-representation of the current face.

        The vertex-representation consists of
        the ``[vertices, rays, lines]`` that face contains.

        INPUT:

        - ``names`` -- if ``True`` returns the names of the ``[vertices, rays, lines]``
          as given on initialization of the :class:`CombinatorialPolyhedron`

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
            # if dual, the vertex-represention corresponds to the coatom-representation
            length = self.set_coatom_repr()
            if names and self._V:
                return tuple(self._V[self.coatom_repr[i]]
                             for i in range(length))
            else:
                return tuple(Integer(self.coatom_repr[i])
                             for i in range(length))
        else:
            # if not dual, the vertex-represention corresponds to the atom-representation
            length = self.set_atom_repr()
            if names and self._V:
                return tuple(self._V[self.atom_repr[i]]
                             for i in range(length))
            else:
                return tuple(Integer(self.atom_repr[i])
                             for i in range(length))

    def length_vertex_repr(self):
        r"""
        Return the length of the :class:`vertex_repr`.

        Might be faster than `len(self.vertex_repr())`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: all(it.length_vertex_repr() == len(it.vertex_repr()) for _ in it)
            True
        """
        if unlikely(self.face is NULL):
            raise ValueError("iterator not set to a face yet")
        if self.dual:
            return Integer(self.set_coatom_repr())
        else:
            return Integer(self.length_atom_repr())

    def facet_repr(self, names=True):
        r"""
        Return the facet-representation of the current face.

        The facet-representation consists of the facets
        that contain the face and of the equalities of the Polyhedron.

        INPUT:

        - ``names`` -- if ``True`` returns the names of the ``[facets, equations]``
          as given on initialization of :class:`CombinatorialPolyhedron`

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
            sage: [next(it) for _ in range(3)]
            [0, 1, 1]
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
            # if not dual, the facet-represention corresponds to the coatom-representation
            length = self.set_coatom_repr()  # fill self.coatom_repr_face
            if names and self._H:
                return tuple(self._H[self.coatom_repr[i]]
                             for i in range(length)) + self._equalities
            else:
                return tuple(Integer(self.coatom_repr[i])
                             for i in range(length))
        else:
            # if dual, the facet-represention corresponds to the atom-representation
            length = self.set_atom_repr()  # fill self.atom_repr_face
            if names and self._H:
                return tuple(self._H[self.atom_repr[i]]
                             for i in range(length)) + self._equalities
            else:
                return tuple(Integer(self.atom_repr[i])
                             for i in range(length))

    def length_facet_repr(self):
        r"""
        Returns the length of the :meth:`facet_repr`.

        Might be faster then ``len(self.facet_repr())``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: all(it.length_facet_repr() == len(it.facet_repr()) for _ in it)
            True
        """
        if unlikely(self.face is NULL):
            raise ValueError("iterator not set to a face yet")
        if not self.dual:
            return Integer(self.set_coatom_repr())
        else:
            return Integer(self.length_atom_repr())

    def ignore_subfaces(self):
        r"""
        :class:`FaceIterator` will not visit any faces of the current face.

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

        # The current face is added to ``visited_all``.
        # This will make the iterator skip those faces.
        # Also, this face will not be added a second time to ``visited_all``,
        # as there are no new faces.
        self.visited_all[self.nr_visited_all[self.current_dimension]] = self.face
        self.nr_visited_all[self.current_dimension] += 1

    def ignore_supfaces(self):
        r"""
        :class:`FaceIterator` will not visit any faces of the current face.

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

        # The current face is added to ``visited_all``.
        # This will make the iterator skip those faces.
        # Also, this face will not be added a second time to ``visited_all``,
        # as there are no new faces.
        self.visited_all[self.nr_visited_all[self.current_dimension]] = self.face
        self.nr_visited_all[self.current_dimension] += 1

    cdef inline int next_face(self) except -1:
        r"""
        Set attribute ``face`` to the next face and return the dimension.

        Will return the dimension of the Polyhedron on failure.

        The function calls :meth:`FaceIterator.next_face_loop` until a new
        face is set or until the iterator is consumed.

        .. NOTE::

            The face_iterator can be prevented from visiting any subfaces
            (or supfaces in dual mode) as in :meth:`FaceIterator.ignore_subfaces`
            and :meth`FaceIterator.ignore_supfaces`.

            Those methods add the current face to ``visited_all`` before
            visiting sub-/supfaces instead of after. One cannot arbitralily
            add faces to ``visited_all``, as visited_all has a maximal length.
        """
        cdef int dim = self.dimension
        while (not self.next_face_loop()) and (self.current_dimension < dim):
            sig_check()
        return self.current_dimension

    cdef inline int next_face_loop(self) except -1:
        r"""
        Set attribute ``face`` to the next face. On success return `1`.
        Otherwise `0`. Needs to be recalled then.

        If ``self.current_dimension == self.dimension``, then the iterator is
        consumed.
        """
        if unlikely(self.current_dimension == self.dimension):
            # The function is not supposed to be called,
            # just prevent it from crashing.
            raise StopIteration

        # Getting ``[faces, nr_faces, nr_visited_all]`` according to dimension.
        cdef uint64_t **faces = self.newfaces[self.current_dimension]
        cdef size_t nr_faces = self.nr_newfaces[self.current_dimension]
        cdef size_t nr_visited_all = self.nr_visited_all[self.current_dimension]

        if (self.request_dimension > -2) and (self.request_dimension != self.current_dimension):
            # If only a specifice dimension was requested (i.e. ``self.request_dimension > 2``),
            # then we will not yield faces in other dimension.
            self.yet_to_visit = 0

        if self.yet_to_visit:
            # Set ``face`` to the next face.
            self.yet_to_visit -= 1
            self.face = faces[self.yet_to_visit]
            return 1

        if self.current_dimension <= self.lowest_dimension:
            # We will not yield the empty face.
            # We will not yield below requested dimension.
            self.current_dimension += 1
            return 0

        if nr_faces <= 1:
            # There will be no more faces from intersections.
            self.current_dimension += 1
            return 0

        # We will visit the last face now.
        self.nr_newfaces[self.current_dimension] -= 1
        nr_faces -= 1

        if not self.first_time[self.current_dimension]:
            # In this case there exists ``faces[nr_faces + 1]``, of which we
            # have visited all faces, but which was not added to
            # ``visited_all`` yet.
            self.visited_all[nr_visited_all] = faces[nr_faces + 1]
            self.nr_visited_all[self.current_dimension] += 1
            nr_visited_all = self.nr_visited_all[self.current_dimension]
        else:
            # Once we have visited all faces of ``faces[nr_faces]``, we want
            # to add it to ``visited_all``.
            self.first_time[self.current_dimension] = False

        # Get the faces of codimension 1 contained in ``faces[nr_faces]``,
        # which we have not yet visited.
        cdef size_t newfacescounter

        sig_on()
        newfacescounter = get_next_level(
            faces, nr_faces + 1, self.maybe_newfaces[self.current_dimension-1],
            self.newfaces[self.current_dimension-1],
            self.visited_all, nr_visited_all, self.face_length)
        sig_off()

        if newfacescounter:
            # ``faces[nr_faces]`` contains new faces.
            # We will visted them on next call, starting with codimension 1.

            # Setting the variables correclty for next call of ``next_face_loop``.
            self.current_dimension -= 1
            self.first_time[self.current_dimension] = True
            self.nr_newfaces[self.current_dimension] = newfacescounter
            self.nr_visited_all[self.current_dimension] = nr_visited_all
            self.yet_to_visit = newfacescounter
            return 0
        else:
            # ``faces[nr_faces]`` contains no new faces.
            # Hence there is no need to add it to ``visited_all``.
            # NOTE:
            #     For the methods ``ignore_subfaces`` and ``ignore_supfaces``
            #     this step needs to be done, as ``faces[nr_faces]`` might
            #     have been added manually to ``visited_all``.
            #     So this step is required to respect boundaries of ``visited_all``.
            self.first_time[self.current_dimension] = True
            return 0

    cdef size_t length_atom_repr(self) except -1:
        r"""
        Calculate the number of atoms in the current face by counting the
        number of set bits.
        """
        if self.face:
            return count_atoms(self.face, self.face_length)

        # The face was not initialized properly.
        raise LookupError("``FaceIterator`` does not point to a face")

    cdef size_t set_coatom_repr(self) except -1:
        r"""
        Set ``coatom_repr`` to be the coatom-representation of the current face.
        Return its length.
        """
        cdef size_t nr_coatoms = self.coatoms.nr_faces
        cdef uint64_t **coatoms = self.coatoms.data
        cdef size_t face_length = self.face_length
        return bit_repr_to_coatom_repr(self.face, coatoms, nr_coatoms,
                                       face_length, self.coatom_repr)

    cdef size_t set_atom_repr(self) except -1:
        r"""
        Set ``atom_repr`` to be the atom-representation of the current face.
        Return its length.
        """
        cdef size_t face_length = self.face_length
        return bit_repr_to_vertex_list(self.face, self.atom_repr, face_length)

cdef class ListOfAllFaces(SageObject):
    r"""
    A class to generate incidences of :class:`CombinatorialPolyhedron`.

    On initialization all faces of the given :class:`CombinatorialPolyhedron`
    are added and sorted (except coatoms). The incidences can be used to
    generate the ``face_lattice``.

    Might generate the faces of the dual Polyhedron for speed.

    INPUT:

    - :class:`CombinatorialPolyhedron`

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

    The faces are recorded with :class:`FaceIterator` in Bit-representation.
    Once created, all level-sets but the coatoms are sorted with merge sort.
    Non-trivial incidences of elements whos rank differs by 1 are determined
    by intersecting with all coatoms. Then each intersection is looked up in
    the sorted level sets.
    """
    def __init__(self, CombinatorialPolyhedron C):
        r"""
        Initialize :class:`ListOfAllFaces`.

        See :class:`ListOfAllFaces`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: C._record_all_faces() # indirect doctests
            sage: C.face_lattice()
            Finite lattice containing 28 elements

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.base.ListOfAllFaces).run()
        """
        self._mem = MemoryAllocator()
        self.dimension = C.dimension()
        self.dual = False
        if C.bitrep_facets.nr_faces > C.bitrep_vertices.nr_faces:
            self.dual = True
        if C._unbounded:
            self.dual = False
        cdef FaceIterator face_iter = C._face_iter(self.dual)
        self.face_length = face_iter.face_length
        self._V = C._V
        self._H = C._H
        self._equalities = C._equalities

        # copy f_vector for later use
        f_vector = C.f_vector()
        self.f_vector = <size_t *> self._mem.allocarray(self.dimension + 2, sizeof(size_t))
        if self.dual:
            for i in range(-1, self.dimension + 1):
                self.f_vector[i+1] = f_vector[-i-2]
        else:
            for i in range(-1, self.dimension + 1):
                self.f_vector[i+1] = f_vector[i+1]

        # face_counter keeps track, if all faces have been added already
        self.face_counter = <size_t *> self._mem.calloc(self.dimension + 2, sizeof(size_t))
        self.face_counter[0] = 1
        self.face_counter[self.dimension + 1] = 1
        if self.dimension > -1:
            # We will obtain the coatoms from ``CombinatorialPolyhedron``.
            self.face_counter[self.dimension] = self.f_vector[self.dimension]

        # Initialize atoms, coatoms, ``atom_repr`` and ``coatom_repr``.
        if self.dimension == 0:
            # In case of the 0-dimensional Polyhedron, we have to fix atoms and coatoms.
            # So far this didn't matter, as we only iterated over proper faces.
            self.atoms = facets_tuple_to_bit_repr_of_vertices(((),), 1)
            self.coatoms = facets_tuple_to_bit_repr_of_facets(((),), 1)
            self.face_length = self.coatoms.face_length
        else:
            self.atoms = face_iter.atoms
            self.coatoms = face_iter.coatoms
        cdef size_t nr_atoms = self.atoms.nr_faces
        self.atom_repr = <size_t *> self._mem.allocarray(self.coatoms.nr_vertices, sizeof(size_t))
        self.coatom_repr = <size_t *> self._mem.allocarray(self.coatoms.nr_faces, sizeof(size_t))

        # Initialize the data for ``faces``:
        cdef ListOfFaces coatoms_mem
        self.faces_mem = tuple(ListOfFaces(self.f_vector[i+1], nr_atoms)
                               for i in range(-1, self.dimension-1))
        if self.dimension > -1:
            # the coatoms
            self.faces_mem += (self.coatoms,)
        self.faces_mem += (ListOfFaces(1, nr_atoms),)  # the full Polyhedron

        # Setting up a pointer to raw data of ``faces``:
        self.faces = <uint64_t ***> self._mem.allocarray(self.dimension + 2, sizeof(uint64_t **))
        cdef ListOfFaces some_list  # assuming a type
        for i in range(self.dimension + 2):
            some_list = self.faces_mem[i]
            self.faces[i] = some_list.data

        if self.dimension != 0:
            # Initialize the empty face.
            # In case ``dimension == 0``, we would overwrite the coatoms.
            vertex_list_to_bit_repr((), self.faces[0][0], self.face_length)
        # Intialize the full polyhedron
        vertex_list_to_bit_repr(tuple(j for j in range(nr_atoms)),
                                self.faces[self.dimension + 1][0],
                                self.face_length)

        # Attributes for iterating over the incidences.
        self.is_incidence_initialized = 0
        cdef ListOfFaces incidence_face_mem = ListOfFaces(1, nr_atoms)
        self.incidence_face = incidence_face_mem.data[0]
        self.faces_mem += (incidence_face_mem,)  # needs to be stored somewhere

        # Adding all faces, using the iterator.
        cdef int d
        if face_iter.current_dimension != self.dimension:
            # If there are proper faces.
            d = face_iter.next_face()
            while (d == self.dimension - 1):
                # We already have the coatoms.
                d = face_iter.next_face()
            while (d < self.dimension):
                self._add_face(d, face_iter.face)
                d = face_iter.next_face()

        # Sorting the faces, except for coatoms.
        self._sort()

    def __reduce__(self):
        r"""
        Override __reduce__ to indicate that pickle/unpickle will not work.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base import ListOfAllFaces
            sage: P = polytopes.simplex()
            sage: C = CombinatorialPolyhedron(P)
            sage: all_faces = ListOfAllFaces(C)
            sage: all_faces1 = loads(all_faces.dumps())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    cdef int _add_face(self, int face_dim, uint64_t *face) except -1:
        r"""
        Add a face to :class:`ListOfAllFaces`.

        This method is used at initialization only.
        """
        cdef size_t counter = self.face_counter[face_dim + 1]
        cdef size_t max_number = self.f_vector[face_dim + 1]
        if unlikely(counter >= max_number):
            raise IOError("trying to add too many faces to ``ListOfAllFaces``")

        # Actually add the face by copying its data.
        memcpy(self.faces[face_dim + 1][counter], face, self.face_length*8)
        self.face_counter[face_dim + 1] += 1

    cdef int _sort(self) except -1:
        r"""
        Sort each list of ``self.faces`` (except for coatoms).

        This method is used on initialization only.
        """
        cdef int dim = self.dimension
        cdef int i
        for i in range(dim + 2):
            if unlikely(self.f_vector[i] != self.face_counter[i]):
                raise ValueError("``ListOfAllFaces`` does not contain all faces")

        for i in range(0, dim):
            # Sort each level set, except for coatoms, full- and empty polyhedron.
            self._sort_one_list(self.faces[i], self.f_vector[i])

    cdef int _sort_one_list(self, uint64_t **faces, size_t nr_faces) except -1:
        r"""
        Sort ``faces`` of length ``nr_faces``.

        See :meth:`sort`.
        """
        cdef MemoryAllocator mem = MemoryAllocator()

        # Merge sort needs a second list of pointers.
        cdef uint64_t **extra_mem = <uint64_t **> mem.allocarray(nr_faces, sizeof(uint64_t *))

        # Sort the faces using merge sort.
        self._sort_one_list_loop(faces, faces, extra_mem, nr_faces)

    cdef int _sort_one_list_loop(
            self, uint64_t **inp, uint64_t **output1,
            uint64_t **output2, size_t nr_faces) except -1:
        r"""
        This is merge sort.

        Sorts ``inp`` and returns it in ``output1``.

        ..WARNING::

            Input is the same as output1 or output2

        See :meth:`sort`.
        """
        if unlikely(nr_faces == 0):
            # Prevent it from crashing.
            # In this case there is nothing to do anyway.
            return 0

        if nr_faces == 1:
            # The final case, where there is only one element.
            output1[0] = inp[0]
            return 0

        cdef size_t middle = nr_faces//2
        cdef size_t len_upper_half = nr_faces - middle

        # Sort the upper and lower half of ``inp`` iteratively into ``output2``.
        self._sort_one_list_loop(inp, output2, output1, middle)
        self._sort_one_list_loop(&(inp[middle]), &(output2[middle]),
                                 &(output1[middle]), len_upper_half)

        # Merge lower and upper half into ``output1``.
        cdef size_t i = 0        # index through lower half
        cdef size_t j = middle   # index through upper half
        cdef size_t counter = 0  # counts how many elements have been "merged" already
        while i < middle and j < nr_faces:
            # Compare the lowest elements of lower and upper half.
            if self.is_smaller(output2[i], output2[j]):
                output1[counter] = output2[i]
                i += 1
                counter += 1
            else:
                output1[counter] = output2[j]
                j += 1
                counter += 1
        if i < middle:
            # Add the remaining elements of lower half.
            while i < middle:
                output1[counter] = output2[i]
                i += 1
                counter += 1
        else:
            # Add the remaining elements of upper half.
            while j < nr_faces:
                output1[counter] = output2[j]
                j += 1
                counter += 1

    cdef inline size_t find_face(self, int dimension, uint64_t *face) except -1:
        r"""
        Return the index of ``face``, if it is of dimension ``dimension``.

        .. NOTE::

            Will give an index no matter if ``face`` is actual of dimension
            ``dimension``. Check the result with :meth:`is_equal`.

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
            ...
            ValueError: cannot find a facet, as those are not sorted
            sage: it.set_request_dimension(1)
            sage: S = set(find_face_from_iterator(it, C) for _ in it)
            sage: S == set(range(36))
            True
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
        cdef uint64_t **faces = self.faces[dimension + 1]

        while (nr_faces > 1):
            # In each iteration step, we will look for ``face`` in
            # ``faces[start:start+nr_faces]``.
            middle = nr_faces//2
            if self.is_smaller(face, faces[middle + start]):
                # If face is in the list, then in the lower half.
                # Look for face in ``faces[start : start + middle]`` in next step.
                nr_faces = middle
            else:
                # If face is in the list, then in the upper half.
                # Look for face in ``faces[start+middle:start+nr_faces]``, i.e.
                # ``faces[start + middle : (start + middle) + nr_faces - middle]``.
                nr_faces -= middle
                start += middle
        return start

    cdef inline bint is_smaller(self, uint64_t *one, uint64_t *two):
        r"""
        Return `1` if ``one`` is smaller than ``two``, otherwise `0`.
        """
        return memcmp(one, two, self.face_length*8) < 0

    cdef inline int is_equal(self, int dimension, size_t index, uint64_t *face) except -1:
        r"""
        Check wether ``face`` is of dimension ``dimension`` with index ``index``.

        This is used to validate the output of :meth:`find_face`.
        """
        if unlikely(dimension < -1 or dimension > self.dimension
                    or index >= self.f_vector[dimension + 1]):
            raise IndexError()
        cdef uint64_t *face2 = self.faces[dimension+1][index]
        cdef size_t i
        return (0 == memcmp(face, face2, self.face_length*8))

    def vertex_repr(self, dimension, index, names=True):
        r"""
        Return the vertex-representation of the face of dimension ``dimension``
        and index ``index``.

        The vertex-representation consists of
        the ``[vertices, rays, lines]`` that face contains.

        INPUT:

        - ``dimension`` -- dimension of the face
        - ``index`` -- index of the face
        - ``names`` -- if ``True`` returns the names of the ``[vertices, rays, lines]``
          as given on initialization of :class:`CombinatorialPolyhedron`

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
            sage: it.set_request_dimension(1)
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
            # if dual, the vertex-represention corresponds to the coatom-representation
            dimension = self.dimension - 1 - dimension  # if dual, the dimensions are reversed
            length = self.set_coatom_repr(dimension, index)
            if names and self._V:
                return tuple(self._V[self.coatom_repr[i]]
                             for i in range(length))
            else:
                return tuple(Integer(self.coatom_repr[i])
                             for i in range(length))
        else:
            # if not dual, the vertex-represention corresponds to the atom-representation
            length = self.set_atom_repr(dimension, index)
            if names and self._V:
                return tuple(self._V[self.atom_repr[i]]
                             for i in range(length))
            else:
                return tuple(Integer(self.atom_repr[i])
                             for i in range(length))

    def facet_repr(self, dimension, index, names=True):
        r"""
        Return the facet-representation of the face of dimension ``dimension``
        and index ``index``.

        The facet-representation consists of the facets
        that contain the face and of the equalities of the Polyhedron.

        INPUT:

        - ``dimension`` -- dimension of the face
        - ``index`` -- index of the face
        - ``names`` -- if ``True`` returns the names of the ``[facets, equations]``
          as given on initialization of :class:`CombinatorialPolyhedron`

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
            sage: it.set_request_dimension(1)
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

            sage: P = Polyhedron(vertices=[[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice_facet_repr(0)
            (An equation (0, 1) x - 1 == 0, An equation (1, 0) x + 0 == 0)
            sage: C.face_lattice_facet_repr(1)
            (An equation (0, 1) x - 1 == 0, An equation (1, 0) x + 0 == 0)
            sage: C.face_lattice_vertex_repr(0)
            ()
            sage: C.face_lattice_vertex_repr(1)
            (A vertex at (0, 1),)

            sage: P = Polyhedron()
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice_facet_repr(0)
            (An equation -1 == 0,)

            sage: P = Polyhedron(lines=[[0,1]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C.face_lattice_facet_repr(0)
            (An equation (1, 0) x + 0 == 0,)
            sage: C.face_lattice_facet_repr(1)
            (An equation (1, 0) x + 0 == 0,)
        """
        cdef size_t length
        if not self.dual:
            # if not dual, the vertex-represention corresponds to the coatom-representation
            length = self.set_coatom_repr(dimension, index)
            if unlikely((self.coatoms.nr_faces == 0 or self.dimension == 0)
                        and names and self._H is not None):
                # in this case the facet does not correspond to a Hrep
                return self._equalities + self._H
            elif names and self._H:
                return tuple(self._H[self.coatom_repr[i]]
                             for i in range(length)) + self._equalities
            else:
                return tuple(Integer(self.coatom_repr[i])
                             for i in range(length))
        else:
            # if dual, the facet-represention corresponds to the atom-representation
            dimension = self.dimension - 1 - dimension  # if dual, the dimensions are reversed
            length = self.set_atom_repr(dimension, index)
            if names and self._H:
                return tuple(self._H[self.atom_repr[i]]
                             for i in range(length)) + self._equalities
            else:
                return tuple(Integer(self.atom_repr[i])
                             for i in range(length))

    cdef size_t set_coatom_repr(self, int dimension, size_t index) except -1:
        r"""
        Set ``atom_repr`` to be the atom-representation of the face
        of dimension ``dimension`` and index ``index``.
        Return its length.

        .. SEEALSO::

            :class:`ListOfAllFaces`.
        """
        if unlikely(dimension < -1 or dimension > self.dimension):
            raise ValueError("no face of dimension %s"%dimension)
        if unlikely(index >= self.f_vector[dimension + 1]):
            raise IndexError("no %s-th face of dimension %s"%(index, dimension))
        if unlikely(self.coatoms.nr_faces == 0):
            return 0

        cdef size_t nr_coatoms = self.f_vector[self.dimension]
        cdef uint64_t **coatoms = self.faces[self.dimension]
        cdef size_t face_length = self.face_length
        cdef uint64_t *face = self.faces[dimension+1][index]
        return bit_repr_to_coatom_repr(face, coatoms, nr_coatoms,
                                       face_length, self.coatom_repr)

    cdef size_t set_atom_repr(self, int dimension, size_t index) except -1:
        r"""
        Set ``atom_repr`` to be the atom-representation of the face
        of dimension ``dimension`` and index ``index``.
        Return its length.

        .. SEEALSO::

            :class:`ListOfAllFaces`.
        """
        if unlikely(dimension < -1 or dimension > self.dimension):
            raise ValueError("no face of dimension %s"%dimension)
        if unlikely(index >= self.f_vector[dimension + 1]):
            raise IndexError("no %s-th face of dimension %s"%(index, dimension))

        cdef size_t face_length = self.face_length
        cdef uint64_t *face = self.faces[dimension+1][index]
        return bit_repr_to_vertex_list(face, self.atom_repr, face_length)

    cdef void incidence_init(self, int dimension_one, int dimension_two):
        r"""
        Initialize the :class:`ListOfAllFaces` to give incidences between
        ``dimension_one`` and ``dimension_two``.

        This will enable :meth:`next_incidence` to give all such incidences.

        Currently only ``dimension_one == dimension_two + 1`` and incidences
        with empty and full Polyhedron are implemented, which suffices for the
        face-lattice.
        """
        cdef size_t i
        if dimension_one == self.dimension:
            # The full polyhedron is incident to every face.
            if dimension_two < -1:
                raise ValueError("no faces of dimension %s"%dimension_two)
            if dimension_two > self.dimension:
                raise ValueError("no faces of dimension %s"%dimension_two)
            self.incidence_dim_one = dimension_one
            self.incidence_dim_two = dimension_two
            self.incidence_counter_one = 0
            self.incidence_counter_two = 0
            self.is_incidence_initialized = 2
            return

        if dimension_two == -1:
            # The empty polyhedron is incident to every face.
            if dimension_one < -1:
                raise ValueError("no faces of dimension %s"%dimension_two)
            if dimension_one > self.dimension:
                raise ValueError("no faces of dimension %s"%dimension_two)
            self.incidence_dim_one = dimension_one
            self.incidence_dim_two = dimension_two
            self.incidence_counter_one = 0
            self.incidence_counter_two = 0
            self.is_incidence_initialized = 3
            return

        if dimension_one != dimension_two + 1:
            raise ValueError("``dimension_one = dimension_two`` + 1 must hold")
            # At the moment, this is not implemented.

        if dimension_one > self.dimension:
            raise ValueError("no faces of dimension %s"%dimension_one)
        if dimension_two < -1:
            raise ValueError("no faces of dimension %s"%dimension_two)

        self.incidence_dim_one = dimension_one
        self.incidence_dim_two = dimension_two
        self.incidence_counter_one = 0
        self.incidence_counter_two = 0
        self.is_incidence_initialized = 1

    cdef inline bint next_incidence(self, size_t *one, size_t *two):
        r"""
        Set ``one[0]`` and ``two[0]`` to be the next incidence. Return ``True``
        unless there are no more incidences, then return `0`.

        After initialization with :meth:`next_incidence`, this method will give
        all incidences of faces of ``dimension_one`` and ``dimension_two``.
        ``one[0]`` will represent the index of a face in ``dimension_one`` and
        ``two[0]`` will represent the index of a face in ``dimension_two``
        according to their order in :class:`ListOfAllFaces`.

        Use :meth:`vertex_repr` and :meth:`facet_repr` to interpret the output.

        ALGORITHM:

        This is the algorithm for non-trivial cases:
        ``0 < self.dimension_two + 1 == self.dimension_one < self.dimension``

        We intersect each face of dimension ``incidence_dim_one`` with each
        coatom. The result will be looked up in the ``incidence_dim_two`` faces.
        """
        cdef bint result = False
        while ((not result)
                and (self.incidence_counter_one < self.f_vector[self.incidence_dim_one + 1])):
            # Calls next_incidence_loop, until it gives a result or
            # until there are no more incidences.
            result = self.next_incidence_loop(one, two)

        return result

    cdef inline bint next_incidence_loop(self, size_t *one, size_t *two):
        r"""
        Set ``one[0]`` and ``two[0]`` to be the next incidence. Return ``True``
        on success and ``False`` otherwise.

        If it returns ``False``, it needs to be called again, unless
        ``self.incidence_counter_one >= self.f_vector[self.incidence_dim_one + 1])``.

        See :meth:`next_incidence`.
        """
        cdef uint64_t **coatoms = self.faces[self.dimension]
        cdef uint64_t *dimension_one_face  # depending on the index ``incidence_counter_one``

        cdef size_t location  # the index the intersection has, if of correct dimension
        cdef int is_it_equal  # checks if face with index ``location`` is intersection

        if self.is_incidence_initialized == 1:
            # The standard case, where
            # ``0 < self.dimension_two + 1 == self.dimension_one < self.dimension``.

            one[0] = self.incidence_counter_one
            dimension_one_face = self.faces[self.incidence_dim_one + 1][self.incidence_counter_one]

            # Get the intersection of ``dimension_one_face`` with the
            # ``self.incidence_counter_two``-th coatom.
            intersection(dimension_one_face, coatoms[self.incidence_counter_two],
                         self.incidence_face, self.face_length)

            # Get the location of the intersection and
            # check, wether it is correct.
            location = self.find_face(self.incidence_dim_two, self.incidence_face)
            two[0] = location
            is_it_equal = self.is_equal(self.incidence_dim_two,
                                        location, self.incidence_face)

            # Set counters for next function call.
            self.incidence_counter_two += 1
            if self.incidence_counter_two == self.f_vector[self.dimension]:
                self.incidence_counter_one += 1
                self.incidence_counter_two = 0
            return is_it_equal

        if self.is_incidence_initialized == 2:
            # the case where ``dimension_one`` is dimension of Polyhedron.
            return self.next_trivial_incidence(one, two)

        if self.is_incidence_initialized == 3:
            # the case where ``dimension_two`` is `-1`.
            return self.next_trivial_incidence2(one, two)

        if self.is_incidence_initialized == 0:
            return 0

    cdef inline bint next_trivial_incidence(self, size_t *one, size_t *two):
        r"""
        Handling the case where ``dimension_one`` is dimension of Polyhedron.

        See :meth:`next_incidence`.
        """
        one[0] = 0
        two[0] = self.incidence_counter_two
        self.incidence_counter_two += 1

        # Once done, raising ``self.incidence_counter_one``, such that
        # :meth:`next_incidence` recognizes that we are done.
        if self.incidence_counter_two >= self.f_vector[self.incidence_dim_two + 1]:
            self.incidence_counter_one += 1

        return (two[0] < self.f_vector[self.incidence_dim_two + 1])

    cdef inline bint next_trivial_incidence2(self, size_t *one, size_t *two):
        r"""
        Handling the case where ``dimension_two`` is `-1`.

        See :meth:`next_incidence`.
        """
        two[0] = 0
        one[0] = self.incidence_counter_one
        self.incidence_counter_one += 1
        return (one[0] < self.f_vector[self.incidence_dim_one + 1])

cdef class CombinatorialPolyhedron(SageObject):
    r"""
    The class of the Combinatorial Type of a Polyehdron, a Polytope.

    INPUT:

    - ``data`` -- an instance of
      :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`

    or

    - ``data`` --  an instance of
      :class:`~sage.geometry.lattice_polytope.LatticePolytopeClass`

    or

    - ``data`` -- an ``incidence_matrix`` as in
      :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.incidence_matrix`

      * ``vertices`` -- a list of ``[vertices, rays, lines]``, if
        the rows in the incidence_matrix should correspond to names

      * ``facets`` -- a list of facets, if
        the columns in the incidence_matrix should correspond to names

      * ``nr_lines`` -- ``None`` for bounded Polyhedra,
        for unbounded Polyhedra, this needs to be set
        to the correct number of lines,
        i.e. the the maximum number of lines with
        linearly independent directions in the Polyehdron

    or

    - ``data`` -- a ``[list, tuple, iterator]`` of facets,
      each facet given as a list of ``[vertices, rays, lines]``
      if the Polyhedron is unbounded, then rays and lines are required
      if the Polyehdron contains no lines, the rays can be thought of as
      the vertices of the facets deleted from a bounded Polyhedron
      see :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`
      on how to use rays and lines

      * ``facets`` -- a list of names of the facets, if
        the facets given should correspond to names

      * ``nr_lines`` -- ``None`` for bounded Polyhedra,
        for unbounded Polyhedra, this needs to be set
        to the correct number of lines,
        i.e. the the maximum number of lines with
        linearly independent directions in the Polyehdron

    or

    - ``data`` -- an Integer, representing the dimension of a Polyhedron equal
      to its affine hull

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

    Input is an integer::

        sage: CombinatorialPolyhedron(-1).f_vector()
        (1,)
        sage: CombinatorialPolyhedron(0).f_vector()
        (1, 1)
        sage: CombinatorialPolyhedron(5).f_vector()
        (1, 0, 0, 0, 0, 0, 1)

    Specifying the number of lines is important::

        sage: P = Polyhedron(ieqs=[[1,-1,0],[1,1,0]])
        sage: C = CombinatorialPolyhedron(P) #this works fine
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices
        sage: data = P.incidence_matrix()
        sage: vert = P.Vrepresentation()

    Incorrect due to missing number of lines::

        sage: C = CombinatorialPolyhedron(data, vertices=vert)
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 3 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C.vertices()
        (A line in the direction (0, 1),
         A line in the direction (0, 1),
         A line in the direction (0, 1))

    Correct usage with number of lines specified::

        sage: C = CombinatorialPolyhedron(data, vertices=vert, nr_lines=1)
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 0 vertices
        sage: C.f_vector()
        (1, 0, 2, 1)
        sage: C.vertices()
        ()

    Initialization from Polyhedron will automatically specify number of lines::

        sage: P = Polyhedron(rays=[[1,0],[0,1]])
        sage: C = CombinatorialPolyhedron(P) # this works fine
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 1 vertices
        sage: data = P.incidence_matrix()
        sage: vert = P.Vrepresentation()

    Incorrect due to missing number of lines::

        sage: C = CombinatorialPolyhedron(data, vertices=vert)
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 3 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C.vertices()
        (A vertex at (0, 0), A vertex at (0, 0), A vertex at (0, 0))

    Correct usage with number of lines specified::

        sage: C = CombinatorialPolyhedron(data, vertices=vert, nr_lines=0)
        sage: C
        Combinatorial Type of a Polyhedron of dimension 2 with 1 vertices
        sage: C.f_vector()
        (1, 1, 2, 1)
        sage: C.vertices()
        (A vertex at (0, 0),)
    """
    def __init__(self, data, vertices=None, facets=None, nr_lines=None):
        r"""
        Initialize :class:`CombinatorialPolyhedron`.

        See :class:`CombinatorialPolyhedron`.

        TESTS::

            sage: C = CombinatorialPolyhedron([[0,1,2],[0,1,3],
            ....: [0,2,3],[1,2,3]])    # indirect doctests

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron).run()
        """
        self._dimension = -2  # a "NULL" value
        self._edges = NULL
        self._ridges = NULL
        self._face_lattice_incidences = NULL
        self._equalities = ()
        self._all_faces = None
        self._mem_tuple = ()

        # ``_length_edges_list`` should not be touched in an instance
        # of :class:`CombinatorialPolyhedron`. This number can be altered,
        # but should probably be a power of `2` (for memory usage).
        # ``_length_edges_list`` shouldn't be too small for speed and
        # shouldn't be too large, as ``ridges``, ``edges`` and ``incidences``
        # each have a memory overhead of
        # ``self._length_edges_list*2*sizeof(size_t *)``.
        self._length_edges_list = 16348

        if is_Polyhedron(data):
            # Input is ``Polyhedron``.
            vertices = data.Vrepresentation()
            facets = tuple(inequality for inequality in data.Hrepresentation())

            if not data.is_compact():
                self._unbounded = True
                self._nr_lines = int(data.n_lines())
            else:
                self._unbounded = False
                self._nr_lines = 0

            data = data.incidence_matrix()
        elif is_LatticePolytope(data):
            # Input is ``LatticePolytope``.
            self._unbounded = False
            self._nr_lines = 0
            vertices = data.vertices()
            self._length_Vrep = len(vertices)
            facets = data.facets()
            self._length_Hrep = len(facets)
            data = tuple(tuple(vert for vert in facet.vertices())
                         for facet in facets)
        else:
            # Input is different from ``Polyhedron`` and ``LatticePolytope``.
            if nr_lines is None:
                # According to input, the Polyhedron is bounded.
                self._unbounded = False
                self._nr_lines = 0
            else:
                # According to input, the Polyhedron is unbounded.
                self._unbounded = True
                self._nr_lines = int(nr_lines)

        if vertices:
            # The vertices have names, which the user might want later.
            self._V = tuple(vertices)
            self._Vinv = {v: i for i,v in enumerate(self._V)}
        else:
            self._V = None
            self._Vinv = None

        if facets:
            # The facets have names, which the user might want later.
            facets = tuple(facets)

            # Remove equalities from facets.
            test = [1] * len(facets)  # 0 if that facet is an equality
            for i in range(len(facets)):
                if hasattr(facets[i], "is_inequality"):
                    # At the moment this test only works for input being
                    if not facets[i].is_inequality():
                        test[i] = 0
            self._H = tuple(facets[i] for i in range(len(facets)) if test[i])

            # Saving the equalities.
            self._equalities = tuple(facets[i] for i in range(len(facets)) if not test[i])
        else:
            self._H = None

        if is_Matrix(data):
            # Input is incidence-matrix or was converted to it.
            self._length_Hrep = data.ncols()
            self._length_Vrep = data.nrows()

            # Initializing the facets in their Bit-representation.
            self.bitrep_facets = incidence_matrix_to_bit_repr_of_facets(data)

            # Initializing the vertices as their Bit-representation.
            self.bitrep_vertices = incidence_matrix_to_bit_repr_of_vertices(data)

            self._nr_facets = self.bitrep_facets.nr_faces

        elif isinstance(data, Integer):
            # To construct a trivial Polyhedron, equal to its affine hull,
            # one can give an Integer as Input.
            if data < -1:
                ValueError("any Polyhedron must have dimension at least -1")
            self._nr_facets = 0
            self._dimension = data

            # Initializing the facets in their Bit-representation.
            self.bitrep_facets = facets_tuple_to_bit_repr_of_facets((), 0)

            # Initializing the vertices as their Bit-representation.
            self.bitrep_vertices = facets_tuple_to_bit_repr_of_vertices((), 0)

        else:
            # Input is a "list" of facets.
            # The facets given by its ``[vertices, rays, lines]``.
            # Actually at least tuple, list, iterator will work.
            if is_iterator(data):
                data = tuple(data)

            if self._V is None:
                # Get the names of the vertices.
                vertices = sorted(set.union(*map(set, data)))
                nr_vertices = len(vertices)
                if vertices != range(len(vertices)):
                    self._V = tuple(vertices)
                    self._Vinv = {v: i for i,v in enumerate(self._V)}
            else:
                # Assuming the user gave as correct names for the vertices
                # and labeled them instead by `0,...,n`.
                nr_vertices = len(self._V)

            self._length_Vrep = nr_vertices

            # Relabel the vertices to be `0,...,n`.
            if self._V is not None:
                def f(v): return self._Vinv[v]
            else:
                def f(v): return int(v)
            facets = tuple(tuple(f(i) for i in j) for j in data)

            self._nr_facets = len(facets)
            self._length_Hrep = len(facets)

            # Initializing the facets in their Bit-representation.
            self.bitrep_facets = facets_tuple_to_bit_repr_of_facets(facets, nr_vertices)

            # Initializing the vertices as their Bit-representation.
            self.bitrep_vertices = facets_tuple_to_bit_repr_of_vertices(facets, nr_vertices)

    def _repr_(self):
        r"""
        Return a description of the Combinatorial Polyhedron.

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
            sage: tup = tuple((it.vertex_repr(), it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(), it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True

            sage: P = polytopes.cyclic_polytope(4,10)
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_iter()
            sage: it1 = C1.face_iter()
            sage: tup = tuple((it.vertex_repr(), it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(), it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0], [0,-1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_iter()
            sage: it1 = C1.face_iter()
            sage: tup = tuple((it.vertex_repr(), it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(), it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True

            sage: P = Polyhedron(rays=[[1,0,0], [-1,0,0],
            ....:                      [0,-1,0], [0,1,0]])
            sage: C = CombinatorialPolyhedron(P)
            sage: C1 = loads(C.dumps())
            sage: it = C.face_iter()
            sage: it1 = C1.face_iter()
            sage: tup = tuple((it.vertex_repr(), it.facet_repr()) for _ in it)
            sage: tup1 = tuple((it1.vertex_repr(), it1.facet_repr()) for _ in it1)
            sage: tup == tup1
            True
        """
        nr_lines = None
        if self._unbounded:
            nr_lines = Integer(self._nr_lines)
        # Give a constructor by list of facets.
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
        Return the dimension of the Polyhedron.

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
            # Dimension not calculated yet.
            if self._nr_facets == 0:
                # The dimension of a trivial Polyhedron is assumed to contain
                # exactly one "vertex" and for each dimension one "line" as in
                # :class:`~sage.geometry.polyhedron.parent.Polyhedron_base`
                self._dimension = self._length_Vrep - 1
            elif self._unbounded or self._nr_facets <= self._length_Vrep:
                self._dimension = calculate_dimension(self.bitrep_facets)
            else:
                # If the Polyhedron has many facets,
                # calculating the dimenion of the dual will be faster.
                # The dual exists, if the Polyhedron is bounded.
                self._dimension = calculate_dimension(self.bitrep_vertices)
        return Integer(self._dimension)

    def nr_vertices(self):
        r"""
        Return the number of vertices.

        Is equivalent to ``len(self.vertices())``.

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
            # This specific trivial Polyhedron needs special attention.
            return Integer(1)
        elif not self._unbounded:
            # In the unbounded case, we need to actually calculated the vertices,
            # the the V-representation contains also ``lines`` and ``rays``.
            return Integer(self._length_Vrep)
        else:
            return len(self.vertices())

    def vertices(self, names=True):
        r"""
        Return the elements in the ``Vrepresentation`` that are vertices.

        In case of an unbounded Polyhedron, there might be lines and
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
            # Handling the case of a trivial Polyhedron of dimension `0`.
            if names and self._V:
                return (self._V[0],)
            else:
                return (Integer(0),)
        dual = False
        if not self._unbounded:
            # In the bounded case, we already know all the vertices.
            # Whereas, in the unbounded case, some of those "vertices" might
            # be ``rays`` or ``lines``.
            dual = True

        # Get all `0`-dimensional faces from :meth:`face_iter`.
        face_iter = self.face_iter(0, dual=dual)
        return tuple(face_iter.vertex_repr(names=names)[0] for _ in face_iter)

    def nr_facets(self):
        r"""
        Return the number of facets.

        Is equivalent to ``len(self.facets())``.

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
            # This trivial Polyhedron needs special attention.
            return Integer(1)
        return Integer(self._nr_facets)

    def facets(self, names=True):
        r"""
        Return the facets as lists of ``[vertices, rays, lines]``.

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
            # Special attention for this trivial case.
            # There is actually one facet, but we have not initialized it.
            return ((),)

        # Get all facets from :meth:`face_iter`.
        face_iter = self.face_iter(self.dimension() - 1, dual=False)
        tup = tuple(face_iter.vertex_repr(names=names) for _ in face_iter)

        # It is important to have the facets in the exact same order as
        # on input, so that pickle/unpickle by :meth:`reduce` works.
        # Every facet knows its index by the facet-representation.
        face_iter = self.face_iter(self.dimension() - 1, dual=False)
        indices = tuple(face_iter.facet_repr(names=False)[0] for _ in face_iter)
        dic = {}
        for i in range(len(tup)):
            dic[indices[i]] = tup[i]
        return tuple(dic[i] for i in range(len(tup)))

    def edges(self, names=True):
        r"""
        Return the edges of the Polyhedron, i.e. the rank 1 faces.

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
            ....:     alarm(0.1)
            ....:     C.edges()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: len(C.edges())
            19900
        """
        cdef size_t len_edge_list = self._length_edges_list
        if self._edges is NULL:
            # Calculate the edges.
            if self._unbounded:
                self._calculate_edges(dual=False)
            elif self._length_Vrep > self._nr_facets*self._nr_facets:
                # This is a wild estimate
                # that in this case it is better not to use the dual.
                self._calculate_edges(dual=False)
            else:
                # In most bounded cases, one should use the dual.
                self._calculate_ridges(dual=True)
        if self._edges is NULL:
            raise ValueError('could not determine edges')

        # The edges are being saved in a list basically.
        # The first entry represents the first vertex of the first edge.
        # The second entry represents the second vertex of that edge.
        # The third entry represents the first vertex of the second edge etc.

        # Possibly there are many edges.
        # Hence, edges are stored in an array of arrays,
        # with each array containing ``len_edge_list`` of edges.

        # Mapping the indices of the vertices to the names, if requested.
        if self._V is not None and names is True:
            def f(size_t i): return self._V[i]
        else:
            def f(size_t i): return Integer(i)

        # Getting the indices of the `i`-th edge.
        def vertex_one(size_t i):
            return f(self._edges[i // len_edge_list][2*(i % len_edge_list)])
        def vertex_two(size_t i):
            return f(self._edges[i // len_edge_list][2*(i % len_edge_list)+1])

        cdef size_t j
        return tuple((vertex_one(j), vertex_two(j)) for j in range(self._nr_edges))

    def edge_graph(self, names=True):
        r"""
        Return the edge graph.

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
        Return the ridges.

        The ridges of a Polyhedron are the faces
        contained in exactly two facets.

        To obtain all faces of codimnesion 1 use
        :meth:`CombinatorialPolyhedron.face_iter` instead.

        The ridges will be given by the facets, they are contained in.

        INPUT:

        - ``add_equalities`` -- if ``True``, then equalities of the Polyhedron
          will be added

        - ``names`` -- if ``False``, then the facets are given by their indices

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
            ....:     alarm(0.1)
            ....:     C.ridges()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: len(C.ridges())
            19900
        """
        cdef size_t len_ridge_list = self._length_edges_list
        if self._ridges is NULL:
            # Calculate the ridges.
            if self._unbounded:
                self._calculate_ridges(dual=False)
            elif self._length_Vrep*self._length_Vrep < self._nr_facets:
                # This is a wild estimate
                # that in this case it is better to use the dual.
                self._calculate_edges(dual=True)
            else:
                # In most bounded cases, one should not use the dual.
                self._calculate_ridges(dual=False)
        if self._ridges is NULL:
            raise ValueError('could not determine ridges')
        nr_ridges = self._nr_ridges

        # The ridges are being saved in a list basically.
        # The first entry represents the first facet of the first ridge.
        # The second entry represents the second facet of that ridge.
        # The third entry represents the first facet of the second ridge etc.

        # Possibly there are many ridges.
        # Hence, ridges are stored in an array of arrays,
        # with each array containing ``len_ridge_list`` of ridges.

        # Mapping the indices of the vertices to the names, if requested.
        if self._H is not None and names is True:
            def f(size_t i): return self._H[i]
        else:
            def f(size_t i): return Integer(i)

        # Getting the indices of the `i`-th ridge.
        def facet_one(size_t i):
            return f(self._ridges[i // len_ridge_list][2*(i % len_ridge_list)])
        def facet_two(size_t i):
            return f(self._ridges[i // len_ridge_list][2*(i % len_ridge_list)+1])

        cdef size_t j
        if add_equalities:
            # Also getting the equalities for each facet.
            return tuple(
                ((self._equalities + (facet_one(i),)),
                 (self._equalities + (facet_two(i),)))
                for i in range(nr_ridges))
        else:
            return tuple((facet_one(i), facet_two(i))
                         for i in range(nr_ridges))

    def ridge_graph(self, names=True):
        r"""
        Return the ridge graph.

        The ridge graph of a Polyhedron consists of
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
        Calculate the ``f_vector`` of the Polyhedron.

        The ``f_vector`` contains the number of faces of dimension `k`
        for each `k` in ``range(-1, self.dimension() + 1)``.

        .. NOTE::

            To obtain edges and/or ridges as well, first do so. This might
            already calculate the ``f_vector``.

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
            sage: N = combinations(range(25),24)
            sage: C = CombinatorialPolyhedron(N)
            sage: try:
            ....:     alarm(0.5)
            ....:     C.f_vector()
            ....: except:
            ....:     print("alarm!")
            ....:
            alarm!
            sage: C.f_vector()  # long time
            (1,
             25,
             300,
             2300,
             12650,
             53130,
             177100,
             480700,
             1081575,
             2042975,
             3268760,
             4457400,
             5200300,
             5200300,
             4457400,
             3268760,
             2042975,
             1081575,
             480700,
             177100,
             53130,
             12650,
             2300,
             300,
             25,
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

        - ``dimension`` -- if specified, then iterate over only this dimension
        - ``dual`` -- if ``True``, iterate starting with the vertices,
          if ``False``, iterate starting with the facets,
          if ``None``, choose ``True`` or ``False`` for speed

        OUTPUT:

        - :class:`FaceIterator`

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
            # Determine the faster way, to iterate through all faces.
            if self._unbounded or self._nr_facets <= self._length_Vrep:
                dual = False
            else:
                dual = True

        face_iter = self._face_iter(int(dual))
        if dimension is not None:
            # Setting the face iterator to return only requested dimension.
            face_iter.set_request_dimension(dimension)
        return face_iter

    cdef FaceIterator _face_iter(self, bint dual):
        r"""
        A method to obtain the FaceIterator as Cython object.

        See :meth:`CombinatorialPolyhedron.face_iter`
        """
        if dual and self._unbounded:
            raise ValueError("cannot iterate over dual of unbounded Polyhedron")
        return FaceIterator(self, dual)

    def face_lattice(self):
        r"""
        Generate the face-lattice.

        OUTPUT:

        - :class:'~sage.combinat.posets.lattices.FiniteLatticePoset'

        .. NOTE::

            Use :meth:`CombinatorialPolyhedron.face_lattice_dimension` to get
            the dimension for each element in the Face Lattice.
            Use :meth:`CombinatorialPolyhedron.face_lattice_vertex_repr` to get
            the vertex representation for each element in the Face Lattice.
            Use :meth:`CombinatorialPolyhedron.face_lattice_facet_repr` to get
            the facet_repr for each element in the Face Lattice.

        .. WARNING::

            The labeling of the face lattice might depend on archicture
            and implementation.
            Relabeling the face lattice with
            :meth:`CombinatorialPolyhedron.face_lattice_vertex_repr` and/or
            :meth:`CombinatorialPolyhedron.face_lattice_facet_repr` will be
            platform independent.

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
        if not self._face_lattice_incidences:
            # Calculate all incidences.
            self._calculate_face_lattice_incidences()
        if self._face_lattice_incidences is NULL:
            raise TypeError("could not determine face lattice")

        cdef size_t **incidences = self._face_lattice_incidences
        cdef size_t nr_incidences = self._nr_face_lattice_incidences
        cdef size_t len_incidence_list = self._length_edges_list

        # The incidences are being saved in a list basically.
        # The first entry represents the first face of the first incidence.
        # The second entry represents the second face of that incidence.
        # The third entry represents the first face of the second incidence etc.

        # Possibly there are many incidences.
        # Hence, incidences are stored in an array of arrays,
        # with each array containing ``len_incidence_list`` of incidences.

        # Getting the indices of the `i`-th incidence.
        def face_one(size_t i):
            return Integer(incidences[i // len_incidence_list][2*(i % len_incidence_list)])
        def face_two(size_t i):
            return Integer(incidences[i // len_incidence_list][2*(i % len_incidence_list)+1])

        # Edges of the face-lattice/Hasse diagram.
        cdef size_t j
        edges = tuple((face_one(j), face_two(j))
                          for j in range(nr_incidences))

        V = tuple(range(sum(self._f_vector)))
        D = DiGraph([V, edges], format='vertices_and_edges')
        return FiniteLatticePoset(D)

    def face_lattice_dimension(self, index):
        r"""
        Return for each element in :meth:`CombinatorialPolyhedron.face_lattice`
        its dimension.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: C = CombinatorialPolyhedron(P)
            sage: F = C.face_lattice()
            sage: def f(i):
            ....:     return (i, C.face_lattice_dimension(i))
            ....:
            sage: G = F.relabel(f)
            sage: set(G._elements)
            {(0, -1),
             (1, 0),
             (2, 0),
             (3, 0),
             (4, 0),
             (5, 0),
             (6, 0),
             (7, 0),
             (8, 0),
             (9, 1),
             (10, 1),
             (11, 1),
             (12, 1),
             (13, 1),
             (14, 1),
             (15, 1),
             (16, 1),
             (17, 1),
             (18, 1),
             (19, 1),
             (20, 1),
             (21, 2),
             (22, 2),
             (23, 2),
             (24, 2),
             (25, 2),
             (26, 2),
             (27, 3)}
        """
        f_vector = self.f_vector()
        dim = self.dimension()

        # Getting the dimension, by considering the following:
        # The level-set of dimension `d` will have indices `k, k+1, ..., k+n-1`,
        # where `n` is the number of faces of dimension `d` ( ``n = f_vector[d + 1]``)
        # and `k` is the number of face of dimension up to `d`, i.e.
        # ``k = sum(f_vector[:d])``.
        return Integer(max(d for d in range(dim+2) if sum(f_vector[:d]) <= index)) - 1

    def face_lattice_vertex_repr(self, index, names=True):
        r"""
        Return for each element in :meth:`CombinatorialPolyhedron.face_lattice`
        its vertex-representation as in :meth:`ListOfAllFaces.vertex_repr` or
        :meth:`FaceIterator.vertex_repr`.

        If ``names`` is set to ``True``, then names of the
        ``[vertices, rays, lines]`` are used.

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
        self._record_all_faces()                            # Initalize ``_all_faces``, if not done yet.
        dim = self.face_lattice_dimension(index)            # Determine dimension to that index.
        newindex = index - sum(self._f_vector[:dim + 1])    # Index in that level-set.

        # Let ``_all_faces`` determine vertex-representation.
        return self._all_faces.vertex_repr(dim, newindex, names=names)

    def face_lattice_facet_repr(self, index, names=True):
        r"""
        Return for each element in :meth:`CombinatorialPolyhedron.face_lattice`
        its facet-representation as in :meth:`ListOfAllFaces.facet_repr` or
        :meth:`FaceIterator.facet_repr`.

        If ``names`` is set to ``True``, then names of the
        ``facets`` are used.

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
        self._record_all_faces()                            # Initalize ``_all_faces``, if not done yet.
        dim = self.face_lattice_dimension(index)            # Determine dimension to that index.
        newindex = index - sum(self._f_vector[:dim + 1])    # Index in that level-set.

        # Let ``_all_faces`` determine facet-representation.
        return self._all_faces.facet_repr(dim, newindex, names=names)

    cdef int _calculate_f_vector(self) except -1:
        r"""
        Calculate the ``f_vector`` of the Polyhedron.

        See :meth:`f_vector`.
        """
        if self._f_vector:
            return 0  # There is no need to recalculate the f_vector.

        cdef bint dual
        if self._unbounded or self._nr_facets <= self._length_Vrep:
            # In this case the non-dual approach is faster..
            dual = False
        else:
            # In this case the dual approach is faster.
            dual = True
        cdef FaceIterator face_iter = self._face_iter(dual)

        cdef int dim = self.dimension()
        cdef int d  # dimension of the current face of the iterator
        cdef MemoryAllocator mem = MemoryAllocator()

        # Initialize ``f_vector``.
        cdef size_t *f_vector = <size_t *> mem.calloc((dim + 2), sizeof(size_t))
        f_vector[0] = 1         # Face iterator will only visit proper faces.
        f_vector[dim + 1] = 1   # Face iterator will only visit proper faces.

        # For each face in the iterator, add `1` to the corresponding entry in
        # ``f_vector``.
        if self._nr_facets > 0 and dim > 0:
            d = face_iter.next_face()
            while (d < dim):
                sig_check()
                f_vector[d+1] += 1
                d = face_iter.next_face()

        # Copy ``f_vector``.
        if dual:
            # We have calculated the ``f_vector`` of the dual.
            # Reverse it:
            self._f_vector = \
                tuple(Integer(f_vector[dim+1-i]) for i in range(dim+2))
        else:
            self._f_vector = tuple(Integer(f_vector[i]) for i in range(dim+2))

    cdef int _calculate_edges(self, dual) except -1:
        r"""
        Calculate the edges of the Polyhedron.

        If ``dual`` is ``True``, calculate the edges of the dual. In this case,
        this will also calculate the ``f_vector``, if unknown.

        See :meth:`CombinatorialPolyhedron.edges` and :meth:`CombinatorialPolyhedron.ridges`.
        """
        if (self._edges is not NULL and not dual) or (self._ridges is not NULL and dual):
            return 0  # There is no need to recalculate.

        cdef MemoryAllocator mem = MemoryAllocator()
        cdef FaceIterator face_iter = self._face_iter(dual)
        cdef size_t len_edge_list = self._length_edges_list
        cdef int dim = self.dimension()
        cdef int d              # dimension of the current face of ``FaceIterator``
        cdef size_t *f_vector   # calculate f_vector, if not done already
        cdef bint is_f_vector   # True if f_vector was calculated previously

        # For each edge we determine its location in ``edges``
        # by ``edges[one][two]``
        cdef size_t **edges = <size_t**> mem.malloc(sizeof(size_t*))
        cdef size_t one, two

        cdef size_t counter = 0         # the number of edges so far
        cdef size_t current_length = 1  # dynamically enlarge **edges

        if self._f_vector:
            is_f_vector = True
        else:
            # in this case we will calculate the f_vector while we're at it
            is_f_vector = False

        if dim == 1:
            # In this case there is an edge, but its not a proper face.
            edges[0] = <size_t *> mem.allocarray(2, sizeof(size_t))
            edges[0][0] = 0
            edges[0][1] = 1
            counter = 1

            # Success, copy the data to ``CombinatorialPolyhedron``.
            if dual:
                # We have actually calculated the ridges.
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
            # Only calculate the edges.

            if not dual:
                face_iter.set_request_dimension(1)
            else:
                # :meth:`FaceIterator.set_request_dimension`
                # requires the dimension of the original Polyhedron
                face_iter.set_request_dimension(dim - 2)

            if self._nr_facets > 0 and dim > 0:
                # If not, there won't even be any edges. Prevent error message.

                while (face_iter.next_face() == 1):

                    # Determine the position in ``edges``.
                    one = counter // len_edge_list
                    two = counter % len_edge_list

                    # Enlarge ``edges`` if needed.
                    if unlikely(two == 0):
                        if unlikely(one + 1 > current_length):
                            # enlarge **edges
                            current_length *= 2
                            edges = <size_t **> mem.reallocarray(edges, current_length, sizeof(size_t*))

                        edges[one] = <size_t *> mem.allocarray(2 * len_edge_list, sizeof(size_t))

                    # Set up face_iter.atom_repr
                    face_iter.set_atom_repr()

                    # Copy the information.
                    edges[one][2*two] = face_iter.atom_repr[0]
                    edges[one][2*two + 1] = face_iter.atom_repr[1]
                    counter += 1

            # Success, copy the data to ``CombinatorialPolyhedron``.
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
            # While doing the edges one might as well do the f-vector.
            f_vector = <size_t *> mem.calloc(dim + 2, sizeof(size_t))
            f_vector[0] = 1         # This is not a proper face.
            f_vector[dim + 1] = 1   # This is not a proper face.

            counter = 0
            if self._nr_facets > 0 and dim > 0:
                # If not, there won't even be any edges. Prevent error message.

                d = face_iter.next_face()
                while (d < dim):
                    f_vector[d+1] += 1

                    if d == 1:
                        # If it is an edge.

                         # Determine the position in ``edges``.
                        one = counter // len_edge_list
                        two = counter % len_edge_list

                        # Enlarge ``edges`` if needed.
                        if unlikely(two == 0):
                            if unlikely(one + 1 > current_length):
                                # enlarge **edges
                                current_length *= 2
                                edges = <size_t **> mem.reallocarray(edges, current_length, sizeof(size_t*))

                            edges[one] = <size_t *> mem.allocarray(2 * len_edge_list, sizeof(size_t))

                        # Set up face_iter.atom_repr
                        face_iter.set_atom_repr()

                        # Copy the information.
                        edges[one][2*two] = face_iter.atom_repr[0]
                        edges[one][2*two + 1] = face_iter.atom_repr[1]
                        counter += 1

                    d = face_iter.next_face()  # Go to next face.

            # Success, copy the data to ``CombinatorialPolyhedron``.
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
        Calculate the ridges of the Polyhedron.

        If ``dual`` is ``True``, calculate the ridges of the polar.

        See :meth:`edges` and :meth:`ridges`.
        """
        if (self._edges is not NULL and dual) or (self._ridges is not NULL and not dual):
            return 0  # There is no need to recalculate.

        cdef MemoryAllocator mem = MemoryAllocator()
        cdef FaceIterator face_iter = self._face_iter(dual)
        cdef size_t len_ridge_list = self._length_edges_list
        cdef int dim = self.dimension()

        # For each ridge we determine its location in ``ridges``
        # by ``ridges[one][two]``.
        cdef size_t **ridges = <size_t**> mem.malloc(sizeof(size_t*))
        cdef size_t one, two

        cdef size_t counter = 0         # the number of ridges so far
        cdef size_t current_length = 1  # dynamically enlarge **ridges

        if dim == 1 and self._nr_facets > 1:
            # In this case there is a ridge, but its not a proper face.
            ridges[0] = <size_t *> mem.allocarray(2, sizeof(size_t))
            ridges[0][0] = 0
            ridges[0][1] = 1
            counter = 1

            # Success, copy the data to ``CombinatorialPolyhedron``.
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
            # :meth:`FaceIterator.set_request_dimension`
            # requires the dimension of the original Polyhedron
            face_iter.set_request_dimension(1)
        else:
            face_iter.set_request_dimension(dim - 2)

        if self._nr_facets > 1 and dim > 0:
            # If not, there won't even be any ridges
            # as intersection of two distince facets.
            # Prevent error message.

            while (face_iter.next_face() == dim - 2):

                # Determine the position in ``ridges``.
                one = counter // len_ridge_list
                two = counter % len_ridge_list

                # Enlarge ``ridges`` if needed.
                if unlikely(two == 0):
                    if unlikely(one + 1 > current_length):
                        # enlarge **ridges
                        current_length *= 2
                        ridges = <size_t **> mem.reallocarray(ridges, current_length, sizeof(size_t*))

                    ridges[one] = <size_t *> mem.allocarray(2 * len_ridge_list, sizeof(size_t))

                # Set up face_iter.coatom_repr
                face_iter.set_coatom_repr()

                # Copy the information.
                ridges[one][2*two] = face_iter.coatom_repr[0]
                ridges[one][2*two + 1] = face_iter.coatom_repr[1]
                counter += 1

        # Success, copy the data to ``CombinatorialPolyhedron``.
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
        Calculate all incidences for the face lattice.

        See :meth:`face_lattice`.
        """
        if self._face_lattice_incidences:
            return 1  # There is no need to recalculate the incidences.

        cdef size_t len_incidence_list = self._length_edges_list
        cdef int dim = self.dimension()
        f_vector = self.f_vector()
        cdef MemoryAllocator mem = MemoryAllocator()
        self._record_all_faces()  # set up ``self._all_faces``
        cdef ListOfAllFaces all_faces = self._all_faces

        # ``all_faces`` will store its incidences in ``first`` and ``second``.
        cdef size_t first = 0, second = 0

        # ``dimension_one`` and ``dimension_two`` will be the dimensions of the
        # incidences, we currently obtain from ``all_faces``.
        # Almost always ``dimension_two = dimension_one - 1``.
        cdef int dimension_one, dimension_two
        cdef int j  # an index for ``range(dimension_two + 1)``

        # The indices of the incidences in ``all_faces`` are levelwise.
        # Hence, we have to add to each index dependent on dimension:

        # For ``dimension_two`` we add:
        cdef size_t already_seen       # = sum(f_vector[j] for j in range(dimension_two + 1))

        # For ``dimension_one`` we add:
        cdef size_t already_seen_next  # = sum(f_vector[j] for j in range(dimension_two + 2))

        # For each incidence we determine its location in ``incidences``
        # by ``incidences[one][two]``.
        cdef size_t **incidences = <size_t**> mem.malloc(sizeof(size_t*))
        cdef size_t one, two

        cdef size_t counter = 0         # the number of incidences so far
        cdef size_t current_length = 1  # dynamically enlarge **incidences

        if all_faces is None:
            raise ValueError("could not determine a list of all faces")

        dimension_one = 0
        if dim > -1:
            while (f_vector[dimension_one + 1] == 0):
                # Taking care of cases, where there might be no faces
                # of dimension 0, 1, etc (``self._nr_lines > 0``).
                dimension_one += 1
            dimension_two = -1

        while (dimension_one < dim + 1):
            already_seen = sum(f_vector[j] for j in range(dimension_two + 1))
            already_seen_next = already_seen + f_vector[dimension_two + 1]

            if all_faces.dual:
                # If ``dual``, then ``all_faces`` has the dimensions reversed.
                all_faces.incidence_init(dim - 1 - dimension_two, dim - 1 - dimension_one)
            else:
                all_faces.incidence_init(dimension_one, dimension_two)

            # Get all incidences for fixed ``[dimension_one, dimension_two]``.
            while all_faces.next_incidence(&second, &first):

                # Determine the position in ``incidences``.
                one = counter // len_incidence_list
                two = counter % len_incidence_list

                # Enlarge ``incidences`` if needed.
                if unlikely(two == 0):
                    if unlikely(one + 1 > current_length):
                        # enlarge **incidences
                        current_length *= 2
                        incidences = <size_t **> mem.reallocarray(incidences, current_length, sizeof(size_t*))

                    incidences[one] = <size_t *> mem.allocarray(2 * len_incidence_list, sizeof(size_t))

                if all_faces.dual:
                    # If ``dual``, then ``second`` and ``first are flipped.
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

            # Increase dimensions.
            dimension_one += 1
            dimension_two = dimension_one - 1

        # Success, copy the data to ``CombinatorialPolyhedron``.
        self._nr_face_lattice_incidences = counter
        sig_block()
        self._mem_tuple += (mem,)
        self._face_lattice_incidences = incidences
        sig_unblock()

    def _record_all_faces(self):
        r"""
        Initialize :class:`ListOfAllFaces` for the Polyhedron.

        Record and sort all faces of the Polyhedron in that class.

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
            return  # Have recorded all faces already.

        self._all_faces = ListOfAllFaces(self)
        if self._all_faces is None:
            raise ValueError("could not determine a list of all faces")
