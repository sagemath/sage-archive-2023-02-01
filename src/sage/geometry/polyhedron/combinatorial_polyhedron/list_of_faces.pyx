# distutils: language = c++

from __future__ import absolute_import, division
from sage.structure.element import is_Matrix

from libc.string cimport memset
from cysignals.signals cimport sig_check, sig_on, sig_off

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef extern from "bit_vector_operations.cc":
    # Any Bit-representation is assumed to be `chunksize`-Bit aligned.
    cdef const size_t chunksize

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
        ....: from sage.rings.integer cimport smallInteger
        ....:
        ....: def vertex_list_to_bit_repr_wrapper(tup):
        ....:     cdef size_t face_length = max(tup)//64 + 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = <uint64_t *> mem.allocarray(face_length, 8)
        ....:     vertex_list_to_bit_repr(tup, output, face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(face_length))
        ....:
        ....: def vertex_list_to_bit_repr_wrong_size(tup):
        ....:     cdef size_t face_length = 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = <uint64_t *> mem.allocarray(face_length, 8)
        ....:     vertex_list_to_bit_repr(tup, output, face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(face_length))
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
        ....: from sage.rings.integer cimport smallInteger
        ....:
        ....: def incidences_to_bit_reprs_wrapper(tup):
        ....:     cdef size_t face_length = (len(tup)-1)//64 + 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = \
        ....:          <uint64_t *> mem.allocarray(face_length, 8)
        ....:     incidences_to_bit_repr(tup, output, face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(face_length))
        ....:
        ....: def incidences_to_bit_reprs_wrong_size(tup):
        ....:     cdef size_t face_length = 1
        ....:     cdef MemoryAllocator mem = MemoryAllocator()
        ....:     cdef uint64_t *output = \
        ....:         <uint64_t *> mem.allocarray(face_length, 8)
        ....:     incidences_to_bit_repr(tup, output, face_length)
        ....:     return tuple(smallInteger(output[i]) for i in range(face_length))
        ....: ''')  # long time

        sage: incidences_to_bit_reprs_wrapper((0,) * 62 + (1,1))  # long time
        (3,)
        sage: incidences_to_bit_reprs_wrapper((0,) * 61 + (1,0,1) +  # long time
        ....:                                 (0,) * 61 + (1,))
        (5, 4)
        sage: incidences_to_bit_reprs_wrapper((1,) * 64)  # long time
        (-1,)
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
        ....: from sage.rings.integer cimport smallInteger
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
        ....:     return tuple(smallInteger(output[i]) for i in range(length))
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
        ....: from sage.rings.integer cimport smallInteger
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
        ....:     return tuple(smallInteger(output[i]) for i in range(length))
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
        ....: from sage.rings.integer cimport smallInteger
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
        ....:     return tuple(smallInteger(output[i]) for i in range(length))
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
        ....: from sage.rings.integer cimport smallInteger
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
        ....:     return tuple(smallInteger(output[i]) for i in range(length))
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
        ....: from sage.rings.integer cimport smallInteger
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
        ....:     return tuple(smallInteger(output[i]) for i in range(outputlength))
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

cdef class ListOfFaces:
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

    cpdef int calculate_dimension(self) except -2:
        r"""
        Calculate the dimension of a Polyhedron by its facets.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.base \
            ....:     import ListOfFaces, facets_tuple_to_bit_repr_of_facets, \
            ....:            facets_tuple_to_bit_repr_of_vertices
            sage: bi_pyr = ((0,1,4), (1,2,4), (2,3,4), (3,0,4),
            ....:           (0,1,5), (1,2,5), (2,3,5), (3,0,5))
            sage: facets = facets_tuple_to_bit_repr_of_facets(bi_pyr, 6)
            sage: vertices = facets_tuple_to_bit_repr_of_vertices(bi_pyr, 6)
            sage: facets.calculate_dimension()
            3
            sage: vertices.calculate_dimension()
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
            ....:     import ListOfFaces, incidence_matrix_to_bit_repr_of_facets, \
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
            ....:     d2 = facets.calculate_dimension()
            ....:     d3 = vertices.calculate_dimension()
            ....:     if not d1 == d2 == d3:
            ....:         print('calculation_dimension() seems to be incorrect')
        """
        if self.nr_faces == 0:
            raise TypeError("at least one face needed")
        return self.calculate_dimension_loop(self.data, self.nr_faces, self.face_length)

    cdef int calculate_dimension_loop(self, uint64_t **faces, size_t nr_faces,
                                      size_t face_length) except -2:
        r"""
        Calculate the dimension of a Polyhedron by its facets.

        INPUT:

        - ``faces`` -- facets in Bit-representation
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
        return self.calculate_dimension_loop(newfaces, new_nr_faces, face_length) + 1
