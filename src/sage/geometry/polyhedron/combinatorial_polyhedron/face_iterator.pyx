r"""
FaceIterator

This module provides a face iterator for polyhedra.

This iterator in principle works on every graded lattice, where
every interval of length two has exactly 4 elements (diamond property).

It also works on unbounded polyhedra, as those satisfy the diamond property,
except for intervals including the empty face.

Terminology in this module:

- Vertices              -- ``[vertices, rays, lines]`` of the polyhedron.
- Facets                -- facets of the polyhedron.
- Coatoms               -- the faces from which all others are constructed in
                           the face iterator. This will be facets or vertices.
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

EXAMPLES:

Construct a face iterator::

    sage: from sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator \
    ....:         import FaceIterator
    sage: P = polytopes.octahedron()
    sage: C = CombinatorialPolyhedron(P)

    sage: FaceIterator(C, False)
    Iterator over the proper faces of a polyhedron of dimension 3
    sage: FaceIterator(C, False, output_dimension=2)
    Iterator over the 2-faces of a polyhedron of dimension 3

Iterator in the non-dual mode starts with facets.
By default the dimension of the current face is given::

    sage: it = FaceIterator(C, False)
    sage: next(it)
    2

Iterator in the dual-mode start with vertices::

    sage: it = FaceIterator(C, True)
    sage: next(it)
    0

Obtain the vertex-representation::

    sage: it = FaceIterator(C, False)
    sage: next(it)
    2
    sage: it.vertex_repr()
    (A vertex at (0, -1, 0), A vertex at (0, 0, -1), A vertex at (1, 0, 0))

    sage: it.length_vertex_repr()
    3

Obtain the facet-representation::

    sage: it = FaceIterator(C, True)
    sage: next(it)
    0
    sage: it.facet_repr()
    (An inequality (-1, -1, 1) x + 1 >= 0,
     An inequality (-1, -1, -1) x + 1 >= 0,
      An inequality (-1, 1, -1) x + 1 >= 0,
       An inequality (-1, 1, 1) x + 1 >= 0)
    sage: it.facet_repr(names=False)
    (4, 5, 6, 7)

    sage: it.length_facet_repr()
    4

In non-dual mode one can ignore all faces contained in the current face::

    sage: it = FaceIterator(C, False)
    sage: next(it)
    2
    sage: it.facet_repr(names=False)
    (7,)
    sage: it.ignore_subfaces()
    sage: [it.facet_repr(names=False) for _ in it]
    [(6,),
    (5,),
    (4,),
    (3,),
    (2,),
    (1,),
    (0,),
    (5, 6),
    (1, 6),
    (0, 1, 5, 6),
    (4, 5),
    (0, 5),
    (0, 3, 4, 5),
    (3, 4),
    (2, 3),
    (0, 3),
    (0, 1, 2, 3),
    (1, 2),
    (0, 1)]

In dual mode one can ignore all faces that contain the current face::

    sage: it = FaceIterator(C, True)
    sage: next(it)
    0
    sage: it.vertex_repr(names=False)
    (5,)
    sage: it.ignore_supfaces()
    sage: [it.vertex_repr(names=False) for _ in it]
    [(4,),
    (3,),
    (2,),
    (1,),
    (0,),
    (3, 4),
    (2, 4),
    (0, 4),
    (0, 3, 4),
    (0, 2, 4),
    (1, 3),
    (0, 3),
    (0, 1, 3),
    (1, 2),
    (0, 2),
    (0, 1, 2),
    (0, 1)]

AUTHOR:

- Jonathan Kliem (2019-04)
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

from sage.rings.integer     cimport smallInteger
from cysignals.signals      cimport sig_check, sig_on, sig_off
from .conversions           cimport bit_repr_to_vertex_list
from .base                  cimport CombinatorialPolyhedron
from .bit_vector_operations cimport get_next_level, count_atoms, bit_repr_to_coatom_repr

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef class FaceIterator(SageObject):
    r"""
    A class to iterate over all faces of a polyhedron.

    Constructs all proper from the facets. In dual mode, constructs all proper
    faces from the vertices. Dual will be faster for less vertices than facets.

    INPUT:

    - ``C`` -- a :class:`CombinatorialPolyhedron`
    - ``dual`` -- if ``True``, then dual polyhedron is used for iteration
      (only possible for bounded Polyhedra)
    - ``output_dimension`` -- if not ``None``, then the FaceIterator will only yield
      faces of this dimension

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
        sage: it.current_face_dimension()
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

    Construct a FaceIterator only yielding dimension `2` faces::

        sage: P = polytopes.permutahedron(5)
        sage: C = CombinatorialPolyhedron(P)
        sage: it = C.face_iter(dimension=2)
        sage: counter = 0
        sage: for _ in it: counter += 1
        sage: print ('permutahedron(5) has', counter,
        ....:        'faces of dimension 2')
        permutahedron(5) has 150 faces of dimension 2
        sage: C.f_vector()
        (1, 120, 240, 150, 30, 1)

    ALGORITHM:

    The algorithm to visit all proper faces exactly once is roughly
    equivalent to::

        faces = [set(facet) for facet in P.facets()]
        face_iterator(faces, [])

        def face_iterator(faces, visited_all):
            # Visit all faces of a polyhedron `P`, except those contained in
            # any of the visited all.

            # Assumes ``faces`` to be excactly those facets of `P`
            # that are not contained in any of the ``visited_all``.

            # Assumes ``visited_all`` to be some list of faces of
            # a polyhedron `P_2`, which contains `P` as one of its faces.

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
    def __init__(self, CombinatorialPolyhedron C, bint dual, output_dimension=None):
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

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator).run()
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
        self._mem = MemoryAllocator()

        # We will not yield the empty face.
        # If there are `n` lines, than there
        # are no faces below dimension `n`.
        # The dimension of the level-sets in the face lattice jumps from `n` to `-1`.
        self.lowest_dimension = self.nr_lines

        if output_dimension is not None:
            if not output_dimension in range(0,self.dimension):
                raise ValueError("``output_dimension`` must be the dimension of proper faces")
            if self.dual:
                # In dual mode, the dimensions are reversed.
                self.output_dimension = self.dimension - 1 - output_dimension
            else:
                self.output_dimension = output_dimension
            self.lowest_dimension = max(self.nr_lines, self.output_dimension)
        else:
            self.output_dimension = -2

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
            Iterator over the proper faces of a polyhedron of dimension 3

            sage: C.face_iter(1)
            Iterator over the 1-faces of a polyhedron of dimension 3
        """
        if self.output_dimension != -2:
            if self.dual:
                # ouput_dimension is stored with respect to the dual
                intended_dimension = self.dimension - 1 - self.output_dimension
            else:
                intended_dimension = self.output_dimension
            output = "Iterator over the {}-faces".format(intended_dimension)
        else:
            output = "Iterator over the proper faces"
        return output + " of a polyhedron of dimension {}".format(self.dimension)

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
        return smallInteger(self.dual*(self.dimension-1-d) + (1-self.dual)*d)

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

    def current_face_dimension(self):
        r"""
        Return the dimension of the current face.

        EXAMPLES::

            sage: P = polytopes.associahedron(['A', 3])
            sage: C = CombinatorialPolyhedron(P)
            sage: it = C.face_iter()
            sage: next(it)
            2
            sage: it.current_face_dimension()
            2
            sage: all(d == it.current_face_dimension() for d in it)
            True
        """
        if unlikely(self.face is NULL):
            raise ValueError("iterator not set to a face yet")
        if unlikely(self.current_dimension == self.dimension):
            raise ValueError("iterator consumed")
        # If ``dual == 0`` return current dimension,
        # if ``dual == 1`` translate current dimension to dual and then return.
        return smallInteger(self.dual*(self.dimension-1-self.current_dimension) +
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
                return tuple(smallInteger(self.coatom_repr[i])
                             for i in range(length))
        else:
            # if not dual, the vertex-represention corresponds to the atom-representation
            length = self.set_atom_repr()
            if names and self._V:
                return tuple(self._V[self.atom_repr[i]]
                             for i in range(length))
            else:
                return tuple(smallInteger(self.atom_repr[i])
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
            return smallInteger(self.set_coatom_repr())
        else:
            return smallInteger(self.length_atom_repr())

    def facet_repr(self, names=True):
        r"""
        Return the facet-representation of the current face.

        The facet-representation consists of the facets
        that contain the face and of the equalities of the polyhedron.

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
                return tuple(smallInteger(self.coatom_repr[i])
                             for i in range(length))
        else:
            # if dual, the facet-represention corresponds to the atom-representation
            length = self.set_atom_repr()  # fill self.atom_repr_face
            if names and self._H:
                return tuple(self._H[self.atom_repr[i]]
                             for i in range(length)) + self._equalities
            else:
                return tuple(smallInteger(self.atom_repr[i])
                             for i in range(length))

    def length_facet_repr(self):
        r"""
        Returns the length of the :meth:`facet_repr`.

        Might be faster than ``len(self.facet_repr())``.

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
            return smallInteger(self.set_coatom_repr())
        else:
            return smallInteger(self.length_atom_repr())

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

        Will return the dimension of the polyhedron on failure.

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

        if (self.output_dimension > -2) and (self.output_dimension != self.current_dimension):
            # If only a specific dimension was requested (i.e. ``self.output_dimension > -2``),
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
        Compute the number of atoms in the current face by counting the
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
