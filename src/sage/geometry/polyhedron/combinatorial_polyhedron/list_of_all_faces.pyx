# distutils: language = c++

from __future__ import absolute_import, division, print_function
from .list_of_faces \
        import facets_tuple_to_bit_repr_of_facets, \
               facets_tuple_to_bit_repr_of_vertices

from sage.rings.integer cimport smallInteger
from libc.string        cimport memcmp, memcpy, memset
from .list_of_faces     cimport vertex_list_to_bit_repr, bit_repr_to_vertex_list
from .base              cimport CombinatorialPolyhedron
from .face_iterator     cimport FaceIterator

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef extern from "bit_vector_operations.cc":
    cdef void intersection(uint64_t *A, uint64_t *B, uint64_t *C,
                           size_t face_length)
#    Return ``A & ~B == 0``.
#    A is not subset of B, iff there is a vertex in A, which is not in B.
#    ``face_length`` is the length of A and B in terms of uint64_t.

    cdef size_t bit_repr_to_coatom_repr(
            uint64_t *face, uint64_t **coatoms, size_t nr_coatoms,
            size_t face_length, size_t *output)
#        Write the coatom-representation of face in output. Return length.
#        ``face_length`` is the length of ``face`` and ``coatoms[i]``
#        in terms of uint64_t.
#        ``nr_coatoms`` length of ``coatoms``.

cdef class ListOfAllFaces:
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
                return tuple(smallInteger(self.coatom_repr[i])
                             for i in range(length))
        else:
            # if not dual, the vertex-represention corresponds to the atom-representation
            length = self.set_atom_repr(dimension, index)
            if names and self._V:
                return tuple(self._V[self.atom_repr[i]]
                             for i in range(length))
            else:
                return tuple(smallInteger(self.atom_repr[i])
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
                return tuple(smallInteger(self.coatom_repr[i])
                             for i in range(length))
        else:
            # if dual, the facet-represention corresponds to the atom-representation
            dimension = self.dimension - 1 - dimension  # if dual, the dimensions are reversed
            length = self.set_atom_repr(dimension, index)
            if names and self._H:
                return tuple(self._H[self.atom_repr[i]]
                             for i in range(length)) + self._equalities
            else:
                return tuple(smallInteger(self.atom_repr[i])
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
