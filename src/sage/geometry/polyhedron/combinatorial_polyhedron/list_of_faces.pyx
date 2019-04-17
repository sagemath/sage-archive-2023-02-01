from __future__ import absolute_import, division
from sage.structure.element import is_Matrix

from cysignals.signals      cimport sig_on, sig_off
from .bit_vector_operations cimport chunksize, get_next_level, count_atoms

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

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

        sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
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
            ....: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
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

            sage: TestSuite(sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces.ListOfFaces).run()
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

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
            ....:     import ListOfFaces
            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import facets_tuple_to_bit_repr_of_facets, \
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

            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces \
            ....:     import ListOfFaces
            sage: from sage.geometry.polyhedron.combinatorial_polyhedron.conversions \
            ....:     import incidence_matrix_to_bit_repr_of_facets, \
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
