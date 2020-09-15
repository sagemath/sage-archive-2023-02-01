"""
Cython methods for lists of faces.
"""
# ****************************************************************************
#       Copyright (C) 2020 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

include "sage/geometry/polyhedron/combinatorial_polyhedron/face.pxi"

from sage.geometry.polyhedron.combinatorial_polyhedron.list_of_faces cimport *
from libc.string                      cimport memset

#############################################################################
# Face List Initalization
#############################################################################

cdef inline void face_list_init(face_list_t faces, size_t n_faces, size_t n_atoms, size_t n_coatoms, MemoryAllocator mem):
    """
    Sets the initial values for a list of faces with given number of faces
    and number of atoms.
    """
    face_list_shallow_init(faces, n_faces, n_atoms, n_coatoms, mem)
    cdef size_t i
    for i in range(n_faces):
        face_init(faces.faces[i], n_atoms, n_coatoms, mem)

cdef inline void face_list_shallow_init(face_list_t faces, size_t n_faces, size_t n_atoms, size_t n_coatoms, MemoryAllocator mem):
    """
    Initialize ``faces`` completely, but only set up memory for the pointers to the faces.
    """
    faces.n_faces = n_faces
    faces.total_n_faces = n_faces
    faces.n_atoms = n_atoms
    faces.n_coatoms = n_coatoms
    faces.faces = <face_t *> mem.allocarray(n_faces, sizeof(face_t))
    faces.is_not_new_face = <bint *> mem.allocarray(n_faces, sizeof(bint))
    faces.polyhedron_is_simple = False

cdef inline int faces_copy(face_list_t dst, face_list_t src) except -1:
    """
    This is a deep copy. All the data for the faces is copied.

    Asserts that ``dst`` is allocated and fits everything.
    """
    assert dst.total_n_faces >= src.n_faces
    assert dst.n_atoms >= src.n_atoms
    assert dst.n_coatoms >= src.n_coatoms

    dst.n_faces = src.n_faces
    dst.polyhedron_is_simple = src.polyhedron_is_simple

    cdef size_t i
    for i in range(src.n_faces):
        face_copy(dst.faces[i], src.faces[i])

cdef inline int faces_shallow_copy(face_list_t dst, face_list_t src) except -1:
    """
    Copy the pointers to the faces.

    Asserts that ``dst`` contains enough space for the pointers.
    """
    assert dst.total_n_faces >= src.n_faces
    dst.n_atoms = src.n_atoms
    dst.n_coatoms = src.n_coatoms
    dst.polyhedron_is_simple = src.polyhedron_is_simple
    dst.n_faces = src.n_faces

    cdef size_t i
    for i in range(src.n_faces):
        dst.faces[i] = src.faces[i]

#############################################################################
# Face Comparison
#############################################################################

# Todo: Find

# Todo Sort

cdef inline bint is_contained_in_one_fused(face_t face, face_list_t faces, algorithm_variant algorithm) nogil:
    """
    Return whether ``face`` is contained in one of ``faces``.
    """
    cdef size_t i
    for i in range(faces.n_faces):
        if face_issubset_fused(face, faces.faces[i], algorithm):
            return True
    return False

cdef inline bint is_not_maximal_fused(face_list_t faces, size_t j, algorithm_variant algorithm) nogil:
    """
    Return whether face ``j`` is not maximal in ``faces``.
    """
    cdef size_t i
    if algorithm_variant is standard:
        for i in range(j):
            if face_issubset_fused(faces.faces[j], faces.faces[i], algorithm):
                return True
        for i in range(j+1, faces.n_faces):
            if face_issubset_fused(faces.faces[j], faces.faces[i], algorithm):
                return True
        return False
    else:
        # For simple polytopes an intersection of facets is of codimension 2,
        # if and only if it contains a coatom.
        return face_isempty(faces.faces[j])

#############################################################################
# Arithmetic
#############################################################################

cdef inline int face_list_intersection_fused(face_list_t dest, face_list_t A, face_t b, algorithm_variant algorithm) nogil except -1:
    """
    Set ``dest`` to be the intersection of each face of ``A`` with ``b``.
    """
    if not dest.total_n_faces >= A.n_faces:
        raise AssertionError
    if not dest.n_atoms >= A.n_atoms:
        raise AssertionError
    dest.n_faces = A.n_faces
    dest.polyhedron_is_simple = A.polyhedron_is_simple

    cdef size_t i
    for i in range(A.n_faces):
        face_intersection_fused(dest.faces[i], A.faces[i], b, algorithm)


cdef inline size_t get_next_level_fused(
        face_list_t faces,
        face_list_t new_faces,
        face_list_t visited_all, algorithm_variant algorithm) nogil except -1:
    """
    Set ``new_faces`` to be the facets of ``faces.faces[face.n_faces-1]``
    that are not contained in a face of ``visited_all``.

    Reduce the number of faces in ``faces`` by one.

    INPUT:

    - ``faces`` -- containing at least one face
    - ``new_faces`` -- needs to be of same size as ``faces``
    - ``visited_all`` -- the faces which have been visited before

    OUTPUT:

    - set ``new_faces`` to point to the new faces

    ALGORITHM:

    To get all facets of ``faces.faces[faces.n_faces-1]``, we would have to:
    - Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    - Add all the intersection of ``visited_all`` with the last face
    - Out of both the inclusion-maximal ones are of codimension one, i.e. facets.

    As we have visited all faces of ``visited_all``, we alter the algorithm
    to not revisit:
    Step 1: Intersect the first ``n_faces-1`` faces of ``faces`` with the last face.
    Step 2: Out of those the inclusion-maximal ones are some of the facets.
            At least we obtain all of those, that we have not already visited.
            Maybe, we get some more.
    Step 3: Only keep those that we have not already visited.
            We obtain exactly the facets of ``faces[n_faces-1]`` that we have
            not visited yet.
    """
    # We keep track, which face in ``new_faces`` is a new face.
    cdef size_t n_faces = faces.n_faces
    cdef bint* is_not_new_face = new_faces.is_not_new_face
    memset(is_not_new_face, 0, n_faces)

    # Step 1:
    n_faces -= 1
    faces.n_faces -= 1;
    face_list_intersection_fused(new_faces, faces, faces.faces[n_faces-1], algorithm)

    cdef size_t j
    for j in range(n_faces):
        if (is_not_maximal_fused(new_faces, j, algorithm) or  # Step 2
                is_contained_in_one_fused(new_faces.faces[j], visited_all, algorithm)):  # Step 3
            is_not_new_face[j] = True

    # Set ``new_faces`` to point to the correct ones.
    cdef size_t n_new_faces = 0
    for j in range(n_faces):
        if is_not_new_face[j]:
            continue
        # It is a new face of codimension 1.
        # Either ``faces.n_new_faces == j`` or ``new_faces.faces[n_new_faces]`` is not
        # a new face.

        swap_faces(new_faces.faces[j], new_faces.faces[n_new_faces])

        n_new_faces += 1

    new_faces.n_faces = n_new_faces
    return n_new_faces

cdef inline size_t get_next_level_1(
        face_list_t faces,
        face_list_t new_faces,
        face_list_t visited_all) nogil except -1:
    if faces.polyhedron_is_simple:
        return get_next_level_fused(faces, new_faces, visited_all, <simple> 0)
    else:
        return get_next_level_fused(faces, new_faces, visited_all, <standard> 0)

cdef inline size_t bit_rep_to_coatom_rep_1(face_t face, face_list_t coatoms, size_t *output):
    """
    Write the coatom-representation of face in output. Return length.
    ``face_length`` is the length of ``face`` and ``coatoms[i]``
    in terms of uint64_t.
    ``n_coatoms`` length of ``coatoms``.
    """
    cdef size_t count_length = 0
    cdef size_t i
    for i in range(coatoms.n_faces):
        if face_issubset(face, coatoms.faces[i]):
            output[count_length] = i
            count_length += 1
    return count_length

cdef inline bint face_list_check_alignment(face_list_t faces):
    """
    Return whether all faces in ``faces`` are aligned correctly.
    """
    cdef size_t i
    for i in range(faces.n_faces):
       if not face_check_alignment(faces.faces[i]):
           return False
    return True
