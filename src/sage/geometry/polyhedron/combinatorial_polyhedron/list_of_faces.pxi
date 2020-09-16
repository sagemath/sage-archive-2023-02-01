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
from cysignals.signals                cimport sig_on, sig_off

#############################################################################
# Face List Initalization
#############################################################################

cdef inline int face_list_init(face_list_t faces, size_t n_faces, size_t n_atoms, size_t n_coatoms, MemoryAllocator mem) except -1:
    """
    Sets the initial values for a list of faces with given number of faces
    and number of atoms.
    """
    face_list_shallow_init(faces, n_faces, n_atoms, n_coatoms, mem)
    cdef size_t i
    for i in range(n_faces):
        face_init(faces.faces[i], n_atoms, n_coatoms, mem)

cdef inline int face_list_shallow_init(face_list_t faces, size_t n_faces, size_t n_atoms, size_t n_coatoms, MemoryAllocator mem) except -1:
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

cdef inline int face_list_copy(face_list_t dst, face_list_t src) except -1:
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

cdef inline int face_list_shallow_copy(face_list_t dst, face_list_t src) except -1:
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

cdef inline int add_face_shallow(face_list_t faces, face_t face) nogil except -1:
    """
    Add a face to faces.
    """
    if not faces.total_n_faces >= faces.n_faces + 1:
        with gil:
            raise AssertionError
    faces.faces[faces.n_faces][0] = face[0]
    faces.n_faces += 1

cdef inline int add_face_deep(face_list_t faces, face_t face) except -1:
    """
    Add a face to faces.
    """
    assert faces.total_n_faces >= faces.n_faces + 1
    face_copy(faces.faces[faces.n_faces], face)
    faces.n_faces += 1

#############################################################################
# Face Comparison
#############################################################################

cdef void sort_faces_list(face_list_t faces):
    r"""
    Sorts faces in place.
    """
    cdef MemoryAllocator mem = MemoryAllocator()

    # Merge sort needs a second list of pointers.
    cdef face_t* extra_mem = <face_t*> mem.allocarray(faces.n_faces, sizeof(face_t))

    # Sort the faces using merge sort.
    _sort_faces_loop(faces.faces, faces.faces, extra_mem, faces.n_faces)
    cdef size_t i
    if faces.n_faces == 0:
        # Nothing to do. Prevent it from crashing.
        return

cdef void _sort_faces_loop(face_t* inp, face_t* out, face_t* extra_mem, size_t n_faces):
    """
    This is merge sort.

    Sorts ``inp`` and returns it in ``out``.

    .. WARNING::

        ``inp`` is the same as ``out`` or ``extra_mem``.

    See :func:`sort_faces`.
    """
    if n_faces == 0:
        # Prevent it from crashing.
        # In this case there is nothing to do anyway.
        return

    if n_faces == 1:
        # The final case, where there is only one element.
        out[0][0] = inp[0][0]
        return

    cdef size_t middle = n_faces//2
    cdef size_t len_upper_half = n_faces - middle

    # Sort the upper and lower half of ``inp`` iteratively into ``output2``.
    _sort_faces_loop(inp, extra_mem, out, middle)
    _sort_faces_loop(inp+middle, extra_mem+middle,
                     out+middle, len_upper_half)

    # Merge lower and upper half into ``output1``.
    cdef size_t i = 0        # index through lower half
    cdef size_t j = middle   # index through upper half
    cdef size_t counter = 0  # counts how many elements have been "merged" already
    cdef int val
    while i < middle and j < n_faces:
        # Compare the lowest elements of lower and upper half.
        val = face_cmp(extra_mem[i], extra_mem[j])
        if val < 0:
            out[counter][0] = extra_mem[i][0]
            i += 1
            counter += 1
        else:
            out[counter][0] = extra_mem[j][0]
            j += 1
            counter += 1
    if i < middle:
        # Add the remaining elements of lower half.
        while i < middle:
            out[counter][0] = extra_mem[i][0]
            i += 1
            counter += 1
    else:
        # Add the remaining elements of upper half.
        while j < n_faces:
            out[counter][0] = extra_mem[j][0]
            j += 1
            counter += 1

cdef inline size_t find_face(face_t face, face_list_t faces):
    r"""
    Return the index of ``face`` in ``faces``.

    Return ``-1`` if the ``face`` is not contained.

    .. NOTE::

        Assumes that ``faces`` are sorted.
    """
    cdef size_t start = 0
    cdef size_t middle
    cdef size_t n_faces = faces.n_faces
    cdef face_t* faces_pt = faces.faces
    cdef int val


    while (n_faces > 1):
        # In each iteration step, we will look for ``face`` in
        # ``faces_pt[start:start+n_faces]``.
        middle = n_faces//2
        val = face_cmp(face, faces_pt[middle + start])
        if val < 0:
            # If face is in the list, then in the lower half.
            # Look for face in ``faces[start : start + middle]`` in next step.
            n_faces = middle
        elif val > 0:
            # If face is in the list, then in the upper half.
            # Look for face in ``faces[start+middle:start+n_faces]``, i.e.
            # ``faces[start + middle : (start + middle) + n_faces - middle]``.
            n_faces -= middle
            start += middle
        else:
            return middle + start
    if face_cmp(face, faces_pt[start]) == 0:
        return start
    else:
        return -1

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
        with gil:
            raise AssertionError
    if not dest.n_atoms >= A.n_atoms:
        with gil:
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
    memset(is_not_new_face, 0, n_faces*sizeof(bint))

    # Step 1:
    n_faces -= 1
    faces.n_faces -= 1;
    face_list_intersection_fused(new_faces, faces, faces.faces[n_faces], algorithm)

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

cdef inline size_t get_next_level(
        face_list_t faces,
        face_list_t new_faces,
        face_list_t visited_all) nogil except -1:

    cdef size_t output
    sig_on()
    if faces.polyhedron_is_simple:
        output = get_next_level_fused(faces, new_faces, visited_all, <simple> 0)
    else:
        output = get_next_level_fused(faces, new_faces, visited_all, <standard> 0)
    sig_off()
    return output

cdef inline size_t bit_rep_to_coatom_rep(face_t face, face_list_t coatoms, size_t *output):
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
