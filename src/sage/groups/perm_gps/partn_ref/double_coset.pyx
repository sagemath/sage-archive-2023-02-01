r"""
Double cosets

This module implements a general algorithm for computing double coset problems
for pairs of objects. The class of objects in question must be some kind
of structure for which an isomorphism is a permutation in $S_n$ for some $n$,
which we call here the order of the object. Given objects $X$ and $Y$,
the program returns an isomorphism in list permutation form if $X \cong Y$, and
a NULL pointer otherwise.

In order to take advantage of the algorithms in this module for a specific kind
of object, one must implement (in Cython) three functions which will be specific
to the kind of objects in question. Pointers to these functions are passed to
the main function of the module, which is \code{double_coset}. For specific
examples of implementations of these functions, see any of the files in
\code{sage.groups.perm_gps.partn_ref} beginning with "refinement." They are:

A. \code{refine_and_return_invariant}:

    Signature:

    \code{int refine_and_return_invariant(PartitionStack *PS, object S, int *cells_to_refine_by, int ctrb_len)}

    This function should split up cells in the partition at the top of the
    partition stack in such a way that any automorphism that respects the
    partition also respects the resulting partition. The array
    cells_to_refine_by is a list of the beginning positions of some cells which
    have been changed since the last refinement. It is not necessary to use
    this in an implementation of this function, but it will affect performance.
    One should consult \code{refinement_graphs} for more details and ideas for
    particular implementations.

    Output:

    An integer $I$ invariant under the orbits of $S_n$.  That is, if
    $\gamma \in S_n$, then
    $$ I(G, PS, cells_to_refine_by) = I( \gamma(G), \gamma(PS), \gamma(cells_to_refine_by) ) .$$


B. \code{compare_structures}:

    Signature:

    \code{int compare_structures(int *gamma_1, int *gamma_2, object S)}

    This function must implement a total ordering on the set of objects of fixed
    order. Return -1 if \code{gamma_1(S) < gamma_2(S)}, 0 if
    \code{gamma_1(S) == gamma_2(S)}, 1 if \code{gamma_1(S) > gamma_2(S)}.

C. \code{all_children_are_equivalent}:

    Signature:

    \code{bint all_children_are_equivalent(PartitionStack *PS, object S)}

    This function must return False unless it is the case that each discrete
    partition finer than the top of the partition stack is equivalent to the
    others under some automorphism of S. The converse need not hold: if this is
    indeed the case, it still may return False. This function is originally used
    as a consequence of Lemma 2.25 in [1].

DOCTEST:
    sage: import sage.groups.perm_gps.partn_ref.double_coset

REFERENCE:

    [1] McKay, Brendan D. Practical Graph Isomorphism. Congressus Numerantium,
        Vol. 30 (1981), pp. 45-87.

    [2] Leon, Jeffrey. Permutation Group Algorithms Based on Partitions, I:
        Theory and Algorithms. J. Symbolic Computation, Vol. 12 (1991), pp.
        533-583.

"""

#*****************************************************************************
#      Copyright (C) 2006 - 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include 'data_structures_pyx.pxi' # includes bitsets

cdef inline int cmp(int a, int b):
    if a < b: return -1
    elif a == b: return 0
    else: return 1

cdef inline bint stacks_are_equivalent(PartitionStack *PS1, PartitionStack *PS2):
    cdef int i, j, depth = min(PS1.depth, PS2.depth)
    for i from 0 <= i < PS1.degree:
        j = cmp(PS1.levels[i], PS2.levels[i])
        if j == 0: continue
        if ( (j < 0 and PS1.levels[i] <= depth and PS2.levels[i] > depth)
            or (j > 0 and PS2.levels[i] <= depth and PS1.levels[i] > depth) ):
            return 0
    return 1

# Functions

cdef int *double_coset(object S1, object S2, int **partition1, int *ordering2,
    int n, bint (*all_children_are_equivalent)(PartitionStack *PS, object S),
    int (*refine_and_return_invariant)\
         (PartitionStack *PS, object S, int *cells_to_refine_by, int ctrb_len),
    int (*compare_structures)(int *gamma_1, int *gamma_2, object S1, object S2) ):
    """
    Traverse the search space for double coset calculation.

    INPUT:
    S1, S2 -- pointers to the structures
    partition1 -- array representing a partition of the points of S1
    ordering2 -- an ordering of the points of S2 representing a second partition
    n -- the number of points (points are assumed to be 0,1,...,n-1)
    all_children_are_equivalent -- pointer to a function
        INPUT:
        PS -- pointer to a partition stack
        S -- pointer to the structure
        OUTPUT:
        bint -- returns True if it can be determined that all refinements below
            the current one will result in an equivalent discrete partition
    refine_and_return_invariant -- pointer to a function
        INPUT:
        PS -- pointer to a partition stack
        S -- pointer to the structure
        alpha -- an array consisting of numbers, which indicate the starting
            positions of the cells to refine against (will likely be modified)
        OUTPUT:
        int -- returns an invariant under application of arbitrary permutations
    compare_structures -- pointer to a function
        INPUT:
        gamma_1, gamma_2 -- (list) permutations of the points of S1 and S2
        S1, S2 -- pointers to the structures
        OUTPUT:
        int -- 0 if gamma_1(S1) = gamma_2(S2), otherwise -1 or 1 (see docs for cmp),
            such that the set of all structures is well-ordered
    OUTPUT:
    If S1 and S2 are isomorphic, a pointer to an integer array representing an
    isomorphism. Otherwise, a NULL pointer.

    """
    cdef PartitionStack *current_ps, *first_ps, *left_ps
    cdef int first_meets_current = -1
    cdef int current_kids_are_same = 1
    cdef int first_kids_are_same

    cdef int *indicators

    cdef OrbitPartition *orbits_of_subgroup
    cdef int subgroup_primary_orbit_size = 0
    cdef int minimal_in_primary_orbit

    cdef bitset_t *fixed_points_of_generators # i.e. fp
    cdef bitset_t *minimal_cell_reps_of_generators # i.e. mcr
    cdef int len_of_fp_and_mcr = 100
    cdef int index_in_fp_and_mcr = -1

    cdef bitset_t *vertices_to_split
    cdef bitset_t vertices_have_been_reduced
    cdef int *permutation, *id_perm, *cells_to_refine_by
    cdef int *vertices_determining_current_stack

    cdef int *output

    cdef int i, j, k
    cdef bint discrete, automorphism, update_label
    cdef bint backtrack, new_vertex, narrow, mem_err = 0

    if n == 0:
        return NULL

    indicators = <int *> sage_malloc(n * sizeof(int))

    fixed_points_of_generators = <bitset_t *> sage_malloc( len_of_fp_and_mcr * sizeof(bitset_t) )
    minimal_cell_reps_of_generators = <bitset_t *> sage_malloc( len_of_fp_and_mcr * sizeof(bitset_t) )

    vertices_to_split = <bitset_t *> sage_malloc( n * sizeof(bitset_t) )
    permutation = <int *> sage_malloc( n * sizeof(int) )
    id_perm = <int *> sage_malloc( n * sizeof(int) )
    cells_to_refine_by = <int *> sage_malloc( n * sizeof(int) )
    vertices_determining_current_stack = <int *> sage_malloc( n * sizeof(int) )

    output = <int *> sage_malloc( n * sizeof(int) )

    current_ps = PS_new(n, 0)
    left_ps = PS_new(n, 0)
    orbits_of_subgroup = OP_new(n)

    # Check for allocation failures:
    if indicators                         is NULL or \
       fixed_points_of_generators         is NULL or \
       minimal_cell_reps_of_generators    is NULL or \
       vertices_to_split                  is NULL or \
       permutation                        is NULL or \
       id_perm                            is NULL or \
       cells_to_refine_by                 is NULL or \
       vertices_determining_current_stack is NULL or \
       current_ps                         is NULL or \
       orbits_of_subgroup                 is NULL or \
       output                             is NULL:
        if indicators                         is not NULL:
            sage_free(indicators)
        if fixed_points_of_generators         is not NULL:
            sage_free(fixed_points_of_generators)
        if minimal_cell_reps_of_generators    is not NULL:
            sage_free(minimal_cell_reps_of_generators)
        if vertices_to_split                  is not NULL:
            sage_free(vertices_to_split)
        if permutation                        is not NULL:
            sage_free(permutation)
        if id_perm                            is not NULL:
            sage_free(id_perm)
        if cells_to_refine_by                 is not NULL:
            sage_free(cells_to_refine_by)
        if vertices_determining_current_stack is not NULL:
            sage_free(vertices_determining_current_stack)
        if current_ps is not NULL:
            PS_dealloc(current_ps)
        if orbits_of_subgroup is not NULL:
            OP_dealloc(orbits_of_subgroup)
        if output is not NULL:
            sage_free(output)
        raise MemoryError

    # Initialize bitsets, checking for allocation failures:
    cdef bint succeeded = 1
    for i from 0 <= i < len_of_fp_and_mcr:
        try:
            bitset_init(fixed_points_of_generators[i], n)
        except MemoryError:
            succeeded = 0
            for j from 0 <= j < i:
                bitset_clear(fixed_points_of_generators[j])
                bitset_clear(minimal_cell_reps_of_generators[j])
            break
        try:
            bitset_init(minimal_cell_reps_of_generators[i], n)
        except MemoryError:
            succeeded = 0
            for j from 0 <= j < i:
                bitset_clear(fixed_points_of_generators[j])
                bitset_clear(minimal_cell_reps_of_generators[j])
            bitset_clear(fixed_points_of_generators[i])
            break
    if succeeded:
        for i from 0 <= i < n:
            try:
                bitset_init(vertices_to_split[i], n)
            except MemoryError:
                succeeded = 0
                for j from 0 <= j < i:
                    bitset_clear(vertices_to_split[j])
                for j from 0 <= j < len_of_fp_and_mcr:
                    bitset_clear(fixed_points_of_generators[j])
                    bitset_clear(minimal_cell_reps_of_generators[j])
                break
    if succeeded:
        try:
            bitset_init(vertices_have_been_reduced, n)
        except MemoryError:
            succeeded = 0
            for j from 0 <= j < n:
                bitset_clear(vertices_to_split[j])
            for j from 0 <= j < len_of_fp_and_mcr:
                bitset_clear(fixed_points_of_generators[j])
                bitset_clear(minimal_cell_reps_of_generators[j])
    if not succeeded:
        sage_free(indicators)
        sage_free(fixed_points_of_generators)
        sage_free(minimal_cell_reps_of_generators)
        sage_free(vertices_to_split)
        sage_free(permutation)
        sage_free(id_perm)
        sage_free(cells_to_refine_by)
        sage_free(vertices_determining_current_stack)
        PS_dealloc(current_ps)
        PS_dealloc(left_ps)
        OP_dealloc(orbits_of_subgroup)
        raise MemoryError

    bitset_zero(vertices_have_been_reduced)

    # set up the identity permutation
    for i from 0 <= i < n:
        id_perm[i] = i
    if ordering2 is NULL:
        ordering2 = id_perm

    # Copy data of partition to left_ps, and
    # ordering of that partition to current_ps.
    i = 0
    j = 0
    while partition1[i] is not NULL:
        k = 0
        while partition1[i][k] != -1:
            left_ps.entries[j+k] = partition1[i][k]
            left_ps.levels[j+k] = n
            current_ps.entries[j+k] = ordering2[j+k]
            current_ps.levels[j+k] = n
            k += 1
        left_ps.levels[j+k-1] = 0
        current_ps.levels[j+k-1] = 0
        PS_move_min_to_front(current_ps, j, j+k-1)
        j += k
        i += 1
    left_ps.levels[j-1] = -1
    current_ps.levels[j-1] = -1
    left_ps.depth = 0
    current_ps.depth = 0
    left_ps.degree = n
    current_ps.degree = n

    # default values of "infinity"
    for i from 0 <= i < n:
        indicators[i] = -1

    cdef bint possible = 1
    cdef bint unknown = 1

    # Our first refinement needs to check every cell of the partition,
    # so cells_to_refine_by needs to be a list of pointers to each cell.
    j = 1
    cells_to_refine_by[0] = 0
    for i from 0 < i < n:
        if left_ps.levels[i-1] == 0:
            cells_to_refine_by[j] = i
            j += 1

    k = refine_and_return_invariant(left_ps, S1, cells_to_refine_by, j)

    j = 1
    cells_to_refine_by[0] = 0
    for i from 0 < i < n:
        if current_ps.levels[i-1] == 0:
            cells_to_refine_by[j] = i
            j += 1
    j = refine_and_return_invariant(current_ps, S2, cells_to_refine_by, j)
    if k != j:
        possible = 0; unknown = 0
    elif not stacks_are_equivalent(left_ps, current_ps):
        possible = 0; unknown = 0
    else:
        PS_move_all_mins_to_front(current_ps)

    first_ps = NULL
    # Refine down to a discrete partition
    while not PS_is_discrete(left_ps) and possible:
        k = PS_first_smallest(left_ps, vertices_to_split[left_ps.depth])
        i = left_ps.depth
        indicators[i] = split_point_and_refine(left_ps, k, S1, refine_and_return_invariant, cells_to_refine_by)
        vertices_determining_current_stack[current_ps.depth] = PS_first_smallest(current_ps, vertices_to_split[current_ps.depth])
        bitset_unset(vertices_have_been_reduced, current_ps.depth)
        while 1:
            j =  split_point_and_refine(current_ps, vertices_determining_current_stack[i], S2, refine_and_return_invariant, cells_to_refine_by)
            if indicators[i] != j:
                possible = 0
            elif not stacks_are_equivalent(left_ps, current_ps):
                possible = 0
            else:
                PS_move_all_mins_to_front(current_ps)
            if not possible:
                j = vertices_determining_current_stack[i] + 1
                j = bitset_next(vertices_to_split[i], j)
                if j == -1:
                    break
                else:
                    possible = 1
                    vertices_determining_current_stack[i] = j
                    current_ps.depth -= 1 # reset for next refinement
            else: break
        if not all_children_are_equivalent(current_ps, S2):
            current_kids_are_same = current_ps.depth + 1

    if possible:
        if compare_structures(left_ps.entries, current_ps.entries, S1, S2) == 0:
            unknown = 0
        else:
            first_meets_current = current_ps.depth
            first_kids_are_same = current_ps.depth
            first_ps = PS_copy(current_ps)
            current_ps.depth -= 1

    # Main loop:
    while possible and unknown and current_ps.depth != -1:

        # I. Search for a new vertex to split, and update subgroup information
        new_vertex = 0
        if current_ps.depth > first_meets_current:
            # If we are not at a node of the first stack, reduce size of
            # vertices_to_split by using the symmetries we already know.
            if not bitset_check(vertices_have_been_reduced, current_ps.depth):
                for i from 0 <= i <= index_in_fp_and_mcr:
                    j = 0
                    while j < current_ps.depth and bitset_check(fixed_points_of_generators[i], vertices_determining_current_stack[j]):
                        j += 1
                    # If each vertex split so far is fixed by generator i,
                    # then remove elements of vertices_to_split which are
                    # not minimal in their orbits under generator i.
                    if j == current_ps.depth:
                        for k from 0 <= k < n:
                            if bitset_check(vertices_to_split[current_ps.depth], k) and not bitset_check(minimal_cell_reps_of_generators[i], k):
                                bitset_flip(vertices_to_split[current_ps.depth], k)
                bitset_flip(vertices_have_been_reduced, current_ps.depth)
            # Look for a new point to split.
            i = vertices_determining_current_stack[current_ps.depth] + 1
            i = bitset_next(vertices_to_split[current_ps.depth], i)
            if i != -1:
                # There is a new point.
                vertices_determining_current_stack[current_ps.depth] = i
                new_vertex = 1
            else:
                # No new point: backtrack.
                current_ps.depth -= 1
        else:
            # If we are at a node of the first stack, the above reduction
            # will not help. Also, we must update information about
            # primary orbits here.
            if current_ps.depth < first_meets_current:
                # If we are done searching under this part of the first
                # stack, then first_meets_current is one higher, and we
                # are looking at a new primary orbit (corresponding to a
                # larger subgroup in the stabilizer chain).
                first_meets_current = current_ps.depth
                for i from 0 <= i < n:
                    if bitset_check(vertices_to_split[current_ps.depth], i):
                        minimal_in_primary_orbit = i
                        break
            while 1:
                i = vertices_determining_current_stack[current_ps.depth]
                # This was the last point to be split here.
                # If it is in the same orbit as minimal_in_primary_orbit,
                # then count it as an element of the primary orbit.
                if OP_find(orbits_of_subgroup, i) == OP_find(orbits_of_subgroup, minimal_in_primary_orbit):
                    subgroup_primary_orbit_size += 1
                # Look for a new point to split.
                i += 1
                i = bitset_next(vertices_to_split[current_ps.depth], i)
                if i != -1:
                    # There is a new point.
                    vertices_determining_current_stack[current_ps.depth] = i
                    if orbits_of_subgroup.mcr[OP_find(orbits_of_subgroup, i)] == i:
                        new_vertex = 1
                        break
                else:
                    # No new point: backtrack.
                    # Note that now, we are backtracking up the first stack.
                    vertices_determining_current_stack[current_ps.depth] = -1
                    # If every choice of point to split gave something in the
                    # (same!) primary orbit, then all children of the first
                    # stack at this point are equivalent.
                    j = 0
                    for i from 0 <= i < n:
                        if bitset_check(vertices_to_split[current_ps.depth], i):
                            j += 1
                    if j == subgroup_primary_orbit_size and first_kids_are_same == current_ps.depth+1:
                        first_kids_are_same = current_ps.depth
                    # Backtrack.
                    subgroup_primary_orbit_size = 0
                    current_ps.depth -= 1
                    break
        if not new_vertex:
            continue

        if current_kids_are_same > current_ps.depth + 1:
            current_kids_are_same = current_ps.depth + 1

        # II. Refine down to a discrete partition, or until
        # we leave the part of the tree we are interested in
        discrete = 0
        while 1:
            i = current_ps.depth
            while 1:
                k = split_point_and_refine(current_ps, vertices_determining_current_stack[i], S2, refine_and_return_invariant, cells_to_refine_by)
                PS_move_all_mins_to_front(current_ps)
                if indicators[i] != k:
                    possible = 0
                elif not stacks_are_equivalent(left_ps, current_ps):
                    possible = 0
                if not possible:
                    j = vertices_determining_current_stack[i] + 1
                    j = bitset_next(vertices_to_split[i], j)
                    if j == -1:
                        break
                    else:
                        possible = 1
                        vertices_determining_current_stack[i] = j
                        current_ps.depth -= 1 # reset for next refinement
                else: break
            if not possible:
                break
            if PS_is_discrete(current_ps):
                break
            vertices_determining_current_stack[current_ps.depth] = PS_first_smallest(current_ps, vertices_to_split[current_ps.depth])
            bitset_unset(vertices_have_been_reduced, current_ps.depth)
            if not all_children_are_equivalent(current_ps, S2):
                current_kids_are_same = current_ps.depth + 1

        # III. Check for automorphisms and isomorphisms
        automorphism = 0
        if possible:
            PS_get_perm_from(current_ps, first_ps, permutation)
            if compare_structures(permutation, id_perm, S2, S2) == 0:
                automorphism = 1
        if not automorphism and possible:
            # if we get here, discrete must be true
            if current_ps.depth != left_ps.depth:
                possible = 0
            elif compare_structures(left_ps.entries, current_ps.entries, S1, S2) == 0:
                unknown = 0
                break # main loop
            else:
                possible = 0
        if automorphism:
            if index_in_fp_and_mcr < len_of_fp_and_mcr - 1:
                index_in_fp_and_mcr += 1
            bitset_zero(fixed_points_of_generators[index_in_fp_and_mcr])
            bitset_zero(minimal_cell_reps_of_generators[index_in_fp_and_mcr])
            for i from 0 <= i < n:
                if permutation[i] == i:
                    bitset_set(fixed_points_of_generators[index_in_fp_and_mcr], i)
                    bitset_set(minimal_cell_reps_of_generators[index_in_fp_and_mcr], i)
                else:
                    bitset_unset(fixed_points_of_generators[index_in_fp_and_mcr], i)
                    k = i
                    j = permutation[i]
                    while j != i:
                        if j < k: k = j
                        j = permutation[j]
                    if k == i:
                        bitset_set(minimal_cell_reps_of_generators[index_in_fp_and_mcr], i)
                    else:
                        bitset_unset(minimal_cell_reps_of_generators[index_in_fp_and_mcr], i)
            current_ps.depth = first_meets_current
            if OP_merge_list_perm(orbits_of_subgroup, permutation): # if permutation made orbits coarser
                if orbits_of_subgroup.mcr[OP_find(orbits_of_subgroup, minimal_in_primary_orbit)] != minimal_in_primary_orbit:
                    continue # main loop
            if bitset_check(vertices_have_been_reduced, current_ps.depth):
                bitset_and(vertices_to_split[current_ps.depth], vertices_to_split[current_ps.depth], minimal_cell_reps_of_generators[index_in_fp_and_mcr])
            continue # main loop
        if not possible:
            possible = 1
            i = current_ps.depth
            current_ps.depth = min(first_kids_are_same-1, current_kids_are_same-1)
            if i == current_kids_are_same:
                continue # main loop
            if index_in_fp_and_mcr < len_of_fp_and_mcr - 1:
                index_in_fp_and_mcr += 1
            bitset_zero(fixed_points_of_generators[index_in_fp_and_mcr])
            bitset_zero(minimal_cell_reps_of_generators[index_in_fp_and_mcr])
            j = current_ps.depth
            current_ps.depth = i # just for mcr and fixed functions...
            for i from 0 <= i < n:
                if PS_is_mcr(current_ps, i):
                    bitset_set(minimal_cell_reps_of_generators[index_in_fp_and_mcr], i)
                    if PS_is_fixed(current_ps, i):
                        bitset_set(fixed_points_of_generators[index_in_fp_and_mcr], i)
            current_ps.depth = j
            if bitset_check(vertices_have_been_reduced, current_ps.depth):
                bitset_and(vertices_to_split[current_ps.depth], vertices_to_split[current_ps.depth], minimal_cell_reps_of_generators[index_in_fp_and_mcr])

    # End of main loop.

    if possible and not unknown:
        for i from 0 <= i < n:
            output[left_ps.entries[i]] = current_ps.entries[i]
    else:
        sage_free(output)
        output = NULL

    # Deallocate:
    for i from 0 <= i < len_of_fp_and_mcr:
        bitset_clear(fixed_points_of_generators[i])
        bitset_clear(minimal_cell_reps_of_generators[i])
    for i from 0 <= i < n:
        bitset_clear(vertices_to_split[i])
    bitset_clear(vertices_have_been_reduced)
    sage_free(indicators)
    sage_free(fixed_points_of_generators)
    sage_free(minimal_cell_reps_of_generators)
    sage_free(vertices_to_split)
    sage_free(permutation)
    sage_free(id_perm)
    sage_free(cells_to_refine_by)
    sage_free(vertices_determining_current_stack)
    PS_dealloc(current_ps)
    PS_dealloc(left_ps)
    OP_dealloc(orbits_of_subgroup)
    if first_ps is not NULL: PS_dealloc(first_ps)

    return output





