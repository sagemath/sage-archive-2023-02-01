r"""
Automorphism groups and canonical labels

This module implements a general algorithm for computing automorphism groups and
canonical labels of objects. The class of objects in question must be some kind
of structure for which an automorphism is a permutation in $S_n$ for some $n$,
which we call here the order of the object. It should be noted that the word
"canonical" in the term canonical label is used loosely. Given an object $X$,
the canonical label $C(X)$ is just an object for which $C(X) = C(Y)$ for any
$Y$ such that $X \cong Y$.

To understand the algorithm better, a few definitions will help. First one
should note that throughout this module, $S_n$ is defined as the set of
permutations of the set $[n] := \{0, 1, ..., n-1\}$. One speaks of
ordered partitions $\Pi = (C_1, ..., C_k)$ of $[n]$. The $C_i$ are disjoint
subsets of $[n]$ whose union is equal to $[n]$, and we call them cells. A
partition $\Pi_1$ is said to be finer than another partition $\Pi_2$ if every
cell of $\Pi_1$ is a subset of some cell of $\Pi_2$ (in this situation we also
say that $\Pi_2$ is coarser than $\Pi_1$). A partition stack
$\nu = (\Pi_0, ..., \Pi_k)$ is a sequence of partitions such that $\Pi_i$ is
finer than $\Pi_{i-1}$ for each $i$ such that $1 \leq i \leq k$. Sometimes these
are called nested partitions. The depth of $\nu$ is defined to be $k$ and the
degree of $\nu$ is defined to be $n$. Partition stacks are implemented as
\code{PartitionStack} (see automorphism_group_canonical_label.pxd for more
information). Finally, we say that a permutation respects the partition $\Pi$ if
its orbits form a partition which is finer than $\Pi$.

In order to take advantage of the algorithms in this module for a specific kind
of object, one must implement (in Cython) three functions which will be specific
to the kind of objects in question. Pointers to these functions are passed to
the main function of the module, which is \code{get_aut_gp_and_can_lab}. For
specific examples of implementations of these functions, see any of the files in
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

    \code{int compare_structures(int *gamma_1, int *gamma_2, object S1, object S2)}

    This function must implement a total ordering on the set of objects of fixed
    order. Return -1 if \code{gamma_1(S1) < gamma_2(S2)}, 0 if
    \code{gamma_1(S1) == gamma_2(S2)}, 1 if \code{gamma_1(S1) > gamma_2(S2)}.

C. \code{all_children_are_equivalent}:

    Signature:

    \code{bint all_children_are_equivalent(PartitionStack *PS, object S)}

    This function must return False unless it is the case that each discrete
    partition finer than the top of the partition stack is equivalent to the
    others under some automorphism of S. The converse need not hold: if this is
    indeed the case, it still may return False. This function is originally used
    as a consequence of Lemma 2.25 in [1].

DOCTEST:
    sage: import sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label

REFERENCE:

    [1] McKay, Brendan D. Practical Graph Isomorphism. Congressus Numerantium,
        Vol. 30 (1981), pp. 45-87.

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

# Functions

cdef bint all_children_are_equivalent_trivial(PartitionStack *PS, object S):
    return 0

cdef int refine_and_return_invariant_trivial(PartitionStack *PS, object S, int *cells_to_refine_by, int ctrb_len):
    return 0

cdef int compare_structures_trivial(int *gamma_1, int *gamma_2, object S1, object S2):
    return 0

def test_get_aut_gp_and_can_lab_trivially(int n=6, partition=[[0,1,2],[3,4],[5]],
    canonical_label=True, base=False):
    """
    sage: tttt = sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label.test_get_aut_gp_and_can_lab_trivially
    sage: tttt()
    12
    sage: tttt(canonical_label=False, base=False)
    12
    sage: tttt(canonical_label=False, base=True)
    12
    sage: tttt(canonical_label=True, base=True)
    12
    sage: tttt(n=0, partition=[])
    1
    sage: tttt(n=0, partition=[], canonical_label=False, base=False)
    1
    sage: tttt(n=0, partition=[], canonical_label=False, base=True)
    1
    sage: tttt(n=0, partition=[], canonical_label=True, base=True)
    1

    """
    cdef int i, j
    cdef aut_gp_and_can_lab_return *output
    cdef Integer I = Integer(0)
    cdef int **part = <int **> sage_malloc((len(partition)+1) * sizeof(int *))
    for i from 0 <= i < len(partition):
        part[i] = <int *> sage_malloc((len(partition[i])+1) * sizeof(int))
        for j from 0 <= j < len(partition[i]):
            part[i][j] = partition[i][j]
        part[i][len(partition[i])] = -1
    part[len(partition)] = NULL
    output = get_aut_gp_and_can_lab([], part, n, &all_children_are_equivalent_trivial, &refine_and_return_invariant_trivial, &compare_structures_trivial, canonical_label, base, True)
    mpz_set(I.value, output.order)
    print I
    for i from 0 <= i < len(partition):
        sage_free(part[i])
    sage_free(part)
    mpz_clear(output.order)
    sage_free(output.generators)
    if base:
        sage_free(output.base)
    if canonical_label:
        sage_free(output.relabeling)
    sage_free(output)

cdef aut_gp_and_can_lab_return *get_aut_gp_and_can_lab(object S, int **partition, int n,
    bint (*all_children_are_equivalent)(PartitionStack *PS, object S),
    int (*refine_and_return_invariant)\
         (PartitionStack *PS, object S, int *cells_to_refine_by, int ctrb_len),
    int (*compare_structures)(int *gamma_1, int *gamma_2, object S1, object S2),
    bint canonical_label, bint base, bint order):
    """
    Traverse the search space for subgroup/canonical label calculation.

    INPUT:
    S -- pointer to the structure
    partition -- array representing a partition of the points
    len_partition -- length of the partition
    n -- the number of points (points are assumed to be 0,1,...,n-1)
    canonical_label -- whether to search for canonical label - if True, return
        the permutation taking S to its canonical label
    base -- if True, return the base for which the given generating set is a
        strong generating set
    order -- if True, return the order of the automorphism group
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
    pointer to a aut_gp_and_can_lab_return struct

    """
    cdef PartitionStack *current_ps, *first_ps, *label_ps
    cdef int first_meets_current = -1
    cdef int label_meets_current
    cdef int current_kids_are_same = 1
    cdef int first_kids_are_same

    cdef int *current_indicators, *first_indicators, *label_indicators
    cdef int first_and_current_indicator_same
    cdef int label_and_current_indicator_same = -1
    cdef int compared_current_and_label_indicators

    cdef OrbitPartition *orbits_of_subgroup
    cdef mpz_t subgroup_size
    cdef int subgroup_primary_orbit_size = 0
    cdef int minimal_in_primary_orbit
    cdef int size_of_generator_array = 2*n*n

    cdef bitset_t *fixed_points_of_generators # i.e. fp
    cdef bitset_t *minimal_cell_reps_of_generators # i.e. mcr
    cdef int len_of_fp_and_mcr = 100
    cdef int index_in_fp_and_mcr = -1

    cdef bitset_t *vertices_to_split
    cdef bitset_t vertices_have_been_reduced
    cdef int *permutation, *id_perm, *cells_to_refine_by
    cdef int *vertices_determining_current_stack

    cdef int i, j, k
    cdef bint discrete, automorphism, update_label
    cdef bint backtrack, new_vertex, narrow, mem_err = 0

    cdef aut_gp_and_can_lab_return *output = <aut_gp_and_can_lab_return *> sage_malloc(sizeof(aut_gp_and_can_lab_return))
    if output is NULL:
        raise MemoryError

    if n == 0:
        output.generators = NULL
        output.num_gens = 0
        output.relabeling = NULL
        output.base = NULL
        output.base_size = 0
        mpz_init_set_si(output.order, 1)
        return output

    # Allocate:
    mpz_init_set_si(subgroup_size, 1)

    output.generators = <int *> sage_malloc( size_of_generator_array * sizeof(int) )
    output.num_gens = 0
    if base:
        output.base = <int *> sage_malloc( n * sizeof(int) )
        output.base_size = 0

    current_indicators = <int *> sage_malloc( n * sizeof(int) )
    first_indicators = <int *> sage_malloc( n * sizeof(int) )

    if canonical_label:
        label_indicators = <int *> sage_malloc( n * sizeof(int) )
        output.relabeling = <int *> sage_malloc( n * sizeof(int) )

    fixed_points_of_generators = <bitset_t *> sage_malloc( len_of_fp_and_mcr * sizeof(bitset_t) )
    minimal_cell_reps_of_generators = <bitset_t *> sage_malloc( len_of_fp_and_mcr * sizeof(bitset_t) )

    vertices_to_split = <bitset_t *> sage_malloc( n * sizeof(bitset_t) )
    permutation = <int *> sage_malloc( n * sizeof(int) )
    id_perm = <int *> sage_malloc( n * sizeof(int) )
    cells_to_refine_by = <int *> sage_malloc( n * sizeof(int) )
    vertices_determining_current_stack = <int *> sage_malloc( n * sizeof(int) )

    current_ps = PS_new(n, 0)
    orbits_of_subgroup = OP_new(n)

    # Check for allocation failures:
    cdef bint succeeded = 1
    if canonical_label:
        if label_indicators is NULL or output.relabeling is NULL:
            succeeded = 0
            if label_indicators is not NULL:
                sage_free(label_indicators)
            if output.relabeling is not NULL:
                sage_free(output.relabeling)
    if base:
        if output.base is NULL:
            succeeded = 0
    if not succeeded                              or \
       output.generators                  is NULL or \
       current_indicators                 is NULL or \
       first_indicators                   is NULL or \
       fixed_points_of_generators         is NULL or \
       minimal_cell_reps_of_generators    is NULL or \
       vertices_to_split                  is NULL or \
       permutation                        is NULL or \
       id_perm                            is NULL or \
       cells_to_refine_by                 is NULL or \
       vertices_determining_current_stack is NULL or \
       current_ps                         is NULL or \
       orbits_of_subgroup                 is NULL:
        if output.generators                  is not NULL:
            sage_free(output.generators)
        if base:
            if output.base                        is not NULL:
                sage_free(output.base)
        if current_indicators                 is not NULL:
            sage_free(current_indicators)
        if first_indicators                   is not NULL:
            sage_free(first_indicators)
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
        sage_free(output)
        mpz_clear(subgroup_size)
        raise MemoryError

    # Initialize bitsets, checking for allocation failures:
    succeeded = 1
    for i from 0 <= i < len_of_fp_and_mcr:
        try:
            bitset_init(fixed_points_of_generators[i], n)
        except MemoryError:
            succeeded = 0
            for j from 0 <= j < i:
                bitset_free(fixed_points_of_generators[j])
                bitset_free(minimal_cell_reps_of_generators[j])
            break
        try:
            bitset_init(minimal_cell_reps_of_generators[i], n)
        except MemoryError:
            succeeded = 0
            for j from 0 <= j < i:
                bitset_free(fixed_points_of_generators[j])
                bitset_free(minimal_cell_reps_of_generators[j])
            bitset_free(fixed_points_of_generators[i])
            break
    if succeeded:
        for i from 0 <= i < n:
            try:
                bitset_init(vertices_to_split[i], n)
            except MemoryError:
                succeeded = 0
                for j from 0 <= j < i:
                    bitset_free(vertices_to_split[j])
                for j from 0 <= j < len_of_fp_and_mcr:
                    bitset_free(fixed_points_of_generators[j])
                    bitset_free(minimal_cell_reps_of_generators[j])
                break
    if succeeded:
        try:
            bitset_init(vertices_have_been_reduced, n)
        except MemoryError:
            succeeded = 0
            for j from 0 <= j < n:
                bitset_free(vertices_to_split[j])
            for j from 0 <= j < len_of_fp_and_mcr:
                bitset_free(fixed_points_of_generators[j])
                bitset_free(minimal_cell_reps_of_generators[j])
    if not succeeded:
        sage_free(output.generators)
        if base:
            sage_free(output.base)
        sage_free(current_indicators)
        sage_free(first_indicators)
        if canonical_label:
            sage_free(output.relabeling)
            sage_free(label_indicators)
        sage_free(fixed_points_of_generators)
        sage_free(minimal_cell_reps_of_generators)
        sage_free(vertices_to_split)
        sage_free(permutation)
        sage_free(id_perm)
        sage_free(cells_to_refine_by)
        sage_free(vertices_determining_current_stack)
        PS_dealloc(current_ps)
        sage_free(output)
        mpz_clear(subgroup_size)
        raise MemoryError

    bitset_zero(vertices_have_been_reduced)

    # Copy data of partition to current_ps.
    i = 0
    j = 0
    while partition[i] is not NULL:
        k = 0
        while partition[i][k] != -1:
            current_ps.entries[j+k] = partition[i][k]
            current_ps.levels[j+k] = n
            k += 1
        current_ps.levels[j+k-1] = 0
        PS_move_min_to_front(current_ps, j, j+k-1)
        j += k
        i += 1
    current_ps.levels[j-1] = -1
    current_ps.depth = 0
    current_ps.degree = n

    # default values of "infinity"
    for i from 0 <= i < n:
        first_indicators[i] = -1
    if canonical_label:
        for i from 0 <= i < n:
            label_indicators[i] = -1

    # set up the identity permutation
    for i from 0 <= i < n:
        id_perm[i] = i

    # Our first refinement needs to check every cell of the partition,
    # so cells_to_refine_by needs to be a list of pointers to each cell.
    j = 1
    cells_to_refine_by[0] = 0
    for i from 0 < i < n:
        if current_ps.levels[i-1] == 0:
            cells_to_refine_by[j] = i
            j += 1
    # Ignore the invariant this time, since we are
    # creating the root of the search tree.
    refine_and_return_invariant(current_ps, S, cells_to_refine_by, j)
    PS_move_all_mins_to_front(current_ps)

    # Refine down to a discrete partition
    while not PS_is_discrete(current_ps):
        if not all_children_are_equivalent(current_ps, S):
            current_kids_are_same = current_ps.depth + 1
        vertices_determining_current_stack[current_ps.depth] = PS_first_smallest(current_ps, vertices_to_split[current_ps.depth])
        bitset_unset(vertices_have_been_reduced, current_ps.depth)
        i = current_ps.depth
        current_indicators[i] = split_point_and_refine(current_ps, vertices_determining_current_stack[i], S, refine_and_return_invariant, cells_to_refine_by)
        PS_move_all_mins_to_front(current_ps)
        first_indicators[i] = current_indicators[i]
        if canonical_label:
            label_indicators[i] = current_indicators[i]
        if base:
            output.base[i] = vertices_determining_current_stack[i]
            output.base_size += 1

    first_meets_current = current_ps.depth
    first_kids_are_same = current_ps.depth
    first_and_current_indicator_same = current_ps.depth
    first_ps = PS_copy(current_ps)
    if canonical_label:
        label_meets_current = current_ps.depth
        label_and_current_indicator_same = current_ps.depth
        compared_current_and_label_indicators = 0
        label_ps = PS_copy(current_ps)
    current_ps.depth -= 1

    # Main loop:
    while current_ps.depth != -1:

        # If necessary, update label_meets_current
        if canonical_label and label_meets_current > current_ps.depth:
            label_meets_current = current_ps.depth

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
                    # Update group size and backtrack.
                    mpz_mul_si(subgroup_size, subgroup_size, subgroup_primary_orbit_size)
                    subgroup_primary_orbit_size = 0
                    current_ps.depth -= 1
                    break
        if not new_vertex:
            continue

        if current_kids_are_same > current_ps.depth + 1:
            current_kids_are_same = current_ps.depth + 1
        if first_and_current_indicator_same > current_ps.depth:
            first_and_current_indicator_same = current_ps.depth
        if canonical_label and label_and_current_indicator_same >= current_ps.depth:
            label_and_current_indicator_same = current_ps.depth
            compared_current_and_label_indicators = 0

        # II. Refine down to a discrete partition, or until
        # we leave the part of the tree we are interested in
        discrete = 0
        while 1:
            i = current_ps.depth
            current_indicators[i] = split_point_and_refine(current_ps, vertices_determining_current_stack[i], S, refine_and_return_invariant, cells_to_refine_by)
            PS_move_all_mins_to_front(current_ps)
            if first_and_current_indicator_same == i and current_indicators[i] == first_indicators[i]:
                first_and_current_indicator_same += 1
            if canonical_label:
                if compared_current_and_label_indicators == 0:
                    if label_indicators[i] == -1:
                        compared_current_and_label_indicators = -1
                    else:
                        compared_current_and_label_indicators = cmp(current_indicators[i], label_indicators[i])
                if label_and_current_indicator_same == i and compared_current_and_label_indicators == 0:
                    label_and_current_indicator_same += 1
                if compared_current_and_label_indicators > 0:
                    label_indicators[i] = current_indicators[i]
            if first_and_current_indicator_same < current_ps.depth and (not canonical_label or compared_current_and_label_indicators < 0):
                break
            if PS_is_discrete(current_ps):
                discrete = 1
                break
            vertices_determining_current_stack[current_ps.depth] = PS_first_smallest(current_ps, vertices_to_split[current_ps.depth])
            bitset_unset(vertices_have_been_reduced, current_ps.depth)
            if not all_children_are_equivalent(current_ps, S):
                current_kids_are_same = current_ps.depth + 1

        # III. Check for automorphisms and labels
        automorphism = 0
        backtrack = 1 - discrete
        if discrete:
            if current_ps.depth == first_and_current_indicator_same:
                PS_get_perm_from(current_ps, first_ps, permutation)
                if compare_structures(permutation, id_perm, S, S) == 0:
                    automorphism = 1
        if not automorphism:
            if (not canonical_label) or (compared_current_and_label_indicators < 0):
                backtrack = 1
            else:
                # if we get here, discrete must be true
                update_label = 0
                if (compared_current_and_label_indicators > 0) or (current_ps.depth < label_ps.depth):
                    update_label = 1
                else:
                    i = compare_structures(current_ps.entries, label_ps.entries, S, S)
                    if i > 0:
                        update_label = 1
                    elif i < 0:
                        backtrack = 1
                    else:
                        PS_get_perm_from(current_ps, label_ps, permutation)
                        automorphism = 1
                if update_label:
                    PS_copy_from_to(current_ps, label_ps)
                    compared_current_and_label_indicators = 0
                    label_meets_current = current_ps.depth
                    label_and_current_indicator_same = current_ps.depth
                    label_indicators[current_ps.depth] = -1
                    backtrack = 1
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
            if OP_merge_list_perm(orbits_of_subgroup, permutation): # if permutation made orbits coarser
                if n*output.num_gens == size_of_generator_array:
                    # must double its size
                    size_of_generator_array *= 2
                    output.generators = <int *> sage_realloc( output.generators, size_of_generator_array )
                    if output.generators is NULL:
                        mem_err = True
                        break # main loop
                j = n*output.num_gens
                for i from 0 <= i < n:
                    output.generators[j + i] = permutation[i]
                output.num_gens += 1
                if orbits_of_subgroup.mcr[OP_find(orbits_of_subgroup, minimal_in_primary_orbit)] != minimal_in_primary_orbit:
                    current_ps.depth = first_meets_current
                    continue # main loop
            if canonical_label:
                current_ps.depth = label_meets_current
            else:
                current_ps.depth = first_meets_current
            if bitset_check(vertices_have_been_reduced, current_ps.depth):
                bitset_and(vertices_to_split[current_ps.depth], vertices_to_split[current_ps.depth], minimal_cell_reps_of_generators[index_in_fp_and_mcr])
            continue # main loop
        if backtrack:
            i = current_ps.depth
            current_ps.depth = min(max(first_kids_are_same - 1, label_and_current_indicator_same), current_kids_are_same - 1)
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

    mpz_init_set(output.order, subgroup_size)
    if canonical_label:
        for i from 0 <= i < n:
            output.relabeling[label_ps.entries[i]] = i

    # Deallocate:
    for i from 0 <= i < len_of_fp_and_mcr:
        bitset_free(fixed_points_of_generators[i])
        bitset_free(minimal_cell_reps_of_generators[i])
    for i from 0 <= i < n:
        bitset_free(vertices_to_split[i])
    bitset_free(vertices_have_been_reduced)
    sage_free(current_indicators)
    sage_free(first_indicators)
    if canonical_label:
        PS_dealloc(label_ps)
        sage_free(label_indicators)
    sage_free(fixed_points_of_generators)
    sage_free(minimal_cell_reps_of_generators)
    sage_free(vertices_to_split)
    sage_free(permutation)
    sage_free(id_perm)
    sage_free(cells_to_refine_by)
    sage_free(vertices_determining_current_stack)
    PS_dealloc(current_ps)
    PS_dealloc(first_ps)
    OP_dealloc(orbits_of_subgroup)
    mpz_clear(subgroup_size)
    if not mem_err:
        return output
    else:
        mpz_clear(output.order)
        sage_free(output.generators)
        if base:
            sage_free(output.base)
        if canonical_label:
            sage_free(output.relabeling)
        sage_free(output)
        output = NULL
        raise MemoryError














