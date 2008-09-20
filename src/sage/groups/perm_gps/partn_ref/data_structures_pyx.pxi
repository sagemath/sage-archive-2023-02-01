r"""
Data structures

This module implements TODO...

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

include '../../../misc/bitset.pxi'

cdef inline int min(int a, int b):
    if a < b: return a
    else: return b

cdef inline int max(int a, int b):
    if a > b: return a
    else: return b

# OrbitPartitions

cdef inline OrbitPartition *OP_new(int n):
    """
    Allocate and return a pointer to a new OrbitPartition of degree n. Returns a
    null pointer in the case of an allocation failure.
    """
    cdef int i
    cdef OrbitPartition *OP = <OrbitPartition *> \
                                sage_malloc(sizeof(OrbitPartition))
    if OP is NULL:
        return OP
    OP.degree = n
    OP.parent = <int *> sage_malloc( n * sizeof(int) )
    OP.rank = <int *> sage_malloc( n * sizeof(int) )
    OP.mcr = <int *> sage_malloc( n * sizeof(int) )
    OP.size = <int *> sage_malloc( n * sizeof(int) )
    if OP.parent is NULL or OP.rank is NULL or \
       OP.mcr is NULL or OP.size is NULL:
        if OP.parent is not NULL: sage_free(OP.parent)
        if OP.rank is not NULL: sage_free(OP.rank)
        if OP.mcr is not NULL: sage_free(OP.mcr)
        if OP.size is not NULL: sage_free(OP.size)
        return NULL
    for i from 0 <= i < n:
        OP.parent[i] = i
        OP.rank[i] = 0
        OP.mcr[i] = i
        OP.size[i] = 1
    return OP

cdef inline int OP_dealloc(OrbitPartition *OP):
    sage_free(OP.parent)
    sage_free(OP.rank)
    sage_free(OP.mcr)
    sage_free(OP.size)
    sage_free(OP)
    return 0

cdef inline int OP_find(OrbitPartition *OP, int n):
    """
    Report the representative ("root") of the cell which contains n.
    """
    if OP.parent[n] == n:
        return n
    else:
        OP.parent[n] = OP_find(OP, OP.parent[n])
        return OP.parent[n]

cdef inline int OP_join(OrbitPartition *OP, int m, int n):
    """
    Join the cells containing m and n, if they are different.
    """
    cdef int m_root = OP_find(OP, m)
    cdef int n_root = OP_find(OP, n)
    if OP.rank[m_root] > OP.rank[n_root]:
        OP.parent[n_root] = m_root
        OP.mcr[m_root] = min(OP.mcr[m_root], OP.mcr[n_root])
        OP.size[m_root] += OP.size[n_root]
    elif OP.rank[m_root] < OP.rank[n_root]:
        OP.parent[m_root] = n_root
        OP.mcr[n_root] = min(OP.mcr[m_root], OP.mcr[n_root])
        OP.size[n_root] += OP.size[m_root]
    elif m_root != n_root:
        OP.parent[n_root] = m_root
        OP.mcr[m_root] = min(OP.mcr[m_root], OP.mcr[n_root])
        OP.size[m_root] += OP.size[n_root]
        OP.rank[m_root] += 1

cdef inline int OP_merge_list_perm(OrbitPartition *OP, int *gamma):
    """
    Joins the cells of OP which intersect the same orbit of gamma.

    INPUT:
        gamma - an integer array representing i -> gamma[i].

    OUTPUT:
        1 - something changed
        0 - orbits of gamma all contained in cells of OP
    """
    cdef int i, i_root, gamma_i_root, changed = 0
    for i from 0 <= i < OP.degree:
        if gamma[i] == i: continue
        i_root = OP_find(OP, i)
        gamma_i_root = OP_find(OP, gamma[i])
        if i_root != gamma_i_root:
            changed = 1
            OP_join(OP, i_root, gamma_i_root)
    return changed

def OP_represent(int n=9, merges=[(0,1),(2,3),(3,4)], perm=[1,2,0,4,3,6,7,5,8]):
    """
    Demonstration and testing.

    DOCTEST:
        sage: from sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label \
        ... import OP_represent
        sage: OP_represent()
        Allocating OrbitPartition...
        Allocation passed.
        Checking that each element reports itself as its root.
        Each element reports itself as its root.
        Merging:
        Merged 0 and 1.
        Merged 2 and 3.
        Merged 3 and 4.
        Done merging.
        Finding:
        0 -> 0, root: size=2, mcr=0, rank=1
        1 -> 0
        2 -> 2, root: size=3, mcr=2, rank=1
        3 -> 2
        4 -> 2
        5 -> 5, root: size=1, mcr=5, rank=0
        6 -> 6, root: size=1, mcr=6, rank=0
        7 -> 7, root: size=1, mcr=7, rank=0
        8 -> 8, root: size=1, mcr=8, rank=0
        Allocating array to test merge_perm.
        Allocation passed.
        Merging permutation: [1, 2, 0, 4, 3, 6, 7, 5, 8]
        Done merging.
        Finding:
        0 -> 0, root: size=5, mcr=0, rank=2
        1 -> 0
        2 -> 0
        3 -> 0
        4 -> 0
        5 -> 5, root: size=3, mcr=5, rank=1
        6 -> 5
        7 -> 5
        8 -> 8, root: size=1, mcr=8, rank=0
        Deallocating OrbitPartition.
        Done.

    """
    cdef int i
    print "Allocating OrbitPartition..."
    cdef OrbitPartition *OP = OP_new(n)
    if OP is NULL:
        print "Allocation failed!"
        return
    print "Allocation passed."
    print "Checking that each element reports itself as its root."
    good = True
    for i from 0 <= i < n:
        if not OP_find(OP, i) == i:
            print "Failed at i = %d!"%i
            good = False
    if good: print "Each element reports itself as its root."
    print "Merging:"
    for i,j in merges:
        OP_join(OP, i, j)
        print "Merged %d and %d."%(i,j)
    print "Done merging."
    print "Finding:"
    for i from 0 <= i < n:
        j = OP_find(OP, i)
        s = "%d -> %d"%(i, j)
        if i == j:
            s += ", root: size=%d, mcr=%d, rank=%d"%\
                   (OP.size[i], OP.mcr[i], OP.rank[i])
        print s
    print "Allocating array to test merge_perm."
    cdef int *gamma = <int *> sage_malloc( n * sizeof(int) )
    if gamma is NULL:
        print "Allocation failed!"
        OP_dealloc(OP)
        return
    print "Allocation passed."
    for i from 0 <= i < n:
        gamma[i] = perm[i]
    print "Merging permutation: %s"%perm
    OP_merge_list_perm(OP, gamma)
    print "Done merging."
    print "Finding:"
    for i from 0 <= i < n:
        j = OP_find(OP, i)
        s = "%d -> %d"%(i, j)
        if i == j:
            s += ", root: size=%d, mcr=%d, rank=%d"%\
                   (OP.size[i], OP.mcr[i], OP.rank[i])
        print s
    print "Deallocating OrbitPartition."
    sage_free(gamma)
    OP_dealloc(OP)
    print "Done."

# PartitionStacks

cdef inline PartitionStack *PS_new(int n, bint unit_partition):
    """
    Allocate and return a pointer to a new PartitionStack of degree n. Returns a
    null pointer in the case of an allocation failure.
    """
    cdef int i
    cdef PartitionStack *PS = <PartitionStack *> \
                                sage_malloc(sizeof(PartitionStack))
    if PS is NULL:
        return PS
    PS.entries = <int *> sage_malloc( n * sizeof(int) )
    PS.levels = <int *> sage_malloc( n * sizeof(int) )
    if PS.entries is NULL or PS.levels is NULL:
        if PS.entries is not NULL: sage_free(PS.entries)
        if PS.levels is not NULL: sage_free(PS.levels)
        return NULL
    PS.depth = 0
    PS.degree = n
    if unit_partition:
        for i from 0 <= i < n-1:
            PS.entries[i] = i
            PS.levels[i] = n
        PS.entries[n-1] = n-1
        PS.levels[n-1] = -1
    return PS

cdef inline PartitionStack *PS_copy(PartitionStack *PS):
    """
    Allocate and return a pointer to a copy of PartitionStack PS. Returns a null
    pointer in the case of an allocation failure.
    """
    cdef int i
    cdef PartitionStack *PS2 = <PartitionStack *> \
                                sage_malloc(sizeof(PartitionStack))
    if PS2 is NULL:
        return PS2
    PS2.entries = <int *> sage_malloc( PS.degree * sizeof(int) )
    PS2.levels = <int *> sage_malloc( PS.degree * sizeof(int) )
    if PS2.entries is NULL or PS2.levels is NULL:
        if PS2.entries is not NULL: sage_free(PS2.entries)
        if PS2.levels is not NULL: sage_free(PS2.levels)
        return NULL
    PS_copy_from_to(PS, PS2)
    return PS2

cdef inline int PS_copy_from_to(PartitionStack *PS, PartitionStack *PS2):
    """
    Copy all data from PS to PS2.
    """
    PS2.depth = PS.depth
    PS2.degree = PS.degree
    for i from 0 <= i < PS.degree:
        PS2.entries[i] = PS.entries[i]
        PS2.levels[i] = PS.levels[i]

cdef inline PartitionStack *PS_from_list(object L):
    """
    Allocate and return a pointer to a PartitionStack representing L. Returns a
    null pointer in the case of an allocation failure.
    """
    cdef int cell, i, num_cells = len(L), cur_start = 0, cur_len, n = 0
    for cell from 0 <= cell < num_cells:
        n += len(L[cell])
    cdef PartitionStack *PS = PS_new(n, 0)
    cell = 0
    if PS is NULL:
        return PS
    while 1:
        cur_len = len(L[cell])
        for i from 0 <= i < cur_len:
            PS.entries[cur_start + i] = L[cell][i]
            PS.levels[cur_start + i] = n
        PS_move_min_to_front(PS, cur_start, cur_start+cur_len-1)
        cur_start += cur_len
        cell += 1
        if cell == num_cells:
            PS.levels[cur_start-1] = -1
            break
        PS.levels[cur_start-1] = 0
    return PS

cdef inline int PS_dealloc(PartitionStack *PS):
    sage_free(PS.entries)
    sage_free(PS.levels)
    sage_free(PS)
    return 0

cdef inline int PS_print(PartitionStack *PS):
    """
    Print a visual representation of PS.
    """
    cdef int i
    for i from 0 <= i < PS.depth + 1:
        PS_print_partition(PS, i)
    print

cdef inline int PS_print_partition(PartitionStack *PS, int k):
    """
    Print the partition at depth k.
    """
    s = '('
    for i from 0 <= i < PS.degree:
        s += str(PS.entries[i])
        if PS.levels[i] <= k:
            s += '|'
        else:
            s += ','
    s = s[:-1] + ')'
    print s

cdef inline bint PS_is_discrete(PartitionStack *PS):
    """
    Returns whether the deepest partition consists only of singleton cells.
    """
    cdef int i
    for i from 0 <= i < PS.degree:
        if PS.levels[i] > PS.depth:
            return 0
    return 1

cdef inline int PS_num_cells(PartitionStack *PS):
    """
    Returns the number of cells.
    """
    cdef int i, ncells = 0
    for i from 0 <= i < PS.degree:
        if PS.levels[i] <= PS.depth:
            ncells += 1
    return ncells

cdef inline int PS_move_min_to_front(PartitionStack *PS, int start, int end):
    """
    Makes sure that the first element of the segment of entries i with
    start <= i <= end is minimal.
    """
    cdef int i, min_loc = start, minimum = PS.entries[start]
    for i from start < i <= end:
        if PS.entries[i] < minimum:
            min_loc = i
            minimum = PS.entries[i]
    if min_loc != start:
        PS.entries[min_loc] = PS.entries[start]
        PS.entries[start] = minimum

cdef inline bint PS_is_mcr(PartitionStack *PS, int m):
    """
    Returns whether PS.elements[m] (not m!) is the smallest element of its cell.
    """
    return m == 0 or PS.levels[m-1] <= PS.depth

cdef inline bint PS_is_fixed(PartitionStack *PS, int m):
    """
    Returns whether PS.elements[m] (not m!) is in a singleton cell, assuming
    PS_is_mcr(PS, m) is already true.
    """
    return PS.levels[m] <= PS.depth

cdef inline int PS_clear(PartitionStack *PS):
    """
    Sets the current partition to the first shallower one, i.e. forgets about
    boundaries between cells that are new to the current level.
    """
    cdef int i, cur_start = 0
    for i from 0 <= i < PS.degree:
        if PS.levels[i] >= PS.depth:
            PS.levels[i] += 1
        if PS.levels[i] < PS.depth:
            PS_move_min_to_front(PS, cur_start, i)
            cur_start = i+1

cdef inline int PS_move_all_mins_to_front(PartitionStack *PS):
    """
    Move minimal cell elements to the front of each cell.
    """
    cdef int i, cur_start = 0
    for i from 0 <= i < PS.degree:
        if PS.levels[i] <= PS.depth:
            PS_move_min_to_front(PS, cur_start, i)
            cur_start = i+1

cdef inline int PS_split_point(PartitionStack *PS, int v):
    """
    Detaches the point v from the cell it is in, putting the singleton cell
    of just v in front. Returns the position where v is now located.
    """
    cdef int i = 0, index_of_v
    while PS.entries[i] != v:
        i += 1
    index_of_v = i
    while PS.levels[i] > PS.depth:
        i += 1
    if (index_of_v == 0 or PS.levels[index_of_v-1] <= PS.depth) \
       and PS.levels[index_of_v] > PS.depth:
        # if v is first (make sure v is not already alone)
        PS_move_min_to_front(PS, index_of_v+1, i)
        PS.levels[index_of_v] = PS.depth
        return index_of_v
    else:
        # If v is not at front, i.e. v is not minimal in its cell,
        # then move_min_to_front is not necessary since v will swap
        # with the first before making its own cell, leaving it at
        # the front of the other.
        i = index_of_v
        while i != 0 and PS.levels[i-1] > PS.depth:
            i -= 1
        PS.entries[index_of_v] = PS.entries[i+1] # move the second element to v
        PS.entries[i+1] = PS.entries[i] # move the first (min) to second
        PS.entries[i] = v # place v first
        PS.levels[i] = PS.depth
        return i

cdef inline int PS_first_smallest(PartitionStack *PS, bitset_t b):
    """
    Find the first occurrence of the smallest cell of size greater than one,
    store its entries to b, and return its minimum element.
    """
    cdef int i = 0, j = 0, location = 0, n = PS.degree
    bitset_zero(b)
    while 1:
        if PS.levels[i] <= PS.depth:
            if i != j and n > i - j + 1:
                n = i - j + 1
                location = j
            j = i + 1
        if PS.levels[i] == -1: break
        i += 1
    # location now points to the beginning of the first, smallest,
    # nontrivial cell
    i = location
    while 1:
        bitset_flip(b, PS.entries[i])
        if PS.levels[i] <= PS.depth: break
        i += 1
    return PS.entries[location]

cdef inline int PS_get_perm_from(PartitionStack *PS1, PartitionStack *PS2, int *gamma):
    """
    Store the permutation determined by PS2[i] -> PS1[i] for each i, where PS[i]
    denotes the entry of the ith cell of the discrete partition PS.
    """
    cdef int i
    for i from 0 <= i < PS1.degree:
        gamma[PS2.entries[i]] = PS1.entries[i]

def PS_represent(partition=[[6],[3,4,8,7],[1,9,5],[0,2]], splits=[6,1,8,2]):
    """
    Demonstration and testing.

    DOCTEST:
        sage: from sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label \
        ... import PS_represent
        sage: PS_represent()
        Allocating PartitionStack...
        Allocation passed:
        (0,1,2,3,4,5,6,7,8,9)
        Checking that entries are in order and correct level.
        Everything seems in order, deallocating.
        Deallocated.
        Creating PartitionStack from partition [[6], [3, 4, 8, 7], [1, 9, 5], [0, 2]].
        PartitionStack's data:
        entries -> [6, 3, 4, 8, 7, 1, 9, 5, 0, 2]
        levels -> [0, 10, 10, 10, 0, 10, 10, 0, 10, -1]
        depth = 0, degree = 10
        (6|3,4,8,7|1,9,5|0,2)
        Checking PS_is_discrete:
        False
        Checking PS_num_cells:
        4
        Checking PS_is_mcr, min cell reps are:
        [6, 3, 1, 0]
        Checking PS_is_fixed, fixed elements are:
        [6]
        Copying PartitionStack:
        (6|3,4,8,7|1,9,5|0,2)
        Checking for consistency.
        Everything is consistent.
        Clearing copy:
        (0,3,4,8,7,1,9,5,6,2)
        Splitting point 6 from original:
        0
        (6|3,4,8,7|1,9,5|0,2)
        Splitting point 1 from original:
        5
        (6|3,4,8,7|1|5,9|0,2)
        Splitting point 8 from original:
        1
        (6|8|3,4,7|1|5,9|0,2)
        Splitting point 2 from original:
        8
        (6|8|3,4,7|1|5,9|2|0)
        Getting permutation from PS2->PS:
        [6, 1, 0, 8, 3, 9, 2, 7, 4, 5]
        Finding first smallest:
        Minimal element is 5, bitset is:
        0000010001
        Deallocating PartitionStacks.
        Done.

    """
    cdef int i, n = sum([len(cell) for cell in partition]), *gamma
    cdef bitset_t b
    print "Allocating PartitionStack..."
    cdef PartitionStack *PS = PS_new(n, 1), *PS2
    if PS is NULL:
        print "Allocation failed!"
        return
    print "Allocation passed:"
    PS_print(PS)
    print "Checking that entries are in order and correct level."
    good = True
    for i from 0 <= i < n-1:
        if not (PS.entries[i] == i and PS.levels[i] == n):
            print "Failed at i = %d!"%i
            print PS.entries[i], PS.levels[i], i, n
            good = False
    if not (PS.entries[n-1] == n-1 and PS.levels[n-1] == -1):
        print "Failed at i = %d!"%(n-1)
        good = False
    if not PS.degree == n or not PS.depth == 0:
        print "Incorrect degree or depth!"
        good = False
    if good: print "Everything seems in order, deallocating."
    PS_dealloc(PS)
    print "Deallocated."
    print "Creating PartitionStack from partition %s."%partition
    PS = PS_from_list(partition)
    print "PartitionStack's data:"
    print "entries -> %s"%[PS.entries[i] for i from 0 <= i < n]
    print "levels -> %s"%[PS.levels[i] for i from 0 <= i < n]
    print "depth = %d, degree = %d"%(PS.depth,PS.degree)
    PS_print(PS)
    print "Checking PS_is_discrete:"
    print "True" if PS_is_discrete(PS) else "False"
    print "Checking PS_num_cells:"
    print PS_num_cells(PS)
    print "Checking PS_is_mcr, min cell reps are:"
    L = [PS.entries[i] for i from 0 <= i < n if PS_is_mcr(PS, i)]
    print L
    print "Checking PS_is_fixed, fixed elements are:"
    print [PS.entries[l] for l in L if PS_is_fixed(PS, l)]
    print "Copying PartitionStack:"
    PS2 = PS_copy(PS)
    PS_print(PS2)
    print "Checking for consistency."
    good = True
    for i from 0 <= i < n:
        if PS.entries[i] != PS2.entries[i] or PS.levels[i] != PS2.levels[i]:
            print "Failed at i = %d!"%i
            good = False
    if PS.degree != PS2.degree or PS.depth != PS2.depth:
        print "Failure with degree or depth!"
        good = False
    if good:
        print "Everything is consistent."
    print "Clearing copy:"
    PS_clear(PS2)
    PS_print(PS2)
    for s in splits:
        print "Splitting point %d from original:"%s
        print PS_split_point(PS, s)
        PS_print(PS)
    print "Getting permutation from PS2->PS:"
    gamma = <int *> sage_malloc(n * sizeof(int))
    PS_get_perm_from(PS, PS2, gamma)
    print [gamma[i] for i from 0 <= i < n]
    sage_free(gamma)
    print "Finding first smallest:"
    bitset_init(b, n)
    i = PS_first_smallest(PS, b)
    print "Minimal element is %d, bitset is:"%i
    print bitset_string(b)
    bitset_clear(b)
    print "Deallocating PartitionStacks."
    PS_dealloc(PS)
    PS_dealloc(PS2)
    print "Done."

# Functions

cdef inline int split_point_and_refine(PartitionStack *PS, int v, object S,
    int (*refine_and_return_invariant)\
         (PartitionStack *PS, object S, int *cells_to_refine_by, int ctrb_len),
    int *cells_to_refine_by):
    """
    Make the partition stack one longer by copying the last partition in the
    stack, split off a given point, and refine. Return the invariant given by
    the refinement function.

    INPUT:
    PS -- the partition stack to refine
    v -- the point to split
    S -- the structure
    refine_and_return_invariant -- the refinement function provided
    cells_to_refine_by -- an array, contents ignored

    """
    PS.depth += 1
    PS_clear(PS)
    cells_to_refine_by[0] = PS_split_point(PS, v)
    return refine_and_return_invariant(PS, S, cells_to_refine_by, 1)

