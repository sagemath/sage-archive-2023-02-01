r"""
Bandwidth

Definition
----------

The bandwidth `bw(M)` of a matrix `M` is the smallest integer `k` such that all
non-zero entries of `M` are at distance `k` from the diagonal. The bandwidth
`bw(G)` of a graph `G` is the minimum bandwidth of the adjacency matrix of `G`,
over all possible relabellings of its vertices.

**Path spanner:** alternatively, the bandwidth measures how tightly a path
represents the distance of a graph `G`. Indeed, if the vertices of `G` can be
ordered as `v_1,...,v_n` in such a way that `k \times d_G(v_i,v_j) \geq |i-j|` then
`bw(G)\leq k`.

    **Proof:** for all `v_i \sim v_j` (i.e. `d_G(v_i,v_j)=1`), the constraint
    ensures that `k\geq |i-j|`, meaning that adjacent vertices are at distance
    at most `k` in the path ordering. That alone is sufficient to ensure that
    `bw(G)\leq k`.

    As a byproduct, we obtain that `k \times d_G(v_i,v_j) \geq |i-j|` in
    general: let `v_{s_0},...,v_{s_i}` be the vertices of a shortest
    `(v_i,v_j)`-path. We have:

    .. MATH::

        k \times d_G(v_i,v_j) &=    k\times d_G(v_i,v_{s_0}) + k\times d_G(v_{s_0},v_{s_1}) + ... + k\times d_G(v_{s_{i-1}},v_{s_i}) + k\times d_G(v_{s_i},v_j)\\
                              &\geq |v_i-v_{s_0}| + |v_{s_0}-v_{s_1}| + ... + |v_{s_{i-1}}-v_{s_i}| + |v_{s_i}-v_j|\\
                              &\geq |v_i-v_j|\\

Satisfiability of a partial assignment
--------------------------------------

Let us suppose that the first `i` vertices `v_1,...,v_i` of `G` have already
been assigned positions `p_1,...,p_i` in an ordering of `V(G)` of bandwidth
`\leq k`. Where can `v_{i+1}` appear ?

Because of the previous definition, `p_{i+1}` must be at distance at most
`k\times d_G(v_1,v_{i+1})` from `p_1`, and in general at distance at most
`k\times d_G(v_j,v_{i+1})` from `p_j`. Each range is an interval of
`\{1,...,n\}\backslash \{p_1,...,p_i\}`, and because the intersection of two
intervals is again an interval we deduce that in order to satisfy all these
constraints simultaneously `p_j` must belong to an interval defined from this
partial assignment.

Applying this rule to all non-assigned vertices, we deduce that each of them
must be assigned to a given interval of `\{1,...,n\}`. Note that this can also
be extended to the already assigned vertices, by saying that `v_j` with `j<i`
must be assigned within the interval `[p_j,p_j]`.

This problem is not always satisfiable, e.g. 5 vertices cannot all be assigned
to the elements of `[10,13]`. This is a matching problem which, because all
admissible sets are intervals, can be solved quickly.

Solving the matching problem
----------------------------

Let `n` points `v_1,...,v_n` be given, along with two functions `m,M:[n]\mapsto
[n]`. Is there an ordering `p_1,...,p_n` of them such that `m(v_i) \leq p_i \leq
M(v_i)` ? This is equivalent to Hall's bipartite matching theorem, and can in
this specific case be solved by the following algorithm:

- Consider all vertices `v` sorted increasingly according to `M(v)`

- For each of them, assign to `v` the smallest position in `[m(v),M(v)]` which
  has not been assigned yet. If there is none, the assignment problem is not
  satisfiable.

Note that the latest operation can be performed with very few bitset operations
(provided that `n<64`).

The algorithm
-------------

This section contains totally subjective choices, that may be changed in the
hope to get better performances.

- Try to find a satisfiable ordering by filling positions, one after the other
  (and not by trying to find each vertex' position)

- Fill the positions in this order: `0,n-1,1,n-2,3,n-3, ...`

.. NOTE::

    There is some symmetry to break as the reverse of a satisfiable ordering is
    also a satisfiable ordering.

Functions
---------
"""
#*****************************************************************************
#          Copyright (C) 2015 Nathann Cohen <nathann.cohen@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from libc.stdint cimport uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, free
from sage.graphs.distances_all_pairs cimport all_pairs_shortest_path_BFS

ctypedef uint32_t index_t

ctypedef struct range_t:
    index_t m
    index_t M

def bandwidth(G, k=None):
    r"""
    Compute the bandwidth of a graph.

    For a definition of the bandwidth of a graph, see the documentation of the
    :mod:`~sage.graphs.graph_decompositions.bandwidth` module.

    INPUT:

    - ``G`` (a graph)

    - ``k`` -- set to an integer value to test whether `bw(G)\leq k`, or to
      ``None`` (default) to compute `bw(G)`.

    OUTPUT:

    When `k` is an integer value, the function returns either ``False`` or a
    pair ``(ordering, adjacency_matrix)``.

    When `k` is equal to ``None``, the function returns a triple ``(bw,
    ordering, adjacency_matrix)``.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.bandwidth import bandwidth
        sage: bandwidth(graphs.PetersenGraph(),3)
        False
        sage: bandwidth(graphs.PetersenGraph())
        (
                                           [0 1 1 0 1 0 0 0 0 0]
                                           [1 0 0 0 0 1 1 0 0 0]
                                           [1 0 0 1 0 0 0 1 0 0]
                                           [0 0 1 0 0 0 1 0 1 0]
                                           [1 0 0 0 0 0 0 0 1 1]
                                           [0 1 0 0 0 0 0 1 1 0]
                                           [0 1 0 1 0 0 0 0 0 1]
                                           [0 0 1 0 0 1 0 0 0 1]
                                           [0 0 0 1 1 1 0 0 0 0]
        5, [0, 4, 5, 8, 1, 9, 3, 7, 6, 2], [0 0 0 0 1 0 1 1 0 0]
        )
        sage: bandwidth(graphs.ChvatalGraph())
        (
                                                   [0 0 1 1 0 1 1 0 0 0 0 0]
                                                   [0 0 0 1 1 1 0 1 0 0 0 0]
                                                   [1 0 0 0 1 0 0 1 1 0 0 0]
                                                   [1 1 0 0 0 0 0 0 1 1 0 0]
                                                   [0 1 1 0 0 0 1 0 0 1 0 0]
                                                   [1 1 0 0 0 0 0 0 0 0 1 1]
                                                   [1 0 0 0 1 0 0 1 0 0 0 1]
                                                   [0 1 1 0 0 0 1 0 0 0 1 0]
                                                   [0 0 1 1 0 0 0 0 0 0 1 1]
                                                   [0 0 0 1 1 0 0 0 0 0 1 1]
                                                   [0 0 0 0 0 1 0 1 1 1 0 0]
        6, [0, 5, 9, 4, 10, 1, 6, 11, 3, 8, 7, 2], [0 0 0 0 0 1 1 0 1 1 0 0]
        )
    """
    # All that this function does is allocate/free the memory for function
    # bandwidth_C

    cdef int n = G.order()

    cdef unsigned short ** d                = <unsigned short **> malloc( n   * sizeof(unsigned short *))
    cdef unsigned short *  distances        = <unsigned short *>  malloc( n*n * sizeof(unsigned short  ))
    cdef index_t *         current          = <index_t *>         malloc( n   * sizeof(index_t))
    cdef index_t *         ordering         = <index_t *>         malloc( n   * sizeof(index_t))
    cdef index_t *         left_to_order    = <index_t *>         malloc( n   * sizeof(index_t))
    cdef index_t *         index_array_tmp  = <index_t *>         malloc( n   * sizeof(index_t))
    cdef range_t *         range_arrays     = <range_t *>         malloc( n*n * sizeof(range_t))
    cdef range_t **        ith_range_array  = <range_t **>        malloc( n   * sizeof(range_t *))
    cdef range_t *         range_array_tmp  = <range_t *>         malloc( n   * sizeof(range_t))

    if (d               is NULL or
        distances       is NULL or
        current         is NULL or
        ordering        is NULL or
        left_to_order   is NULL or
        index_array_tmp is NULL or
        ith_range_array is NULL or
        range_arrays    is NULL or
        range_array_tmp is NULL):

        free(d)
        free(distances)
        free(current)
        free(ordering)
        free(left_to_order)
        free(index_array_tmp)
        free(range_arrays)
        free(ith_range_array)
        free(range_array_tmp)
        raise MemoryError

    cdef int i,j,kk
    all_pairs_shortest_path_BFS(G,NULL,distances,NULL) # compute the distance matrix
    cdef list int_to_vertex = G.vertices()

    # fill d so that d[i][j] works
    for i in range(n):
        d[i] = distances + i*n

    # ith_range_array
    for i in range(n):
        ith_range_array[i] = range_arrays + i*n

    # initialize left_to_order
    for i in range(n):
        left_to_order[i] = i

    if k is None:
        for kk in range((n-1)//G.diameter(),n):
            if bandwidth_C(n,kk,d,current,ordering,left_to_order,index_array_tmp,ith_range_array,range_array_tmp):
                ans = True
                break
    else:
        ans = bool(bandwidth_C(n,k,d,current,ordering,left_to_order,index_array_tmp,ith_range_array,range_array_tmp))

    if ans:
        from sage.matrix.constructor import Matrix
        order = [int_to_vertex[ordering[i]] for i in range(n)]
        M = Matrix([[int(G.has_edge(u,v)) for u in order] for v in order])
        assert all(abs(i-j)<= (kk if k is None else k)
                   for i,j in M.dict())
        ans = (kk, order, M) if k is None else (order,M)

    free(d)
    free(distances)
    free(current)
    free(ordering)
    free(left_to_order)
    free(index_array_tmp)
    free(range_arrays)
    free(ith_range_array)
    free(range_array_tmp)
    return ans

cdef bint bandwidth_C(int n, int k,
                     unsigned short ** d,
                     index_t *         current,         # choice of vertex for the current position
                     index_t *         ordering,        # the actual ordering of vertices
                     index_t *         left_to_order,   # begins with the assigned vertices, ends with the others
                     index_t *         index_array_tmp, # tmp space
                     range_t **        ith_range_array, # array of ranges, for every step of the algorithm
                     range_t *         range_array_tmp):# tmp space

    cdef int i,v
    cdef int pi # the position for which a vertex is being chosen
    cdef int vi # the vertex being tested at position pi
    cdef int radius
    current[0] = -1

    # At first any vertex can be anywhere
    for v in range(n):
        ith_range_array[0][v].m = 0
        ith_range_array[0][v].M = n-1

    i = 0
    while True:
        current[i] += 1

        if current[i] == n: # All choices for this position have been exhausted
            if i == 0:
                return 0
            i = i-1
            left_to_order[i], left_to_order[current[i]] = left_to_order[current[i]], left_to_order[i]
            continue

        pi = (n-1-i//2) if (i%2) else (i//2) # 0, n-1,1,n-2,2,n-3,3, ... that's an ugly 'if'
        vi = left_to_order[current[i]]

        # Wrong choice
        if (ith_range_array[i][vi].m > pi or
            ith_range_array[i][vi].M < pi):
            continue

        # swap
        left_to_order[i], left_to_order[current[i]] = left_to_order[current[i]], left_to_order[i]
        ordering[pi] = vi

        if i == n-1: # We did it !
            return 1

        # build the range array of depth i+1 knowing that vertex vi is at
        # position pi.
        #
        # k*d[v][vi] >= |p_v-p_{vi}|
        for v in range(n):
            radius = k*d[v][vi]
            ith_range_array[i+1][v].m = max(<int> ith_range_array[i][v].m,pi-radius)
            ith_range_array[i+1][v].M = min(<int> ith_range_array[i][v].M,pi+radius)

        # check feasibility at depth i+1
        if is_matching_feasible(n,ith_range_array[i+1],range_array_tmp, index_array_tmp):
            i += 1
            current[i] = i-1
        else:
            # swap back
            left_to_order[i], left_to_order[current[i]] = left_to_order[current[i]], left_to_order[i]

cdef bint is_matching_feasible(int n, range_t * range_array, range_t * range_array_tmp, index_t * index_array_tmp):
    r"""
    Test if the matching is feasible

    INPUT:

    - ``n`` (integer) -- number of points

    - ``range_array`` -- associates to every point an interval in which the
    point must be given a position.

    - ``range_array_tmp`` -- temporary spaces with the same characteristics as
      ``range_array``

    - ``index_array_tmp`` -- temporary space to associate an integer to every
      point.

    OUTPUT:

    The function must return a boolean, and does not change the content of
    ``range_array``.
    """
    # Heuristic: check if some vertex has an empty range, that's an easy 'no'.
    cdef int v,M,m,j
    for v in range(n):
        if range_array[v].M < range_array[v].m:
            #print range_array[v].m, range_array[v].M
            return 0
        index_array_tmp[v] = 0

    # Sort the guys according to increasing value of M in O(n).
    #
    # Step 1: count the occurrences of each M
    for v in range(n):
        index_array_tmp[range_array[v].M] += 1

    # Step 2: sorted table
    for v in range(1,n):
        index_array_tmp[v] += index_array_tmp[v-1]

    for v in range(n):
        M = range_array[v].M
        m = range_array[v].m
        index_array_tmp[M] -= 1
        range_array_tmp[index_array_tmp[M]].M = M
        range_array_tmp[index_array_tmp[M]].m = m

    # Satisfiability. We use index_array_tmp as a bitset, and mark every
    # assigned position.
    for v in range(n):
        index_array_tmp[v] = 0

    for v in range(n):
        for j in range(range_array_tmp[v].m, range_array_tmp[v].M+1):
            if not index_array_tmp[j]:
                index_array_tmp[j] = 1
                break
        else:
            return 0
    return 1
