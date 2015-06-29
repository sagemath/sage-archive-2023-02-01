# distutils: libraries = gmp
r"""
Cutwidth

This module implements several algorithms to compute the cutwidth of a graph and
the corresponding ordering of the vertices. It also implements tests functions
for evaluation the width of a linear ordering (or layout).

Given an ordering
`v_1,\cdots, v_n` of the vertices of `V(G)`, its *cost* is defined as:

.. MATH::

    c(v_1, ..., v_n) = \max_{1\leq i \leq n-1} c'(\{v_1, ..., v_i\})

Where

.. MATH::

    c'(S) = |\{(u,w)\in E(G)\mid u\in S\text{ and }w\in V(G)\backslash S\}|

The *cutwidth* of a graph `G` is equal to the minimum cost of an ordering of its
vertices.


**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`cutwidth` | Return the cutwidth of the graph and the corresponding vertex ordering.
    :meth:`cutwidth_dyn` | Compute the cutwidth of `G` using an exponential time and space algorithm based on dynamic programming
    :meth:`width_of_cut_decomposition` | Return the width of the cut decomposition induced by the linear ordering `L` of the vertices of `G`


Exponential algorithm for cutwidth
----------------------------------

In order to find an optimal ordering of the vertices for the vertex separation,
this algorithm tries to save time by computing the function `c'(S)` **at most
once** once for each of the sets `S\subseteq V(G)`. These values are stored in
an array of size `2^n` where reading the value of `c'(S)` or updating it can be
done in constant time.

Assuming that we can compute the cost of a set `S` and remember it, finding an
optimal ordering is an easy task. Indeed, we can think of the sequence `v_1,
..., v_n` of vertices as a sequence of *sets* `\{v_1\}, \{v_1,v_2\}, ...,
\{v_1,...,v_n\}`, whose cost is precisely `\max c'(\{v_1\}), c'(\{v_1,v_2\}),
... , c'(\{v_1,...,v_n\})`. Hence, when considering the digraph on the `2^n`
sets `S\subseteq V(G)` where there is an arc from `S` to `S'` if `S'=S\cap
\{v\}` for some `v` (that is, if the sets `S` and `S'` can be consecutive in a
sequence), an ordering of the vertices of `G` corresponds to a *path* from
`\emptyset` to `\{v_1,...,v_n\}`. In this setting, checking whether there exists
a ordering of cost less than `k` can be achieved by checking whether there
exists a directed path `\emptyset` to `\{v_1,...,v_n\}` using only sets of cost
less than `k`. This is just a depth-first-search, for each `k`.

**Lazy evaluation of** `c'`

In the previous algorithm, most of the time is actually spent on the computation
of `c'(S)` for each set `S\subseteq V(G)` -- i.e. `2^n` computations of
neighborhoods. This can be seen as a huge waste of time when noticing that it is
useless to know that the value `c'(S)` for a set `S` is less than `k` if all the
paths leading to `S` have a cost greater than `k`. For this reason, the value of
`c'(S)` is computed lazily during the depth-first search. Explanation :

When the depth-first search discovers a set of size less than `k`, the costs of
its out-neighbors (the potential sets that could follow it in the optimal
ordering) are evaluated. When an out-neighbor is found that has a cost smaller
than `k`, the depth-first search continues with this set, which is explored with
the hope that it could lead to a path toward `\{v_1,...,v_n\}`. On the other
hand, if an out-neighbour has a cost larger than `k` it is useless to attempt to
build a cheap sequence going though this set, and the exploration stops
there. This way, a large number of sets will never be evaluated and *a lot* of
computational time is saved this way.

Besides, some improvement is also made by "improving" the values found by
`c'`. Indeed, `c'(S)` is a lower bound on the cost of a sequence containing the
set `S`, but if all out-neighbors of `S` have a cost of `c'(S) + 5` then one
knows that having `S` in a sequence means a total cost of at least `c'(S) +
5`. For this reason, for each set `S` we store the value of `c'(S)`, and replace
it by `\max (c'(S), \min_{\text{next}})` (where `\min_{\text{next}}` is the
minimum of the costs of the out-neighbors of `S`) once the costs of these
out-neighbors have been evaluated by the algorithm.

This algorithm and its implementation are very similar to
:meth:`sage.graphs.graph_decompositions.vertex_separation.vertex_separation_exp`.
The main difference is in the computation of `c'(S)`. See the :mod:`vertex
separation module's documentation
<sage.graphs.graph_decompositions.vertex_separation>` for more details on this
algorithm.

.. NOTE::

    Because of its current implementation, this algorithm only works on graphs
    on strictly less than 32 vertices. This can be changed to 64 if necessary,
    but 32 vertices already require 4GB of memory.


Authors
-------

- David Coudert (2015-06): Initial version


Methods
-------
"""
#*****************************************************************************
#          Copyright (C) 2015 David Coudert <david.coudert@inria.fr>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#      as published by the Free Software Foundation; either version 2 of
#              the License, or (at your option) any later version.
#                        http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/stdsage.pxi'
include 'sage/ext/interrupt.pxi'
include 'sage/ext/cdefs.pxi'
from sage.graphs.graph_decompositions.fast_digraph cimport FastDigraph, popcount32
from sage.graphs.graph_decompositions.vertex_separation import is_valid_ordering
from libc.stdint cimport uint8_t
from sage.ext.memory cimport check_allocarray
from sage.rings.integer_ring import ZZ


################################################################################
# Measurement function of the width of some layout
################################################################################

def width_of_cut_decomposition(G, L):
    r"""
    Returns the width of the cut decomposition induced by the linear ordering
    `L` of the vertices of `G`.

    If `G` is an instance of :mod:`Graph <sage.graphs.graph>`, this function
    returns the width `cw_L(G)` of the cut decomposition induced by the linear
    ordering `L` of the vertices of `G`.

    .. MATH::

        cw_L(G) =  \max_{0\leq i< |V|-1} |\{(u,w)\in E(G)\mid u\in L[:i]\text{ and }w\in V(G)\setminus L[:i]\}|

    INPUT:

    - ``G`` -- a Graph

    - ``L`` -- a linear ordering of the vertices of ``G``

    EXAMPLES:

    Cut decomposition of a Cycle graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: G = graphs.CycleGraph(6)
        sage: L = G.vertices()
        sage: cutwidth.width_of_cut_decomposition(G, L)
        2

    Cut decomposition of a Path graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: P = graphs.PathGraph(6)
        sage: cutwidth.width_of_cut_decomposition(P, [0, 1, 2, 3, 4, 5])
        1
        sage: cutwidth.width_of_cut_decomposition(P, [5, 0, 1, 2, 3, 4])
        2
        sage: cutwidth.width_of_cut_decomposition(P, [0, 2, 4, 1, 3, 5])
        5

    TESTS:

    Giving a wrong linear ordering::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.width_of_cut_decomposition(Graph(), ['a','b'])
        Traceback (most recent call last):
        ...
        ValueError: The input linear vertex ordering L is not valid for G.
    """
    if not is_valid_ordering(G, L):
        raise ValueError("The input linear vertex ordering L is not valid for G.")

    position = {u:i for i,u in enumerate(L)}

    # We count for each position `i` the number of edges going from vertices at
    # positions in `0..i` to vertices at positions in `i+1..n-1`, for each
    # `x\leq i<n-1`.
    cpt = [0]*len(L)
    for u,v in G.edge_iterator(labels=None):
        x,y = position[u],position[v]
        if x>y:
            x,y = y,x
        # Edge (u,v) contributes 1 to the number of edges going from vertices at
        # positions `0..i` to vertices at positions `i+1..n-1` for each `x\leq
        # i<n-1`.
        for i in range(x,y):
            cpt[i] += 1

    # The width of L is the maximum computed value.
    return max(cpt)


################################################################################
# Front end method for cutwidth
################################################################################

def cutwidth(G, algorithm="exponential", cut_off=0):
    r"""
    Return the cutwidth of the graph and the corresponding vertex ordering.

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``algorithm`` -- (default: ``"exponential"``) Specify the algorithm to use
      among

      - ``exponential`` -- Use an exponential time and space algorithm based on
        dynamic programming. This algorithm only works on graphs with strictly
        less than 32 vertices.

    - ``cut_off`` -- (default: 0) This parameter is used to stop the search as
      soon as a solution with width at most ``cut_off`` is found, if any. If
      this bound cannot be reached, the best solution found is returned.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    EXAMPLES:

    Cutwidth of a Complete Graph::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: G = graphs.CompleteGraph(5)
        sage: cw,L = cutwidth(G, algorithm="exponential"); cw
        6
        sage: K = graphs.CompleteGraph(6)
        sage: cw,L = cutwidth(K, algorithm="exponential"); cw
        9
        sage: cw,L = cutwidth(K+K, algorithm="exponential"); cw
        9

    The cutwidth of a `p\times q` Grid Graph with `p\leq q` is `p+1`::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: G = graphs.Grid2dGraph(3,3)
        sage: cw,L = cutwidth(G, algorithm="exponential"); cw
        4
        sage: G = graphs.Grid2dGraph(3,5)
        sage: cw,L = cutwidth(G, algorithm="exponential"); cw
        4

    TESTS:

    Given a wrong algorithm::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: cutwidth(Graph(), algorithm="SuperFast")
        Traceback (most recent call last):
        ...
        ValueError: Algorithm "SuperFast" has not been implemented yet. Please contribute.

    Given anything else than a Graph::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: cutwidth(range(4))
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a Graph.

    Giving a wrong type cut off::

        sage: from sage.graphs.graph_decompositions.cutwidth import cutwidth
        sage: cutwidth(Graph(), cut_off='toto')
        Traceback (most recent call last):
        ...
        ValueError: The specified cut off parameter must be an integer.
    """
    from sage.graphs.graph import Graph

    if not isinstance(G, Graph):
        raise ValueError('The parameter must be a Graph.')

    if not cut_off in ZZ:
        raise ValueError("The specified cut off parameter must be an integer.")

    if not G.is_connected():
        # The graph has several connected components. We solve the problem on
        # each of them and concatenate the partial orderings. The cutwidth is
        # the maximum over all these subgraphs.
        cw, L = 0, []
        for V in G.connected_components():

            if len(V)==1:
                # We can directly add this vertex to the solution
                L.extend(V)

            else:
                # We build the connected subgraph and do a recursive call to get
                # its cutwidth and corresponding ordering
                H = G.subgraph(V)
                cwH,LH = cutwidth(H, algorithm = algorithm,
                                  cut_off      = max(cut_off,cw))

                # We update the cutwidth and ordering
                cw = max(cw, cwH)
                L.extend(LH)

        return cw, L

    # We have a connected graph and we call the desired algorithm
    if algorithm == "exponential":
        return cutwidth_dyn(G, lower_bound=cut_off)

    else:
        raise ValueError('Algorithm "{}" has not been implemented yet. Please contribute.'.format(algorithm))


################################################################################
# Dynamic Programming algorithm for cutwidth
################################################################################
from sage.graphs.graph_decompositions.vertex_separation cimport find_order

def cutwidth_dyn(G, lower_bound=0):
    r"""
    Dynamic programming algorithm for the cutwidth of a Graph.

    This function uses dynamic programming algorithm for determining an optimal
    layout for the cutwidth of `G`. See the :mod:`module's documentation
    <sage.graphs.graph_decompositions.cutwidth>` for more details on this
    method.

    INPUT:

    - ``G`` -- a Graph

    - ``lower_bound`` -- (default: 0) the algorithm returns immediately if it
      finds a solution lower or equal to ``lower_bound`` (in which case it may
      not be optimal).

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    .. NOTE::

        Because of its current implementation, this algorithm only works on
        graphs on strictly less than 32 vertices. This can be changed to 63 if
        necessary, but 32 vertices already require 4GB of memory.

    TESTS:

    Giving anything else than a Graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.cutwidth_dyn([])
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a Graph.

    Giving a too large Graph::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.cutwidth_dyn(graphs.PathGraph(40))
        Traceback (most recent call last):
        ...
        ValueError: The graph should have at most 31 vertices !

    Giving a wrong type lower bound::

        sage: from sage.graphs.graph_decompositions import cutwidth
        sage: cutwidth.cutwidth_dyn(Graph(), lower_bound='toto')
        Traceback (most recent call last):
        ...
        ValueError: The specified lower bound must be an integer.

    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("The parameter must be a Graph.")

    if G.order() >= 32:
        raise ValueError("The graph should have at most 31 vertices !")

    if not lower_bound in ZZ:
        raise ValueError("The specified lower bound must be an integer.")

    cdef FastDigraph g = FastDigraph(G)

    cdef unsigned int mem = 1 << g.n
    cdef uint8_t * neighborhoods = <uint8_t *> check_allocarray(mem, sizeof(uint8_t))

    memset(neighborhoods, <uint8_t> -1, mem)

    cdef int i, k

    try:
        sig_on()
        for k in range(lower_bound, G.size()):
            for i in range(g.n):
                if exists(g, neighborhoods, 0, 0, i, k) <= k:
                    sig_off()
                    order = find_order(g, neighborhoods, k)
                    return k, [g.int_to_vertices[i] for i in order]
        sig_off()
    finally:
        sage_free(neighborhoods)

    order = find_order(g, neighborhoods, k)
    return k, [g.int_to_vertices[i] for i in order]

cdef inline int exists(FastDigraph g, uint8_t * neighborhoods, int S, int cost_S, int v, int k):
    """
    Check whether an ordering with the given cost `k` exists, and updates data
    in the neighborhoods array at the same time. See the module's documentation.

    INPUT:

    - ``g`` -- a FastDiGraph

    - ``neighborhoods`` -- an array of size `2^(g.n)` storing for each subset
      `X\subseteq V` of vertices of the graph the number of edges from `X` to
      `V\setminus X`.

    - ``S`` -- an integer encoding the predecessor subset of vertices (from
      which is issued the current call).

    - ``cost_S`` -- the number of edges from `S` to `V\setminus S`.

    - ``v`` -- a vertex such that the current subset of vertices is
      `current==S\cup\{v\}`.

    - ``k`` -- the maximum admissible cost for a solution.
    """
    cdef int current = S | 1<<v
    # If this is true, it means the set has not been evaluated yet
    if neighborhoods[current] == <uint8_t>-1:
        # The number of edges from `current` to `V\setminus current` is the
        # number of edges from `S` to `V\setminus S`, minus the number of edges
        # from `S` to vertex `v`, plus the number of edges from `v` to
        # `V\setminus (S\cup \{v\})`. This can be computed adding the degree of
        # `c` to `cost_S`, and then removing twice the number of edges from `S`
        # to `v`.
        neighborhoods[current] = cost_S + g.degree[v] - 2*popcount32(S&g.graph[v])

    # If the cost of this set is too high, there is no point in going further.
    # Same thing if the current set is the whole vertex set.
    if neighborhoods[current] > k or (current == (1<<g.n)-1):
        return neighborhoods[current]

    # Minimum of the costs of the outneighbors, initialized with large constant.
    cdef int mini = (<uint8_t> -1)

    cdef int i
    cdef int next_set

    # For each possible extension of the current set witha vertex, check whether
    # there exists a cheap path toward {1..n}, and update the cost.
    for i in range(g.n):
        if (current >> i)&1: # if i in S
            continue

        mini = min(mini, exists(g, neighborhoods, current, neighborhoods[current], i, k))

        # We have found a path !
        if mini <= k:
            return mini

    # Updating the cost of the current set with the minimum of the cost of its
    # outneighbors.
    neighborhoods[current] = mini

    return neighborhoods[current]
