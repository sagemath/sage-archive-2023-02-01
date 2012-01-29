r"""
Vertex separation

This module implements several algorithms to compute the vertex
separation of a digraph and the corresponding ordering of the
vertices. Given an ordering `v_1, ..., v_n` of the vertices of `V(G)`,
its *cost* is defined as:

.. MATH::

    c(v_1, ..., v_n) = \max_{1\leq i \leq n} c'(\{v_1, ..., v_i\})

Where

.. MATH::

    c'(S) = |N^+_G(S)\backslash S|

The *vertex separation* of a digraph `G` is equal to the minimum cost
of an ordering of its vertices.

**Vertex separation and pathwidth**

The vertex separation is defined on a digraph, but one can obtain from a graph
`G` a digraph `D` with the same vertex set, and in which each edge `uv` of `G`
is replaced by two edges `uv` and `vu` in `D`. The vertex separation of `D` is
equal to the pathwidth of `G`, and the corresponding ordering of the vertices of
`D` encodes an optimal path-decomposition of `G`.

This is a result of Kinnersley [Kin92]_ and Bodlaender [Bod98]_.

**References**

.. [Kin92] *The vertex separation number of a graph equals its path-width*,
  Nancy G. Kinnersley,
  Information Processing Letters,
  Volume 42, Issue 6, Pages 345-350,
  24 July 1992

.. [Bod98] *A partial k-arboretum of graphs with bounded treewidth*,
  Hans L. Bodlaender,
  Theoretical Computer Science,
  Volume 209, Issues 1-2, Pages 1-45,
  6 December 1998

**Authors**

    - Nathann Cohen

Exponential algorithm for vertex separation
-------------------------------------------

In order to find an optimal ordering of the vertices for the vertex separation,
this algorithm tries to save time by computing the function `c'(S)` **at most
once** once for each of the sets `S\subseteq V(G)`. These values are stored in
an array of size `2^n` where reading the value of `c'(S)` or updating it can be
done in constant (and small) time.

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
out-neighbors have been evaluated by the algrithm.

.. NOTE::

    Because of its current implementation, this algorithm only works on graphs
    on less than 32 vertices. This can be changed to 64 if necessary, but 32
    vertices already require 4GB of memory.

**Lower bound on the vertex separation**

One can obtain a lower bound on the vertex separation of a graph in exponential
time but *small* memory by computing once the cost of each set `S`. Indeed, the
cost of a sequence `v_1, ..., v_n` corresponding to sets `\{v_1\}, \{v_1,v_2\},
..., \{v_1,...,v_n\}` is

.. MATH::

    \max c'(\{v_1\}),c'(\{v_1,v_2\}),...,c'(\{v_1,...,v_n\})\geq\max c'_1,...,c'_n

where `c_i` is the minimum cost of a set `S` on `i`
vertices. Evaluating the `c_i` can take time (and in particular more
than the previous exact algorithm), but it does not need much memory
to run.

Methods
-------
"""

include '../../ext/stdsage.pxi'
include '../../ext/cdefs.pxi'
include '../../ext/interrupt.pxi'
include 'fast_digraph.pyx'
from libc.stdint cimport uint8_t, int8_t

#*****************************************************************************
#          Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

###############
# Lower Bound #
###############

def lower_bound(G):
    r"""
    Returns a lower bound on the vertex separation of `G`

    INPUT:

    - ``G`` -- a digraph

    OUTPUT:

    A lower bound on the vertex separation of `D` (see the module's
    documentation).

    .. NOTE::

        This method runs in exponential time but has no memory constraint.


    EXAMPLE:

    On a cycle::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: g = digraphs.Circuit(6)
        sage: lower_bound(g)
        1

    TEST:

    Given anything else than a DiGraph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: g = graphs.CycleGraph(5)
        sage: lower_bound(g)
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a DiGraph.
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise ValueError("The parameter must be a DiGraph.")

    if G.order() >= 32:
        raise ValueError("The graph can have at most 31 vertices.")

    cdef FastDigraph FD = FastDigraph(G)
    cdef int * g = FD.graph
    cdef int n = FD.n

    # minimums[i] is means to store the value of c'_{i+1}
    minimums = <uint8_t *> sage_malloc(sizeof(uint8_t)* n)
    cdef unsigned int i

    # They are initialized to n
    for 0<= i< n:
        minimums[i] = n

    cdef uint8_t tmp, tmp_count

    # We go through all sets
    for 1<= i< <unsigned int> (1<<n):
        tmp_count = <uint8_t> popcount32(i)
        tmp = <uint8_t> compute_out_neighborhood_cardinality(FD, i)

        # And update the costs
        minimums[tmp_count-1] = minimum(minimums[tmp_count-1], tmp)

    # We compute the maximum of all those values
    for 1<= i< n:
        minimums[0] = maximum(minimums[0], minimums[i])

    cdef int min = minimums[0]

    sage_free(minimums)

    return min

####################
# Exact algorithms #
####################

def path_decomposition(G, verbose = False):
    r"""
    Returns the pathwidth of the given graph and the ordering of the vertices
    resulting in a corresponding path decomposition.

    INPUT:

    - ``G`` -- a digraph

    - ``verbose`` (boolean) -- whether to display information on the
      computations.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    .. NOTE::

        Because of its current implementation, this algorithm only works on
        graphs on less than 32 vertices. This can be changed to 54 if necessary,
        but 32 vertices already require 4GB of memory.

    EXAMPLE:

    The vertex separation of a circuit is equal to 2::

        sage: from sage.graphs.graph_decompositions.vertex_separation import path_decomposition
        sage: g = graphs.CycleGraph(6)
        sage: path_decomposition(g)
        (2, [0, 1, 2, 3, 4, 5])

    TEST:

    Given anything else than a Graph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import path_decomposition
        sage: g = digraphs.Circuit(6)
        sage: path_decomposition(g)
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a Graph.
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("The parameter must be a Graph.")

    from sage.graphs.digraph import DiGraph
    return vertex_separation(DiGraph(G), verbose = verbose)

def vertex_separation(G, verbose = False):
    r"""
    Returns an optimal ordering of the vertices and its cost for
    vertex-separation.

    INPUT:

    - ``G`` -- a digraph

    - ``verbose`` (boolean) -- whether to display information on the
      computations.

    OUTPUT:

    A pair ``(cost, ordering)`` representing the optimal ordering of the
    vertices and its cost.

    .. NOTE::

        Because of its current implementation, this algorithm only works on
        graphs on less than 32 vertices. This can be changed to 54 if necessary,
        but 32 vertices already require 4GB of memory.

    EXAMPLE:

    The vertex separation of a circuit is equal to 2::

        sage: from sage.graphs.graph_decompositions.vertex_separation import vertex_separation
        sage: g = digraphs.Circuit(6)
        sage: vertex_separation(g)
        (1, [0, 1, 2, 3, 4, 5])

    TEST:

    Given anything else than a DiGraph::

        sage: from sage.graphs.graph_decompositions.vertex_separation import lower_bound
        sage: g = graphs.CycleGraph(5)
        sage: lower_bound(g)
        Traceback (most recent call last):
        ...
        ValueError: The parameter must be a DiGraph.

    Graphs with non-integer vertices::

        sage: D=digraphs.DeBruijn(2,3)
        sage: vertex_separation(D)
        (2, ['000', '001', '100', '010', '101', '011', '110', '111'])
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise ValueError("The parameter must be a DiGraph.")

    if G.order() >= 32:
        raise ValueError("The graph should have at most 31 vertices !")

    cdef FastDigraph g = FastDigraph(G)

    if verbose:
        print "Memory allocation"
        g.print_adjacency_matrix()

    _sig_on

    cdef unsigned int mem = 1 << g.n
    cdef int8_t * neighborhoods = <int8_t *> sage_malloc(mem)

    if neighborhoods == NULL:
        raise MemoryError("Error allocating memory. I just tried to allocate "+str(mem>>10)+"MB, could that be too much ?")

    memset(neighborhoods, <int8_t> -1, mem)

    cdef int i,j , k
    for 0 <= k <g.n:
        if verbose:
            print "Looking for a strategy of cost", str(k)

        if exists(g, neighborhoods, 0, k) <= k:
            break

    if verbose:
        print "... Found !"
        print "Now computing the ordering"

    cdef list order = find_order(g, neighborhoods, k)

    # Relabelling the vertices
    cdef list vertices = G.vertices()
    for i, j in enumerate(order):
        order[i] = vertices[j]

    sage_free(neighborhoods)

    _sig_off

    return k, order

##############################################################################
# Actual algorithm, breadh-first search and updates of the costs of the sets #
##############################################################################

# Check whether an ordering with the given cost exists, and updates data in the
# neighborhoods array at the same time. See the module's documentation

cdef inline int exists(FastDigraph g, int8_t * neighborhoods, int current, int cost):

    # If this is true, it means the set has not been evaluated yet
    if neighborhoods[current] < 0:
        neighborhoods[current] = compute_out_neighborhood_cardinality(g, current)

    # If the cost of this set is too high, there is no point in going further.
    # Same thing if the current set is the whole vertex set.
    if neighborhoods[current] > cost or (current == (1<<g.n)-1):
        return neighborhoods[current]

    # Minimum of the costs of the outneighbors
    cdef int mini = g.n

    cdef int i
    cdef int next_set


    for 0<= i<g.n:
        if (current >> i)&1:
            continue

        # For each of the out-neighbors next_set of current

        next_set = current | 1<<i

        # Check whether there exists a cheap path toward {1..n}, and updated the
        # cost.
        mini = minimum(mini, exists(g, neighborhoods, next_set, cost))

        # We have found a path !
        if mini <= cost:
            return mini

    # Updating the cost of the current set with the minimum of the cost of its
    # outneighbors.
    neighborhoods[current] = mini

    return neighborhoods[current]

# Returns the ordering once we are sure it exists
cdef list find_order(FastDigraph g, int8_t * neighborhoods, int cost):
    cdef list ordering = []
    cdef int current = 0
    cdef int n = g.n
    cdef int i

    while n:
        # We look for n vertices

        for 0<= i<g.n:
            if (current >> i)&1:
                continue

            # Find the next set with small cost (we know it exists)
            next_set = current | 1<<i
            if neighborhoods[next_set] <= cost:
                ordering.append(i)
                current = next_set
                break

        # One less to find
        n -= 1

    return ordering

# Min/Max functions

cdef inline int minimum(int a, int b):
    if a<b:
        return a
    else:
        return b

cdef inline int maximum(int a, int b):
    if a>b:
        return a
    else:
        return b
