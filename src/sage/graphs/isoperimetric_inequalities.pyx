# cython: binding = True
r"""
Isoperimetric inequalities

This module contains various functions to compute isoperimetric numbers
of a graph.

Authors:

- Peleg Michaeli
- Vincent Delecroix
"""
#*****************************************************************************
#       Copyright (C) 2018 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function

from cysignals.signals cimport sig_on, sig_off
from cysignals.memory cimport check_malloc, sig_free

from sage.graphs.base.static_sparse_graph cimport short_digraph, init_short_digraph, free_short_digraph, out_degree
from sage.graphs.graph_decompositions.fast_digraph cimport FastDigraph, compute_out_neighborhood_cardinality, popcount32

from sage.rings.infinity import Infinity
from sage.rings.rational_field import QQ

def cheeger_constant(g):
    r"""
    Return the cheeger constant of the graph.

    The Cheeger constant of a graph `G = (V,E)` is the minimum of `|\partial S|
    / |Vol(S)|` where `Vol(S)` is the sum of degrees of element in `S`,
    `\partial S` is the edge boundary of `S` (number of edges with one end in
    `S` and one end in `V \setminus S`) and the minimum is taken over all
    non-empty subsets `S` of vertices so that `|Vol(S)| \leq |E|`.

    .. SEEALSO::

        Alternative but similar quantities can be obtained via the methods
        :meth:`edge_isoperimetric_number` and :meth:`vertex_isoperimetric_number`.

    EXAMPLES::

        sage: graphs.PetersenGraph().cheeger_constant()
        1/3

    The Cheeger constant of a cycle on `n` vertices is
    `1/\lfloor n/2 \rfloor`::

        sage: [graphs.CycleGraph(k).cheeger_constant() for k in range(2,10)]
        [1, 1, 1/2, 1/2, 1/3, 1/3, 1/4, 1/4]

    The Cheeger constant of a complete graph on `n` vertices is
    `\lceil n/2 \rceil / (n-1)`::

        sage: [graphs.CompleteGraph(k).cheeger_constant() for k in range(2,10)]
        [1, 1, 2/3, 3/4, 3/5, 2/3, 4/7, 5/8]

    For complete bipartite::

        sage: [graphs.CompleteBipartiteGraph(i,j).cheeger_constant() for i in range(2,7) for j in range(2, i)]
        [3/5, 1/2, 3/5, 5/9, 4/7, 5/9, 1/2, 5/9, 1/2, 5/9]

    More examples::

        sage: G = Graph([(0, 1), (0, 3), (0, 8), (1, 4), (1, 6), (2, 4), (2, 7), (2, 9),
        ....:            (3, 6), (3, 8), (4, 9), (5, 6), (5, 7), (5, 8), (7, 9)])
        sage: G.cheeger_constant()
        1/6

        sage: G = Graph([(0, 1), (0, 2), (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (3, 4), (3, 5)])
        sage: G.cheeger_constant()
        1/2

        sage: Graph([[1,2,3,4],[(1,2),(3,4)]]).cheeger_constant()
        0

    TESTS::

        sage: graphs.EmptyGraph().cheeger_constant()
        Traceback (most recent call last):
        ...
        ValueError: Cheeger constant is not defined for the empty graph
    """
    if g.is_directed():
        raise ValueError("Cheeger constant is only defined on non-oriented graph")
    g._scream_if_not_simple()
    if g.num_verts() == 0:
        raise ValueError("Cheeger constant is not defined for the empty graph")
    elif g.num_verts() == 1:
        return Infinity
    elif not g.is_connected():
        return QQ((0,1))

    cdef short_digraph sd         # a copy of the graph g
    cdef int * subgraph           # vertices of the subgraph (stack)
    cdef int * bitsubgraph        # vertices of the subgraph (bit array of +1 (in) or -1 (not in))
    cdef int k = 0                # number of vertices in subgraph
    cdef int vol = 0              # number of edges in the subgraph
    cdef int boundary = 0         # number of edges in the boundary
    cdef int u = 0                # current vertex
    cdef unsigned long bmin = 1   # value of boundary for the min
    cdef unsigned long vmin = 1   # value of the volume for the min
    cdef int i

    init_short_digraph(sd, g)

    subgraph = <int *> check_malloc(sd.n * sizeof(int))
    bitsubgraph = <int *> check_malloc(sd.n * sizeof(int))

    sig_on()
    try:
        for i in range(sd.n):
            bitsubgraph[i] = -1
        while True:
            while True:
                # add vertex u to the subgraph, update the boundary/volume
                # we repeat the operation until we reach the maximum volume
                # or have no more vertex available
                for i in range(sd.neighbors[u+1] - sd.neighbors[u]):
                    boundary -= bitsubgraph[sd.neighbors[u][i]]
                    vol += 1

                bitsubgraph[u] = 1
                subgraph[k] = u
                u += 1
                k += 1

                if vol > sd.m:
                    break

                if boundary * vmin < bmin * vol:
                    bmin = boundary
                    vmin = vol

                if u == sd.n:
                    break

            # backtrack
            k -= 1
            u = subgraph[k]

            bitsubgraph[u] = -1
            for i in range(sd.neighbors[u+1] - sd.neighbors[u]):
                boundary += bitsubgraph[sd.neighbors[u][i]]
                vol -= 1
            u += 1

            if u == sd.n:
                if k == 0:
                    # end of the loop
                    break
                else:
                    # remove one more vertex in order to continue
                    k -= 1
                    u = subgraph[k]

                    bitsubgraph[u] = -1
                    for i in range(sd.neighbors[u+1] - sd.neighbors[u]):
                        boundary += bitsubgraph[sd.neighbors[u][i]]
                        vol -= 1
                    u += 1
        return QQ((bmin, vmin))

    finally:
        free_short_digraph(sd)
        sig_free(subgraph)
        sig_free(bitsubgraph)
        sig_off()

def edge_isoperimetric_number(g):
    r"""
    Return the edge-isoperimetric number of the graph.

    The edge-isoperimetric number of a graph `G = (V,E)` is also sometimes
    called the *isoperimetric number*. It is defined as the minimum of
    `|\partial S| / |S|` where `\partial S` is the edge boundary of `S` (number
    of edges with one end in `S` and one end in `V \setminus S`) and the
    minimum is taken over all subsets of vertices whose cardinality does not
    exceed half the size `|V|` of the graph.

    .. SEEALSO::

        Alternative but similar quantities can be obtained via the methods
        :meth:`cheeger_constant` and :meth:`vertex_isoperimetric_number`.

    EXAMPLES:

    The edge-isoperimetric number of a complete graph on `n` vertices is
    `\lceil n/2 \rceil`::

        sage: [graphs.CompleteGraph(n).edge_isoperimetric_number() for n in range(2,10)]
        [1, 2, 2, 3, 3, 4, 4, 5]

    The edge-isoperimetric constant of a cycle on `n` vertices is
    `2/\lfloor n/2 \rfloor`::

        sage: [graphs.CycleGraph(n).edge_isoperimetric_number() for n in range(2,15)]
        [1, 2, 1, 1, 2/3, 2/3, 1/2, 1/2, 2/5, 2/5, 1/3, 1/3, 2/7]

    In general, for `d`-regular graphs the edge-isoperimetric number is
    `d` times larger than the Cheeger constant of the graph::

        sage: g = graphs.RandomRegular(3, 10)
        sage: g.edge_isoperimetric_number() == g.cheeger_constant() * 3
        True

    And the edge-isoperimetric constant of a disconnected graph is `0`::

        sage: Graph([[1,2,3,4],[(1,2),(3,4)]]).edge_isoperimetric_number()
        0

    TESTS::

        sage: graphs.EmptyGraph().edge_isoperimetric_number()
        Traceback (most recent call last):
        ...
        ValueError: edge-isoperimetric number not defined for the empty graph
    """
    if g.is_directed():
        raise ValueError("edge isoperimetric number is only defined on non-oriented graph")
    g._scream_if_not_simple()
    if g.num_verts() == 0:
        raise ValueError("edge-isoperimetric number not defined for the empty graph")
    elif g.num_verts() == 1:
        return Infinity
    elif not g.is_connected():
        return QQ((0,1))

    cdef short_digraph sd           # a copy of the graph g
    cdef int * subgraph           # vertices of the subgraph (stack)
    cdef int * bitsubgraph        # vertices of the subgraph (bit array of +1 (in) or -1 (not in))
    cdef int k = 0                  # number of vertices in subgraph
    cdef unsigned long vol = 0      # number of edges in the subgraph
    cdef unsigned long boundary = 0 # number of edges in the boundary
    cdef int u = 0                  # current vertex
    cdef int i

    init_short_digraph(sd, g)

    cdef unsigned long bmin = sd.neighbors[1] - sd.neighbors[0] # value of boundary for the min
    cdef unsigned long vmin = 1     # value of the volume for the min

    subgraph = <int *> check_malloc(sd.n * sizeof(int))
    bitsubgraph = <int *> check_malloc(sd.n * sizeof(int))

    sig_on()
    try:
        for i in range(sd.n):
            bitsubgraph[i] = -1

        while True:
            while True:
                # add vertex u to the subgraph, update the boundary/volume
                # we repeat the operation until we reach the maximum volume
                # or have no more vertex available
                for i in range(sd.neighbors[u+1] - sd.neighbors[u]):
                    boundary -= bitsubgraph[sd.neighbors[u][i]]
                vol += 1

                bitsubgraph[u] = 1
                subgraph[k] = u
                u += 1
                k += 1

                if 2*vol > sd.n:
                    break

                if boundary * vmin < bmin * vol:
                    bmin = boundary
                    vmin = vol

                if u == sd.n:
                    break

            # backtrack
            k -= 1
            u = subgraph[k]

            bitsubgraph[u] = -1
            for i in range(sd.neighbors[u+1] - sd.neighbors[u]):
                boundary += bitsubgraph[sd.neighbors[u][i]]
            vol -= 1
            u += 1

            if u == sd.n:
                if k == 0:
                    # end of the loop
                    break
                else:
                    # remove one more vertex in order to continue
                    k -= 1
                    u = subgraph[k]

                    bitsubgraph[u] = -1
                    for i in range(sd.neighbors[u+1] - sd.neighbors[u]):
                        boundary += bitsubgraph[sd.neighbors[u][i]]
                    vol -= 1
                    u += 1

        return QQ((bmin, vmin))

    finally:
        sig_off()
        free_short_digraph(sd)
        sig_free(subgraph)
        sig_free(bitsubgraph)

def vertex_isoperimetric_number(g):
    r"""
    Return the vertex-isoperimetric number of the graph.

    The vertex-isoperimetric number of a graph `G = (V,E)` is also sometimes
    called the *magnifying constant*. It is defined as the minimum of `|N(S)| /
    |S|` where `|N(S)|` is the vertex boundary of `S` and the minimum is taken
    over the subsets `S` of vertices of size at most half of the vertices.

    .. SEEALSO::

        Alternative but similar quantities can be obtained via the methods
        :meth:`cheeger_constant` and :meth:`edge_isoperimetric_number`.

    EXAMPLES:

    The vertex-isoperimetric number of a complete graph on `n` vertices is
    `\lceil n/2 \rceil/\lfloor n/2 \rfloor`::

        sage: [graphs.CompleteGraph(k).vertex_isoperimetric_number() for k in range(2,15)]
        [1, 2, 1, 3/2, 1, 4/3, 1, 5/4, 1, 6/5, 1, 7/6, 1]

    The vertex-isoperimetric constant of a cycle on `n` vertices is
    `2/\lfloor n/2 \rfloor`::

        sage: [graphs.CycleGraph(k).vertex_isoperimetric_number() for k in range(2,15)]
        [1, 2, 1, 1, 2/3, 2/3, 1/2, 1/2, 2/5, 2/5, 1/3, 1/3, 2/7]

    And the vertex-isoperimetric constant of a disconnected graph is `0`::

        sage: Graph([[1,2,3],[(1,2)]]).vertex_isoperimetric_number()
        0

    TESTS::

        sage: graphs.EmptyGraph().vertex_isoperimetric_number()
        Traceback (most recent call last):
        ...
        ValueError: vertex-isoperimetric number not defined for the empty graph

        sage: graphs.CompleteGraph(64).vertex_isoperimetric_number()
        Traceback (most recent call last):
        ...
        ValueError: vertex-isoperimetric number can not be computed on graphs with more than 32 vertices
    """
    if g.is_directed():
        raise ValueError("vertex-isoperimetric number is only defined on non-oriented graph")
    g._scream_if_not_simple()
    if g.num_verts() == 0:
        raise ValueError("vertex-isoperimetric number not defined for the empty graph")
    elif g.num_verts() == 1:
        return Infinity
    elif not g.is_connected():
        return QQ((0,1))
    elif g.num_verts() > 32:
        raise ValueError("vertex-isoperimetric number can not be computed on graphs with more than 32 vertices")

    cdef FastDigraph FG = FastDigraph(g)

    cdef int card, boundary, S
    cdef int mc = FG.n // 2
    cdef int p = FG.n
    cdef int q = 0

    sig_on()

    for S in range(1, 2**FG.n):
        card = popcount32(S)
        if card <= mc:
            boundary = compute_out_neighborhood_cardinality(FG, S)
            if boundary * q < p * card:
                p = boundary
                q = card

    sig_off()
    return QQ((p,q))
