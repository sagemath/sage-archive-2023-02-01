# -*- coding: utf-8 -*-
# cython: binding=True
r"""
Graph traversals.

**This module implements the following graph traversals**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~lex_BFS` | Perform a lexicographic breadth first search (LexBFS) on the graph.
    :meth:`~lex_UP` | Perform a lexicographic UP search (LexUP) on the graph.
    :meth:`~lex_DFS` | Perform a lexicographic depth first search (LexDFS) on the graph.
    :meth:`~lex_DOWN` | Perform a lexicographic DOWN search (LexDOWN) on the graph.

Methods
-------
"""
# ****************************************************************************
# Copyright (C) 2019 Georgios Giapitzakis Tzintanos <giorgosgiapis@mail.com>
#                    David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import collections

from libc.string cimport memset
from sage.ext.memory_allocator cimport MemoryAllocator
from sage.graphs.base.static_sparse_graph cimport short_digraph
from sage.graphs.base.static_sparse_graph cimport init_short_digraph
from sage.graphs.base.static_sparse_graph cimport free_short_digraph
from sage.graphs.base.static_sparse_graph cimport out_degree, has_edge

def lex_BFS(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a lexicographic breadth first search (LexBFS) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it for the first time)

    - ``initial_vertex`` -- (default: ``None``); the first vertex to
      consider

    ALGORITHM:

    This algorithm maintains for each vertex left in the graph a code
    corresponding to the vertices already removed. The vertex of maximal
    code (according to the lexicographic order) is then removed, and the
    codes are updated.

    Time complexity is `O(n+m)` where `n` is the number of vertices and `m` is
    the number of edges.

    See [CK2008]_ for more details on the algorithm.

    .. SEEALSO::

        * :wikipedia:`Lexicographic_breadth-first_search`

    EXAMPLES:

    A Lex BFS is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: len(g.lex_BFS()) == g.order()
        True

    Lex BFS ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_BFS()
        [1, 2, 3, 5, 4, 6]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: G.lex_BFS(initial_vertex=2)
        [2, 3, 1]

    For a Chordal Graph, a reversed Lex BFS is a Perfect Elimination Order::

        sage: g = graphs.PathGraph(3).lexicographic_product(graphs.CompleteGraph(2))
        sage: g.lex_BFS(reverse=True)  # py2
        [(2, 0), (2, 1), (1, 1), (1, 0), (0, 0), (0, 1)]
        sage: g.lex_BFS(reverse=True)  # py3
        [(2, 1), (2, 0), (1, 1), (1, 0), (0, 1), (0, 0)]

    And the vertices at the end of the tree of discovery are, for chordal
    graphs, simplicial vertices (their neighborhood is a complete graph)::

        sage: g = graphs.ClawGraph().lexicographic_product(graphs.CompleteGraph(2))
        sage: v = g.lex_BFS()[-1]
        sage: peo, tree = g.lex_BFS(initial_vertex = v,  tree=True)
        sage: leaves = [v for v in tree if tree.in_degree(v) ==0]
        sage: all(g.subgraph(g.neighbors(v)).is_clique() for v in leaves)
        True

    Different orderings for different traversals::

        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']

    TESTS:

    Lex BFS ordering of a graph on one vertex::

        sage: Graph(1).lex_BFS(tree=True)
        ([0], Digraph on 1 vertex)

    Lex BFS ordering of an empty (di)graph is an empty sequence::

        sage: g = Graph()
        sage: g.lex_BFS()
        []

    Lex BFS ordering of a symmetric digraph should be the same as the Lex BFS
    ordering of the corresponding undirected graph::

        sage: G = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: H = DiGraph(G)
        sage: G.lex_BFS() == H.lex_BFS()
        True

    """
    # Loops and multiple edges are not needed in Lex BFS
    if G.allows_loops() or G.allows_multiple_edges():
        G = G.to_simple(immutable=False)

    cdef int nV = G.order()

    if not nV:
        if tree:
            from sage.graphs.digraph import DiGraph
            g = DiGraph(sparse=True)
            return [], g
        else:
            return []

    # Build adjacency list of G
    cdef list int_to_v = list(G)

    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_v)

    # Perform Lex BFS

    cdef list code = [[] for i in range(nV)]

    def l_func(x):
        return code[x]

    cdef list value = []

    # Initialize the predecessors array
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *pred = <int *>mem.allocarray(nV, sizeof(int))
    memset(pred, -1, nV * sizeof(int))

    cdef set vertices = set(range(nV))

    cdef int source = 0 if initial_vertex is None else int_to_v.index(initial_vertex)
    code[source].append(nV + 1)

    cdef int now = 1, v, int_neighbor
    while vertices:
        v = max(vertices, key=l_func)
        vertices.remove(v)
        for i in range(0, out_degree(sd, v)):
            int_neighbor = sd.neighbors[v][i]
            if int_neighbor in vertices:
                code[int_neighbor].append(nV - now)
                pred[int_neighbor] = v
        value.append(int_to_v[v])
        now += 1

    free_short_digraph(sd)

    if reverse:
        value.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        g = DiGraph(sparse=True)
        g.add_vertices(G)
        edges = [(int_to_v[i], int_to_v[pred[i]]) for i in range(nV) if pred[i] != -1]
        g.add_edges(edges)
        return value, g

    else:
        return value

def lex_UP(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a lexicographic UP search (LexUP) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it for the first time)

    - ``initial_vertex`` -- (default: ``None``); the first vertex to
      consider

    ALGORITHM:

    This algorithm maintains for each vertex left in the graph a code
    corresponding to the vertices already removed. The vertex of maximal
    code (according to the lexicographic order) is then removed, and the
    codes are updated. During the `i`-th iteration of the algorithm `i` is
    appended to the codes of all neighbors of the selected vertex that are left
    in the graph.

    Time complexity is `O(n+m)` where `n` is the number of vertices and `m` is
    the number of edges.

    See [Mil2017]_ for more details on the algorithm.

    .. SEEALSO::

        * :meth:`~sage.graphs.generic_graph.GenericGraph.Lex_BFS` -- perform a
          lexicographic breadth first search (LexBFS) on the graph

    EXAMPLES:

    A Lex UP is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: len(g.lex_UP()) == g.order()
        True

    Lex UP ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_UP()
        [1, 2, 4, 5, 6, 3]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: G.lex_UP(initial_vertex=2)
        [2, 3, 1]

    Different orderings for different traversals::

        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']

    TESTS:

    Lex UP ordering of a graph on one vertex::

        sage: Graph(1).lex_UP(tree=True)
        ([0], Digraph on 1 vertex)

    Lex UP ordering of an empty (di)graph is an empty sequence::

        sage: g = Graph()
        sage: g.lex_UP()
        []

    Lex UP ordering of a symmetric digraph should be the same as the Lex UP
    ordering of the corresponding undirected graph::

        sage: G = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: H = DiGraph(G)
        sage: G.lex_UP() == H.lex_UP()
        True

    """
    # Loops and multiple edges are not needed in Lex UP
    if G.allows_loops() or G.allows_multiple_edges():
        G = G.to_simple(immutable=False)

    cdef int nV = G.order()

    if not nV:
        if tree:
            from sage.graphs.digraph import DiGraph
            g = DiGraph(sparse=True)
            return [], g
        else:
            return []

    # Build adjacency list of G
    cdef list int_to_v = list(G)

    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_v)

    # Perform Lex UP

    cdef list code = [[] for i in range(nV)]

    def l_func(x):
        return code[x]

    cdef list value = []

    # Initialize the predecessors array
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *pred = <int *>mem.allocarray(nV, sizeof(int))
    memset(pred, -1, nV * sizeof(int))

    cdef set vertices = set(range(nV))

    cdef int source = 0 if initial_vertex is None else int_to_v.index(initial_vertex)
    code[source].append(nV + 1)

    cdef int now = 1, v, int_neighbor
    while vertices:
        v = max(vertices, key=l_func)
        vertices.remove(v)
        for i in range(0, out_degree(sd, v)):
            int_neighbor = sd.neighbors[v][i]
            if int_neighbor in vertices:
                code[int_neighbor].append(now)
                pred[int_neighbor] = v
        value.append(int_to_v[v])
        now += 1

    free_short_digraph(sd)

    if reverse:
        value.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        g = DiGraph(sparse=True)
        g.add_vertices(G)
        edges = [(int_to_v[i], int_to_v[pred[i]]) for i in range(nV) if pred[i] != -1]
        g.add_edges(edges)
        return value, g

    else:
        return value

def lex_DFS(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a lexicographic depth first search (LexDFS) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it for the first time)

    - ``initial_vertex`` -- (default: ``None``); the first vertex to
      consider

    ALGORITHM:

    This algorithm maintains for each vertex left in the graph a code
    corresponding to the vertices already removed. The vertex of maximal
    code (according to the lexicographic order) is then removed, and the
    codes are updated. Lex DFS differs from Lex BFS only in the way codes are
    updated after each iteration.

    Time complexity is `O(n+m)` where `n` is the number of vertices and `m` is
    the number of edges.

    See [CK2008]_ for more details on the algorithm.

    EXAMPLES:

    A Lex DFS is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: len(g.lex_DFS()) == g.order()
        True

    Lex DFS ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_DFS()
        [1, 2, 3, 5, 6, 4]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: G.lex_DFS(initial_vertex=2)
        [2, 3, 1]

    Different orderings for different traversals::

        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']

    TESTS:

    Lex DFS ordering of a graph on one vertex::

        sage: Graph(1).lex_DFS(tree=True)
        ([0], Digraph on 1 vertex)

    Lex DFS ordering of an empty (di)graph is an empty sequence::

        sage: g = Graph()
        sage: g.lex_DFS()
        []

    Lex DFS ordering of a symmetric digraph should be the same as the Lex DFS
    ordering of the corresponding undirected graph::

        sage: G = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: H = DiGraph(G)
        sage: G.lex_DFS() == H.lex_DFS()
        True

    """
    # Loops and multiple edges are not needed in Lex DFS
    if G.allows_loops() or G.allows_multiple_edges():
        G = G.to_simple(immutable=False)

    cdef int nV = G.order()

    if not nV:
        if tree:
            from sage.graphs.digraph import DiGraph
            g = DiGraph(sparse=True)
            return [], g
        else:
            return []

    # Build adjacency list of G
    cdef list int_to_v = list(G)

    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_v)

    # Perform Lex DFS

    # We are using deque in order to prepend items in list efficiently
    cdef list code = [collections.deque([]) for i in range(nV)]

    def l_func(x):
        return code[x]

    cdef list value = []

    # initialize the predecessors array
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *pred = <int *>mem.allocarray(nV, sizeof(int))
    memset(pred, -1, nV * sizeof(int))

    cdef set vertices = set(range(nV))

    cdef int source = 0 if initial_vertex is None else int_to_v.index(initial_vertex)
    code[source].appendleft(0)

    cdef int now = 1, v, int_neighbor
    while vertices:
        v = max(vertices, key=l_func)
        vertices.remove(v)
        for i in range(0, out_degree(sd, v)):
            int_neighbor = sd.neighbors[v][i]
            if int_neighbor in vertices:
                code[int_neighbor].appendleft(now)
                pred[int_neighbor] = v
        value.append(int_to_v[v])
        now += 1

    free_short_digraph(sd)

    if reverse:
        value.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        g = DiGraph(sparse=True)
        g.add_vertices(G)
        edges = [(int_to_v[i], int_to_v[pred[i]]) for i in range(nV) if pred[i] != -1]
        g.add_edges(edges)
        return value, g

    else:
        return value

def lex_DOWN(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a lexicographic DOWN search (LexDOWN) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it for the first time)

    - ``initial_vertex`` -- (default: ``None``); the first vertex to
      consider

    ALGORITHM:

    This algorithm maintains for each vertex left in the graph a code
    corresponding to the vertices already removed. The vertex of maximal
    code (according to the lexicographic order) is then removed, and the
    codes are updated. During the `i`-th iteration of the algorithm `n-i` is
    prepended to the codes of all neighbors of the selected vertex that are left
    in the graph.

    Time complexity is `O(n+m)` where `n` is the number of vertices and `m` is
    the number of edges.

    See [Mil2017]_ for more details on the algorithm.

    .. SEEALSO::

        * :meth:`~sage.graphs.generic_graph.GenericGraph.Lex_DFS` -- perform a
          lexicographic breadth depth first search (LexDFS) on the graph

    EXAMPLES:

    A Lex DOWN is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: len(g.lex_DOWN()) == g.order()
        True

    Lex DOWN ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_DOWN()
        [1, 2, 3, 4, 6, 5]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: G.lex_DOWN(initial_vertex=2)
        [2, 3, 1]

    Different orderings for different traversals::

        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']

    TESTS:

    Lex DOWN ordering of a graph on one vertex::

        sage: Graph(1).lex_DOWN(tree=True)
        ([0], Digraph on 1 vertex)

    Lex DOWN ordering of an empty (di)graph is an empty sequence::

        sage: g = Graph()
        sage: g.lex_DOWN()
        []

    Lex DOWN ordering of a symmetric digraph should be the same as the Lex DOWN
    ordering of the corresponding undirected graph::

        sage: G = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: H = DiGraph(G)
        sage: G.lex_DOWN() == H.lex_DOWN()
        True

    """
    # Loops and multiple edges are not needed in Lex DOWN
    if G.allows_loops() or G.allows_multiple_edges():
        G = G.to_simple(immutable=False)

    cdef int nV = G.order()

    if not nV:
        if tree:
            from sage.graphs.digraph import DiGraph
            g = DiGraph(sparse=True)
            return [], g
        else:
            return []

    # Build adjacency list of G
    cdef list int_to_v = list(G)

    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_v)

    # Perform Lex DOWN

    # We are using deque in order to prepend items in list efficiently
    cdef list code = [collections.deque([]) for i in range(nV)]

    def l_func(x):
        return code[x]

    cdef list value = []

    # initialize the predecessors array
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *pred = <int *>mem.allocarray(nV, sizeof(int))
    memset(pred, -1, nV * sizeof(int))

    cdef set vertices = set(range(nV))

    cdef int source = 0 if initial_vertex is None else int_to_v.index(initial_vertex)
    code[source].appendleft(0)

    cdef int now = 1, v, int_neighbor
    while vertices:
        v = max(vertices, key=l_func)
        vertices.remove(v)
        for i in range(0, out_degree(sd, v)):
            int_neighbor = sd.neighbors[v][i]
            if int_neighbor in vertices:
                code[int_neighbor].appendleft(nV - now)
                pred[int_neighbor] = v
        value.append(int_to_v[v])
        now += 1

    free_short_digraph(sd)

    if reverse:
        value.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        g = DiGraph(sparse=True)
        g.add_vertices(G)
        edges = [(int_to_v[i], int_to_v[pred[i]]) for i in range(nV) if pred[i] != -1]
        g.add_edges(edges)
        return value, g

    else:
        return value
