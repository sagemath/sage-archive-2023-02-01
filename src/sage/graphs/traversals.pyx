# -*- coding: utf-8 -*-
"""
GenericGraph traversals implemented in Cython for efficiency.
"""

import collections

import cython

from libc.string cimport memset

from sage.ext.memory_allocator cimport MemoryAllocator
from sage.graphs.base.static_sparse_graph cimport short_digraph
from sage.graphs.base.static_sparse_graph cimport init_short_digraph
from sage.graphs.base.static_sparse_graph cimport free_short_digraph
from sage.graphs.base.static_sparse_graph cimport out_degree, has_edge

def lex_BFS(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a Lex BFS on the graph.

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

    Time complexity is `O(n+m)` where n is the number of vertices and m is the
    number of edges.

    The implementation of the algorithm is described in
    `this <http://www.cs.toronto.edu/~krueger/papers/unified.ps>`_ article
    by Derek G. Corneil and Richard Krueger.

    EXAMPLES:

    A Lex BFS is obviously an ordering of the vertices::

        sage: from sage.graphs.traversals import lex_BFS
        sage: g = graphs.CompleteGraph(6)
        sage: len(lex_BFS(g)) == g.order()
        True

    Lex BFS ordering of the 3-sun graph

        sage: from sage.graphs.traversals import lex_BFS
        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: lex_BFS(g)
        [1, 2, 3, 5, 4, 6]

    The method also works for directed graphs

        sage: from sage.graphs.traversals import lex_BFS
        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: lex_BFS(G, initial_vertex=2)
        [2, 3, 1]

    """
    # Loops and multiple edges are not needed in Lex BFS
    G.remove_multiple_edges()
    G.remove_loops()

    cdef int nV = G.order()

    if nV == 0:
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

    def l(x):
        return code[x]

    cdef list value = []

    # initialize the predecessors array
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *pred = <int *>mem.allocarray(nV, sizeof(int))
    memset(pred, -1, nV * sizeof(int))

    cdef set vertices = set(range(nV))

    cdef int source = 0 if initial_vertex is None else int_to_v.index(initial_vertex)
    code[source].append(nV + 1)

    cdef int now = 1, v, int_neighbor
    while vertices:
        v = max(vertices, key=l)
        vertices.remove(v)
        for i in range (0, out_degree(sd, v)):
            int_neighbor = sd.neighbors[v][i]
            if int_neighbor in vertices:
                code[int_neighbor].append(nV - now)
                pred[int_neighbor] = v
        value.append(int_to_v[v])
        now += 1

    if reverse:
        value.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        g = DiGraph(sparse=True)
        g.add_vertices(G.vertices())
        edges = [(int_to_v[i], int_to_v[pred[i]]) for i in range(nV) if pred[i] != -1]
        g.add_edges(edges)
        return value, g

    else:
        return value

def lex_DFS(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a Lex DFS on the graph.

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

    Time complexity is `O(n+m)` where n is the number of vertices and m is the
    number of edges.

    The implementation of the algorithm is described in
    `this <http://www.cs.toronto.edu/~krueger/papers/unified.ps>`_ article
    by Derek G. Corneil and Richard Krueger.

    EXAMPLES:

    A Lex DFS is obviously an ordering of the vertices::

        sage: from sage.graphs.traversals import lex_DFS
        sage: g = graphs.CompleteGraph(6)
        sage: len(lex_DFS(g)) == g.order()
        True

    Lex DFS ordering of the 3-sun graph

        sage: from sage.graphs.traversals import lex_DFS
        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: lex_DFS(g)
        [1, 2, 3, 5, 6, 4]

    The method also works for directed graphs

        sage: from sage.graphs.traversals import lex_DFS
        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: lex_DFS(G, initial_vertex=2)
        [2, 3, 1]

    """
    # Loops and multiple edges are not needed in Lex DFS
    G.remove_multiple_edges()
    G.remove_loops()

    cdef int nV = G.order()

    if nV == 0:
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

    def l(x):
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
        v = max(vertices, key=l)
        vertices.remove(v)
        for i in range (0, out_degree(sd, v)):
            int_neighbor = sd.neighbors[v][i]
            if int_neighbor in vertices:
                code[int_neighbor].appendleft(now)
                pred[int_neighbor] = v
        value.append(int_to_v[v])
        now += 1

    if reverse:
        value.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        g = DiGraph(sparse=True)
        g.add_vertices(G.vertices())
        edges = [(int_to_v[i], int_to_v[pred[i]]) for i in range(nV) if pred[i] != -1]
        g.add_edges(edges)
        return value, g

    else:
        return value
