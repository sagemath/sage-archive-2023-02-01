r"""
This module implements an algorithm for generating all strong orientations of
an undirected graph.
It is an adaptation of the algorithm published in [CGMRV16]_.

A strong orientation of a graph is an orientation of its edges such that
the obtained digraph is strongly connected (i.e. there exist a directed path
between each pair of vertices).

ALGORITHM:

It runs in `O(m*n)` amortized time, where `m` is the number of edges and
`n` is the number of vertices. The amortized time can be improved to O(m)
with a more involved method.
In order to avoid trivial symetries, the orientation of an arbitrary edge
is fixed before the start of the enumeration process.

NOTE:

Works only for simple graphs (no multiple edges).

AUTHORS:

- Kolja Knauer and Petru Valicov (2017-01-10): initial version

REFERENCE:

.. [CGMRV16] A. Conte, R. Grossi, A. Marino, R. Rizzi, L. Versari,
  "Directing Road Networks by Listing Strong Orientations.",
  Combinatorial Algorithms: Proceedings of 27th International Workshop,
  IWOCA 2016, August 17-19, 2016, pages 83--95
"""

#*****************************************************************************
#       Copyright (C) 2016 YOUR NAME <petru.valicov@lif.univ-mrs.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.graphs.spanning_tree import kruskal
from sage.graphs.digraph import DiGraph

def strong_orientations_iterator(self):
    r"""
    Returns an iterator over all strong orientations of self.

    First preprocesses the graph and generates a spanning tree.
    Then every orientation of the non-tree edges can be extended to at least
    one new strong orientation by orienting properly the edges of the spanning
    tree (this property is proved in [CGMRV16]_). Therefore, this function
    generates all partial orientations of the non-tree edges and then launches
    the generation algorithm described in [CGMRV16]_.

    INPUT:

    - an undirected graph.

    OUTPUT:

    - an iterator which will produce all strong orientations of this graph.

    NOTE:

    In order to avoid symetries an orientation of an arbitrary edge is fixed.

    EXAMPLES::

    A cycle has one possible (non-symetric) strong orientation::

        sage: g = graphs.CycleGraph(4)
        sage: it = g.strong_orientations_iterator()
        sage: len(list(it))
        1

        A tree cannot be strongly oriented::

        sage: g = graphs.RandomTree(100)
        sage: len(list(g.strong_orientations_iterator()))
        0

        Neither can be a disconnected graph::

        sage: g = graphs.CompleteGraph(6)
        sage: g.add_vertex(7)
        sage: len(list(g.strong_orientations_iterator()))
        0

    TESTS:

    The total number of strong orientations of a graph can be counted using
    the Tutte polynomial evaluated at points (0,2)::

        sage: g = graphs.PetersenGraph()
        sage: nr1 = len(list(g.strong_orientations_iterator()))
        sage: nr2 = g.tutte_polynomial()(0,2)
        sage: nr1 == nr2/2 # The Tutte polynomial counts also the symetrical orientations
        True

    """
    # if the graph has a bridge or is disconnected,
    # then it cannot be strongly oriented
    if not self.is_biconnected() :
        return
    
    V = self.vertices()
    Dg = DiGraph([self.vertices(), self.edges()], pos=self.get_pos())

    # compute an arbitrary spanning tree of the undirected graph
    te = kruskal(self)
    treeEdges = [(u,v) for u,v,_ in te]
    A = [edge for edge in self.edges(labels=False) if edge not in treeEdges]
    
    # initialization of the first binary word 00...0
    # corresponding to the current orientation of the non-tree edges
    existingAedges = [0]*len(A)
    
    # Make the edges of the spanning tree doubly oriented
    for e in treeEdges :
        if Dg.has_edge(e) :
            Dg.add_edge(e[1], e[0])
        else :
            Dg.add_edge(e)
    
    # Generate all orientations for non-tree edges (using Gray code)
    # Each of these orientations can be extended to a strong orientation
    # of G by orienting properly the tree-edges
    previousWord = 0
    i = 0
    nr = 2**(len(A)-1)
    while i < nr :
        word = (i >> 1) ^ i
        bitChanged = word ^ previousWord
        
        bit = 0
        while bitChanged > 1 :
            bitChanged >>= 1;
            bit += 1;

        previousWord = word;
        if existingAedges[bit] == 0 :
            Dg.reverse_edge(A[bit])
            existingAedges[bit] = 1
        else :
            Dg.reverse_edge(A[bit][1], A[bit][0])
            existingAedges[bit] = 0
        # launch the algorithm for enumeration of the solutions
        for sol in _strong_orientations_of_a_mixed_graph(Dg, V, treeEdges):
            yield sol
        i = i + 1


def _strong_orientations_of_a_mixed_graph(Dg, V, E):
    r"""
    Helper function for the generation of all strong orientations.

    Generates all strong orientations of a given partially directed graph
    (also called mixed graph). The algorithm finds bound edges i.e undirected
    edges whose orientation is forced and tries all possible orientations for
    the other edges. See [CGMRV16]_ for more details.

    INPUT:

    - Dg -- the mixed graph. The undirected edges are doubly oriented.
    - V -- the set of vertices
    - E -- the set of undirected edges (these edges are oriented in both ways).
      No labels are allowed.

    OUTPUT:

    - an iterator which will produce all strong orientations of the input
      partially directed graph.

    EXAMPLES::

        sage: from sage.graphs.strong_orientations_generator import _strong_orientations_of_a_mixed_graph
        sage: g = graphs.CycleGraph(5)
        sage: Dg = DiGraph(g) # all edges of g will be doubly oriented
        sage: it = _strong_orientations_of_a_mixed_graph(Dg, g.vertices(), g.edges(labels=False))
        sage: len(list(it)) # there are two orientations of this multigraph
        2
    """
    length = len(E)
    i = 0
    boundEdges = []
    while i < length :
        (u,v) = E[i];
        Dg.delete_edge(u,v)
        if not (v in Dg.depth_first_search(u)) :
            del E[i]
            length -= 1
            Dg.add_edge((u,v))
            Dg.delete_edge((v,u))
            boundEdges.append((v,u))
        else :
            Dg.add_edge((u,v))
            Dg.delete_edge((v,u))
            if not (u in Dg.depth_first_search(v)) :
                del E[i]
                length -= 1
                boundEdges.append((u,v))
                Dg.delete_edge(u,v)
            else :
                i += 1
            Dg.add_edge((v,u))
    
    # if true the obtained orientation is strong
    if not E :
        yield Dg.copy()
        
    else :
        (u,v) = E.pop()
        Dg.delete_edge((v,u))
        for orientation in _strong_orientations_of_a_mixed_graph(Dg, V, E):
            yield orientation
        Dg.add_edge((v,u))
        Dg.delete_edge(u,v)
        for orientation in _strong_orientations_of_a_mixed_graph(Dg, V, E):
            yield orientation
        Dg.add_edge(u,v)
        E.append((u,v))
    Dg.add_edges(boundEdges)
    E.extend(boundEdges)
