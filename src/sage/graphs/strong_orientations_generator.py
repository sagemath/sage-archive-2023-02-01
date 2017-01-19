r"""
Returns an iterator over all strong orientations of self.


A strong orientation of a graph is an orientation of its edges such that
that the obtained digraph is strongly connected (i.e. there exist a directed path
between each pair of vertices).

ALGORITHM:

The algorithm is an adaptation of the algorithm published in [CGMRV16]. It runs in `O(m*n)` amortized time
where `m` is the number of edges and `n` is the number of vertices. The amortized time can be improved to O(m) with a more involved algorithm.
In order to avoid trivial symetries, the orientation of an arbitrary edge is fixed before the start of the enumeration process.
 

INPUT:

- the undirected graph.

OUTPUT:

- an iterator which will produce all strong orientations of the input graph.

NOTE:

In order to avoid symetries an orientation of an arbitrary edge is fixed.


AUTHORS:

- Kolja Knauer and Petru Valicov (2017-01-10): initial version

EXAMPLES:

A cycle has one possible (non-symetric) strong orientation::

    sage: g = graphs.CycleGraph(4)
    sage: it = g.all_strong_orientations_iterator()
    sage: len(list(it))
    1

TESTS:
The total number of strong orientations of a graph can be counted using the Tutte polynomial evaluated at points (0,2)::
    sage: g = graphs.PetersenGraph()
    sage: len(list(g.all_strong_orientations_iterator())) == int(g.tutte_polynomial()(0,2)/2) # the Tutte polynomial counts also the symetrical orientations
    True

REFERENCE:

- [CGMRV16] A. Conte, R. Grossi, A. Marino, R. Rizzi, L. Versari, Directing Road Networks by Listing Strong Orientations. 
*Combinatorial Algorithms: Proceedings of 27th International Workshop, IWOCA 2016*, August 17-19, 2016, pages 83--95

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

# MAIN function
# preprocesses the graph and launches the generation algorithm
def all_strong_orientations_iterator(self):
    # if the graph is a forest then it cannot be strongly oriented
    if self.edge_connectivity() <= 1:
        return
    
    V = self.vertices()
    Dg = DiGraph([self.vertices(), self.edges()])

    # compute an arbitrary spanning tree of the undirected graph
    te = kruskal(self)
    treeEdges = [(edge[0],edge[1]) for edge in self.edges() if edge in te]
    A = [edge for edge in self.edges(labels=False) if edge not in treeEdges]
    
    # initialization of the first binary word 00...0
    # corresponding to the current orientation of the non-tree edges
    existingAedges = []
    for i in xrange(0, len(A)):
        existingAedges.append(0)
    
    # Make the edges of the spanning tree double oriented
    for e in treeEdges :
        if Dg.has_edge(e) :
            Dg.add_edge(e[1], e[0])
        else :
            Dg.add_edge(e)
    
    # Generate all orientations for non-tree edges (using Gray code)
    # Each of these orientations can be extended to a strong orientation
    # of G by orienting properly the tree-edges
    previousWord = 0
    for i in xrange(0, 2**(len(A)-1)):
        word = (i >> 1) ^ i
        bitChanged = word ^ previousWord
        
        bit = 0
        while bitChanged > 1 :
            bitChanged >>= 1; bit += 1;

        previousWord = word;
        if existingAedges[bit] == 0 :
            Dg.reverse_edge(A[bit])
            existingAedges[bit] = 1
        else :
            Dg.reverse_edge(A[bit][1], A[bit][0])
            existingAedges[bit] = 0
        # launch the algorithm for enumeration of the solutions
        for orientation in core_generation_algorithm(Dg, V, treeEdges):
            yield orientation

# INPUT : directed graph, the set of vertices, the set of undirected edges
# explores all strong orientations by finding bound edges and trying the other edges
# See [CGMRV16] for more details
def core_generation_algorithm(Dg, V, E):
    length = len(E)
    i = 0
    boundEdges = []
    while i < length :
        (u,v) = E[i];
        Dg.delete_edge((u,v))
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
                Dg.delete_edge((u,v))
            else :
                i += 1
            Dg.add_edge((v,u))
    
    # if true the obtained orientation is strong
    if not E :
        yield Dg
        
    else :
        (u,v) = E.pop()
        Dg.delete_edge((v,u))
        for orientation in core_generation_algorithm(Dg, V, E):
                    yield orientation
        Dg.add_edge((v,u))
        Dg.delete_edge((u,v))
        for orientation in core_generation_algorithm(Dg, V, E):
                    yield orientation
        Dg.add_edge((u,v))
        E.append((u,v))
    Dg.add_edges(boundEdges)
    E.extend(boundEdges)
