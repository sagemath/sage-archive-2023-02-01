r"""
Bipartite Graphs

This module implements bipartite graphs.

AUTHORS:
    -- Robert L. Miller (2008-01-20): initial version

TESTS:

    sage: B = graphs.CompleteBipartiteGraph(7,9)
    sage: loads(dumps(B)) == B
    True

"""

#*****************************************************************************
#         Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from graph import Graph

class BipartiteGraph(Graph):

    def __init__(self, *args, **kwds):
        r"""
        Bipartite graph.

        INPUT:
        1. Empty: creates a zero vertex bipartite graph.

            sage: B = BipartiteGraph()
            sage: type(B)
            <class 'sage.graphs.bipartite_graph.BipartiteGraph'>
            sage: B.order()
            0

        2. From a graph: without any more information, finds a bipartition.

            sage: B = BipartiteGraph( graphs.CycleGraph(4) )
            sage: B = BipartiteGraph( graphs.CycleGraph(5) )
            Traceback (most recent call last):
            ...
            TypeError: Input graph is not bipartite!

        3. From a NetworkX bipartite graph.

            sage: import networkx
            sage: G = graphs.OctahedralGraph()
            sage: N = networkx.cliques.make_clique_bipartite(G.networkx_graph())
            sage: B = BipartiteGraph(N)

        4. From a graph and a partition. Note that if the input graph is not
        bipartite, then Sage will raise an error. However, if one specifies
        check = False, the offending edges are simply deleted (along with those
        vertices not appearing in either list).

            sage: P = graphs.PetersenGraph()
            sage: partition = [range(5), range(5,10)]
            sage: B = BipartiteGraph(P, partition)
            Traceback (most recent call last):
            ...
            TypeError: Input graph is not bipartite with respect to the given partition!

            sage: B = BipartiteGraph(P, partition, check=False)
            sage: B.left
            [0, 1, 2, 3, 4]
            sage: B.show()

        EXAMPLES:
        Test for arbitrary argument handled by Graph class
            sage: B = BipartiteGraph(None)
            sage: B
            Bipartite graph on 0 vertices

        Copy constructor
            sage: G = Graph({0:[5,6], 1:[4,5], 2:[4,6], 3:[4,5,6]})
            sage: B = BipartiteGraph(G)
            sage: B2 = BipartiteGraph(B)
            sage: B == B2
            True
            sage: B3 = BipartiteGraph(G, range(4), range(4,7))
            sage: B3
            Bipartite graph on 7 vertices
            sage: B3 == B2
            True

        Make sure "copy constructor" returns the same partition for no edges
            sage: G = Graph({0:[], 1:[], 2:[]})
            sage: part = (range(2), [2])
            sage: B = BipartiteGraph(G, part)
            sage: B2 = BipartiteGraph(B)
            sage: B == B2
            True
        """
        if len(args) == 0:
            Graph.__init__(self)
            self.left = []; self.right = []
            return
        arg1 = args[0]
        args = args[1:]
        if isinstance(arg1, BipartiteGraph):
            Graph.__init__(self, arg1, *args, **kwds)
            self.left, self.right = arg1.left, arg1.right
        elif isinstance(arg1, Graph) and \
             len(args) > 0 and isinstance(args[0], (list,tuple)) and \
             len(args[0]) == 2 and isinstance(args[0][0], (list,tuple)):
                # Assume that args[0] is a bipartition
                from copy import copy
                left, right = args[0]; left = copy(left); right = copy(right)
                Graph.__init__(self, arg1.subgraph(list(set(left)|set(right))), *args, **kwds)
                if not kwds.has_key('check') or kwds['check']:
                    while len(left) > 0:
                        a = left.pop(0)
                        if len( set( arg1.neighbors(a) ) & set(left) ) != 0:
                            raise TypeError("Input graph is not bipartite with " + \
                             "respect to the given partition!")
                    while len(right) > 0:
                        a = right.pop(0)
                        if len( set( arg1.neighbors(a) ) & set(right) ) != 0:
                            raise TypeError("Input graph is not bipartite with " + \
                             "respect to the given partition!")
                else:
                    while len(left) > 0:
                        a = left.pop(0)
                        a_nbrs = set( arg1.neighbors(a) ) & set(left)
                        if len( a_nbrs ) != 0:
                            self.delete_edges([(a, b) for b in a_nbrs])
                    while len(right) > 0:
                        a = right.pop(0)
                        a_nbrs = set( arg1.neighbors(a) ) & set(right)
                        if len( a_nbrs ) != 0:
                            self.delete_edges([(a, b) for b in a_nbrs])
                self.left, self.right = copy(args[0][0]), copy(args[0][1])
        elif isinstance(arg1, Graph):
            try:
                Graph.__init__(self, arg1, *args, **kwds)
                self.left, self.right = self.bipartite_sets()
                return
            except:
                raise TypeError("Input graph is not bipartite!")
        else:
            import networkx
            Graph.__init__(self, arg1, *args, **kwds)
            if isinstance(arg1, (networkx.XGraph, networkx.Graph)):
                if hasattr(arg1, 'node_type'):
                    # Assume the graph is bipartite
                    self.left = []
                    self.right = []
                    for v in arg1.nodes_iter():
                        if arg1.node_type[v] == 'Bottom':
                            self.left.append(v)
                        elif arg1.node_type[v] == 'Top':
                            self.right.append(v)
                        else:
                            raise TypeError("NetworkX node_type defies bipartite assumtion (is not 'Top' or 'Bottom')")
                else:
                    try:
                        self.left, self.right = self.bipartite_sets()
                    except:
                        raise TypeError("Input graph is not bipartite!")

    def _repr_(self):
        r"""
        Returns a short string representation of self.

        EXAMPLE:
            sage: B = BipartiteGraph(graphs.CycleGraph(16))
            sage: B
            Bipartite cycle graph: graph on 16 vertices

        """
        s = Graph._repr_(self).lower()
        if 'bipartite' in s:
            return s.capitalize()
        else:
            return 'Bipartite ' + s

    def bipartition(self):
        r"""
        Returns the underlying bipartition of the bipartite graph.

        EXAMPLE:
            sage: B = BipartiteGraph( graphs.CycleGraph(4) )
            sage: B.bipartition()
            ([0, 2], [1, 3])

        """
        return (self.left, self.right)

    def project_left(self):
        r"""
        Projects self onto left vertices: edges are 2-paths in the original.

        EXAMPLE:

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: G = B.project_left()
            sage: G.order(), G.size()
            (10, 10)

        """
        G = Graph()
        G.add_vertices(self.left)
        for v in G:
            for u in self.neighbor_iterator(v):
                G.add_edges([(v,w) for w in self.neighbor_iterator(u)])
        return G

    def project_right(self):
        r"""
        Projects self onto right vertices: edges are 2-paths in the original.

        EXAMPLE:

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: G = B.project_right()
            sage: G.order(), G.size()
            (10, 10)

        """
        G = Graph()
        G.add_vertices(self.left)
        for v in G:
            for u in self.neighbor_iterator(v):
                G.add_edges([(v,w) for w in self.neighbor_iterator(u)])
        return G

    def plot(self, *args, **kwds):
        r"""
        Overrides Graph's plot function, to illustrate the bipartite nature.

        EXAMPLE:

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: B.plot()

        """
        if 'pos' not in kwds.keys():
            kwds['pos'] = None
        if kwds['pos'] is None:
            pos = {}
            l_len = len(self.left)
            r_len = len(self.right)
            if l_len == 1:
                pos[self.left[0]] = [-1, 0]
            elif l_len > 1:
                i = 0
                d = 2./(l_len-1)
                for v in self.left:
                    pos[v] = [-1, 1-i*d]
                    i += 1
            if r_len == 1:
                pos[self.right[0]] = [1, 0]
            elif r_len > 1:
                i = 0
                d = 2./(r_len-1)
                for v in self.right:
                    pos[v] = [1, 1-i*d]
                    i += 1
            kwds['pos'] = pos
        return Graph.plot(self, *args, **kwds)

