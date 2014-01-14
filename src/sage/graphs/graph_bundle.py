r"""
Graph Bundles

This module implements graph bundles.

WARNING::

    This module has known bugs.  See the ``plot`` method, for example.

AUTHORS:
    -- Robert L. Miller (2008-01-20): initial version

TESTS::

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

class GraphBundle(Graph):

    def __init__(self, *args, **kwds):
        r"""
        Graph Bundle.

        Note that an instance of the GraphBundle class is also a Graph instance-
        this is the total space of the bundle. The base graph is self.base.

        INPUT:
        1. Empty: creates a zero vertex trivial graph bundle.

            sage: B = GraphBundle()
            sage: type(B)
            <class 'sage.graphs.graph_bundle.GraphBundle'>
            sage: type(B.base)
            <class 'sage.graphs.graph.Graph'>
            sage: B.order()
            0

        2. From a graph and a partition: the partition determines the fibers,
        and we take the quotient map to the base.

            sage: P = graphs.PetersenGraph()
            sage: partition = [range(5), range(5,10)]
            sage: B = GraphBundle(P, partition)
            sage: B.base
            Graph on 2 vertices
            sage: B.base.size()
            1
            sage: B.fiber[0]
            [0, 1, 2, 3, 4]
            sage: B.projection[0]
            0
            sage: B.projection[5]
            1

        """
        if len(args) == 0:
            Graph.__init__(self, sparse=True)
            self.base = Graph()
            self.fiber = {}
            self.projection = {}
            return
        if isinstance(args[0], Graph):
            G = args[0]
            args = args[1:]
            if isinstance(args[0], (list, tuple)):
                if len(args[0])>0 and isinstance(args[0][0], (list, tuple)):
                    # Assume that the second argument is a partition of the vertices
                    from copy import copy
                    self.fiber = {}
                    self.projection = {}
                    partition = args[0]
                    args = args[1:]
                    Graph.__init__(self, G, sparse=True)
                    self.base = Graph(sparse=True)
                    base_size = len(partition)
                    self.base.add_vertices(xrange(base_size))
                    edge_list = []
                    for j in xrange(base_size):
                        par_j = partition[j]
                        self.fiber[j] = copy(par_j)
                        for i in xrange(j):
                            par_i = partition[i]
                            if len(par_i) < len(par_j):
                                if len( set(self.vertex_boundary(par_i)) & set(par_j) ) > 0:
                                    edge_list.append((i, j))
                            else:
                                if len( set(self.vertex_boundary(par_j)) & set(par_i) ) > 0:
                                    edge_list.append((i, j))
                        for v in par_j:
                            self.projection[v] = j
                    self.base.add_edges(edge_list)

    def edge_lift(self, left_vertex, right_vertex):
        r"""
        Returns the edge lift of an edge in the base. This is by definition a
        bipartite graph, so the order the vertices are input determines which
        partition is on the 'left.'

        EXAMPLE:

            sage: P = graphs.PetersenGraph()
            sage: partition = [range(5), range(5,10)]
            sage: B = GraphBundle(P, partition)
            sage: B.edge_lift(0,1)
            Bipartite petersen graph: graph on 10 vertices

        """
        from bipartite_graph import BipartiteGraph
        return BipartiteGraph(self, [self.fiber[left_vertex], self.fiber[right_vertex]], check=False)

    def _repr_(self):
        r"""
        Returns a short string representation of self.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: partition = [range(5), range(5,10)]
            sage: B = GraphBundle(P, partition)
            sage: B
            Graph bundle:
                Total space: petersen graph: graph on 10 vertices
                Base space: graph on 2 vertices

        """
        s_total = Graph._repr_(self).lower()
        s_base = Graph._repr_(self.base).lower()
        s = "Graph bundle:\n"
        s += "    Total space: " + s_total + "\n"
        s += "    Base space: " + s_base
        return s

    def fibers(self):
        r"""
        Returns a list of the fibers of the bundle, as a partition of the total
        space.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: partition = [range(5), range(5,10)]
            sage: B = GraphBundle(P, partition)
            sage: B.fibers()
            [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]

        """
        return self.fiber.values()

    def plot(self, *args, **kwds):
        r"""
        Overrides Graph's plot function, to illustrate the bundle nature.

        EXAMPLE::

            sage: P = graphs.PetersenGraph()
            sage: partition = [range(5), range(5,10)]
            sage: B = GraphBundle(P, partition)
            sage: #B.plot()

          Test disabled due to bug in GraphBundle.__init__().  See trac #8329.

        """
        if 'pos' not in kwds.keys():
            kwds['pos'] = None
        if kwds['pos'] is None:
            import sage.graphs.generic_graph_pyx as generic_graph_pyx
            if 'iterations' not in kwds.keys():
                kwds['iterations'] = 50
            iters = kwds['iterations']
            total_pos = generic_graph_pyx.spring_layout_fast(self, iterations=iters)
            base_pos = generic_graph_pyx.spring_layout_fast(self.base, iterations=iters)
            for v in base_pos.iterkeys():
                for v_tilde in self.fiber[v]:
                    total_pos[v_tilde][0] = base_pos[v][0]
            tot_x = [p[0] for p in total_pos.values()]
            tot_y = [p[1] for p in total_pos.values()]
            bas_x = [p[0] for p in base_pos.values()]
            bas_y = [p[1] for p in base_pos.values()]
            tot_x_min = min(tot_x)
            tot_x_max = max(tot_x)
            tot_y_min = min(tot_y)
            tot_y_max = max(tot_y)
            bas_x_min = min(bas_x)
            bas_x_max = max(bas_x)
            bas_y_min = min(bas_y)
            bas_y_max = max(bas_y)
            if tot_x_max == tot_x_min and tot_y_max == tot_y_min:
                tot_y_max += 1
                tot_y_min -= 1
            elif tot_y_max == tot_y_min:
                delta = (tot_x_max - tot_x_min)/2.0
                tot_y_max += delta
                tot_y_min -= delta
            if bas_x_max == bas_x_min and bas_y_max == bas_y_min:
                bas_y_max += 1
                bas_y_min -= 1
            elif bas_y_max == bas_y_min:
                delta = (bas_x_max - bas_x_min)/2.0
                bas_y_max += delta
                bas_y_min -= delta
            y_diff = (bas_y_max - tot_y_min) + 2*(bas_y_max - bas_y_min)
            pos = {}
            for v in self:
                pos[('t',v)] = [total_pos[v][0], total_pos[v][1] + y_diff]
            for v in self.base:
                pos[('b',v)] = base_pos[v]
            from copy import copy
            G = copy(self)
            rd = {}
            for v in G:
                rd[v] = ('t',v)
            G.relabel(rd)
            B = copy(self.base)
            rd = {}
            for v in B:
                rd[v] = ('b',v)
            B.relabel(rd)
            E = G.disjoint_union(B)
            kwds['pos'] = pos
            from sage.plot.all import arrow
            G = Graph.plot(E, *args, **kwds)
            G += arrow( ((tot_x_max + tot_x_min)/2.0, tot_y_min + y_diff),
                        ((tot_x_max + tot_x_min)/2.0, bas_y_max), axes=False )
            G.axes(False)
            return G
        else:
            return G.plot(self, *args, **kwds)

