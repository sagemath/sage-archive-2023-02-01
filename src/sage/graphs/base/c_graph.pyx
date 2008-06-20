
#*******************************************************************************
#        Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

from graph_backends import GenericGraphBackend
from sage.rings.integer import Integer

cdef class CGraph:
    cdef int add_arc_unsafe(self, int u, int v) except? -1:
        raise NotImplementedError()
    cdef int has_arc_unsafe(self, int u, int v) except? -1:
        raise NotImplementedError()
    cdef int del_arc_unsafe(self, int u, int v) except? -1:
        raise NotImplementedError()
    cdef int out_neighbors_unsafe(self, int u, int *neighbors, int size) except? -2:
        raise NotImplementedError()
    cdef int in_neighbors_unsafe(self, int u, int *neighbors, int size) except? -2:
        raise NotImplementedError()

    def _in_degree(self, int v):
        """
        Return the number of edges coming into v.

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: SparseGraph(7)._in_degree(3)
            0

        """
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Vertex (%d) is not a vertex of the graph."%v)
        return self.in_degrees[v]

    def _out_degree(self, int v):
        """
        Return the number of edges coming out of v.

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: SparseGraph(7)._out_degree(3)
            0

        """
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Vertex (%d) is not a vertex of the graph."%v)
        return self.out_degrees[v]

    def _num_verts(self):
        """
        Return the number of vertices in the (di)graph.

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: SparseGraph(7)._num_verts()
            7

        """
        return self.num_verts

    def _num_arcs(self):
        """
        Return the number of arcs.

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: SparseGraph(7)._num_arcs()
            0

        """
        return self.num_arcs


class CGraphBackend(GenericGraphBackend):

    _cg = None

    def degree(self, v, directed):
        """
        Return the degree of v.

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraphBackend
            sage: B = SparseGraphBackend(7)
            sage: B.degree(3, False)
            0

        """
        if directed:
            return self._cg._in_degree(v) + self._cg._out_degree(v)
        else:
            return self._cg._out_degree(v)

    def has_vertex(self, v):
        """
        Returns whether v is a vertex of self.

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraphBackend
            sage: B = SparseGraphBackend(7)
            sage: B.has_vertex(6)
            True
            sage: B.has_vertex(7)
            False

        """
        if not isinstance(v, (int,Integer)):
            return False
        return (v >= 0 and v < self._cg._num_verts())

    def iterator_nbrs(self, v):
        """
        Returns an iterator over the neighbors of v.

        EXAMPLE:
            sage: P = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: list(P._backend.iterator_nbrs(0))
            [1, 4, 5]

        """
        return iter(set(self.iterator_in_nbrs(v)) | set(self.iterator_out_nbrs(v)))

    def iterator_in_nbrs(self, v):
        """
        Returns an iterator over the incoming neighbors of v.

        EXAMPLE:
            sage: P = DiGraph(graphs.PetersenGraph().to_directed(), implementation='c_graph')
            sage: list(P._backend.iterator_in_nbrs(0))
            [1, 4, 5]

        """
        return iter(self._cg.in_neighbors(v))

    def iterator_out_nbrs(self, v):
        """
        Returns an iterator over the outgoing neighbors of v.

        EXAMPLE:
            sage: P = DiGraph(graphs.PetersenGraph().to_directed(), implementation='c_graph')
            sage: list(P._backend.iterator_out_nbrs(0))
            [1, 4, 5]

        """
        return iter(self._cg.out_neighbors(v))

    def iterator_verts(self, verts):
        """
        Returns an iterator over the vertices of self intersected with verts.

        EXAMPLE:
            sage: P = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: list(P._backend.iterator_verts(P))
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        """
        if verts is None:
            return iter(xrange(self.num_verts()))
        elif isinstance(verts, (int, Integer)):
            if 0 <= verts and verts < self.num_verts():
                return iter([verts])
            else:
                return
        else:
            return iter([i for i in xrange(self.num_verts()) if i in verts])

    def loops(self, new):
        """
        Returns whether loops are allowed in this graph.

        EXAMPLE:
            sage: G = Graph(implementation='c_graph')
            sage: G._backend.loops(None)
            False
            sage: G._backend.loops(True)
            sage: G._backend.loops(None)
            True

        """
        if new is None:
            return self._loops
        if new:
            self._loops = True
        else:
            self._loops = False

    def name(self, new):
        """
        Returns the name of this graph.

        EXAMPLE:
            sage: G = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: G._backend.name(None)
            'Petersen graph'

        """
        if new is None:
            return self._name
        self._name = new

    def num_edges(self, directed):
        """
        Returns the number of edges in self.

        EXAMPLE:
            sage: G = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: G._backend.num_edges(False)
            15

        """
        if directed:
            return self._cg._num_arcs()
        else:
            i = self._cg._num_arcs()
            k = 0
            if self.loops(None):
                if self.multiple_edges():
                    for j in xrange(self.num_verts()):
                        k += len(self.get_edge_label(j, j))
                else:
                    for j in xrange(self.num_verts()):
                        if self.has_edge(j,j,None):
                            k += 1
            i = (i-k)/2
            return i + k

    def num_verts(self):
        """
        Returns the number of vertices of self.

        EXAMPLE:
            sage: G = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: G._backend.num_verts()
            10

        """
        return self._cg._num_verts()

    def relabel(self, perm, directed):
        """
        Relabels the graph according to perm.

        EXAMPLE:
            sage: G = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: G._backend.relabel(range(9,-1,-1), False)
            sage: G.edges()
            [(0, 2, None),
             (0, 3, None),
             (0, 5, None),
             (1, 3, None),
             (1, 4, None),
             (1, 6, None),
             (2, 4, None),
             (2, 7, None),
             (3, 8, None),
             (4, 9, None),
             (5, 6, None),
             (5, 9, None),
             (6, 7, None),
             (7, 8, None),
             (8, 9, None)]

        """
        old_edges = list(self.iterator_edges(xrange(self.num_verts()), True, directed))
        for u, v, l in old_edges:
            self.del_edge(u, v, l, directed)
        for u, v, l in old_edges:
            self.add_edge(perm[u], perm[v], l, directed)


