
#*******************************************************************************
#        Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

from c_graph import CGraphBackend

cdef class DenseGraph(CGraph):
    """
    Implements compiled dense graphs, as an array of packed bits.

    Creating a new dense graph instance:
        G = DenseGraph(int nverts)

    INPUT:
        nverts -- non-negative integer, number of vertices.

    NOTES:
        DenseGraph does not distinguish whether it is directed or not. In fact,
    the datastructure itself is directed. An edge is simply an arc in both
    directions.

    """

    def __new__(self, int nverts):
        cdef int radix = sizeof(unsigned long) << 3
        self.radix_mod_mask = radix - 1
        cdef int i = 0
        while ((<unsigned long>1)<<i) & self.radix_mod_mask:
            i += 1
        self.radix_div_shift = i
        self.num_verts = nverts
        self.num_arcs = 0
#        print ''.join(reversed(['1' if (1<<i)&radix else '0' for i in xrange(sizeof(unsigned long)*8)]))
#        print ''.join(reversed(['1' if (1<<i)&(radix - 1) else '0' for i in xrange(sizeof(unsigned long)*8)]))
        i = nverts >> self.radix_div_shift
        if nverts & self.radix_mod_mask:
            i += 1
        self.num_longs = i
        self.edges = <unsigned long *> sage_malloc(nverts * self.num_longs * sizeof(unsigned long))
        self.in_degrees = <int *> sage_malloc(nverts * sizeof(int))
        self.out_degrees = <int *> sage_malloc(nverts * sizeof(int))
        if not self.edges or not self.in_degrees or not self.out_degrees:
            if self.edges: sage_free(self.edges)
            if self.in_degrees: sage_free(self.in_degrees)
            if self.out_degrees: sage_free(self.out_degrees)
            raise RuntimeError("Failure allocating memory.")
        for i from 0 <= i < self.num_longs * nverts:
            self.edges[i] = 0
        for i from 0 <= i < nverts:
            self.in_degrees[i] = 0
            self.out_degrees[i] = 0

    def __dealloc__(self):
        sage_free(self.edges)
        sage_free(self.in_degrees)
        sage_free(self.out_degrees)

    cdef int add_arc_unsafe(self, int u, int v):
        """
        Adds arc (u, v) to the graph.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        """
        cdef int place = (u * self.num_longs) + (v >> self.radix_div_shift)
        cdef unsigned long word = (<unsigned long>1) << (v & self.radix_mod_mask)
        if not self.edges[place] & word:
            self.in_degrees[v] += 1
            self.out_degrees[u] += 1
            self.num_arcs += 1
            self.edges[place] |= word

    def add_arc(self, int u, int v):
        """
        Adds arc (u, v) to the graph.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        EXAMPLE:
            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(4,7)
            Traceback (most recent call last):
            ...
            RuntimeError: Second vertex (7) is not a vertex of the graph.
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True

        """
        if u < 0 or u >= self.num_verts:
            raise RuntimeError("First vertex (%d) is not a vertex of the graph."%u)
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Second vertex (%d) is not a vertex of the graph."%v)
        self.add_arc_unsafe(u,v)

    cdef int has_arc_unsafe(self, int u, int v):
        """
        Checks whether arc (u, v) is in the graph.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        OUTPUT:
            0 -- False
            1 -- True

        """
        cdef int place = (u * self.num_longs) + (v >> self.radix_div_shift)
        cdef unsigned long word = (<unsigned long>1) << (v & self.radix_mod_mask)
        if self.edges[place] & word:
            return 1
        else:
            return 0

    def has_arc(self, int u, int v):
        """
        Checks whether arc (u, v) is in the graph.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        EXAMPLE:
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1)
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True

        """
        if u < 0 or u >= self.num_verts or v < 0 or v >= self.num_verts:
            return False
        return self.has_arc_unsafe(u,v) == 1

    cdef int del_arc_unsafe(self, int u, int v):
        """
        Deletes the arc from u to v.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        """
        cdef int place = (u * self.num_longs) + (v >> self.radix_div_shift)
        cdef unsigned long word = (<unsigned long>1) << (v & self.radix_mod_mask)
        if self.edges[place] & word:
            self.in_degrees[v] -= 1
            self.out_degrees[u] -= 1
            self.num_arcs -= 1
            self.edges[place] &= ~word

    def del_arc(self, int u, int v):
        """
        Deletes the arc from u to v.

        INPUT:
            u, v -- integers from 0, ..., n-1, where n is the number of vertices

        EXAMPLE:
            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.has_arc(0,1)
            True
            sage: G.del_arc(0,1)
            sage: G.has_arc(0,1)
            False

        """
        if u < 0 or u >= self.num_verts:
            raise RuntimeError("First vertex (%d) is not a vertex of the graph."%u)
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("Second vertex (%d) is not a vertex of the graph."%v)
        self.del_arc_unsafe(u,v)

    cdef int out_neighbors_unsafe(self, int u, int *neighbors, int size):
        """
        Gives all v such that (u, v) is an arc of the graph.

        INPUT:
            u -- integer from 0, ..., n-1, where n is the number of vertices
            neighbors -- must be a pointer to an (allocated) integer array
            size -- the length of the array

        OUTPUT:
            nonnegative integer -- the number of v such that (u, v) is an arc
            -1 -- indicates that the array has been filled with neighbors, but
        there were more

        """
        cdef int place = (u * self.num_longs), num_nbrs = 0
        cdef int i, v = 0
        cdef unsigned long word, data
        for i from 0 <= i < self.num_longs:
            data = self.edges[place + i]
            word = 1
            while word:
                if word & data:
                    if num_nbrs == size:
                        return -1
                    neighbors[num_nbrs] = v
                    num_nbrs += 1
                word = word << 1
                v += 1
        return num_nbrs

    def out_neighbors(self, int u):
        """
        Gives all v such that (u, v) is an arc of the graph.

        INPUT:
            u -- integer from 0, ..., n-1, where n is the number of vertices

        EXAMPLES:
            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(1,2)
            sage: G.add_arc(1,3)
            sage: G.out_neighbors(0)
            [1]
            sage: G.out_neighbors(1)
            [2, 3]

        """
        cdef int i, size, num_nbrs
        if u < 0 or u >= self.num_verts:
            raise RuntimeError("First vertex (%d) is not a vertex of the graph."%u)
        if self.out_degrees[u] == 0:
            return []
        cdef int size = self.out_degrees[u]
        cdef int *neighbors = <int *> sage_malloc(size * sizeof(int))
        if not neighbors:
            raise RuntimeError("Failure allocating memory.")
        num_nbrs = self.out_neighbors_unsafe(u, neighbors, size)
        output = [neighbors[i] for i from 0 <= i < num_nbrs]
        sage_free(neighbors)
        return output

    cdef int in_neighbors_unsafe(self, int v, int *neighbors, int size):
        """
        Gives all u such that (u, v) is an arc of the graph.

        INPUT:
            v -- integer from 0, ..., n-1, where n is the number of vertices
            neighbors -- must be a pointer to an (allocated) integer array
            size -- the length of the array

        OUTPUT:
            nonnegative integer -- the number of u such that (u, v) is an arc
            -1 -- indicates that the array has been filled with neighbors, but
        there were more

        """
        cdef int place = v >> self.radix_div_shift
        cdef unsigned long word = (<unsigned long>1) << (v & self.radix_mod_mask)
        cdef int i, num_nbrs = 0
        for i from 0 <= i < self.num_verts:
            if self.edges[place + i*self.num_longs] & word:
                if num_nbrs == size:
                    return -1
                neighbors[num_nbrs] = i
                num_nbrs += 1
        return num_nbrs

    def in_neighbors(self, int v):
        """
        Gives all u such that (u, v) is an arc of the graph.

        INPUT:
            v -- integer from 0, ..., n-1, where n is the number of vertices

        EXAMPLES:
            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(3,1)
            sage: G.add_arc(1,3)
            sage: G.in_neighbors(1)
            [0, 3]
            sage: G.in_neighbors(3)
            [1]

        """
        cdef int i, size, num_nbrs
        if v < 0 or v >= self.num_verts:
            raise RuntimeError("First vertex (%d) is not a vertex of the graph."%v)
        if self.in_degrees[v] == 0:
            return []
        cdef int size = self.in_degrees[v]
        cdef int *neighbors = <int *> sage_malloc(size * sizeof(int))
        if not neighbors:
            raise RuntimeError("Failure allocating memory.")
        num_nbrs = self.in_neighbors_unsafe(v, neighbors, size)
        output = [neighbors[i] for i from 0 <= i < num_nbrs]
        sage_free(neighbors)
        return output


def random_stress():
    """
    Randomly search for mistakes in the code.

    EXAMPLE:
    No output indicates that no errors were found.
        sage: from sage.graphs.base.dense_graph import random_stress
        sage: for _ in xrange(400):
        ...    random_stress()

    """
    cdef int i, j, k, l, n
    cdef DenseGraph Gnew
    num_verts = 10
    from random import randint
    from sage.graphs.graph import DiGraph
    from sage.misc.misc import uniq
    Gnew = DenseGraph(num_verts)
    Gold = DiGraph(num_verts, loops=True, implementation='networkx')
    for n from 0 <= n < 100:
        i = randint(0,num_verts-1)
        j = randint(0,num_verts-1)
        k = randint(0,num_verts-1)
        if k != 0:
            Gold.add_edge(i,j)
            Gnew.add_arc_unsafe(i,j)
        else:
            Gold.delete_edge(i,j)
            Gnew.del_arc_unsafe(i,j)
    if Gnew.num_arcs != Gold.size():
        raise RuntimeError( "NO" )
    for i from 0 <= i < num_verts:
        if Gnew.out_degrees[i] != Gold.out_degree(i):
            raise RuntimeError( "NO" )
        if Gnew.in_degrees[i] != Gold.in_degree(i):
            raise RuntimeError( "NO" )


class DenseGraphBackend(CGraphBackend):

    def __init__(self, n):
        """
        Initialize a dense graph with n vertices.

        EXAMPLE:
            sage: import sage.graphs.base.dense_graph
            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edge(0,1,None,False)
            sage: list(D.iterator_edges(range(9), True, True))
            [(0, 1, None), (1, 0, None)]

        """
        self._cg = DenseGraph(n)

    def add_edge(self, u, v, l, directed):
        """
        Adds the edge u,v to self.

        EXAMPLE:
            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edge(0,1,None,False)
            sage: list(D.iterator_edges(range(9), True, True))
            [(0, 1, None), (1, 0, None)]

        """
        if directed:
            self._cg.add_arc(u, v)
        else:
            self._cg.add_arc(u, v)
            self._cg.add_arc(v, u)

    def add_edges(self, edges, directed):
        """
        Add edges from a list.

        EXAMPLE:
            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True, True))
            [(0, 1, None),
             (1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]

        """
        for e in edges:
            u,v = e[:2]
            self.add_edge(u,v,None,directed)

    def add_vertex(self, name):
        raise NotImplementedError("Dense graphs do not support adding vertices.")

    def add_vertices(self, vertices):
        raise NotImplementedError("Dense graphs do not support adding vertices.")

    def del_edge(self, u, v, l, directed):
        """
        Delete edge u,v.

        EXAMPLE:
            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True, True))
            [(0, 1, None),
             (1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]
            sage: D.del_edge(0,1,None,True)
            sage: list(D.iterator_edges(range(9), True, True))
            [(1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]

        """
        if v is None:
            u, v = u[:2]
        if directed:
            self._cg.del_arc(u, v)
        else:
            self._cg.del_arc(u, v)
            self._cg.del_arc(v, u)

    def del_vertex(self, v):
        raise NotImplementedError("Dense graphs do not support deleting vertices.")

    def del_vertices(self, vertices):
        raise NotImplementedError("Dense graphs do not support deleting vertices.")

    def get_edge_label(self, u, v):
        """
        Returns the edge label for u,v. Always None, since dense graphs don't
        support edge labels.

        EXAMPLE:
            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True, True))
            [(0, 1, None),
             (1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]
            sage: D.del_edge(0,1,None,True)
            sage: list(D.iterator_edges(range(9), True, True))
            [(1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]
            sage: D.get_edge_label(1,0)

        """
        if not self.has_edge(u, v, None):
            raise RuntimeError("%s, %s not an edge of the graph."%(u, v))
        return None

    def has_edge(self, u, v, l):
        """
        Returns whether this graph has edge u,v.

        EXAMPLE:
            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: D.has_edge(0,1,None)
            True

        """
        return self._cg.has_arc(u, v)

    def iterator_edges(self, vertices, labels, not_directed):
        """
        Iterate over the edges incident to a sequence of vertices.

        INPUT:
            vertices:     a list of vertex labels
            labels:       boolean
            not_directed: boolean

        OUTPUT:
            a generator which yields edges, with or without labels
            depending on the labels parameter.

        DOCTEST:
            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.iterator_edges([],True,True)
            <listiterator object at ...>
        """
        if not_directed and labels:
            return iter([(v,u,None) for v in vertices for u in self._cg.out_neighbors(v)])
        elif not_directed:
            return iter([(v,u) for v in vertices for u in self._cg.out_neighbors(v)])
        elif labels:
            return iter([(v,u,None) for v in vertices for u in self._cg.out_neighbors(v) if u >= v or u not in vertices])
        else:
            return iter([(v,u) for v in vertices for u in self._cg.out_neighbors(v) if u >= v or u not in vertices])

    def iterator_in_edges(self, vertices, labels):
        """
        Iterate over the incoming edges incident to a sequence of vertices.

        INPUT:
            vertices:     a list of vertex labels
            labels:       boolean

        OUTPUT:
            a generator which yields edges, with or without labels
            depending on the labels parameter.

        DOCTEST:
            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.iterator_in_edges([],True)
            <listiterator object at ...>
        """
        if labels:
            return iter([(u,v,None) for v in vertices for u in self._cg.in_neighbors(v)])
        else:
            return iter([(u,v) for v in vertices for u in self._cg.in_neighbors(v)])

    def iterator_out_edges(self, vertices, labels):
        """
        Iterate over the outbound edges incident to a sequence of vertices.

        INPUT:
            vertices:     a list of vertex labels
            labels:       boolean

        OUTPUT:
            a generator which yields edges, with or without labels
            depending on the labels parameter.

        DOCTEST:
            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.iterator_out_edges([],True)
            <listiterator object at ...>
        """
        if labels:
            return iter([(v,u,None) for v in vertices for u in self._cg.out_neighbors(v)])
        else:
            return iter([(v,u) for v in vertices for u in self._cg.out_neighbors(v)])

    def multiple_edges(self, new):
        """
        Get/set whether or not self allows multiple edges.

        INPUT:
            new: boolean or None

        DOCTEST:
            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.multiple_edges(True)
            Traceback (most recent call last):
            ...
            NotImplementedError: Dense graphs do not support multiple edges.
            sage: G.multiple_edges(None)
            False
        """
        if new is None:
            return False
        if new:
            raise NotImplementedError("Dense graphs do not support multiple edges.")

    def set_edge_label(self, u, v, l, directed):
        """
        Label the edge (u,v) by l.

        INPUT:
            u,v:      vertices
            l:        edge label
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.set_edge_label(1,2,'a',True)
            Traceback (most recent call last):
            ...
            NotImplementedError: Dense graphs do not support edge labels.
        """
        raise NotImplementedError("Dense graphs do not support edge labels.")


