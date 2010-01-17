"""

Fast compiled graphs

This implements the base class for sparse and dense graphs in Sage. It is not
intended for use on its own.

Data Structure
--------------

The class ``CGraph`` contains the following variables::

        cdef int num_verts
        cdef int num_arcs
        cdef int *in_degrees
        cdef int *out_degrees
        cdef bitset_t active_vertices

The bitset ``active_vertices`` is a list of all available vertices for use, but
only the ones which are set are considered to actually be in the graph. The
variables ``num_verts`` and ``num_arcs`` are self-explanatory (note that
``num_verts`` is the number of bits set in ``active_vertices``, not the full
length of the bitset). The arrays ``in_degrees`` and ``out_degrees`` are of the
same length as the bitset.

For more about active vertices, see the documentation for the ``realloc``
method.

Methods
-------

"""
#*******************************************************************************
#        Copyright (C) 2008-9 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

include '../../misc/bitset.pxi'

from graph_backends import GenericGraphBackend
from sage.rings.integer import Integer

cdef class CGraph:
    """
    Compiled sparse and dense graphs.
    """

    ###################################
    # Vertex Functions
    ###################################

    cpdef bint has_vertex(self, int n):
        """
        Return whether ``n`` is in self.

        INPUT:
         - ``n`` - integer

        EXAMPLES:

        Upon initialization, a SparseGraph or DenseGraph has the first
        ``nverts`` vertices::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts = 10, expected_degree = 3, extra_vertices = 10)
            sage: S.has_vertex(6)
            True
            sage: S.has_vertex(12)
            False
            sage: S.has_vertex(24)
            False
            sage: S.has_vertex(-19)
            False

        ::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: D = DenseGraph(nverts = 10, extra_vertices = 10)
            sage: D.has_vertex(6)
            True
            sage: D.has_vertex(12)
            False
            sage: D.has_vertex(24)
            False
            sage: D.has_vertex(-19)
            False

        """
        return n >= 0 and n < self.active_vertices.size and bitset_in(self.active_vertices, n)

    cpdef check_vertex(self, int n):
        """
        If ``n`` is not in self, raise an error.

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts = 10, expected_degree = 3, extra_vertices = 10)
            sage: S.check_vertex(4)
            sage: S.check_vertex(12)
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (12) is not a vertex of the graph.
            sage: S.check_vertex(24)
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (24) is not a vertex of the graph.
            sage: S.check_vertex(-19)
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (-19) is not a vertex of the graph.

        ::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: D = DenseGraph(nverts = 10, extra_vertices = 10)
            sage: D.check_vertex(4)
            sage: D.check_vertex(12)
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (12) is not a vertex of the graph.
            sage: D.check_vertex(24)
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (24) is not a vertex of the graph.
            sage: D.check_vertex(-19)
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (-19) is not a vertex of the graph.


        """
        if not self.has_vertex(n):
            raise RuntimeError("Vertex (%d) is not a vertex of the graph."%n)

    cdef int add_vertex_unsafe(self, int k):
        """
        Adds the vertex k to the graph.

        INPUT:

            k -- nonnegative integer, or -1
                -1 -- function will find first available vertex

        OUTPUT:

            -1 -- indicates that no vertex was added because the
                current allocation is already full, or the vertex is out of
                range

            nonnegative integer -- this vertex is now guaranteed to be
                in the graph

        """
        if k == -1:
            k = bitset_first_in_complement(self.active_vertices)
        elif self.active_vertices.size <= k:
            k = -1
        if k != -1:
            if not bitset_in(self.active_vertices, k):
                self.num_verts += 1
            bitset_add(self.active_vertices, k)
        return k

    def add_vertex(self, int k = -1):
        """
        Adds vertex ``k`` to the graph. If ``k == -1``, a new vertex is added
        and the integer used is returned.

        INPUT:

         - ``k`` -- non-negative integer, or -1

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(3, extra_vertices=3)
            sage: G.add_vertex(3)
            3
            sage: G.add_arc(2,5)
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (5) is not a vertex of the graph.
            sage: G.add_arc(1,3)
            sage: G.has_arc(1,3)
            True
            sage: G.has_arc(2,3)
            False

        ::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(3, extra_vertices=3)
            sage: G.add_vertex(3)
            3
            sage: G.add_arc(2,5)
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (5) is not a vertex of the graph.
            sage: G.add_arc(1,3)
            sage: G.has_arc(1,3)
            True
            sage: G.has_arc(2,3)
            False

        """
        if k >= 2*self.active_vertices.size:
            raise RuntimeError("Requested vertex is past twice the allocated "\
            "range: use realloc.")
        if k >= self.active_vertices.size or (k==-1 and self.active_vertices.size == self.num_verts):
            self.realloc(2*self.active_vertices.size)
        return self.add_vertex_unsafe(k)

    cpdef add_vertices(self, object verts):
        """
        Adds vertices from the iterable ``verts``.

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts=4, extra_vertices=4)
            sage: S.verts()
            [0, 1, 2, 3]
            sage: S.add_vertices([3,5,7,9])
            sage: S.verts()
            [0, 1, 2, 3, 5, 7, 9]
            sage: S.realloc(20)
            sage: S.verts()
            [0, 1, 2, 3, 5, 7, 9]

        ::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: D = DenseGraph(nverts=4, extra_vertices=4)
            sage: D.verts()
            [0, 1, 2, 3]
            sage: D.add_vertices([3,5,7,9])
            sage: D.verts()
            [0, 1, 2, 3, 5, 7, 9]
            sage: D.realloc(20)
            sage: D.verts()
            [0, 1, 2, 3, 5, 7, 9]

        """
        cdef int v
        for v in verts:
            self.add_vertex(v)

    cdef int del_vertex_unsafe(self, int v):
        """
        Deletes the vertex v, along with all edges incident to it.

        INPUT:
            v -- non-negative integer

        """
        cdef int size = 0, num_nbrs, i, *neighbors
        if self.in_degrees[v] > size:
            size = self.in_degrees[v]
        if self.out_degrees[v] > size:
            size = self.out_degrees[v]
        if size > 0:
            neighbors = <int *> sage_malloc(size * sizeof(int))
            if not neighbors:
                raise RuntimeError("Failure allocating memory.")
            # delete each arc incident with v
            num_nbrs = self.in_neighbors_unsafe(v, neighbors, size)
            for i from 0 <= i < num_nbrs:
                self.del_arc_unsafe(neighbors[i], v)
            num_nbrs = self.out_neighbors_unsafe(v, neighbors, size)
            for i from 0 <= i < num_nbrs:
                self.del_arc_unsafe(v, neighbors[i])
            sage_free(neighbors)

        self.num_verts -= 1
        bitset_remove(self.active_vertices, v)

    cpdef del_vertex(self, int v):
        """
        Deletes the vertex ``v``, along with all edges incident to it. If ``v``
        is not in self, fails silently.

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(3)
            sage: G.add_arc(0,1)
            sage: G.add_arc(0,2)
            sage: G.add_arc(1,2)
            sage: G.add_arc(2,0)
            sage: G.del_vertex(2)
            sage: for i in range(2):
            ...    for j in range(2):
            ...        if G.has_arc(i,j):
            ...            print i,j
            0 1
            sage: G = SparseGraph(3)
            sage: G.add_arc(0,1)
            sage: G.add_arc(0,2)
            sage: G.add_arc(1,2)
            sage: G.add_arc(2,0)
            sage: G.del_vertex(1)
            sage: for i in xrange(3):
            ...    for j in xrange(3):
            ...        if G.has_arc(i,j):
            ...            print i,j
            0 2
            2 0

        """
        if self.has_vertex(v):
            self.del_vertex_unsafe(v)

    cpdef int current_allocation(self):
        """
        Report the number of vertices allocated.

        INPUT:

         - ``total`` - integer, the total size to make the array

        Returns -1 and fails if reallocation would destroy any active vertices.

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts=4, extra_vertices=4)
            sage: S.current_allocation()
            8
            sage: S.add_vertex(6)
            6
            sage: S.current_allocation()
            8
            sage: S.add_vertex(10)
            10
            sage: S.current_allocation()
            16
            sage: S.add_vertex(40)
            Traceback (most recent call last):
            ...
            RuntimeError: Requested vertex is past twice the allocated range: use realloc.
            sage: S.realloc(50)
            sage: S.add_vertex(40)
            40
            sage: S.current_allocation()
            50
            sage: S.realloc(30)
            -1
            sage: S.current_allocation()
            50
            sage: S.del_vertex(40)
            sage: S.realloc(30)
            sage: S.current_allocation()
            30

        """
        return self.active_vertices.size

    cpdef list verts(self):
        """
        Returns a list of the vertices in self.

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts=4, extra_vertices=4)
            sage: S.verts()
            [0, 1, 2, 3]
            sage: S.add_vertices([3,5,7,9])
            sage: S.verts()
            [0, 1, 2, 3, 5, 7, 9]
            sage: S.realloc(20)
            sage: S.verts()
            [0, 1, 2, 3, 5, 7, 9]

        """
        cdef int i
        return [i for i from 0 <= i < self.active_vertices.size if bitset_in(self.active_vertices, i)]

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

    cpdef realloc(self, int total):
        """
        Reallocate the number of vertices to use, without actually adding any.

        INPUT:

         - ``total`` - integer, the total size to make the array

        Returns -1 and fails if reallocation would destroy any active vertices.

        EXAMPLES:

        First, note that ``realloc`` is implemented for ``SparseGraph`` and
        ``DenseGraph`` differently, and is not implemented at the ``CGraph``
        level::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.realloc(20)
            Traceback (most recent call last):
            ...
            NotImplementedError

        ::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts=4, extra_vertices=4)
            sage: S.current_allocation()
            8
            sage: S.add_vertex(6)
            6
            sage: S.current_allocation()
            8
            sage: S.add_vertex(10)
            10
            sage: S.current_allocation()
            16
            sage: S.add_vertex(40)
            Traceback (most recent call last):
            ...
            RuntimeError: Requested vertex is past twice the allocated range: use realloc.
            sage: S.realloc(50)
            sage: S.add_vertex(40)
            40
            sage: S.current_allocation()
            50
            sage: S.realloc(30)
            -1
            sage: S.current_allocation()
            50
            sage: S.del_vertex(40)
            sage: S.realloc(30)
            sage: S.current_allocation()
            30

        ::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: D = DenseGraph(nverts=4, extra_vertices=4)
            sage: D.current_allocation()
            8
            sage: D.add_vertex(6)
            6
            sage: D.current_allocation()
            8
            sage: D.add_vertex(10)
            10
            sage: D.current_allocation()
            16
            sage: D.add_vertex(40)
            Traceback (most recent call last):
            ...
            RuntimeError: Requested vertex is past twice the allocated range: use realloc.
            sage: D.realloc(50)
            sage: D.add_vertex(40)
            40
            sage: D.current_allocation()
            50
            sage: D.realloc(30)
            -1
            sage: D.current_allocation()
            50
            sage: D.del_vertex(40)
            sage: D.realloc(30)
            sage: D.current_allocation()
            30

        """
        raise NotImplementedError()

    cpdef add_arc(self, int u, int v):
        """
        This function is implemented at the level of sparse and dense graphs::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.add_arc(0,1)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError()
    cpdef bint has_arc(self, int u, int v) except -1:
        """
        This function is implemented at the level of sparse and dense graphs::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.has_arc(0,1)
            Not Implemented!
            False

        """
        # The following is due to a hard to reproduce bug in Cython where except,
        # cpdef, and classes don't play well together:
        print "Not Implemented!"
        # raise NotImplementedError() ... results in:
        # Exception exceptions.NotImplementedError: NotImplementedError() in 'sage.graphs.base.c_graph.CGraph.has_arc' ignored
        # False
    cpdef del_all_arcs(self, int u, int v):
        """
        This function is implemented at the level of sparse and dense graphs::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.del_all_arcs(0,1)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError()
    cpdef list all_arcs(self, int u, int v):
        """
        This function is implemented at the level of sparse and dense graphs::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.all_arcs(0,1)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError()

    cpdef list in_neighbors(self, int v):
        """
        This function is implemented at the level of sparse and dense graphs::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.in_neighbors(0)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError()
    cpdef list out_neighbors(self, int u):
        """
        This function is implemented at the level of sparse and dense graphs::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.out_neighbors(0)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError()


    def _in_degree(self, int v):
        """
        Return the number of edges coming into v.

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: SparseGraph(7)._in_degree(3)
            0

        TEST::

            sage: g = Graph({1: [2, 5], 2: [1, 5, 3, 4], 3: [2, 5], 4: [3], 5: [2, 3]}, implementation='c_graph')
            sage: g._backend.degree(5, False)
            3

        """
        if not self.has_vertex(v):
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
        if not self.has_vertex(v):
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

cdef int get_vertex(object u, dict vertex_ints, dict vertex_labels,
                      CGraph G) except ? -2:
    """
    Returns an int representing the arbitrary hashable vertex u (whether or not
    u is actually in the graph), or -1 if a new association must be made for u
    to be a vertex.
    """
    cdef int u_int
    if u in vertex_ints:
        return vertex_ints[u]
    try:
        u_int = u
    except TypeError:
        return -1
    if u_int < 0 or u_int >= G.active_vertices.size:
        return -1
    if u_int in vertex_labels:
        return -1
    return u_int

cdef object vertex_label(int u_int, dict vertex_ints, dict vertex_labels,
                      CGraph G):
    """
    Returns the object represented by u_int, or None if this does not represent
    a vertex.
    """
    if u_int in vertex_labels:
        return vertex_labels[u_int]
    elif bitset_in(G.active_vertices, u_int):
        return u_int
    else:
        return None

cdef int check_vertex(object u, dict vertex_ints, dict vertex_labels,
                      CGraph G, CGraph G_rev, bint reverse) except ? -1:
    """
    Returns an int representing the arbitrary hashable vertex u, and updates, if
    necessary, the translation dict and list. Adds a vertex if the label is new.
    """
    cdef int u_int = get_vertex(u, vertex_ints, vertex_labels, G)
    if u_int != -1:
        if not bitset_in(G.active_vertices, u_int):
            bitset_add(G.active_vertices, u_int)
            G.num_verts += 1
            if reverse:
                bitset_add(G_rev.active_vertices, u_int)
                G_rev.num_verts += 1
        return u_int
    u_int = bitset_first_in_complement(G.active_vertices)
    if u_int == -1:
        G.realloc(2*G.active_vertices.size)
        if reverse:
            G_rev.realloc(2*G_rev.active_vertices.size)
        return check_vertex(u, vertex_ints, vertex_labels, G, G_rev, reverse)
    vertex_labels[u_int] = u
    vertex_ints[u] = u_int
    G.add_vertex(u_int)
    if reverse:
        G_rev.add_vertex(u_int)
    return u_int

class CGraphBackend(GenericGraphBackend):
    """
    Base class for sparse and dense graph backends.

    ::

        sage: from sage.graphs.base.c_graph import CGraphBackend

    This class is extended by ``SparseGraphBackend`` and ``DenseGraphBackend``,
    which are fully functional backends. This class is mainly just for vertex
    functions, which are the same for both. A ``CGraphBackend`` will not work on
    its own::

        sage: from sage.graphs.base.c_graph import CGraphBackend
        sage: CGB = CGraphBackend()
        sage: CGB.degree(0,True)
        Traceback (most recent call last):
        ...
        AttributeError: 'CGraphBackend' object has no attribute 'vertex_ints'

    The appropriate way to use these backends is via Sage graphs::

        sage: G = Graph(30, implementation="c_graph")
        sage: G.add_edges([(0,1), (0,3), (4,5), (9, 23)])
        sage: G.edges(labels=False)
        [(0, 1), (0, 3), (4, 5), (9, 23)]

    """

    _cg = None
    _cg_rev = None
    _directed = None

    def has_vertex(self, v):
        """
        Returns whether ``v`` is a vertex of self.

        INPUT:

         - ``v`` - any object

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraphBackend
            sage: B = SparseGraphBackend(7)
            sage: B.has_vertex(6)
            True
            sage: B.has_vertex(7)
            False

        """
        cdef v_int = get_vertex(v, self.vertex_ints, self.vertex_labels, self._cg)
        if v_int == -1:
            return False
        if not bitset_in((<CGraph>self._cg).active_vertices, v_int):
            return False
        return True

    def degree(self, v, directed):
        """
        Return the degree of the vertex ``v``.

        INPUT:

         - ``v`` - a vertex of the graph

        EXAMPLE::

            sage: from sage.graphs.base.sparse_graph import SparseGraphBackend
            sage: B = SparseGraphBackend(7)
            sage: B.degree(3, False)
            0

        """
        cdef v_int = get_vertex(v, self.vertex_ints, self.vertex_labels, self._cg)
        if directed:
            return self._cg._in_degree(v_int) + self._cg._out_degree(v_int)
        else:
            return self._cg._out_degree(v_int)

    def add_vertex(self, object name):
        """
        Add a vertex to self.

        INPUT:

         - ``name`` - the vertex to be added (must be hashable)

        EXAMPLE::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_vertex(10)
            sage: D.add_vertex([])
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'list'

        ::

            sage: S = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: S.add_vertex(10)
            sage: S.add_vertex([])
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'list'

        """
        if name is None:
            name = 0
            while name in self.vertex_ints or (
            name not in self.vertex_labels and bitset_in((<CGraph>self._cg).active_vertices, name)):
                name += 1
        check_vertex(name, self.vertex_ints, self.vertex_labels, self._cg, self._cg_rev, (self._directed and self._cg_rev is not None)) # this will add the vertex

    def add_vertices(self, object vertices):
        """
        Add vertices to self.

        INPUT:

         - ``vertices`` - iterator of vertex labels

        EXAMPLE::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(1)
            sage: D.add_vertices([1,2,3])

        ::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(0)
            sage: G.add_vertices([0,1])
            sage: list(G.iterator_verts(None))
            [0, 1]
            sage: list(G.iterator_edges([0,1], True))
            []

        ::

            sage: import sage.graphs.base.dense_graph
            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_vertices([10,11,12])

        """
        cdef object v
        for v in vertices:
            self.add_vertex(v)

    def del_vertex(self, v):
        """
        Delete a vertex in self, failing silently if the vertex is not in the
        graph.

        INPUT:

         - ``v`` - vertex to be deleted

        EXAMPLE::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.del_vertex(0)
            sage: D.has_vertex(0)
            False

        ::

            sage: S = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: S.del_vertex(0)
            sage: S.has_vertex(0)
            False

        """
        if not self.has_vertex(v):
            return
        cdef int v_int = get_vertex(v, self.vertex_ints, self.vertex_labels, self._cg)

        # delete each arc incident with v, and v
        self._cg.del_vertex(v_int)
        if self._cg_rev is not None:
            self._cg_rev.del_vertex(v_int)

        # add v to unused vertices
        if v_int in self.vertex_labels:
            self.vertex_ints.pop(v)
            self.vertex_labels.pop(v_int)

    def del_vertices(self, vertices):
        """
        Delete vertices from an iterable container.

        INPUT:

         - ``vertices`` - iterator of vertex labels

        EXAMPLE::

            sage: import sage.graphs.base.dense_graph
            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.del_vertices([7,8])
            sage: D.has_vertex(7)
            False
            sage: D.has_vertex(6)
            True

        ::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.del_vertices([1,2,3])
            sage: D.has_vertex(1)
            False
            sage: D.has_vertex(0)
            True

        """
        cdef object v
        for v in vertices:
            self.del_vertex(v)

    def iterator_nbrs(self, v):
        """
        Returns an iterator over the neighbors of ``v``.

        INPUT:

         - ``v`` - a vertex

        EXAMPLE::

            sage: P = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: list(P._backend.iterator_nbrs(0))
            [1, 4, 5]

        """
        return iter(set(self.iterator_in_nbrs(v)) | set(self.iterator_out_nbrs(v)))

    def iterator_in_nbrs(self, v):
        """
        Returns an iterator over the incoming neighbors of ``v``.

        INPUT:

         - ``v`` - a vertex

        EXAMPLE::

            sage: P = DiGraph(graphs.PetersenGraph().to_directed(), implementation='c_graph')
            sage: list(P._backend.iterator_in_nbrs(0))
            [1, 4, 5]

        """
        cdef int u_int
        cdef int v_int = get_vertex(v, self.vertex_ints, self.vertex_labels, self._cg)
        if self._cg_rev is not None: # Sparse
            return iter([vertex_label(u_int, self.vertex_ints, self.vertex_labels, self._cg)
                     for u_int in self._cg_rev.out_neighbors(v_int)])
        else: # Dense
            return iter([vertex_label(u_int, self.vertex_ints, self.vertex_labels, self._cg)
                     for u_int in self._cg.in_neighbors(v_int)])

    def iterator_out_nbrs(self, v):
        """
        Returns an iterator over the outgoing neighbors of ``v``.

        INPUT:

         - ``v`` - a vertex

        EXAMPLE::

            sage: P = DiGraph(graphs.PetersenGraph().to_directed(), implementation='c_graph')
            sage: list(P._backend.iterator_out_nbrs(0))
            [1, 4, 5]

        """
        cdef u_int
        cdef int v_int = get_vertex(v, self.vertex_ints, self.vertex_labels, self._cg)
        return iter([vertex_label(u_int, self.vertex_ints, self.vertex_labels, self._cg)
                     for u_int in self._cg.out_neighbors(v_int)])

    def iterator_verts(self, verts):
        """
        Returns an iterator over the vertices of self intersected with ``verts``.

        INPUT:
         - ``verts`` - an iterable container of objects

        EXAMPLE::

            sage: P = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: list(P._backend.iterator_verts(P))
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        """
        cdef int i
        cdef object v
        if verts is None:
            S = set(self.vertex_ints.iterkeys())
            for i from 0 <= i < (<CGraph>self._cg).active_vertices.size:
                if i not in self.vertex_labels and bitset_in((<CGraph>self._cg).active_vertices, i):
                    S.add(i)
            return iter(S)
        is_hashable = False
        try:
            v = hash(verts)
            is_hashable = True
        except:
            pass
        if is_hashable and self.has_vertex(verts):
            return iter([verts])
        else:
            L = []
            for v in verts:
                if self.has_vertex(v):
                    L.append(v)
            return iter(L)

    def loops(self, new):
        """
        Returns whether loops are allowed in this graph.

        INPUT:

         - ``new`` - boolean (to set) or ``None`` (to get)

        EXAMPLE::

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

        INPUT:

         - ``new`` - boolean (to set) or ``None`` (to get)

        EXAMPLE::

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

        INPUT:

         - ``directed`` - whether to count ``(u,v)`` and ``(v,u)`` as one or two edges

        EXAMPLE::

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
                if self.multiple_edges(None):
                    for j in xrange(self.num_verts()):
                        if self.has_edge(j,j,None):
                            k += len(self.get_edge_label(j, j))
                else:
                    for j in xrange(self.num_verts()):
                        if self.has_edge(j,j,None):
                            k += 1
            i = (i-k)/2
            return i + k

    def num_verts(self):
        """
        Returns the number of vertices in self.

        EXAMPLE::

            sage: G = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: G._backend.num_verts()
            10

        """
        return (<CGraph>self._cg).num_verts

    def relabel(self, perm, directed):
        """
        Relabels the graph according to perm.

        INPUT:
         - ``perm`` - anything which represents a permutation as ``v --> perm[v]``, for example a dict or a list
         - ``directed`` - ignored (this is here for compatibility with other backends)

        EXAMPLE::

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
        cdef int i
        cdef object v
        cdef dict new_vx_ints = {}
        cdef dict new_vx_labels = {}
        for v in self.iterator_verts(None):
            i = get_vertex(v, self.vertex_ints, self.vertex_labels, self._cg)
            new_vx_ints[perm[v]] = i
            new_vx_labels[i] = perm[v]
        self.vertex_ints = new_vx_ints
        self.vertex_labels = new_vx_labels


    def shortest_path(self,x,y):
        r"""
        Returns the shortest path between x and y

        EXAMPLE::

            sage: G = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: G.shortest_path(0,1)
            [0, 1]
        """

        if x==y:
            return 0

        # The function being mostly symmetric in x and y
        # their roles are reversed at the end of each loop
        # For this reason is defined, for example,
        # two dictionaries dist_y and dist_x containing the
        # distances to x and y, and a dictionary
        # dist_current and dist_other, pointing toward the
        # previous two, alternatively.
        #
        # Besides, there is another difference in the fact
        # that for directed graphs we are interested in paths
        # leaving x toward y, so we are considering the out_neighbors
        # on x's side, and in_neighbors on y's side

        cdef int x_int = get_vertex(x, self.vertex_ints, self.vertex_labels, self._cg)
        cdef int y_int = get_vertex(y, self.vertex_ints, self.vertex_labels, self._cg)
        cdef int u = 0
        cdef int v = 0
        cdef int w = 0

        # Each vertex knows its predecessors in the search, for each side
        cdef dict pred_x = {}
        cdef dict pred_y = {}
        cdef dict pred_current = pred_x
        cdef dict pred_other = pred_y

        # Stores the distances from x and y
        cdef dict dist_x = {}
        cdef dict dist_y = {}
        cdef dict dist_current = dist_x
        cdef dict dist_other = dist_y
        dist_x[x_int] = 0
        dist_y[y_int] = 0

        # Lists of vertices whose neighbors have not been explored yet
        cdef list next_x = [x_int]
        cdef list next_y = [y_int]
        cdef list next_current = next_x
        cdef list next_other = next_y
        cdef list next_temporary = []
        cdef list neighbors

        cdef list shortest_path = []

        # We are interested in edges leaving x and entering y, so we
        # are dealing with two different "neighbors" functions
        cdef int out = 1

        # As long as the current side (x or y) is not totally explored ...
        while next_current:
            next_temporary = []

            # Take the next vertex in the list, and study all of its neighbors
            # When a new neighbor is found, it is added into a temporary list
            # When all the vertices in the list are tested
            # and next_current is replaced by the temporary list
            #
            # After this, current and other are reversed, and the loop restarts
            for u in next_current:
                if out == 1:
                    neighbors = self._cg.out_neighbors(u)
                elif self._cg_rev is not None: # Sparse
                    neighbors = self._cg_rev.out_neighbors(u)
                else: # Dense
                    neighbors = self._cg.in_neighbors(u)
                for v in neighbors:
                    # If the neihgbor is new, updates the distances and adds to the list
                    if not dist_current.has_key(v):
                        dist_current[v] = dist_current[u] + 1
                        pred_current[v] = u
                        next_current.append(v)

                        # If the new neighbor is already known by the other side ...

                        if dist_other.has_key(v):
                            # build the shortest path and returns in.

                            w = v

                            while w != x_int:
                                shortest_path.append(vertex_label(w, self.vertex_ints, self.vertex_labels, self._cg))
                                w = pred_x[w]

                            shortest_path.append(x)
                            shortest_path.reverse()

                            if v == y_int:
                                return shortest_path

                            w=pred_y[v]
                            while w != y_int:
                                shortest_path.append(vertex_label(w, self.vertex_ints, self.vertex_labels, self._cg))
                                w = pred_y[w]
                            shortest_path.append(y)

                            return shortest_path

            next_current = next_temporary
            pred_current, pred_other = pred_other, pred_current
            dist_current, dist_other = dist_other, dist_current
            next_current, next_other = next_other, next_current
            out = -out

        return []

    def bidirectional_dijkstra(self,x,y):
        r"""
        Returns the shortest path between x and y using
        a bidirectional version of Dijkstra

        EXAMPLE::

            sage: G = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: for (u,v) in G.edges(labels=None):
            ...      G.set_edge_label(u,v,1)
            sage: G.shortest_path(0,1,by_weight=True)
            [0, 1]

        TEST:

        Bugfix from #7673 ::

            sage: G = Graph(implementation="networkx")
            sage: G.add_edges([(0,1,9),(0,2,8),(1,2,7)])
            sage: Gc = G.copy(implementation='c_graph')
            sage: sp = G.shortest_path_length(0,1,by_weight=True)
            sage: spc = Gc.shortest_path_length(0,1,by_weight=True)
            sage: sp == spc
            True
        """

        if x==y:
            return 0

        # ****************** WARNING **********************
        # Use Python to maintain a heap...
        # Rewrite this in Cython as soon as possible !
        # *************************************************
        from heapq import heappush, heappop


        # As for shortest_path, the roles of x and y are symmetric, hence we define
        # dictionaries like pred_current and pred_other, which represent alternatively
        # pred_x or pred_y according to the side studied
        cdef int x_int = get_vertex(x, self.vertex_ints, self.vertex_labels, self._cg)
        cdef int y_int = get_vertex(y, self.vertex_ints, self.vertex_labels, self._cg)
        cdef int u = 0
        cdef int v = 0
        cdef int w = 0
        cdef int pred
        cdef float distance
        cdef float edge_label
        cdef int side
        cdef float f_tmp

        # Each vertex knows its predecessors in the search, for each side
        cdef dict pred_x = {}
        cdef dict pred_y = {}
        cdef dict pred_current
        cdef dict pred_other

        # Stores the distances from x and y
        cdef dict dist_x = {}
        cdef dict dist_y = {}
        cdef dict dist_current
        cdef dict dist_other

        # Lists of vertices who are left to be explored
        # they are represented as 4-uples : (distance, side, predecessor ,name)
        # 1 indicates x's side, -1 indicates y's, the distance being
        # defined relatively
        cdef list queue = [(0,1,x_int,x_int),(0,-1,y_int,y_int)]
        cdef list neighbors

        cdef list shortest_path = []

        # meeting_vertex is a vertex discovered through x and through y
        # which defines ther shortest path found
        # ( of length shortest_path_length )
        cdef int meeting_vertex = -1
        cdef float shortest_path_length

        # As long as the current side (x or y) is not totally explored ...
        while queue:

            (distance, side, pred, v) = heappop(queue)
            if meeting_vertex != -1 and distance > shortest_path_length:
                break

            if side == 1:
                dist_current, dist_other = dist_x, dist_y
                pred_current, pred_other = pred_x, pred_y
            else:
                dist_current, dist_other = dist_y, dist_x
                pred_current, pred_other = pred_y, pred_x

            if not dist_current.has_key(v):
                pred_current[v] = pred
                dist_current[v] = distance

                if dist_other.has_key(v):
                    f_tmp = distance + dist_other[v]
                    if meeting_vertex == -1 or f_tmp<shortest_path_length:
                        meeting_vertex = v
                        shortest_path_length = f_tmp

                if side == 1:
                    neighbors = self._cg.out_neighbors(v)
                elif self._cg_rev is not None: # Sparse
                    neighbors = self._cg_rev.out_neighbors(v)
                else: # Dense
                    neighbors = self._cg.in_neighbors(v)
                for w in neighbors:
                    # If the neihgbor is new, adds its non-found neighbors to the queue
                    if not dist_current.has_key(w):
                        edge_label = self.get_edge_label(v,w) if side == 1 else self.get_edge_label(w,v)
                        heappush(queue,(distance + edge_label,side,v,w))


        # No meeting point has been found
        if meeting_vertex == -1:
            return []
        else:
            # build the shortest path and returns it.
            w = meeting_vertex

            while w != x_int:
                shortest_path.append(vertex_label(w, self.vertex_ints, self.vertex_labels, self._cg))
                w = pred_x[w]

            shortest_path.append(x)
            shortest_path.reverse()

            if meeting_vertex == y_int:
                return shortest_path

            w=pred_y[meeting_vertex]
            while w != y_int:
                shortest_path.append(vertex_label(w, self.vertex_ints, self.vertex_labels, self._cg))
                w = pred_y[w]
            shortest_path.append(y)

            return shortest_path

    def shortest_path_all_vertices(self, v, cutoff = None):
        r"""
        Returns for each vertex `u` a shortest  `v-u` path.

        INPUT:

        - ``v`` -- a vertex
        - ``cutoff`` -- maximal distance. Longer paths will not be returned

        OUTPUT:

        A list which associates to each vertex `u` the shortest path between
        `u` and `v` if there is one.

        NOTE:

        - The weight of edges is not taken into account.

        ALGORITHM:

        This is just a breadth-first search.

        EXAMPLES:

        On the Petersen Graph::

            sage: g = graphs.PetersenGraph()
            sage: paths = g._backend.shortest_path_all_vertices(0)
            sage: all([ len(paths[v]) == 0 or len(paths[v])-1 == g.distance(0,v) for v in g])
            True

        On a disconnected graph ::

            sage: g = 2*graphs.RandomGNP(20,.3)
            sage: paths = g._backend.shortest_path_all_vertices(0)
            sage: all([ (not paths.has_key(v) and g.distance(0,v) == +Infinity) or len(paths[v])-1 == g.distance(0,v) for v in g])
            True
        """
        cdef list current_layer
        cdef list next_layer
        cdef bitset_t seen
        cdef int v_int
        cdef int u_int
        cdef dict distances_int
        cdef dict distance
        cdef int d

        distances = {}
        d = 0

        v_int = get_vertex(v, self.vertex_ints, self.vertex_labels, self._cg)

        bitset_init(seen,(<CGraph>self._cg).active_vertices.size)
        bitset_set_first_n(seen,0)
        bitset_add(seen,v_int)


        current_layer = [(u_int, v_int) for u_int in self._cg.out_neighbors(v_int)]
        next_layer = []
        distances[v] = [v]

        while current_layer:
            if cutoff is not None and d >= cutoff:
                break

            while current_layer:
                v_int, u_int = current_layer.pop()

                if bitset_not_in(seen,v_int):
                    bitset_add(seen,v_int)
                    distances[vertex_label(v_int, self.vertex_ints, self.vertex_labels, self._cg)] = distances[vertex_label(u_int, self.vertex_ints, self.vertex_labels, self._cg)] + [vertex_label(v_int, self.vertex_ints, self.vertex_labels, self._cg)]
                    next_layer.extend([(u_int,v_int) for u_int in self._cg.out_neighbors(v_int)])

            current_layer = next_layer
            next_layer = []
            d += 1

        # If the graph is not connected, vertices which have not been
        # seen should be associated to the empty path

        #for 0 <= v_int < (<CGraph>self._cg).active_vertices.size:
        #    if bitset_in((<CGraph>self._cg).active_vertices, v_int) and not bitset_in(seen, v_int):
        #        distances[vertex_label(v_int, self.vertex_ints, self.vertex_labels, self._cg)] = []

        bitset_free(seen)
        return distances


    def depth_first_search(self, v, reverse=False, ignore_direction=False):
        r"""
        Returns a depth-first search from vertex v.

        INPUT:

        - ``reverse`` -- considers the reversed graph (
          the out_neighbors become in_neighbors ).

          The default value is ``False``.

        - ``ignore_direction`` -- consider the undirected graph.

          The default value is ``False``.

        EXAMPLE::

            sage: G = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: list(G.depth_first_search(0))
            [0, 5, 8, 6, 9, 7, 2, 3, 4, 1]
        """
        return Search_iterator(self, v, direction=-1, reverse=reverse, ignore_direction=ignore_direction)

    def breadth_first_search(self, v, reverse=False, ignore_direction=False):
        r"""
        Returns a breadth-first search from vertex v.

        INPUT:

        - ``reverse`` -- considers the reversed graph (
          the out_neighbors become in_neighbors ).

          The default value is ``False``.

        - ``ignore_direction`` -- consider the undirected graph.

          The default value is ``False``.


        EXAMPLE::

            sage: G = Graph(graphs.PetersenGraph(), implementation='c_graph')
            sage: list(G.breadth_first_search(0))
            [0, 1, 4, 5, 2, 6, 3, 9, 7, 8]
        """
        return Search_iterator(self, v, direction=0, reverse=reverse, ignore_direction=ignore_direction)

    def is_connected(self):
        r"""
        Returns whether the graph is connected.

        EXAMPLE:

        Petersen's graph is connected ::

           sage: DiGraph(graphs.PetersenGraph(),implementation="c_graph").is_connected()
           True

        While the disjoint union of two of them is not::

           sage: DiGraph(2*graphs.PetersenGraph(),implementation="c_graph").is_connected()
           False

        A graph with non-integer vertex labels::
            sage: Graph(graphs.CubeGraph(3), implementation='c_graph').is_connected()
            True

        """

        cdef int v_int = 0
        v_int = bitset_first((<CGraph>self._cg).active_vertices)

        if v_int == -1:
            return True
        v = vertex_label(v_int, self.vertex_ints, self.vertex_labels, self._cg)
        return len(list(self.depth_first_search(v, ignore_direction=True)) ) == (<CGraph>self._cg).num_verts

    def is_strongly_connected(self):
        r"""
        Returns whether the graph is strongly connected.

        EXAMPLE:

        The circuit on 3 vertices is obviously strongly connected ::

            sage: g = DiGraph({ 0 : [1], 1 : [2], 2 : [0]},implementation="c_graph")
            sage: g.is_strongly_connected()
            True

        But a transitive triangle is not::

            sage: g = DiGraph({ 0 : [1,2], 1 : [2]},implementation="c_graph")
            sage: g.is_strongly_connected()
            False
        """
        cdef int v_int = 0

        # Pick one vertex
        v_int = bitset_first((<CGraph>self._cg).active_vertices)

        if v_int == -1:
            return True

        v = vertex_label(v_int, self.vertex_ints, self.vertex_labels, self._cg)

        return (<CGraph>self._cg).num_verts == len(list(self.depth_first_search(v))) and \
            (<CGraph>self._cg).num_verts == len(list(self.depth_first_search(v, reverse=True)))

    def strongly_connected_component_containing_vertex(self, v):
        r"""
        Returns the strongly connected component containing the given vertex

        INPUT:

        - ``v`` -- a vertex

        EXAMPLE:

        The digraph obtained from the PetersenGraph has an unique
        strongly connected component ::

            sage: g = DiGraph(graphs.PetersenGraph())
            sage: g.strongly_connected_component_containing_vertex(0)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        In the Butterfly DiGraph, each vertex is a strongly connected
        component::

            sage: g = digraphs.ButterflyGraph(3)
            sage: all([[v] == g.strongly_connected_component_containing_vertex(v) for v in g])
            True
        """
        cdef int v_int = get_vertex(v, self.vertex_ints, self.vertex_labels, self._cg)
        cdef set a = set(self.depth_first_search(v))
        cdef set b = set(self.depth_first_search(v, reverse=True))
        return list(a & b)



cdef class Search_iterator:
    cdef graph
    cdef int direction
    cdef list stack
    cdef bitset_t seen
    cdef bool test_out
    cdef bool test_in

    def __init__(self, graph, v, direction=0, reverse=False, ignore_direction=False):
        self.graph = graph
        self.direction = direction

        bitset_init(self.seen,(<CGraph>self.graph._cg).active_vertices.size)
        bitset_set_first_n(self.seen,0)

        self.stack = [get_vertex(v, self.graph.vertex_ints, self.graph.vertex_labels, self.graph._cg)]

        if not self.graph.directed:
            ignore_direction = False

        self.test_out = (not reverse) or ignore_direction
        self.test_in = reverse or ignore_direction

    def __iter__(self):
        return self

    def __next__(self):
        cdef int v_int
        cdef int w_int

        while self.stack:
            v_int = self.stack.pop(self.direction)

            if bitset_not_in(self.seen,v_int):
                value = vertex_label(v_int, self.graph.vertex_ints, self.graph.vertex_labels, self.graph._cg)
                bitset_add(self.seen,v_int)

                if self.test_out:
                    self.stack.extend(self.graph._cg.out_neighbors(v_int))
                if self.test_in:
                    self.stack.extend(self.graph._cg.in_neighbors(v_int))

                break
        else:
            bitset_free(self.seen)
            raise StopIteration

        return value




