r"""
Graph Theory

AUTHOR:
    -- Robert L. Miller (2006-10-22): initial version
    -- William Stein (2006-12-05): Editing
    -- Robert L. Miller (2006-01-13): refactoring, adjusting for
        NetworkX-0.33, fixed plotting bugs
    -- Robert L. Miller (2006-01-23): basic tutorial, edge labels, loops,
        multiple edges & arcs

TUTORIAL:

    I. The Basics

        1. Graph Format
            SAGE graphs are actually NetworkX graphs, wrapped in a SAGE class.
        In fact, any graph can produce its underlying NetworkX graph as
        follows:

            sage: G = graphs.PetersenGraph()
            sage: N = G.networkx_graph()
            sage: N  # random output
            <networkx.graph.Graph object at 0xa6bf1d0>

        The NetworkX graph is essentially a dictionary of dictionaries:

            sage: N.adj
            {0: {1: None, 4: None, 5: None}, 1: {0: None, 2: None, 6: None}, 2: {1: None, 3: None, 7: None}, 3: {8: None, 2: None, 4: None}, 4: {0: None, 9: None, 3: None}, 5: {0: None, 8: None, 7: None}, 6: {8: None, 1: None, 9: None}, 7: {9: None, 2: None, 5: None}, 8: {3: None, 5: None, 6: None}, 9: {4: None, 6: None, 7: None}}

        Each dictionary key is a vertex label, and each key in the following
        dictionary is a neighbor of that vertex. In undirected graphs, there
        is reduncancy. For example, the dictionary containing the entry
        1: {2: None} implies it must contain 2: {1: None}. The innermost entry
        of None is related to edge labelling (see section I.3.).

        2. Databases

        For some commonly used graphs to play with, type

        graphs.

        and hit <tab>. Most of these graphs come with their own custom plot, so you
        can see how people usually visualize these graphs. Work is currently in progress
        to include a database of known graphs that can be searched by certain
        parameters.

            sage: G = graphs.PetersenGraph()
            sage: G.plot().save('sage.png')    # or G.show()
            sage: G.degree_histogram()
            [0, 0, 0, 10]
            sage: G.adjacency_matrix()
            [0 1 0 0 1 1 0 0 0 0]
            [1 0 1 0 0 0 1 0 0 0]
            [0 1 0 1 0 0 0 1 0 0]
            [0 0 1 0 1 0 0 0 1 0]
            [1 0 0 1 0 0 0 0 0 1]
            [1 0 0 0 0 0 0 1 1 0]
            [0 1 0 0 0 0 0 0 1 1]
            [0 0 1 0 0 1 0 0 0 1]
            [0 0 0 1 0 1 1 0 0 0]
            [0 0 0 0 1 0 1 1 0 0]

            sage: S = G.random_subgraph(.7)
            sage: S.plot().save('sage.png')    # or S.show()
            sage: S.density()         # random output (depends on choice of random graph)
            0.33333333333333331

        3. Labels

        Each vertex can have any hashable object as a label. These are things like
        strings, numbers, and tuples. Each edge is given a default label of None, but
        if specified, edges can have any label at all. Edges between nodes u and v are
        represented typically as (u, v, l), where l is the label for the edge.

NOTE: Many functions are passed directly on to NetworkX, and in this
case the documentation is based on the NetworkX docs.
"""

#*****************************************************************************
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

## IMPORTANT: Do not import networkx at module scope.  It takes a
## surprisingliy long time to initialize itself.  It's better if it is
## imported in functions, so it only gets started if it is actually
## going to be used.

from random import random

from sage.structure.sage_object import SageObject
from sage.plot.plot import Graphics, GraphicPrimitive_NetworkXGraph

class GenericGraph(SageObject):

    def __contains__(self, vertex):
        """
        Return True if vertex is one of the vertices of this graph, i.e.,
        is equal to one of the vertices.

        INPUT:
            vertex -- an integer

        OUTPUT:
            bool -- True or False

        EXAMPLES:
            sage: g = Graph({0:[1,2,3], 2:[5]}); g
            Simple graph on 5 vertices (no loops, no multiple edges)
            sage: 2 in g
            True
            sage: 10 in g
            False
        """
        return vertex in self._nxg

    def __getitem__(self,vertex):
        """
        G[vertex] returns the neighbors (in & out if digraph) of vertex.
        """
        return self.neighbors(vertex)

    def __iter__(self):
        """
        Return an iterator over the vertices. Allows 'for v in G'.
        """
        return self.vertex_iterator()

    def __len__(self):
        return len(self._nxg.adj)

    def __str__(self):
        if self._nxg.name != "No Name":
            return self._nxg.name
        else: return repr(self)

    def _latex_(self):
        return repr(self)

    def _matrix_(self, R=None):
        """
        EXAMPLES:
            sage: G = graphs.CompleteBipartiteGraph(2,3)
            sage: m = matrix(G); m.parent()
            Full MatrixSpace of 5 by 5 sparse matrices over Integer Ring
            sage: m
            [0 0 1 1 1]
            [0 0 1 1 1]
            [1 1 0 0 0]
            [1 1 0 0 0]
            [1 1 0 0 0]
            sage: factor(m.charpoly())
            (x^2 - 6) * x^3
        """
        if R is None:
            return self.am()
        else:
            return self.am().change_ring(R)

    def networkx_graph(self):
        """
        Creates a NetworkX graph from the SAGE graph.
        """
        return self._nxg.copy()

    def networkx_info(self, vertex=None):
        """
        Returns NetworkX information about the graph or the given node.
        """
        self._nxg.info(vertex)

    ### General properties

    def density(self):
        """
        Returns the density (number of edges divided by number of possible
        edges).
        """
        import networkx
        return networkx.density(self._nxg)

    def order(self):
        """
        Returns the number of vertices.
        """
        return self._nxg.order()

    def size(self):
        """
        Returns the number of edges.
        """
        return self._nxg.size()

    ### Vertex handlers

    def add_vertices(self, vertices):
        """
        Add vertices to the graph from an iterable container of vertices.
        """
        self._nxg.add_nodes_from(vertices)

    def clear(self):
        """
        Empties the graph of vertices and edges, removes name.
        """
        self._nxg.clear()

    def neighbors(self, vertex):
        """
        Return a list of neighbors (in and out if directed) of vertex.
        """
        return list(self.neighbor_iterator(vertex))

    def random_subgraph(self, p, inplace=False, create_using=None):
        """
        Return a random subgraph that contains each vertex with prob. p.
        """
        vertices = []
        p = float(p)
        for v in self:
            if random() < p:
                vertices.append(v)
        return self.subgraph(vertices, inplace, create_using)

    def vertex_iterator(self, vertices=None):
        """
        Returns an iterator over the given vertices. Returns False if not given
        a vertex, sequence, iterator or None. None is equivalent to a list of
        every vertex.
        """
        return self._nxg.prepare_nbunch(vertices)

    def vertices(self):
        """
        Return a list of the vertex keys.
        """
        return self._nxg.nodes()

    ### Constructors

    def am(self):
        """
        Shorter call for adjacency matrix makes life easier.
        """
        return self.adjacency_matrix()

class Graph(GenericGraph):
    """
    Undirected graph.

    EXAMPLES:
        sage: g = Graph({0:[1,2,3], 2:[5]}); g
        Simple graph on 5 vertices (no loops, no multiple edges)
        sage: g.vertices()
        [0, 1, 2, 3, 5]
        sage: g.edges()
        [(0, 1, None), (0, 2, None), (0, 3, None), (2, 5, None)]
        sage: g.plot().save('sage.png')
    """

    ### Note: NetworkX function print_dna not wrapped.

    def __init__(self, data=None, pos=None, loops=False, **kwds):
        """
        Create a graph object.

        INPUT:
            data -- can be any of the following:
                1. A NetworkX graph
                2. A dictionary of dictionaries
                3. A dictionary of lists
                4. A numpy matrix or ndarray
                5. A pygraphviz agraph
                6. A scipy sparse matrix

            pos -- a positioning dictionary: for example, the
            spring layout from NetworkX for the 5-cycle is
            {0: [-0.91679746, 0.88169588],
             1: [ 0.47294849, 1.125     ],
             2: [ 1.125     ,-0.12867615],
             3: [ 0.12743933,-1.125     ],
             4: [-1.125     ,-0.50118505]}
            name -- (must be an explicitly named parameter, i.e.,
                     name="complete") gives the graph a name
            loops -- boolean, whether to allow loops (ignored if data is an instance of
            the Graph class)
            multiedges -- boolean, whether to allow multiple edges (ignored if data is
            an instance of the Graph class)

        EXAMPLES:
        We illustrate the first four input formats (the other two
        involve packages that are currently not standard in SAGE):

        1. A networkx graph:
            sage: import networkx
            sage: g = networkx.Graph({0:[1,2,3], 2:[5]})
            sage: Graph(g)
            Simple graph on 5 vertices (no loops, no multiple edges)

        2. A dictionary of dictionaries:
            sage: g = Graph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}); g
            Simple graph on 5 vertices (no loops, no multiple edges)

        The labels ('x', 'z', 'a', 'out') are labels for edges. For example, 'out' is
        the label for the edge on 2 and 5. Labels can be used as weights, if all the
        labels share some common parent.

        3. A dictionary of lists:
            sage: g = Graph({0:[1,2,3], 2:[5]}); g
            Simple graph on 5 vertices (no loops, no multiple edges)

        4. A numpy matrix or ndarray:
            TODO

        Other examples:
            sage: G = Graph(name="Null graph")
            sage: G
            Null graph: a simple graph on 0 vertices (no loops, no multiple edges)
            sage: P = Graph({0:[1,4,5],1:[0,2,6],2:[1,3,7],3:[2,4,8],4:[0,3,9],5:[0,7,8],6:[1,8,9],7:[2,5,9],8:[3,5,6],9:[4,6,7]}, name="Petersen graph")
            sage: P
            Petersen graph: a simple graph on 10 vertices (no loops, no multiple edges)
        """
        import networkx
        if isinstance(data, Graph):
            self._nxg = data.networkx_graph()
        elif isinstance(data, networkx.Graph):
            self._nxg = networkx.XGraph(data, selfloops=loops, **kwds)
        elif isinstance(data, networkx.XGraph):
            self._nxg = data
        else:
            self._nxg = networkx.XGraph(data, selfloops=loops, **kwds)
            if kwds.has_key('name'):
                self._nxg.name = kwds['name']
        self.__pos = pos

    def _repr_(self):
        if not self._nxg.name is None and not self._nxg.name == "":
            name = self._nxg.name
            name = name + ": a s"
        else:
            name = "S"
        if self.loops():
            loops = "with"
        else:
            loops = "no"
        if self.multiple_edges():
            multi = "with"
        else:
            multi = "no"
        return name + "imple graph on %d vertices (%s loops, %s multiple edges)"%(len(self._nxg.adj),loops,multi)

    def copy(self):
        """
        Creates a copy of the graph.
        """
        G = Graph(self._nxg, name=self._nxg.name)
        return G

    def to_directed(self):
        return DiGraph(self._nxg.to_directed(), pos=self.__pos)

    def to_undirected(self):
        return self.copy()

    ### General properties

    def is_directed(self):
        return False

    def loops(self, new=None):
        """
        Returns whether loops are permitted in the graph.

        INPUT:
        new -- boolean, changes whether loops are permitted in the graph.
        """
        if not new is None:
            if new:
                self._nxg.allow_selfloops()
            else:
                self._nxg.ban_selfloops()
        return self._nxg.selfloops

    def multiple_edges(self, new=None):
        """
        Returns whether multiple edges are permitted in the graph.

        INPUT:
        new -- boolean, changes whether multiple edges are permitted in the graph.
        """
        if not new is None:
            if new:
                self._nxg.allow_multiedges()
            else:
                self._nxg.ban_multiedges()
        return self._nxg.multiedges

    ### Vertex handlers

    def add_vertex(self, name=None):
        """
        Creates an isolated vertex

        INPUT:
        name -- Name of the new vertex. If no name is specified, then the
        vertex will be represented by the least integer not already represen-
        ting a vertex. Name must be an immutable object.
        """
        ### TODO- add doc note about representing other objects as vertices
        ### This will be done when such representation is implemented
        if name is None: # then find an integer to use as a key
            i = 0
            while self._nxg.adj.has_key(i):
                i=i+1
            self._nxg.add_node(i)
        else:
            self._nxg.add_node(name)

    def delete_vertex(self, vertex):
        """
        Deletes vertex, removing all incident edges.
        """
        self._nxg.delete_node(vertex)

    def delete_vertices(self, vertices):
        """
        Remove vertices from the graph taken from an iterable container of
        vertices.
        """
        self._nxg.delete_nodes_from(vertices)

    def has_vertex(self, vertex):
        return self._nxg.has_node(vertex)

    def neighbor_iterator(self, vertex):
        """
        Return an iterator over neighbors of vertex.
        """
        return self._nxg.neighbors_iter(vertex)

    def vertex_boundary(self, vertices1, vertices2=None):
        """
        Returns a list of all vertices in the external boundary of vertices1,
        intersected with vertices2. If vertices2 is None, then vertices2 is the
        complement of vertices1.
        """
        return self._nxg.node_boundary(vertices1, vertices2)

    def loop_vertices(self):
        """
        Returns a list of vertices with loops.
        """
        return self._nxg.nodes_with_selfloops()

    ### Edge Handlers

    def add_edge(self, u, v=None, label=None):
        """
        Adds an edge between u and v.

        INPUT:
        The following forms are all accepted:

        G.add_edge( 1, 2 )
        G.add_edge( (1, 2) )
        G.add_edges( [ (1, 2) ] )
        G.add_edge( 1, 2, 'label' )
        G.add_edge( (1, 2, 'label') )
        G.add_edges( [ (1, 2, 'label') ] )

        WARNING:
        The following intuitive input results in nonintuitive output:
        sage.: G = Graph()
        sage.: G.add_edge((1,2), 'label')
        sage.: G.networkx_graph().adj
        {'label': {(1, 2): None}, (1, 2): {'label': None}}

        Use one of these instead:
        sage.: G = Graph()
        sage.: G.add_edge((1,2), label="label")
        sage.: G.networkx_graph().adj
        {1: {2: 'label'}, 2: {1: 'label'}}

        sage.: G = Graph()
        sage.: G.add_edge(1,2,'label')
        sage.: G.networkx_graph().adj
        {1: {2: 'label'}, 2: {1: 'label'}}
        """
        self._nxg.add_edge(u, v, label)

    def add_edges(self, edges):
        """
        Add edges from an iterable container.
        """
        self._nxg.add_edges_from( edges )

    def delete_edge(self, u, v=None, label=None):
        r"""
        Delete the edge \{u, v\}, return silently if vertices or edge does not
        exist.

        INPUT:
        The following forms are all accepted:

        G.delete_edge( 1, 2 )
        G.delete_edge( (1, 2) )
        G.delete_edges( [ (1, 2) ] )
        G.delete_edge( 1, 2, 'label' )
        G.delete_edge( (1, 2, 'label') )
        G.delete_edges( [ (1, 2, 'label') ] )
        """
        self._nxg.delete_edge(u, v, label)

    def delete_edges(self, edges):
        """
        Delete edges from an iterable container.
        """
        self._nxg.delete_edges_from(edges)

    def delete_multiedge(self, u, v):
        """
        Deletes all edges on u and v.
        """
        self._nxg.delete_multiedge(u, v)

    def edges(self, labels=True):
        """
        Return a list of edges. Each edge is a triple (u,v,l) where u and v are
        vertices and l is a label.

        INPUT:
        label -- if False, each edge is a tuple (u,v) of vertices.
        """
        L = self._nxg.edges()
        if labels:
            return L
        else:
            K = []
            for u,v,l in L:
                K.append((u,v))
            return K

    def edge_boundary(self, vertices1, vertices2=None, labels=True):
        """
        Returns a list of edges (u,v,l) with u in vertices1 and v in vertices2.
        If vertices2 is None, then it is set to the complement of vertices1.

        INPUT:
        label -- if False, each edge is a tuple (u,v) of vertices.
        """
        L = self._nxg.edge_boundary(vertices1, vertices2)
        if labels:
            return L
        else:
            K = []
            for u,v,l in L:
                K.append((u,v))
            return K

    def edge_iterator(self, vertices=None):
        """
        Returns an iterator over the edges incident with any vertex given.
        If vertices is None, then returns an iterator over all edges.
        """
        return self._nxg.edges_iter(vertices)

    def edges_incident(self, vertices=None, labels=True):
        """
        Returns a list of edges incident with any vertex given. If vertex is
        None, returns a list of all edges in graph.

        INPUT:
        label -- if False, each edge is a tuple (u,v) of vertices.
        """
        L = self._nxg.edges(vertices)
        if labels:
            return L
        else:
            K = []
            for u,v,l in L:
                K.append((u,v))
            return K

    def has_edge(self, u, v=None, label=None):
        r"""
        Returns True if \{u, v\} is an edge, False otherwise.

        INPUT:
        The following forms are accepted by NetworkX:

        G.has_edge( 1, 2 )
        G.has_edge( (1, 2) )
        G.has_edge( 1, 2, 'label' )
        """
        return self._nxg.has_edge(u, v)

    def edge_label(self, u, v=None, label=None):
        """
        Returns the label of an edge.
        """
        return self._nxg.get_edge(u,v)

    def remove_multiple_edges(self):
        self._nxg.remove_all_multiedges()

    def remove_loops(self, vertices=None):
        """
        Removes loops on vertices in vertices. If vertices is None, removes all loops.
        """
        if vertices is None:
            self._nxg.remove_all_selfloops()
        else:
            for v in vertices:
                self.delete_multiedge(v,v)

    def loop_edges(self):
        """
        Returns a list of all loops in the graph.
        """
        return self._nxg.selfloop_edges()

    def number_of_loops(self):
        """
        Returns the number of edges that are loops.
        """
        return self._nxg.number_of_selfloops()

    ### Degree functions

    def degree(self, vertices=None, with_labels=False):
        """
        Gives the degree of a vertex or of vertices.

        INPUT:
        If vertices is a single vertex, returns the number of neighbors of
        vertex. If vertices is an iterable container of vertices, returns a
        list of degrees. If vertices is None, same as listing all vertices.

        OUTPUT:
        Single vertex- an integer. Multiple vertices- a list of integers. If
        with_labels is True, then returns a dictionary mapping each vertex to
        its degree.

        EXAMPLE:
        sage: P = graphs.PetersenGraph()
        sage: P.degree(5)
        3

        sage: K = graphs.CompleteGraph(9)
        sage: K.degree()
        [8, 8, 8, 8, 8, 8, 8, 8, 8]
        """
        return self._nxg.degree(vertices, with_labels)

    def degree_histogram(self):
        """
        Returns a list, whose ith entry is the frequency of degree i.
        """
        import networkx
        return networkx.degree_histogram(self._nxg)

    def degree_iterator(self, vertices=None, with_labels=False):
        """
        with_labels=False:
            returns an iterator over degrees.
        with_labels=True:
            returns an iterator over tuples (vertex, degree).
        """
        return self._nxg.degree_iter(vertices, with_labels)

    ### Representations

    def adjacency_matrix(self, sparse=True):
        """
        Returns the adjacency matrix of the digraph. Each vertex is
        represented by its position in the list returned by the vertices()
        function.
        """
        n = len(self._nxg.adj)
        verts = self.vertices()
        D = {}
        for e in self.edge_iterator():
            i,j,l = e
            i = verts.index(i)
            j = verts.index(j)
            D[(i,j)] = 1
            D[(j,i)] = 1
        from sage.rings.integer_mod_ring import IntegerModRing
        from sage.matrix.constructor import matrix
        M = matrix(IntegerModRing(2), n, n, D, sparse=sparse)
        return M

    ### Construction

    def add_cycle(self, vertices):
        """
        Adds a cycle to the graph with the given vertices. If the vertices are
        already present, only the edges are added.

        INPUT:
        vertices -- a list of indices for the vertices of the cycle to be
        added.

        EXAMPLES:
        sage: G = Graph()
        sage: for i in range(10): G.add_vertex(name=i)
        sage.: show(G)
        sage: G.add_cycle(range(20)[10:20])
        sage.: show(G)
        sage: G.add_cycle(range(10))
        sage.: show(G)
        """
        self._nxg.add_cycle(vertices)

    def add_path(self, vertices):
        """
        Adds a cycle to the graph with the given vertices. If the vertices are
        already present, only the edges are added.

        INPUT:
        vertices -- a list of indices for the vertices of the cycle to be
        added.

        EXAMPLES:
        sage: G = Graph()
        sage: for i in range(10): G.add_vertex(name=i)
        sage.: show(G)
        sage: G.add_path(range(20)[10:20])
        sage.: show(G)
        sage: G.add_path(range(10))
        sage.: show(G)
        """
        self._nxg.add_path(vertices)

    def subgraph(self, vertices, inplace=False, create_using=None):
        """
        Returns the subgraph induced by the given vertices.

        INPUT:
        inplace -- Using inplace is True will simply delete the extra vertices
        and edges from the current graph. This will modify the graph, and re-
        turn itself.
        vertices -- Vertices can be a single vertex or an iterable container
        of vertices, e.g. a list, set, graph, file or numeric array.
        create_using -- Can be an existing graph object or a call to a graph
        object, such as create_using=DiGraph().
        """
        if inplace:
            self._nxg = self._nxg.subgraph(vertices, inplace, create_using)
            return self
        else:
            NXG = self._nxg.subgraph(vertices, inplace, create_using)
            return Graph(NXG)

    ### Visualization

    def plot(self, pos=None, vertex_labels=True, node_size=200):
        GG = Graphics()
        if pos is None:
            if self.__pos is None:
                NGP = GraphicPrimitive_NetworkXGraph(self._nxg, pos=None, vertex_labels=vertex_labels, node_size=node_size)
            else:
                NGP = GraphicPrimitive_NetworkXGraph(self._nxg, pos=self.__pos, vertex_labels=vertex_labels, node_size=node_size)
        GG.append(NGP)
        GG.axes(False)
        return GG

    def show(self, pos=None, vertex_labels=True, node_size=200, **kwds):
        """
        INPUT:
            pos -- ??
            with_labels -- bool (default: True)
            node_size -- how big the nodes are
            other named options -- All other options are passed onto
            the show command; e.g., dpi=50 will make a small plot.

        EXAMPLES:
            sage: g = Graph({0:[1,2,3], 2:[5]}); g
            Simple graph on 5 vertices (no loops, no multiple edges)
            sage: g.plot().save('sage.png')
        """
        self.plot(pos=pos, vertex_labels=vertex_labels, node_size=node_size).show(**kwds)

class DiGraph(GenericGraph):
    """
    Directed graph (no loops, multiple edges only with distinct labels, non-hyper).
    """

    def __init__(self, data=None, pos=None, loops=False, **kwds):
        """
        Initialize digraph.

        INPUT:
        data -- can be any of the following:
            1 NetworkX digraph
            2 dictionary of dictionaries
            3 dictionary of lists
            4 numpy matrix or ndarray
            5 pygraphviz agraph
            6 scipy sparse matrix
        pos -- a positioning dictionary: for example, the
        spring layout from NetworkX for the 5-cycle is
            {0: [-0.91679746, 0.88169588],
             1: [ 0.47294849, 1.125     ],
             2: [ 1.125     ,-0.12867615],
             3: [ 0.12743933,-1.125     ],
             4: [-1.125     ,-0.50118505]}
        name -- (in kwds) gives the graph a name
        loops -- boolean, whether to allow loops
        multiedges -- boolean, whether to allow multiple edges

        EXAMPLES:
        needed
        """
        import networkx
        if isinstance(data, DiGraph):
            self._nxg = data.networkx_graph()
        elif isinstance(data, networkx.DiGraph):
            self._nxg = networkx.XDiGraph(data, selfloops=loops, **kwds)
        elif isinstance(data, networkx.XDiGraph):
            self._nxg = data
        else:
            self._nxg = networkx.XDiGraph(data, selfloops=loops, **kwds)
        self.__pos = pos

    def _repr_(self):
        if not self._nxg.name is None and not self._nxg.name == "":
            name = self._nxg.name
            name = name + ": a s"
        else: name = "S"
        return name + "imple directed graph on %d vertices"%len(self._nxg.adj)

    def copy(self):
        """
        Creates a copy of the graph.
        """
        G = DiGraph(self._nxg, name=self._nxg.name)
        return G

    ### General Properties

    def is_directed(self):
        return True

    def loops(self, new=None):
        """
        Returns whether loops are permitted in the digraph.

        INPUT:
        new -- boolean, changes whether loops are permitted in the digraph.
        """
        if not new is None:
            if new:
                self._nxg.allow_selfloops()
            else:
                self._nxg.ban_selfloops()
        return self._nxg.selfloops

    def multiple_arcs(self, new=None):
        """
        Returns whether multiple arcs are permitted in the digraph.

        INPUT:
        new -- boolean, changes whether multiple arcs are permitted in the digraph.
        """
        if not new is None:
            if new:
                self._nxg.allow_multiedges()
            else:
                self._nxg.ban_multiedges()
        return self._nxg.multiedges

    ### Vertex Handlers

    def add_vertex(self, name=None):
        """
        Creates an isolated vertex.

        INPUT:
        n -- Name of the new vertex. If no name is specified, then the vertex
        will be represented by the least integer not already representing a
        vertex. Name must be an immutable object.
        """
        ### TODO- add doc note about representing other objects as vertices
        if name is None: # then find an integer to use as a key
            i = 0
            while self._nxg.succ.has_key(i):
                i=i+1
            self._nxg.add_node(i)
        else:
            self._nxg.add_node(name)

    def delete_vertex(self, vertex):
        """
        Deletes vertex, removing all incident arcs.
        """
        self._nxg.delete_node(vertex)

    def delete_vertices(self, vertices):
        """
        Remove vertices from the digraph taken from an iterable container of
        vertices.
        """
        self._nxg.delete_nodes_from(vertices)

    def neighbor_iterator(self, vertex):
        """
        Return an iterator over neighbors (connected either way) of vertex.
        """
        A = list(self._nxg.pred[vertex].iterkeys())
        B = list(self._nxg.succ[vertex].iterkeys())
        C = []
        for V in A:
            if not V in B:
                C += [V]
        for V in B:
            C += [V]
        return iter(C)

    def loop_vertices(self):
        """
        Returns a list of vertices with loops.
        """
        return self._nxg.nodes_with_selfloops()

    ### Arc Handlers

    def add_arc(self, u, v=None, label=None):
        """
        Adds an arc from u to v.

        INPUT:
        The following forms are all accepted by NetworkX:
        INPUT:
        The following forms are all accepted:

        G.add_arc( 1, 2 )
        G.add_arc( (1, 2) )
        G.add_arcs( [ (1, 2) ] )
        G.add_arc( 1, 2, 'label' )
        G.add_arc( (1, 2, 'label') )
        G.add_arcs( [ (1, 2, 'label') ] )

        WARNING:
        The following intuitive input results in nonintuitive output:
        sage: G = DiGraph()
        sage: G.add_arc((1,2),'label')
        sage: G.networkx_graph().adj           # random output order
        {'label': {}, (1, 2): {'label': None}}

        Use one of these instead:
        sage.: G = DiGraph()
        sage.: G.add_arc((1,2), label="label")
        sage.: G.networkx_graph().adj
        {1: {2: 'label'}, 2: {}}

        sage.: G = DiGraph()
        sage.: G.add_edge(1,2,'label')
        sage.: G.networkx_graph().adj
        {1: {2: 'label'}, 2: {}}
        """
        self._nxg.add_edge(u, v, label)

    def add_arcs(self, arcs):
        """
        Add arcs from an iterable container.
        """
        self._nxg.add_edges_from( arcs )

    def arcs(self, labels=True):
        """
        Return a list of arcs.

        INPUT:
        label -- if False, each edge is a tuple (u,v) of vertices.
        """
        L = self._nxg.edges()
        if labels:
            return L
        else:
            K = []
            for u,v,l in L:
                K.append((u,v))
            return K

    def arc_iterator(self, vertices=None, labels=True):
        """
        Returns an iterator over the arcs pointing out of the given
        set of vertices. If vertices is None, then returns an iterator over
        all arcs.

        INPUT:
        label -- if False, each edge is a tuple (u,v) of vertices.
        """
        L = self._nxg.edges_iter(vertices)
        if labels:
            return L
        else:
            K = []
            for u,v,l in L:
                K.append((u,v))
            return K

    def delete_arc(self, u, v=None, label=None):
        r"""
        Delete the arc from u to v, return silently if vertices or edge does
        not exist.

        INPUT:
        The following forms are all accepted:

        G.delete_arc( 1, 2 )
        G.delete_arc( (1, 2) )
        G.delete_arcs( [ (1, 2) ] )
        G.delete_arc( 1, 2, 'label' )
        G.delete_arc( (1, 2, 'label') )
        G.delete_arcs( [ (1, 2, 'label') ] )
        """
        self._nxg.delete_edge(u, v, label)

    def delete_arcs(self, arcs):
        """
        Delete arcs from an iterable container.
        """
        self._nxg.delete_edges_from(arcs)

    def delete_multiarc(self, u, v):
        """
        Deletes all arcs from u to v.
        """
        self._nxg.delete_multiedge(u, v)

    def incoming_arc_iterator(self, vertices=None, labels=True):
        """
        Return an iterator over all arriving arcs from vertices, or over all
        arcs if vertices is None.

        INPUT:
        label -- if False, each edge is a tuple (u,v) of vertices.
        """
        L = self._nxg.in_edges_iter(vertices)
        if labels:
            return L
        else:
            K = []
            for u,v,l in L:
                K.append((u,v))
            return K

    def incoming_arcs(self, vertices=None, labels=True):
        """
        Returns a list of arcs arriving at vertices.

        INPUT:
        label -- if False, each edge is a tuple (u,v) of vertices.
        """
        L = self._nxg.in_edges(vertices)
        if labels:
            return L
        else:
            K = []
            for u,v,l in L:
                K.append((u,v))
            return K

    def outgoing_arc_iterator(self, vertices=None):
        """
        Return an iterator over all departing arcs from vertices, or over all
        arcs if vertices is None.
        """
        return self._nxg.out_edges_iter(vertices)

    def outgoing_arcs(self, vertices=None, labels=True):
        """
        Returns a list of arcs departing from vertices.

        INPUT:
        label -- if False, each edge is a tuple (u,v) of vertices.
        """
        L = self._nxg.out_edges(vertices)
        if labels:
            return L
        else:
            K = []
            for u,v,l in L:
                K.append((u,v))
            return K

    def predecessor_iterator(self, vertex):
        """
        Returns an iterator over predecessor vertices of vertex.
        """
        return self._nxg.predecessors_iter(vertex)

    def predecessors(self, vertex):
        return list(self.predecessor_iterator(vertex))

    def successor_iterator(self, vertex):
        """
        Returns an iterator over successor vertices of vertex.
        """
        return self._nxg.successors_iter(vertex)

    def successors(self, vertex):
        return list(self.successor_iterator(vertex))

    def arc_label(self, u, v=None, label=None):
        """
        Returns the label of an arc.
        """
        return self._nxg.get_edge(u,v)

    def remove_multiple_arcs(self):
        self._nxg.remove_all_multiedges()

    def remove_loops(self, vertices=None):
        """
        Removes loops on vertices in vertices. If vertices is None, removes all loops.
        """
        if vertices is None:
            self._nxg.remove_all_selfloops()
        else:
            for v in vertices:
                self.delete_multiarc(v,v)

    def loop_arcs(self):
        """
        Returns a list of all loops in the graph.
        """
        return self._nxg.selfloop_edges()

    def number_of_loops(self):
        """
        Returns the number of arcs that are loops.
        """
        return self._nxg.number_of_selfloops()

    ### Degree functions

    def degree_iterator(self, vertices=None, with_labels=False):
        """
        with_labels=False:
            returns an iterator over degrees (in + out).
        with_labels=True:
            returns an iterator over tuples (vertex, degree (in + out) ).
        """
        return self._nxg.degree_iter(vertices, with_labels)

    def degree(self, vertices=None, with_labels=False):
        """
        Gives the degree (in + out) of a vertex or of vertices.

        INPUT:
        If vertices is a single vertex, returns the number of neighbors of
        vertex. If vertices is an iterable container of vertices, returns a
        list of degrees. If vertices is None, same as listing all vertices.

        OUTPUT:
        Single vertex- an integer. Multiple vertices- a list of integers. If
        with_labels is True, then returns a dictionary mapping each vertex to
        its degree.

        EXAMPLE:
        needed
        """
        return self._nxg.degree(vertices, with_labels)

    def in_degree(self, vertices=None, with_labels=False):
        return self._nxg.in_degree(vertices, with_labels)

    def in_degree_iterator(self, vertices=None, with_labels=False):
        """
        Same as degree_iterator, but for in degree.
        """
        return self._nxg.in_degree_iter(vertices, with_labels)

    def out_degree(self, vertices=None, with_labels=False):
        return self._nxg.out_degree(vertices, with_labels)

    def out_degree_iterator(self, vertices=None, with_labels=False):
        """
        Same as degree_iterator, but for out degree.
        """
        return self._nxg.out_degree_iter(vertices, with_labels)

    ### Representations

    def adjacency_matrix(self, sparse=True):
        """
        Returns the adjacency matrix of the digraph. Each vertex is
        represented by its position in the list returned by the vertices()
        function.
        """
        n = len(self._nxg.adj)
        verts = self.vertices()
        D = {}
        for e in self.arc_iterator():
            i,j,l = e
            i = verts.index(i)
            j = verts.index(j)
            D[(i,j)] = 1
        from sage.rings.integer_mod_ring import IntegerModRing
        from sage.matrix.constructor import matrix
        M = matrix(IntegerModRing(2), n, n, D, sparse=sparse)
        return M

    ### Contructors

    def reverse(self):
        """
        Returns a copy of digraph with arcs reversed in direction.
        """
        NXG = self._nxg.reverse()
        G = DiGraph(NXG)
        return G

    def subgraph(self, vertices, inplace=False, create_using=None):
        """
        Returns the subgraph induced by the given vertices.

        INPUT:
        inplace -- Using inplace is True will simply delete the extra vertices
        and edges from the current graph. This will modify the graph, and re-
        turn itself.
        vertices -- Vertices can be a single vertex or an iterable container
        of vertices, e.g. a list, set, graph, file or numeric array.
        create_using -- Can be an existing graph object or a call to a graph
        object, such as create_using=DiGraph().
        """
        if inplace:
            self._nxg = self._nxg.subgraph(vertices, inplace, create_using)
            return self
        else:
            NXG = self._nxg.subgraph(vertices, inplace, create_using)
            return DiGraph(NXG)

    def to_directed(self):
        return self.copy()

    def to_undirected(self):
        return Graph(self._nxg.to_undirected(), pos=self.__pos)

    ### Visualization

    def plot(self, pos=None, vertex_labels=True, node_size=200):
        GG = Graphics()
        if pos is None:
            if self.__pos is None:
                NGP = GraphicPrimitive_NetworkXGraph(self._nxg, pos=None, vertex_labels=vertex_labels, node_size=node_size)
            else:
                NGP = GraphicPrimitive_NetworkXGraph(self._nxg, pos=self.__pos, vertex_labels=vertex_labels, node_size=node_size)
        GG.append(NGP)
        GG.axes(False)
        return GG

    def show(self, pos=None, vertex_labels=True, node_size=200):
        self.plot(pos, vertex_labels, node_size).show()

class Network(GenericGraph):
    """
    Weighted multigraph: directed or undirected, allowing loops (non-hyper).
    """
    pass

### Hypergraphs and Complexes

class HyperGraph(SageObject):
    """
    Edges are simply subsets of the vertex set.
    """
    pass

class GenericComplex(SageObject):
    pass

class SimplicialComplex(GenericComplex):
    pass

class CubicalComplex(GenericComplex):
    pass

