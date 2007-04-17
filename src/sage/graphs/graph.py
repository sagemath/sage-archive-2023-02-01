r"""
Graph Theory

AUTHOR:
    -- Robert L. Miller (2006-10-22): initial version
    -- William Stein (2006-12-05): Editing
    -- Robert L. Miller (2007-01-13): refactoring, adjusting for
        NetworkX-0.33, fixed plotting bugs
                        (2007-01-23): basic tutorial, edge labels, loops,
        multiple edges & arcs
                        (2007-02-07): graph6 and sparse6 formats, matrix input
    -- Emily Kirkmann (2007-02-11): added graph_border option to plot and show
    -- Robert L. Miller (2007-02-12): vertex color-maps, graph boundaries,
        graph6 helper functions in SageX
                        SAGE Days 3 (2007-02-17--21): 3d plotting in Tachyon
                        (2007-02-25): display a partition
                        (2007-02-28): associate arbitrary objects to vertices,
        edge and arc label display (in 2d), edge coloring
                        (2007-03-21): Automorphism group, isomorphism check,
        canonical label

TUTORIAL:

    I. The Basics

        1. Graph Format

            A. The SAGE Graph Class: NetworkX plus

            SAGE graphs are actually NetworkX graphs, wrapped in a SAGE class.
            In fact, any graph can produce its underlying NetworkX graph. For example,

                sage: import networkx
                sage: G = graphs.PetersenGraph()
                sage: N = G.networkx_graph()
                sage: isinstance(N, networkx.graph.Graph)
                True

            The NetworkX graph is essentially a dictionary of dictionaries:

                sage: N.adj
                {0: {1: None, 4: None, 5: None}, 1: {0: None, 2: None, 6: None}, 2: {1: None, 3: None, 7: None}, 3: {8: None, 2: None, 4: None}, 4: {0: None, 9: None, 3: None}, 5: {0: None, 8: None, 7: None}, 6: {8: None, 1: None, 9: None}, 7: {9: None, 2: None, 5: None}, 8: {3: None, 5: None, 6: None}, 9: {4: None, 6: None, 7: None}}

            Each dictionary key is a vertex label, and each key in the following
            dictionary is a neighbor of that vertex. In undirected graphs, there
            is reduncancy: for example, the dictionary containing the entry
            1: {2: None} implies it must contain 2: {1: None}. The innermost entry
            of None is related to edge labelling (see section I.3.).

            B. Supported formats

            SAGE Graphs can be created from a wide range of inputs. A few examples are
            covered here.

                i.   a. NetworkX dictionary format:

                sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], 5: [7, 8], 6: [8,9], 7: [9]}
                sage: G = Graph(d); G
                Graph on 10 vertices
                sage: G.save('sage.png')

                    b. A NetworkX graph:

                sage: K = networkx.complete_bipartite_graph(12,7)
                sage: G = Graph(K)
                sage: G.degree()
                [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 12, 12, 12, 12, 12, 12, 12]

                ii. graph6 or sparse6 format:

                sage: s = ':I`AKGsaOs`cI]Gb~'
                sage: G = Graph(s); G
                Looped multi-graph on 10 vertices
                sage: G.save('sage.png')

                iii. adjacency matrix: In an adjacency matrix, each column and each row represent
                a vertex. If a 1 shows up in row i, column j, there is an edge (i,j).

                sage: M = Matrix([(0,1,0,0,1,1,0,0,0,0),(1,0,1,0,0,0,1,0,0,0),(0,1,0,1,0,0,0,1,0,0),(0,0,1,0,1,0,0,0,1,0),(1,0,0,1,0,0,0,0,0,1),(1,0,0,0,0,0,0,1,1,0),(0,1,0,0,0,0,0,0,1,1),(0,0,1,0,0,1,0,0,0,1),(0,0,0,1,0,1,1,0,0,0),(0,0,0,0,1,0,1,1,0,0)])
                sage: M
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
                sage: G = Graph(M); G
                Graph on 10 vertices
                sage: G.save('sage.png')

                iv. incidence matrix: In an incidence matrix, each row represents a vertex
                and each column reprensents an edge.

                sage: M = Matrix([(-1,0,0,0,1,0,0,0,0,0,-1,0,0,0,0),(1,-1,0,0,0,0,0,0,0,0,0,-1,0,0,0),(0,1,-1,0,0,0,0,0,0,0,0,0,-1,0,0),(0,0,1,-1,0,0,0,0,0,0,0,0,0,-1,0),(0,0,0,1,-1,0,0,0,0,0,0,0,0,0,-1),(0,0,0,0,0,-1,0,0,0,1,1,0,0,0,0),(0,0,0,0,0,0,0,1,-1,0,0,1,0,0,0),(0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0),(0,0,0,0,0,0,0,0,1,-1,0,0,0,1,0),(0,0,0,0,0,0,1,-1,0,0,0,0,0,0,1)])
                sage: M
                [-1  0  0  0  1  0  0  0  0  0 -1  0  0  0  0]
                [ 1 -1  0  0  0  0  0  0  0  0  0 -1  0  0  0]
                [ 0  1 -1  0  0  0  0  0  0  0  0  0 -1  0  0]
                [ 0  0  1 -1  0  0  0  0  0  0  0  0  0 -1  0]
                [ 0  0  0  1 -1  0  0  0  0  0  0  0  0  0 -1]
                [ 0  0  0  0  0 -1  0  0  0  1  1  0  0  0  0]
                [ 0  0  0  0  0  0  0  1 -1  0  0  1  0  0  0]
                [ 0  0  0  0  0  1 -1  0  0  0  0  0  1  0  0]
                [ 0  0  0  0  0  0  0  0  1 -1  0  0  0  1  0]
                [ 0  0  0  0  0  0  1 -1  0  0  0  0  0  0  1]
                sage: G = Graph(M); G
                Graph on 10 vertices
                sage: G.save('sage.png')

        2. Generators

        For some commonly used graphs to play with, type

            sage.: graphs.

        and hit <tab>. Most of these graphs come with their own custom plot, so you
        can see how people usually visualize these graphs.

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

            sage: S = G.subgraph([0,1,2,3])
            sage: S.plot().save('sage.png')    # or S.show()
            sage: S.density()
            1/2

            sage: L = graphs_query.get_list_of_graphs(nodes=7, diameter=5)
            sage.: graphs_list.show_graphs(L)

        3. Labels

        Each vertex can have any hashable object as a label. These are things like
        strings, numbers, and tuples. Each edge is given a default label of None, but
        if specified, edges can have any label at all. Edges between nodes u and v are
        represented typically as (u, v, l), where l is the label for the edge.

        Note that vertex labels themselves cannot be mutable items:

            sage: M = Matrix( [[0,0],[0,0]] )
            sage: G = Graph({ 0 : { M : None } })
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        However, if one wants to define a dictionary, with the same keys and arbitrary objects
        for entries, one can make that association:

            sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(),2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
            sage: d[2]
            Moebius-Kantor Graph: Graph on 16 vertices
            sage: T = graphs.TetrahedralGraph()
            sage: T.vertices()
            [0, 1, 2, 3]
            sage: T.associate(d)
            sage: T.obj(1)
            Flower Snark: Graph on 20 vertices

        4. Database

        There is a database available for searching for graphs that satisfy a certain set
        of parameters, including number of vertices and edges, density, maximum and minimum
        degree, diameter, radius, and connectivity. If you wish to search a database of
        graphs by parameter, type

            sage.: graphs_query.

        and hit tab.

            sage: L = graphs_query.get_list_of_graphs(nodes=7, diameter=5)
            sage.: graphs_list.show_graphs(L)

        6. Visualization

        To see a graph G you are working with, right now there are two main options:

            sage: G = graphs.RandomGNP(15,.3)

        You can view the graph in two dimensions via matplotlib:

            sage.: G.show()

        Or you can view it in three dimensions via Tachyon:

            sage.: G.show3d()

NOTE: Many functions are passed directly on to NetworkX, and in this
case the documentation is based on the NetworkX docs.
"""

#*****************************************************************************
#      Copyright (C) 2006 - 2007 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from random import random
from sage.structure.sage_object import SageObject
from sage.plot.plot import Graphics, GraphicPrimitive_NetworkXGraph
import sage.graphs.graph_fast as graph_fast

class GenericGraph(SageObject):
    """
    Base class for graphs and digraphs.
    """

    def __cmp__(self, other):
        """
        Comparison of self and other. Must be in the same class, have the same
        settings for loops and multiedges, output the same vertex list (in order)
        and the same adjacency matrix.

        EXAMPLES:
            sage: G = graphs.EmptyGraph()
            sage: H = Graph()
            sage: G == H
            True
            sage: G.to_directed() == H.to_directed()
            True
            sage: G = graphs.RandomGNP(8,.9999)
            sage: H = graphs.CompleteGraph(8)
            sage: G == H # (quasi-) random output (most often true)
            True
            sage: G = Graph( {0:[1,2,3,4,5,6,7]} )
            sage: H = Graph( {1:[0], 2:[0], 3:[0], 4:[0], 5:[0], 6:[0], 7:[0]} )
            sage: G == H
            True
            sage: G.loops(True)
            True
            sage: G == H
            False
            sage: G = graphs.RandomGNP(9,.3).to_directed()
            sage: H = graphs.RandomGNP(9,.3).to_directed()
            sage: G == H # random output (most often false)
            False
        """
        if type(self) != type(other):
            return 1
        elif self.loops() != other.loops():
            return 1
        else:
            try:
                if self.multiple_arcs() != other.multiple_arcs():
                    return 1
            except AttributeError:
                if self.multiple_edges() != other.multiple_edges():
                    return 1
        if self.vertices() != other.vertices():
            return 1
        comp = enum(self) - enum(other)
        if comp < 0:
            return -1
        elif comp == 0:
            return 0
        elif comp > 0:
            return 1


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
            Graph on 5 vertices
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
        """
        len(G) returns the number of vertices in G.
        """
        return len(self._nxg.adj)

    def __str__(self):
        """
        str(G) returns the name of the graph, unless it doesn't have one, in
        which case it returns the default refresentation.
        """
        name = self._nxg.name
        if name != "No Name" and (not name is None):
            return self._nxg.name
        else:
            return repr(self)

    def _latex_(self):
        """
        TODO: Apply Craig's algo.
        """
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

        EXAMPLE:
            sage: G = graphs.TetrahedralGraph()
            sage: N = G.networkx_graph()
            sage: type(N)
            <class 'networkx.xgraph.XGraph'>
        """
        return self._nxg.copy()

    def networkx_info(self, vertex=None):
        """
        Returns NetworkX information about the graph or the given node.
        """
        self._nxg.info(vertex)

    def __get_pos__(self):
        return self._pos

    def __set_pos__(self, pos):
        self._pos = pos

    ### General properties

    def name(self, new=None, set_to_none=False):
        """
        Returns the name of the (di)graph.

        INPUT:
        new -- if not None, then this becomes the new name of the (di)graph.
        set_to_none -- if True, removes any name

        EXAMPLE:
            sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], 5: [7, 8], 6: [8,9], 7: [9]}
            sage: G = Graph(d); G
            Graph on 10 vertices
            sage: G.name("Petersen Graph"); G
            'Petersen Graph'
            Petersen Graph: Graph on 10 vertices
            sage: G.name(set_to_none=True); G
            Graph on 10 vertices
        """
        if not new is None:
            if not isinstance(new, str):
                raise TypeError, "New name must be a string."
            self._nxg.name = new
        if set_to_none:
            self._nxg.name = None
        return self._nxg.name

    def loops(self, new=None):
        """
        Returns whether loops are permitted in the graph.

        INPUT:
        new -- boolean, changes whether loops are permitted in the graph.

        EXAMPLE:
            sage: G = Graph(); G
            Graph on 0 vertices
            sage: G.loops(True); G
            True
            Looped graph on 0 vertices

            sage: D = DiGraph(); D
            Digraph on 0 vertices
            sage: D.loops(True); D
            True
            Looped digraph on 0 vertices
        """
        if not new is None:
            if new:
                self._nxg.allow_selfloops()
            else:
                self._nxg.ban_selfloops()
        return self._nxg.selfloops

    def density(self):
        """
        Returns the density (number of edges divided by number of possible
        edges).

        EXAMPLE:
            sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], 5: [7, 8], 6: [8,9], 7: [9]}
            sage: G = Graph(d); G.density()
            1/3
        """
        from sage.rings.rational import Rational
        n = self.order()
        n = (n**2 - n)/2
        return Rational(self.size())/Rational(n)

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

    def has_vertex(self, vertex):
        """
        Indicates whether vertex is a vertex of the graph.

        EXAMPLE:
            sage: graphs.PetersenGraph().has_vertex(99)
            False
        """
        return self._nxg.has_node(vertex)

    def add_vertex(self, name=None):
        """
        Creates an isolated vertex.

        INPUT:
        n -- Name of the new vertex. If no name is specified, then the vertex
        will be represented by the least integer not already representing a
        vertex. Name must be an immutable object.

        EXAMPLES:
            sage: G = Graph(); G.add_vertex(); G
            Graph on 1 vertex

            sage: D = DiGraph(); D.add_vertex(); D
            Digraph on 1 vertex
        """
        if name is None: # then find an integer to use as a key
            i = 0
            while self.has_vertex(i):
                i=i+1
            self._nxg.add_node(i)
        else:
            self._nxg.add_node(name)

    def add_vertices(self, vertices):
        """
        Add vertices to the (di)graph from an iterable container of vertices.

        EXAMPLES:
            sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], 5: [7,8], 6: [8,9], 7: [9]}
            sage: G = Graph(d)
            sage: G.add_vertices([10,11,12])
            sage: G.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
            sage: G.add_vertices(graphs.CycleGraph(25).vertices())
            sage: G.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        """
        self._nxg.add_nodes_from(vertices)

    def delete_vertex(self, vertex):
        """
        Deletes vertex, removing all incident edges(arcs).

        EXAMPLES:
            sage: G = graphs.WheelGraph(9)
            sage: G.delete_vertex(0); G.save('sage.png')

            sage: D = DiGraph({0:[1,2,3,4,5],1:[2],2:[3],3:[4],4:[5],5:[1]})
            sage: D.delete_vertex(0); D
            Digraph on 5 vertices
        """
        self._nxg.delete_node(vertex)

    def delete_vertices(self, vertices):
        """
        Remove vertices from the (di)graph taken from an iterable container of
        vertices.

        EXAMPLE:
            sage: D = DiGraph({0:[1,2,3,4,5],1:[2],2:[3],3:[4],4:[5],5:[1]})
            sage: D.delete_vertices([1,2,3,4,5]); D
            Digraph on 1 vertex
        """
        self._nxg.delete_nodes_from(vertices)

    def get_boundary(self):
        return self.__boundary

    def set_boundary(self, boundary):
        self.__boundary = boundary

    def vertex_boundary(self, vertices1, vertices2=None):
        """
        Returns a list of all vertices in the external boundary of vertices1,
        intersected with vertices2. If vertices2 is None, then vertices2 is the
        complement of vertices1.

        EXAMPLE:
            sage: G = graphs.CubeGraph(4)
            sage: l = ['0111', '0000', '0001', '0011', '0010', '0101', '0100', '1111', '1101', '1011', '1001']
            sage: G.vertex_boundary(['0000', '1111'], l)
            ['0111', '1011', '1101', '0010', '0100', '0001']
        """
        return self._nxg.node_boundary(vertices1, vertices2)

    def associate(self, vertex_dict):
        """
        Associate arbitrary objects with each vertex, via an association dictionary.

        INPUT:
            vertex_dict -- the association dictionary

        EXAMPLES:
            sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), 2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
            sage: d[2]
            Moebius-Kantor Graph: Graph on 16 vertices
            sage: T = graphs.TetrahedralGraph()
            sage: T.vertices()
            [0, 1, 2, 3]
            sage: T.associate(d)
            sage: T.obj(1)
            Flower Snark: Graph on 20 vertices
        """
        for v in self.vertices():
            if not vertex_dict.has_key(v):
                vertex_dict[v] = None
        self._assoc = vertex_dict

    def obj(self, vertex):
        """
        Retrieve the object associated with a given vertex.

        INPUT:
            vertex -- the given vertex

        EXAMPLES:
            sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), 2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
            sage: d[2]
            Moebius-Kantor Graph: Graph on 16 vertices
            sage: T = graphs.TetrahedralGraph()
            sage: T.vertices()
            [0, 1, 2, 3]
            sage: T.associate(d)
            sage: T.obj(1)
            Flower Snark: Graph on 20 vertices
        """
        return self._assoc[vertex]

    def loop_vertices(self):
        """
        Returns a list of vertices with loops.

        EXAMPLE:
            sage: G = Graph({0 : [0], 1: [1,2,3], 2: [3]}, loops=True)
            sage: G.loop_vertices()
            [0, 1]
        """
        return self._nxg.nodes_with_selfloops()

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

    def relabel(self, perm, inplace=True, quick=False):
        r"""
        Uses a dictionary or permutation to relabel the (di)graph.
        If perm is a dictionary, each old vertex v is a key in the
        dictionary, and its new label is d[v]. If perm is a list,
        we think of it as a map i \mapsto perm[i] (only for graphs
        with V = {0,1,...,n-1} ). If perm is a per mutation, the
        permutation is simply applied to the graph, under the
        assumption that V = {0,1,...,n-1} is the vertex set, and
        the permutation acts on the set {1,2,...,n}, where we think
        of n = 0.

        INPUT:
            quick -- if True, simply return the enumeration of the new graph
        without constructing it. Requires that perm is of type list.

        EXAMPLES:
            sage: G = Graph({0:[1],1:[2],2:[]})
            sage: G.am()
            [0 1 0]
            [1 0 1]
            [0 1 0]

        Relabeling using a list:
            sage: G.relabel([0,2,1])
            sage: G.am()
            [0 0 1]
            [0 0 1]
            [1 1 0]

        Relabeling using a dictionary:
            sage: G.relabel({0:2,2:0})
            sage: G.am()
            [0 1 1]
            [1 0 0]
            [1 0 0]

        Relabeling using a SAGE permutation:
            sage: from sage.groups.perm_gps.permgroup import SymmetricGroup
            sage: S = SymmetricGroup(3)
            sage: gamma = S('(3,2)')
            sage: G.relabel(gamma)
            sage: G.am()
            [0 0 1]
            [0 0 1]
            [1 1 0]
        """
        if type(perm) == list:
            if quick:
                n = self.order()
                numbr = 0
                if isinstance(self, Graph):
                    for i,j,l in self.edge_iterator():
                        numbr += 1<<((n-(perm[i]+1))*n + n-(perm[j]+1))
                        numbr += 1<<((n-(perm[j]+1))*n + n-(perm[i]+1))
                elif isinstance(self, DiGraph):
                    for i,j,l in self.arc_iterator():
                        numbr += 1<<((n-(perm[i]+1))*n + n-(perm[j]+1))
                return numbr
            if isinstance(self, Graph):
                oldd = self._nxg.adj
                newd = {}
                for v in oldd.iterkeys():
                    oldtempd = oldd[v]
                    newtempd = {}
                    for w in oldtempd.iterkeys():
                        newtempd[perm[w]] = oldtempd[w]
                    newd[perm[v]] = newtempd
                if inplace:
                    self._nxg.adj = newd
                else:
                    G = self.copy()
                    G._nxg.adj = newd
                    return G
            else: # DiGraph
                oldsucc = self._nxg.succ
                oldpred = self._nxg.pred
                newsucc = {}
                newpred = {}
                for v in oldsucc.iterkeys():
                    oldtempsucc = oldsucc[v]
                    newtempsucc = {}
                    for w in oldtempsucc.iterkeys():
                        newtempsucc[perm[w]] = oldtempsucc[w]
                    newsucc[perm[v]] = newtempsucc
                for v in oldpred.iterkeys():
                    oldtemppred = oldpred[v]
                    newtemppred = {}
                    for w in oldtemppred.iterkeys():
                        newtemppred[perm[w]] = oldtemppred[w]
                    newpred[perm[v]] = newtemppred
                if inplace:
                    self._nxg.adj = newsucc
                    self._nxg.succ = self._nxg.adj
                    self._nxg.pred = newpred
                else:
                    D = self.copy()
                    D._nxg.adj = newsucc
                    D._nxg.succ = D._nxg.adj
                    D._nxg.pred = newpred
                    return D
        from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
        if type(perm) == PermutationGroupElement:
            n = self.order()
            dict = {}
            llist = perm.list()
            for i in range(1,n):
                dict[i] = llist[i-1]%n
            if n > 0:
                dict[0] = llist[n-1]%n
            perm = dict
        if type(perm) == type({}):
            keys = perm.keys()
            verts = self.vertices()
            for v in verts:
                if v not in keys:
                    perm[v] = v
            for v in perm.iterkeys():
                if v in verts:
                    try:
                        ddddd = {perm[v]:0}
                    except TypeError:
                        raise ValueError, "perm dictionary must be of the format {a:a1, b:b1, ...} where a,b,... are vertices and a1,b1,... are hashable"
            if isinstance(self, Graph):
                oldd = self._nxg.adj
                newd = {}
                for v in oldd.iterkeys():
                    oldtempd = oldd[v]
                    newtempd = {}
                    for w in oldtempd.iterkeys():
                        newtempd[perm[w]] = oldtempd[w]
                    newd[perm[v]] = newtempd
                if inplace:
                    self._nxg.adj = newd
                else:
                    G = self.copy()
                    G._nxg.adj = newd
                    return G
            else: # DiGraph
                oldsucc = self._nxg.succ
                oldpred = self._nxg.pred
                newsucc = {}
                newpred = {}
                for v in oldsucc.iterkeys():
                    oldtempsucc = oldsucc[v]
                    newtempsucc = {}
                    for w in oldtempsucc.iterkeys():
                        newtempsucc[perm[w]] = oldtempsucc[w]
                    newsucc[perm[v]] = newtempsucc
                for v in oldpred.iterkeys():
                    oldtemppred = oldpred[v]
                    newtemppred = {}
                    for w in oldtemppred.iterkeys():
                        newtemppred[perm[w]] = oldtemppred[w]
                    newpred[perm[v]] = newtemppred
                if inplace:
                    self._nxg.adj = newsucc
                    self._nxg.succ = self._nxg.adj
                    self._nxg.pred = newpred
                else:
                    D = self.copy()
                    D._nxg.adj = newsucc
                    D._nxg.succ = D._nxg.adj
                    D._nxg.pred = newpred
                    return D

    ### Constructors

    def am(self):
        """
        Shorter call for adjacency matrix makes life easier.
        """
        return self.adjacency_matrix()

    ### Visualization

    def plot(self, pos=None, layout=None, vertex_labels=True, edge_labels=False,
             node_size=200, graph_border=False, color_dict=None, partition=None,
             edge_colors=None, scaling_term=0.05):
        """
        Returns a graphics object representing the (di)graph.

        INPUT:
            pos -- an optional positioning dictionary
            layout -- what kind of layout to use, takes precedence over pos
                'circular' -- plots the graph with vertices evenly distributed on a circle
                'spring' -- uses the traditional spring layout, ignores the graphs current positions
            vertex_labels -- whether to print vertex labels
            edge_labels -- whether to print edge(arc) labels. By default, False, but if True, the result
                of str(l) is printed on the edge for each label l. Labels equal to None are not printed.
            node_size -- size of vertices displayed
            graph_border -- whether to include a box around the graph
            color_dict -- optional dictionary to specify vertex colors: each key is a color recognizable
                by matplotlib, and each corresponding entry is a list of vertices. If a vertex is not listed,
                it looks invisible on the resulting plot (it doesn't get drawn).
            edge_colors -- a dictionary specifying edge colors: each key is a color recognized by
                matplotlib, and each entry is a list of edges.
            partition -- a partition of the vertex set. if specified, plot will show each cell in a different
                color. color_dict takes precedence.
            scaling_term -- default is 0.05. if nodes are getting chopped off, increase; if graph
                is too small, decrease. should be positive, but values much bigger than
                1/8 won't be useful unless the nodes are huge

        EXAMPLES:
            sage: from math import sin, cos, pi
            sage: P = graphs.PetersenGraph()
            sage: d = {'#FF0000':[0,5], '#FF9900':[1,6], '#FFFF00':[2,7], '#00FF00':[3,8], '#0000FF':[4,9]}
            sage: pos_dict = {}
            sage: for i in range(5):
            ...    x = float(cos(pi/2 + ((2*pi)/5)*i))
            ...    y = float(sin(pi/2 + ((2*pi)/5)*i))
            ...    pos_dict[i] = [x,y]
            ...
            sage: for i in range(10)[5:]:
            ...    x = float(0.5*cos(pi/2 + ((2*pi)/5)*i))
            ...    y = float(0.5*sin(pi/2 + ((2*pi)/5)*i))
            ...    pos_dict[i] = [x,y]
            ...
            sage: pl = P.plot(pos=pos_dict, color_dict=d)
            sage: pl.save('sage.png')

            sage: C = graphs.CubeGraph(8)
            sage: P = C.plot(vertex_labels=False, node_size=0, graph_border=True)
            sage: P.save('sage.png')

            sage: G = graphs.HeawoodGraph()
            sage: for u,v,l in G.edges():
            ...    G.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: G.plot(edge_labels=True).save('sage.png')

            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []} )
            sage: for u,v,l in D.arcs():
            ...    D.set_arc_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: D.plot(edge_labels=True, layout='circular').save('sage.png')

            sage: from sage.plot.plot import rainbow
            sage: C = graphs.CubeGraph(5)
            sage: R = rainbow(5)
            sage: edge_colors = {}
            sage: for i in range(5):
            ...    edge_colors[R[i]] = []
            sage: for u,v,l in C.edges():
            ...    for i in range(5):
            ...        if u[i] != v[i]:
            ...            edge_colors[R[i]].append((u,v,l))
            sage: C.plot(vertex_labels=False, node_size=0, edge_colors=edge_colors).save('sage.png')
        """
        from sage.plot.plot import networkx_plot, rainbow
        import networkx
        if color_dict is None and not partition is None:
            l = len(partition)
            R = rainbow(l)
            color_dict = {}
            for i in range(l):
                color_dict[R[i]] = partition[i]
        if pos is None and layout is None:
            if not self._pos is None:
                pos = self._pos
        elif layout == 'circular':
            from math import sin, cos, pi
            n = self.order()
            verts = self.vertices()
            pos = {}
            for i in range(n):
                x = float(cos((pi/2) + ((2*pi)/n)*i))
                y = float(sin((pi/2) + ((2*pi)/n)*i))
                pos[verts[i]] = [x,y]
        elif layout == 'spring':
            pos = None
        if pos is None:
            pos = graph_fast.spring_layout_fast(self)
        else:
            for v in pos:
                for a in range(len(pos[v])):
                    pos[v][a] = float(pos[v][a])
        G = networkx_plot(self._nxg, pos=pos, vertex_labels=vertex_labels, node_size=node_size, color_dict=color_dict, edge_colors=edge_colors, graph_border=graph_border, scaling_term=scaling_term)
        if edge_labels:
            from sage.plot.plot import text
            K = Graphics()
            for u,v,l in self._nxg.edges():
                if not l is None:
                    K += text(str(l), [(pos[u][0] + pos[v][0])/2, (pos[u][1] + pos[v][1])/2])
            K.range(xmin=G.xmin(), xmax=G.xmax(), ymin=G.ymin(), ymax=G.ymax())
            G += K
            G.axes(False)
        return G

    def show(self, pos=None, layout=None, vertex_labels=True, edge_labels=False, node_size=200,
             graph_border=False, color_dict=None, edge_colors=None, partition=None,
             scaling_term=0.05, talk=False, **kwds):
        """
        Shows the (di)graph.

        INPUT:
            pos -- an optional positioning dictionary
            layout -- what kind of layout to use, takes precedence over pos
                'circular' -- plots the graph with vertices evenly distributed on a circle
                'spring' -- uses the traditional spring layout, ignores the graphs current positions
            vertex_labels -- whether to print vertex labels
            edge_labels -- whether to print edge(arc) labels. By default, False, but if True, the result
                of str(l) is printed on the edge for each label l. Labels equal to None are not printed.
            node_size -- size of vertices displayed
            graph_border -- whether to include a box around the graph
            color_dict -- optional dictionary to specify vertex colors: each key is a color recognizable
                by matplotlib, and each corresponding entry is a list of vertices. If a vertex is not listed,
                it looks invisible on the resulting plot (it doesn't get drawn).
            edge_colors -- a dictionary specifying edge colors: each key is a color recognized by
                matplotlib, and each entry is a list of edges.
            partition -- a partition of the vertex set. if specified, plot will show each cell in a different
                color. color_dict takes precedence.
            scaling_term -- default is 0.05. if nodes are getting chopped off, increase; if graph
                is too small, decrease. should be positive, but values much bigger than
                1/8 won't be useful unless the nodes are huge
            talk -- if true, prints large nodes with white backgrounds so that labels are legible on slies

        EXAMPLES:
            sage: from math import sin, cos, pi
            sage: P = graphs.PetersenGraph()
            sage: d = {'#FF0000':[0,5], '#FF9900':[1,6], '#FFFF00':[2,7], '#00FF00':[3,8], '#0000FF':[4,9]}
            sage: pos_dict = {}
            sage: for i in range(5):
            ...    x = float(cos(pi/2 + ((2*pi)/5)*i))
            ...    y = float(sin(pi/2 + ((2*pi)/5)*i))
            ...    pos_dict[i] = [x,y]
            ...
            sage: for i in range(10)[5:]:
            ...    x = float(0.5*cos(pi/2 + ((2*pi)/5)*i))
            ...    y = float(0.5*sin(pi/2 + ((2*pi)/5)*i))
            ...    pos_dict[i] = [x,y]
            ...
            sage: pl = P.plot(pos=pos_dict, color_dict=d)
            sage: pl.save('sage.png')

            sage: C = graphs.CubeGraph(8)
            sage: P = C.plot(vertex_labels=False, node_size=0, graph_border=True)
            sage: P.save('sage.png')

            sage: G = graphs.HeawoodGraph()
            sage: for u,v,l in G.edges():
            ...    G.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: G.plot(edge_labels=True).save('sage.png')

            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []} )
            sage: for u,v,l in D.arcs():
            ...    D.set_arc_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: D.plot(edge_labels=True, layout='circular').save('sage.png')

            sage: from sage.plot.plot import rainbow
            sage: C = graphs.CubeGraph(5)
            sage: R = rainbow(5)
            sage: edge_colors = {}
            sage: for i in range(5):
            ...    edge_colors[R[i]] = []
            sage: for u,v,l in C.edges():
            ...    for i in range(5):
            ...        if u[i] != v[i]:
            ...            edge_colors[R[i]].append((u,v,l))
            sage: C.plot(vertex_labels=False, node_size=0, edge_colors=edge_colors).save('sage.png')
        """
        if talk:
            node_size = 500
            if partition is None:
                color_dict = {'#FFFFFF':self.vertices()}
        self.plot(pos=pos, layout=layout, vertex_labels=vertex_labels, edge_labels=edge_labels, node_size=node_size, color_dict=color_dict, edge_colors=edge_colors, graph_border=graph_border, partition=partition, scaling_term=scaling_term).show(**kwds)

class Graph(GenericGraph):
    r"""
    Undirected graph.

    INPUT:
        data -- can be any of the following:
            1. A NetworkX graph
            2. A dictionary of dictionaries
            3. A dictionary of lists
            4. A numpy matrix or ndarray
            5. A graph6 or sparse6 string
            6. A SAGE adjacency matrix or incidence matrix
            7. A pygraphviz agraph
            8. A scipy sparse matrix

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
        format -- if None, Graph tries to guess- can be several values, including:
            'graph6' -- Brendan McKay's graph6 format, in a string (if the string has
                        multiple graphs, the first graph is taken)
            'sparse6' -- Brendan McKay's sparse6 format, in a string (if the string has
                        multiple graphs, the first graph is taken)
            'adjacency_matrix' -- a square SAGE matrix M, with M[i][j] equal to the number
                        of edges \{i,j\}
            'labeled_adjacency_matrix' -- a square SAGE matrix M, with M[i][j] equal to the
                        label of the single edge \{i,j\}
            'incidence_matrix' -- a SAGE matrix, with one column C for each edge, where
                        if C represents \{i, j\}, C[i] is -1 and C[j] is 1
        boundary -- a list of boundary vertices, if none, graph is considered as a 'graph
                    without boundary'

    EXAMPLES:
    We illustrate the first six input formats (the other two
    involve packages that are currently not standard in SAGE):

    1. A NetworkX graph:
        sage: import networkx
        sage: g = networkx.Graph({0:[1,2,3], 2:[5]})
        sage: Graph(g)
        Graph on 5 vertices

    2. A dictionary of dictionaries:
        sage: g = Graph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}); g
        Graph on 5 vertices

    The labels ('x', 'z', 'a', 'out') are labels for edges. For example, 'out' is
    the label for the edge on 2 and 5. Labels can be used as weights, if all the
    labels share some common parent.

    3. A dictionary of lists:
        sage: g = Graph({0:[1,2,3], 2:[5]}); g
        Graph on 5 vertices

    4. A numpy matrix or ndarray:
        sage: import numpy
        sage: A = numpy.array([[0,1,1],[1,0,1],[1,1,0]])
        sage: Graph(A)
        Graph on 3 vertices

    5. A graph6 or sparse6 string:
    SAGE automatically recognizes whether a string is in graph6 or sage6 format:

        sage: s = ':I`AKGsaOs`cI]Gb~'
        sage: Graph(s)
        Looped multi-graph on 10 vertices

    There are also list functions to take care of lists of graphs:

        sage: s = ':IgMoqoCUOqeb\n:I`AKGsaOs`cI]Gb~\n:I`EDOAEQ?PccSsge\N\n'
        sage: graphs_list.from_sparse6(s)
        [Looped multi-graph on 10 vertices, Looped multi-graph on 10 vertices, Looped multi-graph on 10 vertices]

    6. A SAGE matrix:
    Note: If format is not specified, then SAGE assumes a square matrix is an adjacency
    matrix, and a nonsquare matrix is an incidence matrix.

        A. an adjacency matrix:

        sage: M = graphs.PetersenGraph().am(); M
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
        sage: Graph(M)
        Graph on 10 vertices

        B. an incidence matrix:

        sage: M = Matrix(6, [-1,0,0,0,1, 1,-1,0,0,0, 0,1,-1,0,0, 0,0,1,-1,0, 0,0,0,1,-1, 0,0,0,0,0]); M
        [-1  0  0  0  1]
        [ 1 -1  0  0  0]
        [ 0  1 -1  0  0]
        [ 0  0  1 -1  0]
        [ 0  0  0  1 -1]
        [ 0  0  0  0  0]
        sage: Graph(M)
        Graph on 6 vertices
    """
    def __init__(self, data=None, pos=None, loops=False, format=None, boundary=None, **kwds):
        import networkx
        from sage.structure.element import is_Matrix
        if format is None:
            if isinstance(data, str):
                if data[:10] == ">>graph6<<":
                    data = data[10:]
                    format = 'graph6'
                elif data[:11] == ">>sparse6<<":
                    data = data[11:]
                    format = 'sparse6'
                elif data[0] == ':':
                    format = 'sparse6'
                else:
                    format = 'graph6'
            elif is_Matrix(data):
                if data.is_square():
                    format = 'adjacency_matrix'
                else:
                    format = 'incidence_matrix'
            elif isinstance(data, Graph):
                self._nxg = data.networkx_graph()
            elif isinstance(data, networkx.Graph):
                self._nxg = networkx.XGraph(data, selfloops=loops, **kwds)
            elif isinstance(data, networkx.XGraph):
                self._nxg = data
            else:
                self._nxg = networkx.XGraph(data, selfloops=loops, **kwds)
        if format == 'graph6':
            if not isinstance(data, str):
                raise ValueError, 'If input format is graph6, then data must be a string'
            from sage.rings.integer import Integer
            n = data.find('\n')
            if n == -1:
                n = len(data)
            s = data[:n]
            n, s = graph_fast.N_inverse(s)
            m = graph_fast.R_inverse(s, n)
            d = {}
            k = 0
            for i in range(n):
                d[i] = {}
                for j in range(i):
                    if m[k] == '1':
                        d[i][j] = None
                    k += 1
            self._nxg = networkx.XGraph(d)
        elif format == 'sparse6':
            from sage.rings.arith import ceil, floor
            from sage.misc.functional import log
            n = data.find('\n')
            if n == -1:
                n = len(data)
            s = data[:n]
            n, s = graph_fast.N_inverse(s[1:])
            k = ceil(log(n,2))
            l = [graph_fast.binary(ord(i)-63) for i in s]
            for i in range(len(l)):
                l[i] = '0'* (6-len(l[i])) + l[i]
            bits = ''.join(l)
            b = []
            x = []
            for i in range(floor(len(bits)/(k+1))):
                b.append(int(bits[(k+1)*i:(k+1)*i+1],2))
                x.append(int(bits[(k+1)*i+1:(k+1)*i+k+1],2))
            v = 0
            edges = []
            for i in range(len(b)):
                if b[i] == 1:
                    v += 1
                if x[i] > v:
                    v = x[i]
                else:
                    if v < n:
                        edges.append((x[i],v))
            d = {}
            for i,j in edges:
                if d.has_key(i):
                    if d[i].has_key(j):
                        if d[i][j] is None:
                            d[i][j] = [None,None]
                        else:
                            d[i][j].append(None)
                    d[i][j] = None
                else:
                    d[i] = {j : None}
            for i in [j for j in range(n) if not d.has_key(j)]:
                d[i] = {}
            self._nxg = networkx.XGraph(d, selfloops = True, multiedges = True)
        elif format == 'adjacency_matrix':
            d = {}
            for i in range(data.nrows()):
                d[i] = {}
            self._nxg = networkx.XGraph(d, selfloops = loops, **kwds)
            e = []
            for i,j in data.nonzero_positions():
                if i < j and kwds.get('multiedges',False):
                    e += [(i,j)]*int(data[i][j])
                elif i < j:
                    e.append((i,j))
                elif i == j and loops and kwds.get(multiedges,False):
                    e += [(i,j)]*int(data[i][j])
                elif i == j and loops:
                    e.append((i,j))
            self._nxg.add_edges_from(e)
        elif format == 'weighted_adjacency_matrix':
            d = {}
            for i in range(data.nrows()):
                d[i] = {}
            self._nxg = networkx.XGraph(d, selfloops = loops, **kwds)
            e = []
            for i,j in data.nonzero_positions():
                if i < j:
                    e.append((i,j,data[i][j]))
                elif i == j and loops:
                    e.append((i,j,data[i][j]))
            self._nxg.add_edges_from(e)
        elif format == 'incidence_matrix':
            b = True
            for c in data.columns():
                d = c.dict()
                if not len(d) == 2:
                    b = False
                else:
                    k = d.keys()
                    if not (d[k[0]] == -1 * d[k[1]] and abs(d[k[0]]) == 1):
                        b = False
            if not b:
                raise AttributeError, "Incidence Matrix must have one 1 and one -1 per column."
            else:
                d = {}
                for i in range(data.nrows()):
                    d[i] = {}
                self._nxg = networkx.XGraph(d, selfloops = loops, **kwds)
                e = []
                for c in data.columns():
                    k = c.dict().keys()
                    e.append((k[0],k[1]))
                self._nxg.add_edges_from(e)
        if kwds.has_key('name'):
            self._nxg.name = kwds['name']
        self._pos = pos
        self._boundary = boundary

    def _repr_(self):
        name = ""
        if self.loops():
            name += "looped "
        if self.multiple_edges():
            name += "multi-"
        name += "graph on %d vert"%self.order()
        if self.order() == 1:
            name += "ex"
        else:
            name += "ices"
        name = name.capitalize()
        if not self._nxg.name is None and not self._nxg.name == "":
            name = self._nxg.name + ": " + name
        return name

    def copy(self):
        """
        Creates a copy of the graph.
        """
        G = Graph(self._nxg, name=self._nxg.name, pos=self._pos, loops=self.loops(), boundary=self._boundary)
        return G

    def to_directed(self):
        """
        Returns a directed version of the graph. A single edge becomes two
        arcs, one in each direction.

        EXAMPLE:
            sage: graphs.PetersenGraph().to_directed()
            Digraph on 10 vertices
        """
        return DiGraph(self._nxg.to_directed(), pos=self._pos)

    def to_undirected(self):
        """
        Since the graph is already undirected, simply returns a copy of itself.

        EXAMPLE:
            sage: graphs.PetersenGraph().to_undirected()
            Petersen graph: Graph on 10 vertices
        """
        return self.copy()

    ### General properties

    def is_directed(self):
        """
        Since graph is undirected, returns False.
        """
        return False

    def multiple_edges(self, new=None):
        """
        Returns whether multiple edges are permitted in the graph.

        INPUT:
        new -- boolean, changes whether multiple edges are permitted in the graph.

        EXAMPLE:
            sage: G = Graph(multiedges=True); G
            Multi-graph on 0 vertices
            sage: G.multiple_edges(False); G
            False
            Graph on 0 vertices
        """
        if not new is None:
            if new:
                self._nxg.allow_multiedges()
            else:
                self._nxg.ban_multiedges()
        return self._nxg.multiedges

    ### Vertex handlers

    def neighbor_iterator(self, vertex):
        """
        Return an iterator over neighbors of vertex.

        EXAMPLE:
            sage: G = graphs.CubeGraph(3)
            sage: for i in G.neighbor_iterator('010'):
            ...    print i
            011
            000
            110
        """
        return self._nxg.neighbors_iter(vertex)

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
        sage: G = Graph()
        sage: G.add_edge((1,2), 'label')
        sage: G.networkx_graph().adj           # random output order
        {'label': {(1, 2): None}, (1, 2): {'label': None}}

        Use one of these instead:
        sage: G = Graph()
        sage: G.add_edge((1,2), label="label")
        sage: G.networkx_graph().adj           # random output order
        {1: {2: 'label'}, 2: {1: 'label'}}

        sage: G = Graph()
        sage: G.add_edge(1,2,'label')
        sage: G.networkx_graph().adj           # random output order
        {1: {2: 'label'}, 2: {1: 'label'}}
        """
        self._nxg.add_edge(u, v, label)

    def add_edges(self, edges):
        """
        Add edges from an iterable container.

        EXAMPLE:
            sage: G = graphs.DodecahedralGraph()
            sage: H = Graph()
            sage: H.add_edges( G.edge_iterator() ); H
            Graph on 20 vertices
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

        EXAMPLES:
            sage: G = graphs.CompleteGraph(19)
            sage: G.size()
            171
            sage: G.delete_edge( 1, 2 )
            sage: G.delete_edge( (3, 4) )
            sage: G.delete_edges( [ (5, 6), (7, 8) ] )
            sage: G.delete_edge( 9, 10, 'label' )
            sage: G.delete_edge( (11, 12, 'label') )
            sage: G.delete_edges( [ (13, 14, 'label') ] )
            sage: G.size()
            164
            sage: G.has_edge( (11, 12) )
            False

            Note that even though the edge (11, 12) has no label, it still gets
            deleted: NetworkX does not pay attention to labels here.
        """
        self._nxg.delete_edge(u, v, label)

    def delete_edges(self, edges):
        """
        Delete edges from an iterable container.

        EXAMPLE:
            sage: K12 = graphs.CompleteGraph(12)
            sage: K4 = graphs.CompleteGraph(4)
            sage: K12.size()
            66
            sage: K12.delete_edges(K4.edge_iterator())
            sage: K12.size()
            60
        """
        self._nxg.delete_edges_from(edges)

    def delete_multiedge(self, u, v):
        """
        Deletes all edges on u and v.

        EXAMPLE:
            sage: G = Graph(multiedges=True)
            sage: G.add_edges([(0,1), (0,1), (0,1), (1,2), (2,3)])
            sage: G.edges()
            [(0, 1, None), (0, 1, None), (0, 1, None), (1, 2, None), (2, 3, None)]
            sage: G.delete_multiedge( 0, 1 )
            sage: G.edges()
            [(1, 2, None), (2, 3, None)]
        """
        self._nxg.delete_multiedge(u, v)

    def edges(self, labels=True):
        """
        Return a list of edges. Each edge is a triple (u,v,l) where u and v are
        vertices and l is a label.

        INPUT:
        labels -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLES:
            sage: graphs.DodecahedralGraph().edges()
            [(0, 1, None), (0, 10, None), (0, 19, None), (1, 8, None), (1, 2, None), (2, 3, None), (2, 6, None), (3, 19, None), (3, 4, None), (4, 17, None), (4, 5, None), (5, 6, None), (5, 15, None), (6, 7, None), (7, 8, None), (7, 14, None), (8, 9, None), (9, 10, None), (9, 13, None), (10, 11, None), (11, 12, None), (11, 18, None), (12, 16, None), (12, 13, None), (13, 14, None), (14, 15, None), (15, 16, None), (16, 17, None), (17, 18, None), (18, 19, None)]

            sage: graphs.DodecahedralGraph().edges(labels=False)
            [(0, 1), (0, 10), (0, 19), (1, 8), (1, 2), (2, 3), (2, 6), (3, 19), (3, 4), (4, 17), (4, 5), (5, 6), (5, 15), (6, 7), (7, 8), (7, 14), (8, 9), (9, 10), (9, 13), (10, 11), (11, 12), (11, 18), (12, 16), (12, 13), (13, 14), (14, 15), (15, 16), (16, 17), (17, 18), (18, 19)]
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
        labels -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLE:
            sage: K = graphs.CompleteBipartiteGraph(9,3)
            sage: len(K.edge_boundary( [0,1,2,3,4,5,6,7,8], [9,10,11] ))
            27
            sage: K.size()
            27
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

        INPUT:
        labels -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLE:
            sage: for i in graphs.PetersenGraph().edge_iterator([0]):
            ...    print i
            (0, 1, None)
            (0, 4, None)
            (0, 5, None)
        """
        return self._nxg.edges_iter(vertices)

    def edges_incident(self, vertices=None, labels=True):
        """
        Returns a list of edges incident with any vertex given. If vertex is
        None, returns a list of all edges in graph.

        INPUT:
        label -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLE:
            sage: graphs.PetersenGraph().edges_incident([0,9], labels=False)
            [(0, 1), (0, 4), (0, 5), (9, 4), (9, 6), (9, 7)]
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

        EXAMPLE:
            sage: graphs.EmptyGraph().has_edge(9,2)
            False
        """
        return self._nxg.has_edge(u, v)

    def set_edge_label(self, u, v, l):
        """
        Set the edge label of a given edge.

        INPUT:
            u, v -- the vertices of the edge
            l -- the new label

        EXAMPLE:
            sage: G = graphs.HeawoodGraph()
            sage: for u,v,l in G.edges():
            ...    G.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: G.plot(edge_labels=True).save('sage.png')
        """
        if self.has_edge(u, v):
            self._nxg.adj[u][v] = l
            self._nxg.adj[v][u] = l

    def edge_label(self, u, v=None):
        """
        Returns the label of an edge.

        EXAMPLE:
            sage: G = Graph({0 : {1 : 'edgelabel'}})
            sage: G.edges(labels=False)
            [(0, 1)]
            sage: G.edge_label( 0, 1 )
            'edgelabel'
        """
        return self._nxg.get_edge(u,v)

    def edge_labels(self):
        """
        Returns a list of edge labels.

        EXAMPLE:
            sage: G = Graph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}})
            sage: G.edge_labels()
            ['x', 'z', 'a', 'out']
        """
        labels = []
        for u,v,l in self.edges():
            labels.append(l)
        return labels

    def remove_multiple_edges(self):
        """
        Removes all multiple edges, retaining one edge for each.

        EXAMPLE:
            sage: G = Graph(multiedges=True)
            sage: G.add_edges( [ (0,1), (0,1), (0,1), (0,1), (1,2) ] )
            sage: G.edges(labels=False)
            [(0, 1), (0, 1), (0, 1), (0, 1), (1, 2)]

            sage: G.remove_multiple_edges()
            sage: G.edges(labels=False)
            [(0, 1), (1, 2)]
        """
        self._nxg.remove_all_multiedges()

    def remove_loops(self, vertices=None):
        """
        Removes loops on vertices in vertices. If vertices is None, removes all loops.

        EXAMPLE
            sage: G = Graph(loops=True)
            sage: G.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: G.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: G.remove_loops()
            sage: G.edges(labels=False)
            [(2, 3)]
            sage: G.loops()
            True
        """
        if vertices is None:
            self._nxg.remove_all_selfloops()
        else:
            for v in vertices:
                self.delete_multiedge(v,v)

    def loop_edges(self):
        """
        Returns a list of all loops in the graph.

        EXAMPLE:
            sage: G = Graph(loops=True)
            sage: G.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: G.loop_edges()
            [(0, 0, None), (1, 1, None), (2, 2, None), (3, 3, None)]
        """
        return self._nxg.selfloop_edges()

    def number_of_loops(self):
        """
        Returns the number of edges that are loops.

        EXAMPLE:
            sage: G = Graph(loops=True)
            sage: G.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: G.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: G.number_of_loops()
            4
        """
        return self._nxg.number_of_selfloops()

    ### Degree functions

    def degree(self, vertices=None, labels=False):
        """
        Gives the degree of a vertex or of vertices.

        INPUT:
        vertices -- If vertices is a single vertex, returns the number of
        neighbors of vertex. If vertices is an iterable container of vertices,
        returns a list of degrees. If vertices is None, same as listing all vertices.
        labels -- see OUTPUT

        OUTPUT:
        Single vertex- an integer. Multiple vertices- a list of integers. If
        labels is True, then returns a dictionary mapping each vertex to
        its degree.

        EXAMPLES:
            sage: P = graphs.PetersenGraph()
            sage: P.degree(5)
            3

            sage: K = graphs.CompleteGraph(9)
            sage: K.degree()
            [8, 8, 8, 8, 8, 8, 8, 8, 8]
        """
        return self._nxg.degree(vertices, with_labels=labels)

    def degree_histogram(self):
        """
        Returns a list, whose ith entry is the frequency of degree i.

        EXAMPLE:
            sage: G = graphs.Grid2dGraph(9,12)
            sage: G.degree_histogram()
            [0, 0, 4, 34, 70]
        """
        import networkx
        return networkx.degree_histogram(self._nxg)

    def degree_iterator(self, vertices=None, labels=False):
        """
        INPUT:
        labels=False:
            returns an iterator over degrees.
        labels=True:
            returns an iterator over tuples (vertex, degree).
        vertices -- if specified, restrict to this subset.

        EXAMPLES:
            sage: G = graphs.Grid2dGraph(3,4)
            sage: for i in G.degree_iterator():
            ...    print i
            3
            4
            2
            3
            3
            2
            3
            2
            3
            3
            2
            4
            sage: for i in G.degree_iterator(labels=True):
            ...    print i
            ((0, 1), 3)
            ((1, 2), 4)
            ((0, 0), 2)
            ((2, 1), 3)
            ((0, 2), 3)
            ((2, 0), 2)
            ((1, 3), 3)
            ((2, 3), 2)
            ((2, 2), 3)
            ((1, 0), 3)
            ((0, 3), 2)
            ((1, 1), 4)
        """
        return self._nxg.degree_iter(vertices, with_labels=labels)

    ### Representations

    def adjacency_matrix(self, sparse=True):
        """
        Returns the adjacency matrix of the graph. Each vertex is
        represented by its position in the list returned by the vertices()
        function.

        EXAMPLE:
            sage: G = graphs.CubeGraph(4)
            sage: G.adjacency_matrix()
            [0 1 0 0 0 1 0 1 0 1 0 0 0 0 0 0]
            [1 0 0 0 1 0 1 0 1 0 0 0 0 0 0 0]
            [0 0 0 1 0 1 0 1 0 0 0 0 0 0 0 1]
            [0 0 1 0 1 0 1 0 0 0 0 0 0 0 1 0]
            [0 1 0 1 0 1 0 0 0 0 0 0 0 1 0 0]
            [1 0 1 0 1 0 0 0 0 0 0 0 1 0 0 0]
            [0 1 0 1 0 0 0 1 0 0 0 1 0 0 0 0]
            [1 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0 1 0 1 0 1 0 0]
            [1 0 0 0 0 0 0 0 1 0 1 0 1 0 0 0]
            [0 0 0 0 0 0 0 1 0 1 0 1 0 0 0 1]
            [0 0 0 0 0 0 1 0 1 0 1 0 0 0 1 0]
            [0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 1]
            [0 0 0 0 1 0 0 0 1 0 0 0 1 0 1 0]
            [0 0 0 1 0 0 0 0 0 0 0 1 0 1 0 1]
            [0 0 1 0 0 0 0 0 0 0 1 0 1 0 1 0]

        """
        n = len(self._nxg.adj)
        verts = self.vertices()
        D = {}
        for e in self.edge_iterator():
            i,j,l = e
            i = verts.index(i)
            j = verts.index(j)
            if D.has_key((i,j)) and self.multiple_edges():
                D[(i,j)] += 1
                D[(j,i)] += 1
            else:
                D[(i,j)] = 1
                D[(j,i)] = 1
        from sage.rings.integer_mod_ring import IntegerModRing
        from sage.rings.integer_ring import IntegerRing
        from sage.matrix.constructor import matrix
        if self.multiple_edges:
            R = IntegerRing()
        else:
            R = IntegerModRing(2)
        M = matrix(R, n, n, D, sparse=sparse)
        return M

    def incidence_matrix(self, sparse=True):
        """
        Returns an incidence matrix of the graph. Each row is a vertex, and
        each column is an edge.

        EXAMPLE:
            sage: G = graphs.CubeGraph(3)
            sage: G.incidence_matrix()
            [-1 -1 -1  0  0  0  0  0  0  0  0  0]
            [ 1  0  0 -1 -1  0  0  0  0  0  0  0]
            [ 0  0  0  1  0 -1 -1  0  0  0  0  0]
            [ 0  1  0  0  0  0  1 -1  0  0  0  0]
            [ 0  0  0  0  1  0  0  0 -1 -1  0  0]
            [ 0  0  1  0  0  0  0  0  0  1 -1  0]
            [ 0  0  0  0  0  0  0  1  0  0  1 -1]
            [ 0  0  0  0  0  1  0  0  1  0  0  1]
        """
        from sage.matrix.constructor import matrix
        from copy import copy
        n = len(self._nxg.adj)
        verts = self.vertices()
        d = [0]*n
        cols = []
        for i, j, l in self.edge_iterator():
            col = copy(d)
            i = verts.index(i)
            j = verts.index(j)
            col[i] = -1
            col[j] = 1
            cols.append(col)
        return matrix(cols, sparse=sparse).transpose()

    def __bit_vector(self):
        vertices = self.vertices()
        n = len(vertices)
        nc = int(n*(n - 1))/int(2)
        bit_vector = set()
        for e,f,g in self.edge_iterator():
            c = vertices.index(e)
            d = vertices.index(f)
            a,b = sorted([c,d])
            p = int(b*(b - 1))/int(2) + a
            bit_vector.add(p)
        bit_vector = sorted(bit_vector)
        s = []
        j = 0
        for i in bit_vector:
            s.append( '0'*(i - j) + '1' )
            j = i + 1
        s = "".join(s)
        s += '0'*(nc-len(s))
        return s

    def graph6_string(self):
        """
        Returns the graph6 representation of the graph as an ASCII string. Only valid
        for simple (no loops, multiple edges) graphs on 0 to 262143 vertices.

        EXAMPLE:
            sage: G = graphs.KrackhardtKiteGraph()
            sage: G.graph6_string()
            'IvUqwK@?G'
        """
        n = self.order()
        if n > 262143:
            raise ValueError, 'graph6 format supports graphs on 0 to 262143 vertices only.'
        elif self.loops() or self.multiple_edges():
            raise ValueError, 'graph6 format supports only simple graphs (no loops, no multiple edges)'
        else:
            return graph_fast.N(n) + graph_fast.R(self.__bit_vector())

    def sparse6_string(self):
        """
        Returns the sparse6 representation of the graph as an ASCII string. Only valid
        for undirected graphs on 0 to 262143 vertices, but loops and multiple edges are
        permitted.

        EXAMPLE:
            sage: G = graphs.BullGraph()
            sage: G.sparse6_string()
            ':Da@en'
        """
        n = self.order()
        if n > 262143:
            raise ValueError, 'sparse6 format supports graphs on 0 to 262143 vertices only.'
        else:
            vertices = self.vertices()
            n = len(vertices)
            edges = self.edges(labels=False)
            for e in edges: # replace edge labels with natural numbers (by index in vertices)
                e = (vertices.index(e[0]),vertices.index(e[1]))
            # order edges
            def cmp(x, y):
                if x[1] < y[1]:
                    return -1
                elif x[1] > y[1]:
                    return 1
                elif x[1] == y[1]:
                    if x[0] < y[0]:
                        return -1
                    if x[0] > y[0]:
                        return 1
                    else:
                        return 0
            edges.sort(cmp)

            # encode bit vector
            from sage.rings.arith import ceil
            from sage.misc.functional import log
            k = ceil(log(n,2))
            v = 0
            i = 0
            m = 0
            s = ''
            while m < len(edges):
                if edges[m][1] > v + 1:
                    sp = graph_fast.binary(edges[m][1])
                    sp = '0'*(k-len(sp)) + sp
                    s += '1' + sp
                    v = edges[m][1]
                elif edges[m][1] == v + 1:
                    sp = graph_fast.binary(edges[m][0])
                    sp = '0'*(k-len(sp)) + sp
                    s += '1' + sp
                    v += 1
                    m += 1
                else:
                    sp = graph_fast.binary(edges[m][0])
                    sp = '0'*(k-len(sp)) + sp
                    s += '0' + sp
                    m += 1

            # encode s as a 6-string, as in R(x), but padding with 1's
            # pad on the right to make a multiple of 6
            s = s + ( '1' * ((6 - len(s))%6) )

            # split into groups of 6, and convert numbers to decimal, adding 63
            six_bits = ''
            for i in range(len(s)/6):
                six_bits += chr( int( s[6*i:6*(i+1)], 2) + 63 )
            return ':' + graph_fast.N(n) + six_bits

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
        object, such as create_using=DiGraph(). Must be a NetworkX object.

        EXAMPLES:
            sage: G = graphs.CompleteGraph(9)
            sage: H = G.subgraph([0,1,2]); H
            Graph on 3 vertices
            sage: G
            Complete graph: Graph on 9 vertices
            sage: K = G.subgraph([0,1,2], inplace=True); K
            Subgraph of (Complete graph): Graph on 3 vertices
            sage: G
            Subgraph of (Complete graph): Graph on 3 vertices
            sage: G is K
            True
        """
        if inplace:
            self._nxg = self._nxg.subgraph(vertices, inplace, create_using)
            return self
        else:
            NXG = self._nxg.subgraph(vertices, inplace, create_using)
            return Graph(NXG)

    ### Visualization

    def plot3d(self, bgcolor=(1,1,1), vertex_color=(1,0,0), edge_color=(0,0,0), pos3d=None):
        """
        Plots the graph using Tachyon, and returns a Tachyon object containing
        a representation of the graph.

        INPUT:
            bgcolor
            vertex_color
            edge_color
            pos3d -- a position dictionary for the vertices

        EXAMPLES:
            sage: D = graphs.DodecahedralGraph()
            sage: P3D = D.plot3d()
            sage: P3D.save('sage.png') # long time

            sage: G = graphs.PetersenGraph()
            sage: G.plot3d(vertex_color=(0,0,1)).save('sage.png') # long time

            sage: C = graphs.CubeGraph(4)
            sage: C.plot3d(edge_color=(0,1,0), vertex_color=(1,1,1), bgcolor=(0,0,0)).save('sage.png') # long time
        """
        TT, pos3d = tachyon_vertex_plot(self, bgcolor=bgcolor, vertex_color=vertex_color, pos3d=pos3d)
        TT.texture('edge', ambient=0.1, diffuse=0.9, specular=0.03, opacity=1.0, color=edge_color)
        for u,v,l in self.edges():
            TT.fcylinder( (pos3d[u][0],pos3d[u][1],pos3d[u][2]),\
                          (pos3d[v][0],pos3d[v][1],pos3d[v][2]), .02,'edge')
        return TT

    def show3d(self, bgcolor=(1,1,1), vertex_color=(1,0,0), edge_color=(0,0,0), pos3d=None, **kwds):
        """
        Plots the graph using Tachyon, and shows the resulting plot.

        INPUT:
            bgcolor -- background color
            vertex_color -- vertex color
            edge_color -- edge color
            (pos3d -- currently ignored, pending GSL random point distribution in sphere...)

        EXAMPLES:
            sage: D = graphs.DodecahedralGraph()
            sage: P3D = D.plot3d()
            sage: P3D.save('sage.png') # long time

            sage: G = graphs.PetersenGraph()
            sage: G.plot3d(vertex_color=(0,0,1)).save('sage.png') # long time

            sage: C = graphs.CubeGraph(4)
            sage: C.plot3d(edge_color=(0,1,0), vertex_color=(1,1,1), bgcolor=(0,0,0)).save('sage.png') # long time
        """
        self.plot3d(bgcolor=bgcolor, vertex_color=vertex_color, edge_color=edge_color).show(**kwds)

    ### Connected components

    def is_connected(self):
        """
        Indicates whether the graph is connected. Note that in a graph, path
        connected is equivalent to connected.

        EXAMPLE:
            sage: G = Graph( { 0 : [1, 2], 1 : [2], 3 : [4, 5], 4 : [5] } )
            sage: G.is_connected()
            False
            sage: G.add_edge(0,3)
            sage: G.is_connected()
            True
        """
        import networkx
        return networkx.component.is_connected(self._nxg)

    def connected_components(self):
        """
        Returns a list of lists of vertices, each list representing a
        connected component. The list is ordered from largest to smallest
        component.

        EXAMPLE:
            sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: G.connected_components()
            [[0, 1, 2, 3], [4, 5, 6]]
        """
        import networkx
        return networkx.component.connected_components(self._nxg)

    def connected_components_number(self):
        """
        Returns the number of connected components.

        EXAMPLE:
            sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: G.connected_components_number()
            2
        """
        import networkx
        return networkx.component.number_connected_components(self._nxg)

    def connected_components_subgraphs(self):
        """
        Returns a list of connected components as graph objects.

        EXAMPLE:
            sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: L = G.connected_components_subgraphs()
            sage.: graphs_list.show_graphs(L)
        """
        cc = self.connected_components()
        list = []
        for c in cc:
            list.append(self.subgraph(c, inplace=False))
        return list

    def connected_component_containing_vertex(self, vertex):
        """
        Returns a list of the vertices connected to vertex.

        EXAMPLE:
            sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: G.connected_component_containing_vertex(0)
            [0, 1, 2, 3]
        """
        import networkx
        return networkx.component.node_connected_component(self._nxg, vertex)

    ### Automorphism and isomorphism

    def automorphism_group(self, partition=None, translation=False):
        """
        Returns the largest subgroup of the automorphism group of the graph
        whose orbit partition is finer than the partition given. If no
        partition is given, the unit partition is used and the entire
        automorphism group is given.

        INPUT:
            translation -- if True, then output is the tuple (group, dict),
        where dict is a dictionary translating from keys == vertices to
        entries == elements of {1,2,...,n} (since permutation groups can
        currently only act on positive integers).

        EXAMPLES:
            sage: L = graphs_query.get_list_of_graphs(nodes=4)
            sage.: graphs_list.show_graphs(L)
            sage: for g in L:
            ...    G = g.automorphism_group()
            ...    G.order(), G.gens()
            (24, ((2,3), (1,2), (1,4)))
            (4, ((2,3), (1,4)))
            (2, ((1,2),))
            (8, ((2,3), (1,4), (1,3)(2,4)))
            (6, ((2,3), (1,2)))
            (6, ((1,2), (1,4)))
            (2, ((1,4)(2,3),))
            (2, ((1,2),))
            (8, ((1,3), (1,4)(2,3)))
            (4, ((2,4), (1,3)))
            (24, ((2,3), (1,2), (1,4)))

            sage: C = graphs.CubeGraph(4)
            sage: G = C.automorphism_group()
            sage: M = G.character_table()
            sage: M.determinant()
            712483534798848
            sage: G.order()
            384

            sage: D = graphs.DodecahedralGraph()
            sage: G = D.automorphism_group()
            sage: A5 = AlternatingGroup(5)
            sage: Z2 = CyclicPermutationGroup(2)
            sage: H = A5.direct_product(Z2)[0] #see documentation for direct_product to explain the [0]
            sage: G.is_isomorphic(H)
            True
        """
        if self.multiple_edges():
            raise NotImplementedError, "Search algorithm does not support multiple edges yet."
        else:
            from sage.graphs.graph_isom import search_tree, perm_group_elt
            from sage.groups.perm_gps.permgroup import PermutationGroup
            if partition is None:
                partition = [self.vertices()]
            if translation:
                a,b = search_tree(self, partition, dict=True, lab=False, dig=self.loops())
            else:
                a = search_tree(self, partition, dict=False, lab=False, dig=self.loops())
            a = PermutationGroup([perm_group_elt(aa) for aa in a])
            if translation:
                return a,b
            else:
                return a

    def is_isomorphic(self, other, proof=False):
        """
        Tests for isomorphism between self and other.

        INPUT:
            proof -- if True, then output is (a,b), where a is a boolean and b is either a map or
        None.

        EXAMPLES:
            sage: D = graphs.DodecahedralGraph()
            sage: E = D.copy()
            sage: gamma = SymmetricGroup(20).random_element()
            sage: E.relabel(gamma)
            sage: D.is_isomorphic(E)
            True

            sage: D = graphs.DodecahedralGraph()
            sage: S = SymmetricGroup(20)
            sage: gamma = S.random_element()
            sage: E = D.copy()
            sage: E.relabel(gamma)
            sage: a,b = D.is_isomorphic(E, proof=True); a
            True
            sage: import networkx
            sage: from sage.plot.plot import GraphicsArray
            sage: position_D = networkx.spring_layout(D._nxg)
            sage: position_E = {}
            sage: for vert in position_D:
            ...    position_E[b[vert]] = position_D[vert]
            sage.: GraphicsArray([D.plot(pos=position_D), E.plot(pos=position_E)]).show()
        """
        if self.multiple_edges():
            raise NotImplementedError, "Search algorithm does not support multiple edges yet."
        from sage.graphs.graph_isom import search_tree
        if proof:
            b,a = self.canonical_label(proof=True)
            d,c = other.canonical_label(proof=True)
            map = {}
            cc = c.items()
            for vert in self.vertices():
                for aa,bb in cc:
                    if bb == a[vert]:
                        map[vert] = aa
                        break
            if enum(b) == enum(d):
                return True, map
            else:
                return False, None
        else:
            from sage.graphs.graph_isom import search_tree
            b = self.canonical_label()
            d = other.canonical_label()
            return enum(b) == enum(d)

    def canonical_label(self, partition=None, proof=False):
        """
        Returns the canonical label with respect to the partition. If no
        partition is given, uses the unit partition.

        EXAMPLE:
            sage: D = graphs.DodecahedralGraph()
            sage: E = D.canonical_label(); E
            Dodecahedron: Graph on 20 vertices
            sage: D.canonical_label(proof=True)
            (Dodecahedron: Graph on 20 vertices, {0: 0, 1: 19, 2: 16, 3: 15, 4: 9, 5: 1, 6: 10, 7: 8, 8: 14, 9: 12, 10: 17, 11: 11, 12: 5, 13: 6, 14: 2, 15: 4, 16: 3, 17: 7, 18: 13, 19: 18})
            sage: D.is_isomorphic(E)
            True
        """
        if self.multiple_edges():
            raise NotImplementedError, "Search algorithm does not support multiple edges yet."
        from sage.graphs.graph_isom import search_tree
        if partition is None:
            partition = [self.vertices()]
        if proof:
            a,b,c = search_tree(self, partition, proof=True, dig=self.loops())
            return b,c
        else:
            a,b = search_tree(self, partition, dig=self.loops())
            return b

class DiGraph(GenericGraph):
    """
    Directed graph.

    INPUT:
        data -- can be any of the following:
            1. A NetworkX digraph
            2. A dictionary of dictionaries
            3. A dictionary of lists
            4. A numpy matrix or ndarray
            5. A SAGE adjacency matrix or incidence matrix
            6. pygraphviz agraph
            7. scipy sparse matrix

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
                 the DiGraph class)
        multiedges -- boolean, whether to allow multiple edges (ignored if data is
        an instance of the DiGraph class)
        format -- if None, DiGraph tries to guess- can be several values, including:
            'adjacency_matrix' -- a square SAGE matrix M, with M[i][j] equal to the number
                                  of edges \{i,j\}
            'incidence_matrix' -- a SAGE matrix, with one column C for each edge, where
                                  if C represents \{i, j\}, C[i] is -1 and C[j] is 1
        boundary -- a list of boundary vertices, if none, digraph is considered as a 'digraph
                    without boundary'
    EXAMPLES:
    1. A NetworkX digraph:
        sage: import networkx
        sage: g = networkx.DiGraph({0:[1,2,3], 2:[5]})
        sage: DiGraph(g)
        Digraph on 5 vertices

    2. A dictionary of dictionaries:
        sage: g = DiGraph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}); g
        Digraph on 5 vertices

    The labels ('x', 'z', 'a', 'out') are labels for arcs. For example, 'out' is
    the label for the arc from 2 to 5. Labels can be used as weights, if all the
    labels share some common parent.

    3. A dictionary of lists:
        sage: g = DiGraph({0:[1,2,3], 2:[5]}); g
        Digraph on 5 vertices

    4. A numpy matrix or ndarray:
        sage: import numpy
        sage: A = numpy.array([[0,1,0],[1,0,0],[1,1,0]])
        sage: DiGraph(A)
        Digraph on 3 vertices

    5. A SAGE matrix:
    Note: If format is not specified, then SAGE assumes a square matrix is an adjacency
    matrix, and a nonsquare matrix is an incidence matrix.

        A. an adjacency matrix:

        sage: M = Matrix([[0, 1, 1, 1, 0],[0, 0, 0, 0, 0],[0, 0, 0, 0, 1],[0, 0, 0, 0, 0],[0, 0, 0, 0, 0]]); M
        [0 1 1 1 0]
        [0 0 0 0 0]
        [0 0 0 0 1]
        [0 0 0 0 0]
        [0 0 0 0 0]
        sage: DiGraph(M)
        Digraph on 5 vertices

        B. an incidence matrix:

        sage: M = Matrix(6, [-1,0,0,0,1, 1,-1,0,0,0, 0,1,-1,0,0, 0,0,1,-1,0, 0,0,0,1,-1, 0,0,0,0,0]); M
        [-1  0  0  0  1]
        [ 1 -1  0  0  0]
        [ 0  1 -1  0  0]
        [ 0  0  1 -1  0]
        [ 0  0  0  1 -1]
        [ 0  0  0  0  0]
        sage: DiGraph(M)
        Digraph on 6 vertices
    """

    def __init__(self, data=None, pos=None, loops=False, format=None, boundary=None, **kwds):
        import networkx
        from sage.structure.element import is_Matrix
        if format is None:
            if is_Matrix(data):
                if data.is_square():
                    format = 'adjacency_matrix'
                else:
                    format = 'incidence_matrix'
            elif isinstance(data, DiGraph):
                self._nxg = data.networkx_graph()
            elif isinstance(data, networkx.DiGraph):
                self._nxg = networkx.XDiGraph(data, selfloops=loops, **kwds)
            elif isinstance(data, networkx.XDiGraph):
                self._nxg = data
            else:
                self._nxg = networkx.XDiGraph(data, selfloops=loops, **kwds)
        if format == 'adjacency_matrix':
            d = {}
            for i in range(data.nrows()):
                d[i] = {}
            self._nxg = networkx.XDiGraph(d, selfloops = loops, **kwds)
            e = []
            for i,j in data.nonzero_positions():
                if i == j and loops and kwds.get('multiedges',False):
                    e += [(i,j)]*int(data[i][j])
                elif i == j and loops:
                    e.append((i,j))
                elif not i == j and kwds.get('multiedges',False):
                    e += [(i,j)]*int(data[i][j])
                elif not i == j:
                    e.append((i,j))
            self._nxg.add_edges_from(e)
        elif format == 'weighted_adjacency_matrix':
            d = {}
            for i in range(data.nrows()):
                d[i] = {}
            self._nxg = networkx.XDiGraph(d, selfloops = loops, **kwds)
            e = []
            for i,j in data.nonzero_positions():
                if i != j:
                    e.append((i,j,data[i][j]))
                elif i == j and loops:
                    e.append((i,j,data[i][j]))
            self._nxg.add_edges_from(e)
        elif format == 'incidence_matrix':
            b = True
            for c in data.columns():
                d = c.dict()
                if not len(d) == 2:
                    b = False
                else:
                    k = d.keys()
                    if not d[k[0]] == -1 * d[k[1]]:
                        b = False
            if not b:
                raise AttributeError, "Incidence Matrix must have one 1 and one -1 per column."
            else:
                d = {}
                for i in range(data.nrows()):
                    d[i] = {}
                self._nxg = networkx.XDiGraph(d, selfloops = loops, **kwds)
                e = []
                for c in data.columns():
                    k = c.dict().keys()
                    if c[k[0]] == -1:
                        e.append((k[0],k[1]))
                    else:
                        e.append((k[1],k[0]))
                self._nxg.add_edges_from(e)
        if kwds.has_key('name'):
            self._nxg.name = kwds['name']
        self._pos = pos
        self._boundary = boundary

    def _repr_(self):
        name = ""
        if self.loops():
            name += "looped "
        if self.multiple_arcs():
            name += "multi-"
        name += "digraph on %d vert"%self.order()
        if self.order() == 1:
            name += "ex"
        else:
            name += "ices"
        name = name.capitalize()
        if not self._nxg.name is None and not self._nxg.name == "":
            name = self._nxg.name + ": " + name
        return name

    def copy(self):
        """
        Creates a copy of the graph.
        """
        G = DiGraph(self._nxg, name=self._nxg.name, pos=self._pos, loops=self.loops(), boundary=self._boundary)
        return G

    def to_directed(self):
        """
        Since the graph is already directed, simply returns a copy of itself.

        EXAMPLE:
            sage: DiGraph({0:[1,2,3],4:[5,1]}).to_directed()
            Digraph on 6 vertices
        """
        return self.copy()

    def to_undirected(self):
        """
        Returns an undirected version of the graph. Every arc becomes an edge.

        EXAMPLE:
            sage: D = DiGraph({0:[1,2],1:[0]})
            sage: G = D.to_undirected()
            sage: D.arcs(labels=False)
            [(0, 1), (0, 2), (1, 0)]
            sage: G.edges(labels=False)
            [(0, 1), (0, 2)]
        """
        return Graph(self._nxg.to_undirected(), pos=self._pos)

    ### General Properties

    def is_directed(self):
        """
        Since digraph is directed, returns True.
        """
        return True

    def multiple_arcs(self, new=None):
        """
        Returns whether multiple arcs are permitted in the digraph.

        INPUT:
        new -- boolean, changes whether multiple arcs are permitted in the digraph.

        EXAMPLE:
            sage: D = DiGraph(multiedges=True); D
            Multi-digraph on 0 vertices
            sage: D.multiple_arcs(False); D
            False
            Digraph on 0 vertices
        """
        if not new is None:
            if new:
                self._nxg.allow_multiedges()
            else:
                self._nxg.ban_multiedges()
        return self._nxg.multiedges

    ### Vertex Handlers

    def neighbor_iterator(self, vertex):
        """
        Return an iterator over neighbors (connected either way) of vertex.

        EXAMPLE:
            sage: D = graphs.CubeGraph(3).to_directed()
            sage: for i in D.neighbor_iterator('010'):
            ...    print i
            011
            000
            110
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
        sage: G = DiGraph()
        sage: G.add_arc((1,2), label="label")
        sage: G.networkx_graph().adj           # random output order
        {1: {2: 'label'}, 2: {}}

        sage: G = DiGraph()
        sage: G.add_arc(1,2,'label')
        sage: G.networkx_graph().adj           # random output order
        {1: {2: 'label'}, 2: {}}
        """
        self._nxg.add_edge(u, v, label)

    def add_arcs(self, arcs):
        """
        Add arcs from an iterable container.

        EXAMPLE:
            sage: G = graphs.DodecahedralGraph().to_directed()
            sage: H = DiGraph()
            sage: H.add_arcs( G.arc_iterator() ); H
            Digraph on 20 vertices
        """
        self._nxg.add_edges_from( arcs )

    def delete_arc(self, u, v=None, label=None):
        r"""
        Delete the arc from u to v, return silently if vertices or arc does
        not exist.

        INPUT:
        The following forms are all accepted:

        G.delete_arc( 1, 2 )
        G.delete_arc( (1, 2) )
        G.delete_arcs( [ (1, 2) ] )
        G.delete_arc( 1, 2, 'label' )
        G.delete_arc( (1, 2, 'label') )
        G.delete_arcs( [ (1, 2, 'label') ] )

        EXAMPLES:
            sage: D = graphs.CompleteGraph(19).to_directed()
            sage: D.size()
            342
            sage: D.delete_arc( 1, 2 )
            sage: D.delete_arc( (3, 4) )
            sage: D.delete_arcs( [ (5, 6), (7, 8) ] )
            sage: D.delete_arc( 9, 10, 'label' )
            sage: D.delete_arc( (11, 12, 'label') )
            sage: D.delete_arcs( [ (13, 14, 'label') ] )
            sage: D.size()
            335
            sage: D.has_arc( (11, 12) )
            False

            Note that even though the edge (11, 12) has no label, it still gets
            deleted: NetworkX does not pay attention to labels here.
        """
        self._nxg.delete_edge(u, v, label)

    def delete_arcs(self, arcs):
        """
        Delete arcs from an iterable container.

        EXAMPLE:
            sage: K12 = graphs.CompleteGraph(12).to_directed()
            sage: K4 = graphs.CompleteGraph(4).to_directed()
            sage: K12.size()
            132
            sage: K12.delete_arcs(K4.arc_iterator())
            sage: K12.size()
            120
        """
        self._nxg.delete_edges_from(arcs)

    def delete_multiarc(self, u, v):
        """
        Deletes all arcs from u to v.

        EXAMPLE:
            sage: D = DiGraph(multiedges=True)
            sage: D.add_arcs([(0,1), (0,1), (0,1), (1,0), (1,2), (2,3)])
            sage: D.arcs()
            [(0, 1, None), (0, 1, None), (0, 1, None), (1, 0, None), (1, 2, None), (2, 3, None)]
            sage: D.delete_multiarc( 0, 1 )
            sage: D.arcs()
            [(1, 0, None), (1, 2, None), (2, 3, None)]
        """
        self._nxg.delete_multiedge(u, v)

    def arcs(self, labels=True):
        """
        Return a list of arcs. Each arc is a triple (u,v,l) where the arc is
        from u to v, with label l.

        INPUT:
        labels -- if False, each arc is a tuple (u,v) of vertices.

        EXAMPLES:
            sage: D = graphs.DodecahedralGraph().to_directed()
            sage: D.arcs()
            [(0, 1, None), (0, 10, None), (0, 19, None), (1, 0, None), (1, 8, None), (1, 2, None), (2, 1, None), (2, 3, None), (2, 6, None), (3, 2, None), (3, 19, None), (3, 4, None), (4, 17, None), (4, 3, None), (4, 5, None), (5, 4, None), (5, 6, None), (5, 15, None), (6, 2, None), (6, 5, None), (6, 7, None), (7, 8, None), (7, 6, None), (7, 14, None), (8, 1, None), (8, 7, None), (8, 9, None), (9, 8, None), (9, 10, None), (9, 13, None), (10, 0, None), (10, 9, None), (10, 11, None), (11, 10, None), (11, 12, None), (11, 18, None), (12, 16, None), (12, 11, None), (12, 13, None), (13, 9, None), (13, 12, None), (13, 14, None), (14, 7, None), (14, 13, None), (14, 15, None), (15, 16, None), (15, 5, None), (15, 14, None), (16, 17, None), (16, 12, None), (16, 15, None), (17, 16, None), (17, 18, None), (17, 4, None), (18, 11, None), (18, 17, None), (18, 19, None), (19, 0, None), (19, 18, None), (19, 3, None)]
            sage: D.arcs(labels = False)
            [(0, 1), (0, 10), (0, 19), (1, 0), (1, 8), (1, 2), (2, 1), (2, 3), (2, 6), (3, 2), (3, 19), (3, 4), (4, 17), (4, 3), (4, 5), (5, 4), (5, 6), (5, 15), (6, 2), (6, 5), (6, 7), (7, 8), (7, 6), (7, 14), (8, 1), (8, 7), (8, 9), (9, 8), (9, 10), (9, 13), (10, 0), (10, 9), (10, 11), (11, 10), (11, 12), (11, 18), (12, 16), (12, 11), (12, 13), (13, 9), (13, 12), (13, 14), (14, 7), (14, 13), (14, 15), (15, 16), (15, 5), (15, 14), (16, 17), (16, 12), (16, 15), (17, 16), (17, 18), (17, 4), (18, 11), (18, 17), (18, 19), (19, 0), (19, 18), (19, 3)]
        """
        L = self._nxg.edges()
        if labels:
            return L
        else:
            K = []
            for u,v,l in L:
                K.append((u,v))
            return K

    def arc_boundary(self, vertices1, vertices2=None, labels=True):
        """
        Returns a list of edges (u,v,l) with u in vertices1 and v in vertices2.
        If vertices2 is None, then it is set to the complement of vertices1.

        INPUT:
        labels -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLE:
            sage: K = graphs.CompleteBipartiteGraph(9,3).to_directed()
            sage: len(K.arc_boundary( [0,1,2,3,4,5,6,7,8], [9,10,11] ))
            27
            sage: K.size()
            54

            Note that the arc boundary preserves direction: compare this example to
            the one in edge_boundary in the Graph class.
        """
        L = self._nxg.edge_boundary(vertices1, vertices2)
        if labels:
            return L
        else:
            K = []
            for u,v,l in L:
                K.append((u,v))
            return K

    def arc_iterator(self, vertices=None):
        """
        Returns an iterator over the arcs pointing out of the given
        set of vertices. If vertices is None, then returns an iterator over
        all arcs.

        EXAMPLE:
            sage: D = DiGraph( { 0 : [1,2], 1: [0] } )
            sage: for i in D.arc_iterator([0]):
            ...    print i
            (0, 1, None)
            (0, 2, None)
        """
        return self._nxg.edges_iter(vertices)

    def incoming_arc_iterator(self, vertices=None):
        """
        Return an iterator over all arriving arcs from vertices, or over all
        arcs if vertices is None.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: for a in D.incoming_arc_iterator([0]):
            ...    print a
            (1, 0, None)
            (4, 0, None)
        """
        return self._nxg.in_edges_iter(vertices)

    def incoming_arcs(self, vertices=None, labels=True):
        """
        Returns a list of arcs arriving at vertices.

        INPUT:
        labels -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.incoming_arcs([0])
            [(1, 0, None), (4, 0, None)]
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

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: for a in D.outgoing_arc_iterator([0]):
            ...    print a
            (0, 1, None)
            (0, 2, None)
            (0, 3, None)
        """
        return self._nxg.out_edges_iter(vertices)

    def outgoing_arcs(self, vertices=None, labels=True):
        """
        Returns a list of arcs departing from vertices.

        INPUT:
        labels -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.outgoing_arcs([0])
            [(0, 1, None), (0, 2, None), (0, 3, None)]
        """
        L = self._nxg.out_edges(vertices)
        if labels:
            return L
        else:
            K = []
            for u,v,l in L:
                K.append((u,v))
            return K

    def has_arc(self, u, v=None, label=None):
        """
        Returns True if there is an arc from u to v, False otherwise.

        INPUT:
        The following forms are accepted by NetworkX:

        D.has_arc( 1, 2 )
        D.has_arc( (1, 2) )
        D.has_arc( 1, 2, 'label' )

        EXAMPLE:
            sage: DiGraph().has_arc(9,2)
            False
        """
        return self._nxg.has_edge(u,v)

    def set_arc_label(self, u, v, l):
        """
        Set the arc label of a given arc.

        INPUT:
            u, v -- the vertices (and direction) of the arc
            l -- the new label

        EXAMPLE:
            sage: SD = DiGraph( { 1:[18,2], 2:[5,3], 3:[4,6], 4:[7,2], 5:[4], 6:[13,12], 7:[18,8,10], 8:[6,9,10], 9:[6], 10:[11,13], 11:[12], 12:[13], 13:[17,14], 14:[16,15], 15:[2], 16:[13], 17:[15,13], 18:[13] } )
            sage: SD.set_arc_label(1, 18, 'discrete')
            sage: SD.set_arc_label(4, 7, 'discrete')
            sage: SD.set_arc_label(2, 5, 'h = 0')
            sage: SD.set_arc_label(7, 18, 'h = 0')
            sage: SD.set_arc_label(7, 10, 'aut')
            sage: SD.set_arc_label(8, 10, 'aut')
            sage: SD.set_arc_label(8, 9, 'label')
            sage: SD.set_arc_label(8, 6, 'no label')
            sage: SD.set_arc_label(13, 17, 'k > h')
            sage: SD.set_arc_label(13, 14, 'k = h')
            sage: SD.set_arc_label(17, 15, 'v_k finite')
            sage: SD.set_arc_label(14, 15, 'v_k m.c.r.')
            sage: posn = {1:[ 3,-3],  2:[0,2],  3:[0, 13],  4:[3,9],  5:[3,3],  6:[16, 13], 7:[6,1],  8:[6,6],  9:[6,11], 10:[9,1], 11:[10,6], 12:[13,6], 13:[16,2], 14:[10,-6], 15:[0,-10], 16:[14,-6], 17:[16,-10], 18:[6,-4]}
            sage: SD.plot(pos=posn, node_size=400, color_dict={'#FFFFFF':range(1,19)}, edge_labels=True).save('search_tree.png')
        """
        if self.has_arc(u, v):
            self._nxg.adj[u][v] = l

    def arc_label(self, u, v=None):
        """
        Returns the label of an arc.

        EXAMPLE:
            sage: D = DiGraph({0 : {1 : 'edgelabel'}})
            sage: D.arcs(labels=False)
            [(0, 1)]
            sage: D.arc_label( 0, 1 )
            'edgelabel'
        """
        return self._nxg.get_edge(u,v)

    def arc_labels(self):
        """
        Returns a list of edge labels.

        EXAMPLE:
            sage: G = DiGraph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}})
            sage: G.arc_labels()
            ['x', 'z', 'a', 'out']
        """
        labels = []
        for u,v,l in self.arcs():
            labels.append(l)
        return labels

    def predecessor_iterator(self, vertex):
        """
        Returns an iterator over predecessor vertices of vertex.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: for a in D.predecessor_iterator(0):
            ...    print a
            1
            4
        """
        return self._nxg.predecessors_iter(vertex)

    def predecessors(self, vertex):
        """
        Returns a list of predecessor vertices of vertex.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.predecessors(0)
            [1, 4]
        """
        return list(self.predecessor_iterator(vertex))

    def successor_iterator(self, vertex):
        """
        Returns an iterator over successor vertices of vertex.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: for a in D.successor_iterator(0):
            ...    print a
            1
            2
            3
        """
        return self._nxg.successors_iter(vertex)

    def successors(self, vertex):
        """
        Returns a list of successor vertices of vertex.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.successors(0)
            [1, 2, 3]
        """
        return list(self.successor_iterator(vertex))

    def remove_multiple_arcs(self):
        """
        Removes all multiple arcs, retaining one arc for each.

        EXAMPLE:
            sage: D = DiGraph(multiedges=True)
            sage: D.add_arcs( [ (0,1), (0,1), (0,1), (0,1), (1,2) ] )
            sage: D.arcs(labels=False)
            [(0, 1), (0, 1), (0, 1), (0, 1), (1, 2)]
            sage: D.remove_multiple_arcs()
            sage: D.arcs(labels=False)
            [(0, 1), (1, 2)]
        """
        self._nxg.remove_all_multiedges()

    def remove_loops(self, vertices=None):
        """
        Removes loops on vertices in vertices. If vertices is None, removes all loops.

        EXAMPLE:
            sage: D = DiGraph(loops=True)
            sage: D.add_arcs( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: D.arcs(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: D.remove_loops()
            sage: D.arcs(labels=False)
            [(2, 3)]
            sage: D.loops()
            True
        """
        if vertices is None:
            self._nxg.remove_all_selfloops()
        else:
            for v in vertices:
                self.delete_multiarc(v,v)

    def loop_arcs(self):
        """
        Returns a list of all loops in the graph.

        EXAMPLE:
            sage: D = DiGraph(loops=True)
            sage: D.add_arcs( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: D.loop_arcs()
            [(0, 0, None), (1, 1, None), (2, 2, None), (3, 3, None)]
        """
        return self._nxg.selfloop_edges()

    def number_of_loops(self):
        """
        Returns the number of arcs that are loops.

        EXAMPLE:
            sage: D = DiGraph(loops=True)
            sage: D.add_arcs( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: D.arcs(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: D.number_of_loops()
            4
        """
        return self._nxg.number_of_selfloops()

    ### Degree functions

    def degree(self, vertices=None, labels=False):
        """
        Gives the degree (in + out) of a vertex or of vertices.

        INPUT:
        vertices -- If vertices is a single vertex, returns the number of
        neighbors of vertex. If vertices is an iterable container of vertices,
        returns a list of degrees. If vertices is None, same as listing all vertices.
        labels -- see OUTPUT

        OUTPUT:
        Single vertex- an integer. Multiple vertices- a list of integers. If
        labels is True, then returns a dictionary mapping each vertex to
        its degree.

        EXAMPLES:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.degree(vertices = [0,1,2], labels=True)
            {0: 5, 1: 4, 2: 3}
            sage: D.degree()
            [5, 4, 3, 3, 3, 2]
        """
        return self._nxg.degree(vertices, with_labels=labels)

    def degree_iterator(self, vertices=None, labels=False):
        """
        INPUT:
        labels=False:
            returns an iterator over degrees.
        labels=True:
            returns an iterator over tuples (vertex, degree).
        vertices -- if specified, restrict to this subset.

        EXAMPLE:
            sage: D = graphs.Grid2dGraph(2,4).to_directed()
            sage: for i in D.degree_iterator():
            ...    print i
            6
            6
            4
            4
            6
            4
            4
            6
            sage: for i in D.degree_iterator(labels=True):
            ...    print i
            ((0, 1), 6)
            ((1, 2), 6)
            ((0, 0), 4)
            ((0, 3), 4)
            ((0, 2), 6)
            ((1, 3), 4)
            ((1, 0), 4)
            ((1, 1), 6)
        """
        return self._nxg.degree_iter(vertices, with_labels=labels)

    def in_degree(self, vertices=None, labels=False):
        """
        Same as degree, but for in degree.

        EXAMPLES:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.in_degree(vertices = [0,1,2], labels=True)
            {0: 2, 1: 2, 2: 2}
            sage: D.in_degree()
            [2, 2, 2, 2, 1, 1]
        """
        return self._nxg.in_degree(vertices, with_labels=labels)

    def in_degree_iterator(self, vertices=None, labels=False):
        """
        Same as degree_iterator, but for in degree.

        EXAMPLES:
            sage: D = graphs.Grid2dGraph(2,4).to_directed()
            sage: for i in D.in_degree_iterator():
            ...    print i
            3
            3
            2
            2
            3
            2
            2
            3
            sage: for i in D.in_degree_iterator(labels=True):
            ...    print i
            ((0, 1), 3)
            ((1, 2), 3)
            ((0, 0), 2)
            ((0, 3), 2)
            ((0, 2), 3)
            ((1, 3), 2)
            ((1, 0), 2)
            ((1, 1), 3)
        """
        return self._nxg.in_degree_iter(vertices, with_labels=labels)

    def out_degree(self, vertices=None, labels=False):
        """
        Same as degree, but for out degree.

        EXAMPLES:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.out_degree(vertices = [0,1,2], labels=True)
            {0: 3, 1: 2, 2: 1}
            sage: D.out_degree()
            [3, 2, 1, 1, 2, 1]
        """
        return self._nxg.out_degree(vertices, with_labels=labels)

    def out_degree_iterator(self, vertices=None, labels=False):
        """
        Same as degree_iterator, but for out degree.

        EXAMPLES:
            sage: D = graphs.Grid2dGraph(2,4).to_directed()
            sage: for i in D.out_degree_iterator():
            ...    print i
            3
            3
            2
            2
            3
            2
            2
            3
            sage: for i in D.out_degree_iterator(labels=True):
            ...    print i
            ((0, 1), 3)
            ((1, 2), 3)
            ((0, 0), 2)
            ((0, 3), 2)
            ((0, 2), 3)
            ((1, 3), 2)
            ((1, 0), 2)
            ((1, 1), 3)
        """
        return self._nxg.out_degree_iter(vertices, with_labels=labels)

    ### Representations

    def adjacency_matrix(self, sparse=True):
        """
        Returns the adjacency matrix of the digraph. Each vertex is
        represented by its position in the list returned by the vertices()
        function.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.adjacency_matrix()
            [0 1 1 1 0 0]
            [1 0 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [1 0 0 0 0 1]
            [0 1 0 0 0 0]
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

    def incidence_matrix(self, sparse=True):
        """
        Returns an incidence matrix of the graph. Each row is a vertex, and
        each column is an edge.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.incidence_matrix()
            [-1 -1 -1  1  0  0  0  1  0  0]
            [ 1  0  0 -1 -1  0  0  0  0  1]
            [ 0  1  0  0  1 -1  0  0  0  0]
            [ 0  0  1  0  0  1 -1  0  0  0]
            [ 0  0  0  0  0  0  1 -1 -1  0]
            [ 0  0  0  0  0  0  0  0  1 -1]
        """
        from sage.matrix.constructor import matrix
        from copy import copy
        n = len(self._nxg.adj)
        verts = self.vertices()
        d = [0]*n
        cols = []
        for i, j, l in self.arc_iterator():
            col = copy(d)
            i = verts.index(i)
            j = verts.index(j)
            col[i] = -1
            col[j] = 1
            cols.append(col)
        return matrix(cols, sparse=sparse).transpose()

    ### Contruction

    def reverse(self):
        """
        Returns a copy of digraph with arcs reversed in direction.

        TODO: results in error because of the following NetworkX bug (0.33) - trac #92

        EXAMPLES:
            sage: import networkx
            sage: D = networkx.XDiGraph({ 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] })
            sage: D.reverse()
            Traceback (most recent call last):
            ...
            ValueError: too many values to unpack
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

        EXAMPLES:
            sage: D = graphs.CompleteGraph(9).to_directed()
            sage: H = D.subgraph([0,1,2]); H
            Digraph on 3 vertices
            sage: D
            Digraph on 9 vertices
            sage: K = D.subgraph([0,1,2], inplace=True); K
            Subgraph of (None): Digraph on 3 vertices
            sage: D
            Subgraph of (None): Digraph on 3 vertices
            sage: D is K
            True
        """
        if inplace:
            self._nxg = self._nxg.subgraph(vertices, inplace, create_using)
            return self
        else:
            NXG = self._nxg.subgraph(vertices, inplace, create_using)
            return DiGraph(NXG)

    ### Visualization

    def plot3d(self, bgcolor=(1,1,1), vertex_color=(1,0,0), arc_color=(0,0,0), pos3d=None):
        """
        Plots the graph using Tachyon, and returns a Tachyon object containing
        a representation of the graph.

        INPUT:
            bgcolor
            vertex_color
            arc_color
            pos3d -- a position dictionary for the vertices

        NOTE:
            The weaknesses of the NetworkX spring layout are illustrated even further in the
            case of digraphs: my guess is that digraphs weren't even considered in the authoring
            of this algorithm. The following example illustrates this.

        EXAMPLE:
            # This is a running example

            # A directed version of the dodecahedron
            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []} )

            # If I use an undirected version of my graph, the output is as expected
            sage: import networkx
            sage: pos3d=networkx.spring_layout(graphs.DodecahedralGraph()._nxg, dim=3)
            sage: D.plot3d(pos3d=pos3d).save('sage.png') # long time

            # However, if I use the directed version, everything gets skewed bizarrely:
            sage: D.plot3d().save('sage.png') # long time
        """
        TT, pos3d = tachyon_vertex_plot(self, bgcolor=bgcolor, vertex_color=vertex_color, pos3d=pos3d)
        TT.texture('arc', ambient=0.1, diffuse=0.9, specular=0.03, opacity=1.0, color=arc_color)
        for u,v,l in self.arcs():
            TT.fcylinder( (pos3d[u][0],pos3d[u][1],pos3d[u][2]),\
                          (pos3d[v][0],pos3d[v][1],pos3d[v][2]), .02,'arc')
            TT.fcylinder( (0.25*pos3d[u][0] + 0.75*pos3d[v][0],\
                           0.25*pos3d[u][1] + 0.75*pos3d[v][1],\
                           0.25*pos3d[u][2] + 0.75*pos3d[v][2],),
                          (pos3d[v][0],pos3d[v][1],pos3d[v][2]), .0325,'arc')
        return TT

    def show3d(self, bgcolor=(1,1,1), vertex_color=(1,0,0), edge_color=(0,0,0), pos3d=None, **kwds):
        """
        Plots the graph using Tachyon, and shows the resulting plot.

        INPUT:
            bgcolor -- background color
            vertex_color -- vertex color
            edge_color -- edge color
            (pos3d -- currently ignored, pending GSL random point distribution in sphere...)

        NOTE:
            The weaknesses of the NetworkX spring layout are illustrated even further in the
            case of digraphs: my guess is that digraphs weren't even considered in the authoring
            of this algorithm. The following example illustrates this.

        EXAMPLE:
            # This is a running example

            # A directed version of the dodecahedron
            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []} )

            # If I use an undirected version of my graph, the output is as expected
            sage: import networkx
            sage: pos3d=networkx.spring_layout(graphs.DodecahedralGraph()._nxg, dim=3)
            sage: D.plot3d(pos3d=pos3d).save('sage.png') # long time

            # However, if I use the directed version, everything gets skewed bizarrely:
            sage: D.plot3d().save('sage.png') # long time
        """
        self.plot3d(bgcolor=bgcolor, vertex_color=vertex_color, arc_color=arc_color).show(**kwds)

    ### TODO: Connected components?

    ### Automorphism and isomorphism

    def automorphism_group(self, partition=None, translation=False):
        """
        Returns the largest subgroup of the automorphism group of the digraph
        whose orbit partition is finer than the partition given. If no
        partition is given, the unit partition is used and the entire
        automorphism group is given.

        INPUT:
            translation -- if True, then output is the tuple (group, dict),
        where dict is a dictionary translating from keys == vertices to
        entries == elements of {1,2,...,n} (since permutation groups can
        currently only act on positive integers).

        EXAMPLES:
            TODO
        """
        if self.multiple_arcs():
            raise NotImplementedError, "Search algorithm does not support multiple edges yet."
        else:
            from sage.graphs.graph_isom import search_tree, perm_group_elt
            from sage.groups.perm_gps.permgroup import PermutationGroup
            if partition is None:
                partition = [self.vertices()]
            if translation:
                a,b = search_tree(self, partition, dict=True, lab=False, dig=True)
            else:
                a = search_tree(self, partition, dict=False, lab=False, dig=True)
            a = PermutationGroup([perm_group_elt(aa) for aa in a])
            if translation:
                return a,b
            else:
                return a

    def is_isomorphic(self, other, proof=False):
        """
        Tests for isomorphism between self and other.

        INPUT:
            proof -- if True, then output is (a,b), where a is a boolean and b is either a map or
        None.

        EXAMPLES:
            TODO
        """
        if self.multiple_arcs():
            raise NotImplementedError, "Search algorithm does not support multiple edges yet."
        from sage.graphs.graph_isom import search_tree
        if proof:
            b,a = self.canonical_label(proof=True)
            d,c = other.canonical_label(proof=True)
            if enum(b) == enum(d):
                map = {}
                cc = c.items()
                for vert in self.vertices():
                    for aa,bb in cc:
                        if bb == a[vert]:
                            map[vert] = aa
                            break
                return True, map
            else:
                return False, None
        else:
            from sage.graphs.graph_isom import search_tree
            b = self.canonical_label()
            d = other.canonical_label()
            return enum(b) == enum(d)

    def canonical_label(self, partition=None, proof=False):
        """
        Returns the canonical label with respect to the partition. If no
        partition is given, uses the unit partition.

        EXAMPLE:
            TODO
        """
        if self.multiple_arcs():
            raise NotImplementedError, "Search algorithm does not support multiple edges yet."
        from sage.graphs.graph_isom import search_tree
        if partition is None:
            partition = [self.vertices()]
        if proof:
            a,b,c = search_tree(self, partition, proof=True, dig=True)
            return b,c
        else:
            a,b = search_tree(self, partition, dig=True)
            return b

def tachyon_vertex_plot(g, bgcolor=(1,1,1), vertex_color=(1,0,0), pos3d=None):
    import networkx
    from math import sqrt
    from sage.plot.tachyon import Tachyon
    c = [0,0,0]
    r = []
    verts = g.vertices()
    spring = False
    if pos3d is None:
        pos3d = graph_fast.spring_layout_fast(g, dim=3)
    try:
        for v in verts:
            c[0] += pos3d[v][0]
            c[1] += pos3d[v][1]
            c[2] += pos3d[v][2]
    except KeyError:
        raise KeyError, "Oops! You haven't specified positions for all the vertices."
    order = g.order()
    c[0] = c[0]/order
    c[1] = c[1]/order
    c[2] = c[2]/order
    for v in verts:
        pos3d[v][0] = pos3d[v][0] - c[0]
        pos3d[v][1] = pos3d[v][1] - c[1]
        pos3d[v][2] = pos3d[v][2] - c[2]
        r.append(abs(sqrt((pos3d[v][0])**2 + (pos3d[v][1])**2 + (pos3d[v][2])**2)))
    r = max(r)
    for v in verts:
        pos3d[v][0] = pos3d[v][0]/r
        pos3d[v][1] = pos3d[v][1]/r
        pos3d[v][2] = pos3d[v][2]/r
    TT = Tachyon(camera_center=(1.4,1.4,1.4), antialiasing=13)
    TT.light((4,3,2), 0.02, (1,1,1))
    TT.texture('node', ambient=0.1, diffuse=0.9, specular=0.03, opacity=1.0, color=vertex_color)
    TT.texture('bg', ambient=1, diffuse=1, specular=0, opacity=1.0, color=bgcolor)
    TT.plane((-1.6,-1.6,-1.6), (1.6,1.6,1.6), 'bg')
    for v in verts:
        TT.sphere((pos3d[v][0],pos3d[v][1],pos3d[v][2]), .06, 'node')
    return TT, pos3d

def enum(graph, quick=False):
    """
    Used for isomorphism checking.

    INPUT:
        quick -- now we know that the vertices are 0,1,...,n-1

    EXAMPLES:
        sage: from sage.graphs.graph import enum
        sage: enum(graphs.DodecahedralGraph())
        646827340296833569479885332381965103655612500627043016896502674924517797573929148319427466126170568392555309533861838850L
        sage: enum(graphs.MoebiusKantorGraph())
        29627597595494233374689380190219099810725571659745484382284031717525232288040L
        sage: enum(graphs.FlowerSnark())
        645682215283153372602620320081348424178216159521280462146968720908564261127120716040952785862033320307812724373694972050L
        sage: enum(graphs.CubeGraph(3))
        6100215452666565930L
        sage: enum(graphs.CubeGraph(4))
        31323620658472264895128471376615338141839885567113523525061169966087480352810L
        sage: enum(graphs.CubeGraph(5))
        56178607138625465573345383656463935701397275938329921399526324254684498525419117323217491887221387354861371989089284563861938014744765036177184164647909535771592043875566488828479926184925998575521710064024379281086266290501476331004707336065735087197243607743454550839234461575558930808225081956823877550090L
        sage: enum(graphs.CubeGraph(6))
        17009933328531023098235951265708015080189260525466600242007791872273951170067729430659625711869482140011822425402311004663919203785115296476561677814427201708237805402966561863692388687547518491537427897858240566495945005294876576523289206747123399572439707189803821880345487300688962557172856432472391025950779306221469432919735886988596366979797317084123956762362685536557279604675024249987913439836592296340787741671304722135394212035449285260308821361913500205796919484488876249630521666898413890977354122918711285458724686283296097840711521153201188450783978019001984591992381570913097193343212274205747843852376395748070926193308573472616983062165141386183945049871456376379041631456999916186868438148001405477879591035696239287238767746380404501285533026300096772164676955425088646172718295360584249310479706751274583871684827338312536787740914529353458829503642591918761588296961192261166874864565050490306157300749101788751129640698534818737753110920871293122429238702542726347017441416450649382146313791818349648006634724962025571237208317435310419071153813687071275479812184286929976456778629116002591936357623320676067640749567446551071011889378108453641887998273235139859889734259803684619153716302058849155208478850L
    """
    enumeration = 0
    n = graph.order()
    if quick:
        if isinstance(graph, Graph):
            for i, j, l in graph.edge_iterator():
                enumeration += 1 << ((n-(i+1))*n + n-(j+1))
                enumeration += 1 << ((n-(j+1))*n + n-(i+1))
        elif isinstance(graph, DiGraph):
            for i, j, l in graph.arc_iterator():
                enumeration += 1 << ((n-(i+1))*n + n-(j+1))
        return enumeration
    M = graph.am()
    for i, j in M.nonzero_positions():
        enumeration += 1 << ((n-(i+1))*n + n-(j+1))
    return enumeration




