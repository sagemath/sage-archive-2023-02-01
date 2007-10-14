r"""
Graph Theory

This module implements many graph theoretic operations and concepts.

AUTHOR:
    -- Robert L. Miller (2006-10-22): initial version
    -- William Stein (2006-12-05): Editing
    -- Robert L. Miller (2007-01-13): refactoring, adjusting for
        NetworkX-0.33, fixed plotting bugs
                        (2007-01-23): basic tutorial, edge labels, loops,
                                      multiple edges and arcs
                        (2007-02-07): graph6 and sparse6 formats, matrix input
    -- Emily Kirkmann (2007-02-11): added graph_border option to plot and show
    -- Robert L. Miller (2007-02-12): vertex color-maps, graph boundaries,
        graph6 helper functions in Cython
                        SAGE Days 3 (2007-02-17--21): 3d plotting in Tachyon
                        (2007-02-25): display a partition
                        (2007-02-28): associate arbitrary objects to vertices,
        edge and arc label display (in 2d), edge coloring
                        (2007-03-21): Automorphism group, isomorphism check,
        canonical label
                        (2007-06-07--09): NetworkX function wrapping
    -- Michael W. Hansen (2007-06-09): Topological sort generation
    -- Emily Kirkman, Robert L. Miller SAGE Days 4: Finished wrapping NetworkX
    -- Emily Kirkman (2007-07-21): Genus (including circular planar, all
        embeddings and all planar embeddings), all paths, interior paths
    -- Bobby Moretti (2007-08-12): fixed up plotting of graphs with
       edge colors differentiated by label
    -- Jason Grout (2007-09-25): Added functions, bug fixes, and
       general enhancements

\subsection{Graph Format}

\subsubsection{The SAGE Graph Class: NetworkX plus}

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
            is redundancy: for example, the dictionary containing the entry
            \verb|1: {2: None}| implies it must contain \verb|{2: {1: None}|.
            The innermost entry of \var{None} is related to edge labeling
            (see section \ref{Graph:labels}).

            \subsubsection{Supported formats}


            SAGE Graphs can be created from a wide range of inputs. A few examples are
            covered here.

            \begin{itemize}

                  \item NetworkX dictionary format:

                sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], \
                      5: [7, 8], 6: [8,9], 7: [9]}
                sage: G = Graph(d); G
                Graph on 10 vertices
                sage: G.plot().save('sage.png')    # or G.show()

                    \item A NetworkX graph:

                sage: K = networkx.complete_bipartite_graph(12,7)
                sage: G = Graph(K)
                sage: G.degree()
                [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 12, 12, 12, 12, 12, 12, 12]

                \item graph6 or sparse6 format:

                sage: s = ':I`AKGsaOs`cI]Gb~'
                sage: G = Graph(s); G
                Looped multi-graph on 10 vertices
                sage: G.plot().save('sage.png')    # or G.show()

                \item adjacency matrix In an adjacency matrix, each column and each row represent
                a vertex. If a 1 shows up in row i, column j, there is an edge (i,j).

                sage: M = Matrix([(0,1,0,0,1,1,0,0,0,0),(1,0,1,0,0,0,1,0,0,0), \
                (0,1,0,1,0,0,0,1,0,0), (0,0,1,0,1,0,0,0,1,0),(1,0,0,1,0,0,0,0,0,1), \
                (1,0,0,0,0,0,0,1,1,0), (0,1,0,0,0,0,0,0,1,1),(0,0,1,0,0,1,0,0,0,1), \
                (0,0,0,1,0,1,1,0,0,0), (0,0,0,0,1,0,1,1,0,0)])
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
                sage: G.plot().save('sage.png')    # or G.show()

                \item incidence matrix: In an incidence matrix, each row represents a vertex
                and each column reprensents an edge.

                sage: M = Matrix([(-1,0,0,0,1,0,0,0,0,0,-1,0,0,0,0), \
                (1,-1,0,0,0,0,0,0,0,0,0,-1,0,0,0),(0,1,-1,0,0,0,0,0,0,0,0,0,-1,0,0), \
                (0,0,1,-1,0,0,0,0,0,0,0,0,0,-1,0),(0,0,0,1,-1,0,0,0,0,0,0,0,0,0,-1), \
                (0,0,0,0,0,-1,0,0,0,1,1,0,0,0,0),(0,0,0,0,0,0,0,1,-1,0,0,1,0,0,0), \
                (0,0,0,0,0,1,-1,0,0,0,0,0,1,0,0),(0,0,0,0,0,0,0,0,1,-1,0,0,0,1,0), \
                (0,0,0,0,0,0,1,-1,0,0,0,0,0,0,1)])
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
                sage: G.plot().save('sage.png')    # or G.show()

        \end{itemize}

        \subsection{Generators}

        For some commonly used graphs to play with, type

            sage.: graphs.

        and hit \kbd{tab}. Most of these graphs come with their own custom plot, so you
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

            sage: G = GraphDatabase()
            sage: L = G.get_list(num_vertices=7, diameter=5)
            sage.: graphs_list.show_graphs(L)

            \subsection{Labels}\label{Graph:labels}

        Each vertex can have any hashable object as a label. These are things like
        strings, numbers, and tuples. Each edge is given a default label of \var{None}, but
        if specified, edges can have any label at all. Edges between vertices $u$ and $v$ are
        represented typically as \verb|(u, v, l)|, where \var{l} is the label for the edge.

        Note that vertex labels themselves cannot be mutable items:

            sage: M = Matrix( [[0,0],[0,0]] )
            sage: G = Graph({ 0 : { M : None } })
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        However, if one wants to define a dictionary, with the same keys and arbitrary objects
        for entries, one can make that association:

            sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), \
                  2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
            sage: d[2]
            Moebius-Kantor Graph: Graph on 16 vertices
            sage: T = graphs.TetrahedralGraph()
            sage: T.vertices()
            [0, 1, 2, 3]
            sage: T.associate(d)
            sage: T.obj(1)
            Flower Snark: Graph on 20 vertices

        \subsection{Database}

        There is a database available for searching for graphs that satisfy a certain set
        of parameters, including number of vertices and edges, density, maximum and minimum
        degree, diameter, radius, and connectivity. If you wish to search a database of
        graphs by parameter, type

            sage.: graphs_query.

        and hit \kbd{tab}.

            sage: graphs_query = GraphDatabase()
            sage: L = graphs_query.get_list(num_vertices=7, diameter=5)
            sage.: graphs_list.show_graphs(L)

        \subsection{Visualization}

        To see a graph G you are working with, right now there are two main options.
        You can view the graph in two dimensions via matplotlib with \method{show()}.

            sage: G = graphs.RandomGNP(15,.3)
            sage.: G.show()

        Or you can view it in three dimensions via Tachyon with \method{show3d()}.

            sage.: G.show3d()

            \note{Many functions are passed directly on to NetworkX.
              In these cases, the documentation is based on the
              NetworkX docs.}

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
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ

class GenericGraph(SageObject):
    """
    Base class for graphs and digraphs.

    """

    def __cmp__(self, other):
        """
        Comparison of self and other. Must be in the same class, have the same
        settings for loops and multiedges, output the same vertex list (in order)
        and the same adjacency matrix.

        Note that this is _not_ an isomorphism test.

        Note that the less-than and greater-than value returned here
        doesn't mean much.  The equality test is the useful thing.

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
            sage: G == H
            False
            sage: G = graphs.RandomGNP(9,.3).to_directed()
            sage: H = graphs.RandomGNP(9,.3).to_directed()
            sage: G == H # random output (most often false)
            False

        """
        # If the graphs have different properties, they are not equal.
        if type(self) != type(other):
            return 1
        elif self.loops() != other.loops():
            return 1
        elif self.multiple_edges() != other.multiple_edges():
            return 1

        # If the vertices have different labels, the graphs are not equal.
        if self.vertices() != other.vertices():
            return 1

        # Check that the edges are the same.
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
        Return an iterator over the vertices which allows \code{for v in G} syntax.

        """
        return self.vertex_iterator()

    def __len__(self):
        """
        len(G) returns the number of vertices in G.

        """
        return len(self._nxg.adj)

    def __str__(self):
        """
        str(G) returns the name of the graph, unless the name is the empty string, in
        which case it returns the default representation.

        """
        if self._nxg.name != '':
            return self._nxg.name
        else:
            return repr(self)

    def _latex_(self):
        """
        To include a graph in LaTeX document, see function
        Graph.write_to_eps().

        """
        raise NotImplementedError('To include a graph in LaTeX document, see function Graph.write_to_eps().')

    def _matrix_(self, R=None):
        """

        Returns the adjacency matrix of the graph over the specified ring.

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
            x^3 * (x^2 - 6)
        """
        if R is None:
            return self.am()
        else:
            return self.am().change_ring(R)

    def networkx_graph(self):
        """
        Creates a new NetworkX graph from the SAGE graph.

        EXAMPLE:
            sage: G = graphs.TetrahedralGraph()
            sage: N = G.networkx_graph()
            sage: type(N)
            <class 'networkx.xgraph.XGraph'>

        Note that this returns a copy of the actual internal object,
        not the actual internal networkX object.

            sage: G = graphs.TetrahedralGraph()
            sage: N = G.networkx_graph()
            sage: G._nxg is N
            False

        """
        return self._nxg.copy()

    def networkx_info(self, vertex=None):
        """
        Returns NetworkX information about the graph or the given vertex.

        """
        self._nxg.info(vertex)

    def __get_pos__(self):
        """
        Returns the position dictionary, a dictionary specifying the coordinates
        of each vertex.
        """
        return self._pos

    def __set_pos__(self, pos):
        """
        Sets the position dictionary, a dictionary specifying the
        coordinates of each vertex.
        """
        self._pos = pos

    ### General properties

    def name(self, new=None):
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
            sage: G.name(new=""); G
            Graph on 10 vertices
            sage: G.name()

        """
        if new is None:
            if self._nxg.name != '':
                return self._nxg.name
            else:
                return None

        if not isinstance(new, str):
            raise TypeError, "New name must be a string."

        self._nxg.name = new
        if new != '':
            return self._nxg.name
        else:
            return None

    def loops(self, new=None):
        """
        Returns whether loops are permitted in the graph.

        INPUT:
        new -- boolean, changes whether loops are permitted in the graph.

        EXAMPLE:
            sage: G = Graph(); G
            Graph on 0 vertices
            sage: G.loops(True)

            sage: D = DiGraph(); D
            Digraph on 0 vertices
            sage: D.loops()
            False
            sage: D.loops(True)
            sage: D.loops()
            True

        """
        if new is not None:
            if new:
                self._nxg.allow_selfloops()
            else:
                self._nxg.ban_selfloops()
        else: return self._nxg.selfloops

    def multiple_edges(self, new=None):
        """
        Returns whether multiple edges are permitted in the (di)graph.

        INPUT:
        new -- boolean. If specified, changes whether multiple edges are
        permitted in the graph.

        EXAMPLE:
            sage: G = Graph(multiedges=True); G
            Multi-graph on 0 vertices
            sage: G.multiple_edges(False); G
            Graph on 0 vertices
            sage: D = DiGraph(multiedges=True); D
            Multi-digraph on 0 vertices
            sage: D.multiple_edges(False); D
            Digraph on 0 vertices

        """
        if new is not None:
            if new:
                self._nxg.allow_multiedges()
            else:
                self._nxg.ban_multiedges()
        else: return self._nxg.multiedges

    def density(self):
        """
        Returns the density (number of edges divided by number of possible
        edges).

        In the case of a multigraph, raises an error, since there is an
        infinite number of possible edges.

        EXAMPLE:
            sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], 5: [7, 8], 6: [8,9], 7: [9]}
            sage: G = Graph(d); G.density()
            1/3
            sage: G = Graph({0:[1,2], 1:[0] }); G.density()
            2/3
            sage: G = DiGraph({0:[1,2], 1:[0] }); G.density()
            1/2

        Note that there are more possible edges on a looped graph:
            sage: G.loops(True)
            sage: G.density()
            1/3

        """
        if self.multiple_edges():
            raise TypeError("Density is not well-defined for multigraphs.")
        from sage.rings.rational import Rational
        n = self.order()
        if self.loops():
            if self.is_directed():
                return Rational(self.size())/Rational(n**2)
            else:
                return Rational(self.size())/Rational((n**2 + n)/2)
        else:
            if self.is_directed():
                return Rational(self.size())/Rational((n**2 - n))
            else:
                return Rational(self.size())/Rational((n**2 - n)/2)

    def to_simple(self):
        """
        Returns a simple version of itself (i.e., undirected and loops
        and multiple edges are removed).

        EXAMPLE:
            sage: G = DiGraph(loops=True,multiedges=True)
            sage: G.add_edges( [ (0,0), (1,1), (2,2), (2,3), (2,3), (3,2) ] )
            sage: G.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (2, 3), (3, 2)]
            sage: H=G.to_simple()
            sage: H.edges(labels=False)
            [(2, 3)]
            sage: H.is_directed()
            False
            sage: H.loops()
            False
            sage: H.multiple_edges()
            False

        """
        g=self.to_undirected()
        g.loops(False)
        g.multiple_edges(False)
        return g


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
        Creates an isolated vertex.  If the vertex already exists, then nothing is done.

        INPUT:
        n -- Name of the new vertex. If no name is specified, then the vertex
        will be represented by the least integer not already representing a
        vertex. Name must be an immutable object.

        As it is implemented now, if a graph $G$ has a large number of
        vertices with numeric labels, then G.add_vertex() could
        potentially be slow.

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
        Vertices that already exist in the graph will not be added again.

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
        Deletes vertex, removing all incident edges.  Deleting a non-existant
        vertex will raise an exception.

        EXAMPLES:
            sage: G = graphs.WheelGraph(9)
            sage: G.delete_vertex(0); G.save('sage.png')

            sage: D = DiGraph({0:[1,2,3,4,5],1:[2],2:[3],3:[4],4:[5],5:[1]})
            sage: D.delete_vertex(0); D
            Digraph on 5 vertices
            sage: D.vertices()
            [1, 2, 3, 4, 5]
            sage: D.delete_vertex(0)
            Traceback (most recent call last):
            ...
            NetworkXError: node 0 not in graph

        """
        self._nxg.delete_node(vertex)

    def delete_vertices(self, vertices):
        """
        Remove vertices from the (di)graph taken from an iterable container of
        vertices.  Deleting a non-existant vertex will raise an exception.

        EXAMPLE:
            sage: D = DiGraph({0:[1,2,3,4,5],1:[2],2:[3],3:[4],4:[5],5:[1]})
            sage: D.delete_vertices([1,2,3,4,5]); D
            Digraph on 1 vertex
            sage: D.vertices()
            [0]
            sage: D.delete_vertices([1])
            Traceback (most recent call last):
            ...
            NetworkXError: node 1 not in graph

        """
        self._nxg.delete_nodes_from(vertices)

    def get_boundary(self):
        return self._boundary

    def set_boundary(self, boundary):
        if isinstance(boundary,list):
            self._boundary = boundary

    def vertex_boundary(self, vertices1, vertices2=None):
        """
        Returns a list of all vertices in the external boundary of vertices1,
        intersected with vertices2. If vertices2 is None, then vertices2 is the
        complement of vertices1.  This is much faster if vertices1 is smaller than
        vertices2.

        The external boundary of a set of vertices is the union of the
        neighborhoods of each vertex in the set.  Note that in this
        implementation, since vertices2 defaults to the complement of
        vertices1, if a vertex $v$ has a loop, then vertex_boundary(v)
        will not contain $v$.

        EXAMPLE:
            sage: G = graphs.CubeGraph(4)
            sage: l = ['0111', '0000', '0001', '0011', '0010', '0101', '0100', '1111', '1101', '1011', '1001']
            sage: G.vertex_boundary(['0000', '1111'], l)
            ['0010', '0100', '0001', '0111', '1011', '1101']

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
        try:
            for v in vertex_dict:
                self._assoc[v] = vertex_dict[v]
        except:
            self._assoc = {}
            for v in vertex_dict:
                self._assoc[v] = vertex_dict[v]

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
        try:
            return self._assoc[vertex]
        except:
            return None

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
        Empties the graph of vertices and edges and removes name,
        boundary, associated objects, and position information.

        EXAMPLE:
            sage: G=graphs.CycleGraph(4); G.associate({0:'vertex0'})
            sage: G.order(); G.size()
            4
            4
            sage: len(G._pos)
            4
            sage: G.name()
            'Cycle graph'
            sage: G.obj(0)
            'vertex0'
            sage: G.clear()
            sage: G.order(); G.size()
            0
            0
            sage: len(G._pos)
            0
            sage: G.name()
            sage: G.obj(0)

        """
        self._nxg.clear()
        self._pos=[]
        self._boundary=[]
        self._assoc=None

    def neighbors(self, vertex):
        """
        Return a list of neighbors (in and out if directed) of vertex.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: list(P.neighbor_iterator(3))
            [8, 2, 4]

        """
        return list(self.neighbor_iterator(vertex))

    def random_subgraph(self, p, inplace=False, create_using=None):
        """
        Return a random subgraph that contains each vertex with prob. p.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: P.random_subgraph(.25) # random
            Subgraph of (Petersen graph): Graph on 3 vertices

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

        INPUT:
            vertices -- iterated vertices are these intersected with the
                vertices of the (di)graph

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: for v in P.vertex_iterator():
            ...    print v
            ...
            0
            1
            2
            3
            4
            5
            6
            7
            8
            9

        Note that since the intersection option is available, the
        vertex_iterator() function is sub-optimal, speedwise, but note the
        following optimization:
            sage.: timeit V = P.vertices()
            100000 loops, best of 3: 8.85 [micro]s per loop
            sage.: timeit V = list(P.vertex_iterator())
            100000 loops, best of 3: 5.74 [micro]s per loop
            sage.: timeit V = list(P._nxg.adj.iterkeys())
            100000 loops, best of 3: 3.45 [micro]s per loop

        In other words, if you want a fast vertex iterator, call the dictionary
        directly.

        """
        return self._nxg.prepare_nbunch(vertices)

    def vertices(self, boundary_first=False):
        """
        Return a list of the vertices.

        INPUT:
            boundary_first -- Return the boundary vertices first.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: P.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        Note that the output of the vertices() function is always sorted. This
        is sub-optimal, speedwise, but note the following optimizations:
            sage.: timeit V = P.vertices()
            100000 loops, best of 3: 8.85 [micro]s per loop
            sage.: timeit V = list(P.vertex_iterator())
            100000 loops, best of 3: 5.74 [micro]s per loop
            sage.: timeit V = list(P._nxg.adj.iterkeys())
            100000 loops, best of 3: 3.45 [micro]s per loop

        In other words, if you want a fast vertex iterator, call the dictionary
        directly.

        """
        if not boundary_first:
            return sorted(list(self.vertex_iterator()))

        bdy_verts = []
        int_verts = []
        for v in self.vertex_iterator():
            if v in self._boundary:
                bdy_verts.append(v)
            else:
                int_verts.append(v)
        return sorted(bdy_verts) + sorted(int_verts)

    def relabel(self, perm, quick=False, inplace=True):
        r"""
        Uses a dictionary, list, or permutation to relabel the (di)graph.
        If perm is a dictionary d, each old vertex v is a key in the
        dictionary, and its new label is d[v].

        If perm is a list, we think of it as a map $i \mapsto perm[i]$
        with the assumption that the vertices are $\{0,1,...,n-1\}$.

        If perm is a permutation, the permutation is simply applied to
        the graph, under the assumption that the vertices are
        $\{0,1,...,n-1\}$.  The permutation acts on the set
        $\{1,2,...,n\}$, where we think of $n = 0$.


        INPUT:
            quick -- if True, simply return the enumeration of the new graph
        without constructing it. Requires that perm is of type list.

        EXAMPLES:
            sage: G = graphs.PathGraph(3)
            sage: G.am()
            [0 1 0]
            [1 0 1]
            [0 1 0]

        Relabeling using a dictionary:
            sage: G.relabel({1:2,2:1}, inplace=False).am()
            [0 0 1]
            [0 0 1]
            [1 1 0]

        Relabeling using a list:
            sage: G.relabel([0,2,1], inplace=False).am()
            [0 0 1]
            [0 0 1]
            [1 1 0]

        Relabeling using a SAGE permutation:
            sage: from sage.groups.perm_gps.permgroup_named import SymmetricGroup
            sage: S = SymmetricGroup(3)
            sage: gamma = S('(1,2)')
            sage: G.relabel(gamma, inplace=False).am()
            [0 0 1]
            [0 0 1]
            [1 1 0]

        """
        if not inplace:
            G = self.copy()
            G.relabel(perm, quick, True)
            return G
        if type(perm) == list:
            if quick:
                n = self.order()
                numbr = 0
                if isinstance(self, Graph):
                    for i,j,l in self.edge_iterator():
                        numbr += 1<<((n-(perm[i]+1))*n + n-(perm[j]+1))
                        numbr += 1<<((n-(perm[j]+1))*n + n-(perm[i]+1))
                elif isinstance(self, DiGraph):
                    for i,j,l in self.edge_iterator():
                        numbr += 1<<((n-(perm[i]+1))*n + n-(perm[j]+1))
            if isinstance(self, Graph):
                oldd = self._nxg.adj
                newd = {}
                for v in oldd.iterkeys():
                    oldtempd = oldd[v]
                    newtempd = {}
                    for w in oldtempd.iterkeys():
                        newtempd[perm[w]] = oldtempd[w]
                    newd[perm[v]] = newtempd
                self._nxg.adj = newd
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
                self._nxg.adj = newsucc
                self._nxg.succ = self._nxg.adj
                self._nxg.pred = newpred

        from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
        if type(perm) == PermutationGroupElement:
            n = self.order()
            dict = {}
            llist = perm.list()
            for i in xrange(1,n):
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
                self._nxg.adj = newd
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
                self._nxg.adj = newsucc
                self._nxg.succ = self._nxg.adj
                self._nxg.pred = newpred

    ### Cliques

    def cliques(self):
        """
        Returns the list of maximal cliques.  Each maximal clique is
        represented by a list of vertices.

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.  (See examples below).

        Maximal cliques are the largest complete subgraphs containing a
        given point.  This function is based on Networkx's implementation
        of the Bron and Kerbosch Algorithm, [1].

        REFERENCE:
            [1] Coen Bron and Joep Kerbosch. (1973). Algorithm 457: Finding
                All Cliques of an Undirected Graph. Commun. ACM. v 16. n 9.
                pages 575-577. ACM Press. [Online] Available:
                http://www.ram.org/computing/rambin/rambin.html

        EXAMPLES:
            sage: (graphs.ChvatalGraph()).cliques()
            [[0, 1], [0, 4], [0, 6], [0, 9], [2, 1], [2, 3], [2, 6], [2, 8], [3, 4], [3, 7], [3, 9], [5, 1], [5, 4], [5, 10], [5, 11], [7, 1], [7, 8], [7, 11], [8, 4], [8, 10], [10, 6], [10, 9], [11, 6], [11, 9]]
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage.: D.show(figsize=[2,2])
            sage: D.cliques()
            Traceback (most recent call last):
            ...
            TypeError: Function defined for undirected graphs only.  See documentation.
            sage: D = D.to_undirected()
            sage.: D.show(figsize=[2,2])
            sage: D.cliques()
            [[0, 1, 2], [0, 1, 3]]
        """
        if (self.is_directed()):
            raise TypeError('Function defined for undirected graphs only.  See documentation.')
        else:
            import networkx.cliques
            return networkx.cliques.find_cliques(self._nxg)

    def cliques_get_max_clique_graph(self, name=''):
        """
        Returns a graph constructed with maximal cliques as vertices,
        and edges between maximal cliques with common members in
        the original graph.

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.  (See examples below).

        INPUT:
           name -- The name of the new graph.

        EXAMPLES:
            sage: (graphs.ChvatalGraph()).cliques_get_max_clique_graph()
            Graph on 24 vertices
            sage.: ((graphs.ChvatalGraph()).cliques_get_max_clique_graph()).show(figsize=[2,2], vertex_size=20, vertex_labels=False)
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage.: D.show(figsize=[2,2])
            sage: D.cliques_get_max_clique_graph()
            Traceback (most recent call last):
            ...
            TypeError: Function defined for undirected graphs only.  See documentation.
            sage: D = D.to_undirected()
            sage.: D.show(figsize=[2,2])
            sage: D.cliques_get_max_clique_graph()
            Graph on 2 vertices
            sage.: (D.cliques_get_max_clique_graph()).show(figsize=[2,2])
        """
        if (self.is_directed()):
            raise TypeError('Function defined for undirected graphs only.  See documentation.')
        else:
            import networkx.cliques
            return Graph(networkx.cliques.make_max_clique_graph(self._nxg, name=name, create_using=networkx.xgraph.XGraph()))

    # Add fpos (below) when Bipartite class is wrapped.
    def cliques_get_clique_bipartite(self, **kwds):
        """
        Returns a bipartite graph constructed such that cliques are the
        top vertices and the bottom vertices are retained from the given graph.
        Top and bottom vertices are connected if the bottom vertex belongs to
        the clique represented by a top vertex.

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.  (See examples below).

        EXAMPLES:
            sage: (graphs.ChvatalGraph()).cliques_get_clique_bipartite()
            Graph on 36 vertices
            sage.: ((graphs.ChvatalGraph()).cliques_get_clique_bipartite()).show(figsize=[2,2], vertex_size=20, vertex_labels=False)
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage.: D.show(figsize=[2,2])
            sage: D.cliques_get_clique_bipartite()
            Traceback (most recent call last):
            ...
            TypeError: Function defined for undirected graphs only.  See documentation.
            sage: D = D.to_undirected()
            sage.: D.show(figsize=[2,2])
            sage: D.cliques_get_clique_bipartite()
            Graph on 6 vertices
            sage.: (D.cliques_get_clique_bipartite()).show(figsize=[2,2])
        """
        if (self.is_directed()):
            raise TypeError('Function defined for undirected graphs only.  See documentation.')
        else:
            import networkx.cliques
            return Graph(networkx.cliques.make_clique_bipartite(self._nxg, **kwds))

    # TODO: Also implement project_down and project_up after Bipartite class.

    def clique_number(self, cliques=None):
        """
        Returns the size of the largest clique of the graph (clique number).

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.  (See examples below).

        INPUT:
            -- cliques - list of cliques (if already computed)

        EXAMPLES:
            sage: C = Graph('DJ{')
            sage: C.clique_number()
            4
            sage: E = C.cliques()
            sage: E
            [[4, 1, 2, 3], [4, 0]]
            sage: C.clique_number(cliques=E)
            4
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage.: D.show(figsize=[2,2])
            sage: D.clique_number()
            Traceback (most recent call last):
            ...
            TypeError: Function defined for undirected graphs only.  See documentation.
            sage: D = D.to_undirected()
            sage.: D.show(figsize=[2,2])
            sage: D.clique_number()
            3
        """
        if (self.is_directed()):
            raise TypeError('Function defined for undirected graphs only.  See documentation.')
        else:
            import networkx.cliques
            return networkx.cliques.graph_clique_number(self._nxg, cliques)

    def cliques_vertex_clique_number(self, vertices=None, with_labels=False, cliques=None):
        r"""
        Returns a list of sizes of the largest maximal cliques containing
        each vertex.  (Returns a single value if only one input vertex).

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.  (See examples below).

        INPUT:
            -- vertices - the vertices to inspect (default is entire graph)
            -- with_labels - (boolean) default False returns list as above
                             True returns a dictionary keyed by vertex labels
            -- cliques - list of cliques (if already computed)

        EXAMPLES:
            sage: C = Graph('DJ{')
            sage: C.cliques_vertex_clique_number()
            [2, 4, 4, 4, 4]
            sage: E = C.cliques()
            sage: E
            [[4, 1, 2, 3], [4, 0]]
            sage: C.cliques_vertex_clique_number(cliques=E)
            [2, 4, 4, 4, 4]
            sage: F = graphs.Grid2dGraph(2,3)
            sage: F.cliques_vertex_clique_number(with_labels=True)
            {(0, 1): 2, (1, 2): 2, (0, 0): 2, (1, 1): 2, (1, 0): 2, (0, 2): 2}
            sage: F.cliques_vertex_clique_number(vertices=[(0, 1), (1, 2)])
            [2, 2]
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage.: D.show(figsize=[2,2])
            sage: D.cliques_vertex_clique_number()
            Traceback (most recent call last):
            ...
            TypeError: Function defined for undirected graphs only.  See documentation.
            sage: D = D.to_undirected()
            sage.: D.show(figsize=[2,2])
            sage: D.cliques_vertex_clique_number()
            [3, 3, 3, 3]
        """
        if (self.is_directed()):
            raise TypeError('Function defined for undirected graphs only.  See documentation.')
        else:
            import networkx.cliques
            return networkx.cliques.node_clique_number(self._nxg, vertices, with_labels, cliques)

    def cliques_number_of(self, vertices=None, cliques=None, with_labels=False):
        """
        Returns a list of the number of maximal cliques containing
        each vertex.  (Returns a single value if only one input vertex).

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.  (See examples below).

        INPUT:
            -- vertices - the vertices to inspect (default is entire graph)
            -- with_labels - (boolean) default False returns list as above
                             True returns a dictionary keyed by vertex labels
            -- cliques - list of cliques (if already computed)

        EXAMPLES:
            sage: C = Graph('DJ{')
            sage: C.cliques_number_of()
            [1, 1, 1, 1, 2]
            sage: E = C.cliques()
            sage: E
            [[4, 1, 2, 3], [4, 0]]
            sage: C.cliques_number_of(cliques=E)
            [1, 1, 1, 1, 2]
            sage: F = graphs.Grid2dGraph(2,3)
            sage: F.cliques_number_of(with_labels=True)
            {(0, 1): 3, (1, 2): 2, (0, 0): 2, (1, 1): 3, (1, 0): 2, (0, 2): 2}
            sage: F.cliques_number_of(vertices=[(0, 1), (1, 2)])
            [3, 2]
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage.: D.show(figsize=[2,2])
            sage: D.cliques_number_of()
            Traceback (most recent call last):
            ...
            TypeError: Function defined for undirected graphs only.  See documentation.
            sage: D = D.to_undirected()
            sage.: D.show(figsize=[2,2])
            sage: D.cliques_number_of()
            [2, 2, 1, 1]
        """
        if (self.is_directed()):
            raise TypeError('Function defined for undirected graphs only.  See documentation.')
        else:
            import networkx.cliques
            return networkx.cliques.number_of_cliques(self._nxg, vertices, cliques, with_labels)

    def cliques_containing_vertex(self, vertices=None, cliques=None, with_labels=False):
        """
        Returns the cliques containing each vertex, represented as a list of
        lists.  (Returns a single list if only one input vertex).

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.  (See examples below).

        INPUT:
            -- vertices - the vertices to inspect (default is entire graph)
            -- with_labels - (boolean) default False returns list as above
                             True returns a dictionary keyed by vertex labels
            -- cliques - list of cliques (if already computed)

        EXAMPLES:
            sage: C = Graph('DJ{')
            sage: C.cliques_containing_vertex()
            [[[4, 0]], [[4, 1, 2, 3]], [[4, 1, 2, 3]], [[4, 1, 2, 3]], [[4, 1, 2, 3], [4, 0]]]
            sage: E = C.cliques()
            sage: E
            [[4, 1, 2, 3], [4, 0]]
            sage: C.cliques_containing_vertex(cliques=E)
            [[[4, 0]], [[4, 1, 2, 3]], [[4, 1, 2, 3]], [[4, 1, 2, 3]], [[4, 1, 2, 3], [4, 0]]]
            sage: F = graphs.Grid2dGraph(2,3)
            sage: F.cliques_containing_vertex(with_labels=True)
            {(0, 1): [[(0, 1), (0, 0)], [(0, 1), (0, 2)], [(0, 1), (1, 1)]], (1, 2): [[(1, 2), (0, 2)], [(1, 2), (1, 1)]], (0, 0): [[(0, 1), (0, 0)], [(1, 0), (0, 0)]], (1, 1): [[(0, 1), (1, 1)], [(1, 2), (1, 1)], [(1, 0), (1, 1)]], (1, 0): [[(1, 0), (0, 0)], [(1, 0), (1, 1)]], (0, 2): [[(0, 1), (0, 2)], [(1, 2), (0, 2)]]}
            sage: F.cliques_containing_vertex(vertices=[(0, 1), (1, 2)])
            [[[(0, 1), (0, 0)], [(0, 1), (0, 2)], [(0, 1), (1, 1)]], [[(1, 2), (0, 2)], [(1, 2), (1, 1)]]]
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage.: D.show(figsize=[2,2])
            sage: D.cliques_containing_vertex()
            Traceback (most recent call last):
            ...
            TypeError: Function defined for undirected graphs only.  See documentation.
            sage: D = D.to_undirected()
            sage.: D.show(figsize=[2,2])
            sage: D.cliques_containing_vertex()
            [[[0, 1, 2], [0, 1, 3]], [[0, 1, 2], [0, 1, 3]], [[0, 1, 2]], [[0, 1, 3]]]
        """
        if (self.is_directed()):
            raise TypeError('Function defined for undirected graphs only.  See documentation.')
        else:
            import networkx.cliques
            return networkx.cliques.cliques_containing_node(self._nxg, vertices, cliques, with_labels)

    ### Cluster

    def cluster_triangles(self, nbunch=None, with_labels=False):
        r"""
        Returns the number of triangles for nbunch of vertices as an
        ordered list.

        The clustering coefficient of a graph is the fraction of
        possible triangles that are triangles,
        c_i = triangles_i / (k_i*(k_i-1)/2)
        where k_i is the degree of vertex i, [1].  A coefficient for
        the whole graph is the average of the c_i.  Transitivity is
        the fraction of all possible triangles which are triangles,
        T = 3*triangles/triads, [1].

        INPUT:
            -- nbunch - The vertices to inspect.  If nbunch=None, returns
                data for all vertices in the graph
            -- with_labels - (boolean) default False returns list as above
                             True returns dict keyed by vertex labels.

        REFERENCE:
            [1] Aric Hagberg, Dan Schult and Pieter Swart. NetworkX
                documentation. [Online] Available:
                https://networkx.lanl.gov/reference/networkx/

        EXAMPLES:
            sage: (graphs.FruchtGraph()).cluster_triangles()
            [1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0]
            sage: (graphs.FruchtGraph()).cluster_triangles(with_labels=True)
            {0: 1, 1: 1, 2: 0, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 0, 9: 1, 10: 1, 11: 0}
            sage: (graphs.FruchtGraph()).cluster_triangles(nbunch=[0,1,2])
            [1, 1, 0]
        """
        import networkx
        return networkx.triangles(self._nxg, nbunch, with_labels)

    def clustering_average(self):
        r"""
        Returns the average clustering coefficient.

        The clustering coefficient of a graph is the fraction of
        possible triangles that are triangles,
        c_i = triangles_i / (k_i*(k_i-1)/2)
        where k_i is the degree of vertex i, [1].  A coefficient for
        the whole graph is the average of the c_i.  Transitivity is
        the fraction of all possible triangles which are triangles,
        T = 3*triangles/triads, [1].

        REFERENCE:
            [1] Aric Hagberg, Dan Schult and Pieter Swart. NetworkX
                documentation. [Online] Available:
                https://networkx.lanl.gov/reference/networkx/

        EXAMPLES:
            sage: (graphs.FruchtGraph()).clustering_average()
            0.25
        """
        import networkx
        return networkx.average_clustering(self._nxg)

    def clustering_coeff(self, nbunch=None, with_labels=False, weights=False):
        r"""
        Returns the clustering coefficient for each vertex in nbunch
        as an ordered list.

        The clustering coefficient of a graph is the fraction of
        possible triangles that are triangles,
        c_i = triangles_i / (k_i*(k_i-1)/2)
        where k_i is the degree of vertex i, [1].  A coefficient for
        the whole graph is the average of the c_i.  Transitivity is
        the fraction of all possible triangles which are triangles,
        T = 3*triangles/triads, [1].

        INPUT:
            -- nbunch - the vertices to inspect (default None returns
                        data on all vertices in graph)
            -- with_labels - (boolean) default False returns list as above
                             True returns dict keyed by vertex labels.
            -- weights - default is False.  If both with_labels and weights
                        are True, then returns a clustering coefficient dict
                        and a dict of weights based on degree.  Weights are
                        the fraction of connected triples in the graph that
                        include the keyed vertex.

        REFERENCE:
            [1] Aric Hagberg, Dan Schult and Pieter Swart. NetworkX
                documentation. [Online] Available:
                https://networkx.lanl.gov/reference/networkx/

        EXAMPLES:
            sage: (graphs.FruchtGraph()).clustering_coeff()
            [0.33333333333333331, 0.33333333333333331, 0.0, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.33333333333333331, 0.0, 0.33333333333333331, 0.33333333333333331, 0.0]
            sage: (graphs.FruchtGraph()).clustering_coeff(with_labels=True)
            {0: 0.33333333333333331, 1: 0.33333333333333331, 2: 0.0, 3: 0.33333333333333331, 4: 0.33333333333333331, 5: 0.33333333333333331, 6: 0.33333333333333331, 7: 0.33333333333333331, 8: 0.0, 9: 0.33333333333333331, 10: 0.33333333333333331, 11: 0.0}
            sage: (graphs.FruchtGraph()).clustering_coeff(with_labels=True,weights=True)
            ({0: 0.33333333333333331, 1: 0.33333333333333331, 2: 0.0, 3: 0.33333333333333331, 4: 0.33333333333333331, 5: 0.33333333333333331, 6: 0.33333333333333331, 7: 0.33333333333333331, 8: 0.0, 9: 0.33333333333333331, 10: 0.33333333333333331, 11: 0.0}, {0: 0.083333333333333329, 1: 0.083333333333333329, 2: 0.083333333333333329, 3: 0.083333333333333329, 4: 0.083333333333333329, 5: 0.083333333333333329, 6: 0.083333333333333329, 7: 0.083333333333333329, 8: 0.083333333333333329, 9: 0.083333333333333329, 10: 0.083333333333333329, 11: 0.083333333333333329})
            sage: (graphs.FruchtGraph()).clustering_coeff(nbunch=[0,1,2])
            [0.33333333333333331, 0.33333333333333331, 0.0]
            sage: (graphs.FruchtGraph()).clustering_coeff(nbunch=[0,1,2],with_labels=True,weights=True)
            ({0: 0.33333333333333331, 1: 0.33333333333333331, 2: 0.0}, {0: 0.083333333333333329, 1: 0.083333333333333329, 2: 0.083333333333333329})
        """
        import networkx
        return networkx.clustering(self._nxg, nbunch, with_labels, weights)

    def cluster_transitivity(self):
        r"""
        Returns the transitivity (fraction of transitive triangles)
        of the graph.

        The clustering coefficient of a graph is the fraction of
        possible triangles that are triangles,
        c_i = triangles_i / (k_i*(k_i-1)/2)
        where k_i is the degree of vertex i, [1].  A coefficient for
        the whole graph is the average of the c_i.  Transitivity is
        the fraction of all possible triangles which are triangles,
        T = 3*triangles/triads, [1].

        REFERENCE:
            [1] Aric Hagberg, Dan Schult and Pieter Swart. NetworkX
                documentation. [Online] Available:
                https://networkx.lanl.gov/reference/networkx/

        EXAMPLES:
            sage: (graphs.FruchtGraph()).cluster_transitivity()
            0.25
        """
        import networkx
        return networkx.transitivity(self._nxg)

    ### Cores

    def cores(self, with_labels=False):
        """
        Returns the core number for each vertex in an ordered list.

        'K-cores in graph theory were introduced by Seidman in 1983
        and by Bollobas in 1984 as a method of (destructively) simplifying
        graph topology to aid in analysis and visualization. They have been
        more recently defined as the following by Batagelj et al: given a
        graph G with vertices set V and edges set E, the k-core is computed
        by pruning all the vertices (with their respective edges) with degree
        less than k. That means that if a vertex u has degree d_u, and it has
        n neighbors with degree less than k, then the degree of u becomes d_u - n,
        and it will be also pruned if k > d_u - n.  This operation can be
        useful to filter or to study some properties of the graphs. For
        instance, when you compute the 2-core of graph G, you are cutting
        all the vertices which are in a tree part of graph. (A tree is a
        graph with no loops),' [1].

        INPUT:
            -- with_labels - default False returns list as described above.
                             True returns dict keyed by vertex labels.

        REFERENCE:
            [1] K-core. Wikipedia. (2007). [Online] Available:
                http://en.wikipedia.org/wiki/K-core
            [2] Boris Pittel, Joel Spencer and Nicholas Wormald. Sudden
                Emergence of a Giant k-Core in a Random Graph. (1996).
                J. Combinatorial Theory. Ser B 67. pages 111-151. [Online]
                Available: http://cs.nyu.edu/cs/faculty/spencer/papers/k-core.pdf

        EXAMPLES:
            sage: (graphs.FruchtGraph()).cores()
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
            sage: (graphs.FruchtGraph()).cores(with_labels=True)
            {0: 3, 1: 3, 2: 3, 3: 3, 4: 3, 5: 3, 6: 3, 7: 3, 8: 3, 9: 3, 10: 3, 11: 3}
        """
        import networkx.cores
        return networkx.cores.find_cores(self._nxg, with_labels)

    ### Distance

    def distance(self, u, v):
        """
        Returns the (directed) distance from u to v in the (di)graph, i.e. the
        length of the shortest path from u to v.

        EXAMPLES:
            sage: G = graphs.CycleGraph(9)
            sage: G.distance(0,1)
            1
            sage: G.distance(0,4)
            4
            sage: G.distance(0,5)
            4
            sage: G = Graph( {0:[], 1:[]} )
            sage: G.distance(0,1)
            +Infinity

        """
        return self.shortest_path_length(u, v)

    def eccentricity(self, v=None, dist_dict=None, with_labels=False):
        """
        Return the eccentricity of vertex (or vertices) v.

        The eccentricity of a vertex is the maximum distance to any other
        vertex.

        INPUT:
            v -- either a single vertex or a list of vertices. If it is not
        specified, then it is taken to be all vertices.
            dist_dict -- optional, a dict of dicts of distance.
            with_labels -- Whether to return a list or a dict.

        EXAMPLES:
            sage: G = graphs.KrackhardtKiteGraph()
            sage: G.eccentricity()
            [4, 4, 4, 4, 4, 3, 3, 2, 3, 4]
            sage: G.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: G.eccentricity(7)
            2
            sage: G.eccentricity([7,8,9])
            [3, 4, 2]
            sage: G.eccentricity([7,8,9], with_labels=True) == {8: 3, 9: 4, 7: 2}
            True
            sage: G = Graph( { 0 : [], 1 : [], 2 : [1] } )
            sage: G.eccentricity()
            [+Infinity, +Infinity, +Infinity]
            sage: G = Graph({0:[]})
            sage: G.eccentricity(with_labels=True)
            {0: 0}
            sage: G = Graph({0:[], 1:[]})
            sage: G.eccentricity(with_labels=True)
            {0: +Infinity, 1: +Infinity}

        """
        import networkx
        try:
            return networkx.eccentricity(self._nxg, v, dist_dict, with_labels)
        except networkx.NetworkXError:
            from sage.rings.infinity import Infinity
            e = {}
            if v is None:
                v = self.vertices()
            elif not isinstance(v, list):
                v = [v]
            for u in v:
                e[u] = Infinity
            if with_labels:
                return e
            else:
                if len(e)==1: return e.values()[0] # return single value
                return e.values()

    def radius(self):
        """
        Returns the radius of the (di)graph.

        The radius is defined to be the minimum eccentricity of any vertex,
        where the eccentricity is the maximum distance to any other vertex.

        EXAMPLES:
        The more symmetric a graph is, the smaller (diameter - radius) is.
            sage: G = graphs.BarbellGraph(9, 3)
            sage: G.radius()
            3
            sage: G.diameter()
            6

            sage: G = graphs.OctahedralGraph()
            sage: G.radius()
            2
            sage: G.diameter()
            2

        """
        return min(self.eccentricity())

    def center(self):
        """
        Returns the set of vertices in the center, i.e. whose eccentricity is
        equal to the radius of the (di)graph.

        In other words, the center is the set of vertices achieving the
        minimum eccentricity.

        EXAMPLES:
            sage: G = graphs.DiamondGraph()
            sage: G.center()
            [1, 2]
            sage: P = graphs.PetersenGraph()
            sage: P.subgraph(P.center()) == P
            True
            sage: S = graphs.StarGraph(19)
            sage: S.center()
            [0]
            sage: G = Graph()
            sage: G.center()
            []
            sage: G.add_vertex()
            sage: G.center()
            [0]

        """
        e = self.eccentricity(with_labels=True)
        try:
            r = min(e.values())
        except:
            return []
        return [v for v in e if e[v]==r]

    def diameter(self):
        """
        Returns the largest distance between any two vertices. Returns
        Infinity if the (di)graph is not connected.

        EXAMPLES:
            sage: G = graphs.PetersenGraph()
            sage: G.diameter()
            2
            sage: G = Graph( { 0 : [], 1 : [], 2 : [1] } )
            sage: G.diameter()
            +Infinity

        Although max( {} ) is usually defined as -Infinity, since the diameter
        will never be negative, we define it to be zero:
            sage: G = graphs.EmptyGraph()
            sage: G.diameter()
            0

        """
        e = self.eccentricity()
        if not isinstance(e, list):
            e = [e]
        if len(e) == 0:
            return 0
        return max(e)

    def periphery(self):
        """
        Returns the set of vertices in the periphery, i.e. whose eccentricity
        is equal to the diameter of the (di)graph.

        In other words, the periphery is the set of vertices achieving the
        maximum eccentricity.

        EXAMPLES:
            sage: G = graphs.DiamondGraph()
            sage: G.periphery()
            [0, 3]
            sage: P = graphs.PetersenGraph()
            sage: P.subgraph(P.periphery()) == P
            True
            sage: S = graphs.StarGraph(19)
            sage: S.periphery()
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: G = Graph()
            sage: G.periphery()
            []
            sage: G.add_vertex()
            sage: G.periphery()
            [0]

        """
        e = self.eccentricity(with_labels=True)
        try:
            r = max(e.values())
        except:
            return []
        return [v for v in e if e[v]==r]

    ### Paths

    def shortest_path(self, u, v, by_weight=False, bidirectional=True):
        """
        Returns a list of vertices representing some shortest path from u to
        v: if there is no path from u to v, the list is empty.

        INPUT:
            by_weight -- if False, uses a breadth first search. If True, takes
        edge weightings into account, using Dijkstra's algorithm.
            bidirectional -- if True, the algorithm will expand vertices from
        u and v at the same time, making two spheres of half the usual radius.
        This generally doubles the speed (consider the total volume in each
        case).

        EXAMPLE:
            sage: D = graphs.DodecahedralGraph()
            sage: D.shortest_path(4, 9)
            [4, 17, 16, 12, 13, 9]
            sage: D.shortest_path(5, 5)
            [5]
            sage: D.delete_vertices([9,12,14])
            sage: D.shortest_path(13, 4)
            []
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} } )
            sage: G.plot(edge_labels=True).save('sage.png')
            sage: G.shortest_path(0, 3)
            [0, 4, 3]
            sage: G.shortest_path(0, 3, by_weight=True)
            [0, 1, 2, 3]

        """ #         TODO- multiple edges??
        if u == v: # to avoid a NetworkX bug
            return [u]
        import networkx
        if by_weight:
            if bidirectional:
                try:
                    L = networkx.bidirectional_dijkstra(self._nxg, u, v)[1]
                except:
                    L = False
            else:
                L = networkx.dijkstra_path(self._nxg, u, v)
        else:
            if bidirectional:
                L = networkx.shortest_path(self._nxg, u, v)
            else:
                try:
                    L = networkx.single_source_shortest_path(self._nxg, u)[v]
                except:
                    L = False
        if L:
            return L
        else:
            return []

    def shortest_path_length(self, u, v, by_weight=False,
                                         bidirectional=True,
                                         weight_sum=None):
        """
        Returns the minimal length of paths from u to v: if there is no path
        from u to v, returns Infinity.

        INPUT:
            by_weight -- if False, uses a breadth first search. If True, takes
        edge weightings into account, using Dijkstra's algorithm.
            bidirectional -- if True, the algorithm will expand vertices from
        u and v at the same time, making two spheres of half the usual radius.
        This generally doubles the speed (consider the total volume in each
        case).
            weight_sum -- if False, returns the number of edges in the path.
        If True, returns the sum of the weights of these edges. Default
        behavior is to have the same value as by_weight.

        EXAMPLE:
            sage: D = graphs.DodecahedralGraph()
            sage: D.shortest_path_length(4, 9)
            5
            sage: D.shortest_path_length(5, 5)
            0
            sage: D.delete_vertices([9,12,14])
            sage: D.shortest_path_length(13, 4)
            +Infinity
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} } )
            sage: G.plot(edge_labels=True).save('sage.png')
            sage: G.shortest_path_length(0, 3)
            2
            sage: G.shortest_path_length(0, 3, by_weight=True)
            3

        """
        if weight_sum is None:
            weight_sum = by_weight
        path = self.shortest_path(u, v, by_weight, bidirectional)
        length = len(path) - 1
        if length == -1:
            from sage.rings.infinity import Infinity
            return Infinity
        if weight_sum:
            wt = 0
            for j in range(length):
                wt += self.edge_label(path[j], path[j+1])
            return wt
        else:
            return length

    def shortest_paths(self, u, by_weight=False, cutoff=None):
        """
        Returns a dictionary d of shortest paths d[v] from u to v, for each
        vertex v connected by a path from u.

        INPUT:
            by_weight -- if False, uses a breadth first search. If True, uses
        Dijkstra's algorithm to find the shortest paths by weight.
            cutoff -- integer depth to stop search. Ignored if by_weight is
        True.

        EXAMPLES:
            sage: D = graphs.DodecahedralGraph()
            sage: D.shortest_paths(0)
            {0: [0], 1: [0, 1], 2: [0, 1, 2], 3: [0, 19, 3], 4: [0, 19, 3, 4], 5: [0, 19, 3, 4, 5], 6: [0, 1, 2, 6], 7: [0, 1, 8, 7], 8: [0, 1, 8], 9: [0, 10, 9], 10: [0, 10], 11: [0, 10, 11], 12: [0, 10, 11, 12], 13: [0, 10, 9, 13], 14: [0, 1, 8, 7, 14], 15: [0, 10, 11, 12, 16, 15], 16: [0, 10, 11, 12, 16], 17: [0, 19, 18, 17], 18: [0, 19, 18], 19: [0, 19]}
            sage: D.shortest_paths(0, cutoff=2)
            {0: [0], 1: [0, 1], 2: [0, 1, 2], 3: [0, 19, 3], 8: [0, 1, 8], 9: [0, 10, 9], 10: [0, 10], 11: [0, 10, 11], 18: [0, 19, 18], 19: [0, 19]}
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} } )
            sage: G.plot(edge_labels=True).save('sage.png')
            sage: G.shortest_paths(0, by_weight=True)
            {0: [0], 1: [0, 1], 2: [0, 1, 2], 3: [0, 1, 2, 3], 4: [0, 4]}

        """
        import networkx
        if by_weight:
            return networkx.single_source_dijkstra_path(self._nxg, u)
        else:
            return networkx.single_source_shortest_path(self._nxg, u, cutoff)

    def shortest_path_lengths(self, u, by_weight=False, weight_sums=None):
        """
        Returns a dictionary of shortest path lengths keyed by targets that
        are connected by a path from u.

        INPUT:
            by_weight -- if False, uses a breadth first search. If True, takes
        edge weightings into account, using Dijkstra's algorithm.
            weight_sums -- if False, returns the number of edges in each path.
        If True, returns the sum of the weights of these edges. Default
        behavior is to have the same value as by_weight.

        EXAMPLES:
            sage: D = graphs.DodecahedralGraph()
            sage: D.shortest_path_lengths(0)
            {0: 0, 1: 1, 2: 2, 3: 2, 4: 3, 5: 4, 6: 3, 7: 3, 8: 2, 9: 2, 10: 1, 11: 2, 12: 3, 13: 3, 14: 4, 15: 5, 16: 4, 17: 3, 18: 2, 19: 1}
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} } )
            sage: G.plot(edge_labels=True).save('sage.png')
            sage: G.shortest_path_lengths(0, by_weight=True)
            {0: 0, 1: 1, 2: 2, 3: 3, 4: 2}

        """
        if weight_sums is None:
            weight_sums = by_weight
        paths = self.shortest_paths(u, by_weight)
        if weight_sums:
            weights = {}
            for v in paths:
                wt = 0
                path = paths[v]
                for j in range(len(path) - 1):
                    wt += self.edge_label(path[j], path[j+1])
                weights[v] = wt
            return weights
        else:
            lengths = {}
            for v in paths:
                lengths[v] = len(paths[v]) - 1
            return lengths

    def shortest_path_all_pairs(self):
        """
        Uses the Floyd-Warshall algorithm to find a shortest path for each
        pair of vertices.

        OUTPUT:
            A tuple (dist, pred). They are both dicts of dicts. The first
        indicates the length dist[u][v] of the shortest weighted path from u
        to v. The second is more complicated-- it indicates the predecessor
        pred[u][v] of v in the shortest path from u to v.

        EXAMPLE:
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} } )
            sage: G.plot(edge_labels=True).save('sage.png')
            sage: dist, pred = G.shortest_path_all_pairs()
            sage: dist
            {0: {0: 0, 1: 1, 2: 2, 3: 3, 4: 2}, 1: {0: 1, 1: 0, 2: 1, 3: 2, 4: 3}, 2: {0: 2, 1: 1, 2: 0, 3: 1, 4: 3}, 3: {0: 3, 1: 2, 2: 1, 3: 0, 4: 2}, 4: {0: 2, 1: 3, 2: 3, 3: 2, 4: 0}}
            sage: pred
            {0: {0: None, 1: 0, 2: 1, 3: 2, 4: 0}, 1: {0: 1, 1: None, 2: 1, 3: 2, 4: 0}, 2: {0: 1, 1: 2, 2: None, 3: 2, 4: 3}, 3: {0: 1, 1: 2, 2: 3, 3: None, 4: 3}, 4: {0: 4, 1: 0, 2: 3, 3: 4, 4: None}}
            sage: pred[0]
            {0: None, 1: 0, 2: 1, 3: 2, 4: 0}

        So for example the shortest weighted path from 0 to 3 is obtained as
        follows. The predecessor of 3 is pred[0][3] == 2, the predecessor of 2
        is pred[0][2] == 1, and the predecessor of 1 is pred[0][1] == 0.
        """
        from sage.rings.infinity import Infinity
        dist = {}
        pred = {}
        adj = self._nxg.adj
        verts = self.vertices()
        for u in verts:
            dist[u] = {}
            pred[u] = {}
            for v in verts:
                if adj[u].has_key(v):
                    dist[u][v] = adj[u][v]
                    pred[u][v] = u
                else:
                    dist[u][v] = Infinity
                    pred[u][v] = None
            dist[u][u] = 0

        for w in verts:
            for u in verts:
                for v in verts:
                    if dist[u][v] > dist[u][w] + dist[w][v]:
                        dist[u][v] = dist[u][w] + dist[w][v]
                        pred[u][v] = pred[w][v]

        return dist, pred

    ### Searches

    def breadth_first_search(self, u):
        """
        Returns an iterator over vertices in a breadth-first ordering.

        EXAMPLES:
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} } )
            sage: list(G.breadth_first_search(0))
            [0, 1, 4, 2, 3]
            sage: list(G.depth_first_search(0))
            [0, 4, 3, 2, 1]
            sage: D = DiGraph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} } )
            sage: list(D.breadth_first_search(0))
            [0, 1, 2, 3, 4]
            sage: list(D.depth_first_search(0))
            [0, 1, 2, 3, 4]

        """
        # This function is straight from an old version of networkx
        if self.is_directed():
            neighbors=self.successor_iterator
        else:
            neighbors=self.neighbor_iterator
        # nlist=[u] # list of nodes in a BFS order
        yield u
        seen={} # nodes seen
        queue=[u] # FIFO queue
        seen[u]=True
        while queue:
            v=queue.pop(0)  # this is expensive, should use a faster FIFO queue
            for w in neighbors(v):
                if w not in seen:
                    seen[w]=True
                    queue.append(w)
                    # nlist.append(w)
                    yield w
        # return nlist

    def depth_first_search(self, u):
        """
        Returns an iterator over vertices in a depth-first ordering.

        EXAMPLES:
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} } )
            sage: list(G.breadth_first_search(0))
            [0, 1, 4, 2, 3]
            sage: list(G.depth_first_search(0))
            [0, 4, 3, 2, 1]
            sage: D = DiGraph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} } )
            sage: list(D.breadth_first_search(0))
            [0, 1, 2, 3, 4]
            sage: list(D.depth_first_search(0))
            [0, 1, 2, 3, 4]

        """
        # This function is straight from an old version of networkx
        if self.is_directed():
            neighbors=self.successor_iterator
        else:
            neighbors=self.neighbor_iterator
        # nlist=[] # list of nodes in a DFS preorder
        seen={} # nodes seen
        queue=[u]  # use as LIFO queue
        seen[u]=True
        while queue:
            v=queue.pop()
            # nlist.append(v)
            yield v
            for w in neighbors(v):
                if w not in seen:
                    seen[w]=True
                    queue.append(w)
        # return nlist

    ### Constructors

    def am(self):
        """
        Shorter call for adjacency matrix makes life easier.

        """
        return self.adjacency_matrix()

    def complement(self):
        """
        Returns the complement of the (di)graph.

        The complement of a graph has the same vertices, but exactly those
        edges that are not in the original graph. This is not well defined for
        graphs with loops or multiple edges.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: P.plot().save('sage.png')
            sage: PC = P.complement()
            sage: PC.plot().save('sage.png')

        """
        if self.loops():
            raise TypeError('(Di)Graph complement not well defined for (di)graphs with loops.')
        if self.is_directed():
            if self.multiple_edges():
                raise TypeError('Digraph complement not well defined for graphs with multiple edges.')
            import networkx
            D = DiGraph(networkx.complement(self._nxg), pos=self._pos)
            return D
        else:
            if self.multiple_edges():
                raise TypeError('Graph complement not well defined for graphs with multiple edges.')
            import networkx
            G = Graph(networkx.complement(self._nxg), pos=self._pos)
            return G


    def line_graph(self):
        """
        Returns the line graph of the (di)graph.

        The line graph of an undirected graph G is an undirected graph
        H such that the vertices of H are the edges of G and two
        vertices e and f of H are adjacent if e and f share a common
        vertex in G.  In other words, an edge in H represents a path
        of length 2 in G.

        The line graph of a directed graph G is a directed graph H
        such that the vertices of H are the edges of G and two
        vertices e and f of H are adjacent if e and f share a common
        vertex in G and the terminal vertex of e is the initial vertex
        of f.  In other words, an edge in H represents a (directed)
        path of length 2 in G.


        EXAMPLE:
            sage: g=graphs.CompleteGraph(4)
            sage: h=g.line_graph()
            sage: h.vertices()
            [(0, 1, None),
            (0, 2, None),
            (0, 3, None),
            (1, 2, None),
            (1, 3, None),
            (2, 3, None)]
            sage: h.am()
            [0 1 1 1 1 0]
            [1 0 1 1 0 1]
            [1 1 0 0 1 1]
            [1 1 0 0 1 1]
            [1 0 1 1 0 1]
            [0 1 1 1 1 0]
            sage: g = DiGraph([[1..4],lambda i,j: i<j])
            sage: h = g.line_graph()
            sage: h.vertices()
            [(1, 2, None),
            (1, 3, None),
            (1, 4, None),
            (2, 3, None),
            (2, 4, None),
            (3, 4, None)]
            sage: h.edges()
            [((1, 2, None), (2, 3, None), None),
             ((1, 2, None), (2, 4, None), None),
             ((1, 3, None), (3, 4, None), None),
             ((2, 3, None), (3, 4, None), None)]
        """
        if self.is_directed():
            G=DiGraph()
            G.add_vertices(self.edges())
            for v in self:
                # Connect appropriate incident edges of the vertex v
                G.add_edges([(e,f) for e in self.incoming_edge_iterator(v) \
                             for f in self.edge_iterator(v)])
            return G
        else:
            G=Graph()
            # We must sort the edges' endpoints so that (1,2,None) is
            # seen as the same edge as (2,1,None).
            elist=[(min(i[0:2]),max(i[0:2]),i[2])
                   for i in self.edge_iterator()]
            G.add_vertices(elist)
            for v in self:
                elist=[(min(i[0:2]),max(i[0:2]),i[2])
                       for i in self.edge_iterator(v)]
                G.add_edges([(e, f) for e in elist for f in elist])
            return G



    def disjoint_union(self, other):
        """
        Returns the disjoint union of self and other.

        If there are common vertices to both, they will be renamed.

        EXAMPLE:
            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: D.disjoint_union(P)
            union( Dodecahedron, Petersen graph ): Graph on 30 vertices

        """
        if (self.is_directed() and not other.is_directed()) or (not self.is_directed() and other.is_directed()):
            raise TypeError('Both arguments must be of the same class.')
        repeat = False
        for u in self.vertices():
            for v in other.vertices():
                if u == v:
                    repeat = True
                    break
            if repeat: break
        if repeat:
            rename = ('0,','1,')
        else:
            rename = False
        import networkx
        if self.is_directed():
            return DiGraph(networkx.union(self._nxg, other._nxg, rename=rename))
        else:
            return Graph(networkx.union(self._nxg, other._nxg, rename=rename))

    def union(self, other):
        """
        Returns the union of self and other.

        If there are common vertices to both, they will be renamed.

        EXAMPLE:
            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: D.union(P)
            Graph on 20 vertices

        """
        if (self.is_directed() and not other.is_directed()) or (not self.is_directed() and other.is_directed()):
            raise TypeError('Both arguments must be of the same class.')
        if self.is_directed():
            G = DiGraph()
        else:
            G = Graph()
        G.add_vertices(self.vertices())
        G.add_vertices(other.vertices())
        G.add_edges(self.edges())
        G.add_edges(other.edges())
        return G

    def cartesian_product(self, other):
        """
        Returns the Cartesian product of self and other.

        The Cartesian product of G and H is the graph L with vertex set
        V(L) equal to the Cartesian product of the vertices V(G) and V(H), and
        ((u,v), (w,x)) is an edge iff either
            - (u, w) is an edge of self and v = x, or
            - (v, x) is an edge of other and u = w.

        EXAMPLES:
            sage: Z = graphs.CompleteGraph(2)
            sage: C = graphs.CycleGraph(5)
            sage: P = C.cartesian_product(Z); P
            Graph on 10 vertices
            sage: P.plot().save('sage.png')

            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: C = D.cartesian_product(P); C
            Graph on 200 vertices
            sage: C.plot().save('sage.png')

        """
        if (self.is_directed() and not other.is_directed()) or (not self.is_directed() and other.is_directed()):
            raise TypeError('Both arguments must be of the same class.')
        if self.is_directed():
            G = DiGraph()
        else:
            G = Graph()
        verts = []
        for a in self.vertices():
            for b in other.vertices():
                G.add_vertex((a,b))
                verts.append((a,b))
        for i in range(len(verts)):
            for j in range(i):
                u,v = verts[i]
                w,x = verts[j]
                if (self.has_edge(u, w) and v == x) or (other.has_edge(v, x) and u == w):
                    G.add_edge((u,v), (w,x))
        return G

    def tensor_product(self, other):
        """
        Returns the tensor product of self and other.

        The tensor product of G and H is the graph L with vertex set
        V(L) equal to the Cartesian product of the vertices V(G) and V(H), and
        ((u,v), (w,x)) is an edge iff
            - (u, w) is an edge of self, and
            - (v, x) is an edge of other.

        EXAMPLES:
            sage: Z = graphs.CompleteGraph(2)
            sage: C = graphs.CycleGraph(5)
            sage: T = C.tensor_product(Z); T
            Graph on 10 vertices
            sage: T.plot().save('sage.png')

            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: T = D.tensor_product(P); T
            Graph on 200 vertices
            sage: T.plot().save('sage.png')

        """
        if (self.is_directed() and not other.is_directed()) or (not self.is_directed() and other.is_directed()):
            raise TypeError('Both arguments must be of the same class.')
        if self.is_directed():
            G = DiGraph()
        else:
            G = Graph()
        verts = []
        for a in self.vertices():
            for b in other.vertices():
                G.add_vertex((a,b))
                verts.append((a,b))
        for i in range(len(verts)):
            for j in range(i):
                u,v = verts[i]
                w,x = verts[j]
                if self.has_edge(u, w) and other.has_edge(v, x):
                    G.add_edge((u,v), (w,x))
        return G

    def lexicographic_product(self, other):
        """
        Returns the lexicographic product of self and other.

        The lexicographic product of G and H is the graph L with vertex set
        V(L) equal to the Cartesian product of the vertices V(G) and V(H), and
        ((u,v), (w,x)) is an edge iff
            - (u, w) is an edge of self, or
            - u = w and (v, x) is an edge of other.

        EXAMPLES:
            sage: Z = graphs.CompleteGraph(2)
            sage: C = graphs.CycleGraph(5)
            sage: L = C.lexicographic_product(Z); L
            Graph on 10 vertices
            sage: L.plot().save('sage.png')

            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: L = D.lexicographic_product(P); L
            Graph on 200 vertices
            sage: L.plot().save('sage.png')

        """
        if (self.is_directed() and not other.is_directed()) or (not self.is_directed() and other.is_directed()):
            raise TypeError('Both arguments must be of the same class.')
        if self.is_directed():
            G = DiGraph()
        else:
            G = Graph()
        verts = []
        for a in self.vertices():
            for b in other.vertices():
                G.add_vertex((a,b))
                verts.append((a,b))
        for i in range(len(verts)):
            for j in range(i):
                u,v = verts[i]
                w,x = verts[j]
                if self.has_edge(u, w) or (u == w and other.has_edge(v, x)):
                    G.add_edge((u,v), (w,x))
        return G

    def strong_product(self, other):
        """
        Returns the strong product of self and other.

        The strong product of G and H is the graph L with vertex set
        V(L) equal to the Cartesian product of the vertices V(G) and V(H), and
        ((u,v), (w,x)) is an edge iff either
            - (u, w) is an edge of self and v = x, or
            - (v, x) is an edge of other and u = w, or
            - (u, w) is an edge of self and (v, x) is an edge of
        other.


        EXAMPLES:
            sage: Z = graphs.CompleteGraph(2)
            sage: C = graphs.CycleGraph(5)
            sage: S = C.strong_product(Z); S
            Graph on 10 vertices
            sage: S.plot().save('sage.png')

            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: S = D.strong_product(P); S
            Graph on 200 vertices
            sage: S.plot().save('sage.png')

        """
        if (self.is_directed() and not other.is_directed()) or (not self.is_directed() and other.is_directed()):
            raise TypeError('Both arguments must be of the same class.')
        if self.is_directed():
            G = DiGraph()
        else:
            G = Graph()

        verts = []
        for a in self.vertices():
            for b in other.vertices():
                G.add_vertex((a,b))
                verts.append((a,b))
        for i in range(len(verts)):
            for j in range(i):
                u,v = verts[i]
                w,x = verts[j]
                if (self.has_edge(u, w) and v == x) or \
                   (other.has_edge(v, x) and u == w) or \
                   (self.has_edge(u, w) and other.has_edge(v, x)):
                    G.add_edge((u,v), (w,x))
        return G

    def disjunctive_product(self, other):
        """
        Returns the disjunctive product of self and other.

        The disjunctive product of G and H is the graph L with vertex set
        V(L) equal to the Cartesian product of the vertices V(G) and V(H), and
        ((u,v), (w,x)) is an edge iff either
            - (u, w) is an edge of self, or
            - (v, x) is an edge of other.

        EXAMPLES:
            sage: Z = graphs.CompleteGraph(2)
            sage: D = Z.disjunctive_product(Z); D
            Graph on 4 vertices
            sage: D.plot().save('sage.png')

            sage: C = graphs.CycleGraph(5)
            sage: D = C.disjunctive_product(Z); D
            Graph on 10 vertices
            sage: D.plot().save('sage.png')

        """
        if (self.is_directed() and not other.is_directed()) or (not self.is_directed() and other.is_directed()):
            raise TypeError('Both arguments must be of the same class.')
        if self.is_directed():
            G = DiGraph()
        else:
            G = Graph()
        verts = []
        for a in self.vertices():
            for b in other.vertices():
                G.add_vertex((a,b))
                verts.append((a,b))
        for i in range(len(verts)):
            for j in range(i):
                u,v = verts[i]
                w,x = verts[j]
                if self.has_edge(u, w) or other.has_edge(v, x):
                    G.add_edge((u,v), (w,x))
        return G

    ### Visualization
    def _color_by_label(self, format='hex'):
        """
        Logic for coloring by label (factored out from plot() for use
        in 3d plots, etc)
        """
        from sage.plot.plot import rainbow
        edge_labels = []
        if self.is_directed():
            iterator = self.edge_iterator
        else:
            iterator = self.edge_iterator
        for e in iterator():
            i = 0
            while i < len(edge_labels):
                if not edge_labels[i][0][2] == e[2]:
                    i += 1
                else:
                    edge_labels[i].append(e)
                    break
            if i == len(edge_labels):
                edge_labels.append([e])
        num_labels = len(edge_labels)
        r = rainbow(num_labels, format=format)
        edge_colors = {}
        for i in range(num_labels):
            edge_colors[r[i]] = edge_labels[i]
        return edge_colors

    def plot(self, pos=None, layout=None, vertex_labels=True,
             edge_labels=False, vertex_size=200, graph_border=False,
             vertex_colors=None, partition=None, edge_colors=None,
             scaling_term=0.05, iterations=50,
             color_by_label=False, heights=None):
        """
        Returns a graphics object representing the (di)graph.

        INPUT:
            pos -- an optional positioning dictionary
            layout -- what kind of layout to use, takes precedence over pos
                'circular' -- plots the graph with vertices evenly distributed on a circle
                'spring' -- uses the traditional spring layout, ignores the graphs current positions
            vertex_labels -- whether to print vertex labels
            edge_labels -- whether to print edge labels. By default, False, but if True, the result
                of str(l) is printed on the edge for each label l. Labels equal to None are not printed.
            vertex_size -- size of vertices displayed
            graph_border -- whether to include a box around the graph
            vertex_colors -- optional dictionary to specify vertex colors: each key is a color recognizable
                by matplotlib, and each corresponding entry is a list of vertices. If a vertex is not listed,
                it looks invisible on the resulting plot (it doesn't get drawn).
            edge_colors -- a dictionary specifying edge colors: each key is a color recognized by
                matplotlib, and each entry is a list of edges.
            partition -- a partition of the vertex set. if specified, plot will show each cell in a different
                color. vertex_colors takes precedence.
            scaling_term -- default is 0.05. if vertices are getting chopped off, increase; if graph
                is too small, decrease. should be positive, but values much bigger than
                1/8 won't be useful unless the vertices are huge
            iterations -- how many iterations of the spring layout algorithm to
                go through, if applicable
            color_by_label -- if True, color edges by their labels
            heights -- if specified, this is a dictionary from a set of
                floating point heights to a set of vertices

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
            sage: pl = P.plot(pos=pos_dict, vertex_colors=d)
            sage: pl.save('sage.png')

            sage: C = graphs.CubeGraph(8)
            sage: P = C.plot(vertex_labels=False, vertex_size=0, graph_border=True)
            sage: P.save('sage.png')

            sage: G = graphs.HeawoodGraph()
            sage: for u,v,l in G.edges():
            ...    G.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: G.plot(edge_labels=True).save('sage.png')

            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []} )
            sage: for u,v,l in D.edges():
            ...    D.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
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
            sage: C.plot(vertex_labels=False, vertex_size=0, edge_colors=edge_colors).save('sage.png')

        """
        from sage.plot.plot import networkx_plot
        import networkx
        if vertex_colors is None:
            if partition is not None:
                l = len(partition)
                R = rainbow(l)
                vertex_colors = {}
                for i in range(l):
                    vertex_colors[R[i]] = partition[i]
            elif len(self._boundary) != 0:
                vertex_colors = {}
                bdy_verts = []
                int_verts = []
                for v in self.vertex_iterator():
                    if v in self._boundary:
                        bdy_verts.append(v)
                    else:
                        int_verts.append(v)
                vertex_colors['#ffffff'] = bdy_verts
                vertex_colors['#999999'] = int_verts
        if pos is None and layout is None and heights is None:
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
        elif heights is not None:
            pos = {}
            mmax = max([len(ccc) for ccc in heights.values()])
            dist = (1.0/(mmax+1))
            for height in heights:
                num_xes = len(heights[height])
                if num_xes == 0: continue
                j = (mmax - num_xes)/2.0
                for k in range(num_xes):
                    pos[heights[height][k]] = [ dist * (j+k+1), height ]
        if pos is None:
            pos = graph_fast.spring_layout_fast(self, iterations=iterations)
        else:
            for v in pos:
                for a in range(len(pos[v])):
                    pos[v][a] = float(pos[v][a])

        if color_by_label:
            edge_colors = self._color_by_label()

        G = networkx_plot(self._nxg, pos=pos, vertex_labels=vertex_labels, vertex_size=vertex_size, vertex_colors=vertex_colors, edge_colors=edge_colors, graph_border=graph_border, scaling_term=scaling_term)
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

    def show(self, pos=None, layout=None, vertex_labels=True,
             edge_labels=False, vertex_size=200, graph_border=False,
             vertex_colors=None, edge_colors=None, partition=None,
             scaling_term=0.05, talk=False, iterations=50,
             color_by_label=False, heights=None, **kwds):
        """
        Shows the (di)graph.

        INPUT:
            pos -- an optional positioning dictionary
            layout -- what kind of layout to use, takes precedence over pos
                'circular' -- plots the graph with vertices evenly distributed on a circle
                'spring' -- uses the traditional spring layout, ignores the graphs current positions
            vertex_labels -- whether to print vertex labels
            edge_labels -- whether to print edgeedge labels. By default, False, but if True, the result
                of str(l) is printed on the edge for each label l. Labels equal to None are not printed.
            vertex_size -- size of vertices displayed
            graph_border -- whether to include a box around the graph
            vertex_colors -- optional dictionary to specify vertex colors: each key is a color recognizable
                by matplotlib, and each corresponding entry is a list of vertices. If a vertex is not listed,
                it looks invisible on the resulting plot (it doesn't get drawn).
            edge_colors -- a dictionary specifying edge colors: each key is a color recognized by
                matplotlib, and each entry is a list of edges.
            partition -- a partition of the vertex set. if specified, plot will show each cell in a different
                color. vertex_colors takes precedence.
            scaling_term -- default is 0.05. if vertices are getting chopped off, increase; if graph
                is too small, decrease. should be positive, but values much bigger than
                1/8 won't be useful unless the vertices are huge
            talk -- if true, prints large vertices with white backgrounds so that labels are legible on slies
            iterations -- how many iterations of the spring layout algorithm to
                go through, if applicable
            color_by_label -- if True, color edges by their labels
            heights -- if specified, this is a dictionary from a set of
                floating point heights to a set of vertices

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
            sage: pl = P.plot(pos=pos_dict, vertex_colors=d)
            sage: pl.save('sage.png')

            sage: C = graphs.CubeGraph(8)
            sage: P = C.plot(vertex_labels=False, vertex_size=0, graph_border=True)
            sage: P.save('sage.png')

            sage: G = graphs.HeawoodGraph()
            sage: for u,v,l in G.edges():
            ...    G.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: G.plot(edge_labels=True).save('sage.png')

            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []} )
            sage: for u,v,l in D.edges():
            ...    D.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
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
            sage: C.plot(vertex_labels=False, vertex_size=0, edge_colors=edge_colors).save('sage.png')

        """
        if talk:
            vertex_size = 500
            if partition is None:
                vertex_colors = {'#FFFFFF':self.vertices()}
        self.plot(pos=pos, layout=layout, vertex_labels=vertex_labels,
                  edge_labels=edge_labels, vertex_size=vertex_size,
                  vertex_colors=vertex_colors, edge_colors=edge_colors,
                  graph_border=graph_border, partition=partition,
                  scaling_term=scaling_term, iterations=iterations,
                  color_by_label=color_by_label,
                  heights=heights).show(**kwds)

    def transitive_closure(self):
        r"""
        Modifies a graph to be its transitive closure and returns the
        modified graph.

        The transitive closure of a graph G has an edge (x,y) if and
        only if there is a path between x and y in G.

        The transitive closure of any strongly connected component of
        a graph is a complete graph.  In particular, the transitive
        closure of a connected undirected graph is a complete graph.
        The transitive closure of a directed acyclic graph is a
        directed acyclic graph representing the full partial order.

        EXAMPLES:
            sage: g=graphs.PathGraph(4)
            sage: g.transitive_closure()==graphs.CompleteGraph(4)
            True
            sage: g=DiGraph({0:[1,2], 1:[3], 2:[5,6]})
            sage: g.transitive_closure().edges(labels=False)
            [(0, 1), (0, 2), (0, 3), (0, 5), (0, 6), (1, 3), (2, 5), (2, 6)]

        """
        G = self.copy()
        for v in G:
            # todo optimization opportunity: we are adding edges that
            # are already in the graph and we are adding edges
            # one at a time.
            for e in G.breadth_first_search(v):
                G.add_edge((v,e))
        return G

    def antisymmetric(self):
        r"""
        Returns True if the relation given by the graph is
        antisymmetric and False otherwise.

        A graph represents an antisymmetric relation if there being a
        path from a vertex x to a vertex y implies that there is not a
        path from y to x unless x=y.

        A directed acyclic graph is antisymmetric.  An undirected
        graph is never antisymmetric unless it is just a union of
        isolated vertices.

        sage: graphs.RandomGNP(20,0.5).antisymmetric()
        False
        sage: graphs.RandomDirectedGNR(20,0.5).antisymmetric()
        True

        """
        if not self.is_directed():
            if self.size()-len(self.loop_edges())>0:
                return False
            else:
                return True

        g = self.copy()
        g.multiple_edges(False)
        g.loops(False)
        g = g.transitive_closure()
        gpaths = g.edges(labels=False)
        for e in gpaths:
            if (e[1],e[0]) in gpaths:
                return False
        return True

    def independent_set(self, vertices=None):
        r"""
        Returns True if the set of vertices defines an independent set and False otherwise.  If no vertices are passed, then we assume the whole graph.

        EXAMPLE:
            sage: G = graphs.EmptyGraph()
            sage: G.add_vertices([0..10])
            sage: G.independent_set()
            True
            sage: graphs.PathGraph(3).independent_set([0,2])
            True
            sage: graphs.PathGraph(3).independent_set([0,1])
            False
        """
        return self.subgraph(vertices).to_simple().size()==0

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
        name -- (must be an explicitly named parameter, i.e., name="complete")
            gives the graph a name
        loops -- boolean, whether to allow loops (ignored if data is an
            instance of the Graph class)
        multiedges -- boolean, whether to allow multiple edges (ignored if
            data is an instance of the Graph class)
        format -- if None, Graph tries to guess- can be several values,
            including:
            'graph6' -- Brendan McKay's graph6 format, in a string (if the
                string has multiple graphs, the first graph is taken)
            'sparse6' -- Brendan McKay's sparse6 format, in a string (if the
                string has multiple graphs, the first graph is taken)
            'adjacency_matrix' -- a square SAGE matrix M, with M[i][j] equal
                to the number of edges \{i,j\}
            'weighted_adjacency_matrix' -- a square SAGE matrix M, with M[i][j]
                equal to the weight of the single edge \{i,j\}
            'incidence_matrix' -- a SAGE matrix, with one column C for each
                edge, where if C represents \{i, j\}, C[i] is -1 and C[j] is 1
            'elliptic_curve_congruence' -- data must be an iterable container
                of elliptic curves, and the graph produced has each curve as a
                vertex (it's Cremona label) and an edge E-F labelled p if and
                only if E is congruent to F mod p
        boundary -- a list of boundary vertices, if empty, graph is considered
            as a 'graph without boundary'

    EXAMPLES:
    We illustrate the first six input formats (the other two
    involve packages that are currently not standard in SAGE):

    1. A NetworkX XGraph:
        sage: import networkx
        sage: g = networkx.XGraph({0:[1,2,3], 2:[5]})
        sage: Graph(g)
        Graph on 5 vertices

    In this single case, we do not make a copy of g, but just wrap the actual
    NetworkX passed.  We do this for performance reasons.

        sage: import networkx
        sage: g = networkx.XGraph({0:[1,2,3], 2:[5]})
        sage: G = Graph(g)
        sage: H = Graph(g)
        sage: G._nxg is H._nxg
        True

    2. A NetworkX graph:
        sage: import networkx
        sage: g = networkx.Graph({0:[1,2,3], 2:[5]})
        sage: DiGraph(g)
        Digraph on 5 vertices

    Note that in this case, we copy the networkX structure.

        sage: import networkx
        sage: g = networkx.Graph({0:[1,2,3], 2:[5]})
        sage: G = Graph(g)
        sage: H = Graph(g)
        sage: G._nxg is H._nxg
        False




    3. A dictionary of dictionaries:
        sage: g = Graph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}); g
        Graph on 5 vertices

    The labels ('x', 'z', 'a', 'out') are labels for edges. For example, 'out' is
    the label for the edge on 2 and 5. Labels can be used as weights, if all the
    labels share some common parent.

    4. A dictionary of lists:
        sage: g = Graph({0:[1,2,3], 2:[5]}); g
        Graph on 5 vertices

    5. A list of vertices and a function describing adjacencies.  Note
       that the list of vertices and the function must be enclosed in
       a list (i.e., [list of vertices, function]).

       Construct the Paley graph over GF(13).

        sage: g=Graph([GF(13), lambda i,j: i!=j and (i-j).is_square()])
        sage: g.vertices()
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        sage: g.adjacency_matrix()
        [0 1 0 1 1 0 0 0 0 1 1 0 1]
        [1 0 1 0 1 1 0 0 0 0 1 1 0]
        [0 1 0 1 0 1 1 0 0 0 0 1 1]
        [1 0 1 0 1 0 1 1 0 0 0 0 1]
        [1 1 0 1 0 1 0 1 1 0 0 0 0]
        [0 1 1 0 1 0 1 0 1 1 0 0 0]
        [0 0 1 1 0 1 0 1 0 1 1 0 0]
        [0 0 0 1 1 0 1 0 1 0 1 1 0]
        [0 0 0 0 1 1 0 1 0 1 0 1 1]
        [1 0 0 0 0 1 1 0 1 0 1 0 1]
        [1 1 0 0 0 0 1 1 0 1 0 1 0]
        [0 1 1 0 0 0 0 1 1 0 1 0 1]
        [1 0 1 1 0 0 0 0 1 1 0 1 0]

       Construct the line graph of a complete graph.

        sage: g=graphs.CompleteGraph(4)
        sage: line_graph=Graph([g.edges(labels=false), \
                 lambda i,j: len(set(i).intersection(set(j)))>0])
        sage: line_graph.vertices()
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        sage: line_graph.adjacency_matrix()
        [0 1 1 1 1 0]
        [1 0 1 1 0 1]
        [1 1 0 0 1 1]
        [1 1 0 0 1 1]
        [1 0 1 1 0 1]
        [0 1 1 1 1 0]

    6. A numpy matrix or ndarray:
        sage: import numpy
        sage: A = numpy.array([[0,1,1],[1,0,1],[1,1,0]])
        sage: Graph(A)
        Graph on 3 vertices

    7. A graph6 or sparse6 string:
    SAGE automatically recognizes whether a string is in graph6 or sage6 format:

        sage: s = ':I`AKGsaOs`cI]Gb~'
        sage: Graph(s)
        Looped multi-graph on 10 vertices

    There are also list functions to take care of lists of graphs:

        sage: s = ':IgMoqoCUOqeb\n:I`AKGsaOs`cI]Gb~\n:I`EDOAEQ?PccSsge\N\n'
        sage: graphs_list.from_sparse6(s)
        [Looped multi-graph on 10 vertices, Looped multi-graph on 10 vertices, Looped multi-graph on 10 vertices]

    8. A SAGE matrix:
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
    def __init__(self, data=None, pos=None, loops=False, format=None, boundary=[], **kwds):
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
            elif isinstance(data, networkx.XGraph):
                self._nxg = data
            elif isinstance(data, networkx.Graph):
                self._nxg = networkx.XGraph(data, selfloops=loops, **kwds)
            elif isinstance(data,list) and len(data)>=2 and callable(data[1]):
                # Pass XGraph a dict of lists describing the adjacencies
                self._nxg = networkx.XGraph(dict([[i]+[[j for j in data[0] if data[1](i,j)]] for i in data[0]]), selfloops=loops, **kwds)
            else:
                self._nxg = networkx.XGraph(data, selfloops=loops, **kwds)
        if format == 'graph6':
            if not isinstance(data, str):
                raise ValueError, 'If input format is graph6, then data must be a string.'
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
            from math import ceil, floor
            from sage.misc.functional import log
            n = data.find('\n')
            if n == -1:
                n = len(data)
            s = data[:n]
            n, s = graph_fast.N_inverse(s[1:])
            k = int(ceil(log(n,2)))
            bits = ''.join([graph_fast.binary(ord(i)-63).zfill(6) for i in s])
            b = []
            x = []
            for i in range(int(floor(len(bits)/(k+1)))):
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
        elif format == 'elliptic_curve_congruence':
            from sage.rings.arith import lcm, prime_divisors, prange
            from sage.misc.misc import prod
            self._nxg = networkx.XGraph(None, selfloops=loops, **kwds)
            curves = list(data)
            self.add_vertices( [curve.cremona_label() for curve in curves] )
            for i in range(self.order()):
                for j in range(i):
                    E = curves[i]
                    F = curves[j]
                    M = E.conductor()
                    N = F.conductor()
                    MN = lcm(M, N)
                    p_MN = prime_divisors(MN)
                    lim = prod([(j^(MN.ord(j)) + j^(MN.ord(j)-1)) for j in p_MN])
                    a_E = E.anlist(lim)
                    a_F = F.anlist(lim)
                    l_list = [p for p in prange(lim) if p not in p_MN ]
                    p_edges = l_list
                    for l in l_list:
                        n = a_E[l] - a_F[l]
                        if n != 0:
                            P = prime_divisors(n)
                            p_edges = [p for p in p_edges if p in P]
                    if len(p_edges) > 0:
                        self.add_edge(E.cremona_label(), F.cremona_label(), str(p_edges)[1:-1])
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

        EXAMPLE:
            sage: g=Graph({0:[0,1,1,2]},loops=True,multiedges=True)
            sage: g==g.copy()
            True


        """
        G = Graph(self._nxg.copy(), name=self._nxg.name, pos=self._pos, boundary=self._boundary)
        return G

    def to_directed(self):
        """
        Returns a directed version of the graph. A single edge becomes two
        edges, one in each direction.

        EXAMPLE:
            sage: graphs.PetersenGraph().to_directed()
            Petersen graph: Digraph on 10 vertices

        """
        return DiGraph(self._nxg.to_directed(), name=self._nxg.name, pos=self._pos, boundary=self._boundary)

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

    def edges(self, labels=True, sort=True):
        """
        Return a list of edges. Each edge is a triple (u,v,l) where u
        and v are vertices and l is a label.

        INPUT:
            labels -- (bool; default: True) if False, each edge is a
                      tuple (u,v) of vertices.
            sort -- (bool; default: True) if True, ensure that the list
                    of edges is sorted.

        OUTPUT:
            A list of tuples.  It is safe to change the returned list.

        EXAMPLES:
            sage: graphs.DodecahedralGraph().edges()
            [(0, 1, None), (0, 10, None), (0, 19, None), (1, 2, None), (1, 8, None), (2, 3, None), (2, 6, None), (3, 4, None), (3, 19, None), (4, 5, None), (4, 17, None), (5, 6, None), (5, 15, None), (6, 7, None), (7, 8, None), (7, 14, None), (8, 9, None), (9, 10, None), (9, 13, None), (10, 11, None), (11, 12, None), (11, 18, None), (12, 13, None), (12, 16, None), (13, 14, None), (14, 15, None), (15, 16, None), (16, 17, None), (17, 18, None), (18, 19, None)]

            sage: graphs.DodecahedralGraph().edges(labels=False)
            [(0, 1), (0, 10), (0, 19), (1, 2), (1, 8), (2, 3), (2, 6), (3, 4), (3, 19), (4, 5), (4, 17), (5, 6), (5, 15), (6, 7), (7, 8), (7, 14), (8, 9), (9, 10), (9, 13), (10, 11), (11, 12), (11, 18), (12, 13), (12, 16), (13, 14), (14, 15), (15, 16), (16, 17), (17, 18), (18, 19)]
        """
        L = self._nxg.edges()
        if sort:
            L.sort()
        if labels:
            return L
        else:
            return [(u,v) for u,v,_ in L]


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
            4
            2
            3
            2
            3
            3
            2
            3
            sage: for i in G.degree_iterator(labels=True):
            ...    print i
            ((0, 1), 3)
            ((1, 2), 4)
            ((0, 0), 2)
            ((2, 1), 3)
            ((1, 1), 4)
            ((2, 0), 2)
            ((1, 3), 3)
            ((2, 3), 2)
            ((2, 2), 3)
            ((1, 0), 3)
            ((0, 3), 2)
            ((0, 2), 3)

        """
        return self._nxg.degree_iter(vertices, with_labels=labels)

    ### Centrality

    def centrality_betweenness(self, normalized=True):
        r"""
        Returns the betweenness centrality (fraction of number of shortest
        paths that go through each vertex) as a dictionary keyed by vertices.
        The betweenness is normalized by default to be in range (0,1).  This
        wraps Networkx's implementation of the algorithm described in [1].

        Measures of the centrality of a vertex within a graph determine the
        relative importance of that vertex to its graph.  Vertices that occur
        on more shortest paths between other vertices have higher betweenness
        than vertices that occur on less.

        INPUT:
            normalized -- boolean (default True) - if set to False, result
                          is not normalized.

        REFERENCE:
            [1] Ulrik Brandes. (2003). Faster Evaluation of Shortest-Path
                Based Centrality Indices. [Online] Available:
                http://citeseer.nj.nec.com/brandes00faster.html

        EXAMPLES:
            sage: (graphs.ChvatalGraph()).centrality_betweenness()
            {0: 0.069696969696969688, 1: 0.069696969696969688, 2: 0.060606060606060601, 3: 0.060606060606060601, 4: 0.069696969696969688, 5: 0.069696969696969688, 6: 0.060606060606060601, 7: 0.060606060606060601, 8: 0.060606060606060601, 9: 0.060606060606060601, 10: 0.060606060606060601, 11: 0.060606060606060601}
            sage: (graphs.ChvatalGraph()).centrality_betweenness(normalized=False)
            {0: 7.6666666666666661, 1: 7.6666666666666661, 2: 6.6666666666666661, 3: 6.6666666666666661, 4: 7.6666666666666661, 5: 7.6666666666666661, 6: 6.6666666666666661, 7: 6.6666666666666661, 8: 6.6666666666666661, 9: 6.6666666666666661, 10: 6.6666666666666661, 11: 6.6666666666666661}
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage.: D.show(figsize=[2,2])
            sage: D = D.to_undirected()
            sage.: D.show(figsize=[2,2])
            sage: D.centrality_betweenness()
            {0: 0.16666666666666666, 1: 0.16666666666666666, 2: 0.0, 3: 0.0}

        """
        import networkx
        return networkx.betweenness_centrality(self._nxg, normalized)

    def centrality_degree(self, v=None):
        r"""
        Returns the degree centrality (fraction of vertices connected to) as
        a dictionary of values keyed by vertex.  The degree centrality is
        normalized to be in range (0,1).

        Measures of the centrality of a vertex within a graph determine the
        relative importance of that vertex to its graph.  Degree centrality
        measures the number of links incident upon a vertex.

        INPUT:
            v -- a vertex label (to find degree centrality of only one vertex)

        EXAMPLES:
            sage: (graphs.ChvatalGraph()).centrality_degree()
            {0: 0.36363636363636365, 1: 0.36363636363636365, 2: 0.36363636363636365, 3: 0.36363636363636365, 4: 0.36363636363636365, 5: 0.36363636363636365, 6: 0.36363636363636365, 7: 0.36363636363636365, 8: 0.36363636363636365, 9: 0.36363636363636365, 10: 0.36363636363636365, 11: 0.36363636363636365}
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage.: D.show(figsize=[2,2])
            sage: D = D.to_undirected()
            sage.: D.show(figsize=[2,2])
            sage: D.centrality_degree()
            {0: 1.0, 1: 1.0, 2: 0.66666666666666663, 3: 0.66666666666666663}
            sage: D.centrality_degree(v=1)
            1.0
        """
        import networkx
        return networkx.degree_centrality(self._nxg, v)

    def centrality_closeness(self, v=None):
        r"""
        Returns the closeness centrality (1/average distance to all vertices) as
        a dictionary of values keyed by vertex.  The degree centrality is
        normalized to be in range (0,1).

        Measures of the centrality of a vertex within a graph determine the
        relative importance of that vertex to its graph.  'Closeness centrality
        may be defined as the total graph-theoretic distance of a given vertex
        from all other vertices... Closeness is an inverse measure of centrality
        in that a larger value indicates a less central actor while a smaller
        value indicates a more central actor,' [1].

        INPUT:
            v -- a vertex label (to find degree centrality of only one vertex)

        REFERENCE:
            [1] Stephen P Borgatti. (1995). Centrality and AIDS. [Online]
                Available: http://www.analytictech.com/networks/centaids.htm

        EXAMPLES:
            sage: (graphs.ChvatalGraph()).centrality_closeness()
            {0: 0.61111111111111116, 1: 0.61111111111111116, 2: 0.61111111111111116, 3: 0.61111111111111116, 4: 0.61111111111111116, 5: 0.61111111111111116, 6: 0.61111111111111116, 7: 0.61111111111111116, 8: 0.61111111111111116, 9: 0.61111111111111116, 10: 0.61111111111111116, 11: 0.61111111111111116}
            sage: D = DiGraph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage.: D.show(figsize=[2,2])
            sage: D = D.to_undirected()
            sage.: D.show(figsize=[2,2])
            sage: D.centrality_closeness()
            {0: 1.0, 1: 1.0, 2: 0.75, 3: 0.75}
            sage: D.centrality_closeness(v=1)
            1.0
        """
        import networkx
        return networkx.closeness_centrality(self._nxg, v)

    ### Spectrum

    def spectrum(self, laplacian=False):
        """
        Returns the spectrum of the graph, the eigenvalues of the adjacency
        matrix

        INPUT:
            laplacian -- if True, use the Laplacian matrix instead (see
                self.kirchhoff_matrix())

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: P.spectrum()
	    [-2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0]
            sage: P.spectrum(laplacian=True)   # random low-order bits (at least for first eigenvalue)
	    [-1.41325497305e-16, 2.0, 2.0, 2.0, 2.0, 2.0, 5.0, 5.0, 5.0, 5.0]

        """
        from sage.matrix.constructor import matrix
        from sage.rings.real_double import RDF
        if laplacian:
            M = self.kirchhoff_matrix()
        else:
            M = self.am()
        M = matrix(RDF, M.rows())
        E = M.eigen_left()[0]
        v = [e.real() for e in E]
	v.sort()
	return v

    ### Representations

    def adjacency_matrix(self, sparse=True, boundary_first=False, over_integers=False):
        """
        Returns the adjacency matrix of the graph. Each vertex is
        represented by its position in the list returned by the vertices()
        function.

        If the graph allows multiple edges, then the returned matrix is over
        the integers, otherwise it is over the ring with two elements.

        INPUT:
            sparse -- whether to represent with a sparse matrix
            boundary_first -- whether to represent the boundary vertices in
                the upper left block
            over_integers -- overrides checking multiple edges

        EXAMPLE:
            sage: G = graphs.CubeGraph(4)
            sage: G.adjacency_matrix()
            [0 1 1 0 1 0 0 0 1 0 0 0 0 0 0 0]
            [1 0 0 1 0 1 0 0 0 1 0 0 0 0 0 0]
            [1 0 0 1 0 0 1 0 0 0 1 0 0 0 0 0]
            [0 1 1 0 0 0 0 1 0 0 0 1 0 0 0 0]
            [1 0 0 0 0 1 1 0 0 0 0 0 1 0 0 0]
            [0 1 0 0 1 0 0 1 0 0 0 0 0 1 0 0]
            [0 0 1 0 1 0 0 1 0 0 0 0 0 0 1 0]
            [0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 1]
            [1 0 0 0 0 0 0 0 0 1 1 0 1 0 0 0]
            [0 1 0 0 0 0 0 0 1 0 0 1 0 1 0 0]
            [0 0 1 0 0 0 0 0 1 0 0 1 0 0 1 0]
            [0 0 0 1 0 0 0 0 0 1 1 0 0 0 0 1]
            [0 0 0 0 1 0 0 0 1 0 0 0 0 1 1 0]
            [0 0 0 0 0 1 0 0 0 1 0 0 1 0 0 1]
            [0 0 0 0 0 0 1 0 0 0 1 0 1 0 0 1]
            [0 0 0 0 0 0 0 1 0 0 0 1 0 1 1 0]

        """
        n = len(self._nxg.adj)
        if boundary_first:
            verts = self.vertices(boundary_first=True)
        else:
            verts = self.vertices()
        D = {}
        for i,j,l in self.edge_iterator():
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
        if self.multiple_edges() or over_integers:
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
            [ 0  1  0  0  0  0  1 -1  0  0  0  0]
            [ 0  0  0  1  0 -1 -1  0  0  0  0  0]
            [-1 -1 -1  0  0  0  0  0  0  0  0  0]
            [ 1  0  0 -1 -1  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  1  0  0  1 -1]
            [ 0  0  0  0  0  1  0  0  1  0  0  1]
            [ 0  0  1  0  0  0  0  0  0  1 -1  0]
            [ 0  0  0  0  1  0  0  0 -1 -1  0  0]

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

    def weighted_adjacency_matrix(self, sparse=True, boundary_first=False):
        """
        Returns the weighted adjacency matrix of the graph. Each vertex is
        represented by its position in the list returned by the vertices()
        function.

        EXAMPLES:
            sage: G = Graph()
            sage: G.add_edges([(0,1,1),(1,2,2),(0,2,3),(0,3,4)])
            sage: M = G.weighted_adjacency_matrix(); M
            [0 1 3 4]
            [1 0 2 0]
            [3 2 0 0]
            [4 0 0 0]
            sage: H = Graph(data=M, format='weighted_adjacency_matrix')
            sage: H == G
            True

        """
        if self.multiple_edges():
            raise NotImplementedError, "Don't know how to represent weights for a multigraph."

        n = len(self._nxg.adj)
        if boundary_first:
            verts = self.vertices(boundary_first=True)
        else:
            verts = self.vertices()
        D = {}
        for e in self.edge_iterator():
            i,j,l = e
            i = verts.index(i)
            j = verts.index(j)
            D[(i,j)] = l
            D[(j,i)] = l
        from sage.matrix.constructor import matrix
        M = matrix(D, sparse=sparse)
        return M

    def kirchhoff_matrix(self, weighted=False, boundary_first=False):
        """
        Returns the Kirchhoff matrix (a.k.a. the Laplacian) of the graph.

        The Kirchhoff matrix is defined to be D - M, where D is the diagonal
        degree matrix (each diagonal entry is the degree of the corresponding
        vertex), and M is the adjacency matrix.

        If weighted == True, the weighted adjacency matrix is used for M, and
        the diagonal entries are the row-sums of M.

        AUTHOR:
            Tom Boothby

        EXAMPLES:
            sage: G = Graph()
            sage: G.add_edges([(0,1,1),(1,2,2),(0,2,3),(0,3,4)])
            sage: M = G.kirchhoff_matrix(weighted=True); M
            [ 8 -1 -3 -4]
            [-1  3 -2  0]
            [-3 -2  5  0]
            [-4  0  0  4]
            sage: M = G.kirchhoff_matrix(); M
            [ 3 -1 -1 -1]
            [-1  2 -1  0]
            [-1 -1  2  0]
            [-1  0  0  1]
            sage: G.set_boundary([2,3])
            sage: M = G.kirchhoff_matrix(weighted=True, boundary_first=True); M
            [ 5  0 -3 -2]
            [ 0  4 -4  0]
            [-3 -4  8 -1]
            [-2  0 -1  3]
            sage: M = G.kirchhoff_matrix(boundary_first=True); M
            [ 2  0 -1 -1]
            [ 0  1 -1  0]
            [-1 -1  3 -1]
            [-1  0 -1  2]

        """
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import IntegerRing

        if weighted:
            M = self.weighted_adjacency_matrix(boundary_first=boundary_first)
        else:
            M = self.adjacency_matrix(boundary_first=boundary_first, over_integers=True)
        A = list(-M)
        S = [sum(M[i]) for i in range(M.nrows())]
        for i in range(len(A)):
            A[i][i] = S[i]
        return M.parent()(A)

    def is_circular_planar(self, ordered=True):
        """
        Returns True if a graph with boundary is circular planar, and
        False otherwise.  A graph (with nonempty boundary) is circular
        planar if it has a planar embedding in which all boundary vertices
        can be drawn in order on a disc boundary, with all the interior
        vertices drawn inside the disc.

        Note -- This function assumes that the graph has nonempty
                boundary.  (Circular Planarity has no definition for
                graphs without boundary).
             -- The current version relies on computing the genus of a
                slightly modified graph so it is time-expensive and not
                reasonable to use for graphs with > 12 vertices.
             -- Also since the current version relies on computing the
                genus, it is necessary that the graph be connected in
                order to use Euler's formula.

        INPUT:
            ordered -- whether or not to consider the order of the boundary
                       (set ordered=False to see if there is any possible
                       boundary order that will satisfy circular planarity)

        EXAMPLES:
            sage: g439 = Graph({1:[5,7], 2:[5,6], 3:[6,7], 4:[5,6,7]})
            sage: g439.set_boundary([1,2,3,4])
            sage.: g439.show(figsize=[2,2], vertex_labels=True, vertex_size=175)
            sage: g439.is_circular_planar()
            False
            sage: g439.set_boundary([1,2,3])
            sage: g439.is_circular_planar()
            True

        Order matters:
            sage: K23 = graphs.CompleteBipartiteGraph(2,3)
            sage: K23.set_boundary([0,1,2,3])
            sage: K23.is_circular_planar()
            False
            sage: K23.set_boundary([0,2,1,3]) # Diff Order!
            sage: K23.is_circular_planar()
            True
            sage: K23.is_circular_planar(ordered=False)
            True
        """
        if not self.is_connected():
            raise TypeError("Graph must be connected to use Euler's Formula to compute minimal genus.")
        from sage.rings.infinity import Infinity
        from sage.combinat.all import CyclicPermutationsOfPartition
        from sage.graphs.graph_genus1 import trace_faces, nice_copy

        graph = nice_copy(self)
        boundary = graph.get_boundary()

        extra = 0
        while graph.has_vertex(extra):
            extra=extra+1
        graph.add_vertex(extra)

        for vertex in boundary:
            graph.add_edge(vertex,extra)

        verts = len(graph.vertices())
        edges = len(graph.edges())

        # Construct a list of all rotation systems for graph
        part = []
        for vertex in graph.vertices():
            if vertex != extra:
                part.append(graph.neighbors(vertex))
        if not ordered:
            part.append(graph.neighbors(extra))

        all_perms = []
        for p in CyclicPermutationsOfPartition(part):
            if ordered:
                p.append(boundary)
            all_perms.append(p)

        max_faces = -Infinity
        for p in all_perms:
            t = trace_faces(graph, p)
            num = len(t)
            if num > max_faces:
                max_faces = num
        genus = (2 - verts + edges - max_faces)/2
        if genus == 0: return True
        else: return False

    def genus(self):
        """
        Returns the minimal genus of the graph.  The genus of a compact
        surface is the number of handles it has.  The genus of a graph
        is the minimal genus of the surface it can be embedded into.

        Note -- This function uses Euler's formula and thus it is
                necessary to consider only connected graphs.

        EXAMPLES:
            sage: (graphs.PetersenGraph()).genus()
            1
            sage: (graphs.CubeGraph(3)).genus()
            0
            sage: K23 = graphs.CompleteBipartiteGraph(2,3)
            sage: K23.genus()
            0
            sage: K33 = graphs.CompleteBipartiteGraph(3,3)
            sage: K33.genus()
            1
        """
        if not self.is_connected():
            raise TypeError("Graph must be connected to use Euler's Formula to compute minimal genus.")
        from sage.rings.infinity import Infinity
        from sage.combinat.all import CyclicPermutationsOfPartition
        from sage.graphs.graph_genus1 import trace_faces, nice_copy

        graph = nice_copy(self)

        verts = len(graph.vertices())
        edges = len(graph.edges())

        # Construct a list of all rotation systems for graph
        part = []
        for vertex in graph.vertices():
            part.append(graph.neighbors(vertex))

        all_perms = []
        for p in CyclicPermutationsOfPartition(part):
            all_perms.append(p)

        max_faces = -Infinity
        for p in all_perms:
            t = trace_faces(graph, p)
            faces = len(t)
            if faces > max_faces:
                max_faces = faces
        return (2-verts+edges-max_faces)/2

    def interior_paths(self, start, end):
        """
        Returns an exhaustive list of paths (also lists) through
        only interior vertices from vertex start to vertex end in the
        graph.

        Note -- start and end do not necessarily have to be boundary
                vertices.

        INPUT:
            start -- the vertex of the graph to search for paths from
            end -- the vertex of the graph to search for paths to

        EXAMPLES:
            sage: eg1 = Graph({0:[1,2], 1:[4], 2:[3,4], 4:[5], 5:[6]})
            sage: eg1.all_paths(0,6)
            [[0, 1, 4, 5, 6], [0, 2, 4, 5, 6]]
            sage: eg2 = eg1.copy()
            sage: eg2.set_boundary([0,1,3])
            sage: eg2.interior_paths(0,6)
            [[0, 2, 4, 5, 6]]
            sage: eg2.all_paths(0,6)
            [[0, 1, 4, 5, 6], [0, 2, 4, 5, 6]]
            sage: eg3 = graphs.PetersenGraph()
            sage: eg3.set_boundary([0,1,2,3,4])
            sage: eg3.all_paths(1,4)
            [[1, 0, 4],
             [1, 0, 5, 8, 3, 2, 7, 9, 4],
             [1, 0, 5, 8, 3, 4],
             [1, 0, 5, 8, 6, 9, 4],
             [1, 0, 5, 8, 6, 9, 7, 2, 3, 4],
             [1, 0, 5, 7, 9, 4],
             [1, 0, 5, 7, 9, 6, 8, 3, 4],
             [1, 0, 5, 7, 2, 3, 8, 6, 9, 4],
             [1, 0, 5, 7, 2, 3, 4],
             [1, 2, 3, 8, 5, 0, 4],
             [1, 2, 3, 8, 5, 7, 9, 4],
             [1, 2, 3, 8, 6, 9, 4],
             [1, 2, 3, 8, 6, 9, 7, 5, 0, 4],
             [1, 2, 3, 4],
             [1, 2, 7, 9, 4],
             [1, 2, 7, 9, 6, 8, 3, 4],
             [1, 2, 7, 9, 6, 8, 5, 0, 4],
             [1, 2, 7, 5, 0, 4],
             [1, 2, 7, 5, 8, 3, 4],
             [1, 2, 7, 5, 8, 6, 9, 4],
             [1, 6, 8, 3, 2, 7, 9, 4],
             [1, 6, 8, 3, 2, 7, 5, 0, 4],
             [1, 6, 8, 3, 4],
             [1, 6, 8, 5, 0, 4],
             [1, 6, 8, 5, 7, 9, 4],
             [1, 6, 8, 5, 7, 2, 3, 4],
             [1, 6, 9, 4],
             [1, 6, 9, 7, 2, 3, 8, 5, 0, 4],
             [1, 6, 9, 7, 2, 3, 4],
             [1, 6, 9, 7, 5, 0, 4],
             [1, 6, 9, 7, 5, 8, 3, 4]]
            sage: eg3.interior_paths(1,4)
            [[1, 6, 8, 5, 7, 9, 4], [1, 6, 9, 4]]
        """
        H = self.copy()
        for vertex in self.get_boundary():
            if (vertex != start and vertex != end):
                H.delete_vertex(vertex)
        return H.all_paths(start, end)

    def all_paths(self, start, end):
        """
        Returns a list of all paths (also lists) between a pair of
        vertices (start, end) in the graph.

        EXAMPLES:
            sage: eg1 = Graph({0:[1,2], 1:[4], 2:[3,4], 4:[5], 5:[6]})
            sage: eg1.all_paths(0,6)
            [[0, 1, 4, 5, 6], [0, 2, 4, 5, 6]]
            sage: eg2 = graphs.PetersenGraph()
            sage: eg2.all_paths(1,4)
            [[1, 0, 4],
             [1, 0, 5, 8, 3, 2, 7, 9, 4],
             [1, 0, 5, 8, 3, 4],
             [1, 0, 5, 8, 6, 9, 4],
             [1, 0, 5, 8, 6, 9, 7, 2, 3, 4],
             [1, 0, 5, 7, 9, 4],
             [1, 0, 5, 7, 9, 6, 8, 3, 4],
             [1, 0, 5, 7, 2, 3, 8, 6, 9, 4],
             [1, 0, 5, 7, 2, 3, 4],
             [1, 2, 3, 8, 5, 0, 4],
             [1, 2, 3, 8, 5, 7, 9, 4],
             [1, 2, 3, 8, 6, 9, 4],
             [1, 2, 3, 8, 6, 9, 7, 5, 0, 4],
             [1, 2, 3, 4],
             [1, 2, 7, 9, 4],
             [1, 2, 7, 9, 6, 8, 3, 4],
             [1, 2, 7, 9, 6, 8, 5, 0, 4],
             [1, 2, 7, 5, 0, 4],
             [1, 2, 7, 5, 8, 3, 4],
             [1, 2, 7, 5, 8, 6, 9, 4],
             [1, 6, 8, 3, 2, 7, 9, 4],
             [1, 6, 8, 3, 2, 7, 5, 0, 4],
             [1, 6, 8, 3, 4],
             [1, 6, 8, 5, 0, 4],
             [1, 6, 8, 5, 7, 9, 4],
             [1, 6, 8, 5, 7, 2, 3, 4],
             [1, 6, 9, 4],
             [1, 6, 9, 7, 2, 3, 8, 5, 0, 4],
             [1, 6, 9, 7, 2, 3, 4],
             [1, 6, 9, 7, 5, 0, 4],
             [1, 6, 9, 7, 5, 8, 3, 4]]
        """
        all_paths = []
        paths_helper(start, end, self, all_paths)
        return all_paths

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
            for i in range(len(edges)): # replace edge labels with natural numbers (by index in vertices)
                edges[i] = (vertices.index(edges[i][0]),vertices.index(edges[i][1]))
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
            from math import ceil
            from sage.misc.functional import log
            k = int(ceil(log(n,2)))
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
        sage: G.add_vertices(range(10)); G
        Graph on 10 vertices
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
        sage: G.add_vertices(range(10)); G
        Graph on 10 vertices
        sage.: show(G)
        sage: G.add_path(range(20)[10:20])
        sage.: show(G)
        sage: G.add_path(range(10))
        sage.: show(G)

        """
        self._nxg.add_path(vertices)

    def subgraph(self, vertices=None, inplace=False, create_using=None):
        """
        Returns the subgraph induced by the given vertices.

        INPUT:
        inplace -- Using inplace is True will simply delete the extra vertices
        and edges from the current graph. This will modify the graph, and re-
        turn itself.
        vertices -- Vertices can be a single vertex or an iterable container
        of vertices, e.g. a list, set, graph, file or numeric array.  If not passed, defaults to the entire graph.
        create_using -- Can be an existing graph object or a call to a graph
        object, such as create_using=DiGraph(). Must be a NetworkX object.

        EXAMPLES:
            sage: G = graphs.CompleteGraph(9)
            sage: H = G.subgraph([0,1,2]); H
            Subgraph of (Complete graph): Graph on 3 vertices
            sage: G
            Complete graph: Graph on 9 vertices
            sage: G.subgraph([0,1,2], inplace=True); G
            Subgraph of (Complete graph): Graph on 3 vertices
            sage: G.subgraph()==G
            True

        """
        if inplace:
            self._nxg = self._nxg.subgraph(vertices, inplace, create_using)
        else:
            NXG = self._nxg.subgraph(vertices, inplace, create_using)
            return Graph(NXG)

    ### Visualization

    def write_to_eps(self, filename, iterations=50):
        r"""
        Writes a plot of the graph to filename in eps format.

        It is relatively simple to include this file in a latex document:

        INPUT:
            filename
            iterations -- how many iterations of the spring layout algorithm to
                go through, if applicable

        \code{\\usepackage{graphics}} must appear before the beginning
        of the document, and \code{\\includegraphics {filename.eps}}
        will include it in your latex doc.  Note: you cannot use
        pdflatex to print the resulting document, use TeX and
        Ghostscript or something similar instead.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: P.write_to_eps('sage.eps')
        """
        from sage.graphs.print_graphs import print_graph_eps
        if self._pos is None:
            pos = graph_fast.spring_layout_fast(self, iterations=iterations)
        else:
            pos = self._pos
            keys = pos.keys()
            for v in self.vertices():
                if v not in keys:
                    pos = graph_fast.spring_layout_fast(self, iterations=iterations)
                    break
        xmin = 0.0
        ymin = 0.0
        xmax = -1.0
        ymax = -1.0
        for v in pos:
            x,y = pos[v]
            if (x > xmax):
                xmax = x
            if (x < xmin):
                xmin = x
            if (y > ymax):
                ymax = y
            if (y < ymin):
                ymin = y
        for v in pos:
            pos[v][0] = 1.8*(pos[v][0] - xmin)/(xmax - xmin) - 0.9
            pos[v][1] = 1.8*(pos[v][1] - ymin)/(ymax - ymin) - 0.9
        if filename[-4:] != '.eps':
            filename += '.eps'
        f = open(filename, 'w')
        f.write( print_graph_eps(self.vertices(), self.edge_iterator(), pos) )
        f.close()

    def plot3d(self, bgcolor=(1,1,1),
               vertex_colors=None, vertex_size=0.06,
               edge_colors=None, edge_size=0.02,
               pos3d=None, iterations=50, color_by_label=False, **kwds):
        """
        Plots the graph using Tachyon, and returns a Tachyon object containing
        a representation of the graph.

        INPUT:
            bgcolor -- rgb tuple (default: (1,1,1))
            vertex_size -- float (default: 0.06)
            vertex_colors -- optional dictionary to specify vertex colors:
                each key is a color recognizable by tachyon (rgb tuple
                (default: (1,0,0))), and each corresponding entry is a list of
                vertices. If a vertex is not listed, it looks invisible on the
                resulting plot (it doesn't get drawn).
            edge_colors -- a dictionary specifying edge colors: each key is a
                color recognized by tachyon ( default: (0,0,0) ), and each
                entry is a list of edges.
            edge_size -- float (default: 0.02)
            pos3d -- a position dictionary for the vertices
            iterations -- how many iterations of the spring layout algorithm to
                go through, if applicable
            xres -- resolution
            yres -- resolution
            **kwds -- passed on to the Tachyon command

        EXAMPLES:
            sage: D = graphs.DodecahedralGraph()
            sage: P3D = D.plot3d()
            sage: P3D.save('sage.png') # long time

            sage: G = graphs.PetersenGraph()
            sage: G.plot3d(vertex_colors={(0,0,1):G.vertices()}).save('sage.png') # long time

            sage: C = graphs.CubeGraph(4)
            sage: C.plot3d(edge_colors={(0,1,0):C.edges()}, vertex_colors={(1,1,1):C.vertices()}, bgcolor=(0,0,0)).save('sage.png') # long time

            sage: K = graphs.CompleteGraph(3)
            sage: K.plot3d(edge_colors={(1,0,0):[(0,1,None)], (0,1,0):[(0,2,None)], (0,0,1):[(1,2,None)]}).save('sage.png') # long time

        """
        TT, pos3d = tachyon_vertex_plot(self, bgcolor=bgcolor, vertex_colors=vertex_colors,
                                        vertex_size=vertex_size, pos3d=pos3d, iterations=iterations, **kwds)
        edges = self.edges()

        if color_by_label:
            if edge_colors is  None:
                # do the coloring
                edge_colors = self._color_by_label(format='rgbtuple')

        if edge_colors is None:
            edge_colors = { (0,0,0) : edges }

        i = 0

        for color in edge_colors:
            i += 1
            TT.texture('edge_color_%d'%i, ambient=0.1, diffuse=0.9, specular=0.03, opacity=1.0, color=color)
            for u, v, l in edge_colors[color]:
                TT.fcylinder( (pos3d[u][0],pos3d[u][1],pos3d[u][2]), (pos3d[v][0],pos3d[v][1],pos3d[v][2]), edge_size,'edge_color_%d'%i)

        return TT

    def show3d(self, bgcolor=(1,1,1),
               vertex_colors=None, vertex_size=0.06,
               edge_colors=None, edge_size=0.02,
               pos3d=None, iterations=50, color_by_label=False,
               **kwds):
        """
        Plots the graph using Tachyon, and shows the resulting plot.

        INPUT:
            bgcolor -- rgb tuple (default: (1,1,1))
            vertex_size -- float (default: 0.06)
            vertex_colors -- optional dictionary to specify vertex colors:
                each key is a color recognizable by tachyon (rgb tuple
                (default: (1,0,0))), and each corresponding entry is a list of
                vertices. If a vertex is not listed, it looks invisible on the
                resulting plot (it doesn't get drawn).
            edge_colors -- a dictionary specifying edge colors: each key is a
                color recognized by tachyon ( default: (0,0,0) ), and each
                entry is a list of edges.
            edge_size -- float (default: 0.02)
            pos3d -- a position dictionary for the vertices
            iterations -- how many iterations of the spring layout algorithm to
                go through, if applicable
            xres -- resolution
            yres -- resolution
            **kwds -- passed on to the Tachyon command

        EXAMPLES:
            sage: D = graphs.DodecahedralGraph()
            sage: P3D = D.plot3d()
            sage: P3D.save('sage.png') # long time

            sage: G = graphs.PetersenGraph()
            sage: G.plot3d(vertex_colors={(0,0,1):G.vertices()}).save('sage.png') # long time

            sage: C = graphs.CubeGraph(4)
            sage: C.plot3d(edge_colors={(0,1,0):C.edges()}, vertex_colors={(1,1,1):C.vertices()}, bgcolor=(0,0,0)).save('sage.png') # long time

            sage: K = graphs.CompleteGraph(3)
            sage: K.plot3d(edge_colors={(1,0,0):[(0,1,None)], (0,1,0):[(0,2,None)], (0,0,1):[(1,2,None)]}).save('sage.png') # long time

        """
        self.plot3d(bgcolor=bgcolor, vertex_colors=vertex_colors,
                    edge_colors=edge_colors, vertex_size=vertex_size,
                    edge_size=edge_size, iterations=iterations,
                    color_by_label=color_by_label, **kwds).show()

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

    ### Coloring

    def bipartite_color(self):
        """
        Returns a dictionary with vertices as the keys and the color
        class as the values.  Fails with an error if the graph is not
        bipartite.

        EXAMPLE:

            sage: graphs.CycleGraph(4).bipartite_color()
            {0: 1, 1: 0, 2: 1, 3: 0}
            sage: graphs.CycleGraph(5).bipartite_color()
            Traceback (most recent call last):
            ...
            NetworkXError: graph is not bipartite

        """
        import networkx.generators.bipartite
        return networkx.generators.bipartite.bipartite_color(self._nxg)


    def bipartite_sets(self):
        """
        Returns (X,Y) where X and Y are the nodes in each bipartite
        set of graph G.  Fails with an error if graph is not
        bipartite.

        EXAMPLE:

            sage: graphs.CycleGraph(4).bipartite_sets()
            ([0, 2], [1, 3])
            sage: graphs.CycleGraph(5).bipartite_sets()
            Traceback (most recent call last):
            ...
            NetworkXError: graph is not bipartite

        """
        import networkx.generators.bipartite
        return networkx.generators.bipartite.bipartite_sets(self._nxg)


    def is_bipartite(self):
        """
        Returns True if graph G is bipartite, False if not.

        Traverse the graph G with depth-first-search and color nodes.
        This function uses the corresponding NetworkX function.

        EXAMPLE:

            sage: graphs.CycleGraph(4).is_bipartite()
            True
            sage: graphs.CycleGraph(5).is_bipartite()
            False


        """
        import networkx.generators.bipartite
        return networkx.generators.bipartite.is_bipartite(self._nxg)




    ### Automorphism and isomorphism

    def automorphism_group(self, partition=None, translation=False,
                           verbosity=0):
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
            sage: graphs_query = GraphDatabase()
            sage: L = graphs_query.get_list(num_vertices=4)
            sage.: graphs_list.show_graphs(L)
            sage: for g in L:
            ...    G = g.automorphism_group()
            ...    G.order(), G.gens()
            (24, ((2,3), (1,2), (1,4)))
            (4, ((2,3), (1,4)))
            (2, ((1,2),))
            (8, ((1,2), (1,4)(2,3)))
            (6, ((1,2), (1,4)))
            (6, ((2,3), (1,2)))
            (2, ((1,4)(2,3),))
            (2, ((1,2),))
            (8, ((2,3), (1,4), (1,3)(2,4)))
            (4, ((2,3), (1,4)))
            (24, ((2,3), (1,2), (1,4)))

            sage: C = graphs.CubeGraph(4)
            sage: G = C.automorphism_group()
            sage: M = G.character_table()
            sage: M.determinant()
            -712483534798848
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
                a,b = search_tree(self, partition, dict=True, lab=False, dig=self.loops(), verbosity=verbosity)
            else:
                a = search_tree(self, partition, dict=False, lab=False, dig=self.loops(), verbosity=verbosity)
            if len(a) != 0:
                a = PermutationGroup([perm_group_elt(aa) for aa in a])
            else:
                a = PermutationGroup([[]])
            if translation:
                return a,b
            else:
                return a

    def is_isomorphic(self, other, certify=False, verbosity=0):
        """
        Tests for isomorphism between self and other.

        INPUT:
            certify -- if True, then output is (a,b), where a is a boolean and b is either a map or
        None.

        EXAMPLES:
            sage: from sage.groups.perm_gps.permgroup_named import SymmetricGroup
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
            sage: a,b = D.is_isomorphic(E, certify=True); a
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
        if certify:
            if self.order() != other.order():
                return False, None
            if self.size() != other.size():
                return False, None
            if sorted(list(self.degree_iterator())) != sorted(list(other.degree_iterator())):
                return False, None
            b,a = self.canonical_label(certify=True, verbosity=verbosity)
            d,c = other.canonical_label(certify=True, verbosity=verbosity)
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
            if self.order() != other.order():
                return False
            if self.size() != other.size():
                return False
            if sorted(list(self.degree_iterator())) != sorted(list(other.degree_iterator())):
                return False
            from sage.graphs.graph_isom import search_tree
            b = self.canonical_label(verbosity=verbosity)
            d = other.canonical_label(verbosity=verbosity)
            return enum(b) == enum(d)

    def canonical_label(self, partition=None, certify=False, verbosity=0):
        """
        Returns the canonical label with respect to the partition. If no
        partition is given, uses the unit partition.

        EXAMPLE:
            sage: D = graphs.DodecahedralGraph()
            sage: E = D.canonical_label(); E
            Dodecahedron: Graph on 20 vertices
            sage: D.canonical_label(certify=True)
            (Dodecahedron: Graph on 20 vertices, {0: 0, 1: 19, 2: 16, 3: 15, 4: 9, 5: 1, 6: 10, 7: 8, 8: 14, 9: 12, 10: 17, 11: 11, 12: 5, 13: 6, 14: 2, 15: 4, 16: 3, 17: 7, 18: 13, 19: 18})
            sage: D.is_isomorphic(E)
            True

        """
        if self.multiple_edges():
            raise NotImplementedError, "Search algorithm does not support multiple edges yet."
        from sage.graphs.graph_isom import search_tree
        if partition is None:
            partition = [self.vertices()]
        if certify:
            a,b,c = search_tree(self, partition, certify=True, dig=self.loops(), verbosity=verbosity)
            return b,c
        else:
            a,b = search_tree(self, partition, dig=self.loops(), verbosity=verbosity)
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
    1. A NetworkX XDiGraph:
        sage: import networkx
        sage: g = networkx.XDiGraph({0:[1,2,3], 2:[5]})
        sage: DiGraph(g)
        Digraph on 5 vertices

    In this single case, we do not make a copy of g, but just wrap the actual
    NetworkX passed.  We do this for performance reasons.

        sage: import networkx
        sage: g = networkx.XDiGraph({0:[1,2,3], 2:[5]})
        sage: G = DiGraph(g)
        sage: H = DiGraph(g)
        sage: G._nxg is H._nxg
        True

    2. A NetworkX digraph:
        sage: import networkx
        sage: g = networkx.DiGraph({0:[1,2,3], 2:[5]})
        sage: DiGraph(g)
        Digraph on 5 vertices

    Note that in this case, we copy the networkX structure.

        sage: import networkx
        sage: g = networkx.DiGraph({0:[1,2,3], 2:[5]})
        sage: G = DiGraph(g)
        sage: H = DiGraph(g)
        sage: G._nxg is H._nxg
        False


    3. A dictionary of dictionaries:
        sage: g = DiGraph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}); g
        Digraph on 5 vertices

    The labels ('x', 'z', 'a', 'out') are labels for edges. For example, 'out' is
    the label for the edge from 2 to 5. Labels can be used as weights, if all the
    labels share some common parent.

    4. A dictionary of lists:
        sage: g = DiGraph({0:[1,2,3], 2:[5]}); g
        Digraph on 5 vertices

    5. A list of vertices and a function describing adjacencies.  Note
       that the list of vertices and the function must be enclosed in
       a list (i.e., [list of vertices, function]).

       We construct a graph on the integers 1 through 12 such that
       there is a directed edge from i to j if and only if i divides j.

        sage: g=DiGraph([[1..12],lambda i,j: i!=j and i.divides(j)])
        sage: g.vertices()
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        sage: g.adjacency_matrix()
        [0 1 1 1 1 1 1 1 1 1 1 1]
        [0 0 0 1 0 1 0 1 0 1 0 1]
        [0 0 0 0 0 1 0 0 1 0 0 1]
        [0 0 0 0 0 0 0 1 0 0 0 1]
        [0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 1]
        [0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0]


    6. A numpy matrix or ndarray:
        sage: import numpy
        sage: A = numpy.array([[0,1,0],[1,0,0],[1,1,0]])
        sage: DiGraph(A)
        Digraph on 3 vertices

    7. A SAGE matrix:
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

    def __init__(self, data=None, pos=None, loops=False, format=None, boundary=[], **kwds):
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
            elif isinstance(data, networkx.XDiGraph):
                self._nxg = data
            elif isinstance(data, networkx.DiGraph):
                self._nxg = networkx.XDiGraph(data, selfloops=loops, **kwds)
            elif isinstance(data, str):
                format = 'dig6'
            elif isinstance(data,list) and len(data)>=2 and callable(data[1]):
                # Pass XGraph a dict of lists describing the adjacencies
                self._nxg = networkx.XDiGraph(dict([[i]+[[j for j in data[0] if data[1](i,j)]] for i in data[0]]), selfloops=loops, **kwds)
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
        elif format == 'dig6':
            if not isinstance(data, str):
                raise ValueError, 'If input format is dig6, then data must be a string.'
            n = data.find('\n')
            if n == -1:
                n = len(data)
            s = data[:n]
            n, s = graph_fast.N_inverse(s)
            m = graph_fast.D_inverse(s, n)
            d = {}
            k = 0
            for i in range(n):
                d[i] = {}
                for j in range(n):
                    if m[k] == '1':
                        d[i][j] = None
                    k += 1
            self._nxg = networkx.XDiGraph(d)
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

        EXAMPLE:
            sage: g=DiGraph({0:[0,1,1,2],1:[0,1]},loops=True,multiedges=True)
            sage: g==g.copy()
            True

        """
        G = DiGraph(self._nxg.copy(), name=self._nxg.name, pos=self._pos, boundary=self._boundary)
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
        Returns an undirected version of the graph. Every directed edge becomes an edge.

        EXAMPLE:
            sage: D = DiGraph({0:[1,2],1:[0]})
            sage: G = D.to_undirected()
            sage: D.edges(labels=False)
            [(0, 1), (0, 2), (1, 0)]
            sage: G.edges(labels=False)
            [(0, 1), (0, 2)]

        """
        return Graph(self._nxg.to_undirected(), name=self._nxg.name, pos=self._pos, boundary=self._boundary)

    ### General Properties

    def is_directed(self):
        """
        Since digraph is directed, returns True.

        """
        return True

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
        return iter(set(self._nxg.successors_iter(vertex)) \
                    | set(self._nxg.predecessors_iter(vertex)))

    ### Edge Handlers

    def add_edge(self, u, v=None, label=None):
        """
        Adds an edge from u to v.

        INPUT:
        The following forms are all accepted by NetworkX:
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
        sage: G = DiGraph()
        sage: G.add_edge((1,2),'label')
        sage: G.networkx_graph().adj           # random output order
        {'label': {}, (1, 2): {'label': None}}

        Use one of these instead:
        sage: G = DiGraph()
        sage: G.add_edge((1,2), label="label")
        sage: G.networkx_graph().adj           # random output order
        {1: {2: 'label'}, 2: {}}

        sage: G = DiGraph()
        sage: G.add_edge(1,2,'label')
        sage: G.networkx_graph().adj           # random output order
        {1: {2: 'label'}, 2: {}}

        """
        self._nxg.add_edge(u, v, label)

    def add_edges(self, edges):
        """
        Add edges from an iterable container.

        EXAMPLE:
            sage: G = graphs.DodecahedralGraph().to_directed()
            sage: H = DiGraph()
            sage: H.add_edges( G.edge_iterator() ); H
            Digraph on 20 vertices

        """
        self._nxg.add_edges_from( edges )

    def delete_edge(self, u, v=None, label=None):
        r"""
        Delete the edge from u to v, return silently if vertices or edge does
        not exist.

        INPUT:
        The following forms are all accepted:

        G.delete_edge( 1, 2 )
        G.delete_edge( (1, 2) )
        G.delete_edges( [ (1, 2) ] )
        G.delete_edge( 1, 2, 'label' )
        G.delete_edge( (1, 2, 'label') )
        G.delete_edges( [ (1, 2, 'label') ] )

        EXAMPLES:
            sage: D = graphs.CompleteGraph(19).to_directed()
            sage: D.size()
            342
            sage: D.delete_edge( 1, 2 )
            sage: D.delete_edge( (3, 4) )
            sage: D.delete_edges( [ (5, 6), (7, 8) ] )
            sage: D.delete_edge( 9, 10, 'label' )
            sage: D.delete_edge( (11, 12, 'label') )
            sage: D.delete_edges( [ (13, 14, 'label') ] )
            sage: D.size()
            335
            sage: D.has_edge( (11, 12) )
            False

            Note that even though the edge (11, 12) has no label, it still gets
            deleted: NetworkX does not pay attention to labels here.

        """
        self._nxg.delete_edge(u, v, label)

    def delete_edges(self, edges):
        """
        Delete edges from an iterable container.

        EXAMPLE:
            sage: K12 = graphs.CompleteGraph(12).to_directed()
            sage: K4 = graphs.CompleteGraph(4).to_directed()
            sage: K12.size()
            132
            sage: K12.delete_edges(K4.edge_iterator())
            sage: K12.size()
            120

        """
        self._nxg.delete_edges_from(edges)

    def delete_multiedge(self, u, v):
        """
        Deletes all edges from u to v.

        EXAMPLE:
            sage: D = DiGraph(multiedges=True)
            sage: D.add_edges([(0,1), (0,1), (0,1), (1,0), (1,2), (2,3)])
            sage: D.edges()
            [(0, 1, None), (0, 1, None), (0, 1, None), (1, 0, None), (1, 2, None), (2, 3, None)]
            sage: D.delete_multiedge( 0, 1 )
            sage: D.edges()
            [(1, 0, None), (1, 2, None), (2, 3, None)]

        """
        self._nxg.delete_multiedge(u, v)

    def edges(self, labels=True, sort=True):
        """
        Return a list of edges. Each edge is a triple (u,v,l) where the edge is
        from u to v, with label l.

        INPUT:
            labels -- (bool; default: True) if False, each edge is a
                      tuple (u,v) of vertices.
            sort -- (bool; default: True) if True, ensure that the list
                    of edges is sorted.

        OUTPUT:
            A list of tuples.  It is safe to change the returned list.


        EXAMPLES:
            sage: D = graphs.DodecahedralGraph().to_directed()
            sage: D.edges()
            [(0, 1, None), (0, 10, None), (0, 19, None), (1, 0, None), (1, 2, None), (1, 8, None), (2, 1, None), (2, 3, None), (2, 6, None), (3, 2, None), (3, 4, None), (3, 19, None), (4, 3, None), (4, 5, None), (4, 17, None), (5, 4, None), (5, 6, None), (5, 15, None), (6, 2, None), (6, 5, None), (6, 7, None), (7, 6, None), (7, 8, None), (7, 14, None), (8, 1, None), (8, 7, None), (8, 9, None), (9, 8, None), (9, 10, None), (9, 13, None), (10, 0, None), (10, 9, None), (10, 11, None), (11, 10, None), (11, 12, None), (11, 18, None), (12, 11, None), (12, 13, None), (12, 16, None), (13, 9, None), (13, 12, None), (13, 14, None), (14, 7, None), (14, 13, None), (14, 15, None), (15, 5, None), (15, 14, None), (15, 16, None), (16, 12, None), (16, 15, None), (16, 17, None), (17, 4, None), (17, 16, None), (17, 18, None), (18, 11, None), (18, 17, None), (18, 19, None), (19, 0, None), (19, 3, None), (19, 18, None)]
            sage: D.edges(labels = False)
            [(0, 1), (0, 10), (0, 19), (1, 0), (1, 2), (1, 8), (2, 1), (2, 3), (2, 6), (3, 2), (3, 4), (3, 19), (4, 3), (4, 5), (4, 17), (5, 4), (5, 6), (5, 15), (6, 2), (6, 5), (6, 7), (7, 6), (7, 8), (7, 14), (8, 1), (8, 7), (8, 9), (9, 8), (9, 10), (9, 13), (10, 0), (10, 9), (10, 11), (11, 10), (11, 12), (11, 18), (12, 11), (12, 13), (12, 16), (13, 9), (13, 12), (13, 14), (14, 7), (14, 13), (14, 15), (15, 5), (15, 14), (15, 16), (16, 12), (16, 15), (16, 17), (17, 4), (17, 16), (17, 18), (18, 11), (18, 17), (18, 19), (19, 0), (19, 3), (19, 18)]
        """
        L = self._nxg.edges()
        if sort:
            L.sort()
        if labels:
            return L
        else:
            return [(u,v) for u,v,_ in L]

    def edge_boundary(self, vertices1, vertices2=None, labels=True):
        """
        Returns a list of edges (u,v,l) with u in vertices1 and v in vertices2.
        If vertices2 is None, then it is set to the complement of vertices1.

        INPUT:
        labels -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLE:
            sage: K = graphs.CompleteBipartiteGraph(9,3).to_directed()
            sage: len(K.edge_boundary( [0,1,2,3,4,5,6,7,8], [9,10,11] ))
            27
            sage: K.size()
            54

            Note that the edge boundary preserves direction: compare this example to
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

    def edge_iterator(self, vertices=None):
        """
        Returns an iterator over the edges pointing out of the given
        set of vertices. If vertices is None, then returns an iterator over
        all edges.

        EXAMPLE:
            sage: D = DiGraph( { 0 : [1,2], 1: [0] } )
            sage: for i in D.edge_iterator([0]):
            ...    print i
            (0, 1, None)
            (0, 2, None)

        """
        return self._nxg.edges_iter(vertices)

    def incoming_edge_iterator(self, vertices=None):
        """
        Return an iterator over all arriving edges from vertices, or over all
        edges if vertices is None.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: for a in D.incoming_edge_iterator([0]):
            ...    print a
            (1, 0, None)
            (4, 0, None)

        """
        return self._nxg.in_edges_iter(vertices)

    def incoming_edges(self, vertices=None, labels=True):
        """
        Returns a list of edges arriving at vertices.

        INPUT:
        labels -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.incoming_edges([0])
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

    def outgoing_edge_iterator(self, vertices=None):
        """
        Return an iterator over all departing edges from vertices, or over all
        edges if vertices is None.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: for a in D.outgoing_edge_iterator([0]):
            ...    print a
            (0, 1, None)
            (0, 2, None)
            (0, 3, None)

        """
        return self._nxg.out_edges_iter(vertices)

    def outgoing_edges(self, vertices=None, labels=True):
        """
        Returns a list of edges departing from vertices.

        INPUT:
        labels -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.outgoing_edges([0])
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

    def has_edge(self, u, v=None, label=None):
        """
        Returns True if there is an edge from u to v, False otherwise.

        INPUT:
        The following forms are accepted by NetworkX:

        D.has_edge( 1, 2 )
        D.has_edge( (1, 2) )
        D.has_edge( 1, 2, 'label' )

        EXAMPLE:
            sage: DiGraph().has_edge(9,2)
            False

        """
        return self._nxg.has_edge(u,v)

    def set_edge_label(self, u, v, l):
        """
        Set the edge label of a given edge.

        INPUT:
            u, v -- the vertices (and direction) of the edge
            l -- the new label

        EXAMPLE:
            sage: SD = DiGraph( { 1:[18,2], 2:[5,3], 3:[4,6], 4:[7,2], 5:[4], 6:[13,12], 7:[18,8,10], 8:[6,9,10], 9:[6], 10:[11,13], 11:[12], 12:[13], 13:[17,14], 14:[16,15], 15:[2], 16:[13], 17:[15,13], 18:[13] } )
            sage: SD.set_edge_label(1, 18, 'discrete')
            sage: SD.set_edge_label(4, 7, 'discrete')
            sage: SD.set_edge_label(2, 5, 'h = 0')
            sage: SD.set_edge_label(7, 18, 'h = 0')
            sage: SD.set_edge_label(7, 10, 'aut')
            sage: SD.set_edge_label(8, 10, 'aut')
            sage: SD.set_edge_label(8, 9, 'label')
            sage: SD.set_edge_label(8, 6, 'no label')
            sage: SD.set_edge_label(13, 17, 'k > h')
            sage: SD.set_edge_label(13, 14, 'k = h')
            sage: SD.set_edge_label(17, 15, 'v_k finite')
            sage: SD.set_edge_label(14, 15, 'v_k m.c.r.')
            sage: posn = {1:[ 3,-3],  2:[0,2],  3:[0, 13],  4:[3,9],  5:[3,3],  6:[16, 13], 7:[6,1],  8:[6,6],  9:[6,11], 10:[9,1], 11:[10,6], 12:[13,6], 13:[16,2], 14:[10,-6], 15:[0,-10], 16:[14,-6], 17:[16,-10], 18:[6,-4]}
            sage: SD.plot(pos=posn, vertex_size=400, vertex_colors={'#FFFFFF':range(1,19)}, edge_labels=True).save('search_tree.png')

        """
        if self.has_edge(u, v):
            self._nxg.adj[u][v] = l

    def edge_label(self, u, v=None):
        """
        Returns the label of an edge.

        EXAMPLE:
            sage: D = DiGraph({0 : {1 : 'edgelabel'}})
            sage: D.edges(labels=False)
            [(0, 1)]
            sage: D.edge_label( 0, 1 )
            'edgelabel'

        """
        return self._nxg.get_edge(u,v)

    def edge_labels(self):
        """
        Returns a list of edge labels.

        EXAMPLE:
            sage: G = DiGraph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}})
            sage: G.edge_labels()
            ['x', 'z', 'a', 'out']

        """
        labels = []
        for u,v,l in self.edges():
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

    def remove_multiple_edges(self):
        """
        Removes all multiple edges, retaining one edge for each.

        EXAMPLE:
            sage: D = DiGraph(multiedges=True)
            sage: D.add_edges( [ (0,1), (0,1), (0,1), (0,1), (1,2) ] )
            sage: D.edges(labels=False)
            [(0, 1), (0, 1), (0, 1), (0, 1), (1, 2)]
            sage: D.remove_multiple_edges()
            sage: D.edges(labels=False)
            [(0, 1), (1, 2)]

        """
        self._nxg.remove_all_multiedges()

    def remove_loops(self, vertices=None):
        """
        Removes loops on vertices in vertices. If vertices is None, removes all loops.

        EXAMPLE:
            sage: D = DiGraph(loops=True)
            sage: D.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: D.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: D.remove_loops()
            sage: D.edges(labels=False)
            [(2, 3)]
            sage: D.loops()
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
            sage: D = DiGraph(loops=True)
            sage: D.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: D.loop_edges()
            [(0, 0, None), (1, 1, None), (2, 2, None), (3, 3, None)]

        """
        return self._nxg.selfloop_edges()

    def number_of_loops(self):
        """
        Returns the number of edges that are loops.

        EXAMPLE:
            sage: D = DiGraph(loops=True)
            sage: D.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: D.edges(labels=False)
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
            6
            4
            4
            4
            6
            sage: for i in D.degree_iterator(labels=True):
            ...    print i
            ((0, 1), 6)
            ((1, 2), 6)
            ((0, 0), 4)
            ((0, 2), 6)
            ((1, 3), 4)
            ((1, 0), 4)
            ((0, 3), 4)
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
            3
            2
            2
            2
            3
            sage: for i in D.in_degree_iterator(labels=True):
            ...    print i
            ((0, 1), 3)
            ((1, 2), 3)
            ((0, 0), 2)
            ((0, 2), 3)
            ((1, 3), 2)
            ((1, 0), 2)
            ((0, 3), 2)
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
            3
            2
            2
            2
            3
            sage: for i in D.out_degree_iterator(labels=True):
            ...    print i
            ((0, 1), 3)
            ((1, 2), 3)
            ((0, 0), 2)
            ((0, 2), 3)
            ((1, 3), 2)
            ((1, 0), 2)
            ((0, 3), 2)
            ((1, 1), 3)

        """
        return self._nxg.out_degree_iter(vertices, with_labels=labels)

    ### Representations

    def adjacency_matrix(self, sparse=True):
        """
        Returns the adjacency matrix of the digraph. Each vertex is
        represented by its position in the list returned by the vertices()
        function.

        If the graph allows multiple edges, then the returned matrix is over
        the integers, otherwise it is over the ring with two elements.

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
        for i,j,l in self.edge_iterator():
            i = verts.index(i)
            j = verts.index(j)
            if D.has_key((i,j)) and self.multiple_edges():
                D[(i,j)] += 1
            else:
                D[(i,j)] = 1
        from sage.rings.integer_mod_ring import IntegerModRing
        from sage.rings.integer_ring import IntegerRing
        from sage.matrix.constructor import matrix
        if self.multiple_edges():
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
        for i, j, l in self.edge_iterator():
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
        Returns a copy of digraph with edges reversed in direction.

        EXAMPLES:
            sage: D = DiGraph({ 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] })
            sage: D.reverse()
            Reverse of (): Digraph on 6 vertices

        """
        NXG = self._nxg.reverse()
        G = DiGraph(NXG)
        return G

    def subgraph(self, vertices=None, inplace=False, create_using=None):
        """
        Returns the subgraph induced by the given vertices.

        INPUT:
        inplace -- Using inplace is True will simply delete the extra vertices
        and edges from the current graph.
        vertices -- Vertices can be a single vertex or an iterable container
        of vertices, e.g. a list, set, graph, file or numeric array.  If not passed, defaults to the entire graph.
        create_using -- Can be an existing graph object or a call to a graph
        object, such as create_using=DiGraph().

        EXAMPLES:
            sage: D = graphs.CompleteGraph(9).to_directed()
            sage: H = D.subgraph([0,1,2]); H
            Subgraph of (Complete graph): Digraph on 3 vertices
            sage: D
            Complete graph: Digraph on 9 vertices
            sage: D.subgraph([0,1,2], inplace=True); D
            Subgraph of (Complete graph): Digraph on 3 vertices

        """
        if inplace:
            self._nxg = self._nxg.subgraph(vertices, inplace, create_using)
        else:
            NXG = self._nxg.subgraph(vertices, inplace, create_using)
            return DiGraph(NXG)

    ### Visualization

    def plot3d(self, bgcolor=(1,1,1), vertex_colors=None,
               vertex_size=0.06,
               edge_size=0.02,
               edge_size2=0.0325,
               edge_colors=None, pos3d=None,
               color_by_label=False, **kwds):
        """
        Plots the graph using Tachyon, and returns a Tachyon object containing
        a representation of the graph.

        INPUT:
            bgcolor -- rgb tuple (default: (1,1,1))
            vertex_size -- float (default: 0.06)
            vertex_colors -- optional dictionary to specify vertex colors:
                each key is a color recognizable by tachyon (rgb tuple
                (default: (1,0,0))), and each corresponding entry is a list of
                vertices. If a vertex is not listed, it looks invisible on the
                resulting plot (it doesn't get drawn).
            edge_colors -- a dictionary specifying edge colors: each key is a
                color recognized by tachyon ( default: (0,0,0) ), and each
                entry is a list of edges.
            edge_size -- float (default: 0.02)
            edge_size2 -- float (default: 0.0325)
            pos3d -- a position dictionary for the vertices
            iterations -- how many iterations of the spring layout algorithm to
                go through, if applicable
            xres -- resolution
            yres -- resolution
            **kwds -- passed on to the Tachyon command

        EXAMPLE:
        A directed version of the dodecahedron
            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []} )
            sage: D.plot3d().save('sage.png') # long time

            sage: P = graphs.PetersenGraph().to_directed()
            sage: from sage.plot.plot import rainbow
            sage: edges = P.edges()
            sage: R = rainbow(len(edges), 'rgbtuple')
            sage: edge_colors = {}
            sage: for i in range(len(edges)):
            ...       edge_colors[R[i]] = [edges[i]]
            sage: P.plot3d(edge_colors=edge_colors).save('sage.png') # long time

        """
        TT, pos3d = tachyon_vertex_plot(self, bgcolor=bgcolor, vertex_colors=vertex_colors,
                                        vertex_size=vertex_size, pos3d=pos3d, **kwds)
        edges = self.edges()
        i = 0

        if color_by_label:
            if edge_colors is None:
                edge_colors = self._color_by_label(format='rgbtuple')

        if edge_colors is None:
            edge_colors = { (0,0,0) : edges }

        for color in edge_colors:
            i += 1
            TT.texture('edge_color_%d'%i, ambient=0.1, diffuse=0.9, specular=0.03, opacity=1.0, color=color)
            for u,v,l in edge_colors[color]:
                TT.fcylinder( (pos3d[u][0],pos3d[u][1],pos3d[u][2]),
                              (pos3d[v][0],pos3d[v][1],pos3d[v][2]), edge_size,'edge_color_%d'%i)
                TT.fcylinder( (0.25*pos3d[u][0] + 0.75*pos3d[v][0],
                               0.25*pos3d[u][1] + 0.75*pos3d[v][1],
                               0.25*pos3d[u][2] + 0.75*pos3d[v][2],),
                              (pos3d[v][0],pos3d[v][1],pos3d[v][2]), edge_size2,'edge_color_%d'%i)
        return TT

    def show3d(self, bgcolor=(1,1,1), vertex_colors=None,
               vertex_size=0.06,
               edge_size=0.02,
               edge_size2=0.0325,
               edge_colors=None,
               pos3d=None, color_by_label=False, **kwds):
        """
        Plots the graph using Tachyon, and shows the resulting plot.

        INPUT:
            bgcolor -- rgb tuple (default: (1,1,1))
            vertex_size -- float (default: 0.06)
            vertex_colors -- optional dictionary to specify vertex colors:
                each key is a color recognizable by tachyon (rgb tuple
                (default: (1,0,0))), and each corresponding entry is a list of
                vertices. If a vertex is not listed, it looks invisible on the
                resulting plot (it doesn't get drawn).
            edge_colors -- a dictionary specifying edge colors: each key is a
                color recognized by tachyon ( default: (0,0,0) ), and each
                entry is a list of edges.
            edge_size -- float (default: 0.02)
            edge_size2 -- float (default: 0.0325)
            pos3d -- a position dictionary for the vertices
            iterations -- how many iterations of the spring layout algorithm to
                go through, if applicable
            xres -- resolution
            yres -- resolution
            **kwds -- passed on to the Tachyon command

        EXAMPLE:
        A directed version of the dodecahedron
            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []} )
            sage: D.plot3d().save('sage.png') # long time

            sage: P = graphs.PetersenGraph().to_directed()
            sage: from sage.plot.plot import rainbow
            sage: edges = P.edges()
            sage: R = rainbow(len(edges), 'rgbtuple')
            sage: edge_colors = {}
            sage: for i in range(len(edges)):
            ...       edge_colors[R[i]] = [edges[i]]
            sage: P.plot3d(edge_colors=edge_colors).save('sage.png') # long time

        """

        self.plot3d(bgcolor=bgcolor, vertex_colors=vertex_colors,
                    vertex_size=vertex_size, edge_size=edge_size,
                    edge_size2=edge_size2, edge_colors=edge_colors,
                    color_by_label=color_by_label, **kwds).show()

    ### Connected components

    def is_connected(self):
        """
        Indicates whether the digraph is connected. Note that in a digraph,
        the direction of the edges does not effect whether it is connected or
        not.

        EXAMPLE:
            sage: D = DiGraph( { 0 : [1, 2], 1 : [2], 3 : [4, 5], 4 : [5] } )
            sage: D.is_connected()
            False
            sage: D.add_edge(0,3)
            sage: D.is_connected()
            True

        """
        import networkx
        return networkx.component.is_connected(self._nxg.to_undirected())

    def connected_components(self):
        """
        Returns a list of lists of vertices, each list representing a
        connected component. The list is ordered from largest to smallest
        component.

        EXAMPLE:
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: D.connected_components()
            [[0, 1, 2, 3], [4, 5, 6]]

        """
        import networkx
        return networkx.component.connected_components(self._nxg.to_undirected())

    def connected_components_number(self):
        """
        Returns the number of connected components.

        EXAMPLE:
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: D.connected_components_number()
            2

        """
        import networkx
        return networkx.component.number_connected_components(self._nxg.to_undirected())

    def connected_components_subgraphs(self):
        """
        Returns a list of connected components as graph objects.

        EXAMPLE:
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: L = D.connected_components_subgraphs()
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
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: D.connected_component_containing_vertex(0)
            [0, 1, 2, 3]

        """
        import networkx
        return networkx.component.node_connected_component(self._nxg.to_undirected(), vertex)

    ### Automorphism and isomorphism

    def automorphism_group(self, partition=None, translation=False, verbosity=0):
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

        EXAMPLE:
            sage: D = DiGraph( { 0:[1], 1:[2], 2:[3], 3:[4], 4:[0] } )
            sage: D.automorphism_group()
            Permutation Group with generators [(1,2,3,4,5)]

        """
        if self.multiple_edges():
            raise NotImplementedError, "Search algorithm does not support multiple edges yet."
        else:
            from sage.graphs.graph_isom import search_tree, perm_group_elt
            from sage.groups.perm_gps.permgroup import PermutationGroup
            if partition is None:
                partition = [self.vertices()]
            if translation:
                a,b = search_tree(self, partition, dict=True, lab=False, dig=True, verbosity=verbosity)
            else:
                a = search_tree(self, partition, dict=False, lab=False, dig=True, verbosity=verbosity)
            if len(a) != 0:
                a = PermutationGroup([perm_group_elt(aa) for aa in a])
            else:
                a = PermutationGroup([[]])
            if translation:
                return a,b
            else:
                return a

    def is_isomorphic(self, other, certify=False, verbosity=0):
        """
        Tests for isomorphism between self and other.

        INPUT:
            certify -- if True, then output is (a,b), where a is a boolean and b is either a map or
        None.

        EXAMPLE:
            sage: A = DiGraph( { 0 : [1,2] } )
            sage: B = DiGraph( { 1 : [0,2] } )
            sage: A.is_isomorphic(B, certify=True)
            (True, {0: 1, 1: 0, 2: 2})

        """
        if self.multiple_edges():
            raise NotImplementedError, "Search algorithm does not support multiple edges yet."
        from sage.graphs.graph_isom import search_tree
        if certify:
            if self.order() != other.order():
                return False, None
            if self.size() != other.size():
                return False, None
            if sorted(list(self.in_degree_iterator())) != sorted(list(other.in_degree_iterator())):
                return False, None
            if sorted(list(self.out_degree_iterator())) != sorted(list(other.out_degree_iterator())):
                return False, None
            b,a = self.canonical_label(certify=True, verbosity=verbosity)
            d,c = other.canonical_label(certify=True, verbosity=verbosity)
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
            if self.order() != other.order():
                return False
            if self.size() != other.size():
                return False
            if sorted(list(self.in_degree_iterator())) != sorted(list(other.in_degree_iterator())):
                return False
            if sorted(list(self.out_degree_iterator())) != sorted(list(other.out_degree_iterator())):
                return False
            from sage.graphs.graph_isom import search_tree
            b = self.canonical_label(verbosity=verbosity)
            d = other.canonical_label(verbosity=verbosity)
            return enum(b) == enum(d)

    def canonical_label(self, partition=None, certify=False, verbosity=0):
        """
        Returns the canonical label with respect to the partition. If no
        partition is given, uses the unit partition.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: DP = P.to_directed()
            sage: DP.canonical_label().adjacency_matrix()
            [0 0 0 0 0 0 0 1 1 1]
            [0 0 0 0 1 0 1 0 0 1]
            [0 0 0 1 0 0 1 0 1 0]
            [0 0 1 0 0 1 0 0 0 1]
            [0 1 0 0 0 1 0 0 1 0]
            [0 0 0 1 1 0 0 1 0 0]
            [0 1 1 0 0 0 0 1 0 0]
            [1 0 0 0 0 1 1 0 0 0]
            [1 0 1 0 1 0 0 0 0 0]
            [1 1 0 1 0 0 0 0 0 0]


        """
        if self.multiple_edges():
            raise NotImplementedError, "Search algorithm does not support multiple edges yet."
        from sage.graphs.graph_isom import search_tree
        if partition is None:
            partition = [self.vertices()]
        if certify:
            a,b,c = search_tree(self, partition, certify=True, dig=True, verbosity=verbosity)
            return b,c
        else:
            a,b = search_tree(self, partition, dig=True, verbosity=verbosity)
            return b

    ### DIG6 format

    def __bit_vector(self):
        vertices = self.vertices()
        n = len(vertices)
        ns = n*n
        bit_vector = set()
        for e,f,g in self.edge_iterator():
            c = vertices.index(e)
            d = vertices.index(f)
            p = c*n + d
            bit_vector.add(p)
        bit_vector = sorted(bit_vector)
        s = []
        j = 0
        for i in bit_vector:
            s.append( '0'*(i - j) + '1' )
            j = i + 1
        s = "".join(s)
        s += '0'*(ns-len(s))
        return s

    def dig6_string(self):
        """
        Returns the dig6 representation of the digraph as an ASCII string.
        Valid for single (no multiple edges) digraphs on 0 to 262143 vertices.

        EXAMPLE: TODO

        """
        n = self.order()
        if n > 262143:
            raise ValueError, 'dig6 format supports graphs on 0 to 262143 vertices only.'
        elif self.multiple_edges():
            raise ValueError, 'dig6 format does not support multiple edges.'
        else:
            return graph_fast.N(n) + graph_fast.R(self.__bit_vector())

    ### Directed Acyclic Graphs (DAGs)

    def is_directed_acyclic(self):
        """
        Returns whether the digraph is acyclic or not.

        A directed graph is acyclic if for any vertex v, there is no directed
        path that starts and ends at v. Every directed acyclic graph (dag)
        corresponds to a partial ordering of its vertices, however multiple
        dags may lead to the same partial ordering.

        EXAMPLES:
            sage: D = DiGraph({ 0:[1,2,3], 4:[2,5], 1:[8], 2:[7], 3:[7], 5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10] })
            sage: D.plot(layout='circular').save('dag.png')
            sage: D.is_directed_acyclic()
            True

            sage: D.add_edge(9,7)
            sage: D.is_directed_acyclic()
            True

            sage: D.add_edge(7,4)
            sage: D.is_directed_acyclic()
            False

        """
        import networkx
        return networkx.is_directed_acyclic_graph(self._nxg)

    def topological_sort(self):
        """
        Returns a topological sort of the digraph if it is acyclic, and
        raises a TypeError if the digraph contains a directed cycle.

        A topological sort is an ordering of the vertices of the digraph such
        that each vertex comes before all of its successors. That is, if u
        comes before v in the sort, then there may be a directed path from u
        to v, but there will be no directed path from v to u.

        EXAMPLES:
            sage: D = DiGraph({ 0:[1,2,3], 4:[2,5], 1:[8], 2:[7], 3:[7], 5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10] })
            sage: D.plot(layout='circular').save('dag.png')
            sage: D.topological_sort()
            [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]

            sage: D.add_edge(9,7)
            sage: D.topological_sort()
            [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]

            sage: D.add_edge(7,4)
            sage: D.topological_sort()
            Traceback (most recent call last):
            ...
            TypeError: Digraph is not acyclic-- there is no topological sort.

        NOTE:
        There is a recursive version of this in NetworkX, but it has
        problems:
            sage: import networkx
            sage: D = DiGraph({ 0:[1,2,3], 4:[2,5], 1:[8], 2:[7], 3:[7], 5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10] })
            sage: N = D.networkx_graph()
            sage: networkx.topological_sort(N)
            [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]
            sage: networkx.topological_sort_recursive(N) is None
            True

        """
        import networkx
        S = networkx.topological_sort(self._nxg)
        if S is None:
            raise TypeError('Digraph is not acyclic-- there is no topological sort.')
        else:
            return S

    def topological_sort_generator(self):
        """
        Returns a list of all topological sorts of the digraph if it is
        acyclic, and raises a TypeError if the digraph contains a directed
        cycle.

        A topological sort is an ordering of the vertices of the digraph such
        that each vertex comes before all of its successors. That is, if u
        comes before v in the sort, then there may be a directed path from u
        to v, but there will be no directed path from v to u. See also
        Graph.topological_sort().

        AUTHORS:
            Michael W. Hansen -- original implementation
            Robert L. Miller -- wrapping, documentation

        REFERENCE:
            [1] Pruesse, Gara and Ruskey, Frank. Generating Linear Extensions
                Fast. SIAM J. Comput., Vol. 23 (1994), no. 2, pp. 373-386.

        EXAMPLES:
            sage: D = DiGraph({ 0:[1,2], 1:[3], 2:[3,4] })
            sage: D.plot(layout='circular').save('dag.png')
            sage: D.topological_sort_generator()
            [[0, 1, 2, 3, 4], [0, 1, 2, 4, 3], [0, 2, 1, 3, 4], [0, 2, 1, 4, 3], [0, 2, 4, 1, 3]]

            sage: for sort in D.topological_sort_generator():
            ...       for edge in D.edge_iterator():
            ...           u,v,l = edge
            ...           if sort.index(u) > sort.index(v):
            ...               print "This should never happen."

        """
        from sage.graphs.linearextensions import linearExtensions
        try:
            return linearExtensions(self._nxg)
        except:
            raise TypeError('Digraph is not acyclic-- there is no topological sort (or there was an error in sage/graphs/linearextensions.py).')

def tachyon_vertex_plot(g, bgcolor=(1,1,1),
                        vertex_colors=None,
                        vertex_size=0.06,
                        pos3d=None,
                        iterations=50, **kwds):
    import networkx
    from math import sqrt
    from sage.plot.tachyon import Tachyon

    c = [0,0,0]
    r = []
    verts = g.vertices()

    if vertex_colors is None:
        vertex_colors = { (1,0,0) : verts }
    if pos3d is None:
        pos3d = graph_fast.spring_layout_fast(g, dim=3, iterations=iterations)
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
    if r == 0:
        r = 1
    for v in verts:
        pos3d[v][0] = pos3d[v][0]/r
        pos3d[v][1] = pos3d[v][1]/r
        pos3d[v][2] = pos3d[v][2]/r
    TT = Tachyon(camera_center=(1.4,1.4,1.4), antialiasing=13, **kwds)
    TT.light((4,3,2), 0.02, (1,1,1))
    TT.texture('bg', ambient=1, diffuse=1, specular=0, opacity=1.0, color=bgcolor)
    TT.plane((-1.6,-1.6,-1.6), (1.6,1.6,1.6), 'bg')

    i = 0
    for color in vertex_colors:
        i += 1
        TT.texture('node_color_%d'%i, ambient=0.1, diffuse=0.9,
                   specular=0.03, opacity=1.0, color=color)
        for v in vertex_colors[color]:
            TT.sphere((pos3d[v][0],pos3d[v][1],pos3d[v][2]), vertex_size, 'node_color_%d'%i)

    return TT, pos3d

def enum(graph, quick=False):
    """
    Used for isomorphism checking.

    INPUT:
        quick -- now we know that the vertices are 0,1,...,n-1

    CAUTION:
    Enumeration should not be expected to be the same from one version to the
    next. It depends on what ordering you put on the graph, which will be
    consistent as long as the functions used are deterministic. However, if
    the functions change a little bit, for example when a new version of
    NetworkX comes out, then the enumeration may be different as well. For
    example, in moving from NX 0.33 to 0.34, the default ordering of the
    vertices of the cube graph changed, which is reasonable since they are
    strings of the form '011100100010', not simply integers.

    EXAMPLES:
        sage: from sage.graphs.graph import enum
        sage: enum(graphs.DodecahedralGraph())
        646827340296833569479885332381965103655612500627043016896502674924517797573929148319427466126170568392555309533861838850
        sage: enum(graphs.MoebiusKantorGraph())
        29627597595494233374689380190219099810725571659745484382284031717525232288040
        sage: enum(graphs.FlowerSnark())
        645682215283153372602620320081348424178216159521280462146968720908564261127120716040952785862033320307812724373694972050
        sage: enum(graphs.CubeGraph(3))
        7535809024060107030
        sage: enum(graphs.CubeGraph(4))
        47267715876236163882872165742917649077474356975346093231312192918052414226710
        sage: enum(graphs.CubeGraph(5))
        73383767099440499978977371110767635000712632058761917166935299414577892531740903317276163088753429176806392542119313927407270916017052236205784763916564139194844709702330411308151107850749347729631684377105917135158882718579653914986774438564026130776549642103326266773293337981002726861228331644217045614870
        sage: enum(graphs.CubeGraph(6))
        426330773284506488918634734865759041919956843959533319933690768247529558188799395712215446423075008058591366298704228052382156035290705594805355797326008656345430976382884644704793815838747426479646427610678507924554229761217219227470873068998288136471101733055146437055591054200049489396547818542127638584078729984152123059286362563567564461701218389382825275571923027468257894398406386537450007469620559809790444332966965101712836129760203217607202553686523550484118794498888998087405724741434040120926560135219397307574002437217860019320176154290286982859306143571386395302492049075037862488404582194490395200617034630435370849810881764417871098062828471680512013766342500009882083287831864769598280796427782079253352467753414067084099529009749310841716158019483591570426310985975244360481910644650705334099521468011318423160985507943967440420335063304346970852155833601505006094802527630580485975440894441077533686879443347828349415299434860114778542645481417844030011732317840312124857577317624806403294574388054782407793941193395139576232216820471248926353033028753066607315177945805410490322510890417889986155290813751020615454150372673790818813879164226820716569511273513116613734799963972530877463950648401174889762511716630

    """
    enumeration = 0
    n = graph.order()
    if quick:
        if isinstance(graph, Graph):
            for i, j, l in graph.edge_iterator():
                enumeration += 1 << ((n-(i+1))*n + n-(j+1))
                enumeration += 1 << ((n-(j+1))*n + n-(i+1))
        elif isinstance(graph, DiGraph):
            for i, j, l in graph.edge_iterator():
                enumeration += 1 << ((n-(i+1))*n + n-(j+1))
        return enumeration
    M = graph.am()
    for i, j in M.nonzero_positions():
        enumeration += 1 << ((n-(i+1))*n + n-(j+1))
    return ZZ(enumeration)

def paths_helper(start, end, G, all_paths, p=None):
    """
    The recursive helper for path finding calls.  (i.e.: all_paths
    and interior_paths).  Spawns potential path for each unvisited
    neighbor of current vertex and appends all succesful paths to
    one list.  (Note that paths themselves are lists of vertices).

    INPUT:
        start -- the vertex to start path search at
        end -- the vertex to find a path to
        all_paths -- the list (should initially be empty) to append
                     all successful paths to
        p -- the current path to update (via appending a vertex)
    """

    if p is None:
        # ONLY ONCE
        p = [start]

    plist = []
    # At each vertex, fill list of spawning paths (i.e. all neighbors)
    for i in range(len(G[p[-1]])):
        if G[p[-1]][i] not in p:
            plist.append(p + [G[p[-1]][i]])

    # If path completes, add to list
    if (p[-1] == end):
        all_paths.append(p)

    # Recursion: (look at all neighbors of all neighbors)
    for p in plist:
        paths_helper(start, end, G, all_paths, p)


