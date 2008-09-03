r"""
Graph Theory

This module implements many graph theoretic operations and concepts.

AUTHORS:
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
    -- Robert L. Miller (Sage Days 7): Edge labeled graph isomorphism
    -- Tom Boothby (Sage Days 7): Miscellaneous awesomeness
    -- Tom Boothby (2008-01-09): Added graphviz output

\subsection{Graph Format}

\subsubsection{The \sage Graph Class: NetworkX plus}

            \sage graphs are actually NetworkX graphs, wrapped in a \sage class.
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


            \sage Graphs can be created from a wide range of inputs. A few examples are
            covered here.

            \begin{itemize}

                \item NetworkX dictionary format:

                sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], \
                      5: [7, 8], 6: [8,9], 7: [9]}
                sage: G = Graph(d); G
                Graph on 10 vertices
                sage: G.plot().show()    # or G.show()

                \item A NetworkX graph:

                sage: K = networkx.complete_bipartite_graph(12,7)
                sage: G = Graph(K)
                sage: G.degree()
                [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 12, 12, 12, 12, 12, 12, 12]

                \item graph6 or sparse6 format:

                sage: s = ':I`AKGsaOs`cI]Gb~'
                sage: G = Graph(s, sparse=True); G
                Looped multi-graph on 10 vertices
                sage: G.plot().show()    # or G.show()

                Note that the \verb+\+ character is an escape character in Python, and
                also a character used by graph6 strings:

                sage: G = Graph('Ihe\n@GUA')
                Traceback (most recent call last):
                ...
                RuntimeError: The string (Ihe) seems corrupt: for n = 10, the string is too short.

                In Python, the escaped character \verb+\+ is represented by \verb+\\+:

                sage: G = Graph('Ihe\\n@GUA')
                sage: G.plot().show()    # or G.show()

                \item adjacency matrix: In an adjacency matrix, each column and each row represent
                a vertex. If a 1 shows up in row $i$, column $j$, there is an edge $(i,j)$.

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
                sage: G.plot().show()    # or G.show()

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
                sage: G.plot().show()    # or G.show()

        \end{itemize}

        \subsection{Generators}

        If you wish to iterate through all the isomorphism types of graphs,
        type, for example:

            sage: for g in graphs(4):
            ...     print g.spectrum()
            [0.0, 0.0, 0.0, 0.0]
            ...
            [-1.0, -1.0, -1.0, 3.0]

        For some commonly used graphs to play with, type

            sage: graphs.[tab]          # not tested

        and hit \kbd{tab}. Most of these graphs come with their own custom plot, so you
        can see how people usually visualize these graphs.

            sage: G = graphs.PetersenGraph()
            sage: G.plot().show()    # or G.show()
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
            sage: S.plot().show()    # or S.show()
            sage: S.density()
            1/2

            sage: G = GraphDatabase()
            sage: L = G.get_list(num_vertices=7, diameter=5)
            sage: graphs_list.show_graphs(L)

        \subsection{Labels}\label{Graph:labels}

        Each vertex can have any hashable object as a label. These are things like
        strings, numbers, and tuples. Each edge is given a default label of \code{None}, but
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
            sage: T.set_vertices(d)
            sage: T.get_vertex(1)
            Flower Snark: Graph on 20 vertices

        \subsection{Database}

        There is a database available for searching for graphs that satisfy a certain set
        of parameters, including number of vertices and edges, density, maximum and minimum
        degree, diameter, radius, and connectivity. If you wish to search a database of
        graphs by parameter, type

            sage: graphs_query.[tab]             # not tested

        and hit \kbd{tab}.

            sage: graphs_query = GraphDatabase()
            sage: L = graphs_query.get_list(num_vertices=7, diameter=5)
            sage: graphs_list.show_graphs(L)

        \subsection{Visualization}

        To see a graph $G$ you are working with, right now there are two main options.
        You can view the graph in two dimensions via matplotlib with \method{show()}.

            sage: G = graphs.RandomGNP(15,.3)
            sage: G.show()

        Or you can view it in three dimensions via jmol with \method{show3d()}.

            sage: G.show3d()

"""

#*****************************************************************************
#      Copyright (C) 2006 - 2007 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.prandom import random
from sage.structure.sage_object import SageObject
import sage.graphs.graph_fast as graph_fast
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational import Rational

class GenericGraph(SageObject):
    """
    Base class for graphs and digraphs.

    """

    # Nice defaults for plotting arrays of graphs (see sage.misc.functional.show)
    graphics_array_defaults =  {'layout': 'circular', 'vertex_size':50, 'vertex_labels':False, 'graph_border':True}

    def __add__(self, other_graph):
        """
        Returns the disjoint union of self and other.

        If there are common vertices to both, they will be renamed.

        EXAMPLE:
            sage: G = graphs.CycleGraph(3)
            sage: H = graphs.CycleGraph(4)
            sage: J = G + H; J
            Cycle graph disjoint_union Cycle graph: Graph on 7 vertices
            sage: J.vertices()
            [0, 1, 2, 3, 4, 5, 6]

        """
        if isinstance(other_graph, GenericGraph):
            return self.disjoint_union(other_graph, verbose_relabel=False)

    def __eq__(self, other):
        """
        Comparison of self and other. For equality, must be in the same class,
        have the same settings for loops and multiedges, output the same
        vertex list (in order) and the same adjacency matrix.

        Note that this is _not_ an isomorphism test.

        EXAMPLES:
            sage: G = graphs.EmptyGraph()
            sage: H = Graph()
            sage: G == H
            True
            sage: G.to_directed() == H.to_directed()
            True
            sage: G = graphs.RandomGNP(8,.9999)
            sage: H = graphs.CompleteGraph(8)
            sage: G == H # most often true
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
            sage: G == H # most often false
            False
            sage: G = Graph(multiedges=True)
            sage: G.add_edge(0,1)
            sage: H = G.copy()
            sage: H.add_edge(0,1)
            sage: G == H
            False

        """
        # inputs must be (di)graphs:
        if not isinstance(other, GenericGraph):
            raise TypeError("Cannot compare graph to non-graph (%s)."%str(other))
        g1_is_graph = isinstance(self, Graph)
        g2_is_graph = isinstance(other, Graph)
        if g1_is_graph != g2_is_graph:
            return False
        if self.multiple_edges() != other.multiple_edges():
            return False
        if self.loops() != other.loops():
            return False
        if self.vertices() != other.vertices():
            return False
        verts = self.vertices()
        # Finally, we are prepared to check edges:
        if not self.multiple_edges():
            for i in verts:
                for j in verts:
                    if self.has_edge(i,j) != other.has_edge(i,j):
                        return False
                    if self.has_edge(i,j) and self._weighted and other._weighted:
                        if self.edge_label(i,j) != other.edge_label(i,j):
                            return False
            return True
        else:
            for i in verts:
                for j in verts:
                    if self.has_edge(i, j):
                        edges1 = self.edge_label(i, j)
                    else:
                        edges1 = []
                    if other.has_edge(i, j):
                        edges2 = other.edge_label(i, j)
                    else:
                        edges2 = []
                    if len(edges1) != len(edges2):
                        return False
                    if sorted(edges1) != sorted(edges2) and self._weighted and other._weighted:
                        return False
            return True

    def __hash__(self):
        """
        Since graphs are mutable, they should not be hashable, so we return a type error.

        EXAMPLES:
            sage: hash(Graph())
            Traceback (most recent call last):
            ...
            TypeError: graphs are mutable, and thus not hashable

        """
        raise TypeError, "graphs are mutable, and thus not hashable"

    def __mul__(self, n):
        """
        Returns the sum of a graph with itself n times.

        EXAMPLE:
            sage: G = graphs.CycleGraph(3)
            sage: H = G*3; H
            Cycle graph disjoint_union Cycle graph disjoint_union Cycle graph: Graph on 9 vertices
            sage: H.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8]

        """
        if isinstance(n, (int, long, Integer)):
            if n < 1:
                raise TypeError('Multiplication of a graph and a nonpositive integer is not defined.')
            if n == 1:
                return self.copy()
            return sum([self]*(n-1), self)
        else:
            raise TypeError('Multiplication of a graph and something other than an integer is not defined.')

    def __ne__(self, other):
        """
        Tests for inequality, complement of __eq__.

        EXAMPLES:
            sage: g = Graph()
            sage: g2 = g.copy()
            sage: g == g
            True
            sage: g != g
            False
            sage: g2 == g
            True
            sage: g2 != g
            False
            sage: g is g
            True
            sage: g2 is g
            False

        """
        return (not (self == other))

    def __rmul__(self, n):
        """
        Returns the sum of a graph with itself n times.

        EXAMPLE:
            sage: G = graphs.CycleGraph(3)
            sage: H = int(3)*G; H
            Cycle graph disjoint_union Cycle graph disjoint_union Cycle graph: Graph on 9 vertices
            sage: H.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8]

        """
        return self*n

    def __str__(self):
        """
        str(G) returns the name of the graph, unless the name is the empty string, in
        which case it returns the default representation.

        EXAMPLE:
            sage: G = graphs.PetersenGraph()
            sage: str(G)
            'Petersen graph'

        """
        if self.name():
            return self.name()
        else:
            return repr(self)

    def _bit_vector(self):
        """
        Returns a string representing the edges of the (simple) graph for graph6 and dig6 strings.

        EXAMPLE:
            sage: G = graphs.PetersenGraph()
            sage: G._bit_vector()
            '101001100110000010000001001000010110000010110'
            sage: len([a for a in G._bit_vector() if a == '1'])
            15
            sage: G.num_edges()
            15

        """
        vertices = self.vertices()
        n = len(vertices)
        if self._directed:
            total_length = n*n
            bit = lambda x,y : x*n + y
        else:
            total_length = int(n*(n - 1))/int(2)
            n_ch_2 = lambda b : int(b*(b-1))/int(2)
            bit = lambda x,y : n_ch_2(max([x,y])) + min([x,y])
        bit_vector = set()
        for u,v,_ in self.edge_iterator():
            bit_vector.add(bit(vertices.index(u), vertices.index(v)))
        bit_vector = sorted(bit_vector)
        s = []
        j = 0
        for i in bit_vector:
            s.append( '0'*(i - j) + '1' )
            j = i + 1
        s = "".join(s)
        s += '0'*(total_length-len(s))
        return s

    def _latex_(self):
        """
        To include a graph in LaTeX document, see function
        Graph.write_to_eps().

        EXAMPLE:
            sage: G = graphs.PetersenGraph()
            sage: latex(G)
            Traceback (most recent call last):
            ...
            NotImplementedError: To include a graph in LaTeX document, see function Graph.write_to_eps().
            sage: G._latex_()
            Traceback (most recent call last):
            ...
            NotImplementedError: To include a graph in LaTeX document, see function Graph.write_to_eps().

        """
        raise NotImplementedError('To include a graph in LaTeX document, see function Graph.write_to_eps().')

    def _matrix_(self, R=None):
        """
        Returns the adjacency matrix of the graph over the specified ring.

        EXAMPLES:
            sage: G = graphs.CompleteBipartiteGraph(2,3)
            sage: m = matrix(G); m.parent()
            Full MatrixSpace of 5 by 5 dense matrices over Integer Ring
            sage: m
            [0 0 1 1 1]
            [0 0 1 1 1]
            [1 1 0 0 0]
            [1 1 0 0 0]
            [1 1 0 0 0]
            sage: G._matrix_()
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

    def _repr_(self):
        """
        Return a string representation of self.

        EXAMPLE:
            sage: G = graphs.PetersenGraph()
            sage: G._repr_()
            'Petersen graph: Graph on 10 vertices'

        """
        name = ""
        if self.loops():
            name += "looped "
        if self.multiple_edges():
            name += "multi-"
        if self._directed:
            name += "di"
        name += "graph on %d vert"%self.order()
        if self.order() == 1:
            name += "ex"
        else:
            name += "ices"
        name = name.capitalize()
        if self.name() != '':
            name = self.name() + ": " + name
        return name

    ### Formats

    def copy(self, implementation='networkx', sparse=True):
        """
        Creates a copy of the graph.

        EXAMPLES:
            sage: g=Graph({0:[0,1,1,2]},loops=True,multiedges=True,implementation='networkx')
            sage: g==g.copy()
            True
            sage: g=DiGraph({0:[0,1,1,2],1:[0,1]},loops=True,multiedges=True,implementation='networkx')
            sage: g==g.copy()
            True

        Note that vertex associations are also kept:
            sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), 2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
            sage: T = graphs.TetrahedralGraph()
            sage: T.set_vertices(d)
            sage: T2 = T.copy()
            sage: T2.get_vertex(0)
            Dodecahedron: Graph on 20 vertices

        Notice that the copy is at least as deep as the objects:
            sage: T2.get_vertex(0) is T.get_vertex(0)
            False

        TESTS:
        We make copies of the _pos and _boundary attributes.
            sage: g = graphs.PathGraph(3)
            sage: h = g.copy()
            sage: h._pos is g._pos
            False
            sage: h._boundary is g._boundary
            False

        """
        from copy import copy
        if self._directed:
            G = DiGraph(self, name=self.name(), pos=copy(self._pos), boundary=copy(self._boundary), implementation=implementation, sparse=sparse)
        else:
            G = Graph(self, name=self.name(), pos=copy(self._pos), boundary=copy(self._boundary), implementation=implementation, sparse=sparse)
        if hasattr(self, '_assoc'):
            G._assoc = {}
            for v in self._assoc:
                try:
                    G._assoc[v] = self._assoc[v].copy()
                except:
                    G._assoc[v] = copy(self._assoc[v])
        if hasattr(self, '_embedding'):
            G._embedding = self._embedding
        G._weighted = self._weighted
        return G

    def networkx_graph(self, copy=True):
        """
        Creates a new NetworkX graph from the SAGE graph.

        INPUT:
        copy -- if False, and the underlying implementation is a NetworkX graph,
            then the actual object itself is returned.

        EXAMPLES:
            sage: G = graphs.TetrahedralGraph()
            sage: N = G.networkx_graph()
            sage: type(N)
            <class 'networkx.xgraph.XGraph'>

            sage: G = graphs.TetrahedralGraph()
            sage: G = Graph(G, implementation='networkx')
            sage: N = G.networkx_graph()
            sage: G._backend._nxg is N
            False

            sage: G = Graph(graphs.TetrahedralGraph(), implementation='networkx')
            sage: N = G.networkx_graph(copy=False)
            sage: G._backend._nxg is N
            True

        """
        try:
            if copy:
                return self._backend._nxg.copy()
            else:
                return self._backend._nxg
        except:
            import networkx
            if self._directed:
                class_type = networkx.XDiGraph
            else:
                class_type = networkx.XGraph
            N = class_type(selfloops=self.loops(), multiedges=self.multiple_edges(),
                           name=self.name())
            N.add_nodes_from(self.vertices())
            N.add_edges_from(self.edges())
            return N

    def adjacency_matrix(self, sparse=None, boundary_first=False):
        """
        Returns the adjacency matrix of the (di)graph. Each vertex is
        represented by its position in the list returned by the vertices()
        function.

        The matrix returned is over the integers.  If a different ring
        is desired, use either the change_ring function or the matrix
        function.

        INPUT:
            sparse -- whether to represent with a sparse matrix
            boundary_first -- whether to represent the boundary vertices in
                the upper left block

        EXAMPLES:
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

            sage: matrix(GF(2),G) # matrix over GF(2)
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

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.adjacency_matrix()
            [0 1 1 1 0 0]
            [1 0 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [1 0 0 0 0 1]
            [0 1 0 0 0 0]

        TESTS:
            sage: graphs.CubeGraph(8).adjacency_matrix().parent()
            Full MatrixSpace of 256 by 256 dense matrices over Integer Ring
            sage: graphs.CubeGraph(9).adjacency_matrix().parent()
            Full MatrixSpace of 512 by 512 sparse matrices over Integer Ring


        """
        n = self.order()
        if sparse is None:
            if n <= 256 or self.density() > 0.05:
                sparse=False
            else:
                sparse=True

        verts = self.vertices(boundary_first=boundary_first)
        D = {}
        directed = self._directed
        multiple_edges = self.multiple_edges()
        for i,j,l in self.edge_iterator():
            i = verts.index(i)
            j = verts.index(j)
            if multiple_edges and (i,j) in D:
                D[(i,j)] += 1
                if not directed:
                    D[(j,i)] += 1
            else:
                D[(i,j)] = 1
                if not directed:
                    D[(j,i)] = 1
        from sage.rings.integer_ring import IntegerRing
        from sage.matrix.constructor import matrix
        M = matrix(IntegerRing(), n, n, D, sparse=sparse)
        return M

    am = adjacency_matrix # shorter call makes life easier

    def incidence_matrix(self, sparse=True):
        """
        Returns an incidence matrix of the (di)graph. Each row is a vertex, and
        each column is an edge. Note that in the case of graphs, there is a
        choice of orientation for each edge.

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
        n = self.order()
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
            sage: G = Graph(sparse=True)
            sage: G.add_edges([(0,1,1),(1,2,2),(0,2,3),(0,3,4)])
            sage: M = G.weighted_adjacency_matrix(); M
            [0 1 3 4]
            [1 0 2 0]
            [3 2 0 0]
            [4 0 0 0]
            sage: H = Graph(data=M, format='weighted_adjacency_matrix', sparse=True)
            sage: H == G
            True

        """
        if self.multiple_edges():
            raise NotImplementedError, "Don't know how to represent weights for a multigraph."

        verts = self.vertices(boundary_first=boundary_first)
        D = {}
        if self._directed:
            for i,j,l in self.edge_iterator():
                i = verts.index(i)
                j = verts.index(j)
                D[(i,j)] = l
        else:
            for i,j,l in self.edge_iterator():
                i = verts.index(i)
                j = verts.index(j)
                D[(i,j)] = l
                D[(j,i)] = l
        from sage.matrix.constructor import matrix
        M = matrix(D, sparse=sparse)
        return M

    def kirchhoff_matrix(self, weighted=None, boundary_first=False):
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
            sage: G = Graph(sparse=True)
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
            sage: M = G.laplacian_matrix(boundary_first=True); M
            [ 2  0 -1 -1]
            [ 0  1 -1  0]
            [-1 -1  3 -1]
            [-1  0 -1  2]

        """
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import IntegerRing

        if weighted is None:
            weighted = self._weighted

        if weighted:
            M = self.weighted_adjacency_matrix(boundary_first=boundary_first)
        else:
            M = self.adjacency_matrix(boundary_first=boundary_first)
        A = list(-M)
        S = [sum(M.row(i)) for i in range(M.nrows())]
        for i in range(len(A)):
            A[i][i] = S[i]
        return M.parent()(A)

    laplacian_matrix = kirchhoff_matrix

    ### Attributes

    def get_boundary(self):
        """
        Returns the boundary of the (di)graph.

        EXAMPLE:
            sage: G = graphs.PetersenGraph()
            sage: G.set_boundary([0,1,2,3,4])
            sage: G.get_boundary()
            [0, 1, 2, 3, 4]

        """
        return self._boundary

    def set_boundary(self, boundary):
        """
        Sets the boundary of the (di)graph.

        EXAMPLE:
            sage: G = graphs.PetersenGraph()
            sage: G.set_boundary([0,1,2,3,4])
            sage: G.get_boundary()
            [0, 1, 2, 3, 4]
            sage: G.set_boundary((1..4))
            sage: G.get_boundary()
            [1, 2, 3, 4]
        """
        if isinstance(boundary,list):
            self._boundary = boundary
        else:
            self._boundary = list(boundary)

    def set_embedding(self, embedding):
        """
        Sets a combinatorial embedding dictionary to _embedding attribute.
        Dictionary is organized with vertex labels as keys and a list of
        each vertex's neighbors in clockwise order.

        Dictionary is error-checked for validity.

        INPUT:
            embedding -- a dictionary

        EXAMPLES:
            sage: G = graphs.PetersenGraph()
            sage: G.set_embedding({0: [1, 5, 4], 1: [0, 2, 6], 2: [1, 3, 7], 3: [8, 2, 4], 4: [0, 9, 3], 5: [0, 8, 7], 6: [8, 1, 9], 7: [9, 2, 5], 8: [3, 5, 6], 9: [4, 6, 7]})
            sage: G.set_embedding({'s': [1, 5, 4], 1: [0, 2, 6], 2: [1, 3, 7], 3: [8, 2, 4], 4: [0, 9, 3], 5: [0, 8, 7], 6: [8, 1, 9], 7: [9, 2, 5], 8: [3, 5, 6], 9: [4, 6, 7]})
            Traceback (most recent call last):
            ...
            Exception: embedding is not valid for Petersen graph

        """
        if self.check_embedding_validity(embedding):
            self._embedding = embedding
        else:
            raise Exception('embedding is not valid for %s'%self)

    def get_embedding(self):
        """
        Returns the attribute _embedding if it exists.  _embedding
        is a dictionary organized with vertex labels as keys and a list of
        each vertex's neighbors in clockwise order.

        Error-checked to insure valid embedding is returned.

        EXAMPLES:
            sage: G = graphs.PetersenGraph()
            sage: G.genus()
            1
            sage: G.get_embedding()
            {0: [1, 5, 4],
             1: [0, 2, 6],
             2: [1, 3, 7],
             3: [8, 2, 4],
             4: [0, 9, 3],
             5: [0, 8, 7],
             6: [8, 1, 9],
             7: [9, 2, 5],
             8: [3, 5, 6],
             9: [4, 6, 7]}

        """
        if self.check_embedding_validity():
            return self._embedding
        else:
            raise Exception('%s has been modified and the embedding is no longer valid.'%self)

    def check_embedding_validity(self, embedding=None):
        """
        Checks whether an _embedding attribute is defined on self and if so,
        checks for accuracy.  Returns True if everything is okay, False otherwise.

        If embedding=None will test the attribute _embedding.

        EXAMPLES:
            sage: d = {0: [1, 5, 4], 1: [0, 2, 6], 2: [1, 3, 7], 3: [8, 2, 4], 4: [0, 9, 3], 5: [0, 8, 7], 6: [8, 1, 9], 7: [9, 2, 5], 8: [3, 5, 6], 9: [4, 6, 7]}
            sage: G = graphs.PetersenGraph()
            sage: G.check_embedding_validity(d)
            True

        """
        if embedding is None:
            embedding = getattr(self, '_embedding', None)
        if embedding is None:
            return False
        if len(embedding) != self.order():
            return False
        if self._directed:
            connected = lambda u,v : self.has_edge(u,v) or self.has_edge(v,u)
        else:
            connected = lambda u,v : self.has_edge(u,v)
        for v in embedding:
            if not self.has_vertex(v):
                return False
            if len(embedding[v]) != len(self.neighbors(v)):
                return False
            for u in embedding[v]:
                if not connected(v,u):
                    return False
        return True

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
        if new is False:
            self.remove_loops()
        return self._backend.loops(new)

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
        return self._backend.multiple_edges(new)

    def name(self, new=None):
        """
        Returns the name of the (di)graph.

        INPUT:
        new -- if not None, then this becomes the new name of the (di)graph.
        (if new == '', removes any name)

        EXAMPLE:
            sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], 5: [7, 8], 6: [8,9], 7: [9]}
            sage: G = Graph(d); G
            Graph on 10 vertices
            sage: G.name("Petersen Graph"); G
            Petersen Graph: Graph on 10 vertices
            sage: G.name(new=""); G
            Graph on 10 vertices
            sage: G.name()
            ''

        """
        return self._backend.name(new)

    def get_pos(self):
        """
        Returns the position dictionary, a dictionary specifying the coordinates
        of each vertex.

        EXAMPLES:
        By default, the position of a graph is None:
            sage: G = Graph()
            sage: G.get_pos()
            sage: G.get_pos() is None
            True
            sage: P = G.plot(save_pos=True)
            sage: G.get_pos()
            {}

        Some of the named graphs come with a pre-specified positioning:
            sage: G = graphs.PetersenGraph()
            sage: G.get_pos()
            {0: [..., ...],
             ...
             9: [..., ...]}

        """
        return self._pos

    def set_pos(self, pos):
        """
        Sets the position dictionary, a dictionary specifying the
        coordinates of each vertex.

        EXAMPLE:
        Note that set_pos will allow you to do ridiculous things, which
        will not blow up until plotting:
            sage: G = graphs.PetersenGraph()
            sage: G.get_pos()
            {0: [..., ...],
             ...
             9: [..., ...]}

            sage: G.set_pos('spam')
            sage: P = G.plot()
            Traceback (most recent call last):
            ...
            TypeError: string indices must be integers

        """
        self._pos = pos

    def weighted(self, new=None):
        """
        Returns whether the (di)graph is to be considered as a weighted (di)graph.

        Note that edge weightings can still exist for (di)graphs G where
        G.weighted() is False.

        EXAMPLE:
        Here we have two graphs with different labels, but weighted is False
        for both, so we just check for the presence of edges:
            sage: G = Graph({0:{1:'a'}}, implementation='networkx')
            sage: H = Graph({0:{1:'b'}}, implementation='networkx')
            sage: G == H
            True

        Now one is weighted and the other is not, so the comparison is done as
        if neither is weighted:
            sage: G.weighted(True)
            sage: H.weighted()
            False
            sage: G == H
            True

        However, if both are weighted, then we finally compare 'a' to 'b'.
            sage: H.weighted(True)
            sage: G == H
            False

        """
        if new is not None:
            if new:
                self._weighted = new
        else:
            return self._weighted

    ### Properties

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
        sage: digraphs.RandomDirectedGNR(20,0.5).antisymmetric()
        True

        """
        if not self._directed:
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
            if self._directed:
                return Rational(self.size())/Rational(n**2)
            else:
                return Rational(self.size())/Rational((n**2 + n)/2)
        else:
            if self._directed:
                return Rational(self.size())/Rational((n**2 - n))
            else:
                return Rational(self.size())/Rational((n**2 - n)/2)

    def is_eulerian(self):
        """
        Return true if the graph has an tour that visits each edge exactly once.

        EXAMPLES:
            sage: graphs.CompleteGraph(4).is_eulerian()
            False
            sage: graphs.CycleGraph(4).is_eulerian()
            True
            sage: g = DiGraph({0:[1,2], 1:[2]}); g.is_eulerian()
            False
            sage: g = DiGraph({0:[2], 1:[3], 2:[0,1], 3:[2]}); g.is_eulerian()
            True

        """
        if not self.is_connected():
            return False
        if self._directed:
            for i in self.vertex_iterator():
                # loops don't matter since they count in both the in and out degree.
                if self.in_degree(i) != self.out_degree(i):
                    return False
        else:
            for i in self.degree_iterator():
                # loops don't matter since they add an even number to the degree
                if i % 2 != 0:
                    return False
        return True

    def is_tree(self):
        """
        Return True if the graph is a tree.

        EXAMPLES:
            sage: for g in graphs.trees(6):
            ...     g.is_tree()
            True
            True
            True
            True
            True
            True

        """
        if not self.is_connected():
            return False
        if self.num_verts() != self.num_edges() + 1:
            return False
        return True

    def is_forest(self):
        """
        Return True if the graph is a forest, i.e. a disjoint union of trees.

        EXAMPLE:
            sage: seven_acre_wood = sum(graphs.trees(7), Graph())
            sage: seven_acre_wood.is_forest()
            True

        """
        for g in self.connected_components_subgraphs():
            if not g.is_tree():
                return False
        return True

    def order(self):
        """
        Returns the number of vertices. Note that len(G) returns the number of
        vertices in G also.

        EXAMPLE:
            sage: G = graphs.PetersenGraph()
            sage: G.order()
            10

            sage: G = graphs.TetrahedralGraph()
            sage: len(G)
            4

        """
        return self._backend.num_verts()

    __len__ = order

    num_verts = order

    def size(self):
        """
        Returns the number of edges.

        EXAMPLE:
            sage: G = graphs.PetersenGraph()
            sage: G.size()
            15

        """
        return self._backend.num_edges(self._directed)

    num_edges = size

    ### Planarity

    def is_planar(self, on_embedding=None, kuratowski=False, set_embedding=False, set_pos=False):
        """
        Returns True if the graph is planar, and False otherwise.  This wraps the
        reference implementation provided by John Boyer of the linear time
        planarity algorithm by edge addition due to Boyer Myrvold.  (See reference
        code in graphs.planarity).

        Note -- the argument on_embedding takes precedence over set_embedding.
                This means that only the on_embedding combinatorial embedding
                will be tested for planarity and no _embedding attribute will
                be set as a result of this function call, unless on_embedding
                is None.

        REFERENCE:
            [1] John M. Boyer and Wendy J. Myrvold, On the Cutting Edge:
                Simplified O(n) Planarity by Edge Addition.  Journal of Graph
                Algorithms and Applications, Vol. 8, No. 3, pp. 241-273, 2004.

        INPUT:
            kuratowski -- returns a tuple with boolean as first entry.  If the
                       graph is nonplanar, will return the Kuratowski subgraph
                       or minor as the second tuple entry.  If the graph is
                       planar, returns None as the second entry.
            on_embedding -- the embedding dictionary to test planarity on.  (i.e.:
                       will return True or False only for the given embedding.)
            set_embedding -- whether or not to set the instance field variable that
                       contains a combinatorial embedding (clockwise ordering
                       of neighbors at each vertex).  This value will only be
                       set if a planar embedding is found.  It is stored as a
                       Python dict: {v1: [n1,n2,n3]} where v1 is a vertex and
                       n1,n2,n3 are its neighbors.
            set_pos -- whether or not to set the position dictionary (for
                       plotting) to reflect the combinatorial embedding.  Note
                       that this value will default to False if set_emb is set
                       to False.  Also, the position dictionary will only be
                       updated if a planar embedding is found.

        EXAMPLES:
            sage: g = graphs.CubeGraph(4)
            sage: g.is_planar()
            False

            sage: g = graphs.CircularLadderGraph(4)
            sage: g.is_planar(set_embedding=True)
            True
            sage: g.get_embedding()
            {0: [1, 4, 3],
             1: [2, 5, 0],
             2: [3, 6, 1],
             3: [0, 7, 2],
             4: [0, 5, 7],
             5: [1, 6, 4],
             6: [2, 7, 5],
             7: [4, 6, 3]}

            sage: g = graphs.PetersenGraph()
            sage: (g.is_planar(kuratowski=True))[1].adjacency_matrix()
            [0 1 0 0 0 1 0 0 0 0]
            [1 0 1 0 0 0 1 0 0 0]
            [0 1 0 1 0 0 0 0 0 0]
            [0 0 1 0 1 0 0 0 1 0]
            [0 0 0 1 0 0 0 0 0 1]
            [1 0 0 0 0 0 0 1 1 0]
            [0 1 0 0 0 0 0 0 1 1]
            [0 0 0 0 0 1 0 0 0 1]
            [0 0 0 1 0 1 1 0 0 0]
            [0 0 0 0 1 0 1 1 0 0]

            sage: k43 = graphs.CompleteBipartiteGraph(4,3)
            sage: result = k43.is_planar(kuratowski=True); result
            (False, Graph on 6 vertices)
            sage: result[1].is_isomorphic(graphs.CompleteBipartiteGraph(3,3))
            True
        """
        if on_embedding:
            if self.check_embedding_validity(on_embedding):
                return (0 == self.genus(minimal=False,set_embedding=False,on_embedding=on_embedding))
            else:
                raise Exception('on_embedding is not a valid embedding for %s.'%self)
        else:
            from sage.graphs.planarity import is_planar
            G = self.to_undirected()
            planar = is_planar(G,kuratowski=kuratowski,set_pos=set_pos,set_embedding=set_embedding)
            if kuratowski:
                bool_result = planar[0]
            else:
                bool_result = planar
            if bool_result:
                if set_pos:
                    self._pos = G._pos
                if set_embedding:
                    self._embedding = G._embedding
            return planar

    def is_circular_planar(self, ordered=True, kuratowski=False, set_embedding=False, set_pos=False):
        """
        A graph (with nonempty boundary) is circular planar if it has a
        planar embedding in which all boundary vertices can be drawn in
        order on a disc boundary, with all the interior vertices drawn
        inside the disc.

        Returns True if the graph is circular planar, and False if it is
        not.  If kuratowski is set to True, then this function will return
        a tuple, with boolean first entry and second entry the Kuratowski
        subgraph or minor isolated by the Boyer-Myrvold algorithm.
        Note that this graph might contain a vertex or edges that
        were not in the initial graph.  These would be elements referred to
        below as parts of the wheel and the star, which were added to the
        graph to require that the boundary can be drawn on the boundary of a
        disc, with all other vertices drawn inside (and no edge crossings).
        For more information, refer to reference [2].

        This is a linear time algorithm to test for circular planarity.  It
        relies on the edge-addition planarity algorithm due to Boyer-Myrvold.
        We accomplish linear time for circular planarity by modifying the graph
        before running the general planarity algorithm.

        REFERENCE:
            [1] John M. Boyer and Wendy J. Myrvold, On the Cutting Edge:
                Simplified O(n) Planarity by Edge Addition.  Journal of Graph
                Algorithms and Applications, Vol. 8, No. 3, pp. 241-273, 2004.
            [2] Kirkman, Emily A. O(n) Circular Planarity Testing. [Online]
                Available: soon!

        INPUT:
            ordered -- whether or not to consider the order of the boundary
                       (set ordered=False to see if there is any possible
                       boundary order that will satisfy circular planarity)
            kuratowski -- if set to True, returns a tuple with boolean first
                       entry and the Kuratowski subgraph or minor as the second
                       entry.  See notes above.
            on_embedding -- the embedding dictionary to test planarity on.  (i.e.:
                       will return True or False only for the given embedding.)
            set_embedding -- whether or not to set the instance field variable that
                       contains a combinatorial embedding (clockwise ordering
                       of neighbors at each vertex).  This value will only be
                       set if a circular planar embedding is found.  It is
                       stored as a Python dict: {v1: [n1,n2,n3]} where v1 is
                       a vertex and n1,n2,n3 are its neighbors.
            set_pos -- whether or not to set the position dictionary (for
                       plotting) to reflect the combinatorial embedding.  Note
                       that this value will default to False if set_emb is set
                       to False.  Also, the position dictionary will only be
                       updated if a circular planar embedding is found.

        EXAMPLES:
            sage: g439 = Graph({1:[5,7], 2:[5,6], 3:[6,7], 4:[5,6,7]}, implementation='networkx')
            sage: g439.set_boundary([1,2,3,4])
            sage: g439.show(figsize=[2,2], vertex_labels=True, vertex_size=175)
            sage: g439.is_circular_planar()
            False
            sage: g439.is_circular_planar(kuratowski=True)
            (False, Graph on 7 vertices)
            sage: g439.set_boundary([1,2,3])
            sage: g439.is_circular_planar(set_embedding=True, set_pos=False)
            True
            sage: g439.is_circular_planar(kuratowski=True)
            (True, None)
            sage: g439.get_embedding()
            {1: [7, 5],
             2: [5, 6],
             3: [6, 7],
             4: [7, 6, 5],
             5: [4, 2, 1],
             6: [4, 3, 2],
             7: [3, 4, 1]}

        Order matters:
            sage: K23 = graphs.CompleteBipartiteGraph(2,3)
            sage: K23.set_boundary([0,1,2,3])
            sage: K23.is_circular_planar()
            False
            sage: K23.is_circular_planar(ordered=False)
            True
            sage: K23.set_boundary([0,2,1,3]) # Diff Order!
            sage: K23.is_circular_planar(set_embedding=True)
            True
        """
        from sage.graphs.planarity import is_planar
        graph = self.to_undirected()
        boundary = self.get_boundary()
        if hasattr(graph, '_embedding'):
            del(graph._embedding)

        extra = 0
        while graph.has_vertex(extra):
            extra=extra+1
        graph.add_vertex(extra)

        for vertex in boundary:
            graph.add_edge(vertex,extra)

        extra_edges = []
        if ordered: # WHEEL
            for i in range(len(boundary)-1):
                if not graph.has_edge(boundary[i],boundary[i+1]):
                    graph.add_edge(boundary[i],boundary[i+1])
                    extra_edges.append((boundary[i],boundary[i+1]))
            if not graph.has_edge(boundary[-1],boundary[0]):
                graph.add_edge(boundary[-1],boundary[0])
                extra_edges.append((boundary[-1],boundary[0]))
        # else STAR (empty list of extra edges)

        result = is_planar(graph,kuratowski=kuratowski,set_embedding=set_embedding,circular=True)

        if kuratowski:
            bool_result = result[0]
        else:
            bool_result = result

        if bool_result:
            graph.delete_vertex(extra)
            graph.delete_edges(extra_edges)

            if hasattr(graph,'_embedding'):
                # strip the embedding to fit original graph
                graph._embedding.pop(extra)
                for u,v in extra_edges:
                    graph._embedding[u].pop(graph._embedding[u].index(v))
                    graph._embedding[v].pop(graph._embedding[v].index(u))
                for w in boundary:
                    graph._embedding[w].pop(graph._embedding[w].index(extra))

                if set_embedding:
                    self._embedding = graph._embedding

            if (set_pos and set_embedding):
                self.set_planar_positions()
        return result

    def set_planar_positions(self, set_embedding=False, on_embedding=None, external_face=None, test=False, circular=False):
        """
        Uses Schnyder's algorithm to determine positions for a planar embedding of
        self, raising an error if self is not planar.

        INPUT:
            set_embedding -- if True, sets the combinatorial embedding used (see
        self.get_embedding())
            on_embedding -- dict: provide a combinatorial embedding
            external_face -- ignored
            test -- if True, perform sanity tests along the way
            circular -- ignored

        EXAMPLES:
            sage: g = graphs.PathGraph(10)
            sage: g.set_planar_positions(test=True)
            True
            sage: g = graphs.BalancedTree(3,4)
            sage: g.set_planar_positions(test=True)
            True
            sage: g = graphs.CycleGraph(7)
            sage: g.set_planar_positions(test=True)
            True
            sage: g = graphs.CompleteGraph(5)
            sage: g.set_planar_positions(test=True,set_embedding=True)
            Traceback (most recent call last):
            ...
            Exception: Complete graph is not a planar graph.

        """
        from sage.graphs.schnyder import _triangulate, _normal_label, _realizer, _compute_coordinates

        G = self.to_undirected()
        try:
            G._embedding = self._embedding
        except AttributeError:
            pass
        embedding_copy = None
        if set_embedding:
            if not (G.is_planar(set_embedding=True)):
                raise Exception('%s is not a planar graph.'%self)
            embedding_copy = G._embedding
        else:
            if on_embedding is not None:
                if G.check_embedding_validity(on_embedding):
                    if not (G.is_planar(on_embedding=on_embedding)):
                        raise Exception( 'Provided embedding is not a planar embedding for %s.'%self )
                else:
                    raise Exception('Provided embedding is not a valid embedding for %s. Try putting set_embedding=True.'%self)
            else:
                if hasattr(G,'_embedding'):
                    if G.check_embedding_validity():
                        if not (G.is_planar(on_embedding=G._embedding)):
                            raise Exception('%s has nonplanar _embedding attribute.  Try putting set_embedding=True.'%self)
                        embedding_copy = G._embedding
                    else:
                        raise Exception('Provided embedding is not a valid embedding for %s. Try putting set_embedding=True.'%self)
                else:
                    G.is_planar(set_embedding=True)

        # The following is what was breaking the code.  It is where we were specifying the external
        #       face ahead of time.  This is definitely a TODO:
        #
        # Running is_planar(set_embedding=True) has set attribute self._embedding
        #if external_face is None:
        #    faces = trace_faces( self, self._embedding )
        #    faces.sort(key=len)
        #    external_face = faces[-1]

        #n = len(external_face)
        #other_added_edges = []
        #if n > 3:
        #    v1, v2, v3 = external_face[0][0], external_face[int(n/3)][0], external_face[int(2*n/3)][0]
        #    if not self.has_edge( (v1,v2) ):
        #        self.add_edge( (v1, v2) )
        #        other_added_edges.append( (v1, v2) )
        #    if not self.has_edge( (v2,v3) ):
        #        self.add_edge( (v2, v3) )
        #        other_added_edges.append( (v2, v3) )
        #    if not self.has_edge( (v3,v1) ):
        #        self.add_edge( (v3, v1) )
        #        other_added_edges.append( (v3, v1) )
        #    if not self.is_planar(set_embedding=True): # get new combinatorial embedding (with added edges)
        #        raise Exception('Modified graph %s is not planar.  Try specifying an external face.'%self)

        # Triangulate the graph
        extra_edges = _triangulate( G, G._embedding)

        # Optional error-checking
        if test:
            G.is_planar(set_embedding=True) # to get new embedding
            test_faces = G.trace_faces(G._embedding)
            for face in test_faces:
                if len(face) != 3:
                    raise Exception('BUG: Triangulation returned face: %'%face)

        G.is_planar(set_embedding=True)
        faces = G.trace_faces(G._embedding)
        # Assign a normal label to the graph
        label = _normal_label( G, G._embedding, faces[0] )

        # Get dictionary of tree nodes from the realizer
        tree_nodes = _realizer( G, label)

        # Compute the coordinates and store in position dictionary (attr self._pos)
        _compute_coordinates( G, tree_nodes )
        self._pos = G._pos

        # Delete all the edges added to the graph
        #G.delete_edges( extra_edges )
        #self.delete_edges( other_added_edges )

        if embedding_copy is not None:
            self._embedding = embedding_copy

        if test:    # Optional error-checking, ( looking for edge-crossings O(n^2) ).
            return self.is_drawn_free_of_edge_crossings() # returns true if tests pass
        else:
            return

    def is_drawn_free_of_edge_crossings(self):
        """
        Returns True is the position dictionary for this graph is set
        and that position dictionary gives a planar embedding.

        This simply checks all pairs of edges that don't share a vertex
        to make sure that they don't intersect.

        NOTE:  This function require that _pos attribute is set.  (Returns
        false otherwise.)

        EXAMPLES:
            sage: D = graphs.DodecahedralGraph()
            sage: D.set_planar_positions()
            sage: D.is_drawn_free_of_edge_crossings()
            True

        """
        if self._pos is None:
            return False

        G = self.to_undirected()
        for edge1 in G.edges(labels = False):
            for edge2 in G.edges(labels = False):
                if edge1[0] == edge2[0] or edge1[0] == edge2[1] or edge1[1] == edge2[0] or edge1[1] == edge2[1]:
                    continue
                p1, p2 = self._pos[edge1[0]], self._pos[edge1[1]]
                dy = Rational(p2[1] - p1[1])
                dx = Rational(p2[0] - p1[0])
                q1, q2 = self._pos[edge2[0]], self._pos[edge2[1]]
                db = Rational(q2[1] - q1[1])
                da = Rational(q2[0] - q1[0])
                if(da * dy == db * dx):
                    if dx != 0:
                        t1 = Rational(q1[0] - p1[0])/dx
                        t2 = Rational(q2[0] - p1[0])/dx
                        if (0 <= t1 and t1 <= 1) or (0 <= t2 and t2 <= 1):
                            if p1[1] + t1 * dy == q1[1] or p1[1] + t2 * dy == q2[1]:
                                return False
                    else:
                        t1 = Rational(q1[1] - p1[1])/dy
                        t2 = Rational(q2[1] - p1[1])/dy
                        if (0 <= t1 and t1 <= 1) or (0 <= t2 and t2 <= 1):
                            if p1[0] + t1 * dx == q1[0] or p1[0] + t2 * dx == q2[0]:
                                return False
                else:
                    s = (dx * Rational(q1[1] - p1[1]) + dy * Rational(p1[0] - q1[0])) / (da * dy - db * dx)
                    t = (da * Rational(p1[1] - q1[1]) + db * Rational(q1[0] - p1[0])) / (db * dx - da * dy)

                    if s >= 0 and s <= 1 and t >= 0 and t <= 1:
                        print 'fail on', p1, p2, ' : ',q1, q2
                        print edge1, edge2
                        return False
        return True

    def genus(self, set_embedding=True, on_embedding=None, minimal=True, maximal=False, circular=False, ordered=True):
        """
        Returns the minimal genus of the graph.  The genus of a compact
        surface is the number of handles it has.  The genus of a graph
        is the minimal genus of the surface it can be embedded into.

        Note -- This function uses Euler's formula and thus it is
                necessary to consider only connected graphs.

        INPUT:
            set_embedding (boolean) -- whether or not to store an embedding attribute of the computed
                                       (minimal) genus of the graph.  (Default is True).
            on_embedding (dict) -- a combinatorial embedding to compute the genus of the graph on.
                                       Note that this must be a valid embedding for the graph.  The
                                       dictionary structure is given by:
                                       { vertex1: [neighbor1, neighbor2, neighbor3], vertex2: [neighbor] }
                                       where there is a key for each vertex in the graph and a (clockwise)
                                       ordered list of each vertex's neighbors as values.  on_embedding
                                       takes precedence over a stored _embedding attribute if minimal is
                                       set to False.
                                       Note that as a shortcut, the user can enter on_embedding=True to
                                       compute the genus on the current _embedding attribute.  (see eg's.)
            minimal (boolean) -- whether or not to compute the minimal genus of the graph (i.e., testing
                                       all embeddings).  If minimal is False, then either maximal must be
                                       True or on_embedding must not be None.  If on_embedding is not None,
                                       it will take priority over minimal.  Similarly, if maximal is True,
                                       it will take priority over minimal.
            maximal (boolean) -- whether or not to compute the maximal genus of the graph (i.e., testin
                                       all embeddings).  If maximal is False, then either minimal must be
                                       True or on_embedding must not be None.  If on_embedding is not None,
                                       it will take priority over maximal.  However, maximal takes priority
                                       over the default minimal.
            circular (boolean) -- whether or not to compute the genus preserving a planar embedding of the
                                       boundary.  (Default is False).  If circular is True, on_embedding is
                                       not a valid option.
            ordered (boolean) -- if circular is True, then whether or not the boundary order may be permuted.
                                       (Default is True, which means the boundary order is preserved.)

        EXAMPLES:
            sage: g = graphs.PetersenGraph()
            sage: g.genus() # tests for minimal genus by default
            1
            sage: g.genus(on_embedding=True, maximal=True) # on_embedding overrides minimal and maximal arguments
            1
            sage: g.genus(maximal=True) # setting maximal to True overrides default minimal=True
            3
            sage: g.genus(on_embedding=g.get_embedding()) # can also send a valid combinatorial embedding dict
            3
            sage: (graphs.CubeGraph(3)).genus()
            0
            sage: K23 = graphs.CompleteBipartiteGraph(2,3)
            sage: K23.genus()
            0
            sage: K33 = graphs.CompleteBipartiteGraph(3,3)
            sage: K33.genus()
            1

        Using the circular argument, we can compute the minimal genus preserving a planar, ordered boundary:
            sage: cube = graphs.CubeGraph(3)
            sage: cube.set_boundary(['001','110'])
            sage: cube.genus()
            0
            sage: cube.is_circular_planar()
            False
            sage: cube.genus(circular=True) #long time
            1
            sage: cube.genus(circular=True, maximal=True) #long time
            3
            sage: cube.genus(circular=True, on_embedding=True) #long time
            3
        """
        if not self.is_connected():
            raise TypeError("Graph must be connected to use Euler's Formula to compute minimal genus.")
        from sage.combinat.all import CyclicPermutationsOfPartition

        G = self.to_undirected()
        verts = G.order()
        edges = G.size()

        if circular:
            boundary = G.get_boundary()
            if hasattr(G, '_embedding'):
                del(G._embedding)

            extra = 0
            while G.has_vertex(extra):
                extra=extra+1
            G.add_vertex(extra)
            verts += 1

            for vertex in boundary:
                G.add_edge(vertex,extra)

            extra_edges = []
            if ordered: # WHEEL
                for i in range(len(boundary)-1):
                    if not G.has_edge(boundary[i],boundary[i+1]):
                        G.add_edge(boundary[i],boundary[i+1])
                        extra_edges.append((boundary[i],boundary[i+1]))
                if not G.has_edge(boundary[-1],boundary[0]):
                    G.add_edge(boundary[-1],boundary[0])
                    extra_edges.append((boundary[-1],boundary[0]))
                # else STAR (empty list of extra edges)

            edges = G.size()

        if on_embedding is not None:
            if on_embedding: #i.e., if on_embedding True (returns False if on_embedding is of type dict)
                if not hasattr(self,'_embedding'):
                    raise Exception("Graph must have attribute _embedding set to compute current (embedded) genus.")
                faces = len(self.trace_faces(self._embedding))
                return (2-verts+edges-faces)/2
            else: # compute genus on the provided dict
                faces = len(self.trace_faces(on_embedding))
                return (2-verts+edges-faces)/2
        else: # then compute either maximal or minimal genus of all embeddings
            # Construct an intitial combinatorial embedding for graph
            part = []
            for vertex in G.vertices():
                part.append(G.neighbors(vertex))

            # Iterate through all embeddings
            from sage.rings.infinity import infinity
            max_faces = -1
            min_faces = infinity
            min_embedding = []
            max_embedding = []
            labels = G.vertices()
            for p in CyclicPermutationsOfPartition(part):
                # Make dict of node labels embedding
                comb_emb = {}
                for i in range(len(p)):
                    comb_emb[labels[i]] = p[i]
                t = G.trace_faces(comb_emb)
                faces = len(t)
                if faces > max_faces:
                    max_faces = faces
                    min_embedding = comb_emb
                if faces < min_faces:
                    min_faces = faces
                    max_embedding = comb_emb

            if maximal:
                faces = min_faces
                emb = max_embedding
            else:
                faces = max_faces
                emb = min_embedding

            if set_embedding:
                if not circular:
                    # Make dict of node labels embedding
                    self._embedding = emb
                else: # for circular, we must strip extra vertices and edges from embedding
                    emb.pop(extra)
                    for u,v in extra_edges:
                        emb[u].pop(emb[u].index(v))
                        emb[v].pop(emb[v].index(u))
                    for w in boundary:
                        emb[w].pop(emb[w].index(extra))
                    self._embedding = emb

        return (2-verts+edges-faces)/2

    def trace_faces(self, comb_emb):
        """
        A helper function for finding the genus of a graph.
        Given a graph and a combinatorial embedding (rot_sys),
        this function will compute the faces (returned as a list
        of lists of edges (tuples) of the particular embedding.

        Note -- rot_sys is an ordered list based on the hash order
        of the vertices of graph.  To avoid confusion, it might be
        best to set the rot_sys based on a 'nice_copy' of the graph.

        INPUT:
            comb_emb -- a combinatorial embedding dictionary.  Format:
                    { v1:[v2,v3], v2:[v1], v3:[v1] } (clockwise ordering
                    of neighbors at each vertex.)

        EXAMPLES:
            sage: T = graphs.TetrahedralGraph()
            sage: T.trace_faces({0: [1, 3, 2], 1: [0, 2, 3], 2: [0, 3, 1], 3: [0, 1, 2]})
            [[(0, 1), (1, 2), (2, 0)],
             [(3, 2), (2, 1), (1, 3)],
             [(2, 3), (3, 0), (0, 2)],
             [(0, 3), (3, 1), (1, 0)]]

        """
        from sage.sets.set import Set

        # Establish set of possible edges
        edgeset = Set([])
        for edge in self.to_undirected().edges():
            edgeset = edgeset.union(Set([(edge[0],edge[1]),(edge[1],edge[0])]))

        # Storage for face paths
        faces = []
        path = []
        for edge in edgeset:
            path.append(edge)
            edgeset -= Set([edge])
            break  # (Only one iteration)

        # Trace faces
        while (len(edgeset) > 0):
            neighbors = comb_emb[path[-1][-1]]
            next_node = neighbors[(neighbors.index(path[-1][-2])+1)%(len(neighbors))]
            tup = (path[-1][-1],next_node)
            if tup == path[0]:
                faces.append(path)
                path = []
                for edge in edgeset:
                    path.append(edge)
                    edgeset -= Set([edge])
                    break  # (Only one iteration)
            else:
                path.append(tup)
                edgeset -= Set([tup])
        if (len(path) != 0): faces.append(path)
        return faces

    ### Connectivity

    def is_connected(self):
        """
        Indicates whether the (di)graph is connected. Note that in a graph, path
        connected is equivalent to connected.

        EXAMPLE:
            sage: G = Graph( { 0 : [1, 2], 1 : [2], 3 : [4, 5], 4 : [5] } )
            sage: G.is_connected()
            False
            sage: G.add_edge(0,3)
            sage: G.is_connected()
            True
            sage: D = DiGraph( { 0 : [1, 2], 1 : [2], 3 : [4, 5], 4 : [5] } )
            sage: D.is_connected()
            False
            sage: D.add_edge(0,3)
            sage: D.is_connected()
            True
            sage: D = DiGraph({1:[0], 2:[0]})
            sage: D.is_connected()
            True

        """
        if self.order() == 0:
            return True
        v = self.vertex_iterator().next()
        conn_verts = list(self.breadth_first_search(v, ignore_direction=True))
        return len(conn_verts) == self.num_verts()

    def connected_components(self):
        """
        Returns a list of lists of vertices, each list representing a
        connected component. The list is ordered from largest to smallest
        component.

        EXAMPLE:
            sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: G.connected_components()
            [[0, 1, 2, 3], [4, 5, 6]]
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: D.connected_components()
            [[0, 1, 2, 3], [4, 5, 6]]

        """
        seen = []
        components = []
        for v in self:
            if v not in seen:
                c = list(self.breadth_first_search(v, ignore_direction=True))
                c.sort()
                seen.extend(c)
                components.append(c)
        components.sort(lambda comp1, comp2: cmp(len(comp2), len(comp1)))
        return components

    def connected_components_number(self):
        """
        Returns the number of connected components.

        EXAMPLE:
            sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: G.connected_components_number()
            2
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: D.connected_components_number()
            2

        """
        return len(self.connected_components())

    def connected_components_subgraphs(self):
        """
        Returns a list of connected components as graph objects.

        EXAMPLE:
            sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: L = G.connected_components_subgraphs()
            sage: graphs_list.show_graphs(L)
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: L = D.connected_components_subgraphs()
            sage: graphs_list.show_graphs(L)

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
            sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
            sage: D.connected_component_containing_vertex(0)
            [0, 1, 2, 3]

        """
        c = list(self.breadth_first_search(vertex, ignore_direction=True))
        c.sort()
        return c

    def blocks_and_cut_vertices(self):
        """
        Computes the blocks and cut vertices of the graph. In the case of a
        digraph, this computation is done on the underlying graph.

        A cut vertex is one whose deletion increases the number of connected
        components. A block is a maximal induced subgraph which itself has no
        cut vertices. Two distinct blocks cannot overlap in more than a single
        cut vertex.

        OUTPUT:
        ( B, C ), where B is a list of blocks- each is a list of vertices and
            the blocks are the corresponding induced subgraphs- and C is a list
            of cut vertices.

        EXAMPLES:
            sage: graphs.PetersenGraph().blocks_and_cut_vertices()
            ([[0, 1, 2, 3, 8, 5, 7, 9, 4, 6]], [])
            sage: graphs.PathGraph(6).blocks_and_cut_vertices()
            ([[5, 4], [4, 3], [3, 2], [2, 1], [0, 1]], [4, 3, 2, 1])
            sage: graphs.CycleGraph(7).blocks_and_cut_vertices()
            ([[0, 1, 2, 3, 4, 5, 6]], [])
            sage: graphs.KrackhardtKiteGraph().blocks_and_cut_vertices()
            ([[9, 8], [8, 7], [0, 1, 3, 2, 5, 6, 4, 7]], [8, 7])

        ALGORITHM:
        8.3.8 in [1]. Notice the typo - the stack must also be considered as one
            of the blocks at termination.

        REFERENCE:
            [1] D. Jungnickel, Graphs, Networks and Algorithms, Springer, 2005.
        """
        G = self.to_undirected()
        map_to_ints = G.relabel(return_map=True)
        s = G.vertex_iterator().next()
        nr = [0]*G.num_verts()
        p = [0]*G.num_verts()
        L = [0]*G.num_verts()
        edges = {}
        for u,v in G.edges(labels=False):
            edges[(u,v)] = False
            edges[(v,u)] = False
        i = 1
        v = s
        nr[s] = 1
        L[s] = 1
        C = []
        B = []
        S = [s]
        while True:
            while True:
                u = -1
                for u in G.neighbor_iterator(v):
                    if not edges[(v,u)]: break
                if u == -1: break
                if edges[(v,u)]: break
                edges[(v,u)] = True
                edges[(u,v)] = True
                if nr[u] == 0:
                    p[u] = v
                    i += 1
                    nr[u] = i
                    L[u] = i
                    S.append(u)
                    v = u
                else:
                    L[v] = min([ L[v], nr[u] ])
            if p[v] != s:
                if L[v] < nr[p[v]]:
                    L[p[v]] = min([ L[p[v]], L[v] ])
                else:
                    C.append(p[v])
                    B_k = []
                    while True:
                        u = S.pop(-1)
                        B_k.append(u)
                        if u == v: break
                    B_k.append(p[v])
                    B.append(B_k)
            else:
                if any([ not edges[(s,u)] for u in G.neighbors(s)]):
                    C.append(s)
                B_k = []
                while True:
                    u = S.pop(-1)
                    B_k.append(u)
                    if u == v: break
                B_k.append(s)
                B.append(B_k)
            v = p[v]
            if p[v] == 0 and all([edges[(v,u)] for u in G.neighbors(v)]):
                break
        if len(S) != 0:
            B.append(S)
        return B, C

    ### Vertex handlers

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
            sage: G = Graph(sparse=True); G.add_vertex(); G
            Graph on 1 vertex

            sage: D = DiGraph(sparse=True); D.add_vertex(); D
            Digraph on 1 vertex

        """
        self._backend.add_vertex(name)

    def add_vertices(self, vertices):
        """
        Add vertices to the (di)graph from an iterable container of vertices.
        Vertices that already exist in the graph will not be added again.

        EXAMPLES:
            sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], 5: [7,8], 6: [8,9], 7: [9]}
            sage: G = Graph(d, sparse=True)
            sage: G.add_vertices([10,11,12])
            sage: G.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
            sage: G.add_vertices(graphs.CycleGraph(25).vertices())
            sage: G.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]

        """
        self._backend.add_vertices(vertices)

    def delete_vertex(self, vertex, in_order=False):
        """
        Deletes vertex, removing all incident edges.  Deleting a non-existant
        vertex will raise an exception.

        INPUT:
        in_order -- (default False) If True, this deletes the ith vertex in
            the sorted list of vertices, i.e. G.vertices()[i]

        EXAMPLES:
            sage: G = Graph(graphs.WheelGraph(9), sparse=True)
            sage: G.delete_vertex(0); G.show()

            sage: D = DiGraph({0:[1,2,3,4,5],1:[2],2:[3],3:[4],4:[5],5:[1]}, implementation='networkx')
            sage: D.delete_vertex(0); D
            Digraph on 5 vertices
            sage: D.vertices()
            [1, 2, 3, 4, 5]
            sage: D.delete_vertex(0)
            Traceback (most recent call last):
            ...
            NetworkXError: node 0 not in graph

            sage: G = graphs.CompleteGraph(4).line_graph(labels=False)
            sage: G.vertices()
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: G.delete_vertex(0, in_order=True)
            sage: G.vertices()
            [(0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

        """
        if in_order:
            vertex = self.vertices()[vertex]
        self._backend.del_vertex(vertex)

    def delete_vertices(self, vertices):
        """
        Remove vertices from the (di)graph taken from an iterable container of
        vertices.  Deleting a non-existant vertex will raise an exception.

        EXAMPLE:
            sage: D = DiGraph({0:[1,2,3,4,5],1:[2],2:[3],3:[4],4:[5],5:[1]}, implementation='networkx')
            sage: D.delete_vertices([1,2,3,4,5]); D
            Digraph on 1 vertex
            sage: D.vertices()
            [0]
            sage: D.delete_vertices([1])
            Traceback (most recent call last):
            ...
            NetworkXError: node 1 not in graph

        """
        self._backend.del_vertices(vertices)

    def has_vertex(self, vertex):
        """
        Return True if vertex is one of the vertices of this graph.

        INPUT:
            vertex -- an integer

        OUTPUT:
            bool -- True or False

        EXAMPLES:
            sage: g = Graph({0:[1,2,3], 2:[4]}); g
            Graph on 5 vertices
            sage: 2 in g
            True
            sage: 10 in g
            False
            sage: graphs.PetersenGraph().has_vertex(99)
            False

        """
        return self._backend.has_vertex(vertex)

    __contains__ = has_vertex

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

        In a digraph, the external boundary of a vertex v are those vertices u
        with an arc (v, u).

        EXAMPLE:
            sage: G = graphs.CubeGraph(4)
            sage: l = ['0111', '0000', '0001', '0011', '0010', '0101', '0100', '1111', '1101', '1011', '1001']
            sage: G.vertex_boundary(['0000', '1111'], l)
            ['0111', '0001', '0010', '0100', '1101', '1011']

            sage: D = DiGraph({0:[1,2], 3:[0]})
            sage: D.vertex_boundary([0])
            [1, 2]

        """
        vertices1 = [v for v in vertices1 if v in self]
        output = set()
        if self._directed:
            for v in vertices1:
                output.update(self.successor_iterator(v))
        else:
            for v in vertices1:
                output.update(self.neighbor_iterator(v))
        if vertices2 is not None:
            output.intersection_update(vertices2)
        return list(output)

    def set_vertices(self, vertex_dict):
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
            sage: T.set_vertices(d)
            sage: T.get_vertex(1)
            Flower Snark: Graph on 20 vertices

        """
        try:
            for v in vertex_dict:
                self._assoc[v] = vertex_dict[v]
        except:
            self._assoc = {}
            for v in vertex_dict:
                self._assoc[v] = vertex_dict[v]

    def set_vertex(self, vertex, object):
        """
        Associate an arbitrary object with a vertex.

        INPUT:
            vertex -- which vertex
            object -- object to associate to vertex

        EXAMPLE:
            sage: T = graphs.TetrahedralGraph()
            sage: T.vertices()
            [0, 1, 2, 3]
            sage: T.set_vertex(1, graphs.FlowerSnark())
            sage: T.get_vertex(1)
            Flower Snark: Graph on 20 vertices

        """
        try:
            self._assoc[vertex] = object
        except:
            self._assoc = {}
            self._assoc[vertex] = object

    def get_vertex(self, vertex):
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
            sage: T.set_vertices(d)
            sage: T.get_vertex(1)
            Flower Snark: Graph on 20 vertices

        """
        try:
            return self._assoc[vertex]
        except:
            return None

    def get_vertices(self, verts=None):
        """
        Return a dictionary of the objects associated to each vertex.

        INPUT:
            verts -- iterable container of vertices

        EXAMPLES:
            sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), 2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
            sage: T = graphs.TetrahedralGraph()
            sage: T.set_vertices(d)
            sage: T.get_vertices([1,2])
            {1: Flower Snark: Graph on 20 vertices,
             2: Moebius-Kantor Graph: Graph on 16 vertices}

        """
        if verts is None:
            verts = self.vertices()
        output = {}
        for v in verts:
            try:
                output[v] = self._assoc[v]
            except:
                output[v] = None
        return output

    def loop_vertices(self):
        """
        Returns a list of vertices with loops.

        EXAMPLE:
            sage: G = Graph({0 : [0], 1: [1,2,3], 2: [3]}, loops=True)
            sage: G.loop_vertices()
            [0, 1]

        """
        if self.loops():
            return [v for v in self if self.has_edge(v,v)]
        else:
            return []

    def vertex_iterator(self, vertices=None):
        """
        Returns an iterator over the given vertices. Returns False if not given
        a vertex, sequence, iterator or None. None is equivalent to a list of
        every vertex. Note that \code{for v in G} syntax is allowed.

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
            ...
            8
            9

            sage: G = graphs.TetrahedralGraph()
            sage: for i in G:
            ...    print i
            0
            1
            2
            3

        Note that since the intersection option is available, the
        vertex_iterator() function is sub-optimal, speedwise, but note the
        following optimization:
            sage: timeit V = P.vertices()                   # not tested
            100000 loops, best of 3: 8.85 [micro]s per loop
            sage: timeit V = list(P.vertex_iterator())      # not tested
            100000 loops, best of 3: 5.74 [micro]s per loop
            sage: timeit V = list(P._nxg.adj.iterkeys())    # not tested
            100000 loops, best of 3: 3.45 [micro]s per loop

        In other words, if you want a fast vertex iterator, call the dictionary
        directly.

        """
        return self._backend.iterator_verts(vertices)

    __iter__ = vertex_iterator

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
            sage: D = G.to_directed()
            sage: for i in D.neighbor_iterator('010'):
            ...    print i
            011
            000
            110

            sage: D = DiGraph({0:[1,2], 3:[0]})
            sage: list(D.neighbor_iterator(0))
            [1, 2, 3]

        """
        if self._directed:
            return iter(set(self.successor_iterator(vertex)) \
                    | set(self.predecessor_iterator(vertex)))
        else:
            return self._backend.iterator_nbrs(vertex)

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
            sage: timeit V = P.vertices()                     # not tested
            100000 loops, best of 3: 8.85 [micro]s per loop
            sage: timeit V = list(P.vertex_iterator())        # not tested
            100000 loops, best of 3: 5.74 [micro]s per loop
            sage: timeit V = list(P._nxg.adj.iterkeys())      # not tested
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

    def neighbors(self, vertex):
        """
        Return a list of neighbors (in and out if directed) of vertex.

        G[vertex] also works.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: sorted(P.neighbors(3))
            [2, 4, 8]
            sage: sorted(P[4])
            [0, 3, 9]

        """
        return list(self.neighbor_iterator(vertex))

    __getitem__ = neighbors

    ### Edge handlers

    def add_edge(self, u, v=None, label=None):
        """
        Adds an edge from u and v.

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
            sage: G = Graph(implementation='networkx')
            sage: G.add_edge((1,2), 'label')
            sage: G.networkx_graph().adj           # random output order
            {'label': {(1, 2): None}, (1, 2): {'label': None}}

        Use one of these instead:
            sage: G = Graph(implementation='networkx')
            sage: G.add_edge((1,2), label="label")
            sage: G.networkx_graph().adj           # random output order
            {1: {2: 'label'}, 2: {1: 'label'}}

            sage: G = Graph(implementation='networkx')
            sage: G.add_edge(1,2,'label')
            sage: G.networkx_graph().adj           # random output order
            {1: {2: 'label'}, 2: {1: 'label'}}

        """
        if label is None:
            if v is None:
                try:
                    u, v, label = u
                except:
                    u, v = u
                    label = None
        if not self.loops() and u==v:
            return
        self._backend.add_edge(u, v, label, self._directed)

    def add_edges(self, edges):
        """
        Add edges from an iterable container.

        EXAMPLE:
            sage: G = graphs.DodecahedralGraph()
            sage: H = Graph(implementation='networkx')
            sage: H.add_edges( G.edge_iterator() ); H
            Graph on 20 vertices
            sage: G = graphs.DodecahedralGraph().to_directed()
            sage: H = DiGraph(implementation='networkx')
            sage: H.add_edges( G.edge_iterator() ); H
            Digraph on 20 vertices

        """
        self._backend.add_edges(edges, self._directed)

    def delete_edge(self, u, v=None, label=None):
        r"""
        Delete the edge from u to v, returning silently if vertices or edge does
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

        """
        self._backend.del_edge(u, v, label, self._directed)

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

            sage: K12 = graphs.CompleteGraph(12).to_directed()
            sage: K4 = graphs.CompleteGraph(4).to_directed()
            sage: K12.size()
            132
            sage: K12.delete_edges(K4.edge_iterator())
            sage: K12.size()
            120

        """
        for e in edges:
            self.delete_edge(e)

    def delete_multiedge(self, u, v):
        """
        Deletes all edges from u and v.

        EXAMPLE:
            sage: G = Graph(multiedges=True, implementation='networkx')
            sage: G.add_edges([(0,1), (0,1), (0,1), (1,2), (2,3)])
            sage: G.edges()
            [(0, 1, None), (0, 1, None), (0, 1, None), (1, 2, None), (2, 3, None)]
            sage: G.delete_multiedge( 0, 1 )
            sage: G.edges()
            [(1, 2, None), (2, 3, None)]

            sage: D = DiGraph(multiedges=True)
            sage: D.add_edges([(0,1,1), (0,1,2), (0,1,3), (1,0), (1,2), (2,3)])
            sage: D.edges()
            [(0, 1, 1), (0, 1, 2), (0, 1, 3), (1, 0, None), (1, 2, None), (2, 3, None)]
            sage: D.delete_multiedge( 0, 1 )
            sage: D.edges()
            [(1, 0, None), (1, 2, None), (2, 3, None)]

        """
        if self.multiple_edges():
            for l in self.edge_label(u, v):
                self.delete_edge(u, v, l)
        else:
            self.delete_edge(u, v)

    def set_edge_label(self, u, v, l):
        """
        Set the edge label of a given edge.

        NOTE:
            There can be only one edge from u to v for this to make sense.
        Otherwise, an error is raised.

        INPUT:
            u, v -- the vertices (and direction if digraph) of the edge
            l -- the new label

        EXAMPLES:
            sage: SD = DiGraph( { 1:[18,2], 2:[5,3], 3:[4,6], 4:[7,2], 5:[4], 6:[13,12], 7:[18,8,10], 8:[6,9,10], 9:[6], 10:[11,13], 11:[12], 12:[13], 13:[17,14], 14:[16,15], 15:[2], 16:[13], 17:[15,13], 18:[13] }, implementation='networkx')
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
            sage: SD.plot(pos=posn, vertex_size=400, vertex_colors={'#FFFFFF':range(1,19)}, edge_labels=True).show()

            sage: G = graphs.HeawoodGraph()
            sage: for u,v,l in G.edges():
            ...    G.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: G.edges()
                [(0, 1, '(0,1)'),
                 (0, 5, '(0,5)'),
                 (0, 13, '(0,13)'),
                 ...
                 (11, 12, '(11,12)'),
                 (12, 13, '(12,13)')]

            sage: g = Graph({0: [0,1,1,2]}, loops=True, multiedges=True, implementation='networkx')
            sage: g.set_edge_label(0,0,'test')
            sage: g.edges()
            [(0, 0, 'test'), (0, 1, None), (0, 1, None), (0, 2, None)]
            sage: g.add_edge(0,0,'test2')
            sage: g.set_edge_label(0,0,'test3')
            Traceback (most recent call last):
            ...
            RuntimeError: Cannot set edge label, since there are multiple edges from 0 to 0.

            sage: dg = DiGraph({0 : [1], 1 : [0]}, sparse=True)
            sage: dg.set_edge_label(0,1,5)
            sage: dg.set_edge_label(1,0,9)
            sage: dg.outgoing_edges(1)
            [(1, 0, 9)]
            sage: dg.incoming_edges(1)
            [(0, 1, 5)]
            sage: dg.outgoing_edges(0)
            [(0, 1, 5)]
            sage: dg.incoming_edges(0)
            [(1, 0, 9)]

            sage: G = Graph({0:{1:1}}, implementation='c_graph')
            sage: G.num_edges()
            1
            sage: G.set_edge_label(0,1,1)
            sage: G.num_edges()
            1

        """
        if self.multiple_edges():
            if len(self.edge_label(u, v)) > 1:
                raise RuntimeError("Cannot set edge label, since there are multiple edges from %s to %s."%(u,v))
        self._backend.set_edge_label(u, v, l, self._directed)

    def has_edge(self, u, v=None, label=None):
        r"""
        Returns True if (u, v) is an edge, False otherwise.

        INPUT:
        The following forms are accepted by NetworkX:

        G.has_edge( 1, 2 )
        G.has_edge( (1, 2) )
        G.has_edge( 1, 2, 'label' )

        EXAMPLE:
            sage: graphs.EmptyGraph().has_edge(9,2)
            False
            sage: DiGraph().has_edge(9,2)
            False
            sage: G = Graph(implementation='networkx')
            sage: G.add_edge(0,1,"label")
            sage: G.has_edge(0,1,"different label")
            False
            sage: G.has_edge(0,1,"label")
            True

        """
        return self._backend.has_edge(u, v, label)

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

            sage: D = graphs.DodecahedralGraph().to_directed()
            sage: D.edges()
            [(0, 1, None), (0, 10, None), (0, 19, None), (1, 0, None), (1, 2, None), (1, 8, None), (2, 1, None), (2, 3, None), (2, 6, None), (3, 2, None), (3, 4, None), (3, 19, None), (4, 3, None), (4, 5, None), (4, 17, None), (5, 4, None), (5, 6, None), (5, 15, None), (6, 2, None), (6, 5, None), (6, 7, None), (7, 6, None), (7, 8, None), (7, 14, None), (8, 1, None), (8, 7, None), (8, 9, None), (9, 8, None), (9, 10, None), (9, 13, None), (10, 0, None), (10, 9, None), (10, 11, None), (11, 10, None), (11, 12, None), (11, 18, None), (12, 11, None), (12, 13, None), (12, 16, None), (13, 9, None), (13, 12, None), (13, 14, None), (14, 7, None), (14, 13, None), (14, 15, None), (15, 5, None), (15, 14, None), (15, 16, None), (16, 12, None), (16, 15, None), (16, 17, None), (17, 4, None), (17, 16, None), (17, 18, None), (18, 11, None), (18, 17, None), (18, 19, None), (19, 0, None), (19, 3, None), (19, 18, None)]
            sage: D.edges(labels = False)
            [(0, 1), (0, 10), (0, 19), (1, 0), (1, 2), (1, 8), (2, 1), (2, 3), (2, 6), (3, 2), (3, 4), (3, 19), (4, 3), (4, 5), (4, 17), (5, 4), (5, 6), (5, 15), (6, 2), (6, 5), (6, 7), (7, 6), (7, 8), (7, 14), (8, 1), (8, 7), (8, 9), (9, 8), (9, 10), (9, 13), (10, 0), (10, 9), (10, 11), (11, 10), (11, 12), (11, 18), (12, 11), (12, 13), (12, 16), (13, 9), (13, 12), (13, 14), (14, 7), (14, 13), (14, 15), (15, 5), (15, 14), (15, 16), (16, 12), (16, 15), (16, 17), (17, 4), (17, 16), (17, 18), (18, 11), (18, 17), (18, 19), (19, 0), (19, 3), (19, 18)]

        """
        L = list(self.edge_iterator(labels=labels))
        if sort:
            L.sort()
        return L

    def edge_boundary(self, vertices1, vertices2=None, labels=True):
        """
        Returns a list of edges (u,v,l) with u in vertices1 and v in vertices2.
        If vertices2 is None, then it is set to the complement of vertices1.

        In a digraph, the external boundary of a vertex v are those vertices u
        with an arc (v, u).

        INPUT:
        labels -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLE:
            sage: K = graphs.CompleteBipartiteGraph(9,3)
            sage: len(K.edge_boundary( [0,1,2,3,4,5,6,7,8], [9,10,11] ))
            27
            sage: K.size()
            27

        Note that the edge boundary preserves direction:
            sage: K = graphs.CompleteBipartiteGraph(9,3).to_directed()
            sage: len(K.edge_boundary( [0,1,2,3,4,5,6,7,8], [9,10,11] ))
            27
            sage: K.size()
            54

            sage: D = DiGraph({0:[1,2], 3:[0]})
            sage: D.edge_boundary([0])
            [(0, 1, None), (0, 2, None)]
            sage: D.edge_boundary([0], labels=False)
            [(0, 1), (0, 2)]

        """
        vertices1 = [v for v in vertices1 if v in self]
        output = []
        if self._directed:
            output.extend(self.outgoing_edge_iterator(vertices1))
        else:
            output.extend(self.edge_iterator(vertices1))
        if vertices2 is not None:
            output = [e for e in output if e[1] in vertices2]
        else:
            output = list(output)
        if not labels:
            output = [(u,v) for u,v,l in output]
        return output

    def edge_iterator(self, vertices=None, labels=True, ignore_direction=False):
        """
        Returns an iterator over the edges incident with any vertex given. If
        the graph is directed, iterates over edges going out only. If vertices
        is None, then returns an iterator over all edges. If self is directed,
        returns outgoing edges only.

        INPUT:
        labels -- if False, each edge is a tuple (u,v) of vertices.
        ignore_direction -- (default False) only applies to directed graphs. If
            True, searches across edges in either direction.

        EXAMPLE:
            sage: for i in graphs.PetersenGraph().edge_iterator([0]):
            ...    print i
            (0, 1, None)
            (0, 4, None)
            (0, 5, None)
            sage: D = DiGraph( { 0 : [1,2], 1: [0] } )
            sage: for i in D.edge_iterator([0]):
            ...    print i
            (0, 1, None)
            (0, 2, None)

            sage: G = graphs.TetrahedralGraph()
            sage: list(G.edge_iterator(labels=False))
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

            sage: D = DiGraph({1:[0], 2:[0]})
            sage: list(D.edge_iterator(0))
            []
            sage: list(D.edge_iterator(0, ignore_direction=True))
            [(1, 0, None), (2, 0, None)]

        """
        if vertices is None:
            vertices = self
        elif vertices in self:
            vertices = [vertices]
        else:
            vertices = [v for v in vertices if v in self]
        if ignore_direction and self._directed:
            for v in self._backend.iterator_out_edges(vertices, labels):
                yield v
            for v in self._backend.iterator_in_edges(vertices, labels):
                yield v
        elif self._directed:
            for v in self._backend.iterator_out_edges(vertices, labels):
                yield v
        else:
            for v in self._backend.iterator_edges(vertices, labels, self._directed):
                yield v

    def edges_incident(self, vertices=None, labels=True):
        """
        Returns a list of edges incident with any vertex given. If vertices is
        None, returns a list of all edges in graph. For digraphs, only lists
        outward edges.

        INPUT:
        label -- if False, each edge is a tuple (u,v) of vertices.

        EXAMPLE:
            sage: graphs.PetersenGraph().edges_incident([0,9], labels=False)
            [(0, 1), (0, 4), (0, 5), (9, 4), (9, 6), (9, 7)]
            sage: D = DiGraph({0:[1]})
            sage: D.edges_incident([0])
            [(0, 1, None)]
            sage: D.edges_incident([1])
            []

        """
        if vertices in self:
            vertices = [vertices]
        v = list(self.edge_boundary(vertices, labels=labels))
        v.sort()
        return v

    def edge_label(self, u, v=None):
        """
        Returns the label of an edge. Note that if the graph allows multiple
        edges, then a list of labels on the edge is returned.

        EXAMPLE:
            sage: G = Graph({0 : {1 : 'edgelabel'}}, implementation='networkx')
            sage: G.edges(labels=False)
            [(0, 1)]
            sage: G.edge_label( 0, 1 )
            'edgelabel'
            sage: D = DiGraph({0 : {1 : 'edgelabel'}}, implementation='networkx')
            sage: D.edges(labels=False)
            [(0, 1)]
            sage: D.edge_label( 0, 1 )
            'edgelabel'

            sage: G = Graph(multiedges=True)
            sage: [G.add_edge(0,1,i) for i in range(1,6)]
            [None, None, None, None, None]
            sage: sorted(G.edge_label(0,1))
            [1, 2, 3, 4, 5]

        """
        return self._backend.get_edge_label(u,v)

    def edge_labels(self):
        """
        Returns a list of edge labels.

        EXAMPLE:
            sage: G = Graph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}, implementation='networkx')
            sage: G.edge_labels()
            ['x', 'z', 'a', 'out']
            sage: G = DiGraph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}, implementation='networkx')
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
            sage: G = Graph(multiedges=True, implementation='networkx')
            sage: G.add_edges( [ (0,1), (0,1), (0,1), (0,1), (1,2) ] )
            sage: G.edges(labels=False)
            [(0, 1), (0, 1), (0, 1), (0, 1), (1, 2)]

            sage: G.remove_multiple_edges()
            sage: G.edges(labels=False)
            [(0, 1), (1, 2)]

            sage: D = DiGraph(multiedges=True)
            sage: D.add_edges( [ (0,1,1), (0,1,2), (0,1,3), (0,1,4), (1,2) ] )
            sage: D.edges(labels=False)
            [(0, 1), (0, 1), (0, 1), (0, 1), (1, 2)]
            sage: D.remove_multiple_edges()
            sage: D.edges(labels=False)
            [(0, 1), (1, 2)]

        """
        if self.multiple_edges():
            if self._directed:
                for v in self:
                    for u in self.predecessor_iterator(v):
                        edges = self.edge_boundary([u], [v])
                        if len(edges) > 1:
                            self.delete_edges(edges[1:])
            else:
                for v in self:
                    for u in self.neighbor_iterator(v):
                        edges = self.edge_boundary([v], [u])
                        if len(edges) > 1:
                            self.delete_edges(edges[1:])

    def remove_loops(self, vertices=None):
        """
        Removes loops on vertices in vertices. If vertices is None, removes all loops.

        EXAMPLE
            sage: G = Graph(4, loops=True)
            sage: G.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: G.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: G.remove_loops()
            sage: G.edges(labels=False)
            [(2, 3)]
            sage: G.loops()
            True

            sage: D = DiGraph(4, loops=True)
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
            vertices = self
        for v in vertices:
            if self.has_edge(v,v):
                self.delete_multiedge(v,v)

    def loop_edges(self):
        """
        Returns a list of all loops in the graph.

        EXAMPLE:
            sage: G = Graph(4, loops=True)
            sage: G.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: G.loop_edges()
            [(0, 0, None), (1, 1, None), (2, 2, None), (3, 3, None)]

            sage: D = DiGraph(4, loops=True)
            sage: D.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: D.loop_edges()
            [(0, 0, None), (1, 1, None), (2, 2, None), (3, 3, None)]

            sage: G = Graph(4, loops=True, multiedges=True)
            sage: G.add_edges([(i,i) for i in range(4)])
            sage: G.loop_edges()
            [(0, 0, None), (1, 1, None), (2, 2, None), (3, 3, None)]

        """
        if self.multiple_edges():
            return [(v,v,l) for v in self.loop_vertices() for l in self.edge_label(v,v)]
        else:
            return [(v,v,self.edge_label(v,v)) for v in self.loop_vertices()]

    def number_of_loops(self):
        """
        Returns the number of edges that are loops.

        EXAMPLE:
            sage: G = Graph(4, loops=True)
            sage: G.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: G.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: G.number_of_loops()
            4

            sage: D = DiGraph(4, loops=True)
            sage: D.add_edges( [ (0,0), (1,1), (2,2), (3,3), (2,3) ] )
            sage: D.edges(labels=False)
            [(0, 0), (1, 1), (2, 2), (2, 3), (3, 3)]
            sage: D.number_of_loops()
            4

        """
        return len(self.loop_edges())

    ### Modifications

    def clear(self):
        """
        Empties the graph of vertices and edges and removes name,
        boundary, associated objects, and position information.

        EXAMPLE:
            sage: G=graphs.CycleGraph(4); G.set_vertices({0:'vertex0'})
            sage: G.order(); G.size()
            4
            4
            sage: len(G._pos)
            4
            sage: G.name()
            'Cycle graph'
            sage: G.get_vertex(0)
            'vertex0'
            sage: H = G.copy(implementation='c_graph', sparse=True)
            sage: H.clear()
            sage: H.order(); H.size()
            0
            0
            sage: len(H._pos)
            0
            sage: H.name()
            ''
            sage: H.get_vertex(0)
            sage: H = G.copy(implementation='networkx')
            sage: H.clear()
            sage: H.order(); H.size()
            0
            0
            sage: len(H._pos)
            0
            sage: H.name()
            ''
            sage: H.get_vertex(0)

        """
        self.name('')
        self.delete_vertices(self.vertices())
        self._pos=[]
        self._boundary=[]
        self._assoc=None

    ### Degree functions

    def degree(self, vertices=None, labels=False):
        """
        Gives the degree (in + out for digraphs) of a vertex or of vertices.

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

            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.degree(vertices = [0,1,2], labels=True)
            {0: 5, 1: 4, 2: 3}
            sage: D.degree()
            [5, 4, 3, 3, 3, 2]

        """
        if labels:
            return dict(self.degree_iterator(vertices,labels))
        elif vertices in self and not labels:
            return self.degree_iterator(vertices,labels).next()
        else:
            return list(self.degree_iterator(vertices,labels))

    def degree_histogram(self):
        """
        Returns a list, whose ith entry is the frequency of degree i.

        EXAMPLE:
            sage: G = graphs.Grid2dGraph(9,12)
            sage: G.degree_histogram()
            [0, 0, 4, 34, 70]

            sage: G = graphs.Grid2dGraph(9,12).to_directed()
            sage: G.degree_histogram()
            [0, 0, 0, 0, 4, 0, 34, 0, 70]

        """
        degree_sequence = self.degree()
        dmax = max(degree_sequence) + 1
        frequency = [0]*dmax
        for d in degree_sequence:
            frequency[d] += 1
        return frequency

    def degree_iterator(self, vertices=None, labels=False):
        """
        Returns an iterator over the degrees of the (di)graph. In the case of a
        digraph, the degree is defined as the sum of the in-degree and the
        out-degree, i.e. the total number of edges incident to a given vertex.

        INPUT:
        labels=False: returns an iterator over degrees.
        labels=True: returns an iterator over tuples (vertex, degree).
        vertices -- if specified, restrict to this subset.

        EXAMPLES:
            sage: G = graphs.Grid2dGraph(3,4)
            sage: for i in G.degree_iterator():
            ...    print i
            3
            4
            2
            ...
            2
            3
            sage: for i in G.degree_iterator(labels=True):
            ...    print i
            ((0, 1), 3)
            ((1, 2), 4)
            ((0, 0), 2)
            ...
            ((0, 3), 2)
            ((0, 2), 3)

            sage: D = graphs.Grid2dGraph(2,4).to_directed()
            sage: for i in D.degree_iterator():
            ...    print i
            6
            6
            ...
            4
            6
            sage: for i in D.degree_iterator(labels=True):
            ...    print i
            ((0, 1), 6)
            ((1, 2), 6)
            ...
            ((0, 3), 4)
            ((1, 1), 6)

        """
        if vertices is None:
            vertices = self
        elif vertices in self:
            vertices = [vertices]
        else:
            vertices = [v for v in vertices if v in self]
        if labels:
            filter = lambda v, self: (v, self._backend.degree(v, self._directed))
        else:
            filter = lambda v, self: self._backend.degree(v, self._directed)
        for v in vertices:
            yield filter(v, self)

    ### Substructures

    def subgraph(self, vertices=None, edges=None, inplace=False,
                       vertex_property=None, edge_property=None):
        """
        Returns the subgraph containing the given vertices and edges.
        If either vertices or edges are not specified, they are
        assumed to be all vertices or edges.  If edges are not
        specified, returns the subgraph induced by the vertices.

        INPUT:
        inplace -- Using inplace is True will simply delete the extra vertices
            and edges from the current graph. This will modify the graph.
        vertices -- Vertices can be a single vertex or an iterable container
            of vertices, e.g. a list, set, graph, file or numeric array.
            If not passed, defaults to the entire graph.
        edges -- As with vertices, edges can be a single edge or an iterable
            container of edges (e.g., a list, set, file, numeric array,
            etc.).  If not edges are not specified, then all edges are
            assumed and the returned graph is an induced subgraph.  In
            the case of multiple edges, specifying an edge as (u,v)
            means to keep all edges (u,v), regardless of the label.
        vertex_property -- If specified, this is expected to be a function on
            vertices, which is intersected with the vertices specified, if any are.
        edge_property -- If specified, this is expected to be a function on
            edges, which is intersected with the edges specified, if any are.

        EXAMPLES:
            sage: G = graphs.CompleteGraph(9)
            sage: H = G.subgraph([0,1,2]); H
            Subgraph of (Complete graph): Graph on 3 vertices
            sage: G
            Complete graph: Graph on 9 vertices
            sage: J = G.subgraph(edges=[(0,1)])
            sage: J.edges(labels=False)
            [(0, 1)]
            sage: J.vertices()==G.vertices()
            True
            sage: G.subgraph([0,1,2], inplace=True); G
            Subgraph of (Complete graph): Graph on 3 vertices
            sage: G.subgraph()==G
            True

            sage: D = graphs.CompleteGraph(9).to_directed()
            sage: H = D.subgraph([0,1,2]); H
            Subgraph of (Complete graph): Digraph on 3 vertices
            sage: H = D.subgraph(edges=[(0,1), (0,2)])
            sage: H.edges(labels=False)
            [(0, 1), (0, 2)]
            sage: H.vertices()==D.vertices()
            True
            sage: D
            Complete graph: Digraph on 9 vertices
            sage: D.subgraph([0,1,2], inplace=True); D
            Subgraph of (Complete graph): Digraph on 3 vertices
            sage: D.subgraph()==D
            True

        A more complicated example involving multiple edges and labels.

            sage: G = Graph(multiedges=True, implementation='networkx')
            sage: G.add_edges([(0,1,'a'), (0,1,'b'), (1,0,'c'), (0,2,'d'), (0,2,'e'), (2,0,'f'), (1,2,'g')])
            sage: G.subgraph(edges=[(0,1), (0,2,'d'), (0,2,'not in graph')]).edges()
            [(0, 1, 'a'), (0, 1, 'b'), (0, 1, 'c'), (0, 2, 'd')]
            sage: J = G.subgraph(vertices=[0,1], edges=[(0,1,'a'), (0,2,'c')])
            sage: J.edges()
            [(0, 1, 'a')]
            sage: J.vertices()
            [0, 1]

            sage: D = DiGraph(multiedges=True, implementation='networkx')
            sage: D.add_edges([(0,1,'a'), (0,1,'b'), (1,0,'c'), (0,2,'d'), (0,2,'e'), (2,0,'f'), (1,2,'g')])
            sage: D.subgraph(edges=[(0,1), (0,2,'d'), (0,2,'not in graph')]).edges()
            [(0, 1, 'a'), (0, 1, 'b'), (0, 2, 'd')]
            sage: H = D.subgraph(vertices=[0,1], edges=[(0,1,'a'), (0,2,'c')])
            sage: H.edges()
            [(0, 1, 'a')]
            sage: H.vertices()
            [0, 1]

        Using the property arguments:
            sage: P = graphs.PetersenGraph()
            sage: S = P.subgraph(vertex_property = lambda v : v%2 == 0)
            sage: S.vertices()
            [0, 2, 4, 6, 8]

            sage: C = graphs.CubeGraph(2)
            sage: S = C.subgraph(edge_property=(lambda e: e[0][0] == e[1][0]))
            sage: C.edges()
            [('00', '01', None),
             ('10', '00', None),
             ('11', '01', None),
             ('11', '10', None)]
            sage: S.edges()
            [('00', '01', None), ('11', '10', None)]

        TESTS:
        We should delete unused _pos dictionary entries
            sage: g = graphs.PathGraph(10)
            sage: h = g.subgraph([3..5])
            sage: h._pos.keys()
            [3, 4, 5]
        """
        if vertices is None:
            vertices = []
        elif vertices in self:
            vertices = [v for v in self if v != vertices]
        else:
            vertices = [v for v in self if v not in vertices]
        if vertex_property is not None:
            for v in self:
                if v not in vertices and not vertex_property(v):
                    vertices.append(v)
        if inplace:
            G = self
        else:
            G = self.copy()
        if G._pos is not None:
            for v in vertices:
                del G._pos[v]
        G.name("Subgraph of (%s)"%self.name())
        if edges is None and edge_property is not None:
            G.delete_edges([e for e in self.edge_iterator() if not edge_property(e)])
        elif edges is not None:
            if G._directed:
                edges_graph = G.edges()
                edges_to_keep_labeled = [e for e in edges if len(e)==3]
                edges_to_keep_unlabeled = [e for e in edges if len(e)==2]
            else:
                edges_graph = [sorted(e[0:2])+[e[2]] for e in G.edges()]
                edges_to_keep_labeled = [sorted(e[0:2])+[e[2]] for e in edges if len(e)==3]
                edges_to_keep_unlabeled = [sorted(e) for e in edges if len(e)==2]
            edges_to_delete=[]
            for e in edges_graph:
                if e not in edges_to_keep_labeled and e[0:2] not in edges_to_keep_unlabeled:
                    edges_to_delete.append(tuple(e))
            if edge_property is not None:
                for e in G.edge_iterator():
                    if e not in edges_to_delete and not edge_property(e):
                        edges_to_delete.append(e)
            G.delete_edges(edges_to_delete)
        G.delete_vertices(vertices)
        if not inplace:
            return G

    def random_subgraph(self, p, inplace=False):
        """
        Return a random subgraph that contains each vertex with prob. p.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: P.random_subgraph(.25)
            Subgraph of (Petersen graph): Graph on 4 vertices

        """
        vertices = []
        p = float(p)
        for v in self:
            if random() < p:
                vertices.append(v)
        return self.subgraph(vertices=vertices, inplace=inplace)

    def is_clique(self, vertices=None, directed_clique=False):
        """
        Returns True if the set \code{vertices} is a clique, False if not.  A
        clique is a set of vertices such that there is an edge between
        any two vertices.

        INPUT:
            vertices -- Vertices can be a single vertex or an iterable
            container of vertices, e.g. a list, set, graph, file or
            numeric array.  If not passed, defaults to the entire graph.
            directed_clique -- (default False) If set to False, only
        consider the underlying undirected graph.  If set to True and the
        graph is directed, only return True if all possible edges in
        _both_ directions exist.

        EXAMPLE:
            sage: g = graphs.CompleteGraph(4)
            sage: g.is_clique([1,2,3])
            True
            sage: g.is_clique()
            True
            sage: h = graphs.CycleGraph(4)
            sage: h.is_clique([1,2])
            True
            sage: h.is_clique([1,2,3])
            False
            sage: h.is_clique()
            False
            sage: i = graphs.CompleteGraph(4).to_directed()
            sage: i.delete_edge([0,1])
            sage: i.is_clique()
            True
            sage: i.is_clique(directed_clique=True)
            False

        """
        if directed_clique and self._directed:
            subgraph=self.subgraph(vertices)
            subgraph.loops(False)
            subgraph.multiple_edges(False)
            n=subgraph.order()
            return subgraph.size()==n*(n-1)
        else:
            subgraph=self.subgraph(vertices).to_simple()
            n=subgraph.order()
            return subgraph.size()==n*(n-1)/2

    def is_independent_set(self, vertices=None):
        """
        Returns True if the set \code{vertices} is an independent set, False
        if not.  An independent set is a set of vertices such that
        there is no edge between any two vertices.

        INPUT:
            vertices -- Vertices can be a single vertex or an iterable
            container of vertices, e.g. a list, set, graph, file or
            numeric array.  If not passed, defaults to the entire graph.

        EXAMPLE:
            sage: graphs.CycleGraph(4).is_independent_set([1,3])
            True
            sage: graphs.CycleGraph(4).is_independent_set([1,2,3])
            False

        """
        return self.subgraph(vertices).to_simple().size()==0

    def is_subgraph(self, other):
        """
        Tests whether self is a subgraph of other.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: G = P.subgraph(range(6))
            sage: G.is_subgraph(P)
            True

        """
        self_verts = self.vertices()
        for v in self_verts:
            if v not in other:
                return False
        return other.subgraph(self_verts) == self

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
        return networkx.triangles(self.networkx_graph(copy=False), nbunch, with_labels)

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
        return networkx.average_clustering(self.networkx_graph(copy=False))

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
        return networkx.clustering(self.networkx_graph(copy=False), nbunch, with_labels, weights)

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
        return networkx.transitivity(self.networkx_graph(copy=False))

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
        return networkx.cores.find_cores(self.networkx_graph(copy=False), with_labels)

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
        if v is None:
            v = self.vertices()
        elif not isinstance(v, list):
            v = [v]
        e = {}
        infinite = False
        for u in v:
            if dist_dict is None:
                length = self.shortest_path_lengths(u)
            else:
                length = dist_dict[u]
            if len(length) != self.num_verts():
                infinite = True
                break
            e[u] = max(length.values())
        if infinite:
            from sage.rings.infinity import Infinity
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
            sage: G = Graph(sparse=True)
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

    def girth(self):
        """
        Computes the girth of the graph. For directed graphs, computes the girth
        of the undirected graph.

        The girth is the length of the shortest cycle in the graph. Graphs
        without cycles have infinite girth.

        EXAMPLES:
            sage: graphs.TetrahedralGraph().girth()
            3
            sage: graphs.CubeGraph(3).girth()
            4
            sage: graphs.PetersenGraph().girth()
            5
            sage: graphs.HeawoodGraph().girth()
            6
            sage: graphs.trees(9).next().girth()
            +Infinity

        """
        G = self.to_undirected()
        G.relabel() # vertices are now {0, ..., n-1}
        n = G.num_verts()
        best = n+1
        for i in range(n-2):
            span = [0]*n
            span[i] = 1
            depth = 1
            thisList = [i]
            while 2*depth <= best and 3 < best:
                nextList = []
                for v in thisList:
                    for u in G.neighbors(v):
                        if not span[u]:
                            span[u] = 1
                            nextList.append(u)
                        else:
                            if u in thisList:
                                best = depth*2-1
                                break
                            if u in nextList:
                                best = depth*2
                thisList = nextList
                depth += 1
        if best == n+1:
            from sage.rings.infinity import Infinity
            return Infinity
        return best



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
            sage: G = Graph(sparse=True)
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

    def interior_paths(self, start, end):
        """
        Returns an exhaustive list of paths (also lists) through
        only interior vertices from vertex start to vertex end in the
        (di)graph.

        Note -- start and end do not necessarily have to be boundary
                vertices.

        INPUT:
            start -- the vertex of the graph to search for paths from
            end -- the vertex of the graph to search for paths to

        EXAMPLES:
            sage: eg1 = Graph({0:[1,2], 1:[4], 2:[3,4], 4:[5], 5:[6]})
            sage: sorted(eg1.all_paths(0,6))
            [[0, 1, 4, 5, 6], [0, 2, 4, 5, 6]]
            sage: eg2 = eg1.copy()
            sage: eg2.set_boundary([0,1,3])
            sage: sorted(eg2.interior_paths(0,6))
            [[0, 2, 4, 5, 6]]
            sage: sorted(eg2.all_paths(0,6))
            [[0, 1, 4, 5, 6], [0, 2, 4, 5, 6]]
            sage: eg3 = graphs.PetersenGraph()
            sage: eg3.set_boundary([0,1,2,3,4])
            sage: sorted(eg3.all_paths(1,4))
            [[1, 0, 4],
             [1, 0, 5, 7, 2, 3, 4],
             [1, 0, 5, 7, 2, 3, 8, 6, 9, 4],
             [1, 0, 5, 7, 9, 4],
             [1, 0, 5, 7, 9, 6, 8, 3, 4],
             [1, 0, 5, 8, 3, 2, 7, 9, 4],
             [1, 0, 5, 8, 3, 4],
             [1, 0, 5, 8, 6, 9, 4],
             [1, 0, 5, 8, 6, 9, 7, 2, 3, 4],
             [1, 2, 3, 4],
             [1, 2, 3, 8, 5, 0, 4],
             [1, 2, 3, 8, 5, 7, 9, 4],
             [1, 2, 3, 8, 6, 9, 4],
             [1, 2, 3, 8, 6, 9, 7, 5, 0, 4],
             [1, 2, 7, 5, 0, 4],
             [1, 2, 7, 5, 8, 3, 4],
             [1, 2, 7, 5, 8, 6, 9, 4],
             [1, 2, 7, 9, 4],
             [1, 2, 7, 9, 6, 8, 3, 4],
             [1, 2, 7, 9, 6, 8, 5, 0, 4],
             [1, 6, 8, 3, 2, 7, 5, 0, 4],
             [1, 6, 8, 3, 2, 7, 9, 4],
             [1, 6, 8, 3, 4],
             [1, 6, 8, 5, 0, 4],
             [1, 6, 8, 5, 7, 2, 3, 4],
             [1, 6, 8, 5, 7, 9, 4],
             [1, 6, 9, 4],
             [1, 6, 9, 7, 2, 3, 4],
             [1, 6, 9, 7, 2, 3, 8, 5, 0, 4],
             [1, 6, 9, 7, 5, 0, 4],
             [1, 6, 9, 7, 5, 8, 3, 4]]
            sage: sorted(eg3.interior_paths(1,4))
            [[1, 6, 8, 5, 7, 9, 4], [1, 6, 9, 4]]
            sage: dg = DiGraph({0:[1,3,4], 1:[3], 2:[0,3,4],4:[3]}, boundary=[4])
            sage: sorted(dg.all_paths(0,3))
            [[0, 1, 3], [0, 3], [0, 4, 3]]
            sage: sorted(dg.interior_paths(0,3))
            [[0, 1, 3], [0, 3]]
            sage: ug = dg.to_undirected()
            sage: sorted(ug.all_paths(0,3))
            [[0, 1, 3], [0, 2, 3], [0, 2, 4, 3], [0, 3], [0, 4, 2, 3], [0, 4, 3]]
            sage: sorted(ug.interior_paths(0,3))
            [[0, 1, 3], [0, 2, 3], [0, 3]]
        """
        H = self.copy()
        for vertex in self.get_boundary():
            if (vertex != start and vertex != end):
                H.delete_edges(H.edges_incident(vertex))
        return H.all_paths(start, end)

    def all_paths(self, start, end):
        """
        Returns a list of all paths (also lists) between a pair of
        vertices (start, end) in the (di)graph.

        EXAMPLES:
            sage: eg1 = Graph({0:[1,2], 1:[4], 2:[3,4], 4:[5], 5:[6]})
            sage: eg1.all_paths(0,6)
            [[0, 1, 4, 5, 6], [0, 2, 4, 5, 6]]
            sage: eg2 = graphs.PetersenGraph()
            sage: sorted(eg2.all_paths(1,4))
            [[1, 0, 4],
             [1, 0, 5, 7, 2, 3, 4],
             [1, 0, 5, 7, 2, 3, 8, 6, 9, 4],
             [1, 0, 5, 7, 9, 4],
             [1, 0, 5, 7, 9, 6, 8, 3, 4],
             [1, 0, 5, 8, 3, 2, 7, 9, 4],
             [1, 0, 5, 8, 3, 4],
             [1, 0, 5, 8, 6, 9, 4],
             [1, 0, 5, 8, 6, 9, 7, 2, 3, 4],
             [1, 2, 3, 4],
             [1, 2, 3, 8, 5, 0, 4],
             [1, 2, 3, 8, 5, 7, 9, 4],
             [1, 2, 3, 8, 6, 9, 4],
             [1, 2, 3, 8, 6, 9, 7, 5, 0, 4],
             [1, 2, 7, 5, 0, 4],
             [1, 2, 7, 5, 8, 3, 4],
             [1, 2, 7, 5, 8, 6, 9, 4],
             [1, 2, 7, 9, 4],
             [1, 2, 7, 9, 6, 8, 3, 4],
             [1, 2, 7, 9, 6, 8, 5, 0, 4],
             [1, 6, 8, 3, 2, 7, 5, 0, 4],
             [1, 6, 8, 3, 2, 7, 9, 4],
             [1, 6, 8, 3, 4],
             [1, 6, 8, 5, 0, 4],
             [1, 6, 8, 5, 7, 2, 3, 4],
             [1, 6, 8, 5, 7, 9, 4],
             [1, 6, 9, 4],
             [1, 6, 9, 7, 2, 3, 4],
             [1, 6, 9, 7, 2, 3, 8, 5, 0, 4],
             [1, 6, 9, 7, 5, 0, 4],
             [1, 6, 9, 7, 5, 8, 3, 4]]
            sage: dg = DiGraph({0:[1,3], 1:[3], 2:[0,3]})
            sage: sorted(dg.all_paths(0,3))
            [[0, 1, 3], [0, 3]]
            sage: ug = dg.to_undirected()
            sage: sorted(ug.all_paths(0,3))
            [[0, 1, 3], [0, 2, 3], [0, 3]]
        """
        if self.is_directed():
            iterator=self.successor_iterator
        else:
            iterator=self.neighbor_iterator
        all_paths = []      # list of
        act_path = []       # the current path
        act_path_iter = []  # the neighbor/successor-iterators of the current path
        done = False
        s=start
        while not done:
            if s==end:      # if path completes, add to list
                all_paths.append(act_path+[s])
            else:
                if s not in act_path:   # we want vertices just once in a path
                    act_path.append(s)  # extend current path
                    act_path_iter.append(iterator(s))  # save the state of the neighbor/successor-iterator of the current vertex
            s=None
            while (s is None) and not done:
                try:
                    s=act_path_iter[-1].next()  # try to get the next neighbor/successor, ...
                except (StopIteration):         # ... if there is none ...
                    act_path.pop()              # ... go one step back
                    act_path_iter.pop()
                if len(act_path)==0:            # there is no other vertex ...
                    done = True                 # ... so we are done
        return all_paths


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
            sage: D.delete_edges(D.edges_incident(13))
            sage: D.shortest_path(13, 4)
            []
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse=True )
            sage: G.plot(edge_labels=True).show()
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
                    L = networkx.bidirectional_dijkstra(self.networkx_graph(copy=False), u, v)[1]
                except:
                    L = False
            else:
                L = networkx.dijkstra_path(self.networkx_graph(copy=False), u, v)
        else:
            if bidirectional:
                L = networkx.shortest_path(self.networkx_graph(copy=False), u, v)
            else:
                try:
                    L = networkx.single_source_shortest_path(self.networkx_graph(copy=False), u)[v]
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
            sage: D.delete_edges(D.edges_incident(13))
            sage: D.shortest_path_length(13, 4)
            +Infinity
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse=True )
            sage: G.plot(edge_labels=True).show()
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
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse=True )
            sage: G.plot(edge_labels=True).show()
            sage: G.shortest_paths(0, by_weight=True)
            {0: [0], 1: [0, 1], 2: [0, 1, 2], 3: [0, 1, 2, 3], 4: [0, 4]}

        """
        import networkx
        if by_weight:
            return networkx.single_source_dijkstra_path(self.networkx_graph(copy=False), u)
        else:
            return networkx.single_source_shortest_path(self.networkx_graph(copy=False), u, cutoff)

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
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse=True )
            sage: G.plot(edge_labels=True).show()
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

    def shortest_path_all_pairs(self, by_weight=True, default_weight=1):
        """
        Uses the Floyd-Warshall algorithm to find a shortest weighted
        path for each pair of vertices.

        The weights (labels) on the vertices can be anything that can
        be compared and can be summed.

        INPUT:

            by_weight -- If False, figure distances by the numbers of edges.
            default_weight -- (defaults to 1) The default weight to
            assign edges that don't have a weight (i.e., a label).

        OUTPUT:
            A tuple (dist, pred). They are both dicts of dicts. The first
        indicates the length dist[u][v] of the shortest weighted path from u
        to v. The second is more complicated-- it indicates the predecessor
        pred[u][v] of v in the shortest path from u to v.

        EXAMPLE:
            sage: G = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse=True )
            sage: G.plot(edge_labels=True).show()
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

            sage: G = Graph( { 0: {1:None}, 1: {2:None}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, sparse=True )
            sage: G.shortest_path_all_pairs(by_weight=False)
            ({0: {0: 0, 1: 1, 2: 2, 3: 2, 4: 1},
            1: {0: 1, 1: 0, 2: 1, 3: 2, 4: 2},
            2: {0: 2, 1: 1, 2: 0, 3: 1, 4: 2},
            3: {0: 2, 1: 2, 2: 1, 3: 0, 4: 1},
            4: {0: 1, 1: 2, 2: 2, 3: 1, 4: 0}},
            {0: {0: None, 1: 0, 2: 1, 3: 4, 4: 0},
            1: {0: 1, 1: None, 2: 1, 3: 2, 4: 0},
            2: {0: 1, 1: 2, 2: None, 3: 2, 4: 3},
            3: {0: 4, 1: 2, 2: 3, 3: None, 4: 3},
            4: {0: 4, 1: 0, 2: 3, 3: 4, 4: None}})
            sage: G.shortest_path_all_pairs()
            ({0: {0: 0, 1: 1, 2: 2, 3: 3, 4: 2},
            1: {0: 1, 1: 0, 2: 1, 3: 2, 4: 3},
            2: {0: 2, 1: 1, 2: 0, 3: 1, 4: 3},
            3: {0: 3, 1: 2, 2: 1, 3: 0, 4: 2},
            4: {0: 2, 1: 3, 2: 3, 3: 2, 4: 0}},
            {0: {0: None, 1: 0, 2: 1, 3: 2, 4: 0},
            1: {0: 1, 1: None, 2: 1, 3: 2, 4: 0},
            2: {0: 1, 1: 2, 2: None, 3: 2, 4: 3},
            3: {0: 1, 1: 2, 2: 3, 3: None, 4: 3},
            4: {0: 4, 1: 0, 2: 3, 3: 4, 4: None}})
            sage: G.shortest_path_all_pairs(default_weight=200)
            ({0: {0: 0, 1: 200, 2: 5, 3: 4, 4: 2},
            1: {0: 200, 1: 0, 2: 200, 3: 201, 4: 202},
            2: {0: 5, 1: 200, 2: 0, 3: 1, 4: 3},
            3: {0: 4, 1: 201, 2: 1, 3: 0, 4: 2},
            4: {0: 2, 1: 202, 2: 3, 3: 2, 4: 0}},
            {0: {0: None, 1: 0, 2: 3, 3: 4, 4: 0},
            1: {0: 1, 1: None, 2: 1, 3: 2, 4: 0},
            2: {0: 4, 1: 2, 2: None, 3: 2, 4: 3},
            3: {0: 4, 1: 2, 2: 3, 3: None, 4: 3},
            4: {0: 4, 1: 0, 2: 3, 3: 4, 4: None}})

        """
        from sage.rings.infinity import Infinity
        dist = {}
        pred = {}
        verts = self.vertices()
        for u in verts:
            dist[u] = {}
            pred[u] = {}
            for v in verts:
                if self.has_edge(u, v):
                    if by_weight is False:
                        dist[u][v] = 1
                    elif self.edge_label(u, v) is None:
                        dist[u][v] = default_weight
                    else:
                        dist[u][v] = self.edge_label(u, v)
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

    def breadth_first_search(self, u, ignore_direction=False):
        """
        Returns an iterator over vertices in a breadth-first ordering.

        INPUT:
        u -- vertex at which to start search
        ignore_direction -- (default False) only applies to directed graphs. If
            True, searches across edges in either direction.

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

            sage: D = DiGraph({1:[0], 2:[0]})
            sage: list(D.breadth_first_search(0))
            [0]
            sage: list(D.breadth_first_search(0, ignore_direction=True))
            [0, 1, 2]

        """
        # This function is straight from an old version of networkx
        if not self._directed or ignore_direction:
            neighbors=self.neighbor_iterator
        else:
            neighbors=self.successor_iterator
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

    def depth_first_search(self, u, ignore_direction=False):
        """
        Returns an iterator over vertices in a depth-first ordering.

        INPUT:
        u -- vertex at which to start search
        ignore_direction -- (default False) only applies to directed graphs. If
            True, searches across edges in either direction.

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

            sage: D = DiGraph({1:[0], 2:[0]})
            sage: list(D.depth_first_search(0))
            [0]
            sage: list(D.depth_first_search(0, ignore_direction=True))
            [0, 2, 1]

        """
        # This function is straight from an old version of networkx
        if not self._directed or ignore_direction:
            neighbors=self.neighbor_iterator
        else:
            neighbors=self.successor_iterator
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

    ### Constructors

    def add_cycle(self, vertices):
        """
        Adds a cycle to the graph with the given vertices. If the vertices are
        already present, only the edges are added.

        For digraphs, adds the directed cycle, whose orientation is determined
        by the list. Adds edges (vertices[u], vertices[u+1]) and (vertices[-1],
        vertices[0]).

        INPUT:
        vertices -- a list of indices for the vertices of the cycle to be
        added.

        EXAMPLES:
            sage: G = Graph(implementation='networkx')
            sage: G.add_vertices(range(10)); G
            Graph on 10 vertices
            sage: show(G)
            sage: G.add_cycle(range(20)[10:20])
            sage: show(G)
            sage: G.add_cycle(range(10))
            sage: show(G)

            sage: D = DiGraph(implementation='networkx')
            sage: D.add_cycle(range(4))
            sage: D.edges()
            [(0, 1, None), (1, 2, None), (2, 3, None), (3, 0, None)]

        """
        self.add_path(vertices)
        self.add_edge(vertices[-1], vertices[0])

    def add_path(self, vertices):
        """
        Adds a cycle to the graph with the given vertices. If the vertices are
        already present, only the edges are added.

        For digraphs, adds the directed path vertices[0], ..., vertices[-1].

        INPUT:
        vertices -- a list of indices for the vertices of the cycle to be
        added.

        EXAMPLES:
            sage: G = Graph(implementation='networkx')
            sage: G.add_vertices(range(10)); G
            Graph on 10 vertices
            sage: show(G)
            sage: G.add_path(range(20)[10:20])
            sage: show(G)
            sage: G.add_path(range(10))
            sage: show(G)

            sage: D = DiGraph(sparse=True)
            sage: D.add_path(range(4))
            sage: D.edges()
            [(0, 1, None), (1, 2, None), (2, 3, None)]

        """
        vert1 = vertices[0]
        for v in vertices[1:]:
            self.add_edge(vert1, v)
            vert1 = v

    def complement(self):
        """
        Returns the complement of the (di)graph.

        The complement of a graph has the same vertices, but exactly those
        edges that are not in the original graph. This is not well defined for
        graphs with multiple edges.

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: P.plot().show()
            sage: PC = P.complement()
            sage: PC.plot().show()

            sage: graphs.TetrahedralGraph().complement().size()
            0
            sage: graphs.CycleGraph(4).complement().edges()
            [(0, 2, None), (1, 3, None)]
            sage: graphs.CycleGraph(4).complement()
            complement(Cycle graph): Graph on 4 vertices
            sage: Graph(multiedges=True).complement()
            Traceback (most recent call last):
            ...
            TypeError: Complement not well defined for (di)graphs with multiple edges.

        """
        if self.multiple_edges():
            raise TypeError('Complement not well defined for (di)graphs with multiple edges.')
        G = self.copy()
        G.delete_edges(G.edges())
        G.name('complement(%s)'%self.name())
        for u in self:
            for v in self:
                if not self.has_edge(u,v):
                    G.add_edge(u,v)
        return G

    def line_graph(self, labels=True):
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
            sage: h2=g.line_graph(labels=False)
            sage: h2.vertices()
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: h2.am()==h.am()
            True
            sage: g = DiGraph([[1..4],lambda i,j: i<j], implementation='networkx')
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
        if self._directed:
            G=DiGraph(implementation='networkx')
            G.add_vertices(self.edges(labels=labels))
            for v in self:
                # Connect appropriate incident edges of the vertex v
                G.add_edges([(e,f) for e in self.incoming_edge_iterator(v, labels=labels) \
                             for f in self.outgoing_edge_iterator(v, labels=labels)])
            return G
        else:
            G=Graph(implementation='networkx')
            # We must sort the edges' endpoints so that (1,2,None) is
            # seen as the same edge as (2,1,None).
            if labels:
                elist=[(min(i[0:2]),max(i[0:2]),i[2])
                       for i in self.edge_iterator()]
            else:
                elist=[(min(i),max(i))
                       for i in self.edge_iterator(labels=False)]
            G.add_vertices(elist)
            for v in self:
                if labels:
                    elist=[(min(i[0:2]),max(i[0:2]),i[2])
                           for i in self.edge_iterator(v)]
                else:
                    elist=[(min(i),max(i))
                           for i in self.edge_iterator(v, labels=False)]
                G.add_edges([(e, f) for e in elist for f in elist])
            return G

    def to_simple(self):
        """
        Returns a simple version of itself (i.e., undirected and loops
        and multiple edges are removed).

        EXAMPLE:
            sage: G = DiGraph(loops=True,multiedges=True,sparse=True)
            sage: G.add_edges( [ (0,0), (1,1), (2,2), (2,3,1), (2,3,2), (3,2) ] )
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

    def disjoint_union(self, other, verbose_relabel=True):
        """
        Returns the disjoint union of self and other.

        If the graphs have common vertices, the vertices will be
        renamed to form disjoint sets.

        INPUT:
            verbose_relabel -- (defaults to True) If True and the
            graphs have common vertices, then each vertex v in the
            first graph will be changed to '0,v' and each vertex u in
            the second graph will be changed to '1,u'.  If False, the
            vertices of the first graph and the second graph will be
            relabeled with consecutive integers.

        EXAMPLE:
            sage: G = graphs.CycleGraph(3)
            sage: H = graphs.CycleGraph(4)
            sage: J = G.disjoint_union(H); J
            Cycle graph disjoint_union Cycle graph: Graph on 7 vertices
            sage: J.vertices()
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (1, 3)]
            sage: J = G.disjoint_union(H, verbose_relabel=False); J
            Cycle graph disjoint_union Cycle graph: Graph on 7 vertices
            sage: J.vertices()
            [0, 1, 2, 3, 4, 5, 6]

        If the vertices are already disjoint and verbose_relabel is
        True, then the vertices are not relabeled.
            sage: G=Graph({'a': ['b']}, implementation='networkx')
            sage: G.name("Custom path")
            sage: G.name()
            'Custom path'
            sage: H=graphs.CycleGraph(3)
            sage: J=G.disjoint_union(H); J
            Custom path disjoint_union Cycle graph: Graph on 5 vertices
            sage: J.vertices()
            [0, 1, 2, 'a', 'b']

        """
        if (self._directed and not other._directed) or (not self._directed and other._directed):
            raise TypeError('Both arguments must be of the same class.')

        if not verbose_relabel:
            r_self = {}; r_other = {}; i = 0
            for v in self:
                r_self[v] = i; i += 1
            for v in other:
                r_other[v] = i; i += 1
            G = self.relabel(r_self, inplace=False).union(other.relabel(r_other, inplace=False))
        elif any(u==v for u in self for v in other):
            r_self = dict([[v,(0,v)] for v in self])
            r_other = dict([[v,(1,v)] for v in other])
            G = self.relabel(r_self, inplace=False).union(other.relabel(r_other, inplace=False))
        else:
            G = self.union(other)

        G.name('%s disjoint_union %s'%(self.name(), other.name()))
        return G

    def union(self, other):
        """
        Returns the union of self and other.

        If the graphs have common vertices, the common vertices will
        be identified.

        EXAMPLE:
            sage: G = graphs.CycleGraph(3)
            sage: H = graphs.CycleGraph(4)
            sage: J = G.union(H); J
            Graph on 4 vertices
            sage: J.vertices()
            [0, 1, 2, 3]
            sage: J.edges(labels=False)
            [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)]

        """
        if (self._directed and not other._directed) or (not self._directed and other._directed):
            raise TypeError('Both arguments must be of the same class.')
        if self._directed:
            G = DiGraph()
        else:
            G = Graph(implementation='networkx')
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
            sage: P.plot().show()

            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: C = D.cartesian_product(P); C
            Graph on 200 vertices
            sage: C.plot().show()

        """
        if (self._directed and not other._directed) or (not self._directed and other._directed):
            raise TypeError('Both arguments must be of the same class.')
        if self._directed:
            G = DiGraph()
        else:
            G = Graph(implementation='networkx')
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
        Returns the tensor product, also called the categorical product, of self
        and other.

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
            sage: T.plot().show()

            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: T = D.tensor_product(P); T
            Graph on 200 vertices
            sage: T.plot().show()

        """
        if (self._directed and not other._directed) or (not self._directed and other._directed):
            raise TypeError('Both arguments must be of the same class.')
        if self._directed:
            G = DiGraph()
        else:
            G = Graph(implementation='networkx')
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

    categorical_product = tensor_product

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
            sage: L.plot().show()

            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: L = D.lexicographic_product(P); L
            Graph on 200 vertices
            sage: L.plot().show()

        """
        if (self._directed and not other._directed) or (not self._directed and other._directed):
            raise TypeError('Both arguments must be of the same class.')
        if self._directed:
            G = DiGraph()
        else:
            G = Graph(implementation='networkx')
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
            - (u, w) is an edge of self and (v, x) is an edge of other.
        In other words, the edges of the strong product is the union of the
        edges of the tensor and Cartesian products.

        EXAMPLES:
            sage: Z = graphs.CompleteGraph(2)
            sage: C = graphs.CycleGraph(5)
            sage: S = C.strong_product(Z); S
            Graph on 10 vertices
            sage: S.plot().show()

            sage: D = graphs.DodecahedralGraph()
            sage: P = graphs.PetersenGraph()
            sage: S = D.strong_product(P); S
            Graph on 200 vertices
            sage: S.plot().show()

        """
        if (self._directed and not other._directed) or (not self._directed and other._directed):
            raise TypeError('Both arguments must be of the same class.')
        if self._directed:
            G = DiGraph()
        else:
            G = Graph(implementation='networkx')

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
            sage: D.plot().show()

            sage: C = graphs.CycleGraph(5)
            sage: D = C.disjunctive_product(Z); D
            Graph on 10 vertices
            sage: D.plot().show()

        """
        if (self._directed and not other._directed) or (not self._directed and other._directed):
            raise TypeError('Both arguments must be of the same class.')
        if self._directed:
            G = DiGraph()
        else:
            G = Graph(implementation='networkx')
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

    def transitive_closure(self):
        r"""
        Computes the transitive closure of a graph and returns it.
        The original graph is not modified.

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
            sage: g=DiGraph({0:[1,2], 1:[3], 2:[4,5]})
            sage: g.transitive_closure().edges(labels=False)
            [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 3), (2, 4), (2, 5)]

        """
        G = self.copy()
        for v in G:
            # todo optimization opportunity: we are adding edges that
            # are already in the graph and we are adding edges
            # one at a time.
            for e in G.breadth_first_search(v):
                G.add_edge((v,e))
        return G

    def transitive_reduction(self):
        r"""
        Returns a transitive reduction of a graph.  The original graph
        is not modified.

        A transitive reduction H of G has a path from x to y if and
        only if there was a path from x to y in G.  Deleting any edge
        of H destroys this property.  A transitive reduction is not
        unique in general.  A transitive reduction has the same
        transitive closure as the original graph.

        A transitive reduction of a complete graph is a tree.  A
        transitive reduction of a tree is itself.


        EXAMPLES:
            sage: g=graphs.PathGraph(4)
            sage: g.transitive_reduction()==g
            True
            sage: g=graphs.CompleteGraph(5)
            sage: edges = g.transitive_reduction().edges(); len(edges)
            4
            sage: g=DiGraph({0:[1,2], 1:[2,3,4,5], 2:[4,5]})
            sage: g.transitive_reduction().size()
            5

        """
        from sage.rings.infinity import Infinity
        G = self.copy()
        for e in self.edge_iterator():
            # Try deleting the edge, see if we still have a path
            # between the vertices.
            G.delete_edge(e)
            if G.distance(e[0],e[1])==Infinity:
                # oops, we shouldn't have deleted it
                G.add_edge(e)
        return G

    ### Visualization

    def _color_by_label(self, format='hex'):
        """
        Logic for coloring by label (factored out from plot() for use
        in 3d plots, etc)

        EXAMPLE:
            sage: G = AlternatingGroup(5).cayley_graph()
            sage: G.num_edges()
            120
            sage: G._color_by_label()
            {'#00ffff': [((1,4)(3,5), (1,5,4), (3,4,5)),
             ...],
             '#ff0000': [((1,4)(3,5), (1,5,4,2,3), (1,2,3,4,5)),
             ...]}
        """
        from sage.plot.plot import rainbow
        edge_labels = []
        for e in self.edge_iterator():
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
            scaling_term=0.05, iterations=50, loop_size=.1, talk=False,
            color_by_label=False, heights=None, edge_style=None, save_pos=False,
            tree_root=None, tree_orientation="down"):
        """
        Returns a graphics object representing the (di)graph.

        INPUT:
            pos -- an optional positioning dictionary
            layout -- what kind of layout to use, takes precedence over pos
                'circular' -- plots the graph with vertices evenly distributed
                    on a circle
                'spring' -- uses the traditional spring layout, using the
                    graph's current positions as initial positions
                'tree' -- the (di)graph must be a tree. One can specify the root
                    of the tree using the keyword tree_root, otherwise a root
                    will be selected at random. Then the tree will be plotted in
                    levels, depending on minimum distance for the root.
            vertex_labels -- whether to print vertex labels
            edge_labels -- whether to print edge labels. By default, False,
                but if True, the result of str(l) is printed on the edge for
                each label l. Labels equal to None are not printed (to set edge
                labels, see set_edge_label).
            vertex_size -- size of vertices displayed
            graph_border -- whether to include a box around the graph
            vertex_colors -- optional dictionary to specify vertex colors: each
                key is a color recognizable by matplotlib, and each corresponding
                entry is a list of vertices. If a vertex is not listed, it looks
                invisible on the resulting plot (it doesn't get drawn).
            edge_colors -- a dictionary specifying edge colors: each key is a
                color recognized by matplotlib, and each entry is a list of edges.
            partition -- a partition of the vertex set. if specified, plot will
                show each cell in a different color. vertex_colors takes precedence.
            scaling_term -- default is 0.05. if vertices are getting chopped off,
                increase; if graph is too small, decrease. should be positive, but
                values much bigger than 1/8 won't be useful unless the vertices
                are huge
            talk -- if true, prints large vertices with white backgrounds so that
                labels are legible on slides
            iterations -- how many iterations of the spring layout algorithm to
                go through, if applicable
            color_by_label -- if True, color edges by their labels
            heights -- if specified, this is a dictionary from a set of
                floating point heights to a set of vertices
            edge_style -- keyword arguments passed into the
                edge-drawing routine.  This currently only works for
                directed graphs, since we pass off the undirected graph to
                networkx
            tree_root -- a vertex of the tree to be used as the root for
                the layout="tree" option. If no root is specified, then one
                is chosen at random. Ignored unless layout='tree'.
            tree_orientation -- "up" or "down" (default is "down").
                If "up" (resp., "down"), then the root of the tree will
                appear on the bottom (resp., top) and the tree will grow
                upwards (resp. downwards). Ignored unless layout='tree'.
            save_pos -- save position computed during plotting

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
            sage: pl.show()

            sage: C = graphs.CubeGraph(8)
            sage: P = C.plot(vertex_labels=False, vertex_size=0, graph_border=True)
            sage: P.show()

            sage: G = graphs.HeawoodGraph()
            sage: for u,v,l in G.edges():
            ...    G.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: G.plot(edge_labels=True).show()

            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []}, implementation='networkx' )
            sage: for u,v,l in D.edges():
            ...    D.set_edge_label(u,v,'(' + str(u) + ',' + str(v) + ')')
            sage: D.plot(edge_labels=True, layout='circular').show()

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
            sage: C.plot(vertex_labels=False, vertex_size=0, edge_colors=edge_colors).show()

            sage: D = graphs.DodecahedralGraph()
            sage: Pi = [[6,5,15,14,7],[16,13,8,2,4],[12,17,9,3,1],[0,19,18,10,11]]
            sage: D.show(partition=Pi)

            sage: G = graphs.PetersenGraph()
            sage: G.loops(True)
            sage: G.add_edge(0,0)
            sage: G.show()

            sage: D = DiGraph({0:[0,1], 1:[2], 2:[3]}, loops=True)
            sage: D.show()
            sage: D.show(edge_colors={(0,1,0):[(0,1,None),(1,2,None)],(0,0,0):[(2,3,None)]})

            sage: pos = {0:[0.0, 1.5], 1:[-0.8, 0.3], 2:[-0.6, -0.8], 3:[0.6, -0.8], 4:[0.8, 0.3]}
            sage: g = Graph({0:[1], 1:[2], 2:[3], 3:[4], 4:[0]})
            sage: g.plot(pos=pos, layout='spring', iterations=0)

            sage: G = Graph()
            sage: P = G.plot()
            sage: P.axes()
            False
            sage: G = DiGraph()
            sage: P = G.plot()
            sage: P.axes()
            False

            sage: G = graphs.PetersenGraph()
            sage: G.get_pos()
            {0: [6.12..., 1.0...],
             1: [-0.95..., 0.30...],
             2: [-0.58..., -0.80...],
             3: [0.58..., -0.80...],
             4: [0.95..., 0.30...],
             5: [1.53..., 0.5...],
             6: [-0.47..., 0.15...],
             7: [-0.29..., -0.40...],
             8: [0.29..., -0.40...],
             9: [0.47..., 0.15...]}
            sage: P = G.plot(save_pos=True, layout='spring')
            sage: G.get_pos()
            {0: [-0.39..., 0.06...],
             1: [-0.26..., 0.94...],
             2: [-0.29..., 0.43...],
             3: [-0.40..., -0.70...],
             4: [-0.88..., -0.46...],
             5: [0.75..., -0.1...],
             6: [0.32..., 0.28...],
             7: [0.67..., 0.52...],
             8: [0.44..., -0.72...],
             9: [0.05..., -0.19...]}

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.plot(heights={0:[0], 1:[4,5,1], 2:[2], 3:[3,6]})

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.plot(heights={0:[0], 1:[4,5,1], 2:[2], 3:[3,6]})
            sage: t.set_edge_label(0,1,-7)
            sage: t.set_edge_label(0,5,3)
            sage: t.set_edge_label(0,5,99)
            sage: t.set_edge_label(1,2,1000)
            sage: t.set_edge_label(3,2,'spam')
            sage: t.set_edge_label(2,6,3/2)
            sage: t.set_edge_label(0,4,66)
            sage: t.plot(heights={0:[0], 1:[4,5,1], 2:[2], 3:[3,6]}, edge_labels=True)

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.plot(layout='tree')

            sage: t = DiGraph('JCC???@A??GO??CO??GO??')
            sage: t.plot(layout='tree', tree_root=0, tree_orientation="up")

        """
        if edge_style is None:
            edge_style={}
        if talk:
            vertex_size = 500
            if partition is None:
                vertex_colors = {'#FFFFFF':self.vertices()}
        from sage.plot.plot import networkx_plot, Graphics, rainbow
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
                vertex_colors['#fec7b8'] = int_verts
                vertex_colors['#b3e8ff'] = bdy_verts
            else:
                vertex_colors={'#fec7b8': self.vertices()}
        if pos is None and layout is None:
            if heights is None:
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
        elif layout == 'tree':
            if not self.is_tree():
                raise RuntimeError("Cannot use tree layout on this graph: self.is_tree() returns False.")
            verts = self.vertices()
            if tree_root is None:
                from sage.misc.prandom import randrange
                root = verts[randrange(self.num_verts())]
            else:
                root = tree_root
            # BFS search for heights
            seen = [root]
            queue = [root]
            heights = [-1]*self.num_verts()
            heights[verts.index(root)] = 0
            while queue:
                u = queue.pop(0)
                for v in self.neighbors(u):
                    if v not in seen:
                        seen.append(v)
                        queue.append(v)
                        heights[verts.index(v)] = heights[verts.index(u)] + 1
            if tree_orientation == 'down':
                maxx = max(heights)
                heights = [maxx - heights[i] for i in xrange(self.num_verts())]
            heights_dict = {}
            for v in self:
                if not heights_dict.has_key(heights[verts.index(v)]):
                    heights_dict[heights[verts.index(v)]] = [v]
                else:
                    heights_dict[heights[verts.index(v)]].append(v)
            heights = heights_dict
        if heights is not None and pos is None and self.num_verts() > 0:
            pos = {}
            mmax = max([len(ccc) for ccc in heights.values()])
            ymin = min(heights.keys())
            ymax = max(heights.keys())
            dist = ((ymax-ymin)/(mmax+1.0))
            for height in heights:
                num_xes = len(heights[height])
                if num_xes == 0: continue
                j = (mmax - num_xes)/2.0
                for k in range(num_xes):
                    pos[heights[height][k]] = [ dist * (j+k+1) + random()*(dist*0.03), height ]
        if pos is None or layout == 'spring' or heights is not None:
            pos = graph_fast.spring_layout_fast(self, iterations=iterations, vpos=pos, height=(heights is not None))
        else:
            for v in pos:
                for a in range(len(pos[v])):
                    pos[v][a] = float(pos[v][a])

        if color_by_label:
            edge_colors = self._color_by_label()

        G = networkx_plot(self.networkx_graph(copy=False), pos=pos, vertex_labels=vertex_labels, \
          vertex_size=vertex_size, vertex_colors=vertex_colors, \
          edge_colors=edge_colors, graph_border=graph_border, \
          scaling_term=scaling_term, draw_edges=(not self._directed))
        if save_pos:
            self.set_pos(pos)
        if self._directed:
            from sage.plot.plot import arrow
            P = Graphics()
            from math import sqrt, pi
            vertex_radius = sqrt(vertex_size/pi)
            edge_style.setdefault('rgbcolor', (0,0,0))
            edge_style.setdefault('arrowshorten', vertex_radius*2)
            edge_style.setdefault('arrowsize', 10)

            if edge_colors is None:
                for u,v,_ in self.edge_iterator():
                    if u != v:
                        P += arrow((pos[u][0],pos[u][1]),(pos[v][0],pos[v][1]), **edge_style)
            else:
                for color in edge_colors:
                    for u,v,_ in edge_colors[color]:
                        if u != v:
                            this_edge=edge_style.copy()
                            this_edge['rgbcolor']=color
                            P += arrow((pos[u][0],pos[u][1]),(pos[v][0],pos[v][1]), **this_edge)
            limits = (G.xmin(), G.xmax(), G.ymin(), G.ymax())
            G = P + G
            G.xmin(limits[0])
            G.xmax(limits[1])
            G.ymin(limits[2])
            G.ymax(limits[3])
        if edge_labels:
            from sage.plot.plot import text
            K = Graphics()
            for u,v,l in self.edges():
                if not l is None:
                    K += text(str(l), [(pos[u][0] + pos[v][0])/2, (pos[u][1] + pos[v][1])/2])
            K.axes_range(xmin=G.xmin(), xmax=G.xmax(), ymin=G.ymin(), ymax=G.ymax())
            G += K
        if self.loops():
            from sage.plot.plot import circle
            L = []
            for v in self.loop_vertices():
                L.append(circle((pos[v][0],pos[v][1]-loop_size), loop_size, rgbcolor=(0,0,0)))
            G = sum(L) + G
        G.axes(False)
        return G

    def show(self, **kwds):
        """
        Shows the (di)graph.

        For syntax and lengthy documentation, see G.plot?. Any options not used by
        plot will be passed on to the Graphics.show method.

        EXAMPLE:
            sage: C = graphs.CubeGraph(8)
            sage: P = C.plot(vertex_labels=False, vertex_size=0, graph_border=True)
            sage: P.show()

        """
        kwds.setdefault('figsize', [4,4])
        import inspect
        vars = inspect.getargspec(GenericGraph.plot)[0]
        plot_kwds = {}
        for kwd in vars:
            if kwds.has_key(kwd):
                plot_kwds[kwd] = kwds.pop(kwd)
        self.plot(**plot_kwds).show(**kwds)

    def plot3d(self, bgcolor=(1,1,1), vertex_colors=None, vertex_size=0.06,
                     edge_colors=None, edge_size=0.02, edge_size2=0.0325,
                     pos3d=None, iterations=50, color_by_label=False,
                     engine='jmol', **kwds):
        """
        Plot a graph in three dimensions.

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
            edge_size2 -- float (default: 0.0325), used for Tachyon sleeves
            pos3d -- a position dictionary for the vertices
            iterations -- how many iterations of the spring layout algorithm to
                go through, if applicable
            engine -- which renderer to use. Options:
                'jmol' -- default
                'tachyon'
            xres -- resolution
            yres -- resolution
            **kwds -- passed on to the rendering engine

        EXAMPLES:
            sage: G = graphs.CubeGraph(5)
            sage: G.plot3d(iterations=500, edge_size=None, vertex_size=0.04)

        We plot a fairly complicated Cayley graph:
            sage: A5 = AlternatingGroup(5); A5
            Alternating group of order 5!/2 as a permutation group
            sage: G = A5.cayley_graph()
            sage: G.plot3d(vertex_size=0.03, edge_size=0.01, vertex_colors={(1,1,1):G.vertices()}, bgcolor=(0,0,0), color_by_label=True, iterations=200)

        Some Tachyon examples:
            sage: D = graphs.DodecahedralGraph()
            sage: P3D = D.plot3d(engine='tachyon')
            sage: P3D.show() # long time

            sage: G = graphs.PetersenGraph()
            sage: G.plot3d(engine='tachyon', vertex_colors={(0,0,1):G.vertices()}).show() # long time

            sage: C = graphs.CubeGraph(4)
            sage: C.plot3d(engine='tachyon', edge_colors={(0,1,0):C.edges()}, vertex_colors={(1,1,1):C.vertices()}, bgcolor=(0,0,0)).show() # long time

            sage: K = graphs.CompleteGraph(3)
            sage: K.plot3d(engine='tachyon', edge_colors={(1,0,0):[(0,1,None)], (0,1,0):[(0,2,None)], (0,0,1):[(1,2,None)]}).show() # long time

        A directed version of the dodecahedron
            sage: D = DiGraph( { 0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4], 4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13], 13: [14], 14: [15], 15: [16], 16: [17], 17: [18], 18: [19], 19: []} )
            sage: D.plot3d().show() # long time

            sage: P = graphs.PetersenGraph().to_directed()
            sage: from sage.plot.plot import rainbow
            sage: edges = P.edges()
            sage: R = rainbow(len(edges), 'rgbtuple')
            sage: edge_colors = {}
            sage: for i in range(len(edges)):
            ...       edge_colors[R[i]] = [edges[i]]
            sage: P.plot3d(engine='tachyon', edge_colors=edge_colors).show() # long time

        """
        if engine == 'jmol':
            from sage.plot.plot3d.all import sphere, line3d, arrow3d
            kwds.setdefault('aspect_ratio', [1,1,1])
            verts = self.vertices()

            if vertex_colors is None:
                vertex_colors = { (1,0,0) : verts }
            if pos3d is None:
                pos3d = graph_fast.spring_layout_fast(self, dim=3, iterations=iterations)

            if color_by_label:
                if edge_colors is  None:
                        # do the coloring
                        edge_colors = self._color_by_label(format='rgbtuple')
            elif edge_colors is None:
                edge_colors = { (0,0,0) : self.edges() }

            # by default turn off the frame
            if not kwds.has_key('frame'):
                kwds['frame'] = False
            # by default make the background given by bgcolor
            if not kwds.has_key('background'):
                kwds['background'] = bgcolor
            try:
                graphic = 0
                for color in vertex_colors:
                    for v in vertex_colors[color]:
                        graphic += sphere(center=pos3d[v], size=vertex_size, color=color, **kwds)
                if self._directed:
                    for color in edge_colors:
                        for u, v, l in edge_colors[color]:
                            graphic += arrow3d(pos3d[u], pos3d[v], radius=edge_size, color=color, closed=False, **kwds)

                else:
                    for color in edge_colors:
                        for u, v, l in edge_colors[color]:
                            graphic += line3d([pos3d[u], pos3d[v]], radius=edge_size, color=color, closed=False, **kwds)

                return graphic

            except KeyError:
                raise KeyError, "Oops! You haven't specified positions for all the vertices."

        elif engine == 'tachyon':
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
                if self._directed:
                    for u,v,l in edge_colors[color]:
                        TT.fcylinder( (pos3d[u][0],pos3d[u][1],pos3d[u][2]),
                                      (pos3d[v][0],pos3d[v][1],pos3d[v][2]), edge_size,'edge_color_%d'%i)
                        TT.fcylinder( (0.25*pos3d[u][0] + 0.75*pos3d[v][0],
                                       0.25*pos3d[u][1] + 0.75*pos3d[v][1],
                                       0.25*pos3d[u][2] + 0.75*pos3d[v][2],),
                                      (pos3d[v][0],pos3d[v][1],pos3d[v][2]), edge_size2,'edge_color_%d'%i)
                else:
                    for u, v, l in edge_colors[color]:
                        TT.fcylinder( (pos3d[u][0],pos3d[u][1],pos3d[u][2]), (pos3d[v][0],pos3d[v][1],pos3d[v][2]), edge_size,'edge_color_%d'%i)

            return TT

        else:
            raise TypeError("Rendering engine (%s) not implemented."%engine)

    def show3d(self, bgcolor=(1,1,1), vertex_colors=None, vertex_size=0.06,
                     edge_colors=None, edge_size=0.02, edge_size2=0.0325,
                     pos3d=None, iterations=50, color_by_label=False,
                     engine='jmol', **kwds):
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
            edge_size2 -- float (default: 0.0325), used for Tachyon sleeves
            pos3d -- a position dictionary for the vertices
            iterations -- how many iterations of the spring layout algorithm to
                go through, if applicable
            engine -- which renderer to use. Options:
                'jmol' -- default
                'tachyon'
            xres -- resolution
            yres -- resolution
            **kwds -- passed on to the Tachyon command

        EXAMPLES:
            sage: G = graphs.CubeGraph(5)
            sage: G.show3d(iterations=500, edge_size=None, vertex_size=0.04)

        We plot a fairly complicated Cayley graph:
            sage: A5 = AlternatingGroup(5); A5
            Alternating group of order 5!/2 as a permutation group
            sage: G = A5.cayley_graph()
            sage: G.show3d(vertex_size=0.03, edge_size=0.01, edge_size2=0.02, vertex_colors={(1,1,1):G.vertices()}, bgcolor=(0,0,0), color_by_label=True, iterations=200)

        Some Tachyon examples:
            sage: D = graphs.DodecahedralGraph()
            sage: D.show3d(engine='tachyon') # long time

            sage: G = graphs.PetersenGraph()
            sage: G.show3d(engine='tachyon', vertex_colors={(0,0,1):G.vertices()}) # long time

            sage: C = graphs.CubeGraph(4)
            sage: C.show3d(engine='tachyon', edge_colors={(0,1,0):C.edges()}, vertex_colors={(1,1,1):C.vertices()}, bgcolor=(0,0,0)) # long time

            sage: K = graphs.CompleteGraph(3)
            sage: K.show3d(engine='tachyon', edge_colors={(1,0,0):[(0,1,None)], (0,1,0):[(0,2,None)], (0,0,1):[(1,2,None)]}) # long time

        """
        self.plot3d(bgcolor=bgcolor, vertex_colors=vertex_colors,
                    edge_colors=edge_colors, vertex_size=vertex_size, engine=engine,
                    edge_size=edge_size, iterations=iterations, edge_size2=edge_size2,
                    color_by_label=color_by_label, **kwds).show()

    def _graphviz_string_helper(self, graph_string, edge_string):
        r"""
        Returns a representation in the DOT language, ready to render in graphviz.

        Use \code{graphviz_string} instead.

        INPUT:
            -- graph_string: a string, "graph" for undirected graphs or
               "digraph" for directed graphs.
            -- edge_string: a string, "--" for undirected graphs or "->" for
               directed graphs.

        WARNING:
            Internal function, not for external use!

        REFERENCES:
            http://www.graphviz.org/doc/info/lang.html

        EXAMPLE:
            sage: G = Graph({0:{1:None,2:None}, 1:{0:None,2:None}, 2:{0:None,1:None,3:'foo'}, 3:{2:'foo'}}, implementation='networkx')
            sage: s = G.graphviz_string() # indirect doctest
            sage: s
            'graph {\n"0";"1";"2";"3";\n"0"--"1";"0"--"2";"1"--"2";"2"--"3"[label="foo"];\n}'
        """
        s = '%s {\n' % graph_string
        for v in self.vertex_iterator():
            s+= '"%s";'%v
        s+= '\n'
        for u, v, label in self.edge_iterator():
            if label is None:
                s+= '"%s"%s"%s";' % (u, edge_string, v)
            else:
                s+= '"%s"%s"%s"[label="%s"];' % (u, edge_string, v, label)
        s+= "\n}"
        return s

    def graphviz_string(self):
        r"""
        Returns a representation in the DOT language, ready to render in graphviz.

        EXAMPLES:
            sage: G = Graph({0:{1:None,2:None}, 1:{0:None,2:None}, 2:{0:None,1:None,3:'foo'}, 3:{2:'foo'}}, implementation='networkx')
            sage: s = G.graphviz_string()
            sage: s
            'graph {\n"0";"1";"2";"3";\n"0"--"1";"0"--"2";"1"--"2";"2"--"3"[label="foo"];\n}'
        """
        raise NotImplementedError, "GenericGraph subclasses must override graphviz_string()"

    def graphviz_to_file_named(self, filename):
        r"""
        Write a representation in the DOT language to the named file, ready to
        render in graphviz.

        EXAMPLES:
            sage: G = Graph({0:{1:None,2:None}, 1:{0:None,2:None}, 2:{0:None,1:None,3:'foo'}, 3:{2:'foo'}}, implementation='networkx')
            sage: G.graphviz_to_file_named(os.environ['SAGE_TESTDIR']+'/temp_graphviz')
            sage: open(os.environ['SAGE_TESTDIR']+'/temp_graphviz').read()
            'graph {\n"0";"1";"2";"3";\n"0"--"1";"0"--"2";"1"--"2";"2"--"3"[label="foo"];\n}'
        """
        return open(filename, 'wt').write(self.graphviz_string())

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

            sage: D = P.to_directed()
            sage: D.delete_edge(7,9)
            sage: D.spectrum()
            [-2.0, -2.0, -2.0, -1.7..., 0.8..., 1.0, 1.0, 1.0, 1.0, 2.9...]

        """
        from sage.matrix.constructor import matrix
        from sage.rings.real_double import RDF
        if laplacian:
            M = self.kirchhoff_matrix()
        else:
            M = self.am()
        M = matrix(RDF, M.rows())
        E = M.right_eigenvectors()[0]
        v = [e.real() for e in E]
        v.sort()
        return v

    def characteristic_polynomial(self, var='x', laplacian=False):
        """
        Returns the characteristic polynomial of the adjacency matrix
        of the (di)graph.

        INPUT:
            laplacian -- if True, use the Laplacian matrix instead (see
                self.kirchhoff_matrix())

        EXAMPLE:
            sage: P = graphs.PetersenGraph()
            sage: P.characteristic_polynomial()
            x^10 - 15*x^8 + 75*x^6 - 24*x^5 - 165*x^4 + 120*x^3 + 120*x^2 - 160*x + 48
            sage: P.characteristic_polynomial(laplacian=True)
            x^10 - 30*x^9 + 390*x^8 - 2880*x^7 + 13305*x^6 - 39882*x^5 + 77640*x^4 - 94800*x^3 + 66000*x^2 - 20000*x

        """
        if laplacian:
            return self.kirchhoff_matrix().charpoly(var=var)
        else:
            return self.am().charpoly(var=var)

    def eigenspaces(self, laplacian=False):
        """
        Returns the eigenspaces of the adjacency matrix of the graph.

        INPUT:
            laplacian -- if True, use the Laplacian matrix instead (see
                self.kirchhoff_matrix())

        EXAMPLE:
            sage: C = graphs.CycleGraph(5)
            sage: E = C.eigenspaces()
            sage: E[0][0]
            -1.618...
            sage: E[1][0]  # eigenspace computation is somewhat random
            Vector space of degree 5 and dimension 1 over Real Double Field
            User basis matrix:
            [ 0.632... -0.632... -0.447... -0.013... 0.073...]

            sage: D = C.to_directed()
            sage: F = D.eigenspaces()
            sage: abs(E[0][0] - F[0][0]) < 0.00001
            True

        """
        if laplacian:
            M = self.kirchhoff_matrix()
        else:
            M = self.am()
        from sage.matrix.constructor import matrix
        from sage.rings.real_double import RDF
        M = matrix(RDF, M.rows())
        return M.eigenspaces()

    ### Automorphism and isomorphism

    def relabel(self, perm=None, inplace=True, return_map=False):
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

        If no arguments are provided, the graph is relabeled to be on the
        vertices $\{0,1,...,n-1\}$.

        INPUT:
            inplace -- default is True. If True, modifies the graph and returns
        nothing. If False, returns a relabeled copy of the graph.
            return_map -- default is False. If True, returns the dictionary
        representing the map.

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

        Relabeling to simpler labels:

            sage: G = graphs.CubeGraph(3)
            sage: G.vertices()
            ['000', '001', '010', '011', '100', '101', '110', '111']
            sage: G.relabel()
            sage: G.vertices()
            [0, 1, 2, 3, 4, 5, 6, 7]

            sage: G = graphs.CubeGraph(3)
            sage: expecting = {'000': 0, '001': 1, '010': 2, '011': 3, '100': 4, '101': 5, '110': 6, '111': 7}
            sage: G.relabel(return_map=True) == expecting
            True

        TESTS:
            sage: P = Graph(graphs.PetersenGraph(), sparse=True)
            sage: P.delete_edge([0,1])
            sage: P.add_edge((4,5))
            sage: P.add_edge((2,6))
            sage: P.delete_vertices([0,1])
            sage: P.relabel()

        """
        if perm is None:
            verts = self.vertices() # vertices() returns a sorted list:
            perm = {}; i = 0        # this guarantees consistent relabeling
            for v in verts:
                perm[v] = i
                i += 1
        if not inplace:
            G = self.copy()
            G.relabel(perm)
            return G
        if type(perm) is list:
            perm = dict( [ [i,perm[i]] for i in xrange(len(perm)) ] )
        from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
        if type(perm) is PermutationGroupElement:
            n = self.order()
            ddict = {}
            llist = perm.list()
            for i in xrange(1,n):
                ddict[i] = llist[i-1]%n
            if n > 0:
                ddict[0] = llist[n-1]%n
            perm = ddict
        if type(perm) is not dict:
            raise TypeError("Type of perm is not supported for relabeling.")
        keys = perm.keys()
        verts = self.vertices()
        for v in verts:
            if v not in keys:
                perm[v] = v
        for v in perm.iterkeys():
            if v in verts:
                try:
                    hash(perm[v])
                except TypeError:
                    raise ValueError, "perm dictionary must be of the format {a:a1, b:b1, ...} where a,b,... are vertices and a1,b1,... are hashable"
        self._backend.relabel(perm, self._directed)
        if self._pos:
            new_pos = {}
            for v in perm.iterkeys():
                try:
                    new_pos[perm[v]] = self._pos[v]
                except KeyError:
                    pass
                    # since we allow position dicts to specify
                    # positions for only some of the vertices...
            self._pos = new_pos
        self._boundary = [perm[v] for v in self._boundary]
        if hasattr(self, '_assoc'):
            new_assoc = {}
            for v in self._assoc:
                new_assoc[perm[v]] = self._assoc[v]
            self._assoc = new_assoc
        if hasattr(self, '_embedding'):
            new_embedding = {}
            for v in self._embedding:
                new_embedding[perm[v]] = self._embedding[v]
            self._embedding = new_embedding
        if return_map:
            return perm

    def degree_to_cell(self, vertex, cell):
        """
        Returns the number of edges from vertex to an edge in cell. In the case
        of a digraph, returns a tuple (in_degree, out_degree).

        EXAMPLES:
            sage: G = graphs.CubeGraph(3)
            sage: cell = G.vertices()[:3]
            sage: G.degree_to_cell('011', cell)
            2
            sage: G.degree_to_cell('111', cell)
            0

            sage: D = DiGraph({ 0:[1,2,3], 1:[3,4], 3:[4,5]})
            sage: cell = [0,1,2]
            sage: D.degree_to_cell(5, cell)
            (0, 0)
            sage: D.degree_to_cell(3, cell)
            (2, 0)
            sage: D.degree_to_cell(0, cell)
            (0, 2)

        """
        if self._directed:
            in_neighbors_in_cell = set([a for a,_,_ in self.incoming_edges(vertex)]) & set(cell)
            out_neighbors_in_cell = set([a for _,a,_ in self.outgoing_edges(vertex)]) & set(cell)
            return (len(in_neighbors_in_cell), len(out_neighbors_in_cell))
        else:
            neighbors_in_cell = set(self.neighbors(vertex)) & set(cell)
            return len(neighbors_in_cell)

    def is_equitable(self, partition, quotient_matrix=False):
        """
        Checks whether the given partition is equitable with respect to self.

        A partition is equitable with respect to a graph if for every pair of
        cells C1, C2 of the partition, the number of edges from a vertex of C1
        to C2 is the same, over all vertices in C1.

        INPUT:
            partition -- a list of lists
            quotient_matrix -- (default False) if True, and the partition is
        equitable, returns a matrix over the integers whose rows and columns
        represent cells of the partition, and whose i,j entry is the number of
        vertices in cell j adjacent to each vertex in cell i (since the
        partition is equitable, this is well defined)

        EXAMPLE:
            sage: G = graphs.PetersenGraph()
            sage: G.is_equitable([[0,4],[1,3,5,9],[2,6,8],[7]])
            False
            sage: G.is_equitable([[0,4],[1,3,5,9],[2,6,8,7]])
            True
            sage: G.is_equitable([[0,4],[1,3,5,9],[2,6,8,7]], quotient_matrix=True)
            [1 2 0]
            [1 0 2]
            [0 2 1]

            sage: ss = (graphs.WheelGraph(6)).line_graph(labels=False)
            sage: prt = [[(0, 1)], [(0, 2), (0, 3), (0, 4), (1, 2), (1, 4)], [(2, 3), (3, 4)]]

            sage: ss.is_equitable(prt)
            Traceback (most recent call last):
            ...
            TypeError: Partition ([[(0, 1)], [(0, 2), (0, 3), (0, 4), (1, 2), (1, 4)], [(2, 3), (3, 4)]]) is not valid for this graph: vertices are incorrect.

            sage: ss = (graphs.WheelGraph(5)).line_graph(labels=False)
            sage: ss.is_equitable(prt)
            False

        """
        from sage.misc.flatten import flatten
        from sage.misc.misc import uniq
        if sorted(flatten(partition, max_level=1)) != self.vertices():
            raise TypeError("Partition (%s) is not valid for this graph: vertices are incorrect."%partition)
        if any(len(cell)==0 for cell in partition):
            raise TypeError("Partition (%s) is not valid for this graph: there is a cell of length 0."%partition)
        if quotient_matrix:
            from sage.matrix.constructor import Matrix
            from sage.rings.integer_ring import IntegerRing
            n = len(partition)
            M = Matrix(IntegerRing(), n)
            for i in xrange(n):
                for j in xrange(n):
                    cell_i = partition[i]
                    cell_j = partition[j]
                    degrees = [self.degree_to_cell(u, cell_j) for u in cell_i]
                    if len(uniq(degrees)) > 1:
                        return False
                    if self._directed:
                        M[i, j] = degrees[0][0]
                    else:
                        M[i, j] = degrees[0]
            return M
        else:
            for cell1 in partition:
                for cell2 in partition:
                    degrees = [self.degree_to_cell(u, cell2) for u in cell1]
                    if len(uniq(degrees)) > 1:
                        return False
            return True

    def coarsest_equitable_refinement(self, partition, sparse=False):
        """
        Returns the coarsest partition which is finer than the input partition,
        and equitable with respect to self.

        A partition is equitable with respect to a graph if for every pair of
        cells C1, C2 of the partition, the number of edges from a vertex of C1
        to C2 is the same, over all vertices in C1.

        A partition P1 is finer than P2 (P2 is coarser than P1) if every cell of
        P1 is a subset of a cell of P2.

        INPUT:
            partition -- a list of lists
            sparse -- (default False) whether to use sparse or dense
        representation- for small graphs, use dense for speed

        EXAMPLES:
            sage: G = graphs.PetersenGraph()
            sage: G.coarsest_equitable_refinement([[0],range(1,10)])
            [[0], [2, 3, 6, 7, 8, 9], [1, 4, 5]]
            sage: G = graphs.CubeGraph(3)
            sage: verts = G.vertices()
            sage: Pi = [verts[:1], verts[1:]]
            sage: Pi
            [['000'], ['001', '010', '011', '100', '101', '110', '111']]
            sage: G.coarsest_equitable_refinement(Pi)
            [['000'], ['011', '101', '110'], ['111'], ['001', '010', '100']]

        Note that given an equitable partition, this function returns that
        partition:
            sage: P = graphs.PetersenGraph()
            sage: prt = [[0], [1, 4, 5], [2, 3, 6, 7, 8, 9]]
            sage: P.coarsest_equitable_refinement(prt)
            [[0], [1, 4, 5], [2, 3, 6, 7, 8, 9]]

            sage: ss = (graphs.WheelGraph(6)).line_graph(labels=False)
            sage: prt = [[(0, 1)], [(0, 2), (0, 3), (0, 4), (1, 2), (1, 4)], [(2, 3), (3, 4)]]
            sage: ss.coarsest_equitable_refinement(prt)
            Traceback (most recent call last):
            ...
            TypeError: Partition ([[(0, 1)], [(0, 2), (0, 3), (0, 4), (1, 2), (1, 4)], [(2, 3), (3, 4)]]) is not valid for this graph: vertices are incorrect.

            sage: ss = (graphs.WheelGraph(5)).line_graph(labels=False)
            sage: ss.coarsest_equitable_refinement(prt)
            [[(0, 1)], [(1, 2), (1, 4)], [(0, 3)], [(0, 2), (0, 4)], [(2, 3), (3, 4)]]

        ALGORITHM:
            Brendan D. McKay's Master's Thesis, University of Melbourne, 1976.

        """
        from sage.misc.flatten import flatten
        if sorted(flatten(partition, max_level=1)) != self.vertices():
            raise TypeError("Partition (%s) is not valid for this graph: vertices are incorrect."%partition)
        if any(len(cell)==0 for cell in partition):
            raise TypeError("Partition (%s) is not valid for this graph: there is a cell of length 0."%partition)
        if self.multiple_edges():
            raise TypeError("Refinement function does not support multiple edges.")
        G = self.copy()
        perm_to = G.relabel(return_map=True)
        partition = [[perm_to[b] for b in cell] for cell in partition]
        perm_from = {}
        for v in self:
            perm_from[perm_to[v]] = v
        n = G.num_verts()
        if sparse:
            from sage.graphs.base.sparse_graph import SparseGraph
            CG = SparseGraph(n)
        else:
            from sage.graphs.base.dense_graph import DenseGraph
            CG = DenseGraph(n)
        for i in range(n):
            for j in range(n):
                if G.has_edge(i,j):
                    CG.add_arc(i,j)
        from sage.graphs.graph_isom import PartitionStack
        nu = PartitionStack(partition)
        alpha = []
        for i in xrange(len(partition)):
            alpha.append(sum([len(cell) for cell in partition[:i]]))
        nu.refine(CG, alpha, n, (self._directed or self.loops()), True)
        repr = nu.repr_at_k(1)
        result = [[int(b) for b in a.split(',')] for a in repr[1:-1].split('|')]
        return [[perm_from[b] for b in cell] for cell in result]

    def automorphism_group(self, partition=None, translation=False,
                           verbosity=0, edge_labels=False, order=False,
                           return_group=True, orbits=False):
        """
        Returns the largest subgroup of the automorphism group of the (di)graph
        whose orbit partition is finer than the partition given. If no
        partition is given, the unit partition is used and the entire
        automorphism group is given.

        INPUT:
        translation -- if True, then output includes a a dictionary
            translating from keys == vertices to entries == elements of
            {1,2,...,n} (since permutation groups can currently only act on
            positive integers).
        partition -- default is the unit partition, otherwise computes the
            subgroup of the full automorphism group respecting the partition.
        edge_labels -- default False, otherwise allows only permutations
            respecting edge labels.
        order -- (default False) if True, compute the order of the automorphism
            group
        return_group -- default True
        orbits -- returns the orbits of the group acting on the vertices of the
            graph

        OUTPUT:
            The order of the output is group, translation, order, orbits.
        However, there are options to turn each of these on or off.

        EXAMPLES:
        Graphs:
            sage: graphs_query = GraphDatabase()
            sage: L = graphs_query.get_list(num_vertices=4)
            sage: graphs_list.show_graphs(L)
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

        Multigraphs:
            sage: G = Graph(multiedges=True, implementation='networkx')
            sage: G.add_edge(('a', 'b'))
            sage: G.add_edge(('a', 'b'))
            sage: G.add_edge(('a', 'b'))
            sage: G.automorphism_group()
            Permutation Group with generators [(1,2)]

        Digraphs:
            sage: D = DiGraph( { 0:[1], 1:[2], 2:[3], 3:[4], 4:[0] } )
            sage: D.automorphism_group()
            Permutation Group with generators [(1,2,3,4,5)]

        Edge labeled graphs:
            sage: G = Graph(implementation='networkx')
            sage: G.add_edges( [(0,1,'a'),(1,2,'b'),(2,3,'c'),(3,4,'b'),(4,0,'a')] )
            sage: G.automorphism_group(edge_labels=True)
            Permutation Group with generators [(1,4)(2,3)]

        You can also ask for just the order of the group:
            sage: G = graphs.PetersenGraph()
            sage: G.automorphism_group(return_group=False, order=True)
            120

        Or, just the orbits (recall the Petersen graph is transitive!)
            sage: G = graphs.PetersenGraph()
            sage: G.automorphism_group(return_group=False, orbits=True)
            [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]

        """
        from sage.graphs.graph_isom import search_tree, perm_group_elt
        from sage.groups.perm_gps.permgroup import PermutationGroup
        dig = (self._directed or self.loops())
        if partition is None:
            partition = [self.vertices()]
        if edge_labels:
            G, partition = graph_isom_equivalent_non_edge_labeled_graph(self, partition)
            A = search_tree(G, partition, lab=False, dict_rep=True, dig=dig, verbosity=verbosity, order=order)
            if order:
                a,b,c = A
            else:
                a,b = A
            # b is a translation of the labelings
            acting_vertices = {}
            translation_d = {}
            m = G.order()
            for v in self:
                translation_d[v] = b[('o',v)]
                if b[('o',v)] == m:
                    acting_vertices[v] = 0
                else:
                    acting_vertices[v] = b[('o',v)]
            real_aut_gp = []
            n = self.order()
            for gen in a:
                gen_restr = [0]*n
                for v in self.vertex_iterator():
                    gen_restr[acting_vertices[v]] = gen[acting_vertices[v]]
                if gen_restr not in real_aut_gp:
                    real_aut_gp.append(gen_restr)
            id = range(n)
            if id in real_aut_gp:
                real_aut_gp.remove(id)
            a = real_aut_gp
            b = translation_d
        elif self.multiple_edges():
            G, partition = graph_isom_equivalent_non_multi_graph(self, partition)
            A = search_tree(G, partition, lab=False, dict_rep=True, dig=dig, verbosity=verbosity, order=order)
            if order:
                a,b,c = A
            else:
                a,b = A
            # b is a translation of the labelings
            acting_vertices = {}
            translation_d = {}
            m = G.order()
            for v in self:
                translation_d[v] = b[('o',v)]
                if b[('o',v)] == m:
                    acting_vertices[v] = 0
                else:
                    acting_vertices[v] = b[('o',v)]
            real_aut_gp = []
            n = self.order()
            for gen in a:
                gen_restr = [0]*n
                for v in self.vertex_iterator():
                    gen_restr[acting_vertices[v]] = gen[acting_vertices[v]]
                if gen_restr not in real_aut_gp:
                    real_aut_gp.append(gen_restr)
            id = range(n)
            if id in real_aut_gp:
                real_aut_gp.remove(id)
            a = real_aut_gp
            b = translation_d
        else:
            if translation:
                A = search_tree(self, partition, dict_rep=True, lab=False, dig=dig, verbosity=verbosity, order=order)
                if order:
                    a,b,c = A
                else:
                    a,b = A
            else:
                a = search_tree(self, partition, dict_rep=False, lab=False, dig=dig, verbosity=verbosity, order=order)
                if order:
                    a,c = a
        output = []
        if return_group:
            if len(a) != 0:
                output.append(PermutationGroup([perm_group_elt(aa) for aa in a]))
            else:
                output.append(PermutationGroup([[]]))
        if translation:
            output.append(b)
        if order:
            output.append(c)
        if orbits:
            from sage.graphs.graph_isom import OrbitPartition
            O = OrbitPartition(self.num_verts())
            for gen in a:
                O.incorporate_permutation(gen)
            orbit_dict = {}
            for i in range(self.num_verts()):
                j = O.find(i)
                if orbit_dict.has_key(j):
                    orbit_dict[j].append(i)
                else:
                    orbit_dict[j] = [i]
            output.append(orbit_dict.values())
        # A Python switch statement!
        return { 0: None,
                 1: output[0],
                 2: tuple(output),
                 3: tuple(output),
                 4: tuple(output)
               }[len(output)]

    def is_vertex_transitive(self, partition=None, verbosity=0,
                           edge_labels=False, order=False,
                           return_group=True, orbits=False):
        """
        Returns whether the automorphism group of self is transitive within
        the partition provided, by default the unit partition of the vertices
        of self (thus by default tests for vertex transitivity in the usual
        sense).

        EXAMPLE:
            sage: G = Graph({0:[1],1:[2]})
            sage: G.is_vertex_transitive()
            False
            sage: P = graphs.PetersenGraph()
            sage: P.is_vertex_transitive()
            True
            sage: D = graphs.DodecahedralGraph()
            sage: D.is_vertex_transitive()
            True
            sage: R = graphs.RandomGNP(2000, .01)
            sage: R.is_vertex_transitive()
            False

        """
        if partition is None:
            partition = [self.vertices()]
        new_partition = self.automorphism_group(partition,
                          verbosity=verbosity, edge_labels=edge_labels,
                          order=False, return_group=False, orbits=True)
        for cell in partition:
            for new_cell in new_partition:
                if cell[0] in new_cell:
                    if any([c not in new_cell for c in cell[1:]]):
                        return False
        return True

    def is_isomorphic(self, other, certify=False, verbosity=0, edge_labels=False):
        """
        Tests for isomorphism between self and other.

        INPUT:
            certify -- if True, then output is (a,b), where a is a boolean and b
                is either a map or None.
            edge_labels -- default False, otherwise allows only permutations
                respecting edge labels.

        EXAMPLES:
        Graphs:
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
            sage: from sage.plot.plot import GraphicsArray
            sage: from sage.graphs.graph_fast import spring_layout_fast
            sage: position_D = spring_layout_fast(D)
            sage: position_E = {}
            sage: for vert in position_D:
            ...    position_E[b[vert]] = position_D[vert]
            sage: GraphicsArray([D.plot(pos=position_D), E.plot(pos=position_E)]).show()

        Multigraphs:
            sage: G = Graph(multiedges=True)
            sage: G.add_edge((0,1,1))
            sage: G.add_edge((0,1,2))
            sage: G.add_edge((0,1,3))
            sage: G.add_edge((0,1,4))
            sage: H = Graph(multiedges=True, implementation='networkx')
            sage: H.add_edge((3,4))
            sage: H.add_edge((3,4))
            sage: H.add_edge((3,4))
            sage: H.add_edge((3,4))
            sage: G.is_isomorphic(H)
            True

        Digraphs:
            sage: A = DiGraph( { 0 : [1,2] } )
            sage: B = DiGraph( { 1 : [0,2] } )
            sage: A.is_isomorphic(B, certify=True)
            (True, {0: 1, 1: 0, 2: 2})

        Edge labeled graphs:
            sage: G = Graph(implementation='networkx')
            sage: G.add_edges( [(0,1,'a'),(1,2,'b'),(2,3,'c'),(3,4,'b'),(4,0,'a')] )
            sage: H = G.relabel([1,2,3,4,0], inplace=False)
            sage: G.is_isomorphic(H, edge_labels=True)
            True

        """
        if certify:
            if self._directed != other._directed:
                return False, None
            if self.order() != other.order():
                return False, None
            if self.size() != other.size():
                return False, None
            if self._directed:
                if sorted(list(self.in_degree_iterator())) != sorted(list(other.in_degree_iterator())):
                    return False, None
                if sorted(list(self.out_degree_iterator())) != sorted(list(other.out_degree_iterator())):
                    return False, None
            else:
                if sorted(list(self.degree_iterator())) != sorted(list(other.degree_iterator())):
                    return False, None
            b,a = self.canonical_label(certify=True, verbosity=verbosity, edge_labels=edge_labels)
            d,c = other.canonical_label(certify=True, verbosity=verbosity, edge_labels=edge_labels)
            if b == d:
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
            if self._directed != other._directed:
                return False
            if self.order() != other.order():
                return False
            if self.size() != other.size():
                return False
            if self._directed:
                if sorted(list(self.in_degree_iterator())) != sorted(list(other.in_degree_iterator())):
                    return False
                if sorted(list(self.out_degree_iterator())) != sorted(list(other.out_degree_iterator())):
                    return False
            else:
                if sorted(list(self.degree_iterator())) != sorted(list(other.degree_iterator())):
                    return False
            b = self.canonical_label(verbosity=verbosity, edge_labels=edge_labels)
            d = other.canonical_label(verbosity=verbosity, edge_labels=edge_labels)
            return b == d

    def canonical_label(self, partition=None, certify=False, verbosity=0, edge_labels=False):
        """
        Returns the canonical label with respect to the partition. If no
        partition is given, uses the unit partition.

        INPUT:
            partition -- if given, the canonical label with respect to this
                partition will be computed. The default is the unit partition.
            certify -- if True, a dictionary mapping from the (di)graph to its
                canonical label will be given.
            verbosity -- gets passed to nice: prints helpful output.
            edge_labels -- default False, otherwise allows only permutations
                respecting edge labels.

        EXAMPLE:
            sage: D = graphs.DodecahedralGraph()
            sage: E = D.canonical_label(); E
            Dodecahedron: Graph on 20 vertices
            sage: D.canonical_label(certify=True)
            (Dodecahedron: Graph on 20 vertices, {0: 0, 1: 19, 2: 16, 3: 15, 4: 9, 5: 1, 6: 10, 7: 8, 8: 14, 9: 12, 10: 17, 11: 11, 12: 5, 13: 6, 14: 2, 15: 4, 16: 3, 17: 7, 18: 13, 19: 18})
            sage: D.is_isomorphic(E)
            True

        Multigraphs:
            sage: G = Graph(multiedges=True)
            sage: G.add_edge((0,1))
            sage: G.add_edge((0,1))
            sage: G.add_edge((0,1))
            sage: G.canonical_label()
            Multi-graph on 2 vertices

        Digraphs:
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

        Edge labeled graphs:
            sage: G = Graph(implementation='networkx')
            sage: G.add_edges( [(0,1,'a'),(1,2,'b'),(2,3,'c'),(3,4,'b'),(4,0,'a')] )
            sage: G.canonical_label(edge_labels=True)
            Graph on 5 vertices

        """
        from sage.graphs.graph_isom import search_tree
        dig = (self.loops() or self._directed)
        if partition is None:
            partition = [self.vertices()]
        if edge_labels:
            G, partition = graph_isom_equivalent_non_edge_labeled_graph(self, partition)
            a,b,c = search_tree(G, partition, certify=True, dig=dig, verbosity=verbosity)
            # c is a permutation to the canonical label of G, which depends only on isomorphism class of self.
            H = self.copy()
            relabeling = {}
            for v in H:
                relabeling[v] = c[('o',v)]
            H.relabel(relabeling)
            if certify:
                return H, relabeling
            else:
                return H
        if self.multiple_edges():
            G, partition = graph_isom_equivalent_non_multi_graph(self, partition)
            a,b,c = search_tree(G, partition, certify=True, dig=dig, verbosity=verbosity)
            # c is a permutation to the canonical label of G, which depends only on isomorphism class of self.
            H = self.copy()
            relabeling = {}
            for v in H:
                relabeling[v] = c[('o',v)]
            H.relabel(relabeling)
            if certify:
                return H, relabeling
            else:
                return H
        else:
            if certify:
                a,b,c = search_tree(self, partition, certify=True, dig=dig, verbosity=verbosity)
                return b,c
            else:
                a,b = search_tree(self, partition, dig=dig, verbosity=verbosity)
                return b


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
        weighted -- whether graph thinks of itself as weighted or not. See
            self.weighted()
        format -- if None, Graph tries to guess- can be several values,
            including:
            'graph6' -- Brendan McKay's graph6 format, in a string (if the
                string has multiple graphs, the first graph is taken)
            'sparse6' -- Brendan McKay's sparse6 format, in a string (if the
                string has multiple graphs, the first graph is taken)
            'adjacency_matrix' -- a square SAGE matrix M, with M[i,j] equal
                to the number of edges \{i,j\}
            'weighted_adjacency_matrix' -- a square SAGE matrix M, with M[i,j]
                equal to the weight of the single edge \{i,j\}. Given this
                format, weighted is ignored (assumed True).
            'incidence_matrix' -- a SAGE matrix, with one column C for each
                edge, where if C represents \{i, j\}, C[i] is -1 and C[j] is 1
            'elliptic_curve_congruence' -- data must be an iterable container
                of elliptic curves, and the graph produced has each curve as a
                vertex (it's Cremona label) and an edge E-F labelled p if and
                only if E is congruent to F mod p
        boundary -- a list of boundary vertices, if empty, graph is considered
            as a 'graph without boundary'
        implementation -- what to use as a backend for the graph. Currently, the
            options are either 'networkx' or 'c_graph'
        sparse -- only for implementation == 'c_graph'. Whether to use sparse or
            dense graphs as backend. Note that currently dense graphs do not have
            edge labels, nor can they be multigraphs
        vertex_labels -- only for implementation == 'c_graph'. Whether to allow
            any object as a vertex (slower), or only the integers 0, ..., n-1,
            where n is the number of vertices.

    EXAMPLES:
    We illustrate the first six input formats (the other two
    involve packages that are currently not standard in SAGE):

    1. A NetworkX XGraph:
        sage: import networkx
        sage: g = networkx.XGraph({0:[1,2,3], 2:[4]})
        sage: Graph(g)
        Graph on 5 vertices

    In this single case, we do not make a copy of g, but just wrap the actual
    NetworkX passed- the Sage graph becomes a wrapper for the NetworkX graph.

        sage: import networkx
        sage: g = networkx.XGraph({0:[1,2,3], 2:[4]})
        sage: G = Graph(g, implementation='networkx')
        sage: H = Graph(g, implementation='networkx')
        sage: G._backend._nxg is H._backend._nxg
        True

    2. A NetworkX graph:
        sage: import networkx
        sage: g = networkx.Graph({0:[1,2,3], 2:[4]})
        sage: DiGraph(g)
        Digraph on 5 vertices

    Note that in this case, we copy the networkX structure.

        sage: import networkx
        sage: g = networkx.Graph({0:[1,2,3], 2:[4]})
        sage: G = Graph(g, implementation='networkx')
        sage: H = Graph(g, implementation='networkx')
        sage: G._backend._nxg is H._backend._nxg
        False

    3. A dictionary of dictionaries:
        sage: g = Graph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}, implementation='networkx'); g
        Graph on 5 vertices

    The labels ('x', 'z', 'a', 'out') are labels for edges. For example, 'out' is
    the label for the edge on 2 and 5. Labels can be used as weights, if all the
    labels share some common parent.

    4. A dictionary of lists:
        sage: g = Graph({0:[1,2,3], 2:[4]}); g
        Graph on 5 vertices

    5. A list of vertices and a function describing adjacencies.  Note
       that the list of vertices and the function must be enclosed in
       a list (i.e., [list of vertices, function]).

       Construct the Paley graph over GF(13).

        sage: g=Graph([GF(13), lambda i,j: i!=j and (i-j).is_square()], implementation='networkx')
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
                 lambda i,j: len(set(i).intersection(set(j)))>0], implementation='networkx')
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
        sage: Graph(A, implementation='networkx')
        Graph on 3 vertices

    7. A graph6 or sparse6 string:
    Sage automatically recognizes whether a string is in graph6 or sparse6 format:

        sage: s = ':I`AKGsaOs`cI]Gb~'
        sage: Graph(s, sparse=True)
        Looped multi-graph on 10 vertices

        sage: G = Graph('G?????')
        sage: G = Graph("G'?G?C")
        Traceback (most recent call last):
        ...
        RuntimeError: The string seems corrupt: valid characters are
        ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
        sage: G = Graph('G??????')
        Traceback (most recent call last):
        ...
        RuntimeError: The string (G??????) seems corrupt: for n = 8, the string is too long.

        sage: G = Graph(":I'AKGsaOs`cI]Gb~")
        Traceback (most recent call last):
        ...
        RuntimeError: The string seems corrupt: valid characters are
        ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

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

        sage: Graph(matrix([[1,2],[3,4]]),loops=True)
        Looped graph on 2 vertices

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
    _directed = False

    def __init__(self, data=None, pos=None, loops=False, format=None,
                 boundary=[], weighted=False, implementation='networkx',
                 sparse=True, vertex_labels=True, **kwds):
        """
        TEST:
            sage: G = Graph()
            sage: loads(dumps(G)) == G
            True

        """
        if implementation == 'networkx':
            import networkx
            from sage.graphs.base.graph_backends import NetworkXGraphBackend
        elif implementation == 'c_graph':
            from sage.graphs.base.sparse_graph import SparseGraphBackend
            from sage.graphs.base.dense_graph import DenseGraphBackend
            if kwds.get('multiedges', False) or weighted:
                sparse = True
            CGB = SparseGraphBackend if sparse else DenseGraphBackend
        else:
            raise NotImplementedError()
        from sage.structure.element import is_Matrix
        self._weighted = weighted
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
                if implementation == 'networkx':
                    self._backend = NetworkXGraphBackend(data.networkx_graph())
                elif implementation == 'c_graph':
                    if not sparse and data.multiple_edges():
                        raise TypeError("Dense graphs do not support multiple edges.")
                    self._backend = CGB(data.num_verts())
                    if data.vertices() != range(data.num_verts()):
                        verts = data.vertices()
                        self._backend.relabel(dict([[i, verts[i]] for i in xrange(data.num_verts())]), False)
                    for u,v,l in data.edge_iterator():
                        self._backend.add_edge(u,v,l,False)
                self.loops(data.loops())
                self.multiple_edges(data.multiple_edges())
                self.name(data.name())
            elif isinstance(data,list) and len(data)>=2 and callable(data[1]):
                if implementation == 'networkx':
                    # Pass XGraph a dict of lists describing the adjacencies
                    self._backend = NetworkXGraphBackend(networkx.XGraph(dict([[i]+[[j for j in data[0] if data[1](i,j)]] for i in data[0]]), selfloops=loops, **kwds))
                elif implementation == 'c_graph':
                    self._backend = CGB(len(data[0]))
                    if sorted(data[0]) != range(len(data[0])):
                        verts = data[0]
                        self._backend.relabel(dict([[i, verts[i]] for i in xrange(len(data[0]))]), False)
                    for u in xrange(len(data[0])):
                        for v in xrange(u+1):
                            uu,vv = data[0][u], data[0][v]
                            if data[1](uu,vv):
                                self.add_edge(uu,vv)
                    self.loops(loops)
            elif isinstance(data, dict):
                if implementation == 'networkx':
                    self._backend = NetworkXGraphBackend(networkx.XGraph(data, selfloops=loops, **kwds))
                elif implementation == 'c_graph':
                    keys = data.keys()
                    for u in data:
                        if isinstance(data[u], (list,dict)):
                            for v in data[u]:
                                if v not in keys:
                                    keys.append(v)
                        else:
                            raise TypeError("Dict must be in NetworkX format. (was: %s)"%data)
                    self._backend = CGB(len(keys))
                    self.loops(loops)
                    self.multiple_edges(kwds.get('multiedges', False))
                    if sorted(keys) != range(len(keys)):
                        verts = keys
                        self._backend.relabel(dict([[i, verts[i]] for i in xrange(len(keys))]), False)
                    for u in data:
                        if isinstance(data[u], dict):
                            for v in data[u].keys():
                                if data[u][v] is None and not self.has_edge(u,v):
                                    self.add_edge(u, v)
                                elif data[u][v] is not None and not self.has_edge(u,v,data[u][v]):
                                    self.add_edge(u, v, data[u][v])
                        elif isinstance(data[u], list):
                            for v in data[u]:
                                if not self.has_edge(u,v):
                                    self.add_edge(u, v)
            elif hasattr(data,'adj'):
                import networkx
                if isinstance(data, (networkx.Graph, networkx.XGraph)):
                    if implementation == 'networkx':
                        if isinstance(data, networkx.XGraph):
                            self._backend = NetworkXGraphBackend(data)
                        else:
                            self._backend = NetworkXGraphBackend(networkx.XGraph(data, selfloops=loops, **kwds))
                    elif implementation == 'c_graph':
                        keys = data.adj.keys()
                        for u in data:
                            for v in data[u]:
                                if v not in keys:
                                    keys.append(v)
                        if not sparse and isinstance(data, networkx.XGraph) and data.multiedges:
                            raise TypeError("Dense graphs do not support multiple edges.")
                        self._backend = CGB(data.order())
                        self.loops(loops)
                        self.multiple_edges(kwds.get('multiedges', False))
                        if sorted(keys) != range(data.order()):
                            self._backend.relabel(dict([[i, keys[i]] for i in xrange(len(keys))]), False)
                        if isinstance(data, networkx.XGraph):
                            for u,v,l in data.edges_iter():
                                self.add_edge(u,v,l)
                        else:
                            for u,v in data.edges_iter():
                                self.add_edge(u,v)
                else:
                    raise NotImplementedError()
            else:
                if implementation == 'networkx':
                    if isinstance(data, (int, Integer)):
                        self._backend = NetworkXGraphBackend(networkx.XGraph(selfloops=loops, **kwds))
                        self.add_vertices(xrange(data))
                    else:
                        self._backend = NetworkXGraphBackend(networkx.XGraph(data, selfloops=loops, **kwds))
                elif implementation == 'c_graph':
                    if data is None:
                        self._backend = CGB(0)
                        self.loops(loops)
                        self.multiple_edges(kwds.get('multiedges',False))
                        self.name(kwds.get('name',''))
                    elif isinstance(data, (int, Integer)):
                        self._backend = CGB(data)
                        self.loops(loops)
                        self.multiple_edges(kwds.get('multiedges',False))
                        self.name(kwds.get('name',''))
                    else:
                        raise NotImplementedError("Datatype not suppported for C graphs.")

        if format == 'graph6':
            if not isinstance(data, str):
                raise ValueError, 'If input format is graph6, then data must be a string.'
            n = data.find('\n')
            if n == -1:
                n = len(data)
            ss = data[:n]
            n, s = graph_fast.N_inverse(ss)
            m = graph_fast.R_inverse(s, n)
            expected = n*(n-1)/2 + (6 - n*(n-1)/2)%6
            if len(m) > expected:
                raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too long."%(ss,n))
            elif len(m) < expected:
                raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too short."%(ss,n))
            if implementation == 'networkx':
                self._backend = NetworkXGraphBackend(networkx.XGraph())
                self.add_vertices(xrange(n))
            else:
                self._backend = CGB(n)
            k = 0
            for i in xrange(n):
                for j in xrange(i):
                    if m[k] == '1':
                        self.add_edge(i, j)
                    k += 1
            self.loops(False)
            self.multiple_edges(False)
        elif format == 'sparse6':
            from math import ceil, floor
            from sage.misc.functional import log
            n = data.find('\n')
            if n == -1:
                n = len(data)
            s = data[:n]
            n, s = graph_fast.N_inverse(s[1:])
            if n == 0:
                edges = []
            else:
                k = int(ceil(log(n,2)))
                ords = [ord(i) for i in s]
                if any(o > 126 or o < 63 for o in ords):
                    raise RuntimeError("The string seems corrupt: valid characters are \n" + ''.join([chr(i) for i in xrange(63,127)]))
                bits = ''.join([graph_fast.binary(o-63).zfill(6) for o in ords])
                b = []
                x = []
                for i in xrange(int(floor(len(bits)/(k+1)))):
                    b.append(int(bits[(k+1)*i:(k+1)*i+1],2))
                    x.append(int(bits[(k+1)*i+1:(k+1)*i+k+1],2))
                v = 0
                edges = []
                for i in xrange(len(b)):
                    if b[i] == 1:
                        v += 1
                    if x[i] > v:
                        v = x[i]
                    else:
                        if v < n:
                            edges.append((x[i],v))
            if implementation == 'networkx':
                self._backend = NetworkXGraphBackend(networkx.XGraph(selfloops = True, multiedges = True))
                self.add_vertices(xrange(n))
            else:
                self._backend = CGB(n)
                self.loops(True)
                self.multiple_edges(True)
            for i,j in edges:
                self.add_edge(i,j)
        elif format == 'adjacency_matrix':
            if implementation == 'networkx':
                self._backend = NetworkXGraphBackend(networkx.XGraph(selfloops = loops, **kwds))
                self.add_vertices(xrange(data.nrows()))
            else:
                self._backend = CGB(data.nrows())
                self.loops(loops)
            e = []
            for i,j in data.nonzero_positions():
                if i < j and kwds.get('multiedges',False):
                    e += [(i,j)]*int(data[i][j])
                elif i < j:
                    e.append((i,j))
                elif i == j and loops and kwds.get('multiedges',False):
                    e += [(i,j)]*int(data[i][j])
                elif i == j and loops:
                    e.append((i,j))
            self.add_edges(e)
        elif format == 'weighted_adjacency_matrix':
            self._weighted = True
            if implementation == 'networkx':
                self._backend = NetworkXGraphBackend(networkx.XGraph(selfloops = loops, **kwds))
                self.add_vertices(xrange(data.nrows()))
            else:
                self._backend = CGB(data.nrows())
                self.loops(loops)
            e = []
            for i,j in data.nonzero_positions():
                if i < j:
                    e.append((i,j,data[i,j]))
                elif i == j and loops:
                    e.append((i,j,data[i,j]))
            self.add_edges(e)
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
            if implementation == 'networkx':
                self._backend = NetworkXGraphBackend(networkx.XGraph(selfloops = loops, **kwds))
                self.add_vertices(xrange(data.nrows()))
            else:
                self._backend = CGB(data.nrows())
            e = []
            for c in data.columns():
                k = c.dict().keys()
                e.append((k[0],k[1]))
            self.add_edges(e)
        elif format == 'elliptic_curve_congruence':
            from sage.rings.arith import lcm, prime_divisors, prange
            from sage.misc.misc import prod
            if implementation == 'networkx':
                self._backend = NetworkXGraphBackend(networkx.XGraph(None, selfloops=loops, **kwds))
            else:
                raise NotImplementedError()
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
        self._pos = pos
        self._boundary = boundary
        self.name( kwds.get('name', None) )

    ### Formats

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
            return graph_fast.N(n) + graph_fast.R(self._bit_vector())

    def sparse6_string(self):
        """
        Returns the sparse6 representation of the graph as an ASCII string. Only valid
        for undirected graphs on 0 to 262143 vertices, but loops and multiple edges are
        permitted.

        EXAMPLE:
            sage: G = graphs.BullGraph()
            sage: G.sparse6_string()
            ':Da@en'

            sage: G = Graph()
            sage: G.sparse6_string()
            ':?'

            sage: G = Graph(loops=True, multiedges=True)
            sage: Graph(':?') == G
            True

        """
        n = self.order()
        if n == 0:
            return ':?'
        if n > 262143:
            raise ValueError, 'sparse6 format supports graphs on 0 to 262143 vertices only.'
        else:
            vertices = self.vertices()
            n = len(vertices)
            edges = self.edges(labels=False)
            for i in range(len(edges)): # replace edge labels with natural numbers (by index in vertices)
                edges[i] = (vertices.index(edges[i][0]),vertices.index(edges[i][1]))
            # order edges
            edges.sort(compare_edges)

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

    ### Attributes

    def is_directed(self):
        """
        Since graph is undirected, returns False.

        EXAMPLE:
            sage: Graph().is_directed()
            False
        """
        return False

    ### Properties

    def eulerian_circuit(self, return_vertices=False, labels=True):
        """
        Return a list of edges forming an eulerian circuit if one
        exists.  Otherwise return False.

        This is implemented using Fleury's algorithm.  This could be
        extended to find eulerian paths too (check for existence and
        make sure you start on an odd-degree vertex if one exists).

        INPUT:
            return_vertices -- optionally provide a list of vertices
                for the path
            labels -- whether to return edges with labels (3-tuples)

        OUTPUT:
            either ([edges], [vertices]) or [edges] of an Eulerian circuit

        EXAMPLES:

            sage: g=graphs.CycleGraph(5);
            sage: g.eulerian_circuit()
            [(0, 1, None), (1, 2, None), (2, 3, None), (3, 4, None), (4, 0, None)]
            sage: g.eulerian_circuit(labels=False)
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
            sage: g = graphs.CompleteGraph(7)
            sage: edges, vertices = g.eulerian_circuit(return_vertices=True)
            sage: vertices
            [0, 1, 2, 0, 3, 1, 4, 0, 5, 1, 6, 2, 3, 4, 2, 5, 3, 6, 4, 5, 6, 0]
            sage: graphs.CompleteGraph(4).eulerian_circuit()
            False

        """
        if not self.is_eulerian():
            return False

        edge_list = []
        vertex_list = []
        g = self.copy()
        # Get first vertex
        v = g.vertex_iterator().next()
        vertex_list.append(v)
        while g.size()>0:
            for e in g.edges_incident(v, labels=labels):
                g.delete_edge(e)
                if g.is_connected():
                    break
                else:
                    g.add_edge(e)
            else:
                # Our only choice is a cut edge
                g.delete_edge(e)
                g.delete_vertex(v)

            # the following code is here so that we don't rely the
            # order of vertices in the edge tuple.
            v = e[1] if v==e[0] else e[0]
            edge_list.append(e)
            vertex_list.append(v)

        if return_vertices:
            return edge_list, vertex_list
        else:
            return edge_list

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
        try:
            self.bipartite_color()
            return True
        except:
            return False

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
            RuntimeError: Graph is not bipartite.

        """
        # Straight from the NetworkX source:
        color = {}
        for u in self:
            if u in color:
                continue
            queue = [u]
            color[u] = 1
            while queue:
                v = queue.pop()
                c = 1-color[v]
                for w in self.neighbors(v):
                    if w in color:
                        if color[w] == color[v]:
                            raise RuntimeError("Graph is not bipartite.")
                    else:
                        color[w] = c
                        queue.append(w)
        return color

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
            RuntimeError: Graph is not bipartite.

        """
        color = self.bipartite_color()
        left = [v for v in color if color[v] == 1]
        right = [v for v in color if color[v] == 0]
        return (left, right)

    def chromatic_polynomial(self):
        """
        Returns the chromatic polynomial of the graph G.

        EXAMPLES:
            sage: G = Graph({0:[1,2,3],1:[2]})
            sage: factor(G.chromatic_polynomial())
            (x - 2) * x * (x - 1)^2

            sage: g = graphs.trees(5).next()
            sage: g.chromatic_polynomial().factor()
            x * (x - 1)^4

            sage: seven_acre_wood = sum(graphs.trees(7), Graph())
            sage: seven_acre_wood.chromatic_polynomial()
            x^77 - 66*x^76 ... - 2515943049305400*x^60 ... - 66*x^12 + x^11

            sage: for i in range(2,7):
            ...     graphs.CompleteGraph(i).chromatic_polynomial().factor()
            (x - 1) * x
            (x - 2) * (x - 1) * x
            (x - 3) * (x - 2) * (x - 1) * x
            (x - 4) * (x - 3) * (x - 2) * (x - 1) * x
            (x - 5) * (x - 4) * (x - 3) * (x - 2) * (x - 1) * x

        """
        from sage.graphs.chrompoly import chromatic_polynomial
        return chromatic_polynomial(self)

    def chromatic_number(self):
        """
        Returns the minimal number of colors needed to color the
        vertices of the graph G.

        EXAMPLES:
            sage: G = Graph({0:[1,2,3],1:[2]})
            sage: G.chromatic_number()
            3
        """
        from sage.graphs.graph_coloring import chromatic_number
        return chromatic_number(self)

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
            sage: D.show(figsize=[2,2])
            sage: D = D.to_undirected()
            sage: D.show(figsize=[2,2])
            sage: D.centrality_betweenness()
            {0: 0.16666666666666666, 1: 0.16666666666666666, 2: 0.0, 3: 0.0}

        """
        import networkx
        return networkx.betweenness_centrality(self.networkx_graph(copy=False), normalized)

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
            sage: D.show(figsize=[2,2])
            sage: D = D.to_undirected()
            sage: D.show(figsize=[2,2])
            sage: D.centrality_degree()
            {0: 1.0, 1: 1.0, 2: 0.66666666666666663, 3: 0.66666666666666663}
            sage: D.centrality_degree(v=1)
            1.0
        """
        import networkx
        return networkx.degree_centrality(self.networkx_graph(copy=False), v)

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
            sage: D.show(figsize=[2,2])
            sage: D = D.to_undirected()
            sage: D.show(figsize=[2,2])
            sage: D.centrality_closeness()
            {0: 1.0, 1: 1.0, 2: 0.75, 3: 0.75}
            sage: D.centrality_closeness(v=1)
            1.0
        """
        import networkx
        return networkx.closeness_centrality(self.networkx_graph(copy=False), v)

    ### Constructors

    def to_directed(self, implementation='networkx'):
        """
        Returns a directed version of the graph. A single edge becomes two
        edges, one in each direction.

        EXAMPLE:
            sage: graphs.PetersenGraph().to_directed()
            Petersen graph: Digraph on 10 vertices

        """
        D = DiGraph(name=self.name(), pos=self._pos, boundary=self._boundary,
                    multiedges=self.multiple_edges(), implementation=implementation)
        D.name(self.name())
        D.add_vertices(self.vertex_iterator())
        for u,v,l in self.edge_iterator():
            D.add_edge(u,v,l)
            D.add_edge(v,u,l)
        if hasattr(self, '_embedding'):
            D._embedding = self._embedding
        D._weighted = self._weighted
        return D

    def to_undirected(self):
        """
        Since the graph is already undirected, simply returns a copy of itself.

        EXAMPLE:
            sage: graphs.PetersenGraph().to_undirected()
            Petersen graph: Graph on 10 vertices

        """
        return self.copy()

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
            sage: P.write_to_eps(tmp_dir() + 'sage.eps')
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

    def graphviz_string(self):
       r"""
       Returns a representation in the DOT language, ready to render in graphviz.

       REFERENCES:
           http://www.graphviz.org/doc/info/lang.html

       EXAMPLE:
           sage: G = Graph({0:{1:None,2:None}, 1:{0:None,2:None}, 2:{0:None,1:None,3:'foo'}, 3:{2:'foo'}}, implementation='networkx')
           sage: s = G.graphviz_string()
           sage: s
           'graph {\n"0";"1";"2";"3";\n"0"--"1";"0"--"2";"1"--"2";"2"--"3"[label="foo"];\n}'
       """
       return self._graphviz_string_helper("graph", "--") # edge_string is "--" for undirected graphs

    ### Cliques

    def cliques(self):
        """
        Returns the list of maximal cliques.  Each maximal clique is
        represented by a list of vertices.

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.

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
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques()
            [[0, 1, 2], [0, 1, 3]]
        """
        import networkx.cliques
        return networkx.cliques.find_cliques(self.networkx_graph(copy=False))

    def cliques_get_max_clique_graph(self, name=''):
        """
        Returns a graph constructed with maximal cliques as vertices,
        and edges between maximal cliques with common members in
        the original graph.

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.

        INPUT:
           name -- The name of the new graph.

        EXAMPLES:
            sage: (graphs.ChvatalGraph()).cliques_get_max_clique_graph()
            Graph on 24 vertices
            sage: ((graphs.ChvatalGraph()).cliques_get_max_clique_graph()).show(figsize=[2,2], vertex_size=20, vertex_labels=False)
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_get_max_clique_graph()
            Graph on 2 vertices
            sage: (G.cliques_get_max_clique_graph()).show(figsize=[2,2])
        """
        import networkx.cliques
        return Graph(networkx.cliques.make_max_clique_graph(self.networkx_graph(copy=False), name=name, create_using=networkx.xgraph.XGraph()), implementation='networkx')

    def cliques_get_clique_bipartite(self, **kwds):
        """
        Returns a bipartite graph constructed such that cliques are the
        right vertices and the left vertices are retained from the given graph.
        Right and left vertices are connected if the bottom vertex belongs to
        the clique represented by a top vertex.

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.

        EXAMPLES:
            sage: (graphs.ChvatalGraph()).cliques_get_clique_bipartite()
            Bipartite graph on 36 vertices
            sage: ((graphs.ChvatalGraph()).cliques_get_clique_bipartite()).show(figsize=[2,2], vertex_size=20, vertex_labels=False)
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_get_clique_bipartite()
            Bipartite graph on 6 vertices
            sage: (G.cliques_get_clique_bipartite()).show(figsize=[2,2])

        """
        import networkx.cliques
        from bipartite_graph import BipartiteGraph
        return BipartiteGraph(networkx.cliques.make_clique_bipartite(self.networkx_graph(copy=False), **kwds))

    def clique_number(self, cliques=None):
        """
        Returns the size of the largest clique of the graph (clique number).

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.

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
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.clique_number()
            3
        """
        import networkx.cliques
        return networkx.cliques.graph_clique_number(self.networkx_graph(copy=False), cliques)

    def cliques_vertex_clique_number(self, vertices=None, with_labels=False, cliques=None):
        r"""
        Returns a list of sizes of the largest maximal cliques containing
        each vertex.  (Returns a single value if only one input vertex).

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.

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
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_vertex_clique_number()
            [3, 3, 3, 3]
        """
        import networkx.cliques
        return networkx.cliques.node_clique_number(self.networkx_graph(copy=False), vertices, with_labels, cliques)

    def cliques_number_of(self, vertices=None, cliques=None, with_labels=False):
        """
        Returns a list of the number of maximal cliques containing
        each vertex.  (Returns a single value if only one input vertex).

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.

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
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_number_of()
            [2, 2, 1, 1]
        """
        import networkx.cliques
        return networkx.cliques.number_of_cliques(self.networkx_graph(copy=False), vertices, cliques, with_labels)

    def cliques_containing_vertex(self, vertices=None, cliques=None, with_labels=False):
        """
        Returns the cliques containing each vertex, represented as a list of
        lists.  (Returns a single list if only one input vertex).

        Currently only implemented for undirected graphs.  Use to_undirected
        to convert a digraph to an undirected graph.

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
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_containing_vertex()
            [[[0, 1, 2], [0, 1, 3]], [[0, 1, 2], [0, 1, 3]], [[0, 1, 2]], [[0, 1, 3]]]
        """
        import networkx.cliques
        return networkx.cliques.cliques_containing_node(self.networkx_graph(copy=False), vertices, cliques, with_labels)

    ### Miscellaneous

    def min_spanning_tree(self, weight_function=lambda e: 1,
                          algorithm='Kruskal',
                          starting_vertex=None ):
        """
        Returns the edges of a minimum spanning tree, if one exists,
        otherwise returns False.

        INPUT:
            weight_function -- A function that takes an edge and
                returns a numeric weight.  Defaults to assigning each edge
                a weight of 1.
            algorithm -- Three variants of algorithms are implemented:
                'Kruskal', 'Prim fringe', and 'Prim edge' (the last
                two are variants of Prim's algorithm).  Defaults to
                'Kruskal'.  Currently, 'Prim fringe' ignores the
                labels on the edges.
            starting_vertex -- The vertex with which to start Prim's
                algorithm.

        OUTPUT:
            the edges of a minimum spanning tree.

        EXAMPLES:
            sage: g=graphs.CompleteGraph(5)
            sage: len(g.min_spanning_tree())
            4
            sage: weight = lambda e: 1/( (e[0]+1)*(e[1]+1) )
            sage: g.min_spanning_tree(weight_function=weight)
            [(3, 4, None), (2, 4, None), (1, 4, None), (0, 4, None)]
            sage: g.min_spanning_tree(algorithm='Prim edge', starting_vertex=2, weight_function=weight)
            [(2, 4, None), (3, 4, None), (1, 3, None), (0, 4, None)]
            sage: g.min_spanning_tree(algorithm='Prim fringe', starting_vertex=2, weight_function=weight)
            [(4, 2), (3, 4), (1, 4), (0, 4)]


        """
        if self.is_connected()==False:
            return False

        if algorithm=='Kruskal':
            # Kruskal's algorithm
            edges=[]
            sorted_edges_iterator=iter(sorted(self.edges(), key=weight_function))
            union_find = dict([(v,None) for v in self.vertex_iterator()])
            while len(edges) < self.order()-1:
                # get next edge
                e=sorted_edges_iterator.next()
                components=[]
                for start_v in e[0:2]:
                    v=start_v
                    children=[]

                    # Find the component a vertex lives in.
                    while union_find[v] != None:
                        children.append(v)
                        v=union_find[v]

                    # Compress the paths as much as we can for
                    # efficiency reasons.
                    for child in children:
                        union_find[child]=v

                    components.append(v)

                if components[0]!=components[1]:
                    # put in edge
                    edges.append(e)

                    # Union the components by making one the parent of the
                    # other.
                    union_find[components[0]]=components[1]

            return edges

        elif algorithm=='Prim fringe':
            if starting_vertex is None:
                v = self.vertex_iterator().next()
            else:
                v = starting_vertex
            tree=set([v])
            edges=[]

            # initialize fringe_list with v's neighbors.  fringe_list
            # contains fringe_vertex: (vertex_in_tree, weight) for each
            # fringe vertex
            fringe_list=dict([u,(v,weight_function((v,u)))] for u in self[v])

            for i in xrange(self.order()-1):
                # Find the smallest-weight fringe vertex
                v=min(fringe_list,key=lambda x: fringe_list[x][1])
                edges.append((v,fringe_list[v][0]))
                tree.add(v)
                fringe_list.pop(v)

                # Update fringe list
                for neighbor in [u for u in self[v] if u not in tree]:
                    w=weight_function((v,neighbor))
                    if neighbor not in fringe_list or \
                           (neighbor in fringe_list and fringe_list[neighbor][1]>w):
                        fringe_list[neighbor]=(v,weight_function((v,neighbor)))
            return edges

        elif algorithm=='Prim edge':
            if starting_vertex is None:
                v = self.vertex_iterator().next()
            else:
                v = starting_vertex
            sorted_edges=sorted(self.edges(), key=weight_function)
            tree=set([v])
            edges=[]

            for i in xrange(self.order()-1):
                # Find a minimum-weight edge connecting a vertex in
                # the tree to something outside the tree.  Remove the
                # edges between tree vertices for efficiency.

                for i in xrange(len(sorted_edges)):
                    e=sorted_edges[i]
                    v0,v1=e[0],e[1]
                    if v0 in tree:
                        if v1 not in tree:
                            edges.append(e)
                            sorted_edges[i:i+1]=[]
                            tree.add(v1)
                            break
                        else:
                            sorted_edges[i:i+1]=[]
                    elif v1 in tree:
                        edges.append(e)
                        sorted_edges[i:i+1]=[]
                        tree.add(v0)
                        break
            return edges
        else:
            raise NotImplementedError, "Minimum Spanning Tree algorithm '%s' is not implemented."%algorithm


class DiGraph(GenericGraph):
    """
    Directed graph.

    INPUT:
        data -- can be any of the following:
            1. A NetworkX digraph
            2. A dictionary of dictionaries
            3. A dictionary of lists
            4. A numpy matrix or ndarray
            5. A Sage adjacency matrix or incidence matrix
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
        weighted -- whether digraph thinks of itself as weighted or not. See
            self.weighted()
        format -- if None, DiGraph tries to guess- can be several values, including:
            'adjacency_matrix' -- a square Sage matrix M, with M[i,j] equal to the number
                                  of edges \{i,j\}
            'incidence_matrix' -- a Sage matrix, with one column C for each edge, where
                                  if C represents \{i, j\}, C[i] is -1 and C[j] is 1
            'weighted_adjacency_matrix' -- a square Sage matrix M, with M[i,j]
                equal to the weight of the single edge \{i,j\}. Given this
                format, weighted is ignored (assumed True).
        boundary -- a list of boundary vertices, if none, digraph is considered as a 'digraph
                    without boundary'
        implementation -- what to use as a backend for the graph. Currently, the
            options are either 'networkx' or 'c_graph'
        sparse -- only for implementation == 'c_graph'. Whether to use sparse or
            dense graphs as backend. Note that currently dense graphs do not have
            edge labels, nor can they be multigraphs
        vertex_labels -- only for implementation == 'c_graph'. Whether to allow
            any object as a vertex (slower), or only the integers 0, ..., n-1,
            where n is the number of vertices.

    EXAMPLES:
    1. A NetworkX XDiGraph:
        sage: import networkx
        sage: g = networkx.XDiGraph({0:[1,2,3], 2:[4]})
        sage: DiGraph(g)
        Digraph on 5 vertices

    In this single case, we do not make a copy of g, but just wrap the actual
    NetworkX passed.  We do this for performance reasons.

        sage: import networkx
        sage: g = networkx.XDiGraph({0:[1,2,3], 2:[4]})
        sage: G = DiGraph(g, implementation='networkx')
        sage: H = DiGraph(g, implementation='networkx')
        sage: G._backend._nxg is H._backend._nxg
        True

    2. A NetworkX digraph:
        sage: import networkx
        sage: g = networkx.DiGraph({0:[1,2,3], 2:[4]})
        sage: DiGraph(g)
        Digraph on 5 vertices

    Note that in this case, we copy the networkX structure.

        sage: import networkx
        sage: g = networkx.DiGraph({0:[1,2,3], 2:[4]})
        sage: G = DiGraph(g, implementation='networkx')
        sage: H = DiGraph(g, implementation='networkx')
        sage: G._backend._nxg is H._backend._nxg
        False


    3. A dictionary of dictionaries:
        sage: g = DiGraph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}, implementation='networkx'); g
        Digraph on 5 vertices

    The labels ('x', 'z', 'a', 'out') are labels for edges. For example, 'out' is
    the label for the edge from 2 to 5. Labels can be used as weights, if all the
    labels share some common parent.

    4. A dictionary of lists:
        sage: g = DiGraph({0:[1,2,3], 2:[4]}); g
        Digraph on 5 vertices

    5. A list of vertices and a function describing adjacencies.  Note
       that the list of vertices and the function must be enclosed in
       a list (i.e., [list of vertices, function]).

       We construct a graph on the integers 1 through 12 such that
       there is a directed edge from i to j if and only if i divides j.

        sage: g=DiGraph([[1..12],lambda i,j: i!=j and i.divides(j)], implementation='networkx')
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
        sage: DiGraph(A, implementation='networkx')
        Digraph on 3 vertices

    7. A \sage matrix:
    Note: If format is not specified, then \sage assumes a square matrix is an adjacency
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

    8. A c_graph implemented DiGraph can be constructed from a networkx
    implemented DiGraph if its vertex set is equal to range(n):

        sage: D = DiGraph({0:[1],1:[2],2:[0]}, implementation="networkx")
        sage: E = DiGraph(D,implementation="c_graph")
        sage: D == E
        True

    9. A dig6 string:
    Sage automatically recognizes whether a string is in dig6 format, which
    is a directed version of graph6:

        sage: D = DiGraph('IRAaDCIIOWEOKcPWAo')
        sage: D
        Digraph on 10 vertices

        sage: D = DiGraph('IRAaDCIIOEOKcPWAo')
        Traceback (most recent call last):
        ...
        RuntimeError: The string (IRAaDCIIOEOKcPWAo) seems corrupt: for n = 10, the string is too short.

        sage: D = DiGraph("IRAaDCI'OWEOKcPWAo")
        Traceback (most recent call last):
        ...
        RuntimeError: The string seems corrupt: valid characters are
        ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

    """
    _directed = True

    def __init__(self, data=None, pos=None, loops=False, format=None,
                 boundary=[], weighted=False, implementation='networkx',
                 sparse=True, vertex_labels=True, **kwds):
        """
        TEST:
            sage: D = DiGraph()
            sage: loads(dumps(D)) == D
            True
        """
        if implementation == 'networkx':
            import networkx
            from sage.graphs.base.graph_backends import NetworkXGraphBackend
        elif implementation == 'c_graph':
            from sage.graphs.base.sparse_graph import SparseGraphBackend
            from sage.graphs.base.dense_graph import DenseGraphBackend
            if kwds.get('multiedges', False) or weighted:
                sparse = True
            CGB = SparseGraphBackend if sparse else DenseGraphBackend
        else:
            raise NotImplementedError()
        from sage.structure.element import is_Matrix
        self._weighted = weighted
        if format is None:
            if is_Matrix(data):
                if data.is_square():
                    format = 'adjacency_matrix'
                else:
                    format = 'incidence_matrix'
            elif isinstance(data, DiGraph):
                if implementation == 'networkx':
                    self._backend = NetworkXGraphBackend(data.networkx_graph())
                    self.loops(data.loops())
                    self.multiple_edges(data.multiple_edges())
                    self.name(data.name())
                elif implementation == 'c_graph':
                    if not sparse and data.multiple_edges():
                        raise TypeError("Dense graphs do not support multiple edges.")
                    self._backend = CGB(data.num_verts())
                    self.loops(data.loops())
                    self.multiple_edges(data.multiple_edges())
                    self.name(data.name())
                    if data.vertices() != range(data.num_verts()):
                        verts = data.vertices()
                        self._backend.relabel(dict([[i, verts[i]] for i in xrange(data.num_verts())]), False)
                    for u,v,l in data.edge_iterator():
                        self._backend.add_edge(u,v,l,True)
            elif hasattr(data, 'adj'):
                import networkx
                if isinstance(data, (networkx.XGraph, networkx.Graph)) and not isinstance(data, (networkx.XDiGraph, networkx.DiGraph)):
                    data = data.to_directed()
                if isinstance(data, (networkx.XDiGraph, networkx.DiGraph)):
                    if implementation == 'networkx':
                        if isinstance(data, networkx.XDiGraph):
                            self._backend = NetworkXGraphBackend(data)
                        else:
                            self._backend = NetworkXGraphBackend(networkx.XDiGraph(data, selfloops=loops, **kwds))
                    elif implementation == 'c_graph':
                        keys = data.adj.keys()
                        for u in data:
                            for v in data[u]:
                                if v not in keys:
                                    keys.append(v)
                        if not sparse and isinstance(data, networkx.XDiGraph) and data.multiedges:
                            raise TypeError("Dense graphs do not support multiple edges.")
                        self._backend = CGB(data.order())
                        self.loops(loops)
                        self.multiple_edges(kwds.get('multiedges', False))
                        if sorted(keys) != range(data.order()):
                            self._backend.relabel(dict([[i, keys[i]] for i in xrange(len(keys))]), False)
                        if isinstance(data, networkx.XDiGraph):
                            for u,v,l in data.edges_iter():
                                self.add_edge(u,v,l)
                        else:
                            for u,v in data.edges_iter():
                                self.add_edge(u,v)
                else:
                    raise NotImplementedError()
            elif isinstance(data, str):
                format = 'dig6'
            elif isinstance(data,list) and len(data)>=2 and callable(data[1]):
                if implementation == 'networkx':
                    # Pass XGraph a dict of lists describing the adjacencies
                    self._backend = NetworkXGraphBackend(networkx.XDiGraph(dict([[i]+[[j for j in data[0] if data[1](i,j)]] for i in data[0]]), selfloops=loops, **kwds))
                elif implementation == 'c_graph':
                    self._backend = CGB(len(data[0]))
                    if sorted(data[0]) != range(len(data[0])):
                        verts = data[0]
                        self._backend.relabel(dict([[i, verts[i]] for i in xrange(len(data[0]))]), False)
                    for u in xrange(len(data[0])):
                        for v in xrange(len(data[0])):
                            uu,vv = data[0][u], data[0][v]
                            if data[1](uu,vv):
                                self.add_edge(uu,vv)
                    self.loops(loops)
            elif isinstance(data, dict):
                if implementation == 'networkx':
                    self._backend = NetworkXGraphBackend(networkx.XDiGraph(data, selfloops=loops, **kwds))
                elif implementation == 'c_graph':
                    keys = data.keys()
                    for u in data:
                        if isinstance(data[u], (list,dict)):
                            for v in data[u]:
                                if v not in keys:
                                    keys.append(v)
                        else:
                            raise TypeError("Dict must be in NetworkX format. (was: %s)"%data)
                    self._backend = CGB(len(keys))
                    self.loops(loops)
                    self.multiple_edges(kwds.get('multiedges', False))
                    if sorted(keys) != range(len(keys)):
                        verts = keys
                        self._backend.relabel(dict([[i, verts[i]] for i in xrange(len(keys))]), False)
                    for u in data:
                        if isinstance(data[u], dict):
                            for v in data[u].keys():
                                if data[u][v] is None and not self.has_edge(u,v):
                                    self.add_edge(u, v)
                                elif data[u][v] is not None and not self.has_edge(u,v,data[u][v]):
                                    self.add_edge(u, v, data[u][v])
                        elif isinstance(data[u], list):
                            for v in data[u]:
                                if not self.has_edge(u,v):
                                    self.add_edge(u, v)
            elif isinstance(data, (int, Integer)):
                if implementation == 'networkx':
                    self._backend = NetworkXGraphBackend(networkx.XDiGraph(selfloops=loops, **kwds))
                    self.add_vertices(xrange(data))
                elif implementation == 'c_graph':
                    self._backend = CGB(data)
                    self.loops(loops)
                    self.multiple_edges(kwds.get('multiedges',False))
                    self.name(kwds.get('name',''))
            else:
                if implementation == 'networkx':
                    self._backend = NetworkXGraphBackend(networkx.XDiGraph(data, selfloops=loops, **kwds))
                elif implementation == 'c_graph':
                    if data is None:
                        self._backend = CGB(0)
                        self.loops(loops)
                        self.multiple_edges(kwds.get('multiedges',False))
                    else:
                        raise NotImplementedError("Datatype not suppported for C graphs.")
        if format == 'adjacency_matrix':
            if implementation == 'networkx':
                self._backend = NetworkXGraphBackend(networkx.XDiGraph(selfloops = loops, **kwds))
                self.add_vertices(xrange(data.nrows()))
            else:
                self._backend = CGB(data.nrows())
                self.loops(loops)
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
            self.add_edges(e)
        elif format == 'weighted_adjacency_matrix':
            self._weighted = True
            if implementation == 'networkx':
                self._backend = NetworkXGraphBackend(networkx.XDiGraph(d, selfloops = loops, **kwds))
                self.add_vertices(xrange(data.nrows()))
            else:
                self._backend = CGB(data.nrows())
                self.loops(loops)
            e = []
            for i,j in data.nonzero_positions():
                if i != j:
                    e.append((i,j,data[i][j]))
                elif i == j and loops:
                    e.append((i,j,data[i][j]))
            self.add_edges(e)
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
            if implementation == 'networkx':
                self._backend = NetworkXGraphBackend(networkx.XDiGraph(selfloops = loops, **kwds))
                self.add_vertices(xrange(data.nrows()))
            else:
                self._backend = CGB(data.nrows())
            e = []
            for c in data.columns():
                k = c.dict().keys()
                if c[k[0]] == -1:
                    e.append((k[0],k[1]))
                else:
                    e.append((k[1],k[0]))
            self.add_edges(e)
        elif format == 'dig6':
            if not isinstance(data, str):
                raise ValueError, 'If input format is dig6, then data must be a string.'
            n = data.find('\n')
            if n == -1:
                n = len(data)
            ss = data[:n]
            n, s = graph_fast.N_inverse(ss)
            m = graph_fast.D_inverse(s, n)
            expected = n**2
            if len(m) > expected:
                raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too long."%(ss,n))
            elif len(m) < expected:
                raise RuntimeError("The string (%s) seems corrupt: for n = %d, the string is too short."%(ss,n))
            if implementation == 'networkx':
                self._backend = NetworkXGraphBackend(networkx.XDiGraph())
                self.add_vertices(xrange(n))
            else:
                self._backend = CGB(n)
            k = 0
            for i in range(n):
                for j in range(n):
                    if m[k] == '1':
                        self.add_edge(i, j)
                    k += 1
            self.loops(False)
        self._pos = pos
        self._boundary = boundary
        self.name( kwds.get('name', None) )

    ### Formats

    def dig6_string(self):
        """
        Returns the dig6 representation of the digraph as an ASCII string.
        Valid for single (no multiple edges) digraphs on 0 to 262143 vertices.

        EXAMPLE:
            sage: D = DiGraph()
            sage: D.dig6_string()
            '?'
            sage: D.add_edge(0,1)
            sage: D.dig6_string()
            'AO'

        """
        n = self.order()
        if n > 262143:
            raise ValueError, 'dig6 format supports graphs on 0 to 262143 vertices only.'
        elif self.multiple_edges():
            raise ValueError, 'dig6 format does not support multiple edges.'
        else:
            return graph_fast.N(n) + graph_fast.R(self._bit_vector())

    ### Attributes

    def is_directed(self):
        """
        Since digraph is directed, returns True.

        EXAMPLE:
            sage: DiGraph().is_directed()
            True
        """
        return True

    ### Properties

    def is_directed_acyclic(self):
        """
        Returns whether the digraph is acyclic or not.

        A directed graph is acyclic if for any vertex v, there is no directed
        path that starts and ends at v. Every directed acyclic graph (dag)
        corresponds to a partial ordering of its vertices, however multiple
        dags may lead to the same partial ordering.

        EXAMPLES:
            sage: D = DiGraph({ 0:[1,2,3], 4:[2,5], 1:[8], 2:[7], 3:[7], 5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10] })
            sage: D.plot(layout='circular').show()
            sage: D.is_directed_acyclic()
            True

            sage: D.add_edge(9,7)
            sage: D.is_directed_acyclic()
            True

            sage: D.add_edge(7,4)
            sage: D.is_directed_acyclic()
            False

        """
        try:
            S = self.topological_sort()
            return True
        except:
            return False

    def to_directed(self):
        """
        Since the graph is already directed, simply returns a copy of itself.

        EXAMPLE:
            sage: DiGraph({0:[1,2,3],4:[5,1]}).to_directed()
            Digraph on 6 vertices

        """
        return self.copy()

    def to_undirected(self, implementation='networkx'):
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
        G = Graph(name=self.name(), pos=self._pos, boundary=self._boundary,
                  multiedges=self.multiple_edges(), loops=self.loops(), implementation=implementation)
        G.name(self.name())
        G.add_vertices(self.vertex_iterator())
        G.add_edges(self.edge_iterator())
        if hasattr(self, '_embedding'):
            G._embedding = self._embedding
        G._weighted = self._weighted
        return G

    ### Edge Handlers

    def incoming_edge_iterator(self, vertices=None, labels=True):
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
        if vertices is None:
            vertices = self
        elif vertices in self:
            vertices = [vertices]
        else:
            vertices = [v for v in vertices if v in self]
        return self._backend.iterator_in_edges(vertices, labels)

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
        return list(self.incoming_edge_iterator(vertices, labels=labels))

    def outgoing_edge_iterator(self, vertices=None, labels=True):
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
        if vertices is None:
            vertices = self
        elif vertices in self:
            vertices = [vertices]
        else:
            vertices = [v for v in vertices if v in self]
        return self._backend.iterator_out_edges(vertices, labels)

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
        return list(self.outgoing_edge_iterator(vertices, labels=labels))

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
        return self._backend.iterator_in_nbrs(vertex)

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
        return self._backend.iterator_out_nbrs(vertex)

    def successors(self, vertex):
        """
        Returns a list of successor vertices of vertex.

        EXAMPLE:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.successors(0)
            [1, 2, 3]

        """
        return list(self.successor_iterator(vertex))

    ### Degree functions

    def in_degree(self, vertices=None, labels=False):
        """
        Same as degree, but for in degree.

        EXAMPLES:
            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.in_degree(vertices = [0,1,2], labels=True)
            {0: 2, 1: 2, 2: 2}
            sage: D.in_degree()
            [2, 2, 2, 2, 1, 1]
            sage: G = graphs.PetersenGraph().to_directed()
            sage: G.in_degree(0)
            3

        """
        if vertices in self:
            for d in self.in_degree_iterator(vertices):
                return d # (weird, but works: only happends once!)
        elif labels:
            di = {}
            for v, d in self.in_degree_iterator(vertices, labels=labels):
                di[v] = d
            return di
        else:
            return list(self.in_degree_iterator(vertices, labels=labels))

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
        if vertices is None:
            vertices = self.vertex_iterator()
        elif vertices in self:
            d = 0
            for e in self.incoming_edge_iterator(vertices):
                d += 1
            yield d
            return
        if labels:
            for v in vertices:
                d = 0
                for e in self.incoming_edge_iterator(v):
                    d += 1
                yield (v, d)
        else:
            for v in vertices:
                d = 0
                for e in self.incoming_edge_iterator(v):
                    d += 1
                yield d

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
        if vertices in self:
            for d in self.out_degree_iterator(vertices):
                return d # (weird, but works: only happends once!)
        elif labels:
            di = {}
            for v, d in self.out_degree_iterator(vertices, labels=labels):
                di[v] = d
            return di
        else:
            return list(self.out_degree_iterator(vertices, labels=labels))

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
        if vertices is None:
            vertices = self.vertex_iterator()
        elif vertices in self:
            d = 0
            for e in self.outgoing_edge_iterator(vertices):
                d += 1
            yield d
            return
        if labels:
            for v in vertices:
                d = 0
                for e in self.outgoing_edge_iterator(v):
                    d += 1
                yield (v, d)
        else:
            for v in vertices:
                d = 0
                for e in self.outgoing_edge_iterator(v):
                    d += 1
                yield d

    ### Construction

    def reverse(self):
        """
        Returns a copy of digraph with edges reversed in direction.

        EXAMPLES:
            sage: D = DiGraph({ 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] })
            sage: D.reverse()
            Reverse of (): Digraph on 6 vertices

        """
        H = DiGraph(multiedges=self.multiple_edges(), loops=self.loops(), implementation='networkx')
        H.add_vertices(self)
        H.add_edges( [ (v,u,d) for (u,v,d) in self.edge_iterator() ] )
        name = self.name()
        if name is None:
            name = ''
        H.name("Reverse of (%s)"%name)
        return H

    ### Directed Acyclic Graphs (DAGs)

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
            sage: D.plot(layout='circular').show()
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
        S = networkx.topological_sort(self.networkx_graph(copy=False))
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
            sage: D.plot(layout='circular').show()
            sage: D.topological_sort_generator()
            [[0, 1, 2, 3, 4], [0, 1, 2, 4, 3], [0, 2, 1, 3, 4], [0, 2, 1, 4, 3], [0, 2, 4, 1, 3]]

            sage: for sort in D.topological_sort_generator():
            ...       for edge in D.edge_iterator():
            ...           u,v,l = edge
            ...           if sort.index(u) > sort.index(v):
            ...               print "This should never happen."

        """
        from sage.graphs.linearextensions import LinearExtensions
        try:
            return LinearExtensions(self).list()
        except TypeError:
            raise TypeError('Digraph is not acyclic-- there is no topological sort (or there was an error in sage/graphs/linearextensions.py).')

    ### Visualization

    def graphviz_string(self):
       r"""
       Returns a representation in the DOT language, ready to render in graphviz.

       REFERENCES:
           http://www.graphviz.org/doc/info/lang.html

       EXAMPLE:
           sage: G = DiGraph({0:{1:None,2:None}, 1:{2:None}, 2:{3:'foo'}, 3:{}}, implementation='networkx')
           sage: s = G.graphviz_string(); s
           'digraph {\n"0";"1";"2";"3";\n"0"->"1";"0"->"2";"1"->"2";"2"->"3"[label="foo"];\n}'
       """
       return self._graphviz_string_helper("digraph", "->") # edge_string is "->" for directed graphs





def tachyon_vertex_plot(g, bgcolor=(1,1,1),
                        vertex_colors=None,
                        vertex_size=0.06,
                        pos3d=None,
                        iterations=50, **kwds):
    """
    Helper function for plotting graphs in 3d with Tachyon. Returns a plot containing
    only the vertices, as well as the 3d position dictionary used for the plot.

    EXAMPLE:
        sage: G = graphs.TetrahedralGraph()
        sage: from sage.graphs.graph import tachyon_vertex_plot
        sage: T,p = tachyon_vertex_plot(G)
        sage: type(T)
        <class 'sage.plot.tachyon.Tachyon'>
        sage: type(p)
        <type 'dict'>

    """
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


def graph_isom_equivalent_non_multi_graph(g, partition):
    r"""
    Helper function for canonical labeling of multi-(di)graphs.

    The idea for this function is that the main algorithm for computing isomorphism
    of graphs does not allow multiple edges. Instead of making some very difficult
    changes to that, we can simply modify the multigraph into a non-multi graph
    that carries essentially the same information. For each pair of vertices
    $\{u,v\}$, if there is at most one edge between $u$ and $v$, we do nothing,
    but if there are more than one, we split each edge into two, introducing a new
    vertex. These vertices really represent edges, so we keep them in their own
    part of a partition -- to distinguish them from genuine vertices. Then the
    canonical label and automorphism group is computed, and in the end, we strip
    off the parts of the generators that describe how these new vertices move,
    and we have the automorphism group of the original multi-graph. Similarly, by
    putting the additional vertices in their own cell of the partition, we guarantee
    that the relabeling leading to a canonical label moves genuine vertices amongst
    themselves, and hence the canonical label is still well-defined, when we forget
    about the additional vertices.


    EXAMPLE:
        sage: from sage.graphs.graph import graph_isom_equivalent_non_multi_graph
        sage: G = Graph(multiedges=True)
        sage: G.add_edge((0,1,1))
        sage: G.add_edge((0,1,2))
        sage: G.add_edge((0,1,3))
        sage: graph_isom_equivalent_non_multi_graph(G, [[0,1]])
        (Graph on 5 vertices, [[('o', 0), ('o', 1)], [('x', 0), ('x', 1), ('x', 2)]])

    """
    if g._directed:
        G = DiGraph(loops=g.loops())
    else:
        G = Graph(loops=g.loops(), implementation='networkx')
    G.add_vertices([('o', v) for v in g.vertices()]) # 'o' for original
    if g._directed:
        edges_with_multiplicity = [[u,v] for u,v,_ in g.edge_iterator()]
    else:
        edges_with_multiplicity = [sorted([u,v]) for u,v,_ in g.edge_iterator()]
    index = 0
    while len(edges_with_multiplicity) > 0:
        [u,v] = edges_with_multiplicity.pop(0)
        m = edges_with_multiplicity.count([u,v]) + 1
        if m == 1:
            G.add_edge((('o',u),('o',v)))
        else:
            for _ in xrange(m):
                G.add_edges([(('o',u), ('x', index)), (('x', index), ('o',v))]) # 'x' for extra
                index += 1
            edges_with_multiplicity = [e for e in edges_with_multiplicity if e != [u,v]]
    new_partition = [[('o',v) for v in cell] for cell in partition] + [[('x',i) for i in xrange(index)]]
    return G, new_partition


def graph_isom_equivalent_non_edge_labeled_graph(g, partition):
    """
    Helper function for canonical labeling of edge labeled (di)graphs.

    Translates to a bipartite incidence-structure type graph appropriate for
    computing canonical labels of edge labeled graphs. Note that this is actually
    computationally equivalent to implementing a change on an inner loop of the
    main algorithm-- namely making the refinement procedure sort for each label.

    If the graph is a multigraph, it is translated to a non-multigraph, where each
    edge is labeled with a dictionary describing how many edges of each label were
    originally there. Then in either case we are working on a graph without multiple
    edges. At this point, we create another (bipartite) graph, whose left vertices
    are the original vertices of the graph, and whose right vertices represent the
    edges. We partition the left vertices as they were originally, and the right
    vertices by common labels: only automorphisms taking edges to like-labeled edges
    are allowed, and this additional partition information enforces this on the
    bipartite graph.


    EXAMPLE:
        sage: G = Graph(multiedges=True, implementation='networkx')
        sage: G.add_edges([(0,1,i) for i in range(10)])
        sage: G.add_edge(1,2,'string')
        sage: G.add_edge(2,3)
        sage: from sage.graphs.graph import graph_isom_equivalent_non_edge_labeled_graph
        sage: graph_isom_equivalent_non_edge_labeled_graph(G, [G.vertices()])
        (Graph on 7 vertices,
         [[('o', 0), ('o', 1), ('o', 2), ('o', 3)], [('x', 0)], [('x', 1)], [('x', 2)]])

    """
    if g.multiple_edges():
        if g._directed:
            G = DiGraph(loops=g.loops())
        else:
            G = Graph(loops=g.loops(), implementation='networkx')
        G.add_vertices(g.vertices())
        for u,v,l in g.edge_iterator():
            if not G.has_edge(u,v):
                G.add_edge(u,v,[[l,1]])
            else:
                label_list = G.edge_label(u,v)
                seen_label = False
                for i in range(len(label_list)):
                    if label_list[i][0] == l:
                        label_list[i][1] += 1
                        seen_label = True
                        break
                if not seen_label:
                    label_list.append([l,1])
        g = G
    edge_partition = []
    if g._directed:
        G = DiGraph(loops=g.loops())
    else:
        G = Graph(loops=g.loops(), implementation='networkx')
    G.add_vertices([('o', v) for v in g.vertices()]) # 'o' for original
    index = 0
    for u,v,l in g.edge_iterator():
        if len([a for a in edge_partition if a[0] == l]) == 0:
            edge_partition.append([l, [index]])
        else:
            i = 0
            while edge_partition[i][0] != l:
                i += 1
            edge_partition[i][1].append(index)
        G.add_edges([(('o',u), ('x', index)), (('x', index), ('o',v))]) # 'x' for extra
        index += 1
    new_partition = [[('o',v) for v in cell] for cell in partition] + [[('x',v) for v in a[1]] for a in edge_partition]
    return G, new_partition



def compare_edges(x, y):
    """
    Compare edge x to edge y, return -1 if x < y, 1 if x > y, else 0.

    EXAMPLE:
        sage: G = graphs.PetersenGraph()
        sage: E = G.edges()
        sage: from sage.graphs.graph import compare_edges
        sage: compare_edges(E[0], E[2])
        -1
        sage: compare_edges(E[0], E[1])
        -1
        sage: compare_edges(E[0], E[0])
        0
        sage: compare_edges(E[1], E[0])
        1

    """
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


