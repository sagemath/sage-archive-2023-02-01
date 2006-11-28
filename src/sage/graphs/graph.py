r"""
Abstract base classes for graphs; wrapper for NetworkX.

Many functions are passed directly to NetworkX, and in this case the documen-
tation is paraphrased from the NX docs.

AUTHOR:
    - Robert L. Miller (2006-10-22): initial version

\section{Tutorial}

EXAMPLES:
TODO
"""

#*****************************************************************************
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

import networkx as NX     # the LANL library for graph theory
from sage.structure.sage_object import SageObject
from sage.plot.plot import Graphics, GraphicPrimitive_NetworkXGraph
from sage.matrix.constructor import matrix
from sage.rings.integer_mod_ring import IntegerModRing
from random import random

class GenericGraph(SageObject):
    pass

class SimpleGraph(GenericGraph):

    def __getitem__(self,vertex):
        """
        G[vertex] returns the neighbors (in & out if digraph) of vertex.
        """
        return self.neighbors(vertex)

    ### Vertex handlers

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

class Graph(SimpleGraph):
    """
    Undirected, simple graph (no loops, no multiple edges, non-hyper).
    """

    ### Note: NetworkX function print_dna not wrapped.

    def __contains__(self, vertex):
        """
        Allows 'v in G'.
        """
        return vertex in self.__nxg

    def __init__(self, data=None, pos=None, **kwds):
        """
        Initialize graph.

        INPUT:
        data -- can be any of the following:
            1 NetworkX graph
            2 dictionary of dictionaries
            3 dictionary of lists
            4 numpy matrix or ndarray
            5 pygraphviz agraph
            6 scipy sparse matrix
        pos -- an optional positioning dictionary: for example, the
        spring layout from NetworkX for the 5-cycle is
            {   0: [-0.91679746, 0.88169588,],
                1: [ 0.47294849, 1.125     ,],
                2: [ 1.125     ,-0.12867615,],
                3: [ 0.12743933,-1.125     ,],
                4: [-1.125     ,-0.50118505,]   }
        name -- (in kwds) gives the graph a name

        EXAMPLES:
        sage: G = Graph(name="Null graph")
        sage: G
        Null graph: a simple graph on 0 vertices

        sage: P = Graph({0:[1,4,5],1:[0,2,6],2:[1,3,7],3:[2,4,8],4:[0,3,9],5:[0,7,8],6:[1,8,9],7:[2,5,9],8:[3,5,6],9:[4,6,7]},name="Petersen graph")
        sage: P
        Petersen graph: a simple graph on 10 vertices

        """
        if isinstance(data, Graph):
            self.__nxg = data.networkx_graph()
        elif isinstance(data, NX.Graph):
            self.__nxg = data
        else:
            self.__nxg = NX.Graph(data, **kwds)
        ### NOTE: Name bug in NetworkX supposedly fixed, so check if the fol-
        ### lowing fix is unnecessary
        self.__nxg.name=kwds.get("name","No Name")
        ### end fix
        self.__pos = pos

    def __iter__(self):
        """
        Return an iterator over the vertices. Allows 'for v in G'.
        """
        return self.vertex_iterator()

    def __len__(self):
        return len(self.__nxg.adj)

    def __str__(self):
        if self.__nxg.name != "No Name":
            return self.__nxg.name
        else: return repr(self)

    def _latex_(self):
        return repr(self)

    def _matrix_(self):
        return self.am()

    def _repr_(self):
        if self.__nxg.name != "No Name":
            name = self.__nxg.name
            name = name + ": a s"
        else: name = "S"
        return name + "imple graph on %d vertices"%len(self.__nxg.adj)

    def clear(self):
        """
        Empties the graph of vertices and edges, removes name.
        """
        self.__nxg.clear()

    def copy(self):
        """
        Creates a copy of the graph.
        """
        G = Graph(self.__nxg, name=self.__nxg.name)
        return G

    def networkx_graph(self):
        """
        Creates a NetworkX graph from the SAGE graph.
        """
        return self.__nxg.copy()

    def to_directed(self):
        pass # NOT YET IMPLEMENTED
        ### NX: to_directed

    def to_undirected(self):
        return self,copy()

    ### General properties

    def density(self):
        """
        Returns the density.
        """
        return NX.density(self.__nxg)

    def is_directed(self):
        return False

    def networkx_info(self, vertex=None):
        """
        Returns NetworkX information about the graph or the given node.
        """
        self.__nxg.info(vertex)

    def order(self):
        """
        Returns the number of vertices.
        """
        return self.__nxg.order()

    def size(self):
        """
        Returns the number of edges.
        """
        return self.__nxg.size()

    ### Vertex handlers

    def add_vertex(self, name=None):
        """
        Creates an isolated vertex

        INPUT:
        n -- Name of the new vertex. If no name is specified, then the vertex
        will be represented by the least integer not already representing a
        vertex. Name must be an immutable object.
        """
        ### TODO- add doc note about representing other objects as vertices
        ### This will be done when such representation is implemented
        if name is None: # then find an integer to use as a key
            i = 0
            while self.__nxg.adj.has_key(i):
                    i=i+1
            self.__nxg.add_node(i)
        else:
            self.__nxg.add_node(name)

    def add_vertices(self, vertices):
        """
        Add vertices to the graph from an iterable container of vertices.
        """
        self.__nxg.add_nodes_from(vertices)

    def delete_vertex(self, vertex):
        """
        Deletes vertex, removing all incident edges.
        """
        self.__nxg.delete_node(vertex)

    def delete_vertices(self, vertices):
        """
        Remove vertices from the graph taken from an iterable container of
        vertices.
        """
        self.__nxg.delete_nodes_from(vertices)

    def has_vertex(self, vertex):
        return self.__nxg.has_node(vertex)

    def neighbor_iterator(self, vertex):
        """
        Return an iterator over neighbors of vertex.
        """
        return self.__nxg.neighbors_iter(vertex)

    def vertex_boundary(self, vertices1, vertices2=None):
        """
        Returns a list of all vertices in the external boundary of vertices1,
        intersected with vertices2. If vertices2 is None, then vertices2 is the
        complement of vertices1.
        """
        return self.__nxg.node_boundary(vertices1, vertices2)

    def vertex_iterator(self, vertices=None):
        """
        Returns an iterator over the given vertices. Returns False if not given
        a vertex, sequence, iterator or None. None is equivalent to a list of
        every vertex.
        """
        return self.__nxg.prepare_nbunch(vertices)

    def vertices(self):
        """
        Return a list of the vertex keys.
        """
        return self.__nxg.nodes()

    ### Edge Handlers

    def add_edge(self, u, v=None):
        """
        Adds an edge between u and v.

        INPUT:
        The following forms are all accepted by NetworkX:

        G.add_edge( 1, 2 )
        G.add_edge( (1, 2) )
        G.add_edges(  [ (1, 2) ]  )
        """
        self.__nxg.add_edge(u, v)

    def add_edges(self, edges):
        """
        Add edges from an iterable container.
        """
        self.__nxg.add_edges_from( edges )

    def delete_edge(self, u, v=None):
        r"""
        Delete the edge \{u, v\}, return silently if vertices or edge does not
        exist.
        """
        self.__nxg.delete_edge(u, v)

    def delete_edges(self, edges):
        """
        Delete edges from an iterable container.
        """
        self.__nxg.delete_edges_from(edges)

    def edges(self):
        """
        Return a list of edges.
        """
        return self.__nxg.edges()

    def edge_boundary(self, vertices1, vertices2=None):
        r"""
        Returns a list of edges \{u, v\} with u in vertices1 and v in vertices2.
        If vertices2 is None, then it is set to the complement of vertices1.
        """
        return self.__nxg.edge_boundary(vertices1, vertices2)

    def edge_iterator(self, vertices=None):
        """
        Returns an iterator over the edges incident with any vertex given.
        If vertices is None, then returns an iterator over all edges.
        """
        return self.__nxg.edges_iter(vertices)

    def edges_incident(self, vertices=None):
        """
        Returns a list of edges incident with any vertex given. If vertex is
        None, returns a list of all edges in graph.
        """
        return self.__nxg.edges(vertices)

    def has_edge(self, u, v=None):
        r"""
        Returns True if \{u, v\} is an edge, False otherwise.

        INPUT:
        The following forms are accepted by NetworkX:

        G.has_edge( 1, 2 )
        G.has_edge( (1, 2) )
        """
        return self.__nxg.has_edge(u, v)

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
        return self.__nxg.degree(vertices, with_labels)

    def degree_histogram(self):
        """
        Returns a list, whose ith entry is the frequency of degree i.
        """
        return NX.degree_histogram(self.__nxg)

    def degree_iterator(self, vertices=None, with_labels=False):
        """
        with_labels=False:
            returns an iterator over degrees.
        with_labels=True:
            returns an iterator over tuples (vertex, degree).
        """
        return self.__nxg.degree_iter(vertices, with_labels)

    ### Representations

    def adjacency_matrix(self, sparse=True):
        """
        Returns the adjacency matrix of the digraph. Each vertex is represented by its
        position in the list returned by the vertices() function.
        """
        n = len(self.__nxg.adj)
        verts = self.vertices()
        D = {}
        for e in self.edge_iterator():
            i,j = e
            i = verts.index(i)
            j = verts.index(j)
            D[(i,j)] = 1
            D[(j,i)] = 1
        M = matrix(IntegerModRing(2), n, n, D, sparse=sparse)
        return M

    def am(self):
        """
        Shorter call for adjacency matrix makes life easier.
        """
        return self.adjacency_matrix()

    ### Construction

    def add_cycle(self, vertices):
        self.__nxg.add_cycle(vertices)

    def add_path(self, vertices):
        self.__nxg.add_path(vertices)

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
            self.__nxg = self.__nxg.subgraph(vertices, inplace, create_using)
            return self
        else:
            NXG = self.__nxg.subgraph(vertices, inplace, create_using)
            return Graph(NXG, name=NXG.name) ### this bug is supposedly fixed

    ### Visualization

    def plot(self, pos=None,
                   with_labels=True,
                   node_size=200):
        GG = Graphics()
        if pos is None:
            if self.__pos is None:
                NGP = GraphicPrimitive_NetworkXGraph(self.__nxg, pos=None, with_labels=with_labels, node_size=node_size)
            else:
                NGP = GraphicPrimitive_NetworkXGraph(self.__nxg, pos=self.__pos, with_labels=with_labels, node_size=node_size)
        GG.append(NGP)
        GG.axes(False)
        return GG

    def show(self, pos=None,
                   with_labels=True,
                   node_size=200):
        self.plot(pos, with_labels, node_size).show()

class DiGraph(SimpleGraph):
    """
    Directed, simple graph (no loops, no multiple edges, non-hyper).
    """

    def __contains__(self, vertex):
        """
        Allows 'v in G'
        """
        return vertex in self.__nxg

    def __init__(self, data=None, pos=None, **kwds):
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
        pos -- an optional positioning dictionary: for example, the
        spring layout from NetworkX for the 5-cycle is
            {   0: [-0.91679746, 0.88169588,],
                1: [ 0.47294849, 1.125     ,],
                2: [ 1.125     ,-0.12867615,],
                3: [ 0.12743933,-1.125     ,],
                4: [-1.125     ,-0.50118505,]   }
        name -- (in kwds) gives the graph a name

        EXAMPLES:
        needed
        """
        if isinstance(data, DiGraph):
            self.__nxg = data.networkx_graph()
        elif isinstance(data, NX.DiGraph):
            self.__nxg = data
        else:
            self.__nxg = NX.DiGraph(data, **kwds)
        ### NOTE: Name bug in NetworkX supposedly fixed, so check if the fol-
        ### lowing fix is unnecessary
        self.__nxg.name=kwds.get("name","No Name")
        ### end fix
        self.__pos = pos

    def __iter__(self):
        """
        Return an iterator over the vertices. Allows 'for v in G'.
        """
        return self.vertex_iterator()

    def clear(self):
        """
        Empties the graph of vertices and edges, removes name.
        """
        self.__nxg.clear()

    def copy(self):
        """
        Creates a copy of the graph.
        """
        G = DiGraph(self.__nxg, name=self.__nxg.name)
        return G

    def is_directed(self):
        return True

    def networkx_graph(self):
        return self.__nxg.copy()

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
            while self.__nxg.succ.has_key(i):
                    i=i+1
            self.__nxg.add_node(i)
        else:
            self.__nxg.add_node(name)

    def add_vertices(self, vertices):
        """
        Add vertices to the graph from an iterable container of vertices.
        """
        self.__nxg.add_nodes_from(vertices)

    def delete_vertex(self, vertex):
        """
        Deletes vertex, removing all incident edges.
        """
        self.__nxg.delete_node(vertex)

    def delete_vertices(self, vertices):
        """
        Remove vertices from the graph taken from an iterable container of
        vertices.
        """
        self.__nxg.delete_nodes_from(vertices)

    def neighbor_iterator(self, vertex):
        """
        Return an iterator over neighbors (connected either way) of vertex.
        """
        A = list(self.__nxg.pred[vertex].iterkeys())
        B = list(self.__nxg.succ[vertex].iterkeys())
        C = []
        for V in A:
            if not V in B:
                C += [V]
        for V in B:
            C += [V]
        return iter(C)

    def vertex_iterator(self, vertices=None):
        """
        Returns an iterator over the given vertices. Returns False if not given
        a vertex, sequence, iterator or None. None is equivalent to a list of
        every vertex.
        """
        return self.__nxg.prepare_nbunch(vertices)

    def vertices(self):
        """
        Return a list of the vertex keys.
        """
        return self.__nxg.nodes()

    ### Arc Handlers

    def add_arc(self, u, v=None):
        """
        Adds an arc from u to v.

        INPUT:
        The following forms are all accepted by NetworkX:

        G.add_edge( 1, 2 )
        G.add_edge( (1, 2) )
        G.add_edges(  [ (1, 2) ]  )
        """
        self.__nxg.add_edge(u, v)

    def add_arcs(self, arcs):
        """
        Add arcs from an iterable container.
        """
        self.__nxg.add_edges_from( arcs )

    def arcs(self):
        """
        Return a list of arcs.
        """
        return self.__nxg.edges()

    def arc_iterator(self, vertices=None):
        """
        Returns an iterator over the arcs pointing out of the given set of vertices.
        If vertices is None, then returns an iterator over all arcs.
        """
        return self.__nxg.edges_iter(vertices)

    def delete_arc(self, u, v=None):
        r"""
        Delete the arc from u to v, return silently if vertices or edge does
        not exist.
        """
        self.__nxg.delete_edge(u, v)

    def delete_arcs(self, arcs):
        """
        Delete arcs from an iterable container.
        """
        self.__nxg.delete_edges_from(edges)

    def incoming_arc_iterator(self, vertices=None):
        """
        Return an iterator over all arriving arcs from vertices, or over all
        arcs if vertices is None.
        """
        return self.__nxg.in_edges_iter(vertices)

    def incoming_arcs(self, vertices=None):
        """
        Returns a list of arcs arriving at vertices.
        """
        return self.__nxg.in_edges(vertices)

    def outgoing_arc_iterator(self, vertices=None):
        """
        Return an iterator over all departing arcs from vertices, or over all
        arcs if vertices is None.
        """
        return self.__nxg.out_edges_iter(vertices)

    def outgoing_arcs(self, vertices=None):
        """
        Returns a list of arcs departing from vertices.
        """
        return self.__nxg.out_edges(vertices)

    def predecessor_iterator(self, vertex):
        """
        Returns an iterator over predecessor vertices of vertex.
        """
        return self.__nxg.predecessors_iter(vertex)

    def predecessors(self, vertex):
        return list(self.predecessor_iterator(vertex))

    def successor_iterator(self, vertex):
        """
        Returns an iterator over successor vertices of vertex.
        """
        return self.__nxg.successors_iter(vertex)

    def successors(self, vertex):
        return list(self.successor_iterator(vertex))

    ### Degree functions

    def degree_iterator(self, vertices=None, with_labels=False):
        """
        with_labels=False:
            returns an iterator over degrees (in + out).
        with_labels=True:
            returns an iterator over tuples (vertex, degree (in + out) ).
        """
        return self.__nxg.degree_iter(vertices, with_labels)

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
        return self.__nxg.degree(vertices, with_labels)

    def in_degree(self, vertices=None, with_labels=False):
        return self.__nxg.in_degree(vertices, with_labels)

    def in_degree_iterator(self, vertices=None, with_labels=False):
        """
        Same as degree_iterator, but for in degree.
        """
        return self.__nxg.in_degree_iter(vertices, with_labels)

    def out_degree(self, vertices=None, with_labels=False):
        return self.__nxg.out_degree(vertices, with_labels)

    def out_degree_iterator(self, vertices=None, with_labels=False):
        """
        Same as degree_iterator, but for out degree.
        """
        return self.__nxg.out_degree_iter(vertices, with_labels)

    ### Representations

    def adjacency_matrix(self, sparse=True):
        """
        Returns the adjacency matrix of the digraph. Each vertex is represented by its
        position in the list returned by the vertices() function.
        """
        n = len(self.__nxg.adj)
        verts = self.vertices()
        D = {}
        for e in self.arc_iterator():
            i,j = e
            i = verts.index(i)
            j = verts.index(j)
            D[(i,j)] = 1
        M = matrix(IntegerModRing(2), n, n, D, sparse=sparse)
        return M

    def am(self):
        """
        Shorter call for adjacency matrix makes life easier.
        """
        return self.adjacency_matrix()

    ### Contructors

    def reverse(self):
        """
        Returns a copy of digraph with arcs reversed in direction.
        """
        NXG = self.__nxg.reverse()
        G = DiGraph(NXG, name=NXG.name) # name thing: check bug
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
            self.__nxg = self.__nxg.subgraph(vertices, inplace, create_using)
            return self
        else:
            NXG = self.__nxg.subgraph(vertices, inplace, create_using)
            return DiGraph(NXG, name=NXG.name) ### this bug is supposedly fixed

    def to_directed(self):
        return self.copy()

    def to_undirected(self):
        NXG = self.__nxg.to_undirected()
        return DiGraph(NXG, name=NXG.name)

    ### Visualization

    def plot(self, pos=None,
                   with_labels=True,
                   node_size=200):
        GG = Graphics()
        if pos is None:
            if self.__pos is None:
                NGP = GraphicPrimitive_NetworkXGraph(self.__nxg, pos=None, with_labels=with_labels, node_size=node_size)
            else:
                NGP = GraphicPrimitive_NetworkXGraph(self.__nxg, pos=self.__pos, with_labels=with_labels, node_size=node_size)
        GG.append(NGP)
        GG.axes(False)
        return GG

    def show(self, pos=None,
                   with_labels=True,
                   node_size=200):
        self.plot(pos, with_labels, node_size).show()

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
