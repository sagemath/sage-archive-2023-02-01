"""
Implements various backends for Sage graphs.

"""



#*******************************************************************************
#        Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer

#from random import random
#from sage.plot.plot import Graphics, GraphicPrimitive_NetworkXGraph
#import sage.graphs.graph_fast as graph_fast
#from sage.rings.integer import Integer
#from sage.rings.integer_ring import ZZ
#from sage.graphs.graph_coloring import chromatic_number, chromatic_polynomial
#from sage.rings.rational import Rational

class GenericGraphBackend(SageObject):
    """
    A generic wrapper for the backend of a graph.  Various graph classes use
    extensions of this class.  Note, this graph has a number of placeholder
    functions, so the doctests are rather silly.

    DOCTEST:
        sage: import sage.graphs.base.graph_backends

    """
    _loops = False
    _multiple_edges = False
    _name = ''
    def add_edge(self, u, v, l, directed):
        """
        Add an edge (u,v) to self, with label l.  If directed is True, this is
        interpreted as an arc from u to v.

        INPUT:
            u,v:      vertices
            l:        edge label
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.add_edge(1,2,'a',True)
            Traceback (most recent call last):
            ...
            NotImplementedError
         """
        raise NotImplementedError()
    def add_edges(self, edges, directed):
        """
        Add a sequence of edges to self.  If directed is True, these are
        interpreted as arcs.

        INPUT:
            edges:    iterator
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.add_edges([],True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def add_vertex(self, name):
        """
        Add a labelled vertex to self.

        INPUT:
            name: vertex label

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.add_vertex(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def add_vertices(self, vertices):
        """
        Add labelled vertices to self.

        INPUT:
            vertices: iterator of vertex labels

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.add_vertices([1,2,3])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def degree(self, v, directed):
        """
        Returns the total number of vertices incident to v.

        INPUT:
            v:       a vertex label
            directed: boolean
        OUTPUT:
            degree of v

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.degree(1, False)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def del_edge(self, u, v, l, directed):
        """
        Deletes the edge (u,v) with label l.

        INPUT:
            u,v:      vertices
            l:        edge label
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.del_edge(1,2,'a',True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def del_vertex(self, v):
        """
        Delete a labelled vertex in self.

        INPUT:
            v: vertex label

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.del_vertex(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def del_vertices(self, vertices):
        """
        Delete labelled vertices in self.

        INPUT:
            vertices: iterator of vertex labels

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.del_vertices([1,2,3])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def get_edge_label(self, u, v):
        """
        Returns the edge label of (u,v).

        INPUT:
            u,v: vertex labels

        OUTPUT:
            label of (u,v)

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.get_edge_label(1,2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def has_edge(self, u, v, l):
        """
        True if self has an edge (u,v) with label l.

        INPUT:
            u,v: vertex labels
            l: label

        OUTPUT:
            boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.has_edge(1,2,'a')
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def has_vertex(self, v):
        """
        True if self has a vertex with label v.

        INPUT:
            v: vertex label

        OUTPUT:
            boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.has_vertex(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
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
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_edges([],True,True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
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
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_in_edges([],True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
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
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_out_edges([],True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def iterator_nbrs(self, v):
        """
        Iterate over the vertices adjacent to v.

        INPUT:
            v: vertex label

        OUTPUT:
            a generator which yields vertex labels

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_nbrs(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def iterator_in_nbrs(self, v):
        """
        Iterate over the vertices u such that the edge (u,v) is in self
        (that is, predecessors of v).

        INPUT:
            v: vertex label

        OUTPUT:
            a generator which yields vertex labels

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_in_nbrs(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def iterator_out_nbrs(self, v):
        """
        Iterate over the vertices u such that the edge (v,u) is in self
        (that is, successors of v).

        INPUT:
            v: vertex label

        OUTPUT:
            a generator which yields vertex labels

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_out_nbrs(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def iterator_verts(self, verts):
        """
        Iterate over the vertices v with labels in verts.

        INPUT:
            vertex: vertex labels

        OUTPUT:
            a generator which yields vertices

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_verts(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def loops(self, new):
        """
        Get/set whether or not self allows loops.

        INPUT:
            new: boolean or None

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.loops(True)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: G.loops(None)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def multiple_edges(self, new):
        """
        Get/set whether or not self allows multiple edges.

        INPUT:
            new: boolean or None

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.multiple_edges(True)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: G.multiple_edges(None)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def name(self, new):
        """
        Get/set name of self.

        INPUT:
            new: string or None

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.name("A Generic Graph")
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: G.name(None)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def num_edges(self, directed):
        """
        The number of edges in self

        INPUT:
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.num_edges(True)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: G.num_edges(False)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def num_verts(self):
        """
        The number of vertices in self

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.num_verts()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def relabel(self, perm, directed):
        """
        Relabel the vertices of self by a permutation.
        INPUT:
            perm:     permutation
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.relabel([],False)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def set_edge_label(self, u, v, l, directed):
        """
        Label the edge (u,v) by l.

        INPUT:
            u,v:      vertices
            l:        edge label
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.set_edge_label(1,2,'a',True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()

class NetworkXGraphBackend(GenericGraphBackend):
    """
    A wrapper for NetworkX as the backend of a graph.

    DOCTEST:
        sage: import sage.graphs.base.graph_backends

    """

    _nxg = None

    def __init__(self, N=None):
        """
        Initialize the backend with NetworkX graph N.

        DOCTEST:
        """
        if N is None:
            import networkx
            N = networkx.XGraph()
        self._nxg = N

    def add_edge(self, u, v, l, directed):
        """
        Add an edge (u,v) to self, with label l.  If directed is True, this is
        interpreted as an arc from u to v.

        INPUT:
            u,v:      vertices
            l:        edge label
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.add_edge(1,2,'a',True)
         """
        self._nxg.add_edge(u, v, l)

    def add_edges(self, edges, directed):
        """
        Add a sequence of edges to self.  If directed is True, these are
        interpreted as arcs.

        INPUT:
            edges:    iterator
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.add_edges([],True)
        """
        for e in edges:
            try:
                u,v,l = e
            except:
                u,v = e
                l = None
            self.add_edge(u,v,l,directed)

    def add_vertex(self, name):
        """
        Add a labelled vertex to self.

        INPUT:
            name: vertex label

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.add_vertex(0)
        """
        if name is None: # then find an integer to use as a key
            i = 0
            while self.has_vertex(i):
                i=i+1
            name = i
        self._nxg.add_node(name)

    def add_vertices(self, vertices):
        """
        Add labelled vertices to self.

        INPUT:
            vertices: iterator of vertex labels

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.add_vertices([1,2,3])
        """
        for v in vertices:
            self.add_vertex(v)

    def degree(self, v, directed):
        """
        Returns the total number of vertices incident to v.

        INPUT:
            v:       a vertex label
            directed: boolean
        OUTPUT:
            degree of v

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.add_vertices(range(3))
            sage: G.degree(1, False)
            0
        """
        return self._nxg.degree(v)

    def del_edge(self, u, v, l, directed):
        """
        Deletes the edge (u,v) with label l.

        INPUT:
            u,v:      vertices
            l:        edge label
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.del_edge(1,2,'a',True)
        """
        self._nxg.delete_edge(u, v, l)

    def del_vertex(self, v):
        """
        Delete a labelled vertex in self.

        INPUT:
            v: vertex label

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.del_vertex(0)
            Traceback (most recent call last):
            ...
            NetworkXError: node 0 not in graph
        """
        self._nxg.delete_node(v)

    def del_vertices(self, vertices):
        """
        Delete labelled vertices in self.

        INPUT:
            vertices: iterator of vertex labels

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.del_vertices([1,2,3])
            Traceback (most recent call last):
            ...
            NetworkXError: node 1 not in graph
        """
        for v in vertices:
            self.del_vertex(v)

    def get_edge_label(self, u, v):
        """
        Returns the edge label of (u,v).

        INPUT:
            u,v: vertex labels

        OUTPUT:
            label of (u,v)

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.get_edge_label(1,2)
            Traceback (most recent call last):
            ...
            NetworkXError: Edge (1,2) requested via get_edge does not exist.
        """
        return self._nxg.get_edge(u, v)

    def has_edge(self, u, v, l):
        """
        True if self has an edge (u,v) with label l.

        INPUT:
            u,v: vertex labels
            l: label

        OUTPUT:
            boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.has_edge(1,2,'a')
            False
        """
        return self._nxg.has_edge(u, v, l)

    def has_vertex(self, v):
        """
        True if self has a vertex with label v.

        INPUT:
            v: vertex label

        OUTPUT:
            boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.has_vertex(0)
            False
        """
        return self._nxg.has_node(v)

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
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.iterator_edges([],True,True)
            <generator object at ...>
        """
        if labels:
            filter = lambda u,v,l: (u,v,l)
        else:
            filter = lambda u,v,l: (u,v)
        if not_directed:
            for u,v,l in self._nxg.edges_iter():
                if u in vertices or v in vertices:
                    yield filter(u,v,l)
        else:
            for u,v,l in self._nxg.edges_iter(vertices):
                yield filter(u,v,l)

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
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.iterator_in_edges([],True)
            <generator object at ...>
        """
        if labels:
            for e in self._nxg.in_edges_iter(vertices):
                yield e
        else:
            for u,v,l in self._nxg.in_edges_iter(vertices):
                yield (u,v)

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
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.iterator_out_edges([],True)
            <generator object at ...>
        """
        if labels:
            for e in self._nxg.out_edges_iter(vertices):
                yield e
        else:
            for u,v,l in self._nxg.out_edges_iter(vertices):
                yield (u,v)

    def iterator_nbrs(self, v):
        """
        Iterate over the vertices adjacent to v.

        INPUT:
            v: vertex label

        OUTPUT:
            a generator which yields vertex labels

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.iterator_nbrs(0)
            <generator object at ...>
        """
        return self._nxg.neighbors_iter(v)

    def iterator_in_nbrs(self, v):
        """
        Iterate over the vertices u such that the edge (u,v) is in self
        (that is, predecessors of v).

        INPUT:
            v: vertex label

        OUTPUT:
            a generator which yields vertex labels

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.iterator_in_nbrs(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'XGraph' object has no attribute 'predecessors_iter'
        """
        return self._nxg.predecessors_iter(v)

    def iterator_out_nbrs(self, v):
        """
        Iterate over the vertices u such that the edge (v,u) is in self
        (that is, successors of v).

        INPUT:
            v: vertex label

        OUTPUT:
            a generator which yields vertex labels

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.iterator_out_nbrs(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'XGraph' object has no attribute 'successors_iter'
        """
        return self._nxg.successors_iter(v)

    def iterator_verts(self, verts):
        """
        Iterate over the vertices v with labels in verts.

        INPUT:
            vertex: vertex labels

        OUTPUT:
            a generator which yields vertices

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.iterator_verts(0)
            <listiterator object at ...>
        """
        return iter(self._nxg.prepare_nbunch(verts))

    def loops(self, new):
        """
        Get/set whether or not self allows loops.

        INPUT:
            new: boolean or None

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.loops(True)
            sage: G.loops(None)
            True
        """
        if new is None:
            return self._nxg.selfloops
        if new:
            self._nxg.allow_selfloops()
        else:
            self._nxg.ban_selfloops()

    def multiple_edges(self, new):
        """
        Get/set whether or not self allows multiple edges.

        INPUT:
            new: boolean or None

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.multiple_edges(True)
            sage: G.multiple_edges(None)
            True
        """
        if new is None:
            return self._nxg.multiedges
        if new:
            self._nxg.allow_multiedges()
        else:
            self._nxg.ban_multiedges()

    def name(self, new):
        """
        Get/set name of self.

        INPUT:
            new: string or None

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.name("A NetworkX Graph")
            sage: G.name(None)
            'A NetworkX Graph'
        """
        if new is None:
            return self._nxg.name
        self._nxg.name = new

    def num_edges(self, directed):
        """
        The number of edges in self

        INPUT:
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.num_edges(True)
            0
            sage: G.num_edges(False)
            0
        """
        return self._nxg.size()

    def num_verts(self):
        """
        The number of vertices in self

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.num_verts()
            0
        """
        return self._nxg.order()

    def relabel(self, perm, directed):
        """
        Relabel the vertices of self by a permutation.
        INPUT:
            perm:     permutation
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.relabel([],False)
        """
        if directed:
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
        else:
            oldd = self._nxg.adj
            newd = {}
            for v in oldd.iterkeys():
                oldtempd = oldd[v]
                newtempd = {}
                for w in oldtempd.iterkeys():
                    newtempd[perm[w]] = oldtempd[w]
                newd[perm[v]] = newtempd
            self._nxg.adj = newd

    def set_edge_label(self, u, v, l, directed):
        """
        Label the edge (u,v) by l.

        INPUT:
            u,v:      vertices
            l:        edge label
            directed: boolean

        DOCTEST:
            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.set_edge_label(1,2,'a',True)
        """
        if not self.has_edge(u, v, None):
            return
        if self.multiple_edges(None):
            if directed:
                self._nxg.succ[u][v] = [l]
                self._nxg.pred[v][u] = [l]
            else:
                self._nxg.adj[u][v] = [l]
                self._nxg.adj[v][u] = [l]
        else:
            if directed:
                self._nxg.succ[u][v] = l
                self._nxg.pred[v][u] = l
            else:
                self._nxg.adj[u][v] = l
                self._nxg.adj[v][u] = l







