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
from sage.graphs.base.all import SparseGraph, DenseGraph

#from random import random
#from sage.plot.plot import Graphics, GraphicPrimitive_NetworkXGraph
#import sage.graphs.graph_fast as graph_fast
#from sage.rings.integer import Integer
#from sage.rings.integer_ring import ZZ
#from sage.graphs.graph_coloring import chromatic_number, chromatic_polynomial
#from sage.rings.rational import Rational

class GenericGraphBackend(SageObject):
    _loops = False
    _multiple_edges = False
    _name = ''
    def add_edge(self, *args, **kwds):
        raise NotImplementedError()
    def add_vertex(self, *args, **kwds):
        raise NotImplementedError()
    def degree(self, *args, **kwds):
        raise NotImplementedError()
    def del_edge(self, *args, **kwds):
        raise NotImplementedError()
    def del_vertex(self, *args, **kwds):
        raise NotImplementedError()
    def get_edge_label(self, *args, **kwds):
        raise NotImplementedError()
    def has_edge(self, *args, **kwds):
        raise NotImplementedError()
    def has_vertex(self, *args, **kwds):
        raise NotImplementedError()
    def iterator_edges(self, *args, **kwds):
        raise NotImplementedError()
    def iterator_in_edges(self, *args, **kwds):
        raise NotImplementedError()
    def iterator_out_edges(self, *args, **kwds):
        raise NotImplementedError()
    def iterator_nbrs(self, *args, **kwds):
        raise NotImplementedError()
    def iterator_in_nbrs(self, *args, **kwds):
        raise NotImplementedError()
    def iterator_out_nbrs(self, *args, **kwds):
        raise NotImplementedError()
    def iterator_verts(self, *args, **kwds):
        raise NotImplementedError()
    def loops(self, *args, **kwds):
        raise NotImplementedError()
    def multiple_edges(self, *args, **kwds):
        raise NotImplementedError()
    def name(self, *args, **kwds):
        raise NotImplementedError()
    def num_edges(self, *args, **kwds):
        raise NotImplementedError()
    def num_verts(self, *args, **kwds):
        raise NotImplementedError()
    def relabel(self, *args, **kwds):
        raise NotImplementedError()
    def set_edge_label(self, *args, **kwds):
        raise NotImplementedError()

class NetworkXGraphBackend(GenericGraphBackend):
    _nxg = None

    def add_edge(self, *args, **kwds):
        self._nxg.add_edge(args[0], args[1], args[2])

    def add_vertex(self, *args, **kwds):
        name = args[0]
        if name is None: # then find an integer to use as a key
            i = 0
            while self.has_vertex(i):
                i=i+1
            name = i
        self._nxg.add_node(name)

    def degree(self, *args, **kwds):
        return self._nxg.degree(args[0])

    def del_edge(self, *args, **kwds):
        self._nxg.delete_edge(args[0], args[1], args[2])

    def del_vertex(self, *args, **kwds):
        self._nxg.delete_node(args[0])

    def get_edge_label(self, *args, **kwds):
        return self._nxg.get_edge(args[0],args[1])

    def has_edge(self, *args, **kwds):
        return self._nxg.has_edge(args[0], args[1], args[2])

    def has_vertex(self, *args, **kwds):
        return self._nxg.has_node(args[0])

    def iterator_edges(self, *args, **kwds):
        vertices, labels, ignore_direction = args
        if labels:
            filter = lambda u,v,l: (u,v,l)
        else:
            filter = lambda u,v,l: (u,v)
        if ignore_direction:
            for u,v,l in self._nxg.edges_iter():
                if u in vertices or v in vertices:
                    yield filter(u,v,l)
        else:
            for u,v,l in self._nxg.edges_iter(vertices):
                yield filter(u,v,l)

    def iterator_in_edges(self, *args, **kwds):
        if args[1]:
            for e in self._nxg.in_edges_iter(args[0]):
                yield e
        else:
            for u,v,l in self._nxg.in_edges_iter(args[0]):
                yield (u,v)

    def iterator_out_edges(self, *args, **kwds):
        if args[1]:
            for e in self._nxg.out_edges_iter(args[0]):
                yield e
        else:
            for u,v,l in self._nxg.out_edges_iter(args[0]):
                yield (u,v)

    def iterator_nbrs(self, *args, **kwds):
        return self._nxg.neighbors_iter(args[0])

    def iterator_in_nbrs(self, *args, **kwds):
        return self._nxg.predecessors_iter(args[0])

    def iterator_out_nbrs(self, *args, **kwds):
        return self._nxg.successors_iter(args[0])

    def iterator_verts(self, *args, **kwds):
        return self._nxg.prepare_nbunch(args[0])

    def loops(self, *args, **kwds):
        if args[0] is None:
            return self._nxg.selfloops
        if args[0]:
            self._nxg.allow_selfloops()
        else:
            self._nxg.ban_selfloops()

    def multiple_edges(self, *args, **kwds):
        if args[0] is None:
            return self._nxg.multiedges
        if args[0]:
            self._nxg.allow_multiedges()
        else:
            self._nxg.ban_multiedges()

    def name(self, *args, **kwds):
        if args[0] is None:
            return self._nxg.name
        self._nxg.name = args[0]

    def num_edges(self, *args, **kwds):
        return self._nxg.size()

    def num_verts(self, *args, **kwds):
        return self._nxg.order()

    def relabel(self, *args, **kwds):
        perm, directed = args[0], args[1]
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

    def set_edge_label(self, *args, **kwds):
        u, v, l, directed = args
        if not self.has_edge(u, v, None):
            return
        if self.multiple_edges(None):
            if len(self.get_edge_label(u, v)) > 1:
                raise RuntimeError("Cannot set edge label, since there are multiple edges from %s to %s."%(u,v))
            if directed:
                self._nxg.adj[u][v] = [l]
            else:
                self._nxg.adj[u][v] = [l]
                self._nxg.adj[v][u] = [l]
        else:
            if directed:
                self._nxg.adj[u][v] = l
            else:
                self._nxg.adj[u][v] = l
                self._nxg.adj[v][u] = l

class DenseGraphBackend(GenericGraphBackend):
    pass

class SparseGraphBackend(GenericGraphBackend):
    pass
