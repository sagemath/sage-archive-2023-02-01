r"""
Graphic Matroids

Theory
======

Let `G = (V,E)` be a graph and let `C` be the collection of the edge sets
of cycles in `G`. The corresponding graphic matroid `M(G)` has ground set `E`
and circuits `C`.
"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2017 Zachary Gershkoff <zgersh2@lsu.edu>
#       Copyright (C) 2017 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from .matroid import Matroid

from sage.graphs.graph import Graph
from copy import copy, deepcopy

#I'll put this here for now but I suspect it belongs in another file.
def contract_edge(G, e):
    """
    Contract an element of a graphic matroid.
    """
    G.allow_multiple_edges(True)
    G.allow_loops(True)

    if e not in G.edges():
        raise ValueError("The specified edge is not in the graph.")

    G.delete_edge(e)

    #If e was a loop, stop there. Otherwise, merge the vertices.
    if not e[0] == e[1]:
        # merge_vertices() loses multiedges, so we put them on as loops afterwards
        edge_label_list = []
        for edge in G.edges():
            if ((edge[0] == e[0] and edge[1] == e[1]) or
                (edge[0] == e[1] and edge[1] == e[0])):
                edge_label_list.append(edge[2])

        G.merge_vertices([e[0], e[1]])

        for edge_label in edge_label_list:
            G.add_edge(e[0], e[0], edge_label)

class GraphicMatroid(Matroid):
    """
    The graphic matroid class.

    INPUT:
    G -- a SageMath graph
    groundset (optional) -- a list in 1-1 correspondence with G.edges()
    """

    def __init__(self, G, groundset = None):


        if groundset is None:
            #Try to construct a ground set based on the edge labels.
            #If that fails, use range() to come up with a groundset.
            groundset = G.edge_labels()

        groundset_set = frozenset(groundset)

        #if the provided ground set is incomplete, it gets overwriten
        if len(groundset_set) != G.num_edges():
            groundset = range(G.num_edges())
            groundset_set = frozenset(groundset)

        self._groundset = groundset_set

        #Map ground set elements to graph edges:
        self._groundset_edge_map = {x: y for (x,y) in zip(groundset, G.edges(labels=False))}
        #Construct a graph and assign edge labels corresponding to the ground set
        edge_list = []
        for i, e in enumerate(G.edges()):
            edge_list.append((e[0],e[1],groundset[i]))
        self._G = Graph(edge_list, loops=True, multiedges=True)

    def _rank(self, X):
        """
        Return the rank of a set `X`.

        This method does no checking on `X`, and
        `X` may be assumed to have the same interface as `frozenset`.

        INPUT:

        - `X` -- an object with Python's `frozenset` interface.

        OUTPUT:

        The rank of `X` in the matroid.
        """
        #we get ValueErrors if loops and multiedges are not enabled
        self._H = Graph(loops=True, multiedges=True)

        for x in X:
            self._H.add_edge(self._groundset_edge_map[x])
            if not self._H.is_forest():
                self._H.delete_edge(self._groundset_edge_map[x])

        return len(self._H.edges())

    def groundset(self):
        """
        Returns the ground set of the matroid as a frozenset.
        """
        return self._groundset

    def graph(self):
        """
        Returns a graph that has a cycle matroid equal to the matroid.
        """
        return copy(self._G)

    def _repr_(self):
        """
        Returns a string representation of the matroid.
        """
        self._mrank = str(self._rank(self._groundset))
        self._elts = str(len(self._groundset))

        return "A graphic matroid of rank " + self._mrank + " on " + self._elts + " elements."

    def _minor(self, contractions=frozenset([]), deletions=frozenset([])):
        """
        Return a minor.
        INPUT:
        contractions -- frozenset, subset of self.groundset() to be contracted
        deletions -- frozenset, subset of self.groundset() to be deleted

        Assumptions: contractions are independent, deletions are coindependent,
        contractions and deletions are disjoint.
        """
        g = self.graph()

        #self._new_groundset_edge_map = copy(self._groundset_edge_map)
        cont_edges = self._groundset_to_edges(contractions)
        del_edges = self._groundset_to_edges(deletions)
        #for e in cont_edges:
        while cont_edges != []:
            e = cont_edges.pop()
            contract_edge(g, e)
            del_edges = self._update_edges(e, del_edges)
            cont_edges = self._update_edges(e, cont_edges)

        g.delete_edges(del_edges)
        #for x in deletions:
            #e = (self._new_groundset_edge_map[x][0], self._new_groundset_edge_map[x][1], x)

            #self._new_G.delete_edge(self._edge)

        return GraphicMatroid(deepcopy(g))

    def _update_edges(self, edge, edges):
        """
        After a contraction, updates a list of edges to exclude the vertex
        that was removed.
        """
        v0 = edge[0]
        v1 = edge[1]
        new_edges = []
        for e in edges:
            if e[0] == v1:
                if v0 <= e[1]:
                    e = (v0, e[1], e[2])
                else:
                    e = (e[1], v0, e[2])
            if e[1] == v1:
                if v0 >= e[0]:
                    e = (e[0], v0, e[2])
                else:
                    e = (v0, e[0], e[2])
            new_edges.append(e)
        return new_edges

    def _has_minor(self, N, certificate = False):
        """
        Checks if the matroid has a minor isomoprhic to M(H).
        """
        # The graph minor algorithm is faster but it doesn't make sense
        # to use it if M(H) is not 3-connected, because of all the possible
        # Whitney switches or 1-sums that will give the same matroid.
        if isinstance(N,GraphicMatroid) and N.is_3connected():
            # Graph.minor() does not work with multigraphs
            G = self.graph()
            G.allow_loops(False)
            G.allow_multiple_edges(False)
            H = N.graph()
            H.allow_loops(False)
            H.allow_multiple_edges(False)

            try:
                # Graph.minor() returns a certificate if there is one
                # and a ValueError if there isn't.
                cert = G.minor(H)
                if certificate:
                    return cert
                else:
                    return True
            except ValueError:
                return False
        else:
            # otherwise use the default method for abstract matroids
            return Matroid._has_minor(self,N)

    def _groundset_to_edges(self, X):
        """
        Given a subset of the ground set, this will return the corresponding
        set of edges.
        """
        edge_list = []
        for x in X:
            v0 = self._groundset_edge_map[x][0]
            v1 = self._groundset_edge_map[x][1]
            edge_list.append((v0, v1, x))
        return edge_list

    def _subgraph_from_set(self,X):
        """
        Returns the subgraph induced by the edges corresponding to the elements of X.
        """
        edge_list = self._groundset_to_edges(X)
        return Graph(edge_list, loops=True, multiedges=True)

    def _corank(self, X):
        """
        Returns the corank of the set X in the matroid.
        """
        components = self._G.connected_components_number()
        g = self.graph()
        g.delete_edges(self._groundset_to_edges(X))
        return (len(X) - (g.connected_components_number() - components))

    def _is_independent(self, X):
        """
        Tests if the set is an independent set of the matroid.
        """
        g = self._subgraph_from_set(X)
        return g.is_forest()

    def _is_circuit(self, X):
        """
        Tests if the given set is a circuit.
        """
        g = self._subgraph_from_set(X)
        return g.is_cycle()

    def _closure(self, X):
        """
        Returns the closure of a set.
        """
        X = set(X)
        Y = self.groundset().difference(X)
        edgelist = self._groundset_to_edges(Y)
        g = self._subgraph_from_set(X)
        V = g.vertices()
        components = g.connected_components_number()
        for e in edgelist:
            # an edge is in the closure iff both its vertices are
            # in the induced subgraph, and the edge doesn't connect components
            if e[0] in V and e[1] in V:
                g.add_edge(e)
                if g.connected_components_number() >= components:
                    X.add(e[2])
                else:
                    g.delete_edge(e)
        return frozenset(X)

    def _max_independent(self,X):
        """
        Returns a maximum independent subset of a set.
        """
        res = set()
        g = self.graph()
        edgelist = self._groundset_to_edges(X)
        #for e in edgelist:
        while edgelist != []:
            e = edgelist.pop()
            if e not in g.loops():
                contract_edge(g,e)
                edgelist = self._update_edges(e, edgelist)
                res.add(e[2])
        return frozenset(res)

    def _max_coindependent(self,X):
        """
        Returns a maximum coindependent subset of a set.
        """
        res = set()
        components = self._G.connected_components_number()
        g = self.graph()
        edgelist = self._groundset_to_edges(X)
        for e in edgelist:
            g.delete_edge(e)
            if g.connected_components_number() > components:
                g.add_edge(e)
            else:
                res.add(e[2])
        return frozenset(res)
