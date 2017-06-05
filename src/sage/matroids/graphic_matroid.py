r"""
Graphic Matroids

Theory
======

Let `G = (V,E)` be a graph and let `C` be the collection of the edge sets of cycles in `G`. The corresponding graphic matroid `M(G)` has ground set `E` and circuits `C`.
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
from .utilities import sanitize_contractions_deletions, setprint_s

#I'll put this here for now but I suspect it belongs in another file.
def contract_edge(G, edgelist):
    """
    Contract an edge or a list of edges.
    """
    G.allow_multiple_edges(True)
    G.allow_loops(True)
    if edgelist in G.edges():
        edgelist = [edgelist]

    for e in edgelist:
        if e not in G.edges():
            raise ValueError("The specified edge is not in the graph.")
    G.delete_edges(edgelist)
    #If e was a loop, stop there. Otherwise, merge the vertices.
    if e[0] == e[1]:
        pass
    else:
        # merge_vertices() loses multiedges, so we put them on as loops afterwards
        edge_label_list = []
        for edge in G.edges():
            if (edge[0] == e[0] and edge[1] == e[1]) or (edge[0] == e[1] and edge[1] == e[0]):
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
        self._H = Graph()

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

        self._new_G = copy(self._G)

        self._new_groundset_edge_map = copy(self._groundset_edge_map)
        edge_list = []
        for x in contractions:
            edge_list.append((self._new_groundset_edge_map[x][0],
                              self._new_groundset_edge_map[x][1], x))
            #Putting the label on the edge make sure the correct element is removed.
            #Without this the result would be isomorphic but the labels would be wrong.

            contract_edge(self._new_G, edge_list)

            #Since this changes vertex labels, I need to update the groundset edge map.
            for key in self._new_groundset_edge_map.keys():
                if self._new_groundset_edge_map[key][0] == self._edge[1]:
                    self._new_groundset_edge_map[key] = (self._edge[0],
                                                         self._new_groundset_edge_map[key][1])
                if self._new_groundset_edge_map[key][1] == self._edge[1]:
                    self._new_groundset_edge_map[key] = (self._new_groundset_edge_map[key][0],
                                                         self._edge[0])

        edge_list = []
        for x in deletions:
            edge_list.append((self._new_groundset_edge_map[x][0],
                              self._new_groundset_edge_map[x][1], x))
            self._new_G.delete_edges(edge_list)

        return GraphicMatroid(deepcopy(self._new_G))