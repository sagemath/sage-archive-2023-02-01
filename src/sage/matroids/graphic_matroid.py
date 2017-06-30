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
from .utilities import newlabel, contract_edge, update_edges
from itertools import combinations
import random

class GraphicMatroid(Matroid):
    """
    The graphic matroid class.

    INPUT:
    G -- a SageMath graph
    groundset (optional) -- a list in 1-1 correspondence with G.edges()

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: edgelist = [('a','b'),('b','c'),('c','d'),('d','a')]
        sage: M = GraphicMatroid(Graph(edgelist))
        sage: M.graph().edges()
        [(0, 1, 0), (0, 3, 1), (1, 2, 2), (2, 3, 3)]
        sage: M.is_isomorphic(Matroid(graphs.CycleGraph(4)))
        True

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

        # Force vertices to be integers
        vertex_to_integer_map = {vertex: integer for (vertex, integer) in zip(
            G.vertices(), range(len(G.vertices())))}

        # Construct a graph and assign edge labels corresponding to the ground set
        edge_list = []
        for i, e in enumerate(G.edges()):
            edge_list.append((vertex_to_integer_map[e[0]],
                vertex_to_integer_map[e[1]],groundset[i]))
        self._G = Graph(edge_list, loops=True, multiedges=True)

        # The graph should be connected to make computations easier
        while self._G.connected_components_number() > 1:
            # This will give a list containing lists of vertices
            comps = self._G.connected_components()
            # Choose random vertices to avoid the graphs being plotted
            # as bouquets
            v1 = random.choice(comps[0])
            v2 = random.choice(comps[1])
            self._G.add_edge((v1, v2, None))
            contract_edge(self._G, (v1, v2, None))

        # Map ground set elements to graph edges:
        # The the edge labels should already be the elements.
        self._groundset_edge_map = ({l: (u, v) for
            (u, v, l) in self._G.edges()})

    #COPYING, LOADING, SAVING

    def __copy__(self):
        """
        Create a shallow copy.
        """
        return GraphicMatroid(self._G)
        # something about name wrangling?



    def _rank(self, X):
        """
        Return the rank of a set `X`.

        This method does no checking on `X`, and
        `X` may be assumed to have the same interface as `frozenset`.

        INPUT:

        - `X` -- an object with Python's `frozenset` interface.

        OUTPUT:

        The rank of `X` in the matroid.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: edgelist = [(0,0,0), (0,1,1), (0,2,2), (0,3,3), (1,2,4), (1,3,5)]
            sage: M = GraphicMatroid(Graph(edgelist, loops=True, multiedges=True))
            sage: M.rank([0])
            0
            sage: M.rank([1,2])
            2
            sage: M.rank([1,2,4])
            2
            sage: M.rank(M.groundset()) == M.full_rank()
            True
            sage: edgelist = [(0,0,0), (1,2,1), (1,2,2), (2,3,3)]
            sage: M = GraphicMatroid(Graph(edgelist, loops=True, multiedges=True))
            sage: M.rank(M.groundset())
            2
            sage: M.rank([0,3])
            1

        """
        edges = self.groundset_to_edges(X)
        vertices = set([u for (u, v, l) in edges]).union(
            set([v for (u, v, l) in edges]))
        # This counts components:
        from sage.sets.disjoint_set import DisjointSet
        DS_vertices = DisjointSet(vertices)
        for (u, v, l) in edges:
            DS_vertices.union(u,v)
        return (len(vertices) - DS_vertices.number_of_subsets())

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
            del_edges = update_edges(e, del_edges)
            cont_edges = update_edges(e, cont_edges)

        g.delete_edges(del_edges)
        #for x in deletions:
            #e = (self._new_groundset_edge_map[x][0], self._new_groundset_edge_map[x][1], x)

            #self._new_G.delete_edge(self._edge)

        return GraphicMatroid(deepcopy(g))

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

    def groundset_to_edges(self, X):
        """
        Given a subset of the ground set, this will return the corresponding
        set of edges.
        """
        for x in X:
            if x not in self._groundset:
                raise ValueError("input must be a subset of the ground set")
        return self._groundset_to_edges(X)

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
                edgelist = update_edges(e, edgelist)
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

    def _circuit(self, X):
        """
        Returns a minimal dependent subset.
        """
        circuit = set()
        g = self._subgraph_from_set(X)
        for e in g.edges():
            g.delete_edge(e)
            if g.is_forest():
                g.add_edge(e)
                circuit.add(e[2])
        return frozenset(circuit)

    def _coclosure(self, X):
        """
        Returns the coclosure of a set.
        """
        g = self.graph()
        g.delete_edges(self._groundset_to_edges(X))
        components = g.connected_components_number()
        X = set(X)
        Y = self.groundset().difference(X)
        for e in self._groundset_to_edges(Y):
            g.delete_edge(e)
            if g.connected_components_number() > components:
                X.add(e[2])
            g.add_edge(e)
        return frozenset(X)

    def _is_coindependent(self, X):
        """
        Tests if input is coindependent.
        """
        g = self.graph()
        components = g.connected_components_number()
        g.delete_edges(self._groundset_to_edges(X))
        if g.connected_components_number() == components:
            return True
        else:
            return False

    def _cocircuit(self, X):
        """
        Returns a minimal codependent subset.
        """
        cocircuit = set()
        codependent_edges = []
        g = self.graph()
        edges = self._groundset_to_edges(X)
        components = g.connected_components_number()
        for e in edges:
            codependent_edges.append(e)
            g.delete_edge(e)
            if g.connected_components_number() > components:
                break
        for e in codependent_edges:
            g.add_edge(e)
            # if that repaired the components, then e is part of the cocircuit
            if g.connected_components_number() == components:
                cocircuit.add(e[2])
            g.delete_edge(e)
        return frozenset(cocircuit)

    def _is_cocircuit(self, X):
        """
        Tests if the input is a cocircuit.
        """
        edges = self._groundset_to_edges(X)
        g = self.graph()
        components = g.connected_components_number()
        g.delete_edges(edges)
        # This should have made exactly 1 more component
        if not g.connected_components_number() == (components + 1):
            return False
        for e in edges:
            g.add_edge(e)
            # Every time an edge is added, it should repair the components
            if not g.connected_components_number() == components:
                return False
            g.delete_edge(e)
        return True

    def _is_closed(self, X):
        """
        Test if input is a closed set.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        Boolean.
        """
        # Take the set of vertices of the edges corresponding to the elements,
        # and check if there are other edges incident with two of those vertices.
        # Also, the must not be loops outside of X.
        vertex_set = set()
        edge_list = self._groundset_to_edges(X)
        if not set(self._G.loops()).issubset(set(edge_list)):
            return False

        Y = self.groundset().difference(X)
        edge_list2 = self._groundset_to_edges(Y)
        for e in edge_list:
            vertex_set.add(e[0])
            vertex_set.add(e[1])
        for e in edge_list2:
            if e[0] in vertex_set and e[1] in vertex_set:
                return False
        return True

    def _is_isomorphic(self, other, certificate = False):
        """
        Test if ``self`` is isomorphic to ``other``.

        INPUT:

        - ``other`` -- A matroid.
        - ``certificate`` -- Boolean

        OUTPUT:

        - If ``certificate`` is ``False``, Boolean.
        - If ``certificate`` is ``True``, a tuple containing a boolean and a dictionary
          giving the isomorphism or None.
        """
        if isinstance(other,GraphicMatroid) and other.is_3connected():
            # Graph.is_isomorphic() supports multigraphs
            # This could be made faster by using self._G instead of self.graph()
            G = self.graph()
            H = other.graph()
            return G.is_isomorphic(H, certificate=certificate)
        else:
            return Matroid._is_isomorphic(self, other, certificate=certificate)

    def _isomorphism(self, other):
        """
        Return isomorphism from ``self`` to ``other``, if such an isomorphism exists.
        """
        # TODO: If M is M(K_4) and N is M(W_3), then M._isomorphism(N) = True
        # but N._isomorphism(M) gives ImportError: No module named basis_matroid
        # from basis_exchange_matroid.pyx
        if isinstance(other,GraphicMatroid) and other.is_3connected():
            G = self.graph()
            H = other.graph()
            return G.is_isomorphic(H, certificate=True)[1]
        else:
            return Matroid._isomorphism(self, other)

    def graphic_extension(self, u, v=None, element=None):
        """
        Returns a graphic matroid with a specified element added.
        """
        # TODO: Coloops should be handled in a more logical way
        # Perhaps make a new method to connect/disconnect matroid components
        # in the graph
        if element is None:
            element = newlabel(self.groundset())
        elif element in self.groundset():
            raise ValueError("cannot extend by element already in groundset")
        # If u or v are not already vertices, the graph package will
        # make them into vertices
        # If v is None, make a loop at u, not a coloop
        # Since this is extension, not coextension.
        if v is None:
            v = u
        G = self.graph()
        G.add_edge(u,v,element)
        return GraphicMatroid(deepcopy(G))

    def graphic_extensions(self, element=None, vertices=None):
        """
        Returns an iterable containing the graphic extensions of self
        by an element.

        INPUT:

        - ``element`` -- (optional) The name of the newly added element in
          each extension.
        - ``vertices`` -- (optional) A set of vertices over which the extension
          may be taken. If not given, will use all vertices.

        OUTPUT:

        An iterable containing matroids.

        .. NOTE::

            The extension by a loop will always occur.
            The extension by a coloop will never occur.
        """
        matroid_list = []
        G = self.graph()
        if element is None:
            element = newlabel(self.groundset())
        else:
            if element in self.groundset():
                raise ValueError("cannot extend by element already in groundset")
        if vertices is None:
            vertices = self._G.vertices()

        # First extend by a loop, then consider every pair of vertices.
        # We'll make the loop a new component.
        v = G.add_vertex()
        G.add_edge(v,v,element)
        matroid_list.append(GraphicMatroid(deepcopy(G)))
        G.delete_vertex(v)

        pairs = combinations(vertices, 2)
        for p in pairs:
            G.add_edge(p[0], p[1], element)
            matroid_list.append(GraphicMatroid(deepcopy(G)))
            G.delete_edge(p[0], p[1], element)
        return iter(matroid_list)

    def graphic_coextension(self, u, v=None, X=None, element=None):
        """
        Returns a matroid coextended by a new element.

        INPUT:

        - ``u`` -- The vertex to be split. If u is not a vertex of the
          matroid's graph, then the new element will be a coloop.
        - ``v`` -- (optional) The name of the new vertex resulting
          from the splitting.
        - ``X`` -- (optional) A list of the matroid elements
          corresponding to
          edges of ``u`` that move to ``v`` after splitting.
        - ``element`` -- (optional) The name of the newly added element.

        OUTPUT:

        An instance of GraphicMatroid extended by the new element.

        .. NOTE::

            A loop on ``u`` will stay a loop unless it is in ``X``.
        """
        if element is None:
            element = newlabel(self.groundset())
        else:
            if element in self.groundset():
                raise ValueError("cannot extend by element already in groundset")
        # To prevent an error for iterating over None:
        if X is None:
            X = []

        G = self.graph()
        vertices = G.vertices()
        if u not in vertices:
            vertices.append(u)
            G.add_vertex(u)
        edgelist = self._groundset_to_edges(X)
        edges_on_u = G.edges_incident(u)
        if u == v or v in vertices:
            raise ValueError("v must be a distinct vertex")
        elif v is None:
            v = G.add_vertex()
        for e in edgelist:
            if e not in edges_on_u:
                # if e is a loop, put it on u and v
                # otherwise raise an error
                if e[0] == e[1]:
                    G.add_edge(u, v, e[2])
                    G.delete_edge(e)
                else:
                    raise ValueError("the edges are not all incident with u")

            elif e[0] == u:
                G.add_edge(v, e[1], e[2])
            elif e[1] == u:
                G.add_edge(e[0], v, e[2])
            G.delete_edge(e)
        G.add_edge(u, v, element)

        return GraphicMatroid(deepcopy(G))

    def twist(self, X):
        """
        Performs a Whitney twist on the graph if `X` is part of a
        2-separation.

        INPUT:

        - ``X`` - The set of elements to be twisted with respect
          to the rest of the matroid. The connectivity of ``X`` must be 1,
          and deletion of the edges corresponding to the elements of ``X``
          must not create new nontrivial components.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: edgelist = [(0,1,0), (1,2,1), (1,2,2), (2,3,3), (2,3,4), (2,3,5), (3,0,6)]
            sage: M = GraphicMatroid(Graph(edgelist, multiedges=True))
            sage: M1 = M.twist([0,1,2]); M1.graph().edges()
            [(0, 1, 1), (0, 1, 2), (0, 3, 6), (1, 2, 0), (2, 3, 3), (2, 3, 4), (2, 3, 5)]
            sage: M2 = M.twist([0,1,3])
            Traceback (most recent call last):
            ...
            ValueError: the input must display a 2-separation and not a 1-separation


        """
        # We require two things:
        # The connectivity of X is less than 2,
        # and the graph M is connected if we delete X
        if not set(X).issubset(self.groundset()):
            raise ValueError("X must be a subset of the ground set")
        connectivity = self.connectivity(X)
        if connectivity != 1:
            raise ValueError("the input must display a 2-separation "
                + "and not a 1-separation")
        G = self.graph()
        X_edges = self.groundset_to_edges(X)
        G.delete_edges(X_edges)
        isolated_vertices = [v for v in G if G.degree(v) == 0]
        G.delete_vertices(isolated_vertices)
        if not G.is_connected():
            raise ValueError("the input must be a separation of the graph")

        # Determine the vertices
        X_vertices = set([e[0] for e in X_edges]).union(
            set([e[1] for e in X_edges]))
        Y_edges = self.groundset_to_edges(self.groundset().difference(set(X)))
        Y_vertices = set([e[0] for e in Y_edges]).union(
            set([e[1] for e in Y_edges]))
        vertices = X_vertices.intersection(Y_vertices)
        a = list(vertices)[0]
        b = list(vertices)[1]

        edges = [(u, v, l) for (u, v, l) in X_edges if (
            u in vertices or v in vertices)]
        G = self.graph()
        for (u, v, l) in edges:
            G.delete_edge(u, v, l)
            if u == a:
                u = b
            elif u == b:
                u = a
            if v == a:
                v = b
            elif v == b:
                v = a
            G.add_edge(u, v, l)
        return GraphicMatroid(G)


    def one_sum(self, X, u, v):
        """
        Similar to above, except rearranging blocks
        """
        pass
