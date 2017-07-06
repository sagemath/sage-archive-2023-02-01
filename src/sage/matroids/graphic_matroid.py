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
from .utilities import newlabel, split_vertex
from itertools import combinations
from sage.rings.integer import Integer

class GraphicMatroid(Matroid):
    """
    The graphic matroid class.

    INPUT:
    - ``G`` -- a SageMath graph
    - ``groundset`` -- (optional) a list in 1-1 correspondence with
      ``G.edges()``

    EXAMPLES::

        sage: from sage.matroids.advanced import *
        sage: M = GraphicMatroid(graphs.BullGraph()); M
        Graphic matroid of rank 4 on 5 elements.
        sage: N = GraphicMatroid(graphs.CompleteBipartiteGraph(3,3)); N
        Graphic matroid of rank 5 on 9 elements.
    """

    def __init__(self, G, groundset = None):
        """
        See class definition for full documentation.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
            sage: G1 = graphs.CycleGraph(3); G2 = graphs.DiamondGraph()
            sage: G = G1.disjoint_union(G2)
            sage: M = GraphicMatroid(G)
            sage: M
            Graphic matroid of rank 5 on 8 elements.
            sage: M.graph()
            Looped multi-graph on 6 vertices
            sage: M.graph().is_connected()
            True
            sage: M.is_connected()
            False

        """

        if groundset is None:
            #Try to construct a ground set based on the edge labels.
            #If that fails, use range() to come up with a groundset.
            groundset = G.edge_labels()

        groundset_set = frozenset(groundset)

        # if the provided ground set is incomplete, it gets overwriten
        # invalidate `None` as label
        if len(groundset_set) != G.num_edges() or None in groundset_set:
            groundset = range(G.num_edges())
            groundset_set = frozenset(groundset)

        self._groundset = groundset_set

        # Map vertices on input graph to vertices in self._G
        self._vertex_map = {v: v for v in G.vertices()}
        comps = G.connected_components()
        while len(comps) > 1:
            comp = comps.pop()
            v1 = comps[-1][-1]
            v2 = comp[0]
            self._vertex_map[v2] = v1
            comps[-1].extend(comp)

        # Construct a graph and assign edge labels corresponding to the ground set
        edge_list = []
        for i, e in enumerate(G.edges()):
            edge_list.append((self._vertex_map[e[0]],
                self._vertex_map[e[1]], groundset[i]))
        self._G = Graph(edge_list, loops=True, multiedges=True, weighted=True,
            data_structure='static_sparse')
        # Map ground set elements to graph edges:
        # The the edge labels should already be the elements.
        self._groundset_edge_map = ({l: (u, v) for
            (u, v, l) in self._G.edges()})

    #COPYING, LOADING, SAVING

    def __copy__(self):
        """
        Create a shallow copy.

        EXAMPLES::

            sage: M = Matroid(graphs.RandomGNP(5,.5))
            sage: N = copy(M)
            sage: M == N
            True
        """
        N = GraphicMatroid(self._G)
        if getattr(self, '__custom_name') is not None:  # because of name wrangling, this is not caught by the default copy
            N.rename(getattr(self, '__custom_name'))
        return N


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
            sage: M.rank(M.groundset())
            3
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

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph())
            sage: M.groundset()
            frozenset({0, 1, 2, 3, 4})
            sage: G = graphs.CompleteGraph(3).disjoint_union(graphs.CompleteGraph(4))
            sage: M = Matroid(G); M.groundset()
            frozenset({0, 1, 2, 3, 4, 5, 6, 7, 8})
            sage: M = Matroid(Graph([(0, 1, 'a'), (0, 2, 'b'), (0, 3, 'c')]))
            sage: M.groundset()
            frozenset({'a', 'b', 'c'})
        """
        return self._groundset

    def graph(self):
        """
        Return a graph that has a cycle matroid equal to the matroid.

        The graph will always have loops and multiedges enabled.

        EXAMPLES::

            sage: M = Matroid(Graph([(0, 1, 'a'), (0, 2, 'b'), (0, 3, 'c')]))
            sage: M.groundset()
            frozenset({'a', 'b', 'c'})
            sage: M.graph().edges()
            [(0, 1, 'a'), (0, 2, 'b'), (0, 3, 'c')]
            sage: M = Matroid(graphs.CompleteGraph(5))
            sage: M.graph()
            Looped multi-graph on 5 vertices
        """
        # Return a mutable graph
        return self._G.copy(data_structure='sparse')

    def vertex_map(self):
        """
        Return a dictionary mapping the input vertices to the current vertices.

        The graph for the matroid is alway connected. If the constructor is
        given a graph with multiple components, it will connect them. The
        Python dictionary given by this method has the vertices from the
        input graph as keys, and the corresponding vertex label after any
        merging as values.

        EXAMPLES::

            sage: G = Graph([(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5),
            ....: (6, 7), (6, 8), (7, 8), (8, 8), (7, 8)], multiedges=True, loops=True)
            sage: M = Matroid(G)
            sage: M.graph().edges()
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (2, 4, 3),
             (2, 5, 4),
             (4, 5, 5),
             (5, 7, 6),
             (5, 8, 7),
             (7, 8, 8),
             (7, 8, 9),
             (8, 8, 10)]
            sage: M.vertex_map()
            {0: 0, 1: 1, 2: 2, 3: 2, 4: 4, 5: 5, 6: 5, 7: 7, 8: 8}
        """
        return copy(self._vertex_map)

    def _repr_(self):
        """
        Returns a string representation of the matroid.

        EXAMPLES::

            sage: M = Matroid(graphs.CompleteGraph(5))
            sage: M
            Graphic matroid of rank 4 on 10 elements
            sage: G = Graph([(0, 0), (0, 1), (0, 2), (1, 1), (2, 2)], loops=True)
            sage: M = Matroid(G)
            sage: M
            Graphic matroid of rank 2 on 5 elements
        """
        self._mrank = str(self._rank(self._groundset))
        self._elts = str(len(self._groundset))

        return "Graphic matroid of rank " + self._mrank + " on " + self._elts + " elements"

    def _minor(self, contractions=frozenset([]), deletions=frozenset([])):
        """
        Return a minor.

        INPUT:

        - ``contractions`` -- frozenset, subset of ``self.groundset()`` to be contracted
        -  ``deletions`` -- frozenset, subset of ``self.groundset()`` to be deleted

        Assumptions: contractions are independent, deletions are coindependent,
        contractions and deletions are disjoint.

        EXAMPLES::

            sage: M = Matroid(graphs.CompleteGraph(5))
            sage: M._minor(deletions=frozenset([0,1,2]))
            Graphic matroid of rank 4 on 7 elements.
            sage: M._minor(contractions=frozenset([0,1,2]))
            Graphic matroid of rank 1 on 7 elements.
            sage: M = Matroid(graphs.PetersenGraph()); M.groundset()
            frozenset({0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14})
            sage: N = M._minor(deletions = frozenset([0, 3, 5, 9]), contractions =
            ....: frozenset([1, 2, 11])); N
            Graphic matroid of rank 6 on 8 elements.
        """
        g = self.graph()
        cont_edges = self._groundset_to_edges(contractions)
        del_edges = self._groundset_to_edges(deletions)
        # deletions first so contractions don't mess up the vertices
        g.delete_edges(del_edges)
        g.contract_edges(cont_edges)

        return GraphicMatroid(g)

    def _has_minor(self, N, certificate = False):
        """
        Check if the matroid has a minor isomoprhic to M(H).

        INPUT:

        - ``N`` - A matroid.
        - ``certificate`` - (default: ``False``)
        """
        # This method needs work...

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
        Return a list of edges corresponding to a set of ground set elements.

        INPUT:

        - ``X`` -- a subset of the ground set.

        OUTPUT:

        A list of graph edges.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph()); M.groundset()
            frozenset({0, 1, 2, 3, 4})
            sage: M.groundset_to_edges([2,3,4])
            [(1, 2, 2), (1, 3, 3), (2, 3, 4)]
            sage: M.groundset_to_edges([2,3,4,5])
            Traceback (most recent call last):
            ...
            ValueError: input must be a subset of the ground set
        """
        for x in X:
            if x not in self._groundset:
                raise ValueError("input must be a subset of the ground set")
        return self._groundset_to_edges(X)

    def _groundset_to_edges(self, X):
        """
        Return a list of edges corresponding to a set of ground set elements.

        INPUT:

        - ``X`` -- a subset of the ground set.

        OUTPUT:

        A list of graph edges.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph()); M.groundset()
            frozenset({0, 1, 2, 3, 4})
            sage: M._groundset_to_edges([2,3,4])
            [(1, 2, 2), (1, 3, 3), (2, 3, 4)]
        """
        edge_list = []
        for x in X:
            v0 = self._groundset_edge_map[x][0]
            v1 = self._groundset_edge_map[x][1]
            edge_list.append((v0, v1, x))
        return edge_list

    def subgraph_from_set(self, X):
        """
        Return the subgraph corresponding to M restricted to X.

        INPUT:

        - ``X`` -- a subset of the ground set.

        OUTPUT:

        A SageMath graph.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph()); M.groundset()
            frozenset({0, 1, 2, 3, 4})
            sage: M.subgraph_from_set([0,1,2])
            Looped multi-graph on 3 vertices
            sage: M.subgraph_from_set([3,4,5])
            Traceback (most recent call last):
            ...
            ValueError: input must be a subset of the ground set
        """
        for x in X:
            if x not in self._groundset:
                raise ValueError("input must be a subset of the ground set")
        return self._subgraph_from_set(X)

    def _subgraph_from_set(self, X):
        """
        Return the subgraph corresponding to M restricted to X.

        INPUT:

        - ``X`` -- a subset of the ground set.

        OUTPUT:

        A SageMath graph.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph()); M.groundset()
            frozenset({0, 1, 2, 3, 4})
            sage: M._subgraph_from_set([0,1,2])
            Looped multi-graph on 3 vertices
        """
        edge_list = self._groundset_to_edges(X)
        return Graph(edge_list, loops=True, multiedges=True)

    def _corank(self, X):
        """
        Return the corank of the set `X` in the matroid.

        Internal version that does no input checking.

        INPUT:

        - ``X`` -- An iterable container of ground set elements.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: M = Matroid(graphs.CompleteBipartiteGraph(3,3))
            sage: M._corank([0,1,2])
            2
            sage: M._corank([1,2,3])
            3
        """
        g = self.graph()
        g.delete_edges(self._groundset_to_edges(X))
        return (len(X) - (g.connected_components_number() - Integer(1)))

    def _is_independent(self, X):
        """
        Test if a set is an independent set of the matroid.

        INPUT:

        - ``X`` -- An iterable container of ground set elements.

        OUTPUT:

        Boolean.

        EXAMPLES::

        sage: M = Matroid(graphs.DiamondGraph())
        sage: M._is_independent([0,1,2])
        False
        sage: M._is_independent([0,1,3])
        True
        """
        g = self._subgraph_from_set(X)
        return g.is_forest()

    def _is_circuit(self, X):
        """
        Test if input is a circuit.

        INPUT:

        - ``X`` -- An iterable container of ground set elements.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph())
            sage: M._is_circuit([0,1,2])
            True
            sage: M._is_circuit([0,1,2,3])
            False
            sage: M._is_circuit([0,1,3])
            False
        """
        g = self._subgraph_from_set(X)
        return g.is_cycle()

    def _closure(self, X):
        """
        Return the closure of a set.

        INPUT:

        - ``X`` -- An iterable container of ground set elements.

        OUTPUT:

        ``frozenset`` instance containing a subset of the ground set.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph())
            sage: M._closure([0])
            frozenset({0})
            sage: M._closure([0,1])
            frozenset({0, 1, 2})
            sage: M._closure(M.groundset())
            frozenset({0, 1, 2, 3, 4})
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

    def _max_independent(self, X):
        """
        Compute a maximal independent subset.

        INPUT:

        - ``X`` -- An object with Python's ``frozenset`` interface containing
          a subset of ``self.groundset()``.

        OUTPUT:

        ``frozenset`` instance containing a subset of the ground set.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph())
            sage: M._max_independent(M.groundset())
            frozenset({1, 3, 4})
            sage: M._max_independent(frozenset([0,1,2]))
            frozenset({1, 2})
            sage: M._max_independent(frozenset([3,4]))
            frozenset({3, 4})
            sage: N = M.graphic_extension(0, element='a')
            sage: N._max_independent(frozenset(['a']))
            frozenset()
        """
        if self._is_independent(X):
            return X

        edgelist = self._groundset_to_edges(X)
        G = self._subgraph_from_set(X)
        ind_list = []
        for (u, v, l) in edgelist:
            if G.is_cut_edge(u, v, l):
                ind_list.append(l)
            else:
                G.delete_edge(u, v, l)
        return frozenset(ind_list)

    def _max_coindependent(self, X):
        """
        Compute a maximal coindependent subset.

        INPUT:

        - ``X`` -- An iterable container of ground set elements.

        OUTPUT:

        ``frozenset`` instance containing a subset of the ground set.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph())
            sage: M._max_coindependent(M.groundset())
            frozenset({0, 2})
            sage: M._max_coindependent([2,3,4])
            frozenset({2, 3})
            sage: N = M.graphic_extension(0, element='a')
            sage: N.max_coindependent([0,1,2,'a'])
            frozenset({0, 2, 'a'})
        """
        res = set()
        g = self.graph()
        edgelist = self._groundset_to_edges(X)
        for e in edgelist:
            g.delete_edge(e)
            if g.connected_components_number() > 1:
                g.add_edge(e)
            else:
                res.add(e[2])
        return frozenset(res)

    def _circuit(self, X):
        """
        Return a minimal dependent subset.

        INPUT:

        - ``X`` -- An iterable container of ground set elements.

        OUTPUT:

        ``frozenset`` instance containing a subset of the groundset.
        A ``ValueError`` is raised if the set contains no circuit.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph())
            sage: M._circuit(M.groundset())
            frozenset({2, 3, 4})
            sage: N = Matroid(graphs.CompleteBipartiteGraph(3,3))
            sage: N._circuit([0, 1, 2, 6, 7, 8])
            frozenset({1, 2, 7, 8})
            sage: N._circuit([0,1,2])
            Traceback (most recent call last):
            ...
            ValueError: no circuit in independent set
        """
        circuit = set()
        g = self._subgraph_from_set(X)
        if g.is_forest():
            raise ValueError("no circuit in independent set")
        for e in g.edges():
            g.delete_edge(e)
            if g.is_forest():
                g.add_edge(e)
                circuit.add(e[2])
        return frozenset(circuit)

    def _coclosure(self, X):
        """
        Return the coclosure of a set.

        INPUT:

        - ``X`` -- An iterable container of ground set elements.
        OUTPUT:

        ``frozenset`` instance containing a subset of the groundset.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph())
            sage: M._coclosure([0])
            frozenset({0, 1})
            sage: M._coclosure([0,1])
            frozenset({0, 1})
            sage: N = M.graphic_extension(0, element='a')
            sage: N._coclosure([3])
            frozenset({3, 4})
            sage: N = M.graphic_coextension(0, element='a')
            sage: N._coclosure([3])
            frozenset({3, 4, 'a'})
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
        Test if input is coindependent.

        INPUT:

        - ``X`` -- An iterable container of ground set elements.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph())
            sage: M._is_coindependent([0])
            True
            sage: M._is_coindependent([0,1])
            False
            sage: N = M.graphic_coextension(3, element='a')
            sage: N._is_coindependent([0,'a'])
            False
            sage: N1 = N.graphic_extension(3, element='b')
            sage: N1._is_coindependent([0,'b'])
            True
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
        Return a graphic matroid extended by a new element.

        INPUT:

        - ``u`` -- a vertex in the matroid's graph.
        - ``v`` -- (default: ``None``) another vertex. If not specified or if ``v`` is
          ``u``, then, the new element will be a loop.
        - ``element`` -- (default: ``None``) the label of the new element. If
          not specified, a new label will be generated automatically.

        OUTPUT:

        A GraphicMatroid with the specified element added.

        EXAMPLES::

            sage: M = Matroid(graphs.CompleteGraph(4))
            sage: M1 = M.graphic_extension(0,1,'a'); M1
            Graphic matroid of rank 3 on 7 elements.
            sage: M1.graph().edges()
            [(0, 1, 0), (0, 1, 'a'), (0, 2, 1), (0, 3, 2), (1, 2, 3), (1, 3, 4), (2, 3, 5)]
            sage: M2 = M1.graphic_extension(3); M2
            Graphic matroid of rank 3 on 8 elements.

        """
        # __init()__ forces the graph to be connected, so this should
        # never make a coloop
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
        return GraphicMatroid(G)

    def graphic_extensions(self, element=None, vertices=None):
        """
        Return an iterable containing the graphic extensions.

        INPUT:

        - ``element`` -- (optional) The name of the newly added element in
          each extension.
        - ``vertices`` -- (optional) A set of vertices over which the extension
          may be taken. If not given, will use all vertices.

        OUTPUT:

        An iterable containing instances of GraphicMatroid.

        .. NOTE::

            The extension by a loop will always occur.
            The extension by a coloop will never occur.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph())
            sage: I = M.graphic_extensions('a')
            sage: for N in I:
            ....:     N.graph().edges()
            [(0, 0, 'a'), (0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 3, 4)]
            [(0, 1, 0), (0, 1, 'a'), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 3, 4)]
            [(0, 1, 0), (0, 2, 1), (0, 2, 'a'), (1, 2, 2), (1, 3, 3), (2, 3, 4)]
            [(0, 1, 0), (0, 2, 1), (0, 3, 'a'), (1, 2, 2), (1, 3, 3), (2, 3, 4)]
            [(0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 2, 'a'), (1, 3, 3), (2, 3, 4)]
            [(0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 3, 3), (1, 3, 'a'), (2, 3, 4)]
            [(0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 3, 4), (2, 3, 'a')]

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
        # Put the loop on the first vertex.
        G.add_edge(vertices[0], vertices[0], element)
        matroid_list.append(GraphicMatroid(G))
        G.delete_edge(vertices[0], vertices[0], element)

        pairs = combinations(vertices, 2)
        for p in pairs:
            G.add_edge(p[0], p[1], element)
            matroid_list.append(GraphicMatroid(G))
            G.delete_edge(p[0], p[1], element)
        return iter(matroid_list)

    def graphic_coextension(self, u, X=None, element=None):
        """
        Return a matroid coextended by a new element.

        INPUT:

        - ``u`` -- The vertex to be split. If u is not a vertex of the
          matroid's graph, then the new element will be a coloop.
        - ``X`` -- (optional) A list of the matroid elements corresponding to
          edges of ``u`` that move to the new vertex after splitting.
        - ``element`` -- (optional) The name of the newly added element.

        OUTPUT:

        An instance of GraphicMatroid extended by the new element.

        .. NOTE::

            A loop on ``u`` will stay a loop unless it is in ``X``.

        EXAMPLES::

            sage: M = Matroid(graphs.WheelGraph(5))
            sage: M1 = M.graphic_coextension(0, X=[1,2], element='a')
            sage: M1.graph().edges()
            [(0, 1, 0),
             (0, 4, 3),
             (0, 5, 'a'),
             (1, 2, 4),
             (1, 4, 5),
             (2, 3, 6),
             (2, 5, 1),
             (3, 4, 7),
             (3, 5, 2)]

        TESTS::

            sage: M = Matroid(graphs.CycleGraph(3))
            sage: M = M.graphic_extension(0, element='a')
            sage: M.graph().edges()
            [(0, 0, 'a'), (0, 1, 0), (0, 2, 1), (1, 2, 2)]
            sage: M1 = M.graphic_coextension(0, X=[1], element='b')
            sage: M1.graph().edges()
            [(0, 0, 'a'), (0, 1, 0), (0, 3, 'b'), (1, 2, 2), (2, 3, 1)]
            sage: M2 = M.graphic_coextension(0, X=[1, 'a'], element='b')
            sage: M2.graph().edges()
            [(0, 1, 0), (0, 3, 'a'), (0, 3, 'b'), (1, 2, 2), (2, 3, 1)]

        ::

            sage: M = Matroid(graphs.CycleGraph(3))
            sage: M = M.graphic_coextension(u=4, element='a')
            sage: M.graph()
            Looped multi-graph on 4 vertices
            sage: M.graph().loops()
            []
            sage: M = M.graphic_coextension(u=4, element='a')
            Traceback (most recent call last):
            ...
            ValueError: cannot extend by element already in ground set

        """
        if element is None:
            element = newlabel(self.groundset())
        else:
            if element in self.groundset():
                raise ValueError("cannot extend by element already in ground set")
        # To prevent an error for iterating over None:
        if X is None:
            X = []

        G = self.graph()
        vertices = G.vertices()
        if u not in vertices:
            vertices.append(u)
            G.add_vertex(u)
        edgelist = self._groundset_to_edges(X)
        v = G.add_vertex()

        split_vertex(G = G, u = u, v = v, edges=edgelist)
        G.add_edge(u, v, element)

        return GraphicMatroid(G)

    def twist(self, X):
        """
        Performs a Whitney twist on the graph.

        `X` must be part of a 2-separation.
        The connectivity of ``X`` must be 1, and the subgraph induced by X must
        intersect the subgraph induced by the rest of the elements on exactly
        two vertices.

        INPUT:

        - ``X`` -- The set of elements to be twisted with respect
          to the rest of the matroid.

        OUTPUT:

        An instance of ``GraphicMatroid`` isomorphic to this matroid but with
        a graph that is not necessarily isomorphic.

        EXAMPLES::

            sage: edgelist = [(0,1,0), (1,2,1), (1,2,2), (2,3,3), (2,3,4), (2,3,5), (3,0,6)]
            sage: M = Matroid(Graph(edgelist, multiedges=True))
            sage: M1 = M.twist([0,1,2]); M1.graph().edges()
            [(0, 1, 1), (0, 1, 2), (0, 3, 6), (1, 2, 0), (2, 3, 3), (2, 3, 4), (2, 3, 5)]
            sage: M2 = M.twist([0,1,3])
            Traceback (most recent call last):
            ...
            ValueError: the input must display a 2-separation that is not a 1-separation

        TESTS::

            sage: edgedict = {0: [1, 2], 1: [2, 3], 2: [3], 3: [4, 5], 4: [5]}
            sage: M = Matroid(Graph(edgedict))
            sage: M.graph().edges()
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (1, 3, 3),
             (2, 3, 4),
             (3, 4, 5),
             (3, 5, 6),
             (4, 5, 7)]
            sage: M1 = M.twist([0, 1]); M1.graph().edges()
            [(0, 1, 1),
             (0, 2, 0),
             (1, 2, 2),
             (1, 3, 3),
             (2, 3, 4),
             (3, 4, 5),
             (3, 5, 6),
             (4, 5, 7)]
            sage: M2 = M1.twist([0, 1, 2]); M2.graph().edges()
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (1, 3, 3),
             (2, 3, 4),
             (3, 4, 5),
             (3, 5, 6),
             (4, 5, 7)]
            sage: M1 == M
            False
            sage: M2 == M
            True
            sage: M2.twist([3, 4])
            Traceback (most recent call last):
            ...
            ValueError: too many vertices in the intersection
        """
        # We require two things:
        # (1) The connectivity of X is 1,
        # (2) X intersects the rest of the graph on 2 vertices
        if not set(X).issubset(self.groundset()):
            raise ValueError("X must be a subset of the ground set")
        connectivity = self.connectivity(X)
        if connectivity != 1:
            raise ValueError("the input must display a 2-separation "
                + "that is not a 1-separation")
        G = self.graph()

        # Determine the vertices
        X_edges = self.groundset_to_edges(X)
        X_vertices = set([e[0] for e in X_edges]).union(
            set([e[1] for e in X_edges]))
        Y_edges = self.groundset_to_edges(self.groundset().difference(set(X)))
        Y_vertices = set([e[0] for e in Y_edges]).union(
            set([e[1] for e in Y_edges]))
        vertices = X_vertices.intersection(Y_vertices)
        if len(vertices) != 2:
            raise ValueError("too many vertices in the intersection")
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
        Arrange matroid components in the graph.

        The matroid's graph must be connected even if the matroid is not
        connected, but if there are multiple matroid components, the user may
        choose how they are arranged in the graph. This method will take the
        block of the graph that represents `X` and attach it by vertex `u` to
        another vertex `v` in the graph.

        INPUT:

        - ``X`` -- A subset of the ground set.
        - ``u`` -- A vertex spanned by the edges of the elements in ``X``.
        - ``v`` -- A vertex spanned by the edges of the elements not in ``X``.

        EXAMPLES::

            sage: edgedict = {0:[1, 2], 1:[2, 3], 2:[3], 3:[4, 5], 6:[4, 5]}
            sage: M = Matroid(Graph(edgedict))
            sage: M.graph().edges()
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (1, 3, 3),
             (2, 3, 4),
             (3, 4, 5),
             (3, 5, 6),
             (4, 6, 7),
             (5, 6, 8)]
            sage: M1 = M.one_sum(u = 3, v = 1, X = [5, 6, 7, 8])
            sage: M1.graph().edges()
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (1, 3, 3),
             (1, 4, 5),
             (1, 5, 6),
             (2, 3, 4),
             (4, 6, 7),
             (5, 6, 8)]
            sage: M2 = M.one_sum(u = 4, v = 3, X = [5, 6, 7, 8])
            sage: M2.graph().edges()
            [(0, 1, 0),
             (0, 2, 1),
             (1, 2, 2),
             (1, 3, 3),
             (2, 3, 4),
             (3, 6, 7),
             (3, 7, 5),
             (5, 6, 8),
             (5, 7, 6)]


        TESTS::

            sage: M = Matroid(graphs.CompleteGraph(4))
            sage: M.one_sum(u = 1, v = 2, X = [0,1])
            Traceback (most recent call last):
            ...
            ValueError: the input must display a 1-separation

        ::

            sage: M = Matroid(graphs.BullGraph())
            sage: M.graph().edges()
            [(0, 1, 0), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 4, 4)]
            sage: M1 = M.one_sum(u = 3, v = 0, X = [3,4])
            Traceback (most recent call last):
            ...
            ValueError: too many vertices in the intersection

            sage: M1 = M.one_sum(u = 3, v = 2, X = [3])
            sage: M1.graph().edges()
            [(0, 1, 0), (0, 2, 1), (1, 2, 2), (2, 4, 4), (2, 5, 3)]

            sage: M2 = M1.one_sum(u = 5, v = 0, X = [3,4])
            sage: M2.graph().edges()
            [(0, 1, 0), (0, 2, 1), (0, 3, 3), (1, 2, 2), (3, 4, 4)]

            sage: M = Matroid(graphs.BullGraph())
            sage: M.one_sum(u = 0, v = 1, X = [3])
            Traceback (most recent call last):
            ...
            ValueError: first vertex must be spanned by the input

            sage: M.one_sum(u = 1, v = 3, X = [3])
            Traceback (most recent call last):
            ...
            ValueError: second vertex must be spanned by the rest of the graph
        """
        # We require two things:
        # (1) The connectivity of X is 0,
        # (2) X intersects the rest of the graph on 1
        if not set(X).issubset(self.groundset()):
            raise ValueError("X must be a subset of the ground set")
        connectivity = self.connectivity(X)
        if connectivity != 0:
            raise ValueError("the input must display a 1-separation")
        G = self.graph()
        if u not in G or v not in G:
            raise ValueError("the vertices must already be in the graph")

        # Determine the vertex
        X_edges = self.groundset_to_edges(X)
        X_vertices = set([e[0] for e in X_edges]).union(
            set([e[1] for e in X_edges]))
        if u not in X_vertices:
            raise ValueError("first vertex must be spanned by the input")
        Y_edges = self.groundset_to_edges(self.groundset().difference(set(X)))
        Y_vertices = set([e[0] for e in Y_edges]).union(
            set([e[1] for e in Y_edges]))
        if v not in Y_vertices:
            raise ValueError("second vertex must be spanned by " +
                "the rest of the graph")
        vertices = X_vertices.intersection(Y_vertices)
        if len(vertices) != 1:
            raise ValueError("too many vertices in the intersection")
        a = vertices.pop()
        b = G.add_vertex()

        edgeset = set(X_edges).intersection(set(G.edges_incident(a)))
        split_vertex(G, u = a, v = b, edges = edgeset)
        # If u was the cut vertex, u is now detached from our component
        # so we merge the new vertex. Otherwise we can merge u
        if u == a:
            G.merge_vertices([v, b])
        else:
            G.merge_vertices([v, u])

        return GraphicMatroid(G)

    # Comparison:

    def _vertex_stars(self):
        """
        Computes the set of edge labels around each vertex.

        Internal method for hashing purposes.

        INPUT:

        None.

        OUTPUT:

        A ``frozenset`` of ``frozenset``s containing the edge labels around
        each vertex.

        EXAMPLES::

            sage: M = Matroid(graphs.DiamondGraph())
            sage: M._vertex_stars()
            frozenset({frozenset({0, 2, 3}),
                       frozenset({1, 2, 4}),
                       frozenset({3, 4}),
                       frozenset({0, 1})})

            sage: N = Matroid(graphs.BullGraph()); N._vertex_stars()
            frozenset({frozenset({0, 2, 3}),
                       frozenset({4}),
                       frozenset({1, 2, 4}),
                       frozenset({3}),
                       frozenset({0, 1})})

        """
        star_list = []
        for v in self._G.vertices():
            star = [l for (u, v, l) in self._G.edges_incident(v)]
            star_list.append(frozenset(star))
        return frozenset(star_list)

    def __hash__(self):
        r"""
        Return an invariant of the matroid.

        This function is called when matroids are added to a set. It is very
        desirable to override it so it can distinguish matroids on the same
        groundset, which is a very typical use case!

        .. WARNING::

            This method is linked to __richcmp__ (in Cython) and __cmp__ or
            __eq__/__ne__ (in Python). If you override one, you should (and in
            Cython: MUST) override the other!

        EXAMPLES::

            sage: M = Matroid(graphs.CompleteGraph(3))
            sage: N = Matroid(graphs.CycleGraph(3))
            sage: O = Matroid(graphs.ButterflyGraph())
            sage: hash(M) == hash(N)
            True
            sage: hash(O) == hash(N)
            False
            sage: P = Matroid(Graph([(0, 1, 'a'), (0, 2, 'b'), (1, 2, 'c')]))
            sage: hash(P) == hash(M)
            False
        """
        return hash(self._vertex_stars())

    def __eq__(self, other):
        """
        Compare two matroids.

        INPUT:

        - ``other`` -- A matroid.

        OUTPUT:

        ``True`` if ``self`` and ``other`` have the same graph; ``False``
        otherwise.

        EXAMPLES::

            sage: M = Matroid(graphs.CompleteGraph(3))
            sage: N = Matroid(graphs.CycleGraph(3))
            sage: O = Matroid(graphs.ButterflyGraph())
            sage: P = Matroid(Graph([(0, 1, 'a'), (0, 2, 'b'), (1, 2, 'c')]))
            sage: M == N
            True
            sage: M == O
            False
            sage: M == P
            False

        """
        # Graph.__eq__() will ignore edge labels unless we turn on weighted()
        # This will be done in __init__()
        if not isinstance(other, GraphicMatroid):
            return False
        return (self._G == other._G)

    def __ne__(self, other):
        """
        Compare two matroids.

        INPUT:

        - ``other`` -- A matroid.

        OUTPUT:

        ``False`` if ``self`` and ``other`` have the same graph; ``True``
        otherwise.

        EXAMPLES::

            sage: M = Matroid(graphs.CycleGraph(4))
            sage: N = Matroid(graphs.CompleteBipartiteGraph(2,2))
            sage: O = Matroid(graphs.PetersenGraph())
            sage: M != N
            True
            sage: M.equals(N)
            True
            sage: M != O
            True

        """
        return (not self == other)
