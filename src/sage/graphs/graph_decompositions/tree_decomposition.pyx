# cython: binding=True
r"""
Tree decompositions

This module implements tree-decomposition methods.

A tree-decomposition of a graph `G = (V, E)` is a pair `(X, T)`, where `X=\{X_1,
X_2, \ldots, X_t\}` is a familly of subsets of `V`, usually called *bags*, and
`T` is a tree of order `t` whose nodes are the subsets `X_i` satisfying the
following properties:

1. The union of all sets `X_i` equals `V`. That is, each vertex of the graph `G`
  is associated with at least one tree node.
2. For every edge `(v, w)` in the graph, there is a subset `X_i` that contains
  both `v` and `w`. That is, each edge of the graph `G` appears in a tree node.
3. The nodes associated with vertex `v \in V` form a connected subtree of
  `T`. That is, if `X_i` and `X_j` both contain a vertex `v \in V`, then all
  nodes `X_k` of the tree in the (unique) path between `X_i` and `X_j` contain
  `v` as well, and we have `X_i \cap X_j \subseteq X_k`.

The *width* of a tree decomposition is the size of the largest set `X_i` minus
one, i.e., `\max_{X_i \in X} |X_i| - 1`, and the *treewidth* `tw(G)` of a graph
`G` is the minimum width among all possible tree decompositions of `G`. Observe
that, the size of the largest set is diminished by one in order to make the
treewidth of a tree equal to one.

.. SEEALSO::

    - :wikipedia:`Tree_decomposition`
    - :wikipedia:`Treewidth`


**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`treewidth` | Compute the tree-width of `G` (and provide a decomposition).
    :meth:`is_valid_tree_decomposition` | Check whether `T` is a valid tree-decomposition for `G`.
    :meth:`reduced_tree_decomposition(T)` | Return a reduced tree-decomposition of `T`.


.. TODO:

    - Ensure that we use decomposition by clique separators as preprocessing for
      treewidth
    - Add method to return a *nice* tree decomposition
    - Add methods to compute treelength and examples in the module documentation
      of the difference between treewidth and treelength
    - Approximation of treelength based on
      :meth:`~sage.graphs.graph.Graph.lex_M`
    - Approximation of treelength based on BFS Layering
    - upgrade tdlib to 0.9.0 :trac:`30813` 


Methods
-------
"""
# ****************************************************************************
#       Copyright (C) 2020 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.sets.set import Set
from sage.misc.cachefunc import cached_function
from itertools import combinations
from itertools import chain
from sage.features import PythonModule
from sage.sets.disjoint_set import DisjointSet

def is_valid_tree_decomposition(G, T):
    r"""
    Check whether `T` is a valid tree-decomposition for `G`.

    INPUT:

    - ``G`` -- a sage Graph

    - ``T`` -- a tree decomposition, i.e., a tree whose vertices are the bags
      (subsets of vertices) of the decomposition

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import is_valid_tree_decomposition
        sage: K = graphs.CompleteGraph(4)
        sage: T = Graph()
        sage: T.add_vertex(Set(K))
        sage: is_valid_tree_decomposition(K, T)
        True

        sage: G = graphs.RandomGNP(10, .2)
        sage: T = G.treewidth(certificate=True)
        sage: is_valid_tree_decomposition(G, T)
        True

    The union of the bags is the set of vertices of `G`::

        sage: G = graphs.PathGraph(4)
        sage: T = G.treewidth(certificate=True)
        sage: _ = G.add_vertex()
        sage: is_valid_tree_decomposition(G, T)
        False

    Each edge of `G` is contained in a bag::

        sage: G = graphs.PathGraph(4)
        sage: T = G.treewidth(certificate=True)
        sage: G.add_edge(0, 3)
        sage: is_valid_tree_decomposition(G, T)
        False

    The bags containing a vertex `v` form a subtree of `T`::

        sage: G = graphs.PathGraph(4)
        sage: X1, X2, X3 = Set([0, 1]), Set([1, 2]), Set([2, 3])
        sage: T = Graph([(X1, X3), (X3, X2)])
        sage: is_valid_tree_decomposition(G, T)
        False

    TESTS:

    Check that both parameters are sage graphs::

        sage: is_valid_tree_decomposition("foo", Graph())
        Traceback (most recent call last):
        ...
        ValueError: the first parameter must be a sage Graph
        sage: is_valid_tree_decomposition(Graph(), "foo")
        Traceback (most recent call last):
        ...
        ValueError: the second parameter must be a sage Graph

    Check that `T` is a tree::

        sage: is_valid_tree_decomposition(Graph(), Graph(2))
        Traceback (most recent call last):
        ...
        ValueError: the second parameter must be a tree

    The vertices of `T` must be iterables::

        sage: is_valid_tree_decomposition(Graph(1), Graph(1))
        Traceback (most recent call last):
        ...
        ValueError: the vertices of T must be iterables

    Small cases::

        sage: is_valid_tree_decomposition(Graph(), Graph())
        True
        sage: is_valid_tree_decomposition(Graph(1), Graph())
        False
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the first parameter must be a sage Graph")
    if not isinstance(T, Graph):
        raise ValueError("the second parameter must be a sage Graph")
    if not T:
        return not G
    elif not T.is_tree():
        raise ValueError("the second parameter must be a tree")

    for X in T:
        try:
            _ = list(X)
        except:
            raise ValueError("the vertices of T must be iterables")

    # 1. The union of the bags equals V
    if set(G) != set(chain(*T)):
        return False

    # 2. Each edge of G is contained in a bag
    vertex_to_bags = {u: set() for u in G}
    for Xi in T:
        for u in Xi:
            vertex_to_bags[u].add(Xi)

    for u, v in G.edge_iterator(labels=False):
        if all(v not in Xi for Xi in vertex_to_bags[u]):
            return False

    # 3. The bags containing a vertex v form a connected subset of T
    for X in vertex_to_bags.values():
        D = DisjointSet(X)
        for Xi in X:
            for Xj in T.neighbor_iterator(Xi):
                if Xj in X:
                    D.union(Xi, Xj)
        if D.number_of_subsets() > 1:
            return False

    return True

def reduced_tree_decomposition(T):
    r"""
    Return a reduced tree-decomposition of `T`.

    We merge all edges between two sets `S` and `S'` where `S` is a subset of
    `S'`. To do so, we use a simple union-find data structure to record merge
    operations and the good sets.

    .. WARNING::

        This method assumes that the vertices of the input tree `T` are hashable
        and have attribute ``issuperset``, e.g., ``frozenset`` or
        :class:`~sage.sets.set.Set_object_enumerated_with_category`.

    INPUT:

    - ``T`` -- a tree-decomposition

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import reduced_tree_decomposition
        sage: from sage.graphs.graph_decompositions.tree_decomposition import is_valid_tree_decomposition
        sage: G = graphs.PathGraph(3)
        sage: T = Graph()
        sage: T.add_path([Set([0]), Set([0, 1]), Set([1]), Set([1, 2]),  Set([2])])
        sage: T.order()
        5
        sage: is_valid_tree_decomposition(G, T)
        True
        sage: T2 = reduced_tree_decomposition(T)
        sage: is_valid_tree_decomposition(G, T2)
        True
        sage: T2.order()
        2

    TESTS::

        sage: G = graphs.PathGraph(3)
        sage: T = G.treewidth(certificate=True)
        sage: is_valid_tree_decomposition(G, T)
        True
        sage: T == reduced_tree_decomposition(T)
        True
    """
    def get_ancestor(ancestor, u):
        if ancestor[u] == u:
            return u
        ancestor[u] = get_ancestor(ancestor, ancestor[u])
        return ancestor[u]

    ancestor = {u: u for u in T}
    for u, v in T.edge_iterator(labels=False):
        u = get_ancestor(ancestor, u)
        v = get_ancestor(ancestor, v)
        if u == v:
            continue
        elif u.issuperset(v):
            ancestor[v] = u
        elif v.issuperset(u):
            ancestor[u] = v

    from sage.graphs.graph import Graph
    H = Graph(multiedges=False, name="Reduced tree-decomposition of {}".format(T.name()))
    for u, v in T.edge_iterator(labels=False):
        u = get_ancestor(ancestor, u)
        v = get_ancestor(ancestor, v)
        if u != v:
            H.add_edge(u, v)
    return H


def treewidth(self, k=None, certificate=False, algorithm=None):
    r"""
    Computes the tree-width of `G` (and provides a decomposition)

    INPUT:

    - ``k`` -- integer (default: ``None``); indicates the width to be
      considered. When ``k`` is an integer, the method checks that the graph has
      treewidth `\leq k`. If ``k`` is ``None`` (default), the method computes
      the optimal tree-width.

    - ``certificate`` -- boolean (default: ``False``); whether to return the
      tree-decomposition itself.

    - ``algorithm`` -- whether to use ``"sage"`` or ``"tdlib"`` (requires the
      installation of the 'tdlib' package). The default behaviour is to use
      'tdlib' if it is available, and Sage's own algorithm when it is not.

    OUTPUT:

        ``g.treewidth()`` returns the treewidth of ``g``. When ``k`` is
        specified, it returns ``False`` when no tree-decomposition of width
        `\leq k` exists or ``True`` otherwise. When ``certificate=True``, the
        tree-decomposition is also returned.

    ALGORITHM:

        This function virtually explores the graph of all pairs
        ``(vertex_cut,cc)``, where ``vertex_cut`` is a vertex cut of the graph
        of cardinality `\leq k+1`, and ``connected_component`` is a connected
        component of the graph induced by ``G-vertex_cut``.

        We deduce that the pair ``(vertex_cut,cc)`` is feasible with tree-width
        `k` if ``cc`` is empty, or if a vertex ``v`` from ``vertex_cut`` can be
        replaced with a vertex from ``cc``, such that the pair
        ``(vertex_cut+v,cc-v)`` is feasible.

    .. NOTE::

        The implementation would be much faster if ``cc``, the argument of the
        recursive function, was a bitset. It would also be very nice to not copy
        the graph in order to compute connected components, for this is really a
        waste of time.

    .. SEEALSO::

        :meth:`~sage.graphs.graph_decompositions.vertex_separation.path_decomposition`
        computes the pathwidth of a graph. See also the
        :mod:`~sage.graphs.graph_decompositions.vertex_separation` module.

    EXAMPLES:

    The PetersenGraph has treewidth 4::

        sage: graphs.PetersenGraph().treewidth()
        4
        sage: graphs.PetersenGraph().treewidth(certificate=True)
        Tree decomposition: Graph on 6 vertices

    The treewidth of a 2d grid is its smallest side::

        sage: graphs.Grid2dGraph(2,5).treewidth()
        2
        sage: graphs.Grid2dGraph(3,5).treewidth()
        3

    TESTS::

        sage: g = graphs.PathGraph(3)
        sage: g.treewidth()
        1
        sage: g = 2*graphs.PathGraph(3)
        sage: g.treewidth()
        1
        sage: g.treewidth(certificate=True)
        Tree decomposition: Graph on 4 vertices
        sage: g.treewidth(2)
        True
        sage: g.treewidth(1)
        True
        sage: Graph(1).treewidth()
        0
        sage: Graph(0).treewidth()
        -1
        sage: g = graphs.PetersenGraph()
        sage: g.treewidth(k=2)
        False
        sage: g.treewidth(k=6)
        True
        sage: g.treewidth(certificate=True).is_tree()
        True
        sage: g.treewidth(k=3, certificate=True)
        False
        sage: T = g.treewidth(k=4,certificate=True)
        sage: T
        Tree decomposition: Graph on 6 vertices
        sage: from sage.graphs.graph_decompositions.tree_decomposition import is_valid_tree_decomposition
        sage: is_valid_tree_decomposition(g, T)
        True

    All edges do appear (:trac:`17893`)::

        sage: from itertools import combinations
        sage: g = graphs.PathGraph(10)
        sage: td = g.treewidth(certificate=True)
        sage: for bag in td:
        ....:    g.delete_edges(list(combinations(bag,2)))
        sage: g.size()
        0

    :trac:`19358`::

        sage: g = Graph()
        sage: for i in range(3):
        ....:     for j in range(2):
        ....:         g.add_path([i,(i,j),(i+1)%3])
        sage: g.treewidth()
        2

    The decomposition is a tree (:trac:`23546`)::

        sage: g = Graph({0:[1,2], 3:[4,5]})
        sage: t = g.treewidth(certificate=True)
        sage: t.is_tree()
        True
        sage: vertices = set()
        sage: for s in t.vertices():
        ....:     vertices = vertices.union(s)
        sage: vertices == set(g.vertices())
        True

    Check that the use of atoms and clique separators is correct
    (:trac:`30993`)::

        sage: g = 2 * graphs.Grid2dGraph(2, 3)
        sage: g.treewidth(algorithm='sage')
        2

    Trivially true::

        sage: graphs.PetersenGraph().treewidth(k=35)
        True
        sage: graphs.PetersenGraph().treewidth(k=35,certificate=True)
        Tree decomposition: Graph on 1 vertex

    Bad input::

        sage: graphs.PetersenGraph().treewidth(k=-3)
        Traceback (most recent call last):
        ...
        ValueError: k(=-3) must be a nonnegative integer
    """
    g = self

    # Check Input
    if algorithm is None or algorithm == 'tdlib':
        try:
            import sage.graphs.graph_decompositions.tdlib as tdlib
            tdlib_found = True
        except ImportError:
            tdlib_found = False

    elif algorithm != "sage":
        raise ValueError("'algorithm' must be equal to 'tdlib', 'sage', or None")

    if algorithm is None:
        if tdlib_found:
            algorithm = 'tdlib'
        else:
            algorithm = 'sage'

    if k is not None and k < 0:
        raise ValueError("k(={}) must be a nonnegative integer".format(k))

    # Stupid cases
    from sage.graphs.graph import Graph
    if not g.order():
        if certificate:
            return Graph()
        elif k is None:
            return -1
        else:
            return True

    if k is not None and k >= g.order() - 1:
        if certificate:
            return Graph({Set(g): []}, name="Tree decomposition")
        return True

    # TDLIB
    if algorithm == 'tdlib':
        if not tdlib_found:
            from sage.features import FeatureNotPresentError
            raise FeatureNotPresentError(PythonModule('sage.graphs.graph_decompositions.tdlib',
                                                      spkg='tdlib'))

        T = tdlib.treedecomposition_exact(g, -1 if k is None else k)
        width = tdlib.get_width(T)

        if certificate:
            return T if (k is None or width <= k) else False
        elif k is None:
            return width
        else:
            return width <= k

    # The treewidth of a graph is the maximum over its atoms. So, we decompose
    # the graph by clique minimal separators, compute the treewidth of each of
    # its atoms, and combine the results.
    # This decomposition also deals with disconnected cases.
    atoms, cliques = g.atoms_and_clique_separators()
    if cliques:
        if not certificate:
            if k is None:
                return max(g.subgraph(a).treewidth() for a in atoms)
            elif max(len(c) for c in cliques) - 1 > k:
                return False
            else:
                return all(g.subgraph(a).treewidth(k) for a in atoms)
        else:
            # We compute the tree decomposition of each atom
            T = [g.subgraph(a).treewidth(certificate=True) for a in atoms]
            # and merge the resulting trees
            while len(T) > 1:
                A = T.pop()
                B = T.pop()
                C = cliques.pop()
                if not C:
                    # g is not connected and so some cliques are empty
                    # We connect the trees by an edge
                    A.add_edge(next(A.vertex_iterator()), next(B.vertex_iterator()))
                else:
                    # We search for a vertex in A and B containing clique C
                    ua, ub = None, None
                    for u in A:
                        if u.issuperset(C):
                            ua = u
                    for u in B:
                        if u.issuperset(C):
                            ub = u
                    if ua and ub:
                        A.add_edge(ua, ub)
                    else:
                        # This should never happen
                        raise RuntimeError("something goes wrong."
                                           "Please report the issue to"
                                           "sage-devel@googlegroups.com")

                # We merge A and B
                A.add_edges(B.edges())
                T.append(A)

            return T[0]

    # Forcing k to be defined
    if k is None:
        for i in range(max(0, g.clique_number() - 1, min(g.degree())), g.order()+1):
            ans = g.treewidth(k=i, certificate=certificate)
            if ans:
                return ans if certificate else i

    # This is the recursion described in the method's documentation. All
    # computations are cached, and depends on the pair ``cut,
    # connected_component`` only.
    #
    # It returns either a boolean or the corresponding tree-decomposition, as a
    # list of edges between vertex cuts (as it is done for the complete
    # tree-decomposition at the end of the main function.
    @cached_function
    def rec(cut, cc):
        # Easy cases
        if len(cut) > k:
            return False
        if len(cc) + len(cut) <= k + 1:
            return [(cut, cut.union(cc))] if certificate else True

        # We explore all possible extensions of the cut
        for v in cc:

            # New cuts and connected components, with v respectively added
            # and removed
            cutv = cut.union([v])
            ccv = cc.difference([v])

            # The values returned by the recursive calls.
            sons = []

            # Removing v may have disconnected cc. We iterate on its
            # connected components
            for cci in g.subgraph(ccv).connected_components(sort=False):

                # The recursive subcalls. We remove on-the-fly the vertices
                # from the cut which play no role in separating the
                # connected component from the rest of the graph.
                reduced_cut = frozenset([x for x in cutv
                                         if any(xx in cci for xx in g.neighbor_iterator(x))])
                son = rec(reduced_cut, frozenset(cci))
                if not son:
                    break

                if certificate:
                    sons.extend(son)
                    sons.append((cut, cutv))
                    sons.append((cutv, reduced_cut))

            # Weird Python syntax which is useful once in a lifetime : if break
            # was never called in the loop above, we return "sons".
            else:
                return sons if certificate else True

        return False

    # Main call to rec function, i.e. rec({v}, V-{v})
    V = list(g)
    v = frozenset([V.pop()])
    TD = rec(v, frozenset(V))

    if TD is False:
        return False

    if not certificate:
        return True

    # Building the Tree-Decomposition graph. Its vertices are cuts of the
    # decomposition, and there is an edge from a cut C1 to a cut C2 if C2 is an
    # immediate subcall of C1
    G = Graph()
    G.add_edges(((Set(x), Set(y)) for x,y in TD), loops=False)

    # The Tree-Decomposition may contain a lot of useless nodes.
    # We merge all edges between two sets S,S' where S is a subset of S'
    G = reduced_tree_decomposition(G)

    G.name("Tree decomposition")
    return G
