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
        sage: graphs.PetersenGraph().treewidth(k=2)
        False
        sage: graphs.PetersenGraph().treewidth(k=6)
        True
        sage: graphs.PetersenGraph().treewidth(certificate=True).is_tree()
        True
        sage: graphs.PetersenGraph().treewidth(k=3,certificate=True)
        False
        sage: graphs.PetersenGraph().treewidth(k=4,certificate=True)
        Tree decomposition: Graph on 6 vertices

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

    Trivially true::

        sage: graphs.PetersenGraph().treewidth(k=35)
        True
        sage: graphs.PetersenGraph().treewidth(k=35,certificate=True)
        Tree decomposition: Graph on 1 vertex

    Bad input:

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

    # Disconnected cases
    if not g.is_connected():
        if not certificate:
            if k is None:
                return max(cc.treewidth() for cc in g.connected_components_subgraphs())
            else:
                return all(cc.treewidth(k) for cc in g.connected_components_subgraphs())
        else:
            T = [cc.treewidth(certificate=True) for cc in g.connected_components_subgraphs()]
            tree = Graph([list(chain(*T)),
                          list(chain(*[t.edges(labels=False, sort=False) for t in T]))],
                         format='vertices_and_edges', name="Tree decomposition")
            v = next(T[0].vertex_iterator())
            for t in T[1:]:
                tree.add_edge(next(t.vertex_iterator()),v)
            return tree

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
    G = Graph(name="Tree decomposition")
    G.add_edges(((Set(x), Set(y)) for x,y in TD), loops=False)

    # The Tree-Decomposition contains a lot of useless nodes.
    #
    # We merge all edges between two sets S,S' where S is a subset of S'
    changed = True
    while changed:
        changed = False
        for v in G.vertices(sort=False):
            for u in G.neighbor_iterator(v):
                if u.issuperset(v):
                    G.merge_vertices([u, v])  # the new vertex is named 'u'
                    changed = True
                    break

    return G
