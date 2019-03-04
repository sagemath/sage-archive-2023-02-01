r"""
Functions related to domination.

This file contains varied functions related to domination as well as an
implementation of the algorithm described in _Enumerating minimal
dominating sets in $K_t$-free graphs and variants_ by Marthe Bonamy,
Oscar Defrain, Marc Heinrich, Michał Pilipczuk, and Jean-Florent Raymond
(https://arxiv.org/abs/1810.00789 ) to enumerate the minimal dominating
sets of a graph.

EXAMPLES::

to be added

AUTHORS:

- Jean-Florent Raymond (2019-03-04): initial version
"""

# ****************************************************************************
#       Copyright (C) 2019 Jean-Florent Raymond  <raymond@tu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.all import *


def closed_neighbor_iterator(self, vertex):
    r'''
    Return an iterator over the closed neighbors of ``vertex``.

    INPUT:
    - ``vertex`` -- the vertex of ``self``, the closed neighborhood of
    which we want to iterate over.

    EXAMPLES::

        sage: g = graphs.CubeGraph(3)
        sage: for i in g.closed_neighbor_iterator('010'):
        ....:     print(i)
        010
        011
        000
        110

    ::

        sage: g = Graph(3, loops = True)
        sage: g.add_edge(0,1)
        sage: g.add_edge(0,0)
        sage: list(g.closed_neighbor_iterator(0))
        [0, 1]
        sage: list(g.closed_neighbor_iterator(2))
        [2]

    TESTS::

        sage: G = graphs.CubeGraph(3)
        sage: list(G.closed_neighbor_iterator('013'))
        LookupError: vertex (013) is not a vertex of the graph
    '''

    if not self.has_vertex(vertex):
        raise LookupError(
            'vertex ({0}) is not a vertex of the graph'.format(vertex))

    if not self.has_edge(vertex, vertex):
        yield vertex

    for neighbor in self.neighbor_iterator(vertex):
        yield neighbor
    return


def neighbors_of_a_set(self, vertex_subset):
    r'''
    Return the set of neighbors of an iterable of vertices.

    INPUT:

    - ``vertex_subset`` -- iterable; contains the vertices of ``self``
    for which we want to enumerate neighbors

    OUTPUT:

    An iterator over vertices of ``self`` that does not belong to
    ``vertex_subset`` and have a neighbor in this set.

    EXAMPLES::

        sage: g = graphs.PathGraph(5)
        sage: g.neighbors_of_a_set({2,3})
        {1, 4}
    '''

    vertex_subset_list = list(vertex_subset)
    if vertex_subset_list:
        neigh = set.union(*(set(self.neighbor_iterator(v))
                            for v in vertex_subset_list))
        neigh.difference_update(vertex_subset_list)
        return neigh
    return set()


def closed_neighbors_of_a_set(self, vertex_subset):
    r'''
    Return the closed neighborhood of an iterable of vertices.

    INPUT:

    - ``vertex_subset`` -- an iterable of vertices of self

    OUTPUT:

    An iterator over vertices of ``self`` that belong to
    ``vertex_subset`` or have a neighbor in this set.

    EXAMPLES::

        sage: G = graphs.PathGraph(5)
        sage: G.neighbors_of_a_set({2,3})
        {1, 2, 3, 4}
    '''

    vertex_subset_list = list(vertex_subset)
    if vertex_subset_list:
        return set.union(*(set(self.closed_neighbor_iterator(v))
                           for v in vertex_subset_list))
    return set()


def private_neighbors(self, vertex, vertex_subset):
    r'''
    Return the private neighbors of a vertex with repect to an iterable.

    INPUT:
    
    - ``vertex`` -- a vertex of ``self``
    - ``vertex_subset`` -- an iterable of vertices of self

    OUTPUT:

    Return the closed neighbors of ``vertex`` that are not closed
    neighbors of an other vertex of ``vertex_subset``.

    EXAMPLES::

        sage: g = graphs.PathGraph(5)
        sage: list(g.private_neighbors(1, [1, 3, 4]))
        [1, 0]

        sage: list(g.private_neighbors(1, [3, 4]))
        [1, 0]

        sage: list(g.private_neighbors(1, [3, 4, 0]))
        []

    '''
    
    closed_neighborhood_vs = set(self.closed_neighbors_of_a_set(
        u for u in vertex_subset if u!=vertex))
    
    
    return (neighbor for neighbor in self.closed_neighbor_iterator(vertex)
            if neighbor not in closed_neighborhood_vs)


def is_dominating(self, p_dominating, p_dominated=None):
    r'''
    Return whether the first set dominates the second one.

    We say that as set `D` of vertices of a graph `G`dominated a set `S`
    if every vertex of `S` belongs to `D` or is adjacent to a vertex of
    `D`.  Also, `D` is a dominating set of `G` if it dominates `V(G)`.

    INPUT:
    
    - ``p_dominating`` -- an iterable of vertices of ``self``
    - ``vertex_subset`` -- (default: `None`) an iterable of vertices of
    ``self``

    OUTPUT:

    Return whether ``p_dominating`` dominates ``p_dominated``.
    If ``p_dominated`` is set to ``None``, returns whether
    ``p_dominating`` dominates ``self``.

    EXAMPLES::

        sage: g = graphs.CycleGraph(5)
        sage: g.is_dominating([0,1], [4, 2])
        True

        sage: g.is_dominating([0,1])
        False

    TESTS::

        sage: g.is_dominating([0,1], {2, 42})
        LookupError: vertex (42) is not a vertex of the graph
    '''
    if p_dominated is None:
        sp_dominated = set(self.vertex_iterator())
    else:
        sp_dominated = set(p_dominated)

        # using a set to check membership repeatedly faster
        all_vertices = set(self.vertex_iterator()) 

        # check that sp_dominated is a subset of self.vertices
        try:
            bad_boy = next(v for v in sp_dominated if v not in all_vertices)
            raise LookupError(
                'vertex ({0}) is not a vertex of the graph'.format(bad_boy))
        except StopIteration:
            pass

    actually_dominated = self.closed_neighbors_of_a_set(p_dominating)
    return sp_dominated <= actually_dominated


def is_redundant(self, dom, focus=None):
    r'''
    Return whether a vertex iterable has redundant vertices.

    Let `G` be a graph and `D` be a subset of its vertices. A vertex `v`
    of `D` is said to be redundant in `S` if every closed neighbors of
    `v` that belongs to `S` is dominated by `D \ {v}`.

    INPUT:
    
    - ``dom`` -- an iterable of vertices of ``self``
    - ``focus`` -- (default: `None`) an iterable of vertices of
    ``self``

    OUTPUT:

    Return whether ``dom`` has a redundant vertex in ``focus``.
    When called with value ``None`` for ``focus``, the function is ran
    with ``focus`` being equal to all vertices of ``self``.

    EXAMPLES::

        sage: G = graphs.CubeGraph(3)
        sage: G.is_redundant(['000', '101'], ['011'])
        True
        sage: G.is_redundant(['000', '101'])
        False
    '''

    if focus is None:
        focus = self.vertices()

    # dominator[v] (for v in focus) will contain the list of vertices
    # of dom that are adjacent to v
    dominator = dict()

    # For every x in dom, has_private[x] will be True
    # iff x has a private neighbor in focus
    has_private = dict()

    # Counts the number of vertices of dom with a private neighbor in S:
    irredundant_ds = 0

    for v in focus:
        dominator[v] = []  # Initialization

    for x in dom:
        has_private[x] = False      # Initialization
        for v in self.closed_neighbor_iterator(x):
            if v in focus:
                # x dominates all its closed neighbors:
                dominator[v].append(x)

    # Now we can compute has_private[]:
    for v in focus:
        # Here we do care neither about vertices dominated by more than
        # one vertex of dom (they are not private neighbor of anybody)
        # nor about vertices not dominated by dom (idem).
        if len(dominator[v]) == 1:
            # If v is dominated only by one vertex x,
            # then v is a private neighbor of x
            if not has_private[dominator[v][0]]:
                # If it was not know yet that x has a private
                has_private[dominator[v][0]] = True
                # We found a new useful vertex
                irredundant_ds = irredundant_ds + 1

    # True iff each elt of dom is irredundant
    return len(dom) != irredundant_ds


# Adding the functions defined above to the Graph class

Graph.closed_neighbor_iterator = closed_neighbor_iterator
Graph.neighbors_of_a_set = neighbors_of_a_set
Graph.closed_neighbors_of_a_set = closed_neighbors_of_a_set
Graph.private_neighbors = private_neighbors

Graph.is_dominating = is_dominating
Graph.is_redundant = is_redundant


# The code below is specific to the enumeration of minimal dominating sets

def _parent(G, dom, V_prev):
    r'''
    Return an subset of dom that is irredundant in V_prev.

    For internal use.

    INPUT:
    
    - ``dom`` -- an iterable of vertices of ``self``
    - ``V_prev`` -- an iterable of vertices of ``self``

    OUTPUT:

    Return the list obtained from ``dom`` by iteratively removing those
    vertices of mininum index that have no private neighbor in V_prev.

    EXAMPLES:

        sage: G = graphs.PathGraph(4)
        sage: G.add_vertices([4, 5])
        sage: G.add_edges([(4, 1), (5, 2)])
        sage: _parent(G, [0, 2, 4, 5], [1, 2])
        [4, 5]
        sage: _parent(G, [0, 2, 4, 5], [1, 3])
        [2]

    .. WARNING:
    
    We assume that vertices are sortable (i.e. they can be compared).
    '''

    # The list where we search vertices
    D_start = sorted(dom, reverse=True)

    # The list to be output at the end, that we construct:
    D_end = []

    while D_start:
        v = D_start.pop() # element of min index
        priv = set(G.closed_neighbor_iterator(v))
        # We remove the vertices already dominated
        # by other vertices of (D_end union D_start)
        priv.difference_update(*(G.closed_neighbor_iterator(u)
                                 for u in D_start if u != v))
        priv.difference_update(*(G.closed_neighbor_iterator(u)
                                 for u in D_end if u != v))
        # Now priv is the private neighborhood of v
        # in G wrt D_start + D_end
        if priv.intersection(V_prev) != set():
            # if v has a private in V_prev, we keep it
            D_end.append(v)

    return D_end


def _peel(G, A):
    r'''
    Return a peeling of a bicolored graph.

    For internal use.
    Given a graph `G` and a subset `A` of its vertices, a peeling
    of `(G,A)` is a list $[(u_0, V_0), \dots, (u_{p+1}, V_{p+1})]$ such
    that $u_0$ and $u_{p+1}$ are `None`, $V_0$ is the empty set,
    $V_{p+1} = V(G)$, $V_p = A$ and for every $i \in $\{1, \dots, p\}$,
    $V_{i-1} = V_i \setminus N[v_i]$, for some $u_i\in V_i$.

    INPUT:
    
    - ``G`` -- a graph
    - ``A`` -- an iterable of vertices of ``G``

    OUTPUT:

    A peeling of `(G,A)`.

    EXAMPLES:

        sage: G = Graph(10); _peel(G, range(5))
        [(None, set()),
        (4, {4}),
        (3, {3, 4}),
        (2, {2, 3, 4}),
        (1, {1, 2, 3, 4}),
        (0, {0, 1, 2, 3, 4}),
        (None, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])]

        sage: G = graphs.PathGraph(10); _peel(G, (i for i in range(10) if i%2==0))
        [(None, set()),
        (8, {8}),
        (6, {6, 8}),
        (4, {4, 6, 8}),
        (2, {2, 4, 6, 8}),
        (0, {0, 2, 4, 6, 8}),
        (None, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])]

    '''

    Acomp = set(G.vertices())
    Acomp.difference_update(A)  # Acomp  = V - A

    peeling = [(None, G.vertices())]
    H = copy(G)
    H.delete_vertices(Acomp)
    while H.order() > 0:
        ui = next(H.vertex_iterator())  # pick some vertex of H
        Vi = set(H.vertex_iterator())
        peeling.append((ui, Vi))
        H.delete_vertices(H.closed_neighbor_iterator(ui))
    peeling.append((None, set()))
    peeling.reverse()
    return peeling


def _cand_ext_enum(G, dom, u_next, V_next):
    r'''
    Return an iterator over the candidate extensions of ``dom`` to ``V_next``.

    For internal use. (See description below.)
    
    INPUT:

    - `G` -- a graph
    - `dom` -- an iterable over some vertices of `G`
    - `u_next` -- a vertex of `G`
    - `V_next` -- an iterable over some vertices
    
    OUTPUT:

    An iterator over all sets `X` that dominate
    $N(u_next) \cap V_next \setminus N[dom]$ and are inclusion-wise
    minimal (hereafter called candidate extensions of ``dom`` to ``V_next``).
    '''

    def _aux_with_rep(G, dom, u_next, V_next):
        # Auxilliary routine. Does the same as _cand_ext_enum,
        # except that the same set may be output several times
        # (up to |G| times).
        #
        # In order to later remove duplicates, we here output pairs
        # (dom,i) where dom is the output candidate extension and i counts
        # how many elements were already output.
        #

        # S = neighbors of u_next in V_next that are not yet dominated:
        S = set.intersection(set(G.neighbor_iterator(u_next)), set(
            V_next)) - G.closed_neighbors_of_a_set(dom)

        # True iff u_next is dominated by dom:
        u_next_dom_by_dom = any(
            v in dom for v in G.closed_neighbor_iterator(u_next))

        if u_next_dom_by_dom:
            # In this case, u_next is already dominated by dom,
            # so only S has to be dominated.
            # We enumerate the minimal DSs of the bicolored graph G(S):

            cand_ext_index = 0
            for X in minimal_dominating_sets(G, list(S)):
                yield (X, cand_ext_index)
                cand_ext_index += 1

        elif not S:
            # In this case, only u_next has to be dominated
            cand_ext_index = 0
            for w in G.closed_neighbor_iterator(u_next):
                # Notice that the case w = u_next is included
                yield ({w}, cand_ext_index)
                cand_ext_index += 1

        else:
            # In this case, both u_next and S have to be dominated

            S_plus = copy(S)
            S_plus.add(u_next)  # S_plus = S + {u_next}

            yield ({u_next}, 0)  # The trivial extension
            # Start from 1 because we already output the 0-th elt:
            cand_ext_index = 1

            for w in G.neighbor_iterator(u_next):
                # Note that w never belongs to dom,
                # as we are not in the first case of the if statement

                S_minus = set.difference(
                    S, set(G.closed_neighbor_iterator(w)))
                # S_minus: vertices of S that still need to be
                # dominated, assuming w is included in the DS

                for Q in minimal_dominating_sets(G, S_minus):
                    sQ = set(Q)
                    NQ = G.closed_neighbors_of_a_set(sQ)
                    Nw_minus = set.intersection(
                        set(G.closed_neighbor_iterator(w)), S_plus)
                    if not NQ >= Nw_minus:
                        # If Nw_minus is not included in i.e. if w has
                        # a private neighbor in V_next wrt Q + {w}:
                        sQ.add(w)
                        yield (sQ, cand_ext_index)
                        cand_ext_index += 1
    #
    # End of aux_with_rep routine

    # Here we use aux_with_rep twice to enumerate the minimal
    # dominating sets while avoiding repeated outputs
    for (X, i) in _aux_with_rep(G, dom, u_next, V_next):
        for (Y, j) in _aux_with_rep(G, dom, u_next, V_next):
            if j >= i:
                # This is the first time we meet X: we output it
                yield X
                break
            elif set(Y) == set(X):
                # X has already been output in the past: we ignore it
                break


def minimal_dominating_sets(self, vertices_to_dominate=None):
    r'''
    Return an iterator over the minimal dominating sets of the graph.

    INPUT:

    - `self` -- a graph
    - `vertices_to_dominate` -- vertex iterable or None (default: `None`)

    OUTPUT:

    An iterator over the inclusion-minimal sets of vertices of `self`
    that dominate `vertices_to_dominate`.

    ALGORITHM:

    The algorithm is described in https://arxiv.org/abs/1810.00789 .

    EXAMPLES:

        sage: G = Graph()
        sage: for ds in minimal_dominating_sets(G):
               .. print(ds)
        sage: set([])

        sage: G = graphs.ButterflyGraph()
        sage: sorted(list(minimal_dominating_sets(G)))
        sage: [{0, 1}, {1, 3}, {0, 2}, {2, 3}, {4}]
        sage: sorted(list(minimal_dominating_sets(G, [0,3])))
        [{0}, {3}, {4}]
        sage: sorted(list(minimal_dominating_sets(G, [4])))
        [{4}, {0}, {1}, {2}, {3}]

    TESTS:
        
        sage: G = graphs.PetersenGraph()
        sage: for d in minimal_dominating_sets(G):
        ....:     print(d)
        {0, 2, 6}
        {1, 2, 6, 7, 9}
        {1, 2, 3, 6, 8}
        {2, 4, 5, 6}
        {0, 1, 5, 6, 8}
        {0, 1, 4, 6, 9}
        {0, 3, 6, 7}
        {3, 5, 6}
        {4, 6, 7}
        {0, 1, 2, 3, 4}
        {1, 4, 5}
        {1, 3, 7}
        {0, 1, 2, 5, 7}
        {0, 8, 2, 9}
        {8, 1, 9}
        {8, 2, 4}
        {0, 3, 4, 5, 8}
        {0, 8, 7}
        {8, 1, 4, 7}
        {2, 3, 5, 7, 8}
        {5, 6, 7, 8, 9}
        {0, 9, 3}
        {1, 3, 5, 9}
        {2, 3, 4, 7, 9}
        {3, 4, 6, 8, 9}
        {9, 2, 5}
        {0, 4, 5, 7, 9}

    '''

    def tree_search(G, plng, dom, i):
        # Internal routine
        # Recursively generates the leaves descendant of (dom,i)
        # G: graph; plng: peeling

        if i == len(plng) - 2:
            # we reached a leaf, i.e. dom is a minimal DS of vertices_to_dominate
            # '-2' because the last cell of plng is used
            # for the vertices of G - vertices_to_dominate
            yield dom
            return

        u_next, V_next = plng[i+1]

        if G.is_dominating(dom, V_next):  # if dom dominates V_{i+1}
            # then dom is its unique extension: we recurse on it
            for L in tree_search(G, plng, dom, i + 1):
                yield L
            return

        # For every candidate extension
        for can_ext in _cand_ext_enum(G, dom, u_next, V_next):

            # We complete dom with can_ext -> canD
            canD = set.union(set(can_ext), set(dom))

            if (not G.is_redundant(canD, V_next)) and set(dom) == set(_parent(G, canD, plng[i][1])):
                # If canD is a legitimate child of dom and is not
                # redundant, we recurse on it:
                for Di in tree_search(G, plng, canD, i + 1):
                    yield Di
    ##
    # end of tree-search routine

    if vertices_to_dominate is None:
        vertices_to_dominate = self.vertices()

    elif not vertices_to_dominate:  # base case: vertices_to_dominate is empty
        yield set()  # the empty set/list is the only minimal DS of the empty set of vertex
        return

    peeling = _peel(self, vertices_to_dominate)

    for dom in tree_search(self, peeling, set(), 0):
        # we generate the leaves of the search tree that are descendant of (empty set, 0)
        yield dom
