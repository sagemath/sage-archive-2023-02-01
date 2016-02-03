r"""
Partial cubes

The code in this module is originally from the PADS library by David Eppstein,
which is available at http://www.ics.uci.edu/~eppstein/PADS/ under the MIT
license. The algorithm for partial cube recognition in quadratic time has been
described in [Eppstein2008]_.

REFERENCE:

.. [Eppstein2008] David Eppstein,
  "Recognizing partial cubes in quadratic time",
  J. Graph Algorithms and Applications 15 (2): 269-293, 2011.
  Available at http://arxiv.org/abs/0705.1025

Functions
---------
"""

def breadth_first_level_search(G, start, ignore_direction=False, neighbors=None):
    r"""
    Generate a sequence of dictionaries, each mapping the vertices from level i
    to a set of their neighbours at level i+1.

    Originally written by D. Eppstein for the PADS library
    (http://www.ics.uci.edu/~eppstein/PADS/).

    INPUT:

    - ``G`` -- a graph to perform the search on.

    - ``start`` -- vertex or list of vertices from which to start the traversal.

    - ``ignore_direction`` -- (default ``False``) only applies to directed
      graphs. If ``True``, searches across edges in either direction.

    - ``neighbors`` -- a function giving the neighbors of a vertex.
      The function should take a vertex and return a list of
      vertices.  For a graph, ``neighbors`` is by default the
      :meth:`~GenericGraph.neighbor_iterator` function of the graph. For a
      digraph, the ``neighbors`` function defaults to the
      :meth:`~DiGraph.neighbor_out_iterator` function of the graph.

    EXAMPLE::

        sage: H = graphs.HeawoodGraph()
        sage: list(sage.graphs.partial_cube.breadth_first_level_search(H, 0))
        [{0: {1, 5, 13}},
         {1: {2, 10}, 5: {4, 6}, 13: {8, 12}},
         {2: {3, 7}, 4: {3, 9}, 6: {7, 11}, 8: {7, 9}, 10: {9, 11}, 12: {3, 11}},
         {3: set(), 7: set(), 9: set(), 11: set()}]

    """
    if neighbors is None:
        if not G._directed or ignore_direction:
            neighbors = G.neighbor_iterator
        else:
            neighbors = G.neighbor_out_iterator
    visited = set()
    if isinstance(start, list):
        currentLevel = start
    else:
        currentLevel = [start]
    while currentLevel:
        visited.update(currentLevel)
        nextLevel = set()
        levelGraph = {v: set() for v in currentLevel}
        for v in currentLevel:
            for w in neighbors(v):
                if w not in visited:
                    levelGraph[v].add(w)
                    nextLevel.add(w)
        yield levelGraph
        currentLevel = nextLevel

def depth_first_traversal(G, start, ignore_direction=False,
                          neighbors=None):
    r"""
    Generate a sequence of triples (v,w,edgetype) for DFS of graph G.

    Originally written by D. Eppstein for the PADS library
    (http://www.ics.uci.edu/~eppstein/PADS/).

    INPUT:

    - ``G`` -- a graph to perform the search on.

    - ``start`` -- vertex or list of vertices from which to start the traversal.

    - ``ignore_direction`` -- (default False) only applies to directed graphs.
      If True, searches across edges in either direction.

    - ``neighbors`` -- a function giving the neighbors of a vertex.
      The function should take a vertex and return a list of
      vertices.  For a graph, ``neighbors`` is by default the
      :meth:`~GenericGraph.neighbor_iterator` function of the graph. For a
      digraph, the ``neighbors`` function defaults to the
      :meth:`~DiGraph.neighbor_out_iterator` function of the graph.

    OUTPUT:

    - a generator of triples ``(v,w,edgetype)``, where ``edgetype`` is ``True``
      if the algorithm is progressing via the edge ``vw``, or ``False`` if the
      algorithm is backtracking via the edge ``wv``.

    EXAMPLE::

        sage: H = graphs.HeawoodGraph()
        sage: t = list(sage.graphs.partial_cube.depth_first_traversal(H, 0))
        sage: len(t)
        26

    """
    if neighbors is None:
        if not G._directed or ignore_direction:
            neighbors=G.neighbor_iterator
        else:
            neighbors=G.neighbor_out_iterator
    else:
        neighbors = lambda v: iter(neighbors(v))
    seen=set([])
    if not isinstance(start, list):
        start = [start]

    for v in start:
        if v in seen:
            continue
        seen.add(v)
        stack = [(v, neighbors(v))]
        while stack:
            parent, children = stack[-1]
            try:
                child = next(children)
                if child not in seen:
                    yield (parent, child, True)
                    seen.add(child)
                    stack.append((child, neighbors(child)))
            except StopIteration:
                stack.pop()
                if stack:
                    yield (stack[-1][0], parent, False)

def is_partial_cube(G, certificate=False):
    r"""
    Test whether the given graph is a partial cube.

    A partial cube is a graph that can be isometrically embedded into a
    hypercube, i.e., its vertices can be labelled with (0,1)-vectors of some
    fixed length such that the distance between any two vertices in the graph
    equals the Hamming distance of their labels.

    Originally written by D. Eppstein for the PADS library
    (http://www.ics.uci.edu/~eppstein/PADS/), see also [Eppstein2008]_.
    The algorithm runs in O(n^2) time, where n is the number of vertices.

    INPUT:

    - ``certificate`` (boolean; ``False``) -- The function returns ``True``
      or ``False`` according to the graph, when ``certificate = False``. When
      ``certificate = True`` and the graph is a partial cube, the function
      returns ``(True, mapping)``, where ``mapping`` is an isometric mapping of
      the vertices of the graph to the vertices of a hypercube ((0, 1)-strings
      of a fixed length). When ``certificate = True`` and the graph is not a
      partial cube, ``(False, None)`` is returned.

    EXAMPLES:

    The Petersen graph is not a partial cube::

        sage: g = graphs.PetersenGraph()
        sage: g.is_partial_cube()
        False

    All prisms are partial cubes::

        sage: g = graphs.CycleGraph(10).cartesian_product(graphs.CompleteGraph(2))
        sage: g.is_partial_cube()
        True

    TESTS:

    The returned mapping is an isometric embedding into a hypercube::

        sage: g = graphs.DesarguesGraph()
        sage: _, m = g.is_partial_cube(certificate=True)
        sage: m # random
        {0: '00000',
         1: '00001',
         2: '00011',
         3: '01011',
         4: '11011',
         5: '11111',
         6: '11110',
         7: '11100',
         8: '10100',
         9: '00100',
         10: '01000',
         11: '10001',
         12: '00111',
         13: '01010',
         14: '11001',
         15: '10111',
         16: '01110',
         17: '11000',
         18: '10101',
         19: '00110'}
        sage: all(all(g.distance(u, v) == len([i for i in range(len(m[u])) if m[u][i] != m[v][i]]) for v in m) for u in m)
        True

    A graph without vertices is trivially a partial cube::

        sage: Graph().is_partial_cube(certificate = True)
        (True, {})

    """
    G._scream_if_not_simple()

    if G.order() == 0:
        if certificate:
            return (True, {})
        else:
            return True

    if certificate:
        fail = (False, None)
    else:
        fail = False

    if not G.is_connected():
        return fail
    n = G.order()

    # Initial sanity check: are there few enough edges?
    # Needed so that we don't try to use union-find on a dense
    # graph and incur superquadratic runtimes.
    if 1 << (2*G.size()//n) > n:
        return fail

    # Set up data structures for algorithm:
    # - CG: contracted graph at current stage of algorithm
    # - UF: union find data structure representing known edge equivalences
    # - NL: limit on number of remaining available labels
    from sage.graphs.digraph import DiGraph
    from sage.graphs.graph import Graph
    from sage.sets.disjoint_set import DisjointSet
    CG = DiGraph({v: {w: (v, w) for w in G[v]} for v in G})
    UF = DisjointSet(CG.edges(labels = False))
    NL = n-1

    # Main contraction loop in place of the original algorithm's recursion
    while CG.order() > 1:
        if not Graph(CG).is_bipartite():
            return fail

        # Find max degree vertex in CG, and update label limit
        deg, root = max((len(CG[v]), v) for v in CG)
        if deg > NL:
            return fail
        NL -= deg

        # Set up bitvectors on vertices
        bitvec = {v:0 for v in CG}
        neighbors = {}
        for i, neighbor in enumerate(CG[root]):
            bitvec[neighbor] = 1 << i
            neighbors[1 << i] = neighbor

        # Breadth first search to propagate bitvectors to the rest of the graph
        for LG in breadth_first_level_search(CG, root):
            for v in LG:
                for w in LG[v]:
                    bitvec[w] |= bitvec[v]

        # Make graph of labeled edges and union them together
        labeled = Graph([CG.vertices(), []])
        for v, w in CG.edge_iterator(labels = False):
            diff = bitvec[v]^bitvec[w]
            if not diff or bitvec[w] &~ bitvec[v] == 0:
                continue    # zero edge or wrong direction
            if diff not in neighbors:
                return fail
            neighbor = neighbors[diff]
            UF.union(CG.edge_label(v, w), CG.edge_label(root, neighbor))
            UF.union(CG.edge_label(w, v), CG.edge_label(neighbor, root))
            labeled.add_edge(v, w)

        # Map vertices to components of labeled-edge graph
        component = {}
        for i, SCC in enumerate(labeled.connected_components()):
            for v in SCC:
                component[v] = i

        # generate new compressed subgraph
        NG = DiGraph(labeled.connected_components_number())
        for v, w, t in CG.edge_iterator():
            if bitvec[v] == bitvec[w]:
                vi = component[v]
                wi = component[w]
                if vi == wi:
                    return fail
                if wi in NG.neighbors_out(vi):
                    UF.union(NG.edge_label(vi, wi), t)
                else:
                    NG.add_edge(vi, wi, t)
        CG = NG

    # Make a digraph with edges labeled by the equivalence classes in UF
    g = DiGraph({v: {w: UF.find((v, w)) for w in G[v]} for v in G})

    # Check that no two edges on a single vertex have the same label
    action = {v: {} for v in g}
    reverse = {}
    for v, w, t in g.edge_iterator():
        if t in action[v]:
            return fail
        action[v][t] = w
        rt = g.edge_label(w, v)
        if t not in reverse:
            reverse[t] = rt
            reverse[rt] = t
    current = initialState = next(g.vertex_iterator())

    # Find list of tokens that lead to the initial state
    activeTokens = set()
    for LG in breadth_first_level_search(g, initialState):
        for v in LG:
            for w in LG[v]:
                activeTokens.add(g.edge_label(w, v))
    for t in activeTokens:
        if reverse[t] in activeTokens:
            return fail
    activeTokens = list(activeTokens)

    # Rest of data structure: point from states to list and list to states
    activeForState = {v: -1 for v in g}
    statesForPos = [[] for i in activeTokens]

    def scan(v):
        """Find the next token that is effective for v."""
        a = next(i for i in range(activeForState[v]+1, len(activeTokens))
                 if activeTokens[i] is not None
                    and activeTokens[i] in action[v])
        activeForState[v] = a
        statesForPos[a].append(v)

    # Initialize isometric embedding into a hypercube
    if certificate:
        dim = 0
        tokmap = {}
        for t in reverse:
            if t not in tokmap:
                tokmap[t] = tokmap[reverse[t]] = 1 << dim
                dim += 1
        embed = {initialState: 0}

    # Set initial active states
    for v in g:
        if v != current:
            try:
                scan(v)
            except StopIteration:
                return fail

    # Traverse the graph, maintaining active tokens
    for prev, current, fwd in depth_first_traversal(g, initialState):
        if not fwd:
            prev, current = current, prev
        elif certificate:
            embed[current] = embed[prev] ^ tokmap[g.edge_label(prev, current)]

        # Add token to end of list, point to it from old state
        activeTokens.append(g.edge_label(prev, current))
        activeForState[prev] = len(activeTokens) - 1
        statesForPos.append([prev])

        # Inactivate reverse token, find new token for its states
        activeTokens[activeForState[current]] = None
        for v in statesForPos[activeForState[current]]:
            if v != current:
                try:
                    scan(v)
                except StopIteration:
                    return fail

    # All checks passed, return the result
    if certificate:
        format = "{0:0%db}" % dim
        return (True, {v: format.format(l) for v, l in embed.items()})
    else:
        return True
