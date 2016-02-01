def breadth_first_level_search(G, start, ignore_direction=False,
                               distance=None, neighbors=None):
    from sage.rings.semirings.non_negative_integer_semiring import NN
    if (distance is not None and distance not in NN):
        raise ValueError("distance must be a non-negative integer, not {0}".format(distance))
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
    G._scream_if_not_simple()
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

    g = DiGraph({v: {w: UF.find((v, w)) for w in G[v]} for v in G})
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

    # find list of tokens that lead to the initial state
    activeTokens = set()
    for LG in breadth_first_level_search(g, initialState):
        for v in LG:
            for w in LG[v]:
                activeTokens.add(g.edge_label(w, v))
    for t in activeTokens:
        if reverse[t] in activeTokens:
            return fail
    activeTokens = list(activeTokens)

    # rest of data structure: point from states to list and list to states
    activeForState = {v: -1 for v in g}
    statesForPos = [[] for i in activeTokens]

    def scan(v):
        """Find the next token that is effective for v."""
        try:
            a = next(i for i in range(activeForState[v]+1,
                                      len(activeTokens))
                     if activeTokens[i] is not None
                     and activeTokens[i] in action[v])
            activeForState[v] = a
            statesForPos[a].append(v)
        except StopIteration:
            return fail

    if certificate:
        dim = 0
        tokmap = {}
        for t in reverse:
            if t not in tokmap:
                tokmap[t] = tokmap[reverse[t]] = 1 << dim
                dim += 1
        embed = {initialState: 0}

    # set initial active states
    for v in g:
        if v != current:
            scan(v)

    # traverse the graph, maintaining active tokens
    for prev, current, fwd in depth_first_traversal(g, initialState):
        if not fwd:
            prev, current = current, prev
        elif certificate:
            embed[current] = embed[prev] ^ tokmap[g.edge_label(prev, current)]

        # add token to end of list, point to it from old state
        activeTokens.append(g.edge_label(prev, current))
        activeForState[prev] = len(activeTokens) - 1
        statesForPos.append([prev])

        # inactivate reverse token, find new token for its states
        activeTokens[activeForState[current]] = None
        for v in statesForPos[activeForState[current]]:
            if v != current:
                scan(v)

    if certificate:
        format = "{0:0%db}" % dim
        return (True, {v: format.format(l) for v, l in embed.items()})
    else:
        return True
