r"""
A module for computing the genus and possible embeddings of a graph.

AUTHOR:
    -- Emily A. Kirkman (2007-07-21): initial version

"""

#*****************************************************************************
#                     Copyright (C) 2007 Emily A. Kirkman
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

def nice_copy(graph):
    """
    Creates a 'nice' copy of the graph (with vertices labeled
    0 through n-1).  Also copies the boundary to be in the
    same order with (possibly) new vertex labels.

    INPUT:
        graph -- the graph to make a nice copy of

    EXAMPLES:
        sage: from sage.graphs import graph_genus1
        sage: G = graphs.PetersenGraph()
        sage: G.set_boundary([0,1,2,3,4])
        sage: H = graph_genus1.nice_copy(G)
        sage: H == G
        True
        sage: J = Graph({'alpha':['beta', 'epsilon'], 'gamma':['beta', 'epsilon']})
        sage: K = graph_genus1.nice_copy(J)
        sage: K == J
        False
        sage: J.set_boundary(['beta','alpha'])
        sage: L = graph_genus1.nice_copy(J)
        sage: L.get_boundary()
        [1, 0]
    """
    boundary = graph.get_boundary()
    graph = graph.copy()

    newboundary = []
    relabeldict = {}
    i = 0
    for v in graph.vertices():
        relabeldict[v] = i
        i += 1
    if boundary is not None:
        for v in boundary:
            newboundary.append(relabeldict[v])
    graph.relabel(relabeldict)
    graph.set_boundary(newboundary)
    return graph

def trace_faces(graph, rot_sys):
    """
    A helper function for finding the genus of a graph.
    Given a graph and a combinatorial embedding (rot_sys),
    this function will compute the faces (returned as a list
    of lists of edges (tuples) of the particular embedding.

    Note -- rot_sys is an ordered list based on the hash order
    of the vertices of graph.  To avoid confusion, it might be
    best to set the rot_sys based on a 'nice_copy' of the graph.

    INPUT:
        graph -- a graph to compute the faces of
        rot_sys -- the combinatorial embedding of graph (an
                   ordered list by vertex hash order)

    EXAMPLES:
        sage: from sage.graphs import graph_genus1
        sage: J = Graph({'alpha':['beta', 'epsilon'], 'gamma':['beta', 'epsilon']})
        sage: J.set_boundary(['beta','alpha'])
        sage: K = graph_genus1.nice_copy(J)
        sage: rot = []
        sage: for node in K.vertices():
        ...     rot.append(K[node])
        sage: rot
        [[1, 2], [0, 3], [0, 3], [1, 2]]
        sage: graph_genus1.trace_faces(K,rot)
        [[(0, 1), (1, 3), (3, 2), (2, 0)], [(1, 0), (0, 2), (2, 3), (3, 1)]]
    """
    from sage.sets.set import Set

    # Make dict of node labels embedding
    comb_emb = {}
    labels = graph.vertices()
    for i in range(len(rot_sys)):
        comb_emb[labels[i]] = rot_sys[i]

    # Establish set of possible edges
    edgeset = Set([])
    for edge in graph.edges():
        edgeset = edgeset.union(Set([(edge[0],edge[1]),(edge[1],edge[0])]))

    # Storage for face paths
    faces = []
    path = []
    for edge in edgeset:
        path.append(edge)
        edgeset -= Set([edge])
        break  # (Only one iteration)

    # Trace faces
    while (len(edgeset) > 0):
        neighbors = comb_emb[path[-1][-1]]
        next_node = neighbors[(neighbors.index(path[-1][-2])+1)%(len(neighbors))]
        tup = (path[-1][-1],next_node)
        if tup == path[0]:
            faces.append(path)
            path = []
            for edge in edgeset:
                path.append(edge)
                edgeset -= Set([edge])
                break  # (Only one iteration)
        else:
            path.append(tup)
            edgeset -= Set([tup])

    if (len(path) != 0): faces.append(path)
    return faces

def all_embeddings(graph):
    """
    Returns a list of tuples, one for each possible embedding.
    The tuples have the minimal genus of the particular
    embedding as the first entry.  The second entry is a list
    of lists of edges (tuples) that represent the embedding
    via face traces.  Each inner list represents one face in
    the embedding.

    Note -- returns list of tuples:
        (genus,[list of lists representing face traces])

    INPUT:
        graph -- the graph to find all possible embeddings of

    EXAMPLES:
        sage: from sage.graphs import graph_genus1
        sage: J = Graph({'alpha':['beta', 'epsilon'], 'gamma':['beta', 'epsilon']})
        sage: J.set_boundary(['beta','alpha'])
        sage: graph_genus1.all_embeddings(J)
        [(0, [[(0, 1), (1, 3), (3, 2), (2, 0)], [(1, 0), (0, 2), (2, 3), (3, 1)]])]
        sage: K23 = graphs.CompleteBipartiteGraph(2,3)
        sage: graph_genus1.all_embeddings(K23)
        [(1,
          [[(1, 2),
            (2, 0),
            (0, 3),
            (3, 1),
            (1, 4),
            (4, 0),
            (0, 2),
            (2, 1),
            (1, 3),
            (3, 0),
            (0, 4),
            (4, 1)]]),
         (0,
          [[(1, 2), (2, 0), (0, 4), (4, 1)],
           [(1, 3), (3, 0), (0, 2), (2, 1)],
           [(0, 3), (3, 1), (1, 4), (4, 0)]]),
         (0,
          [[(1, 2), (2, 0), (0, 3), (3, 1)],
           [(1, 3), (3, 0), (0, 4), (4, 1)],
           [(1, 4), (4, 0), (0, 2), (2, 1)]]),
         (1,
          [[(1, 2),
            (2, 0),
            (0, 4),
            (4, 1),
            (1, 3),
            (3, 0),
            (0, 2),
            (2, 1),
            (1, 4),
            (4, 0),
            (0, 3),
            (3, 1)]])]
    """
    from sage.combinat.all import CyclicPermutationsOfPartition

    graph = nice_copy(graph)

    verts = len(graph.vertices())
    edges = len(graph.edges())

    # Construct a list of all rotation systems for graph
    part = []
    for node in graph.vertices():
        part.append(graph.neighbors(node))
    all_perms = []
    for p in CyclicPermutationsOfPartition(part):
        all_perms.append(p)

    embeddings = []
    for p in all_perms:
        t = trace_faces(graph, p)
        faces = len(t)

        g = (2-verts+edges-faces)/2
        embeddings.append((g,t))

    return embeddings

def planar_embeddings(graph):
    """
    Returns a list of lists of lists of edges (tuples), where
    each inner inner list of edges represents the tracing of
    one face in the embedding, and each inner list of lists
    represents one planar embedding.  The list is an
    exhaustive list of all embeddings of graph with genus=0.
    Returns an empty list if there are no planar embeddings
    of the graph.

    INPUT:
        graph -- the graph to find all planar embeddings of

    EXAMPLES:
        sage: from sage.graphs import graph_genus1
        sage: K5 = graphs.CompleteGraph(5)
        sage: graph_genus1.planar_embeddings(K5)
        []

        sage: K23 = graphs.CompleteBipartiteGraph(2,3)
        sage: graph_genus1.planar_embeddings(K23)
        [[[(1, 2), (2, 0), (0, 4), (4, 1)],
          [(1, 3), (3, 0), (0, 2), (2, 1)],
          [(0, 3), (3, 1), (1, 4), (4, 0)]],
         [[(1, 2), (2, 0), (0, 3), (3, 1)],
          [(1, 3), (3, 0), (0, 4), (4, 1)],
          [(1, 4), (4, 0), (0, 2), (2, 1)]]]
        sage: graph_genus1.all_embeddings(K23)
        [(1,
          [[(1, 2),
            (2, 0),
            (0, 3),
            (3, 1),
            (1, 4),
            (4, 0),
            (0, 2),
            (2, 1),
            (1, 3),
            (3, 0),
            (0, 4),
            (4, 1)]]),
         (0,
          [[(1, 2), (2, 0), (0, 4), (4, 1)],
           [(1, 3), (3, 0), (0, 2), (2, 1)],
           [(0, 3), (3, 1), (1, 4), (4, 0)]]),
         (0,
          [[(1, 2), (2, 0), (0, 3), (3, 1)],
           [(1, 3), (3, 0), (0, 4), (4, 1)],
           [(1, 4), (4, 0), (0, 2), (2, 1)]]),
         (1,
          [[(1, 2),
            (2, 0),
            (0, 4),
            (4, 1),
            (1, 3),
            (3, 0),
            (0, 2),
            (2, 1),
            (1, 4),
            (4, 0),
            (0, 3),
            (3, 1)]])]

        sage: g = graphs.CycleGraph(9)
        sage: graph_genus1.planar_embeddings(g)
        [[[(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 0), (0, 1)],
          [(5, 4), (4, 3), (3, 2), (2, 1), (1, 0), (0, 8), (8, 7), (7, 6), (6, 5)]]]
        sage: graph_genus1.all_embeddings(g)
        [(0,
          [[(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 0), (0, 1)],
           [(5, 4), (4, 3), (3, 2), (2, 1), (1, 0), (0, 8), (8, 7), (7, 6), (6, 5)]])]
    """
    from sage.combinat.all import CyclicPermutationsOfPartition

    graph = nice_copy(graph)

    verts = len(graph.vertices())
    edges = len(graph.edges())

    # Construct a list of all rotation systems for graph
    part = []
    for node in graph.vertices():
        part.append(graph.neighbors(node))
    all_perms = []
    for p in CyclicPermutationsOfPartition(part):
        all_perms.append(p)

    # planar embeddings
    plan_emb = []

    for p in all_perms:
        t = trace_faces(graph, p)
        faces = len(t)

        # return planar embeddings
        if ((2-verts+edges-faces)/2 == 0):
            plan_emb.append(t)

    return plan_emb
