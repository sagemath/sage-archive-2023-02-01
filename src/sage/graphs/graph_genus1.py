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

def nice_copy(g):
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
    boundary = g.get_boundary()
    graph = g.copy()

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
    if hasattr(g,'__embedding__'):
        emb = g.__embedding__
        for vertex in emb:
            for nbr in emb[vertex]:
                nbr = relabeldict[nbr]
            vertex = relabeldict[nbr]
        graph.__embedding__ = emb
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

