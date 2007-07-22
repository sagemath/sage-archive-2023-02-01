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
    Returns list of tuples:
    (genus,[list of lists representing face traces])
    """

    graph = nice_copy(graph)

    verts = len(graph.vertices())
    edges = len(graph.edges())

    # Construct a list of all rotation systems for graph
    part = []
    for node in graph.vertices():
        part.append(graph.neighbors(node))
    all_perms = []
    for p in part_perm(part):
        all_perms.append(p)

    embeddings = []
    for p in all_perms:
        t = trace_faces(graph, p)
        faces = len(t)

        g = (2-verts+edges-faces)/2
        embeddings.append((g,t))

    return embeddings

def planar_embeddings(graph):

    graph = nice_copy(graph)

    verts = len(graph.vertices())
    edges = len(graph.edges())

    # Construct a list of all rotation systems for graph
    part = []
    for node in graph.vertices():
        part.append(graph.neighbors(node))
    all_perms = []
    for p in part_perm(part):
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