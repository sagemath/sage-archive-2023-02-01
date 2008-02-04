cdef extern from "planarity/graph.h":
    # TODO: point out to Robert how much shit I was able to erase already
    ctypedef struct graphNode:
        int v
        int link[2]
    ctypedef graphNode * graphNodeP

    ctypedef struct BM_graph:
        graphNodeP G
        int N
    ctypedef BM_graph * graphP

    cdef int OK, EMBEDFLAGS_PLANAR, NONPLANAR, NOTOK

    cdef graphP gp_New()
    cdef void gp_Free(graphP *pGraph)
    cdef int gp_InitGraph(graphP theGraph, int N)
    cdef int gp_AddEdge(graphP theGraph, int u, int ulink, int v, int vlink)
    cdef int gp_Embed(graphP theGraph, int embedFlags)
    cdef int gp_SortVertices(graphP theGraph)

def is_planar(g, set_pos=True, set_emb=True, circular=False):
    # create to and from mappings to relabel vertices to the set {0,...,n-1}
    listto = g.vertices()
    ffrom = {}
    for vvv in listto:
        ffrom[vvv] = listto.index(vvv)
    to = {}
    for i from 0 <= i < len(listto):
        to[i] = listto[i]
    g.relabel(ffrom)

    cdef graphP theGraph
    theGraph = gp_New()
    cdef int status
    status = gp_InitGraph(theGraph, g.order())
    if status != OK:
        raise RuntimeError("status does not equal ok.")
    for u, v, _ in g.edge_iterator():
        gp_AddEdge(theGraph, u, 0, v, 0)
    status = gp_Embed(theGraph, EMBEDFLAGS_PLANAR)
    gp_SortVertices(theGraph)

    # use to and from mappings to relabel vertices back from the set {0,...,n-1}
    g.relabel(to)

    if status == NOTOK:
        raise RuntimeError("not ok.")
    elif status == NONPLANAR:
        # TODO: Kuratowski subgraph isolator
        gp_Free(&theGraph)
        return False
    else:
        if not circular:
            if set_emb:
                # TODO: explain to self why the following line doesn't compile
                #cdef int i
                emb_dict = {}
                for i in range(theGraph.N):
                #for i from 0 <= i < theGraph.N:
                    linked_list = []
                    j = theGraph.G[i].link[1]
                    while j >= theGraph.N:
                        linked_list.append(to[theGraph.G[j].v])
                        j = theGraph.G[j].link[1]
                    emb_dict[to[i]] = linked_list
                g.__embedding__ = emb_dict
            if set_pos:
                schnyder(g,emb_dict)
        else:
            if set_emb:
                # Take counter-clockwise embedding if circular planar test
                # Also, pos must be set after removing extra vertex and edges
                #
                # TODO: explain to self why the following line doesn't compile
                #cdef int i
                emb_dict = {}
                for i in range(theGraph.N):
                #for i from 0 <= i < theGraph.N:
                    linked_list = []
                    j = theGraph.G[i].link[0]
                    while j >= theGraph.N:
                        linked_list.append(to[theGraph.G[j].v])
                        j = theGraph.G[j].link[0]
                    emb_dict[to[i]] = linked_list
                g.__embedding__ = emb_dict
        gp_Free(&theGraph)
        return True

def schnyder(g, emb_dict, circular=False):
    # TODO, because everyone loves a pretty picture of a graph
    from sage.graphs.graph_genus1 import trace_faces
    from sage.sets.set import Set
    # TODO: ignoring label for now, should only need to chnage this once.
    # g.relabel(ffrom)
    faces = trace_faces(g, emb_dict.values())
    """
    ext_face = faces[0] # actually should pick longest one
    if circular:
        #TODO
        # want to set external face as one with all boundary nodes in order on it
        pass

    triang_faces = []
    extra_edges = []
    for face in faces:
        if len(face) > 3:
            # use a Set to guarantee each extra edge is used in both directions

            tris = [face[2*i:2*i+2] for i in range(len(face)/2)]
            for tri in tris:
                connect = (tri[0][0],tri[-1][-1])
                tri.append(connect)
                triang_faces.append(tri)
                extra_edges.append(connect)
                g.add_edge(connect)
        else:
            triang_faces.append(face)
    """
    pass





## list from c array:
# cdef int * array
# ........
## *without* a cdef blah List ....
# (or cdef object List)
# List = [array[i] for i from 0 <= i < length_of_array]




