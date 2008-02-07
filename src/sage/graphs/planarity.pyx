from sage.sets.set import Set
from graph_genus1 import trace_faces

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

def triangulate(g, comb_emb):
    cdef i

    faces = trace_faces(g, comb_emb)
#    edges_used = Set([])
#    edges_free = Set([])

    for face in faces:
        print face

    print 'starting to loop:'

    # first check if graphs is one of
    # single vertex
    # two vertices
    # o--o--o
    # not connected?

    for face in faces:
        print face
        new_face = []
        if len(face) < 3:
            print 'oops'
            return
        if len(face) == 3:
            continue
        elif len(face) == 4:
            new_face = (face[1][1], face[0][0])
            if g.has_edge(new_face):
                new_face = (face[2][1], face[1][0])
            g.add_edge(new_face)
        else:
            N = len(face)
            i = 0
            while i < N-1:
                new_edge = (face[i+1][1], face[i][0])
                if g.has_edge(new_edge) or new_edge[0] == new_edge[1]:
                    new_face.append(face[i])
                    if i == N - 2:
                        break
                    i = i + 1
                    continue
                    #new_edge = (face[i+1][1], face[i][0])

                g.add_edge(new_edge)
                new_face.append((new_edge[1], new_edge[0]))
                i = i + 2
            if i != N:
                new_face.append(face[-1])
            faces.append(new_face)

def normal_label(g, comb_emb, external_face):
    print 'using external face', external_face
    contracted = []
    contractible = []

    v1 = external_face[0][0]
    v1_neighbors = Set(g.neighbors(v1))

    neighbor_count = {}
    for v in g.vertices():
        neighbor_count[v] = len(v1_neighbors.intersection( Set(g.neighbors(v))))

    for v in v1_neighbors:
        if v == external_face[1][0] or v == external_face[2][0]:
            continue
        if neighbor_count[v] == 2:
            contractible.append(v)

    while g.order() > 3:
        try:
            v = contractible.pop()
        except:
            print 'fart'
            break
        print 'going to contract', v
        v_neighbors = Set(g.neighbors(v))
        contracted.append( (v, v_neighbors, v_neighbors - v1_neighbors - Set([v1])) )
        g.delete_vertex(v)
        g.show()
        v1_neighbors -= Set([v])
        for w in v_neighbors - v1_neighbors - Set([v1]):
            print 'adding edge:', v1, w
            g.add_edge( (v1, w) )
            g.show()
        if g.order() == 3:
            break
        v1_neighbors += v_neighbors - Set([v1])
        for w in v_neighbors - Set([v1]):
            print 'neighbors of', w, ':', g.neighbors(w)
            print 'neighbors of', v1,':', v1_neighbors, g.neighbors(v1)
            new_neighbor_count = len(v1_neighbors.intersection( Set(g.neighbors(w))))
            print w, 'neighbor count is now', new_neighbor_count
            if new_neighbor_count != neighbor_count[w]:
                if new_neighbor_count == 2:
                    contractible.append(w)
                    print 'declaring', w, 'contractible'
                elif neighbor_count[w] == 2:
                    contractible.remove(w)
                    print 'declaring', w, 'no longer contractible'
            neighbor_count[w] = new_neighbor_count

    while len(contracted) > 0:
        v, new_neighbors, neighbors_to_delete = contracted.pop()
        for w in new_neighbors:
            g.add_edge((v,w))
        for w in neighbors_to_delete:
            g.delete_edge((v1,w))
        g.show()

def schnyder(g, emb_dict, circular=False):
    # TODO
    pass



