from sage.sets.set import Set
from graph_genus1 import trace_faces
from graph import DiGraph

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
    """
    Given a connected graph g with at least 3 vertices and a planar combinatorial
    embedding comb_emb of g, modify g in place to form a graph whose faces are
    all triangles, and return the set of newly created edges.
    """

    #print 'starting to loop:'

    # first make sure that the graph has at least 3 vertices, and that it is connected
    if g.order() < 3:
        raise NotImplementedError("triangulate() only accepts graphs with more than 2 vertices as input.")
    if g.is_connected == False:
        raise NotImplementedError("triangulate() only knows how to handle connected graphs.")

    if g.order() == 3 and len(g.edges()) == 2:                # if g is o--o--o
        vertex_list = g.vertices()
        if len(g.neighbors(vertex_list[0]) == 2):                     # figure our which of the vertices already has two neighbors
            new_edge = (vertex_list[1], vertex_list[2])    # and connect the other two together.
        elif len(g.neighbors(vertex_list[1])) == 2:
            new_edge = (vertex_list[0], vertex_list[2])
        else:
            new_edge = (vertex_list[0], vertex_list[1])

        g.add_edge(new_edge)
        return [new_edge]

    # At this point we know that the graph is connected, has at least 3 vertices, and
    # that it is not the graph o--o--o. This is where the real work starts.

    cdef i                                  # i is just an index variable that we will use for looping.

    faces = trace_faces(g, comb_emb)        # We start by finding all of the faces of this embedding.

    edges_added = []                        # The list of edges that we add to the graph.
                                            # This will be returned at the end.


    # The simple way to triangular a face is to just pick a vertex and draw
    # an edge from that vertex to every other vertex in the face. Think that this
    # might ultimately result in graphs that don't look very nice when we draw them
    # we have decided on a different strategy.

    for face in faces:
        #print face
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
            edges_added.append(new_face)
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
                edges_added.append(new_edge)
                new_face.append((new_edge[1], new_edge[0]))
                i = i + 2
            if i != N:
                new_face.append(face[-1])
            faces.append(new_face)

    return edges_added

def normal_label(g, comb_emb, external_face):
    #print 'using external face', external_face
    contracted = []
    contractible = []

    labels = {}

    external_vertices = [external_face[0][0], external_face[1][0], external_face[2][0]]
    external_vertices.sort()
    v1,v2,v3 = external_vertices
    v1_neighbors = Set(g.neighbors(v1))

    neighbor_count = {}
    for v in g.vertices():
        neighbor_count[v] = 0
    for v in g.neighbors(v1):
        neighbor_count[v] = len(v1_neighbors.intersection( Set(g.neighbors(v))))


    for v in v1_neighbors:
        if v in [v1,v2,v3]:
            continue
        if neighbor_count[v] == 2:
            contractible.append(v)

    # contraction phase:

    while g.order() > 3:
        try:
            v = contractible.pop()
        except:
            print 'fart'

            break
        #print 'going to contract', v
        v_neighbors = Set(g.neighbors(v))
        contracted.append( (v, v_neighbors, v_neighbors - v1_neighbors - Set([v1])) )
        g.delete_vertex(v)
        #g.show()
        v1_neighbors -= Set([v])
        for w in v_neighbors - v1_neighbors - Set([v1]):
            #print 'adding edge:', v1, w
            g.add_edge( (v1, w) )
            #g.show()
        if g.order() == 3:
            break
        v1_neighbors += v_neighbors - Set([v1])
        contractible = []
        for w in g.neighbors(v1):
            if(len(v1_neighbors.intersection( Set(g.neighbors(w))))) == 2 and w not in [v1, v2, v3]:
                contractible.append(w)

#        for w in v_neighbors - Set([v1]):
#            print 'neighbors of', w, ':', g.neighbors(w)
#            print 'neighbors of', v1,':', v1_neighbors, g.neighbors(v1)
#            new_neighbor_count = len(v1_neighbors.intersection( Set(g.neighbors(w))))
#            print neighbor_count
#            print w, 'neighbor count is now', new_neighbor_count
#            if new_neighbor_count != neighbor_count[w] and (w not in [v1,v2,v3]):
#                if new_neighbor_count == 2:
#                    print 'declaring', w, 'contractible'
#                    contractible.append(w)
#                elif neighbor_count[w] == 2:
#           print 'declaring', w, 'no longer contractible'
#                    contractible.remove(w)
#            neighbor_count[w] = new_neighbor_count

    # expansion phase:

    v1, v2, v3 = g.vertices() #always in sorted order

    labels[v1] = {(v2,v3):1}
    labels[v2] = {(v1,v3):2}
    labels[v3] = {(v1,v2):3}

    while len(contracted) > 0:
        v, new_neighbors, neighbors_to_delete = contracted.pop()
        #print 'going to add back vertex', v
        #print v, new_neighbors, neighbors_to_delete
        labels[v] = {}

        for w in neighbors_to_delete:
            g.delete_edge((v1,w))

        if len(neighbors_to_delete) == 0:
            # we are adding v into the face new_neighbors
            w1, w2, w3 = sorted(new_neighbors)
            #print '3 neighbors:', w1, w2, w3
            labels[v] = {(w1, w2): labels[w3].pop((w1,w2)), (w2,w3) : labels[w1].pop((w2,w3)), (w1,w3) : labels[w2].pop((w1,w3))}
            labels[w1][tuple(sorted((w2,v)))] = labels[v][(w2,w3)]
            labels[w1][tuple(sorted((w3,v)))] = labels[v][(w2,w3)]

            labels[w2][tuple(sorted((w1,v)))] = labels[v][(w1,w3)]
            labels[w2][tuple(sorted((w3,v)))] = labels[v][(w1,w3)]

            labels[w3][tuple(sorted((w1,v)))] = labels[v][(w1,w2)]
            labels[w3][tuple(sorted((w2,v)))] = labels[v][(w1,w2)]
        else:
            new_neighbors_set = Set(new_neighbors)
            angles_out_of_v1 = set()
            vertices_in_order = []
            l = []
            for angle in labels[v1].keys():
                if len(Set(angle).intersection(new_neighbors_set)) == 2:
                    angles_out_of_v1.add(angle)
                    l = l + list(angle)
            # find a unique element in l
            l.sort()
            i = 0
            while i < len(l):
                if l[i] == l[i+1]:
                    i = i + 2
                else:
                    break

            angle_set = Set(angles_out_of_v1)

            vertices_in_order.append(l[i])
            while len(angles_out_of_v1) > 0:
                for angle in angles_out_of_v1:
                    if vertices_in_order[-1] in angle:
                        break
                if angle[0] == vertices_in_order[-1]:
                    vertices_in_order.append(angle[1])
                else:
                    vertices_in_order.append(angle[0])
                angles_out_of_v1.remove(angle)

            w = vertices_in_order

            # is w[0] a 2 or a 3?
            top_label = labels[w[0]][tuple(sorted((v1, w[1])))]
            if top_label == 3:
                bottom_label = 2
            else:
                bottom_label = 3
            i = 0
            while i < len(w) - 1:
                labels[v][ tuple(sorted((w[i],w[i+1]))) ] = 1
                labels[w[i]][ tuple(sorted( (w[i+1], v) )) ] = top_label
                labels[w[i+1]][ tuple(sorted( (w[i], v) )) ] = bottom_label
                i = i + 1

            labels[v][tuple(sorted( (v1, w[0])))] = bottom_label
            labels[v][tuple(sorted( (v1, w[-1])))] = top_label


            labels[w[0]][tuple(sorted((v1,v)))] = top_label
            labels[w[-1]][tuple(sorted((v1,v)))] = bottom_label
            labels[v1][tuple(sorted( (w[0],v) ))] = 1
            labels[v1][tuple(sorted( (w[-1],v) ))] = 1

            #delete all the extra labels

            for angle in angle_set:
                labels[v1].pop( angle )

            labels[w[0]].pop( tuple (sorted( (v1, w[1]) ) ))
            labels[w[-1]].pop( tuple (sorted( (v1, w[-2]) ) ))

            i = 1
            while i < len(w) - 1:
                labels[w[i]].pop(tuple(sorted( (v1, w[i+1]))))
                labels[w[i]].pop(tuple(sorted( (v1, w[i-1]))))
                i = i + 1

        #for labeling in labels:
            #print labeling, labels[labeling]

        for w in new_neighbors:
            g.add_edge((v,w))
        #g.show()

    return labels, (v1, v2, v3)

def realizer(g, x):
    normal_labeling, (v1, v2, v3) = x
    """
    Given a graph g and a normal labeling return the realizer
    of that graph and normal labeling.
    """

    realizer = DiGraph()

    tree_nodes = {}
    for v in g:
        tree_nodes[v] = [TreeNode(label = v, children = []), TreeNode(label = v, children = []), TreeNode(label = v, children = [])]


    for v in g:
        ones = []
        twos = []
        threes = []
        l = [ones,twos,threes]
        for angle, value in normal_labeling[v].items():
            l[value - 1] += list(angle)

        ones.sort()
        twos.sort()
        threes.sort()

        i = 0
        while i < len(ones) - 1:
            if ones[i] == ones[i+1]:
                realizer.add_edge( (ones[i], v), label=1)
    #            print 'adding edge in tree from', v, 'to', ones[i]
                tree_nodes[v][0].append_child(tree_nodes[ones[i]][0])
    #            tree_nodes[v][0].pretty_print()
                i = i + 1
            i = i + 1
        i = 0
        while i < len(twos) - 1:
            if twos[i] == twos[i+1]:
                realizer.add_edge( (twos[i], v), label=2)
                tree_nodes[v][1].append_child(tree_nodes[twos[i]][1])
                i = i + 1
            i = i + 1
        i = 0
        while i < len(threes) - 1:
            if threes[i] == threes[i+1]:
                realizer.add_edge( (threes[i], v), label=3)
                tree_nodes[v][2].append_child(tree_nodes[threes[i]][2])
                i = i + 1
            i = i + 1

    compute_coordinates(realizer, (tree_nodes, (v1, v2, v3)))

    realizer.show(edge_labels=True)


#    tree_nodes[v1.label][0].pretty_print()
#    tree_nodes[v2.label][1].pretty_print()
#    tree_nodes[v3.label][2].pretty_print()

    return tree_nodes, (v1, v2, v3)

def compute_coordinates(g, x):

    #print "ENTERING COMPUTE COORDINATES"

    tree_nodes, (v1, v2, v3) = x
    # find the roots of each tree:
#    t1, t2, t3 = tree_nodes[g.vertices()[0]]
    t1, t2, t3 = tree_nodes[v1][0], tree_nodes[v2][1], tree_nodes[v3][2]
    #while t1.parent is not None:
    #    t1 = t1.parent
    #while t2.parent is not None:
    #    t2 = t2.parent
    #while t3.parent is not None:
    #    t3 = t3.parent

    #print "\tUSING ROOT VERTICES: ", t1, t2, t3


    t1.compute_number_of_descendants()
    t2.compute_number_of_descendants()
    t3.compute_number_of_descendants()

    t1.compute_depth_of_self_and_children()
    t2.compute_depth_of_self_and_children()
    t3.compute_depth_of_self_and_children()

    #print "TREE 1:"

    #t1.pretty_print()

    #print "TREE 2:"

    #t2.pretty_print()

    #print "TREE 3:"

    #t3.pretty_print()

    coordinates = {}
    coordinates[t1.label] = [2 * g.order() - 5, 0]
    coordinates[t2.label] = [0, 2 * g.order() - 5]
    coordinates[t3.label] = [0, 0]

    coordinates[t1.label] = [g.order() - 2, 1]
    coordinates[t2.label] = [0, g.order() - 2]
    coordinates[t3.label] = [1, 0]

    for v in g.vertices():
        if v not in [t1.label,t2.label,t3.label]:
            #print "Computing coordinates for ", v
            r = list((0,0,0))
            for i in [0,1,2]:
                #print "\t Computing size of region", i
                p = tree_nodes[v][(i + 1) % 3]
#                label_list = []
                #print "\t\t Tracing up tree ", (i + 1) % 3
                while p is not None:
#                    label_list.append(p.label)
                    q = tree_nodes[p.label][i].number_of_descendants
                    #print "\t\t\t Adding", q
                    r[i] += q
                    p = p.parent
#                print label_list
                p = tree_nodes[v][(i - 1) % 3]
#                label_list = []
                #print "\t\t Tracing up tree ", (i - 1) % 3
                while p is not None:
#                    label_list.append(p.label)
                    q = tree_nodes[p.label][i].number_of_descendants# - tree_nodes[v][i].number_of_descendants
                    #print "\t\t\t Adding", q
                    r[i] += q
                    p = p.parent
#                print label_list

                q = tree_nodes[v][i].number_of_descendants
                #print "\t\t Subtracting, ", q
                r[i] -= q
                #print "\t\t Subtracting, ", q
                q = tree_nodes[v][(i-1)%3].depth
                r[i] -= q
#            for i in [0,1,2]:
#                print "Adjusting vertex counts to get triangle counts:"
#                r[i] = 2 * r[i] - 5
#                print "\t Triangle counts:", r

            #print v, r, sum(r)
            if sum(r) != g.order() - 1:
                raise Exception("fuck")

            coordinates[v] = r[:-1]

    #for key, value in coordinates.items():
        #print key, value

    g.set_pos(coordinates)
    #g.show()

class TreeNode():
    def __init__(self, parent = None, children = None, label = None):
        if children is None:
            children = []
        self.parent = parent
        self.children = children
        self.label = label
        self.number_of_descendants = 1


    def compute_number_of_descendants(self):
        n = 1
        for child in self.children:
            n += child.compute_number_of_descendants()
        self.number_of_descendants = n
        return n

    def compute_depth_of_self_and_children(self):
        if self.parent is None:
            self.depth = 1
        else:
            self.depth = self.parent.depth + 1
        for child in self.children:
            child.compute_depth_of_self_and_children()

    def append_child(self, child):
        if child in self.children:
            return
        self.children.append(child)
        child.parent = self

    def pretty_print(self, d = 0):
    #    if d > 1:
    #        return
        print ' ' * 4 * d , self.label, " : num descendants:", self.number_of_descendants
        for child in self.children:
            child.pretty_print(d + 1)

def schnyder(g, emb_dict, circular=False):
    # TODO
    pass



