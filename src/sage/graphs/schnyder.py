"""
Schnyder's Algorithm for straight-line planar embeddings

A module for computing the (x,y) coordinates for a straight-line planar
embedding of any connected planar graph with at least three vertices.  Uses
Walter Schnyder's Algorithm.

AUTHORS:

- Jonathan Bober, Emily Kirkman (2008-02-09) --  initial version

REFERENCE:

.. [1] Schnyder, Walter. Embedding Planar Graphs on the Grid.
       Proc. 1st Annual ACM-SIAM Symposium on Discrete Algorithms,
       San Francisco (1994), pp. 138-147.
"""
#*****************************************************************************
#      Copyright (C) 2008 Jonathan Bober and Emily Kirkman
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from sage.sets.set import Set
from all import DiGraph

def _triangulate(g, comb_emb):
    """
    Helper function to schnyder method for computing coordinates in the plane to
    plot a planar graph with no edge crossings.

    Given a connected graph g with at least 3 vertices and a planar combinatorial
    embedding comb_emb of g, modify g in place to form a graph whose faces are
    all triangles, and return the set of newly created edges.

    The simple way to triangulate a face is to just pick a vertex and draw
    an edge from that vertex to every other vertex in the face. Think that this
    might ultimately result in graphs that don't look very nice when we draw them
    so we have decided on a different strategy.  After handling special cases,
    we add an edge to connect every third vertex on each face.  If this edge is
    a repeat, we keep the first edge of the face and retry the process at the next
    edge in the face list.  (By removing the special cases we guarantee that this
    method will work on one of these attempts.)

    INPUT:

    - g -- the graph to triangulate
    - ``comb_emb`` -- a planar combinatorial embedding of g

    OUTPUT:

    A list of edges that are added to the graph (in place)

    EXAMPLES::

        sage: from sage.graphs.schnyder import _triangulate
        sage: g = Graph(graphs.CycleGraph(4))
        sage: g.is_planar(set_embedding=True)
        True
        sage: _triangulate(g, g._embedding)
        [(2, 0), (1, 3)]

        sage: g = graphs.PathGraph(3)
        sage: g.is_planar(set_embedding=True)
        True
        sage: _triangulate(g, g._embedding)
        [(0, 2)]
    """
    # first make sure that the graph has at least 3 vertices, and that it is connected
    if g.order() < 3:
        raise ValueError("A Graph with less than 3 vertices doesn't have any triangulation.")
    if g.is_connected() == False:
        raise NotImplementedError("_triangulate() only knows how to handle connected graphs.")

    if g.order() == 3 and len(g.edges()) == 2:             # if g is o--o--o
        vertex_list = g.vertices()
        if len(g.neighbors(vertex_list[0])) == 2:          # figure out which of the vertices already has two neighbors
            new_edge = (vertex_list[1], vertex_list[2])    # and connect the other two together.
        elif len(g.neighbors(vertex_list[1])) == 2:
            new_edge = (vertex_list[0], vertex_list[2])
        else:
            new_edge = (vertex_list[0], vertex_list[1])

        g.add_edge(new_edge)
        return [new_edge]

    # At this point we know that the graph is connected, has at least 3 vertices, and
    # that it is not the graph o--o--o. This is where the real work starts.

    faces = g.faces(comb_emb)        # We start by finding all of the faces of this embedding.

    edges_added = []                        # The list of edges that we add to the graph.
                                            # This will be returned at the end.


    for face in faces:
        new_face = []
        if len(face) < 3:
            raise RuntimeError('Triangulate method created face %s with < 3 edges.'%face)
        if len(face) == 3:
            continue # This face is already triangulated
        elif len(face) == 4:  # In this special case just add diagonal edge to square
            new_face = (face[1][1], face[0][0])
            if g.has_edge(new_face):
                new_face = (face[2][1], face[1][0])
            g.add_edge(new_face)
            edges_added.append(new_face)
        else:
            N = len(face)
            i = 0
            while i < N-1:
                new_edge = (face[i+1][1], face[i][0])  # new_edge is from third vertex in face to first
                if g.has_edge(new_edge) or new_edge[0] == new_edge[1]: # check for repeats
                    new_face.append(face[i]) # if repeated, keep first edge in face instead
                    if i == N - 2:           # if we are two from the end, found a triangle already
                        break
                    i = i + 1
                    continue

                g.add_edge(new_edge)
                edges_added.append(new_edge)
                new_face.append((new_edge[1], new_edge[0]))
                i = i + 2
            if i != N:
                new_face.append(face[-1])
            faces.append(new_face)

    return edges_added

def _normal_label(g, comb_emb, external_face):
    r"""
    Helper function to schnyder method for computing coordinates in
    the plane to plot a planar graph with no edge crossings.

    Constructs a normal labelling of a triangular graph g, given the
    planar combinatorial embedding of g and a designated external
    face.  Returns labels dictionary.  The normal label is constructed
    by first contracting the graph down to its external face, then
    expanding the graph back to the original while simultaneously
    adding angle labels.

    INPUT:

    - g -- the graph to find the normal labeling of (g must be triangulated)
    - ``comb_emb`` -- a planar combinatorial embedding of g
    - ``external_face`` -- the list of three edges in the external face of g

    OUTPUT:

    x -- tuple with entries

        x[0] = dict of dicts of normal labeling for each vertex of g and each
        adjacent neighbors u,v (u < v) of vertex:

        { vertex : { (u,v): angel_label } }

        x[1] = (v1,v2,v3) tuple of the three vertices of the external face.

    EXAMPLES::

        sage: from sage.graphs.schnyder import _triangulate, _normal_label, _realizer
        sage: g = Graph(graphs.CycleGraph(7))
        sage: g.is_planar(set_embedding=True)
        True
        sage: faces = g.faces(g._embedding)
        sage: _triangulate(g, g._embedding)
        [(2, 0), (4, 2), (6, 4), (5, 0), (3, 5), (1, 3), (4, 0), (3, 0)]
        sage: tn = _normal_label(g, g._embedding, faces[0])
        sage: _realizer(g, tn)
        ({0: [<sage.graphs.schnyder.TreeNode instance at ...>]},
         (0, 1, 2))
    """
    contracted = []
    contractible = []

    labels = {}

    external_vertices = sorted([external_face[0][0], external_face[1][0], external_face[2][0]])
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
        except Exception:
            raise RuntimeError('Contractible list is empty but graph still has %d vertices.  (Expected 3.)'%g.order())

            break
        # going to contract v
        v_neighbors = Set(g.neighbors(v))
        contracted.append( (v, v_neighbors, v_neighbors - v1_neighbors - Set([v1])) )
        g.delete_vertex(v)
        v1_neighbors -= Set([v])
        for w in v_neighbors - v1_neighbors - Set([v1]):
            # adding edge (v1, w)
            g.add_edge( (v1, w) )
        if g.order() == 3:
            break
        v1_neighbors += v_neighbors - Set([v1])
        contractible = []
        for w in g.neighbors(v1):
            if(len(v1_neighbors.intersection( Set(g.neighbors(w))))) == 2 and w not in [v1, v2, v3]:
                contractible.append(w)


    # expansion phase:

    v1, v2, v3 = g.vertices() #always in sorted order

    labels[v1] = {(v2,v3):1}
    labels[v2] = {(v1,v3):2}
    labels[v3] = {(v1,v2):3}

    while len(contracted) > 0:
        v, new_neighbors, neighbors_to_delete = contracted.pop()
        # going to add back vertex v
        labels[v] = {}

        for w in neighbors_to_delete:
            g.delete_edge((v1,w))

        if len(neighbors_to_delete) == 0:
            # we are adding v into the face new_neighbors
            w1, w2, w3 = sorted(new_neighbors)

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

        for w in new_neighbors:
            g.add_edge((v,w))

    return labels, (v1,v2,v3)

def _realizer(g, x, example=False):
    """
    Given a triangulated graph g and a normal labeling constructs the
    realizer and returns a dictionary of three trees determined by the
    realizer, each spanning all interior vertices and rooted at one of
    the three external vertices.

    A realizer is a directed graph with edge labels that span all interior
    vertices from each external vertex.  It is determined by giving direction
    to the edges that have the same angle label on both sides at a vertex.
    (Thus the direction actually points to the parent in the tree.)  The
    edge label is set as whatever the matching angle label is.  Then from
    any interior vertex, following the directed edges by label will
    give a path to each of the three external vertices.

    INPUT:

    - g -- the graph to compute the realizer of
    - x -- tuple with entries

        x[0] = dict of dicts representing a normal labeling of g.  For
        each vertex of g and each adjacent neighbors u,v (u < v) of
        vertex:  { vertex : { (u,v): angle_label } }

        x[1] = (v1, v2, v3) tuple of the three external vertices (also
        the roots of each tree)

    OUTPUT:

    - x -- tuple with entries

        x[0] = dict of lists of TreeNodes:

        { root_vertex : [ list of all TreeNodes under root_vertex ] }

        x[1] = (v1,v2,v3) tuple of the three external vertices (also the
        roots of each tree)

    EXAMPLES::

        sage: from sage.graphs.schnyder import _triangulate, _normal_label, _realizer
        sage: g = Graph(graphs.CycleGraph(7))
        sage: g.is_planar(set_embedding=True)
        True
        sage: faces = g.faces(g._embedding)
        sage: _triangulate(g, g._embedding)
        [(2, 0), (4, 2), (6, 4), (5, 0), (3, 5), (1, 3), (4, 0), (3, 0)]
        sage: tn = _normal_label(g, g._embedding, faces[0])
        sage: _realizer(g, tn)
        ({0: [<sage.graphs.schnyder.TreeNode instance at ...>]},
         (0, 1, 2))

    """
    normal_labeling, (v1, v2, v3) = x
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
                tree_nodes[v][0].append_child(tree_nodes[ones[i]][0])
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

    _compute_coordinates(realizer, (tree_nodes, (v1, v2, v3)))

    if example:
        realizer.show(talk=True, edge_labels=True)

    return tree_nodes, (v1,v2,v3)

def _compute_coordinates(g, x):
    r"""
    Given a triangulated graph g with a dict of trees given by the
    realizer and tuple of the external vertices, we compute the
    coordinates of a planar geometric embedding in the grid.

    The coordinates will be set to the ``_pos`` attribute of g.

    INPUT:

    - g -- the graph to compute the coordinates of
    - x -- tuple with entries

        x[0] = dict of tree nodes for the three trees with each external
        vertex as root:

        { root_vertex : [ list of all TreeNodes under root_vertex ] }

        x[1] = (v1, v2, v3) tuple of the three external vertices (also
        the roots of each tree)

    EXAMPLES::

        sage: from sage.graphs.schnyder import _triangulate, _normal_label, _realizer, _compute_coordinates
        sage: g = Graph(graphs.CycleGraph(7))
        sage: g.is_planar(set_embedding=True)
        True
        sage: faces = g.faces(g._embedding)
        sage: _triangulate(g, g._embedding)
        [(2, 0), (4, 2), (6, 4), (5, 0), (3, 5), (1, 3), (4, 0), (3, 0)]
        sage: tn = _normal_label(g, g._embedding, faces[0])
        sage: r = _realizer(g, tn)
        sage: _compute_coordinates(g,r)
        sage: g.get_pos()
        {0: [5, 1], 1: [0, 5], 2: [1, 0], 3: [1, 4], 4: [2, 1], 5: [2, 3], 6: [3, 2]}
    """

    tree_nodes, (v1, v2, v3) = x
    # find the roots of each tree:
    t1, t2, t3 = tree_nodes[v1][0], tree_nodes[v2][1], tree_nodes[v3][2]

    # Compute the number of descendants and depth of each node in
    # each tree.
    t1.compute_number_of_descendants()
    t2.compute_number_of_descendants()
    t3.compute_number_of_descendants()

    t1.compute_depth_of_self_and_children()
    t2.compute_depth_of_self_and_children()
    t3.compute_depth_of_self_and_children()

    coordinates = {} # the dict to pass to g.set_pos()

    # Setting coordinates for external vertices
    coordinates[t1.label] = [g.order() - 2, 1]
    coordinates[t2.label] = [0, g.order() - 2]
    coordinates[t3.label] = [1, 0]

    for v in g.vertices():
        if v not in [t1.label,t2.label,t3.label]:
            # Computing coordinates for v
            r = list((0,0,0))

            for i in [0,1,2]:
                # Computing size of region i:

                # Tracing up tree (i + 1) % 3
                p = tree_nodes[v][(i + 1) % 3]
                while p is not None:
                    q = tree_nodes[p.label][i].number_of_descendants
                    # Adding number of descendants from Tree i nodes with
                    # labels on path up tree (i + 1) % 3
                    r[i] += q
                    p = p.parent

                # Tracing up tree (i - 1) % 3
                p = tree_nodes[v][(i - 1) % 3]
                while p is not None:
                    q = tree_nodes[p.label][i].number_of_descendants
                    # Adding number of descendants from Tree i nodes with
                    # labels on path up tree (i - 1) % 3
                    r[i] += q
                    p = p.parent

                q = tree_nodes[v][i].number_of_descendants
                # Subtracting
                r[i] -= q

                # Subtracting
                q = tree_nodes[v][(i-1)%3].depth
                r[i] -= q

            if sum(r) != g.order() - 1:
                raise RuntimeError("Computing coordinates failed: vertex %s's coordinates sum to %s.  Expected %s"%(v,sum(r),g.order()-1))

            coordinates[v] = r[:-1]

    g.set_pos(coordinates) # Setting _pos attribute to store coordinates

class TreeNode():
    """
    A class to represent each node in the trees used by ``_realizer`` and
    ``_compute_coordinates`` when finding a planar geometric embedding in
    the grid.

    Each tree node is doubly linked to its parent and children.

    INPUT:

    - ``parent`` -- the parent TreeNode of ``self``
    - ``children`` -- a list of TreeNode children of ``self``
    - ``label`` -- the associated realizer vertex label

    EXAMPLES::

        sage: from sage.graphs.schnyder import TreeNode
        sage: tn = TreeNode(label=5)
        sage: tn2 = TreeNode(label=2,parent=tn)
        sage: tn3 = TreeNode(label=3)
        sage: tn.append_child(tn3)
        sage: tn.compute_number_of_descendants()
        2
        sage: tn.number_of_descendants
        2
        sage: tn3.number_of_descendants
        1
        sage: tn.compute_depth_of_self_and_children()
        sage: tn3.depth
        2
    """
    def __init__(self, parent = None, children = None, label = None):
        """
        INPUT:

        - ``parent`` -- the parent TreeNode of ``self``
        - ``children`` -- a list of TreeNode children of ``self``
        - ``label`` -- the associated realizer vertex label

        EXAMPLE::

            sage: from sage.graphs.schnyder import TreeNode
            sage: tn = TreeNode(label=5)
            sage: tn2 = TreeNode(label=2,parent=tn)
            sage: tn3 = TreeNode(label=3)
            sage: tn.append_child(tn3)
            sage: tn.compute_number_of_descendants()
            2
            sage: tn.number_of_descendants
            2
            sage: tn3.number_of_descendants
            1
            sage: tn.compute_depth_of_self_and_children()
            sage: tn3.depth
            2
        """
        if children is None:
            children = []
        self.parent = parent
        self.children = children
        self.label = label
        self.number_of_descendants = 1

    def compute_number_of_descendants(self):
        """
        Computes the number of descendants of self and all descendants.

        For each TreeNode, sets result as attribute self.number_of_descendants

        EXAMPLES::

            sage: from sage.graphs.schnyder import TreeNode
            sage: tn = TreeNode(label=5)
            sage: tn2 = TreeNode(label=2,parent=tn)
            sage: tn3 = TreeNode(label=3)
            sage: tn.append_child(tn3)
            sage: tn.compute_number_of_descendants()
            2
            sage: tn.number_of_descendants
            2
            sage: tn3.number_of_descendants
            1
            sage: tn.compute_depth_of_self_and_children()
            sage: tn3.depth
            2

        """
        n = 1
        for child in self.children:
            n += child.compute_number_of_descendants()
        self.number_of_descendants = n
        return n

    def compute_depth_of_self_and_children(self):
        """
        Computes the depth of self and all descendants.

        For each TreeNode, sets result as attribute self.depth

        EXAMPLES::

            sage: from sage.graphs.schnyder import TreeNode
            sage: tn = TreeNode(label=5)
            sage: tn2 = TreeNode(label=2,parent=tn)
            sage: tn3 = TreeNode(label=3)
            sage: tn.append_child(tn3)
            sage: tn.compute_number_of_descendants()
            2
            sage: tn.number_of_descendants
            2
            sage: tn3.number_of_descendants
            1
            sage: tn.compute_depth_of_self_and_children()
            sage: tn3.depth
            2
        """
        if self.parent is None:
            self.depth = 1
        else:
            self.depth = self.parent.depth + 1
        for child in self.children:
            child.compute_depth_of_self_and_children()

    def append_child(self, child):
        """
        Add a child to list of children.

        EXAMPLES::

            sage: from sage.graphs.schnyder import TreeNode
            sage: tn = TreeNode(label=5)
            sage: tn2 = TreeNode(label=2,parent=tn)
            sage: tn3 = TreeNode(label=3)
            sage: tn.append_child(tn3)
            sage: tn.compute_number_of_descendants()
            2
            sage: tn.number_of_descendants
            2
            sage: tn3.number_of_descendants
            1
            sage: tn.compute_depth_of_self_and_children()
            sage: tn3.depth
            2
        """
        if child in self.children:
            return
        self.children.append(child)
        child.parent = self
