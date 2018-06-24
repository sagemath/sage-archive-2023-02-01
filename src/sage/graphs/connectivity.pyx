# cython: binding=True
r"""
Connectivity related functions

This module implements the connectivity based functions for graphs and digraphs.
The methods in this module are also available as part of GenericGraph, DiGraph 
or Graph classes as aliases, and these methods can be accessed through this module
or as class methods.
Here is what the module can do:

**For both directed and undirected graphs:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`is_connected` | Test whether the (di)graph is connected.
    :meth:`connected_components` | Return the list of connected components
    :meth:`connected_components_number` | Return the number of connected components.
    :meth:`connected_components_subgraphs` | Return a list of connected components as graph objects.
    :meth:`connected_component_containing_vertex` | Return a list of the vertices connected to vertex.
    :meth:`connected_components_sizes` | Return the sizes of the connected components as a list.
    :meth:`blocks_and_cut_vertices` | Compute the blocks and cut vertices of the graph.
    :meth:`blocks_and_cuts_tree` | Compute the blocks-and-cuts tree of the graph.
    :meth:`is_cut_edge` | Return True if the input edge is a cut-edge or a bridge.
    :meth:`is_cut_vertex` | Return True if the input vertex is a cut-vertex.
    :meth:`edge_connectivity` | Return the edge connectivity of the graph.
    :meth:`vertex_connectivity` | Return the vertex connectivity of the graph.

**For DiGraph:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`is_strongly_connected` | Returns whether the current ``DiGraph`` is strongly connected.
    :meth:`strongly_connected_components_digraph` | Returns the digraph of the strongly connected components
    :meth:`strongly_connected_components_subgraphs` | Returns the strongly connected components as a list of subgraphs.
    :meth:`strongly_connected_component_containing_vertex` | Returns the strongly connected component containing a given vertex.
    :meth:`strong_articulation_points` | Return the strong articulation points of this digraph.

**For undirected graphs:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`bridges` | Returns a list of the bridges (or cut edges) of given undirected graph.

Methods
-------
"""
from __future__ import absolute_import
from sage.rings.integer import Integer

def is_connected(G):
    """
    Indicates whether the (di)graph is connected. Note that in a graph,
    path connected is equivalent to connected.

    INPUT:

    - ``G`` (generic_graph) - the input graph.

    .. SEEALSO::

        - :meth:`~Graph.is_biconnected`

    EXAMPLES::

        sage: from sage.graphs.connectivity import is_connected
        sage: G = Graph( { 0 : [1, 2], 1 : [2], 3 : [4, 5], 4 : [5] } )
        sage: is_connected(G)
        False
        sage: G.is_connected()
        False
        sage: G.add_edge(0,3)
        sage: is_connected(G)
        True
        sage: D = DiGraph( { 0 : [1, 2], 1 : [2], 3 : [4, 5], 4 : [5] } )
        sage: is_connected(D)
        False
        sage: D.add_edge(0,3)
        sage: is_connected(D)
        True
        sage: D = DiGraph({1:[0], 2:[0]})
        sage: is_connected(D)
        True

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import is_connected
        sage: is_connected('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if G.order() == 0:
        return True

    try:
        return G._backend.is_connected()
    except AttributeError:
        v = next(G.vertex_iterator())
        conn_verts = list(G.depth_first_search(v, ignore_direction=True))
        return len(conn_verts) == G.num_verts()


def connected_components(G):
    """
    Returns the list of connected components.

    Returns a list of lists of vertices, each list representing a
    connected component. The list is ordered from largest to smallest
    component.

    INPUT:

    - ``G`` (generic_graph) - the input graph.

    EXAMPLES::

        sage: from sage.graphs.connectivity import connected_components
        sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: connected_components(G)
        [[0, 1, 2, 3], [4, 5, 6]]
        sage: G.connected_components()
        [[0, 1, 2, 3], [4, 5, 6]]
        sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: connected_components(D)
        [[0, 1, 2, 3], [4, 5, 6]]

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_components
        sage: connected_components('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    seen = set()
    components = []
    for v in G:
        if v not in seen:
            c = connected_component_containing_vertex(G, v)
            seen.update(c)
            components.append(c)
    components.sort(key=lambda comp: -len(comp))
    return components


def connected_components_number(G):
    """
    Returns the number of connected components.

    INPUT:

    - ``G`` (generic_graph) - the input graph.

    EXAMPLES::

        sage: from sage.graphs.connectivity import connected_components_number
        sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: connected_components_number(G)
        2
        sage: G.connected_components_number()
        2
        sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: connected_components_number(D)
        2

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_components_number
        sage: connected_components_number('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    return len(connected_components(G))


def connected_components_subgraphs(G):
    """
    Returns a list of connected components as graph objects.

    EXAMPLES::

        sage: from sage.graphs.connectivity import connected_components_subgraphs
        sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: L = connected_components_subgraphs(G)
        sage: graphs_list.show_graphs(L)
        sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: L = connected_components_subgraphs(D)
        sage: graphs_list.show_graphs(L)
        sage: L = D.connected_components_subgraphs()
        sage: graphs_list.show_graphs(L)

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_components_subgraphs
        sage: connected_components_subgraphs('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    cc = connected_components(G)
    list = []
    for c in cc:
        list.append(G.subgraph(c, inplace=False))
    return list


def connected_component_containing_vertex(G, vertex):
    """
    Returns a list of the vertices connected to vertex.

    INPUT:

    - ``G`` (generic_graph) - the input graph.

    - ``v`` - the vertex to search for.

    EXAMPLES::

        sage: from sage.graphs.connectivity import connected_component_containing_vertex
        sage: G = Graph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: connected_component_containing_vertex(G,0)
        [0, 1, 2, 3]
        sage: G.connected_component_containing_vertex(0)
        [0, 1, 2, 3]
        sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: connected_component_containing_vertex(D,0)
        [0, 1, 2, 3]

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_component_containing_vertex
        sage: connected_component_containing_vertex('I am not a graph',0)
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    try:
        c = list(G._backend.depth_first_search(vertex, ignore_direction=True))
    except AttributeError:
        c = list(G.depth_first_search(vertex, ignore_direction=True))

    c.sort()
    return c


def connected_components_sizes(G):
    """
    Return the sizes of the connected components as a list.

    The list is sorted from largest to lower values.

    EXAMPLES::

        sage: from sage.graphs.connectivity import connected_components_sizes
        sage: for x in graphs(3):    print(connected_components_sizes(x))
        [1, 1, 1]
        [2, 1]
        [3]
        [3]
        sage: for x in graphs(3):    print(x.connected_components_sizes())
        [1, 1, 1]
        [2, 1]
        [3]
        [3]

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_components_sizes
        sage: connected_components_sizes('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    return sorted((len(cc) for cc in connected_components(G)),reverse=True)


def blocks_and_cut_vertices(G, algorithm="Tarjan_Boost"):
    """
    Computes the blocks and cut vertices of the graph.

    In the case of a digraph, this computation is done on the underlying
    graph.

    A cut vertex is one whose deletion increases the number of connected
    components. A block is a maximal induced subgraph which itself has no
    cut vertices. Two distinct blocks cannot overlap in more than a single
    cut vertex.

    INPUT:

    - ``algorithm`` -- The algorithm to use in computing the blocks
        and cut vertices of ``G``. The following algorithms are supported:

      - ``"Tarjan_Boost"`` (default) -- Tarjan's algorithm
        (Boost implementation).

      - ``"Tarjan_Sage"`` -- Tarjan's algorithm
        (Sage implementation).

    OUTPUT: ``(B, C)``, where ``B`` is a list of blocks - each is a list of
    vertices and the blocks are the corresponding induced subgraphs - and
    ``C`` is a list of cut vertices.

    ALGORITHM:

      We implement the algorithm proposed by Tarjan in [Tarjan72]_. The
      original version is recursive. We emulate the recursion using a stack.

    .. SEEALSO::

        - :meth:`blocks_and_cuts_tree`
        - :func:`sage.graphs.base.boost_graph.blocks_and_cut_vertices`
        - :meth:`~Graph.is_biconnected`
        - :meth:`~Graph.bridges`

    EXAMPLES:

    We construct a trivial example of a graph with one cut vertex::

        sage: from sage.graphs.connectivity import blocks_and_cut_vertices
        sage: rings = graphs.CycleGraph(10)
        sage: rings.merge_vertices([0, 5])
        sage: blocks_and_cut_vertices(rings)
        ([[0, 1, 4, 2, 3], [0, 6, 9, 7, 8]], [0])
        sage: rings.blocks_and_cut_vertices()
        ([[0, 1, 4, 2, 3], [0, 6, 9, 7, 8]], [0])
        sage: blocks_and_cut_vertices(rings, algorithm="Tarjan_Sage")
        ([[0, 1, 2, 3, 4], [0, 6, 7, 8, 9]], [0])

    The Petersen graph is biconnected, hence has no cut vertices::

        sage: blocks_and_cut_vertices(graphs.PetersenGraph())
        ([[0, 1, 4, 5, 2, 6, 3, 7, 8, 9]], [])

    Decomposing paths to pairs::

        sage: g = graphs.PathGraph(4) + graphs.PathGraph(5)
        sage: blocks_and_cut_vertices(g)
        ([[2, 3], [1, 2], [0, 1], [7, 8], [6, 7], [5, 6], [4, 5]], [1, 2, 5, 6, 7])

    A disconnected graph::

        sage: g = Graph({1:{2:28, 3:10}, 2:{1:10, 3:16}, 4:{}, 5:{6:3, 7:10, 8:4}})
        sage: blocks_and_cut_vertices(g)
        ([[1, 2, 3], [5, 6], [5, 7], [5, 8], [4]], [5])

    TESTS::

        sage: blocks_and_cut_vertices(Graph(0))
        ([], [])
        sage: blocks_and_cut_vertices(Graph(1))
        ([[0]], [])
        sage: blocks_and_cut_vertices(Graph(2))
        ([[0], [1]], [])

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_components_sizes
        sage: connected_components_sizes('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if (algorithm=="Tarjan_Boost"):
        from sage.graphs.base.boost_graph import blocks_and_cut_vertices
        return blocks_and_cut_vertices(G)

    if (algorithm!="Tarjan_Sage"):
        raise NotImplementedError("Blocks and cut vertices algorithm '%s' is not implemented." % algorithm)

    # If algorithm is "Tarjan_Sage"
    blocks = []
    cut_vertices = set()

    # We iterate over all vertices to ensure that we visit each connected
    # component of the graph
    seen = set()
    for start in G.vertex_iterator():
        if start in seen:
            continue

        # Special case of an isolated vertex
        if not G.degree(start):
            blocks.append([start])
            seen.add(start)
            continue

        # Each vertex is numbered with an integer from 1...|V(G)|,
        # corresponding to the order in which it is discovered during the
        # DFS.
        number = {}
        num = 1

        # Associates to each vertex v the smallest number of a vertex that
        # can be reached from v in the orientation of the graph that the
        # algorithm creates.
        low_point = {}

        # Associates to each vertex an iterator over its neighbors
        neighbors = {}

        stack = [start]
        edge_stack = []
        start_already_seen = False

        while stack:
            v = stack[-1]
            seen.add(v)

            # The first time we meet v
            if not v in number:
                # We number the vertices in the order they are reached
                # during DFS
                number[v] = num
                neighbors[v] = G.neighbor_iterator(v)
                low_point[v] = num
                num += 1

            try:
                # We consider the next of its neighbors
                w = next(neighbors[v])

                # If we never met w before, we remember the direction of
                # edge vw, and add w to the stack.
                if not w in number:
                    edge_stack.append( (v,w) )
                    stack.append(w)

                # If w is an ancestor of v in the DFS tree, we remember the
                # direction of edge vw
                elif number[w]<number[v]:
                    edge_stack.append( (v,w) )
                    low_point[v] = min(low_point[v], number[w])

            # We went through all of v's neighbors
            except StopIteration:
                # We trackback, so w takes the value of v and we pop the
                # stack
                w = stack.pop()

                # Test termination of the algorithm
                if not stack:
                    break

                v = stack[-1]

                # Propagating the information : low_point[v] indicates the
                # smallest vertex (the vertex x with smallest number[x])
                # that can be reached from v
                low_point[v] = min(low_point[v], low_point[w])

                # The situation in which there is no path from w to an
                # ancestor of v : we have identified a new biconnected
                # component
                if low_point[w] >= number[v]:
                    new_block = set()
                    nw = number[w]
                    u1,u2 = edge_stack.pop()
                    while number[u1] >= nw:
                        new_block.add(u1)
                        u1,u2 = edge_stack.pop()
                    new_block.add(u1)
                    blocks.append(sorted(list(new_block)))

                    # We update the set of cut vertices.
                    #
                    # If v is start, then we add it only if it belongs to
                    # several blocks.
                    if (not v is start) or start_already_seen:
                        cut_vertices.add(v)
                    else:
                        start_already_seen = True

    return blocks,sorted(list(cut_vertices))


def blocks_and_cuts_tree(G):
    """
    Returns the blocks-and-cuts tree of ``self``.

    This new graph has two different kinds of vertices, some representing
    the blocks (type B) and some other the cut vertices of the graph
    ``self`` (type C).

    There is an edge between a vertex `u` of type B and a vertex `v` of type
    C if the cut-vertex corresponding to `v` is in the block corresponding
    to `u`.

    The resulting graph is a tree, with the additional characteristic
    property that the distance between two leaves is even. When ``self`` is
    not connected, the resulting graph is a forest.

    When ``self`` is biconnected, the tree is reduced to a single node of
    type `B`.

    We referred to [HarPri]_ and [Gallai]_ for blocks and cuts tree.

    .. SEEALSO::

        - :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cut_vertices`
        - :meth:`~Graph.is_biconnected`

    EXAMPLES::

        sage: from sage.graphs.connectivity import blocks_and_cuts_tree
        sage: T = blocks_and_cuts_tree(graphs.KrackhardtKiteGraph()); T
        Graph on 5 vertices
        sage: T.is_isomorphic(graphs.PathGraph(5))
        True
        sage: from sage.graphs.connectivity import blocks_and_cuts_tree
        sage: T = graphs.KrackhardtKiteGraph().blocks_and_cuts_tree(); T
        Graph on 5 vertices

    The distance between two leaves is even::

        sage: T = blocks_and_cuts_tree(graphs.RandomTree(40))
        sage: T.is_tree()
        True
        sage: leaves = [v for v in T if T.degree(v) == 1]
        sage: all(T.distance(u,v) % 2 == 0 for u in leaves for v in leaves)
        True

    The tree of a biconnected graph has a single vertex, of type `B`::

        sage: T = blocks_and_cuts_tree(graphs.PetersenGraph())
        sage: T.vertices()
        [('B', (0, 1, 4, 5, 2, 6, 3, 7, 8, 9))]

    TESTS:

    When ``self`` is not connected, the resulting graph is a forest (:trac:`24163`)::

        sage: from sage.graphs.connectivity import blocks_and_cuts_tree
        sage: T = blocks_and_cuts_tree(Graph(2))
        sage: T.is_forest()
        True

    If ``G`` is not a Sage graph, an error is raised::

        sage: blocks_and_cuts_tree('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    from sage.graphs.graph import Graph
    B, C = G.blocks_and_cut_vertices()
    B = map(tuple, B)
    g = Graph()
    for bloc in B:
        g.add_vertex(('B', bloc))
        for c in bloc:
            if c in C:
                g.add_edge(('B', bloc), ('C', c))
    return g

def is_cut_edge(G, u, v=None, label=None):
    """
    Returns True if the input edge is a cut-edge or a bridge.

    A cut edge (or bridge) is an edge that when removed increases
    the number of connected components.  This function works with
    simple graphs as well as graphs with loops and multiedges.  In
    a digraph, a cut edge is an edge that when removed increases
    the number of (weakly) connected components.

    INPUT: The following forms are accepted

    - is_cut_edge(G, 1, 2 )

    - is_cut_edge(G, (1, 2) )

    - is_cut_edge(G, 1, 2, 'label' )

    - is_cut_edge(G, (1, 2, 'label') )

    OUTPUT:

    - Returns True if (u,v) is a cut edge, False otherwise

    EXAMPLES::

        sage: from sage.graphs.connectivity import is_cut_edge
        sage: G = graphs.CompleteGraph(4)
        sage: is_cut_edge(G,0,2)
        False
        sage: G.is_cut_edge(0,2)
        False

        sage: G = graphs.CompleteGraph(4)
        sage: G.add_edge((0,5,'silly'))
        sage: is_cut_edge(G,(0,5,'silly'))
        True

        sage: G = Graph([[0,1],[0,2],[3,4],[4,5],[3,5]])
        sage: is_cut_edge(G,(0,1))
        True

        sage: G = Graph([[0,1],[0,2],[1,1]], loops = True)
        sage: is_cut_edge(G,(1,1))
        False

        sage: G = digraphs.Circuit(5)
        sage: is_cut_edge(G,(0,1))
        False

        sage: G = graphs.CompleteGraph(6)
        sage: is_cut_edge(G,(0,7))
        Traceback (most recent call last):
        ...
        ValueError: edge not in graph

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: is_cut_edge('I am not a graph',0)
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if label is None:
        if v is None:
            try:
                u, v, label = u
            except ValueError:
                u, v = u
                label = None

    if not G.has_edge(u,v):
        raise ValueError('edge not in graph')

    # If edge (u,v) is a pending edge, it is also a cut-edge
    if G.degree(u) == 1 or G.degree(v) == 1:
        return True
    elif G.allows_multiple_edges():
        # If we have two or more edges between u and v, it is not a cut-edge
        if len([(uu,vv) for uu,vv,ll in G.edges_incident(u) if uu == v or vv == v]) > 1:
            return False

    g = G.copy(immutable=False) if G.is_immutable() else G
    g.delete_edge(u,v,label)
    if g.is_directed():
        # (u,v) is a cut-edge if u is not in the connected
        # component containing v of self-(u,v)
        sol = not u in connected_component_containing_vertex(g,v)
    else:
        # (u,v) is a cut-edge if there is no path from u to v in
        # self-(u,v)
        sol = not g.distance(u,v) < g.order()

    g.add_edge(u,v,label)
    return sol


def is_cut_vertex(G, u, weak=False):
    r"""
    Returns True if the input vertex is a cut-vertex.

    A vertex is a cut-vertex if its removal from the (di)graph increases the
    number of (strongly) connected components. Isolated vertices or leafs
    are not cut-vertices. This function works with simple graphs as well as
    graphs with loops and multiple edges.

    INPUT:

    - ``G`` (generic_graph) - the input graph.

    - ``u`` -- a vertex

    - ``weak`` -- (default: ``False``) boolean set to `True` if the
      connectivity of directed graphs is to be taken in the weak sense, that
      is ignoring edges orientations.

    OUTPUT:

    Returns True if ``u`` is a cut-vertex, and False otherwise.

    EXAMPLES:

    Giving a LollipopGraph(4,2), that is a complete graph with 4 vertices with a pending edge::

        sage: from sage.graphs.connectivity import is_cut_vertex
        sage: G = graphs.LollipopGraph(4,2)
        sage: is_cut_vertex(G,0)
        False
        sage: is_cut_vertex(G,3)
        True
        sage: G.is_cut_vertex(3)
        True

    Comparing the weak and strong connectivity of a digraph::

        sage: from sage.graphs.connectivity import is_strongly_connected
        sage: D = digraphs.Circuit(6)
        sage: is_strongly_connected(D)
        True
        sage: is_cut_vertex(D,2)
        True
        sage: is_cut_vertex(D, 2, weak=True)
        False

    Giving a vertex that is not in the graph::

        sage: G = graphs.CompleteGraph(6)
        sage: is_cut_vertex(G,7)
        Traceback (most recent call last):
        ...
        ValueError: The input vertex is not in the vertex set.

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: is_cut_vertex('I am not a graph',0)
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if not u in G:
        raise ValueError('The input vertex is not in the vertex set.')

    # Initialization
    if not G.is_directed() or weak:
        # Weak connectivity

        if G.degree(u) < 2:
            # An isolated or a leaf vertex is not a cut vertex
            return False

        neighbors_func = [G.neighbor_iterator]
        start = next(G.neighbor_iterator(u))
        CC = set(G.vertex_iterator())

    else:
        # Strong connectivity for digraphs

        if G.out_degree(u) == 0 or G.in_degree(u) == 0:
            # A vertex without in or out neighbors is not a cut vertex
            return False

        # We consider only the strongly connected component containing u
        CC = set(strongly_connected_component_containing_vertex(G, u))

        # We perform two DFS starting from an out neighbor of u and avoiding
        # u. The first DFS follows the edges directions, and the second is
        # in the reverse order. If both allow to reach all neighbors of u,
        # then u is not a cut vertex
        neighbors_func = [G.neighbor_out_iterator, G.neighbor_in_iterator]
        start = next(G.neighbor_out_iterator(u))

    CC.discard(u)
    CC.discard(start)
    for neighbors in neighbors_func:

        # We perform a DFS starting from a neighbor of u and avoiding u
        queue = [start]
        seen = set(queue)
        targets = set(G.neighbor_iterator(u))&CC
        targets.discard(start)
        while queue and targets:
            v = queue.pop()
            for w in neighbors(v):
                if not w in seen and w in CC:
                    seen.add(w)
                    queue.append(w)
                    targets.discard(w)

        # If some neighbors cannot be reached, u is a cut vertex.
        if targets:
            return True

    return False


def edge_connectivity(G,
                      value_only = True,
                      implementation = None,
                      use_edge_labels = False,
                      vertices = False,
                      solver = None,
                      verbose = 0):
    r"""
    Returns the edge connectivity of the graph.

    For more information, see the
    `Wikipedia article on connectivity
    <http://en.wikipedia.org/wiki/Connectivity_(graph_theory)>`_.

    .. NOTE::

        When the graph is a directed graph, this method actually computes
        the *strong* connectivity, (i.e. a directed graph is strongly
        `k`-connected if there are `k` disjoint paths between any two
        vertices `u, v`). If you do not want to consider strong
        connectivity, the best is probably to convert your ``DiGraph``
        object to a ``Graph`` object, and compute the connectivity of this
        other graph.

    INPUT:

    - ``G`` (generic_graph) - the input graph.

    - ``value_only`` -- boolean (default: ``True``)

      - When set to ``True`` (default), only the value is returned.

      - When set to ``False``, both the value and a minimum edge cut
        are returned.

    - ``implementation`` -- selects an implementation:

      - When set to ``None`` (default): selects the best implementation
        available.

      - When set to ``"boost"``, we use the Boost graph library (which is
        much more efficient). It is not available when ``edge_labels=True``,
        and it is unreliable for directed graphs (see :trac:`18753`).

      - When set to ``"Sage"``, we use Sage's implementation.

    - ``use_edge_labels`` -- boolean (default: ``False``)

      - When set to ``True``, computes a weighted minimum cut
        where each edge has a weight defined by its label. (If
        an edge has no label, `1` is assumed.). Implies
        ``boost`` = ``False``.

      - When set to ``False``, each edge has weight `1`.

    - ``vertices`` -- boolean (default: ``False``)

      - When set to ``True``, also returns the two sets of
        vertices that are disconnected by the cut. Implies
        ``value_only=False``.

    - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
      solver to be used (ignored if ``implementation='boost'``). If set to
      ``None``, the default one is used. For more information on LP solvers
      and which default solver is used, see the method
      :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
      of the class
      :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: ``0``). Sets the level of
      verbosity. Set to 0 by default, which means quiet.

    EXAMPLES:

    A basic application on the PappusGraph::

        sage: from sage.graphs.connectivity import edge_connectivity
        sage: g = graphs.PappusGraph()
        sage: edge_connectivity(g)
        3
        sage: g.edge_connectivity()
        3

    The edge connectivity of a complete graph ( and of a random graph )
    is its minimum degree, and one of the two parts of the bipartition
    is reduced to only one vertex. The cutedges isomorphic to a
    Star graph::

        sage: g = graphs.CompleteGraph(5)
        sage: [ value, edges, [ setA, setB ]] = edge_connectivity(g,vertices=True)
        sage: value
        4
        sage: len(setA) == 1 or len(setB) == 1
        True
        sage: cut = Graph()
        sage: cut.add_edges(edges)
        sage: cut.is_isomorphic(graphs.StarGraph(4))
        True

    Even if obviously in any graph we know that the edge connectivity
    is less than the minimum degree of the graph::

        sage: g = graphs.RandomGNP(10,.3)
        sage: min(g.degree()) >= edge_connectivity(g)
        True

    If we build a tree then assign to its edges a random value, the
    minimum cut will be the edge with minimum value::

        sage: g = graphs.RandomGNP(15,.5)
        sage: tree = Graph()
        sage: tree.add_edges(g.min_spanning_tree())
        sage: for u,v in tree.edge_iterator(labels=None):
        ....:      tree.set_edge_label(u,v,random())
        sage: minimum = min([l for u,v,l in tree.edge_iterator()])
        sage: [value, [(u,v,l)]] = edge_connectivity(tree, value_only=False, use_edge_labels=True)
        sage: l == minimum
        True

    When ``value_only = True`` and ``implementation="sage"``, this function is
    optimized for small connectivity values and does not need to build a
    linear program.

    It is the case for graphs which are not connected ::

        sage: g = 2 * graphs.PetersenGraph()
        sage: edge_connectivity(g, implementation="sage")
        0.0

    For directed graphs, the strong connectivity is tested
    through the dedicated function ::

        sage: g = digraphs.ButterflyGraph(3)
        sage: edge_connectivity(g, implementation="sage")
        0.0

    We check that the result with Boost is the same as the result without
    Boost ::

        sage: g = graphs.RandomGNP(15,.3)
        sage: edge_connectivity(g) == edge_connectivity(g, implementation="sage")
        True

    Boost interface also works with directed graphs ::

        sage: edge_connectivity(digraphs.Circuit(10), implementation = "boost", vertices = True)
        [1, [(0, 1)], [{0}, {1, 2, 3, 4, 5, 6, 7, 8, 9}]]

    However, the Boost algorithm is not reliable if the input is directed
    (see :trac:`18753`)::

        sage: g = digraphs.Path(3)
        sage: edge_connectivity(g)
        0.0
        sage: edge_connectivity(g, implementation="boost")
        1
        sage: g.add_edge(1,0)
        sage: edge_connectivity(g)
        0.0
        sage: edge_connectivity(g, implementation="boost")
        0

    TESTS:

    Checking that the two implementations agree::

        sage: for i in range(10):
        ....:     g = graphs.RandomGNP(30,0.3)
        ....:     e1 = edge_connectivity(g, implementation="boost")
        ....:     e2 = edge_connectivity(g, implementation="sage")
        ....:     assert (e1 == e2)

    Disconnected graphs and ``vertices=True``::

        sage: g = graphs.PetersenGraph()
        sage: edge_connectivity((2*g), vertices=True)
        [0, [], [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]]]

    If ``G`` is not a Sage graph, an error is raised::

        sage: edge_connectivity('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    G._scream_if_not_simple(allow_loops=True)
    g=G

    if vertices:
        value_only=False

    if implementation is None:
        if use_edge_labels or g.is_directed():
            implementation = "sage"
        else:
            implementation = "boost"

    implementation = implementation.lower()
    if implementation not in ["boost", "sage"]:
        raise ValueError("'implementation' must be set to 'boost', 'sage' or None.")
    elif implementation=="boost" and use_edge_labels:
        raise ValueError("The Boost implementation is currently not able to handle edge labels")

    # Otherwise, an error is created
    if g.num_edges() == 0 or g.num_verts() == 0:
        if value_only:
            return 0
        elif vertices:
            return [0,[],[[],[]]]
        else:
            return [0,[]]

    if implementation == "boost":
        from sage.graphs.base.boost_graph import edge_connectivity

        [obj, edges] = edge_connectivity(g)

        if value_only:
            return obj

        val = [obj, edges]

        if vertices and not obj == 0:
            H = G.copy()
            H.delete_edges(edges)

            if H.is_directed():
                a = set(H.breadth_first_search([x for x,y in edges]))
                b = set(H).difference(a)
                val.append([a,b])
            else:
                val.append(connected_components(H))
        elif vertices:
            val.append(connected_components(G))

        return val

    if use_edge_labels:
        from sage.rings.real_mpfr import RR
        weight=lambda x: x if x in RR else 1
    else:
        weight=lambda x: 1


    # Better methods for small connectivity tests,
    # when one is not interested in cuts...
    if value_only and not use_edge_labels:

        if G.is_directed():
            if not is_strongly_connected(G):
                return 0.0

        else:
            if not is_connected(G):
                return 0.0

            h = G.strong_orientation()
            if not is_strongly_connected(h):
                return 1.0


    if g.is_directed():
        reorder_edge = lambda x,y : (x,y)
    else:
        reorder_edge = lambda x,y : (x,y) if x<= y else (y,x)

    from sage.numerical.mip import MixedIntegerLinearProgram

    p = MixedIntegerLinearProgram(maximization=False, solver=solver)

    in_set = p.new_variable(binary = True)
    in_cut = p.new_variable(binary = True)

    # A vertex has to be in some set
    for v in g:
        p.add_constraint(in_set[0,v]+in_set[1,v],max=1,min=1)

    # There is no empty set
    p.add_constraint(p.sum([in_set[1,v] for v in g]),min=1)
    p.add_constraint(p.sum([in_set[0,v] for v in g]),min=1)

    if g.is_directed():
        # There is no edge from set 0 to set 1 which
        # is not in the cut
        for (u,v) in g.edge_iterator(labels=None):
            p.add_constraint(in_set[0,u] + in_set[1,v] - in_cut[(u,v)], max = 1)
    else:

        # Two adjacent vertices are in different sets if and only if
        # the edge between them is in the cut
        for (u,v) in g.edge_iterator(labels=None):
            p.add_constraint(in_set[0,u]+in_set[1,v]-in_cut[reorder_edge(u,v)],max=1)
            p.add_constraint(in_set[1,u]+in_set[0,v]-in_cut[reorder_edge(u,v)],max=1)

    p.set_objective(p.sum([weight(l ) * in_cut[reorder_edge(u,v)] for (u,v,l) in g.edge_iterator()]))

    obj = p.solve(objective_only=value_only, log=verbose)

    if use_edge_labels is False:
        obj = Integer(round(obj))

    if value_only:
        return obj

    else:
        val = [obj]

        in_cut = p.get_values(in_cut)
        in_set = p.get_values(in_set)

        edges = []
        for (u,v,l) in g.edge_iterator():
            if in_cut[reorder_edge(u,v)] == 1:
                edges.append((u,v,l))

        val.append(edges)

        if vertices:
            a = []
            b = []
            for v in g:
                if in_set[0,v] == 1:
                    a.append(v)
                else:
                    b.append(v)
            val.append([a,b])

        return val

def vertex_connectivity(G, value_only=True, sets=False, k=None, solver=None, verbose=0):
    r"""
    Return the vertex connectivity of the graph.

    For more information, see :wikipedia:`Connectivity_(graph_theory)` and
    :wikipedia:`K-vertex-connected_graph`.

    .. NOTE::

        * When the graph is directed, this method actually computes the
          *strong* connectivity, (i.e. a directed graph is strongly
          `k`-connected if there are `k` vertex disjoint paths between any
          two vertices `u, v`). If you do not want to consider strong
          connectivity, the best is probably to convert your ``DiGraph``
          object to a ``Graph`` object, and compute the connectivity of this
          other graph.

        * By convention, a complete graph on `n` vertices is `n-1`
          connected. In this case, no certificate can be given as there is
          no pair of vertices split by a cut of order `k-1`. For this
          reason, the certificates returned in this situation are empty.

    INPUT:

    - ``G`` (generic_graph) - the input graph.

    - ``value_only`` -- boolean (default: ``True``)

      - When set to ``True`` (default), only the value is returned.

      - When set to ``False`` , both the value and a minimum vertex cut are
        returned.

    - ``sets`` -- boolean (default: ``False``)

      - When set to ``True``, also returns the two sets of vertices that
        are disconnected by the cut.  Implies ``value_only=False``

    - ``k`` -- integer (default: ``None``) When specified, check if the
      vertex connectivity of the (di)graph is larger or equal to `k`. The
      method thus outputs a boolean only.

    - ``solver`` -- (default: ``None``) Specify a Linear Program (LP) solver
      to be used. If set to ``None``, the default one is used. For more
      information on LP solvers, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.  Use method
      :meth:`sage.numerical.backends.generic_backend.default_mip_solver` to
      know which default solver is used or to set the default solver.

    - ``verbose`` -- integer (default: ``0``). Sets the level of
      verbosity. Set to 0 by default, which means quiet.

    EXAMPLES:

    A basic application on a ``PappusGraph``::

       sage: from sage.graphs.connectivity import vertex_connectivity
       sage: g=graphs.PappusGraph()
       sage: vertex_connectivity(g)
       3
       sage: g.vertex_connectivity()
       3

    In a grid, the vertex connectivity is equal to the minimum degree, in
    which case one of the two sets is of cardinality `1`::

       sage: g = graphs.GridGraph([ 3,3 ])
       sage: [value, cut, [ setA, setB ]] = vertex_connectivity(g, sets=True)
       sage: len(setA) == 1 or len(setB) == 1
       True

    A vertex cut in a tree is any internal vertex::

       sage: tree = graphs.RandomTree(15)
       sage: val, [cut_vertex] = vertex_connectivity(tree, value_only=False)
       sage: tree.degree(cut_vertex) > 1
       True

    When ``value_only = True``, this function is optimized for small
    connectivity values and does not need to build a linear program.

    It is the case for connected graphs which are not connected::

       sage: g = 2 * graphs.PetersenGraph()
       sage: vertex_connectivity(g)
       0

    Or if they are just 1-connected::

       sage: g = graphs.PathGraph(10)
       sage: vertex_connectivity(g)
       1

    For directed graphs, the strong connectivity is tested
    through the dedicated function::

       sage: g = digraphs.ButterflyGraph(3)
       sage: vertex_connectivity(g)
       0

    A complete graph on `10` vertices is `9`-connected::

       sage: g = graphs.CompleteGraph(10)
       sage: vertex_connectivity(g)
       9

    A complete digraph on `10` vertices is `9`-connected::

       sage: g = DiGraph(graphs.CompleteGraph(10))
       sage: vertex_connectivity(g)
       9

    When parameter ``k`` is set, we only check for the existence of a
    vertex cut of order at least ``k``::

       sage: g = graphs.PappusGraph()
       sage: vertex_connectivity(g, k=3)
       True
       sage: vertex_connectivity(g, k=4)
       False

    TESTS:

    Giving negative value to parameter ``k``::

       sage: g = graphs.PappusGraph()
       sage: vertex_connectivity(g, k=-1)
       Traceback (most recent call last):
       ...
       ValueError: parameter k must be strictly positive

    The empty graph has vertex connectivity 0, is considered connected but
    not biconnected. The empty digraph is considered strongly connected::

       sage: from sage.graphs.connectivity import is_strongly_connected
       sage: from sage.graphs.connectivity import is_connected
       sage: empty = Graph()
       sage: vertex_connectivity(empty)
       0
       sage: vertex_connectivity(empty, k=1) == is_connected(empty)
       True
       sage: vertex_connectivity(Graph(), k=2) == empty.is_biconnected()
       True
       sage: vertex_connectivity(DiGraph(), k=1) == is_strongly_connected(DiGraph())
       True

    If ``G`` is not a Sage graph, an error is raised::

        sage: vertex_connectivity('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph

    Complete Graph with loops or multiple edges (:trac:`25589`)::

        sage: G = Graph([(0, 1), (0, 1)], multiedges=True)
        sage: G.vertex_connectivity()
        1
        sage: G = graphs.CompleteGraph(4)
        sage: G.allow_loops(True)
        sage: G.add_edge(0, 0)
        sage: G.vertex_connectivity(value_only=False, verbose=1)
        (3, [])
        sage: G.allow_multiple_edges(True)
        sage: G.add_edge(0, 1)
        sage: G.vertex_connectivity(value_only=False, verbose=1)
        (3, [])
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    g = G

    if k is not None:
        if k < 1:
            raise ValueError("parameter k must be strictly positive")
        if g.order() == 0:
            # We follow the convention of is_connected, is_biconnected and
            # is_strongly_connected
            return k == 1
        if (g.is_directed() and k > min(min(g.in_degree()), min(g.out_degree()))) \
           or (not g.is_directed() and (k > min(g.degree()))):
            return False
        value_only = True
        sets = False

    elif sets:
        value_only = False

    # When the graph is complete, the MILP below is infeasible.
    if (g.is_clique(directed_clique=g.is_directed()) \
        or (not g.is_directed() and g.to_simple().is_clique())):
        if k is not None:
            return g.order() > k
        if value_only:
            return max(g.order()-1, 0)
        elif not sets:
            return max(g.order()-1, 0), []
        else:
            return max(g.order()-1, 0), [], [[], []]

    if value_only:
        if G.is_directed():
            if not is_strongly_connected(G):
                return 0 if k is None else False

        else:
            if not is_connected(G):
                return 0 if k is None else False

            if len(G.blocks_and_cut_vertices()[0]) > 1:
                return 1 if k is None else (k == 1)

        if k == 1:
            # We know that the (di)graph is (strongly) connected
            return True

    from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

    p = MixedIntegerLinearProgram(maximization=False, solver=solver)

    # Sets 0 and 2 are "real" sets while set 1 represents the cut
    in_set = p.new_variable(binary=True)

    # A vertex has to be in some set
    for v in g:
        p.add_constraint(in_set[0, v] + in_set[1, v] + in_set[2, v], max=1, min=1)

    # There is no empty set
    p.add_constraint(p.sum(in_set[0, v] for v in g), min=1)
    p.add_constraint(p.sum(in_set[2, v] for v in g), min=1)

    if g.is_directed():
        # There is no edge from set 0 to set 1 which is not in the cut
        for u, v in g.edge_iterator(labels=None):
            p.add_constraint(in_set[0, u] + in_set[2, v], max=1)
    else:
        # Two adjacent vertices are in different sets if and only if
        # the edge between them is in the cut
        for u, v in g.edge_iterator(labels=None):
            p.add_constraint(in_set[0, u] + in_set[2, v], max=1)
            p.add_constraint(in_set[2, u] + in_set[0, v], max=1)

    if k is not None:
        # To check if the vertex connectivity is at least k, we check if
        # there exists a cut of order at most k-1. If the ILP is infeasible,
        # the vertex connectivity is >= k.
        p.add_constraint(p.sum(in_set[1, v] for v in g) <= k-1)
        try:
            p.solve(objective_only=True, log=verbose)
            return False
        except MIPSolverException:
            return True

    else:
        p.set_objective(p.sum(in_set[1, v] for v in g))

    if value_only:
        return Integer(round(p.solve(objective_only=True, log=verbose)))

    val = Integer(round(p.solve(log=verbose)))

    in_set = p.get_values(in_set)

    cut = []
    a = []
    b = []

    for v in g:
        if in_set[0, v]:
            a.append(v)
        elif in_set[1, v]:
            cut.append(v)
        else:
            b.append(v)

    if sets:
        return val, cut, [a, b]

    return val, cut


def is_strongly_connected(G):
    r"""
    Returns whether the current ``DiGraph`` is strongly connected.

    EXAMPLES:

    The circuit is obviously strongly connected ::

        sage: from sage.graphs.connectivity import is_strongly_connected
        sage: g = digraphs.Circuit(5)
        sage: is_strongly_connected(g)
        True
        sage: g.is_strongly_connected()
        True

    But a transitive triangle is not::

        sage: g = DiGraph({ 0 : [1,2], 1 : [2]})
        sage: is_strongly_connected(g)
        False

    TESTS:

    If ``G`` is not a Sage DiGraph, an error is raised::

        sage: is_strongly_connected('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage DiGraph
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise TypeError("the input must be a Sage DiGraph")

    if G.order()==1:
        return True

    try:
        return G._backend.is_strongly_connected()

    except AttributeError:
        return len(G.strongly_connected_components()) == 1


def strongly_connected_components_digraph(G, keep_labels = False):
    r"""
    Returns the digraph of the strongly connected components

    INPUT:

    - ``G`` (DiGraph) - the input graph.

    - ``keep_labels`` -- boolean (default: False)

    The digraph of the strongly connected components of a graph `G` has
    a vertex per strongly connected component included in `G`. There
    is an edge from a component `C_1` to a component `C_2` if there is
    an edge from one to the other in `G`.

    EXAMPLES:

    Such a digraph is always acyclic ::

        sage: from sage.graphs.connectivity import strongly_connected_components_digraph
        sage: g = digraphs.RandomDirectedGNP(15,.1)
        sage: scc_digraph = strongly_connected_components_digraph(g)
        sage: scc_digraph.is_directed_acyclic()
        True
        sage: scc_digraph = g.strongly_connected_components_digraph()
        sage: scc_digraph.is_directed_acyclic()
        True

    The vertices of the digraph of strongly connected components are
    exactly the strongly connected components::

        sage: g = digraphs.ButterflyGraph(2)
        sage: scc_digraph = strongly_connected_components_digraph(g)
        sage: g.is_directed_acyclic()
        True
        sage: all([ Set(scc) in scc_digraph.vertices() for scc in g.strongly_connected_components()])
        True

    The following digraph has three strongly connected components,
    and the digraph of those is a chain::

        sage: g = DiGraph({0:{1:"01", 2: "02", 3: "03"}, 1: {2: "12"}, 2:{1: "21", 3: "23"}})
        sage: scc_digraph = strongly_connected_components_digraph(g)
        sage: scc_digraph.vertices()
        [{0}, {3}, {1, 2}]
        sage: scc_digraph.edges()
        [({0}, {1, 2}, None), ({0}, {3}, None), ({1, 2}, {3}, None)]

    By default, the labels are discarded, and the result has no
    loops nor multiple edges. If ``keep_labels`` is ``True``, then
    the labels are kept, and the result is a multi digraph,
    possibly with multiple edges and loops. However, edges in the
    result with same source, target, and label are not duplicated
    (see the edges from 0 to the strongly connected component
    `\{1,2\}` below)::

        sage: g = DiGraph({0:{1:"0-12", 2: "0-12", 3: "0-3"}, 1: {2: "1-2", 3: "1-3"}, 2:{1: "2-1", 3: "2-3"}})
        sage: scc_digraph = strongly_connected_components_digraph(g, keep_labels = True)
        sage: scc_digraph.vertices()
        [{0}, {3}, {1, 2}]
        sage: scc_digraph.edges()
        [({0}, {1, 2}, '0-12'),
         ({0}, {3}, '0-3'),
         ({1, 2}, {1, 2}, '1-2'),
         ({1, 2}, {1, 2}, '2-1'),
         ({1, 2}, {3}, '1-3'),
         ({1, 2}, {3}, '2-3')]

    TESTS:

    If ``G`` is not a Sage DiGraph, an error is raised::

        sage: strongly_connected_components_digraph('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage DiGraph
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise TypeError("the input must be a Sage DiGraph")

    from sage.sets.set import Set

    scc = G.strongly_connected_components()
    scc_set = [Set(_) for _ in scc]

    d = {}
    for i,c in enumerate(scc):
        for v in c:
            d[v] = i

    if keep_labels:
        g = DiGraph(multiedges=True, loops=True)
        g.add_vertices(range(len(scc)))

        g.add_edges( set((d[u], d[v], label) for (u,v,label) in G.edges() ) )
        g.relabel(scc_set, inplace=True)

    else:
        g = DiGraph(multiedges=False, loops=False)
        g.add_vertices(range(len(scc)))

        g.add_edges(((d[u], d[v]) for u, v in G.edges(labels=False)), loops=False)
        g.relabel(scc_set, inplace=True)

    return g


def strongly_connected_components_subgraphs(G):
    r"""
    Returns the strongly connected components as a list of subgraphs.

    EXAMPLES:

    In the symmetric digraph of a graph, the strongly connected components are the connected
    components::

        sage: from sage.graphs.connectivity import strongly_connected_components_subgraphs
        sage: g = graphs.PetersenGraph()
        sage: d = DiGraph(g)
        sage: strongly_connected_components_subgraphs(d)
        [Subgraph of (Petersen graph): Digraph on 10 vertices]
        sage: d.strongly_connected_components_subgraphs()
        [Subgraph of (Petersen graph): Digraph on 10 vertices]

    TESTS:

    If ``G`` is not a Sage DiGraph, an error is raised::

        sage: strongly_connected_components_subgraphs('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage DiGraph
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise TypeError("the input must be a Sage DiGraph")

    return [G.subgraph(_) for _ in G.strongly_connected_components()]


def strongly_connected_component_containing_vertex(G, v):
    """
    Returns the strongly connected component containing a given vertex

    INPUT:

    - ``G`` (DiGraph) - the input graph.

    - ``v`` -- a vertex

    EXAMPLES:

    In the symmetric digraph of a graph, the strongly connected components are the connected
    components::

        sage: from sage.graphs.connectivity import strongly_connected_component_containing_vertex
        sage: g = graphs.PetersenGraph()
        sage: d = DiGraph(g)
        sage: strongly_connected_component_containing_vertex(d,0)
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: d.strongly_connected_component_containing_vertex(0)
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    TESTS:

    If ``G`` is not a Sage DiGraph, an error is raised::

        sage: strongly_connected_component_containing_vertex('I am not a graph',0)
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage DiGraph
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise TypeError("the input must be a Sage DiGraph")

    if G.order()==1:
        return [v]

    try:
        return G._backend.strongly_connected_component_containing_vertex(v)

    except AttributeError:
        raise AttributeError("This function is only defined for C graphs.")


def strong_articulation_points(G):
    r"""
    Return the strong articulation points of this digraph.

    A vertex is a strong articulation point if its deletion increases the
    number of strongly connected components. This method implements the
    algorithm described in [ILS2012]_. The time complexity is dominated by
    the time complexity of the immediate dominators finding algorithm.

    OUTPUT: The list of strong articulation points.

    EXAMPLES:

    Two cliques sharing a vertex::

        sage: from sage.graphs.connectivity import strong_articulation_points
        sage: D = digraphs.Complete(4)
        sage: D.add_clique([3, 4, 5, 6])
        sage: strong_articulation_points(D)
        [3]
        sage: D.strong_articulation_points()
        [3]

    Two cliques connected by some arcs::

        sage: D = digraphs.Complete(4) * 2
        sage: D.add_edges([(0, 4), (7, 3)])
        sage: sorted( strong_articulation_points(D) )
        [0, 3, 4, 7]
        sage: D.add_edge(1, 5)
        sage: sorted( strong_articulation_points(D) )
        [3, 7]
        sage: D.add_edge(6, 2)
        sage: strong_articulation_points(D)
        []

    .. SEEALSO::

        - :meth:`~sage.graphs.digraph.DiGraph.strongly_connected_components`
        - :meth:`~sage.graphs.base.boost_graph.dominator_tree`

    TESTS:

    All strong articulation points are found::

        sage: from sage.graphs.connectivity import strong_articulation_points
        sage: def sap_naive(G):
        ....:     nscc = len(G.strongly_connected_components())
        ....:     S = []
        ....:     for u in G:
        ....:         H = copy(G)
        ....:         H.delete_vertex(u)
        ....:         if len(H.strongly_connected_components()) > nscc:
        ....:             S.append(u)
        ....:     return S
        sage: D = digraphs.RandomDirectedGNP(20, 0.1)
        sage: X = sap_naive(D)
        sage: SAP = strong_articulation_points(D)
        sage: set(X) == set(SAP)
        True

    Trivial cases::

        sage: strong_articulation_points(DiGraph())
        []
        sage: strong_articulation_points(DiGraph(1))
        []
        sage: strong_articulation_points(DiGraph(2))
        []

    If ``G`` is not a Sage DiGraph, an error is raised::

        sage: strong_articulation_points('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage DiGraph
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise TypeError("the input must be a Sage DiGraph")

    # The method is applied on each strongly connected component
    if is_strongly_connected(G):
        # Make a mutable copy of self
        L = [ DiGraph( [(u, v) for u, v in G.edge_iterator(labels=0) if u != v],
                           data_structure='sparse', immutable=False) ]
    else:
        # Get the list of strongly connected components of self as mutable
        # subgraphs
        L = [ G.subgraph(scc, immutable=False) for scc in G.strongly_connected_components() ]

    SAP = []
    for g in L:
        n = g.order()
        if n <= 1:
            continue
        if n == 2:
            SAP.extend( g.vertices() )
            continue

        # 1. Choose arbitrarily a vertex r, and test whether r is a strong
        # articulation point.
        r = next(g.vertex_iterator())
        E = g.incoming_edges(r) + g.outgoing_edges(r)
        g.delete_vertex(r)
        if not is_strongly_connected(g):
            SAP.append(r)
        g.add_edges(E)

        # 2. Compute the set of non-trivial immediate dominators in g
        Dr = set( g.dominator_tree(r, return_dict=True).values() )

        # 3. Compute the set of non-trivial immediate dominators in the
        # reverse digraph
        DRr = set( g.dominator_tree(r, return_dict=True, reverse=True).values() )

        # 4. Store D(r) + DR(r) - r
        SAP.extend( Dr.union(DRr).difference([r, None]) )

    return SAP

def bridges(G, labels=True):
    r"""
    Returns a list of the bridges (or cut edges).

    A bridge is an edge whose deletion disconnects the undirected graph.
    A disconnected graph has no bridge.

    INPUT:

    - ``labels`` -- (default: ``True``) if ``False``, each bridge is a tuple
      `(u, v)` of vertices

    EXAMPLES::

        sage: from sage.graphs.connectivity import bridges
        sage: from sage.graphs.connectivity import is_connected
        sage: g = 2*graphs.PetersenGraph()
        sage: g.add_edge(1,10)
        sage: is_connected(g)
        True
        sage: bridges(g)
        [(1, 10, None)]
        sage: g.bridges()
        [(1, 10, None)]


    TESTS:

    Ticket :trac:`23817` is solved::

        sage: G = Graph()
        sage: G.add_edge(0, 1)
        sage: bridges(G)
        [(0, 1, None)]
        sage: G.allow_loops(True)
        sage: G.add_edge(0, 0)
        sage: G.add_edge(1, 1)
        sage: bridges(G)
        [(0, 1, None)]

    If ``G`` is not a Sage Graph, an error is raised::

        sage: bridges('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be an Undirected Sage graph
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise TypeError("the input must be an Undirected Sage graph")

    # Small graphs and disconnected graphs have no bridge
    if G.order() < 2 or not is_connected(G):
        return []

    B,C = G.blocks_and_cut_vertices()

    # A block of size 2 is a bridge, unless the vertices are connected with
    # multiple edges.
    ME = set(G.multiple_edges(labels=False))
    my_bridges = []
    for b in B:
        if len(b) == 2 and not tuple(b) in ME:
            if labels:
                my_bridges.append((b[0], b[1], G.edge_label(b[0], b[1])))
            else:
                my_bridges.append(tuple(b))

    return my_bridges

# ==============================================================================
# Methods for finding 3-vertex-connected components and building SPQR-tree
# ==============================================================================

def cleave(G, cut_vertices=None, virtual_edges=True):
    r"""
    Return the connected subgraphs separated by the input vertex cut.

    Given a connected (multi)graph `G` and a vertex cut `X`, this method
    computes the list of subgraphs of `G` induced by each connected component
    `c` of `G\setminus X` plus `X`, i.e., `G[c\cup X]`.

    INPUT:

    - ``G`` -- a Graph.

    - ``cut_vertices`` -- (default: ``None``) a set of vertices representing a
      vertex cut of ``G``. If not vertex cut is given, the method will compute
      one via a call to :meth:`~sage.graphs.connectivity.vertex_connectivity`.

    - ``virtual_edges`` -- boolean (default: ``True``); whether to add virtual
      edges to the sides of the cut or not. A virtual edge is an edge between a
      pair of vertices of the cut that are not connected by an edge in ``G``.

    OUTPUT: A triple `(S, C, f)`, where

    - `S` is a list of the graphs that are sides of the vertex cut.

    - `C` is the graph of the cocycles. For each pair of vertices of the cut, it
      has one copy of each edge connecting them in ``G`` per sides of the cut
      plus one extra copy. Furthermore, when ``virtual_edges == True``, if a
      pair of vertices of the cut is not connected by an edge in ``G``, then it
      has one virtual edge between them per sides of the cut.

    - `f` is the complement of the subgraph of ``G`` induced but the vertex
      cut. Hence, its vertex set is the vertex cut, and its edge set is the set
      of virtual edges (i.e., edges between pairs of vertices of the cut that
      are not connected by an edge in ``G``).

    EXAMPLES:

    If there is an edge between cut vertices::

        sage: from sage.graphs.connectivity import cleave
        sage: G = Graph(2)
        sage: for _ in range(3):
        ....:     G.add_clique([0, 1, G.add_vertex(), G.add_vertex()])
        sage: S,C,f = cleave(G, cut_vertices=[0, 1])
        sage: [g.order() for g in S]
        [4, 4, 4]
        sage: C.order(), C.size()
        (2, 4)
        sage: f.vertices()
        [0, 1]

    If cut vertices doesn't have edge between them::

        sage: G.delete_edge(0, 1)
        sage: S,C,f = cleave(G, cut_vertices=[0, 1])
        sage: [g.order() for g in S]
        [4, 4, 4]
        sage: C.order(), C.size()
        (2, 3)
        sage: f.vertices(), f.edges()
        ([0, 1], [(0, 1, None)])

    If `G` is a biconnected multigraph::

        sage: G = graphs.CompleteBipartiteGraph(2,3)
        sage: G.add_edge(2, 3)
        sage: G.allow_multiple_edges(True)
        sage: G.add_edges(G.edges())
        sage: G.add_edges([(0, 1), (0, 1), (0, 1)])
        sage: S,C,f = cleave(G, cut_vertices=[0, 1])
        sage: for g in S:
        ....:     print(g.edges(labels=0))
        [(0, 1), (0, 1), (0, 1), (0, 2), (0, 2), (0, 3), (0, 3), (1, 2), (1, 2), (1, 3), (1, 3), (2, 3), (2, 3)]
        [(0, 1), (0, 1), (0, 1), (0, 4), (0, 4), (1, 4), (1, 4)]

    TESTS::

        sage: cleave(Graph(2))
        Traceback (most recent call last):
        ...
        ValueError: this method is designed for connected graphs only
        sage: cleave(graphs.PathGraph(3), cut_vertices=[5])
        Traceback (most recent call last):
        ...
        ValueError: vertex 5 is not a vertex of the input graph
        sage: cleave(Graph())
        Traceback (most recent call last):
        ...
        ValueError: the input graph has no vertex cut
        sage: cleave(graphs.CompleteGraph(5))
        Traceback (most recent call last):
        ...
        ValueError: the input graph has no vertex cut
        sage: cleave(graphs.CompleteGraph(5), cut_vertices=[3, 4])
        Traceback (most recent call last):
        ...
        ValueError: the set cut_vertices is not a vertex cut of the graph
    """
    if not G.is_connected():
        raise ValueError("this method is designed for connected graphs only")

    # If a vertex cut is given, we check that it is valid. Otherwise, we compute
    # a small vertex cut
    if cut_vertices:
        for u in cut_vertices:
            if not u in G:
                raise ValueError("vertex {} is not a vertex of the input graph".format(u))
        cut_vertices = list(cut_vertices)
    else:
        cut_size,cut_vertices = G.vertex_connectivity(value_only=False)
        if not cut_vertices:
            # Typical example is a clique
            raise ValueError("the input graph has no vertex cut")

    H = G.copy()
    H.delete_vertices(cut_vertices)
    CC = H.connected_components()
    if len(CC) == 1:
        raise ValueError("the set cut_vertices is not a vertex cut of the graph")

    # We identify the virtual edges, i.e., pairs of vertices of the vertex cut
    # that are not connected by an edge in G
    from sage.graphs.graph import Graph
    K = G.subgraph(cut_vertices)
    if virtual_edges:
        if K.allows_multiple_edges():
            virtual_cut_graph = K.to_simple().complement()
        else:
            virtual_cut_graph = K.complement()
    else:
        virtual_cut_graph = Graph()

    # We now build the graphs in each side of the cut, including the vertices
    # from the vertex cut
    cut_sides = []
    for comp in CC:
        h = G.subgraph(comp + cut_vertices)
        if virtual_edges:
            h.add_edges(virtual_cut_graph.edges())
        cut_sides.append(h)

    # We build the cocycles for re-assembly. For each edge between a pair of
    # vertices of the cut in the original graph G, a bond with one edge more
    # than the number of cut sides is needed. For pairs of vertices of the cut
    # that are not connected by an edge in G, a bond with one edge per cut side
    # is needed.
    cocycles = Graph(multiedges=True)
    if K.size():
        cocycles.add_edges(K.edges() * (len(cut_sides) + 1))
    if virtual_edges and virtual_cut_graph:
        cocycles.add_edges(virtual_cut_graph.edges() * len(cut_sides))

    return cut_sides, cocycles, virtual_cut_graph

def spqr_tree(G):
    r"""
    Decomposes a graph into cycles, cocycles, and 3-connected blocks summed over
    cocycles.

    Splits a two-connected graph at each two-vertex separation, giving a unique
    decomposition into 3-connected blocks, cycles, and cocycles. The cocycles
    are dipole graphs with one edge per real edge between the included vertices
    and one additional (virtual) edge per connected component resulting from
    deletion of the vertices in the cut.

    OUTPUT: `(R,S,P,SPR,tricomp,SPQR-tree)` where `R` is a list of rigid
    (three-connected) graphs, `S` a list of series graphs(cycles), `P` a list of
    parallel classes(cocycle graphs), `SPR` the so-called three-block into which
    a two-connected graph uniquely decomposes, tricomp a list of all
    triconnected components which is generated from `R` and `S` blocks and
    `SQPR-tree` a tree whose vertices are labeled with the vertices of
    three-blocks in the decomposition and the block's type.

    EXAMPLES::

        sage: from sage.graphs.connectivity import spqr_tree
        sage: G = Graph({0:[1,2],3:[1,2],4:[1,2],5:[1,2],6:[1,2]})
        sage: R,S,P,SPR,tricomp,SPQR_tree = spqr_tree(G)
        sage: R
        []
        sage: C3 = graphs.CycleGraph(3)
        sage: all(h.is_isomorphic(C3) for h in S)
        True

        sage: G = Graph(2)
        sage: for i in range(3):
        ....:     G.add_clique([0, 1, G.add_vertex(), G.add_vertex()])
        sage: R,S,P,SPR,tricomp,SPQR_tree = spqr_tree(G)
        sage: CG4 = graphs.CompleteGraph(4)
        sage: all(h.is_isomorphic(CG4) for h in R)
        True

        sage: G = graphs.CompleteBipartiteGraph(2,3)
        sage: G.add_edge(2, 3)
        sage: G.allow_multiple_edges(True)
        sage: G.add_edges(G.edges())
        sage: G.add_edges([[0,1],[0,1],[0,1]])
        sage: R,S,P,SPR,tricomp,SPQR_tree = spqr_tree(G)
        sage: tricomp
        [Subgraph of (): Multi-graph on 3 vertices,
        Subgraph of (Complete bipartite graph): Graph on 4 vertices]

    TESTS::

    sage: G = Graph({0:[1,2,3,4,5,6],1:[2,3],2:[3],4:[5,6],5:[6]})
    sage: spqr_tree(G)
    Traceback (most recent call last):
    ...
    ValueError: generation of SPQR trees is only implemented for 2-connected graphs
    """
    from sage.graphs.graph import Graph
    cut_size, cut_vertices = G.vertex_connectivity(value_only=False)

    if cut_size < 2:
        raise ValueError("generation of SPQR trees is only implemented for 2-connected graphs")
    elif cut_size > 2:
        return [G],[],[],[G],[G],Graph({('R',tuple(G.vertices())):[]})
    elif G.is_cycle():
        return [],[G],[],[G],[G],Graph({('S',tuple(G.vertices())):[]})

    R_blocks = []
    two_blocks = []
    cocycles = []
    cuts = []
    cycles = []
    SG = G.copy()

    # Split_multiple_edge Algorithm, it will make SG as simple graph.  If G is a
    # multigraph, simplify while recording virtual edges that will be needed to
    # 2-sum with cocycles during reconstruction. After this step, next step will
    # be finding S and R, blocks,
    if SG.has_multiple_edges():
        mults = sorted(SG.multiple_edges(labels=False))
        for i in range(len(mults) - 1):
            if mults[i] == mults[i + 1]:
                cocycles.append(frozenset(mults[i]))
                SG.delete_edge(mults[i])
            else:
                cuts.append(frozenset(mults[i]))
        cuts.append(frozenset(mults[-1]))

    # If G simplifies to a cycle.  Otherwise, seed the list of blocks to be
    # 2-split with G
    if SG.is_cycle():
        cycles.extend([frozenset(e) for e in SG.edge_iterator(labels=False)])
    else:
        SG.allow_multiple_edges(False)
        two_blocks.append((SG, cut_vertices))


    # Each minor of G in two_blocks has a 2-vertex cut; we split at this cut and
    # check each side for S or R or a 2-vertex cut
    while two_blocks:
        B,B_cut = two_blocks.pop()
        # B will be always simple graph.
        S, C, f = cleave(B, cut_vertices=B_cut)
        for K in S:
            if K.is_cycle():
                # Add all the edges of K to cycles list.
                cycles.extend(frozenset(e) for e in K.edge_iterator(labels=False))
            else:
                K_cut_size,K_cut_vertices = K.vertex_connectivity(value_only=False)
                if K_cut_size == 2:
                    two_blocks.append((K, K_cut_vertices))
                else:
                    R_blocks.append(K)
        # add a P block if cleave says we need it; mark virtual edge
        cocycles.extend([frozenset(e) for e in C.edge_iterator(labels=False)])
        for e in f.edge_iterator(labels=False):
            fe = frozenset(e)
            if not fe in cuts:
                cuts.append(fe)

    # Cycles may or may not(if G is a cycle with order > 3) triangulated; We
    # must undo this for S-blocks Virtual edges to be used in cycle assembly
    # from triangles will be multiple edges; check to be sure they are not
    # already in a P-block, for then we cannot 2-sum cycles at this edge for a
    # S-block Reconstruct S-blocks from smaller polygons, where cocycles permit.

    cycles_graph = Graph(cycles, multiedges=True)
    for e in cycles_graph.edges(labels=False):
        if cycles_graph.subgraph(edges=[e]).has_multiple_edges():
            f = frozenset(e)
            if not f in cocycles:
                cuts.remove(f)
                while cycles_graph.has_edge(e):
                    cycles_graph.delete_edge(e)

    # Things which are getting deleted from cycles_graph store it in tmp.
    tmp = []
    polygons = []
    for component in cycles_graph.connected_components_subgraphs():
        for block in component.blocks_and_cut_vertices()[0]:
            tmp.append(cycles_graph.subgraph(vertices=block))

    # After above step, there is no use of cycles_graph.
    del cycles_graph

    # tmp will be seeded with the union of cycles; each multiple edge
    # corresponds to a sequence of 2-sums of polygons; when cycles have no more
    # multiple edges they are S-blocks and are removed from list

    while tmp:
        block = tmp.pop()
        #print('b'," Bi-Connected:",block.is_biconnected()," Is-Cycle:",block.is_cycle()," Multi-Edge:",block.has_multiple_edges())
        if block.is_cycle():
            polygons.append(block)
        elif block.has_multiple_edges():
            f = block.multiple_edges()[0]
            SS = block.subgraph(vertices=[v for v in block.vertices() if v not in [f[0],f[1]]])
            for SK in SS.connected_components():
                SKS = block.subgraph(vertices=SK+[f[0],f[1]])
                while SKS.has_edge(f):
                    SKS.delete_edge(f)
                SKS.add_edge(f)
                #print("SKS",SKS.is_biconnected())
                tmp.append(SKS)
        else:
            R_blocks.append(block)


    # as with blocks_and_cuts_tree, we name vertices of the edge-sum tree with a
    # type (P = cocycle, S = cycle or R = three-block) and the list of G's
    # vertices in this block, noting only cocycles have multiple edges

    SPR = R_blocks + polygons
    Treeverts = []
    Tree = Graph()
    for T in SPR:
        if len(T) == 2:
            # It's a virtual edge component.
            block_type = 'P'
        elif T.is_cycle():
            block_type = 'S'
        else:
            block_type = 'R'
        Treeverts.append((block_type,tuple(T.vertices())))
    if len(Treeverts) == 1:
        Tree.add_vertex(Treeverts[0])

    P_blocks = []
    cocycles_graph = Graph(cocycles, multiedges=True)

    for cut in cuts:
        u,v = cut
        if cut in cocycles:
            # Updating the P-Block
            P_blocks.append(cocycles_graph.subgraph(vertices=[u, v]))
            for i in range(len(SPR)):
                if SPR[i].has_edge(u, v):
                    Tree.add_edge(('P',(u, v)),Treeverts[i])
        else:
            for i in range(len(SPR)):
                if SPR[i].has_edge(u, v):
                    for j in range(i + 1, len(SPR)):
                        if SPR[j].has_edge(u, v):
                            Tree.add_edge(Treeverts[i],Treeverts[j])

    thcomps = []
    for component in polygons:
        if component.order()==3 and component.size()==3:
            thcomps.append(component)
    thcomps += R_blocks
    #print(process.get_memory_info()[0])
    return R_blocks, polygons, P_blocks, R_blocks + polygons + P_blocks, thcomps, Tree
