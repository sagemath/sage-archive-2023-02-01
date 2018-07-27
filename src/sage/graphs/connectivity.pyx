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
    :meth:`cleave` | Return the connected subgraphs separated by the input vertex cut.
    :meth:`spqr_tree` | Return a SPQR-tree representing the triconnected components of the graph.
    :meth:`spqr_tree_to_graph` | Return the graph represented by the SPQR-tree `T`.

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
      vertex cut of ``G``. If no vertex cut is given, the method will compute
      one via a call to :meth:`~sage.graphs.connectivity.vertex_connectivity`.

    - ``virtual_edges`` -- boolean (default: ``True``); whether to add virtual
      edges to the sides of the cut or not. A virtual edge is an edge between a
      pair of vertices of the cut that are not connected by an edge in ``G``.

    OUTPUT: A triple `(S, C, f)`, where

    - `S` is a list of the graphs that are sides of the vertex cut.

    - `C` is the graph of the cocycles. For each pair of vertices of the cut,
      if there exists an edge between them, `C` has one copy of each edge
      connecting them in ``G`` per sides of the cut plus one extra copy.
      Furthermore, when ``virtual_edges == True``, if a pair of vertices of the
      cut is not connected by an edge in ``G``, then it has one virtual edge
      between them per sides of the cut.

    - `f` is the complement of the subgraph of ``G`` induced by the vertex
      cut. Hence, its vertex set is the vertex cut, and its edge set is the set
      of virtual edges (i.e., edges between pairs of vertices of the cut that
      are not connected by an edge in ``G``). When ``virtual_edges == False``,
      the edge set is empty.

    EXAMPLES:

    If there is an edge between cut vertices::

        sage: from sage.graphs.connectivity import cleave
        sage: G = Graph(2)
        sage: for _ in range(3):
        ....:     G.add_clique([0, 1, G.add_vertex(), G.add_vertex()])
        sage: S1,C1,f1 = cleave(G, cut_vertices=[0, 1])
        sage: [g.order() for g in S1]
        [4, 4, 4]
        sage: C1.order(), C1.size()
        (2, 4)
        sage: f1.vertices(), f1.edges()
        ([0, 1], [])

    If ``virtual_edges == False`` and there is an edge between cut vertices::

        sage: G.subgraph([0, 1]).complement() == Graph([[0, 1], []])
        True
        sage: S2,C2,f2 = cleave(G, cut_vertices=[0, 1], virtual_edges = False)
        sage: (S1 == S2, C1 == C2, f1 == f2)
        (True, True, True)

    If cut vertices doesn't have edge between them::

        sage: G.delete_edge(0, 1)
        sage: S1,C1,f1 = cleave(G, cut_vertices=[0, 1])
        sage: [g.order() for g in S1]
        [4, 4, 4]
        sage: C1.order(), C1.size()
        (2, 3)
        sage: f1.vertices(), f1.edges()
        ([0, 1], [(0, 1, None)])

    If ``virtual_edges == False`` and the cut vertices are not connected by an
    edge::

        sage: G.subgraph([0, 1]).complement() == Graph([[0, 1], []])
        False
        sage: S2,C2,f2 = cleave(G, cut_vertices=[0, 1], virtual_edges = False)
        sage: [g.order() for g in S2]
        [4, 4, 4]
        sage: C2.order(), C2.size()
        (2, 0)
        sage: f2.vertices(), f2.edges()
        ([0, 1], [])
        sage: (S1 == S2, C1 == C2, f1 == f2)
        (False, False, False)

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
        virtual_cut_graph = Graph([cut_vertices, []])

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
    cocycles = Graph([cut_vertices, []], multiedges=True)
    if K.size():
        cocycles.add_edges(K.edges() * (len(cut_sides) + 1))
    if virtual_edges and virtual_cut_graph:
        cocycles.add_edges(virtual_cut_graph.edges() * len(cut_sides))

    return cut_sides, cocycles, virtual_cut_graph

def spqr_tree(G):
    r"""
    Return an SPQR-tree representing the triconnected components of the graph.

    An SPQR-tree is a tree data structure used to represent the triconnected
    components of a biconnected (multi)graph and the 2-vertex cuts separating
    them. A node of a SPQR-tree, and the graph associated with it, can be one of
    the following four types:

    - ``S`` -- the associated graph is a cycle with at least three vertices.
      ``S`` stands for ``series``.

    - ``P`` -- the associated graph is a dipole graph, a multigraph with two
      vertices and three or more edges. ``P`` stands for ``parallel``.

    - ``Q`` -- the associated graph has a single real edge. This trivial case is
      necessary to handle the graph that has only one edge.

    - ``R`` -- the associated graph is a 3-connected graph that is not a cycle
      or dipole. ``R`` stands for ``rigid``.

    This method decomposes a biconnected graph into cycles, cocycles, and
    3-connected blocks summed over cocycles, and arranges them as a SPQR-tree.
    More precisely, it splits the graph at each of its 2-vertex cuts, giving a
    unique decomposition into 3-connected blocks, cycles and cocycles. The
    cocycles are dipole graphs with one edge per real edge between the included
    vertices and one additional (virtual) edge per connected component resulting
    from deletion of the vertices in the cut. See :wikipedia:`SPQR_tree`.

    OUTPUT: ``SPQR-tree`` a tree whose vertices are labeled with the block's type
    and the subgraph of three-blocks in the decomposition.

    EXAMPLES::

        sage: from sage.graphs.connectivity import spqr_tree
        sage: G = Graph(2)
        sage: for i in range(3):
        ....:     G.add_clique([0, 1, G.add_vertex(), G.add_vertex()])
        sage: Tree = spqr_tree(G)
        sage: Tree.order()
        4
        sage: K4 = graphs.CompleteGraph(4)
        sage: all(u[1].is_isomorphic(K4) for u in Tree.vertices() if u[0] == 'R')
        True
        sage: from sage.graphs.connectivity import spqr_tree_to_graph
        sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
        True

        sage: G = Graph(2)
        sage: for i in range(3):
        ....:     G.add_path([0, G.add_vertex(), G.add_vertex(), 1])
        sage: Tree = spqr_tree(G)
        sage: Tree.order()
        4
        sage: C4 = graphs.CycleGraph(4)
        sage: all(u[1].is_isomorphic(C4) for u in Tree.vertices() if u[0] == 'S')
        True
        sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
        True

        sage: G.allow_multiple_edges(True)
        sage: G.add_edges(G.edges())
        sage: Tree = spqr_tree(G)
        sage: Tree.order()
        13
        sage: all(u[1].is_isomorphic(C4) for u in Tree.vertices() if u[0] == 'S')
        True
        sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
        True

        sage: G = graphs.CycleGraph(6)
        sage: Tree = spqr_tree(G)
        sage: Tree.order()
        1
        sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
        True
        sage: G.add_edge(0, 3)
        sage: Tree = spqr_tree(G)
        sage: Tree.order()
        3
        sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
        True

    TESTS::

        sage: G = graphs.PathGraph(4)
        sage: spqr_tree(G)
        Traceback (most recent call last):
        ...
        ValueError: generation of SPQR-trees is only implemented for 2-connected graphs

        sage: G = Graph([(0, 0)], loops=True)
        sage: spqr_tree(G)
        Traceback (most recent call last):
        ...
        ValueError: generation of SPQR-trees is only implemented for graphs without loops
    """
    from sage.graphs.graph import Graph
    from collections import Counter

    if G.has_loops():
        raise ValueError("generation of SPQR-trees is only implemented for graphs without loops")

    cut_size, cut_vertices = G.vertex_connectivity(value_only=False)

    if cut_size < 2:
        raise ValueError("generation of SPQR-trees is only implemented for 2-connected graphs")
    elif cut_size > 2:
        return Graph({('R', Graph(G, immutable=True)):[]}, name='SPQR-tree of {}'.format(G.name()))

    # Split_multiple_edge Algorithm. If the input graph has multiple edges, we
    # make SG a simple graph while recording virtual edges that will be needed
    # to 2-sum with cocycles during reconstruction. After this step, next steps
    # will be finding S and R blocks.
    if G.has_multiple_edges():
        SG = G.to_simple()
        counter_multiedges = Counter([frozenset(e) for e in G.multiple_edges(labels=False)])
    else:
        SG = G
        counter_multiedges = Counter()

    # If SG simplifies to a cycle, we return the corresponding SPQR tree
    if SG.is_cycle():
        T = Graph(name='SPQR-tree of {}'.format(G.name()))
        S_node = ('S', Graph(SG, immutable=True))
        T.add_vertex(S_node)
        for e,num in counter_multiedges.items():
            P_node = ('P', Graph([e] * (num + 1), multiedges=True, immutable=True))
            T.add_edge(S_node, P_node)
        return T

    # We now search for 2-vertex cuts. If a cut is found, we split the graph and
    # check each side for S or R blocks or another 2-vertex cut
    R_blocks = []
    virtual_edges = set()
    two_blocks = [(SG, cut_vertices)]
    cycles_graph = Graph(multiedges=True)
    cocycles_count = Counter()

    while two_blocks:
        B,B_cut = two_blocks.pop()
        # B will be always simple graph.
        S, C, f = cleave(B, cut_vertices=B_cut)
        for K in S:
            if K.is_cycle():
                # Add all the edges of K to cycles list.
                cycles_graph.add_edges(e for e in K.edge_iterator(labels=False))
            else:
                K_cut_size,K_cut_vertices = K.vertex_connectivity(value_only=False)
                if K_cut_size == 2:
                    # The graph has a 2-vertex cut. We add it to the stack
                    two_blocks.append((K, K_cut_vertices))
                else:
                    # The graph is 3-vertex connected
                    R_blocks.append(('R', Graph(K, immutable=True)))
        # Store edges of cocycle (P block) and virtual edges (if any)
        cocycles_count.update(frozenset(e) for e in C.edge_iterator(labels=False))
        virtual_edges.update(frozenset(e) for e in f.edge_iterator(labels=False))


    # Cycles of order > 3 may have been triangulated; We undo this to reduce
    # the number of S-blocks. We start removing edges of the triangulation.
    count = Counter([frozenset(e) for e in cycles_graph.multiple_edges(labels=False)])
    for e,num in count.items():
        if num == 2 and e in virtual_edges:
            virtual_edges.discard(e)
            for _ in range(num):
                cycles_graph.delete_edge(e)
                cocycles_count[e] -= 1

    # We then extract cycles to form S_blocks
    S_blocks = []
    if not cycles_graph:
        tmp = []
    elif cycles_graph.is_connected():
        tmp = [cycles_graph]
    else:
        tmp = list(cycles_graph.connected_components_subgraphs())
    while tmp:
        block = tmp.pop()
        if block.is_cycle():
            S_blocks.append(('S', Graph(block, immutable=True)))
        elif block.has_multiple_edges():
            cut = block.multiple_edges(labels=False)[0]
            S,C,f = cleave(block, cut_vertices=cut)
            # Remove the cut edge from `S block` components if present
            for h in S:
                while h.has_edge(cut):
                    h.delete_edge(cut)
                h.add_edge(cut)
                tmp.append(h)
        else:
            # This should never happen
            raise ValueError("something goes wrong")


    # We now build the SPQR tree
    Tree = Graph(name='SPQR tree of {}'.format(G.name()))
    SR_blocks = S_blocks + R_blocks
    Tree.add_vertices(SR_blocks)
    for e,num in cocycles_count.items():
        if num:
            P_block = ('P', Graph([e] * (num + max(0, counter_multiedges[e] - 1)), multiedges=True, immutable=True))
            for block in SR_blocks:
                # Note: here we use a try...except statement since the immutable
                # graph backend raises an error if an end vertex of the edge is
                # not in the graph.
                try:
                    if block[1].has_edge(e):
                        Tree.add_edge(block, P_block)
                except:
                    continue

    # We finally add P blocks to account for multiple edges of the input graph
    # that are not involved in any separator of the graph
    for e,num in counter_multiedges.items():
        if not cocycles_count[e]:
            P_block = ('P', Graph([e] * (num + 1), multiedges=True, immutable=True))
            for block in SR_blocks:
                try:
                    if block[1].has_edge(e):
                        Tree.add_edge(block, P_block)
                        break
                except:
                    continue

    return Tree

def spqr_tree_to_graph(T):
    r"""
    Return the graph represented by the SPQR-tree `T`.

    The main purpose of this method is to test :meth:`spqr_tree`.

    INPUT:

    - ``T`` -- a SPQR tree as returned by :meth:`spqr_tree`.

    OUTPUT: a (multi) graph

    EXAMPLES:

    :wikipedia:`SPQR_tree` reference paper example::

        sage: from sage.graphs.connectivity import spqr_tree
        sage: from sage.graphs.connectivity import spqr_tree_to_graph
        sage: G = Graph([(1, 2), (1, 4), (1, 8), (1, 12), (3, 4), (2, 3),
        ....: (2, 13), (3, 13), (4, 5), (4, 7), (5, 6), (5, 8), (5, 7), (6, 7),
        ....: (8, 11), (8, 9), (8, 12), (9, 10), (9, 11), (9, 12), (10, 12)])
        sage: T = spqr_tree(G)
        sage: H = spqr_tree_to_graph(T)
        sage: H.is_isomorphic(G)
        True

    A small multigraph ::

        sage: G = Graph([(0, 2), (0, 2), (1, 3), (2, 3)], multiedges=True)
        sage: for i in range(3):
        ....:     G.add_clique([0, 1, G.add_vertex(), G.add_vertex()])
        sage: for i in range(3):
        ....:     G.add_clique([2, 3, G.add_vertex(), G.add_vertex()])
        sage: T = spqr_tree(G)
        sage: H = spqr_tree_to_graph(T)
        sage: H.is_isomorphic(G)
        True

    TESTS::

        sage: H = spqr_tree_to_graph(Graph())
        sage: H.is_isomorphic(Graph())
        True
    """
    from sage.graphs.graph import Graph
    from collections import Counter

    count_G = Counter()
    for t,g in T.vertex_iterator():
        if not t == 'P':
            count_G.update(g.edge_iterator(labels=False))

    for t,g in T.vertex_iterator():
        if t == 'P':
            count_g = Counter(g.edge_iterator(labels=False))
            for e,num in count_g.items():
                count_G[e] = abs(count_g[e] - count_G[e])

    G = Graph(multiedges=True)
    for e,num in count_G.items():
        for _ in range(num):
            G.add_edge(e)

    return G

class LinkedListNode:
    """
    Helper class for ``Triconnectivity``.
    Node in a linked list.
    Has pointers to its previous node and next node.
    If this node is the `head` of the linked list, reference to the linked list
    object is stored in `listobj`.
    """
    prev = None
    next = None
    data = None #edge or int
    listobj = None
    def set_data(self, e):
        self.data = e
    def get_data(self):
        return self.data
    def set_obj(self, l):
        self.listobj = l
    def clear_obj(self):
        self.listobj = None
    def replace(self, node):
        """
        Replace self node with ``node`` in the corresponding linked list.
        """
        if self.prev == None and self.next == None:
            self.listobj.set_head(node)
        elif self.prev == None:
            self.listobj.head = node
            node.next = self.next
            node.listobj = self.listobj
        elif self.next == None:
            self.prev.next = node
            node.prev = self.prev
        else:
            self.prev.next = node
            self.next.prev = node
            node.prev = self.prev
            node.next = self.next

class LinkedList:
    """
    A helper class for ``Triconnectivity``.
    A linked list with a head and a tail pointer
    """
    head = None
    tail = None
    length = 0
    def remove(self, node):
        """
        Remove the node ``node`` from the linked list.
        """
        if node.prev == None and node.next == None:
            self.head = None
            self.tail = None
        elif node.prev == None: # node is head
            self.head = node.next
            node.next.prev = None
            node.next.set_obj(self)
        elif node.next == None: #node is tail
            node.prev.next = None
            self.tail = node.prev
        else:
            node.prev.next = node.next
            node.next.prev = node.prev
        self.length -= 1
    def set_head(self, h):
        """
        Set the node ``h`` as the head of the linked list.
        """
        self.head = h
        self.tail = h
        self.length = 1
        h.set_obj(self)
    def append(self, node):
        """
        Append the node ``node`` to the linked list.
        """
        if self.head == None:
            self.set_head(node)
        else:
            self.tail.next = node
            node.prev = self.tail
            self.tail = node
            self.length += 1
    def get_head(self):
        return self.head
    def get_length(self):
        return self.length
    def replace(self, node1, node2):
        """
        Replace the node ``node1`` with ``node2`` in the linked list.
        """
        if node1.prev == None and node1.next == None:
            self.head = node2
            self.tail = node2
        elif node1.prev == None: # head has to be replaced
            node1.next.prev = node2
            node2.next = node1.next
            self.head = node2
        elif node1.next == None: # tail has to be replaced
            node1.prev.next = node2
            node2.prev = node1.prev
            self.tail = node2
        else:
            node1.prev.next = node2
            node1.next.prev = node2
            node2.prev = node1.prev
            node2.next = node1.next
    def push_front(self, node):
        """
        Add node ``node`` to the beginning of the linked list.
        """
        if self.head == None:
            self.head = node
            self.tail = node
            node.set_obj(self)
        else:
            self.head.clear_obj()
            self.head.prev = node
            node.next = self.head
            self.head = node
            node.set_obj(self)
        self.length += 1
    def to_string(self):
        temp = self.head
        s = ""
        while temp:
            s += "  " + str(temp.get_data())
            temp = temp.next
        return s
    def concatenate(self, lst2):
        """
        Concatenates lst2 to self.
        Makes lst2 empty.
        """
        self.tail.next = lst2.head
        lst2.head.prev = self.tail
        self.tail = lst2.tail
        self.length += lst2.length
        lst2.head = None
        lst2.length = 0

class Component:
    """
    A helper class for ``Triconnectivity``.
    A connected component.
    `edge_list` contains the list of edges belonging to the component.
    `component_type` stores the type of the component.
        - 0 if bond.
        - 1 if polygon.
        - 2 is triconnected component.
    """
    edge_list = LinkedList()
    component_type = 0  #bond = 0, polygon = 1, triconnected = 2
    def __init__(self, edge_list, type_c):
        """
        `edge_list` is a list of edges to be added to the component.
        `type_c` is the type of the component.
        """
        self.edge_list = LinkedList()
        for e in edge_list:
            e_node = LinkedListNode()
            e_node.set_data(e)
            self.edge_list.append(e_node)
        self.component_type = type_c
    def add_edge(self, e):
        e_node = LinkedListNode()
        e_node.set_data(e)
        self.edge_list.append(e_node)
    def finish_tric_or_poly(self, e):
        """
        Edge `e` is the last edge to be added to the component.
        Classify the component as a polygon or triconnected component
        depending on the number of edges belonging to it.
        """
        e_node = LinkedListNode()
        e_node.set_data(e)
        self.edge_list.append(e_node)
        if self.edge_list.get_length() >= 4:
            self.component_type = 2
        else:
            self.component_type = 1
    def __str__(self):
        """
        Function for printing the component.
        """
        if self.component_type == 0:
            type_str = "Bond: "
        elif self.component_type == 1:
            type_str =  "Polygon: "
        else:
            type_str = "Triconnected: "
        return type_str + self.edge_list.to_string()
    def get_edge_list(self):
        """
        Return a list of edges belonging to the component.
        """
        e_list = []
        e_node = self.edge_list.get_head()
        while e_node:
            e_list.append(e_node.get_data())
            e_node = e_node.next
        return e_list

from sage.graphs.base.sparse_graph cimport SparseGraph

class Triconnectivity:
    """
    This module implements the algorithm for finding the triconnected
    components of a biconnected graph.
    A biconnected graph is a graph where deletion of any one vertex does
    not disconnect the graph.

    INPUT:

    - ``G`` -- The input graph.

    - ``check`` (default: ``True``) -- Boolean to indicate whether ``G``
        needs to be tested for biconnectivity.

    OUTPUT:

    No output, the triconnected components are printed.
    The triconnected components are stored in `comp_list_new and `comp_type`.
    `comp_list_new` is a list of components, with `comp_list_new[i]` contains
    the list of edges in the $i^{th}$ component. `comp_type[i]` stores the type
    of the $i^{th}$ component - 1 for bond, 2 for polygon, 3 for triconnected
    component. The output can be accessed through these variables.

    ALGORITHM:

      We implement the algorithm proposed by Tarjan in [Tarjan72]_. The
      original version is recursive. We emulate the recursion using a stack.

    ALGORITHM::
    We implement the algorithm proposed by Hopcroft and Tarjan in
    [Hopcroft1973]_ and later corrected by Gutwenger and Mutzel in
    [Gut2001]_.

    .. SEEALSO::

        - :meth:`~Graph.is_biconnected`

    EXAMPLES::

        An example from [Hopcroft1973]_:

        sage: from sage.graphs.connectivity import Triconnectivity
        sage: G = Graph()
        sage: G.add_edges([(1,2),(1,4),(1,8),(1,12),(1,13),(2,3),(2,13),(3,4)])
        sage: G.add_edges([(3,13),(4,5),(4,7),(5,6),(5,7),(5,8),(6,7),(8,9),(8,11)])
        sage: G.add_edges([(8,12),(9,10),(9,11),(9,12),(10,11),(10,12)])
        sage: tric = Triconnectivity(G)
        Triconnected:  [(8, 9, None), (9, 10, None), (10, 11, None), (9, 11, None), (8, 11, None), (10, 12, None), (9, 12, None), (8, 12, 'newVEdge0')]
        Bond:  [(8, 12, None), (8, 12, 'newVEdge0'), (8, 12, 'newVEdge1')]
        Polygon:  [(8, 12, 'newVEdge1'), (1, 12, None), (8, 1, 'newVEdge2')]
        Bond:  [(1, 8, None), (8, 1, 'newVEdge2'), (8, 1, 'newVEdge3')]
        Polygon:  [(5, 8, None), (8, 1, 'newVEdge3'), (4, 5, 'newVEdge8'), (4, 1, 'newVEdge9')]
        Polygon:  [(5, 6, None), (6, 7, None), (5, 7, 'newVEdge5')]
        Bond:  [(5, 7, None), (5, 7, 'newVEdge5'), (5, 7, 'newVEdge6')]
        Polygon:  [(5, 7, 'newVEdge6'), (4, 7, None), (5, 4, 'newVEdge7')]
        Bond:  [(5, 4, 'newVEdge7'), (4, 5, 'newVEdge8'), (4, 5, None)]
        Bond:  [(1, 4, None), (4, 1, 'newVEdge9'), (4, 1, 'newVEdge10')]
        Polygon:  [(3, 4, None), (4, 1, 'newVEdge10'), (3, 1, 'newVEdge11')]
        Triconnected:  [(1, 2, None), (2, 3, None), (3, 1, 'newVEdge11'), (3, 13, None), (2, 13, None), (1, 13, None)]

        An example from [Gut2001]_

        sage: G = Graph()
        sage: G.add_edges([(1,2),(1,4),(2,3),(2,5),(3,4),(3,5),(4,5),(4,6),(5,7),(5,8)])
        sage: G.add_edges([(5,14),(6,8),(7,14),(8,9),(8,10),(8,11),(8,12),(9,10),(10,13)])
        sage: G.add_edges([(10,14),(10,15),(10,16),(11,12),(11,13),(12,13),(14,15),(14,16),(15,16)])
        sage: tric = Triconnectivity(G)
        Polygon:  [(6, 8, None), (4, 6, None), (5, 8, 'newVEdge12'), (5, 4, 'newVEdge13')]
        Polygon:  [(8, 9, None), (9, 10, None), (8, 10, 'newVEdge1')]
        Bond:  [(8, 10, 'newVEdge1'), (8, 10, None), (8, 10, 'newVEdge4'), (10, 8, 'newVEdge5')]
        Triconnected:  [(8, 11, None), (11, 12, None), (8, 12, None), (12, 13, None), (11, 13, None), (8, 13, 'newVEdge3')]
        Polygon:  [(8, 13, 'newVEdge3'), (10, 13, None), (8, 10, 'newVEdge4')]
        Triconnected:  [(10, 15, None), (14, 15, None), (15, 16, None), (10, 16, None), (14, 16, None), (10, 14, 'newVEdge6')]
        Bond:  [(10, 14, 'newVEdge6'), (14, 10, 'newVEdge7'), (10, 14, None)]
        Polygon:  [(14, 10, 'newVEdge7'), (10, 8, 'newVEdge5'), (5, 14, 'newVEdge10'), (5, 8, 'newVEdge11')]
        Polygon:  [(5, 7, None), (7, 14, None), (5, 14, 'newVEdge9')]
        Bond:  [(5, 14, None), (5, 14, 'newVEdge9'), (5, 14, 'newVEdge10')]
        Bond:  [(5, 8, None), (5, 8, 'newVEdge11'), (5, 8, 'newVEdge12')]
        Bond:  [(5, 4, 'newVEdge13'), (4, 5, 'newVEdge14'), (4, 5, None)]
        Triconnected:  [(2, 3, None), (3, 4, None), (4, 5, 'newVEdge14'), (3, 5, None), (2, 5, None), (2, 4, 'newVEdge15')]
        Polygon:  [(1, 2, None), (2, 4, 'newVEdge15'), (1, 4, None)]

        An example with multi-edges and accessing the triconnected components:

        sage: G = Graph()
        sage: G.allow_multiple_edges(True)
        sage: G.add_edges([(1,2),(2,3),(3,4),(4,5),(1,5),(1,5),(2,3)])
        sage: tric = Triconnectivity(G)
        Bond:  [(1, 5, None), (1, 5, None), (1, 5, 'newVEdge0')]
        Bond:  [(2, 3, None), (2, 3, None), (2, 3, 'newVEdge1')]
        Polygon:  [(4, 5, None), (1, 5, 'newVEdge0'), (3, 4, None), (1, 2,
        None), (2, 3, 'newVEdge1')]
        sage: tric.comp_list_new
        [[(1, 5, None), (1, 5, None), (1, 5, 'newVEdge0')],
        [(2, 3, None), (2, 3, None), (2, 3, 'newVEdge1')],
        [(4, 5, None), (1, 5, 'newVEdge0'), (3, 4, None), (1, 2, None), (2, 3, 'newVEdge1')]]
        sage: tric.comp_type
        [0, 0, 1]

        An example of a triconnected graph:

        sage: G2 = Graph()
        sage: G2.allow_multiple_edges(True)
        sage: G2.add_edges([('a','b'),('a','c'),('a','d'),('b','c'),('b','d'),('c','d')])
        sage: tric2 = Triconnectivity(G2)
        Triconnected:  [('a', 'b', None), ('b', 'c', None), ('a', 'c', None), ('c', 'd', None), ('b', 'd', None), ('a', 'd', None)]

    TESTS::

        A disconnected graph:

        sage: from sage.graphs.connectivity import Triconnectivity
        sage: G = Graph([(1,2),(3,5)])
        sage: tric = Triconnectivity(G)
        Traceback (most recent call last):
        ...
        ValueError: Graph is disconnected

        A graph with a cut vertex:

        sage: from sage.graphs.connectivity import Triconnectivity
        sage: G = Graph([(1,2),(1,3),(2,3),(3,4),(3,5),(4,5)])
        sage: tric = Triconnectivity(G)
        Traceback (most recent call last):
        ...
        ValueError: Graph has a cut vertex

    """
    def __init__(self, G, check=True):
        from sage.graphs.graph import Graph
        # graph_copy is a SparseGraph of the input graph `G`
        # We relabel the edges with increasing numbers to be able to
        # distinguish between multi-edges
        self.graph_copy = Graph(multiedges=True)

        # Add all the vertices first
        # there is a possibility of isolated vertices
        for v in G.vertex_iterator():
            self.graph_copy.add_vertex(v)

        edges = G.edges()
        # dict to map new edges with the old edges
        self.edge_label_dict = {}
        for i in range(len(edges)):
            newEdge = tuple([edges[i][0], edges[i][1], i])
            self.graph_copy.add_edge(newEdge)
            self.edge_label_dict[newEdge] = edges[i]

        # type SparseGraph
        self.graph_copy = self.graph_copy.copy(implementation='c_graph')

        # mapping of vertices to integers in c_graph
        self.vertex_to_int = self.graph_copy.relabel(inplace=True, return_map=True)

        # mapping of integers to original vertices
        self.int_to_vertex = dict([(v,k) for k,v in self.vertex_to_int.items()])
        self.n = self.graph_copy.order() # number of vertices
        self.m = self.graph_copy.size() # number of edges

        # status of each edge: unseen=0, tree=1, frond=2, removed=3
        self.edge_status = dict((e, 0) for e in self.graph_copy.edges())

        # Edges of the graph which are in the reverse direction in palm tree
        self.reverse_edges = set()
        self.dfs_number = [0 for i in range(self.n)] # DFS number of vertex i

        # Linked list of fronds entering vertex i in the order they are visited
        self.highpt = [LinkedList() for i in range(self.n)]

        # A dictionary whose key is an edge e, value is a pointer to element in
        # self.highpt containing the edge e. Used in the `path_search` function.
        self.in_high = dict((e, None) for e in self.graph_copy.edges())

        # Translates DFS number of a vertex to its new number
        self.old_to_new = [0 for i in range(self.n+1)]
        self.newnum = [0 for i in range(self.n)] # new number of vertex i
        self.node_at = [0 for i in range(self.n+1)] # node at dfs number of i
        self.lowpt1 = [None for i in range(self.n)] # lowpt1 number of vertex i
        self.lowpt2 = [None for i in range(self.n)] # lowpt2 number of vertex i

        # i^th value contains a LinkedList of incident edges of vertex i
        self.adj = [LinkedList() for i in range(self.n)]

        # A dictionary whose key is an edge, value is a pointer to element in
        # self.adj containing the edge. Used in the `path_search` function.
        self.in_adj = {}
        self.nd = [None for i in range(self.n)] # number of descendants of vertex i

        # Parent vertex of vertex i in the palm tree
        self.parent = [None for i in range(self.n)]
        self.degree = [None for i in range(self.n)] # Degree of vertex i
        self.tree_arc = [None for i in range(self.n)] # Tree arc entering the vertex i
        self.vertex_at = [1 for i in range(self.n)] # vertex with DFS number of i
        self.dfs_counter = 0
        self.components_list = [] # list of components of `graph_copy`
        self.graph_copy_adjacency = [[] for i in range(self.n)] # Stores adjacency list

        # Dictionary of (e, True/False) to denote if edge e starts a path
        self.starts_path = dict((e, False) for e in self.graph_copy.edges())

        self.is_biconnected = True # Boolean to store if the graph is biconnected or not
        self.cut_vertex = None # If graph is not biconnected

        # Label used for virtual edges, incremented at every new virtual edge
        self.virtual_edge_num = 0

        self.new_path = False # Boolean used to store if new path is started

        # Stacks used in `path_search` function
        self.e_stack = []
        self.t_stack_h = [None for i in range(2*self.m + 1)]
        self.t_stack_a = [None for i in range(2*self.m + 1)]
        self.t_stack_b = [None for i in range(2*self.m + 1)]
        self.t_stack_top = 0
        self.t_stack_a[self.t_stack_top] = -1

        # The final triconnected components are stored
        self.comp_list_new = [] # i^th entry is list of edges in i^th component
        self.comp_type = [] # i^th entry is type of i^th component


        # Trivial cases
        if self.n < 2:
            raise ValueError("Graph is not biconnected")
        if self.n == 2:
            if self.m < 3:
                raise ValueError("Graph is not biconnected")
            comp = Component([], 0)
            for e in self.graph_copy.edges():
                comp.add_edge(e)
            self.components_list.append(comp)
            return

        # Triconnectivity algorithm
        self.__split_multi_egdes()

        # Build adjacency list
        for e in self.graph_copy.edges():
            self.graph_copy_adjacency[e[0]].append(e)
            self.graph_copy_adjacency[e[1]].append(e)

        self.dfs_counter = 0 # Initialisation for dfs1()
        self.start_vertex = 0 # Initialisation for dfs1()
        self.cut_vertex = self.__dfs1(self.start_vertex, check=check)

        if check:
            # If graph is disconnected
            if self.dfs_counter < self.n:
                self.is_biconnected = False
                raise ValueError("Graph is disconnected")

            # If graph has a cut vertex
            if self.cut_vertex != None:
                self.cut_vertex = self.int_to_vertex[self.cut_vertex]
                self.is_biconnected = False
                raise ValueError("Graph has a cut vertex")

        # Identify reversed edges to reflect the palm tree arcs and fronds
        for e in self.graph_copy.edges():
            up = (self.dfs_number[e[1]] - self.dfs_number[e[0]]) > 0
            if (up and self.edge_status[e]==2) or (not up and self.edge_status[e]==1):
                # Add edge to the set reverse_edges
                self.reverse_edges.add(e)

        self.__build_acceptable_adj_struct()
        self.__dfs2()

        self.__path_search(self.start_vertex)

        # last split component
        c = Component([],0)
        while self.e_stack:
            c.add_edge(self.__estack_pop())
        c.component_type = 2 if c.edge_list.get_length() > 4 else 1
        self.components_list.append(c)
        c = None

        self.__assemble_triconnected_components()
        self.print_triconnected_components()

    def __tstack_push(self, h, a, b):
        """
        Push `(h,a,b)` triple on Tstack
        """
        self.t_stack_top += 1
        self.t_stack_h[self.t_stack_top] = h
        self.t_stack_a[self.t_stack_top] = a
        self.t_stack_b[self.t_stack_top] = b

    def __tstack_push_eos(self):
        """
        Push end-of-stack marker on Tstack
        """
        self.t_stack_top += 1
        self.t_stack_a[self.t_stack_top] = -1

    def __tstack_not_eos(self):
        """
        Return true iff end-of-stack marker is not on top of Tstack
        """
        return self.t_stack_a[self.t_stack_top] != -1

    def __estack_pop(self):
        """
        Pop from estack and return the popped element
        """
        e = self.e_stack[-1]
        self.e_stack = self.e_stack[0:-1]
        return e

    def __new_component(self, edges=[], type_c=0):
        """
        Create a new component, add `edges` to it.
        type_c = 0 for bond, 1 for polygon, 2 for triconnected component
        """
        c = Component(edges, type_c)
        self.components_list.append(c)
        return c

    def __high(self, v):
        """
        Return the high(v) value, which is the first value in highpt list of `v`
        """
        head = self.highpt[v].head
        if head == None:
            return 0
        else:
            return head.data

    def __del_high(self, e):
        if e in self.in_high:
            it = self.in_high[e]
            if it:
                if e in self.reverse_edges:
                    v = e[0]
                else:
                    v = e[1]
                self.highpt[v].remove(it)

    def __split_multi_egdes(self):
        """
        Iterate through all the edges, and constructs bonds wherever
        multiedges are present.

        If there are `k` multiple edges between `u` and `v`, then `k+1`
        edges will be added to a new component (one of them is a virtual edge),
        all the `k` edges are removed from the graph and a virtual edge is
        between `u` and `v` is added to the graph. Instead of deleting edges
        from the graph, they are instead marked as removed, i.e., 3 in
        the dictionary `edge_status`.

        No return value. Update the `components_list` and `graph_copy`.
        `graph_copy` will become simple graph after this function.

        """
        comp = []
        if self.graph_copy.has_multiple_edges():
            sorted_edges = sorted(self.graph_copy.edges())
            for i in range(len(sorted_edges) - 1):

                # Find multi edges and add to component and delete from graph
                if (sorted_edges[i][0] == sorted_edges[i + 1][0]) and \
                   (sorted_edges[i][1] == sorted_edges[i + 1][1]):
                    self.edge_status[sorted_edges[i]] = 3 # edge removed
                    comp.append(sorted_edges[i])
                else:
                    if comp:
                        comp.append(sorted_edges[i])
                        self.edge_status[sorted_edges[i]] = 3 # edge removed

                        # Add virtual edge to graph_copy
                        newVEdge = tuple([sorted_edges[i][0], sorted_edges[i][1], "newVEdge"+str(self.virtual_edge_num)])
                        self.graph_copy.add_edge(newVEdge)
                        self.virtual_edge_num += 1

                        # mark unseen for newVEdge
                        self.edge_status[newVEdge] = 0

                        comp.append(newVEdge)
                        self.__new_component(comp)
                    comp = []
            if comp:
                comp.append(sorted_edges[i+1])
                self.edge_status[sorted_edges[i+1]] = 3 # edge removed

                # Add virtual edge to graph_copy
                newVEdge = tuple([sorted_edges[i+1][0], sorted_edges[i+1][1], "newVEdge"+str(self.virtual_edge_num)])
                self.graph_copy.add_edge(newVEdge)
                self.virtual_edge_num += 1
                self.edge_status[newVEdge] = 0

                comp.append(newVEdge)
                self.__new_component(comp)

    def __dfs1(self, v, u=None, check=True):
        """
        This function builds the palm tree of the graph using a dfs traversal.
        Also populates the lists lowpt1, lowpt2, nd, parent, and dfs_number.
        It updates the dict `edge_status` to reflect palm tree arcs and fronds.

        Input::

        - ``v`` -- The start vertex for DFS.

        - ``u`` -- The parent vertex of ``v`` in the palm tree.

        - ``check`` -- if True, the graph is tested for biconnectivity. If the
            graph has a cut vertex, the cut vertex is returned. If set to False,
            the graph is assumed to be biconnected, function returns None.

        Output::

        - If ``check`` is set to True, and a cut vertex is found, the cut vertex
          is returned. If no cut vertex is found, return value is None.
          If ``check`` is set to False, ``None`` is returned.
        """
        first_son = None # For testing biconnectivity
        s1 = None # Storing the cut vertex, if there is one
        self.dfs_counter += 1
        self.dfs_number[v] = self.dfs_counter
        self.parent[v] = u
        self.degree[v] = self.graph_copy.degree(v)
        self.lowpt1[v] = self.lowpt2[v] = self.dfs_number[v]
        self.nd[v] = 1
        for e in self.graph_copy_adjacency[v]:
            if self.edge_status[e]:
                continue

            w = e[0] if e[0] != v else e[1] # Opposite vertex of edge e
            if self.dfs_number[w] == 0:
                self.edge_status[e] = 1 # tree edge
                if first_son is None:
                    first_son = w
                self.tree_arc[w] = e
                s1 = self.__dfs1(w, v, check)

                if check:
                    # Check for cut vertex.
                    # The situation in which there is no path from w to an
                    # ancestor of v : we have identified a cut vertex
                    if (self.lowpt1[w] >= self.dfs_number[v]) and (w != first_son or u != None):
                        s1 = v

                # Calculate the `lowpt1` and `lowpt2` values.
                # `lowpt1` is the smallest vertex (the vertex x with smallest
                # dfs_number[x]) that can be reached from v.
                # `lowpt2` is the next smallest vertex that can be reached from v.
                if self.lowpt1[w] < self.lowpt1[v]:
                        self.lowpt2[v] = min(self.lowpt1[v], self.lowpt2[w])
                        self.lowpt1[v] = self.lowpt1[w]

                elif self.lowpt1[w] == self.lowpt1[v]:
                    self.lowpt2[v] = min(self.lowpt2[v], self.lowpt2[w])

                else:
                    self.lowpt2[v] = min(self.lowpt2[v], self.lowpt1[w])

                self.nd[v] += self.nd[w]

            else:
                self.edge_status[e] = 2 #frond
                if self.dfs_number[w] < self.lowpt1[v]:
                    self.lowpt2[v] = self.lowpt1[v]
                    self.lowpt1[v] = self.dfs_number[w]
                elif self.dfs_number[w] > self.lowpt1[v]:
                    self.lowpt2[v] = min(self.lowpt2[v], self.dfs_number[w])

        return s1 # s1 is None is graph does not have a cut vertex


    def __build_acceptable_adj_struct(self):
        """
        Builds the adjacency lists for each vertex with certain properties
        of the ordering, using the ``lowpt1`` and ``lowpt2`` values.

        The list ``adj`` and the dictionary ``in_adj`` are populated.

        `phi` values of each edge are calculated using the `lowpt` values of
        incident vertices. The edges are then sorted by the `phi` values and
        added to adjacency list.
        """
        max = 3*self.n + 2
        bucket = [[] for i in range(max+1)]

        for e in self.graph_copy.edges():
            edge_type = self.edge_status[e]
            if edge_type == 3: # If edge status is `removed`, go to next edge
                continue

            # compute phi value
            # bucket sort adjacency list by phi values
            if e in self.reverse_edges:
                if edge_type==1: # tree arc
                    if self.lowpt2[e[0]] < self.dfs_number[e[1]]:
                        phi = 3*self.lowpt1[e[0]]
                    else:
                        phi = 3*self.lowpt1[e[0]] + 2
                else: # tree frond
                    phi = 3*self.dfs_number[e[0]]+1
            else:
                if edge_type==1: # tree arc
                    if self.lowpt2[e[1]] < self.dfs_number[e[0]]:
                        phi = 3*self.lowpt1[e[1]]
                    else:
                        phi = 3*self.lowpt1[e[1]] + 2
                else: # tree frond
                    phi = 3*self.dfs_number[e[1]]+1

            bucket[phi].append(e)

        # Populate `adj` and `in_adj` with the sorted edges
        for i in range(1,max+1):
            for e in bucket[i]:
                node = LinkedListNode()
                node.set_data(e)
                if e in self.reverse_edges:
                    self.adj[e[1]].append(node)
                    self.in_adj[e] = node
                else:
                    self.adj[e[0]].append(node)
                    self.in_adj[e] = node

    def __path_finder(self, v):
        """
        This function is a helper function for `dfs2` function.
        Calculate `newnum[v]` and identify the edges which start a new path.
        """
        self.newnum[v] = self.dfs_counter - self.nd[v] + 1
        e_node = self.adj[v].get_head()
        while e_node:
            e = e_node.get_data()
            e_node = e_node.next
            w = e[1] if e[0] == v else e[0] # opposite vertex of e
            if self.new_path:
                self.new_path = False
                self.starts_path[e] = True
            if self.edge_status[e] == 1: # tree arc
                self.__path_finder(w)
                self.dfs_counter -= 1
            else:
                # Identified a new frond that enters `w`. Add to `highpt[w]`.
                highpt_node = LinkedListNode()
                highpt_node.set_data(self.newnum[v])
                self.highpt[w].append(highpt_node)
                self.in_high[e] = highpt_node
                self.new_path = True


    def __dfs2(self):
        """
        Update the values of lowpt1 and lowpt2 lists with the help of
        new numbering obtained from `Path Finder` funciton.
        Populate `highpt` values.
        """
        self.in_high = dict((e, None) for e in self.graph_copy.edges())
        self.dfs_counter = self.n
        self.newnum = [0 for i in range(self.n)]
        self.starts_path = dict((e, False) for e in self.graph_copy.edges())

        self.new_path = True

        # We call the pathFinder function with the start vertex
        self.__path_finder(self.start_vertex)

        # Update `old_to_new` values with the calculated `newnum` values
        for v in self.graph_copy.vertices():
            self.old_to_new[self.dfs_number[v]] = self.newnum[v]

        # Update lowpt values according to `newnum` values.
        for v in self.graph_copy.vertices():
            self.node_at[self.newnum[v]] = v
            self.lowpt1[v] = self.old_to_new[self.lowpt1[v]]
            self.lowpt2[v] = self.old_to_new[self.lowpt2[v]]

    def __path_search(self, v):
        """
        Find the separation pairs and construct the split components.
        Check for type-1 and type-2 separation pairs, and construct
        the split components while also creating new virtual edges wherever
        required.
        """
        y = 0
        vnum = self.newnum[v]
        outv = self.adj[v].get_length()
        e_node = self.adj[v].get_head()
        while e_node:
            e = e_node.get_data()
            it = e_node

            if e in self.reverse_edges:
                w = e[0] # target
            else:
                w = e[1]
            wnum = self.newnum[w]
            if self.edge_status[e] == 1: # e is a tree arc
                if self.starts_path[e]: # if a new path starts at edge e
                    y = 0
                    # Pop all (h,a,b) from tstack where a > lowpt1[w]
                    if self.t_stack_a[self.t_stack_top] > self.lowpt1[w]:
                        while self.t_stack_a[self.t_stack_top] > self.lowpt1[w]:
                            y = max(y, self.t_stack_h[self.t_stack_top])
                            b = self.t_stack_b[self.t_stack_top]
                            self.t_stack_top -= 1
                        self.__tstack_push(y, self.lowpt1[w], b)

                    else:
                        self.__tstack_push(wnum + self.nd[w] - 1, self.lowpt1[w], vnum)
                    self.__tstack_push_eos()

                self.__path_search(w)

                self.e_stack.append(self.tree_arc[w])

                temp_node = self.adj[w].get_head()
                temp = temp_node.get_data()
                if temp in self.reverse_edges:
                    temp_target = temp[0]
                else:
                    temp_target = temp[1]

                # Type-2 separation pair check
                # while v is not the start_vertex
                while vnum != 1 and ((self.t_stack_a[self.t_stack_top] == vnum) or \
                        (self.degree[w] == 2 and self.newnum[temp_target] > wnum)):

                    a = self.t_stack_a[self.t_stack_top]
                    b = self.t_stack_b[self.t_stack_top]
                    e_virt = None
                    if a == vnum and self.parent[self.node_at[b]] == self.node_at[a]:
                        self.t_stack_top -= 1

                    else:
                        e_ab = None
                        if self.degree[w] == 2 and self.newnum[temp_target] > wnum:
                            # found type-2 separation pair - (v, temp_target)
                            e1 = self.__estack_pop()
                            e2 = self.__estack_pop()
                            self.adj[w].remove(self.in_adj[e2])

                            if e2 in self.reverse_edges:
                                x = e2[0] # target
                            else:
                                x = e2[1] # target

                            e_virt = tuple([v, x, "newVEdge"+str(self.virtual_edge_num)])
                            self.graph_copy.add_edge(e_virt)
                            self.virtual_edge_num += 1
                            self.degree[v] -= 1
                            self.degree[x] -= 1

                            if e2 in self.reverse_edges:
                                e2_source = e2[1] # target
                            else:
                                e2_source = e2[0]
                            if e2_source != w:
                                raise ValueError("Graph is not biconnected")

                            comp = Component([e1, e2, e_virt], 1)
                            self.components_list.append(comp)
                            comp = None

                            if self.e_stack:
                                e1 = self.e_stack[-1]
                                if e1 in self.reverse_edges:
                                    if e1[1] == x and e1[0] == v:
                                        e_ab = self.__estack_pop()
                                        self.adj[x].remove(self.in_adj[e_ab])
                                        self.__del_high(e_ab)
                                else:
                                    if e1[0] == x and e1[1] == v:
                                        e_ab = self.__estack_pop()
                                        self.adj[x].remove(self.in_adj[e_ab])
                                        self.__del_high(e_ab)

                        else: # found type-2 separation pair - (self.node_at[a], self.node_at[b])
                            h = self.t_stack_h[self.t_stack_top]
                            self.t_stack_top -= 1

                            comp = Component([],0)
                            while True:
                                xy = self.e_stack[-1]
                                if xy in self.reverse_edges:
                                    x = xy[1]
                                    xy_target = xy[0]
                                else:
                                    x = xy[0]
                                    xy_target = xy[1]
                                if not (a <= self.newnum[x] and self.newnum[x] <= h and \
                                    a <= self.newnum[xy_target] and self.newnum[xy_target] <= h):
                                    break
                                if (self.newnum[x] == a and self.newnum[xy_target] == b) or \
                                    (self.newnum[xy_target] == a and self.newnum[x] == b):
                                    e_ab = self.__estack_pop()
                                    if e_ab in self.reverse_edges:
                                        e_ab_source = e_ab[1] # source
                                    else:
                                        e_ab_source = e_ab[0] # source
                                    self.adj[e_ab_source].remove(self.in_adj[e_ab])
                                    self.__del_high(e_ab)

                                else:
                                    eh = self.__estack_pop()
                                    if eh in self.reverse_edges:
                                        eh_source = eh[1]
                                    else:
                                        eh_source = eh[0]
                                    if it != self.in_adj[eh]:
                                        self.adj[eh_source].remove(self.in_adj[eh])
                                        self.__del_high(eh)

                                    comp.add_edge(eh)
                                    self.degree[x] -= 1
                                    self.degree[xy_target] -= 1

                            e_virt = tuple([self.node_at[a], self.node_at[b], "newVEdge"+str(self.virtual_edge_num)])
                            self.graph_copy.add_edge(e_virt)
                            self.virtual_edge_num += 1
                            comp.finish_tric_or_poly(e_virt)
                            self.components_list.append(comp)
                            comp = None
                            x = self.node_at[b]

                        if e_ab is not None:
                            comp = Component([e_ab, e_virt], type_c=0)
                            e_virt = tuple([v, x, "newVEdge"+str(self.virtual_edge_num)])
                            self.graph_copy.add_edge(e_virt)
                            self.virtual_edge_num += 1
                            comp.add_edge(e_virt)
                            self.degree[x] -= 1
                            self.degree[v] -= 1
                            self.components_list.append(comp)
                            comp = None

                        self.e_stack.append(e_virt)
                        e_virt_node = LinkedListNode()
                        e_virt_node.set_data(e_virt)
                        # Replace `it` node with `e_virt_node`
                        it.replace(e_virt_node)
                        it = e_virt_node

                        self.in_adj[e_virt] = it
                        self.degree[x] += 1
                        self.degree[v] += 1
                        self.parent[x] = v
                        self.tree_arc[x] = e_virt
                        self.edge_status[e_virt] = 1
                        w = x
                        wnum = self.newnum[w]

                # start type-1 check
                if self.lowpt2[w] >= vnum and self.lowpt1[w] < vnum and \
                    (self.parent[v] != self.start_vertex or outv >= 2):
                    # type-1 separation pair - (self.node_at[self.lowpt1[w]], v)

                    # Create a new component and add edges to it
                    comp = Component([],0)
                    if not self.e_stack:
                        raise ValueError("stack is empty")
                    while self.e_stack:
                        xy = self.e_stack[-1]
                        if xy in self.reverse_edges:
                            xx = self.newnum[xy[1]] #source
                            y = self.newnum[xy[0]] #target
                        else:
                            xx = self.newnum[xy[0]] #source
                            y = self.newnum[xy[1]] #target

                        if not ((wnum <= xx and  xx < wnum + self.nd[w]) or \
                            (wnum <= y and y < wnum + self.nd[w])):
                            break

                        comp.add_edge(self.__estack_pop())
                        self.__del_high(xy)
                        self.degree[self.node_at[xx]] -= 1
                        self.degree[self.node_at[y]] -= 1

                    e_virt = tuple([v, self.node_at[self.lowpt1[w]], "newVEdge"+str(self.virtual_edge_num)])
                    self.graph_copy.add_edge(e_virt) # Add virtual edge to graph
                    self.virtual_edge_num += 1
                    comp.finish_tric_or_poly(e_virt) # Add virtual edge to component
                    self.components_list.append(comp)
                    comp = None

                    if (xx == vnum and y == self.lowpt1[w]) or \
                        (y == vnum and xx == self.lowpt1[w]):
                        comp_bond = Component([],type_c = 0) # new triple bond
                        eh = self.__estack_pop()
                        if self.in_adj[eh] != it:
                            if eh in self.reverse_edges:
                                self.adj[eh[1]].remove(self.in_adj[eh])
                            else:
                                self.adj[eh[0]].remove(self.in_adj[eh])

                        comp_bond.add_edge(eh)
                        comp_bond.add_edge(e_virt)
                        e_virt = tuple([v, self.node_at[self.lowpt1[w]], "newVEdge"+str(self.virtual_edge_num)])
                        self.graph_copy.add_edge(e_virt)
                        self.virtual_edge_num += 1
                        comp_bond.add_edge(e_virt)
                        self.in_high[e_virt] = self.in_high[eh]
                        self.degree[v] -= 1
                        self.degree[self.node_at[self.lowpt1[w]]] -= 1

                        self.components_list.append(comp_bond)
                        comp_bond = None
                    if self.node_at[self.lowpt1[w]] != self.parent[v]:
                        self.e_stack.append(e_virt)

                        e_virt_node = LinkedListNode()
                        e_virt_node.set_data(e_virt)
                        # replace `it` node with `e_virt_node`
                        it.replace(e_virt_node)
                        it = e_virt_node

                        self.in_adj[e_virt] = it
                        if not e_virt in self.in_high and self.__high(self.node_at[self.lowpt1[w]]) < vnum:
                            vnum_node = LinkedListNode()
                            vnum_node.set_data(vnum)
                            self.highpt[self.node_at[self.lowpt1[w]]].push_front(vnum_node)
                            self.in_high[e_virt] = vnum_node

                        self.degree[v] += 1
                        self.degree[self.node_at[self.lowpt1[w]]] += 1

                    else:
                        self.adj[v].remove(it)
                        comp_bond = Component([e_virt], type_c=0)
                        e_virt = tuple([self.node_at[self.lowpt1[w]], v, "newVEdge"+str(self.virtual_edge_num)])
                        self.graph_copy.add_edge(e_virt)
                        self.virtual_edge_num += 1
                        comp_bond.add_edge(e_virt)

                        eh = self.tree_arc[v];
                        comp_bond.add_edge(eh)

                        self.components_list.append(comp_bond)
                        comp_bond = None

                        self.tree_arc[v] = e_virt
                        self.edge_status[e_virt] = 1
                        self.in_adj[e_virt] = self.in_adj[eh]
                        e_virt_node = LinkedListNode()
                        e_virt_node.set_data(e_virt)
                        self.in_adj[eh] = e_virt_node
                        # end type-1 search

                # if an path starts at edge e, empty the tstack.
                if self.starts_path[e]:
                    while self.__tstack_not_eos():
                        self.t_stack_top -= 1
                    self.t_stack_top -= 1

                while self.__tstack_not_eos() and self.t_stack_b[self.t_stack_top] != vnum \
                    and self.__high(v) > self.t_stack_h[self.t_stack_top]:
                    self.t_stack_top -= 1

                outv -= 1

            else: # e is a frond
                if self.starts_path[e]:
                    y = 0
                    # pop all (h,a,b) from tstack where a > w
                    if self.t_stack_a[self.t_stack_top] > wnum:
                        while self.t_stack_a[self.t_stack_top] > wnum:
                            y = max(y, self.t_stack_h[self.t_stack_top])
                            b = self.t_stack_b[self.t_stack_top]
                            self.t_stack_top -= 1
                        self.__tstack_push(y, wnum, b)

                    else:
                        self.__tstack_push(vnum, wnum, vnum)
                self.e_stack.append(e) # add (v,w) to ESTACK

            # Go to next edge in adjacency list
            e_node = e_node.next

    def __assemble_triconnected_components(self):
        """
        Iterate through all the split components built by `path_finder` and
        merges two bonds or two polygons that share an edge for contructing the
        final triconnected components.
        Subsequently, convert the edges in triconnected components into original
        vertices and edges. The triconnected components are stored in
        `self.comp_list_new` and `self.comp_type`.
        """
        comp1 = {} # The index of first component that an edge belongs to
        comp2 = {} # The index of second component that an edge belongs to
        item1 = {} # Pointer to the edge node in component1
        item2 = {} # Pointer to the edge node in component2
        num_components = len(self.components_list)
        visited = [False for i in range(num_components)]

        # For each edge, we populate the comp1, comp2, item1 and item2 values
        for i in range(num_components): # for each component
            e_node = self.components_list[i].edge_list.get_head()
            while e_node: # for each edge
                e = e_node.get_data()
                if e not in item1:
                    comp1[e] = i
                    item1[e] = e_node
                else:
                    comp2[e] = i
                    item2[e] = e_node

                e_node = e_node.next

        # For each edge in a component, if the edge is a virtual edge, merge
        # the two components the edge belongs to
        for i in range(num_components):
            c1 = self.components_list[i]
            c1_type = c1.component_type
            l1 = c1.edge_list
            visited[i] = True

            if l1.get_length() == 0:
                continue

            if c1_type == 0 or c1_type == 1:
                e_node = self.components_list[i].edge_list.get_head()
                # Iterate through each edge in the component
                while e_node:
                    e = e_node.get_data()
                    e_node_next = e_node.next
                    # The label of a virtual edge is a string
                    if not isinstance(e[2], str):
                        e_node = e_node_next
                        continue

                    j = comp1[e]
                    if visited[j]:
                        j = comp2[e]
                        if visited[j]:
                            e_node = e_node_next
                            continue
                        e_node2 = item2[e]
                    else:
                        e_node2 = item1[e]

                    c2 = self.components_list[j]

                    # If the two components are not the same type, do not merge
                    if (c1_type != c2.component_type):
                        e_node = e_node_next # Go to next edge
                        continue

                    visited[j] = True
                    l2 = c2.edge_list

                    # Remove the corresponding virtual edges in both the components
                    # and merge the components
                    l2.remove(e_node2)
                    l1.concatenate(l2)

                    # if `e_node_next` was empty, after merging two components,
                    # more edges are added to the component.
                    if not e_node_next:
                        e_node_next = e_node.next # Go to next edge

                    l1.remove(e_node)

                    e_node = e_node_next

        # Convert connected components into original graph vertices and edges
        self.comp_list_new = []
        self.comp_type = []
        for i in range(len(self.components_list)):
            if self.components_list[i].edge_list.get_length() > 0:
                e_list = self.components_list[i].get_edge_list()
                e_list_new = []
                # For each edge, get the original source, target and label
                for e in e_list:
                    source = self.int_to_vertex[e[0]]
                    target = self.int_to_vertex[e[1]]
                    if isinstance(e[2], str):
                        label = e[2]
                    else:
                        label = self.edge_label_dict[(source, target,e[2])][2]
                    e_list_new.append(tuple([source, target, label]))
                # Add the component data to `comp_list_new` and `comp_type`
                self.comp_type.append(self.components_list[i].component_type)
                self.comp_list_new.append(e_list_new)

    def print_triconnected_components(self):
        """
        Print all the triconnected components along with the type of
        each component.
        """
        for i in range(len(self.comp_list_new)):
            if self.comp_type[i] == 0:
                print "Bond: ",
            elif self.comp_type[i] == 1:
                print "Polygon: ",
            else:
                print "Triconnected: ",
            print self.comp_list_new[i]

