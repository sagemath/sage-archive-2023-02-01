from libcpp cimport bool
from libcpp cimport string
from sage.graphs.graph import Graph
def FurerGadget(int k, prefix=None):
    r"""
    Return a Furer gadget of order ``k`` and their coloring.
    
    Construct the Furer gadget described for the first time in [CFI1992]_,
    a graph composed by a middle layer of `2^(k-1)` nodes and two sets of nodes `(a_0, ... , a_{k-1})` and `(b_0, ... , b_{k-1})`.
    Each node in the middle is connected to one of either `a_i` or `b_i`, for each i in [0,k[.
    To read about the complete construction, see the paper mentioned above.
    The returned coloring colors the middle section with one color, and then each pair `(a_i, b_i)` with another color.
    Since this method is mainly used to create Furer gadgets for the Cai-Furer-Immerman construction, returning gadgets that don't
    always have the same vertex labels is important, that's why there is a parameter to manually set a prefix to be appended to each vertex label.
    
    INPUT:

    - ``k`` -- The order of the returned FurerGadget

    - ``prefix`` -- Prefix to be appended to each vertex label, so as to individualise the returned Furer gadget
    
    OUTPUT:
    
    - ``G`` -- The Furer gadget of order ``k``
    
    - ``coloring`` -- A list of list of vertices, representing the partition induced by the coloring on the vertices of ``G``
    """
    from itertools import repeat as rep, chain, combinations
    cdef int i
    cdef tuple t
    cdef set s
    G = Graph()
    if prefix is not None:
        for i in range(k):
            G.add_vertex((prefix, i, 'a'))
            G.add_vertex((prefix, i, 'b'))
    else:
        for i in range(k):
            G.add_vertex((i, 'a'))
            G.add_vertex((i, 'b'))
    powerset = chain.from_iterable(combinations(range(k), r) for r in range(k+1) if r % 2 == 0)
    if prefix is not None:
        for t in powerset:
            s = set(t)
            for i in range(k):
                if i in s:
                    G.add_edge((prefix,t), (prefix,i,'a'))
                else: 
                    G.add_edge((prefix,t), (prefix,i,'b'))
    else:
        for t in powerset:
            s = set(t)
            for i in range(k):
                if i in s:
                    G.add_edge(t, (i,'a'))
                else: 
                    G.add_edge(t, (i,'b'))
    cdef list partition = []
    if prefix is not None:
        for i in range(k):
            partition.append([(prefix, i, 'a'),(prefix, i, 'b')])
    else:
        for i in range(k):
            partition.append([(i, 'a'),(i, 'b')])
    partition.append(powerset)
    return G, partition

def CaiFurerImmermanGraph(G, bool twisted=False):
    r"""
    Return the a Cai-Furer-Immerman graph from ``G``, possibly a twisted one, and a partition of its nodes.

    A Cai-Furer-Immerman graph from/on ``G`` is a graph created by applying the transformation
    described in [CFI1992]_ on a graph ``G``, that is substituting every vertex v in ``G`` with a
    Furer gadget F(v) of order d equal to the degree of the vertex, and then substituting every edge (v,u) in ``G`` 
    with a pair of edges, one connecting the two "a" nodes of F(v) and F(u) and the other their two "b" nodes.
    The returned coloring of the vertices is made by the union of the colorings of each single Furer gadget,
    individualised for each vertex of ``G``.
    To understand better what these "a" and "b" nodes are, see the documentation on  Furer gadgets.
    
    Furthermore, this method can apply what is described in the paper mentioned above as a "twist" on an edge,
    that is taking only one of the pairs of edges introduced in the new graph and swap two of their extremes,
    making each edge go from an "a" node to a "b" node. This is only doable if the original graph G is connected.
    
    A CaiFurerImmerman graph on a graph with no balanced vertex separators smaller than s and its twisted version
    cannot be distinguished by k-WL for any k < s.

    INPUT:

    - ``G`` -- An undirected graph on which to construct the Cai-Furer-Immerman graph

    - ``twisted`` -- A boolean indicating if the version to construct is a twisted one or not
    
    OUTPUT:
    
    - ``H`` -- The Cai-Furer-Immerman graph on ``G``
    
    - ``coloring`` -- A list of list of vertices, representing the partition induced by the coloring on the vertices of ``H``
    """
    cdef bool isConnected = G.is_connected()
    newG = Graph()
    cdef list total_partition = []
    cdef dict edge_index = {}
    cdef list p
    cdef int i,j
    cdef str s
    for v in G:
        Fk, p = FurerGadget(G.degree(v), v)
        total_partition += p
        newG=newG.union(Fk)
        edge_index[v] = 0
    for v,u in G.edges(labels=False):
        i = edge_index[v]
        edge_index[v] += 1
        j = edge_index[u]
        edge_index[u] += 1
        edge_va = (v, (i, 'a'))
        edge_vb = (v, (i, 'b'))
        edge_ua = (u, (j, 'a'))
        edge_ub = (u, (j, 'b'))
        if isConnected and twisted:
            temp = edge_ua
            edge_ua = edge_ub
            edge_ua = temp
            isConnected = False
        newG.add_edge(edge_va, edge_ua)
        newG.add_edge(edge_vb, edge_ub)
    if(twisted and G.is_connected()):
        s = " twisted"
    else:
        s = ""
    newG.name("CaiFurerImmerman" + s + " graph constructed from a " + G.name())
    return newG, total_partition

def EgawaGraph(int p, int s):
    r"""
    Returns the Egawa graph with parameters `p`, `s`.

    Egawa graphs are a peculiar family of graphs devised by Yoshimi Egawa in 
    [Ega1981]_ .
    The Shrikhande graph is a special case of this family of graphs, with parameters
    `(1,0)`.
    All the graphs in this family are not recognizable by 1-WL (Weisfeiler Lehamn algorithm of the first order) and 2-WL,
    that is their orbits are not correctly returned by k-WL for k lower than 3.
    
    Furthermore, all the graphs in this family are distance-regular, but they are
    not distance-transitive if `p \neq 0`.
    
    The Egawa graph with parameters `(0, s)` is isomorphic to the Hamming graph
    with parameters `(s, 4)`, when the underlying set of the Hamming graph is `[0,1,2,3]`
    
    INPUT:

    - ``p`` -- power to which the graph named `Y` in the reference provided above will be raised

    - ``s`` -- power to which the graph named `X` in the reference provided above will be raised
    
    OUTPUT:
    
    - ``G`` -- The Egawa graph with parameters (p,s)

    EXAMPLES:

    Every Egawa graph is distance regular.  ::

        sage: g = graphs.EgawaGraph(1, 2)
        sage: g.is_distance_regular()
        True

    An Egawa graph with parameters (0,s) is isomorphic to the Hamming graph
    with parameters (s, 4).  ::

        sage: g = graphs.EgawaGraph(0, 4)
        sage: g.is_isomorphic(graphs.HammingGraph(4,4))
        True
    """
    from sage.graphs.generators.basic import CompleteGraph
    from itertools import product, chain, repeat
    g = Graph(name="Egawa Graph with parameters " + str(p) + "," + str(s), multiedges=False)
    X = CompleteGraph(4)
    Y = Graph('O?Wse@UgqqT_LUebWkbT_')
    cdef tuple prefix, suffix, u, v
    cdef int i, el
    g.add_vertices(product(*chain(repeat(Y, p), repeat(X,s))))
    for v in g:
        for i in range(p):
            prefix = v[:i]
            suffix = v[i+1:]
            for el in Y.neighbor_iterator(v[i]):
                u = prefix + (el,) + suffix
                g.add_edge(v,u)
        for i in range(p, s+p):
            prefix = v[:i]
            suffix = v[i+1:]
            for el in X:
                if el == v[i]: continue
                u = prefix + (el,) + suffix
                g.add_edge(v,u)
    return g

def HammingGraph(int n, int q, X=None):
    r"""
    Returns the Hamming graph with parameters ``n``, ``q`` over set ``X``.

    Hamming graphs are graphs over the cartesian product of n copies of ``X``,
    where `q = |X|`, where the vertices, labelled with the corresponding tuple in `X^n`,
    are connected if the Hamming distance between their labels is 1.
    All Hamming graphs are regular, vertex-transitive and distance-regular.
    
    Hamming graphs with parameters `(1,q)` represent the complete graph with
    q vertices over the set ``X``.
    
    INPUT:

    - ``n`` -- power to which ``X`` will be raised to provide vertices for the Hamming graph

    - ``q`` -- cardinality of ``X``

    - ``X`` -- list of labels representing the vertices of the underlying graph the Hamming graph will be based on;
               if None (or left unused), the list `[0, ... , q-1]` will be used
    
    OUTPUT:
    
    - ``G`` -- The Hamming graph with parameters `(n,q,X)`

    EXAMPLES:

    Every Hamming graph is distance-regular, regular and vertex-transitive.  ::

        sage: g = graphs.HammingGraph(3, 7)
        sage: g.is_distance_regular()
        True
        sage: g.is_regular() 
        True
        sage: g.is_vertex_transitive()
        True

    A Hamming graph with parameters (1,q) is isomorphic to the Complete graph
    with parameter q.  ::

        sage: g = graphs.HammingGraph(1, 23)
        sage: g.is_isomorphic(graphs.CompleteGraph(23))
        True
    
    If a parameter ``q`` is provided which is not equal to ``X``'s cardinality, an exception is raised. ::
    
        sage: X = ['a','b','c','d','e']
        sage: g = graphs.HammingGraph(2, 3, X)
        Traceback (most recent call last):
        ...
        ValueError: q must be the cardinality of X
        
    REFERENCES:
    
    For a more accurate description, see the following wikipedia page:
    :wikipedia:`Hamming_graph`
    """
    from itertools import product, repeat
    cdef int i
    cdef tuple v, prefix, suffix, u
    if not X:
        X = list(range(q))
    if q != len(X):
        raise ValueError("q must be the cardinality of X")
    g = Graph(name="Hamming Graph with parameters " + str(n) + "," + str(q), multiedges=False)
    g.add_vertices(product(*repeat(X, n)))
    for v in g:
        for i in range(n):
            prefix = v[:i]
            suffix = v[i+1:]
            for el in X:
                if el == v[i]: continue
                u = prefix + (el,) + suffix
                g.add_edge(v,u)
    return g
