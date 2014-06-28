# -*- coding: utf-8 -*-
r"""
Families of graphs

The methods defined here appear in :mod:`sage.graphs.graph_generators`.

AUTHORS:

- David Coudert (2012) Ringed Trees

"""

###########################################################################
#
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################

# import from Sage library
from sage.graphs.graph import Graph
from sage.graphs import graph
from math import sin, cos, pi

def JohnsonGraph(n, k):
    r"""
    Returns the Johnson graph with parameters `n, k`.

    Johnson graphs are a special class of undirected graphs defined from systems
    of sets. The vertices of the Johnson graph `J(n,k)` are the `k`-element
    subsets of an `n`-element set; two vertices are adjacent when they meet in a
    `(k-1)`-element set. For more information about Johnson graphs, see the
    corresponding :wikipedia:`Wikipedia page <Johnson_graph>`.

    EXAMPLES:

    The Johnson graph is a Hamiltonian graph.  ::

        sage: g = graphs.JohnsonGraph(7, 3)
        sage: g.is_hamiltonian()
        True

    Every Johnson graph is vertex transitive.  ::

        sage: g = graphs.JohnsonGraph(6, 4)
        sage: g.is_vertex_transitive()
        True

    The complement of the Johnson graph `J(n,2)` is isomorphic to the Knesser
    Graph `K(n,2)`.  In paritcular the complement of `J(5,2)` is isomorphic to
    the Petersen graph.  ::

        sage: g = graphs.JohnsonGraph(5,2)
        sage: g.complement().is_isomorphic(graphs.PetersenGraph())
        True
    """

    g = Graph(name="Johnson graph with parameters "+str(n)+","+str(k))
    from sage.combinat.subset import Set, Subsets

    S = Set(range(n))
    g.add_vertices(Subsets(S, k))

    for sub in Subsets(S, k-1):
        elem_left = S - sub
        for i in elem_left:
            for j in elem_left:
                if j <= i:
                    continue
                g.add_edge(sub+Set([i]),sub+Set([j]))

    return g


def KneserGraph(n,k):
    r"""
    Returns the Kneser Graph with parameters `n, k`.

    The Kneser Graph with parameters `n,k` is the graph
    whose vertices are the `k`-subsets of `[0,1,\dots,n-1]`, and such
    that two vertices are adjacent if their corresponding sets
    are disjoint.

    For example, the Petersen Graph can be defined
    as the Kneser Graph with parameters `5,2`.

    EXAMPLE::

        sage: KG=graphs.KneserGraph(5,2)
        sage: print KG.vertices()
        [{4, 5}, {1, 3}, {2, 5}, {2, 3}, {3, 4}, {3, 5}, {1, 4}, {1, 5}, {1, 2}, {2, 4}]
        sage: P=graphs.PetersenGraph()
        sage: P.is_isomorphic(KG)
        True

    TESTS::

        sage: KG=graphs.KneserGraph(0,0)
        Traceback (most recent call last):
        ...
        ValueError: Parameter n should be a strictly positive integer
        sage: KG=graphs.KneserGraph(5,6)
        Traceback (most recent call last):
        ...
        ValueError: Parameter k should be a strictly positive integer inferior to n
    """

    if not n>0:
        raise ValueError("Parameter n should be a strictly positive integer")
    if not (k>0 and k<=n):
        raise ValueError("Parameter k should be a strictly positive integer inferior to n")

    g = Graph(name="Kneser graph with parameters {},{}".format(n,k))

    from sage.combinat.subset import Subsets
    S = Subsets(n,k)
    if k>n/2:
        g.add_vertices(S)

    s0 = S.underlying_set()    # {1,2,...,n}
    for s in S:
        for t in Subsets(s0.difference(s), k):
            g.add_edge(s,t)

    return g

def BalancedTree(r, h):
    r"""
    Returns the perfectly balanced tree of height `h \geq 1`,
    whose root has degree `r \geq 2`.

    The number of vertices of this graph is
    `1 + r + r^2 + \cdots + r^h`, that is,
    `\frac{r^{h+1} - 1}{r - 1}`. The number of edges is one
    less than the number of vertices.

    INPUT:

    - ``r`` -- positive integer `\geq 2`. The degree of the root node.

    - ``h`` -- positive integer `\geq 1`. The height of the balanced tree.

    OUTPUT:

    The perfectly balanced tree of height `h \geq 1` and whose root has
    degree `r \geq 2`. A ``NetworkXError`` is returned if `r < 2` or
    `h < 1`.

    ALGORITHM:

    Uses `NetworkX <http://networkx.lanl.gov>`_.

    EXAMPLES:

    A balanced tree whose root node has degree `r = 2`, and of height
    `h = 1`, has order 3 and size 2::

        sage: G = graphs.BalancedTree(2, 1); G
        Balanced tree: Graph on 3 vertices
        sage: G.order(); G.size()
        3
        2
        sage: r = 2; h = 1
        sage: v = 1 + r
        sage: v; v - 1
        3
        2

    Plot a balanced tree of height 5, whose root node has degree `r = 3`::

        sage: G = graphs.BalancedTree(3, 5)
        sage: G.show()   # long time

    A tree is bipartite. If its vertex set is finite, then it is planar. ::

        sage: r = randint(2, 5); h = randint(1, 7)
        sage: T = graphs.BalancedTree(r, h)
        sage: T.is_bipartite()
        True
        sage: T.is_planar()
        True
        sage: v = (r^(h + 1) - 1) / (r - 1)
        sage: T.order() == v
        True
        sage: T.size() == v - 1
        True

    TESTS:

     Normally we would only consider balanced trees whose root node
     has degree `r \geq 2`, but the construction degenerates
     gracefully::

        sage: graphs.BalancedTree(1, 10)
        Balanced tree: Graph on 2 vertices

        sage: graphs.BalancedTree(-1, 10)
        Balanced tree: Graph on 1 vertex

    Similarly, we usually want the tree must have height `h \geq 1`
    but the algorithm also degenerates gracefully here::

        sage: graphs.BalancedTree(3, 0)
        Balanced tree: Graph on 1 vertex

        sage: graphs.BalancedTree(5, -2)
        Balanced tree: Graph on 0 vertices

        sage: graphs.BalancedTree(-2,-2)
        Balanced tree: Graph on 0 vertices
    """
    import networkx
    return Graph(networkx.balanced_tree(r, h), name="Balanced tree")

def BarbellGraph(n1, n2):
    r"""
    Returns a barbell graph with ``2*n1 + n2`` nodes. The argument ``n1``
    must be greater than or equal to 2.

    A barbell graph is a basic structure that consists of a path graph
    of order ``n2`` connecting two complete graphs of order ``n1`` each.

    This constructor depends on `NetworkX <http://networkx.lanl.gov>`_
    numeric labels. In this case, the ``n1``-th node connects to the
    path graph from one complete graph and the ``n1 + n2 + 1``-th node
    connects to the path graph from the other complete graph.

    INPUT:

    - ``n1`` -- integer `\geq 2`. The order of each of the two
      complete graphs.

    - ``n2`` -- nonnegative integer. The order of the path graph
      connecting the two complete graphs.

    OUTPUT:

    A barbell graph of order ``2*n1 + n2``. A ``ValueError`` is
    returned if ``n1 < 2`` or ``n2 < 0``.

    ALGORITHM:

    Uses `NetworkX <http://networkx.lanl.gov>`_.

    PLOTTING:

    Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, each barbell
    graph will be displayed with the two complete graphs in the
    lower-left and upper-right corners, with the path graph connecting
    diagonally between the two. Thus the ``n1``-th node will be drawn at a
    45 degree angle from the horizontal right center of the first
    complete graph, and the ``n1 + n2 + 1``-th node will be drawn 45
    degrees below the left horizontal center of the second complete graph.

    EXAMPLES:

    Construct and show a barbell graph ``Bar = 4``, ``Bells = 9``::

        sage: g = graphs.BarbellGraph(9, 4); g
        Barbell graph: Graph on 22 vertices
        sage: g.show() # long time

    An ``n1 >= 2``, ``n2 >= 0`` barbell graph has order ``2*n1 + n2``. It
    has the complete graph on ``n1`` vertices as a subgraph. It also has
    the path graph on ``n2`` vertices as a subgraph. ::

        sage: n1 = randint(2, 2*10^2)
        sage: n2 = randint(0, 2*10^2)
        sage: g = graphs.BarbellGraph(n1, n2)
        sage: v = 2*n1 + n2
        sage: g.order() == v
        True
        sage: K_n1 = graphs.CompleteGraph(n1)
        sage: P_n2 = graphs.PathGraph(n2)
        sage: s_K = g.subgraph_search(K_n1, induced=True)
        sage: s_P = g.subgraph_search(P_n2, induced=True)
        sage: K_n1.is_isomorphic(s_K)
        True
        sage: P_n2.is_isomorphic(s_P)
        True

    Create several barbell graphs in a Sage graphics array::

        sage: g = []
        sage: j = []
        sage: for i in range(6):
        ...       k = graphs.BarbellGraph(i + 2, 4)
        ...       g.append(k)
        ...
        sage: for i in range(2):
        ...       n = []
        ...       for m in range(3):
        ...           n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ...       j.append(n)
        ...
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time

    TESTS:

    The input ``n1`` must be `\geq 2`::

        sage: graphs.BarbellGraph(1, randint(0, 10^6))
        Traceback (most recent call last):
        ...
        ValueError: Invalid graph description, n1 should be >= 2
        sage: graphs.BarbellGraph(randint(-10^6, 1), randint(0, 10^6))
        Traceback (most recent call last):
        ...
        ValueError: Invalid graph description, n1 should be >= 2

    The input ``n2`` must be `\geq 0`::

        sage: graphs.BarbellGraph(randint(2, 10^6), -1)
        Traceback (most recent call last):
        ...
        ValueError: Invalid graph description, n2 should be >= 0
        sage: graphs.BarbellGraph(randint(2, 10^6), randint(-10^6, -1))
        Traceback (most recent call last):
        ...
        ValueError: Invalid graph description, n2 should be >= 0
        sage: graphs.BarbellGraph(randint(-10^6, 1), randint(-10^6, -1))
        Traceback (most recent call last):
        ...
        ValueError: Invalid graph description, n1 should be >= 2
    """
    # sanity checks
    if n1 < 2:
        raise ValueError("Invalid graph description, n1 should be >= 2")
    if n2 < 0:
        raise ValueError("Invalid graph description, n2 should be >= 0")

    pos_dict = {}

    for i in range(n1):
        x = float(cos((pi / 4) - ((2 * pi) / n1) * i) - (n2 / 2) - 1)
        y = float(sin((pi / 4) - ((2 * pi) / n1) * i) - (n2 / 2) - 1)
        j = n1 - 1 - i
        pos_dict[j] = (x, y)
    for i in range(n1, n1 + n2):
        x = float(i - n1 - (n2 / 2) + 1)
        y = float(i - n1 - (n2 / 2) + 1)
        pos_dict[i] = (x, y)
    for i in range(n1 + n2, (2 * n1) + n2):
        x = float(
            cos((5 * (pi / 4)) + ((2 * pi) / n1) * (i - n1 - n2))
            + (n2 / 2) + 2)
        y = float(
            sin((5 * (pi / 4)) + ((2 * pi) / n1) * (i - n1 - n2))
            + (n2 / 2) + 2)
        pos_dict[i] = (x, y)

    import networkx
    G = networkx.barbell_graph(n1, n2)
    return Graph(G, pos=pos_dict, name="Barbell graph")

def BubbleSortGraph(n):
    r"""
    Returns the bubble sort graph `B(n)`.

    The vertices of the bubble sort graph are the set of permutations on
    `n` symbols. Two vertices are adjacent if one can be obtained from the
    other by swapping the labels in the `i`-th and `(i+1)`-th positions for
    `1 \leq i \leq n-1`. In total, `B(n)` has order `n!`. Thus, the order
    of `B(n)` increases according to `f(n) = n!`.

    INPUT:

    - ``n`` -- positive integer. The number of symbols to permute.

    OUTPUT:

    The bubble sort graph `B(n)` on `n` symbols. If `n < 1`, a
    ``ValueError`` is returned.

    EXAMPLES::

        sage: g = graphs.BubbleSortGraph(4); g
        Bubble sort: Graph on 24 vertices
        sage: g.plot() # long time

    The bubble sort graph on `n = 1` symbol is the trivial graph `K_1`::

        sage: graphs.BubbleSortGraph(1)
        Bubble sort: Graph on 1 vertex

    If `n \geq 1`, then the order of `B(n)` is `n!`::

        sage: n = randint(1, 8)
        sage: g = graphs.BubbleSortGraph(n)
        sage: g.order() == factorial(n)
        True

    TESTS:

    Input ``n`` must be positive::

        sage: graphs.BubbleSortGraph(0)
        Traceback (most recent call last):
        ...
        ValueError: Invalid number of symbols to permute, n should be >= 1
        sage: graphs.BubbleSortGraph(randint(-10^6, 0))
        Traceback (most recent call last):
        ...
        ValueError: Invalid number of symbols to permute, n should be >= 1

    AUTHORS:

    - Michael Yurko (2009-09-01)
    """
    # sanity checks
    if n < 1:
        raise ValueError(
            "Invalid number of symbols to permute, n should be >= 1")
    if n == 1:
        from sage.graphs.generators.basic import CompleteGraph
        return Graph(CompleteGraph(n), name="Bubble sort")
    from sage.combinat.permutation import Permutations
    #create set from which to permute
    label_set = [str(i) for i in xrange(1, n + 1)]
    d = {}
    #iterate through all vertices
    for v in Permutations(label_set):
        v = list(v) # So we can easily mutate it
        tmp_dict = {}
        #add all adjacencies
        for i in xrange(n - 1):
            #swap entries
            v[i], v[i + 1] = v[i + 1], v[i]
            #add new vertex
            new_vert = ''.join(v)
            tmp_dict[new_vert] = None
            #swap back
            v[i], v[i + 1] = v[i + 1], v[i]
        #add adjacency dict
        d[''.join(v)] = tmp_dict
    return Graph(d, name="Bubble sort")

def CirculantGraph(n, adjacency):
    r"""
    Returns a circulant graph with n nodes.

    A circulant graph has the property that the vertex `i` is connected
    with the vertices `i+j` and `i-j` for each j in adj.

    INPUT:


    -  ``n`` - number of vertices in the graph

    -  ``adjacency`` - the list of j values


    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, each circulant
    graph will be displayed with the first (0) node at the top, with
    the rest following in a counterclockwise manner.

    Filling the position dictionary in advance adds O(n) to the
    constructor.

    .. SEEALSO::

        * :meth:`sage.graphs.generic_graph.GenericGraph.is_circulant`
          -- checks whether a (di)graph is circulant, and/or returns
          all possible sets of parameters.

    EXAMPLES: Compare plotting using the predefined layout and
    networkx::

        sage: import networkx
        sage: n = networkx.cycle_graph(23)
        sage: spring23 = Graph(n)
        sage: posdict23 = graphs.CirculantGraph(23,2)
        sage: spring23.show() # long time
        sage: posdict23.show() # long time

    We next view many cycle graphs as a Sage graphics array. First we
    use the ``CirculantGraph`` constructor, which fills in
    the position dictionary::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ...    k = graphs.CirculantGraph(i+3,i)
        ...    g.append(k)
        ...
        sage: for i in range(3):
        ...    n = []
        ...    for m in range(3):
        ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ...    j.append(n)
        ...
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time

    Compare to plotting with the spring-layout algorithm::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ...    spr = networkx.cycle_graph(i+3)
        ...    k = Graph(spr)
        ...    g.append(k)
        ...
        sage: for i in range(3):
        ...    n = []
        ...    for m in range(3):
        ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ...    j.append(n)
        ...
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time

    Passing a 1 into adjacency should give the cycle.

    ::

        sage: graphs.CirculantGraph(6,1)==graphs.CycleGraph(6)
        True
        sage: graphs.CirculantGraph(7,[1,3]).edges(labels=false)
        [(0, 1),
        (0, 3),
        (0, 4),
        (0, 6),
        (1, 2),
        (1, 4),
        (1, 5),
        (2, 3),
        (2, 5),
        (2, 6),
        (3, 4),
        (3, 6),
        (4, 5),
        (5, 6)]
    """
    from sage.graphs.graph_plot import _circle_embedding

    if not isinstance(adjacency,list):
        adjacency=[adjacency]

    G = Graph(n, name="Circulant graph ("+str(adjacency)+")")
    _circle_embedding(G, range(n))

    for v in G:
        G.add_edges([(v,(v+j)%n) for j in adjacency])

    return G

def CubeGraph(n):
    r"""
    Returns the hypercube in `n` dimensions.

    The hypercube in `n` dimension is build upon the binary
    strings on `n` bits, two of them being adjacent if
    they differ in exactly one bit. Hence, the distance
    between two vertices in the hypercube is the Hamming
    distance.

    EXAMPLES:

    The distance between `0100110` and `1011010` is
    `5`, as expected ::

        sage: g = graphs.CubeGraph(7)
        sage: g.distance('0100110','1011010')
        5

    Plot several `n`-cubes in a Sage Graphics Array ::

        sage: g = []
        sage: j = []
        sage: for i in range(6):
        ...    k = graphs.CubeGraph(i+1)
        ...    g.append(k)
        ...
        sage: for i in range(2):
        ...    n = []
        ...    for m in range(3):
        ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ...    j.append(n)
        ...
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show(figsize=[6,4]) # long time

    Use the plot options to display larger `n`-cubes

    ::

        sage: g = graphs.CubeGraph(9)
        sage: g.show(figsize=[12,12],vertex_labels=False, vertex_size=20) # long time

    AUTHORS:

    - Robert Miller
    """
    theta = float(pi/n)

    d = {'':[]}
    dn={}
    p = {'':(float(0),float(0))}
    pn={}

    # construct recursively the adjacency dict and the positions
    for i in xrange(n):
        ci = float(cos(i*theta))
        si = float(sin(i*theta))
        for v,e in d.iteritems():
            v0 = v+'0'
            v1 = v+'1'
            l0 = [v1]
            l1 = [v0]
            for m in e:
                l0.append(m+'0')
                l1.append(m+'1')
            dn[v0] = l0
            dn[v1] = l1
            x,y = p[v]
            pn[v0] = (x, y)
            pn[v1] = (x+ci, y+si)
        d,dn = dn,{}
        p,pn = pn,{}

    # construct the graph
    r = Graph(name="%d-Cube"%n)
    r.add_vertices(d.keys())
    for u,L in d.iteritems():
        for v in L:
            r.add_edge(u,v)
    r.set_pos(p)

    return r

def DorogovtsevGoltsevMendesGraph(n):
    """
    Construct the n-th generation of the Dorogovtsev-Goltsev-Mendes
    graph.

    EXAMPLE::

        sage: G = graphs.DorogovtsevGoltsevMendesGraph(8)
        sage: G.size()
        6561

    REFERENCE:

    - [1] Dorogovtsev, S. N., Goltsev, A. V., and Mendes, J.
      F. F., Pseudofractal scale-free web, Phys. Rev. E 066122
      (2002).
    """
    import networkx
    return Graph(networkx.dorogovtsev_goltsev_mendes_graph(n),\
           name="Dorogovtsev-Goltsev-Mendes Graph, %d-th generation"%n)

def FoldedCubeGraph(n):
    r"""
    Returns the folded cube graph of order `2^{n-1}`.

    The folded cube graph on `2^{n-1}` vertices can be obtained from a cube
    graph on `2^n` vertices by merging together opposed
    vertices. Alternatively, it can be obtained from a cube graph on
    `2^{n-1}` vertices by adding an edge between opposed vertices. This
    second construction is the one produced by this method.

    For more information on folded cube graphs, see the corresponding
    :wikipedia:`Wikipedia page <Folded_cube_graph>`.

    EXAMPLES:

    The folded cube graph of order five is the Clebsch graph::

        sage: fc = graphs.FoldedCubeGraph(5)
        sage: clebsch = graphs.ClebschGraph()
        sage: fc.is_isomorphic(clebsch)
        True
    """

    if n < 1:
        raise ValueError("The value of n must be at least 2")

    g = CubeGraph(n-1)
    g.name("Folded Cube Graph")

    # Complementing the binary word
    def complement(x):
        x = x.replace('0','a')
        x = x.replace('1','0')
        x = x.replace('a','1')
        return x

    for x in g:
        if x[0] == '0':
            g.add_edge(x,complement(x))

    return g


def FriendshipGraph(n):
    r"""
    Returns the friendship graph `F_n`.

    The friendship graph is also known as the Dutch windmill graph. Let
    `C_3` be the cycle graph on 3 vertices. Then `F_n` is constructed by
    joining `n \geq 1` copies of `C_3` at a common vertex. If `n = 1`,
    then `F_1` is isomorphic to `C_3` (the triangle graph). If `n = 2`,
    then `F_2` is the butterfly graph, otherwise known as the bowtie
    graph. For more information, see this
    `Wikipedia article on the friendship graph <http://en.wikipedia.org/wiki/Friendship_graph>`_.

    INPUT:

    - ``n`` -- positive integer; the number of copies of `C_3` to use in
      constructing `F_n`.

    OUTPUT:

    - The friendship graph `F_n` obtained from `n` copies of the cycle
      graph `C_3`.

    .. seealso::

        - :meth:`GraphGenerators.ButterflyGraph`

    EXAMPLES:

    The first few friendship graphs. ::

        sage: A = []; B = []
        sage: for i in range(9):
        ...       g = graphs.FriendshipGraph(i + 1)
        ...       A.append(g)
        sage: for i in range(3):
        ...       n = []
        ...       for j in range(3):
        ...           n.append(A[3*i + j].plot(vertex_size=20, vertex_labels=False))
        ...       B.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(B)
        sage: G.show()  # long time

    For `n = 1`, the friendship graph `F_1` is isomorphic to the cycle
    graph `C_3`, whose visual representation is a triangle. ::

        sage: G = graphs.FriendshipGraph(1); G
        Friendship graph: Graph on 3 vertices
        sage: G.show()  # long time
        sage: G.is_isomorphic(graphs.CycleGraph(3))
        True

    For `n = 2`, the friendship graph `F_2` is isomorphic to the
    butterfly graph, otherwise known as the bowtie graph. ::

        sage: G = graphs.FriendshipGraph(2); G
        Friendship graph: Graph on 5 vertices
        sage: G.is_isomorphic(graphs.ButterflyGraph())
        True

    If `n \geq 1`, then the friendship graph `F_n` has `2n + 1` vertices
    and `3n` edges. It has radius 1, diameter 2, girth 3, and
    chromatic number 3. Furthermore, `F_n` is planar and Eulerian. ::

        sage: n = randint(1, 10^3)
        sage: G = graphs.FriendshipGraph(n)
        sage: G.order() == 2*n + 1
        True
        sage: G.size() == 3*n
        True
        sage: G.radius()
        1
        sage: G.diameter()
        2
        sage: G.girth()
        3
        sage: G.chromatic_number()
        3
        sage: G.is_planar()
        True
        sage: G.is_eulerian()
        True

    TESTS:

    The input ``n`` must be a positive integer. ::

        sage: graphs.FriendshipGraph(randint(-10^5, 0))
        Traceback (most recent call last):
        ...
        ValueError: n must be a positive integer
    """
    # sanity checks
    if n < 1:
        raise ValueError("n must be a positive integer")
    # construct the friendship graph
    if n == 1:
        from sage.graphs.generators.basic import CycleGraph
        G = CycleGraph(3)
        G.name("Friendship graph")
        return G
    # build the edge and position dictionaries
    from sage.functions.trig import cos, sin
    from sage.rings.real_mpfr import RR
    from sage.symbolic.constants import pi
    N = 2*n + 1           # order of F_n
    d = (2*pi) / (N - 1)  # angle between external nodes
    edge_dict = {}
    pos_dict = {}
    for i in range(N - 2):
        if i & 1:  # odd numbered node
            edge_dict.setdefault(i, [i + 1, N - 1])
        else:      # even numbered node
            edge_dict.setdefault(i, [N - 1])
        pos_dict.setdefault(i, [RR(cos(i*d)), RR(sin(i*d))])
    edge_dict.setdefault(N - 2, [0, N - 1])
    pos_dict.setdefault(N - 2, [RR(cos(d * (N-2))), RR(sin(d * (N-2)))])
    pos_dict.setdefault(N - 1, [0, 0])
    return Graph(edge_dict, pos=pos_dict, name="Friendship graph")

def FuzzyBallGraph(partition, q):
    r"""
    Construct a Fuzzy Ball graph with the integer partition
    ``partition`` and ``q`` extra vertices.

    Let `q` be an integer and let `m_1,m_2,...,m_k` be a set of positive
    integers.  Let `n=q+m_1+...+m_k`.  The Fuzzy Ball graph with partition
    `m_1,m_2,...,m_k` and `q` extra vertices is the graph constructed from the
    graph `G=K_n` by attaching, for each `i=1,2,...,k`, a new vertex `a_i` to
    `m_i` distinct vertices of `G`.

    For given positive integers `k` and `m` and nonnegative
    integer `q`, the set of graphs ``FuzzyBallGraph(p, q)`` for
    all partitions `p` of `m` with `k` parts are cospectral with
    respect to the normalized Laplacian.

    EXAMPLES::

        sage: graphs.FuzzyBallGraph([3,1],2).adjacency_matrix()
        [0 1 1 1 1 1 1 0]
        [1 0 1 1 1 1 1 0]
        [1 1 0 1 1 1 1 0]
        [1 1 1 0 1 1 0 1]
        [1 1 1 1 0 1 0 0]
        [1 1 1 1 1 0 0 0]
        [1 1 1 0 0 0 0 0]
        [0 0 0 1 0 0 0 0]


    Pick positive integers `m` and `k` and a nonnegative integer `q`.
    All the FuzzyBallGraphs constructed from partitions of `m` with
    `k` parts should be cospectral with respect to the normalized
    Laplacian::

        sage: m=4; q=2; k=2
        sage: g_list=[graphs.FuzzyBallGraph(p,q) for p in Partitions(m, length=k)]
        sage: set([g.laplacian_matrix(normalized=True).charpoly() for g in g_list])  # long time (7s on sage.math, 2011)
        set([x^8 - 8*x^7 + 4079/150*x^6 - 68689/1350*x^5 + 610783/10800*x^4 - 120877/3240*x^3 + 1351/100*x^2 - 931/450*x])
    """
    from sage.graphs.generators.basic import CompleteGraph
    if len(partition)<1:
        raise ValueError("partition must be a nonempty list of positive integers")
    n=q+sum(partition)
    g=CompleteGraph(n)
    curr_vertex=0
    for e,p in enumerate(partition):
        g.add_edges([(curr_vertex+i, 'a{0}'.format(e+1)) for i in range(p)])
        curr_vertex+=p
    return g

def FibonacciTree(n):
    r"""
    Returns the graph of the Fibonacci Tree `F_{i}` of order `n`.
    `F_{i}` is recursively defined as the a tree with a root vertex
    and two attached child trees `F_{i-1}` and `F_{i-2}`, where
    `F_{1}` is just one vertex and `F_{0}` is empty.

    INPUT:

    - ``n`` - the recursion depth of the Fibonacci Tree

    EXAMPLES::

        sage: g = graphs.FibonacciTree(3)
        sage: g.is_tree()
        True

    ::

        sage: l1 = [ len(graphs.FibonacciTree(_)) + 1 for _ in range(6) ]
        sage: l2 = list(fibonacci_sequence(2,8))
        sage: l1 == l2
        True

    AUTHORS:

    - Harald Schilly and Yann Laigle-Chapuy (2010-03-25)
    """
    T = Graph(name="Fibonacci-Tree-%d"%n)
    if n == 1: T.add_vertex(0)
    if n < 2: return T

    from sage.combinat.combinat import fibonacci_sequence
    F = list(fibonacci_sequence(n + 2))
    s = 1.618 ** (n / 1.618 - 1.618)
    pos = {}

    def fib(level, node, y):
        pos[node] = (node, y)
        if level < 2: return
        level -= 1
        y -= s
        diff = F[level]
        T.add_edge(node, node - diff)
        if level == 1: # only one child
            pos[node - diff] = (node, y)
            return
        T.add_edge(node, node + diff)
        fib(level, node - diff, y)
        fib(level - 1, node + diff, y)

    T.add_vertices(xrange(sum(F[:-1])))
    fib(n, F[n + 1] - 1, 0)
    T.set_pos(pos)

    return T

def GeneralizedPetersenGraph(n,k):
    r"""
    Returns a generalized Petersen graph with `2n` nodes. The variables
    `n`, `k` are integers such that `n>2` and `0<k\leq\lfloor(n-1)`/`2\rfloor`

    For `k=1` the result is a graph isomorphic to the circular ladder graph
    with the same `n`. The regular Petersen Graph has `n=5` and `k=2`.
    Other named graphs that can be described using this notation include
    the Desargues graph and the Moebius-Kantor graph.

    INPUT:

    - ``n`` - the number of nodes is `2*n`.

    - ``k`` - integer `0<k\leq\lfloor(n-1)`/`2\rfloor`. Decides
      how inner vertices are connected.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the generalized
    Petersen graphs are displayed as an inner and outer cycle pair, with
    the first n nodes drawn on the outer circle. The first (0) node is
    drawn at the top of the outer-circle, moving counterclockwise after that.
    The inner circle is drawn with the (n)th node at the top, then
    counterclockwise as well.

    EXAMPLES: For `k=1` the resulting graph will be isomorphic to a circular
    ladder graph. ::

        sage: g = graphs.GeneralizedPetersenGraph(13,1)
        sage: g2 = graphs.CircularLadderGraph(13)
        sage: g.is_isomorphic(g2)
        True

    The Desargues graph::

        sage: g = graphs.GeneralizedPetersenGraph(10,3)
        sage: g.girth()
        6
        sage: g.is_bipartite()
        True

    AUTHORS:

    - Anders Jonsson (2009-10-15)
    """
    if (n < 3):
            raise ValueError("n must be larger than 2")
    if (k < 1 or k>((n-1)/2)):
            raise ValueError("k must be in 1<= k <=floor((n-1)/2)")
    pos_dict = {}
    G = Graph()
    for i in range(n):
        x = float(cos((pi/2) + ((2*pi)/n)*i))
        y = float(sin((pi/2) + ((2*pi)/n)*i))
        pos_dict[i] = (x,y)
    for i in range(n, 2*n):
        x = float(0.5*cos((pi/2) + ((2*pi)/n)*i))
        y = float(0.5*sin((pi/2) + ((2*pi)/n)*i))
        pos_dict[i] = (x,y)
    for i in range(n):
        G.add_edge(i, (i+1) % n)
        G.add_edge(i, i+n)
        G.add_edge(i+n, n + (i+k) % n)
    return Graph(G, pos=pos_dict, name="Generalized Petersen graph (n="+str(n)+",k="+str(k)+")")

def HararyGraph( k, n ):
    r"""
    Returns the Harary graph on `n` vertices and connectivity `k`, where
    `2 \leq k < n`.

    A `k`-connected graph `G` on `n` vertices requires the minimum degree
    `\delta(G)\geq k`, so the minimum number of edges `G` should have is
    `\lceil kn/2\rceil`. Harary graphs achieve this lower bound, that is,
    Harary graphs are minimal `k`-connected graphs on `n` vertices.

    The construction provided uses the method CirculantGraph.  For more
    details, see the book D. B. West, Introduction to Graph Theory, 2nd
    Edition, Prentice Hall, 2001, p. 150--151; or the `MathWorld article on
    Harary graphs <http://mathworld.wolfram.com/HararyGraph.html>`_.

    EXAMPLES:

    Harary graphs `H_{k,n}`::

        sage: h = graphs.HararyGraph(5,9); h
        Harary graph 5, 9: Graph on 9 vertices
        sage: h.order()
        9
        sage: h.size()
        23
        sage: h.vertex_connectivity()
        5

    TESTS:

    Connectivity of some Harary graphs::

        sage: n=10
        sage: for k in range(2,n):
        ...       g = graphs.HararyGraph(k,n)
        ...       if k != g.vertex_connectivity():
        ...          print "Connectivity of Harary graphs not satisfied."
    """
    if k < 2:
        raise ValueError("Connectivity parameter k should be at least 2.")
    if k >= n:
        raise ValueError("Number of vertices n should be greater than k.")

    if k%2 == 0:
        G = CirculantGraph( n, range(1,k//2+1) )
    else:
        if n%2 == 0:
            G = CirculantGraph( n, range(1,(k-1)//2+1) )
            for i in range(n):
                G.add_edge( i, (i + n//2)%n )
        else:
            G = HararyGraph( k-1, n )
            for i in range((n-1)//2 + 1):
                G.add_edge( i, (i + (n-1)//2)%n )
    G.name('Harary graph {0}, {1}'.format(k,n))
    return G

def HyperStarGraph(n,k):
    r"""
    Returns the hyper-star graph HS(n,k).

    The vertices of the hyper-star graph are the set of binary strings
    of length n which contain k 1s. Two vertices, u and v, are adjacent
    only if u can be obtained from v by swapping the first bit with a
    different symbol in another position.

    INPUT:

    -  ``n``

    -  ``k``

    EXAMPLES::

        sage: g = graphs.HyperStarGraph(6,3)
        sage: g.plot() # long time

    REFERENCES:

    - Lee, Hyeong-Ok, Jong-Seok Kim, Eunseuk Oh, and Hyeong-Seok Lim.
      "Hyper-Star Graph: A New Interconnection Network Improving the
      Network Cost of the Hypercube." In Proceedings of the First EurAsian
      Conference on Information and Communication Technology, 858-865.
      Springer-Verlag, 2002.

    AUTHORS:

    - Michael Yurko (2009-09-01)
    """
    from sage.combinat.combination import Combinations
    # dictionary associating the positions of the 1s to the corresponding
    # string: e.g. if n=6 and k=3, comb_to_str([0,1,4])=='110010'
    comb_to_str={}
    for c in Combinations(n,k):
        L = ['0']*n
        for i in c:
            L[i]='1'
        comb_to_str[tuple(c)] = ''.join(L)

    g = Graph(name="HS(%d,%d)"%(n,k))
    g.add_vertices(comb_to_str.values())

    for c in Combinations(range(1,n),k): # 0 is not in c
        L = []
        u = comb_to_str[tuple(c)]
        # switch 0 with the 1s
        for i in xrange(len(c)):
            v = tuple([0]+c[:i]+c[i+1:])
            g.add_edge( u , comb_to_str[v] )

    return g

def LCFGraph(n, shift_list, repeats):
    """
    Returns the cubic graph specified in LCF notation.

    LCF (Lederberg-Coxeter-Fruchte) notation is a concise way of
    describing cubic Hamiltonian graphs. The way a graph is constructed
    is as follows. Since there is a Hamiltonian cycle, we first create
    a cycle on n nodes. The variable shift_list = [s_0, s_1, ...,
    s_k-1] describes edges to be created by the following scheme: for
    each i, connect vertex i to vertex (i + s_i). Then, repeats
    specifies the number of times to repeat this process, where on the
    jth repeat we connect vertex (i + j\*len(shift_list)) to vertex (
    i + j\*len(shift_list) + s_i).

    INPUT:


    -  ``n`` - the number of nodes.

    -  ``shift_list`` - a list of integer shifts mod n.

    -  ``repeats`` - the number of times to repeat the
       process.


    EXAMPLES::

        sage: G = graphs.LCFGraph(4, [2,-2], 2)
        sage: G.is_isomorphic(graphs.TetrahedralGraph())
        True

    ::

        sage: G = graphs.LCFGraph(20, [10,7,4,-4,-7,10,-4,7,-7,4], 2)
        sage: G.is_isomorphic(graphs.DodecahedralGraph())
        True

    ::

        sage: G = graphs.LCFGraph(14, [5,-5], 7)
        sage: G.is_isomorphic(graphs.HeawoodGraph())
        True

    The largest cubic nonplanar graph of diameter three::

        sage: G = graphs.LCFGraph(20, [-10,-7,-5,4,7,-10,-7,-4,5,7,-10,-7,6,-5,7,-10,-7,5,-6,7], 1)
        sage: G.degree()
        [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
        sage: G.diameter()
        3
        sage: G.show()  # long time

    PLOTTING: LCF Graphs are plotted as an n-cycle with edges in the
    middle, as described above.

    REFERENCES:

    - [1] Frucht, R. "A Canonical Representation of Trivalent
      Hamiltonian Graphs." J. Graph Th. 1, 45-60, 1976.

    - [2] Grunbaum, B.  Convex Polytope es. New York: Wiley,
      pp. 362-364, 1967.

    - [3] Lederberg, J. 'DENDRAL-64: A System for Computer
      Construction, Enumeration and Notation of Organic Molecules
      as Tree Structures and Cyclic Graphs. Part II. Topology of
      Cyclic Graphs.' Interim Report to the National Aeronautics
      and Space Administration. Grant NsG 81-60. December 15,
      1965.  http://profiles.nlm.nih.gov/BB/A/B/I/U/_/bbabiu.pdf.
    """
    import networkx
    pos_dict = {}
    for i in range(n):
        x = float(cos(pi/2 + ((2*pi)/n)*i))
        y = float(sin(pi/2 + ((2*pi)/n)*i))
        pos_dict[i] = [x,y]
    return Graph(networkx.LCF_graph(n, shift_list, repeats),\
                 pos=pos_dict, name="LCF Graph")

def MycielskiGraph(k=1, relabel=True):
    r"""
    Returns the `k`-th Mycielski Graph.

    The graph `M_k` is triangle-free and has chromatic number
    equal to `k`. These graphs show, constructively, that there
    are triangle-free graphs with arbitrarily high chromatic
    number.

    The Mycielski graphs are built recursively starting with
    `M_0`, an empty graph; `M_1`, a single vertex graph; and `M_2`
    is the graph `K_2`.  `M_{k+1}` is then built from `M_k`
    as follows:

    If the vertices of `M_k` are `v_1,\ldots,v_n`, then the
    vertices of `M_{k+1}` are
    `v_1,\ldots,v_n,w_1,\ldots,w_n,z`. Vertices `v_1,\ldots,v_n`
    induce a copy of `M_k`. Vertices `w_1,\ldots,w_n` are an
    independent set. Vertex `z` is adjacent to all the
    `w_i`-vertices. Finally, vertex `w_i` is adjacent to vertex
    `v_j` iff `v_i` is adjacent to `v_j`.

    INPUT:

    - ``k`` Number of steps in the construction process.

    - ``relabel`` Relabel the vertices so their names are the integers
      ``range(n)`` where ``n`` is the number of vertices in the graph.

    EXAMPLE:

    The Mycielski graph `M_k` is triangle-free and has chromatic
    number equal to `k`. ::

        sage: g = graphs.MycielskiGraph(5)
        sage: g.is_triangle_free()
        True
        sage: g.chromatic_number()
        5

    The graphs `M_4` is (isomorphic to) the Grotzsch graph. ::

        sage: g = graphs.MycielskiGraph(4)
        sage: g.is_isomorphic(graphs.GrotzschGraph())
        True

    REFERENCES:

    -  [1] Weisstein, Eric W. "Mycielski Graph."
       From MathWorld--A Wolfram Web Resource.
       http://mathworld.wolfram.com/MycielskiGraph.html

    """
    g = Graph()
    g.name("Mycielski Graph " + str(k))

    if k<0:
        raise ValueError("parameter k must be a nonnegative integer")

    if k == 0:
        return g

    if k == 1:
        g.add_vertex(0)
        return g

    if k == 2:
        g.add_edge(0,1)
        return g

    g0 = MycielskiGraph(k-1)
    g = MycielskiStep(g0)
    g.name("Mycielski Graph " + str(k))
    if relabel: g.relabel()

    return g

def MycielskiStep(g):
    r"""
    Perform one iteration of the Mycielski construction.

    See the documentation for ``MycielskiGraph`` which uses this
    method. We expose it to all users in case they may find it
    useful.

    EXAMPLE. One iteration of the Mycielski step applied to the
    5-cycle yields a graph isomorphic to the Grotzsch graph ::

        sage: g = graphs.CycleGraph(5)
        sage: h = graphs.MycielskiStep(g)
        sage: h.is_isomorphic(graphs.GrotzschGraph())
        True
    """

    # Make a copy of the input graph g
    gg = g.copy()

    # rename a vertex v of gg as (1,v)
    renamer = dict( [ (v, (1,v)) for v in g.vertices() ] )
    gg.relabel(renamer)

    # add the w vertices to gg as (2,v)
    wlist = [ (2,v) for v in g.vertices() ]
    gg.add_vertices(wlist)

    # add the z vertex as (0,0)
    gg.add_vertex((0,0))

    # add the edges from z to w_i
    gg.add_edges( [ ( (0,0) , (2,v) ) for v in g.vertices() ] )

    # make the v_i w_j edges
    for v in g.vertices():
        gg.add_edges( [ ((1,v),(2,vv)) for vv in g.neighbors(v) ] )

    return gg

def NKStarGraph(n,k):
    r"""
    Returns the (n,k)-star graph.

    The vertices of the (n,k)-star graph are the set of all arrangements of
    n symbols into labels of length k. There are two adjacency rules for
    the (n,k)-star graph. First, two vertices are adjacent if one can be
    obtained from the other by swapping the first symbol with another
    symbol. Second, two vertices are adjacent if one can be obtained from
    the other by swapping the first symbol with an external symbol (a
    symbol not used in the original label).

    INPUT:

    -  ``n``

    -  ``k``

    EXAMPLES::

        sage: g = graphs.NKStarGraph(4,2)
        sage: g.plot() # long time

    REFERENCES:

    - Wei-Kuo, Chiang, and Chen Rong-Jaye. "The (n, k)-star graph: A
      generalized star graph." Information Processing Letters 56,
      no. 5 (December 8, 1995): 259-264.

    AUTHORS:

    - Michael Yurko (2009-09-01)
    """
    from sage.combinat.permutation import Arrangements
    #set from which to permute
    set = [str(i) for i in xrange(1,n+1)]
    #create dict
    d = {}
    for v in Arrangements(set,k):
        v = list(v) # So we can easily mutate it
        tmp_dict = {}
        #add edges of dimension i
        for i in xrange(1,k):
            #swap 0th and ith element
            v[0], v[i] = v[i], v[0]
            #convert to str and add to list
            vert = "".join(v)
            tmp_dict[vert] = None
            #swap back
            v[0], v[i] = v[i], v[0]
        #add other edges
        tmp_bit = v[0]
        for i in set:
            #check if external
            if not (i in v):
                v[0] = i
                #add edge
                vert = "".join(v)
                tmp_dict[vert] = None
            v[0] = tmp_bit
        d["".join(v)] = tmp_dict
    return Graph(d, name="(%d,%d)-star"%(n,k))

def NStarGraph(n):
    r"""
    Returns the n-star graph.

    The vertices of the n-star graph are the set of permutations on n
    symbols. There is an edge between two vertices if their labels differ
    only in the first and one other position.

    INPUT:

    -  ``n``

    EXAMPLES::

        sage: g = graphs.NStarGraph(4)
        sage: g.plot() # long time

    REFERENCES:

    - S.B. Akers, D. Horel and B. Krishnamurthy, The star graph: An
      attractive alternative to the previous n-cube. In: Proc. Internat.
      Conf. on Parallel Processing (1987), pp. 393--400.

    AUTHORS:

    - Michael Yurko (2009-09-01)
    """
    from sage.combinat.permutation import Permutations
    #set from which to permute
    set = [str(i) for i in xrange(1,n+1)]
    #create dictionary of lists
    #vertices are adjacent if the first element
    #is swapped with the ith element
    d = {}
    for v in Permutations(set):
        v = list(v) # So we can easily mutate it
        tmp_dict = {}
        for i in xrange(1,n):
            if v[0] != v[i]:
                #swap 0th and ith element
                v[0], v[i] = v[i], v[0]
                #convert to str and add to list
                vert = "".join(v)
                tmp_dict[vert] = None
                #swap back
                v[0], v[i] = v[i], v[0]
        d["".join(v)] = tmp_dict
    return Graph(d, name = "%d-star"%n)

def OddGraph(n):
    r"""
    Returns the Odd Graph with parameter `n`.

    The Odd Graph with parameter `n` is defined as the
    Kneser Graph with parameters `2n-1,n-1`.
    Equivalently, the Odd Graph is the graph whose vertices
    are the `n-1`-subsets of `[0,1,\dots,2(n-1)]`, and such
    that two vertices are adjacent if their corresponding sets
    are disjoint.

    For example, the Petersen Graph can be defined
    as the Odd Graph with parameter `3`.

    EXAMPLE::

        sage: OG=graphs.OddGraph(3)
        sage: print OG.vertices()
        [{4, 5}, {1, 3}, {2, 5}, {2, 3}, {3, 4}, {3, 5}, {1, 4}, {1, 5}, {1, 2}, {2, 4}]
        sage: P=graphs.PetersenGraph()
        sage: P.is_isomorphic(OG)
        True

    TESTS::

        sage: KG=graphs.OddGraph(1)
        Traceback (most recent call last):
        ...
        ValueError: Parameter n should be an integer strictly greater than 1
    """

    if not n>1:
        raise ValueError("Parameter n should be an integer strictly greater than 1")
    g = KneserGraph(2*n-1,n-1)
    g.name("Odd Graph with parameter %s" % n)
    return g

def PaleyGraph(q):
    r"""
    Paley graph with `q` vertices

    Parameter `q` must be the power of a prime number and congruent
    to 1 mod 4.

    EXAMPLES::

        sage: G=graphs.PaleyGraph(9);G
        Paley graph with parameter 9: Graph on 9 vertices
        sage: G.is_regular()
        True

    A Paley graph is always self-complementary::

        sage: G.complement().is_isomorphic(G)
        True
    """
    from sage.rings.finite_rings.integer_mod import mod
    from sage.rings.finite_rings.constructor import FiniteField
    assert q.is_prime_power(), "Parameter q must be a prime power"
    assert mod(q,4)==1, "Parameter q must be congruent to 1 mod 4"
    g = Graph([FiniteField(q,'a'), lambda i,j: (i-j).is_square()],
    loops=False, name = "Paley graph with parameter %d"%q)
    return g

def HanoiTowerGraph(pegs, disks, labels=True, positions=True):
    r"""
    Returns the graph whose vertices are the states of the
    Tower of Hanoi puzzle, with edges representing legal moves between states.

    INPUT:

    - ``pegs`` - the number of pegs in the puzzle, 2 or greater
    - ``disks`` - the number of disks in the puzzle, 1 or greater
    - ``labels`` - default: ``True``, if ``True`` the graph contains
      more meaningful labels, see explanation below.  For large instances,
      turn off labels for much faster creation of the graph.
    - ``positions`` - default: ``True``, if ``True`` the graph contains
      layout information.  This creates a planar layout for the case
      of three pegs.  For large instances, turn off layout information
      for much faster creation of the graph.

    OUTPUT:

    The Tower of Hanoi puzzle has a certain number of identical pegs
    and a certain number of disks, each of a different radius.
    Initially the disks are all on a single peg, arranged
    in order of their radii, with the largest on the bottom.

    The goal of the puzzle is to move the disks to any other peg,
    arranged in the same order.  The one constraint is that the
    disks resident on any one peg must always be arranged with larger
    radii lower down.

    The vertices of this graph represent all the possible states
    of this puzzle.  Each state of the puzzle is a tuple with length
    equal to the number of disks, ordered by largest disk first.
    The entry of the tuple is the peg where that disk resides.
    Since disks on a given peg must go down in size as we go
    up the peg, this totally describes the state of the puzzle.

    For example ``(2,0,0)`` means the large disk is on peg 2, the
    medium disk is on peg 0, and the small disk is on peg 0
    (and we know the small disk must be above the medium disk).
    We encode these tuples as integers with a base equal to
    the number of pegs, and low-order digits to the right.

    Two vertices are adjacent if we can change the puzzle from
    one state to the other by moving a single disk.  For example,
    ``(2,0,0)`` is adjacent to ``(2,0,1)`` since we can move
    the small disk off peg 0 and onto (the empty) peg 1.
    So the solution to a 3-disk puzzle (with at least
    two pegs) can be expressed by the shortest path between
    ``(0,0,0)`` and ``(1,1,1)``.  For more on this representation
    of the graph, or its properties, see [ARETT-DOREE]_.

    For greatest speed we create graphs with integer vertices,
    where we encode the tuples as integers with a base equal
    to the number of pegs, and low-order digits to the right.
    So for example, in a 3-peg puzzle with 5 disks, the
    state ``(1,2,0,1,1)`` is encoded as
    `1\ast 3^4 + 2\ast 3^3 + 0\ast 3^2 + 1\ast 3^1 + 1\ast 3^0 = 139`.

    For smaller graphs, the labels that are the tuples are informative,
    but slow down creation of the graph.  Likewise computing layout
    information also incurs a significant speed penalty. For maximum
    speed, turn off labels and layout and decode the
    vertices explicitly as needed.  The
    :meth:`sage.rings.integer.Integer.digits`
    with the ``padsto`` option is a quick way to do this, though you
    may want to reverse the list that is output.

    PLOTTING:

    The layout computed when ``positions = True`` will
    look especially good for the three-peg case, when the graph is known
    to be planar.  Except for two small cases on 4 pegs, the graph is
    otherwise not planar, and likely there is a better way to layout
    the vertices.

    EXAMPLES:

    A classic puzzle uses 3 pegs.  We solve the 5 disk puzzle using
    integer labels and report the minimum number of moves required.
    Note that `3^5-1` is the state where all 5 disks
    are on peg 2. ::

        sage: H = graphs.HanoiTowerGraph(3, 5, labels=False, positions=False)
        sage: H.distance(0, 3^5-1)
        31

    A slightly larger instance. ::

        sage: H = graphs.HanoiTowerGraph(4, 6, labels=False, positions=False)
        sage: H.num_verts()
        4096
        sage: H.distance(0, 4^6-1)
        17

    For a small graph, labels and layout information can be useful.
    Here we explicitly list a solution as a list of states. ::

        sage: H = graphs.HanoiTowerGraph(3, 3, labels=True, positions=True)
        sage: H.shortest_path((0,0,0), (1,1,1))
        [(0, 0, 0), (0, 0, 1), (0, 2, 1), (0, 2, 2), (1, 2, 2), (1, 2, 0), (1, 1, 0), (1, 1, 1)]

    Some facts about this graph with `p` pegs and `d` disks:

    - only automorphisms are the "obvious" ones - renumber the pegs.
    - chromatic number is less than or equal to `p`
    - independence number is `p^{d-1}`

    ::

        sage: H = graphs.HanoiTowerGraph(3,4,labels=False,positions=False)
        sage: H.automorphism_group().is_isomorphic(SymmetricGroup(3))
        True
        sage: H.chromatic_number()
        3
        sage: len(H.independent_set()) == 3^(4-1)
        True

    TESTS:

    It is an error to have just one peg (or less). ::

        sage: graphs.HanoiTowerGraph(1, 5)
        Traceback (most recent call last):
        ...
        ValueError: Pegs for Tower of Hanoi graph should be two or greater (not 1)

    It is an error to have zero disks (or less). ::

        sage: graphs.HanoiTowerGraph(2, 0)
        Traceback (most recent call last):
        ...
        ValueError: Disks for Tower of Hanoi graph should be one or greater (not 0)

    .. rubric:: Citations

    .. [ARETT-DOREE] Arett, Danielle and Doree, Suzanne
       "Coloring and counting on the Hanoi graphs"
       Mathematics Magazine, Volume 83, Number 3, June 2010, pages 200-9


    AUTHOR:

    - Rob Beezer, (2009-12-26), with assistance from Su Doree

    """

    # sanitize input
    from sage.rings.all import Integer
    pegs = Integer(pegs)
    if pegs < 2:
        raise ValueError("Pegs for Tower of Hanoi graph should be two or greater (not %d)" % pegs)
    disks = Integer(disks)
    if disks < 1:
        raise ValueError("Disks for Tower of Hanoi graph should be one or greater (not %d)" % disks)

    # Each state of the puzzle is a tuple with length
    # equal to the number of disks, ordered by largest disk first
    # The entry of the tuple is the peg where that disk resides
    # Since disks on a given peg must go down in size as we go
    # up the peg, this totally describes the puzzle
    # We encode these tuples as integers with a base equal to
    # the number of pegs, and low-order digits to the right

    # complete graph on number of pegs when just a single disk
    edges = [[i,j] for i in range(pegs) for j in range(i+1,pegs)]

    nverts = 1
    for d in range(2, disks+1):
        prevedges = edges      # remember subgraph to build from
        nverts = pegs*nverts   # pegs^(d-1)
        edges = []

        # Take an edge, change its two states in the same way by adding
        # a large disk to the bottom of the same peg in each state
        # This is accomplished by adding a multiple of pegs^(d-1)
        for p in range(pegs):
            largedisk = p*nverts
            for anedge in prevedges:
                edges.append([anedge[0]+largedisk, anedge[1]+largedisk])

        # Two new states may only differ in the large disk
        # being the only disk on two different pegs, thus
        # otherwise being a common state with one less disk
        # We construct all such pairs of new states and add as edges
        from sage.combinat.subset import Subsets
        for state in range(nverts):
            emptypegs = range(pegs)
            reduced_state = state
            for i in range(d-1):
                apeg = reduced_state % pegs
                if apeg in emptypegs:
                    emptypegs.remove(apeg)
                reduced_state = reduced_state//pegs
            for freea, freeb in Subsets(emptypegs, 2):
                edges.append([freea*nverts+state,freeb*nverts+state])

    H = Graph({}, loops=False, multiedges=False)
    H.add_edges(edges)


    # Making labels and/or computing positions can take a long time,
    # relative to just constructing the edges on integer vertices.
    # We try to minimize coercion overhead, but need Sage
    # Integers in order to use digits() for labels.
    # Getting the digits with custom code was no faster.
    # Layouts are circular (symmetric on the number of pegs)
    # radiating outward to the number of disks (radius)
    # Algorithm uses some combination of alternate
    # clockwise/counterclockwise placements, which
    # works well for three pegs (planar layout)
    #
    from sage.functions.trig import sin, cos, csc
    if labels or positions:
        mapping = {}
        pos = {}
        a = Integer(-1)
        one = Integer(1)
        if positions:
            radius_multiplier = 1 + csc(pi/pegs)
            sine = []; cosine = []
            for i in range(pegs):
                angle = 2*i*pi/float(pegs)
                sine.append(sin(angle))
                cosine.append(cos(angle))
        for i in range(pegs**disks):
            a += one
            state = a.digits(base=pegs, padto=disks)
            if labels:
                state.reverse()
                mapping[i] = tuple(state)
                state.reverse()
            if positions:
                locx = 0.0; locy = 0.0
                radius = 1.0
                parity = -1.0
                for index in range(disks):
                    p = state[index]
                    radius *= radius_multiplier
                    parity *= -1.0
                    locx_temp = cosine[p]*locx - parity*sine[p]*locy + radius*cosine[p]
                    locy_temp = parity*sine[p]*locx + cosine[p]*locy - radius*parity*sine[p]
                    locx = locx_temp
                    locy = locy_temp
                pos[i] = (locx,locy)
        # set positions, then relabel (not vice versa)
        if positions:
            H.set_pos(pos)
        if labels:
            H.relabel(mapping)

    return H

def line_graph_forbidden_subgraphs():
    r"""
    Returns the 9 forbidden subgraphs of a line graph.

    `Wikipedia article on the line graphs
    <http://en.wikipedia.org/wiki/Line_graph>`_

    The graphs are returned in the ordering given by the Wikipedia
    drawing, read from left to right and from top to bottom.

    EXAMPLE::

        sage: graphs.line_graph_forbidden_subgraphs()
        [Claw graph: Graph on 4 vertices,
        Graph on 6 vertices,
        Graph on 6 vertices,
        Graph on 5 vertices,
        Graph on 6 vertices,
        Graph on 6 vertices,
        Graph on 6 vertices,
        Graph on 6 vertices,
        Graph on 5 vertices]

    """
    from sage.graphs.all import Graph
    from sage.graphs.generators.basic import ClawGraph
    graphs = [ClawGraph()]

    graphs.append(Graph({
                0: [1, 2, 3],
                1: [2, 3],
                4: [2],
                5: [3]
                }))

    graphs.append(Graph({
                0: [1, 2, 3, 4],
                1: [2, 3, 4],
                3: [4],
                2: [5]
                }))

    graphs.append(Graph({
                0: [1, 2, 3],
                1: [2, 3],
                4: [2, 3]
                }))

    graphs.append(Graph({
                0: [1, 2, 3],
                1: [2, 3],
                4: [2],
                5: [3, 4]
                }))

    graphs.append(Graph({
                0: [1, 2, 3, 4],
                1: [2, 3, 4],
                3: [4],
                5: [2, 0, 1]
                }))

    graphs.append(Graph({
                5: [0, 1, 2, 3, 4],
                0: [1, 4],
                2: [1, 3],
                3: [4]
                }))

    graphs.append(Graph({
                1: [0, 2, 3, 4],
                3: [0, 4],
                2: [4, 5],
                4: [5]
                }))

    graphs.append(Graph({
                0: [1, 2, 3],
                1: [2, 3, 4],
                2: [3, 4],
                3: [4]
                }))

    return graphs


def petersen_family(generate=False):
    r"""
    Returns the Petersen family

    The Petersen family is a collection of 7 graphs which are the forbidden
    minors of the linklessly embeddable graphs. For more information see the
    :wikipedia:`Petersen_family`.

    INPUT:

    - ``generate`` (boolean) -- whether to generate the family from the
      `\Delta-Y` transformations. When set to ``False`` (default) a hardcoded
      version of the graphs (with a prettier layout) is returned.

    EXAMPLE::

        sage: graphs.petersen_family()
        [Petersen graph: Graph on 10 vertices,
         Complete graph: Graph on 6 vertices,
         Multipartite Graph with set sizes [3, 3, 1]: Graph on 7 vertices,
         Graph on 8 vertices,
         Graph on 9 vertices,
         Graph on 7 vertices,
         Graph on 8 vertices]

    The two different inputs generate the same graphs::

        sage: F1 = graphs.petersen_family(generate=False)
        sage: F2 = graphs.petersen_family(generate=True)
        sage: F1 = [g.canonical_label().graph6_string() for g in F1]
        sage: F2 = [g.canonical_label().graph6_string() for g in F2]
        sage: set(F1) == set(F2)
        True
    """
    from sage.graphs.generators.smallgraphs import PetersenGraph
    if not generate:
        from sage.graphs.generators.basic import CompleteGraph, \
             CompleteBipartiteGraph, CompleteMultipartiteGraph
        from sage.graphs.graph_plot import _circle_embedding
        l = [PetersenGraph(), CompleteGraph(6),
             CompleteMultipartiteGraph([3, 3, 1])]
        g = CompleteBipartiteGraph(4, 4)
        g.delete_edge(0, 4)
        g.name("")
        l.append(g)
        g = Graph('HKN?Yeb')
        _circle_embedding(g, [1, 2, 4, 3, 0, 5])
        _circle_embedding(g, [6, 7, 8], radius=.6, shift=1.25)
        l.append(g)
        g = Graph('Fs\\zw')
        _circle_embedding(g, [1, 2, 3])
        _circle_embedding(g, [4, 5, 6], radius=.7)
        g.get_pos()[0] = (0, 0)
        l.append(g)
        g = Graph('GYQ[p{')
        _circle_embedding(g, [1, 4, 6, 0, 5, 7, 3], shift=0.25)
        g.get_pos()[2] = (0, 0)
        l.append(g)
        return l

    def DeltaYTrans(G, triangle):
        """
        Apply a Delta-Y transformation to a given triangle of G.
        """
        a, b, c = triangle
        G = G.copy()
        G.delete_edges([(a, b), (b, c), (c, a)])
        v = G.order()
        G.add_edges([(a, v), (b, v), (c, v)])
        return G.canonical_label()

    def YDeltaTrans(G, v):
        """
        Apply a Y-Delta transformation to a given vertex v of G.
        """
        G = G.copy()
        a, b, c = G.neighbors(v)
        G.delete_vertex(v)
        G.add_cycle([a, b, c])
        return G.canonical_label()

    # We start from the Petersen Graph, and apply Y-Delta transform
    # for as long as we generate new graphs.
    P = PetersenGraph()

    l = set([])
    l_new = [P.canonical_label().graph6_string()]

    while l_new:
        g = l_new.pop(0)
        if g in l:
            continue
        l.add(g)
        g = Graph(g)
        # All possible Delta-Y transforms
        for t in g.subgraph_search_iterator(Graph({1: [2, 3], 2: [3]})):
            l_new.append(DeltaYTrans(g, t).graph6_string())
        # All possible Y-Delta transforms
        for v in g:
            if g.degree(v) == 3:
                l_new.append(YDeltaTrans(g, v).graph6_string())

    return [Graph(x) for x in l]


def WheelGraph(n):
    """
    Returns a Wheel graph with n nodes.

    A Wheel graph is a basic structure where one node is connected to
    all other nodes and those (outer) nodes are connected cyclically.

    This constructor depends on NetworkX numeric labels.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, each wheel
    graph will be displayed with the first (0) node in the center, the
    second node at the top, and the rest following in a
    counterclockwise manner.

    With the wheel graph, we see that it doesn't take a very large n at
    all for the spring-layout to give a counter-intuitive display. (See
    Graphics Array examples below).

    EXAMPLES: We view many wheel graphs with a Sage Graphics Array,
    first with this constructor (i.e., the position dictionary
    filled)::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ...    k = graphs.WheelGraph(i+3)
        ...    g.append(k)
        ...
        sage: for i in range(3):
        ...    n = []
        ...    for m in range(3):
        ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ...    j.append(n)
        ...
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time

    Next, using the spring-layout algorithm::

        sage: import networkx
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ...    spr = networkx.wheel_graph(i+3)
        ...    k = Graph(spr)
        ...    g.append(k)
        ...
        sage: for i in range(3):
        ...    n = []
        ...    for m in range(3):
        ...        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ...    j.append(n)
        ...
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time

    Compare the plotting::

        sage: n = networkx.wheel_graph(23)
        sage: spring23 = Graph(n)
        sage: posdict23 = graphs.WheelGraph(23)
        sage: spring23.show() # long time
        sage: posdict23.show() # long time
    """
    pos_dict = {}
    pos_dict[0] = (0,0)
    for i in range(1,n):
        x = float(cos((pi/2) + ((2*pi)/(n-1))*(i-1)))
        y = float(sin((pi/2) + ((2*pi)/(n-1))*(i-1)))
        pos_dict[i] = (x,y)
    import networkx
    G = networkx.wheel_graph(n)
    return Graph(G, pos=pos_dict, name="Wheel graph")

def trees(vertices):
    r"""
    Returns a generator of the distinct trees on a fixed number of vertices.

    INPUT:

    -  ``vertices`` - the size of the trees created.

    OUTPUT:

    A generator which creates an exhaustive, duplicate-free listing
    of the connected free (unlabeled) trees with ``vertices`` number
    of vertices.  A tree is a graph with no cycles.

    ALGORITHM:

    Uses an algorithm that generates each new tree
    in constant time.  See the documentation for, and implementation
    of, the :mod:`sage.graphs.trees` module, including a citation.

    EXAMPLES:

    We create an iterator, then loop over its elements. ::

        sage: tree_iterator = graphs.trees(7)
        sage: for T in tree_iterator:
        ...     print T.degree_sequence()
        [2, 2, 2, 2, 2, 1, 1]
        [3, 2, 2, 2, 1, 1, 1]
        [3, 2, 2, 2, 1, 1, 1]
        [4, 2, 2, 1, 1, 1, 1]
        [3, 3, 2, 1, 1, 1, 1]
        [3, 3, 2, 1, 1, 1, 1]
        [4, 3, 1, 1, 1, 1, 1]
        [3, 2, 2, 2, 1, 1, 1]
        [4, 2, 2, 1, 1, 1, 1]
        [5, 2, 1, 1, 1, 1, 1]
        [6, 1, 1, 1, 1, 1, 1]

    The number of trees on the first few vertex counts.
    This is sequence A000055 in Sloane's OEIS. ::

        sage: [len(list(graphs.trees(i))) for i in range(0, 15)]
        [1, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159]
    """
    from sage.graphs.trees import TreeIterator
    return iter(TreeIterator(vertices))

def RingedTree(k, vertex_labels = True):
    r"""
    Return the ringed tree on k-levels.

    A ringed tree of level `k` is a binary tree with `k` levels (counting
    the root as a level), in which all vertices at the same level are connected
    by a ring.

    More precisely, in each layer of the binary tree (i.e. a layer is the set of
    vertices `[2^i...2^{i+1}-1]`) two vertices `u,v` are adjacent if `u=v+1` or
    if `u=2^i` and `v=`2^{i+1}-1`.

    Ringed trees are defined in [CFHM12]_.

    INPUT:

    - ``k`` -- the number of levels of the ringed tree.

    - ``vertex_labels`` (boolean) -- whether to label vertices as binary words
      (default) or as integers.

    EXAMPLE::

        sage: G = graphs.RingedTree(5)
        sage: P = G.plot(vertex_labels=False, vertex_size=10)
        sage: P.show() # long time
        sage: G.vertices()
        ['', '0', '00', '000', '0000', '0001', '001', '0010', '0011', '01',
         '010', '0100', '0101', '011', '0110', '0111', '1', '10', '100',
         '1000', '1001', '101', '1010', '1011', '11', '110', '1100', '1101',
         '111', '1110', '1111']

    TEST::

        sage: G = graphs.RingedTree(-1)
        Traceback (most recent call last):
        ...
        ValueError: The number of levels must be >= 1.
        sage: G = graphs.RingedTree(5, vertex_labels = False)
        sage: G.vertices()
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]

    REFERENCES:

    .. [CFHM12] On the Hyperbolicity of Small-World and Tree-Like Random Graphs
      Wei Chen, Wenjie Fang, Guangda Hu, Michael W. Mahoney
      http://arxiv.org/abs/1201.1717
    """
    if k<1:
        raise ValueError('The number of levels must be >= 1.')

    from sage.graphs.graph_plot import _circle_embedding

    # Creating the Balanced tree, which contains most edges already
    g = BalancedTree(2,k-1)
    g.name('Ringed Tree on '+str(k)+' levels')

    # We consider edges layer by layer
    for i in range(1,k):
        vertices = range(2**(i)-1,2**(i+1)-1)

        # Add the missing edges
        g.add_cycle(vertices)

        # And set the vertices' positions
        radius = i if i <= 1 else 1.5**i
        shift = -2**(i-2)+.5 if i > 1 else 0
        _circle_embedding(g, vertices, radius = radius, shift = shift)

    # Specific position for the central vertex
    g.get_pos()[0] = (0,0.2)

    # Relabel vertices as binary words
    if not vertex_labels:
        return g

    vertices = ['']
    for i in range(k-1):
        for j in range(2**(i)-1,2**(i+1)-1):
            v = vertices[j]
            vertices.append(v+'0')
            vertices.append(v+'1')

    g.relabel(vertices)

    return g

def SymplecticGraph(d,q):
    r"""
    Returns the Symplectic graph `Sp(d,q)`

    The Symplectic Graph `Sp(d,q)` is built from a projective space of dimension
    `d-1` over a field `F_q`, and a symplectic form `f`. Two vertices `u,v` are
    made adjacent if `f(u,v)=0`.

    See the `page on symplectic graphs on Andries Brouwer's website
    <http://www.win.tue.nl/~aeb/graphs/Sp.html>`_.

    INPUT:

    - ``d,q`` (integers) -- note that only even values of `d` are accepted by
      the function.

    EXAMPLES::

        sage: g = graphs.SymplecticGraph(6,2)
        sage: g.is_strongly_regular(parameters=True)
        (63, 30, 13, 15)
        sage: set(g.spectrum()) == {-5, 3, 30}
        True
    """
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.modules.free_module import VectorSpace
    from sage.schemes.projective.projective_space import ProjectiveSpace
    from sage.matrix.constructor import identity_matrix, block_matrix, zero_matrix

    if d < 1 or d%2 != 0:
        raise ValueError("d must be even and greater than 2")

    F = FiniteField(q,"x")
    M = block_matrix(F, 2, 2,
                     [zero_matrix(F,d/2),
                      identity_matrix(F,d/2),
                      -identity_matrix(F,d/2),
                      zero_matrix(F,d/2)])

    V = VectorSpace(F,d)
    PV = list(ProjectiveSpace(d-1,F))
    G = Graph([map(tuple,PV), lambda x,y:V(x)*(M*V(y)) == 0], loops = False)
    G.name("Symplectic Graph Sp("+str(d)+","+str(q)+")")
    G.relabel()
    return G

def AffineOrthogonalPolarGraph(d,q,sign="+"):
    r"""
    Returns the affine polar graph `VO^+(d,q),VO^-(d,q)` or `VO(d,q)`.

    Affine Polar graphs are built from a `d`-dimensional vector space over
    `F_q`, and a quadratic form which is hyperbolic, elliptic or parabolic
    according to the value of ``sign``.

    Note that `VO^+(d,q),VO^-(d,q)` are strongly regular graphs, while `VO(d,q)`
    is not.

    For more information on Affine Polar graphs, see `Affine Polar
    Graphs page of Andries Brouwer's website
    <http://www.win.tue.nl/~aeb/graphs/VO.html>`_.

    INPUT:

    - ``d`` (integer) -- ``d`` must be even if ``sign != None``, and odd
      otherwise.

    - ``q`` (integer) -- a power of a prime number, as `F_q` must exist.

    - ``sign`` -- must be qual to ``"+"``, ``"-"``, or ``None`` to compute
      (respectively) `VO^+(d,q),VO^-(d,q)` or `VO(d,q)`. By default
      ``sign="+"``.

    .. NOTE::

        The graph `VO^\epsilon(d,q)` is the graph induced by the
        non-neighbors of a vertex in an :meth:`Orthogonal Polar Graph
        <OrthogonalPolarGraph>` `O^\epsilon(d+2,q)`.

    EXAMPLES:

    The :meth:`Brouwer-Haemers graph <BrouwerHaemersGraph>` is isomorphic to
    `VO^-(4,3)`::

        sage: g = graphs.AffineOrthogonalPolarGraph(4,3,"-")
        sage: g.is_isomorphic(graphs.BrouwerHaemersGraph())
        True

    Some examples from `Brouwer's table or strongly regular graphs
    <http://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`_::

        sage: g = graphs.AffineOrthogonalPolarGraph(6,2,"-"); g
        Affine Polar Graph VO^-(6,2): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 27, 10, 12)
        sage: g = graphs.AffineOrthogonalPolarGraph(6,2,"+"); g
        Affine Polar Graph VO^+(6,2): Graph on 64 vertices
        sage: g.is_strongly_regular(parameters=True)
        (64, 35, 18, 20)

    When ``sign is None``::

        sage: g = graphs.AffineOrthogonalPolarGraph(5,2,None); g
        Affine Polar Graph VO^-(5,2): Graph on 32 vertices
        sage: g.is_strongly_regular(parameters=True)
        False
        sage: g.is_regular()
        True
        sage: g.is_vertex_transitive()
        True
    """
    if sign in ["+","-"]:
        s = 1 if sign == "+" else -1
        if d%2 == 1:
            raise ValueError("d must be even when sign!=None")
    else:
        if d%2 == 0:
            raise ValueError("d must be odd when sign==None")
        s = 0

    from sage.interfaces.gap import gap
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.modules.free_module import VectorSpace
    from sage.matrix.constructor import Matrix
    from sage.libs.gap.libgap import libgap
    from itertools import combinations

    M = Matrix(libgap.InvariantQuadraticForm(libgap.GeneralOrthogonalGroup(s,d,q))['matrix'])
    F = libgap.GF(q).sage()
    V = list(VectorSpace(F,d))

    G = Graph()
    G.add_vertices(map(tuple,V))
    for x,y in combinations(V,2):
        if not (x-y)*M*(x-y):
            G.add_edge(tuple(x),tuple(y))

    G.name("Affine Polar Graph VO^"+str('+' if s == 1 else '-')+"("+str(d)+","+str(q)+")")
    G.relabel()
    return G

def OrthogonalPolarGraph(m, q, sign="+"):
    r"""
    Returns the Orthogonal Polar Graph `O^{\epsilon}(m,q)`.

    For more information on Orthogonal Polar graphs, see see the `page of
    Andries Brouwer's website <http://www.win.tue.nl/~aeb/graphs/srghub.html>`_.

    INPUT:

    - ``m,q`` (integers) -- `q` must be a prime power.

    - ``sign`` -- ``"+"`` or ``"-"`` if `m` is even, ``"+"`` (default)
      otherwise.

    EXAMPLES::

        sage: G = graphs.OrthogonalPolarGraph(6,3,"+"); G
        Orthogonal Polar Graph O^+(6, 3): Graph on 130 vertices
        sage: G.is_strongly_regular(parameters=True)
        (130, 48, 20, 16)
        sage: G = graphs.OrthogonalPolarGraph(6,3,"-"); G
        Orthogonal Polar Graph O^-(6, 3): Graph on 112 vertices
        sage: G.is_strongly_regular(parameters=True)
        (112, 30, 2, 10)
        sage: G = graphs.OrthogonalPolarGraph(5,3); G
        Orthogonal Polar Graph O(5, 3): Graph on 40 vertices
        sage: G.is_strongly_regular(parameters=True)
        (40, 12, 2, 4)

    TESTS::

        sage: G = graphs.OrthogonalPolarGraph(4,3,"")
        Traceback (most recent call last):
        ...
        ValueError: sign must be equal to either '-' or '+' when m is even
        sage: G = graphs.OrthogonalPolarGraph(5,3,"-")
        Traceback (most recent call last):
        ...
        ValueError: sign must be equal to either '' or '+' when m is odd
    """
    from sage.schemes.projective.projective_space import ProjectiveSpace
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.modules.free_module_element import free_module_element as vector
    from sage.matrix.constructor import Matrix
    from sage.libs.gap.libgap import libgap
    from itertools import combinations

    if m % 2 == 0:
        if sign != "+" and sign != "-":
            raise ValueError("sign must be equal to either '-' or '+' when "
                             "m is even")
    else:
        if sign != "" and sign != "+":
            raise ValueError("sign must be equal to either '' or '+' when "
                             "m is odd")
        sign = ""

    e = {'+': 1,
         '-': -1,
         '' : 0}[sign]

    M = Matrix(libgap.InvariantQuadraticForm(libgap.GeneralOrthogonalGroup(e,m,q))['matrix'])
    Fq = libgap.GF(q).sage()
    PG = ProjectiveSpace(m - 1, Fq)
    m_over_two = m // 2

    def F(x):
        return x*M*x

    V = [x for x in PG if F(vector(x)) == 0]

    G = Graph([V,lambda x,y:F(vector(x)-vector(y))==0],loops=False)

    G.relabel()
    G.name("Orthogonal Polar Graph O" + ("^" + sign if sign else "") + str((m, q)))
    return G
