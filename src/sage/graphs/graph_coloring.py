"""
Graph Coloring

AUTHORS:

- Tom Boothby (2008-02-21): Initial version
- Carlo Hamalainen (2009-03-28): minor change: switch to C++ DLX solver
- Nathann Cohen (2009-10-24): Coloring methods using linear programming
"""

#*****************************************************************************
#           Copyright (C) 2008 Tom Boothby <boothby@u.washington.edu>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.matrices.dlxcpp import DLXCPP
from sage.all import Matrix, vector, QQ
from sage.plot.colors import rainbow
from chrompoly import chromatic_polynomial
from graph_generators import GraphGenerators

def all_graph_colorings(G,n,count_only=False):
    r"""
    Computes all `n`-colorings of the graph `G` by casting the graph
    coloring problem into an exact cover problem, and passing this
    into an implementation of the Dancing Links algorithm described
    by Knuth (who attributes the idea to Hitotumatu and Noshita).

    The construction works as follows. Columns:

    * The first `|V|` columns correspond to a vertex -- a `1` in this
      column indicates that that vertex has a color.

    * After those `|V|` columns, we add `n*|E|` columns -- a `1` in
      these columns indicate that a particular edge is
      incident to a vertex with a certain color.

    Rows:

    * For each vertex, add `n` rows; one for each color `c`.  Place
      a `1` in the column corresponding to the vertex, and a `1`
      in the appropriate column for each edge incident to the
      vertex, indicating that that edge is incident to the
      color `c`.

    * If `n > 2`, the above construction cannot be exactly covered
      since each edge will be incident to only two vertices
      (and hence two colors) - so we add `n*|E|` rows, each one
      containing a `1` for each of the `n*|E|` columns.  These
      get added to the cover solutions "for free" during the
      backtracking.

    Note that this construction results in `n*|V| + 2*n*|E| + n*|E|`
    entries in the matrix.  The Dancing Links algorithm uses a
    sparse representation, so if the graph is simple, `|E| \leq |V|^2`
    and `n <= |V|`, this construction runs in `O(|V|^3)` time.
    Back-conversion to a coloring solution is a simple scan of the
    solutions, which will contain `|V| + (n-2)*|E|` entries,  so
    runs in `O(|V|^3)` time also.  For most graphs, the conversion
    will be much faster -- for example, a planar graph will be
    transformed for `4`-coloring in linear time since `|E| = O(|V|)`.

    REFERENCES:

    http://www-cs-staff.stanford.edu/~uno/papers/dancing-color.ps.gz

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import all_graph_colorings
        sage: G = Graph({0:[1,2,3],1:[2]})
        sage: n = 0
        sage: for C in all_graph_colorings(G,3):
        ...       parts = [C[k] for k in C]
        ...       for P in parts:
        ...           l = len(P)
        ...           for i in range(l):
        ...               for j in range(i+1,l):
        ...                   if G.has_edge(P[i],P[j]):
        ...                       raise RuntimeError, "Coloring Failed."
        ...       n+=1
        sage: print "G has %s 3-colorings."%n
        G has 12 3-colorings.

    TESTS::

        sage: G = Graph({0:[1,2,3],1:[2]})
        sage: for C in all_graph_colorings(G,0): print C
        sage: for C in all_graph_colorings(G,-1): print C
        Traceback (most recent call last):
        ...
        ValueError: n must be non-negative.
    """

    if n == 0: return
    if n < 0: raise ValueError, "n must be non-negative."

    V = G.vertices()
    E = G.edges()

    nV=len(V)
    nE=len(E)

    ones = []
    N = xrange(n)
    Vd= {}
    colormap = {}
    k = 0
    for i in range(nV):
        v = V[i]
        Vd[v] = i
        for c in N:
            ones.append([k, [i]])
            colormap[k] = (v,c)
            k+=1

    kk = nV
    for e in E:
        for c in N:
            v0 = n*Vd[e[0]]+c
            v1 = n*Vd[e[1]]+c
            ones[v0][1].append(kk+c)
            ones[v1][1].append(kk+c)
        kk+=n

    if n > 2:
        for i in range(n*nE):
            ones.append([k+i, [nV+i]])

    colors = rainbow(n)

    for i in range(len(ones)): ones[i] = ones[i][1]

    try:
        for a in DLXCPP(ones):
            if count_only:
                yield 1
                continue
            coloring = {}
            for x in a:
                if colormap.has_key(x):
                    v,c = colormap[x]
                    if coloring.has_key(colors[c]):
                        coloring[colors[c]].append(v)
                    else:
                        coloring[colors[c]] = [v]
            yield coloring
    except RuntimeError:
        raise RuntimeError, "Too much recursion!  Graph coloring failed."

def first_coloring(G, n=0, hex_colors=False):
    r"""
    Given a graph, and optionally a natural number `n`, returns
    the first coloring we find with at least `n` colors.

    INPUT:

    - ``hex_colors`` -- (default: ``False``) when set to ``True``, the
      partition returned is a dictionary whose keys are colors and whose
      values are the color classes (ideal for plotting).

    -  ``n`` -- The minimal number of colors to try.

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import first_coloring
        sage: G = Graph({0: [1, 2, 3], 1: [2]})
        sage: first_coloring(G, 3)
        [[1, 3], [0], [2]]
    """
    o = G.order()
    for m in xrange(n, o + 1):
        for C in all_graph_colorings(G, m):
            if hex_colors:
                return C
            else:
                return C.values()

def number_of_n_colorings(G,n):
    r"""
    Given a graph `G` and a natural number `n`, returns the number of
    `n`-colorings of the graph.

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import number_of_n_colorings
        sage: G = Graph({0:[1,2,3],1:[2]})
        sage: number_of_n_colorings(G,3)
        12
    """
    #Take care of the stupid stuff
    if n == 1:
        return int(len(G.edges()) == 0)
    if n < 1:
        if n == 0:
            return int(len(G.vertices()) == 0)
        else:
            #negative colors?? what does that even mean?
            return 0

    m = 0
    for C in all_graph_colorings(G,n,count_only=True):
        m+=1
    return m

def numbers_of_colorings(G):
    r"""
    Returns the number of `n`-colorings of the graph `G` for `n` from
    `0` to `|V|`.

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import numbers_of_colorings
        sage: G = Graph({0:[1,2,3],1:[2]})
        sage: numbers_of_colorings(G)
        [0, 0, 0, 12, 72]
    """
    o = G.order()
    return [number_of_n_colorings(G,i) for i in range(0,o+1)]

def chromatic_number(G):
    r"""
    Returns the minimal number of colors needed to color the
    vertices of the graph `G`.

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import chromatic_number
        sage: G = Graph({0:[1,2,3],1:[2]})
        sage: chromatic_number(G)
        3

        sage: G = graphs.PetersenGraph()
        sage: G.chromatic_number()
        3
    """
    o = G.order()
    if o == 0:
        return 0
    if len(G.edges()) == 0:
        return 1
    elif G.is_bipartite(): #can we do it in linear time?
        return 2
    else: #counting cliques is faster than our brute-force method...
        m = G.clique_number()
    if m >= o-1: #marginal improvement... if there's an o-1 clique and not an o clique, don't waste our time coloring.
        return m
    for n in range(m,o+1):
        for C in all_graph_colorings(G,n):
            return n

def vertex_coloring(g, k=None, value_only=False, hex_colors=False, log=0):
    r"""
    Computes the chromatic number of the given graph or tests its
    `k`-colorability. See http://en.wikipedia.org/wiki/Graph_coloring for
    further details on graph coloring.

    INPUT:

    - ``g`` -- a graph.

    - ``k`` -- (default: ``None``) tests whether the graph is `k`-colorable.
      The function returns a partition of the vertex set in `k` independent
      sets if possible and ``False`` otherwise.

    - ``value_only`` -- (default: ``False``):

      - When set to ``True``, only the chromatic number is returned.

      - When set to ``False`` (default), a partition of the vertex set into
        independent sets is returned if possible.

    - ``hex_colors`` -- (default: ``False``) when set to ``True``, the
      partition returned is a dictionary whose keys are colors and whose
      values are the color classes (ideal for plotting).

    - ``log`` -- (default: ``0``) as vertex-coloring is an `NP`-complete
      problem, this function may take some time depending on the graph.
      Use ``log`` to define the level of verbosity you want from the
      linear program solver. By default ``log=0``, meaning that there will
      be no message printed by the solver.

    OUTPUT:

    - If ``k=None`` and ``value_only=False``, then return a partition of the
      vertex set into the minimum possible of independent sets.

    - If ``k=None`` and ``value_only=True``, return the chromatic number.

    - If ``k`` is set and ``value_only=None``, return ``False`` if the
      graph is not `k`-colorable, and a partition of the vertex set into
      `k` independent sets otherwise.

    - If ``k`` is set and ``value_only=True``, test whether the graph is
      `k`-colorable, and return ``True`` or ``False`` accordingly.

    EXAMPLE::

       sage: from sage.graphs.graph_coloring import vertex_coloring
       sage: g = graphs.PetersenGraph()
       sage: vertex_coloring(g, value_only=True) # optional - requires GLPK or CBC
       3
    """
    from sage.numerical.mip import MixedIntegerLinearProgram
    from sage.plot.colors import rainbow

    # If k==None, tries to find an optimal coloring
    if k is None:
        # No need to start a linear program if the graph is an
        # independent set or bipartite.
        # - Independent set
        if g.size() == 0:
            if value_only:
                return 1
            elif hex_colors:
                return {rainbow(1)[0]: g.vertices()}
            else:
                return g.vertices()
        # - Bipartite set
        if g.is_bipartite():
            if value_only:
                return 2
            elif hex_colors:
                return dict(zip(rainbow(2), g.bipartite_sets()))
            else:
                return g.bipartite_sets()
        # - No need to try any k smaller than the maximum clique in the graph
        # - max, because the graph could be triangle-free.
        k = max(3, g.clique_number())
        while True:
            # tries to color the graph, increasing k each time it fails.
            tmp = vertex_coloring(g, k=k, value_only=value_only,
                                  hex_colors=hex_colors, log=log)
            if tmp is not False:
                if value_only:
                    return k
                else:
                    return tmp
            k += 1
    else:
        # Is the graph empty?
        # If the graph is empty, something should be returned.
        # This is not so stupid, as the graph could be emptied
        # by the test of degeneracy.
        if g.order() == 0:
            if value_only:
                return True
            elif hex_colors:
                return dict([(color, []) for color in rainbow(k)])
            else:
                return [[] for i in xrange(k)]
        # Is the graph connected?
        # This is not so stupid, as the graph could be disconnected
        # by the test of degeneracy (as previously).
        if not g.is_connected():
            if value_only:
                for component in g.connected_components():
                    tmp = vertex_coloring(g.subgraph(component), k=k,
                                          value_only=value_only,
                                          hex_colors=hex_colors,
                                          log=log)
                    if tmp is False:
                        return False
                return True
            colorings = []
            for component in g.connected_components():
                tmp = vertex_coloring(g.subgraph(component), k=k,
                                      value_only=value_only,
                                      hex_colors=False, log=log)
                if tmp is False:
                    return False
                colorings.append(tmp)
            value = [[] for color in xrange(k)]
            for color in xrange(k):
                for component in colorings:
                    value[color].extend(component[color])
            if hex_colors:
                return dict(zip(rainbow(k), value))
            else:
                return value

        # Degeneracy
        # Vertices whose degree is less than k are of no importance in
        # the coloring.
        if min(g.degree()) < k:
            vertices = set(g.vertices())
            deg = []
            tmp = [v for v in vertices if g.degree(v) < k]
            while len(tmp) > 0:
                v = tmp.pop(0)
                neighbors = list(set(g.neighbors(v)) & vertices)
                if v in vertices and len(neighbors) < k:
                    vertices.remove(v)
                    tmp.extend(neighbors)
                    deg.append(v)
            if value_only:
                return vertex_coloring(g.subgraph(list(vertices)), k=k,
                                       value_only=value_only,
                                       hex_colors=hex_colors,
                                       log=log)
            value = vertex_coloring(g.subgraph(list(vertices)), k=k,
                                    value_only=value_only,
                                    hex_colors=False,
                                    log=log)
            if value is False:
                return False
            while len(deg) > 0:
                for classe in value:
                    if len(list(set(classe) & set(g.neighbors(deg[-1])))) == 0:
                        classe.append(deg[-1])
                        deg.pop(-1)
                        break
            if hex_colors:
                return dict(zip(rainbow(k), value))
            else:
                return value

        p = MixedIntegerLinearProgram(maximization=True)
        color = p.new_variable(dim=2)
        # a vertex has exactly one color
        [p.add_constraint(sum([color[v][i] for i in xrange(k)]), min=1, max=1)
             for v in g.vertices()]
        # adjacent vertices have different colors
        [p.add_constraint(color[u][i] + color[v][i], max=1)
             for (u, v) in g.edge_iterator(labels=None)
                 for i in xrange(k)]
        # anything is good as an objective value as long as it is satisfiable
        p.add_constraint(color[g.vertex_iterator().next()][0],  max=1, min=1)
        p.set_objective(color[g.vertex_iterator().next()][0])
        p.set_binary(color)
        from sage.numerical.mip import MIPSolverException
        try:
            if value_only:
                p.solve(objective_only=True, log=log)
                return True
            else:
                chi = p.solve(log=log)
        except MIPSolverException:
            return False

        color = p.get_values(color)
        # builds the color classes
        classes = [[] for i in xrange(k)]
        [classes[i].append(v)
             for i in xrange(k)
                 for v in g.vertices()
                     if color[v][i] == 1]
        if hex_colors:
            return dict(zip(rainbow(len(classes)), classes))
        else:
            return classes

def edge_coloring(g, value_only=False, vizing=False, hex_colors=False, log=0):
    r"""
    Properly colors the edges of a graph. See the URL
    http://en.wikipedia.org/wiki/Edge_coloring for further details on
    edge coloring.

    INPUT:

    - ``g`` -- a graph.

    - ``value_only`` -- (default: ``False``):

      - When set to ``True``, only the chromatic index is returned.

      - When set to ``False``, a partition of the edge set into
        matchings is returned if possible.

    - ``vizing`` -- (default: ``False``):

      - When set to ``True``, tries to find a `\Delta + 1`-edge-coloring,
        where `\Delta` is equal to the maximum degree in the graph.

      - When set to ``False``, tries to find a `\Delta`-edge-coloring,
        where `\Delta` is equal to the maximum degree in the graph. If
        impossible, tries to find and returns a `\Delta + 1`-edge-coloring.
        This implies that ``value_only=False``.

    - ``hex_colors`` -- (default: ``False``) when set to ``True``, the
      partition returned is a dictionary whose keys are colors and whose
      values are the color classes (ideal for plotting).

    - ``log`` -- (default: ``0``) as edge-coloring is an `NP`-complete
      problem, this function may take some time depending on the graph. Use
      ``log`` to define the level of verbosity you wantfrom the linear
      program solver. By default ``log=0``, meaning that there will be no
      message printed by the solver.

    OUTPUT:

    In the following, `\Delta` is equal to the maximum degree in the graph
    ``g``.

    - If ``vizing=True`` and ``value_only=False``, return a partition of
      the edge set into `\Delta + 1` matchings.

    - If ``vizing=False`` and ``value_only=True``, return the chromatic index.

    - If ``vizing=False`` and ``value_only=False``, return a partition of
      the edge set into the minimum number of matchings.

    - If ``vizing=True`` and ``value_only=True``, should return something,
      but mainly you are just trying to compute the maximum degree of the
      graph, and this is not the easiest way. By Vizing's theorem, a graph
      has a chromatic index equal to `\Delta` or to `\Delta + 1`.

    EXAMPLE::

       sage: from sage.graphs.graph_coloring import edge_coloring
       sage: g = graphs.PetersenGraph()
       sage: edge_coloring(g, value_only=True) # optional - requires GLPK or CBC
       4

    Complete graphs are colored using the linear-time round-robin coloring::

       sage: from sage.graphs.graph_coloring import edge_coloring
       sage: len(edge_coloring(graphs.CompleteGraph(20)))
       19
    """
    from sage.numerical.mip import MixedIntegerLinearProgram
    from sage.plot.colors import rainbow
    from sage.numerical.mip import MIPSolverException

    if g.is_clique():
        if value_only:
            return g.order() if g.order() % 2 == 0 else g.order() + 1
        vertices = g.vertices()
        r = round_robin(g.order())
        classes = [[] for v in g]
        if g.order() % 2 == 0 and not vizing:
            classes.pop()
        for (u, v, c) in r.edge_iterator():
            classes[c].append((vertices[u], vertices[v]))
        if hex_colors:
            return dict(zip(rainbow(len(classes)), classes))
        else:
            return classes

    p = MixedIntegerLinearProgram(maximization=True)
    color = p.new_variable(dim=2)
    obj = {}
    k = max(g.degree())
    # reorders the edge if necessary...
    R = lambda x: x if (x[0] <= x[1]) else (x[1], x[0], x[2])
    # Vizing's coloring uses Delta + 1 colors
    if vizing:
        value_only = False
        k += 1
    #  A vertex can not have two incident edges with the same color.
    [p.add_constraint(
            sum([color[R(e)][i] for e in g.edges_incident(v)]), max=1)
                for v in g.vertex_iterator()
                    for i in xrange(k)]
    # an edge must have a color
    [p.add_constraint(sum([color[R(e)][i] for i in xrange(k)]), max=1, min=1)
         for e in g.edge_iterator()]
    # anything is good as an objective value as long as it is satisfiable
    e = g.edge_iterator().next()
    p.set_objective(color[R(e)][0])
    p.set_binary(color)
    try:
        if value_only:
            p.solve(objective_only=True, log=log)
        else:
            chi = p.solve(log=log)
    except:
        if value_only:
            return k + 1
        else:
            # if the coloring with Delta colors fails, tries Delta + 1
            return edge_coloring(g, vizing=True, hex_colors=hex_colors, log=log)
    if value_only:
        return k
    # Builds the color classes
    color = p.get_values(color)
    classes = [[] for i in xrange(k)]
    [classes[i].append(e)
         for e in g.edge_iterator()
             for i in xrange(k)
                 if color[R(e)][i] == 1]
    # if needed, builds a dictionary from the color classes adding colors
    if hex_colors:
        return dict(zip(rainbow(len(classes)), classes))
    else:
        return classes

def round_robin(n):
    r"""
    Computes a round-robin coloring of the complete graph on `n` vertices.

    A round-robin coloring of the complete graph `G` on `2n` vertices
    (`V = [0, \dots, 2n - 1]`) is a proper coloring of its edges such that
    the edges with color `i` are all the `(i + j, i - j)` plus the
    edge `(2n - 1, i)`.

    If `n` is odd, one obtain a round-robin coloring of the complete graph
    through the round-robin coloring of the graph with `n + 1` vertices.

    INPUT:

    - ``n`` -- the number of vertices in the complete graph.

    OUTPUT:

    - A ``CompleteGraph`` with labelled edges such that the label of each
      edge is its color.

    EXAMPLES::

        sage: from sage.graphs.graph_coloring import round_robin
        sage: round_robin(3).edges()
        [(0, 1, 2), (0, 2, 1), (1, 2, 0)]

    ::

        sage: round_robin(4).edges()
        [(0, 1, 2), (0, 2, 1), (0, 3, 0), (1, 2, 0), (1, 3, 1), (2, 3, 2)]


    For higher orders, the coloring is still proper and uses the expected
    number of colors.

    ::

        sage: g = round_robin(9)
        sage: sum([Set([e[2] for e in g.edges_incident(v)]).cardinality() for v in g]) == 2*g.size()
        True
        sage: Set([e[2] for e in g.edge_iterator()]).cardinality()
        9

    ::

        sage: g = round_robin(10)
        sage: sum([Set([e[2] for e in g.edges_incident(v)]).cardinality() for v in g]) == 2*g.size()
        True
        sage: Set([e[2] for e in g.edge_iterator()]).cardinality()
        9
    """
    if n <= 1:
        raise ValueError("There must be at least two vertices in the graph.")
    mod = lambda x, y: x - y*(x // y)
    if n % 2 == 0:
        g = GraphGenerators().CompleteGraph(n)
        for i in xrange(n - 1):
            g.set_edge_label(n - 1, i, i)
            for j in xrange(1, (n - 1) // 2 + 1):
                g.set_edge_label(mod(i - j, n - 1), mod(i + j, n - 1), i)
        return g
    else:
        g = round_robin(n + 1)
        g.delete_vertex(n)
        return g

class Test:
    r"""
    This class performs randomized testing for all_graph_colorings.
    Since everything else in this file is derived from
    all_graph_colorings, this is a pretty good randomized tester for
    the entire file.  Note that for a graph `G`, ``G.chromatic_polynomial()``
    uses an entirely different algorithm, so we provide a good,
    independent test.
    """

    def random(self,tests = 1000):
        r"""
        Calls ``self.random_all_graph_colorings()``.  In the future, if
        other methods are added, it should call them, too.

        TESTS::

            sage: from sage.graphs.graph_coloring import Test
            sage: Test().random(1)
        """
        self.random_all_graph_colorings(tests)

    def random_all_graph_colorings(self,tests = 1000):
        r"""
        Verifies the results of ``all_graph_colorings()`` in three ways:

        #. all colorings are unique

        #. number of m-colorings is `P(m)` (where `P` is the chromatic
           polynomial of the graph being tested)

        #. colorings are valid -- that is, that no two vertices of
           the same color share an edge.

        TESTS::

            sage: from sage.graphs.graph_coloring import Test
            sage: Test().random_all_graph_colorings(1)
        """
        from sage.all import Set

        G = GraphGenerators().RandomGNP(10,.5)
        Q = G.chromatic_polynomial()
        N = G.chromatic_number()
        m = N

        S = Set([])

        for C in all_graph_colorings(G, m):
            parts = [C[k] for k in C]
            for P in parts:
                l = len(P)
                for i in range(l):
                    for j in range(i+1,l):
                        if G.has_edge(P[i],P[j]):
                            raise RuntimeError, "Coloring Failed."

            #make the dict into a set for quick uniqueness checking
            S+= Set([Set([(k,tuple(C[k])) for k in C])])

        if len(S) != Q(m):
            raise RuntimeError, "Incorrect number of unique colorings!"
