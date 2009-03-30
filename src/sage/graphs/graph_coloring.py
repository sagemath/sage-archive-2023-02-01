"""
Graph Coloring Functions

AUTHORS:
    -- Tom Boothby   (2008-02-21): Initial version
    -- Carlo Hamalainen (2009-03-28): minor change: switch to C++ DLX solver
"""

#*****************************************************************************
#           Copyright (C) 2008 Tom Boothby <boothby@u.washington.edu>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.matrices.dlxcpp import DLXCPP
from sage.all import Matrix, vector, QQ
from sage.plot.plot import rainbow
from chrompoly import chromatic_polynomial
from graph_generators import GraphGenerators

def all_graph_colorings(G,n,count_only=False):
    """
    Computes all n-colorings of the graph G by casting the graph
    coloring problem into an exact cover problem, and passing this
    into an implementation of the Dancing Links algorithm described
    by Knuth (who attributes the idea to Hitotumatu and Noshita).

    The construction works as follows:
     (columns)
      * The first |V| columns correspond to a vertex -- a 1 in this
           column indicates that that vertex has a color.
      * After those |V| columns, we add n*|E| columns -- a 1 in
           these columns indicate that a particular edge is
           incident to a vertex with a certain color.

     (rows)
      * For each vertex, add n rows; one for each color c.  Place
           a 1 in the column corresponding to the vertex, and a 1
           in the appropriate column foreach edge incident to the
           vertex, indicating that that edge is incident to the
           color c.
      * If n > 2, the above construction cannot be exactly covered
           since each edge will be incident to only two vertices
           (and hence two colors) - so we add n*|E| rows, each one
           containing a 1 for each of the n*|E| columns.  These
           get added to the cover solutions "for free" during the
           backtracking.

    Note that this construction results in n*|V| + 2*n*|E| + n*|E|
    entries in the matrix.  The Dancing Links algorithm uses a
    sparse representation, so if the graph is simple, |E| <= |V|^2
    and n <= |V|, this construction runs in O(|V|^3) time.
    Back-conversion to a coloring solution is a simple scan of the
    solutions, which will contain |V| + (n-2)*|E| entries,  so
    runs in O(|V|^3) time also.  For most graphs, the conversion
    will be much faster -- for example, a planar graph will be
    transformed for 4-coloring in linear time since |E| = O(|V|).

    REFERENCES:
        http://www-cs-staff.stanford.edu/~uno/papers/dancing-color.ps.gz

    EXAMPLES:
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

def first_coloring(G,n=0):
    """
    Given a graph, and optionally a natural number n, returns
    the first coloring we find with at least n colors.

    EXAMPLES:
        sage: from sage.graphs.graph_coloring import first_coloring
        sage: G = Graph({0:[1,2,3],1:[2]})
        sage: first_coloring(G,3)
        {'#00ff00': [1, 3], '#ff0000': [0], '#0000ff': [2]}
    """
    o = G.order()
    for m in range(n,o+1):
        for C in all_graph_colorings(G,m):
            return C

def number_of_n_colorings(G,n):
    """
    Given a graph G and a natural number n, returns the number of
    n-colorings of the graph.

    EXAMPLES:
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
    """
    Returns the number of n-colorings of the graph G for n from
    0 to |V|.

    EXAMPLES:
        sage: from sage.graphs.graph_coloring import numbers_of_colorings
        sage: G = Graph({0:[1,2,3],1:[2]})
        sage: numbers_of_colorings(G)
        [0, 0, 0, 12, 72]
    """
    o = G.order()
    return [number_of_n_colorings(G,i) for i in range(0,o+1)]

def chromatic_number(G):
    """
    Returns the minimal number of colors needed to color the
    vertices of the graph G.

    EXAMPLES:
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
        m = max([len(c) for c in G.cliques()])
    if m >= o-1: #marginal improvement... if there's an o-1 clique and not an o clique, don't waste our time coloring.
        return m
    for n in range(m,o+1):
        for C in all_graph_colorings(G,n):
            return n

class Test:
    """
    This class performs randomized testing for all_graph_colorings.
    Since everything else in this file is derived from
    all_graph_colorings, this is a pretty good randomized tester for
    the entire file.  Note that for a graph G, G.chromatic_polynomial()
    uses an entirely different algorithm, so we provide a good,
    independent test.
    """

    def random(self,tests = 1000):
        """
        Calls self.random_all_graph_colorings().  In the future, if
        other methods are added, it should call them, too.

        TESTS:
            sage: from sage.graphs.graph_coloring import Test
            sage: Test().random(1)
        """
        self.random_all_graph_colorings(tests)

    def random_all_graph_colorings(self,tests = 1000):
        """
        Verifies the results of all_graph_colorings in three ways:
            1) all colorings are unique
            2) number of m-colorings is P(m) (where P is the chromatic
               polynomial of the graph being tested)
            3) colorings are valid -- that is, that no two vertices of
               the same color share an edge.

        TESTS:
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
