r"""
Hyperbolicity

**Definition** :

    The hyperbolicity `\delta` of a graph `G` has been defined by Gromov
    [Gromov87]_ as follows (we give here the so-called 4-points condition):

      Let `a, b, c, d` be vertices of the graph, let `S_1`, `S_2` and `S_3` be
      defined by

      .. MATH::

          S_1 = dist(a, b) + dist(b, c)\\
          S_2 = dist(a, c) + dist(b, d)\\
          S_3 = dist(a, d) + dist(b, c)\\

      and let `M_1` and `M_2` be the two largest values among `S_1`, `S_2`, and
      `S_3`. We define `hyp(a, b, c, d) = M_1 - M_2`, and the hyperbolicity
      `\delta(G)` of the graph is the maximum of `hyp` over all possible
      4-tuples `(a, b, c, d)` divided by 2. That is, the graph is said
      `\delta`-hyperbolic when

      .. MATH::

          \delta(G) = \frac{1}{2}\max_{a,b,c,d\in V(G)}hyp(a, b, c, d)

      (note that `hyp(a, b, c, d)=0` whenever two elements among `a,b,c,d` are equal)

**Some known results** :

    - Trees and cliques are `0`-hyperbolic

    - `n\times n` grids are `n-1`-hyperbolic

    - Cycles are approximately `n/4`-hyperbolic

    - Chordal graphs are `\leq 1`-hyperbolic

    Besides, the hyperbolicity of a graph is the maximum over all its
    biconnected components.

**Algorithms and complexity** :

    The time complexity of the naive implementation (i.e. testing all 4-tuples)
    is `O( n^4 )`, and an algorithm with time complexity `O(n^{3.69})` has been
    proposed in [FIV12]_.  This remains very long for large-scale graphs, and
    much harder to implement.


    An improvement over the naive algorithm has been proposed in [CCL12]_, and
    is implemented in the current module. Like the naive algorithm, it has
    complexity `O(n^4)` but behaves much better in practice. It uses the
    following fact :

      Assume that `S_1 = dist(a, b) + dist(c, d)` is the largest among
      `S_1,S_2,S_3`. We have

      .. MATH::

        S_2 + S_3 =& dist(a, c) + dist(b, d) + dist(a, d) + dist(b, c)\\
        =& [ dist(a, c) + dist(b, c) ] + [ dist(a, d) + dist(b, d)]\\
        \geq &dist(a,b) + dist(a,b)\\
        \geq &2dist(a,b)\\

      Now, since `S_1` is the largest sum, we have

      .. MATH::

        hyp(a, b, c, d) =& S_1 - \max\{S_2, S_3\}\\
        \leq& S_1 - \frac{S_2+ S_3}{2}\\
        \leq& S_1 - dist(a, b)\\
        =& dist(c, d)\\

      We obtain similarly that `hyp(a, b, c, d) \leq dist(a, b)`.  Consequently,
      in the implementation, we ensure that `S_1` is larger than `S_2` and `S_3`
      using an ordering of the pairs by decreasing lengths. Furthermore, we use
      the best value `h` found so far to cut exploration.

    The worst case time complexity of this algorithm is `O(n^4)`, but it
    performs very well in practice since it cuts the search space.  This
    algorithm can be turned into an approximation algorithm since at any step of
    its execution we maintain an upper and a lower bound. We can thus stop
    execution as soon as a multiplicative approximation factor or an additive
    one is proven.

TODO:

- Add exact methods for the hyperbolicity of chordal graphs

- Add method for partitioning the graph with clique separators

**This module contains the following functions**

At Python level :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~hyperbolicity` | Return the hyperbolicity of the graph or an approximation of this value.
    :meth:`~hyperbolicity_distribution` | Return the hyperbolicity distribution of the graph or a sampling of it.

REFERENCES:

.. [CCL12] N. Cohen, D. Coudert, and A. Lancin. Exact and approximate algorithms
   for computing the hyperbolicity of large-scale graphs.  Research Report
   RR-8074, Sep. 2012. [`<http://hal.inria.fr/hal-00735481>`_].

.. [FIV12] H. Fournier, A. Ismail, and A. Vigneron. Computing the Gromov
   hyperbolicity of a discrete metric space. ArXiv, Tech. Rep. arXiv:1210.3323,
   Oct. 2012. [`<http://arxiv.org/abs/1210.3323>`_].

.. [Gromov87] M. Gromov. Hyperbolic groups. Essays in Group Theory, 8:75--263,
   1987.

AUTHORS:

- David Coudert (2012): initial version, exact and approximate algorithm,
  distribution, sampling


Methods
-------
"""

###############################################################################
#           Copyright (C) 2012 David Coudert <david.coudert@inria.fr>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###############################################################################

# imports
from sage.graphs.graph import Graph
from sage.graphs.distances_all_pairs cimport c_distances_all_pairs
from sage.rings.arith import binomial
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from sage.functions.other import floor
from sage.data_structures.bitset import Bitset
from libc.stdint cimport uint16_t, uint32_t, uint64_t
include "sage/ext/stdsage.pxi"


# Defining a pair of vertices as a C struct
ctypedef struct pair:
    uint16_t s
    uint16_t t


######################################################################
# Speedup functions
######################################################################

def _my_subgraph(G, vertices, relabel=False, return_map=False):
    r"""
    Return the subgraph containing the given vertices

    This method considers only the connectivity. Therefore, edge labels are
    ignored as well as any other decoration of the graph (vertex position,
    etc.).

    If ``relabel`` is ``True``, the vertices of the new graph are relabeled with
    integers in the range '0\cdots |vertices|-1'. The relabeling map is returned
    if ``return_map`` is also ``True``.

    TESTS:

    Giving anything else than a Graph::

        sage: from sage.graphs.hyperbolicity import _my_subgraph as mysub
        sage: mysub([],[])
        Traceback (most recent call last):
        ...
        ValueError: The input parameter must be a Graph.

    Subgraph of a PetersenGraph::

        sage: from sage.graphs.hyperbolicity import _my_subgraph as mysub
        sage: H = mysub(graphs.PetersenGraph(), [0,2,4,6])
        sage: H.edges(labels=None)
        [(0, 4)]
        sage: H.vertices()
        [0, 2, 4, 6]
    """
    if not isinstance(G,Graph):
        raise ValueError("The input parameter must be a Graph.")
    H = Graph()
    if not vertices:
        return (H,{}) if (relabel and return_map) else H

    if relabel:
        map = dict(zip(iter(vertices),xrange(len(vertices))))
    else:
        map = dict(zip(iter(vertices),iter(vertices)))

    B = {}
    for v in G.vertex_iterator():
        B[v] = False
    for v in vertices:
        B[v] = True
        H.add_vertex(map[v])

    for u in vertices:
        for v in G.neighbor_iterator(u):
            if B[v]:
                H.add_edge(map[u],map[v])

    return (H,map) if (relabel and return_map) else H


######################################################################
# Building blocks
######################################################################

cdef inline int __hyp__(unsigned short ** distances, int a, int b, int c, int d):
    """
    Return the hyperbolicity of the given 4-tuple.
    """
    cdef int S1, S2, S3, h
    S1 = distances[a][b] + distances[c][d]
    S2 = distances[a][c] + distances[b][d]
    S3 = distances[a][d] + distances[b][c]
    if S1 >= S2:
        if S2 > S3:
            h = S1-S2
        else:
            h = abs(S1-S3)
    else:
        if S1 > S3:
            h = S2-S1
        else:
            h = abs(S2-S3)
    return h

######################################################################
# Basic algorithm for the hyperbolicity
######################################################################

cdef tuple __hyperbolicity_basic_algorithm__(int N,
                                             unsigned short **  distances,
                                             verbose = False):
    """
    Returns **twice** the hyperbolicity of a graph, and a certificate.

    This method implements the basic algorithm for computing the hyperbolicity
    of a graph, that is iterating over all 4-tuples. See the module's
    documentation for more details.

    INPUTS:

    - ``N`` -- number of vertices of the graph.

    - ``distances`` -- path distance matrix (see the distance_all_pairs module).

    - ``verbose`` -- (default: ``False``) is boolean. Set to True to display
      some information during execution.

    OUTPUTS:

    This function returns a tuple ( h, certificate ), where:

    - ``h`` -- the maximum computed value over all 4-tuples, and so is twice the
      hyperbolicity of the graph. If no such 4-tuple is found, -1 is returned.

    - ``certificate`` -- 4-tuple of vertices maximizing the value `h`. If no
      such 4-tuple is found, the empty list [] is returned.
    """
    cdef int a, b, c, d, S1, S2, S3, hh, h_LB
    cdef list certificate

    h_LB = -1
    for 0 <= a < N-3:
        for a < b < N-2:
            for b < c < N-1:
                for c < d < N:

                    # We compute the hyperbolicity of the 4-tuple
                    hh = __hyp__(distances, a, b, c, d)

                    # We compare the value with previously known bound
                    if hh > h_LB:
                        h_LB = hh
                        certificate = [a, b, c, d]

                        if verbose:
                            print 'New lower bound:',ZZ(hh)/2

    # Last, we return the computed value and the certificate
    if h_LB != -1:
        return ( h_LB, certificate )
    else:
        return ( -1, [] )


######################################################################
# Decomposition methods
######################################################################

def elimination_ordering_of_simplicial_vertices(G, max_degree=4, verbose=False):
    r"""
    Return an elimination ordering of simplicial vertices.

    An elimination ordering of simplicial vertices is an elimination ordering of
    the vertices of the graphs such that the induced subgraph of their neighbors
    is a clique. More precisely, as long as the graph has a vertex ``u`` such
    that the induced subgraph of its neighbors is a clique, we remove ``u`` from
    the graph, add it to the elimination ordering (list of vertices), and
    repeat. This method is inspired from the decomposition of a graph by
    clique-separators.

    INPUTS:

    - ``G`` -- a Graph

    - ``max_degree`` -- (default: 4) maximum degree of the vertices to consider.
      The running time of this method depends on the value of this parameter.

    - ``verbose`` -- (default: ``False``) is boolean set to ``True`` to display
      some information during execution.

    OUTPUT:

    - ``elim`` -- A ordered list of vertices such that vertex ``elim[i]`` is
      removed before vertex ``elim[i+1]``.

    TESTS:

    Giving anything else than a Graph::

        sage: from sage.graphs.hyperbolicity import elimination_ordering_of_simplicial_vertices
        sage: elimination_ordering_of_simplicial_vertices([])
        Traceback (most recent call last):
        ...
        ValueError: The input parameter must be a Graph.

    Giving two small bounds on degree::

        sage: from sage.graphs.hyperbolicity import elimination_ordering_of_simplicial_vertices
        sage: elimination_ordering_of_simplicial_vertices(Graph(), max_degree=0)
        Traceback (most recent call last):
        ...
        ValueError: The parameter max_degree must be > 0.

    Giving a graph built from a bipartite graph plus an edge::

        sage: G = graphs.CompleteBipartiteGraph(2,10)
        sage: G.add_edge(0,1)
        sage: from sage.graphs.hyperbolicity import elimination_ordering_of_simplicial_vertices
        sage: elimination_ordering_of_simplicial_vertices(G)
        [2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 11]
        sage: elimination_ordering_of_simplicial_vertices(G,max_degree=1)
        []
    """
    if verbose:
        print 'Entering elimination_ordering_of_simplicial_vertices'

    if not isinstance(G,Graph):
        raise ValueError("The input parameter must be a Graph.")
    elif max_degree < 1:
        raise ValueError("The parameter max_degree must be > 0.")

    # We make a copy of the graph. We use a NetworkX graph since modifications
    # are a bit faster this way.
    import networkx
    ggnx = networkx.empty_graph()
    for u,v in G.edge_iterator(labels=None):
        ggnx.add_edge(u,v)

    from sage.combinat.combination import Combinations
    cdef list elim = []
    cdef set L = set()

    # We identify vertices of degree at most max_degree
    for u,d in ggnx.degree_iter():
        if d<=max_degree:
            L.add(u)

    while L:
        # We pick up a vertex and check if the induced subgraph of its neighbors
        # is a clique. If True, we record it, remove it from the graph, and
        # update the list of vertices of degree at most max_degree.
        u = L.pop()
        X = ggnx.neighbors(u)
        if all(ggnx.has_edge(v,w) for v,w in Combinations(X,2).list()):
            elim.append(u)
            ggnx.remove_node(u)
            for v,d in ggnx.degree_iter(X):
                if d<=max_degree:
                    L.add(v)

    if verbose:
        print 'Propose to eliminate',len(elim),'of the',G.num_verts(),'vertices'
        print 'End elimination_ordering_of_simplicial_vertices'

    return elim


######################################################################
# Greedy dominating set
######################################################################

def _greedy_dominating_set(H):
    r"""
    Returns a greedy approximation of a dominating set
    """
    V = sorted([(d,u) for u,d in H.degree_iterator(None,True)],reverse=True)
    DOM = []
    seen = set()
    for _,u in V:
        if not u in seen:
            seen.add(u)
            DOM.append(u)
            seen.update(H.neighbor_iterator(u))
    return DOM

######################################################################
# Compute the hyperbolicity using a path decreasing length ordering
######################################################################

cdef inline _invert_cells(pair * tab, uint32_t idxa, uint32_t idxb):
    cdef pair tmp = tab[idxa]
    tab[idxa] = tab[idxb]
    tab[idxb] = tmp


cdef _order_pairs_according_elimination(elim, int N, pair *pairs, uint32_t nb_pairs, uint32_t *last_pair):
    r"""
    Re-order the pairs of vertices according the set of eliminated vertices.

    We put pairs of vertices with an extremity in ``elim`` at the end of the
    array of pairs.  We record the positions of the first pair of vertices in
    ``elim``.  If ``elim`` is empty, the ordering is unchanged and we set
    ``last_pair=nb_pairs``.
    """
    cdef uint32_t j, jmax

    jmax = nb_pairs-1
    if nb_pairs<=1 or not elim:
        last_pair[0] = nb_pairs
    else:
        B = Bitset(iter(elim))
        j = 0
        while j<jmax:

            if pairs[j].s in B or pairs[j].t in B:
                while (pairs[jmax].s in B or pairs[jmax].t in B) and (j<jmax):
                    jmax -= 1

                if j<jmax:
                    _invert_cells(pairs, j, jmax)

                if jmax>0:
                    jmax -= 1

            else: # This pair is at a correct position.
                j += 1

        # We record the position of the first pair of vertices in elim
        last_pair[0] = jmax+1

cdef tuple __hyperbolicity__(int N,
                             unsigned short **  distances,
                             int D,
                             int h_LB,
                             float approximation_factor,
                             float additive_gap,
                             elim,
                             verbose = False):
    """
    Return the hyperbolicity of a graph.

    This method implements the exact and the approximate algorithms proposed in
    [CCL12]_. See the module's documentation for more details.

    INPUTS:

    - ``N`` -- number of vertices of the graph

    - ``distances`` -- path distance matrix

    - ``D`` -- diameter of the graph

    - ``h_LB`` -- lower bound on the hyperbolicity

    - ``approximation_factor`` -- When the approximation factor is set to some
      value larger than 1.0, the function stop computations as soon as the ratio
      between the upper bound and the best found solution is less than the
      approximation factor. When the approximation factor is 1.0, the problem is
      solved optimaly.

     - ``additive_gap`` -- When sets to a positive number, the function stop
       computations as soon as the difference between the upper bound and the
       best found solution is less than additive gap. When the gap is 0.0, the
       problem is solved optimaly.

    - ``elim`` -- list of vertices that should not be considered during
      computations.

    - ``verbose`` -- (default: ``False``) is boolean set to ``True`` to display
      some information during execution

    OUTPUTS:

    This function returns a tuple ( h, certificate, h_UB ), where:

    - ``h`` -- is an integer. When 4-tuples with hyperbolicity larger or equal
     to `h_LB are found, h is the maximum computed value and so twice the
     hyperbolicity of the graph. If no such 4-tuple is found, it returns -1.

    - ``certificate`` -- is a list of vertices. When 4-tuples with hyperbolicity
      larger that h_LB are found, certificate is the list of the 4 vertices for
      which the maximum value (and so the hyperbolicity of the graph) has been
      computed. If no such 4-tuple is found, it returns the empty list [].

    - ``h_UB`` -- is an integer equal to the proven upper bound for `h`. When
      ``h == h_UB``, the returned solution is optimal.
    """
    cdef int i, j, l, l1, l2, h, hh, h_UB, a, b, c, d, S1, S2, S3
    cdef uint32_t x, y
    cdef dict distr = {}
    cdef list certificate = []

    # ==> Allocates and fills nb_pairs_of_length
    #
    # nb_pairs_of_length[d] is the number of pairs of vertices at distance d
    cdef uint32_t * nb_pairs_of_length = <uint32_t *>sage_calloc(D+1,sizeof(uint32_t))
    if nb_pairs_of_length == NULL:
        sage_free(nb_pairs_of_length)
        raise MemoryError

    for 0 <= i < N:
        for i+1 <= j < N:
            nb_pairs_of_length[ distances[i][j] ] += 1

    # ==> Allocates pairs_of_length
    #
    # pairs_of_length[d] is the list of pairs of vertices at distance d
    cdef pair ** pairs_of_length = <pair **>sage_malloc(sizeof(pair *)*(D+1))
    if pairs_of_length == NULL:
        sage_free(nb_pairs_of_length)
        raise MemoryError

    # ==> Allocates cpt_pairs
    #
    # (temporary variable used to fill pairs_of_length)
    cdef uint32_t * cpt_pairs = <uint32_t *>sage_calloc(D+1,sizeof(uint32_t))
    if cpt_pairs == NULL:
        sage_free(nb_pairs_of_length)
        sage_free(pairs_of_length)
        raise MemoryError

    # ==> Allocates pairs_of_length[d] for all d
    for 1 <= i <= D:
        pairs_of_length[i] = <pair *>sage_malloc(sizeof(pair)*nb_pairs_of_length[i])

        if nb_pairs_of_length[i] > 0 and pairs_of_length[i] == NULL:
            while i>1:
                i -= 1
                sage_free(pairs_of_length[i])
            sage_free(nb_pairs_of_length)
            sage_free(pairs_of_length)
            sage_free(cpt_pairs)
            raise MemoryError

    # ==> Fills pairs_of_length[d] for all d
    for 0 <= i < N:
        for i+1 <= j < N:
            l = distances[i][j]
            pairs_of_length[ l ][ cpt_pairs[ l ] ].s = i
            pairs_of_length[ l ][ cpt_pairs[ l ] ].t = j
            cpt_pairs[ l ] += 1

    sage_free(cpt_pairs)

    if verbose:
        print "Current 2 connected component has %d vertices and diameter %d" %(N,D)
        print "Paths length distribution:", [ (l, nb_pairs_of_length[l]) for l in range(1, D+1) ]

    # ==> Allocates last_pair
    #
    # We change the ordering of the pairs to avoid considering pairs touching
    # eliminated vertices (i.e., vertices in elim). We store in last_pair[l] the
    # index of the last pair of length l to consider.
    cdef uint32_t * last_pair = <uint32_t *>sage_malloc(sizeof(uint32_t)*(D+1))
    if last_pair == NULL:
        for 1 <= i <= D:
            sage_free(pairs_of_length[i])
        sage_free(nb_pairs_of_length)
        sage_free(pairs_of_length)
        sage_free(cpt_pairs)
        raise MemoryError

    for 1 <= l <= D:
        _order_pairs_according_elimination(elim, N, pairs_of_length[l], nb_pairs_of_length[l], last_pair+l)

    approximation_factor = min(approximation_factor, D)
    additive_gap = min(additive_gap, D)

    # We create the list of triples (sum,length1,length2) sorted in
    # co-lexicographic order: decreasing by sum, decreasing by length2,
    # decreasing length1. This is to ensure a valid ordering for S1, to avoid
    # some tests, and to ease computation of bounds.
    cdef list triples = []
    for i in range(D,0,-1):
        for j in range(D,i-1,-1):
            triples.append((i+j,j,i))

    # We use some short-cut variables for efficiency
    cdef pair * pairs_of_length_l1
    cdef pair * pairs_of_length_l2
    cdef uint32_t nb_pairs_of_length_l1, nb_pairs_of_length_l2
    h = h_LB
    h_UB = D
    cdef int STOP = 0
    for S1, l1, l2 in triples:

        # If we cannot improve further, we stop
        #
        # See the module's documentation for an proof that this cut is
        # valid. Remember that the triples are sorted in a specific order.
        if l2 <= h:
            h_UB = h
            break

        if h_UB > l2:
            h_UB = l2

            if verbose:
                print "New upper bound:",ZZ(h_UB)/2

            # Termination if required approximation is found
            if certificate and ( (h_UB <= h*approximation_factor) or (h_UB-h <= additive_gap) ):
                break

        pairs_of_length_l1 = pairs_of_length[l1]
        nb_pairs_of_length_l1 = last_pair[l1]
        x = 0
        while x < nb_pairs_of_length_l1:
            a = pairs_of_length_l1[x].s
            b = pairs_of_length_l1[x].t

            # If we cannot improve further, we stop
            #
            # See the module's documentation for an proof that this cut is
            # valid.
            if l2 <= h:
                STOP = 1
                break

            # We do not want to test pairs of pairs twice if l1 == l2
            elif l1 == l2:
                y = x+1
            else:
                y = 0

            pairs_of_length_l2 = pairs_of_length[l2]
            nb_pairs_of_length_l2 = last_pair[l2]
            while y < nb_pairs_of_length_l2:
                c = pairs_of_length_l2[y].s
                d = pairs_of_length_l2[y].t

                # If two points are equal, the value will be 0, so we skip the
                # test.
                if a == c or a == d or b == c or b == d:
                    y += 1
                    continue

                # We compute the hyperbolicity of the 4-tuple. We have S1 = l1 +
                # l2, and the order in which pairs are visited allow us to claim
                # that S1 = max( S1, S2, S3 ). If at some point S1 is not the
                # maximum value, the order ensures that the maximum value has
                # previously been checked.
                S2 = distances[a][c] + distances[b][d]
                S3 = distances[a][d] + distances[b][c]
                if S2 > S3:
                    hh = S1 - S2
                else:
                    hh = S1 - S3

                if h < hh or not certificate:
                    # We update current bound on the hyperbolicity and the
                    # search space
                    h = hh
                    certificate = [a, b, c, d]

                    if verbose:
                        print "New lower bound:",ZZ(hh)/2

                    # Termination if required approximation is found
                    if (h_UB <= h*approximation_factor) or (h_UB-h <= additive_gap):
                        STOP = 1
                        break

                # We go for the next pair c-d
                y += 1
                # We cut current exploration if we know we can not improve lower bound
                #
                # See the module's documentation for an proof that this cut is
                # valid.
                if l2 <= h:
                    STOP = 1
                    h_UB = h
                    break

            if STOP:
                break

            # We go for the next pair a-b
            x += 1

        if STOP:
            break

    # We now free the memory
    sage_free(nb_pairs_of_length)
    for 1 <= i <= D:
        sage_free(pairs_of_length[i])
    sage_free(pairs_of_length)
    sage_free(last_pair)

    # Last, we return the computed value and the certificate
    if len(certificate) == 0:
        return ( -1, [], h_UB )
    else:
        return (h, certificate, h_UB)


def hyperbolicity(G, algorithm='cuts', approximation_factor=None, additive_gap=None, verbose = False):
    r"""
    Return the hyperbolicity of the graph or an approximation of this value.

    The hyperbolicity of a graph has been defined by Gromov [Gromov87]_ as
    follows: Let `a, b, c, d` be vertices of the graph, let `S_1 = dist(a, b) +
    dist(b, c)`, `S_2 = dist(a, c) + dist(b, d)`, and `S_3 = dist(a, d) +
    dist(b, c)`, and let `M_1` and `M_2` be the two largest values among `S_1`,
    `S_2`, and `S_3`. We have `hyp(a, b, c, d) = |M_1 - M_2|`, and the
    hyperbolicity of the graph is the maximum over all possible 4-tuples `(a, b,
    c, d)` divided by 2. The worst case time complexity is in `O( n^4 )`.

    See the documentation of :mod:`sage.graphs.hyperbolicity` for more
    information.

    INPUT:

    - ``G`` -- a Graph

    - ``algorithm`` -- (default: ``'cuts'``) specifies the algorithm to use
      among:

          - ``'basic'`` is an exhaustive algorithm considering all possible
            4-tuples and so have time complexity in `O(n^4)`.

          - ``'cuts'`` is an exact algorithm proposed in [CCL12_]. It considers
            the 4-tuples in an ordering allowing to cut the search space as soon
            as a new lower bound is found (see the module's documentation). This
            algorithm can be turned into a approximation algorithm.

          - ``'cuts+'`` is an additive constant approximation algorithm. It
            proceeds by first removing the simplicial vertices and then applying
            the ``'cuts'`` algorithm on the remaining graph, as documented in
            [CCL12_]. By default, the additive constant of the approximation is
            one. This value can be increased by setting the ``additive_gap`` to
            the desired value, provide ``additive_gap \geq 1``. In some cases,
            the returned result is proven optimal. However, this algorithm
            *cannot* be used to compute an approximation with multiplicative
            factor, and so the ``approximation_factor`` parameter is just
            ignored here.

          - ``'dom'`` is an approximation with additive constant four. It
            computes the hyperbolicity of the vertices of a dominating set of
            the graph. This is sometimes slower than ``'cuts'`` and sometimes
            faster. Try it to know if it is interesting for you.
            The ``additive_gap`` and ``approximation_factor`` parameters cannot
            be used in combination with this method and so are ignored.

    - ``approximation_factor`` -- (default: None) When the approximation factor
      is set to some value (larger than 1.0), the function stop computations as
      soon as the ratio between the upper bound and the best found solution is
      less than the approximation factor. When the approximation factor is 1.0,
      the problem is solved optimaly. This parameter is used only when the
      chosen algorithm is ``'cuts'``.

    - ``additive_gap`` -- (default: None) When sets to a positive number, the
      function stop computations as soon as the difference between the upper
      bound and the best found solution is less than additive gap. When the gap
      is 0.0, the problem is solved optimaly. This parameter is used only when
      the chosen algorithm is ``'cuts'`` or ``'cuts+'``. The parameter must be
      ``\geq 1`` when used with ``'cuts+'``.

    - ``verbose`` -- (default: ``False``) is a boolean set to True to display
      some information during execution: new upper and lower bounds, etc.

    OUTPUT:

    This function returns the tuple ( delta, certificate, delta_UB ), where:

    - ``delta`` -- the hyperbolicity of the graph (half-integer value).

    - ``certificate`` -- is the list of the 4 vertices for which the maximum
      value has been computed, and so the hyperbolicity of the graph.

    - ``delta_UB`` -- is an upper bound for ``delta``. When ``delta ==
      delta_UB``, the returned solution is optimal. Otherwise, the approximation
      factor if ``delta_UB/delta``.

    EXAMPLES:

    Hyperbolicity of a `3\times 3` grid::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = graphs.GridGraph([3,3])
        sage: hyperbolicity(G,algorithm='cuts')
        (2, [(0, 0), (0, 2), (2, 0), (2, 2)], 2)
        sage: hyperbolicity(G,algorithm='basic')
        (2, [(0, 0), (0, 2), (2, 0), (2, 2)], 2)

    Hyperbolicity of a PetersenGraph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = graphs.PetersenGraph()
        sage: hyperbolicity(G,algorithm='cuts')
        (1/2, [0, 1, 2, 3], 1/2)
        sage: hyperbolicity(G,algorithm='basic')
        (1/2, [0, 1, 2, 3], 1/2)
        sage: hyperbolicity(G,algorithm='dom')
        (0, [0, 2, 8, 9], 1)

    Asking for an approximation::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = graphs.GridGraph([2,10])
        sage: hyperbolicity(G,algorithm='cuts', approximation_factor=1.5)
        (1, [(0, 0), (0, 9), (1, 0), (1, 9)], 3/2)
        sage: hyperbolicity(G,algorithm='cuts', approximation_factor=4)
        (1, [(0, 0), (0, 9), (1, 0), (1, 9)], 4)
        sage: hyperbolicity(G,algorithm='cuts', additive_gap=2)
        (1, [(0, 0), (0, 9), (1, 0), (1, 9)], 3)
        sage: hyperbolicity(G,algorithm='cuts+')
        (1, [(0, 0), (0, 9), (1, 0), (1, 9)], 2)
        sage: hyperbolicity(G,algorithm='cuts+', additive_gap=2)
        (1, [(0, 0), (0, 9), (1, 0), (1, 9)], 3)
        sage: hyperbolicity(G,algorithm='dom')
        (1, [(0, 1), (0, 9), (1, 0), (1, 8)], 5)

    Comparison of results::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: for i in xrange(10): # long time
        ...       G = graphs.RandomBarabasiAlbert(100,2)
        ...       d1,_,_ = hyperbolicity(G,algorithm='basic')
        ...       d2,_,_ = hyperbolicity(G,algorithm='cuts')
        ...       d3,_,_ = hyperbolicity(G,algorithm='cuts+')
        ...       l3,_,u3 = hyperbolicity(G,approximation_factor=2)
        ...       if d1!=d2 or d1<d3 or l3>d1 or u3<d1:
        ...          print "That's not good!"

    TESTS:

    Giving anything else than a Graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: hyperbolicity([])
        Traceback (most recent call last):
        ...
        ValueError: The input parameter must be a Graph.

    Giving a non connected graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = Graph([(0,1),(2,3)])
        sage: hyperbolicity(G)
        Traceback (most recent call last):
        ...
        ValueError: The input Graph must be connected.

    Giving wrong approximation factor::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = graphs.PetersenGraph()
        sage: hyperbolicity(G,algorithm='cuts', approximation_factor=0.1)
        Traceback (most recent call last):
        ...
        ValueError: The approximation factor must be >= 1.0.

    Giving negative additive gap::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = Graph()
        sage: hyperbolicity(G,algorithm='cuts', additive_gap=-1)
        Traceback (most recent call last):
        ...
        ValueError: The additive gap must be >= 0 when using the 'cuts' algorithm.

    Asking for an unknown algorithm::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = Graph()
        sage: hyperbolicity(G,algorithm='tip top')
        Traceback (most recent call last):
        ...
        ValueError: Algorithm 'tip top' not yet implemented. Please contribute.
    """
    if not isinstance(G,Graph):
        raise ValueError("The input parameter must be a Graph.")
    if not algorithm in ['basic', 'cuts', 'cuts+', 'dom']:
        raise ValueError("Algorithm '%s' not yet implemented. Please contribute." %(algorithm))
    if approximation_factor is None:
        approximation_factor = 1.0
    elif algorithm=='cuts' and (not approximation_factor in RR or approximation_factor < 1.0):
        raise ValueError("The approximation factor must be >= 1.0.")
    elif algorithm!='cuts':
        print "The approximation_factor is ignored when using the '%s' algorithm." %(algorithm)
    if additive_gap is None:
        additive_gap = 1.0 if algorithm=='cuts+' else 0.0
    elif algorithm=='cuts' or algorithm=='cuts+':
        if not additive_gap in RR:
            raise ValueError("The additive gap must be a real positive number.")
        elif algorithm=='cuts' and additive_gap < 0.0:
            raise ValueError("The additive gap must be >= 0 when using the '%s' algorithm." %(algorithm))
        elif algorithm=='cuts+' and additive_gap < 1.0:
            raise ValueError("The additive gap must be >= 1 when using the '%s' algorithm." %(algorithm))
    else:
        print "The additive_gap is ignored when using the '%s' algorithm." %(algorithm)

    # The hyperbolicity is defined on connected graphs
    if not G.is_connected():
        raise ValueError("The input Graph must be connected.")

    # The hyperbolicity of some classes of graphs is known. If it is easy and
    # fast to test that a graph belongs to one of these classes, we do it.
    if G.num_verts() <= 3:
        # The hyperbolicity of a graph with 3 vertices is 0.
        # The certificate is the set of vertices.
        return 0, G.vertices(), 0

    elif G.num_verts() == G.num_edges() + 1:
        # G is a tree
        # Any set of 4 vertices is a valid certificate
        return 0, G.vertices()[:4], 0

    elif G.is_clique():
        # Any set of 4 vertices is a valid certificate
        return 0, G.vertices()[:4], 0

    cdef unsigned short * _distances_
    cdef unsigned short ** distances
    cdef int i, j, k, N, hyp, hyp_UB, hh, hh_UB, D
    cdef dict distr = {}
    cdef list certificate = []
    cdef list certif
    cdef dict mymap, myinvmap

    # We search for the largest 2-connected component. Indeed, the hyperbolicity
    # of a graph is the maximum over its 2-connected components.
    B,C = G.blocks_and_cut_vertices()

    if verbose:
        # we compute the distribution of size of the blocks
        for V in B:
            i = len(V)
            if i in distr:
                distr[ i ] += 1
            else:
                distr[ i ] = 1
        print "Graph with %d blocks" %(len(B))
        print "Blocks size distribution:", distr

    hyp = 0
    for V in B:

        # The hyperbolicity of a graph with 3 vertices is 0, and a graph cannot
        # have hyperbolicity larger than N/2. So we consider only larger
        # 2-connected subgraphs.
        if len(V) > max( 3, 2*hyp) :

            # We build the subgraph and we relabel the vertices to ensure
            # integer vertex names in the range [0..N-1] since the
            # c_distances_all_pairs uses integer labels in the range [0..N-1].
            H, mymap = _my_subgraph(G, V, relabel=True, return_map=True)
            N = H.num_verts()

            # We test if the block belongs to a graph class with known
            # hyperbolicity and a fast test.
            if H.is_clique():
                continue

            # We compute the distances and store the results in a 2D array, and
            # the diameter
            _distances_ = c_distances_all_pairs(H)
            distances = <unsigned short **>sage_malloc(sizeof(unsigned short *)*N)
            if distances == NULL:
                sage_free(_distances_)
                raise MemoryError

            D = 0
            for 0 <= i < N:
                distances[i] = _distances_+i*N
                for 0 <= j < N:
                    if distances[i][j] > D:
                        D = distances[i][j]

            # We call the cython function for computing the hyperbolicity with
            # the required parameters.
            if algorithm == 'cuts':
                hh, certif, hh_UB = __hyperbolicity__(N, distances, D, hyp, approximation_factor, 2*additive_gap, [], verbose)
                hh_UB = min( hh_UB, D)

            elif algorithm == 'cuts+':
                # We compute the elimination ordering of simplicial vertices of H
                elim = elimination_ordering_of_simplicial_vertices(H, max(2,floor(N**(1/2.0))), verbose)
                if len(elim) > N-4:
                    # We know that this component has hyperbolicity <= 1
                    certif = H.vertices()[:4]
                    hh = __hyp__(distances,certif[0],certif[1],certif[2],certif[3])
                    hh_UB = 2
                else:
                    hh, certif, hh_UB = __hyperbolicity__(N, distances, D, hyp, 1.0, 2*additive_gap-2, elim, verbose)
                    hh_UB = min( hh_UB+2, D)

            elif algorithm == 'dom':
                # Computes a dominating set DOM of H, and computes the
                # hyperbolicity considering only vertices in DOM
                DOM = _greedy_dominating_set(H)
                elim = [u for u in H.vertex_iterator() if not u in DOM]
                # We need at least 4 vertices
                while len(DOM)<4:
                    DOM.append(elim.pop())
                hh, certif, hh_UB = __hyperbolicity__(N, distances, D, hyp, 1.0, 0.0, elim, verbose)
                hh_UB = min( hh+8, D)

            elif algorithm == 'basic':
                hh, certif = __hyperbolicity_basic_algorithm__(N, distances, verbose)
                hh_UB = hh

            # We test if the new computed value improves upon previous value.
            if hh > hyp or (hh==hyp and not certificate):
                hyp = hh
                hyp_UB = hh_UB

                # We construct the inverse mapping of the relabeling of the
                # vertices of the subgraph
                myinvmap = dict([(mymap[x],x) for x in mymap])

                # We then construct the correct certificate
                certificate = [myinvmap[u] for u in certif]

            # We now release the memory
            sage_free(distances)
            sage_free(_distances_)

    # Last, we return the computed value and the certificate
    return  ZZ(hyp)/2, sorted(certificate), ZZ(hyp_UB)/2


######################################################################
# Distribution of the hyperbolicity of 4-tuples
######################################################################

cdef dict __hyperbolicity_distribution__(int N, unsigned short ** distances):
    """
    Return the distribution of the hyperbolicity of the 4-tuples of the graph.

    The hyperbolicity of a graph has been defined by Gromov [Gromov87]_ as
    follows: Let `a, b, c, d` be vertices of the graph, let `S_1 = dist(a, b) +
    dist(b, c)`, `S_2 = dist(a, c) + dist(b, d)`, and `S_3 = dist(a, d) +
    dist(b, c)`, and let `M_1` and `M_2` be the two largest values among `S_1`,
    `S_2`, and `S_3`. We have `hyp(a, b, c, d) = |M_1 - M_2|`, and the
    hyperbolicity of the graph is the maximum over all possible 4-tuples `(a, b,
    c, d)` divided by 2.

    The computation of the hyperbolicity of each 4-tuple, and so the
    hyperbolicity distribution, takes time in `O( n^4 )`.

    We use ``unsigned long int`` on 64 bits, so ``uint64_t``, to count the
    number of 4-tuples of given hyperbolicity. So we cannot exceed `2^64-1`.
    This value should be sufficient for most users.

    INPUT:

    - ``N`` -- number of vertices of the graph (and side of the matrix)

    - ``distances`` -- matrix of distances in the graph

    OUTPUT:

    - ``hdict`` -- A dictionnary such that hdict[i] is the number of 4-tuples of
      hyperbolicity i among the considered 4-tuples.
    """
    # We initialize the table of hyperbolicity. We use an array of unsigned long
    # int instead of a dictionnary since it is much faster.
    cdef int i

    cdef uint64_t * hdistr = <uint64_t *>sage_calloc(N+1,sizeof(uint64_t))
    if hdistr == NULL:
        raise MemoryError

    # We now compute the hyperbolicity of each 4-tuple
    cdef int a, b, c, d
    for 0 <= a < N-3:
        for a < b < N-2:
            for b < c < N-1:
                for c < d < N:
                    hdistr[ __hyp__(distances, a, b, c, d) ] += 1

    # We prepare the dictionnary of hyperbolicity distribution to return
    Nchoose4 = binomial(N,4)
    cdef dict hdict = {ZZ(i)/2: (ZZ(hdistr[i])/Nchoose4) for 0 <= i <= N if hdistr[i] > 0}

    sage_free(hdistr)

    return hdict


# We use this trick since it is way faster than using the sage randint function.
cdef extern from "stdlib.h":
    long c_libc_random "random"()
    void c_libc_srandom "srandom"(unsigned int seed)

cdef dict __hyperbolicity_sampling__(int N, unsigned short ** distances, uint64_t sampling_size):
    """
    Return a sampling of the hyperbolicity distribution of the graph.

    The hyperbolicity of a graph has been defined by Gromov [Gromov87]_ as
    follows: Let `a, b, c, d` be vertices of the graph, let `S_1 = dist(a, b) +
    dist(b, c)`, `S_2 = dist(a, c) + dist(b, d)`, and `S_3 = dist(a, d) +
    dist(b, c)`, and let `M_1` and `M_2` be the two largest values among `S_1`,
    `S_2`, and `S_3`. We have `hyp(a, b, c, d) = |M_1 - M_2|`, and the
    hyperbolicity of the graph is the maximum over all possible 4-tuples `(a, b,
    c, d)` divided by 2.

    We use ``unsigned long int`` on 64 bits, so ``uint64_t``, to count the
    number of 4-tuples of given hyperbolicity. So we cannot exceed `2^64-1`.
    This value should be sufficient for most users.

    INPUT:

    - ``N`` -- number of vertices of the graph (and side of the matrix)

    - ``distances`` -- matrix of distances in the graph

    - ``sampling_size`` -- number of 4-tuples considered. Default value is 1000.

    OUTPUT:

    - ``hdict`` -- A dictionnary such that hdict[i] is the number of 4-tuples of
                hyperbolicity i among the considered 4-tuples.
    """
    cdef int i, a, b, c, d
    cdef uint64_t j

    if N < 4:
        raise ValueError("N must be at least 4")

    # We initialize the table of hyperbolicity. We use an array of unsigned long
    # int instead of a dictionnary since it is much faster.
    cdef uint64_t * hdistr = <uint64_t *>sage_calloc(N+1,sizeof(uint64_t))
    if hdistr == NULL:
        raise MemoryError

    # We now compute the hyperbolicity of each quadruple
    for 0 <= j < sampling_size:
        a = c_libc_random() % N
        b = c_libc_random() % N
        c = c_libc_random() % N
        d = c_libc_random() % N
        while a == b:
            b = c_libc_random() % N
        while a == c or b == c:
            c = c_libc_random() % N
        while a == d or b == d or c == d:
            d = c_libc_random() % N

        hdistr[ __hyp__(distances, a, b, c, d) ] += 1

    # We prepare the dictionnary of hyperbolicity distribution from sampling
    cdef dict hdict = dict( [ (ZZ(i)/2, ZZ(hdistr[i])/ZZ(sampling_size)) for 0 <= i <= N if hdistr[i] > 0 ] )

    sage_free(hdistr)

    return hdict


def hyperbolicity_distribution(G, algorithm='sampling', sampling_size=10**6):
    r"""
    Return the hyperbolicity distribution of the graph or a sampling of it.

    The hyperbolicity of a graph has been defined by Gromov [Gromov87]_ as
    follows: Let `a, b, c, d` be vertices of the graph, let `S_1 = dist(a, b) +
    dist(b, c)`, `S_2 = dist(a, c) + dist(b, d)`, and `S_3 = dist(a, d) +
    dist(b, c)`, and let `M_1` and `M_2` be the two largest values among `S_1`,
    `S_2`, and `S_3`. We have `hyp(a, b, c, d) = |M_1 - M_2|`, and the
    hyperbolicity of the graph is the maximum over all possible 4-tuples `(a, b,
    c, d)` divided by 2.

    The computation of the hyperbolicity of each 4-tuple, and so the
    hyperbolicity distribution, takes time in `O( n^4 )`.

    INPUT:

    - ``G`` -- a Graph.

    - ``algorithm`` -- (default: 'sampling') When algorithm is 'sampling', it
      returns the distribution of the hyperbolicity over a sample of
      ``sampling_size`` 4-tuples. When algorithm is 'exact', it computes the
      distribution of the hyperbolicity over all 4-tuples. Be aware that the
      computation time can be HUGE.

    - ``sampling_size`` -- (default: `10^6`) number of 4-tuples considered in
      the sampling. Used only when ``algorithm == 'sampling'``.

    OUTPUT:

    - ``hdict`` -- A dictionnary such that hdict[i] is the number of 4-tuples of
      hyperbolicity i.

    EXAMPLES:

    Exact hyperbolicity distribution of the Petersen Graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity_distribution
        sage: G = graphs.PetersenGraph()
        sage: hyperbolicity_distribution(G,algorithm='exact')
        {0: 3/7, 1/2: 4/7}

    Exact hyperbolicity distribution of a `3\times 3` grid::

        sage: from sage.graphs.hyperbolicity import hyperbolicity_distribution
        sage: G = graphs.GridGraph([3,3])
        sage: hyperbolicity_distribution(G,algorithm='exact')
        {0: 11/18, 1: 8/21, 2: 1/126}

    TESTS:

    Giving anythin else than a Graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity_distribution
        sage: hyperbolicity_distribution([])
        Traceback (most recent call last):
        ...
        ValueError: The input parameter must be a Graph.

    Giving a non connected graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity_distribution
        sage: G = Graph([(0,1),(2,3)])
        sage: hyperbolicity_distribution(G)
        Traceback (most recent call last):
        ...
        ValueError: The input Graph must be connected.
    """
    if not isinstance(G,Graph):
        raise ValueError("The input parameter must be a Graph.")
    # The hyperbolicity is defined on connected graphs
    if not G.is_connected():
        raise ValueError("The input Graph must be connected.")

    # The hyperbolicity distribution of some classes of graphs is known. If it
    # is easy and fast to test that a graph belongs to one of these classes, we
    # do it.
    if (G.num_verts()==G.num_edges()+1) or G.is_clique():
        return {0: sampling_size if algorithm=='sampling' else binomial(G.num_verts(),4)}

    cdef int N = G.num_verts()
    cdef int i, j
    cdef unsigned short ** distances
    cdef unsigned short * _distances_
    cdef dict hdict

    # We compute the all pairs shortest path and store the result in a 2D array
    # for faster access.
    H = G.relabel( inplace = False )
    _distances_ = c_distances_all_pairs(H)
    distances = <unsigned short **>sage_malloc(sizeof(unsigned short *)*N)
    if distances == NULL:
        sage_free(_distances_)
        raise MemoryError

    for 0 <= i < N:
        distances[i] = _distances_+i*N

    if algorithm == 'exact':
        hdict = __hyperbolicity_distribution__(N, distances)
    elif algorithm == 'sampling':
        hdict = __hyperbolicity_sampling__(N, distances, sampling_size)
    else:
        raise ValueError("Algorithm '%s' not yet implemented. Please contribute." %(algorithm))

    # We release memory
    sage_free(distances)
    sage_free(_distances_)

    return hdict
