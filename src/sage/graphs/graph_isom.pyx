"""
N.I.C.E. - Nice (as in open source) Isomorphism Check Engine

Automorphism group computation and isomorphism checking for graphs.

This is an open source implementation of Brendan McKay's algorithm for graph
automorphism and isomorphism. McKay released a C version of his algorithm,
named nauty (No AUTomorphisms, Yes?) under a license that is not GPL
compatible. Although the program is open source, reading the source disallows
anyone from recreating anything similar and releasing it under the GPL. Also,
many people have complained that the code is difficult to understand. The
first main goal of NICE was to produce a genuinely open graph isomorphism
program, which has been accomplished. The second goal is for this code to be
understandable, so that computed results can be trusted and further derived
work will be possible.

To determine the isomorphism type of a graph, it is convenient to define a
canonical label for each isomorphism class- essentially an equivalence class
representative. Loosely (albeit incorrectly), the canonical label is defined
by enumerating all labeled graphs, then picking the maximal one in each
isomorphism class. The NICE algorithm is essentially a backtrack search. It
searches through the rooted tree of partition nests (where each partition is
equitable) for implicit and explicit automorphisms, and uses this information
to eliminate large parts of the tree from further searching. Since the leaves
of the search tree are all discrete ordered partitions, each one of these
corresponds to an ordering of the vertices of the graph, i.e. another member
of the isomorphism class. Once the algorithm has finished searching the tree,
it will know which leaf corresponds to the canonical label. In the process,
generators for the automorphism group are also produced.

AUTHORS:
    Robert L. Miller -- (2007-03-20) initial version
    Tom Boothby -- (2007-03-20) help with indicator function
    Robert L. Miller -- (2007-04-07--30) optimizations

REFERENCE:
    [1] McKay, Brendan D. Practical Graph Isomorphism. Congressus Numerantium,
        Vol. 30 (1981), pp. 45-87.

NOTE:
    Often we assume that G is a graph on vertices {0,1,...,n-1}.
"""

#*****************************************************************************
#      Copyright (C) 2006 - 2007 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include '../ext/python_mem.pxi'
from sage.graphs.graph import Graph, DiGraph
from sage.misc.misc import prod, cputime
from copy import copy
from sage.rings.infinity import Infinity
from sage.graphs.graph import enum
from sage.sets.set import Set
from sage.rings.integer import Integer

def finer(Pi1, Pi2):
    """
    Returns True if Pi1 is (weakly) finer than Pi2: i.e. each cell of the
    partition Pi1 is contained in some cell of Pi2.

    EXAMPLES:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import finer
        sage: finer( [[0,1,2],[3,4,5]] , [[0,1],[2,3],[4,5]] )
        False
        sage: finer( [[0,1],[2,3],[4,5],[6,7]] , [[0,1,2,3],[4,5,6,7]] )
        True
        sage: finer( [[0,1],[2,3],[4,5],[6,7]] , [[0,1],[2,3],[4,5],[6,7]] )
        True
        sage: finer([[1,2],[0]] , [[1],[2],[0]])
        False
    """
    for p in Pi1:
        cell = None
        i = 0
        while cell is None:
            if p[0] in Pi2[i]:
                cell = Pi2[i]
            i += 1
        for i in range(1,len(p)):
            if not p[i] in cell:
                return False
    return True

cdef class OrderedPartition:

    def __new__(self, data):
        cdef int i, j, k
        cdef OrderedPartition _data
        if isinstance(data, (list, tuple)):
            self.length = len(data)
            self.data = <int **> PyMem_Malloc( self.length * sizeof(int *) )
            self.sizes = <int *> PyMem_Malloc( self.length * sizeof(int) )
            try:
                for i from 0 <= i < self.length:
                    self.sizes[i] = len(data[i])
                    self.data[i] = <int *> PyMem_Malloc( self.sizes[i] * sizeof(int) )
                    for j from 0 <= j < self.sizes[i]:
                        self.data[i][j] = data[i][j]
            except:
                raise TypeError('Error with input data.')
        elif isinstance(data, OrderedPartition):
            _data = data
            self.length = _data.length
            self.data = <int **> PyMem_Malloc( self.length * sizeof(int *) )
            self.sizes = <int *> PyMem_Malloc( self.length * sizeof(int) )
            for i from 0 <= i < self.length:
                self.sizes[i] = _data.sizes[i]

                #self.sizes[i] = data.size(i)
                #self.data[i] = <int *> PyMem_Malloc( self.sizes[i] * sizeof(int) )
                #for j from 0 <= j < self.sizes[i]:
                #    self.data[i][j] = data[i][j]

    def __dealloc__(self):
        return
        cdef int i, j
        for i from 0 <= i < self.length:
            PyMem_Free(self.data[i])
        PyMem_Free(self.sizes)
        PyMem_Free(self.data)

    def __repr__(self):
        s = "("
        cdef int i, j
        for i from 0 <= i < self.length:
            s += "{"
            for j from 0 <= j < self.sizes[i] - 1:
                s += str(self.data[i][j]) + ","
            if self.sizes[i] > 0:
                s += str(self.data[i][self.sizes[i] - 1])
            s += "},"
        s = s[:-1] + ")"
        return s

    def __getitem__(self, n):
        cdef int i
        return [ self.data[n][i] for i from 0 <= i < self.sizes[n] ]

    def is_discrete(self):
        cdef int i
        for i from 0 <= i < self.length:
            if self.sizes[i] > 1:
                return False
        return True

    cdef public int size(self, int n):
        return self.sizes[n]

def vee(alpha, beta):
    """
    Return the finest partition coarser than alpha and beta.

    EXAMPLES:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import vee
        sage: vee([[0]],[[0]])
        [[0]]
        sage: vee([[0,1]], [[0,1]])
        [[0, 1]]
        sage: vee([[1,2],[0]], [[1],[2],[0]])
        [[1, 2], [0]]
        sage: vee([[0,1,2],[3,4],[5,6,7,8],[9]], [[0],[1],[2,3],[4,5],[6],[7],[8,9]])
        [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]
    """
    p = []
    n = 3
    alpha = copy(alpha)
    beta = copy(beta)
    cell = None
    i = 0
    while (not finer(alpha, beta)) and 0 < n:
        for v in beta[i]:
            for a in alpha:
                if v in a:
                    for w in a:
                        if w not in beta[i]:
                            for j in range(len(beta)):
                                if w in beta[j]:
                                    beta[i] += beta[j]
                                    beta.pop(j)
                                    if i > j: i -= 1
                                    break
                    break
        #done with beta[i]
        i += 1
        n-=1
    return beta

def min_cell_reps(Pi):
    """
    Returns the minimum cell representatives of the partition Pi: one element
    from each cell, minimal in each cell.

    EXAMPLES:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import min_cell_reps
        sage: min_cell_reps( [[0,1],[2,3],[4,5],[6,7]] )
        [0, 2, 4, 6]
        sage: min_cell_reps( [[0,1,2,3],[4,5,6],[7]] )
        [0, 4, 7]
    """
    l = []
    for p in Pi:
        l.append( min(p) )
    return l

def fix(Pi):
    """
    Returns a list of all the elements which live in trivial cells: if
    the partition Pi represents the set of orbits of an action, then
    fix(Pi) is the subset consisting of elements fixed by the action.

    EXAMPLES:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import fix
        sage: fix( [[0],[1,2,3,4],[5,6],[7],[8]] )
        [0, 7, 8]
    """
    l = []
    for cell in Pi:
        if len(cell) == 1:
            l.append(cell[0])
    return l

def orbit_partition(gamma, list_perm=False):
    r"""
    Assuming that G is a graph on vertices {0,1,...,n-1}, and gamma is an
    element of SymmetricGroup(n), returns the partition of the vertex set
    determined by the orbits of gamma, considered as action on the set
    {1,2,...,n} where we take 0 = n. In other words, returns the partition
    determined by a cyclic representation of gamma.

    INPUT:
        list_perm -- if True, assumes gamma is a list representing the map
    i \mapsto gamma[i].

    EXAMPLES:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import orbit_partition
        sage: G = graphs.PetersenGraph()
        sage: S = SymmetricGroup(10)
        sage: gamma = S('(10,1,2,3,4)(5,6,7)(8,9)')
        sage: orbit_partition(gamma)
        [[1, 2, 3, 4, 0], [5, 6, 7], [8, 9]]
        sage: gamma = S('(10,5)(1,6)(2,7)(3,8)(4,9)')
        sage: orbit_partition(gamma)
        [[1, 6], [2, 7], [3, 8], [4, 9], [5, 0]]
    """
    if list_perm:
        n = len(gamma)
        seen = [1] + [0]*(n-1)
        i = 0
        p = 0
        partition = [[0]]
        while sum(seen) < n:
            if gamma[i] != partition[p][0]:
                partition[p].append(gamma[i])
                i = gamma[i]
                seen[i] = 1
            else:
                i = min([j for j in range(n) if seen[j] == 0])
                partition.append([i])
                p += 1
                seen[i] = 1
        return partition
    else:
        n = len(gamma.list())
        l = []
        for i in range(1,n+1):
            orb = gamma.orbit(i)
            if orb not in l: l.append(orb)
        for i in l:
            for j in range(len(i)):
                if i[j] == n:
                    i[j] = 0
        return l

def sat225(Pi, n):
    """
    Returns true if Pi satisfies the conditions of Lemma 2.25 in [1].
    """
    m = 0
    for p in Pi:
        if len(p) > 1:
            m += 1
    # Pi has m nontrivial cells
    if n <= len(Pi) + 4: return True
    elif n == len(Pi) + m: return True
    elif n == len(Pi) + m + 1: return True
    else: return False

cdef degree(int** G, v, W):
    """
    Returns the number of edges from vertices in W to v, i.e. the degree of v
    with respect to W. W is usually a cell in a partition, but needs only be
    a subset of the vertex set.
    """
    cdef int i = 0
    for u in W:
        if G[u][v]:
            i += 1
    return i

cdef degree_inv(int** G, v, W):
    """
    Returns the number of edges from v to vertices in W, i.e. the out-degree of v
    with respect to W. W is usually a cell in a partition, but needs only be
    a subset of the vertex set.
    """
    cdef int i = 0
    for u in W:
        if G[v][u]:
            i += 1
    return i

def is_discrete(Pi):
    """
    Returns true iff every cell in the partition Pi is of size 1.

    EXAMPLES:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import is_discrete
        sage: is_discrete( [ [0,1,2,3], [4], [5], [6] ] )
        False
        sage: is_discrete( [ [0], [1], [2], [3], [4], [5], [6] ] )
        True
    """
    cdef int i
    for i from 0 <= i < len(Pi):
        if len(Pi[i]) != 1:
            return False
    return True

#def is_equitable(G, Pi):
#    """
#    A partition Pi of the vertex set of G is said to be equitable if for any
#    two cells V1 and V2 of Pi, and for any two vertices v1 and v2 in V1, we
#    have degree(v_1, V_2) == degree(v_2, V_2).
#
#    Educational only: not used in the main algorithm. If the partition is
#    equitable, returns True. If not, then returns False, along with the
#    counterexample: False, v1, v2, V1, V2.
#
#    EXAMPLES:
#        sage: import sage.graphs.graph_isom
#        sage: from sage.graphs.graph_isom import is_equitable
#        sage: D  = graphs.DodecahedralGraph()
#        sage: Pi1 = [[0],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]]
#        sage: Pi2 = [[0],[1,10,19],[11,9,18,8,3,2],[12,13,7,17,4,6],[5,14,16],[15]]
#        sage: is_equitable(D, Pi1)
#        (False, 2, 1, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], [0])
#        sage: is_equitable(D, Pi2)
#        True
#
#    To see what is going on:
#        sage.: D.show(partition=Pi1)
#        sage.: D.show(partition=Pi2)
#    """
#    for i in range(len(Pi)):
#        for j in range(i):
#            for v_1 in range(len(Pi[i])):
#                for v_2 in range(v_1):
#                    if degree(G, Pi[i][v_1], Pi[j]) != degree(G, Pi[i][v_2], Pi[j]):
#                        return False, Pi[i][v_1], Pi[i][v_2], Pi[i], Pi[j]
#    return True

def replace_in(Pi, k, L):
    """
    Replaces cell k of the partition Pi with the partition of that cell, L.

    EXAMPLE:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import replace_in
        sage: replace_in([[0,1],[2,3],[4,5,6,7]], 2, [[4],[5],[6,7]] )
        [[0, 1], [2, 3], [4], [5], [6, 7]]
    """
    PiNew = []
    for j in Pi[:k]: PiNew.append(j)
    for j in L: PiNew.append(j)
    for j in Pi[k+1:]: PiNew.append(j)
    return PiNew

cdef sort_by_degree(int** G, A, B, dig):
    """
    Assuming A and B are subsets of the vertex set of G, returns an ordered
    partition of A such that degree(x,B) < degree(y,B) iff x occurs in an
    earlier cell than y.
    """
    ddict = {}
    for a in A:
        dd = degree(G, a, B)
        try:
            ddict[dd].append(a)
        except:
            ddict[dd] = [a]
    if dig:
        ee = []
        for part in ddict.values():
            edict = {}
            for a in part:
                dd = degree_inv(G, a, B)
                try:
                    edict[dd].append(a)
                except:
                    edict[dd] = [a]
            ee += edict.values()
        return ee
    else:
        return ddict.values()

cdef refine(int** G, Pi, alpha, dig):
    """
    The key refinement procedure. Given a graph G, a partition Pi of the
    vertex set, and a collection alpha of disjoint subsets of the vertex set,
    returns a refinement R of the partition Pi, so that each cell C of R has
    the same degree to each component of alpha, i.e. degree(c, W) is the same
    for each W in alpha and each c in a fixed C. The order of the refined
    partition also matters: each time a cell is split, the subcell of maximal
    size takes its place, and the rest go to the end of the list alpha.

    It is a theorem (2.6 in [1]) that refine(G, Pi, Pi) is always the unique
    coarsest equitable partition that is finer than Pi. Further (2.7 in [1]),
    if Pi_prime is an equitable partition coarser than Pi, and if alpha is
    chosen from cells of Pi such that for any W in Pi_prime, X is a subset of
    W for at most one X in Pi - alpha, then refine(G, Pi, Pi) is always the
    unique coarsest equitable partition that is finer than Pi.
    """
    cdef int i, j, q
    alpha = copy(alpha) # we don't want to change the user's alpha
    M = len(alpha)
    PiT = copy(Pi)
    m = 0
    while (not is_discrete(PiT)) and m < M:
        W = alpha[m]
        m += 1
        k = 0
        r = len(PiT)
        while k < r:
            X = sort_by_degree(G, PiT[k], W, dig)
            s = len(X)
            if s != 1:
                t = 0
                L = [ len(X[q]) for q from 0 <= q < s ]
                t = L.index( max( L ) ) # relies on "index" returning smallest index
                for j from m <= j < M:
                    if Set(PiT[k]) == Set(alpha[j]):
                        alpha[j] = X[t]
                        break # you've found the one, so stop constructing Sets
                for i from 0 <= i < t:
                    alpha.append(X[i])
                for i from t+1 <= i < s:
                    alpha.append(X[i])
                M += s - 1
                PiT = replace_in(PiT, k, X)
            k += 1
    return PiT

def which(Pi, v):
    """
    Returns the index i such that v is in Pi[i].
    """
    for i in range(len(Pi)):
        if v in Pi[i]:
            return i

def comp(Pi, v):
    """
    Refines the partition Pi by replacing the cell containing v with a cell
    containing only v, followed by the rest of the cell.
    """
    i = which(Pi, v)
    if len(Pi[i]) == 1:
        return Pi
    else:
        L = []
        for vv in Pi[i]:
            if not vv == v:
                L.append(vv)
        return replace_in(Pi, i, [[v], L])

cdef perp(int** G, Pi, v, dig):
    """
    Refines the partition Pi by cutting out a vertex, then using refine,
    comparing against only that vertex.
    """
    return refine(G, comp(Pi, v), [[v]], dig)

cdef partition_nest(int** G, Pi, V, dig):
    """
    Given a sequence of vertices V = (v_1,...,v_{m-1}) of a graph G, and a
    partition Pi, the partition nest derived from G, Pi, and V is defined to
    be the sequence of partitions (Pi_1,...,Pi_m) where:

    Pi_1 := refine(G, Pi, Pi)

    Pi_k := perp(G, Pi_{k-1}, v_{k-1})
    for 2 <= k <= m.

    C.f. 2.9 in [1].
    """
    L = [refine(G, Pi, Pi, dig)]
    for i in range(len(V)):
        L.append(perp(G, L[i], V[i], dig))
    return L

def first_smallest_non_trivial(Pi):
    """
    As the name suggests, returns the first smallest nontrivial cell of the
    ordered partition Pi.
    """
    l = []
    for p in Pi:
        if len(p) != 1:
            l.append(len(p))
    m = min(l)
    for i in range(len(Pi)):
        if len(Pi[i]) == m:
            return Pi[i]

cdef indicator(int** G, Pi, V, k, n):
    """
    Takes a labelled graph, an ordered partition, and a partition nest, and
    outputs an integer, which is the same under any fixed permutation.

    AUTHORS:
        Tom Boothby -- sum then product method
        Robert Miller -- vertex degree check
    """
    V = V[:k]
    LL = []
    for partition in V:
        a = len(partition)
        for k in range(a):
            LL.append( len(partition[k])*(1 + \
                sum(  [  degree(G, partition[k][0], partition[i]) for i in range(a)  ]  ) ) )
    return prod(LL)

def get_permutation(eta, nu, list_perm=False):
    r"""
    Given two terminal nodes of the search tree, eta and nu, each last
    partition is discrete, and the order of the partition determines a
    permutation gamma such that gamma(eta) = nu. Returns the partition gamma.

    INPUT:
        list_perm -- if True, returns a list L representing the map i \mapsto
    L[i].

    EXAMPLE:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import get_permutation

    The following is an example from searching the nest tree of the dodecahedron with unit partition
        sage: eta = [ [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]], [[0], [5, 14, 16], [15], [4, 6, 7, 12, 13, 17], [2, 3, 8, 9, 11, 18], [1, 10, 19]], [[0], [5], [14, 16], [15], [12, 13], [7, 17], [4, 6], [9, 11], [8, 18], [2, 3], [10], [1, 19]], [[0], [5], [14], [16], [15], [12], [13], [17], [7], [4], [6], [11], [9], [18], [8], [3], [2], [10], [19], [1]] ]
        sage: nu = [ [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]], [[1], [12, 15, 17], [16], [4, 5, 11, 13, 14, 18], [3, 6, 7, 9, 10, 19], [0, 2, 8]], [[1], [12], [15, 17], [16], [4, 5], [14, 18], [11, 13], [3, 6], [7, 19], [9, 10], [2], [0, 8]], [[1], [12], [15], [17], [16], [4], [5], [18], [14], [11], [13], [3], [6], [19], [7], [10], [9], [2], [0], [8]] ]
        sage: get_permutation(eta, nu)
        (1,8,7,14,15,16,17,18,19,20)(12,4,11,3,10,2,9,6,13,5)
        sage: get_permutation(eta, nu, list_perm=True)
        [1, 8, 9, 10, 11, 12, 13, 14, 7, 6, 2, 3, 4, 5, 15, 16, 17, 18, 19, 0]
    """
    a = nu[-1]
    b = eta[-1]
    n = len(b)
    if list_perm:
        gamma = []
        i = 0
        while len(gamma) < n:
            for j in range(n):
                if b[j][0] == i:
                    gamma.append(a[j][0])
                    i += 1
                    break
        return gamma
    else:
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        S = SymmetricGroup(n)
        gamma = []
        for i in range(len(b)):
            if b[i][0] != a[i][0]:
                gamma.append([b[i][0],a[i][0]])
        i = 0
        while i < len(gamma):
            if gamma[i][0] == gamma[i][-1]:
                i += 1
            else:
                for j in range(i+1,len(gamma)):
                    if gamma[i][-1] == gamma[j][0]:
                        gamma[i] = gamma[i] + gamma[j][1:]
                        gamma.pop(j)
                        break
        for i in range(len(gamma)):
            gamma[i] = gamma[i][1:]
        if len(gamma) == 0:
            gamma = S('()')
        else:
            gamma = S(str(gamma)[1:-1].replace('[','(').replace(']',')').\
                   replace(' ', '').replace('(0','('+str(n)).replace(',0',','+str(n)))
        return gamma

def term_pnest_graph(G, nu, enumer=False):
    """
    BDM's G(nu): returns the graph G, relabeled in the order found in
    nu[m], where m is the first index corresponding to a discrete partition.
    Assumes nu is a terminal partition nest in T(G, Pi).
    """
    n = G.order()
    d = {}
    m = 0
    while not is_discrete(nu[m]):
        m += 1
    for i in range(n):
        d[nu[m][i][0]] = i
    if enumer:
        # we know that the vertex set is {0,...,n-1}...
        numbr = 0
        if isinstance(G, Graph):
            for i,j,l in G.edge_iterator():
                numbr += 1<<((n-(d[i]+1))*n + n-(d[j]+1))
                numbr += 1<<((n-(d[j]+1))*n + n-(d[i]+1))
        elif isinstance(G, DiGraph):
            for i,j,l in G.arc_iterator():
                numbr += 1<<((n-(d[i]+1))*n + n-(d[j]+1))
        return numbr
    else:
        ord = nu[m]
        H = G.copy()
        H.relabel(d)
        return H

def search_tree(G, Pi, lab=True, dig=False, dict=False, proof=False, verbosity=0):
    """
    Assumes that the vertex set of G is {0,1,...,n-1} for some n.

    Note that this conflicts with the SymmetricGroup we are using to represent
    automorphisms. The solution is to let the group act on the set
    {1,2,...,n}, under the assumption n = 0.

    INPUT:
        lab--       if True, return the canonical label in addition to the auto-
    morphism group.
        dig--       if True, does not use Lemma 2.25 in [1], and the algorithm is
    valid for digraphs and graphs with loops.
        dict--      if True, explain which vertices are which elements of the set
    {1,2,...,n} in the representation of the automorphism group.
        proof--     if True, return the automorphism between G and its canonical
    label. Forces lab=True.
        verbosity-- 0 - print nothing
                    1 - display state trace
                    2 - with timings
                    3 -

    STATE DIAGRAM:
        sage: SD = DiGraph( { 1:[18,2], 2:[5,3], 3:[4,6], 4:[7,2], 5:[4], 6:[13,12], 7:[18,8,10], 8:[6,9,10], 9:[6], 10:[11,13], 11:[12], 12:[13], 13:[17,14], 14:[16,15], 15:[2], 16:[13], 17:[15,13], 18:[13] } )
        sage: SD.set_arc_label(1, 18, 'discrete')
        sage: SD.set_arc_label(4, 7, 'discrete')
        sage: SD.set_arc_label(2, 5, 'h = 0')
        sage: SD.set_arc_label(7, 18, 'h = 0')
        sage: SD.set_arc_label(7, 10, 'aut')
        sage: SD.set_arc_label(8, 10, 'aut')
        sage: SD.set_arc_label(8, 9, 'label')
        sage: SD.set_arc_label(8, 6, 'no label')
        sage: SD.set_arc_label(13, 17, 'k > h')
        sage: SD.set_arc_label(13, 14, 'k = h')
        sage: SD.set_arc_label(17, 15, 'v_k finite')
        sage: SD.set_arc_label(14, 15, 'v_k m.c.r.')
        sage: posn = {1:[ 3,-3],  2:[0,2],  3:[0, 13],  4:[3,9],  5:[3,3],  6:[16, 13], 7:[6,1],  8:[6,6],  9:[6,11], 10:[9,1], 11:[10,6], 12:[13,6], 13:[16,2], 14:[10,-6], 15:[0,-10], 16:[14,-6], 17:[16,-10], 18:[6,-4]}
        sage: SD.plot(pos=posn, vertex_size=400, vertex_colors={'#FFFFFF':range(1,19)}, edge_labels=True).save('search_tree.png')

    EXAMPLES:
        sage: from sage.groups.perm_gps.permgroup import PermutationGroup
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import search_tree, perm_group_elt
        sage: from sage.graphs.graph import enum

        sage: G = graphs.DodecahedralGraph()
        sage: Pi=[range(20)]
        sage: a,b = search_tree(G, Pi)
        sage: print a, enum(b)
        [[0, 19, 3, 2, 6, 5, 4, 17, 18, 11, 10, 9, 13, 12, 16, 15, 14, 7, 8, 1], [0, 1, 8, 9, 13, 14, 7, 6, 2, 3, 19, 18, 17, 4, 5, 15, 16, 12, 11, 10], [0, 19, 18, 11, 12, 16, 17, 4, 3, 2, 1, 8, 7, 6, 5, 15, 14, 13, 9, 10], [1, 8, 9, 10, 11, 12, 13, 14, 7, 6, 2, 3, 4, 5, 15, 16, 17, 18, 19, 0]] 17318942212009113839976787462421724338461987195898671092180383421848885858584973127639899792828728124797968735273000
        sage: c = search_tree(G, Pi, lab=False)
        sage: print c
        [[0, 19, 3, 2, 6, 5, 4, 17, 18, 11, 10, 9, 13, 12, 16, 15, 14, 7, 8, 1], [0, 1, 8, 9, 13, 14, 7, 6, 2, 3, 19, 18, 17, 4, 5, 15, 16, 12, 11, 10], [0, 19, 18, 11, 12, 16, 17, 4, 3, 2, 1, 8, 7, 6, 5, 15, 14, 13, 9, 10], [1, 8, 9, 10, 11, 12, 13, 14, 7, 6, 2, 3, 4, 5, 15, 16, 17, 18, 19, 0]]
        sage: DodecAut = PermutationGroup([perm_group_elt(aa) for aa in a])
        sage: DodecAut.character_table()
        [                     1                      1                      1                      1                      1                      1                      1                      1                      1                      1]
        [                     1                     -1                      1                      1                     -1                      1                     -1                      1                     -1                     -1]
        [                     3                     -1                      0                     -1  zeta5^3 + zeta5^2 + 1     -zeta5^3 - zeta5^2                      0  zeta5^3 + zeta5^2 + 1     -zeta5^3 - zeta5^2                      3]
        [                     3                     -1                      0                     -1     -zeta5^3 - zeta5^2  zeta5^3 + zeta5^2 + 1                      0     -zeta5^3 - zeta5^2  zeta5^3 + zeta5^2 + 1                      3]
        [                     3                      1                      0                     -1 -zeta5^3 - zeta5^2 - 1     -zeta5^3 - zeta5^2                      0  zeta5^3 + zeta5^2 + 1      zeta5^3 + zeta5^2                     -3]
        [                     3                      1                      0                     -1      zeta5^3 + zeta5^2  zeta5^3 + zeta5^2 + 1                      0     -zeta5^3 - zeta5^2 -zeta5^3 - zeta5^2 - 1                     -3]
        [                     4                      0                      1                      0                     -1                     -1                      1                     -1                     -1                      4]
        [                     4                      0                      1                      0                      1                     -1                     -1                     -1                      1                     -4]
        [                     5                      1                     -1                      1                      0                      0                     -1                      0                      0                      5]
        [                     5                     -1                     -1                      1                      0                      0                      1                      0                      0                     -5]
        sage: DodecAut2 = PermutationGroup([perm_group_elt(cc) for cc in c])
        sage: DodecAut2.character_table()
        [                     1                      1                      1                      1                      1                      1                      1                      1                      1                      1]
        [                     1                     -1                      1                      1                     -1                      1                     -1                      1                     -1                     -1]
        [                     3                     -1                      0                     -1  zeta5^3 + zeta5^2 + 1     -zeta5^3 - zeta5^2                      0  zeta5^3 + zeta5^2 + 1     -zeta5^3 - zeta5^2                      3]
        [                     3                     -1                      0                     -1     -zeta5^3 - zeta5^2  zeta5^3 + zeta5^2 + 1                      0     -zeta5^3 - zeta5^2  zeta5^3 + zeta5^2 + 1                      3]
        [                     3                      1                      0                     -1 -zeta5^3 - zeta5^2 - 1     -zeta5^3 - zeta5^2                      0  zeta5^3 + zeta5^2 + 1      zeta5^3 + zeta5^2                     -3]
        [                     3                      1                      0                     -1      zeta5^3 + zeta5^2  zeta5^3 + zeta5^2 + 1                      0     -zeta5^3 - zeta5^2 -zeta5^3 - zeta5^2 - 1                     -3]
        [                     4                      0                      1                      0                     -1                     -1                      1                     -1                     -1                      4]
        [                     4                      0                      1                      0                      1                     -1                     -1                     -1                      1                     -4]
        [                     5                      1                     -1                      1                      0                      0                     -1                      0                      0                      5]
        [                     5                     -1                     -1                      1                      0                      0                      1                      0                      0                     -5]

        sage: G = graphs.PetersenGraph()
        sage: Pi=[range(10)]
        sage: a,b = search_tree(G, Pi)
        sage: print a, enum(b)
        [[0, 1, 2, 7, 5, 4, 6, 3, 9, 8], [0, 1, 6, 8, 5, 4, 2, 9, 3, 7], [0, 4, 3, 8, 5, 1, 9, 2, 6, 7], [1, 0, 4, 9, 6, 2, 5, 3, 7, 8]] 8716441511243809436161868448
        sage: c = search_tree(G, Pi, lab=False)
        sage: PAut = PermutationGroup([perm_group_elt(aa) for aa in a])
        sage: PAut.character_table()
        [ 1  1  1  1  1  1  1]
        [ 1 -1  1 -1  1 -1  1]
        [ 4 -2  0  1  1  0 -1]
        [ 4  2  0 -1  1  0 -1]
        [ 5  1  1  1 -1 -1  0]
        [ 5 -1  1 -1 -1  1  0]
        [ 6  0 -2  0  0  0  1]
        sage: PAut = PermutationGroup([perm_group_elt(cc) for cc in c])
        sage: PAut.character_table()
        [ 1  1  1  1  1  1  1]
        [ 1 -1  1 -1  1 -1  1]
        [ 4 -2  0  1  1  0 -1]
        [ 4  2  0 -1  1  0 -1]
        [ 5  1  1  1 -1 -1  0]
        [ 5 -1  1 -1 -1  1  0]
        [ 6  0 -2  0  0  0  1]

        sage: G = graphs.CubeGraph(3)
        sage: Pi = []
        sage: for i in range(8):
        ...    b = Integer(i).binary()
        ...    Pi.append(b.zfill(3))
        ...
        sage: Pi = [Pi]
        sage: a,b = search_tree(G, Pi)
        sage: print a, enum(b)
        [[0, 3, 2, 1, 6, 5, 4, 7], [0, 1, 4, 5, 2, 3, 6, 7], [0, 3, 6, 5, 2, 1, 4, 7], [1, 0, 3, 2, 5, 4, 7, 6]] 520239721777506480
        sage: c = search_tree(G, Pi, lab=False)

        sage: PermutationGroup([perm_group_elt(aa) for aa in a]).order()
        48
        sage: PermutationGroup([perm_group_elt(cc) for cc in c]).order()
        48
        sage: DodecAut.order()
        120
        sage: PAut.order()
        120

        sage: D = graphs.DodecahedralGraph()
        sage: a,b,c = search_tree(D, [range(20)], proof=True)
        sage: from sage.plot.plot import GraphicsArray
        sage: import networkx
        sage: position_D = networkx.spring_layout(D._nxg)
        sage: position_b = {}
        sage: for vert in position_D:
        ...    position_b[c[vert]] = position_D[vert]
        sage.: GraphicsArray([D.plot(pos=position_D), b.plot(pos=position_b)]).show()
        sage: c
        {0: 0, 1: 19, 2: 16, 3: 15, 4: 9, 5: 1, 6: 10, 7: 8, 8: 14, 9: 12, 10: 17, 11: 11, 12: 5, 13: 6, 14: 2, 15: 4, 16: 3, 17: 7, 18: 13, 19: 18}

    BENCHMARKS:
    The following examples are given to check modifications to the algorithm
    for optimization-- use sage -t -long to check all the cases.

        sage: G = Graph({0:[]})
        sage: Pi = [[0]]
        sage: a,b = search_tree(G, Pi)
        sage: print a, enum(b)
        [] 0
        sage: a,b = search_tree(G, Pi, dig=True)
        sage: print a, enum(b)
        [] 0
        sage: search_tree(G, Pi, lab=False)
        []

        sage: from sage.graphs.graph_isom import all_labeled_graphs, all_ordered_partitions

        sage: graph2 = all_labeled_graphs(2)
        sage: part2 = all_ordered_partitions(range(2))
        sage: for G in graph2:
        ...    for Pi in part2:
        ...        a,b = search_tree(G, Pi)
        ...        c,d = search_tree(G, Pi, dig=True)
        ...        e = search_tree(G, Pi, lab=False)
        ...        a = str(a); b = str(enum(b)); c = str(c); d = str(enum(d)); e = str(e)
        ...        print a.ljust(15), b.ljust(5), c.ljust(15), d.ljust(5), e.ljust(15)
        []              0     []              0     []
        []              0     []              0     []
        [[1, 0]]        0     [[1, 0]]        0     [[1, 0]]
        [[1, 0]]        0     [[1, 0]]        0     [[1, 0]]
        []              6     []              6     []
        []              6     []              6     []
        [[1, 0]]        6     [[1, 0]]        6     [[1, 0]]
        [[1, 0]]        6     [[1, 0]]        6     [[1, 0]]

        sage: graph3 = all_labeled_graphs(3)
        sage: part3 = all_ordered_partitions(range(3))
        sage: for G in graph3:               # long time (~30 secs)
        ...    for Pi in part3:
        ...        a,b = search_tree(G, Pi)
        ...        c,d = search_tree(G, Pi, dig=True)
        ...        e = search_tree(G, Pi, lab=False)
        ...        a = str(a); b = str(enum(b)); c = str(c); d = str(enum(d)); e = str(e)
        ...        print a.ljust(15), b.ljust(5), c.ljust(15), d.ljust(5), e.ljust(15)
        []              0     []              0     []
        []              0     []              0     []
        [[0, 2, 1]]     0     [[0, 2, 1]]     0     [[0, 2, 1]]
        [[0, 2, 1]]     0     [[0, 2, 1]]     0     [[0, 2, 1]]
        []              0     []              0     []
        []              0     []              0     []
        [[2, 1, 0]]     0     [[2, 1, 0]]     0     [[2, 1, 0]]
        [[2, 1, 0]]     0     [[2, 1, 0]]     0     [[2, 1, 0]]
        []              0     []              0     []
        []              0     []              0     []
        [[1, 0, 2]]     0     [[1, 0, 2]]     0     [[1, 0, 2]]
        [[1, 0, 2]]     0     [[1, 0, 2]]     0     [[1, 0, 2]]
        [[1, 0, 2]]     0     [[1, 0, 2]]     0     [[1, 0, 2]]
        [[2, 1, 0]]     0     [[2, 1, 0]]     0     [[2, 1, 0]]
        [[1, 0, 2]]     0     [[1, 0, 2]]     0     [[1, 0, 2]]
        [[0, 2, 1]]     0     [[0, 2, 1]]     0     [[0, 2, 1]]
        [[2, 1, 0]]     0     [[2, 1, 0]]     0     [[2, 1, 0]]
        [[0, 2, 1]]     0     [[0, 2, 1]]     0     [[0, 2, 1]]
        [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]] 0     [[0, 2, 1], [1, 0, 2]]
        []              10    []              10    []
        []              10    []              10    []
        [[0, 2, 1]]     10    [[0, 2, 1]]     10    [[0, 2, 1]]
        [[0, 2, 1]]     10    [[0, 2, 1]]     10    [[0, 2, 1]]
        []              68    []              68    []
        []              160   []              160   []
        []              68    []              68    []
        []              68    []              68    []
        []              68    []              68    []
        []              160   []              160   []
        []              68    []              68    []
        []              68    []              68    []
        []              10    []              10    []
        []              10    []              10    []
        []              10    []              10    []
        [[0, 2, 1]]     160   [[0, 2, 1]]     160   [[0, 2, 1]]
        []              10    []              10    []
        [[0, 2, 1]]     160   [[0, 2, 1]]     160   [[0, 2, 1]]
        [[0, 2, 1]]     10    [[0, 2, 1]]     10    [[0, 2, 1]]
        [[0, 2, 1]]     10    [[0, 2, 1]]     10    [[0, 2, 1]]
        [[0, 2, 1]]     10    [[0, 2, 1]]     10    [[0, 2, 1]]
        [[0, 2, 1]]     10    [[0, 2, 1]]     10    [[0, 2, 1]]
        [[0, 2, 1]]     10    [[0, 2, 1]]     10    [[0, 2, 1]]
        [[0, 2, 1]]     10    [[0, 2, 1]]     10    [[0, 2, 1]]
        []              68    []              68    []
        []              160   []              160   []
        []              68    []              68    []
        []              68    []              68    []
        []              10    []              10    []
        []              10    []              10    []
        [[2, 1, 0]]     10    [[2, 1, 0]]     10    [[2, 1, 0]]
        [[2, 1, 0]]     10    [[2, 1, 0]]     10    [[2, 1, 0]]
        []              160   []              160   []
        []              68    []              68    []
        []              68    []              68    []
        []              68    []              68    []
        []              10    []              10    []
        [[2, 1, 0]]     160   [[2, 1, 0]]     160   [[2, 1, 0]]
        []              10    []              10    []
        []              10    []              10    []
        [[2, 1, 0]]     160   [[2, 1, 0]]     160   [[2, 1, 0]]
        []              10    []              10    []
        [[2, 1, 0]]     10    [[2, 1, 0]]     10    [[2, 1, 0]]
        [[2, 1, 0]]     10    [[2, 1, 0]]     10    [[2, 1, 0]]
        [[2, 1, 0]]     10    [[2, 1, 0]]     10    [[2, 1, 0]]
        [[2, 1, 0]]     10    [[2, 1, 0]]     10    [[2, 1, 0]]
        [[2, 1, 0]]     10    [[2, 1, 0]]     10    [[2, 1, 0]]
        [[2, 1, 0]]     10    [[2, 1, 0]]     10    [[2, 1, 0]]
        []              78    []              78    []
        []              170   []              170   []
        []              78    []              78    []
        []              78    []              78    []
        []              78    []              78    []
        []              170   []              170   []
        []              78    []              78    []
        []              78    []              78    []
        []              228   []              228   []
        []              228   []              228   []
        [[1, 0, 2]]     228   [[1, 0, 2]]     228   [[1, 0, 2]]
        [[1, 0, 2]]     228   [[1, 0, 2]]     228   [[1, 0, 2]]
        [[1, 0, 2]]     78    [[1, 0, 2]]     78    [[1, 0, 2]]
        []              170   []              170   []
        [[1, 0, 2]]     78    [[1, 0, 2]]     78    [[1, 0, 2]]
        []              170   []              170   []
        []              170   []              170   []
        []              170   []              170   []
        [[1, 0, 2]]     78    [[1, 0, 2]]     78    [[1, 0, 2]]
        [[1, 0, 2]]     78    [[1, 0, 2]]     78    [[1, 0, 2]]
        [[1, 0, 2]]     78    [[1, 0, 2]]     78    [[1, 0, 2]]
        [[1, 0, 2]]     78    [[1, 0, 2]]     78    [[1, 0, 2]]
        [[1, 0, 2]]     78    [[1, 0, 2]]     78    [[1, 0, 2]]
        [[1, 0, 2]]     78    [[1, 0, 2]]     78    [[1, 0, 2]]
        []              160   []              160   []
        []              68    []              68    []
        []              68    []              68    []
        []              68    []              68    []
        []              160   []              160   []
        []              68    []              68    []
        []              68    []              68    []
        []              68    []              68    []
        []              10    []              10    []
        []              10    []              10    []
        [[1, 0, 2]]     10    [[1, 0, 2]]     10    [[1, 0, 2]]
        [[1, 0, 2]]     10    [[1, 0, 2]]     10    [[1, 0, 2]]
        [[1, 0, 2]]     160   [[1, 0, 2]]     160   [[1, 0, 2]]
        []              10    []              10    []
        [[1, 0, 2]]     160   [[1, 0, 2]]     160   [[1, 0, 2]]
        []              10    []              10    []
        []              10    []              10    []
        []              10    []              10    []
        [[1, 0, 2]]     10    [[1, 0, 2]]     10    [[1, 0, 2]]
        [[1, 0, 2]]     10    [[1, 0, 2]]     10    [[1, 0, 2]]
        [[1, 0, 2]]     10    [[1, 0, 2]]     10    [[1, 0, 2]]
        [[1, 0, 2]]     10    [[1, 0, 2]]     10    [[1, 0, 2]]
        [[1, 0, 2]]     10    [[1, 0, 2]]     10    [[1, 0, 2]]
        [[1, 0, 2]]     10    [[1, 0, 2]]     10    [[1, 0, 2]]
        []              170   []              170   []
        []              78    []              78    []
        []              78    []              78    []
        []              78    []              78    []
        []              228   []              228   []
        []              228   []              228   []
        [[2, 1, 0]]     228   [[2, 1, 0]]     228   [[2, 1, 0]]
        [[2, 1, 0]]     228   [[2, 1, 0]]     228   [[2, 1, 0]]
        []              78    []              78    []
        []              170   []              170   []
        []              78    []              78    []
        []              78    []              78    []
        []              170   []              170   []
        [[2, 1, 0]]     78    [[2, 1, 0]]     78    [[2, 1, 0]]
        []              170   []              170   []
        []              170   []              170   []
        [[2, 1, 0]]     78    [[2, 1, 0]]     78    [[2, 1, 0]]
        []              170   []              170   []
        [[2, 1, 0]]     78    [[2, 1, 0]]     78    [[2, 1, 0]]
        [[2, 1, 0]]     78    [[2, 1, 0]]     78    [[2, 1, 0]]
        [[2, 1, 0]]     78    [[2, 1, 0]]     78    [[2, 1, 0]]
        [[2, 1, 0]]     78    [[2, 1, 0]]     78    [[2, 1, 0]]
        [[2, 1, 0]]     78    [[2, 1, 0]]     78    [[2, 1, 0]]
        [[2, 1, 0]]     78    [[2, 1, 0]]     78    [[2, 1, 0]]
        []              228   []              228   []
        []              228   []              228   []
        [[0, 2, 1]]     228   [[0, 2, 1]]     228   [[0, 2, 1]]
        [[0, 2, 1]]     228   [[0, 2, 1]]     228   [[0, 2, 1]]
        []              170   []              170   []
        []              78    []              78    []
        []              78    []              78    []
        []              78    []              78    []
        []              170   []              170   []
        []              78    []              78    []
        []              78    []              78    []
        []              78    []              78    []
        []              170   []              170   []
        []              170   []              170   []
        []              170   []              170   []
        [[0, 2, 1]]     78    [[0, 2, 1]]     78    [[0, 2, 1]]
        []              170   []              170   []
        [[0, 2, 1]]     78    [[0, 2, 1]]     78    [[0, 2, 1]]
        [[0, 2, 1]]     78    [[0, 2, 1]]     78    [[0, 2, 1]]
        [[0, 2, 1]]     78    [[0, 2, 1]]     78    [[0, 2, 1]]
        [[0, 2, 1]]     78    [[0, 2, 1]]     78    [[0, 2, 1]]
        [[0, 2, 1]]     78    [[0, 2, 1]]     78    [[0, 2, 1]]
        [[0, 2, 1]]     78    [[0, 2, 1]]     78    [[0, 2, 1]]
        [[0, 2, 1]]     78    [[0, 2, 1]]     78    [[0, 2, 1]]
        []              238   []              238   []
        []              238   []              238   []
        [[0, 2, 1]]     238   [[0, 2, 1]]     238   [[0, 2, 1]]
        [[0, 2, 1]]     238   [[0, 2, 1]]     238   [[0, 2, 1]]
        []              238   []              238   []
        []              238   []              238   []
        [[2, 1, 0]]     238   [[2, 1, 0]]     238   [[2, 1, 0]]
        [[2, 1, 0]]     238   [[2, 1, 0]]     238   [[2, 1, 0]]
        []              238   []              238   []
        []              238   []              238   []
        [[1, 0, 2]]     238   [[1, 0, 2]]     238   [[1, 0, 2]]
        [[1, 0, 2]]     238   [[1, 0, 2]]     238   [[1, 0, 2]]
        [[1, 0, 2]]     238   [[1, 0, 2]]     238   [[1, 0, 2]]
        [[2, 1, 0]]     238   [[2, 1, 0]]     238   [[2, 1, 0]]
        [[1, 0, 2]]     238   [[1, 0, 2]]     238   [[1, 0, 2]]
        [[0, 2, 1]]     238   [[0, 2, 1]]     238   [[0, 2, 1]]
        [[2, 1, 0]]     238   [[2, 1, 0]]     238   [[2, 1, 0]]
        [[0, 2, 1]]     238   [[0, 2, 1]]     238   [[0, 2, 1]]
        [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]] 238   [[0, 2, 1], [1, 0, 2]]

        sage: C = graphs.CubeGraph(1)
        sage: gens = search_tree(C, [C.vertices()], lab=False)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order()
        2
        sage: C = graphs.CubeGraph(2)
        sage: gens = search_tree(C, [C.vertices()], lab=False)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order()
        8
        sage: C = graphs.CubeGraph(3)
        sage: gens = search_tree(C, [C.vertices()], lab=False)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order()
        48
        sage: C = graphs.CubeGraph(4)
        sage: gens = search_tree(C, [C.vertices()], lab=False)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order()
        384
        sage: C = graphs.CubeGraph(5)
        sage: gens = search_tree(C, [C.vertices()], lab=False)  # long time (~8 secs)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order()                    # long time
        3840
        sage: C = graphs.CubeGraph(6)
        sage: gens = search_tree(C, [C.vertices()], lab=False)  # long time (~50 secs)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order()                    # long time
        46080
    """
    n = G.order()
    if n == 0:
        if lab:
            H = copy(G)
        if dict:
            ddd = {}
        if proof:
            dd = {}
            if dict:
                return [[]], ddd, H, dd
            else:
                return [[]], H, dd
        if lab and dict:
            return [[]], ddd, H
        elif lab:
            return [[]], H
        elif dict:
            return [[]], ddd
        else:
            return [[]]

    Pi = copy(Pi)
    if proof:
        lab=True
    cdef int k, j, i, size, h, hh, index, l, hb, hzf
    cdef int hzb = -389

    #create to and from mappings to relabel vertices
    listto = G.vertices()
    ffrom = {}
    for v in listto:
        ffrom[v] = listto.index(v)
    to = {}
    for i from 0 <= i < len(listto):
        to[i] = listto[i]
    G.relabel(ffrom)
    Pi2 = []
    for cell in Pi:
        newcell = []
        for c in cell:
            newcell.append(ffrom[c])
        Pi2.append(newcell)
    Pi = Pi2

    # create the dense boolean matrix
    cdef int** M = <int**> PyMem_Malloc(n * sizeof(int*))
    if not M:
        raise MemoryError, "error allocating Boolean matrix"
    for k from 0 <= k < n:
        M[k] = <int*> PyMem_Malloc(n * sizeof(int))
        if not M[k]:
            for j from 0 <= j < k:
                PyMem_Free(M[j])
                PyMem_Free(M)
            raise MemoryError, "error allocating Boolean matrix"
        for j from 0 <= j < n:
            M[k][j] = 0
    if isinstance(G, Graph):
        for ii, j, la in G.edge_iterator():
            M[ii][j] = 1
            M[j][ii] = 1
    elif isinstance(G, DiGraph):
        for ii, j, la in G.arc_iterator():
            M[ii][j] = 1

    #begin BDM's algorithm:
    L = 100
    state = 1
    W = {}
    v = {}
    Lambda = {}
    nu = {}
    Phi = {}
    Omega = {}
    e = {}
    zf = {}
    zb = {}
    Theta = []
    output = []
    if verbosity > 1:
        t = cputime()
    if verbosity > 2:
        rho = None
    if verbosity > 4:
        eta = None

    while not state is None:
        if verbosity > 0:
            print '-----'
        if verbosity > 2:
            print 'k: ' + str(k)
            print 'eta: ' + str(eta)
            print 'nu: ' + str(nu)
            print 'rho: ' + str(rho)
        if verbosity > 4:
            print 'e: ' + str(e)
            print 'hh: ' + str(hh)
            print 'hb: ' + str(hb)
            print 'h: ' + str(h)
            print 'zb: ' + str(zb)
            print 'hzb: ' + str(hzb)
            print 'hzf: ' + str(hzf)
            print 'qzb: ' + str(qzb)
            print 'v: ' + str(v)
            print 'tvh: ' + str(tvh)
            print 'W: ' + str(W)
            print 'Lambda: ' + str(Lambda)
            print 'Theta: ' + str(Theta)
        if state == 1:
            if verbosity > 0: print 'state: 1'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            size = 1
            k = 1
            h = 0
            hzb = 0
            index = 0
            l = 0
            Theta = [[i] for i from 0 <= i < n]
            nu[1] = refine(M, Pi, Pi, dig)
            hh = 2
            if not dig:
                if sat225(nu[1], n): hh = 1
            if is_discrete(nu[1]): state = 18
            else:
                W[1] = first_smallest_non_trivial(nu[1])
                v[1] = min(W[1])
                Lambda[1] = 0
                e[1] = 0
                state = 2
        elif state == 2:
            if verbosity > 0: print 'state: 2'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            k += 1
            nu[k] = perp(M, nu[k-1], v[k-1], dig)
            Lambda[k] = indicator(M, Pi, nu.values(), k, n)
            if h == 0: state = 5
            else:
                if hzf == k-1 and Lambda[k] == zf[k]:
                    hzf = k
                if not lab:
                    state = 3
                else:
                    try:
                        qzb = Lambda[k] - zb[k]
                    except:
                        qzb = -1
                    if hzb == k-1 and qzb == 0: hzb = k
                    if qzb > 0: zb[k] = Lambda[k]
                    state = 3
        elif state == 3:
            if verbosity > 0: print 'state: 3'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            if hzf <= k or (lab and qzb >= 0): state = 4 ##changed hzb to hzf, == to <=
            else: state = 6
        elif state == 4:
            if verbosity > 0: print 'state: 4'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            if is_discrete(nu[k]): state = 7
            else:
                W[k] = first_smallest_non_trivial(nu[k])
                v[k] = min(W[k])
                if dig or not sat225(nu[k], n): hh = k+1
                e[k] = 0
                state = 2
        elif state == 5:
            if verbosity > 0: print 'state: 5'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            zf[k] = Lambda[k]
            zb[k] = Lambda[k]
            state = 4
        elif state == 6:
            if verbosity > 0: print 'state: 6'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            kprime = k
            k = min([ hh-1, max(ht-1,hzb) ])
            if k == 0: k = 1 # not in BDM, broke at G = Graph({0:[], 1:[]}), Pi = [[0,1]], lab=False
            if kprime == hh: state = 13
            else:
                l = min([l+1,L])
                Omega[l] = min_cell_reps(nu[hh]) # changed Lambda to Omega
                Phi[l] = fix(nu[hh])
                state = 12
        elif state == 7:
            if verbosity > 0: print 'state: 7'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            if h == 0: state = 18
            elif k < hzf: state = 8 ## BDM had !=, broke at G = Graph({0:[],1:[],2:[]}), Pi = [[0,1,2]]
            else:
                gamma = get_permutation(eta.values(), nu.values(), list_perm=True)
                if verbosity > 3: print gamma
                if enum(G, quick=True) == G.relabel(gamma, inplace=False, quick=True): # if G^gamma == G:
                    state = 10
                else:
                    state = 8
        elif state == 8:
            if verbosity > 0: print 'state: 8'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            if (not lab) or (qzb < 0): state = 6
            elif (qzb > 0) or (k < len(rho)): state = 9
            elif (term_pnest_graph(G, nu.values(), enumer=True) > term_pnest_graph(G, rho.values(), enumer=True)): state = 9
            elif (term_pnest_graph(G, nu.values(), enumer=True) < term_pnest_graph(G, rho.values(), enumer=True)): state = 6
            else:
                gamma = get_permutation(nu.values(), rho.values(), list_perm=True)
                if verbosity > 3: print gamma
                state = 10
        elif state == 9:
            if verbosity > 0: print 'state: 9'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            rho = copy(nu)
            qzb = 0
            hb = k
            hzb = k
            zb[k+1] = Infinity
            state = 6
        elif state == 10:
            if verbosity > 0: print 'state: 10'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            l = min([l+1,L])
            Omega[l] = min_cell_reps(orbit_partition(gamma, list_perm=True))
            Phi[l] = fix(orbit_partition(gamma, list_perm=True))
            if finer( orbit_partition(gamma, list_perm=True), Theta ):
                state = 11
            else:
                Theta = vee( orbit_partition(gamma, list_perm=True), Theta )
                output.append(gamma)
                if tvc in min_cell_reps(Theta) and lab: ## added "and lab"
                    state = 11
                else:
                    k = h
                    state = 13
        elif state == 11:
            if verbosity > 0: print 'state: 11'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            k = hb
            state = 12
        elif state == 12:
            if verbosity > 0: print 'state: 12'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            if e[k] == 1:
                llll = []
                for j from 0 <= j < len(W[k]):
                    if W[k][j] in Omega[l]:
                        llll.append(W[k][j])
                W[k] = llll
            state = 13
        elif state == 13:
            if verbosity > 0: print 'state: 13'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            if k == 0: state = None
            else:
                if k > h: state = 17
                elif k == h: state = 14
                else:
                    h = k
                    tvc = min(W[k])
                    tvh = tvc
                state = 14
        elif state == 14:
            if verbosity > 0: print 'state: 14'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            for cell in Theta:
                if v[k] in cell:
                    if tvh in cell:
                        index += 1
                    else: break
            VVV = []
            for j from 0 <= j < len(W[k]):
                if W[k][j] > v[k]:
                    VVV.append(W[k][j])
            if len(VVV) != 0:
                v[k] = min(VVV)
            else:
                v[k] = Infinity
            if v[k] == Infinity: state = 16
            elif v[k] not in min_cell_reps(Theta): state = 14
            else: state = 15
        elif state == 15:
            if verbosity > 0: print 'state: 15'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            hh = min(hh,k+1)
            hzf = min(hzf,k)
            if not lab or hb < k: state = 2 # changed hzb to hb
            else:
                hb = k # changed hzb to hb
                qzb = 0
                state = 2
        elif state == 16:
            if verbosity > 0: print 'state: 16'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            if len(W[k]) == index and ht == k+1: ht = k
            size = size*index
            index = 0
            k -= 1
            state = 13
        elif state == 17:
            if verbosity > 0: print 'state: 17'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            if e[k] == 0:
                li = W[k]
                for i from 1 <= i <= l:
                    boo = True
                    for j from 1 <= j < k:
                        if v[j] not in Phi[i]:
                            boo = False
                            break
                    if boo:
                        li = [v for v in li if v in Omega[i]]
                    W[k] = li
            e[k] = 1
            VVV = []
            for j from 0 <= j < len(W[k]):
                if W[k][j] > v[k]:
                    VVV.append(W[k][j])
            if len(VVV) != 0:
                v[k] = min(VVV)
            else:
                v[k] = Infinity
            if v[k] != Infinity: state = 15
            k -= 1
            state = 13
        elif state == 18:
            if verbosity > 0: print 'state: 18'
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            h = k
            ht = k
            hzf = k
            zf[k+1] = Infinity
            eta = copy(nu)
            k -= 1
            if not lab: state = 13
            else:
                rho = copy(nu)
                hzb = k ## BDM had k+1
                hb = k ## BDM had k+1
                zb[k+2] = Infinity
                qzb = 0
                state = 13

    for k from 0 <= k < n:
        PyMem_Free(M[k])
    PyMem_Free(M)

    if lab:
        H = term_pnest_graph(G, rho.values())
    G.relabel(to)
    if dict:
        ddd = {}
        for v in G.vertices():
            if ffrom[v] != 0:
                ddd[v] = ffrom[v]
            else:
                ddd[v] = n
    if proof:
        proofpart = rho.values()[-1]
        dd = {}
        for cell in proofpart:
            dd[cell[0]] = proofpart.index(cell)
        if dict:
            return output, ddd, H, dd
        else:
            return output, H, dd
    if lab and dict:
        return output, ddd, H
    elif lab:
        return output, H
    elif dict:
        return output, ddd
    else:
        return output

def perm_group_elt(lperm):
    """
    Given a list permutation of the set {0, 1, ..., n-1},
    returns the corresponding PermutationGroupElement where
    we take 0 = n.
    """
    from sage.groups.perm_gps.permgroup_named import SymmetricGroup
    n = len(lperm)
    S = SymmetricGroup(n)
    Part = orbit_partition(lperm, list_perm=True)
    gens = []
    for z in Part:
        if len(z) > 1:
            if 0 in z:
                zed = z.index(0)
                generator = z[:zed] + [n] + z[zed+1:]
                gens.append(tuple(generator))
            else:
                gens.append(tuple(z))
    E = S(gens)
    return E

# Benchmarking functions

def all_labeled_graphs(n):
    """
    Returns all labeled graphs on n vertices {0,1,...,n-1}. Used in
    classifying isomorphism types (naive approach), and more importantly
    in benchmarking the search algorithm.

    EXAMPLE:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import search_tree, all_labeled_graphs
        sage: from sage.graphs.graph import enum
        sage: Glist = {}
        sage: Giso  = {}
        sage: for n in range(1,5): # long time (~9 secs)
        ...    Glist[n] = all_labeled_graphs(n)
        ...    Giso[n] = []
        ...    for g in Glist[n]:
        ...        a, b = search_tree(g, [range(n)])
        ...        inn = False
        ...        for gi in Giso[n]:
        ...            if enum(b) == enum(gi):
        ...                inn = True
        ...        if not inn:
        ...            Giso[n].append(b)
        sage: for n in Giso: # long time (depends on previous)
        ...    print n, len(Giso[n])
        1 1
        2 2
        3 4
        4 11
        sage.: graphs_list.show_graphs(Giso[4])
    """
    TE = []
    for i in range(n):
        for j in range(i):
            TE.append((i, j))
    m = len(TE)
    Glist= []
    for i in range(2**m):
        G = Graph()
        G.add_vertices(range(n))
        b = Integer(i).binary()
        b = '0'*(m-len(b)) + b
        for i in range(m):
            if int(b[i]):
                G.add_edge(TE[i])
        Glist.append(G)
    return Glist

def kpow(listy, k):
    """
    Returns the subset of the power set of listy consisting of subsets of size
    k. Used in all_ordered_partitions.
    """
    list = []
    if k > 1:
        for L in kpow(listy, k-1):
            for a in listy:
                if not a in L:
                    list.append([a] + L)
    if k == 1:
        for i in listy:
            list.append([i])
    return list

def all_ordered_partitions(listy):
    """
    Returns all ordered partitions of the set {0,1,...,n-1}. Used in
    benchmarking the search algorithm.
    """
    L = []
    for i in range(1,len(listy)+1):
        for cell in kpow(listy, i):
            list_remainder = [x for x in listy if x not in cell]
            remainder_partitions = all_ordered_partitions(list_remainder)
            for remainder in remainder_partitions:
                L.append( [cell] + remainder )
    if len(listy) == 0:
        return [[]]
    else:
        return L

def all_labeled_digraphs_with_loops(n):
    """
    Returns all labeled digraphs on n vertices {0,1,...,n-1}. Used in
    classifying isomorphism types (naive approach), and more importantly
    in benchmarking the search algorithm.

    EXAMPLE:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import search_tree, all_labeled_digraphs_with_loops
        sage: from sage.graphs.graph import enum
        sage: Glist = {}
        sage: Giso  = {}
        sage.: for n in range(1,4): # long time (~130 secs)
        ...    Glist[n] = all_labeled_digraphs_with_loops(n)
        ...    Giso[n] = []
        ...    for g in Glist[n]:
        ...        a, b = search_tree(g, [range(n)], dig=True)
        ...        inn = False
        ...        for gi in Giso[n]:
        ...            if enum(b) == enum(gi):
        ...                inn = True
        ...        if not inn:
        ...            Giso[n].append(b)
        sage.: for n in Giso: # long time (depends on previous)
        ...    print n, len(Giso[n])
        1 2
        2 10
        3 127
    """
    TE = []
    for i in range(n):
        for j in range(n):
            TE.append((i, j))
    m = len(TE)
    Glist= []
    for i in range(2**m):
        G = DiGraph(loops=True)
        G.add_vertices(range(n))
        b = Integer(i).binary()
        b = '0'*(m-len(b)) + b
        for j in range(m):
            if int(b[j]):
                G.add_arc(TE[j])
        Glist.append(G)
    return Glist

def all_labeled_digraphs(n):
    """
    EXAMPLES:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import search_tree, all_labeled_digraphs
        sage: from sage.graphs.graph import enum
        sage: Glist = {}
        sage: Giso  = {}
        sage: for n in range(1,4): # long time (~130 secs)
        ...       Glist[n] = all_labeled_digraphs(n)
        ...       Giso[n] = []
        ...       for g in Glist[n]:
        ...           a, b = search_tree(g, [range(n)], dig=True)
        ...           inn = False
        ...           for gi in Giso[n]:
        ...               if enum(b) == enum(gi):
        ...                   inn = True
        ...           if not inn:
        ...               Giso[n].append(b)
        sage: for n in Giso:          # long time
        ...       print n, len(Giso[n])   # long time
        1 1
        2 3
        3 16
    """
    TE = []
    for i in range(n):
        for j in range(n):
            if i != j: TE.append((i, j))
    m = len(TE)
    Glist= []
    for i in range(2**m):
        G = DiGraph(loops=True)
        G.add_vertices(range(n))
        b = Integer(i).binary()
        b = '0'*(m-len(b)) + b
        for j in range(m):
            if int(b[j]):
                G.add_arc(TE[j])
        Glist.append(G)
    return Glist
