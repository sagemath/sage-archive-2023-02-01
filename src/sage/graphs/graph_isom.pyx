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
                        (2007-07-07--14) PartitionStack and OrbitPartition
    Tom Boothby -- (2007-07-14) datastructure advice
    Robert L. Miller -- (2007-07-16--20) bug fixes

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

include '../ext/cdefs.pxi'
include '../ext/python_mem.pxi'
include '../ext/stdsage.pxi'

from sage.graphs.graph import Graph, DiGraph
from sage.misc.misc import cputime
from sage.rings.integer import Integer

cdef class OrbitPartition:
    """
    TODO: documentation

    EXAMPLES:
        sage: from sage.graphs.graph_isom import OrbitPartition
        sage: K = OrbitPartition(20)
        sage: K.find(7)
        7
        sage: K.union_find(7, 12)
        sage: K.find(12)
        7
        sage: J = OrbitPartition(20)
        sage: J.is_finer_than(K, 20)
        True
        sage: K.is_finer_than(J, 20)
        False

        sage: from sage.graphs.graph_isom import OrbitPartition
        sage: Theta1 = OrbitPartition(10)
        sage: Theta2 = OrbitPartition(10)
        sage: Theta1.union_find(0,1)
        sage: Theta1.union_find(2,3)
        sage: Theta1.union_find(3,4)
        sage: Theta1.union_find(5,6)
        sage: Theta1.union_find(8,9)
        sage: Theta2.vee_with(Theta1, 10)
        sage: for i in range(10):
        ...       print i, Theta2.find(i)
        0 0
        1 0
        2 2
        3 2
        4 2
        5 5
        6 5
        7 7
        8 8
        9 8

    """

    def __new__(self, int n):
        cdef int k
        self.elements = <int *> sage_malloc( n * sizeof(int) )
        if not self.elements:
            raise MemoryError("Error allocating memory.")
        self.sizes = <int *> sage_malloc( n * sizeof(int) )
        if not self.sizes:
            sage_free(self.elements)
            raise MemoryError("Error allocating memory.")
        for k from 0 <= k < n:
            self.elements[k] = -1
            self.sizes[k] = 1

    def __dealloc__(self):
        sage_free(self.elements)
        sage_free(self.sizes)

    def find(self, x):
        return self._find(x)

    cdef int _find(self, int x):
        if self.elements[x] == -1:
            return x
        self.elements[x] = self._find(self.elements[x])
        return self.elements[x]

    def union_find(self, a, b):
        self._union_find(a, b)

    cdef void _union_find(self, int a, int b):
        cdef int aRoot, bRoot
        aRoot = self._find(a)
        bRoot = self._find(b)
        self._union_roots(aRoot, bRoot)

    def union_roots(self, a, b):
        self._union_roots(a, b)

    cdef void _union_roots(self, int a, int b):
        if a < b:
            self.elements[b] = a
            self.sizes[b] += self.sizes[a]
        elif a > b:
            self.elements[a] = b
            self.sizes[a] += self.sizes[b]

    def is_finer_than(self, other, n):
        return self._is_finer_than(other, n) == 1

    cdef int _is_finer_than(self, OrbitPartition other, int n):
        cdef int i
        for i from 0 <= i < n:
            if self.elements[i] != -1 and other.find(self.find(i)) != other.find(i):
                return 0
        return 1

    def vee_with(self, other, n):
        self._vee_with(other, n)

    cdef void _vee_with(self, OrbitPartition other, int n):
        cdef int i
        for i from 0 <= i < n:
            if self.elements[i] == -1:
                self._union_roots(i, self.find(other.find(i)))

    cdef int _is_min_cell_rep(self, int i):
        if self.elements[i] == -1:
            return 1
        return 0

    cdef int _is_fixed(self, int i):
        if self.elements[i] == -1 and self.sizes[i] == 1:
            return 1
        return 0

cdef OrbitPartition _orbit_partition_from_list_perm(int *gamma, int n):
    cdef int i
    cdef OrbitPartition O
    O = OrbitPartition(n)
    for i from 0 <= i < n:
        if i != gamma[i]:
            O._union_find(i, gamma[i])
    return O

cdef class PartitionStack:
    """
    TODO: documentation

    EXAMPLES:

        sage: from sage.graphs.graph_isom import PartitionStack
        sage: P = PartitionStack([range(9, -1, -1)])
        sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 1, 10)
        0
        sage: P.sort_by_function(0, [2,1,2,1], 2, 10)
        0
        sage: P.sort_by_function(4, [2,1,2,1], 3, 10)
        4
        sage: P.sort_by_function(0, [0,1], 4, 10)
        0
        sage: P.sort_by_function(2, [1,0], 5, 10)
        2
        sage: P.sort_by_function(4, [1,0], 6, 10)
        4
        sage: P.sort_by_function(6, [1,0], 7, 10)
        6
        sage: P
        ({5,9,7,1,6,2,8,0,4,3})
        ({5,9,7,1},{6,2,8,0},{4},{3})
        ({5,9},{7,1},{6,2,8,0},{4},{3})
        ({5,9},{7,1},{6,2},{8,0},{4},{3})
        ({5},{9},{7,1},{6,2},{8,0},{4},{3})
        ({5},{9},{7},{1},{6,2},{8,0},{4},{3})
        ({5},{9},{7},{1},{6},{2},{8,0},{4},{3})
        ({5},{9},{7},{1},{6},{2},{8},{0},{4},{3})
        ({5},{9},{7},{1},{6},{2},{8},{0},{4},{3})
        ({5},{9},{7},{1},{6},{2},{8},{0},{4},{3})
        sage: P.is_discrete(7)
        1
        sage: P.is_discrete(6)
        0

        sage: M = graphs.PetersenGraph().am()
        sage: MM = []
        sage: for i in range(10):
        ...     MM.append([])
        ...     for j in range(10):
        ...         MM[i].append(M[i][j])
        sage: P = PartitionStack(10)
        sage: P.split_vertex(0, 1)
        sage: P.refine_by_square_matrix(MM, 1, [0], 10, 0)
        sage: P
        ({0,2,3,6,7,8,9,1,4,5})
        ({0},{2,3,6,7,8,9},{1,4,5})
        ({0},{2,3,6,7,8,9},{1,4,5})
        ({0},{2,3,6,7,8,9},{1,4,5})
        ({0},{2,3,6,7,8,9},{1,4,5})
        ({0},{2,3,6,7,8,9},{1,4,5})
        ({0},{2,3,6,7,8,9},{1,4,5})
        ({0},{2,3,6,7,8,9},{1,4,5})
        ({0},{2,3,6,7,8,9},{1,4,5})
        ({0},{2,3,6,7,8,9},{1,4,5})
        sage: P.split_vertex(1, 2)
        sage: P.refine_by_square_matrix(MM, 2, [7], 10, 0)
        sage: P
        ({0,3,7,8,9,2,6,1,4,5})
        ({0},{3,7,8,9,2,6},{1,4,5})
        ({0},{3,7,8,9},{2,6},{1},{4,5})
        ({0},{3,7,8,9},{2,6},{1},{4,5})
        ({0},{3,7,8,9},{2,6},{1},{4,5})
        ({0},{3,7,8,9},{2,6},{1},{4,5})
        ({0},{3,7,8,9},{2,6},{1},{4,5})
        ({0},{3,7,8,9},{2,6},{1},{4,5})
        ({0},{3,7,8,9},{2,6},{1},{4,5})
        ({0},{3,7,8,9},{2,6},{1},{4,5})


    """
    def __new__(self, data):
        cdef int j, k, n
        cdef PartitionStack _data
        try:
            n = int(data)
            self.entries = <int *> sage_malloc( n * sizeof(int) )
            if not self.entries:
                raise MemoryError("Error allocating memory.")
            self.levels = <int *> sage_malloc( n * sizeof(int) )
            if not self.levels:
                sage_free(self.entries)
                raise MemoryError("Error allocating memory.")
            for k from 0 <= k < n-1:
                self.entries[k] = k
                self.levels[k] = n
            self.entries[n-1] = n-1
            self.levels[n-1] = -1
        except:
            if isinstance(data, list):
                n = sum([len(datum) for datum in data])
                self.entries = <int *> sage_malloc( n * sizeof(int) )
                if not self.entries:
                    raise MemoryError("Error allocating memory.")
                self.levels = <int *> sage_malloc( n * sizeof(int) )
                if not self.levels:
                    sage_free(self.entries)
                    raise MemoryError("Error allocating memory.")
                j = 0
                k = 0
                for cell in data:
                    for entry in cell:
                        self.entries[j] = entry
                        self.levels[j] = n
                        j += 1
                    self.levels[j-1] = 0
                    self._percolate(k, j-1)
                    k = j
                self.levels[j-1] = -1
            elif isinstance(data, PartitionStack):
                _data = data
                j = 0
                while _data.levels[j] != -1: j += 1
                n = j + 1
                self.entries = <int *> sage_malloc( n * sizeof(int) )
                if not self.entries:
                    raise MemoryError("Error allocating memory.")
                self.levels = <int *> sage_malloc( n * sizeof(int) )
                if not self.levels:
                    sage_free(self.entries)
                    raise MemoryError("Error allocating memory.")
                for k from 0 <= k < n:
                    self.entries[k] = _data.entries[k]
                    self.levels[k] = _data.levels[k]
            else:
                raise ValueError("Input must be an int, a list of lists, or a PartitionStack.")

    def __dealloc__(self):
        sage_free(self.entries)
        sage_free(self.levels)

    def __repr__(self):
        k = 0
        s = ''
        while k == 0 or self.levels[k-1] != -1:
            s += '({'
            i = 0
            while i == 0 or self.levels[i-1] != -1:
                s += str(self.entries[i])
                if self.levels[i] <= k:
                    s += '},{'
                else:
                    s += ','
                i += 1
            s = s[:-2] + ')\n'
            k += 1
        return s

    def is_discrete(self, k):
        return self._is_discrete(k)

    cdef int _is_discrete(self, int k):
        cdef int i = 0
        while True:
            if self.levels[i] > k:
                return 0
            if self.levels[i] == -1: break
            i += 1
        return 1

    def num_cells(self, k):
        return self._num_cells(k)

    cdef int _num_cells(self, int k):
        cdef int i = 0, j = 1
        while self.levels[i] != -1:
        #for i from 0 <= i < n-1:
            if self.levels[i] <= k:
                j += 1
            i += 1
        return j

    def sat_225(self, k, n):
        return self._sat_225(k, n) == 1

    cdef int _sat_225(self, int k, int n):
        cdef int i, in_cell = 0
        cdef int nontrivial_cells = 0
        cdef int total_cells = self._num_cells(k)
        if n <= total_cells + 4:
            return 1
        for i from 0 <= i < n-1:
            if self.levels[i] <= k:
                if in_cell:
                    nontrivial_cells += 1
                in_cell = 0
            else:
                in_cell = 1
        if in_cell:
            nontrivial_cells += 1
        if n == total_cells + nontrivial_cells:
            return 1
        if n == total_cells + nontrivial_cells + 1:
            return 1
        return 0

    cdef int _is_min_cell_rep(self, int i, int k):
        return i == 0 or self.levels[i-1] <= k

    cdef int _is_fixed(self, int i, int k):
        """
        Assuming you already know it is a minimum cell representative.
        """
        return self.levels[i] <= k

    def split_vertex(self, v, k):
        """
        Splits the cell in self(k) containing v, putting new cells in place
        in self(k).
        """
        self._split_vertex(v, k)

    cdef int _split_vertex(self, int v, int k):
        cdef int i = 0, j
        while self.entries[i] != v:
            i += 1
        j = i
        while self.levels[i] > k:
            i += 1
        if j == 0 or self.levels[j-1] <= k:
            self._percolate(j+1, i)
        else:
            while j != 0 and self.levels[j-1] > k:
                self.entries[j] = self.entries[j-1]
                j -= 1
            self.entries[j] = v
        self.levels[j] = k
        return j

    def percolate(self, start, end):
        self._percolate(start, end)

    cdef void _percolate(self, int start, int end):
        cdef int i, temp
        for i from end >= i > start:
            if self.entries[i] < self.entries[i-1]:
                temp = self.entries[i]
                self.entries[i] = self.entries[i-1]
                self.entries[i-1] = temp

    def sort_by_function(self, start, degrees, k, n):
        cdef int i
        cdef int *degs = <int *> sage_malloc( 3 * n * sizeof(int) )
        if not degs:
            raise MemoryError("Couldn't allocate...")
        for i from 0 <= i < len(degrees):
            degs[i] = degrees[i]
        return self._sort_by_function(start, degs, k, n)
        sage_free(degs)

    cdef int _sort_by_function(self, int start, int *degrees, int k, int n):
        cdef int i, j, m = 2*n, max, max_location
        cdef int *counts = degrees + n, *output = degrees + 2*n
#        print '|'.join(['%02d'%self.entries[iii] for iii in range(n)])
#        print '|'.join(['%02d'%self.levels[iii] for iii in range(n)])
#        print '|'.join(['%02d'%degrees[iii] for iii in range(n)])
#        print '|'.join(['%02d'%counts[iii] for iii in range(n)])
#        print '|'.join(['%02d'%output[iii] for iii in range(n)])

        for i from 0 <= i < n:
            counts[i] = 0
        i = 0
        while self.levels[i+start] > k:
            counts[degrees[i]] += 1
            i += 1
        counts[degrees[i]] += 1

        # i+start is the right endpoint of the cell now
        max = counts[0]
        max_location = 0
        for j from 0 < j < n:
            if counts[j] > max:
                max = counts[j]
                max_location = j
            counts[j] += counts[j - 1]

        for j from i >= j >= 0:
            counts[degrees[j]] -= 1
            output[counts[degrees[j]]] = self.entries[start+j]

        max_location = counts[max_location]+start

        for j from 0 <= j <= i:
            self.entries[start+j] = output[j]

        j = 1
        while j < n and counts[j] <= i:
            if counts[j] > 0:
                self.levels[start + counts[j] - 1] = k
            self._percolate(start + counts[j-1], start + counts[j] - 1)
            j += 1

        return max_location

    def clear(self, k):
        self._clear(k)

    cdef void _clear(self, int k):
        cdef int i = 0, j = 0
        while self.levels[i] != -1:
            if self.levels[i] >= k:
                self.levels[i] += 1
            if self.levels[i] < k:
                self._percolate(j, i)
                j = i + 1
            i+=1

    def refine_by_square_matrix(self, G_matrix, k, alpha, n, dig):
        cdef int *_alpha, i, j
        cdef int **G
        _alpha = <int *> sage_malloc( 4 * n * sizeof(int) )
        if not _alpha:
            raise MemoryError("Memory!")
        G = <int **> sage_malloc( n * sizeof(int*) )
        if not G:
            sage_free(_alpha)
            raise MemoryError("Memory!")
        for i from 0 <= i < n:
            G[i] = <int *> sage_malloc( n * sizeof(int) )
            if not G[i]:
                for j from 0 <= j < i:
                    sage_free(G[j])
                sage_free(G)
                sage_free(_alpha)
                raise MemoryError("Memory!")
        for i from 0 <= i < n:
            for j from 0 <= j < n:
                G[i][j] = G_matrix[i][j]
        for i from 0 <= i < len(alpha):
            _alpha[i] = alpha[i]
        _alpha[len(alpha)] = -1
        self._refine_by_square_matrix(k, _alpha, n, G, dig)
        sage_free(_alpha)
        for i from 0 <= i < n:
            sage_free(G[i])
        sage_free(G)

    cdef int _refine_by_square_matrix(self, int k, int *alpha, int n, int **G, int dig):
        cdef int m = 0, j # - m iterates through alpha, the indicator cells
                          # - j iterates through the cells of the partition
        cdef int i, t, s, r # local variables:
                # - s plays a double role: outer role indicates whether
                # splitting the cell is necessary, inner role is as an index
                # for augmenting _alpha
                # - i, r iterators
                # - t: holds the first largest subcell from sort function
        cdef int invariant = 1
            # as described in [1], an indicator function Lambda(G, Pi, nu) is
            # needed to differentiate nonisomorphic branches on the search
            # tree. The condition is simply that this invariant not depend
            # on a simultaneous relabeling of the graph G, the root partition
            # Pi, and the partition nest nu. Since the function will execute
            # exactly the same way regardless of the labelling, anything that
            # does not depend on self.entries goes... at least, anything cheap
        cdef int *degrees = alpha + n # alpha assumed to be length 4*n for
                                      # extra scratch space
        while not self._is_discrete(k) and alpha[m] != -1:
            invariant += 1
            j = 0
            while j < n: # j still points at a valid cell
                invariant += 50
#                print ' '
#                print '|'.join(['%02d'%self.entries[iii] for iii in range(n)])
#                print '|'.join(['%02d'%self.levels[iii] for iii in range(n)])
#                print '|'.join(['%02d'%alpha[iii] for iii in range(n)])
#                print '|'.join(['%02d'%degrees[iii] for iii in range(n)])
#                print 'j =', j
#                print 'm =', m
                i = j; s = 0
                while True:
                    degrees[i-j] = self._degree_square_matrix(G, i, alpha[m], k)
                    if degrees[i-j] != degrees[0]: s = 1
                    i += 1
                    if self.levels[i-1] <= k: break
#                print '|'.join(['%02d'%degrees[iii] for iii in range(n)])
                # now: j points to this cell,
                #      i points to the next cell (before refinement)
                if s:
                    invariant += 10
                    t = self._sort_by_function(j, degrees, k, n)
                    # t now points to the first largest subcell
                    invariant += t + degrees[i - j - 1]
                    s = m
                    while alpha[s] != -1:
                        if alpha[s] == j: alpha[s] = t
                        s += 1
                    r = j
                    while True:
                        if r == 0 or self.levels[r-1] == k:
                            if r != t:
                                alpha[s] = r
                                s += 1
                        r += 1
                        if r >= i: break
                    alpha[s] = -1
                    while self.levels[j] > k:
                        j += 1
                    j += 1
                    invariant += (i - j)
                else: j = i
            if not dig: m += 1; continue
            # if we are looking at a digraph, also compute
            # the reverse degrees and sort by them
            j = 0
            while j < n: # j still points at a valid cell
                invariant += 20
#                print ' '
#                print '|'.join(['%02d'%self.entries[iii] for iii in range(n)])
#                print '|'.join(['%02d'%self.levels[iii] for iii in range(n)])
#                print '|'.join(['%02d'%alpha[iii] for iii in range(n)])
#                print '|'.join(['%02d'%degrees[iii] for iii in range(n)])
#                print 'j =', j
#                print 'm =', m
                i = j; s = 0
                while True:
                    degrees[i-j] = self._degree_inv_square_matrix(G, i, alpha[m], k)
                    if degrees[i-j] != degrees[0]: s = 1
                    i += 1
                    if self.levels[i-1] <= k: break
                # now: j points to this cell,
                #      i points to the next cell (before refinement)
                if s:
                    invariant += 7
                    t = self._sort_by_function(j, degrees, k, n)
                    # t now points to the first largest subcell
                    invariant += t + degrees[i - j - 1]
                    s = m
                    while alpha[s] != -1:
                        if alpha[s] == j: alpha[s] = t
                        s += 1
                    r = j
                    while True:
                        if r == 0 or self.levels[r-1] == k:
                            if r != t:
                                alpha[s] = r
                                s += 1
                        r += 1
                        if r >= i: break
                    alpha[s] = -1
                    while self.levels[j] > k:
                        j += 1
                    j += 1
                    invariant += (i - j)
                else: j = i
            m += 1
        return invariant

    def degree_square_matrix(self, G, v, W, k):
        cdef int i, j, n = len(G)
        cdef int **GG = <int **> sage_malloc( n * sizeof(int*) )
        if not GG:
            raise MemoryError("Memory!")
        for i from 0 <= i < n:
            GG[i] = <int *> sage_malloc( n * sizeof(int) )
            if not GG[i]:
                for j from 0 <= j < i:
                    sage_free(GG[j])
                sage_free(GG)
                raise MemoryError("Memory!")
        for i from 0 <= i < n:
            for j from 0 <= j < n:
                GG[i][j] = G[i][j]
        j = self._degree_square_matrix(GG, v, W, k)
        for i from 0 <= i < n:
            sage_free(GG[i])
        sage_free(GG)
        return j

    cdef int _degree_square_matrix(self, int** G, int v, int W, int k):
        """
        G is a square matrix, and W points to the beginning of a cell in the
        k-th part of the stack.
        """
        cdef int i = 0
        v = self.entries[v]
        while True:
            if G[self.entries[W]][v]:
                i += 1
            if self.levels[W] > k: W += 1
            else: break
        return i

    cdef int _degree_inv_square_matrix(self, int** G, int v, int W, int k):
        """
        G is a square matrix, and W points to the beginning of a cell in the
        k-th part of the stack.
        """
        cdef int i = 0
        v = self.entries[v]
        while True:
            if G[v][self.entries[W]]:
                i += 1
            if self.levels[W] > k: W += 1
            else: break
        return i

    cdef int _first_smallest_nontrivial(self, int *W, int k, int n):
        cdef int i = 0, j = 0, location = 0, min = n
        while True:
            W[i] = 0
            if self.levels[i] <= k:
                if i != j and n > i - j + 1:
                    n = i - j + 1
                    location = j
                j = i + 1
            if self.levels[i] == -1: break
            i += 1
        # location now points to the beginning of the first, smallest,
        # nontrivial cell
        while True:
            if min > self.entries[location]:
                min = self.entries[location]
            W[self.entries[location]] = 1
            if self.levels[location] <= k: break
            location += 1
        return min

    cdef void _get_permutation_from(self, PartitionStack zeta, int *gamma):
        cdef int i = 0

        while True:
            gamma[zeta.entries[i]] = self.entries[i]
            i += 1
            if self.levels[i-1] == -1: break

# (TODO)
# Important note: the enumeration should be kept abstract, and only comparison
# functions should be written. This takes up too much memory and time. Simply
# iterate starting with the most significant digit in the matrix, and return
# as soon a contradiction is encountered.

    cdef _enumerate_graph_from_discrete(self, int **G, int n):
        cdef int i, j
        enumeration = Integer(0)
        for i from 0 <= i < n:
            for j from 0 <= j < n:
                if G[i][j]:
                    enumeration += Integer(2)**((n-(self.entries[i]+1))*n + n-(self.entries[j]+1))
        return enumeration

cdef _enumerate_graph_with_permutation(int **G, int n, int *gamma):
    cdef int i, j
    enumeration = Integer(0)
    for i from 0 <= i < n:
        for j from 0 <= j < n:
            if G[i][j]:
                enumeration += Integer(2)**((n-(gamma[i]+1))*n + n-(gamma[j]+1))
    return enumeration

cdef _enumerate_graph(int **G, int n):
    cdef int i, j # enumeration = 0
    enumeration = Integer(0)
    for i from 0 <= i < n:
        for j from 0 <= j < n:
            if G[i][j]:
                enumeration += Integer(2)**((n-(i+1))*n + n-(j+1))
    return enumeration

def _term_pnest_graph(G, PartitionStack nu):
    """
    BDM's G(nu): returns the graph G, relabeled in the order found in
    nu[m], where m is the first index corresponding to a discrete partition.
    Assumes nu is a terminal partition nest in T(G, Pi).
    """
    cdef int i, n
    n = G.order()
    d = {}
    for i from 0 <= i < n:
        d[nu.entries[i]] = i
    H = G.copy()
    H.relabel(d)
    return H

def search_tree(G, Pi, lab=True, dig=False, dict=False, certify=False, verbosity=0):
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
        certify--     if True, return the automorphism between G and its canonical
    label. Forces lab=True.
        verbosity-- 0 - print nothing
                    1 - display state trace
                    2 - with timings
                    3 - display partition nests
                    4 - display orbit partition

    STATE DIAGRAM:
        sage: SD = DiGraph( { 1:[18,2], 2:[5,3], 3:[4,6], 4:[7,2], 5:[4], 6:[13,12], 7:[18,8,10], 8:[6,9,10], 9:[6], 10:[11,13], 11:[12], 12:[13], 13:[17,14], 14:[16,15], 15:[2], 16:[13], 17:[15,13], 18:[13] } )
        sage: SD.set_edge_label(1, 18, 'discrete')
        sage: SD.set_edge_label(4, 7, 'discrete')
        sage: SD.set_edge_label(2, 5, 'h = 0')
        sage: SD.set_edge_label(7, 18, 'h = 0')
        sage: SD.set_edge_label(7, 10, 'aut')
        sage: SD.set_edge_label(8, 10, 'aut')
        sage: SD.set_edge_label(8, 9, 'label')
        sage: SD.set_edge_label(8, 6, 'no label')
        sage: SD.set_edge_label(13, 17, 'k > h')
        sage: SD.set_edge_label(13, 14, 'k = h')
        sage: SD.set_edge_label(17, 15, 'v_k finite')
        sage: SD.set_edge_label(14, 15, 'v_k m.c.r.')
        sage: posn = {1:[ 3,-3],  2:[0,2],  3:[0, 13],  4:[3,9],  5:[3,3],  6:[16, 13], 7:[6,1],  8:[6,6],  9:[6,11], 10:[9,1], 11:[10,6], 12:[13,6], 13:[16,2], 14:[10,-6], 15:[0,-10], 16:[14,-6], 17:[16,-10], 18:[6,-4]}
        sage: SD.plot(pos=posn, vertex_size=400, vertex_colors={'#FFFFFF':range(1,19)}, edge_labels=True).save('search_tree.png')

    EXAMPLES:
        sage: import sage.graphs.graph_isom
        sage: from sage.graphs.graph_isom import search_tree
        sage: from sage.graphs.graph import enum
        sage: from sage.groups.perm_gps.permgroup import PermutationGroup # long time
        sage: from sage.graphs.graph_isom import perm_group_elt # long time

        sage: G = graphs.DodecahedralGraph()
        sage: Pi=[range(20)]
        sage: a,b = search_tree(G, Pi)
        sage: print a, enum(b)
        [[0, 19, 3, 2, 6, 5, 4, 17, 18, 11, 10, 9, 13, 12, 16, 15, 14, 7, 8, 1], [0, 1, 8, 9, 13, 14, 7, 6, 2, 3, 19, 18, 17, 4, 5, 15, 16, 12, 11, 10], [1, 8, 9, 10, 11, 12, 13, 14, 7, 6, 2, 3, 4, 5, 15, 16, 17, 18, 19, 0], [2, 1, 0, 19, 18, 11, 10, 9, 8, 7, 6, 5, 15, 14, 13, 12, 16, 17, 4, 3]] 17318942212009113839976787462421724338461987195898671092180383421848885858584973127639899792828728124797968735273000
        sage: c = search_tree(G, Pi, lab=False)
        sage: print c
        [[0, 19, 3, 2, 6, 5, 4, 17, 18, 11, 10, 9, 13, 12, 16, 15, 14, 7, 8, 1], [0, 1, 8, 9, 13, 14, 7, 6, 2, 3, 19, 18, 17, 4, 5, 15, 16, 12, 11, 10], [1, 8, 9, 10, 11, 12, 13, 14, 7, 6, 2, 3, 4, 5, 15, 16, 17, 18, 19, 0], [2, 1, 0, 19, 18, 11, 10, 9, 8, 7, 6, 5, 15, 14, 13, 12, 16, 17, 4, 3]]
        sage: DodecAut = PermutationGroup([perm_group_elt(aa) for aa in a]) # long time
        sage: DodecAut.character_table() # long time
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
        sage: DodecAut2 = PermutationGroup([perm_group_elt(cc) for cc in c]) # long time
        sage: DodecAut2.character_table() # long time
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
        [[0, 1, 2, 7, 5, 4, 6, 3, 9, 8], [0, 1, 6, 8, 5, 4, 2, 9, 3, 7], [0, 4, 3, 8, 5, 1, 9, 2, 6, 7], [1, 0, 4, 9, 6, 2, 5, 3, 7, 8], [2, 1, 0, 5, 7, 3, 6, 4, 8, 9]] 8715233764864019919698297664
        sage: c = search_tree(G, Pi, lab=False)
        sage: PAut = PermutationGroup([perm_group_elt(aa) for aa in a]) # long time
        sage: PAut.character_table() # long time
        [ 1  1  1  1  1  1  1]
        [ 1 -1  1 -1  1 -1  1]
        [ 4 -2  0  1  1  0 -1]
        [ 4  2  0 -1  1  0 -1]
        [ 5  1  1  1 -1 -1  0]
        [ 5 -1  1 -1 -1  1  0]
        [ 6  0 -2  0  0  0  1]
        sage: PAut = PermutationGroup([perm_group_elt(cc) for cc in c]) # long time
        sage: PAut.character_table() # long time
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
        [[0, 2, 1, 3, 4, 6, 5, 7], [0, 1, 4, 5, 2, 3, 6, 7], [1, 0, 3, 2, 5, 4, 7, 6]] 520239721777506480
        sage: c = search_tree(G, Pi, lab=False)

        sage: PermutationGroup([perm_group_elt(aa) for aa in a]).order() # long time
        48
        sage: PermutationGroup([perm_group_elt(cc) for cc in c]).order() # long time
        48
        sage: DodecAut.order() # long time
        120
        sage: PAut.order() # long time
        120

        sage: D = graphs.DodecahedralGraph()
        sage: a,b,c = search_tree(D, [range(20)], certify=True)
        sage: from sage.plot.plot import GraphicsArray # long time
        sage: import networkx # long time
        sage: position_D = networkx.spring_layout(D._nxg) # long time
        sage: position_b = {} # long time
        sage: for vert in position_D: # long time
        ...    position_b[c[vert]] = position_D[vert]
        sage: GraphicsArray([D.plot(pos=position_D), b.plot(pos=position_b)]).save('sage.png') # long time
        sage: c
        {0: 0, 1: 19, 2: 16, 3: 15, 4: 9, 5: 1, 6: 10, 7: 8, 8: 14, 9: 12, 10: 17, 11: 11, 12: 5, 13: 6, 14: 2, 15: 4, 16: 3, 17: 7, 18: 13, 19: 18}

    BENCHMARKS:
    The following examples are given to check modifications to the algorithm
    for optimization.

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
        sage: for G in graph3:
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
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order() # long time
        2
        sage: C = graphs.CubeGraph(2)
        sage: gens = search_tree(C, [C.vertices()], lab=False)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order() # long time
        8
        sage: C = graphs.CubeGraph(3)
        sage: gens = search_tree(C, [C.vertices()], lab=False)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order() # long time
        48
        sage: C = graphs.CubeGraph(4)
        sage: gens = search_tree(C, [C.vertices()], lab=False)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order() # long time
        384
        sage: C = graphs.CubeGraph(5)
        sage: gens = search_tree(C, [C.vertices()], lab=False)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order() # long time
        3840
        sage: C = graphs.CubeGraph(6)
        sage: gens = search_tree(C, [C.vertices()], lab=False)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order() # long time
        46080
    """
    cdef int i, j, # local variables

    cdef OrbitPartition Theta, OP
    cdef int index = 0, size = 1 # see Theorem 2.33 in [1]

    cdef int L = 100 # memory limit for storing values from fix and mcr
    cdef int **Phi # stores results from fix
    cdef int **Omega # stores results from mcr
    cdef int l = -1 # current index for storing values from fix and mcr-
                    # we start at -1 so that when we increment first,
                    # the first place we write to is 0.
    cdef int **_W # which vertices are relevant in light of the above

    cdef PartitionStack _nu, _zeta, _rho
    cdef int k_rho # the number of partitions in rho
    cdef int k = 0 # the number of partitions in nu
    cdef int h = -1 # longest common ancestor of zeta and nu:
                    # zeta[h] == nu[h], zeta[h+1] != nu[h+1]
    cdef int hb     # longest common ancestor of rho and nu:
                    # rho[hb] == nu[hb], rho[hb+1] != nu[hb+1]
    cdef int hh = 1 # the height of the oldest ancestor of nu
                    # satisfying Lemma 2.25 in [1]
    cdef int ht # smallest such that all descendants of zeta[ht] are equivalent

    cdef mpz_t *Lambda_mpz, *zf_mpz, *zb_mpz # for tracking indicator values
    # zf and zb are indicator vectors remembering Lambda[k] for zeta and rho,
    # respectively
    cdef int hzf      # the max height for which Lambda and zf agree
    cdef int hzb = -1 # the max height for which Lambda and zb agree

    cdef int **M # for the square adjacency matrix
    cdef int *_gamma # for storing permutations
    cdef int *_alpha # for storing pointers to cells of nu[k]:
                     # allocated to be length 4*n for scratch (see functions
                     # _sort_by_function and _refine_by_square_matrix)
    cdef int *_v # list of vertices determining nu
    cdef int *_e # 0 or 1, see states 12 and 17
    cdef int state # keeps track of place in algorithm
    cdef int _dig, tvc, tvh, n = G.order()

    # trivial case
    if n == 0:
        if lab:
            H = G.copy()
        if dict:
            ddd = {}
        if certify:
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

    if certify:
        lab=True

    # create to and from mappings to relabel vertices
    listto = G.vertices()
    ffrom = {}
    for vvv in listto:
        ffrom[vvv] = listto.index(vvv)
    to = {}
    for i from 0 <= i < len(listto):
        to[i] = listto[i]
    G.relabel(ffrom)
    Pi2 = []
    for cell in Pi:
        Pi2.append([ffrom[c] for c in cell])
    Pi = Pi2

    # allocate pointers
    _W = <int **> sage_malloc( 2 * (n + L) * sizeof(int *) )
    if not _W:
        raise MemoryError("Error allocating memory. Perhaps you are out?")
    M = _W + n
    Phi = _W + 2*n
    Omega = _W + 2*n + L

    # allocate pointers for GMP ints
    Lambda_mpz = <mpz_t *> sage_malloc( 3 * (n+2) * sizeof(mpz_t) )
    if not Lambda_mpz:
        sage_free(_W)
        sage_free(_gamma)
        raise MemoryError("Error allocating memory. Perhaps you are out?")
    zf_mpz = Lambda_mpz + n + 2
    zb_mpz = Lambda_mpz + 2*n + 4

    # allocate arrays
    _gamma = <int *> sage_malloc( n * ( 2 * (n + L) + 7 ) * sizeof(int) )
    if not _gamma:
        sage_free(_W)
        raise MemoryError("Error allocating memory. Perhaps you are out?")
    _alpha = _gamma + n*( 2*(n + L) + 1 )
    _v = _alpha + 4*n
    _e = _v + n
    for i from 0 <= i < n:
        _W[i] = _gamma + n + n*i
        M[i] = _gamma + n*( 1 + n + 2*L + i )

    # allocate GMP ints
    for i from 0 <= i < n+2:
        mpz_init(Lambda_mpz[i])
        mpz_init_set_si(zf_mpz[i], -1) # correspond to default values of
        mpz_init_set_si(zb_mpz[i], -1) # "infinity"
    for i from 0 <= i < L:
        Phi[i] = _gamma + n*( 1 + n + i )
        Omega[i] = _gamma + n*( 1 + n + i + L )

    # create the dense boolean matrix
    for i from 0 <= i < n:
        for j from 0 <= j < n:
            M[i][j] = 0
            _W[i][j] = 0
    if isinstance(G, Graph):
        for i, j, la in G.edge_iterator():
            M[i][j] = 1
            M[j][i] = 1
    elif isinstance(G, DiGraph):
        for i, j, la in G.edge_iterator():
            M[i][j] = 1

    # set up the rest of the variables
    _nu = PartitionStack(Pi)
    Theta = OrbitPartition(n)
    G_enum = _enumerate_graph(M, n)
    _output = []
    if dig: _dig = 1
    else: _dig = 0

    if verbosity > 1:
        t = cputime()
    if verbosity > 2:
        _rho = PartitionStack(n)
        _zeta = PartitionStack(n)
    state = 1
    while state != -1:
        if verbosity > 0:
            print '-----'
            print 'state:', state
            print '_nu'
            print [_nu.entries[iii] for iii in range(n)]
            print [_nu.levels[iii] for iii in range(n)]
            if verbosity > 1:
                t = cputime(t)
                print 'time:', t
            if verbosity > 2:
                print 'k: ' + str(k)
                print '_zeta:'
                print [_zeta.entries[iii] for iii in range(n)]
                print [_zeta.levels[iii] for iii in range(n)]
                print '_rho'
                print [_rho.entries[iii] for iii in range(n)]
                print [_rho.levels[iii] for iii in range(n)]
            if verbosity > 3:
                Thetarep = []
                for i from 0 <= i < n:
                    j = Theta._find(i)
                    didit = False
                    for celll in Thetarep:
                        if celll[0] == j:
                            celll.append(i)
                            didit = True
                    if not didit:
                        Thetarep.append([j])
                print 'Theta: ', str(Thetarep)

        if state == 1: # Entry point to algorithm
            # get alpha to point to cells of nu
            j = 1
            _alpha[0] = 0
            for i from 0 < i < n:
                if _nu.levels[i-1] == 0:
                    _alpha[j] = i
                    j += 1
            _alpha[j] = -1

            # "nu[0] := R(G, Pi, Pi)"
            _nu._refine_by_square_matrix(k, _alpha, n, M, _dig)

            if not _dig:
                if _nu._sat_225(k, n): hh = k
            if _nu._is_discrete(k): state = 18; continue

            # store the first smallest nontrivial cell in W[k], and set v[k]
            # equal to its minimum element
            _v[k] = _nu._first_smallest_nontrivial(_W[k], k, n)
            mpz_set_ui(Lambda_mpz[k], 0)
            _e[k] = 0 # see state 12, and 17
            state = 2

        elif state == 2: # Move down the search tree one level by refining nu
            k += 1

            # "nu[k] := nu[k-1] perp v[k-1]"
            _nu._clear(k)

            _alpha[0] = _nu._split_vertex(_v[k-1], k)
            _alpha[1] = -1

            i = _nu._refine_by_square_matrix(k, _alpha, n, M, _dig)

            # add one, then multiply by the invariant
            mpz_add_ui(Lambda_mpz[k], Lambda_mpz[k-1], 1)
            mpz_mul_si(Lambda_mpz[k], Lambda_mpz[k], i)

            # only if this is the first time moving down the search tree:
            if h == -1: state = 5; continue

            # update hzf
            if hzf == k-1 and mpz_cmp(Lambda_mpz[k], zf_mpz[k]) == 0: hzf = k
            if not lab: state = 3; continue

            # "qzb := cmp(Lambda[k], zb[k])"
            if mpz_cmp_si(zb_mpz[k], -1) == 0: # if "zb[k] == oo"
                qzb = -1
            else:
                qzb = mpz_cmp( Lambda_mpz[k], zb_mpz[k] )
            # update hzb
            if hzb == k-1 and qzb == 0: hzb = k

            # if Lambda[k] > zb[k], then zb[k] := Lambda[k]
            # (zb keeps track of the indicator invariants corresponding to
            # rho, the closest canonical leaf so far seen- if Lambda is
            # bigger, then rho must be about to change
            if qzb > 0: mpz_set(zb_mpz[k], Lambda_mpz[k])
            state = 3

        elif state == 3: # attempt to rule out automorphisms while moving down
                         # the tree
            if hzf <= k or (lab and qzb >= 0): # changed hzb to hzf, == to <=
                state = 4
            else: state = 6
            # if k > hzf, then we know that nu currently does not look like
            # zeta, the first terminal node encountered. Then if we are not
            # looking for a canonical label, there is no reason to continue.
            # However, if we are looking for one, and qzb < 0, i.e.
            # Lambda[k] < zb[k], then the indicator is not maximal, and we
            # can't reach a canonical leaf.

        elif state == 4: # at this point we have -not- ruled out the presence
                         # of automorphisms
            if _nu._is_discrete(k): state = 7; continue

            # store the first smallest nontrivial cell in W[k], and set v[k]
            # equal to its minimum element
            _v[k] = _nu._first_smallest_nontrivial(_W[k], k, n)

            if _dig or not _nu._sat_225(k, n): hh = k + 1
            _e[k] = 0 # see state 12, and 17
            state = 2 # continue down the tree

        elif state == 5: # alternative to 3: since we have not yet gotten
                         # zeta, there are no automorphisms to rule out.
                         # instead we record Lambda to zf and zb
                         # (see state 3)
            mpz_set(zf_mpz[k], Lambda_mpz[k])
            mpz_set(zb_mpz[k], Lambda_mpz[k])
            state = 4

        elif state == 6: # at this stage, there is no reason to continue
                         # downward, and an automorphism has not been
                         # discovered
            j = k

            # return to the longest ancestor nu[i] of nu that could have a
            # descendant equivalent to zeta or could improve on rho.
            # All terminal nodes descending from nu[hh] are known to be
            # equivalent, so i < hh. Also, if i > hzb, none of the
            # descendants of nu[i] can improve rho, since the indicator is
            # off (Lambda(nu) < Lambda(rho)). If i >= ht, then no descendant
            # of nu[i] is equivalent to zeta (see [1, p67]).
            if ht-1 > hzb:
                if ht-1 < hh-1:
                    k = ht-1
                else:
                    k = hh-1
            else:
                if hzb < hh-1:
                    k = hzb
                else:
                    k = hh-1

            # TODO: investigate the following line
            if k == -1: k = 0 # not in BDM, broke at G = Graph({0:[], 1:[]}), Pi = [[0,1]], lab=False

            if j == hh: state = 13; continue
            # recall hh: the height of the oldest ancestor of zeta for which
            # Lemma 2.25 is satsified, which implies that all terminal nodes
            # descended from there are equivalent (or simply k if 2.25 does
            # not apply). If we are looking at such a node, then the partition
            # at nu[hh] can be used for later pruning, so we store its fixed
            # set and a set of representatives of its cells
            if l < L-1: l += 1
            for i from 0 <= i < n:
                Omega[l][i] = 0 # changed Lambda to Omega
                Phi[l][i] = 0
                if _nu._is_min_cell_rep(i, hh):
                    Omega[l][i] = 1
                    if _nu._is_fixed(i, hh):
                        Phi[l][i] = 1

            state = 12

        elif state == 7: # we have just arrived at a terminal node of the
                         # search tree T(G, Pi)
            # if this is the first terminal node, go directly to 18, to
            # process zeta
            if h == -1: state = 18; continue

            # hzf is the extremal height of ancestors of both nu and zeta,
            # so if k < hzf, nu is not equivalent to zeta, i.e. there is no
            # automorphism to discover.
            # TODO: investigate why, in practice, the same does not seem to be
            # true for hzf < k... BDM had !=, not <, and this broke at
            # G = Graph({0:[],1:[],2:[]}), Pi = [[0,1,2]]
            if k < hzf: state = 8; continue

            # get the permutation corresponding to this terminal node
            _nu._get_permutation_from(_zeta, _gamma)

            if verbosity > 3:
                print 'automorphism discovered:'
                print [_gamma[iii] for iii in range(n)]

            # if G^gamma == G, the permutation is an automorphism, goto 10
            if G_enum == _enumerate_graph_with_permutation(M, n, _gamma):
                state = 10
            else:
                state = 8

        elif state == 8: # we have just ruled out the presence of automorphism
                         # and have not yet considered whether nu improves on
                         # rho
            # if we are not searching for a canonical label, there is nothing
            # to do here
            if (not lab) or (qzb < 0):
                state = 6; continue

            # if Lambda[k] > zb[k] or nu is shorter than rho, then we have
            # found an improvement for rho
            if (qzb > 0) or (k < k_rho):
                state = 9; continue

            # if G(nu) > G(rho), goto 9
            # if G(nu) < G(rho), goto 6
            # if G(nu) == G(rho), get the automorphism and goto 10
            m1 = _nu._enumerate_graph_from_discrete(M, n)
            m2 = _rho._enumerate_graph_from_discrete(M, n)

            if m1 > m2:
                state = 9; continue
            if m1 < m2:
                state = 6; continue

            _rho._get_permutation_from(_nu, _gamma)
            if verbosity > 3:
                print 'automorphism discovered:'
                print [_gamma[iii] for iii in range(n)]
            state = 10

        elif state == 9: # entering this state, nu is a best-so-far guess at
                         # the canonical label
            _rho = PartitionStack(_nu)
            k_rho = k

            qzb = 0
            hb = k
            hzb = k

            # set zb[k+1] = Infinity
            mpz_set_si(zb_mpz[k+1], -1)
            state = 6

        elif state == 10: # we have an automorphism to process
            # increment l
            if l < L - 1:
                l += 1

            # retrieve the orbit partition, and record the relevant
            # information
            # TODO: this step could be optimized. The variable OP is not
            # really necessary
            OP = _orbit_partition_from_list_perm(_gamma, n)
            for i from 0 <= i < n:
                Omega[l][i] = OP._is_min_cell_rep(i)
                Phi[l][i] = OP._is_fixed(i)

            # if each orbit of gamma is part of an orbit in Theta, then the
            # automorphism is already in the span of those we have seen
            if OP._is_finer_than(Theta, n):
                state = 11
                continue
            # otherwise, incorporate this into Theta
            Theta._vee_with(OP, n)

            # record the automorphism
            _output.append([ Integer(_gamma[i]) for i from 0 <= i < n ])

            # The variable tvc was set to be the minimum element of W[k]
            # the last time we were at state 13 and at a node descending to
            # zeta. If this is a minimal cell representative of Theta and
            # we are searching for a canonical label, goto state 11, i.e.
            # backtrack to the common ancestor of rho and nu, then goto state
            # 12, i.e. consider whether we still need to search downward from
            # there. TODO: explain why
            if Theta.elements[tvc] == -1 and lab: ## added "and lab"
                state = 11
                continue
            k = h
            state = 13

        elif state == 11: # if we are searching for a label, backtrack to the
                          # common ancestor of nu and rho
            k = hb
            state = 12

        elif state == 12: # we are looking at a branch we may have to continue
                          # to search downward on
            # e keeps track of the motion through the search tree. It is set to
            # 1 when you have just finished coming up the search tree, and are
            # at a node in the tree for which there may be more branches left
            # to explore. In this case, intersect W[k] with Omega[l], since
            # there may be an automorphism mapping one element of W[k] to
            # another, hence only one must be investigated.
            if _e[k] == 1:
                for j from 0 <= j < n:
                    if _W[k][j] and not Omega[l][j]:
                        _W[k][j] = 0
            state = 13

        elif state == 13: # hub state
            if k == -1:
                # the algorithm has finished
                state = -1; continue
            if k > h:
                # if we are not at a node of zeta
                state = 17; continue
            if k == h:
                # if we are at a node of zeta, then state 14 can rule out
                # vertices to consider
                state = 14; continue

            # thus, it must be that k < h, and this means we are done
            # searching underneath zeta[k+1], so now, k is the new longest
            # ancestor of nu and zeta:
            h = k

            # set tvc and tvh to the minimum cell representative of W[k]
            # (see states 10 and 14)
            for i from 0 <= i < n:
                if _W[k][i]:
                    tvc = i
                    break
            tvh = tvc
            state = 14

        elif state == 14: # iterate v[k] through W[k] until a minimum cell rep
                          # of Theta is found
            # The variable tvh was set to be the minimum element of W[k]
            # the last time we were at state 13 and at a node descending to
            # zeta. If this is in the same cell of Theta as v[k], increment
            # index (see Theorem 2.33 in [1])
            if Theta._find(_v[k]) == Theta._find(tvh):
                index += 1

            # find the next v[k] in W[k]
            i = _v[k] + 1
            while i < n and not _W[k][i]:
                i += 1
            if i < n:
                _v[k] = i
            else:
                # there is no new vertex to consider at this level
                _v[k] = -1
                state = 16
                continue

            # if the new v[k] is not a minimum cell representative of Theta,
            # then we already considered that rep., and that subtree was
            # isomorphic to the one corresponding to v[k]
            if Theta.elements[_v[k]] != -1: state = 14
            else:
                # otherwise, we do have a vertex to consider
                state = 15

        elif state == 15: # we have a new vertex, v[k], that we must split on
            # hh is smallest such that nu[hh] satisfies Lemma 2.25. If it is
            # larger than k+1, it must be modified, since we are changing that
            # part
            if k + 1 < hh:
                hh = k + 1
            # hzf is maximal such that indicators line up for nu and zeta
            if k < hzf:
                hzf = k
            if not lab or hb < k: # changed hzb to hb
                # in either case there is no need to update hb, which is the
                # length of the common ancestor of nu and rho
                state = 2; continue
            hb = k # changed hzb to hb
            qzb = 0
            state = 2

        elif state == 16: # backtrack one level in the search tree, recording
                          # information relevant to Theorem 2.33
            j = 0
            for i from 0 <= i < n:
                if _W[k][i]: j += 1
            if j == index and ht == k+1: ht = k
            size = size*index
            index = 0
            k -= 1
            state = 13

        elif state == 17: # you have just finished coming up the search tree,
                          # and must now consider going back down.
            if _e[k] == 0:
                # intersect W[k] with each Omega[i] such that {v_0,...,v_(k-1)}
                # is contained in Phi[i]
                for i from 0 <= i <= l:
                    # check if {v_0,...,v_(k-1)} is contained in Phi[i]
                    # i.e. fixed pointwise by the automorphisms so far seen
                    j = 0
                    while j < k and Phi[i][_v[j]]:
                        j += 1
                    # if so, only check the minimal orbit representatives
                    if j == k:
                        for j from 0 <= j < n:
                            if _W[k][j] and not Omega[i][j]:
                                _W[k][j] = 0
            _e[k] = 1 # see state 12

            # see if there is a relevant vertex to split on:
            i = _v[k] + 1
            while i < n and not _W[k][i]:
                i += 1
            if i < n:
                _v[k] = i
                state = 15
                continue
            else:
                _v[k] = -1

            # otherwise backtrack one level
            k -= 1
            state = 13

        elif state == 18: # The first time we encounter a terminal node, we
                          # come straight here to set up zeta. This is a one-
                          # time state.
            # initialize counters for zeta:
            h = k # zeta[h] == nu[h]
            ht = k # nodes descended from zeta[ht] are all equivalent
            hzf = k # max such that indicators for zeta and nu agree

            _zeta = PartitionStack(_nu)

            k -= 1
            if not lab: state = 13; continue

            _rho = PartitionStack(_nu)

            # initialize counters for rho:
            k_rho = k # number of partitions in rho
            hzb = k # max such that indicators for rho and nu agree - BDM had k+1
            hb = k # rho[hb] == nu[hb] - BDM had k+1

            qzb = 0 # Lambda[k] == zb[k], so...
            state = 13

    # deallocate the MP integers
    for i from 0 <= i < n:
        mpz_clear(Lambda_mpz[i])
        mpz_clear(zf_mpz[i])
        mpz_clear(zb_mpz[i])

    sage_free(_W)
    sage_free(Lambda_mpz)
    sage_free(_gamma)
    #mpz_clear(G_enum[0])

    if lab:
        H = _term_pnest_graph(G, _rho)
    G.relabel(to)
    if dict:
        ddd = {}
        for v in G.vertices():
            if ffrom[v] != 0:
                ddd[v] = ffrom[v]
            else:
                ddd[v] = n

    if certify:
        dd = {}
        for i from 0 <= i < n:
            dd[_rho.entries[i]] = i
            # NOTE - this should take the relabeling into account!
        if dict:
            return _output, ddd, H, dd
        else:
            return _output, H, dd
    if lab and dict:
        return _output, ddd, H
    elif lab:
        return _output, H
    elif dict:
        return _output, ddd
    else:
        return _output

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
        sage: for n in range(1,5):
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
        sage: for n in Giso:
        ...    print n, len(Giso[n])
        1 1
        2 2
        3 4
        4 11
        sage: n = 5 # long time
        sage: Glist[n] = all_labeled_graphs(n) # long time
        sage: Giso[n] = [] # long time
        sage: for g in Glist[5]: # long time
        ...    a, b = search_tree(g, [range(n)])
        ...    inn = False
        ...    for gi in Giso[n]:
        ...        if enum(b) == enum(gi):
        ...            inn = True
        ...    if not inn:
        ...        Giso[n].append(b)
        sage: print n, len(Giso[n]) # long time
        5 34
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
        for LL in kpow(listy, k-1):
            for a in listy:
                if not a in LL:
                    list.append([a] + LL)
    if k == 1:
        for i in listy:
            list.append([i])
    return list

def all_ordered_partitions(listy):
    """
    Returns all ordered partitions of the set {0,1,...,n-1}. Used in
    benchmarking the search algorithm.
    """
    LL = []
    for i in range(1,len(listy)+1):
        for cell in kpow(listy, i):
            list_remainder = [x for x in listy if x not in cell]
            remainder_partitions = all_ordered_partitions(list_remainder)
            for remainder in remainder_partitions:
                LL.append( [cell] + remainder )
    if len(listy) == 0:
        return [[]]
    else:
        return LL

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
        sage.: for n in range(1,4):
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
        sage.: for n in Giso:
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
                G.add_edge(TE[j])
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
        sage: for n in range(1,4):
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
        sage: for n in Giso:
        ...       print n, len(Giso[n])
        1 1
        2 3
        3 16
        sage.: n = 4 # long time (4 minutes)
        sage.: Glist[n] = all_labeled_digraphs(n) # long time
        sage.: Giso[n] = [] # long time
        sage.: for g in Glist[n]: # long time
        ...       a, b = search_tree(g, [range(n)], dig=True)
        ...       inn = False
        ...       for gi in Giso[n]:
        ...           if enum(b) == enum(gi):
        ...               inn = True
        ...       if not inn:
        ...           Giso[n].append(b)
        sage.: print n, len(Giso[n]) # long time
        4 218

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
                G.add_edge(TE[j])
        Glist.append(G)
    return Glist

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

def number_of_graphs(n, j = None):
    graph_list = []
    n = int(n)
    l = 2**((n*(n-1))/2)
    print 'Computing canonical labels for %d labeled graphs.'%l
    k = 0
    l = l/10
    if l > 100: l = 100
    for g in all_labeled_graphs(n):
        if k%l == 0:
            print k
        k += 1
        g = g.canonical_label()
        if g not in graph_list:
            graph_list.append(g)
    return len(graph_list)

