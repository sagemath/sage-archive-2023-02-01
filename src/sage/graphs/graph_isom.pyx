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

NOTES:
    1. Often we assume that G is a graph on vertices {0,1,...,n-1}.
    2. There is no s == loads(dumps(s)) type test since none of the classes
        defined here are meant to be instantiated for longer than the algorithm
        runs (i.e. pickling is not relevant here).
"""

#*****************************************************************************
#      Copyright (C) 2006 - 2007 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from sage.graphs.graph import GenericGraph, Graph, DiGraph
from sage.misc.misc import cputime
from sage.rings.integer cimport Integer

cdef class OrbitPartition:
    """
    An OrbitPartition is simply a partition which keeps track of the orbits of
    the part of the automorphism group so far discovered. Essentially a
    union-find datastructure.

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
        """
        Returns an element of the cell which depends only on the cell.

        EXAMPLE:
            sage: from sage.graphs.graph_isom import OrbitPartition
            sage: K = OrbitPartition(20)

        0 and 1 begin in different cells:
            sage: K.find(0)
            0
            sage: K.find(1)
            1

        Now we put them in the same cell:
            sage: K.union_find(0,1)
            sage: K.find(0)
            0
            sage: K.find(1)
            0

        """
        return self._find(x)

    cdef int _find(self, int x):
        if self.elements[x] == -1:
            return x
        self.elements[x] = self._find(self.elements[x])
        return self.elements[x]

    def union_find(self, a, b):
        """
        Merges the cells containing a and b.

        EXAMPLE:
            sage: from sage.graphs.graph_isom import OrbitPartition
            sage: K = OrbitPartition(20)

        0 and 1 begin in different cells:
            sage: K.find(0)
            0
            sage: K.find(1)
            1

        Now we put them in the same cell:
            sage: K.union_find(0,1)
            sage: K.find(0)
            0
            sage: K.find(1)
            0

        """
        self._union_find(a, b)

    cdef int _union_find(self, int a, int b):
        cdef int aRoot, bRoot
        aRoot = self._find(a)
        bRoot = self._find(b)
        self._union_roots(aRoot, bRoot)
        return aRoot != bRoot

    cdef void _union_roots(self, int a, int b):
        if a < b:
            self.elements[b] = a
            self.sizes[b] += self.sizes[a]
        elif a > b:
            self.elements[a] = b
            self.sizes[a] += self.sizes[b]

    def is_finer_than(self, other, n):
        """
        Partition P is finer than partition Q if every cell of P is a subset of
        a cell of Q.

        EXAMPLE:
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

        """
        return self._is_finer_than(other, n) == 1

    cdef int _is_finer_than(self, OrbitPartition other, int n):
        cdef int i
        for i from 0 <= i < n:
            if self.elements[i] != -1 and other.find(self.find(i)) != other.find(i):
                return 0
        return 1

    def vee_with(self, other, n):
        """
        Merges the minimal number of cells such that other is finer than self.

        EXAMPLE:
            sage: from sage.graphs.graph_isom import OrbitPartition
            sage: K = OrbitPartition(20)
            sage: K.union_find(7, 12)
            sage: J = OrbitPartition(20)
            sage: J.is_finer_than(K, 20)
            True
            sage: K.is_finer_than(J, 20)
            False
            sage: J.vee_with(K, 20)
            sage: K.is_finer_than(J, 20)
            True

        """
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

    def incorporate_permutation(self, gamma):
        """
        Unions the cells of self which contain common elements of some orbit of
        gamma.

        INPUT:
        gamma -- a permutation, in list notation

        EXAMPLE:
            sage: from sage.graphs.graph_isom import OrbitPartition
            sage: O = OrbitPartition(9)
            sage: O.incorporate_permutation([0,1,3,2,5,6,7,4,8])
            sage: for i in range(9):
            ...    print i, O.find(i)
            0 0
            1 1
            2 2
            3 2
            4 4
            5 4
            6 4
            7 4
            8 8

        """
        cdef int k, n = len(gamma)
        cdef int *_gamma = <int *> sage_malloc( n * sizeof(int) )
        if not _gamma:
            raise MemoryError("Error allocating memory.")
        for k from 0 <= k < n:
            _gamma[k] = gamma[k]
        self._incorporate_permutation(_gamma, n)
        sage_free(_gamma)

    cdef int _incorporate_permutation(self, int *gamma, int n):
        cdef int i, ch = 0, k
        for i from 0 <= i < n:
            k = self._union_find(i, gamma[i])
            if (not ch) and k:
                ch = 1
        return ch

cdef OrbitPartition orbit_partition_from_list_perm(int *gamma, int n):
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
        sage: from sage.graphs.base.sparse_graph import SparseGraph
        sage: P = PartitionStack([range(9, -1, -1)])
        sage: P.set_k(1)
        sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 10)
        0
        sage: P.set_k(2)
        sage: P.sort_by_function(0, [2,1,2,1], 10)
        0
        sage: P.set_k(3)
        sage: P.sort_by_function(4, [2,1,2,1], 10)
        4
        sage: P.set_k(4)
        sage: P.sort_by_function(0, [0,1], 10)
        0
        sage: P.set_k(5)
        sage: P.sort_by_function(2, [1,0], 10)
        2
        sage: P.set_k(6)
        sage: P.sort_by_function(4, [1,0], 10)
        4
        sage: P.set_k(7)
        sage: P.sort_by_function(6, [1,0], 10)
        6
        sage: P
        (5,9,7,1,6,2,8,0,4,3)
        (5,9,7,1|6,2,8,0|4|3)
        (5,9|7,1|6,2,8,0|4|3)
        (5,9|7,1|6,2|8,0|4|3)
        (5|9|7,1|6,2|8,0|4|3)
        (5|9|7|1|6,2|8,0|4|3)
        (5|9|7|1|6|2|8,0|4|3)
        (5|9|7|1|6|2|8|0|4|3)
        sage: P.is_discrete()
        1
        sage: P.set_k(6)
        sage: P.is_discrete()
        0

        sage: G = SparseGraph(10)
        sage: for i,j,_ in graphs.PetersenGraph().edge_iterator():
        ...    G.add_arc(i,j)
        ...    G.add_arc(j,i)
        sage: P = PartitionStack(10)
        sage: P.set_k(1)
        sage: P.split_vertex(0)
        sage: P.refine(G, [0], 10, 0, 1)
        sage: P
        (0,2,3,6,7,8,9,1,4,5)
        (0|2,3,6,7,8,9|1,4,5)
        sage: P.set_k(2)
        sage: P.split_vertex(1)
        sage: P.refine(G, [7], 10, 0, 1)
        sage: P
        (0,3,7,8,9,2,6,1,4,5)
        (0|3,7,8,9,2,6|1,4,5)
        (0|3,7,8,9|2,6|1|4,5)

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
            self.k = 0
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
                self.k = 0
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
                self.k = _data.k
            else:
                raise ValueError("Input must be an int, a list of lists, or a PartitionStack.")

    def __dealloc__(self):
        sage_free(self.entries)
        sage_free(self.levels)

    def __repr__(self):
        """
        EXAMPLE:
            sage: from sage.graphs.graph_isom import PartitionStack
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P.set_k(1)
            sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 10)
            0
            sage: P.set_k(2)
            sage: P.sort_by_function(0, [2,1,2,1], 10)
            0
            sage: P.set_k(3)
            sage: P.sort_by_function(4, [2,1,2,1], 10)
            4
            sage: P.set_k(4)
            sage: P.sort_by_function(0, [0,1], 10)
            0
            sage: P.set_k(5)
            sage: P.sort_by_function(2, [1,0], 10)
            2
            sage: P.set_k(6)
            sage: P.sort_by_function(4, [1,0], 10)
            4
            sage: P.set_k(7)
            sage: P.sort_by_function(6, [1,0], 10)
            6
            sage: P
            (5,9,7,1,6,2,8,0,4,3)
            (5,9,7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2|8,0|4|3)
            (5|9|7,1|6,2|8,0|4|3)
            (5|9|7|1|6,2|8,0|4|3)
            (5|9|7|1|6|2|8,0|4|3)
            (5|9|7|1|6|2|8|0|4|3)

        """
        k = 0
        s = ''
        while (k == 0 or self.levels[k-1] != -1) and k <= self.k:
            s += self.repr_at_k(k) + '\n'
            k += 1
        return s

    def repr_at_k(self, k):
        """
        Return the k-th line of the representation of self, i.e. the k-th
        partition in the stack.

        EXAMPLE:
            sage: from sage.graphs.graph_isom import PartitionStack
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P.set_k(1)
            sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 10)
            0
            sage: P.set_k(2)
            sage: P.sort_by_function(0, [2,1,2,1], 10)
            0
            sage: P.set_k(3)
            sage: P.sort_by_function(4, [2,1,2,1], 10)
            4
            sage: P.set_k(4)
            sage: P.sort_by_function(0, [0,1], 10)
            0
            sage: P.set_k(5)
            sage: P.sort_by_function(2, [1,0], 10)
            2
            sage: P.set_k(6)
            sage: P.sort_by_function(4, [1,0], 10)
            4
            sage: P.set_k(7)
            sage: P.sort_by_function(6, [1,0], 10)
            6
            sage: P
            (5,9,7,1,6,2,8,0,4,3)
            (5,9,7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2|8,0|4|3)
            (5|9|7,1|6,2|8,0|4|3)
            (5|9|7|1|6,2|8,0|4|3)
            (5|9|7|1|6|2|8,0|4|3)
            (5|9|7|1|6|2|8|0|4|3)

            sage: P.repr_at_k(0)
            '(5,9,7,1,6,2,8,0,4,3)'
            sage: P.repr_at_k(1)
            '(5,9,7,1|6,2,8,0|4|3)'

        """
        s = '('
        i = 0
        while i == 0 or self.levels[i-1] != -1:
            s += str(self.entries[i])
            if self.levels[i] <= k:
                s += '|'
            else:
                s += ','
            i += 1
        s = s[:-1] + ')'
        return s

    def set_k(self, k):
        """
        Sets self.k, the index of the finest partition.

            sage: from sage.graphs.graph_isom import PartitionStack
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P.set_k(1)
            sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 10)
            0
            sage: P.set_k(2)
            sage: P.sort_by_function(0, [2,1,2,1], 10)
            0
            sage: P.set_k(3)
            sage: P.sort_by_function(4, [2,1,2,1], 10)
            4
            sage: P.set_k(4)
            sage: P.sort_by_function(0, [0,1], 10)
            0
            sage: P.set_k(5)
            sage: P.sort_by_function(2, [1,0], 10)
            2
            sage: P.set_k(6)
            sage: P.sort_by_function(4, [1,0], 10)
            4
            sage: P.set_k(7)
            sage: P.sort_by_function(6, [1,0], 10)
            6
            sage: P
            (5,9,7,1,6,2,8,0,4,3)
            (5,9,7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2|8,0|4|3)
            (5|9|7,1|6,2|8,0|4|3)
            (5|9|7|1|6,2|8,0|4|3)
            (5|9|7|1|6|2|8,0|4|3)
            (5|9|7|1|6|2|8|0|4|3)

            sage: P.set_k(2)
            sage: P
            (5,9,7,1,6,2,8,0,4,3)
            (5,9,7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2,8,0|4|3)

        """
        self.k = k

    def is_discrete(self):
        """
        Returns whether the partition consists of only singletons.

        EXAMPLE:
            sage: from sage.graphs.graph_isom import PartitionStack
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P.set_k(1)
            sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 10)
            0
            sage: P.set_k(2)
            sage: P.sort_by_function(0, [2,1,2,1], 10)
            0
            sage: P.set_k(3)
            sage: P.sort_by_function(4, [2,1,2,1], 10)
            4
            sage: P.set_k(4)
            sage: P.sort_by_function(0, [0,1], 10)
            0
            sage: P.set_k(5)
            sage: P.sort_by_function(2, [1,0], 10)
            2
            sage: P.set_k(6)
            sage: P.sort_by_function(4, [1,0], 10)
            4
            sage: P.set_k(7)
            sage: P.sort_by_function(6, [1,0], 10)
            6
            sage: P
            (5,9,7,1,6,2,8,0,4,3)
            (5,9,7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2|8,0|4|3)
            (5|9|7,1|6,2|8,0|4|3)
            (5|9|7|1|6,2|8,0|4|3)
            (5|9|7|1|6|2|8,0|4|3)
            (5|9|7|1|6|2|8|0|4|3)
            sage: P.is_discrete()
            True
            sage: P.set_k(2)
            sage: P
            (5,9,7,1,6,2,8,0,4,3)
            (5,9,7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2,8,0|4|3)
            sage: P.is_discrete()
            False

        """
        return (self._is_discrete() == 1)

    cdef int _is_discrete(self):
        cdef int i = 0
        while True:
            if self.levels[i] > self.k:
                return 0
            if self.levels[i] == -1: break
            i += 1
        return 1

    def num_cells(self):
        """
        Return the number of cells in the finest partition.

        EXAMPLE:
            sage: from sage.graphs.graph_isom import PartitionStack
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P.set_k(1)
            sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 10)
            0
            sage: P
            (1,9,7,5,0,2,8,6,4,3)
            (1,9,7,5|0,2,8,6|4|3)
            sage: P.num_cells()
            4
            sage: P.set_k(2)
            sage: P.sort_by_function(0, [2,1,2,1], 10)
            0
            sage: P
            (5,9,1,7,0,2,8,6,4,3)
            (5,9,1,7|0,2,8,6|4|3)
            (5,9|1,7|0,2,8,6|4|3)
            sage: P.num_cells()
            5

        """
        return self._num_cells()

    cdef int _num_cells(self):
        cdef int i = 0, j = 1
        while self.levels[i] != -1:
        #for i from 0 <= i < n-1:
            if self.levels[i] <= self.k:
                j += 1
            i += 1
        return j

    def sat_225(self, n):
        """
        Whether the finest partition satisfies the hypotheses of Lemma 2.25 in
        [1].

        EXAMPLE:
            sage: from sage.graphs.graph_isom import PartitionStack
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P
            (0,9,8,7,6,5,4,3,2,1)
            sage: P.sat_225(10)
            False
            sage: P.set_k(1)
            sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 10)
            0
            sage: P
            (1,9,7,5,0,2,8,6,4,3)
            (1,9,7,5|0,2,8,6|4|3)
            sage: P.sat_225(10)
            False
            sage: P.set_k(2)
            sage: P.sort_by_function(0, [2,1,2,1], 10)
            0
            sage: P
            (5,9,1,7,0,2,8,6,4,3)
            (5,9,1,7|0,2,8,6|4|3)
            (5,9|1,7|0,2,8,6|4|3)
            sage: P.sat_225(10)
            False
            sage: P.set_k(3)
            sage: P.sort_by_function(4, [2,1,2,1], 10)
            4
            sage: P
            (5,9,1,7,2,6,0,8,4,3)
            (5,9,1,7|2,6,0,8|4|3)
            (5,9|1,7|2,6,0,8|4|3)
            (5,9|1,7|2,6|0,8|4|3)
            sage: P.sat_225(10)
            True

        """
        return self._sat_225(n) == 1

    cdef int _sat_225(self, int n):
        cdef int i, in_cell = 0
        cdef int nontrivial_cells = 0
        cdef int total_cells = self._num_cells()
        if n <= total_cells + 4:
            return 1
        for i from 0 <= i < n-1:
            if self.levels[i] <= self.k:
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

    def split_vertex(self, v):
        """
        Splits the cell in self(k) containing v, putting new cells in place
        in self(k).

        EXAMPLE:
            sage: from sage.graphs.graph_isom import PartitionStack
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P
            (0,9,8,7,6,5,4,3,2,1)
            sage: P.split_vertex(2)
            sage: P
            (2|0,9,8,7,6,5,4,3,1)

        """
        self._split_vertex(v)

    cdef int _split_vertex(self, int v):
        cdef int i = 0, j
        while self.entries[i] != v:
            i += 1
        j = i
        while self.levels[i] > self.k:
            i += 1
        if j == 0 or self.levels[j-1] <= self.k:
            self._percolate(j+1, i)
        else:
            while j != 0 and self.levels[j-1] > self.k:
                self.entries[j] = self.entries[j-1]
                j -= 1
            self.entries[j] = v
        self.levels[j] = self.k
        return j

    def percolate(self, start, end):
        """
        Perform one round of bubble sort, moving the smallest element to the
        front.

        EXAMPLE:
            sage: from sage.graphs.graph_isom import PartitionStack
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P
            (0,9,8,7,6,5,4,3,2,1)
            sage: P.percolate(2,7)
            sage: P
            (0,9,3,8,7,6,5,4,2,1)

        """
        self._percolate(start, end)

    cdef void _percolate(self, int start, int end):
        cdef int i, temp
        for i from end >= i > start:
            if self.entries[i] < self.entries[i-1]:
                temp = self.entries[i]
                self.entries[i] = self.entries[i-1]
                self.entries[i-1] = temp

    def sort_by_function(self, start, degrees, n):
        """
        Sort the cell starting at start using a counting sort, where degrees is
        the function giving the sort. Result is the cell is subdivided into
        cells which have elements all of the same 'degree,' in order.

        EXAMPLE:
            sage: from sage.graphs.graph_isom import PartitionStack
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P.set_k(1)
            sage: P
            (0,9,8,7,6,5,4,3,2,1)
            (0,9,8,7,6,5,4,3,2,1)
            sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 10)
            0
            sage: P
            (1,9,7,5,0,2,8,6,4,3)
            (1,9,7,5|0,2,8,6|4|3)

        """
        cdef int i
        cdef int *degs = <int *> sage_malloc( ( 3 * n + 1 ) * sizeof(int) )
        if not degs:
            raise MemoryError("Couldn't allocate...")
        for i from 0 <= i < len(degrees):
            degs[i] = degrees[i]
        X = self._sort_by_function(start, degs, n)
        sage_free(degs)
        return X

    cdef int _sort_by_function(self, int start, int *degrees, int n):
        cdef int i, j, m = 2*n, max, max_location
        cdef int *counts = degrees + n, *output = degrees + 2*n + 1
#        print '|'.join(['%02d'%self.entries[iii] for iii in range(n)])
#        print '|'.join(['%02d'%self.levels[iii] for iii in range(n)])
#        print '|'.join(['%02d'%degrees[iii] for iii in range(n)])
#        print '|'.join(['%02d'%counts[iii] for iii in range(n)])
#        print '|'.join(['%02d'%output[iii] for iii in range(n)])

        for i from 0 <= i <= n:
            counts[i] = 0
        i = 0
        while self.levels[i+start] > self.k:
            counts[degrees[i]] += 1
            i += 1
        counts[degrees[i]] += 1

        # i+start is the right endpoint of the cell now
        max = counts[0]
        max_location = 0
        for j from 0 < j <= n:
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
        while j <= n and counts[j] <= i:
            if counts[j] > 0:
                self.levels[start + counts[j] - 1] = self.k
            self._percolate(start + counts[j-1], start + counts[j] - 1)
            j += 1

        return max_location

    def clear(self):
        """
        Merges all cells in the partition stack.

        EXAMPLE:
            sage: from sage.graphs.graph_isom import PartitionStack
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P.set_k(1)
            sage: P
            (0,9,8,7,6,5,4,3,2,1)
            (0,9,8,7,6,5,4,3,2,1)
            sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 10)
            0
            sage: P
            (1,9,7,5,0,2,8,6,4,3)
            (1,9,7,5|0,2,8,6|4|3)
            sage: P
            (1,9,7,5,0,2,8,6,4,3)
            (1,9,7,5|0,2,8,6|4|3)
            sage: P.clear()
            sage: P
            (1,9,7,5,0,2,8,6,4,3)
            (1,9,7,5,0,2,8,6,4,3)

        """
        self._clear()

    cdef void _clear(self):
        cdef int i = 0, j = 0
        while self.levels[i] != -1:
            if self.levels[i] >= self.k:
                self.levels[i] += 1
            if self.levels[i] < self.k:
                self._percolate(j, i)
                j = i + 1
            i+=1

    def refine(self, CGraph G, alpha, n, dig, uif, test=False):
        """
        Implementation of Algorithm 2.5 in [1].

        EXAMPLE:
            sage: from sage.graphs.graph_isom import PartitionStack
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P.set_k(1)
            sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 10)
            0
            sage: P.set_k(2)
            sage: P.sort_by_function(0, [2,1,2,1], 10)
            0
            sage: P.set_k(3)
            sage: P.sort_by_function(4, [2,1,2,1], 10)
            4
            sage: P.set_k(4)
            sage: P.sort_by_function(0, [0,1], 10)
            0
            sage: P.set_k(5)
            sage: P.sort_by_function(2, [1,0], 10)
            2
            sage: P.set_k(6)
            sage: P.sort_by_function(4, [1,0], 10)
            4
            sage: P.set_k(7)
            sage: P.sort_by_function(6, [1,0], 10)
            6
            sage: P
            (5,9,7,1,6,2,8,0,4,3)
            (5,9,7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2,8,0|4|3)
            (5,9|7,1|6,2|8,0|4|3)
            (5|9|7,1|6,2|8,0|4|3)
            (5|9|7|1|6,2|8,0|4|3)
            (5|9|7|1|6|2|8,0|4|3)
            (5|9|7|1|6|2|8|0|4|3)
            sage: P.is_discrete()
            1
            sage: P.set_k(6)
            sage: P.is_discrete()
            0

            sage: G = SparseGraph(10)
            sage: for i,j,_ in graphs.PetersenGraph().edge_iterator():
            ...    G.add_arc(i,j)
            ...    G.add_arc(j,i)
            sage: P = PartitionStack(10)
            sage: P.set_k(1)
            sage: P.split_vertex(0)
            sage: P.refine(G, [0], 10, 0, 1)
            sage: P
            (0,2,3,6,7,8,9,1,4,5)
            (0|2,3,6,7,8,9|1,4,5)
            sage: P.set_k(2)
            sage: P.split_vertex(1)
            sage: P.refine(G, [7], 10, 0, 1)
            sage: P
            (0,3,7,8,9,2,6,1,4,5)
            (0|3,7,8,9,2,6|1,4,5)
            (0|3,7,8,9|2,6|1|4,5)

        """
        cdef int *_alpha, i, j
        _alpha = <int *> sage_malloc( ( 4 * n + 1 )* sizeof(int) )
        if not _alpha:
            raise MemoryError("Memory!")
        for i from 0 <= i < len(alpha):
            _alpha[i] = alpha[i]
        _alpha[len(alpha)] = -1
        if test:
            self.test_refine(_alpha, n, G, dig, uif)
        else:
            self._refine(_alpha, n, G, dig, uif)
        sage_free(_alpha)

    cdef int test_refine(self, int *alpha, int n, CGraph g, int dig, int uif) except? -1:
        cdef int i, j, result
        initial_partition = [] # this includes the vertex just split out...
        i = 0
        cell = []
        while i < n:
            cell.append(self.entries[i])
            while self.levels[i] > self.k:
                i += 1
                cell.append(self.entries[i])
            i += 1
            initial_partition.append(cell)
            cell = []
        #
        result = self._refine(alpha, n, g, dig, uif)
        #
        terminal_partition = []
        i = 0
        cell = []
        while i < n:
            cell.append(self.entries[i])
            while self.levels[i] > self.k:
                i += 1
                cell.append(self.entries[i])
            i += 1
            terminal_partition.append(cell)
            cell = []
        #
        if dig:
            G = DiGraph(n, loops=True)
        else:
            G = Graph(n)
        for i from 0 <= i < n:
            for j from 0 <= j < n:
                if g.has_arc_unsafe(i, j):
                    G.add_edge(i,j)
        verify_partition_refinement(G, initial_partition, terminal_partition)
        return result

    cdef int _refine(self, int *alpha, int n, CGraph G, int dig, int uif):
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
        cdef int *degrees = alpha + n # alpha assumed to be length 4*n + 1 for
                                      # extra scratch space
        while not self._is_discrete() and alpha[m] != -1:
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
                    degrees[i-j] = self._degree(G, i, alpha[m])
                    if degrees[i-j] != degrees[0]: s = 1
                    i += 1
                    if self.levels[i-1] <= self.k: break
#                print '|'.join(['%02d'%degrees[iii] for iii in range(n)])
                # now: j points to this cell,
                #      i points to the next cell (before refinement)
                if s:
                    invariant += 10
                    t = self._sort_by_function(j, degrees, n)
                    # t now points to the first largest subcell
                    invariant += t
                    s = m
                    while alpha[s] != -1:
                        if alpha[s] == j:
                            alpha[s] = t
                            break
                        s += 1
                    while alpha[s] != -1: s += 1
                    r = j
                    while True:
                        if r == j or self.levels[r-1] == self.k:
                            if r != t:
                                alpha[s] = r
                                s += 1
                        r += 1
                        if r >= i: break
                    alpha[s] = -1
                    invariant += self._degree(G, i-1, alpha[m])
                    invariant += (i - j)
                    j = i
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
                    degrees[i-j] = self._degree_inv(G, i, alpha[m])
                    if degrees[i-j] != degrees[0]: s = 1
                    i += 1
                    if self.levels[i-1] <= self.k: break
                # now: j points to this cell,
                #      i points to the next cell (before refinement)
                if s:
                    invariant += 7
                    t = self._sort_by_function(j, degrees, n)
                    # t now points to the first largest subcell
                    invariant += t
                    s = m
                    while alpha[s] != -1:
                        if alpha[s] == j:
                            alpha[s] = t
                            break
                        s += 1
                    while alpha[s] != -1: s += 1
                    r = j
                    while True:
                        if r == j or self.levels[r-1] == self.k:
                            if r != t:
                                alpha[s] = r
                                s += 1
                        r += 1
                        if r >= i: break
                    alpha[s] = -1
                    invariant += self._degree(G, i-1, alpha[m])
                    invariant += (i - j)
                    j = i
                else: j = i
            m += 1
        if uif:
            return invariant
        else:
            return 0

    def degree(self, CGraph G, v, W):
        """
        Returns the number of edges in G from self.entries[v] to a vertex in W.

        EXAMPLE:
            sage: from sage.graphs.graph_isom import PartitionStack
            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: P = PartitionStack([range(9, -1, -1)])
            sage: P
            (0,9,8,7,6,5,4,3,2,1)
            sage: G = SparseGraph(10)
            sage: G.add_arc(2,9)
            sage: G.add_arc(3,9)
            sage: G.add_arc(4,9)
            sage: P.degree(G, 1, 0)
            3

        """
        cdef int j
        j = self._degree(G, v, W)
        return j

    cdef int _degree(self, CGraph G, int v, int W):
        """
        G is a CGraph, and W points to the beginning of a cell in the
        k-th part of the stack.
        """
        cdef int i = 0
        v = self.entries[v]
        while True:
            if G.has_arc_unsafe(self.entries[W], v):
                i += 1
            if self.levels[W] > self.k: W += 1
            else: break
        return i

    cdef int _degree_inv(self, CGraph G, int v, int W):
        """
        G is a CGraph, and W points to the beginning of a cell in the
        k-th part of the stack.
        """
        cdef int i = 0
        v = self.entries[v]
        while True:
            if G.has_arc_unsafe(v, self.entries[W]):
                i += 1
            if self.levels[W] > self.k: W += 1
            else: break
        return i

    cdef int _first_smallest_nontrivial(self, int *W, int n):
        cdef int i = 0, j = 0, location = 0, min = n
        while True:
            W[i] = 0
            if self.levels[i] <= self.k:
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
            if self.levels[location] <= self.k: break
            location += 1
        return min

    cdef void _get_permutation_from(self, PartitionStack zeta, int *gamma):
        cdef int i = 0

        while True:
            gamma[zeta.entries[i]] = self.entries[i]
            i += 1
            if self.levels[i-1] == -1: break

    cdef int _compare_with(self, CGraph G, int n, PartitionStack other):
        cdef int i, j
        for i from 0 <= i < n:
            for j from 0 <= j < n:
                if G.has_arc_unsafe(self.entries[i], self.entries[j]):
                    if not G.has_arc_unsafe(other.entries[i], other.entries[j]):
                        return 1
                elif G.has_arc_unsafe(other.entries[i], other.entries[j]):
                    return -1
        return 0

cdef int _is_automorphism(CGraph G, int n, int *gamma):
    cdef int i, j
    for i from 0 <= i < n:
        for j from 0 <= j < n:
            if G.has_arc_unsafe(i, j):
                if not G.has_arc_unsafe(gamma[i], gamma[j]):
                    return 0
    return 1

def _term_pnest_graph(G, PartitionStack nu):
    """
    BDM's G(nu): returns the graph G, relabeled in the order found in
    nu[m], where m is the first index corresponding to a discrete partition.
    Assumes nu is a terminal partition nest in T(G, Pi).

    EXAMPLE:
        sage: from sage.graphs.graph_isom import PartitionStack
        sage: from sage.graphs.base.sparse_graph import SparseGraph
        sage: P = PartitionStack([range(9, -1, -1)])
        sage: P.set_k(1)
        sage: P.sort_by_function(0, [2,1,2,1,2,1,3,4,2,1], 10)
        0
        sage: P.set_k(2)
        sage: P.sort_by_function(0, [2,1,2,1], 10)
        0
        sage: P.set_k(3)
        sage: P.sort_by_function(4, [2,1,2,1], 10)
        4
        sage: P.set_k(4)
        sage: P.sort_by_function(0, [0,1], 10)
        0
        sage: P.set_k(5)
        sage: P.sort_by_function(2, [1,0], 10)
        2
        sage: P.set_k(6)
        sage: P.sort_by_function(4, [1,0], 10)
        4
        sage: P.set_k(7)
        sage: P.sort_by_function(6, [1,0], 10)
        6
        sage: P
        (5,9,7,1,6,2,8,0,4,3)
        (5,9,7,1|6,2,8,0|4|3)
        (5,9|7,1|6,2,8,0|4|3)
        (5,9|7,1|6,2|8,0|4|3)
        (5|9|7,1|6,2|8,0|4|3)
        (5|9|7|1|6,2|8,0|4|3)
        (5|9|7|1|6|2|8,0|4|3)
        (5|9|7|1|6|2|8|0|4|3)
        sage: from sage.graphs.graph_isom import _term_pnest_graph
        sage: _term_pnest_graph(graphs.PetersenGraph(), P).edges(labels=False)
        [(0, 2), (0, 6), (0, 7), (1, 2), (1, 4), (1, 8), (2, 5), (3, 4), (3, 5), (3, 7), (4, 6), (5, 9), (6, 9), (7, 8), (8, 9)]

    """
    cdef int i, j, n
    cdef CGraph M
    if isinstance(G, GenericGraph):
        n = G.order()
        H = G.copy()
    else: # G is a CGraph
        M = G
        n = M.num_verts
        if isinstance(G, SparseGraph):
            H = SparseGraph(n)
        else:
            H = DenseGraph(n)
    d = {}
    for i from 0 <= i < n:
        d[nu.entries[i]] = i
    if isinstance(G, GenericGraph):
        H.relabel(d)
    else:
        for i from 0 <= i < n:
            for j in G.out_neighbors(i):
                H.add_arc(d[i],d[j])
    return H

def search_tree(G, Pi, lab=True, dig=False, dict_rep=False, certify=False,
                verbosity=0, use_indicator_function=True, sparse=False,
                base=False, order=False):
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
        dict_rep--  if True, explain which vertices are which elements of the set
    {1,2,...,n} in the representation of the automorphism group.
        certify--     if True, return the relabeling from G to its canonical
    label. Forces lab=True.
        verbosity-- 0 - print nothing
                    1 - display state trace
                    2 - with timings
                    3 - display partition nests
                    4 - display orbit partition
                    5 - plot the part of the tree traversed during search
        use_indicator_function -- option to turn off indicator function
    (False -> slower)
        sparse -- whether to use sparse or dense representation of the graph
    (ignored if G is already a CGraph - see sage.graphs.base)
        base -- whether to return the first sequence of split vertices (used in
    computing the order of the group)
        order -- whether to return the order of the automorphism group

    STATE DIAGRAM:
        sage: SD = DiGraph( { 1:[18,2], 2:[5,3], 3:[4,6], 4:[7,2], 5:[4], 6:[13,12], 7:[18,8,10], 8:[6,9,10], 9:[6], 10:[11,13], 11:[12], 12:[13], 13:[17,14], 14:[16,15], 15:[2], 16:[13], 17:[15,13], 18:[13] }, implementation='networkx' )
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
        sage: SD.plot(pos=posn, vertex_size=400, vertex_colors={'#FFFFFF':range(1,19)}, edge_labels=True)

    NOTE:
        There is a function, called test_refine, that has the same signature as
    _refine. It calls _refine, then checks to make sure the output is sane. To
    use this, simply add 'test' to the two places this algorithm calls the
    function (states 1 and 2).

    EXAMPLES:
    The following example is due to Chris Godsil:
        sage: HS = graphs.HoffmanSingletonGraph()
        sage: clqs = (HS.complement()).cliques()
        sage: alqs = [Set(c) for c in clqs if len(c) == 15]
        sage: Y = Graph([alqs, lambda s,t: len(s.intersection(t))==0], implementation='networkx')
        sage: Y0,Y1 = Y.connected_components_subgraphs()
        sage: Y0.is_isomorphic(Y1)
        True
        sage: Y0.is_isomorphic(HS)
        True

        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree
        sage: from sage.graphs.base.sparse_graph import SparseGraph
        sage: from sage.graphs.base.dense_graph import DenseGraph
        sage: from sage.groups.perm_gps.permgroup import PermutationGroup
        sage: from sage.graphs.graph_isom import perm_group_elt

        sage: G = graphs.DodecahedralGraph()
        sage: GD = DenseGraph(20)
        sage: GS = SparseGraph(20)
        sage: for i,j,_ in G.edge_iterator():
        ...    GD.add_arc(i,j); GD.add_arc(j,i)
        ...    GS.add_arc(i,j); GS.add_arc(j,i)
        sage: Pi=[range(20)]
        sage: a,b = search_tree(G, Pi)
        sage: asp,bsp = search_tree(GS, Pi)
        sage: ade,bde = search_tree(GD, Pi)
        sage: bsg = Graph(implementation='networkx')
        sage: bdg = Graph(implementation='networkx')
        sage: for i in range(20):
        ...    for j in range(20):
        ...        if bsp.has_arc(i,j):
        ...            bsg.add_edge(i,j)
        ...        if bde.has_arc(i,j):
        ...            bdg.add_edge(i,j)
        sage: print a, b.graph6_string()
        [[0, 19, 3, 2, 6, 5, 4, 17, 18, 11, 10, 9, 13, 12, 16, 15, 14, 7, 8, 1], [0, 1, 8, 9, 13, 14, 7, 6, 2, 3, 19, 18, 17, 4, 5, 15, 16, 12, 11, 10], [1, 8, 9, 10, 11, 12, 13, 14, 7, 6, 2, 3, 4, 5, 15, 16, 17, 18, 19, 0]] S?[PG__OQ@?_?_?P?CO?_?AE?EC?Ac?@O
        sage: a == asp
        True
        sage: a == ade
        True
        sage: b == bsg
        True
        sage: b == bdg
        True
        sage: c = search_tree(G, Pi, lab=False)
        sage: print c
        [[0, 19, 3, 2, 6, 5, 4, 17, 18, 11, 10, 9, 13, 12, 16, 15, 14, 7, 8, 1], [0, 1, 8, 9, 13, 14, 7, 6, 2, 3, 19, 18, 17, 4, 5, 15, 16, 12, 11, 10], [1, 8, 9, 10, 11, 12, 13, 14, 7, 6, 2, 3, 4, 5, 15, 16, 17, 18, 19, 0]]
        sage: DodecAut = PermutationGroup([perm_group_elt(aa) for aa in a])
        sage: DodecAut.character_table()
        [                     1                      1                      1                      1                      1                      1                      1                      1                      1                      1]
        [                     1                     -1                      1                      1                     -1                      1                     -1                      1                     -1                     -1]
        [                     3                     -1                      0                     -1     -zeta5^3 - zeta5^2  zeta5^3 + zeta5^2 + 1                      0     -zeta5^3 - zeta5^2  zeta5^3 + zeta5^2 + 1                      3]
        [                     3                     -1                      0                     -1  zeta5^3 + zeta5^2 + 1     -zeta5^3 - zeta5^2                      0  zeta5^3 + zeta5^2 + 1     -zeta5^3 - zeta5^2                      3]
        [                     3                      1                      0                     -1      zeta5^3 + zeta5^2  zeta5^3 + zeta5^2 + 1                      0     -zeta5^3 - zeta5^2 -zeta5^3 - zeta5^2 - 1                     -3]
        [                     3                      1                      0                     -1 -zeta5^3 - zeta5^2 - 1     -zeta5^3 - zeta5^2                      0  zeta5^3 + zeta5^2 + 1      zeta5^3 + zeta5^2                     -3]
        [                     4                      0                      1                      0                     -1                     -1                      1                     -1                     -1                      4]
        [                     4                      0                      1                      0                      1                     -1                     -1                     -1                      1                     -4]
        [                     5                      1                     -1                      1                      0                      0                     -1                      0                      0                      5]
        [                     5                     -1                     -1                      1                      0                      0                      1                      0                      0                     -5]
        sage: DodecAut2 = PermutationGroup([perm_group_elt(cc) for cc in c])
        sage: DodecAut2.character_table()
        [                     1                      1                      1                      1                      1                      1                      1                      1                      1                      1]
        [                     1                     -1                      1                      1                     -1                      1                     -1                      1                     -1                     -1]
        [                     3                     -1                      0                     -1     -zeta5^3 - zeta5^2  zeta5^3 + zeta5^2 + 1                      0     -zeta5^3 - zeta5^2  zeta5^3 + zeta5^2 + 1                      3]
        [                     3                     -1                      0                     -1  zeta5^3 + zeta5^2 + 1     -zeta5^3 - zeta5^2                      0  zeta5^3 + zeta5^2 + 1     -zeta5^3 - zeta5^2                      3]
        [                     3                      1                      0                     -1      zeta5^3 + zeta5^2  zeta5^3 + zeta5^2 + 1                      0     -zeta5^3 - zeta5^2 -zeta5^3 - zeta5^2 - 1                     -3]
        [                     3                      1                      0                     -1 -zeta5^3 - zeta5^2 - 1     -zeta5^3 - zeta5^2                      0  zeta5^3 + zeta5^2 + 1      zeta5^3 + zeta5^2                     -3]
        [                     4                      0                      1                      0                     -1                     -1                      1                     -1                     -1                      4]
        [                     4                      0                      1                      0                      1                     -1                     -1                     -1                      1                     -4]
        [                     5                      1                     -1                      1                      0                      0                     -1                      0                      0                      5]
        [                     5                     -1                     -1                      1                      0                      0                      1                      0                      0                     -5]



        sage: G = graphs.PetersenGraph()
        sage: Pi=[range(10)]
        sage: a,b = search_tree(G, Pi)
        sage: print a, b.graph6_string()
        [[0, 1, 2, 7, 5, 4, 6, 3, 9, 8], [0, 1, 6, 8, 5, 4, 2, 9, 3, 7], [0, 4, 3, 8, 5, 1, 9, 2, 6, 7], [1, 0, 4, 9, 6, 2, 5, 3, 7, 8]] I@OZCMgs?
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
        sage: print a, b.graph6_string()
        [[0, 2, 1, 3, 4, 6, 5, 7], [0, 1, 4, 5, 2, 3, 6, 7], [1, 0, 3, 2, 5, 4, 7, 6]] GIQ\T_
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
        sage: a,b,c = search_tree(D, [range(20)], certify=True)
        sage: from sage.plot.plot import GraphicsArray
        sage: from sage.graphs.graph_fast import spring_layout_fast
        sage: position_D = spring_layout_fast(D)
        sage: position_b = {}
        sage: for vert in position_D:
        ...    position_b[c[vert]] = position_D[vert]
        sage: graphics_array([D.plot(pos=position_D), b.plot(pos=position_b)]).show()
        sage: c
        {0: 0, 1: 19, 2: 16, 3: 15, 4: 9, 5: 1, 6: 10, 7: 8, 8: 14, 9: 12, 10: 17, 11: 11, 12: 5, 13: 6, 14: 2, 15: 4, 16: 3, 17: 7, 18: 13, 19: 18}

    BENCHMARKS:
    The following examples are given to check modifications to the algorithm
    for optimization.

        sage: G = Graph({0:[]})
        sage: Pi = [[0]]
        sage: a,b = search_tree(G, Pi)
        sage: print a, b.graph6_string()
        [] @
        sage: a,b = search_tree(G, Pi, dig=True)
        sage: print a, b.graph6_string()
        [] @
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
        ...        a = str(a); b = b.graph6_string(); c = str(c); d = d.graph6_string(); e = str(e)
        ...        print a.ljust(15), b.ljust(5), c.ljust(15), d.ljust(5), e.ljust(15)
        []              A?    []              A?    []
        []              A?    []              A?    []
        [[1, 0]]        A?    [[1, 0]]        A?    [[1, 0]]
        [[1, 0]]        A?    [[1, 0]]        A?    [[1, 0]]
        []              A_    []              A_    []
        []              A_    []              A_    []
        [[1, 0]]        A_    [[1, 0]]        A_    [[1, 0]]
        [[1, 0]]        A_    [[1, 0]]        A_    [[1, 0]]

        sage: graph3 = all_labeled_graphs(3)
        sage: part3 = all_ordered_partitions(range(3))
        sage: for G in graph3:
        ...    for Pi in part3:
        ...        a,b = search_tree(G, Pi)
        ...        c,d = search_tree(G, Pi, dig=True)
        ...        e = search_tree(G, Pi, lab=False)
        ...        a = str(a); b = b.graph6_string(); c = str(c); d = d.graph6_string(); e = str(e)
        ...        print a.ljust(15), b.ljust(5), c.ljust(15), d.ljust(5), e.ljust(15)
        []              B?    []              B?    []
        []              B?    []              B?    []
        [[0, 2, 1]]     B?    [[0, 2, 1]]     B?    [[0, 2, 1]]
        [[0, 2, 1]]     B?    [[0, 2, 1]]     B?    [[0, 2, 1]]
        []              B?    []              B?    []
        []              B?    []              B?    []
        [[2, 1, 0]]     B?    [[2, 1, 0]]     B?    [[2, 1, 0]]
        [[2, 1, 0]]     B?    [[2, 1, 0]]     B?    [[2, 1, 0]]
        []              B?    []              B?    []
        []              B?    []              B?    []
        [[1, 0, 2]]     B?    [[1, 0, 2]]     B?    [[1, 0, 2]]
        [[1, 0, 2]]     B?    [[1, 0, 2]]     B?    [[1, 0, 2]]
        [[1, 0, 2]]     B?    [[1, 0, 2]]     B?    [[1, 0, 2]]
        [[2, 1, 0]]     B?    [[2, 1, 0]]     B?    [[2, 1, 0]]
        [[1, 0, 2]]     B?    [[1, 0, 2]]     B?    [[1, 0, 2]]
        [[0, 2, 1]]     B?    [[0, 2, 1]]     B?    [[0, 2, 1]]
        [[2, 1, 0]]     B?    [[2, 1, 0]]     B?    [[2, 1, 0]]
        [[0, 2, 1]]     B?    [[0, 2, 1]]     B?    [[0, 2, 1]]
        [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]] B?    [[0, 2, 1], [1, 0, 2]]
        []              BG    []              BG    []
        []              BG    []              BG    []
        [[0, 2, 1]]     BG    [[0, 2, 1]]     BG    [[0, 2, 1]]
        [[0, 2, 1]]     BG    [[0, 2, 1]]     BG    [[0, 2, 1]]
        []              BO    []              BO    []
        []              B_    []              B_    []
        []              BO    []              BO    []
        []              BO    []              BO    []
        []              BO    []              BO    []
        []              B_    []              B_    []
        []              BO    []              BO    []
        []              BO    []              BO    []
        []              BG    []              BG    []
        []              BG    []              BG    []
        []              BG    []              BG    []
        [[0, 2, 1]]     B_    [[0, 2, 1]]     B_    [[0, 2, 1]]
        []              BG    []              BG    []
        [[0, 2, 1]]     B_    [[0, 2, 1]]     B_    [[0, 2, 1]]
        [[0, 2, 1]]     BG    [[0, 2, 1]]     BG    [[0, 2, 1]]
        [[0, 2, 1]]     BG    [[0, 2, 1]]     BG    [[0, 2, 1]]
        [[0, 2, 1]]     BG    [[0, 2, 1]]     BG    [[0, 2, 1]]
        [[0, 2, 1]]     BG    [[0, 2, 1]]     BG    [[0, 2, 1]]
        [[0, 2, 1]]     BG    [[0, 2, 1]]     BG    [[0, 2, 1]]
        [[0, 2, 1]]     BG    [[0, 2, 1]]     BG    [[0, 2, 1]]
        []              BO    []              BO    []
        []              B_    []              B_    []
        []              BO    []              BO    []
        []              BO    []              BO    []
        []              BG    []              BG    []
        []              BG    []              BG    []
        [[2, 1, 0]]     BG    [[2, 1, 0]]     BG    [[2, 1, 0]]
        [[2, 1, 0]]     BG    [[2, 1, 0]]     BG    [[2, 1, 0]]
        []              B_    []              B_    []
        []              BO    []              BO    []
        []              BO    []              BO    []
        []              BO    []              BO    []
        []              BG    []              BG    []
        [[2, 1, 0]]     B_    [[2, 1, 0]]     B_    [[2, 1, 0]]
        []              BG    []              BG    []
        []              BG    []              BG    []
        [[2, 1, 0]]     B_    [[2, 1, 0]]     B_    [[2, 1, 0]]
        []              BG    []              BG    []
        [[2, 1, 0]]     BG    [[2, 1, 0]]     BG    [[2, 1, 0]]
        [[2, 1, 0]]     BG    [[2, 1, 0]]     BG    [[2, 1, 0]]
        [[2, 1, 0]]     BG    [[2, 1, 0]]     BG    [[2, 1, 0]]
        [[2, 1, 0]]     BG    [[2, 1, 0]]     BG    [[2, 1, 0]]
        [[2, 1, 0]]     BG    [[2, 1, 0]]     BG    [[2, 1, 0]]
        [[2, 1, 0]]     BG    [[2, 1, 0]]     BG    [[2, 1, 0]]
        []              BW    []              BW    []
        []              Bg    []              Bg    []
        []              BW    []              BW    []
        []              BW    []              BW    []
        []              BW    []              BW    []
        []              Bg    []              Bg    []
        []              BW    []              BW    []
        []              BW    []              BW    []
        []              Bo    []              Bo    []
        []              Bo    []              Bo    []
        [[1, 0, 2]]     Bo    [[1, 0, 2]]     Bo    [[1, 0, 2]]
        [[1, 0, 2]]     Bo    [[1, 0, 2]]     Bo    [[1, 0, 2]]
        [[1, 0, 2]]     BW    [[1, 0, 2]]     BW    [[1, 0, 2]]
        []              Bg    []              Bg    []
        [[1, 0, 2]]     BW    [[1, 0, 2]]     BW    [[1, 0, 2]]
        []              Bg    []              Bg    []
        []              Bg    []              Bg    []
        []              Bg    []              Bg    []
        [[1, 0, 2]]     BW    [[1, 0, 2]]     BW    [[1, 0, 2]]
        [[1, 0, 2]]     BW    [[1, 0, 2]]     BW    [[1, 0, 2]]
        [[1, 0, 2]]     BW    [[1, 0, 2]]     BW    [[1, 0, 2]]
        [[1, 0, 2]]     BW    [[1, 0, 2]]     BW    [[1, 0, 2]]
        [[1, 0, 2]]     BW    [[1, 0, 2]]     BW    [[1, 0, 2]]
        [[1, 0, 2]]     BW    [[1, 0, 2]]     BW    [[1, 0, 2]]
        []              B_    []              B_    []
        []              BO    []              BO    []
        []              BO    []              BO    []
        []              BO    []              BO    []
        []              B_    []              B_    []
        []              BO    []              BO    []
        []              BO    []              BO    []
        []              BO    []              BO    []
        []              BG    []              BG    []
        []              BG    []              BG    []
        [[1, 0, 2]]     BG    [[1, 0, 2]]     BG    [[1, 0, 2]]
        [[1, 0, 2]]     BG    [[1, 0, 2]]     BG    [[1, 0, 2]]
        [[1, 0, 2]]     B_    [[1, 0, 2]]     B_    [[1, 0, 2]]
        []              BG    []              BG    []
        [[1, 0, 2]]     B_    [[1, 0, 2]]     B_    [[1, 0, 2]]
        []              BG    []              BG    []
        []              BG    []              BG    []
        []              BG    []              BG    []
        [[1, 0, 2]]     BG    [[1, 0, 2]]     BG    [[1, 0, 2]]
        [[1, 0, 2]]     BG    [[1, 0, 2]]     BG    [[1, 0, 2]]
        [[1, 0, 2]]     BG    [[1, 0, 2]]     BG    [[1, 0, 2]]
        [[1, 0, 2]]     BG    [[1, 0, 2]]     BG    [[1, 0, 2]]
        [[1, 0, 2]]     BG    [[1, 0, 2]]     BG    [[1, 0, 2]]
        [[1, 0, 2]]     BG    [[1, 0, 2]]     BG    [[1, 0, 2]]
        []              Bg    []              Bg    []
        []              BW    []              BW    []
        []              BW    []              BW    []
        []              BW    []              BW    []
        []              Bo    []              Bo    []
        []              Bo    []              Bo    []
        [[2, 1, 0]]     Bo    [[2, 1, 0]]     Bo    [[2, 1, 0]]
        [[2, 1, 0]]     Bo    [[2, 1, 0]]     Bo    [[2, 1, 0]]
        []              BW    []              BW    []
        []              Bg    []              Bg    []
        []              BW    []              BW    []
        []              BW    []              BW    []
        []              Bg    []              Bg    []
        [[2, 1, 0]]     BW    [[2, 1, 0]]     BW    [[2, 1, 0]]
        []              Bg    []              Bg    []
        []              Bg    []              Bg    []
        [[2, 1, 0]]     BW    [[2, 1, 0]]     BW    [[2, 1, 0]]
        []              Bg    []              Bg    []
        [[2, 1, 0]]     BW    [[2, 1, 0]]     BW    [[2, 1, 0]]
        [[2, 1, 0]]     BW    [[2, 1, 0]]     BW    [[2, 1, 0]]
        [[2, 1, 0]]     BW    [[2, 1, 0]]     BW    [[2, 1, 0]]
        [[2, 1, 0]]     BW    [[2, 1, 0]]     BW    [[2, 1, 0]]
        [[2, 1, 0]]     BW    [[2, 1, 0]]     BW    [[2, 1, 0]]
        [[2, 1, 0]]     BW    [[2, 1, 0]]     BW    [[2, 1, 0]]
        []              Bo    []              Bo    []
        []              Bo    []              Bo    []
        [[0, 2, 1]]     Bo    [[0, 2, 1]]     Bo    [[0, 2, 1]]
        [[0, 2, 1]]     Bo    [[0, 2, 1]]     Bo    [[0, 2, 1]]
        []              Bg    []              Bg    []
        []              BW    []              BW    []
        []              BW    []              BW    []
        []              BW    []              BW    []
        []              Bg    []              Bg    []
        []              BW    []              BW    []
        []              BW    []              BW    []
        []              BW    []              BW    []
        []              Bg    []              Bg    []
        []              Bg    []              Bg    []
        []              Bg    []              Bg    []
        [[0, 2, 1]]     BW    [[0, 2, 1]]     BW    [[0, 2, 1]]
        []              Bg    []              Bg    []
        [[0, 2, 1]]     BW    [[0, 2, 1]]     BW    [[0, 2, 1]]
        [[0, 2, 1]]     BW    [[0, 2, 1]]     BW    [[0, 2, 1]]
        [[0, 2, 1]]     BW    [[0, 2, 1]]     BW    [[0, 2, 1]]
        [[0, 2, 1]]     BW    [[0, 2, 1]]     BW    [[0, 2, 1]]
        [[0, 2, 1]]     BW    [[0, 2, 1]]     BW    [[0, 2, 1]]
        [[0, 2, 1]]     BW    [[0, 2, 1]]     BW    [[0, 2, 1]]
        [[0, 2, 1]]     BW    [[0, 2, 1]]     BW    [[0, 2, 1]]
        []              Bw    []              Bw    []
        []              Bw    []              Bw    []
        [[0, 2, 1]]     Bw    [[0, 2, 1]]     Bw    [[0, 2, 1]]
        [[0, 2, 1]]     Bw    [[0, 2, 1]]     Bw    [[0, 2, 1]]
        []              Bw    []              Bw    []
        []              Bw    []              Bw    []
        [[2, 1, 0]]     Bw    [[2, 1, 0]]     Bw    [[2, 1, 0]]
        [[2, 1, 0]]     Bw    [[2, 1, 0]]     Bw    [[2, 1, 0]]
        []              Bw    []              Bw    []
        []              Bw    []              Bw    []
        [[1, 0, 2]]     Bw    [[1, 0, 2]]     Bw    [[1, 0, 2]]
        [[1, 0, 2]]     Bw    [[1, 0, 2]]     Bw    [[1, 0, 2]]
        [[1, 0, 2]]     Bw    [[1, 0, 2]]     Bw    [[1, 0, 2]]
        [[2, 1, 0]]     Bw    [[2, 1, 0]]     Bw    [[2, 1, 0]]
        [[1, 0, 2]]     Bw    [[1, 0, 2]]     Bw    [[1, 0, 2]]
        [[0, 2, 1]]     Bw    [[0, 2, 1]]     Bw    [[0, 2, 1]]
        [[2, 1, 0]]     Bw    [[2, 1, 0]]     Bw    [[2, 1, 0]]
        [[0, 2, 1]]     Bw    [[0, 2, 1]]     Bw    [[0, 2, 1]]
        [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]]
        [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]] Bw    [[0, 2, 1], [1, 0, 2]]

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
        sage: gens = search_tree(C, [C.vertices()], lab=False)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order()
        3840
        sage: C = graphs.CubeGraph(6)
        sage: gens = search_tree(C, [C.vertices()], lab=False)
        sage: PermutationGroup([perm_group_elt(aa) for aa in gens]).order()
        46080

    One can also turn off the indicator function (note- this will take longer)
        sage: D1 = DiGraph({0:[2],2:[0],1:[1]}, loops=True)
        sage: D2 = DiGraph({1:[2],2:[1],0:[0]}, loops=True)
        sage: a,b = search_tree(D1, [D1.vertices()], use_indicator_function=False)
        sage: c,d = search_tree(D2, [D2.vertices()], use_indicator_function=False)
        sage: b==d
        True

    Previously a bug, now the output is correct:
        sage: G = Graph('^????????????????????{??N??@w??FaGa?PCO@CP?AGa?_QO?Q@G?CcA??cc????Bo????{????F_')
        sage: perm = {3:15, 15:3}
        sage: H = G.relabel(perm, inplace=False)
        sage: G.canonical_label() == H.canonical_label()
        True

    Another former bug:
        sage: Graph('Fll^G').canonical_label()
        Graph on 7 vertices

        sage: g = Graph(21)
        sage: g.automorphism_group(return_group=False, order=True)
        51090942171709440000


    """
    cdef int i, j, m # local variables
    cdef int uif = 1 if use_indicator_function else 0
    cdef int _base = 1 if base else 0

    cdef OrbitPartition Theta
    cdef int index = 0 # see Theorem 2.33 in [1]
    size = Integer(1)

    cdef int L = 100 # memory limit for storing values from fix and mcr:
                     # Phi and Omega store specific information about some
                     # of the automorphisms we already know about, and they
                     # are arrays of length L
    cdef int **Phi # stores the fixed point sets of each automorphism
    cdef int **Omega # stores the minimal elements of each cell of the
                     # orbit partition
    cdef int l = -1 # current index for storing values in Phi and Omega-
                    # we start at -1 so that when we increment first,
                    # the first place we write to is 0.
    cdef int **W # for each k, W[k] is a list of the vertices to be searched
                  # down from the current partition nest, at k
                  # Phi and Omega are ultimately used to make the size of W
                  # as small as possible

    cdef PartitionStack nu, zeta, rho
    cdef int k_rho # the number of partitions in rho
    cdef int h = -1 # longest common ancestor of zeta and nu:
                    # zeta[h] == nu[h], zeta[h+1] != nu[h+1]
    cdef int hb     # longest common ancestor of rho and nu:
                    # rho[hb] == nu[hb], rho[hb+1] != nu[hb+1]
    cdef int hh = 1 # the height of the oldest ancestor of nu
                    # satisfying Lemma 2.25 in [1]
    cdef int ht # smallest such that all descendants of zeta[ht]
                # are known to be equivalent

    cdef mpz_t *Lambda_mpz, *zf_mpz, *zb_mpz # for tracking indicator values
    # zf and zb are indicator vectors remembering Lambda[k] for zeta and rho,
    # respectively
    cdef int hzf      # the max height for which Lambda and zf agree
    cdef int hzb = -1 # the max height for which Lambda and zb agree

    cdef CGraph M
    cdef int *gamma # for storing permutations
    cdef int *alpha # for storing pointers to cells of nu[k]:
                     # allocated to be length 4*n + 1 for scratch (see functions
                     # _sort_by_function and _refine)
    cdef int *v # list of vertices determining nu
    cdef int *e # 0 or 1, see states 12 and 17
    cdef int state # keeps track of place in algorithm
    cdef int _dig, tvh, n

    if isinstance(G, GenericGraph):
        n = G.order()
    elif isinstance(G, CGraph):
        M = G
        n = M.num_verts
    else:
        raise TypeError("G must be a Sage graph.")

    # trivial case
    if n == 0:
        if not (lab or dict_rep or certify or base or order):
            return [[]]
        output_tuple = [[[]]]
        if dict_rep:
            output_tuple.append({})
        if lab:
            output_tuple.append(G.copy())
        if certify:
            output_tuple.append({})
        if base:
            output_tuple.append([])
        if order:
            output_tuple.append(1)
        return tuple(output_tuple)

    # allocate int pointers
    W = <int **> sage_malloc( n * sizeof(int *) )
    Phi = <int **> sage_malloc( L * sizeof(int *) )
    Omega = <int **> sage_malloc( L * sizeof(int *) )

    # allocate GMP int pointers
    Lambda_mpz = <mpz_t *> sage_malloc( (n+2) * sizeof(mpz_t) )
    zf_mpz = <mpz_t *> sage_malloc( (n+2) * sizeof(mpz_t) )
    zb_mpz = <mpz_t *> sage_malloc( (n+2) * sizeof(mpz_t) )

    # check for memory errors
    if not (W and Phi and Omega and Lambda_mpz and zf_mpz and zb_mpz):
        sage_free(Lambda_mpz)
        sage_free(zf_mpz)
        sage_free(zb_mpz)
        sage_free(W)
        sage_free(Phi)
        sage_free(Omega)
        raise MemoryError("Error allocating memory.")

    # allocate int arrays
    gamma = <int *> sage_malloc( n * sizeof(int) )
    W[0] = <int *> sage_malloc( (n*n) * sizeof(int) )
    Phi[0] = <int *> sage_malloc( (L*n) * sizeof(int) )
    Omega[0] = <int *> sage_malloc( (L*n) * sizeof(int) )
    alpha = <int *> sage_malloc( (4*n + 1) * sizeof(int) )
    v = <int *> sage_malloc( n * sizeof(int) )
    e = <int *> sage_malloc( n * sizeof(int) )

    # check for memory errors
    if not (gamma and W[0] and Phi[0] and Omega[0] and alpha and v and e):
        sage_free(gamma)
        sage_free(W[0])
        sage_free(Phi[0])
        sage_free(Omega[0])
        sage_free(alpha)
        sage_free(v)
        sage_free(e)
        sage_free(Lambda_mpz)
        sage_free(zf_mpz)
        sage_free(zb_mpz)
        sage_free(W)
        sage_free(Phi)
        sage_free(Omega)
        raise MemoryError("Error allocating memory.")

    # setup double index arrays
    for i from 0 < i < n:
        W[i] = W[0] + n*i
    for i from 0 < i < L:
        Phi[i] = Phi[0] + n*i
    for i from 0 < i < L:
        Omega[i] = Omega[0] + n*i

    # allocate GMP ints
    for i from 0 <= i < n+2:
        mpz_init(Lambda_mpz[i])
        mpz_init_set_si(zf_mpz[i], -1) # correspond to default values of
        mpz_init_set_si(zb_mpz[i], -1) # "infinity"
        # Note that there is a potential memory leak here - if a particular
        # mpz fails to allocate, this is not checked for

    if isinstance(G, GenericGraph):
        # relabel vertices to the set {0,...,n-1}
        G = G.copy()
        ffrom = G.relabel(return_map=True)
        to = {}
        for vvv in ffrom.iterkeys():
            to[ffrom[vvv]] = vvv
        Pi2 = []
        for cell in Pi:
            Pi2.append([ffrom[c] for c in cell])
        Pi = Pi2
        if sparse:
            M = SparseGraph(n)
        else:
            M = DenseGraph(n)
        if isinstance(G, Graph):
            for i, j, la in G.edge_iterator():
                M.add_arc_unsafe(i,j)
                M.add_arc_unsafe(j,i)
        elif isinstance(G, DiGraph):
            for i, j, la in G.edge_iterator():
                M.add_arc_unsafe(i,j)

    # initialize W
    for i from 0 <= i < n:
        for j from 0 <= j < n:
            W[i][j] = 0

    # set up the rest of the variables
    nu = PartitionStack(Pi)
    Theta = OrbitPartition(n)
    output = []
    if dig: _dig = 1
    else: _dig = 0
    if certify:
        lab=True
    if _base:
        base_set = []

    if verbosity > 1:
        t = cputime()
    if verbosity > 2:
        rho = PartitionStack(n)
        zeta = PartitionStack(n)
    state = 1
    while state != -1:
        if verbosity > 0:
            print '-----'
            print 'k: ' + str(nu.k)
            print 'k_rho: ' + str(k_rho)
            print 'hh', hh
            print 'nu:'
            print nu
            if verbosity >= 2:
                t = cputime(t)
                print 'time:', t
            if verbosity >= 3:
                print 'zeta:'
                print [zeta.entries[iii] for iii in range(n)]
                print [zeta.levels[iii] for iii in range(n)]
                print 'rho'
                print [rho.entries[iii] for iii in range(n)]
                print [rho.levels[iii] for iii in range(n)]
            if verbosity >= 4:
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
            if verbosity >= 5:
                if state == 1:
                    verbose_first_time = True
                    verbose_just_refined = False
                elif verbose_first_time:
                    verbose_first_time = False
                    # here we have just gone through step 1, and must now begin
                    # to record information about the tree
                    ST_vis = DiGraph()
                    ST_vis_heights = {0:[nu.repr_at_k(0)]}
                    ST_vis.add_vertex(nu.repr_at_k(0))
                    ST_vis_current_vertex = nu.repr_at_k(0)
                    ST_vis_current_level = 0
                    #ST_vis.show(vertex_size=0)
                if state == 2:
                    verbose_just_refined = True
                elif verbose_just_refined:
                    verbose_just_refined = False
                    # here we have gone through step 2, and must record the
                    # refinement just made
                    while ST_vis_current_level > nu.k-1:
                        ST_vis_current_vertex = ST_vis.predecessors(ST_vis_current_vertex)[0]
                        ST_vis_current_level -= 1
                    if ST_vis_heights.has_key(nu.k):
                        ST_vis_heights[nu.k].append(nu.repr_at_k(nu.k))
                    else:
                        ST_vis_heights[nu.k] = [nu.repr_at_k(nu.k)]
                    ST_vis.add_edge(ST_vis_current_vertex, nu.repr_at_k(nu.k), '%d,%d'%(ST_vis_splitvert, ST_vis_invariant))
                    ST_vis_current_vertex = nu.repr_at_k(nu.k)
                    ST_vis_current_level += 1
                if state == 13 and nu.k == -1:
                    ST_vis_new_heights = {}
                    for ST_vis_k in ST_vis_heights:
                        ST_vis_new_heights[-ST_vis_k] = ST_vis_heights[ST_vis_k]
                    ST_vis.show(vertex_size=0, heights=ST_vis_new_heights, figsize=[30,10], edge_labels=True, edge_colors={(.6,.6,.6):ST_vis.edges()})
            print '-----'
            print 'state:', state


        if state == 1: # Entry point to algorithm
            # get alpha to point to cells of nu
            j = 1
            alpha[0] = 0
            for i from 0 < i < n:
                if nu.levels[i-1] == 0:
                    alpha[j] = i
                    j += 1
            alpha[j] = -1

            # "nu[0] := R(G, Pi, Pi)"
            nu._refine(alpha, n, M, _dig, uif)

            if not _dig:
                if nu._sat_225(n): hh = nu.k
            if nu._is_discrete(): state = 18; continue

            # store the first smallest nontrivial cell in W[k], and set v[k]
            # equal to its minimum element
            v[nu.k] = nu._first_smallest_nontrivial(W[nu.k], n)
            mpz_set_ui(Lambda_mpz[nu.k], 0)
            e[nu.k] = 0 # see state 12, and 17
            state = 2

        elif state == 2: # Move down the search tree one level by refining nu
            nu.k += 1

            # "nu[k] := nu[k-1] perp v[k-1]"
            nu._clear()
            alpha[0] = nu._split_vertex(v[nu.k-1])
            alpha[1] = -1
            i = nu._refine(alpha, n, M, _dig, uif)
            if verbosity >= 5:
                ST_vis_invariant = int(i)
                ST_vis_splitvert = int(v[nu.k-1])

            mpz_set_si(Lambda_mpz[nu.k], i)

            # only if this is the first time moving down the search tree:
            if h == -1: state = 5; continue

            # update hzf
            if hzf == nu.k-1 and mpz_cmp(Lambda_mpz[nu.k], zf_mpz[nu.k]) == 0: hzf = nu.k
            if not lab: state = 3; continue

            # "qzb := cmp(Lambda[k], zb[k])"
            if qzb == 0:
                if mpz_cmp_si(zb_mpz[nu.k], -1) == 0: # if "zb[k] == oo"
                    qzb = -1
                else:
                    qzb = mpz_cmp( Lambda_mpz[nu.k], zb_mpz[nu.k] )
            # update hzb
            if hzb == nu.k-1 and qzb == 0: hzb = nu.k

            # if Lambda[k] > zb[k], then zb[k] := Lambda[k]
            # (zb keeps track of the indicator invariants corresponding to
            # rho, the closest canonical leaf so far seen- if Lambda is
            # bigger, then rho must be about to change
            if qzb > 0: mpz_set(zb_mpz[nu.k], Lambda_mpz[nu.k])
            state = 3

        elif state == 3: # attempt to rule out automorphisms while moving down
                         # the tree
            if hzf <= nu.k or (lab and qzb >= 0): # changed hzb to hzf, == to <=
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
            if nu._is_discrete(): state = 7; continue

            # store the first smallest nontrivial cell in W[k], and set v[k]
            # equal to its minimum element
            v[nu.k] = nu._first_smallest_nontrivial(W[nu.k], n)

            if _dig or not nu._sat_225(n): hh = nu.k + 1
            e[nu.k] = 0 # see state 12, and 17
            state = 2 # continue down the tree

        elif state == 5: # alternative to 3: since we have not yet gotten
                         # zeta, there are no automorphisms to rule out.
                         # instead we record Lambda to zf and zb
                         # (see state 3)
            if _base:
                base_set.append(v[nu.k-1])
            mpz_set(zf_mpz[nu.k], Lambda_mpz[nu.k])
            mpz_set(zb_mpz[nu.k], Lambda_mpz[nu.k])
            state = 4

        elif state == 6: # at this stage, there is no reason to continue
                         # downward, and an automorphism has not been
                         # discovered
            j = nu.k

            # return to the longest ancestor nu[i] of nu that could have a
            # descendant equivalent to zeta or could improve on rho.
            # All terminal nodes descending from nu[hh] are known to be
            # equivalent, so i < hh. Also, if i > hzb, none of the
            # descendants of nu[i] can improve rho, since the indicator is
            # off (Lambda(nu) < Lambda(rho)). If i >= ht, then no descendant
            # of nu[i] is equivalent to zeta (see [1, p67]).
            if ht-1 > hzb:
                if ht-1 < hh-1:
                    nu.k = ht-1
                else:
                    nu.k = hh-1
            else:
                if hzb < hh-1:
                    nu.k = hzb
                else:
                    nu.k = hh-1

            # TODO: investigate the following line
            if nu.k == -1: nu.k = 0 # not in BDM, broke at G = Graph({0:[], 1:[]}), Pi = [[0,1]], lab=False

            if lab:
                if hb > nu.k: # update hb since we are backtracking (not in [1])
                    hb = nu.k # recall hb is the longest common ancestor of rho and nu

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
                if nu._is_min_cell_rep(i, hh):
                    Omega[l][i] = 1
                    if nu._is_fixed(i, hh):
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
            if nu.k < hzf: state = 8; continue

            # get the permutation corresponding to this terminal node
            nu._get_permutation_from(zeta, gamma)

            if verbosity > 3:
                print 'checking for automorphism:'
                print [gamma[iii] for iii in range(n)]

            # if G^gamma == G, goto 10
            if _is_automorphism(M, n, gamma):
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
            if (qzb > 0) or (nu.k < k_rho):
                state = 9; continue

            # if G(nu) > G(rho) (returns 1), goto 9
            # if G(nu) < G(rho) (returns -1), goto 6
            # if G(nu) == G(rho) (returns 0), get the automorphism and goto 10
            m = nu._compare_with(M, n, rho)

            if m > 0:
                state = 9; continue
            if m < 0:
                state = 6; continue

            rho._get_permutation_from(nu, gamma)
            if verbosity > 3:
                print 'automorphism discovered:'
                print [gamma[iii] for iii in range(n)]
            state = 10

        elif state == 9: # entering this state, nu is a best-so-far guess at
                         # the canonical label
            rho = PartitionStack(nu)
            k_rho = nu.k

            qzb = 0
            hb = nu.k
            hzb = nu.k

            # set zb[k+1] = Infinity
            mpz_set_si(zb_mpz[nu.k+1], -1)
            state = 6

        elif state == 10: # we have an automorphism to process
            # increment l
            if l < L - 1:
                l += 1

            for i from 0 <= i < n:
                if gamma[i] == i:
                    Phi[l][i] = 1
                    Omega[l][i] = 1
                else:
                    Phi[l][i] = 0
                    m = i
                    j = gamma[i]
                    while j != i:
                        if j < m: m = j
                        j = gamma[j]
                    if m == i:
                        Omega[l][i] = 1
                    else:
                        Omega[l][i] = 0
            m = Theta._incorporate_permutation(gamma, n)
            # if each orbit of gamma is part of an orbit in Theta, then the
            # automorphism is already in the span of those we have seen
            if not m:
                state = 11
                continue

            # record the automorphism
            output.append([ Integer(gamma[i]) for i from 0 <= i < n ])

            # The variable tvh represents the minimum element of W[k],
            # the last time we were at state 13 and backtracking up
            # zeta. If this is not still a minimal cell representative of Theta,
            # then we need to immediately backtrack to the place where it was
            # defined on a part of zeta, since the rest of the tree is now
            # equivalent. Otherwise, proceed to 11 and 12 before moving back to
            # the hub state.
            if Theta.elements[tvh] == -1:
                state = 11
                continue
            nu.k = h
            state = 13

        elif state == 11: # we have just discovered an automorphism,
                          # but tvh is still a minimal representative for its
                          # orbit in Theta. Therefore we cannot backtrack all
                          # the way to where zeta meets nu. Instead we just use
                          # indicator values to determine where to backtrack.
            if lab:
                nu.k = hb
            else:
                nu.k = h
            state = 12

        elif state == 12:
            # e keeps track of the whether W[k] has been thinned out by Omega
            # and Phi. It is set to 1 when you have just finished coming up the
            # search tree, and have intersected W[k] with Omega[i], for the
            # appropriate i < l, but since there may be an automorphism mapping
            # one element of W[k] to another still, we thin out W[k] again.
            # (see state 17) Coming from 11, this is an explicit automorphism.
            # Coming from 6, this is an implicit automorphism.
            if e[nu.k] == 1:
                for j from 0 <= j < n:
                    if W[nu.k][j] and not Omega[l][j]:
                        W[nu.k][j] = 0
            state = 13

        elif state == 13: # hub state
            if nu.k == -1:
                # the algorithm has finished
                state = -1; continue
            if nu.k > h:
                # if we are not at a node of zeta
                state = 17; continue
            if nu.k == h:
                # if we are at a node of zeta, then we have not yet backtracked
                # UP zeta, so skip the rest of 13
                state = 14; continue

            # thus, it must be that k < h, and this means we are done
            # searching underneath zeta[k+1], so now, k is the new longest
            # ancestor of nu and zeta:
            h = nu.k

            # set tvh to the minimum cell representative of W[k]
            # (see states 10 and 14)
            for i from 0 <= i < n:
                if W[nu.k][i]:
                    tvh = i
                    break
            state = 14

        elif state == 14: # iterate v[k] through W[k] until a minimum cell rep
                          # of Theta is found:
                          # this state gets hit only when we are looking for a
                          # new split off of zeta
            # The variable tvh was set to be the minimum element of W[k]
            # the last time we were at state 13 and backtracking up
            # zeta. If this is in the same cell of Theta as v[k], increment
            # index (see Theorem 2.33 in [1])
            if Theta._find(v[nu.k]) == Theta._find(tvh):
                index += 1

            # find the next v[k] in W[k]
            i = v[nu.k] + 1
            while i < n and not W[nu.k][i]:
                i += 1
            if i < n:
                v[nu.k] = i
            else:
                # there is no new vertex to consider at this level
                v[nu.k] = -1
                state = 16
                continue

            # if the new v[k] is not a minimum cell representative of Theta,
            # then we already considered that rep., and that subtree was
            # isomorphic to the one corresponding to v[k]
            if Theta.elements[v[nu.k]] != -1: state = 14
            else:
                # otherwise, we do have a vertex to consider
                state = 15

        elif state == 15: # we have a new vertex, v[k], that we must split on
            # hh is smallest such that nu[hh] satisfies Lemma 2.25. If it is
            # larger than k+1, it must be modified, since we are changing that
            # part
            if nu.k + 1 < hh:
                hh = nu.k + 1
            # hzf is maximal such that indicators line up for nu and zeta
            if nu.k < hzf:
                hzf = nu.k
            if not lab or hzb < nu.k:
                # in either case there is no need to update hzb, which is the
                # length for which nu and rho have the same indicators
                state = 2; continue
            hzb = nu.k
            qzb = 0
            state = 2

        elif state == 16: # backtrack one level in the search tree, recording
                          # information relevant to Theorem 2.33
            j = 0
            for i from 0 <= i < n:
                if W[nu.k][i]: j += 1
            if j == index and ht == nu.k+1: ht = nu.k
            size *= index
            index = 0
            nu.k -= 1

            if lab:
                if hb > nu.k: # update hb since we are backtracking (not in [1]):
                    hb = nu.k # recall hb is the longest common ancestor of rho and nu

            state = 13

        elif state == 17: # you have just finished coming up the search tree,
                          # and must now consider going back down.
            if e[nu.k] == 0:
                # intersect W[k] with each Omega[i] such that {v_0,...,v_(k-1)}
                # is contained in Phi[i]
                for i from 0 <= i <= l:
                    # check if {v_0,...,v_(k-1)} is contained in Phi[i]
                    # i.e. fixed pointwise by the automorphisms so far seen
                    j = 0
                    while j < nu.k and Phi[i][v[j]]:
                        j += 1
                    # if so, only check the minimal orbit representatives
                    if j == nu.k:
                        for j from 0 <= j < n:
                            if W[nu.k][j] and not Omega[i][j]:
                                W[nu.k][j] = 0
            e[nu.k] = 1 # see state 12

            # see if there is a relevant vertex to split on:
            i = v[nu.k] + 1
            while i < n and not W[nu.k][i]:
                i += 1
            if i < n:
                v[nu.k] = i
                state = 15
                continue
            else:
                v[nu.k] = -1

            # otherwise backtrack one level
            nu.k -= 1
            state = 13

        elif state == 18: # The first time we encounter a terminal node, we
                          # come straight here to set up zeta. This is a one-
                          # time state.
            # initialize counters for zeta:
            h = nu.k # zeta[h] == nu[h]
            ht = nu.k # nodes descended from zeta[ht] are all equivalent
            hzf = nu.k # max such that indicators for zeta and nu agree

            zeta = PartitionStack(nu)

            nu.k -= 1
            if not lab: state = 13; continue

            rho = PartitionStack(nu)

            # initialize counters for rho:
            k_rho = nu.k + 1 # number of partitions in rho
            hzb = nu.k # max such that indicators for rho and nu agree - BDM had k+1
            hb = nu.k # rho[hb] == nu[hb] - BDM had k+1

            qzb = 0 # Lambda[k] == zb[k], so...
            state = 13

    # free the GMP ints
    for i from 0 <= i < n+2:
        mpz_clear(Lambda_mpz[i])
        mpz_clear(zf_mpz[i])
        mpz_clear(zb_mpz[i])

    # free int arrays
    sage_free(gamma)
    sage_free(W[0])
    sage_free(Phi[0])
    sage_free(Omega[0])
    sage_free(alpha)
    sage_free(v)
    sage_free(e)

    # free GMP int pointers
    sage_free(Lambda_mpz)
    sage_free(zf_mpz)
    sage_free(zb_mpz)

    # free int pointers
    sage_free(W)
    sage_free(Phi)
    sage_free(Omega)

    # prepare output
    if not (lab or dict_rep or certify or base or order):
        return output
    output_tuple = [output]
    if dict_rep:
        if isinstance(G, GenericGraph):
            ddd = {}
            for vvv in ffrom.iterkeys(): # v is a C variable
                if ffrom[vvv] != 0:
                    ddd[vvv] = ffrom[vvv]
                else:
                    ddd[vvv] = n
        else: # G is a CGraph
            ddd = {}
            for i from 0 <= i < n:
                ddd[i] = i
        output_tuple.append(ddd)
    if lab:
        H = _term_pnest_graph(G, rho)
        output_tuple.append(H)
    if certify:
        dd = {}
        for i from 0 <= i < n:
            dd[to[rho.entries[i]]] = i
        output_tuple.append(dd)
    if base:
        output_tuple.append(base_set)
    if order:
        output_tuple.append(size)
    return tuple(output_tuple)

# Benchmarking functions

def all_labeled_graphs(n):
    """
    Returns all labeled graphs on n vertices {0,1,...,n-1}. Used in
    classifying isomorphism types (naive approach), and more importantly
    in benchmarking the search algorithm.

    EXAMPLE:
        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree
        sage: from sage.graphs.graph_isom import all_labeled_graphs
        sage: Glist = {}
        sage: Giso  = {}
        sage: for n in range(1,5):
        ...    Glist[n] = all_labeled_graphs(n)
        ...    Giso[n] = []
        ...    for g in Glist[n]:
        ...        a, b = search_tree(g, [range(n)])
        ...        inn = False
        ...        for gi in Giso[n]:
        ...            if b == gi:
        ...                inn = True
        ...        if not inn:
        ...            Giso[n].append(b)
        sage: for n in Giso:
        ...    print n, len(Giso[n])
        1 1
        2 2
        3 4
        4 11
        sage: n = 5
        sage: Glist[n] = all_labeled_graphs(n)
        sage: Giso[n] = []
        sage: for g in Glist[5]:
        ...    a, b = search_tree(g, [range(n)])
        ...    inn = False
        ...    for gi in Giso[n]:
        ...        if b == gi:
        ...            inn = True
        ...    if not inn:
        ...        Giso[n].append(b)
        sage: print n, len(Giso[n]) # long time
        5 34
        sage: graphs_list.show_graphs(Giso[4])
    """
    TE = []
    for i in range(n):
        for j in range(i):
            TE.append((i, j))
    m = len(TE)
    Glist= []
    for i in range(2**m):
        G = Graph(n)
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

    EXAMPLE:
        sage: from sage.graphs.graph_isom import kpow
        sage: kpow(['a', 1, {}], 2)
        [[1, 'a'], [{}, 'a'], ['a', 1], [{}, 1], ['a', {}], [1, {}]]

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

    EXAMPLE:
        sage: from sage.graphs.graph_isom import all_ordered_partitions
        sage: all_ordered_partitions(['a', 1, {}])
        [[['a'], [1], [{}]],
         [['a'], [{}], [1]],
         [['a'], [{}, 1]],
         [['a'], [1, {}]],
         [[1], ['a'], [{}]],
         [[1], [{}], ['a']],
         [[1], [{}, 'a']],
         [[1], ['a', {}]],
         [[{}], ['a'], [1]],
         [[{}], [1], ['a']],
         [[{}], [1, 'a']],
         [[{}], ['a', 1]],
         [[1, 'a'], [{}]],
         [[{}, 'a'], [1]],
         [['a', 1], [{}]],
         [[{}, 1], ['a']],
         [['a', {}], [1]],
         [[1, {}], ['a']],
         [[{}, 1, 'a']],
         [[1, {}, 'a']],
         [[{}, 'a', 1]],
         [['a', {}, 1]],
         [[1, 'a', {}]],
         [['a', 1, {}]]]

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
        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree
        sage: from sage.graphs.graph_isom import all_labeled_digraphs_with_loops
        sage: Glist = {}
        sage: Giso  = {}
        sage: for n in range(1,4):
        ...    Glist[n] = all_labeled_digraphs_with_loops(n)
        ...    Giso[n] = []
        ...    for g in Glist[n]:
        ...        a, b = search_tree(g, [range(n)], dig=True)
        ...        inn = False
        ...        for gi in Giso[n]:
        ...            if b == gi:
        ...                inn = True
        ...        if not inn:
        ...            Giso[n].append(b)
        sage: for n in Giso:
        ...    print n, len(Giso[n])
        1 2
        2 10
        3 104
    """
    TE = []
    for i in range(n):
        for j in range(n):
            TE.append((i, j))
    m = len(TE)
    Glist= []
    for i in range(2**m):
        G = DiGraph(n, loops=True)
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
        sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree
        sage: from sage.graphs.graph_isom import all_labeled_digraphs
        sage: Glist = {}
        sage: Giso  = {}
        sage: for n in range(1,4):
        ...       Glist[n] = all_labeled_digraphs(n)
        ...       Giso[n] = []
        ...       for g in Glist[n]:
        ...           a, b = search_tree(g, [range(n)], dig=True)
        ...           inn = False
        ...           for gi in Giso[n]:
        ...               if b == gi:
        ...                   inn = True
        ...           if not inn:
        ...               Giso[n].append(b)
        sage: for n in Giso:
        ...       print n, len(Giso[n])
        1 1
        2 3
        3 16
    """

##         sage: n = 4 # long time (4 minutes)
##         sage: Glist[n] = all_labeled_digraphs(n) # long time
##         sage: Giso[n] = [] # long time
##         sage: for g in Glist[n]: # long time
##         ...       a, b = search_tree(g, [range(n)], dig=True)
##         ...       inn = False
##         ...       for gi in Giso[n]:
##         ...           if enum(b) == enum(gi):
##         ...               inn = True
##         ...       if not inn:
##         ...           Giso[n].append(b)
##         sage: print n, len(Giso[n]) # long time
##         4 218
    TE = []
    for i in range(n):
        for j in range(n):
            if i != j: TE.append((i, j))
    m = len(TE)
    Glist= []
    for i in range(2**m):
        G = DiGraph(n, loops=True)
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

    EXAMPLE:
        sage: from sage.graphs.graph_isom import perm_group_elt
        sage: perm_group_elt([0,2,1])
        (1,2)
        sage: perm_group_elt([1,2,0])
        (1,2,3)

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
        list_perm -- if True, assumes \var{gamma} is a list representing the map
    $i \mapsto \var{gamma}[i]$.

    EXAMPLES:
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

def verify_partition_refinement(G, initial_partition, terminal_partition):
    """
    Verify that the refinement is correct.

    EXAMPLE:
        sage: from sage.graphs.graph_isom import PartitionStack
        sage: from sage.graphs.base.sparse_graph import SparseGraph
        sage: G = SparseGraph(10)
        sage: for i,j,_ in graphs.PetersenGraph().edge_iterator():
        ...    G.add_arc(i,j)
        ...    G.add_arc(j,i)
        sage: P = PartitionStack(10)
        sage: P.set_k(1)
        sage: P.split_vertex(0)
        sage: P.refine(G, [0], 10, 0, 1)
        sage: P
        (0,2,3,6,7,8,9,1,4,5)
        (0|2,3,6,7,8,9|1,4,5)
        sage: P.set_k(2)
        sage: P.split_vertex(1)

    Note that this line implicitly tests the function verify_partition_refinement:
        sage: P.refine(G, [7], 10, 0, 1, test=True)
        sage: P
        (0,3,7,8,9,2,6,1,4,5)
        (0|3,7,8,9,2,6|1,4,5)
        (0|3,7,8,9|2,6|1|4,5)

    """
    if not G.is_equitable(terminal_partition):
        raise RuntimeError("Resulting partition is not equitable!!!!!!!!!\n"+\
        str(initial_partition) + "\n" + \
        str(terminal_partition) + "\n" + \
        str(G.am()))

