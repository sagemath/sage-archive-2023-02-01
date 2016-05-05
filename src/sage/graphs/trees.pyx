r"""
Generation of trees

This is an implementation of the algorithm for generating trees with `n` vertices
(up to isomorphism) in constant time per tree described in [WRIGHT-ETAL]_.

AUTHORS:

- Ryan Dingman (2009-04-16): initial version

REFERENCES:

    .. [WRIGHT-ETAL] Wright, Robert Alan; Richmond, Bruce; Odlyzko, Andrew; McKay, Brendan D.
       Constant time generation of free trees. SIAM J. Comput. 15 (1986), no. 2,
       540--548.
"""

cdef extern from "limits.h":
    cdef int INT_MAX

include "cysignals/memory.pxi"

# from networkx import MultiGraph

from sage.graphs.graph import Graph
from sage.graphs.base.sparse_graph cimport SparseGraph
from sage.graphs.base.sparse_graph cimport SparseGraphBackend

cdef class TreeIterator:
    r"""
    This class iterates over all trees with n vertices (up to isomorphism).

    EXAMPLES::

        sage: from sage.graphs.trees import TreeIterator
        sage: def check_trees(n):
        ...       trees = []
        ...       for t in TreeIterator(n):
        ...           if t.is_tree() == False:
        ...               return False
        ...           if t.num_verts() != n:
        ...               return False
        ...           if t.num_edges() != n - 1:
        ...               return False
        ...           for tree in trees:
        ...               if tree.is_isomorphic(t) == True:
        ...                   return False
        ...           trees.append(t)
        ...       return True
        sage: print check_trees(10)
        True

    ::

        sage: from sage.graphs.trees import TreeIterator
        sage: count = 0
        sage: for t in TreeIterator(15):
        ...       count += 1
        sage: print count
        7741
    """

    def __init__(self, int vertices):
        r"""
        Initializes an iterator over all trees with `n` vertices.

        EXAMPLES::

            sage: from sage.graphs.trees import TreeIterator
            sage: t = TreeIterator(100) # indirect doctest
            sage: print t
            Iterator over all trees with 100 vertices
        """
        self.vertices = vertices
        self.l = NULL
        self.current_level_sequence = NULL
        self.first_time = 1

    def __dealloc__(self):
        r"""
        EXAMPLES::

            sage: from sage.graphs.trees import TreeIterator
            sage: t = TreeIterator(100)
            sage: t = None # indirect doctest
        """
        if self.l != NULL:
            sig_free(self.l)
            self.l = NULL
        if self.current_level_sequence != NULL:
            sig_free(self.current_level_sequence)
            self.current_level_sequence = NULL

    def __str__(self):
        r"""
        EXAMPLES::

            sage: from sage.graphs.trees import TreeIterator
            sage: t = TreeIterator(100)
            sage: print t # indirect doctest
            Iterator over all trees with 100 vertices
        """
        return "Iterator over all trees with %s vertices"%(self.vertices)

    def __iter__(self):
        r"""
        Returns an iterator over all the trees with `n` vertices.

        EXAMPLES::

            sage: from sage.graphs.trees import TreeIterator
            sage: t = TreeIterator(4)
            sage: list(iter(t))
            [Graph on 4 vertices, Graph on 4 vertices]
        """
        return self

    def __next__(self):
        r"""
        Returns the next tree with `n` vertices

        EXAMPLES::

            sage: from sage.graphs.trees import TreeIterator
            sage: T = TreeIterator(5)
            sage: [t for t in T] # indirect doctest
            [Graph on 5 vertices, Graph on 5 vertices, Graph on 5 vertices]


        TESTS:

        This used to be broken for trees with no vertices
        and was fixed in :trac:`13719` ::

            sage: from sage.graphs.trees import TreeIterator
            sage: T = TreeIterator(0)
            sage: [t for t in T] # indirect doctest
            [Graph on 0 vertices]
        """

        if not self.first_time and self.q == 0:
            raise StopIteration

        if self.first_time == 1:
            if self.vertices == 0:
                self.first_time = 0
                self.q = 0
            else:
                self.l = <int *>sig_malloc(self.vertices * sizeof(int))
                self.current_level_sequence = <int *>sig_malloc(self.vertices * sizeof(int))

                if self.l == NULL or self.current_level_sequence == NULL:
                    raise MemoryError

                self.generate_first_level_sequence()
                self.first_time = 0
        else:
            self.generate_next_level_sequence()

        cdef int i
        cdef int vertex1
        cdef int vertex2
        cdef object G

#        from networkx import MultiGraph
#        G = Graph(self.vertices)
#        cdef object XG = G._backend._nxg
#
#        for i from 2 <= i <= self.vertices:
#            vertex1 = i - 1
#            vertex2 = self.current_level_sequence[i - 1] - 1
#            XG.add_edge(vertex1, vertex2)
#
#        return G

        # Currently, c_graph does not have all the same functionality as networkx.
        # Until it does, we can't generate graphs using the c_graph backend even
        # though it is twice as fast (for our purposes) as networkx.

        G = Graph(self.vertices, implementation='c_graph', sparse=True)
        cdef SparseGraph SG = (<SparseGraphBackend?> G._backend)._cg

        for i from 2 <= i <= self.vertices:
            vertex1 = i - 1
            vertex2 = self.current_level_sequence[i - 1] - 1
            SG.add_arc_unsafe(vertex1, vertex2)
            SG.add_arc_unsafe(vertex2, vertex1)

        return G

    cdef int generate_first_level_sequence(self):
        r"""
        Generates the level sequence representing the first tree with `n` vertices
        """
        cdef int i
        cdef int k

        k = (self.vertices / 2) + 1

        if self.vertices == 4:
            self.p = 3
        else:
            self.p = self.vertices
        self.q = self.vertices - 1
        self.h1 = k
        self.h2 = self.vertices
        if self.vertices % 2 == 0:
            self.c = self.vertices + 1
        else:
            self.c = INT_MAX # oo

        self.r = k

        for i from 1 <= i <= k:
            self.l[i - 1] = i
        for i from k < i <= self.vertices:
            self.l[i - 1] = i - k + 1
        for i from 0 <= i < self.vertices:
            self.current_level_sequence[i] = i
        if self.vertices > 2:
            self.current_level_sequence[k] = 1
        if self.vertices <= 3:
            self.q = 0

        return 0

    cdef int generate_next_level_sequence(self):
        r"""
        Generates the level sequence representing the next tree with `n` vertices
        """
        cdef int i
        cdef int fixit = 0

        cdef int needr = 0
        cdef int needc = 0
        cdef int needh2 = 0

        cdef int n = self.vertices
        cdef int p = self.p
        cdef int q = self.q
        cdef int h1 = self.h1
        cdef int h2 = self.h2
        cdef int c = self.c
        cdef int r = self.r
        cdef int *l = self.l
        cdef int *w = self.current_level_sequence

        if c == n + 1 or p == h2 and (l[h1 - 1] == l[h2 - 1] + 1 and n - h2 > r - h1 or l[h1 - 1] == l[h2 - 1] and n - h2 + 1 < r - h1):
            if (l[r - 1] > 3):
                p = r
                q = w[r - 1]
                if h1 == r:
                    h1 = h1 - 1
                fixit = 1
            else:
                p = r
                r = r - 1
                q = 2

        if p <= h1:
            h1 = p - 1
        if p <= r:
           needr = 1
        elif p <= h2:
           needh2 = 1
        elif l[h2 - 1] == l[h1 - 1] - 1 and n - h2 == r - h1:
            if p <= c:
                needc = 1
        else:
            c = INT_MAX

        cdef int oldp = p
        cdef int delta = q - p
        cdef int oldlq = l[q - 1]
        cdef int oldwq = w[q - 1]
        p = INT_MAX

        for i from oldp <= i <= n:
            l[i - 1] = l[i - 1 + delta]
            if l[i - 1] == 2:
                w[i - 1] = 1
            else:
                p = i
                if l[i - 1] == oldlq:
                    q = oldwq
                else:
                    q = w[i - 1 + delta] - delta
                w[i - 1] = q
            if needr == 1 and l[i - 1] == 2:
                needr = 0
                needh2 = 1
                r = i - 1
            if needh2 == 1 and l[i - 1] <= l[i - 2] and i > r + 1:
                needh2 = 0
                h2 = i - 1
                if l[h2 - 1] == l[h1 - 1] - 1 and n - h2 == r - h1:
                    needc = 1
                else:
                    c = INT_MAX
            if needc == 1:
                if l[i - 1] != l[h1 - h2 + i - 1] - 1:
                    needc = 0
                    c = i
                else:
                    c = i + 1

        if fixit == 1:
            r = n - h1 + 1
            for i from r < i <= n:
                l[i - 1] = i - r + 1
                w[i - 1] = i - 1
            w[r] = 1
            h2 = n
            p = n
            q = p - 1
            c = INT_MAX
        else:
            if p == INT_MAX:
                if l[oldp - 2] != 2:
                    p = oldp - 1
                else:
                    p = oldp - 2
                q = w[p - 1]
            if needh2 == 1:
                h2 = n
                if l[h2 - 1] == l[h1 - 1] - 1 and h1 == r:
                    c = n + 1
                else:
                    c = INT_MAX

        self.p = p
        self.q = q
        self.h1 = h1
        self.h2 = h2
        self.c = c
        self.r = r
        self.l = l
        self.current_level_sequence = w

        return 0
