from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement
from collections import deque

cdef class Iterator(object):
    """
    Iterator class for reflection groups.
    """
    cdef int n
    cdef list S
    cdef list is_positive_root
    cdef str algorithm
    cdef bint tracking_words

    def __init__(self, W, algorithm="depth", tracking_words=True):
        """
        Initalize ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.reflection_group_c import Iterator
            sage: W = ReflectionGroup(["B", 4])
            sage: I = Iterator(W)
            sage: TestSuite(I).run(skip="_test_pickling")
        """
        self.S = W.gens()
        self.n = len(W._index_set)
        self.is_positive_root = W._is_positive_root
        self.tracking_words = tracking_words

        # "breadth" is 1.5 - 2x slower than "depth" since it uses
        # a deque with popleft instead of a list with pop
        if algorithm not in ["depth","breadth"]:
            raise ValueError('The algorithm (="%s") must be either "depth" or "breadth"')
        self.algorithm = algorithm

    cdef list succ(self, PermutationGroupElement u, int first):
        cdef PermutationGroupElement u1, si
        cdef int i, j
        cdef list successors = []

        for i in range(first):
            si = <PermutationGroupElement>(self.S[i])
            u1 = <PermutationGroupElement>(si._mul_(u))
            if self.test(u1, i):
                successors.append((u1, i))
        for i in range(first+1, self.n):
            if self.is_positive_root[u.perm[i]+1]:
                si = <PermutationGroupElement>(self.S[i])
                u1 = <PermutationGroupElement>(si._mul_(u))
                if self.test(u1, i):
                    successors.append((u1, i))
        return successors

    cdef list succ_words(self, PermutationGroupElement u, list word, int first):
        cdef PermutationGroupElement u1, si
        cdef int i, j
        cdef list successors = []
        cdef list word_new

        for i in range(first):
            si = <PermutationGroupElement>(self.S[i])
            u1 = <PermutationGroupElement>(si._mul_(u))
            if self.test(u1, i):
                # try to use word+[i] and the reversed
                word_new = [i]+word
                u1._reduced_word = word_new
                successors.append((u1, word_new, i))
        for i in range(first+1, self.n):
            if self.is_positive_root[u.perm[i]+1]:
                si = <PermutationGroupElement>(self.S[i])
                u1 = <PermutationGroupElement>(si._mul_(u))
                if self.test(u1, i):
                    word_new = [i]+word
                    u1._reduced_word = word_new
                    successors.append((u1, word_new, i))
        return successors

    cdef bint test(self, PermutationGroupElement u1, int i):
        cdef int j

        for j in range(i):
            if not self.is_positive_root[u1.perm[j]+1]:
                return False
        return True

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.root_system.reflection_group_c import Iterator
            sage: W = ReflectionGroup(["B", 4])
            sage: I = Iterator(W)
            sage: len(list(I)) == W.cardinality()
            True
        """
        # the breadth search iterator is ~2x slower as it
        # uses a deque with popleft 
        if self.algorithm == "depth":
            if self.tracking_words:
                return self.iter_words_depth()
            else:
                return self.iter_depth()
        elif self.algorithm == "breadth":
            if self.tracking_words:
                return self.iter_words_breadth()
            else:
                return self.iter_breadth()

    def iter_depth(self):
        """
        Iterate over ``self`` using depth-first-search.

        EXAMPLES::

            sage: from sage.combinat.root_system.reflection_group_c import Iterator
            sage: W = ReflectionGroup(['B', 2])
            sage: I = Iterator(W)
            sage: list(I.iter_depth())
            [(),
             (1,3)(2,6)(5,7),
             (1,5)(2,4)(6,8),
             (1,3,5,7)(2,8,6,4),
             (1,7)(3,5)(4,8),
             (1,7,5,3)(2,4,6,8),
             (2,8)(3,7)(4,6),
             (1,5)(2,6)(3,7)(4,8)]
        """
        cdef tuple node
        cdef list cur = [(self.S[0].parent().one(), -1)]
        cdef PermutationGroupElement u
        cdef int first
        cdef list L = []

        while True:
            if not cur:
                if not L:
                    return
                cur = L.pop()
                continue

            u, first = cur.pop()
            yield u
            L.append(self.succ(u, first))

    def iter_words_depth(self):
        """
        Iterate over ``self`` using depth-first-search and setting
        the reduced word.

        EXAMPLES::

            sage: from sage.combinat.root_system.reflection_group_c import Iterator
            sage: W = ReflectionGroup(['B', 2])
            sage: I = Iterator(W)
            sage: for w in I.iter_words_depth(): w._reduced_word
            []
            [1]
            [0]
            [1, 0]
            [0, 1, 0]
            [0, 1]
            [1, 0, 1]
            [0, 1, 0, 1]
        """
        cdef tuple node
        cdef list cur, word

        cdef PermutationGroupElement u
        cdef int first
        cdef list L = []

        one = self.S[0].parent().one()
        one._reduced_word = []
        cur = [(one, list(), -1)]

        while True:
            if not cur:
                if not L:
                    return
                cur = L.pop()
                continue

            u, word, first = cur.pop()
            yield u
            L.append(self.succ_words(u, word, first))

    def iter_breadth(self):
        """
        Iterate over ``self`` using breadth-first-search.

        EXAMPLES::

            sage: from sage.combinat.root_system.reflection_group_c import Iterator
            sage: W = ReflectionGroup(['B', 2])
            sage: I = Iterator(W)
            sage: list(I.iter_breadth())
            [(),
             (1,3)(2,6)(5,7),
             (1,5)(2,4)(6,8),
             (1,7,5,3)(2,4,6,8),
             (1,3,5,7)(2,8,6,4),
             (2,8)(3,7)(4,6),
             (1,7)(3,5)(4,8),
             (1,5)(2,6)(3,7)(4,8)]
        """
        cdef tuple node
        cdef list cur = [(self.S[0].parent().one(), -1)]
        cdef PermutationGroupElement u
        cdef int first
        L = deque()

        while True:
            if not cur:
                if not L:
                    return
                cur = L.popleft()
                continue

            u, first = cur.pop()
            yield u
            L.append(self.succ(u, first))

    def iter_words_breadth(self):
        """
        Iterate over ``self`` using breadth-first-search and setting
        the reduced word.

        EXAMPLES::

            sage: from sage.combinat.root_system.reflection_group_c import Iterator
            sage: W = ReflectionGroup(['B', 2])
            sage: I = Iterator(W)
            sage: for w in I.iter_words_breadth(): w._reduced_word
            []
            [1]
            [0]
            [0, 1]
            [1, 0]
            [1, 0, 1]
            [0, 1, 0]
            [0, 1, 0, 1]
        """
        cdef tuple node
        cdef list cur, word
        cdef PermutationGroupElement u
        cdef int first
        L = deque()

        one = self.S[0].parent().one()
        one._reduced_word = []
        cur = [(one, list(), -1)]

        while True:
            if not cur:
                if not L:
                    return
                cur = L.popleft()
                continue

            u, word, first = cur.pop()
            yield u
            L.append(self.succ_words(u, word, first))

