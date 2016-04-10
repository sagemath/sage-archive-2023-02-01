r"""
This contains a few time-critial auxillary cython functions for
finite complex or real reflection groups.
"""
#*****************************************************************************
#       Copyright (C) 2011-2016 Christian Stump <christian.stump at gmail.com>
#                     2016 Travis Scrimshaw
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement
from collections import deque

cdef class Iterator(object):
    """
    Iterator class for reflection groups.
    """
    cdef int n
    cdef int N # number of refections/positive roots
    cdef tuple S
    cdef str algorithm
    cdef bint tracking_words

    def __init__(self, W, int N, str algorithm="depth", bint tracking_words=True):
        """
        Initalize ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.reflection_group_c import Iterator
            sage: W = ReflectionGroup(["B", 4])
            sage: I = Iterator(W, W.nr_reflections())
            sage: TestSuite(I).run(skip="_test_pickling")
        """
        self.S = tuple( W.simple_reflection(W._index_set_inverse[i])
                        for i in range(len(W._index_set)) )
        self.n = len(W._index_set)
        self.N = N
        self.tracking_words = tracking_words

        # "breadth" is 1.5x slower than "depth" since it uses
        # a deque with popleft instead of a list with pop
        if algorithm not in ["depth","breadth"]:
            raise ValueError('the algorithm (="%s") must be either "depth" or "breadth"')
        self.algorithm = algorithm

    cdef list succ(self, PermutationGroupElement u, int first):
        cdef PermutationGroupElement u1, si
        cdef int i, j
        cdef list successors = []
        cdef tuple S = self.S
        cdef int N = self.N

        for i in range(first):
            si = <PermutationGroupElement>(S[i])
            if self.test(u, si, i):
                successors.append((si._mul_(u), i))
        for i in range(first+1, self.n):
            if u.perm[i] < N:
                si = <PermutationGroupElement>(S[i])
                if self.test(u, si, i):
                    successors.append((si._mul_(u), i))
        return successors

    cdef list succ_words(self, PermutationGroupElement u, list word, int first):
        cdef PermutationGroupElement u1, si
        cdef int i, j
        cdef list successors = []
        cdef list word_new
        cdef tuple S = self.S
        cdef int N = self.N

        for i in range(first):
            si = <PermutationGroupElement>(S[i])
            if self.test(u, si, i):
                u1 = <PermutationGroupElement>(si._mul_(u))
                # try to use word+[i] and the reversed
                word_new = [i] + word
                u1._reduced_word = word_new
                successors.append((u1, word_new, i))
        for i in range(first+1, self.n):
            if u.perm[i] < self.N:
                si = <PermutationGroupElement>(S[i])
                if self.test(u, si, i):
                    u1 = <PermutationGroupElement>(si._mul_(u))
                    word_new = [i] + word
                    u1._reduced_word = word_new
                    successors.append((u1, word_new, i))
        return successors

    cdef bint test(self, PermutationGroupElement u, PermutationGroupElement si, int i):
        cdef int j

        for j in range(i):
            if u.perm[si.perm[j]] >= self.N:
                return False
        return True

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.root_system.reflection_group_c import Iterator
            sage: W = ReflectionGroup(["B", 4])
            sage: I = Iterator(W, W.nr_reflections())
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
            sage: I = Iterator(W, W.nr_reflections())
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
            sage: I = Iterator(W, W.nr_reflections())
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
            sage: I = Iterator(W, W.nr_reflections())
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
            sage: I = Iterator(W, W.nr_reflections())
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

