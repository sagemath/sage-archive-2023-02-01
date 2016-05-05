r"""
This contains a few time-critial auxillary cython functions for
finite complex or real reflection groups.
"""
#*****************************************************************************
#       Copyright (C) 2011-2016 Christian Stump <christian.stump at gmail.com>
#                     2016 Travis Scrimshaw <tscrimsh at umn.edu>
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
    cdef list noncom

    cdef list noncom_letters(self):
        """
        Return a list ``L`` of lists such that ...

        WARNING:: This is not used as it slows down the computation.

        EXAMPLES::

            sage: from sage.combinat.root_system.reflection_group_c import Iterator
            sage: W = ReflectionGroup(["B", 4])               # optional - gap3
            sage: I = Iterator(W, W.number_of_reflections())  # optional - gap3
            sage: TestSuite(I).run(skip="_test_pickling")     # optional - gap3
        """
        cdef tuple S = self.S
        cdef int n = len(S)
        cdef list noncom = []
        cdef list noncom_i
        for i in range(n):
            si = S[i]
            noncom_i = []
            for j in range(i+1,n):
                sj = S[j]
                if si._mul_(sj) == sj._mul_(si):
                    pass
                else:
                    noncom_i.append(j)
            noncom.append(noncom_i)
        noncom.append(range(n))
        return noncom

    def __init__(self, W, int N, str algorithm="depth", bint tracking_words=True):
        """
        Initalize ``self``.

        EXAMPLES::

            sage: from sage.combinat.root_system.reflection_group_c import Iterator
            sage: W = ReflectionGroup(["B", 4])               # optional - gap3
            sage: I = Iterator(W, W.number_of_reflections())  # optional - gap3
            sage: TestSuite(I).run(skip="_test_pickling")     # optional - gap3
        """
        self.S = tuple(W.simple_reflections())
        self.n = len(W._index_set)
        self.N = N
        self.tracking_words = tracking_words

        # "breadth" is 1.5x slower than "depth" since it uses
        # a deque with popleft instead of a list with pop
        if algorithm not in ["depth","breadth"]:
            raise ValueError('the algorithm (="%s") must be either "depth" or "breadth"')
        self.algorithm = algorithm

#        self.noncom = self.noncom_letters()

    cdef list succ(self, PermutationGroupElement u, int first):
        cdef PermutationGroupElement u1, si
        cdef int i
        cdef list successors = []
        cdef tuple S = self.S
        cdef int N = self.N
#        cdef list nc = self.noncom[first]

        for i in range(first):
            si = <PermutationGroupElement>(S[i])
            if self.test(u, si, i):
                successors.append((si._mul_(u), i))
        for i in range(first+1,self.n):
#        for i in nc:
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
        cdef int j, sij
        cdef int N = self.N
        cdef int* siperm = si.perm
        cdef int* uperm = u.perm

        for j in range(i):
            sij = siperm[j]
            if uperm[sij] >= N:
                return False
        return True

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.root_system.reflection_group_c import Iterator
            sage: W = ReflectionGroup(["B", 4])               # optional - gap3
            sage: I = Iterator(W, W.number_of_reflections())  # optional - gap3
            sage: len(list(I)) == W.cardinality()             # optional - gap3
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
            sage: W = ReflectionGroup(['B', 2])               # optional - gap3
            sage: I = Iterator(W, W.number_of_reflections())  # optional - gap3
            sage: list(I.iter_depth())                        # optional - gap3
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
            sage: W = ReflectionGroup(['B', 2])               # optional - gap3
            sage: I = Iterator(W, W.number_of_reflections())  # optional - gap3
            sage: for w in I.iter_words_depth(): w._reduced_word  # optional - gap3
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
            sage: W = ReflectionGroup(['B', 2])               # optional - gap3
            sage: I = Iterator(W, W.number_of_reflections())  # optional - gap3
            sage: list(I.iter_breadth())                      # optional - gap3
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
            sage: W = ReflectionGroup(['B', 2])               # optional - gap3
            sage: I = Iterator(W, W.number_of_reflections())  # optional - gap3
            sage: for w in I.iter_words_breadth(): w._reduced_word  # optional - gap3
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

def iterator_tracking_words(W):
    r"""
    Return an iterator through the elements of ``self`` together
    with the words in the simple generators.

    The iterator is a breadth first search through the graph of the
    elements of the group with generators.

    EXAMPLES::

        sage: from sage.combinat.root_system.reflection_group_c import iterator_tracking_words
        sage: W = ReflectionGroup(4)                        # optional - gap3
        sage: for w in iterator_tracking_words(W): w        # optional - gap3
        ((), [])
        ((1,3,9)(2,4,7)(5,10,18)(6,11,16)(8,12,19)(13,15,20)(14,17,21)(22,23,24), [0])
        ((1,5,13)(2,6,10)(3,7,14)(4,8,15)(9,16,22)(11,12,17)(18,19,23)(20,21,24), [1])
        ((1,9,3)(2,7,4)(5,18,10)(6,16,11)(8,19,12)(13,20,15)(14,21,17)(22,24,23), [0, 0])
        ((1,7,6,12,23,20)(2,8,17,24,9,5)(3,16,10,19,15,21)(4,14,11,22,18,13), [0, 1])
        ((1,10,4,12,21,22)(2,11,19,24,13,3)(5,15,7,17,16,23)(6,18,8,20,14,9), [1, 0])
        ((1,13,5)(2,10,6)(3,14,7)(4,15,8)(9,22,16)(11,17,12)(18,23,19)(20,24,21), [1, 1])
        ((1,16,12,15)(2,14,24,18)(3,5,19,17)(4,6,22,20)(7,8,23,9)(10,13,21,11), [0, 0, 1])
        ((1,2,12,24)(3,6,19,20)(4,17,22,5)(7,11,23,13)(8,21,9,10)(14,16,18,15), [0, 1, 0])
        ((1,14,12,18)(2,15,24,16)(3,22,19,4)(5,6,17,20)(7,10,23,21)(8,11,9,13), [0, 1, 1])
        ((1,18,12,14)(2,16,24,15)(3,4,19,22)(5,20,17,6)(7,21,23,10)(8,13,9,11), [1, 0, 0])
        ((1,15,12,16)(2,18,24,14)(3,17,19,5)(4,20,22,6)(7,9,23,8)(10,11,21,13), [1, 1, 0])
        ((1,6,23)(2,17,9)(3,10,15)(4,11,18)(5,8,24)(7,12,20)(13,14,22)(16,19,21), [0, 0, 1, 0])
        ((1,22,21,12,4,10)(2,3,13,24,19,11)(5,23,16,17,7,15)(6,9,14,20,8,18), [0, 0, 1, 1])
        ((1,4,21)(2,19,13)(3,11,24)(5,7,16)(6,8,14)(9,18,20)(10,12,22)(15,17,23), [0, 1, 0, 0])
        ((1,17,13,12,5,11)(2,20,10,24,6,21)(3,23,14,19,7,18)(4,9,15,22,8,16), [0, 1, 1, 0])
        ((1,19,9,12,3,8)(2,22,7,24,4,23)(5,21,18,17,10,14)(6,13,16,20,11,15), [1, 0, 0, 1])
        ((1,20,23,12,6,7)(2,5,9,24,17,8)(3,21,15,19,10,16)(4,13,18,22,11,14), [1, 1, 0, 0])
        ((1,11,5,12,13,17)(2,21,6,24,10,20)(3,18,7,19,14,23)(4,16,8,22,15,9), [0, 0, 1, 0, 0])
        ((1,23,6)(2,9,17)(3,15,10)(4,18,11)(5,24,8)(7,20,12)(13,22,14)(16,21,19), [0, 0, 1, 1, 0])
        ((1,8,3,12,9,19)(2,23,4,24,7,22)(5,14,10,17,18,21)(6,15,11,20,16,13), [0, 1, 0, 0, 1])
        ((1,21,4)(2,13,19)(3,24,11)(5,16,7)(6,14,8)(9,20,18)(10,22,12)(15,23,17), [0, 1, 1, 0, 0])
        ((1,12)(2,24)(3,19)(4,22)(5,17)(6,20)(7,23)(8,9)(10,21)(11,13)(14,18)(15,16), [0, 0, 1, 0, 0, 1])
        ((1,24,12,2)(3,20,19,6)(4,5,22,17)(7,13,23,11)(8,10,9,21)(14,15,18,16), [0, 0, 1, 1, 0, 0])
    """
    cdef tuple I = tuple(W.simple_reflections())
    cdef list index_list = range(len(I))

    cdef list level_set_cur = [(W.one(), [])]
    cdef set level_set_old = set([ W.one() ])
    cdef list word
    cdef PermutationGroupElement x, y

    while level_set_cur:
        level_set_new = []
        for x, word in level_set_cur:
            yield x, word
            for i in index_list:
                y = x._mul_(I[i])
                if y not in level_set_old:
                    level_set_old.add(y)
                    level_set_new.append((y, word+[i]))
        level_set_cur = level_set_new

cdef bint has_descent(PermutationGroupElement w, int i, int N):
    return w.perm[i] >= N

cdef int first_descent(PermutationGroupElement w, int n, int N):
    cdef int i
    for i in range(n):
        if has_descent(w,i,N):
            return i
    return -1

cpdef list reduced_word_c(W,w):
    r"""
    Computes a reduced word for the element `w` in the
    reflection group `W` in the positions ``range(n)``.

    EXAMPLES::

        sage: from sage.combinat.root_system.reflection_group_c import reduced_word_c
        sage: W = ReflectionGroup(['B',2])                  # optional - gap3
        sage: [ reduced_word_c(W,w) for w in W ]            # optional - gap3
        [[], [1], [0], [0, 1], [1, 0], [1, 0, 1], [0, 1, 0], [0, 1, 0, 1]]
    """
    cdef tuple S = tuple(W.simple_reflections())
    cdef int n = len(S)
    cdef int N = W._number_of_reflections
    cdef int fdes = 0
    cdef list word = []

    while fdes != -1:
        fdes = first_descent(w,n,N)
        si = S[fdes]
        w = si._mul_(w)
        word.append(fdes)
    return word[:-1]
