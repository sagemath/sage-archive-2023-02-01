from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement
from collections import deque

cdef class Iterator(object):
    cdef int n
    cdef list S
    cdef list is_positive_root
    cdef str algorithm
    cdef bint tracking_words

    def __init__(self, W, algorithm="depth", tracking_words=True):
        self.S = W.gens()
        self.n = len(W._index_set)
        self.is_positive_root = W._is_positive_root
        self.tracking_words = tracking_words

        # "breadth" is 1.5 - 2x slower than "depth" since it uses
        # a deque with popleft instead of a list with pop
        if algorithm not in ["depth","breadth"]:
            raise ValueError('The algorithm (="%s") must be either "depth" or "breadth"')
        self.algorithm = algorithm

    cpdef list succ(self, PermutationGroupElement u, int first):
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

    cpdef list succ_words(self, PermutationGroupElement u, list word, int first):
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

    cpdef bint test(self, PermutationGroupElement u1, int i):
        cdef int j

        for j in range(i):
            if not self.is_positive_root[u1.perm[j]+1]:
                return False
        return True

    def __iter__(self):
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
        cdef tuple node
        cdef list cur = [(self.S[0].parent().one(), -1)]
        cdef PermutationGroupElement u
        cdef int first
        cdef list L = []

        #print "doing depth"

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
        cdef tuple node
        cdef list cur, word
        
        cdef PermutationGroupElement u
        cdef int first
        cdef list L = []

        #print "doing word depth"

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
        cdef tuple node
        cdef list cur = [(self.S[0].parent().one(), -1)]
        cdef PermutationGroupElement u
        cdef int first
        L = deque()

        #print "doing breadth"

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
        cdef tuple node
        cdef list cur, word
        cdef PermutationGroupElement u
        cdef int first
        L = deque()

        one = self.S[0].parent().one()
        one._reduced_word = []
        cur = [(one, list(), -1)]

        #print "doing word breadth"

        while True:
            if not cur:
                if not L:
                    return
                cur = L.popleft()
                continue

            u, word, first = cur.pop()
            yield u
            L.append(self.succ_words(u, word, first))
