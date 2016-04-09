from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement
from collections import deque

cdef class Iterator(object):
    cdef int n
    cdef list S
    cdef list is_positive_root
    cdef str algorithm

    def __init__(self, W, algorithm="depth"):
        self.S = W.gens()
        self.n = len(W._index_set)
        self.is_positive_root = W._is_positive_root

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
            return self.iter_depth()
        elif self.algorithm == "breadth":
            return self.iter_breadth()

    def iter_depth(self):
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

    def iter_breadth(self):
        cdef tuple node
        cdef list cur = [(self.S[0].parent().one(), -1)]
        cdef PermutationGroupElement u
        cdef int first
        # Using a deque with popleft is ~2x slower
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
