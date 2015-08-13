from matroid cimport Matroid
cdef class UnionMatroid(Matroid):
    r"""
    Union Matroids.

    When `A` is an `r` times `E` matrix, the linear matroid `M[A]` has ground
    set `E` and set of independent sets

        `I(A) =\{F \subseteq E :` the columns of `A` indexed by `F` are linearly independent `\}`

    The simplest way to create a LinearMatroid is by giving only a matrix `A`.
    Then, the groundset defaults to ``range(A.ncols())``. Any iterable object
    ``E`` can be given as a groundset. If ``E`` is a list, then ``E[i]`` will
    label the `i`-th column of `A`. Another possibility is to specify a
    *reduced* matrix `B`, to create the matroid induced by `A = [ I B ]`.

    INPUT:

    - ``matroids`` -- a iterator of matroids.

    OUTPUT:

    A ``UnionMatroid`` instance, it's a matroid union of all matroids in ``matroids``.
    """
    def __init__(self, matroids):
        self.matroids = list(matroids)
        E = set()
        for M in self.matroids:
            E.update(M.groundset())
        self._groundset = frozenset(E)

    cpdef groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT:

        A set.
        [NEEDUPDATE]
        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
        """
        return self._groundset

    cpdef _rank(self, X):
        sum_matroid = SumMatroid(self.matroids.delete(self.matroids.groundset()-X))
        d = {}
        for (i,x) in sum_matroid.groundset():
            if not x in d:
                d[x]=set()
            d[x].add(i)
        part_matroid = PartitionMatroid([[(i,x) for i in d[x]] for x in d])
        return len(self.sum_matroid.intersection(part_matroid))

cdef class SumMatroid(Matroid):
    def __init__(self, summands):
        self.summands = list(summands)
        E = set()
        for i in range(len(self.summands)):
            g = self.summands[i].groundset()
            E.update(zip([i]*len(g),g))
        self._groundset = frozenset(E)

    cpdef groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT:

        A set.

        [NEEDUPDATE]

        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
        """
        return self._groundset

    cpdef _rank(self, X):
        partition = {}
        for (i,x) in X:
            if not i in partition:
                partition[i] = set()
            partition[i].add(x)
        rk = 0
        for i, Xi in partition.iteritems():
            rk+= self.summands[i]._rank(Xi)
        return rk

cdef class PartitionMatroid(Matroid):
    def __init__(self, partition):
        self.p = list(partition)
        E = set()
        for P in partition:
            E.update(P)
        self._groundset = frozenset(E)

    cpdef groundset(self):
        """
        Return the groundset of the matroid.

        The groundset is the set of elements that comprise the matroid.

        OUTPUT:

        A set.
        [NEEDUPDATE]
        EXAMPLES::

            sage: M = matroids.named_matroids.Pappus()
            sage: sorted(M.groundset())
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
        """
        return self._groundset

    cpdef _rank(self, X):
        r"""
        Return the rank of a set ``X``.

        This method does no checking on ``X``, and ``X`` may be assumed to
        have the same interface as ``frozenset``.

        INPUT:

        - ``X`` -- an object with Python's ``frozenset`` interface.

        OUTPUT:

        Integer.

        [NEEDUPDATE]

        EXAMPLES::

            sage: M = sage.matroids.matroid.Matroid()
            sage: M._rank([0, 1, 2])
            Traceback (most recent call last):
            ...
            NotImplementedError: subclasses need to implement this.

        """
        X2 = set(X)
        used_indices = set()
        rk = 0
        while len(X2) > 0:
            e = X2.pop()
            for i in range(len(self.p)):
                if e in self.p[i]:
                    if i not in used_indices:
                        used_indices.add(i)
                        rk = rk + 1
                    break
        return rk