from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.combinat.permutation import Permutation


class BaxterPermutations(UniqueRepresentation, Parent):
    r"""
    The combinatorial class of Baxter permutations. A Baxter permutation
    is a permutation avoiding the generalized permutation patterns
    `2-41-3` and `3-14-2`. In other words, a permutation `\sigma` is a
    Baxter permutation if for any subword `u := u_1u_2u_3u_4` of `\sigma`
    such that the letters `u_2` and `u_3` are adjacent in `\sigma`, the
    standardized version of `u` is nor `2413` neither `3142`.

    INPUT:

    - ``n`` -- (default: ``None``) a non negative integer, the size of
    the permutations.

    OUTPUT:

    Returns the combinatorial class of the Baxter permutations of size ``n``
    if ``n`` is not ``None``. Otherwise, returns the combinatorial class
    of all Baxter permutations.

    EXAMPLES::

        sage: BaxterPermutations(5)
        Baxter permutations of size 5
        sage: BaxterPermutations()
        Baxter permutations
    """
    @staticmethod
    def __classcall_private__(classe, n=None):
        if n is None:
            return BaxterPermutations_all()
        return BaxterPermutations_size(n)


class BaxterPermutations_size(BaxterPermutations):
    r"""
    The enumerated set of Baxter permutations of a given size.

    EXAMPLES::

        sage: from sage.combinat.baxter_permutations import BaxterPermutations_size
        sage: BaxterPermutations_size(5)
        Baxter permutations of size 5
    """

    def __init__(self, n):
        self.element_class = type(Permutation([]))
        self._n = n
        from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
        super(BaxterPermutations, self).__init__(category=FiniteEnumeratedSets())

    def _element_constructor_(self, x):
        return self.element_class(x)

    def _repr_(self):
        """
        Return a string representation of ``self``

        EXAMPLES::
        """
        return "Baxter permutations of size %s" % self._n

    def __contains__(self, x):
        r"""
        INPUT:

        - ``x`` -- a list or tuple or permutation.

        OUTPUT:

        Returns ``True`` if and only if ``x`` is a Baxter permutation of
        size ``self._n``.

        EXAMPLES::

            sage: Permutation([2, 1, 4, 3]) in BaxterPermutations(4)
            True
            sage: Permutation([2, 1, 4, 3]) in BaxterPermutations(5)
            False
            sage: Permutation([3, 1, 4, 2]) in BaxterPermutations(4)
            False
            sage: [len([p for p in Permutations(n) if Permutation(list(p)) in BaxterPermutations(n)]) for n in range(7)]
            [1, 1, 2, 6, 22, 92, 422]
            sage: sorted([p for p in Permutations(6) if Permutation(list(p)) in BaxterPermutations(6)]) == sorted(BaxterPermutations(6).list())
            True
        """
        if not isinstance(x, self.element_class) or len(x) != self._n:
            return False
        for i in range(1, len(x) - 1):
            a = x[i]
            b = x[i + 1]
            if a < b:  # Hunting pattern 3-14-2.
                max_l = 0
                for j in range(i):
                    if x[j] > a and x[j] < b and x[j] > max_l:
                        max_l = x[j]
                min_r = len(x) + 1
                for j in range(i + 2, len(x)):
                    if x[j] > a and x[j] < b and x[j] < min_r:
                        min_r = x[j]
                if max_l != 0 and min_r != 0 and max_l > min_r:
                    return False
            else:  # Hunting pattern 2-41-3.
                min_l = len(x) + 1
                for j in range(i):
                    if x[j] < a and x[j] > b and x[j] < min_l:
                        min_l = x[j]
                max_r = 0
                for j in range(i + 2, len(x)):
                    if x[j] < a and x[j] > b and x[j] > max_r:
                        max_r = x[j]
                if min_l != 0 and max_r != 0 and min_l < max_r:
                    return False
        return True

    def __iter__(self):
        r"""
        Efficient generation of Baxter permutations.

        OUTPUT:

        An iterator over the Baxter permutations of size ``self._n``.

        EXAMPLES::

            sage: BaxterPermutations(4).list()
            [[4, 3, 2, 1], [3, 4, 2, 1], [3, 2, 4, 1], [3, 2, 1, 4], [2, 4, 3, 1],
             [4, 2, 3, 1], [2, 3, 4, 1], [2, 3, 1, 4], [2, 1, 4, 3], [4, 2, 1, 3],
             [2, 1, 3, 4], [1, 4, 3, 2], [4, 1, 3, 2], [1, 3, 4, 2], [1, 3, 2, 4],
             [4, 3, 1, 2], [3, 4, 1, 2], [3, 1, 2, 4], [1, 2, 4, 3], [1, 4, 2, 3],
             [4, 1, 2, 3], [1, 2, 3, 4]]
            sage: [len(BaxterPermutations(n)) for n in xrange(9)]
            [1, 1, 2, 6, 22, 92, 422, 2074, 10754]

        ALGORITHM::

        The algorithm using generating trees described in [BBF08]_ is used.
        The idea is that all Baxter permutations of size `n + 1` can be
        obtained by inserting the letter `n + 1` either just before a left
        to right maximum or just after an right to left maximum of a Baxter
        permutation of size `n`.

        REFERENCES::

        .. [BBF08] N. Bonichon, M. Bousquet-Melou, E. Fusy.
           Baxter permutations and plane bipolar orientations.
           Seminaire Lotharingien de combinatoire 61A, article B61Ah, 2008.
        """
        if self._n == 0:
            yield Permutation([])
        elif self._n == 1:
            yield Permutation([1])
        else:
            for b in BaxterPermutations(self._n - 1):
                # Left to right maxima.
                for i in [self._n - 2 - i for i in b.reverse().saliances()]:
                    yield Permutation(b[:i] + [self._n] + b[i:])
                # Right to left maxima.
                for i in b.saliances():
                    yield Permutation(b[:i + 1] + [self._n] + b[i + 1:])

    def _an_element_(self):
        """
        """
        return self.first()

    def cardinality(self):
        r"""
        The number of Baxter permutations.

        OUTPUT:

        Returns the number of Baxter permutations of size ``self._n``.

        EXAMPLES::

            sage: [BaxterPermutations(n).cardinality() for n in xrange(13)]
            [1, 1, 2, 6, 22, 92, 422, 2074, 10754, 58202, 326240, 1882960, 11140560]
        """
        if self._n == 0:
            return 1
        from sage.rings.arith import binomial
        return sum((binomial(self._n + 1, k) *
                    binomial(self._n + 1, k + 1) *
                    binomial(self._n + 1, k + 2)) /
                   ((self._n + 1) * binomial(self._n + 1, 2))
                   for k in xrange(self._n))


class BaxterPermutations_all(DisjointUnionEnumeratedSets, BaxterPermutations):
    r"""
    The enumerated set of all Baxter permutations.

    EXAMPLES::

        sage: from sage.combinat.baxter_permutations import BaxterPermutations_all
        sage: BaxterPermutations_all()
        Baxter permutations
    """
    def __init__(self, n=None):
        self.element_class = type(Permutation([]))
        from sage.categories.examples.infinite_enumerated_sets import NonNegativeIntegers
        from sage.sets.family import Family
        DisjointUnionEnumeratedSets.__init__(self,
                                             Family(NonNegativeIntegers(),
                                                    BaxterPermutations_size),
                                             facade=False, keepkey=False)

    def _repr_(self):
        """
        """
        return "Baxter permutations"

    def __contains__(self, x):
        r"""
        INPUT:

        - ``x`` -- any object.

        OUTPUT:

        Returns ``True`` if and only if ``x`` is a Baxter permutation.

        EXAMPLES::

            sage: Permutation([4, 2, 1, 7, 3, 8, 5, 6]) in BaxterPermutations()
            False
            sage: Permutation([4, 3, 6, 9, 7, 5, 1, 2, 8]) in BaxterPermutations()
            True
        """
        if not isinstance(x, self.element_class):
            return False
        return x in BaxterPermutations(len(x))

    def to_pair_of_twin_binary_trees(self, p):
        r"""
        Computes a bijection between Baxter permutations of size ``self._n``
        and the set of pairs of twin binary trees with ``self._n`` nodes.

        INPUT:

        - ``p`` -- a Baxter permutation.

        OUTPUT:

        Returns the pair of twin binary trees `(T_L, T_R)` where `T_L`
        (resp. `T_R`) is obtained by inserting the letters of ``p`` from
        left to right (resp. right to left) following the the binary search
        tree insertion algorithm.

        EXAMPLES::

            sage: BaxterPermutations().to_pair_of_twin_binary_trees(Permutation([]))
            (., .)
            sage: BaxterPermutations().to_pair_of_twin_binary_trees(Permutation([1, 2, 3]))
            (1[., 2[., 3[., .]]], 3[2[1[., .], .], .])
            sage: BaxterPermutations().to_pair_of_twin_binary_trees(Permutation([3, 4, 1, 2]))
            (3[1[., 2[., .]], 4[., .]], 2[1[., .], 4[3[., .], .]])
        """
        from sage.combinat.binary_tree import LabelledBinaryTree
        left = LabelledBinaryTree(None)
        right = LabelledBinaryTree(None)
        for a in p:
            left = left.binary_search_insert(a)
        for a in p.reverse():
            right = right.binary_search_insert(a)
        return (left, right)
