r"""
Integer vectors modulo the action of a permutation group
"""
#*****************************************************************************
#    Copyright (C) 2010-12 Nicolas Borie <nicolas.borie at math dot u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#              The full text of the GPL is available at:
#                    http://www.gnu.org/licenses/
#*****************************************************************************

from itertools import imap
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.semirings.all import NN

from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

from sage.structure.list_clone import ClonableIntArray
from sage.combinat.backtrack import SearchForest

from sage.combinat.enumeration_mod_permgroup import is_canonical, orbit, canonical_children, canonical_representative_of_orbit_of

from sage.combinat.integer_vector import IntegerVectors

class IntegerVectorsModPermutationGroup(UniqueRepresentation):
    r"""
    Returns an enumerated set containing integer vectors which are
    maximal in their orbit under the action of the permutation group
    ``G`` for the lexicographic order.

    In Sage, a permutation group `G` is viewed as a subgroup of the
    symmetric group `S_n` of degree `n` and `n` is said to be the degree
    of `G`.  Any integer vector `v` is said to be canonical if it
    is maximal in its orbit under the action of `G`. `v` is
    canonical if and only if

    .. math::

        v = \max_{\text{lex order}} \{g \cdot v | g \in G \}

    The action of `G` is on position. This means for example that the
    simple transposition `s_1 = (1, 2)` swaps the first and the second entries
    of any integer vector `v = [a_1, a_2, a_3, \dots , a_n]`

    .. math::

        s_1 \cdot v = [a_2, a_1, a_3, \dots , a_n]

    This functions returns a parent which contains a single integer
    vector by orbit under the action of the permutation group `G`. The
    approach chosen here is to keep the maximal integer vector for the
    lexicographic order in each orbit. Such maximal vector will be
    called canonical integer vector under the action of the
    permutation group `G`.

    INPUT:

    - ``G`` - a permutation group
    - ``sum`` - (default: None) - a nonnegative integer
    - ``max_part`` - (default: None) - a nonnegative integer setting the
      maximum of entries of elements
    - ``sgs`` - (default: None) - a strong generating system of the
      group `G`. If you do not provide it, it will be calculated at the
      creation of the parent

    OUTPUT:

    - If ``sum`` and ``max_part`` are None, it returns the infinite enumerated
      set of all integer vectors (list of integers) maximal in their orbit for
      the lexicographic order.

    - If ``sum`` is an integer, it returns a finite enumerated set containing
      all integer vectors maximal in their orbit for the lexicographic order
      and whose entries sum to ``sum``.

    EXAMPLES:

    Here is the set enumerating integer vectors modulo the action of the cyclic
    group of `3` elements::

        sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3)]]))
        sage: I.category()
        Join of Category of infinite enumerated sets and Category of quotients of sets
        sage: I.cardinality()
        +Infinity
        sage: I.list()
        Traceback (most recent call last):
        ...
        NotImplementedError: infinite list
        sage: p = iter(I)
        sage: for i in range(10): p.next()
        [0, 0, 0]
        [1, 0, 0]
        [2, 0, 0]
        [1, 1, 0]
        [3, 0, 0]
        [2, 1, 0]
        [2, 0, 1]
        [1, 1, 1]
        [4, 0, 0]
        [3, 1, 0]

    The method
    :meth:`~sage.combinat.integer_vectors_mod_permgroup.IntegerVectorsModPermutationGroup_All.is_canonical`
    tests if any integer vector is maximal in its orbit. This method
    is also used in the containment test::

        sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
        sage: I.is_canonical([5,2,0,4])
        True
        sage: I.is_canonical([5,0,6,4])
        False
        sage: I.is_canonical([1,1,1,1])
        True
        sage: [2,3,1,0] in I
        False
        sage: [5,0,5,0] in I
        True
        sage: 'Bla' in I
        False
        sage: I.is_canonical('bla')
        Traceback (most recent call last):
        ...
        AssertionError: bla should be a list or a integer vector

    If you give a value to the extra argument ``sum``, the set returned
    will be a finite set containing only canonical vectors whose entries
    sum to ``sum``.::

        sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3)]]), sum=6)
        sage: I.cardinality()
        10
        sage: I.list()
        [[6, 0, 0], [5, 1, 0], [5, 0, 1], [4, 2, 0], [4, 1, 1],
         [4, 0, 2], [3, 3, 0], [3, 2, 1], [3, 1, 2], [2, 2, 2]]
        sage: I.category()
        Join of Category of finite enumerated sets and Category of subquotients of finite sets and Category of quotients of sets

    To get the orbit of any integer vector `v` under the action of the group,
    use the method :meth:`~sage.combinat.integer_vectors_mod_permgroup.IntegerVectorsModPermutationGroup_All.orbit`;
    we convert the returned set of vectors into a list in increasing lexicographic order,
    to get a reproducible test::

        sage: from sage.combinat.enumeration_mod_permgroup import lex_cmp
        sage: sorted(I.orbit([6,0,0]), cmp=lex_cmp)
        [[0, 0, 6], [0, 6, 0], [6, 0, 0]]
        sage: sorted(I.orbit([5,1,0]), cmp=lex_cmp)
        [[0, 5, 1], [1, 0, 5], [5, 1, 0]]
        sage: I.orbit([2,2,2])
        set([[2, 2, 2]])

    TESTS:

    Let us check that canonical integer vectors of the symmetric group
    are just sorted list of integers::

        sage: I = IntegerVectorsModPermutationGroup(SymmetricGroup(5)) # long time
        sage: p = iter(I) # long time
        sage: for i in range(100): # long time
        ...       v = list(p.next())
        ...       assert sorted(v, reverse=True) == v

    We now check that there is as much of canonical vectors under the
    symmetric group `S_n` whose entries sum to `d` than partitions of
    `d` of at most `n` parts::

        sage: I = IntegerVectorsModPermutationGroup(SymmetricGroup(5)) # long time
        sage: for i in range(10): # long time
        ...       d1 = I.subset(i).cardinality()
        ...       d2 = Partitions(i, max_length=5).cardinality()
        ...       print d1
        ...       assert d1 == d2
        1
        1
        2
        3
        5
        7
        10
        13
        18
        23

    We present a last corner case: trivial groups. For the trivial
    group ``G`` acting on a list of length `n`, all integer vectors of
    length `n` are canonical::

        sage: G = PermutationGroup([[(6,)]]) # long time
        sage: G.cardinality() # long time
        1
        sage: I = IntegerVectorsModPermutationGroup(G) # long time
        sage: for i in range(10): # long time
        ...       d1 = I.subset(i).cardinality()
        ...       d2 = IntegerVectors(i,6).cardinality()
        ...       print d1
        ...       assert d1 == d2
        1
        6
        21
        56
        126
        252
        462
        792
        1287
        2002
    """
    @staticmethod
    def __classcall__(cls, G, sum=None, max_part=None, sgs=None):
        r"""
        TESTS::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3)]]))
            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3)]]), None)
            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3)]]), 2)
            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3)]]), -2)
            Traceback (most recent call last):
            ...
            ValueError: Value -2 in not in Non negative integer semiring.
            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3)]]), 8, max_part=5)
        """
        if sum is None and max_part is None:
            return IntegerVectorsModPermutationGroup_All(G, sgs=sgs)
        else:
            if sum is not None:
                assert (sum == NN(sum))
            if max_part is not None:
                assert (max_part == NN(max_part))
            return IntegerVectorsModPermutationGroup_with_constraints(G, sum, max_part, sgs=sgs)

class IntegerVectorsModPermutationGroup_All(UniqueRepresentation, SearchForest):
    r"""
    A class for integer vectors enumerated up to the action of a
    permutation group.

    A Sage permutation group is viewed as a subgroup of the symmetric
    group `S_n` for a certain `n`. This group has a natural action by
    position on vectors of length `n`. This class implements a set
    which keeps a single vector for each orbit. We say that a vector
    is canonical if it is the maximum in its orbit under the action of
    the permutation group for the lexicographic order.

    EXAMPLES::

        sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
        sage: I
        Integer vectors of length 4 enumerated up to the action of Permutation Group with generators [(1,2,3,4)]
        sage: I.cardinality()
        +Infinity
        sage: TestSuite(I).run()
        sage: it = iter(I)
        sage: [it.next(), it.next(), it.next(), it.next(), it.next()]
        [[0, 0, 0, 0],
         [1, 0, 0, 0],
         [2, 0, 0, 0],
         [1, 1, 0, 0],
         [1, 0, 1, 0]]
        sage: x = it.next(); x
        [3, 0, 0, 0]
        sage: I.first()
        [0, 0, 0, 0]

    TESTS::

        sage: TestSuite(I).run()
    """
    def __init__(self, G, sgs=None):
        """
        TESTS::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: I
            Integer vectors of length 4 enumerated up to the action of Permutation Group with generators [(1,2,3,4)]
            sage: I.category()
            Join of Category of infinite enumerated sets and Category of quotients of sets
            sage: TestSuite(I).run()
        """
        SearchForest.__init__(self, algorithm = 'breadth', category = InfiniteEnumeratedSets().Quotients())
        self._permgroup = G
        self.n = G.degree()

        # self.sgs: strong_generating_system
        if sgs is None:
            self._sgs = G.strong_generating_system()
        else:
            self._sgs = map(lambda x: list(x), list(sgs))

    def _repr_(self):
        """
        TESTS::

            sage: IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3)]]))
            Integer vectors of length 3 enumerated up to the action of Permutation Group with generators [(1,2,3)]
        """
        return "Integer vectors of length %s enumerated up to the action of %s"%(str(self.n), self._permgroup.__repr__())

    def ambient(self):
        r"""
        Return the ambient space from which ``self`` is a quotient.

        EXAMPLES::

            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: S.ambient()
            Integer vectors
        """
        # TODO: Fix me once 'IntegerVectors(length=bla)' will return
        # the integer vectors of length bla
        return IntegerVectors(length=self.n)

    def lift(self, elt):
        r"""
        Lift the element ``elt`` inside the ambient space from which ``self`` is a quotient.

        EXAMPLES::

            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: v = S.lift(S([4,3,0,1])); v
            [4, 3, 0, 1]
            sage: type(v)
            <type 'list'>
        """
        # TODO: For now, Sage integer vectors are just python list.
        # Once Integer vectors will have an element class, update this
        # code properly
        return list(elt)

    def retract(self, elt):
        r"""
        Return the canonical representative of the orbit of the
        integer ``elt`` under the action of the permutation group
        defining ``self``.

        If the element ``elt`` is already maximal in its orbit for
        the lexicographic order, ``elt`` is thus the good
        representative for its orbit.

        EXAMPLES::

            sage: [0,0,0,0] in IntegerVectors(length=4)
            True
            sage: [1,0,0,0] in IntegerVectors(length=4)
            True
            sage: [0,1,0,0] in IntegerVectors(length=4)
            True
            sage: [1,0,1,0] in IntegerVectors(length=4)
            True
            sage: [0,1,0,1] in IntegerVectors(length=4)
            True
            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: S.retract([0,0,0,0])
            [0, 0, 0, 0]
            sage: S.retract([1,0,0,0])
            [1, 0, 0, 0]
            sage: S.retract([0,1,0,0])
            [1, 0, 0, 0]
            sage: S.retract([1,0,1,0])
            [1, 0, 1, 0]
            sage: S.retract([0,1,0,1])
            [1, 0, 1, 0]
        """
        # TODO: Once Sage integer vector will have a data structure
        # based on ClonableIntArray, remove the conversion intarray
        assert len(elt) == self.n, "%s is a quotient set of %s"%(self, self.ambient())
        intarray = self.element_class(self, elt, check=False)
        return self.element_class(self, canonical_representative_of_orbit_of(self._sgs, intarray), check=False)

    def roots(self):
        r"""
        Returns the root of generation of ``self``. This method is
        required to build the tree structure of ``self`` which
        inherits from the class :class:`~sage.combinat.backtrack.SearchForest`.

        EXAMPLES::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: I.roots()
            [[0, 0, 0, 0]]
        """
        return [self.element_class(self, self.n*[0,], check=False)]

    def children(self, x):
        r"""
        Returns the list of children of the element ``x``. This method
        is required to build the tree structure of ``self`` which
        inherits from the class :class:`~sage.combinat.backtrack.SearchForest`.

        EXAMPLES::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: I.children(I([2,1,0,0], check=False))
            [[2, 2, 0, 0], [2, 1, 1, 0], [2, 1, 0, 1]]
        """
        return canonical_children(self._sgs, x, -1)

    def permutation_group(self):
        r"""
        Returns the permutation group given to define ``self``.

        EXAMPLES::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: I.permutation_group()
            Permutation Group with generators [(1,2,3,4)]
        """
        return self._permgroup

    def is_canonical(self, v, check=True):
        r"""
        Returns ``True`` if the integer list ``v`` is maximal in its
        orbit under the action of the permutation group given to
        define ``self``.  Such integer vectors are said to be
        canonical. A vector `v` is canonical if and only if

        .. math::

            v = \max_{\text{lex order}} \{g \cdot v | g \in G \}

        EXAMPLES::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: I.is_canonical([4,3,2,1])
            True
            sage: I.is_canonical([4,0,0,1])
            True
            sage: I.is_canonical([4,0,3,3])
            True
            sage: I.is_canonical([4,0,4,4])
            False
        """
        if check:
            assert isinstance(v, (ClonableIntArray, list)), '%s should be a list or a integer vector'%v
            assert (self.n == len(v)), '%s should be of length %s'%(v, self.n)
            for p in v:
                assert (p == NN(p)), 'Elements of %s should be integers'%s
        return is_canonical(self._sgs, self.element_class(self, list(v), check=False))

    def __contains__(self, v):
        """
        EXAMPLES::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: [2,2,0,0] in I
            True
            sage: [2,0,1,0] in I
            True
            sage: [2,0,0,1] in I
            True
            sage: [2,0,0,2] in I
            False
            sage: [2,0,0,2,12] in I
            False
        """
        try:
            return self.is_canonical(self.element_class(self, list(v), check=False), check=False)
        except Exception:
            return False

    def __call__(self, v, check=True):
        r"""
        Returns an element of ``self`` constructed from ``v`` if
        possible.

        TESTS::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: I([3,2,1,0])
            [3, 2, 1, 0]
        """
        try:
            if v.parent() is self:
                return v
            else:
                raise ValueError('%s shoud be a Python list of integer'%(v))
        except Exception:
            return self.element_class(self, list(v), check=check)

    def orbit(self, v):
        r"""
        Returns the orbit of the integer vector ``v`` under the action of the
        permutation group defining ``self``. The result is a set.

        EXAMPLES:

        In order to get reproducible doctests, we convert the returned sets
        into lists in increasing lexicographic order::

            sage: from sage.combinat.enumeration_mod_permgroup import lex_cmp
            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: sorted(I.orbit([2,2,0,0]), cmp=lex_cmp)
            [[0, 0, 2, 2], [0, 2, 2, 0], [2, 0, 0, 2], [2, 2, 0, 0]]
            sage: sorted(I.orbit([2,1,0,0]), cmp=lex_cmp)
            [[0, 0, 2, 1], [0, 2, 1, 0], [1, 0, 0, 2], [2, 1, 0, 0]]
            sage: sorted(I.orbit([2,0,1,0]), cmp=lex_cmp)
            [[0, 1, 0, 2], [0, 2, 0, 1], [1, 0, 2, 0], [2, 0, 1, 0]]
            sage: sorted(I.orbit([2,0,2,0]), cmp=lex_cmp)
            [[0, 2, 0, 2], [2, 0, 2, 0]]
            sage: I.orbit([1,1,1,1])
            set([[1, 1, 1, 1]])
        """
        assert isinstance(v, (list, ClonableIntArray)), '%s shoud be a Python list or an element of %s'%(v, self)
        try:
            if v.parent() is self:
                return orbit(self._sgs, v)
            raise TypeError
        except Exception:
            return orbit(self._sgs, self.element_class(self, v, check=False))

    def subset(self, sum=None, max_part=None):
        r"""
        Returns the subset of ``self`` containing integer vectors
        whose entries sum to ``sum``.

        EXAMPLES::

            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: S.subset(4)
            Integer vectors of length 4 and of sum 4 enumerated up to
            the action of Permutation Group with generators
            [(1,2,3,4)]
        """
        return IntegerVectorsModPermutationGroup_with_constraints(self.permutation_group(), sum, max_part)

    class Element(ClonableIntArray):
        r"""
        Element class for the set of integer vectors of given sum enumerated modulo
        the action of a permutation group. These vector are clonable lists of integers
        which must check conditions comming form the parent appearing in the method
        :meth:`~sage.structure.list_clone.ClonableIntArray.check`.

        TESTS::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: v = I.element_class(I, [4,3,2,1]); v
            [4, 3, 2, 1]
            sage: TestSuite(v).run()
            sage: I.element_class(I, [4,3,2,5])
            Traceback (most recent call last):
            ...
            AssertionError
        """
        def check(self):
            r"""
            Checks that ``self`` verify the invariants needed for
            living in ``self.parent()``.

            EXAMPLES::

                sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
                sage: v = I.an_element()
                sage: v.check()
                sage: w = I([0,4,0,0], check=False); w
                [0, 4, 0, 0]
                sage: w.check()
                Traceback (most recent call last):
                ...
                AssertionError
            """
            assert self.parent().is_canonical(self)


class IntegerVectorsModPermutationGroup_with_constraints(UniqueRepresentation, SearchForest):
    r"""
    This class models finite enumerated sets of integer vectors with
    constraint enumerated up to the action of a permutation group.
    Integer vectors are enumerated modulo the action of the
    permutation group. To implement that, we keep a single integer
    vector by orbit under the action of the permutation
    group. Elements chosen are vectors maximal in their orbit for the
    lexicographic order.

    For more information see :class:`IntegerVectorsModPermutationGroup`.

    EXAMPLES::

        sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), max_part=1)
        sage: I.list()
        [[0, 0, 0, 0], [1, 0, 0, 0], [1, 1, 0, 0], [1, 0, 1, 0], [1, 1, 1, 0], [1, 1, 1, 1]]
        sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), sum=6, max_part=4)
        sage: I.list()
        [[4, 2, 0, 0], [4, 1, 1, 0], [4, 1, 0, 1], [4, 0, 2, 0], [4, 0, 1, 1],
         [4, 0, 0, 2], [3, 3, 0, 0], [3, 2, 1, 0], [3, 2, 0, 1], [3, 1, 2, 0],
         [3, 1, 1, 1], [3, 1, 0, 2], [3, 0, 3, 0], [3, 0, 2, 1], [3, 0, 1, 2],
         [2, 2, 2, 0], [2, 2, 1, 1], [2, 1, 2, 1]]

    Here is the enumeration of unlabeled graphs over 5 vertices::

        sage: G = IntegerVectorsModPermutationGroup(TransitiveGroup(10,12), max_part=1) # optional
        sage: G.cardinality() # optional
        34

    TESTS::

        sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]),4)
        sage: TestSuite(I).run()
    """
    def __init__(self, G, d, max_part, sgs=None):
        r"""
        TESTS::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), 6, max_part=4)
        """
        SearchForest.__init__(self, algorithm = 'breadth', category = (FiniteEnumeratedSets(), FiniteEnumeratedSets().Quotients()))
        self._permgroup = G
        self.n = G.degree()
        self._sum = d
        if max_part is None:
            self._max_part = -1
        else:
            self._max_part = max_part

        # self.sgs: strong_generating_system
        if sgs is None:
            self._sgs = G.strong_generating_system()
        else:
            self._sgs = map(lambda x: list(x), list(sgs))

    def _repr_(self):
        r"""
        TESTS::

            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]])); S
            Integer vectors of length 4 enumerated up to the action of Permutation Group with generators [(1,2,3,4)]
            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), 6); S
            Integer vectors of length 4 and of sum 6 enumerated up to the action of Permutation Group with generators [(1,2,3,4)]
            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), 6, max_part=4); S
            Vectors of length 4 and of sum 6 whose entries is in {0, ..., 4} enumerated up to the action of Permutation Group with generators [(1,2,3,4)]
            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), max_part=4); S
            Integer vectors of length 4 whose entries is in {0, ..., 4} enumerated up to the action of Permutation Group with generators [(1,2,3,4)]
        """
        if self._sum is not None:
            if self._max_part >= 0:
                return "Vectors of length %s and of sum %s whose entries is in {0, ..., %s} enumerated up to the action of %s"%(self.n, self._sum, self._max_part, self.permutation_group())
            else:
                return "Integer vectors of length %s and of sum %s enumerated up to the action of %s"%(self.n, self._sum, self.permutation_group())
        else:
            return "Integer vectors of length %s whose entries is in {0, ..., %s} enumerated up to the action of %s"%(self.n, self._max_part, self.permutation_group())

    def roots(self):
        r"""
        Returns the root of generation of ``self``.This method is
        required to build the tree structure of ``self`` which
        inherits from the class
        :class:`~sage.combinat.backtrack.SearchForest`.

        EXAMPLES::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: I.roots()
            [[0, 0, 0, 0]]
        """
        return [self.element_class(self, self.n*[0,], check=False)]

    def children(self, x):
        r"""
        Returns the list of children of the element ``x``. This method
        is required to build the tree structure of ``self`` which
        inherits from the class
        :class:`~sage.combinat.backtrack.SearchForest`.

        EXAMPLES::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]))
            sage: I.children(I([2,1,0,0], check=False))
            [[2, 2, 0, 0], [2, 1, 1, 0], [2, 1, 0, 1]]
        """
        return canonical_children(self._sgs, x, -1)

    def permutation_group(self):
        r"""
        Returns the permutation group given to define ``self``.

        EXAMPLES::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3)]]), 5)
            sage: I.permutation_group()
            Permutation Group with generators [(1,2,3)]
        """
        return self._permgroup

    def __contains__(self, v):
        r"""
        TESTS::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]),6)
            sage: [6,0,0,0] in I
            True
            sage: [5,0,1,0] in I
            True
            sage: [0,5,1,0] in I
            False
            sage: [3,0,1,3] in I
            False
            sage: [3,3,1,0] in I
            False
        """
        try:
            return (self(v)).parent() is self
        except Exception:
            return False

    def __call__(self, v, check=True):
        r"""
        Make `v` an element living in ``self``.

        TESTS::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), 4)
            sage: v = I([2,1,0,1]); v
            [2, 1, 0, 1]
            sage: v.parent()
            Integer vectors of length 4 and of sum 4 enumerated up to
            the action of Permutation Group with generators
            [(1,2,3,4)]
        """
        try:
            if v.parent() is self:
                return v
            else:
                raise ValueError('%s shoud be a Python list of integer'%(v))
        except Exception:
            return self.element_class(self, list(v), check=check)

    def __iter__(self):
        r"""
        TESTS::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]),4)
            sage: for i in I: i
            [4, 0, 0, 0]
            [3, 1, 0, 0]
            [3, 0, 1, 0]
            [3, 0, 0, 1]
            [2, 2, 0, 0]
            [2, 1, 1, 0]
            [2, 1, 0, 1]
            [2, 0, 2, 0]
            [2, 0, 1, 1]
            [1, 1, 1, 1]
            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), sum=7, max_part=3)
            sage: for i in I: i
            [3, 3, 1, 0]
            [3, 3, 0, 1]
            [3, 2, 2, 0]
            [3, 2, 1, 1]
            [3, 2, 0, 2]
            [3, 1, 3, 0]
            [3, 1, 2, 1]
            [3, 1, 1, 2]
            [3, 0, 2, 2]
            [2, 2, 2, 1]
        """
        if self._max_part < 0:
            return self.elements_of_depth_iterator(self._sum)
        else:
            SF = SearchForest((self([0]*(self.n), check=False),),
                              lambda x : map(lambda y : self(y, check=False), canonical_children(self._sgs, x, self._max_part)),
                              algorithm = 'breadth')
            if self._sum is None:
                return iter(SF)
            else:
                return SF.elements_of_depth_iterator(self._sum)

    def is_canonical(self, v, check=True):
        r"""
        Returns ``True`` if the integer list ``v`` is maximal in its
        orbit under the action of the permutation group given to
        define ``self``.  Such integer vectors are said to be
        canonical. A vector `v` is canonical if and only if

        .. math::

            v = \max_{\text{lex order}} \{g \cdot v | g \in G \}

        EXAMPLES::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), max_part=3)
            sage: I.is_canonical([3,0,0,0])
            True
            sage: I.is_canonical([1,0,2,0])
            False
            sage: I.is_canonical([2,0,1,0])
            True
        """
        if check:
            assert isinstance(v, (ClonableIntArray, list)), '%s should be a list or a integer vector'%v
            assert (self.n == len(v)), '%s should be of length %s'%(v, self.n)
            for p in v:
                assert (p == NN(p)), 'Elements of %s should be integers'%s
        return is_canonical(self._sgs, self.element_class(self, list(v), check=False))

    def ambient(self):
        r"""
        Return the ambient space from which ``self`` is a quotient.

        EXAMPLES::

            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), 6); S.ambient()
            Integer vectors that sum to 6
            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), 6, max_part=12); S.ambient()
            Integer vectors that sum to 6 with constraints: max_part=12

        .. todo::

         Integer vectors should accept ``max_part`` as a single argument, and the following should change::

            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), max_part=12); S.ambient()
            Integer vectors

        """
        if self._sum is not None:
            if self._max_part <= -1:
                return IntegerVectors(n=self._sum)
            else:
                return IntegerVectors(n=self._sum, max_part=self._max_part)
        else:
            # Fix me once max_part should be accepted as a single
            # argument for integer vectors
            return IntegerVectors(max_part=self._max_part)

    def lift(self, elt):
        r"""
        Lift the element ``elt`` inside the ambient space from which ``self`` is a quotient.

        EXAMPLES::

            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), max_part=1)
            sage: v = S.lift([1,0,1,0]); v
            [1, 0, 1, 0]
            sage: v in IntegerVectors(max_part=1)
            True
            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), sum=6)
            sage: v = S.lift(S.list()[5]); v
            [4, 1, 1, 0]
            sage: v in IntegerVectors(n=6)
            True
        """
        # TODO: For now, Sage integer vectors are just python list.
        # Once Integer vectors will have an element class, update this
        # code properly
        return list(elt)

    def retract(self, elt):
        r"""
        Return the canonical representative of the orbit of the
        integer ``elt`` under the action of the permutation group
        defining ``self``.

        If the element ``elt`` is already maximal in its orbits for
        the lexicographic order, ``elt`` is thus the good
        representative for its orbit.

        EXAMPLES::

            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), sum=2, max_part=1)
            sage: S.retract([1,1,0,0])
            [1, 1, 0, 0]
            sage: S.retract([1,0,1,0])
            [1, 0, 1, 0]
            sage: S.retract([1,0,0,1])
            [1, 1, 0, 0]
            sage: S.retract([0,1,1,0])
            [1, 1, 0, 0]
            sage: S.retract([0,1,0,1])
            [1, 0, 1, 0]
            sage: S.retract([0,0,1,1])
            [1, 1, 0, 0]
        """
        # TODO: Once Sage integer vector will have a data structure
        # based on ClonableIntArray, remove the conversion intarray
        assert len(elt) == self.n, "%s is a quotient set of %s"%(self, self.ambient())
        if self._sum is not None:
            assert sum(elt) == self._sum, "%s is a quotient set of %s"%(self, self.ambient())
        if self._max_part >= 0:
            assert max(elt) <= self._max_part, "%s is a quotient set of %s"%(self, self.ambient())
        intarray = self.element_class(self, elt, check=False)
        return self.element_class(self, canonical_representative_of_orbit_of(self._sgs, intarray), check=False)

    def an_element(self):
        r"""
        Returns an element of ``self`` or raises an EmptySetError when
        ``self`` is empty.

        EXAMPLES::

            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), sum=0, max_part=1); S.an_element()
            [0, 0, 0, 0]
            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), sum=1, max_part=1); S.an_element()
            [1, 0, 0, 0]
            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), sum=2, max_part=1); S.an_element()
            [1, 1, 0, 0]
            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), sum=3, max_part=1); S.an_element()
            [1, 1, 1, 0]
            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), sum=4, max_part=1); S.an_element()
            [1, 1, 1, 1]
            sage: S = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), sum=5, max_part=1); S.an_element()
            Traceback (most recent call last):
            ...
            EmptySetError
        """
        if self._max_part < 0:
            return self([self._sum]+(self.n-1)*[0], check=False)
        else:
            try:
                v = iter(self)
                return v.next()
            except StopIteration:
                from sage.categories.sets_cat import EmptySetError
                raise EmptySetError

    def orbit(self, v):
        r"""
        Returns the orbit of the vector ``v`` under the action of the
        permutation group defining ``self``. The result is a set.

        INPUT:

        - ``v`` - an element of ``self`` or any list of length the
          degree of the permutation group.

        EXAMPLES:

        We convert the result in a list in increasing lexicographic
        order, to get a reproducible doctest::

            sage: from sage.combinat.enumeration_mod_permgroup import lex_cmp
            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]),4)
            sage: I.orbit([1,1,1,1])
            set([[1, 1, 1, 1]])
            sage: sorted(I.orbit([3,0,0,1]), cmp=lex_cmp)
            [[0, 0, 1, 3], [0, 1, 3, 0], [1, 3, 0, 0], [3, 0, 0, 1]]
        """
        assert isinstance(v, (list, ClonableIntArray)), '%s shoud be a Python list or an element of %s'%(v, self)
        try:
            if v.parent() is self:
                return orbit(self._sgs, v)
        except Exception:
            return orbit(self._sgs, self.element_class(self, v, check=False))

    class Element(ClonableIntArray):
        r"""
        Element class for the set of integer vectors with constraints enumerated
        modulo the action of a permutation group. These vectors are clonable lists
        of integers which must check conditions comming form the parent as in
        the method :meth:`~sage.combinat.integer_vectors_mod_permgroup.IntegerVectorsModPermutationGroup_with_constraints.Element.check`.

        TESTS::

            sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), 4)
            sage: v = I.element_class(I, [3,1,0,0]); v
            [3, 1, 0, 0]
            sage: TestSuite(v).run()
            sage: v = I.element_class(I, [3,2,0,0])
            Traceback (most recent call last):
            ...
            AssertionError: [3, 2, 0, 0] should be a integer vector of sum 4
        """
        def check(self):
            r"""
            Checks that ``self`` meets the constraints of being an element of ``self.parent()``.

            EXAMPLES::

                sage: I = IntegerVectorsModPermutationGroup(PermutationGroup([[(1,2,3,4)]]), 4)
                sage: v = I.an_element()
                sage: v.check()
                sage: w = I([0,4,0,0], check=False); w
                [0, 4, 0, 0]
                sage: w.check()
                Traceback (most recent call last):
                ...
                AssertionError
            """
            if self.parent()._sum is not None:
                assert sum(self) == self.parent()._sum, '%s should be a integer vector of sum %s'%(self, self.parent()._sum)
            if self.parent()._max_part >= 0:
                assert max(self) <= self.parent()._max_part, 'Entries of %s must be inferiors to %s'%(self, self.parent()._max_part)
            assert self.parent().is_canonical(self)
