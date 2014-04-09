r"""
Subsets

The combinatorial class of the subsets of a finite set. The set can
be given as a list or a Set or else as an integer `n` which encodes the set
`\{1,2,...,n\}`. See :class:`Subsets` for more information and examples.

AUTHORS:

- Mike Hansen: initial version

- Florent Hivert (2009/02/06): doc improvements + new methods

"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.sets.set import Set
from sage.rings.arith import binomial
from sage.rings.integer import Integer
import sage.combinat.subword as subword
import sage.combinat.choose_nk as choose_nk
import sage.misc.prandom as rnd
import __builtin__
import itertools
from combinat import CombinatorialClass
from sage.sets.set import Set_object_enumerated
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets


def Subsets(s, k=None, submultiset=False):
    """
    Returns the combinatorial class of the subsets of the finite set
    s. The set can be given as a list, Set or any iterable convertible
    to a set. It can alternatively be given a non-negative integer `n`
    which encode the set `\{1,2,\dots,n\}` (i.e. the Sage
    ``range(1,s+1)``).

    A second optional parameter k can be given. In this case, Subsets returns
    the combinatorial class of subsets of s of size k.

    Finally the option ``submultiset`` allows one to deal with sets with
    repeated elements usually called multisets.

    EXAMPLES::

        sage: S = Subsets([1, 2, 3]); S
        Subsets of {1, 2, 3}
        sage: S.cardinality()
        8
        sage: S.first()
        {}
        sage: S.last()
        {1, 2, 3}
        sage: S.random_element()
        {2}
        sage: S.list()
        [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]

    Here is the same example where the set is given as an integer::

        sage: S = Subsets(3)
        sage: S.list()
        [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]

    We demonstrate various the effect of the various options::

        sage: S = Subsets(3, 2); S
        Subsets of {1, 2, 3} of size 2
        sage: S.list()
        [{1, 2}, {1, 3}, {2, 3}]

        sage: S = Subsets([1, 2, 2], submultiset=True); S
        SubMultiset of [1, 2, 2]
        sage: S.list()
        [[], [1], [2], [1, 2], [2, 2], [1, 2, 2]]

        sage: S = Subsets([1, 2, 2, 3], 3, submultiset=True); S
        SubMultiset of [1, 2, 2, 3] of size 3
        sage: S.list()
        [[1, 2, 2], [1, 2, 3], [2, 2, 3]]

        sage: S = Subsets(['a','b','a','b'], 2, submultiset=True); S.list()
        [['a', 'a'], ['a', 'b'], ['b', 'b']]
    """
    if k is not None:
        k=Integer(k)

    if isinstance(s, (int, Integer)):
        if s < 0:
            raise ValueError, "s must be non-negative"
        s = Set(range(1,s+1))

#    if len(Set(s)) != len(s):
#        multi = True

    if k is None:
        if submultiset:
            return SubMultiset_s(s)
        else:
            return Subsets_s(s)
    else:
        if submultiset:
            return SubMultiset_sk(s, k)
        else:
            return Subsets_sk(s, k)





class Subsets_s(CombinatorialClass):
    def __init__(self, s):
        """
        TESTS::

            sage: s = Subsets(Set([1]))
            sage: e = s.first()
            sage: isinstance(e, s.element_class)
            True

        In the following "_test_elements" is temporarily disabled
        until :class:`sage.sets.set.Set_object_enumerated` objects
        pass the category tests::

            sage: S = Subsets([1,2,3])
            sage: TestSuite(S).run(skip=["_test_elements"])

            sage: S = sage.sets.set.Set_object_enumerated([1,2])
            sage: TestSuite(S).run()         # todo: not implemented
        """
        CombinatorialClass.__init__(self, category=FiniteEnumeratedSets())
        self.s = Set(s)

    def __repr__(self):
        """
        TESTS::

            sage: repr(Subsets([1,2,3]))
            'Subsets of {1, 2, 3}'
        """
        return "Subsets of %s"%self.s

    def __contains__(self, value):
        """
        TESTS::

            sage: S = Subsets([1,2,3])
            sage: Set([1,2]) in S
            True
            sage: Set([1,4]) in S
            False
            sage: Set([]) in S
            True
        """
        value = Set(value)
        for v in value:
            if not v in self.s:
                return False
        return True

    def cardinality(self):
        r"""
        Returns the number of subsets of the set s.

        This is given by `2^{|s|}`.

        EXAMPLES::

            sage: Subsets(Set([1,2,3])).cardinality()
            8
            sage: Subsets([1,2,3,3]).cardinality()
            8
            sage: Subsets(3).cardinality()
            8
        """
        return Integer(2**len(self.s))

    def first(self):
        """
        Returns the first subset of s. Since we aren't restricted to
        subsets of a certain size, this is always the empty set.

        EXAMPLES::

            sage: Subsets([1,2,3]).first()
            {}
            sage: Subsets(3).first()
            {}
        """
        return Set([])

    def last(self):
        """
        Returns the last subset of s. Since we aren't restricted to subsets
        of a certain size, this is always the set s itself.

        EXAMPLES::

            sage: Subsets([1,2,3]).last()
            {1, 2, 3}
            sage: Subsets(3).last()
            {1, 2, 3}
        """
        return self.s


    def __iter__(self):
        """
        Iterates through the subsets of s.

        EXAMPLES::

            sage: [sub for sub in Subsets(Set([1,2,3]))]
            [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]
            sage: [sub for sub in Subsets(3)]
            [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]
            sage: [sub for sub in Subsets([1,2,3,3])]
            [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]

        """
        lset = __builtin__.list(self.s)
        #We use the iterator for the subwords of range(len(self.s))
        ind_set = lambda index_list: Set([lset[i] for i in index_list])
        it = itertools.imap(ind_set, subword.Subwords(range(len(lset))))
        for sub in it:
            yield sub

    def random_element(self):
        """
        Returns a random element of the class of subsets of s (in other
        words, a random subset of s).

        EXAMPLES::

            sage: Subsets(3).random_element()
            {2}
            sage: Subsets([4,5,6]).random_element()
            {5}
        """
        lset = __builtin__.list(self.s)
        n = len(self.s)
        return Set(filter(lambda x: rnd.randint(0,1), lset))

    def rank(self, sub):
        """
        Returns the rank of sub as a subset of s.

        EXAMPLES::

            sage: Subsets(3).rank([])
            0
            sage: Subsets(3).rank([1,2])
            4
            sage: Subsets(3).rank([1,2,3])
            7
            sage: Subsets(3).rank([2,3,4]) == None
            True
        """
        subset = Set(sub)
        lset = __builtin__.list(self.s)
        lsubset = __builtin__.list(subset)

        try:
            index_list = sorted(map(lambda x: lset.index(x), lsubset))
        except ValueError:
            return None

        n = len(self.s)
        r = 0

        for i in range(len(index_list)):
            r += binomial(n,i)
        return r + choose_nk.rank(index_list,n)

    def unrank(self, r):
        """
        Returns the subset of s that has rank k.

        EXAMPLES::

            sage: Subsets(3).unrank(0)
            {}
            sage: Subsets([2,4,5]).unrank(1)
            {2}
        """

        lset = __builtin__.list(self.s)
        n = len(lset)

        if r >= self.cardinality() or r < 0:
            return None
        else:
            for k in range(n+1):
                bin = binomial(n,k)
                if r >= bin:
                    r = r - bin
                else:
                    return Set([lset[i] for i in choose_nk.from_rank(r, n, k)])

    def _an_element_(self):
        """
        Returns an example of subset.

        EXAMPLES::

            sage: Subsets(0)._an_element_()
            {}
            sage: Subsets(3)._an_element_()
            {1, 2}
            sage: Subsets([2,4,5])._an_element_()
            {2, 4}
        """
        return self.unrank(self.cardinality() // 2)

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: S3 = Subsets(3); S3([1,2])
            {1, 2}
            sage: S3([0,1,2])
            Traceback (most recent call last):
            ...
            ValueError: [0, 1, 2] not in Subsets of {1, 2, 3}
        """
        return Set(x)

    element_class = Set_object_enumerated


class Subsets_sk(CombinatorialClass):
    def __init__(self, s, k):
        """
        TESTS::

            sage: s = Subsets(Set([1]))
            sage: e = s.first()
            sage: isinstance(e, s.element_class)
            True

        In the following "_test_elements" is temporarily disabled
        until :class:`sage.sets.set.Set_object_enumerated` objects
        pass the category tests::

            sage: S = Subsets(3,2)
            sage: TestSuite(S).run(skip=["_test_elements"])
        """
        CombinatorialClass.__init__(self, category=FiniteEnumeratedSets())
        self.s = Set(s)
        self.k = k

    def __repr__(self):
        """
        TESTS::

            sage: repr(Subsets(3,2))
            'Subsets of {1, 2, 3} of size 2'
        """
        return "Subsets of %s of size %s"%(self.s, self.k)

    def __contains__(self, value):
        """
        TESTS:
            sage: S = Subsets([1,2,3], 2)
            sage: Set([1,2]) in S
            True
            sage: Set([1,4]) in S
            False
            sage: Set([]) in S
            False
        """
        value = Set(value)
        if len(value) != self.k:
            return False
        for v in value:
            if not v in self.s:
                return False
        return True

    def cardinality(self):
        """
        EXAMPLES::

            sage: Subsets(Set([1,2,3]), 2).cardinality()
            3
            sage: Subsets([1,2,3,3], 2).cardinality()
            3
            sage: Subsets([1,2,3], 1).cardinality()
            3
            sage: Subsets([1,2,3], 3).cardinality()
            1
            sage: Subsets([1,2,3], 0).cardinality()
            1
            sage: Subsets([1,2,3], 4).cardinality()
            0
            sage: Subsets(3,2).cardinality()
            3
            sage: Subsets(3,4).cardinality()
            0
        """
        if self.k not in range(len(self.s)+1):
            return 0
        else:
            return binomial(len(self.s),self.k)


    def first(self):
        """
        Returns the first subset of s of size k.

        EXAMPLES::

            sage: Subsets(Set([1,2,3]), 2).first()
            {1, 2}
            sage: Subsets([1,2,3,3], 2).first()
            {1, 2}
            sage: Subsets(3,2).first()
            {1, 2}
            sage: Subsets(3,4).first()
        """
        if self.k not in range(len(self.s)+1):
            return None
        else:
            return Set(__builtin__.list(self.s)[:self.k])



    def last(self):
        """
        Returns the last subset of s of size k.

        EXAMPLES::

            sage: Subsets(Set([1,2,3]), 2).last()
            {2, 3}
            sage: Subsets([1,2,3,3], 2).last()
            {2, 3}
            sage: Subsets(3,2).last()
            {2, 3}
            sage: Subsets(3,4).last()
        """
        if self.k not in range(len(self.s)+1):
            return None
        else:
            return Set(__builtin__.list(self.s)[-self.k:])




    def __iter__(self):
        """
        Iterates through the subsets of s of size k.

        EXAMPLES::

            sage: [sub for sub in Subsets(Set([1,2,3]), 2)]
            [{1, 2}, {1, 3}, {2, 3}]
            sage: [sub for sub in Subsets([1,2,3,3], 2)]
            [{1, 2}, {1, 3}, {2, 3}]
            sage: [sub for sub in Subsets(3,2)]
            [{1, 2}, {1, 3}, {2, 3}]
        """
        if self.k not in range(len(self.s)+1):
            return

        lset = __builtin__.list(self.s)
        #We use the iterator for the subwords of range(len(self.s))
        ind_set = lambda index_list: Set([lset[i] for i in index_list])
        for sub in choose_nk.ChooseNK(len(lset),self.k):
            yield ind_set(sub)



    def random_element(self):
        """
        Returns a random element of the class of subsets of s of size k (in
        other words, a random subset of s of size k).

        EXAMPLES::

            sage: Subsets(3, 2).random_element()
            {1, 2}
            sage: Subsets(3,4).random_element() is None
            True
        """
        lset = __builtin__.list(self.s)
        n = len(self.s)

        if self.k not in range(len(self.s)+1):
            return None
        else:
            return Set([lset[i] for i in choose_nk.ChooseNK(n, self.k).random_element()])

    def rank(self, sub):
        """
        Returns the rank of sub as a subset of s of size k.

        EXAMPLES::

            sage: Subsets(3,2).rank([1,2])
            0
            sage: Subsets([2,3,4],2).rank([3,4])
            2
            sage: Subsets([2,3,4],2).rank([2])
            sage: Subsets([2,3,4],4).rank([2,3,4,5])
        """
        subset = Set(sub)
        lset = __builtin__.list(self.s)
        lsubset = __builtin__.list(subset)

        try:
            index_list = sorted(map(lambda x: lset.index(x), lsubset))
        except ValueError:
            return None

        n = len(self.s)
        r = 0

        if self.k not in range(len(self.s)+1):
            return None
        elif self.k != len(subset):
            return None
        else:
            return choose_nk.rank(index_list,n)


    def unrank(self, r):
        """
        Returns the subset of s that has rank k.

        EXAMPLES::

            sage: Subsets(3,2).unrank(0)
            {1, 2}
            sage: Subsets([2,4,5],2).unrank(0)
            {2, 4}
        """

        lset = __builtin__.list(self.s)
        n = len(lset)

        if self.k not in range(len(self.s)+1):
            return None
        elif r >= self.cardinality() or r < 0:
            return None
        else:
            return Set([lset[i] for i in choose_nk.from_rank(r, n, self.k)])

    def _an_element_(self):
        """
        Returns an example of subset.

        EXAMPLES::

            sage: Subsets(0,0)._an_element_()
            {}
            sage: Subsets(3,2)._an_element_()
            {1, 3}
            sage: Subsets([2,4,5],2)._an_element_()
            {2, 5}
        """
        return self.unrank(self.cardinality() // 2)

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: S32 = Subsets(3,2); S32([1,2])
            {1, 2}
            sage: S32([0,1,2])
            Traceback (most recent call last):
            ...
            ValueError: [0, 1, 2] not in Subsets of {1, 2, 3} of size 2
        """
        return Set(x)

    element_class = Set_object_enumerated


class SubMultiset_s(CombinatorialClass):
    """
    The combinatorial class of the sub multisets of s.

    EXAMPLES::

        sage: S = Subsets([1,2,2,3], submultiset=True)
        sage: S._s
        [1, 2, 2, 3]

    The positions of the unique elements in s are stored in::

        sage: S._indices
        [0, 1, 3]

    and their multiplicities in::

        sage: S._multiplicities
        [1, 2, 1]
        sage: Subsets([1,2,3,3], submultiset=True).cardinality()
        12
        sage: TestSuite(S).run()
    """
    def __init__(self, s):
        """
        Constructs the combinatorial class of the sub multisets of s.

        EXAMPLES::

            sage: S = Subsets([1,2,2,3], submultiset=True)
            sage: Subsets([1,2,3,3], submultiset=True).cardinality()
            12
        """
        CombinatorialClass.__init__(self, category=FiniteEnumeratedSets())

        s = sorted(list(s))
        indices = list(sorted(Set([s.index(a) for a in s])))
        multiplicities = [len([a for a in s if a == s[i]])
                          for i in indices]

        self._s = s
        self._indices = indices
        self._multiplicities = multiplicities

    def __repr__(self):
        """
        TESTS::

            sage: S = Subsets([1, 2, 2, 3], submultiset=True); S
            SubMultiset of [1, 2, 2, 3]
        """
        return "SubMultiset of %s"%self._s

    def __contains__(self, s):
        """
        TESTS::

            sage: S = Subsets([1,2,2,3], submultiset=True)
            sage: [] in S
            True
            sage: [1, 2, 2] in S
            True
            sage: all(i in S for i in S)
            True
            sage: [1, 2, 2, 2] in S
            False
            sage: [1, 3, 2, 2] in S
            True
            sage: [4] in S
            False
        """
        return sorted(s) in subword.Subwords(self._s)

    def __iter__(self):
        """
        Iterates through the subsets of the multiset ``self._s``.  Note
        that each subset is represented by a list of its elements rather than
        a set since we can have multiplicities (no multiset data structure yet
        in sage).

        EXAMPLES::

            sage: S = Subsets([1,2,2,3], submultiset=True)
            sage: S.list()
            [[],
             [1],
             [2],
             [3],
             [1, 2],
             [1, 3],
             [2, 2],
             [2, 3],
             [1, 2, 2],
             [1, 2, 3],
             [2, 2, 3],
             [1, 2, 2, 3]]

        """
        for k in range(len(self._s)+1):
            for s in SubMultiset_sk(self._s, k):
                yield s


class SubMultiset_sk(SubMultiset_s):
    """
    The combinatorial class of the subsets of size k of a multiset s.  Note
    that each subset is represented by a list of the elements rather than a
    set since we can have multiplicities (no multiset data structure yet in
    sage).

    EXAMPLES::

        sage: S = Subsets([1,2,3,3],2, submultiset=True)
        sage: S._k
        2
        sage: S.cardinality()
        4
        sage: S.first()
        [1, 2]
        sage: S.last()
        [3, 3]
        sage: [sub for sub in S]
        [[1, 2], [1, 3], [2, 3], [3, 3]]
        sage: TestSuite(S).run()
        """
    def __init__(self, s, k):
        """
        TEST::

            sage: S = Subsets([1,2,3,3],2, submultiset=True)
            sage: [sub for sub in S]
            [[1, 2], [1, 3], [2, 3], [3, 3]]
        """
        SubMultiset_s.__init__(self, s)
        self._k = k

    def __repr__(self):
        """
        TESTS::

            sage: S = Subsets([1, 2, 2, 3], 3, submultiset=True); S
            SubMultiset of [1, 2, 2, 3] of size 3
        """
        return "SubMultiset of %s of size %s"%(self._s, self._k)

    def __contains__(self, s):
        """
        TESTS::

            sage: S = Subsets([1,2,2,3], 2, submultiset=True)
            sage: [] in S
            False
            sage: [1, 2, 2] in S
            False
            sage: all(i in S for i in S)
            True
            sage: [2, 2] in S
            True
            sage: [1, 3] in S
            True
            sage: [4] in S
            False
            sage: [3, 3] in S
            False
        """
        return sorted(s) in subword.Subwords(self._s, self._k)

    def __iter__(self):
        """
        Iterates through the subsets of size ``self._k`` of the multiset
        ``self._s``. Note that each subset is represented by a list of the
        elements rather than a set since we can have multiplicities (no
        multiset data structure yet in sage).

        EXAMPLES::

            sage: S = Subsets([1,2,2,3],2, submultiset=True)
            sage: S.list()
            [[1, 2], [1, 3], [2, 2], [2, 3]]
        """
        from sage.combinat.integer_vector import IntegerVectors
        for iv in IntegerVectors(self._k, len(self._indices), outer=self._multiplicities):
            yield sum([ [self._s[self._indices[i]]]*iv[i] for i in range(len(iv))], [])

