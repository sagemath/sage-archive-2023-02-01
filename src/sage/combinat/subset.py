r"""
Subsets

The set of subsets of a finite set. The set can be given as a list or a Set
or else as an integer `n` which encodes the set `\{1,2,...,n\}`.
See :class:`Subsets` for more information and examples.

AUTHORS:

- Mike Hansen: initial version

- Florent Hivert (2009/02/06): doc improvements + new methods
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2014 Vincent Delecroix <20100.delecroix@gmail.com>,
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

import sage.misc.prandom as rnd
import itertools

from sage.categories.sets_cat import EmptySetError, Sets
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

from sage.structure.parent import Parent
from sage.structure.element import Element

from sage.sets.set import Set, Set_object_enumerated
from sage.rings.arith import binomial
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
import choose_nk

ZZ_0 = ZZ.zero()

def Subsets(s, k=None, submultiset=False):
    """
    Returns the combinatorial class of the subsets of the finite set ``s``. The
    set can be given as a list, Set or any iterable convertible to a set. It can
    alternatively be given a non-negative integer `n` which encode the set
    `\{1,2,\dots,n\}` (i.e. the Sage ``range(1,s+1)``).

    A second optional parameter ``k`` can be given. In this case, Subsets
    returns the combinatorial class of subsets of ``s`` of size ``k``.

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
        sage: S.random_element()  # random
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


    And it is possible to play with subsets of subsets::

        sage: S = Subsets(3)
        sage: S2 = Subsets(S); S2
        Subsets of Subsets of {1, 2, 3}
        sage: S2.cardinality()
        256
        sage: it = iter(S2)
        sage: [it.next() for _ in xrange(8)]
        [{}, {{}}, {{1}}, {{2}}, {{3}}, {{1, 2}},  {{1, 3}}, {{2, 3}}]
        sage: S2.random_element()     # random
        {{2}, {1, 2, 3}, {}}
        sage: [S2.unrank(k) for k in xrange(256)] == S2.list()
        True

        sage: S3 = Subsets(S2)
        sage: S3.cardinality()
        115792089237316195423570985008687907853269984665640564039457584007913129639936
        sage: S3.unrank(14123091480)
        {{{1, 3}, {1, 2, 3}, {2}, {1}},
         {{2}, {1, 2, 3}, {}, {1, 2}},
         {},
         {{2}, {1, 2, 3}, {}, {3}, {1, 2}},
         {{1, 2, 3}, {}, {1}}, {{2}, {2, 3}, {}, {1, 2}}}

        sage: T = Subsets(S2, 10)
        sage: T.cardinality()
        278826214642518400
        sage: T.unrank(1441231049)
        {{{3}, {1, 2}, {}, {2, 3}, {1}, {1, 3}, ..., {{2, 3}, {}}, {{}}}
    """
    if k is not None:
        k = Integer(k)

    if isinstance(s, (int, Integer)):
        if s < 0:
            raise ValueError("s must be non-negative")
        from sage.sets.integer_range import IntegerRange
        s = IntegerRange(1,s+1)

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

class Subsets_s(Parent):
    r"""
    Subsets of a given set.

    EXAMPLES::

        sage: S = Subsets(4); S
        Subsets of {1, 2, 3, 4}
        sage: S.cardinality()
        16
        sage: Subsets(4).list()
        [{}, {1}, {2}, {3}, {4},
         {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4},
         {1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4},
         {1, 2, 3, 4}]

        sage: S = Subsets(Subsets(Subsets(GF(3)))); S
        Subsets of Subsets of Subsets of Finite Field of size 3
        sage: S.cardinality()
        115792089237316195423570985008687907853269984665640564039457584007913129639936
        sage: S.unrank(3149254230)
        {{{1, 2}, {0, 1, 2}, {0, 2}, {0, 1}},
         {{1, 2}, {}, {0, 2}, {1}, {0, 1, 2}, {2}},
         {{1, 2}, {0}}, {{1, 2}, {0, 1}, {0, 1, 2}, {1}},
         {{0, 2}, {1}}}
    """
    # TODO: Set_object_enumerated does not inherit from Element... so we set
    # directly element_class as Set_object_enumerated
    # (see also below the failed test in __init__)
    element_class = Set_object_enumerated

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
        Parent.__init__(self, category=EnumeratedSets().Finite())
        if s not in EnumeratedSets():
            from sage.misc.misc import uniq
            from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
            s = list(s)
            us = uniq(s)
            if len(us) == len(s):
                s = FiniteEnumeratedSet(s)
            else:
                s = FiniteEnumeratedSet(us)
        self._s  = s

    @property
    def _ls(self):
        r"""
        The list of elements of the underlying set.

        We try as much as possible to *not* use it.

        TESTS::

            sage: S = Subsets([1,2,3,4])
            sage: S._ls
            [1, 2, 3, 4]
        """
        return self._s.list()

    def underlying_set(self):
        r"""
        Return the set of elements.

        EXAMPLES::

            sage: Subsets(GF(13)).underlying_set()
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}
        """
        return self.element_class(self._s)

    def __eq__(self, other):
        r"""
        Equality test

        TESTS::

            sage: Subsets([0,1,2]) == Subsets([1,2,3])
            False
            sage: Subsets([0,1,2]) == Subsets([0,1,2])
            True
            sage: Subsets([0,1,2]) == Subsets([0,1,2],2)
            False
        """
        if self.__class__ != other.__class__:
            return False
        return self._s == other._s

    def __ne__(self, other):
        r"""
        Difference test

        TESTS::

            sage: Subsets([0,1,2]) != Subsets([1,2,3])
            True
            sage: Subsets([0,1,2]) != Subsets([0,1,2])
            False
            sage: Subsets([0,1,2]) != Subsets([0,1,2],2)
            True
        """
        return not self.__eq__(other)

    def _repr_(self):
        """
        TESTS::

            sage: repr(Subsets([1,2,3])) #indirect doctest
            'Subsets of {1, 2, 3}'
        """
        return "Subsets of {}".format(self._s)

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
        if value not in Sets():
            return False
        return all(v in self._s for v in value)

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
        return Integer(1) << self._s.cardinality()

    __len__ = cardinality

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
        return self.element_class([])

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
        return self.element_class(self._s)

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
        k = ZZ_0
        while k <= self._s.cardinality():
            for ss in Subsets_sk(self._s, k)._fast_iterator():
                yield self.element_class(ss)
            k += 1

    def random_element(self):
        """
        Returns a random element of the class of subsets of s (in other
        words, a random subset of s).

        EXAMPLES::

            sage: Subsets(3).random_element()           # random
            {2}
            sage: Subsets([4,5,6]).random_element()     # random
            {5}

            sage: S = Subsets(Subsets(Subsets([0,1,2])))
            sage: S.cardinality()
            115792089237316195423570985008687907853269984665640564039457584007913129639936
            sage: s = S.random_element()
            sage: s     # random
            {{{1, 2}, {2}, {0}, {1}}, {{1, 2}, {0, 1, 2}, {0, 2}, {0}, {0, 1}}, ..., {{1, 2}, {2}, {1}}, {{2}, {0, 2}, {}, {1}}}
            sage: s in S
            True
        """
        k = ZZ.random_element(0, self.cardinality())
        return self.unrank(k)

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
            sage: Subsets(3).rank([2,3,4])
            Traceback (most recent call last):
            ...
            ValueError: {2, 3, 4} is not a subset of {1, 2, 3}
        """
        if sub not in Sets():
            ssub = Set(sub)
            if len(sub) != len(ssub):
                raise ValueError("repeated elements in {}".format(sub))
            sub = ssub

        try:
            index_list = sorted(self._s.rank(x) for x in sub)
        except (ValueError,IndexError):
            raise ValueError("{} is not a subset of {}".format(
                    Set(sub), self._s))

        n = self._s.cardinality()
        r = sum(binomial(n,i) for i in xrange(len(index_list)))
        return r + choose_nk.rank(index_list,n)

    def unrank(self, r):
        """
        Returns the subset of s that has rank k.

        EXAMPLES::

            sage: Subsets(3).unrank(0)
            {}
            sage: Subsets([2,4,5]).unrank(1)
            {2}
            sage: Subsets([1,2,3]).unrank(257)
            Traceback (most recent call last):
            ...
            IndexError: index out of range

        """
        r = Integer(r)
        if r >= self.cardinality() or r < 0:
            raise IndexError("index out of range")
        else:
            k = ZZ_0
            n = self._s.cardinality()
            bin = Integer(1)
            while r >= bin:
                r -= bin
                k += 1
                bin = binomial(n,k)
            return self.element_class([self._s.unrank(i) for i in choose_nk.from_rank(r, n, k)])

    def __call__(self, el):
        r"""
        Workaround for returning non elements.

        See the extensive documentation in
        :meth:`sage.sets.finite_enumerated_set.FiniteEnumeratedSet.__call__`.

        TESTS::

            sage: Subsets(['a','b','c'])(['a','b'])  # indirect doctest
            {'a', 'b'}
        """
        if not isinstance(el, Element):
            return self._element_constructor_(el)
        else:
            return Parent.__call__(self, el)

    def _element_constructor_(self,X):
        """
        TESTS::

            sage: S3 = Subsets(3); S3([1,2]) #indirect doctest
            {1, 2}
            sage: S3([0,1,2])
            Traceback (most recent call last):
            ...
            ValueError: {0, 1, 2} not in Subsets of {1, 2, 3}
        """
        e = self.element_class(X)
        if e not in self:
            raise ValueError("{} not in {}".format(e,self))
        return e

    def an_element(self):
        """
        Returns an example of subset.

        EXAMPLES::

            sage: Subsets(0).an_element()
            {}
            sage: Subsets(3).an_element()
            {1, 2}
            sage: Subsets([2,4,5]).an_element()
            {2, 4}
        """
        return self.unrank(self.cardinality() // 2)

class Subsets_sk(Subsets_s):
    r"""
    Subsets of fixed size of a set.

    EXAMPLES::

        sage: S = Subsets([0,1,2,5,7], 3); S
        Subsets of {0, 1, 2, 5, 7} of size 3
        sage: S.cardinality()
        10
        sage: S.first(), S.last()
        ({0, 1, 2}, {2, 5, 7})
        sage: S.random_element()  # random
        {0, 5, 7}
        sage: S([0,2,7])
        {0, 2, 7}
        sage: S([0,3,5])
        Traceback (most recent call last):
        ...
        ValueError: {0, 3, 5} not in Subsets of {0, 1, 2, 5, 7} of size 3
        sage: S([0])
        Traceback (most recent call last):
        ...
        ValueError: {0} not in Subsets of {0, 1, 2, 5, 7} of size 3
    """
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
        Subsets_s.__init__(self, s)
        self._k = Integer(k)
        if self._k < 0:
            raise ValueError("the integer k (={}) should be non-negative".format(k))

    def _repr_(self):
        """
        TESTS::

            sage: repr(Subsets(3,2)) #indirect doctest
            'Subsets of {1, 2, 3} of size 2'
        """
        return Subsets_s._repr_(self) + " of size {}".format(self._k)

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
        return len(value) == self._k and Subsets_s.__contains__(self,value)

    def __eq__(self, other):
        r"""
        Equality test

        TESTS::

            sage: Subsets(5,3) == Subsets(5,3)
            True
            sage: Subsets(4,2) == Subsets(5,2) or Subsets(4,2) == Subsets(4,3)
            False
        """
        if self.__class__ != other.__class__:
            return False
        return self._s == other._s and self._k == other._k

    def __ne__(self, other):
        r"""
        Difference test

        TESTS::

            sage: Subsets(5,3) != Subsets(5,3)
            False
            sage: Subsets(4,2) != Subsets(5,2) and Subsets(4,2) != Subsets(4,3)
            True
        """
        return not self.__eq__(other)

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
        if self._k > self._s.cardinality():
            return ZZ_0
        return binomial(self._s.cardinality(), self._k)

    __len__ = cardinality

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
            Traceback (most recent call last):
            ...
            EmptySetError
        """
        if self._k < 0 or self._k > self._s.cardinality():
            raise EmptySetError
        else:
            return self.element_class(list(itertools.islice(self._s, self._k)))

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
            Traceback (most recent call last):
            ...
            EmptySetError
        """
        if self._k > self._s.cardinality():
            raise EmptySetError
        else:
            return self.element_class([i for i in itertools.islice(reversed(self._s),self._k)])

    def _fast_iterator(self):
        r"""
        Iterate through the subsets of size k if s.

        Beware that this function yield tuples and not sets. If you need sets
        use __iter__

        EXAMPLES::

            sage: list(Subsets(range(3), 2)._fast_iterator())
            [(0, 1), (0, 2), (1, 2)]
        """
        return itertools.combinations(self._s, self._k)

    def __iter__(self):
        """
        Iterates through the subsets of s of size k.

        EXAMPLES::

            sage: Subsets(Set([1,2,3]), 2).list()
            [{1, 2}, {1, 3}, {2, 3}]
            sage: Subsets([1,2,3,3], 2).list()
            [{1, 2}, {1, 3}, {2, 3}]
            sage: Subsets(3,2).list()
            [{1, 2}, {1, 3}, {2, 3}]
            sage: Subsets(3,3).list()
            [{1, 2, 3}]
        """
        return itertools.imap(self.element_class, self._fast_iterator())

    def random_element(self):
        """
        Returns a random element of the class of subsets of s of size k (in
        other words, a random subset of s of size k).

        EXAMPLES::

            sage: Subsets(3, 2).random_element()
            {1, 2}
            sage: Subsets(3,4).random_element()
            Traceback (most recent call last):
            ...
            EmptySetError
        """
        lset = self._ls

        if self._k > len(lset):
            raise EmptySetError
        else:
            return self.element_class(rnd.sample(lset, self._k))

    def rank(self, sub):
        """
        Returns the rank of sub as a subset of s of size k.

        EXAMPLES::

            sage: Subsets(3,2).rank([1,2])
            0
            sage: Subsets([2,3,4],2).rank([3,4])
            2
            sage: Subsets([2,3,4],2).rank([2])
            Traceback (most recent call last):
            ...
            ValueError: {2} is not a subset of length 2 of {2, 3, 4}
            sage: Subsets([2,3,4],4).rank([2,3,4,5])
            Traceback (most recent call last):
            ...
            ValueError: {2, 3, 4, 5} is not a subset of length 4 of {2, 3, 4}
        """
        sub = Set(sub)
        n = self._s.cardinality()

        if self._k != sub.cardinality() or self._k > n:
            raise ValueError("{} is not a subset of length {} of {}".format(
                    sub, self._k, self._s))

        try:
            index_list = sorted(self._s.rank(x) for x in sub)
        except ValueError:
            raise ValueError("{} is not a subset of length {} of {}".format(
                    sub, self._k, self._s))

        return choose_nk.rank(index_list, n)

    def unrank(self, r):
        """
        Returns the subset of s that has rank k.

        EXAMPLES::

            sage: Subsets(3,2).unrank(0)
            {1, 2}
            sage: Subsets([2,4,5],2).unrank(0)
            {2, 4}
            sage: Subsets([1,2,8],3).unrank(42)
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        lset = self._ls
        n = len(lset)

        if self._k > n or r >= self.cardinality() or r < 0:
            raise IndexError("index out of range")
        else:
            return self.element_class([lset[i] for i in choose_nk.from_rank(r, n, self._k)])

    def an_element(self):
        """
        Returns an example of subset.

        EXAMPLES::

            sage: Subsets(0,0).an_element()
            {}
            sage: Subsets(3,2).an_element()
            {1, 3}
            sage: Subsets([2,4,5],2).an_element()
            {2, 5}
        """
        return self.unrank(self.cardinality() // 2)

def dict_to_list(d):
    r"""
    Return a list whose elements are the elements of i of d repeated with
    multiplicity d[i].

    EXAMPLES::

        sage: from sage.combinat.subset import dict_to_list
        sage: dict_to_list({'a':1, 'b':3})
        ['a', 'b', 'b', 'b']
    """
    l = []
    for i,j in d.iteritems():
        l.extend([i]*j)
    return l

def list_to_dict(l):
    r"""
    Return a dictionnary whose keys are the elements of l and values are the
    multiplicity they appear in l.

    EXAMPLES::

        sage: from sage.combinat.subset import list_to_dict
        sage: list_to_dict(['a', 'b', 'b', 'b'])
        {'a': 1, 'b': 3}
    """
    d = {}
    for elt in l:
        if elt not in d:
            d[elt] = 0
        d[elt] += 1
    return d

class SubMultiset_s(Parent):
    """
    The combinatorial class of the sub multisets of s.

    EXAMPLES::

        sage: S = Subsets([1,2,2,3], submultiset=True)
        sage: S.cardinality()
        12
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
        sage: S.first()
        []
        sage: S.last()
        [1, 2, 2, 3]
    """
    # TODO: list does not inherit from Element... so we set
    # directly element_class as list
    element_class = list

    def __init__(self, s):
        """
        Constructs the combinatorial class of the sub multisets of s.

        EXAMPLES::

            sage: S = Subsets([1,2,2,3], submultiset=True)
            sage: Subsets([1,2,3,3], submultiset=True).cardinality()
            12
            sage: TestSuite(S).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())

        self._d = s
        if not isinstance(s, dict):
            self._d = list_to_dict(s)

    def _repr_(self):
        """
        TESTS::

            sage: S = Subsets([1, 2, 2, 3], submultiset=True); S #indirect doctest
            SubMultiset of [1, 2, 2, 3]
        """
        return "SubMultiset of {}".format(dict_to_list(self._d))

    def __eq__(self, other):
        r"""
        TESTS::

            sage: Subsets([1,2,2,3], submultiset=True) == Subsets([1,2,2,3], submultiset=True)
            True
            sage: Subsets([1,2,2,3], submultiset=True) == Subsets([1,2,3,3], submultiset=True)
            False
        """
        if self.__class__ != other.__class__:
            return False
        return self._d == other._d

    def __ne__(self, other):
        r"""
        TESTS::

            sage: Subsets([1,2,2,3], submultiset=True) != Subsets([1,2,2,3], submultiset=True)
            False
            sage: Subsets([1,2,2,3], submultiset=True) != Subsets([1,2,3,3], submultiset=True)
            True
        """
        return not self.__eq__(other)

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
        dd = {}
        for elt in s:
            if elt in dd:
                dd[elt] += 1
                if dd[elt] > self._d[elt]:
                    return False
            elif elt not in self._d:
                return False
            else:
                dd[elt] = 1
        return True

    def cardinality(self):
        r"""
        Return the cardinality of self

        EXAMPLES::

            sage: S = Subsets([1,1,2,3],submultiset=True)
            sage: S.cardinality()
            12
            sage: len(S.list())
            12

            sage: S = Subsets([1,1,2,2,3],submultiset=True)
            sage: S.cardinality()
            18
            sage: len(S.list())
            18

            sage: S = Subsets([1,1,1,2,2,3],submultiset=True)
            sage: S.cardinality()
            24
            sage: len(S.list())
            24
        """
        from sage.all import prod
        return Integer(prod(k+1 for k in self._d.values()))

    def random_element(self):
        r"""
        Return a random element of self with uniform law

        EXAMPLES::

            sage: S = Subsets([1,1,2,3], submultiset=True)
            sage: S.random_element()
            [2]
        """
        l = []
        for i in self._d:
            l.extend([i]*rnd.randint(0,self._d[i]))
        return l

    def generating_serie(self,variable='x'):
        r"""
        Return the serie (here a polynom) associated to the counting of the
        element of self weighted by the number of element they contain.

        EXAMPLES::

            sage: Subsets([1,1],submultiset=True).generating_serie()
            x^2 + x + 1
            sage: Subsets([1,1,2,3],submultiset=True).generating_serie()
            x^4 + 3*x^3 + 4*x^2 + 3*x + 1
            sage: Subsets([1,1,1,2,2,3,3,4],submultiset=True).generating_serie()
            x^8 + 4*x^7 + 9*x^6 + 14*x^5 + 16*x^4 + 14*x^3 + 9*x^2 + 4*x + 1

            sage: S = Subsets([1,1,1,2,2,3,3,4],submultiset=True)
            sage: S.cardinality()
            72
            sage: sum(S.generating_serie())
            72
        """
        from sage.all import prod
        R = ZZ[variable]
        return prod(R([1]*(n+1)) for n in self._d.values())

    def __iter__(self):
        """
        Iterates through the subsets of ``self``.  Note that each subset is
        represented by a list of its elements rather than a set since we can
        have multiplicities (no multiset data structure yet in sage).

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
        for k in range(sum(self._d.values())+1):
            for s in SubMultiset_sk(self._d, k):
                yield s

    def __call__(self, el):
        r"""
        Workaround for returning non elements.

        See the extensive documentation in
        :meth:`sage.sets.finite_enumerated_set.FiniteEnumeratedSet.__call__`.

        TESTS::

            sage: Subsets(['a','b','b','c'], submultiset=True)(['a','b'])  # indirect doctest
            ['a', 'b']
        """
        if not isinstance(el, Element):
            return self._element_constructor_(el)
        else:
            return Parent.__call__(self, el)

    def _element_constructor_(self,X):
        """
        TESTS::

            sage: S = Subsets(['a','b','b','c'], submultiset=True)
            sage: S(['d'])
            Traceback (most recent call last):
            ...
            ValueError: ['d'] not in SubMultiset of ['a', 'c', 'b', 'b']
        """
        e = self.element_class(X)
        if e not in self:
            raise ValueError("{} not in {}".format(e,self))
        return e



class SubMultiset_sk(SubMultiset_s):
    """
    The combinatorial class of the subsets of size k of a multiset s.  Note
    that each subset is represented by a list of the elements rather than a
    set since we can have multiplicities (no multiset data structure yet in
    sage).

    EXAMPLES::

        sage: S = Subsets([1,2,3,3],2,submultiset=True)
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
        """
    def __init__(self, s, k):
        """
        TESTS::

            sage: S = Subsets([1,2,3,3],2,submultiset=True)
            sage: [sub for sub in S]
            [[1, 2], [1, 3], [2, 3], [3, 3]]
            sage: TestSuite(S).run()
        """
        SubMultiset_s.__init__(self, s)
        self._l = dict_to_list(self._d)
        self._k = k

    def __eq__(self, other):
        r"""
        TESTS::

            sage: Subsets([1,2,2,3], submultiset=True) == Subsets([1,2,2,3], submultiset=True)
            True
            sage: Subsets([1,2,2,3], submultiset=True) == Subsets([1,2,3,3], submultiset=True)
            False
        """
        if self.__class__ != other.__class__:
            return False
        return self._d == other._d and self._k == other._k

    def generating_serie(self,variable='x'):
        r"""
        Return the serie (this case a polynom) associated to the counting of the
        element of self weighted by the number of element they contains

        EXAMPLES::

            sage: x = ZZ['x'].gen()
            sage: l = [1,1,1,1,2,2,3]
            sage: for k in xrange(len(l)):
            ....:    S = Subsets(l,k,submultiset=True)
            ....:    print S.generating_serie(x) == S.cardinality()*x**k
            True
            True
            True
            True
            True
            True
            True
        """
        x = ZZ[variable].gen()
        P = SubMultiset_s.generating_serie(self)
        return P[self._k] * (x**self._k)

    def cardinality(self):
        r"""
        Return the cardinality of self

        EXAMPLES::

            sage: S = Subsets([1,2,2,3,3,3],4,submultiset=True)
            sage: S.cardinality()
            5
            sage: len(list(S))
            5

            sage: S = Subsets([1,2,2,3,3,3],3,submultiset=True)
            sage: S.cardinality()
            6
            sage: len(list(S))
            6
        """
        return Integer(sum(1 for _ in self))

    def _repr_(self):
        """
        TESTS::

            sage: S = Subsets([1, 2, 2, 3], 3, submultiset=True)
            sage: repr(S) #indirect doctest
            'SubMultiset of [1, 2, 2, 3] of size 3'
        """
        return "{} of size {}".format(SubMultiset_s._repr_(self), self._k)

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
        return len(s) == self._k and SubMultiset_s.__contains__(self, s)

    def random_element(self):
        r"""
        Return a random submultiset of given length

        EXAMPLES::

            sage: Subsets(7,3).random_element()
            {1, 4, 7}
            sage: Subsets(7,5).random_element()
            {1, 3, 4, 5, 7}
        """
        return rnd.sample(self._l, self._k)

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
        elts = self._d.keys()
        for iv in IntegerVectors(self._k, len(self._d), outer=self._d.values()):
            yield sum([[elts[i]] * iv[i] for i in range(len(iv))], [])
