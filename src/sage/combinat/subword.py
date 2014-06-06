r"""
Subwords

A subword of a word `w` is a word obtained by deleting the letters at some
(non necessarily adjacent) positions in `w`. It is not to be confused with the
notion of factor where one keeps adjacent positions in `w`. Sometimes it is
useful to allow repeated uses of the same letter of `w` in a "generalized"
subword. We call this a subword with repetitions.

For example:

- "bnjr" is a subword of the word "bonjour" but not a factor;

- "njo" is both a factor and a subword of the word "bonjour";

- "nr" is a subword of "bonjour";

- "rn" is not a subword of "bonjour";

- "nnu" is not a subword of "bonjour";

- "nnu" is a subword with repetitions of "bonjour";

A word can be given either as a string, as a list or as a tuple.


As repetition can occur in the initial word, the subwords of a given words is
not a set in general but an enumerated multiset!

.. TODO::

    - implement subwords with repetitions

    - implement the category of EnumeratedMultiset and inheritate from
      when needed (i.e. the initial word has repeated letters)

AUTHORS:

- Mike Hansen: initial version

- Florent Hivert (2009/02/06): doc improvements + new methods + bug fixes

- Vincent Delecroix (2011/10/03): link to itertools for faster generation,
  documentation, random generation, improvements
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


import sage.rings.arith as arith
import sage.misc.prandom as prandom
from sage.rings.integer import Integer
import itertools
from combinat import CombinatorialClass
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet


def Subwords(w, k=None):
    """
    Returns the set of subwords of w. The word w can be given by either a
    string, a list or a tuple.

    If k is specified, then it returns the combinatorial class of
    subwords of w of length k.

    EXAMPLES::

        sage: S = Subwords(['a','b','c']); S
        Subwords of ['a', 'b', 'c']
        sage: S.first()
        []
        sage: S.last()
        ['a', 'b', 'c']
        sage: S.list()
        [[], ['a'], ['b'], ['c'], ['a', 'b'], ['a', 'c'], ['b', 'c'], ['a', 'b', 'c']]

    The same example using string::

        sage: S = Subwords('abc'); S
        Subwords of abc
        sage: S.first()
        ''
        sage: S.last()
        'abc'
        sage: S.list()
        ['', 'a', 'b', 'c', 'ab', 'ac', 'bc', 'abc']

    The same example using tuple::

        sage: S = Subwords((1,2,3)); S
        Subwords of (1, 2, 3)
        sage: S.first()
        ()
        sage: S.last()
        (1, 2, 3)
        sage: S.list()
        [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]

    Using word with specified length::

        sage: S = Subwords(['a','b','c'], 2); S
        Subwords of ['a', 'b', 'c'] of length 2
        sage: S.list()
        [['a', 'b'], ['a', 'c'], ['b', 'c']]
    """
    datatype = type(w)  # 'datatype' is the type of w
    if datatype not in [str, list, tuple]:
        raise ValueError("datatype should be str, list or tuple")

    build = datatype
    # 'build' is a method to build an element with the same type as
    # the one of w.
    if datatype == str:
        build = lambda x: ''.join(x)

    if k is None:
        return Subwords_w(w, build)
    if not isinstance(k, (int, Integer)):
        raise ValueError("k should be an integer")
    if k < 0 or k > len(w):
        return FiniteEnumeratedSet([])
    return Subwords_wk(w, k, build)


class Subwords_w(CombinatorialClass):
    r"""
    Subwords of a given word
    """
    def __init__(self, w, build):
        """
        TESTS::

            sage: S = Subwords([1,2,3])
            sage: S == loads(dumps(S))
            True
            sage: TestSuite(S).run()
        """
        CombinatorialClass.__init__(self)
        self._w = w   # the word
        self._build = build  # how to build an element with same type as w

    def __reduce__(self):
        r"""
        Pickle (how to construct back the object)

        TESTS::

            sage: S = Subwords((1,2,3))
            sage: S == loads(dumps(S))
            True
            sage: S = Subwords('123')
            sage: S == loads(dumps(S))
            True
            sage: S = Subwords(('a',(1,2,3),('a','b'),'ir'))
            sage: S == loads(dumps(S))
            True
        """
        return (Subwords, (self._w,))

    def __repr__(self):
        """
        TESTS::

            sage: repr(Subwords([1,2,3])) # indirect doctest
            'Subwords of [1, 2, 3]'
        """
        return "Subwords of %s" % str(self._w)

    def __contains__(self, w):
        """
        TESTS::

            sage: [] in Subwords([1,2,3,4,3,4,4])
            True
            sage: [2,3,3,4] in Subwords([1,2,3,4,3,4,4])
            True
            sage: [5,5,3] in Subwords([1,3,3,5,4,5,3,5])
            True
            sage: [3,5,5,3] in Subwords([1,3,3,5,4,5,3,5])
            True
            sage: [3,5,5,3,4] in Subwords([1,3,3,5,4,5,3,5])
            False
            sage: [2,3,3,4] in Subwords([1,2,3,4,3,4,4])
            True
            sage: [2,3,3,1] in Subwords([1,2,3,4,3,4,4])
            False
        """
        if smallest_positions(self._w, w) != False:
            return True
        return False

    def cardinality(self):
        """
        EXAMPLES::

            sage: Subwords([1,2,3]).cardinality()
            8
        """
        return Integer(2)**len(self._w)

    def first(self):
        """
        EXAMPLES::

            sage: Subwords([1,2,3]).first()
            []
            sage: Subwords((1,2,3)).first()
            ()
            sage: Subwords('123').first()
            ''
        """
        return self._build([])

    def last(self):
        """
        EXAMPLES::

            sage: Subwords([1,2,3]).last()
            [1, 2, 3]
            sage: Subwords((1,2,3)).last()
            (1, 2, 3)
            sage: Subwords('123').last()
            '123'
        """
        return self._w

    def random_element(self):
        r"""
        Return a random subword with uniform law

        EXAMPLES::

            sage: S1 = Subwords([1,2,3,2,1,3])
            sage: S2 = Subwords([4,6,6,6,7,4,5,5])
            sage: for i in xrange(100):
            ...     w = S1.random_element()
            ...     if w in S2:
            ...         assert(w == [])
            sage: for i in xrange(100):
            ...     w = S2.random_element()
            ...     if w in S1:
            ...         assert(w == [])
        """
        return self._build(elt for elt in self._w if prandom.randint(0,1))

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: Subwords([1,2,3]).list()
            [[], [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]
            sage: Subwords((1,2,3)).list()
            [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
            sage: Subwords('123').list()
            ['', '1', '2', '3', '12', '13', '23', '123']
        """
        return itertools.chain(*[Subwords_wk(self._w,i,self._build) for i in range(0,len(self._w)+1)])

class Subwords_wk(Subwords_w):
    r"""
    Subwords with fixed length of a given word
    """
    def __init__(self, w, k, build):
        """
        TESTS::

            sage: S = Subwords([1,2,3],2)
            sage: S == loads(dumps(S))
            True
            sage: TestSuite(S).run()
        """
        Subwords_w.__init__(self,w,build)
        self._k = k

    def __reduce__(self):
        r"""
        Pickle (how to construct back the object)

        TESTS::

            sage: S = Subwords('abc',2)
            sage: S == loads(dumps(S))
            True
            sage: S = Subwords(('a',1,'45',(1,2)))
            sage: S == loads(dumps(S))
            True
        """
        return (Subwords,(self._w,self._k))

    def __repr__(self):
        """
        TESTS::

            sage: repr(Subwords([1,2,3],2))  # indirect doctest
            'Subwords of [1, 2, 3] of length 2'
        """
        return "%s of length %s" % (Subwords_w.__repr__(self), self._k)

    def __contains__(self, w):
        """
        TESTS::

            sage: [] in Subwords([1, 3, 3, 5, 4, 5, 3, 5],0)
            True
            sage: [2,3,3,4] in Subwords([1,2,3,4,3,4,4],4)
            True
            sage: [2,3,3,4] in Subwords([1,2,3,4,3,4,4],3)
            False
            sage: [5,5,3] in Subwords([1,3,3,5,4,5,3,5],3)
            True
            sage: [5,5,3] in Subwords([1,3,3,5,4,5,3,5],4)
            False
        """
        return len(w) == self._k and Subwords_w.__contains__(self,w)

    def cardinality(self):
        r"""
        Returns the number of subwords of w of length k.

        EXAMPLES::

            sage: Subwords([1,2,3], 2).cardinality()
            3
        """
        return arith.binomial(len(self._w),self._k)

    def first(self):
        r"""
        EXAMPLES::

            sage: Subwords([1,2,3],2).first()
            [1, 2]
            sage: Subwords([1,2,3],0).first()
            []
            sage: Subwords((1,2,3),2).first()
            (1, 2)
            sage: Subwords((1,2,3),0).first()
            ()
            sage: Subwords('123',2).first()
            '12'
            sage: Subwords('123',0).first()
            ''
        """
        return self._w[:self._k]

    def last(self):
        r"""
        EXAMPLES::

            sage: Subwords([1,2,3],2).last()
            [2, 3]
            sage: Subwords([1,2,3],0).last()
            []
            sage: Subwords((1,2,3),2).last()
            (2, 3)
            sage: Subwords((1,2,3),0).last()
            ()
            sage: Subwords('123',2).last()
            '23'
            sage: Subwords('123',0).last()
            ''

        TESTS::

            sage: Subwords('123', 0).last()  # trac 10534
            ''
        """
        if self._k:
            return self._w[-self._k:]
        return self.first()

    def random_element(self):
        r"""
        Return a random subword of given length with uniform law

        EXAMPLES::

            sage: S1 = Subwords([1,2,3,2,1],3)
            sage: S2 = Subwords([4,4,5,5,4,5,4,4],3)
            sage: for i in xrange(100):
            ...     w = S1.random_element()
            ...     if w in S2:
            ...         assert(w == [])
            sage: for i in xrange(100):
            ...     w = S2.random_element()
            ...     if w in S1:
            ...         assert(w == [])
        """
        sample = prandom.sample(self._w, self._k)
        if type(self._w) == list:
            return sample
        return self._build(sample)

    def __iter__(self):
        """
        EXAMPLES::

            sage: Subwords([1,2,3],2).list()
            [[1, 2], [1, 3], [2, 3]]
            sage: Subwords([1,2,3],0).list()
            [[]]
            sage: Subwords((1,2,3),2).list()
            [(1, 2), (1, 3), (2, 3)]
            sage: Subwords((1,2,3),0).list()
            [()]
            sage: Subwords('abc',2).list()
            ['ab', 'ac', 'bc']
            sage: Subwords('abc',0).list()
            ['']
        """
        if self._k > len(self._w):
            return iter([])
        iterator = itertools.combinations(self._w, self._k)
        if type(self._w) == tuple:
            return iterator
        else:
            return itertools.imap(self._build, iterator)

def smallest_positions(word, subword, pos = 0):
    """
    Returns the smallest positions for which subword appears as a
    subword of word. If pos is specified, then it returns the positions
    of the first appearance of subword starting at pos.

    If subword is not found in word, then it returns False.

    EXAMPLES::

        sage: sage.combinat.subword.smallest_positions([1,2,3,4], [2,4])
        [1, 3]
        sage: sage.combinat.subword.smallest_positions([1,2,3,4,4], [2,4])
        [1, 3]
        sage: sage.combinat.subword.smallest_positions([1,2,3,3,4,4], [3,4])
        [2, 4]
        sage: sage.combinat.subword.smallest_positions([1,2,3,3,4,4], [3,4],2)
        [2, 4]
        sage: sage.combinat.subword.smallest_positions([1,2,3,3,4,4], [3,4],3)
        [3, 4]
        sage: sage.combinat.subword.smallest_positions([1,2,3,4], [2,3])
        [1, 2]
        sage: sage.combinat.subword.smallest_positions([1,2,3,4], [5,5])
        False
        sage: sage.combinat.subword.smallest_positions([1,3,3,4,5],[3,5])
        [1, 4]
        sage: sage.combinat.subword.smallest_positions([1,3,3,5,4,5,3,5],[3,5,3])
        [1, 3, 6]
        sage: sage.combinat.subword.smallest_positions([1,3,3,5,4,5,3,5],[3,5,3],2)
        [2, 3, 6]
        sage: sage.combinat.subword.smallest_positions([1,2,3,4,3,4,4],[2,3,3,1])
        False
        sage: sage.combinat.subword.smallest_positions([1,3,3,5,4,5,3,5],[3,5,3],3)
        False

    TEST:

    We check for :trac:`5534`::

        sage: w = ["a", "b", "c", "d"]; ww = ["b", "d"]
        sage: x = sage.combinat.subword.smallest_positions(w, ww); ww
        ['b', 'd']
    """
    pos -= 1
    res = [None] * len(subword)
    for i in xrange(len(subword)):
        for j in xrange(pos+1, len(word)+1):
            if j == len(word):
                return False
            if word[j] == subword[i]:
                pos = j
                break
        if pos != j:
            return False
        res[i] = pos

    return res

