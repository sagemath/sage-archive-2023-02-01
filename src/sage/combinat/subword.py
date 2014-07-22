r"""
Subwords

A subword of a word $w$ is a word obtained by deleting the letters at some
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

A word can be given either as a string or as a list.

AUTHORS:

- Mike Hansen: initial version

- Florent Hivert (2009/02/06): doc improvements + new methods + bug fixes


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

import sage.combinat.combination as combination
from sage.rings.arith import factorial
import itertools
from combinat import CombinatorialClass


def Subwords(w, k=None):
    """
    Returns the combinatorial class of subwords of w. The word w can be given
    by either a string or a list.

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

    ::

        sage: S = Subwords(['a','b','c'], 2); S
        Subwords of ['a', 'b', 'c'] of length 2
        sage: S.list()
        [['a', 'b'], ['a', 'c'], ['b', 'c']]
    """
    if k == None:
        return Subwords_w(w)
    else:
        if k not in range(0, len(w)+1):
            raise ValueError("k must be between 0 and %s"%len(w))
        else:
            return Subwords_wk(w,k)


class Subwords_w(CombinatorialClass):
    def __init__(self, w):
        """
        TESTS::

            sage: S = Subwords([1,2,3])
            sage: S == loads(dumps(S))
            True
        """
        self.w = w

    def __repr__(self):
        """
        TESTS::

            sage: repr(Subwords([1,2,3]))
            'Subwords of [1, 2, 3]'
        """
        return "Subwords of %s"%self.w

    def __contains__(self, w):
        """
        TESTS::

            sage: [] in Subwords([1,2,3,4,3,4,4])
            True
            sage: [2,3,3,4] in Subwords([1,2,3,4,3,4,4])
            True
            sage: [5, 5, 3] in Subwords([1, 3, 3, 5, 4, 5, 3, 5])
            True
            sage: [3, 5, 5, 3] in Subwords([1, 3, 3, 5, 4, 5, 3, 5])
            True
            sage: [3, 5, 5, 3, 4] in Subwords([1, 3, 3, 5, 4, 5, 3, 5])
            False
            sage: [2,3,3,4] in Subwords([1,2,3,4,3,4,4])
            True
            sage: [2,3,3,1] in Subwords([1,2,3,4,3,4,4])
            False
        """
        if smallest_positions(self.w, w) != False:
            return True
        return False

    def cardinality(self):
        """
        EXAMPLES::

            sage: Subwords([1,2,3]).cardinality()
            8
        """
        return 2**len(self.w)

    def first(self):
        """
        EXAMPLES::

            sage: Subwords([1,2,3]).first()
            []
        """
        return []

    def last(self):
        """
        EXAMPLES::

            sage: Subwords([1,2,3]).last()
            [1, 2, 3]
        """
        return self.w

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: [sw for sw in Subwords([1,2,3])]
            [[], [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]
        """
        #Case 1: recursively build a generator for all the subwords
        #        length k and concatenate them together
        w = self.w
        return itertools.chain(*[Subwords_wk(w,i) for i in range(0,len(w)+1)])


class Subwords_wk(CombinatorialClass):
    def __init__(self, w, k):
        """
        TESTS::

            sage: S = Subwords([1,2,3],2)
            sage: S == loads(dumps(S))
            True
        """
        self.w = w
        self.k = k

    def __repr__(self):
        """
        TESTS::

            sage: repr(Subwords([1,2,3],2))
            'Subwords of [1, 2, 3] of length 2'
        """
        return "Subwords of %s of length %s"%(self.w, self.k)

    def __contains__(self, w):
        """
        TESTS::

            sage: [] in Subwords([1, 3, 3, 5, 4, 5, 3, 5],0)
            True
            sage: [2,3,3,4] in Subwords([1,2,3,4,3,4,4],4)
            True
            sage: [2,3,3,4] in Subwords([1,2,3,4,3,4,4],3)
            False
            sage: [5, 5, 3] in Subwords([1, 3, 3, 5, 4, 5, 3, 5],3)
            True
            sage: [5, 5, 3] in Subwords([1, 3, 3, 5, 4, 5, 3, 5],4)
            False
        """
        if len(w) != self.k:
            return False
        if smallest_positions(self.w, w) != False:
            return True
        return False

    def cardinality(self):
        r"""
        Returns the number of subwords of w of length k.

        EXAMPLES::

            sage: Subwords([1,2,3], 2).cardinality()
            3
        """
        w = self.w
        k = self.k
        return factorial(len(w))/(factorial(k)*factorial(len(w)-k))


    def first(self):
        r"""
        EXAMPLES::

            sage: Subwords([1,2,3],2).first()
            [1, 2]
            sage: Subwords([1,2,3],0).first()
            []
        """
        return self.w[:self.k]

    def last(self):
        r"""
        EXAMPLES::

            sage: Subwords([1,2,3],2).last()
            [2, 3]
        """

        return self.w[-self.k:]



    def __iter__(self):
        """
        EXAMPLES::

            sage: [sw for sw in Subwords([1,2,3],2)]
            [[1, 2], [1, 3], [2, 3]]
            sage: [sw for sw in Subwords([1,2,3],0)]
            [[]]
        """
        w = self.w
        k = self.k
        #Case 1: k == 0
        if k == 0:
            return itertools.repeat([],1)

        #Case 2: build a generator for the subwords of length k
        gen = iter(combination.Combinations(range(len(w)), k))
        return itertools.imap(lambda subword: [w[x] for x in subword], gen)


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
    for i in range(len(subword)):
        for j in range(pos+1, len(word)+1):
            if j == len(word):
                return False
            if word[j] == subword[i]:
                pos = j
                break
        if pos != j:
            return False
        res[i] = pos

    return res

