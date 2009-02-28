r"""
Subwords
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
from combinat import CombinatorialClass, CombinatorialObject

#A subword of a word w is a word obtaining by deleting the letters at
#some of the positions in w.

#Example:
#    - [b,n,j,r] is a subword of [b,o,n,j,o,u,r]
#    - [n,r] is a subword of [b,o,n,j,o,u,r]
#    - [r,n] is not a subword of [b,o,n,j,o,u,r]
#    - [n,n,u] is not a subword of [b,o,n,j,o,u,r]
#    - [n,n,u] is a subword with repetitions of [b,o,n,j,o,u,r]


def Subwords(w, k=None):
    """
    Returns the combinatorial class of subwords of w.

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
            raise ValueError, "k must be between 0 and %s"%len(w)
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

    def count(self):
        """
        EXAMPLES::

            sage: Subwords([1,2,3]).count()
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

    def iterator(self):
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


    def count(self):
        r"""
        Returns the number of subwords of w of length k.

        EXAMPLES::

            sage: Subwords([1,2,3], 2).count()
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



    def iterator(self):
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
        gen = combination.Combinations(range(len(w)), k).iterator()
        return itertools.imap(lambda subword: [w[x] for x in subword], gen)


def smallest_positions(word, subword, pos = 0):
    """
    Returns the smallest positions for which subword appears as a
    subword of word. If pos is specified, then it returns the positions
    of the first appearance of subword starting at pos.

    If subword is not found in word, then it returns False

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
    """
    pos -= 1
    for i in range(len(subword)):
        for j in range(pos+1, len(word)):
            if word[j] == subword[i]:
                pos = j
                break
        if pos == -1:
            return False
        subword[i] = pos

    return subword

