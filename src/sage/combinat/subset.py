r"""
Subsets
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
import sage.rings.integer
import sage.combinat.subword as subword
import sage.combinat.choose_nk as choose_nk
import sage.misc.prandom as rnd
import __builtin__
import itertools
from combinat import CombinatorialClass, CombinatorialObject



def Subsets(s, k=None):
    """
    Returns the combinatorial class of subsets of s.

    If s is a non-negative integer, it returns the subsets of
    range(1,s+1).

    If k is specified, it returns the subsets of s of size
    k.

    EXAMPLES:
        sage: S = Subsets([1,2,3]); S
        Subsets of {1, 2, 3}
        sage: S.count()
        8
        sage: S.first()
        {}
        sage: S.last()
        {1, 2, 3}
        sage: S.random()
        {2}
        sage: S.list()
        [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]

        sage: S = Subsets(3)
        sage: S.list()
        [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]

        sage: S = Subsets(3,2); S
        Subsets of {1, 2, 3} of size 2
        sage: S.list()
        [{1, 2}, {1, 3}, {2, 3}]

    """
    if isinstance(s, (int, sage.rings.integer.Integer)):
        if s < 0:
            raise ValueError, "s must be non-negative"
        set = Set(range(1,s+1))
    else:
        set = Set(s)

    if k == None:
        return Subsets_s(set)
    else:
        if isinstance(k, (int, sage.rings.integer.Integer)):
            return Subsets_sk(set, k)
        else:
            raise TypeError, "k must be an integer"




class Subsets_s(CombinatorialClass):
    def __init__(self, s):
        """
        TESTS:
            sage: S = Subsets([1,2,3])
            sage: S == loads(dumps(S))
            True
        """
        self.s = s

    def __repr__(self):
        """
        TESTS:
            sage: repr(Subsets([1,2,3]))
            'Subsets of {1, 2, 3}'
        """
        return "Subsets of %s"%self.s

    def count(self):
        r"""
        Returns the number of subsets of the set s.

        This is given by $2^{|s|}$.

        EXAMPLES:
            sage: Subsets(Set([1,2,3])).count()
            8
            sage: Subsets([1,2,3,3]).count()
            8
            sage: Subsets(3).count()
            8
        """
        return 2**len(self.s)

    def first(self):
        """
        Returns the first subset of s.  Since we aren't restricted
        to subsets of a certain size, this is always the empty
        set.

        EXAMPLES:
            sage: Subsets([1,2,3]).first()
            {}
            sage: Subsets(3).first()
            {}
        """
        return Set([])

    def last(self):
        """
        Returns the last subset of s.  Since we aren't restricted
        to subsets of a certain size, this is always the set s
        itself.

        EXAMPLES:
            sage: Subsets([1,2,3]).last()
            {1, 2, 3}
            sage: Subsets(3).last()
            {1, 2, 3}
        """
        return self.s


    def iterator(self):
        """
        An iterator for all the subsets of s.

        EXAMPLES:
            sage: [sub for sub in Subsets(Set([1,2,3]))]
            [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]
            sage: [sub for sub in Subsets([1,2,3,3])]
            [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]
            sage: [sub for sub in Subsets(3)]
            [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]
        """
        lset = __builtin__.list(self.s)
        #We use the iterator for the subwords of range(len(self.s))
        ind_set = lambda index_list: Set([lset[i] for i in index_list])
        it = itertools.imap(ind_set, subword.Subwords(range(len(lset))))
        for sub in it:
            yield sub

    def random(self):
        """
        Returns a random subset of s.

        EXAMPLES:
            sage: Subsets(3).random()
            {2}
            sage: Subsets([4,5,6]).random()
            {5}
        """
        lset = __builtin__.list(self.s)
        n = len(self.s)
        return Set(filter(lambda x: rnd.randint(0,1), lset))

    def rank(self, sub):
        """
        Returns the rank of sub as a subset of s.

        EXAMPLES:
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
            index_list = map(lambda x: lset.index(x), lsubset)
            index_list.sort()
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

        EXAMPLES:
            sage: Subsets(3).unrank(0)
            {}
            sage: Subsets([2,4,5]).unrank(1)
            {2}
            sage: s = Subsets([2,4,5])
        """

        lset = __builtin__.list(self.s)
        n = len(lset)

        if r >= self.count() or r < 0:
            return None
        else:
            for k in range(n+1):
                bin = binomial(n,k)
                if r >= bin:
                    r = r - bin
                else:
                    return Set([lset[i] for i in choose_nk.from_rank(r, n, k)])


class Subsets_sk(CombinatorialClass):
    def __init__(self, s, k):
        """
        TESTS:
            sage: S = Subsets(3,2)
            sage: S == loads(dumps(S))
            True
        """
        self.s = s
        self.k = k

    def __repr__(self):
        """
        TESTS:
           sage: repr(Subsets(3,2))
           'Subsets of {1, 2, 3} of size 2'
        """
        return "Subsets of %s of size %s"%(self.s, self.k)

    def count(self):
        """
        EXAMPLES:
            sage: Subsets([1,2,3,3], 2).count()
            3
            sage: Subsets(Set([1,2,3]), 2).count()
            3
            sage: Subsets([1,2,3], 1).count()
            3
            sage: Subsets([1,2,3], 3).count()
            1
            sage: Subsets([1,2,3], 0).count()
            1
            sage: Subsets([1,2,3], 4).count()
            0
            sage: Subsets(3,2).count()
            3
            sage: Subsets(3,4).count()
            0
        """
        if self.k not in range(len(self.s)+1):
            return 0
        else:
            return binomial(len(self.s),self.k)


    def first(self):
        """
        Returns the first subset of s of size k.

        EXAMPLES:
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

        EXAMPLES:
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




    def iterator(self):
        """
        An iterator for all the subsets of s of size k.

        EXAMPLES:
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



    def random(self):
        """
        Returns a random subset of s of size k.

        EXAMPLES:
            sage: Subsets(3, 2).random()
            {1, 2}
            sage: Subsets(3,4).random() is None
            True
        """
        lset = __builtin__.list(self.s)
        n = len(self.s)

        if self.k not in range(len(self.s)+1):
            return None
        else:
            return Set([lset[i] for i in choose_nk.ChooseNK(n, self.k).random()])

    def rank(self, sub):
        """
        Returns the rank of sub as a subset of s of size k.

        EXAMPLES:
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
            index_list = map(lambda x: lset.index(x), lsubset)
            index_list.sort()
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

        EXAMPLES:
            sage: Subsets(3,2).unrank(0)
            {1, 2}
            sage: Subsets([2,4,5],2).unrank(0)
            {2, 4}
        """

        lset = __builtin__.list(self.s)
        n = len(lset)

        if self.k not in range(len(self.s)+1):
            return None
        elif r >= self.count() or r < 0:
            return None
        else:
            return Set([lset[i] for i in choose_nk.from_rank(r, n, self.k)])


