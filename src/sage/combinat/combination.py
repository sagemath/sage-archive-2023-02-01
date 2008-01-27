r"""
Combinations
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

from sage.interfaces.all import gap
from sage.rings.all import QQ, RR, ZZ, Integer, binomial
from combinat import CombinatorialObject, CombinatorialClass
from choose_nk import rank, from_rank

def Combinations(mset, k=None):
    """
    Returns the combinatorial class of combinations of mset. If
    k is specified, then it returns the the combintorial class
    of combinations of mset of size k.

    The combinatorial classes correctly handle the cases where
    mset has duplicate elements.

    EXAMPLES:
        sage: C = Combinations(range(4)); C
        Combinations of [0, 1, 2, 3]
        sage: C.list()
        [[],
         [0],
         [1],
         [2],
         [3],
         [0, 1],
         [0, 2],
         [0, 3],
         [1, 2],
         [1, 3],
         [2, 3],
         [0, 1, 2],
         [0, 1, 3],
         [0, 2, 3],
         [1, 2, 3],
         [0, 1, 2, 3]]
         sage: C.count()
         16

         sage: C2 = Combinations(range(4),2); C2
         Combinations of [0, 1, 2, 3] of length 2
         sage: C2.list()
         [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
         sage: C2.count()
         6

         sage: Combinations([1,2,2,3]).list()
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



    #Check to see if everything in mset is unique
    if isinstance(mset, (int, Integer)):
        mset = range(mset)
    else:
        mset = list(mset)

    d = {}
    for i in mset:
        d[mset.index(i)] = 1

    if len(d) == len(mset):
        if k is None:
            return Combinations_set(mset)
        else:
            return Combinations_setk(mset,k)
    else:
        if k is None:
            return Combinations_mset(mset)
        else:
            return Combinations_msetk(mset,k)

class Combinations_mset(CombinatorialClass):
    def __init__(self, mset):
        """
        TESTS:
            sage: C = Combinations(range(4))
            sage: C == loads(dumps(C))
            True
        """
        self.mset = mset

    def __repr__(self):
        """
        TESTS:
            sage: repr(Combinations(range(4)))
            'Combinations of [0, 1, 2, 3]'
        """
        return "Combinations of %s"%self.mset

    def iterator(self):
        """
        TESTS:
            sage: Combinations(['a','a','b']).list() #indirect test
            [[], ['a'], ['b'], ['a', 'a'], ['a', 'b'], ['a', 'a', 'b']]
        """
        for k in range(len(self.mset)+1):
            for comb in Combinations_msetk(self.mset, k):
                yield comb

    def count(self):
        """
        TESTS:
            sage: Combinations([1,2,3]).count()
            8
            sage: Combinations(['a','a','b']).count()
            6
        """
        c = 0
        for k in range(len(self.mset)+1):
            c +=  Combinations_msetk(self.mset, k).count()
        return c

class Combinations_set(Combinations_mset):
    def iterator(self):
        """
        EXAMPLES:
            sage: Combinations([1,2,3]).list() #indirect test
            [[], [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]
        """
        for k in range(len(self.mset)+1):
            for comb in Combinations_setk(self.mset, k):
                yield comb


    def unrank(self, r):
        """
        EXAMPLES:
            sage: c = Combinations([1,2,3])
            sage: c.list() == map(c.unrank, range(c.count()))
            True
        """
        k = 0
        n = len(self.mset)
        b = binomial(n, k)
        while r >= b:
            r -= b
            k += 1
            b = binomial(n,k)

        return map(lambda i: self.mset[i], from_rank(r, n, k))


    def rank(self, x):
        """
        EXAMPLES:
            sage: c = Combinations([1,2,3])
            sage: range(c.count()) == map(c.rank, c)
            True
        """
        x = map(self.mset.index, x)
        r = 0
        n = len(self.mset)
        for i in range(len(x)):
            r += binomial(n, i)
        r += rank(x, n)
        return r

class Combinations_msetk(CombinatorialClass):
    def __init__(self, mset, k):
        """
        TESTS:
            sage: C = Combinations([1,2,3],2)
            sage: C == loads(dumps(C))
            True
        """
        self.mset = mset
        self.k = k

    def __repr__(self):
        """
        TESTS:
            sage: repr(Combinations([1,2,2,3],2))
            'Combinations of [1, 2, 2, 3] of length 2'
        """
        return "Combinations of %s of length %s"%(self.mset, self.k)

    def list(self):
        """
        Wraps GAP's Combinations.
        EXAMPLES:
            sage: Combinations(['a','a','b'],2).list()
            [['a', 'a'], ['a', 'b']]
        """
        def label(x):
            return self.mset[x]

        items = map(self.mset.index, self.mset)
        ans=eval(gap.eval("Combinations(%s,%s)"%(items,ZZ(self.k))).replace("\n",""))

        return map(lambda x: map(label, x), ans)


    def count(self):
        """
        Returns the size of combinations(mset,k).
        IMPLEMENTATION: Wraps GAP's NrCombinations.

        EXAMPLES:
            sage: mset = [1,1,2,3,4,4,5]
            sage: Combinations(mset,2).count()
            12
        """
        items = map(self.mset.index, self.mset)
        return ZZ(gap.eval("NrCombinations(%s,%s)"%(items,ZZ(self.k))))



class Combinations_setk(Combinations_msetk):
    def _iterator(self, items, len_items,  n):
        for i in range(len_items):
            v = items[i:i+1]
            if n == 1:
                yield v
            else:
                rest = items[i+1:]
                for c in self._iterator(rest, len_items-i-1,  n-1):
                    yield v + c

    def _iterator_zero(self):
        yield []

    def iterator(self):
        """
        Posted by Raymond Hettinger, 2006/03/23, to the Python Cookbook:
        http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/474124

        EXAMPLES:
            sage: Combinations([1,2,3,4,5],3).list()
            [[1, 2, 3],
             [1, 2, 4],
             [1, 2, 5],
             [1, 3, 4],
             [1, 3, 5],
             [1, 4, 5],
             [2, 3, 4],
             [2, 3, 5],
             [2, 4, 5],
             [3, 4, 5]]
        """
        if self.k == 0:
            return self._iterator_zero()
        else:
            return self._iterator(self.mset, len(self.mset), self.k)


    def list(self):
        """
        EXAMPLES:
            sage: Combinations([1,2,3,4,5],3).list()
            [[1, 2, 3],
             [1, 2, 4],
             [1, 2, 5],
             [1, 3, 4],
             [1, 3, 5],
             [1, 4, 5],
             [2, 3, 4],
             [2, 3, 5],
             [2, 4, 5],
             [3, 4, 5]]
        """
        return list(self)


    def unrank(self, r):
        """
        EXAMPLES:
            sage: c = Combinations([1,2,3], 2)
            sage: c.list() == map(c.unrank, range(c.count()))
            True
        """
        return map(lambda i: self.mset[i], from_rank(r, len(self.mset), self.k))


    def rank(self, x):
        """
        EXAMPLES:
            sage: c = Combinations([1,2,3], 2)
            sage: range(c.count()) == map(c.rank, c.list())
            True
        """
        x = map(self.mset.index, x)
        return rank(x, len(self.mset))
