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
from sage.rings.all import ZZ, Integer
from sage.rings.arith import binomial
from combinat import CombinatorialClass
from choose_nk import rank, from_rank
from integer_vector import IntegerVectors
from sage.misc.misc import uniq

def Combinations(mset, k=None):
    """
    Returns the combinatorial class of combinations of mset. If k is
    specified, then it returns the combinatorial class of combinations
    of mset of size k.

    The combinatorial classes correctly handle the cases where mset has
    duplicate elements.

    EXAMPLES::

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
         sage: C.cardinality()
         16

    ::

        sage: C2 = Combinations(range(4),2); C2
        Combinations of [0, 1, 2, 3] of length 2
        sage: C2.list()
        [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
        sage: C2.cardinality()
        6

    ::

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
        TESTS::

            sage: C = Combinations(range(4))
            sage: C == loads(dumps(C))
            True
        """
        self.mset = mset

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: c = Combinations(range(4))
            sage: all( i in c for i in c )
            True
            sage: [3,4] in c
            False
            sage: [0,0] in c
            False
        """
        try:
            x = list(x)
        except TypeError:
            return False

        return all(i in self.mset for i in x) and \
               len(uniq(x)) == len(x)


    def __repr__(self):
        """
        TESTS::

            sage: repr(Combinations(range(4)))
            'Combinations of [0, 1, 2, 3]'
        """
        return "Combinations of %s"%self.mset

    def __iter__(self):
        """
        TESTS::

            sage: Combinations(['a','a','b']).list() #indirect doctest
            [[], ['a'], ['b'], ['a', 'a'], ['a', 'b'], ['a', 'a', 'b']]
        """
        for k in range(len(self.mset)+1):
            for comb in Combinations_msetk(self.mset, k):
                yield comb

    def cardinality(self):
        """
        TESTS::

            sage: Combinations([1,2,3]).cardinality()
            8
            sage: Combinations(['a','a','b']).cardinality()
            6
        """
        c = 0
        for k in range(len(self.mset)+1):
            c +=  Combinations_msetk(self.mset, k).cardinality()
        return c

class Combinations_set(Combinations_mset):
    def __iter__(self):
        """
        EXAMPLES::

            sage: Combinations([1,2,3]).list() #indirect doctest
            [[], [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]
        """
        for k in range(len(self.mset)+1):
            for comb in Combinations_setk(self.mset, k):
                yield comb


    def unrank(self, r):
        """
        EXAMPLES::

            sage: c = Combinations([1,2,3])
            sage: c.list() == map(c.unrank, range(c.cardinality()))
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
        EXAMPLES::

            sage: c = Combinations([1,2,3])
            sage: range(c.cardinality()) == map(c.rank, c)
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
        TESTS::

            sage: C = Combinations([1,2,3],2)
            sage: C == loads(dumps(C))
            True
        """
        self.mset = mset
        self.k = k

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: c = Combinations(range(4),2)
            sage: all( i in c for i in c )
            True
            sage: [0,1] in c
            True
            sage: [0,1,2] in c
            False
            sage: [3,4] in c
            False
            sage: [0,0] in c
            False
        """
        try:
            x = list(x)
        except TypeError:
            return False
        return x in Combinations_mset(self.mset) and len(x) == self.k


    def __repr__(self):
        """
        TESTS::

            sage: repr(Combinations([1,2,2,3],2))
            'Combinations of [1, 2, 2, 3] of length 2'
        """
        return "Combinations of %s of length %s"%(self.mset, self.k)

    def __iter__(self):
        """
        EXAMPLES::

            sage: Combinations(['a','a','b'],2).list() # indirect doctest
            [['a', 'a'], ['a', 'b']]
        """
        items = map(self.mset.index, self.mset)
        indices = uniq(sorted(items))
        counts = [0]*len(indices)
        for i in items:
            counts[indices.index(i)] += 1
        for iv in IntegerVectors(self.k, len(indices), outer=counts):
            yield sum([[self.mset[indices[i]]]*iv[i] for i in range(len(indices))],[])

    def cardinality(self):
        """
        Returns the size of combinations(mset,k). IMPLEMENTATION: Wraps
        GAP's NrCombinations.

        EXAMPLES::

            sage: mset = [1,1,2,3,4,4,5]
            sage: Combinations(mset,2).cardinality()
            12
        """
        items = map(self.mset.index, self.mset)
        return ZZ(gap.eval("NrCombinations(%s,%s)"%(items,ZZ(self.k))))



class Combinations_setk(Combinations_msetk):
    def _iterator(self, items, len_items,  n):
        """
        An iterator for all the n-combinations of items.

        EXAMPLES::

            sage: it = Combinations([1,2,3,4],3)._iterator([1,2,3,4],4,3)
            sage: list(it)
            [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
        """
        for i in range(len_items):
            v = items[i:i+1]
            if n == 1:
                yield v
            else:
                rest = items[i+1:]
                for c in self._iterator(rest, len_items-i-1,  n-1):
                    yield v + c

    def _iterator_zero(self):
        """
        An iterator which just returns the empty list.

        EXAMPLES::

            sage: it = Combinations([1,2,3,4,5],3)._iterator_zero()
            sage: list(it)
            [[]]
        """
        yield []

    def __iter__(self):
        r"""
        Posted by Raymond Hettinger, 2006/03/23, to the Python Cookbook:
        http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/474124

        EXAMPLES::

            sage: Combinations([1,2,3,4,5],3).list() # indirect doctest
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
        EXAMPLES::

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
        EXAMPLES::

            sage: c = Combinations([1,2,3], 2)
            sage: c.list() == map(c.unrank, range(c.cardinality()))
            True
        """
        return map(lambda i: self.mset[i], from_rank(r, len(self.mset), self.k))


    def rank(self, x):
        """
        EXAMPLES::

            sage: c = Combinations([1,2,3], 2)
            sage: range(c.cardinality()) == map(c.rank, c.list())
            True
        """
        x = map(self.mset.index, x)
        return rank(x, len(self.mset))
