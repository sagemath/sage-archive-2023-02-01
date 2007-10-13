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
from sage.rings.all import QQ, RR, ZZ
from combinat import CombinatorialObject, CombinatorialClass

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
            sage: Combinations([1,2,3]).list() #indirect test
            [[], [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]
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
            sage: repr(Combinations([1,2,3],2))
            'Combinations of [1, 2, 3] of length 2'
        """
        return "Combinations of %s of length %s"%(self.mset, self.k)

    def list(self):
        """
        Wraps GAP's Combinations.
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
