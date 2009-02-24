"""
Low-level splits
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

from sage.rings.arith import binomial
import sage.misc.prandom as rnd
import sage.combinat.choose_nk as choose_nk
from combinat import CombinatorialClass

def SplitNK(n, k):
    """
    Returns the combinatorial class of splits of a the set range(n)
    into a set of size k and a set of size n-k.

    EXAMPLES::

        sage: from sage.combinat.split_nk import SplitNK
        sage: S = SplitNK(5,2); S
        Splits of {0, ..., 4} into a set of size 2 and one of size 3
        sage: S.first()
        [[0, 1], [2, 3, 4]]
        sage: S.last()
        [[3, 4], [0, 1, 2]]
        sage: S.list()
        [[[0, 1], [2, 3, 4]],
         [[0, 2], [1, 3, 4]],
         [[0, 3], [1, 2, 4]],
         [[0, 4], [1, 2, 3]],
         [[1, 2], [0, 3, 4]],
         [[1, 3], [0, 2, 4]],
         [[1, 4], [0, 2, 3]],
         [[2, 3], [0, 1, 4]],
         [[2, 4], [0, 1, 3]],
         [[3, 4], [0, 1, 2]]]
    """
    return SplitNK_nk(n,k)

class SplitNK_nk(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS::

            sage: from sage.combinat.split_nk import SplitNK
            sage: S = SplitNK(5,2)
            sage: S == loads(dumps(S))
            True
        """
        self._n = n
        self._k = k

    def __repr__(self):
        """
        TESTS::

            sage: from sage.combinat.split_nk import SplitNK
            sage: repr(SplitNK(5,2))
            'Splits of {0, ..., 4} into a set of size 2 and one of size 3'
        """
        return "Splits of {0, ..., %s} into a set of size %s and one of size %s"%(self._n-1, self._k, self._n-self._k)

    def count(self):
        """
        Returns the number of choices of set partitions of range(n) into a
        set of size k and a set of size n-k.

        EXAMPLES::

            sage: from sage.combinat.split_nk import SplitNK
            sage: SplitNK(5,2).count()
            10
        """
        return binomial(self._n,self._k)

    def iterator(self):
        """
        An iterator for all set partitions of range(n) into a set of size k
        and a set of size n-k in lexicographic order.

        EXAMPLES::

            sage: from sage.combinat.split_nk import SplitNK
            sage: [c for c in SplitNK(5,2)]
            [[[0, 1], [2, 3, 4]],
             [[0, 2], [1, 3, 4]],
             [[0, 3], [1, 2, 4]],
             [[0, 4], [1, 2, 3]],
             [[1, 2], [0, 3, 4]],
             [[1, 3], [0, 2, 4]],
             [[1, 4], [0, 2, 3]],
             [[2, 3], [0, 1, 4]],
             [[2, 4], [0, 1, 3]],
             [[3, 4], [0, 1, 2]]]
        """
        range_n = range(self._n)
        for kset in choose_nk.ChooseNK(self._n,self._k):
            yield [ kset, filter(lambda x: x not in kset, range_n) ]


    def random_element(self):
        """
        Returns a random set partition of range(n) into a set of size k and
        a set of size n-k.

        EXAMPLES::

            sage: from sage.combinat.split_nk import SplitNK
            sage: SplitNK(5,2).random_element()
            [[0, 2], [1, 3, 4]]
        """
        r = rnd.sample(xrange(self._n),self._n)
        r0 = r[:self._k]
        r1 = r[self._k:]
        r0.sort()
        r1.sort()
        return [ r0, r1 ]
