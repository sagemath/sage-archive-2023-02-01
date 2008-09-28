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
from combinat import CombinatorialClass
from sage.rings.arith import factorial
import sage.misc.prandom as rnd
from sage.combinat.misc import DoublyLinkedList


class PermutationsNK(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS:
            sage: from sage.combinat.permutation_nk import PermutationsNK
            sage: a = PermutationsNK(3,2)
            sage: a == loads(dumps(a))
            True
        """
        self._n = n
        self._k = k

    def count(self):
        """
        Returns the number of permutations of k things from a list
        of n things.

        EXAMPLES:
            sage: from sage.combinat.permutation_nk import PermutationsNK
            sage: PermutationsNK(3,2).count()
            6
            sage: PermutationsNK(5,4).count()
            120
        """
        n, k = self._n, self._k
        return factorial(n)/factorial(n-k)

    def iterator(self):
        """
        An iterator for all permutations of k thinkgs from range(n).

        EXAMPLES:
            sage: from sage.combinat.permutation_nk import PermutationsNK
            sage: [ p for p in PermutationsNK(3,2)]
            [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]]
            sage: len(PermutationsNK(5,4).list())
            120
            sage: [1, 2, 2, 0] in PermutationsNK(5,4).list()
            False
        """
        n, k = self._n, self._k
        if k == 0:
            yield []
            return

        if k>n:
            return


        range_n = range(n)
        available = range(n)
        available = DoublyLinkedList(range_n)


        L = range(k)
        L[-1] = 'begin'
        for i in range(k-1):
            available.hide(i)

        finished = False
        while not finished:
            L[-1] = available.next_value[L[-1]]
            if L[-1] != 'end':
                yield L[:]
                continue

            finished = True
            for i in reversed(range(k-1)):
                value = L[i]
                available.unhide(value)
                value = available.next_value[value]
                if value != 'end':
                    L[i] = value
                    available.hide(value)
                    value = 'begin'
                    for i in range(i+1, k-1):
                        L[i] = value = available.next_value[value]
                        available.hide(value)
                    L[-1] = available.next_value[value]
                    yield L[:]
                    finished = False
                    break

        return

    def random_element(self):
        """
        Returns a random permutation of k things from range(n).

        EXAMPLES:
            sage: from sage.combinat.permutation_nk import PermutationsNK
            sage: PermutationsNK(3,2).random_element()
            [0, 1]
        """
        n, k = self._n, self._k
        rng = range(n)
        r = rnd.sample(rng, k)

        return r
