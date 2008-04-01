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
        """
        n, k = self._n, self._k
        if k == 0:
            yield []
            return

        if k>n:
            return


        range_n = range(n)
        available = range(n)
        dll = DoublyLinkedList(range_n)


        L = range(k)
        L[-1] = 'begin'
        for i in range(k-1):
            dll.hide(i)

        finished = False
        while not finished:
            L[-1] = dll.next_value[L[-1]]
            if L[-1] != 'end':
                yield L[:]
                continue

            finished = True
            for i in reversed(range(k-1)):
                value = L[i]
                dll.unhide(value)
                value = dll.next_value[value]
                if value != 'end':
                    L[i] = value
                    dll.hide(value)
                    value = 'begin'
                    for j in reversed(range(i+1, k-1)):
                        value = dll.next_value[value]
                        L[j] = value
                        dll.hide(value)
                    L[-1] = dll.next_value[value]
                    yield L[:]
                    finished = False
                    break

        return

    def random(self):
        """
        Returns a random permutation of k things from range(n).

        EXAMPLES:
            sage: from sage.combinat.permutation_nk import PermutationsNK
            sage: PermutationsNK(3,2).random()
            [0, 1]
        """
        n, k = self._n, self._k
        rng = range(n)
        r = rnd.sample(rng, k)

        return r
