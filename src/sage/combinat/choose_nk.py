"""
Low-level Combinations
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
import sage.misc.prandom as rnd
from sage.rings.arith import binomial
from sage.misc.misc import uniq
from combinat import CombinatorialClass

class ChooseNK(CombinatorialClass):
    """
    Low level combinatorial class of all possible choices of k elements out of
    range(n) without repetitions.

    This is a low-level combinatorial class, with a simplistic interface by
    design. It aim at speed. It's element are returned as plain list of python
    int.

    EXAMPLES::

        sage: from sage.combinat.choose_nk import ChooseNK
        sage: c = ChooseNK(4,2)
        sage: c.first()
        [0, 1]
        sage: c.list()
        [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
        sage: type(c.list()[1])
        <type 'list'>
        sage: type(c.list()[1][1])
        <type 'int'>
    """
    def __init__(self, n, k):
        """
        TESTS::

            sage: from sage.combinat.choose_nk import ChooseNK
            sage: c52 = ChooseNK(5,2)
            sage: c52 == loads(dumps(c52))
            True
        """
        self._n = n
        self._k = k

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.choose_nk import ChooseNK
            sage: c52 = ChooseNK(5,2)
            sage: [0,1] in c52
            True
            sage: [1,1] in c52
            False
            sage: [0,1,3] in c52
            False
        """
        try:
            x = list(x)
        except TypeError:
            return False

        r = range(len(x))
        return all(i in r for i in x) and len(uniq(x)) == self._k


    def cardinality(self):
        """
        Returns the number of choices of k things from a list of n things.

        EXAMPLES::

            sage: from sage.combinat.choose_nk import ChooseNK
            sage: ChooseNK(3,2).cardinality()
            3
            sage: ChooseNK(5,2).cardinality()
            10
        """
        return binomial(self._n, self._k)

    def __iter__(self):
        """
        An iterator for all choices of k things from range(n).

        EXAMPLES::

            sage: from sage.combinat.choose_nk import ChooseNK
            sage: [c for c in ChooseNK(5,2)]
            [[0, 1],
             [0, 2],
             [0, 3],
             [0, 4],
             [1, 2],
             [1, 3],
             [1, 4],
             [2, 3],
             [2, 4],
             [3, 4]]
        """
        k = self._k
        n = self._n
        dif = 1
        if k == 0:
            yield []
            return

        if n < 1+(k-1)*dif:
            return
        else:
            subword = [ i*dif for i in range(k) ]

        yield subword[:]
        finished = False

        while not finished:
            #Find the biggest element that can be increased
            if subword[-1] < n-1:
                subword[-1] += 1
                yield subword[:]
                continue

            finished = True
            for i in reversed(range(k-1)):
                if subword[i]+dif < subword[i+1]:
                    subword[i] += 1
                    #Reset the bigger elements
                    for j in range(1,k-i):
                        subword[i+j] = subword[i]+j*dif
                    yield subword[:]
                    finished = False
                    break

        return


    def random_element(self):
        """
        Returns a random choice of k things from range(n).

        EXAMPLES::

            sage: from sage.combinat.choose_nk import ChooseNK
            sage: ChooseNK(5,2).random_element()
            [0, 2]
        """
        r = sorted(rnd.sample(xrange(self._n),self._k))
        return r

    def unrank(self, r):
        """
        EXAMPLES::

            sage: from sage.combinat.choose_nk import ChooseNK
            sage: c52 = ChooseNK(5,2)
            sage: c52.list() == map(c52.unrank, range(c52.cardinality()))
            True
        """
        if r < 0 or r >= self.cardinality():
            raise ValueError, "rank must be between 0 and %s (inclusive)"%(self.cardinality()-1)
        return from_rank(r, self._n, self._k)

    def rank(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.choose_nk import ChooseNK
            sage: c52 = ChooseNK(5,2)
            sage: range(c52.cardinality()) == map(c52.rank, c52)
            True
        """
        if len(x) != self._k:
            return

        return rank(x, self._n)


def rank(comb, n):
    """
    Returns the rank of comb in the subsets of range(n) of size k.

    The algorithm used is based on combinadics and James McCaffrey's
    MSDN article. See: http://en.wikipedia.org/wiki/Combinadic

    EXAMPLES::

        sage: import sage.combinat.choose_nk as choose_nk
        sage: choose_nk.rank([], 3)
        0
        sage: choose_nk.rank([0], 3)
        0
        sage: choose_nk.rank([1], 3)
        1
        sage: choose_nk.rank([2], 3)
        2
        sage: choose_nk.rank([0,1], 3)
        0
        sage: choose_nk.rank([0,2], 3)
        1
        sage: choose_nk.rank([1,2], 3)
        2
        sage: choose_nk.rank([0,1,2], 3)
        0
    """

    k = len(comb)
    if k > n:
        raise ValueError, "len(comb) must be <= n"

    #Generate the combinadic from the
    #combination
    w = [0]*k
    for i in range(k):
        w[i] = (n-1) - comb[i]

    #Calculate the integer that is the dual of
    #the lexicographic index of the combination
    r = k
    t = 0
    for i in range(k):
        t += binomial(w[i],r)
        r -= 1

    return binomial(n,k)-t-1



def _comb_largest(a,b,x):
    """
    Returns the largest w < a such that binomial(w,b) <= x.

    EXAMPLES::

        sage: from sage.combinat.choose_nk import _comb_largest
        sage: _comb_largest(6,3,10)
        5
        sage: _comb_largest(6,3,5)
        4
    """
    w = a - 1

    while binomial(w,b) > x:
        w -= 1

    return w

def from_rank(r, n, k):
    """
    Returns the combination of rank r in the subsets of range(n) of
    size k when listed in lexicographic order.

    The algorithm used is based on combinadics and James McCaffrey's
    MSDN article. See: http://en.wikipedia.org/wiki/Combinadic

    EXAMPLES::

        sage: import sage.combinat.choose_nk as choose_nk
        sage: choose_nk.from_rank(0,3,0)
        []
        sage: choose_nk.from_rank(0,3,1)
        [0]
        sage: choose_nk.from_rank(1,3,1)
        [1]
        sage: choose_nk.from_rank(2,3,1)
        [2]
        sage: choose_nk.from_rank(0,3,2)
        [0, 1]
        sage: choose_nk.from_rank(1,3,2)
        [0, 2]
        sage: choose_nk.from_rank(2,3,2)
        [1, 2]
        sage: choose_nk.from_rank(0,3,3)
        [0, 1, 2]
    """
    if k < 0:
        raise ValueError, "k must be > 0"
    if k > n:
        raise ValueError, "k must be <= n"

    a = n
    b = k
    x = binomial(n,k) - 1 - r # x is the 'dual' of m
    comb = [None]*k

    for i in range(k):
        comb[i] = _comb_largest(a,b,x)
        x = x - binomial(comb[i], b)
        a = comb[i]
        b = b -1

    for i in range(k):
        comb[i] = (n-1)-comb[i]

    return comb

