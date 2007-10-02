from sage.rings.arith import binomial
import random as rnd
from combinat import CombinatorialClass

def ChooseNK(n, k):
    return ChooseNK_nk(n,k)

class ChooseNK_nk(CombinatorialClass):
    def __init__(self, n, k):
        self.n = n
        self.k = k


    def count(self):
        """
        Returns the number of choices of k things from a list
        of n things.

        EXAMPLES:
            sage: from sage.combinat.choose_nk import ChooseNK
            sage: ChooseNK(3,2).count()
            3
            sage: ChooseNK(5,2).count()
            10
        """
        return binomial(self.n, self.k)

    def iterator(self):
        """
        An iterator for all choies of k thinkgs from range(n).

        EXAMPLES:
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
        k = self.k
        n = self.n
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


    def random(self):
        """
        Returns a random choice of k things from range(n).

        EXAMPLES:
            sage: from sage.combinat.choose_nk import ChooseNK
            sage: ChooseNK(5,2).random() #random
            [0,3]
        """
        r = rnd.sample(xrange(self.n),self.k)
        r.sort()
        return r


def rank(comb, n):
    """
    Returns the rank of comb in the subsets of range(n) of
    size k.

    The algorithm used is based on combinadics and James
    McCaffrey's MSDN article.
    See: http://en.wikipedia.org/wiki/Combinadic

    EXAMPLES:
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
    Helper function for from_rank.
    """
    w = a - 1

    while binomial(w,b) > x:
        w -= 1

    return w

def from_rank(r, n, k):
    """
    Returns the combination of rank r in the subsets of range(n)
    of size k when listed in lexicographic order.

    The algorithm used is based on combinadics and James
    McCaffrey's MSDN article.
    See: http://en.wikipedia.org/wiki/Combinadic


    EXAMPLES:
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

