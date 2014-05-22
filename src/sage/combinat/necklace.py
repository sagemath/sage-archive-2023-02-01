"""
Necklaces

The algorithm used in this file comes from

- Sawada, Joe.  "A fast algorithm to generate necklaces with fixed content", Source
  Theoretical Computer Science archive Volume 301 , Issue 1-3 (May
  2003)
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
from sage.combinat.composition import Composition
from combinat import CombinatorialClass
from sage.rings.arith import euler_phi,factorial, divisors, gcd
from sage.rings.integer import Integer
from sage.misc.misc import prod
from sage.combinat.misc import DoublyLinkedList

def Necklaces(e):
    """
    Returns the combinatorial class of necklaces with evaluation e.

    EXAMPLES::

        sage: Necklaces([2,1,1])
        Necklaces with evaluation [2, 1, 1]
        sage: Necklaces([2,1,1]).cardinality()
        3
        sage: Necklaces([2,1,1]).first()
        [1, 1, 2, 3]
        sage: Necklaces([2,1,1]).last()
        [1, 2, 1, 3]
        sage: Necklaces([2,1,1]).list()
        [[1, 1, 2, 3], [1, 1, 3, 2], [1, 2, 1, 3]]
    """
    return Necklaces_evaluation(e)

class Necklaces_evaluation(CombinatorialClass):
    def __init__(self, e):
        """
        TESTS::

            sage: N = Necklaces([2,2,2])
            sage: N == loads(dumps(N))
            True
        """
        if isinstance(e, Composition):
            self.e = e
        else:
            self.e = Composition(e)


    def __repr__(self):
        """
        TESTS::

            sage: repr(Necklaces([2,1,1]))
            'Necklaces with evaluation [2, 1, 1]'
        """
        return "Necklaces with evaluation %s"%self.e

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [2,1,2,1] in Necklaces([2,2])
            False
            sage: [1,2,1,2] in Necklaces([2,2])
            True
            sage: [1,1,2,2] in Necklaces([2,2])
            True
            sage: all([ n in Necklaces([2,1,3,1]) for n in Necklaces([2,1,3,1])])
            True
        """
        xl = list(x)
        n = sum(self.e)
        e = [0]*len(self.e)
        if len(xl) != n:
            return False

        #Check to make sure xl is a list of integers
        for i in xl:
            if not isinstance(i, (int, Integer)):
                return False
            if i <= 0:
                return False
            if i > n:
                return False
            e[i-1] += 1

        #Check to make sure the evaluation is the same
        if e != self.e:
            return False

        #Check to make sure that x is lexicographically less
        #than all of its cyclic shifts
        cyclic_shift = xl[:]
        for i in range(n - 1):
            cyclic_shift = cyclic_shift[1:] + cyclic_shift[:1]
            if cyclic_shift < xl:
                return False

        return True

    def cardinality(self):
        """
        Returns the number of integer necklaces with the evaluation e.

        EXAMPLES::

            sage: Necklaces([]).cardinality()
            0
            sage: Necklaces([2,2]).cardinality()
            2
            sage: Necklaces([2,3,2]).cardinality()
            30

        Check to make sure that the count matches up with the number of
        Lyndon words generated.

        ::

            sage: comps = [[],[2,2],[3,2,7],[4,2]]+Compositions(4).list()
            sage: ns = [ Necklaces(comp) for comp in comps]
            sage: all( [ n.cardinality() == len(n.list()) for n in ns] )
            True
        """
        evaluation = self.e
        le = list(evaluation)
        if len(le) == 0:
            return 0

        n = sum(le)

        return sum([euler_phi(j)*factorial(n/j) / prod([factorial(ni/j) for ni in evaluation]) for j in divisors(gcd(le))])/n


    def __iter__(self):
        """
        An iterator for the integer necklaces with evaluation e.

        EXAMPLES::

            sage: Necklaces([]).list()    #indirect test
            []
            sage: Necklaces([1]).list()   #indirect test
            [[1]]
            sage: Necklaces([2]).list()   #indirect test
            [[1, 1]]
            sage: Necklaces([3]).list()   #indirect test
            [[1, 1, 1]]
            sage: Necklaces([3,3]).list() #indirect test
            [[1, 1, 1, 2, 2, 2],
             [1, 1, 2, 1, 2, 2],
             [1, 1, 2, 2, 1, 2],
             [1, 2, 1, 2, 1, 2]]
            sage: Necklaces([2,1,3]).list() #indirect test
            [[1, 1, 2, 3, 3, 3],
             [1, 1, 3, 2, 3, 3],
             [1, 1, 3, 3, 2, 3],
             [1, 1, 3, 3, 3, 2],
             [1, 2, 1, 3, 3, 3],
             [1, 2, 3, 1, 3, 3],
             [1, 2, 3, 3, 1, 3],
             [1, 3, 1, 3, 2, 3],
             [1, 3, 1, 3, 3, 2],
             [1, 3, 2, 1, 3, 3]]
        """
        if self.e == []:
            return
        for z in _sfc(self.e):
            yield [x+1 for x in z]



##############################
#Fast Fixed Content Algorithm#
##############################
def _ffc(content, equality=False):
    """
    EXAMPLES::

        sage: from sage.combinat.necklace import _ffc
        sage: list(_ffc([3,3])) #necklaces
        [[0, 1, 0, 1, 0, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 0, 1, 1, 1]]
        sage: list(_ffc([3,3], equality=True)) #Lyndon words
        [[0, 0, 1, 1, 0, 1], [0, 0, 1, 0, 1, 1], [0, 0, 0, 1, 1, 1]]
    """
    e = list(content)
    a = [len(e)-1]*sum(e)
    r = [0] * sum(e)
    a[0] = 0
    e[0] -= 1
    k = len(e)

    rng_k = range(k)
    rng_k.reverse()
    dll = DoublyLinkedList(rng_k)
    if e[0] == 0:
        dll.hide(0)

    for x in _fast_fixed_content(a, e, 2, 1, k, r, 2, dll, equality=equality):
        yield x

def _fast_fixed_content(a, content, t, p, k, r, s, dll, equality=False):
    """
    EXAMPLES::

        sage: from sage.combinat.necklace import _fast_fixed_content
        sage: from sage.combinat.misc import DoublyLinkedList
        sage: e = [3,3]
        sage: a = [len(e)-1]*sum(e)
        sage: r = [0]*sum(e)
        sage: a[0] = 0
        sage: e[0] -= 1
        sage: k = len(e)
        sage: dll = DoublyLinkedList(list(reversed(range(k))))
        sage: if e[0] == 0: dll.hide(0)
        sage: list(_fast_fixed_content(a,e,2,1,k,r,2,dll))
        [[0, 1, 0, 1, 0, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 0, 1, 1, 1]]
        sage: list(_fast_fixed_content(a,e,2,1,k,r,2,dll,True))
        [[0, 0, 1, 1, 0, 1], [0, 0, 1, 0, 1, 1], [0, 0, 0, 1, 1, 1]]
    """
    n = len(a)
    if content[k-1] == n - t + 1:
        if content[k-1] == r[t-p-1]:
            if equality:
                if n == p:
                    yield a
            else:
                if n % p == 0:
                    yield a
        elif content[k-1] > r[t-p-1]:
            yield a
    elif content[0] != n-t+1:
        j = dll.head()
        sp = s
        while j != 'end' and j >= a[t-p-1]:
            #print s, j
            r[s-1] = t-s
            a[t-1] = j
            content[j] -= 1

            if content[j] == 0:
                dll.hide(j)

            if j != k-1:
                sp = t+1

            if j == a[t-p-1]:
                for x in _fast_fixed_content(a[:], content, t+1, p+0, k, r, sp, dll, equality=equality):
                    yield x
            else:
                for x in _fast_fixed_content(a[:], content, t+1, t+0, k, r, sp, dll, equality=equality):
                    yield x

            if content[j] == 0:
                dll.unhide(j)

            content[j] += 1
            j = dll.next(j)
        a[t-1] = k-1
    return


################################
# List Fixed Content Algorithm #
################################
def _lfc(content, equality=False):
    """
    EXAMPLES::

        sage: from sage.combinat.necklace import _lfc
        sage: list(_lfc([3,3])) #necklaces
        [[0, 1, 0, 1, 0, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 0, 1, 1, 1]]
        sage: list(_lfc([3,3], equality=True)) #Lyndon words
        [[0, 0, 1, 1, 0, 1], [0, 0, 1, 0, 1, 1], [0, 0, 0, 1, 1, 1]]
    """
    content = list(content)
    a = [0]*sum(content)
    content[0] -= 1
    k = len(content)

    rng_k = range(k)
    rng_k.reverse()
    dll = DoublyLinkedList(rng_k)

    if content[0] == 0:
        dll.hide(0)

    for z in _list_fixed_content(a, content, 2, 1, k, dll, equality=equality):
        yield z

def _list_fixed_content(a, content, t, p, k, dll, equality=False):
    """
    EXAMPLES::

        sage: from sage.combinat.necklace import _list_fixed_content
        sage: from sage.combinat.misc import DoublyLinkedList
        sage: e = [3,3]
        sage: a = [0]*sum(e)
        sage: e[0] -= 1
        sage: k = len(e)
        sage: dll = DoublyLinkedList(list(reversed(range(k))))
        sage: if e[0] == 0: dll.hide(0)
        sage: list(_list_fixed_content(a,e,2,1,k,dll))
        [[0, 1, 0, 1, 0, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 0, 1, 1, 1]]
        sage: list(_list_fixed_content(a,e,2,1,k,dll,True))
        [[0, 0, 1, 1, 0, 1], [0, 0, 1, 0, 1, 1], [0, 0, 0, 1, 1, 1]]
    """
    n = len(a)
    if t > n:
        if equality:
            if n == p:
                yield a
        else:
            if n % p == 0:
                yield a
    else:
        j = dll.head()
        while j != 'end' and j >= a[t-p-1]:
            a[t-1] = j
            content[j] -= 1

            if content[j] == 0:
                dll.hide(j)

            if j == a[t-p-1]:
                for z in _list_fixed_content(a[:], content[:], t+1, p+0, k, dll, equality=equality):
                    yield z
            else:
                for z in _list_fixed_content(a[:], content[:], t+1, t+0, k, dll, equality=equality):
                    yield z

            if content[j] == 0:
                dll.unhide(j)

            content[j] += 1
            j = dll.next(j)



################################
#Simple Fixed Content Algorithm#
################################
def _sfc(content, equality=False):
    """
    This function sets things up and calls _simple_fixed_content.

    EXAMPLES::

        sage: from sage.combinat.necklace import _sfc
        sage: list(_sfc([3,3])) #necklaces
        [[0, 0, 0, 1, 1, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 1, 0, 1, 0, 1]]
        sage: list(_sfc([3,3], equality=True)) #Lyndon words
        [[0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 1, 1], [0, 0, 1, 1, 0, 1]]
    """
    content = list(content)
    a = [0]*sum(content)
    content[0] -= 1
    k = len(content)
    return _simple_fixed_content(a, content, 2, 1, k, equality=equality)

def _simple_fixed_content(a, content, t, p, k, equality=False):
    """
    EXAMPLES::

        sage: from sage.combinat.necklace import _simple_fixed_content
        sage: content = [3,3]
        sage: a = [0]*sum(content)
        sage: content[0] -= 1
        sage: k = len(content); k
        2
        sage: list(_simple_fixed_content(a, content, 2, 1, k))
        [[0, 0, 0, 1, 1, 1],
         [0, 0, 1, 0, 1, 1],
         [0, 0, 1, 1, 0, 1],
         [0, 1, 0, 1, 0, 1]]
        sage: list(_simple_fixed_content(a, content, 2, 1, k, True))
        [[0, 0, 0, 1, 1, 1], [0, 0, 1, 0, 1, 1], [0, 0, 1, 1, 0, 1]]
    """
    n = len(a)
    if t > n:
        if equality:
            if n == p:
                yield a
        else:
            if n % p == 0:
                yield a
    else:
        r = range(a[t-p-1],k)
        for j in r:
            if content[j] > 0:
                a[t-1] = j
                content[j] -= 1
                if j == a[t-p-1]:
                    for z in _simple_fixed_content(a[:], content, t+1, p+0, k, equality=equality):
                        yield z
                else:
                    for z in _simple_fixed_content(a[:], content, t+1, t+0, k, equality=equality):
                        yield z
                content[j] += 1


def _lyn(w):
    """
    Returns the length of the longest prefix of w that is a Lyndon
    word.

    EXAMPLES::

        sage: import sage.combinat.necklace as necklace
        sage: necklace._lyn([0,1,1,0,0,1,2])
        3
        sage: necklace._lyn([0,0,0,1])
        4
        sage: necklace._lyn([2,1,0,0,2,2,1])
        1
    """

    p = 1
    k = max(w)+1
    for i in range(1, len(w)):
        b = w[i]
        a = w[:i]
        if b < a[i-p] or b > k-1:
            return p
        elif b == a[i-p]:
            pass
        else:
            p = i+1
    return p




