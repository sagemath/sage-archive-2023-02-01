"""
Lyndon words


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

from combinat import CombinatorialClass, CombinatorialObject
from sage.combinat.composition import Composition, Compositions
from sage.rings.all import euler_phi,factorial, divisors, gcd, moebius, Integer
from sage.misc.misc import prod
from sage.combinat.misc import DoublyLinkedList
import __builtin__
import itertools
import necklace
from integer_vector import IntegerVectors
import word

def LyndonWords(e, k=None):
    """
    Returns the combinatorial class of Lyndon words.

    EXAMPLES:
      If e is an integer, then e specifies the length of the
      alphabet; k must also be specified in this case.

        sage: LW = LyndonWords(3,3); LW
        Lyndon words from an alphabet of size 3 of length 3
        sage: LW.first()
        [1, 1, 2]
        sage: LW.last()
        [2, 3, 3]
        sage: LW.random()
        [1, 1, 2]
        sage: LW.count()
        8

      If e is a (weak) composition, then it returns the class of
      Lyndon words that have evaluation e.

        sage: LyndonWords([2, 0, 1]).list()
        [[1, 1, 3]]
        sage: LyndonWords([2, 0, 1, 0, 1]).list()
        [[1, 1, 3, 5], [1, 1, 5, 3], [1, 3, 1, 5]]
        sage: LyndonWords([2, 1, 1]).list()
        [[1, 1, 2, 3], [1, 1, 3, 2], [1, 2, 1, 3]]
    """
    if isinstance(e, (int, Integer)):
        if e > 0:
            if not isinstance(k, (int, Integer)):
                raise TypeError, "k must be a non-negative integer"
            if k < 0:
                raise TypeError, "k must be a non-negative integer"
            return LyndonWords_nk(e, k)
    elif e in Compositions():
        return LyndonWords_evaluation(e)

    raise TypeError, "e must be a positive integer or a composition"

class LyndonWords_evaluation(CombinatorialClass):
    def __init__(self, e):
        """
        TESTS:
            sage: LW21 = LyndonWords([2,1]); LW21
            Lyndon words with evaluation [2, 1]
            sage: LW21 == loads(dumps(LW21))
            True
        """
        self.e = Composition(e)

    def __repr__(self):
        """
        TESTS:
            sage: repr(LyndonWords([2,1,1]))
            'Lyndon words with evaluation [2, 1, 1]'
        """
        return "Lyndon words with evaluation %s"%self.e

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: [1,2,1,2] in LyndonWords([2,2])
            False
            sage: [1,1,2,2] in LyndonWords([2,2])
            True
            sage: all([ lw in LyndonWords([2,1,3,1]) for lw in LyndonWords([2,1,3,1])])
            True
        """
        if x not in necklace.Necklaces(self.e):
            return False

        #Check to make sure that x is aperiodic
        xl = list(x)
        n = sum(self.e)
        cyclic_shift = xl[:]
        for i in range(n - 1):
            cyclic_shift = cyclic_shift[1:] + cyclic_shift[:1]
            if cyclic_shift == xl:
                return False

        return True

    def count(self):
        """
        Returns the number of Lyndon words with the
        evaluation e.

        EXAMPLES:
            sage: LyndonWords([]).count()
            0
            sage: LyndonWords([2,2]).count()
            1
            sage: LyndonWords([2,3,2]).count()
            30

          Check to make sure that the count matches up with the
          number of Lyndon words generated.

            sage: comps = [[],[2,2],[3,2,7],[4,2]]+Compositions(4).list()
            sage: lws = [ LyndonWords(comp) for comp in comps]
            sage: all( [ lw.count() == len(lw.list()) for lw in lws] )
            True

        """
        evaluation = self.e
        le = __builtin__.list(evaluation)
        if len(evaluation) == 0:
            return 0

        n = sum(evaluation)

        return sum([moebius(j)*factorial(n/j) / prod([factorial(ni/j) for ni in evaluation]) for j in divisors(gcd(le))])/n


    def iterator(self):
        """
        An iterator for the Lyndon words with evaluation e.

        EXAMPLES:
            sage: LyndonWords([1]).list()    #indirect test
            [[1]]
            sage: LyndonWords([2]).list()    #indirect test
            []
            sage: LyndonWords([3]).list()    #indirect test
            []
            sage: LyndonWords([3,1]).list()  #indirect test
            [[1, 1, 1, 2]]
            sage: LyndonWords([2,2]).list()  #indirect test
            [[1, 1, 2, 2]]
            sage: LyndonWords([1,3]).list()  #indirect test
            [[1, 2, 2, 2]]
            sage: LyndonWords([3,3]).list()  #indirect test
            [[1, 1, 1, 2, 2, 2], [1, 1, 2, 1, 2, 2], [1, 1, 2, 2, 1, 2]]
            sage: LyndonWords([4,3]).list()  #indirect test
            [[1, 1, 1, 1, 2, 2, 2],
             [1, 1, 1, 2, 1, 2, 2],
             [1, 1, 1, 2, 2, 1, 2],
             [1, 1, 2, 1, 1, 2, 2],
             [1, 1, 2, 1, 2, 1, 2]]

        """
        if self.e == []:
            return

        for z in necklace._sfc(self.e, equality=True):
            yield map(_add_one, z)


def _add_one(x):
    return x+1

class LyndonWords_nk(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS:
            sage: LW23 = LyndonWords(2,3); LW23
            Lyndon words from an alphabet of size 2 of length 3
            sage: LW23== loads(dumps(LW23))
            True
        """
        self.n = Integer(n)
        self.k = Integer(k)

    def __repr__(self):
        """
        TESTS:
            sage: repr(LyndonWords(2, 3))
            'Lyndon words from an alphabet of size 2 of length 3'
        """
        return "Lyndon words from an alphabet of size %s of length %s"%(self.n, self.k)

    def __contains__(self, x):
        """
        TESTS:
            sage: LW33 = LyndonWords(3,3)
            sage: all([lw in LW33 for lw in LW33])
            True
        """

        #Get the content of x and create
        c = filter(lambda x: x!=0, word.evaluation(x))

        #Change x to the "standard form"
        d = {}
        for i in x:
            d[i] = 1
        temp = [i for i in sorted(d.keys())]
        new_x = map( lambda i: temp.index(i)+1, x )

        #Check to see if it is in the Lyndon
        return new_x in LyndonWords_evaluation(c)

    def count(self):
        """
        TESTS:
            sage: [ LyndonWords(3,i).count() for i in range(1, 11) ]
            [3, 3, 8, 18, 48, 116, 312, 810, 2184, 5880]
        """
        if self.k == 0:
            return 1
        else:
            s = 0
            for d in divisors(self.k):
                s += moebius(d)*(self.n**(self.k/d))
        return s/self.k

    def iterator(self):
        """
        TESTS:
           sage: LyndonWords(3,3).list()
           [[1, 1, 2],
            [1, 1, 3],
            [1, 2, 2],
            [1, 2, 3],
            [1, 3, 2],
            [1, 3, 3],
            [2, 2, 3],
            [2, 3, 3]]
        """
        for c in IntegerVectors(self.k, self.n):
            cf = filter(lambda x: x != 0, c)
            nonzero_indices = []
            for i in range(len(c)):
                if c[i] != 0:
                    nonzero_indices.append(i)
            for lw in LyndonWords_evaluation(cf):
                yield map(lambda x: nonzero_indices[x-1]+1, lw)


def StandardBracketedLyndonWords(n, k):
    """
    Returns the combinatorial class of standard bracketed Lyndon
    words from [1, ..., n] of length k.  These are in one to one
    correspondence with the Lyndon words and form a basis for
    the subspace of degree k of the free Lie algebra of rank n.

    EXAMPLES:
        sage: SBLW33 = StandardBracketedLyndonWords(3,3); SBLW33
        Standard bracketed Lyndon words from an alphabet of size 3 of length 3
        sage: SBLW33.first()
        [1, [1, 2]]
        sage: SBLW33.last()
        [[2, 3], 3]
        sage: SBLW33.count()
        8
        sage: SBLW33.random()
        [1, [1, 2]]
    """
    return StandardBracketedLyndonWords_nk(n,k)

class StandardBracketedLyndonWords_nk(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS:
            sage: SBLW = StandardBracketedLyndonWords(3, 2)
            sage: SBLW == loads(dumps(SBLW))
            True
        """
        self.n = Integer(n)
        self.k = Integer(k)

    def __repr__(self):
        """
        TESTS:
            sage: repr(StandardBracketedLyndonWords(3, 3))
            'Standard bracketed Lyndon words from an alphabet of size 3 of length 3'
        """
        return "Standard bracketed Lyndon words from an alphabet of size %s of length %s"%(self.n, self.k)

    def count(self):
        """
        EXAMPLES:
            sage: StandardBracketedLyndonWords(3, 3).count()
            8
            sage: StandardBracketedLyndonWords(3, 4).count()
            18
        """
        return LyndonWords_nk(self.n,self.k).count()

    def iterator(self):
        """
        EXAMPLES:
            sage: StandardBracketedLyndonWords(3, 3).list()
            [[1, [1, 2]],
             [1, [1, 3]],
             [[1, 2], 2],
             [1, [2, 3]],
             [[1, 3], 2],
             [[1, 3], 3],
             [2, [2, 3]],
             [[2, 3], 3]]
        """
        for lw in LyndonWords_nk(self.n, self.k):
            yield standard_bracketing(lw)


def standard_bracketing(lw):
    """
    Returns the standard bracketing of a Lyndon word lw.

    EXAMPLES:
        sage: import sage.combinat.lyndon_word as lyndon_word
        sage: map( lyndon_word.standard_bracketing, LyndonWords(3,3) )
        [[1, [1, 2]],
         [1, [1, 3]],
         [[1, 2], 2],
         [1, [2, 3]],
         [[1, 3], 2],
         [[1, 3], 3],
         [2, [2, 3]],
         [[2, 3], 3]]
    """
    n = max(lw)
    k = len(lw)

    if len(lw) == 1:
        return lw[0]

    for i in range(1,len(lw)):
        if lw[i:] in LyndonWords(n,k):
            return [ standard_bracketing( lw[:i] ), standard_bracketing(lw[i:]) ]

