"""
Weighted Integer Vectors

.. WARNING::

   The list(self) function in this file used the :class:`Permutation_class` class improperly, returning
   the list of, generally speaking, invalid permutations (repeated entries, including 0).
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

from combinat import CombinatorialClass
from __builtin__ import list as builtinlist
from sage.rings.integer import Integer
from sage.combinat.words.word import Word
from permutation import Permutation_class

def WeightedIntegerVectors(n, weight):
    """
    Returns the combinatorial class of integer vectors of n weighted by
    weight, that is, the nonnegative integer vectors `(v_1,\\dots,v_{length(weight)})`
    satisfying `\\sum_i v_i weight[i]==n`.

    EXAMPLES::

        sage: WeightedIntegerVectors(8, [1,1,2])
        Integer vectors of 8 weighted by [1, 1, 2]
        sage: WeightedIntegerVectors(8, [1,1,2]).first()
        [0, 0, 4]
        sage: WeightedIntegerVectors(8, [1,1,2]).last()
        [8, 0, 0]
        sage: WeightedIntegerVectors(8, [1,1,2]).cardinality()
        25
        sage: WeightedIntegerVectors(8, [1,1,2]).random_element()
        [1, 1, 3]
    """
    return WeightedIntegerVectors_nweight(n, weight)

class WeightedIntegerVectors_nweight(CombinatorialClass):
    def __init__(self, n, weight):
        """
        TESTS::

            sage: WIV = WeightedIntegerVectors(8, [1,1,2])
            sage: WIV == loads(dumps(WIV))
            True
        """
        self.n = n
        self.weight = weight

    def __repr__(self):
        """
        TESTS::

            sage: repr(WeightedIntegerVectors(8, [1,1,2]))
            'Integer vectors of 8 weighted by [1, 1, 2]'
        """
        return "Integer vectors of %s weighted by %s"%(self.n, self.weight)

    def __contains__(self, x):
        """
        TESTS::

            sage: [] in WeightedIntegerVectors(0, [])
            True
            sage: [] in WeightedIntegerVectors(1, [])
            False
            sage: [3,0,0] in WeightedIntegerVectors(6, [2,1,1])
            True
            sage: [1] in WeightedIntegerVectors(1, [1])
            True
            sage: [1] in WeightedIntegerVectors(2, [2])
            True
            sage: [2] in WeightedIntegerVectors(4, [2])
            True
            sage: [2, 0] in WeightedIntegerVectors(4, [2, 2])
            True
            sage: [2, 1] in WeightedIntegerVectors(4, [2, 2])
            False
            sage: [2, 1] in WeightedIntegerVectors(6, [2, 2])
            True
            sage: [2, 1, 0] in WeightedIntegerVectors(6, [2, 2])
            False
            sage: [0] in WeightedIntegerVectors(0, [])
            False
        """
        if not isinstance(x, builtinlist):
            return False
        if len(self.weight) != len(x):
            return False
        s = 0
        for i in range(len(x)):
            if not isinstance(x[i], (int, Integer)):
                return False
            s += x[i]*self.weight[i]
        if s != self.n:
            return False

        return True

    def _recfun(self, n, l):
        """
        EXAMPLES::

            sage: w = WeightedIntegerVectors(3, [2,1,1])
            sage: w._recfun(3, [1,1,2])
            [[0, 1, 1], [1, 0, 1], [0, 3, 0], [1, 2, 0], [2, 1, 0], [3, 0, 0]]
        """
        result = []
        w = l[-1]
        l = l[:-1]
        if l == []:
            d = int(n) / int(w)
            if n%w == 0:
                return [[d]]
            else:
                return [] #bad branch...

        for d in range(int(n)/int(w), -1, -1):
            result += [ x + [d] for x in self._recfun(n-d*w, l) ]

        return result

    def list(self):
        """
        TESTS::

            sage: WeightedIntegerVectors(7, [2,2]).list()
            []
            sage: WeightedIntegerVectors(3, [2,1,1]).list()
            [[1, 0, 1], [1, 1, 0], [0, 0, 3], [0, 1, 2], [0, 2, 1], [0, 3, 0]]

        ::

            sage: ivw = [ WeightedIntegerVectors(k, [1,1,1]) for k in range(11) ]
            sage: iv  = [ IntegerVectors(k, 3) for k in range(11) ]
            sage: all( [ sorted(iv[k].list()) == sorted(ivw[k].list()) for k in range(11) ] )
            True

        ::

            sage: ivw = [ WeightedIntegerVectors(k, [2,3,7]) for k in range(11) ]
            sage: all( [ i.cardinality() == len(i.list()) for i in ivw] )
            True
        """
        if len(self.weight) == 0:
            if self.n == 0:
                return [[]]
            else:
                return []

        perm = Word(self.weight).standard_permutation()
        l = [x for x in sorted(self.weight)]
        return [perm.action(_) for _ in self._recfun(self.n,l)]
