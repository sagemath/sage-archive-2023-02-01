r"""
Signed Compositions
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
import composition
import cartesian_product
from sage.rings.all import binomial, Integer
import __builtin__

def SignedCompositions(n):
    """
    Returns the combinatorial class of signed compositions of
    n.

    EXAMPLES:
        sage: SC3 = SignedCompositions(3); SC3
        Signed compositions of 3
        sage: SC3.count()
        18
        sage: len(SC3.list())
        18
        sage: SC3.first()
        [1, 1, 1]
        sage: SC3.last()
        [-3]
        sage: SC3.random()
        [1, -1, 1]
        sage: SC3.list()
        [[1, 1, 1],
         [1, 1, -1],
         [1, -1, 1],
         [1, -1, -1],
         [-1, 1, 1],
         [-1, 1, -1],
         [-1, -1, 1],
         [-1, -1, -1],
         [1, 2],
         [1, -2],
         [-1, 2],
         [-1, -2],
         [2, 1],
         [2, -1],
         [-2, 1],
         [-2, -1],
         [3],
         [-3]]
    """
    return SignedCompositions_n(n)

class SignedCompositions_n(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: SC3 = SignedCompositions(3)
            sage: SC3 == loads(dumps(SC3))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS:
            sage: repr(SignedCompositions(3))
            'Signed compositions of 3'
        """
        return "Signed compositions of %s"%self.n

    def __contains__(self, x):
        """
        TESTS:
            sage: [] in SignedCompositions(0)
            True
            sage: [0] in SignedCompositions(0)
            False
            sage: [2,1,3] in SignedCompositions(6)
            True
            sage: [-2, 1, -3] in SignedCompositions(6)
            True
        """
        if x == []:
            return True

        if isinstance(x, __builtin__.list):
            for i in range(len(x)):
                if not isinstance(x[i], (int, Integer)):
                    return False
                if x[i] == 0:
                    return False
                return True
        else:
            return False

        return sum([abs(i) for i in x]) == self.n

    def count(self):
        """
        TESTS:
            sage: SC4 = SignedCompositions(4)
            sage: SC4.count() == len(SC4.list())
            True
            sage: SignedCompositions(3).count()
            18
        """
        return sum([ binomial(self.n-1, i-1)*2**(i) for i in range(1, self.n+1)])

    def iterator(self):
        """
        TESTS:
            sage: SignedCompositions(0).list()   #indirect test
            [[]]
            sage: SignedCompositions(1).list()   #indirect test
            [[1], [-1]]
            sage: SignedCompositions(2).list()   #indirect test
            [[1, 1], [1, -1], [-1, 1], [-1, -1], [2], [-2]]
        """
        for comp in composition.Compositions(self.n):
            l = len(comp)
            a = [[1,-1] for i in range(l)]
            for sign in cartesian_product.CartesianProduct(*a):
                yield [ sign[i]*comp[i] for i in range(l)]

