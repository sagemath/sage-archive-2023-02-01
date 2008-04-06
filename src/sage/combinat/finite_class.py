r"""
Finite combinatorial classes
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

def FiniteCombinatorialClass(l):
    """
    Returns the combinatorial class with elements in l.

    EXAMPLES:
        sage: F = FiniteCombinatorialClass([1,2,3])
        sage: F.list()
        [1, 2, 3]
        sage: F.count()
        3
        sage: F.random()
        1
        sage: F.first()
        1
        sage: F.last()
        3
    """
    if not isinstance(l, list):
        l = list(l)
    return FiniteCombinatorialClass_l(l)

class FiniteCombinatorialClass_l(CombinatorialClass):
    def __init__(self, l):
        """
        TESTS:
            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: F == loads(dumps(F))
            True
        """
        self.l = l

    def __repr__(self):
        """
        TESTS:
            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: repr(F)
            'Combinatorial class with elements in [1, 2, 3]'
        """
        return "Combinatorial class with elements in %s"%self.l

    def __contains__(self, x):
        """
        TESTS:
            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: 1 in F
            True
            sage: 2 in F
            True
            sage: 4 in F
            False
            sage: ZZ in F
            False
        """
        return x in self.l

    def list(self):
        """
        TESTS:
            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: F.list()
            [1, 2, 3]
        """
        return self.l

    def count(self):
        """
        TESTS:
            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: F.count()
            3
        """
        return len(self.l)

