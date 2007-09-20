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

from combinatorial_algebra import CombinatorialAlgebra
import permutation

def SymmetricGroupAlgebra(R,n):
    """
    Returns the symmetric group algebra of order n over R.

    EXAMPLES:
        sage: QS3 = SymmetricGroupAlgebra(QQ, 3); QS3
        Symmetric group algebra of order 3 over Rational Field
        sage: QS3(1)
        [1, 2, 3]
        sage: QS3(2)
        2*[1, 2, 3]
        sage: basis = [QS3(p) for p in Permutations(3)]
        sage: a = sum(basis); a
        [3, 1, 2] + [1, 2, 3] + [2, 3, 1] + [2, 1, 3] + [3, 2, 1] + [1, 3, 2]
        sage: a^2
        6*[3, 1, 2] + 6*[1, 2, 3] + 6*[2, 3, 1] + 6*[2, 1, 3] + 6*[3, 2, 1] + 6*[1, 3, 2]
        sage: a^2 == 6*a
        True
        sage: b = QS3([3, 1, 2])
        sage: b
        [3, 1, 2]
        sage: b*a
        [3, 1, 2] + [1, 2, 3] + [2, 3, 1] + [2, 1, 3] + [3, 2, 1] + [1, 3, 2]
        sage: b*a == a
        True
    """
    return SymmetricGroupAlgebra_n(R,n)

class SymmetricGroupAlgebra_n(CombinatorialAlgebra):
    def __init__(self, R, n):
        """
        TESTS:
            #sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            #sage: QS3 == loads(dumps(QS3))
            #True
        """
        self.n = n
        self._combinatorial_class = permutation.Permutations(n)
        self._name = "Symmetric group algebra of order %s"%self.n
        self._one = permutation.Permutation(range(1,n+1))
        self._prefix = ""
        CombinatorialAlgebra.__init__(self, R)

    def _multiply_basis(self, left, right):
        return left * right
