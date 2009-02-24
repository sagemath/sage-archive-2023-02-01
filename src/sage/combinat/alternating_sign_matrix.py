r"""
Alternating Sign Matrices
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
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.all import ZZ, factorial
from sage.sets.set import Set
from sage.misc.misc import prod
import copy

def from_contre_tableau(comps):
    """
    Returns an alternating sign matrix from a contretableaux.

    TESTS::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: asm.from_contre_tableau([[1, 2, 3], [1, 2], [1]])
        [0 0 1]
        [0 1 0]
        [1 0 0]
        sage: asm.from_contre_tableau([[1, 2, 3], [2, 3], [3]])
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """
    n = len(comps)
    MS = MatrixSpace(ZZ, n)
    M = [ [0 for _ in range(n)] for _ in range(n) ]

    previous_set = Set([])

    for col in range(n-1, -1, -1):
        s = Set( comps[col] )
        for x in s - previous_set:
            M[x-1][col] = 1

        for x in previous_set - s:
            M[x-1][col] = -1

        previous_set = s

    return MS(M)

def AlternatingSignMatrices(n):
    r"""
    Returns the combinatorial class of alternating sign matrices of
    size n.

    EXAMPLES::

        sage: a2 = AlternatingSignMatrices(2); a2
        Alternating sign matrices of size 2
        sage: for a in a2: print a, "-\n"
        [0 1]
        [1 0]
        -
        [1 0]
        [0 1]
        -
    """
    return AlternatingSignMatrices_n(n)

class AlternatingSignMatrices_n(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS::

            sage: a2 = AlternatingSignMatrices(2); a2
            Alternating sign matrices of size 2
            sage: a2 == loads(dumps(a2))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS::

            sage: repr(AlternatingSignMatrices(2))
            'Alternating sign matrices of size 2'
        """
        return "Alternating sign matrices of size %s"%self.n

    def count(self):
        """
        TESTS::

            sage: [ AlternatingSignMatrices(n).count() for n in range(0, 11)]
            [1, 1, 2, 7, 42, 429, 7436, 218348, 10850216, 911835460, 129534272700]

        ::

            sage: asms = [ AlternatingSignMatrices(n) for n in range(6) ]
            sage: all( [ asm.count() == len(asm.list()) for asm in asms] )
            True
        """
        return prod( [ factorial(3*k+1)/factorial(self.n+k) for k in range(self.n)] )

    def iterator(self):
        """
        TESTS::

            sage: AlternatingSignMatrices(0).list()
            [[]]
            sage: AlternatingSignMatrices(1).list()
            [[1]]
            sage: map(list, AlternatingSignMatrices(2).list())
            [[(0, 1), (1, 0)], [(1, 0), (0, 1)]]
            sage: map(list, AlternatingSignMatrices(3).list())
            [[(0, 0, 1), (0, 1, 0), (1, 0, 0)],
             [(0, 1, 0), (0, 0, 1), (1, 0, 0)],
             [(0, 0, 1), (1, 0, 0), (0, 1, 0)],
             [(0, 1, 0), (1, -1, 1), (0, 1, 0)],
             [(0, 1, 0), (1, 0, 0), (0, 0, 1)],
             [(1, 0, 0), (0, 0, 1), (0, 1, 0)],
             [(1, 0, 0), (0, 1, 0), (0, 0, 1)]]
        """
        for z in ContreTableaux(self.n):
            yield from_contre_tableau(z)


def ContreTableaux(n):
    """
    Returns the combinatorial class of contre tableaux of size n.

    EXAMPLES::

        sage: ct4 = ContreTableaux(4); ct4
        Contre tableaux of size 4
        sage: ct4.count()
        42
        sage: ct4.first()
        [[1, 2, 3, 4], [1, 2, 3], [1, 2], [1]]
        sage: ct4.last()
        [[1, 2, 3, 4], [2, 3, 4], [3, 4], [4]]
        sage: ct4.random_element()
        [[1, 2, 3, 4], [1, 2, 3], [1, 3], [3]]
    """
    return ContreTableaux_n(n)

class ContreTableaux_n(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS::

            sage: ct2 = ContreTableaux(2); ct2
            Contre tableaux of size 2
            sage: ct2 == loads(dumps(ct2))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS::

            sage: repr(ContreTableaux(2))
            'Contre tableaux of size 2'
        """
        return "Contre tableaux of size %s"%self.n

    def count(self):
        """
        TESTS::

            sage: [ ContreTableaux(n).count() for n in range(0, 11)]
            [1, 1, 2, 7, 42, 429, 7436, 218348, 10850216, 911835460, 129534272700]
        """
        return prod( [ factorial(3*k+1)/factorial(self.n+k) for k in range(self.n)] )

    def _iterator_rec(self, i):
        """
        TESTS::

            sage: c = ContreTableaux(2)
            sage: list(c._iterator_rec(0))
            [[]]
            sage: list(c._iterator_rec(1))
            [[[1, 2]]]
            sage: list(c._iterator_rec(2))
            [[[1, 2], [1]], [[1, 2], [2]]]
        """
        if i == 0:
            yield []
        elif i == 1:
            yield [range(1, self.n+1)]
        else:
            for columns in self._iterator_rec(i-1):
                previous_column = columns[-1]
                for column in _next_column_iterator(previous_column, len(previous_column)-1):
                    yield columns + [ column ]

    def iterator(self):
        """
        TESTS::

            sage: ContreTableaux(0).list()     #indirect test
            [[]]
            sage: ContreTableaux(1).list()     #indirect test
            [[[1]]]
            sage: ContreTableaux(2).list()     #indirect test
            [[[1, 2], [1]], [[1, 2], [2]]]
            sage: ContreTableaux(3).list()     #indirect test
            [[[1, 2, 3], [1, 2], [1]],
             [[1, 2, 3], [1, 2], [2]],
             [[1, 2, 3], [1, 3], [1]],
             [[1, 2, 3], [1, 3], [2]],
             [[1, 2, 3], [1, 3], [3]],
             [[1, 2, 3], [2, 3], [2]],
             [[1, 2, 3], [2, 3], [3]]]
        """

        for z in self._iterator_rec(self.n):
            yield z



def _next_column_iterator(previous_column, height, i = None):
    """
    Returns a generator for all columbs of height height properly
    filled from row 1 to i

    TESTS::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: list(asm._next_column_iterator([1], 0))
        [[]]
        sage: list(asm._next_column_iterator([1,5],1))
        [[1], [2], [3], [4], [5]]
        sage: list(asm._next_column_iterator([1,4,5],2))
        [[1, 4], [1, 5], [2, 4], [2, 5], [3, 4], [3, 5], [4, 5]]
    """
    if i is None:
        i = height
    if i == 0:
        yield [-1]*height
    else:
        for column in _next_column_iterator(previous_column, height, i-1):
            min_value = previous_column[i-1]
            if i > 1:
                min_value = max(min_value, column[i-2]+1)
            for value in range(min_value, previous_column[i]+1):
                c = copy.copy(column)
                c[i-1] = value
                yield c

def _previous_column_iterator(column, height, max_value):
    """
    TESTS::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: list(asm._previous_column_iterator([2,3], 3, 4))
        [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
    """
    new_column = [1] + column + [ max_value ] * (height - len(column))
    return _next_column_iterator(new_column, height)

def TruncatedStaircases(n, last_column):
    """
    Returns the combinatorial class of truncated staircases of size n
    with last column last_column.

    EXAMPLES::

        sage: t4 = TruncatedStaircases(4, [2,3]); t4
        Truncated staircases of size 4 with last column [2, 3]
        sage: t4.count()
        4
        sage: t4.first()
        [[4, 3, 2, 1], [3, 2, 1], [3, 2]]
        sage: t4.list()
        [[[4, 3, 2, 1], [3, 2, 1], [3, 2]], [[4, 3, 2, 1], [4, 2, 1], [3, 2]], [[4, 3, 2, 1], [4, 3, 1], [3, 2]], [[4, 3, 2, 1], [4, 3, 2], [3, 2]]]
    """
    return TruncatedStaircases_nlastcolumn(n, last_column)

class TruncatedStaircases_nlastcolumn(CombinatorialClass):
    def __init__(self, n, last_column):
        """
        TESTS::

            sage: t4 = TruncatedStaircases(4, [2,3]); t4
            Truncated staircases of size 4 with last column [2, 3]
            sage: t4 == loads(dumps(t4))
            True
        """
        self.n = n
        self.last_column = last_column

    def __repr__(self):
        """
        TESTS::

            sage: repr(TruncatedStaircases(4, [2,3]))
            'Truncated staircases of size 4 with last column [2, 3]'
        """
        return "Truncated staircases of size %s with last column %s"%(self.n, self.last_column)

    def _iterator_rec(self, i):
        """
        TESTS::

            sage: t = TruncatedStaircases(3, [2,3])
            sage: list(t._iterator_rec(1))
            []
            sage: list(t._iterator_rec(2))
            [[[2, 3]]]
            sage: list(t._iterator_rec(3))
            [[[1, 2, 3], [2, 3]]]
        """
        if i < len(self.last_column):
            return
        elif i == len(self.last_column):
            yield [self.last_column]
        else:
            for columns in self._iterator_rec(i-1):
                previous_column = columns[0]
                for column in _previous_column_iterator(previous_column, len(previous_column)+1, self.n):
                    yield [column] + columns

    def iterator(self):
        """
        EXAMPLES:::

            sage: TruncatedStaircases(4, [2,3]).list() #indirect test
            [[[4, 3, 2, 1], [3, 2, 1], [3, 2]], [[4, 3, 2, 1], [4, 2, 1], [3, 2]], [[4, 3, 2, 1], [4, 3, 1], [3, 2]], [[4, 3, 2, 1], [4, 3, 2], [3, 2]]]
        """
        for z in self._iterator_rec(self.n):
            yield map(lambda x: list(reversed(x)), z)


