r"""
Alternating Sign Matrices
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2012 Pierre Cagne <pierre.cagne@ens.fr>,
#                          Luis Serrano <luisgui.serrano@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import itertools
import copy
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.parent import Parent
from sage.rings.integer import Integer
from sage.combinat.posets.lattices import LatticePoset
from sage.misc.flatten import flatten
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.all import ZZ, factorial
from sage.misc.misc import prod
from sage.sets.set import Set

def from_monotone_triangle(triangle):
    r"""
    Return an alternating sign matrix from a monotone triangle.

    EXAMPLES::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: asm.from_monotone_triangle([[1, 2, 3], [1, 2], [1]])
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: asm.from_monotone_triangle([[1, 2, 3], [2, 3], [3]])
        [0 0 1]
        [0 1 0]
        [1 0 0]
    """
    n = len(triangle)
    matrix = []

    prev = [0]*n
    for line in triangle[::-1]:
        v = [1 if j+1 in line else 0 for j in range(n)]
        row = [a-b for (a, b) in zip(v, prev)]
        matrix.append(row)
        prev = v

    return MatrixSpace(ZZ,n)(matrix)


def to_monotone_triangle(matrix):
    r"""
    Return a monotone triangle from an alternating sign matrix.

    EXAMPLES::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: asm.to_monotone_triangle(Matrix(3,[[1, 0, 0],[0, 1, 0],[0, 0, 1]]))
        [[1, 2, 3], [1, 2], [1]]
        sage: M = Matrix(3,[[0,1, 0],[1, -1, 1],[0, 1, 0]])
        sage: asm.to_monotone_triangle(M)
        [[1, 2, 3], [1, 3], [2]]
        sage: asm.from_monotone_triangle(asm.to_monotone_triangle(M)) == M
        True
    """
    n = matrix.nrows()
    triangle = []
    prev = [0]*n
    for row in matrix:
        add_row=[a+b for (a,b) in itertools.izip(row, prev)]
        line = [i+1 for (i,val) in enumerate(add_row) if val==1]
        triangle.append(line)
        prev = add_row
    return triangle[::-1]



class AlternatingSignMatrices(Parent):
    r"""
    Factory class for the combinatorial structure of alternating sign matrices of size `n`.

    An alternating sign matrix is a square matrix of 0s, 1s and -1s
    such that the sum of each row and column is 1 and the non zero
    entries in each row and column alternate in sign.

    EXAMPLES:

    This will create an instance to manipulate the alternating sign
    matrices of size 3::

        sage: A = AlternatingSignMatrices(3)
        sage: A
        Alternating sign matrices of size 3
        sage: A.cardinality()
        7

    Notably, this implementation allows to make a lattice of it::

        sage: L = A.lattice()
        sage: L
        Finite lattice containing 7 elements
        sage: L.category()
        Category of facade finite lattice posets

    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, n, **kwds):
        r"""
        Factory pattern.

        Check properties on arguments, then call the appropriate class.

        EXAMPLES::

           sage: A = AlternatingSignMatrices(4)
           sage: type(A)
           <class 'sage.combinat.alternating_sign_matrix.AlternatingSignMatrices_n'>

        """
        assert(isinstance(n, (int, Integer)))
        return AlternatingSignMatrices_n(n, **kwds)

class AlternatingSignMatrices_n(AlternatingSignMatrices):
    r"""
    Specialization class for alternating sign matrices of size `n`
    when `n` is an integer.

    INPUT:

    - `n` -- an integer, the size of the matrices.

    - ``use_monotone_triangle`` (default: ``True``) -- if ``True``, the generation
      of the matrices uses monotone triangles, else it will use the earlier and
      now obsolete contre-tableaux implementation;
      must be ``True`` to generate a lattice (with the ``lattice`` method)
    """
    def __init__(self, n, use_monotone_triangles=True):
        r"""
        __init__ ; stock arguments

        TESTS::

            sage: import sage.combinat.alternating_sign_matrix as asm
            sage: A = asm.AlternatingSignMatrices_n(4)
            sage: A == loads(dumps(A))
            True
            sage: A == asm.AlternatingSignMatrices_n(4, use_monotone_triangles=False)
            False

        """
        self._n = n
        self._umt = use_monotone_triangles

    def __repr__(self):
        r"""
        String representation of the class.

        TESTS::

            sage: A = AlternatingSignMatrices(4);A
            Alternating sign matrices of size 4
        """
        return "Alternating sign matrices of size %s" % self._n

    def size(self):
        r"""
        Return the size of the matrices in ``self``.

        TESTS::

            sage: A = AlternatingSignMatrices(4)
            sage: A.size()
            4

        """
        return self._n

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: [AlternatingSignMatrices(n).cardinality() for n in range(0, 11)]
            [1, 1, 2, 7, 42, 429, 7436, 218348, 10850216, 911835460, 129534272700]

        """
        return prod( [ factorial(3*k+1)/factorial(self._n+k)
                       for k in range(self._n)] )

    def __iter__(self):
        r"""
        Iterator on the alternating sign matrices of size `n`.

        If defined using ``use_monotone_triangles``, this iterator
        will use the iteration on the monotone triangles. Else, it
        will use the iteration on contre-tableaux.

        TESTS::

            sage: A = AlternatingSignMatrices(4)
            sage: len(list(A))
            42
        """
        if self._umt:
            for t in MonotoneTriangles(self._n):
                yield from_monotone_triangle(t)
        else:
            for c in ContreTableaux(self._n):
                yield from_contre_tableau(c)

    def __eq__(self, other):
        r"""
        Define equality by having the same attributes.

        TESTS::

            sage: A = AlternatingSignMatrices(4)
            sage: A == loads(dumps(A))
            True

        """
        return (self._n == other._n) and (self._umt == other._umt)

    def _lattice_initializer(self):
        r"""
        Return a 2-tuple to use in argument of ``LatticePoset``.

        For more details about the cover relations, see
        ``MonotoneTriangles``. Notice that the returned matrices are
        made immutable to ensure their hashability required by
        ``LatticePoset``.

        EXAMPLES:

        Proof of the lattice property for alternating sign matrices of
        size 3::

            sage: A = AlternatingSignMatrices(3)
            sage: P = Poset(A._lattice_initializer())
            sage: P.is_lattice()
            True
        """
        assert(self._umt)
        (mts, rels) = MonotoneTriangles(self._n)._lattice_initializer()
        bij = dict((t, from_monotone_triangle(t)) for t in mts)
        for m in bij.itervalues(): m.set_immutable()

        asms, rels = bij.itervalues(), [(bij[a], bij[b]) for (a,b) in rels]
        return (asms, rels)

    def cover_relations(self):
        r"""
        Iterate on the cover relations between the alternating sign
        matrices.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: for (a,b) in A.cover_relations():
            ...     eval('a, b')
            ...
            (
            [1 0 0]  [0 1 0]
            [0 1 0]  [1 0 0]
            [0 0 1], [0 0 1]
            )
            (
            [1 0 0]  [1 0 0]
            [0 1 0]  [0 0 1]
            [0 0 1], [0 1 0]
            )
            (
            [0 1 0]  [ 0  1  0]
            [1 0 0]  [ 1 -1  1]
            [0 0 1], [ 0  1  0]
            )
            (
            [1 0 0]  [ 0  1  0]
            [0 0 1]  [ 1 -1  1]
            [0 1 0], [ 0  1  0]
            )
            (
            [ 0  1  0]  [0 0 1]
            [ 1 -1  1]  [1 0 0]
            [ 0  1  0], [0 1 0]
            )
            (
            [ 0  1  0]  [0 1 0]
            [ 1 -1  1]  [0 0 1]
            [ 0  1  0], [1 0 0]
            )
            (
            [0 0 1]  [0 0 1]
            [1 0 0]  [0 1 0]
            [0 1 0], [1 0 0]
            )
            (
            [0 1 0]  [0 0 1]
            [0 0 1]  [0 1 0]
            [1 0 0], [1 0 0]
            )

        """
        (_, rels) = self._lattice_initializer()
        return (_ for _ in rels)

    def lattice(self):
        r"""
        Return the lattice of the alternating sign matrices of size
        `n`, created by ``LatticePoset``.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: L = A.lattice()
            sage: L
            Finite lattice containing 7 elements

        """
        return LatticePoset(self._lattice_initializer(), cover_relations=True)



class MonotoneTriangles(Parent):
    r"""
    Factory class for the combinatorial structure of monotone
    triangles with `n` rows.

    A monotone triangle is a number triangle `(a_{i,j})_{1 \leq i \leq
    n , 1 \leq j \leq i}` on `\{1, \dots, n\}` such that :

    - `a_{i,j} < a_{i,j+1}`

    - `a_{i+1,j} < a_{i,j} \leq a_{i+1,j+1}`

    This notably requires that the bottom column is ``[1,...,n]``.

    EXAMPLES:

    This represents the monotone triangles with base ``[1,2,3]``::

        sage: M = MonotoneTriangles(3)
        sage: M
        Monotone triangles with 3 rows
        sage: M.cardinality()
        7

    The monotone triangles are a lattice::

        sage: M.lattice()
        Finite lattice containing 7 elements

    """
    __metaclass__ = ClasscallMetaclass
    @staticmethod
    def __classcall_private__(cls, n, **kwds):
        r"""
        Factory pattern.

        Check properties on arguments, then call the appropriate
        class.

        EXAMPLES::

            sage: M = MonotoneTriangles(4)
            sage: type(M)
            <class 'sage.combinat.alternating_sign_matrix.MonotoneTriangles_n'>

        """
        assert(isinstance(n, (int, Integer)))
        return MonotoneTriangles_n(n, **kwds)


class MonotoneTriangles_n(MonotoneTriangles):
    r"""
    Specialization class for monotone triangles with `n` rows when
    ``n`` is an integer.

    INPUT:

    - ``n`` -- integer of type ``int`` or ``Integer``, the number of rows of the monotone triangles in ``self``.

    """
    def __init__(self, n):
        r"""
        __init__ ; stock arguments

        TESTS::

            sage: import sage.combinat.alternating_sign_matrix as asm
            sage: M = asm.MonotoneTriangles_n(4)
            sage: M == loads(dumps(M))
            True

        """
        self._n = n

    def __repr__(self):
        r"""
        String representation.

        TESTS::

            sage: M = MonotoneTriangles(4)
            sage: M
            Monotone triangles with 4 rows

        """
        return "Monotone triangles with %s rows" % self._n

    def __iter__(self):
        r"""
        Iterate on the monotone triangles with `n` rows.

        TESTS::

            sage: M = MonotoneTriangles(4)
            sage: len(list(M))
            42

        """
        for z in _triangular_arrangements_base(range(1, self._n + 1)):
            yield z

    def __eq__(self, other):
        r"""
        TESTS::

            sage: M = MonotoneTriangles(4)
            sage: M == loads(dumps(M))
            True

        """
        return self._n == other._n

    def cardinality(self):
        r"""
        Cardinality of ``self``.

        EXAMPLES::

            sage: M = MonotoneTriangles(4)
            sage: M.cardinality()
            42

        """
        return prod( [ factorial(3*k+1)/factorial(self._n+k)
                       for k in range(self._n)] )

    def _lattice_initializer(self):
        r"""
        Return a 2-tuple to use in argument of ``LatticePoset``.

        This couple is composed by the set of the monotone triangles
        with `n` rows and the cover relations. Specializing this
        function allows to generate the monotone triangles just once,
        and so to speed up the computation in comparison of
        ``(list(self), self.cover_relations())``. Notice that the
        function also switch the representation of monotone triangles
        from list of list to tuple of tuple in order to make them
        hashable (required to make a poset with them).

        EXAMPLES::

            sage: M = MonotoneTriangles(3)
            sage: P = Poset(M._lattice_initializer())
            sage: P.is_lattice()
            True
        """
        # get a list of the elements and switch to a tuple
        # representation
        set_ = list(self)
        set_ = map(lambda x: tuple(map(tuple, x)), set_)

        return (set_, [(a,b) for a in set_ for b in set_ if _is_a_cover(a,b)])

    def cover_relations(self):
        r"""
        Iterate on the cover relations in the set of monotone triangles
        with `n` rows.

        EXAMPLES::

            sage: M = MonotoneTriangles(3)
            sage: for (a,b) in M.cover_relations():
            ...     eval('a, b')
            ...
            ([[1, 2, 3], [1, 2], [1]], [[1, 2, 3], [1, 2], [2]])
            ([[1, 2, 3], [1, 2], [1]], [[1, 2, 3], [1, 3], [1]])
            ([[1, 2, 3], [1, 2], [2]], [[1, 2, 3], [1, 3], [2]])
            ([[1, 2, 3], [1, 3], [1]], [[1, 2, 3], [1, 3], [2]])
            ([[1, 2, 3], [1, 3], [2]], [[1, 2, 3], [1, 3], [3]])
            ([[1, 2, 3], [1, 3], [2]], [[1, 2, 3], [2, 3], [2]])
            ([[1, 2, 3], [1, 3], [3]], [[1, 2, 3], [2, 3], [3]])
            ([[1, 2, 3], [2, 3], [2]], [[1, 2, 3], [2, 3], [3]])
        """
        set_ = list(self)
        return ((a,b) for a in set_ for b in set_ if _is_a_cover(a,b))

    def lattice(self):
        r"""
        Return the lattice of the monotone triangles with `n` rows.

        EXAMPLES::

            sage: M = MonotoneTriangles(3)
            sage: P = M.lattice()
            sage: P
            Finite lattice containing 7 elements

        """
        return LatticePoset(self._lattice_initializer(), cover_relations=True)


def _is_a_cover(mt0, mt1):
    r"""
    Define the cover relations.

    Return ``True`` if and only if the second argument is a cover of
    the first one.

    EXAMPLES::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: asm._is_a_cover([[1,2,3],[1,2],[1]], [[1,2,3],[1,3],[1]])
        True
        sage: asm._is_a_cover([[1,2,3],[1,3],[2]], [[1,2,3],[1,2],[1]])
        False

    """
    diffs = 0
    for (a,b) in itertools.izip(flatten(mt0), flatten(mt1)):
        if a != b:
            if a+1 == b:
                diffs += 1
            else:
                return False
        if diffs > 1:
            return False
    return diffs == 1


def _top_rows(row):
    r"""
    Determine the rows that can be on top of ``row``.

    EXAMPLES::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: asm._top_rows([1,3,4])
        [[1, 3], [1, 4], [2, 3], [2, 4], [3, 4]]
    """
    assert(len(row) >= 2)
    [inf, sup] = row[-2:]
    suffs = [[i] for i in xrange(inf, sup+1)]
    if len(row) == 2: return suffs
    sup = inf
    for inf in row[-3::-1]:
        tmp = []
        for i in xrange(inf, sup):
            tmp += [[i] + x for x in suffs]
        tmp += [[sup] + x for x in suffs if x[0] != sup]
        suffs = tmp
        sup = inf
    return suffs


def _triangular_arrangements_base(row):
    r"""
    Determine all triangular arrangements with base ``row``.

    Notice that used with ``range(1, n+1)``, it returns the monotone
    triangles with ``n`` rows.

    EXAMPLES::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: asm._triangular_arrangements_base([1,4])
        [[[1, 4], [1]], [[1, 4], [2]], [[1, 4], [3]], [[1, 4], [4]]]
    """
    arrs = [[row]]
    size = len(row)
    while size > 1:
        tmp = []
        for t in arrs:
            tmp += [t + [x] for x in _top_rows(t[-1])]
        arrs = tmp
        size -= 1
    return arrs


# Here are the previous implementations of the combinatorial structure
# of the alternating sign matrices. Please, consider it obsolete and
# tend to use the monotone triangles instead.

def from_contre_tableau(comps):
    r"""
    Returns an alternating sign matrix from a contre-tableau.

    EXAMPLES::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: asm.from_contre_tableau([[1, 2, 3], [1, 2], [1]])
        doctest:...: DeprecationWarning: You can use from_monotone_triangle instead.
        See http://trac.sagemath.org/12930 for details.
        [0 0 1]
        [0 1 0]
        [1 0 0]
        sage: asm.from_contre_tableau([[1, 2, 3], [2, 3], [3]])
        [1 0 0]
        [0 1 0]
        [0 0 1]
    """
    from sage.misc.superseded import deprecation
    deprecation(12930, 'You can use from_monotone_triangle instead.')
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


class ContreTableaux(Parent):
    """
    Factory class for the combinatorial class of contre tableaux of size `n`.

    EXAMPLES::

        sage: ct4 = ContreTableaux(4); ct4
        Contre tableaux of size 4
        sage: ct4.cardinality()
        42
    """
    __metaclass__ = ClasscallMetaclass
    @staticmethod
    def __classcall_private__(cls, n, **kwds):
        r"""
        Factory pattern.

        Check properties on arguments, then call the appropriate class.

        EXAMPLES::

            sage: C = ContreTableaux(4)
            sage: type(C)
            <class 'sage.combinat.alternating_sign_matrix.ContreTableaux_n'>

        """
        assert(isinstance(n, (int, Integer)))
        return ContreTableaux_n(n, **kwds)


class ContreTableaux_n(ContreTableaux):
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

    def __eq__(self, other):
        """
        TESTS::

            sage: C = ContreTableaux(4)
            sage: C == loads(dumps(C))
            True

        """
        return self.n == other.n

    def cardinality(self):
        """
        EXAMPLES::

            sage: [ ContreTableaux(n).cardinality() for n in range(0, 11)]
            [1, 1, 2, 7, 42, 429, 7436, 218348, 10850216, 911835460, 129534272700]
        """
        return prod( [ factorial(3*k+1)/factorial(self.n+k) for k in range(self.n)] )

    def _iterator_rec(self, i):
        """
        EXAMPLES::

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

    def __iter__(self):
        """
        EXAMPLES::

            sage: list(ContreTableaux(0))
            [[]]
            sage: list(ContreTableaux(1))
            [[[1]]]
            sage: list(ContreTableaux(2))
            [[[1, 2], [1]], [[1, 2], [2]]]
            sage: list(ContreTableaux(3))
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
    Returns a generator for all columns of height height properly
    filled from row 1 to ``i``

    EXAMPLES::

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
    EXAMPLES::

        sage: import sage.combinat.alternating_sign_matrix as asm
        sage: list(asm._previous_column_iterator([2,3], 3, 4))
        [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
    """
    new_column = [1] + column + [ max_value ] * (height - len(column))
    return _next_column_iterator(new_column, height)


class TruncatedStaircases(Parent):
    """
    Factory class for the combinatorial class of truncated staircases
    of size ``n`` with last column ``last_column``.

    EXAMPLES::

        sage: t4 = TruncatedStaircases(4, [2,3]); t4
        Truncated staircases of size 4 with last column [2, 3]
        sage: t4.cardinality()
        4
    """
    __metaclass__ = ClasscallMetaclass
    @staticmethod
    def __classcall_private__(cls, n, last_column, **kwds):
        r"""
        Factory pattern.

        Check properties on arguments, then call the appropriate class.

        TESTS::

            sage: T = TruncatedStaircases(4, [2,3])
            sage: type(T)
            <class 'sage.combinat.alternating_sign_matrix.TruncatedStaircases_nlastcolumn'>

        """
        assert(isinstance(n, (int, Integer)))
        return TruncatedStaircases_nlastcolumn(n, last_column, **kwds)


class TruncatedStaircases_nlastcolumn(TruncatedStaircases):
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
        EXAMPLES::

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

    def __iter__(self):
        """
        EXAMPLES:::

            sage: list(TruncatedStaircases(4, [2,3]))
            [[[4, 3, 2, 1], [3, 2, 1], [3, 2]], [[4, 3, 2, 1], [4, 2, 1], [3, 2]], [[4, 3, 2, 1], [4, 3, 1], [3, 2]], [[4, 3, 2, 1], [4, 3, 2], [3, 2]]]
        """
        for z in self._iterator_rec(self.n):
            yield map(lambda x: list(reversed(x)), z)

    def __eq__(self, other):
        r"""
        TESTS::

            sage: T = TruncatedStaircases(4, [2,3])
            sage: T == loads(dumps(T))
            True
        """
        return ((self.n == other.n) and
                (self.last_column == other.last_column))

    def cardinality(self):
        r"""
        EXAMPLES::

            sage: T = TruncatedStaircases(4, [2,3])
            sage: T.cardinality()
            4
        """
        c = 0
        for _ in self:
            c += 1
        return c
