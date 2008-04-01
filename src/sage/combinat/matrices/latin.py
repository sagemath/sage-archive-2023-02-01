r"""
Latin squares

A {\it latin square} of order $n$ is an $n \times n$ array such that
each symbol $s \in \{ 0, 1, \dots, n-1\}$ appears precisely once in each
row, and precisely once in each column. A {\it partial latin square} of
order $n$ is an $n \times n$ array such that
each symbol $s \in \{ 0, 1, \dots, n-1\}$ appears at most once in each
row, and at most once in each column. A latin square $L$ is a
{\it completion} of a partial latin square $P$ if $P \subseteq L$. If
$P$ completes to just $L$ then $P$ {\it has unique completion}.

A {\it latin bitrade\/} $(T_1,\, T_2)$ is a pair of partial
latin squares such that:
\begin{enumerate}

\item $\{ (i,\,j) \mid (i,\,j,\,k) \in T_1 \mbox{ for some symbol $k$} \}
\newline
= \{ (i,\,j) \mid (i,\,j,\,k') \in T_2 \mbox{ for some symbol $k'$} \}$;

\item for each $(i,\,j,\,k) \in T_1$ and $(i,\,j,\,k') \in T_2$, $k
\neq k'$;

\item the symbols appearing in row $i$ of
$T_1$
are the same as those of row $i$ of
$T_2$;
the symbols appearing in column $j$ of
$T_1$
are the same as those of column $j$ of
$T_2$.

\end{enumerate}

Intuitively speaking, a bitrade gives the difference between two latin
squares, so if $(T_1,\, T_2)$ is a bitrade
for the pair of latin squares $(L_1,\, L_2)$, then
$L1 = (L2 \setminus T_1) \cup T_2$
and
$L2 = (L1 \setminus T_2) \cup T_1$.

This file contains
\begin{enumerate}

\item
LatinSquare class definition;

\item
some named latin squares (back circulant, forward circulant, abelian
$2$-group);

\item
functions is\_partial\_latin\_square and is\_latin\_square to test if a
LatinSquare object satisfies the definition of a latin square or partial
latin square, respectively;

\item tests for completion and unique completion (these use the C++
implementation of Knuth's dancing links algorithm to solve the problem
as a instance of $0-1$ matrix exact cover);

\item functions for calculating the $\tau_i$ representation of a bitrade
and the genus of the associated hypermap embedding;

\item Markov chain of Jacobson and Matthews (1996) for generating latin
squares uniformly at random (provides a generator interface);

\item a few examples of $\tau_i$ representations of bitrades constructed
from the action of a group on itself by right multiplication, functions
for converting to a pair of LatinSquare objects.

\end{enumerate}

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: B = back_circulant(5)
        sage: print B
        [0 1 2 3 4]
        [1 2 3 4 0]
        [2 3 4 0 1]
        [3 4 0 1 2]
        [4 0 1 2 3]
        sage: B.is_latin_square()
        True
        sage: B[0, 1] = 0
        sage: B.is_latin_square()
        False

        sage: (a, b, c, G) = alternating_group_bitrade_generators(1)
        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: print T1
        [ 0 -1  3  1]
        [-1  1  0  2]
        [ 1  3  2 -1]
        [ 2  0 -1  3]
        sage: print T2
        [ 1 -1  0  3]
        [-1  0  2  1]
        [ 2  1  3 -1]
        [ 0  3 -1  2]
        sage: print T1.nr_filled_cells()
        12
        sage: print genus(T1, T2)
        1

 \section{To do}

\begin{enumerate}

\item latin squares with symbols from a ring instead of the integers
$\{ 0, 1, \dots, n-1 \}$.

\item isotopism testing of latin squares and bitrades via graph
isomorphism (nauty?).

\item Combinatorial constructions for bitrades.

\end{enumerate}

AUTHORS:
    - Carlo Hamalainen (2008-03-23): initial version

TESTS:
    sage: L = abelian_2group(3)
    sage: L == loads(dumps(L))
    True
"""

#*****************************************************************************
#       Copyright (C) 2008 Carlo Hamalainen <carlo.hamalainen@gmail.com>,
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

import sets
import copy

from sage.matrix.all import matrix
from sage.rings.all import ZZ
from sage.rings.all import Integer
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.interfaces.gap import GapElement
from sage.combinat.permutation import Permutation
from sage.interfaces.gap import gap
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.rings.arith import is_prime
from sage.rings.finite_field import FiniteField

#load "dancing_links.spyx"
#load "dancing_links.sage"

from dlxcpp import DLXCPP

class LatinSquare:
    def __init__(self, *args):
        """
        Latin squares.

        This class implements a latin square of order n with
        rows and columns indexed by the set {0, 1, ..., n-1}
        and symbols from the same set. The underlying latin square
        is a matrix(ZZ, n, n). If L is a latin square, then the cell
        at row r, column c is empty if and only if L[r, c] < 0. In this
        way we allow partial latin squares and can speak of completions
        to latin squares, etc.

        There are two ways to declare a latin square:

        Empty latin square of order n:

            sage: n = 3
            sage: L = LatinSquare(n)
            sage: L
            [-1 -1 -1]
            [-1 -1 -1]
            [-1 -1 -1]

        Latin square from a matrix:

            sage: M = matrix(ZZ, [[0, 1], [2, 3]])
            sage: print LatinSquare(M)
            [0 1]
            [2 3]
        """

        if len(args) == 1 and (type(args[0]) == Integer or type(args[0]) == int):
            self.square = matrix(ZZ, args[0], args[0])
            self.clear_cells()
        elif len(args) == 2 and (type(args[0]) == Integer or type(args[0]) == int) and (type(args[1]) == Integer or type(args[1]) == int):
            self.square = matrix(ZZ, args[0], args[1])
            self.clear_cells()
        elif len(args) == 1 and type(args[0]) == Matrix_integer_dense:
            self.square = args[0]
        else:
            raise NotImplemented

    def dumps(self):
        """
        Since the latin square class doesn't hold any other private
        variables we just call dumps on self.square:

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(2) == loads(dumps(back_circulant(2)))
            True
        """

        return dumps(self.square)

    def __str__(self):
        """
        The string representation of a latin square is the same as the
        underlying matrix.

        EXAMPLES:
            sage: print LatinSquare(matrix(ZZ, [[0, 1], [2, 3]])).__str__()
            [0 1]
            [2 3]
        """
        return self.square.__str__()

    def __repr__(self):
        """
        The representation of a latin square is the same as the
        underlying matrix.

        EXAMPLES:
            sage: print LatinSquare(matrix(ZZ, [[0, 1], [2, 3]])).__repr__()
            [0 1]
            [2 3]
        """
        return self.square.__str__()
        return self.square.__repr__()

    def __getitem__(self, rc):
        """
        If L is a LatinSquare then this method
        allows us to evaluate L[r, c].

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(3)
            sage: print B[1, 1]
            2
        """

        r = rc[0]
        c = rc[1]

        return self.square[r, c]

    def __setitem__(self, rc, val):
        """
        If L is a LatinSquare then this method
        allows us to set L[r, c].

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(3)
            sage: B[1, 1] = 10
            sage: print B[1, 1]
            10
        """

        r = rc[0]
        c = rc[1]

        self.square[r, c] = val

    def set_immutable(self):
        """
        A latin square is immutable if the underlying matrix is
        immutable.

        EXAMPLES:
            sage: L = LatinSquare(matrix(ZZ, [[0, 1], [2, 3]]))
            sage: L.set_immutable()
            sage: print {L : 0}   # this would fail without set_immutable()
            {[0 1]
            [2 3]: 0}
        """

        self.square.set_immutable()

    def __hash__(self):
        """
        The hash of a latin square is precisely the hash of the
        underlying matrix.

        EXAMPLES:
            sage: L = LatinSquare(matrix(ZZ, [[0, 1], [2, 3]]))
            sage: L.set_immutable()
            sage: print L.__hash__()
            12
        """

        return self.square.__hash__()

    def __eq__(self, Q):
        """
        Two latin squares are equal if the underlying
        matrices are equal.

        EXAMPLES:
            sage: A = LatinSquare(matrix(ZZ, [[0, 1], [2, 3]]))
            sage: B = LatinSquare(matrix(ZZ, [[0, 4], [2, 3]]))
            sage: print A == B
            False
            sage: B[0, 1] = 1
            sage: print A == B
            True
        """

        return self.square == Q.square

    def copy(self):
        """
        To copy a latin square we must copy the underlying matrix.

        EXAMPLES:
            sage: A = LatinSquare(matrix(ZZ, [[0, 1], [2, 3]]))
            sage: B = A.copy()
            sage: print A
            [0 1]
            [2 3]
        """

        C = LatinSquare(self.square.nrows(), self.square.ncols())
        C.square = self.square.copy()
        return C

    def clear_cells(self):
        """
        Mark every cell in self as being empty.

        EXAMPLES:
            sage: A = LatinSquare(matrix(ZZ, [[0, 1], [2, 3]]))
            sage: A.clear_cells()
            sage: print A
            [-1 -1]
            [-1 -1]
        """

        for r in range(self.square.nrows()):
            for c in range(self.square.ncols()):
                self.square[r, c] = -1;

    def nrows(self):
        """
        Number of rows in the latin square.

        EXAMPLES:
            sage: print LatinSquare(3).nrows()
            3
        """

        return self.square.nrows()

    def ncols(self):
        """
        Number of columns in the latin square.

        EXAMPLES:
            sage: print LatinSquare(3).ncols()
            3
        """
        return self.square.ncols()

    def row(self, x):
        """
        Returns row x of the latin square.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(3).row(0)
            (0, 1, 2)
        """

        return self.square.row(x)

    def column(self, x):
        """
        Returns column x of the latin square.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(3).column(0)
            (0, 1, 2)
        """
        return self.square.column(x)

    def list(self):
        """
        Convert the latin square into a list, in a row-wise manner.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(3).list()
            [0, 1, 2, 1, 2, 0, 2, 0, 1]
        """

        return self.square.list()

    def nr_filled_cells(self):
        """
        Returns the number of filled cells (i.e. cells with a positive
        value) in the partial latin square self.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: LatinSquare(matrix([[0, -1], [-1, 0]])).nr_filled_cells()
            2
        """

        s = 0
        for r in range(self.nrows()):
            for c in range(self.ncols()):
                if self[r, c] >= 0: s += 1
        return s

    def apply_isotopism(self, row_perm, col_perm, sym_perm):
        """
        An isotopism is a permutation of the rows, columns, and symbols of a
        partial latin square self. Use isotopism() to convert a tuple (indexed
        from 0) to a Permutation object.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(5)
            sage: print B
            [0 1 2 3 4]
            [1 2 3 4 0]
            [2 3 4 0 1]
            [3 4 0 1 2]
            [4 0 1 2 3]
            sage: alpha = isotopism((0,1,2,3,4))
            sage: beta  = isotopism((1,0,2,3,4))
            sage: gamma = isotopism((2,1,0,3,4))
            sage: print B.apply_isotopism(alpha, beta, gamma)
            [3 4 2 0 1]
            [0 2 3 1 4]
            [1 3 0 4 2]
            [4 0 1 2 3]
            [2 1 4 3 0]
        """

        #Q = matrix(ZZ, self.nrows(), self.ncols())
        Q = LatinSquare(self.nrows(), self.ncols())

        for r in range(self.nrows()):
            for c in range(self.ncols()):
                try:
                    s2 = sym_perm[self[r, c]] - 1
                except IndexError:
                    s2 = self[r, c]  # we must be leaving the symbol fixed?

                Q[row_perm[r]-1, col_perm[c]-1] = s2

        return Q

    def filled_cells_map(self):
        """
        Number the filled cells of self with integers from {1, 2, 3, ...}.

        INPUT:
            self --    Partial latin square self (empty cells have negative
                    values)

        OUTPUT:
            A dictionary cells_map where cells_map[(i,j)] = m means that
            (i,j) is the m-th filled cell in P, while cells_map[m] = (i,j).

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: (a, b, c, G) = alternating_group_bitrade_generators(1)
            sage: (T1, T2) = bitrade_from_group(a, b, c, G)
            sage: print T1.filled_cells_map()
            {1: (0, 0), 2: (0, 2), 3: (0, 3), 4: (1, 1), 5: (1, 2), 6: (1, 3), 7: (2, 0), 8: (2, 1), 9: (2, 2), 10: (3, 0), 11: (3, 1), 12: (3, 3), (2, 1): 8, (1, 3): 6, (0, 3): 3, (1, 2): 5, (3, 3): 12, (3, 0): 10, (2, 2): 9, (1, 1): 4, (0, 0): 1, (3, 1): 11, (2, 0): 7, (0, 2): 2}
        """

        cells_map = {}
        k = 1

        for r in range(self.nrows()):
            for c in range(self.ncols()):
                e = self[r, c]

                if e < 0: continue

                cells_map[ (r,c) ] = k
                cells_map[k] = (r,c)

                k += 1

        return cells_map

    def top_left_empty_cell(self):
        """
        Returns the least [r, c] such that self[r, c] is an empty cell. If
        all cells are filled then we return None.

        INPUT:
            self --    LatinSquare

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(5)
            sage: B[3, 4] = -1
            sage: print B.top_left_empty_cell()
            [3, 4]
        """

        for r in range(self.nrows()):
            for c in range(self.ncols()):
                if self[r, c] < 0:
                    return [r, c]

        return None

    def is_partial_latin_square(self):
        """
        self is a partial latin square if it is an n by n matrix, and each
        symbol in [0, 1, ..., n-1] appears at most once in each row, and at
        most once in each column.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: LatinSquare(4).is_partial_latin_square()
            True
            sage: back_circulant(3).gcs().is_partial_latin_square()
            True
            sage: back_circulant(6).is_partial_latin_square()
            True
        """

        assert self.nrows() == self.ncols()

        n = self.nrows()

        for r in range(n):
            vals_in_row = {}

            for c in range(n):
                e = self[r, c]

                if e < 0: continue

                # Entry out of range 0, 1, ..., n-1:
                if e >= n: return False

                # Entry has already appeared in this row:
                if vals_in_row.has_key(e): return False

                vals_in_row[e] = True

        for c in range(n):
            vals_in_col = {}

            for r in range(n):
                e = self[r, c]

                if e < 0: continue

                # Entry out of range 0, 1, ..., n-1:
                if e >= n: return False

                # Entry has already appeared in this column:
                if vals_in_col.has_key(e): return False

                vals_in_col[e] = True

        return True

    def is_latin_square(self):
        """
        self is a latin square if it is an n by n matrix, and each
        symbol in [0, 1, ..., n-1] appears exactly once in each row, and
        exactly once in each column.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: abelian_2group(4).is_latin_square()
            True

            sage: forward_circulant(7).is_latin_square()
            True
        """

        # We don't allow latin rectangles:
        if self.nrows() != self.ncols():
            return False

        # Every cell must be filled:
        if len(filter(lambda x: x >= 0, self.list())) != self.nrows()*self.ncols():
            return False

        # By necessity self must be a partial latin square:
        if not self.is_partial_latin_square():
            return False

        return True

    def permissable_values(self, r, c):
        """
        Find all values that do not appear in row r and column c of
        the latin square self. If self[r, c] is filled then we return the
        empty list.

        INPUT:
            self -- LatinSquare
            r --    int; row of the latin square
            c --    int; column of the latin square

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: L = back_circulant(5)
            sage: L[0, 0] = -1
            sage: print L.permissable_values(0, 0)
            [0]
        """

        if self[r, c] >= 0:
            return []

        assert self.nrows() == self.ncols()

        n = self.nrows()

        vals = {}
        for e in range(n):
            vals[e] = True

        for i in range(n):
            if self[i, c] >= 0:
                del vals[ self[i, c] ]

        for j in range(n):
            if self[r, j] >= 0:
                try:
                    del vals[ self[r, j] ]
                except KeyError:
                    # We may have already removed a symbol
                    # in the previous for-loop.
                    pass

        return vals.keys()

    def random_empty_cell(self):
        """
        Find an empty cell of self, uniformly at random.

        INPUT:
            self --    LatinSquare

        RETURNS:
            [r, c] -- cell such that self[r, c] is empty, or
            returns None if self is a (full) latin square.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: P = back_circulant(2)
            sage: P[1,1] = -1
            sage: P.random_empty_cell()
            [1, 1]
        """

        cells = {}

        for r in range(self.nrows()):
            for c in range(self.ncols()):
                if self[r, c] < 0:
                    cells[ (r,c) ] = True

        cells = cells.keys()

        if len(cells) == 0: return None

        rc = cells[ ZZ.random_element(len(cells)) ]

        return [rc[0], rc[1]]

    def is_uniquely_completable(self):
        """
        Returns True if the partial latin square self
        has exactly one completion to a latin square. This is
        just a wrapper for the current best-known algorithm,
        Dancing Links by Knuth. See dancing_links.spyx

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(4).gcs().is_uniquely_completable()
            True

            sage: G = abelian_2group(3).gcs()
            sage: G.is_uniquely_completable()
            True

            sage: G[0, 0] = -1
            sage: G.is_uniquely_completable()
            False
        """

        return self.dlxcpp_has_unique_completion()

    def is_completable(self):
        """
        Returns True if the partial latin square
        can be completed to a latin square.

        EXAMPLES:

        The following partial latin square has no completion
        because there is nowhere that we can place the symbol 0
        in the third row:

            sage: B = LatinSquare(3)

            sage: B[0, 0] = 0
            sage: B[1, 1] = 0
            sage: B[2, 2] = 1

            sage: print B
            [ 0 -1 -1]
            [-1  0 -1]
            [-1 -1  1]

            sage: print B.is_completable()
            False

            sage: B[2, 2] = 0
            sage: print B.is_completable()
            True

        """

        return len(dlxcpp_find_completions(self, nr_to_find = 1)) > 0

    def gcs(self):
        """
        A greedy critical set of a latin square self is found by
        successively removing elements in a row-wise (bottom-up)
        manner, checking for unique completion at each step.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: A = abelian_2group(3)
            sage: G = A.gcs()
            sage: print A
            [0 1 2 3 4 5 6 7]
            [1 0 3 2 5 4 7 6]
            [2 3 0 1 6 7 4 5]
            [3 2 1 0 7 6 5 4]
            [4 5 6 7 0 1 2 3]
            [5 4 7 6 1 0 3 2]
            [6 7 4 5 2 3 0 1]
            [7 6 5 4 3 2 1 0]
            sage: print G
            [ 0  1  2  3  4  5  6 -1]
            [ 1  0  3  2  5  4 -1 -1]
            [ 2  3  0  1  6 -1  4 -1]
            [ 3  2  1  0 -1 -1 -1 -1]
            [ 4  5  6 -1  0  1  2 -1]
            [ 5  4 -1 -1  1  0 -1 -1]
            [ 6 -1  4 -1  2 -1  0 -1]
            [-1 -1 -1 -1 -1 -1 -1 -1]
        """

        n = self.nrows()

        G = self.copy()

        for r in range(n-1, -1, -1):
            for c in range(n-1, -1, -1):
                e = G[r, c]
                G[r, c] = -1

                if not G.dlxcpp_has_unique_completion():
                    G[r, c] = e

        return G

    def dlxcpp_has_unique_completion(self):
        """
        Check if the partial latin square self of order n can be
        embedded in precisely one latin square of order n.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: back_circulant(2).dlxcpp_has_unique_completion()
            True
            sage: P = LatinSquare(2)
            sage: P.dlxcpp_has_unique_completion()
            False
            sage: P[0, 0] = 0
            sage: P.dlxcpp_has_unique_completion()
            True
        """
        return len(dlxcpp_find_completions(self, nr_to_find = 2)) == 1

    def vals_in_row(self, r):
        """
        Returns a dictionary with key e if and only if row r of
        self has the symbol e.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(3)
            sage: B[0, 0] = -1
            sage: print back_circulant(3).vals_in_row(0)
            {0: True, 1: True, 2: True}
        """

        n = self.ncols()
        vals_in_row = {}

        for c in range(n):
            e = self[r, c]
            if e >= 0: vals_in_row[e] = True

        return vals_in_row

    def vals_in_col(self, c):
        """
        Returns a dictionary with key e if and only if column c of
        self has the symbol e.

        EXAMPLES:
            sage: from sage.combinat.matrices.latin import *
            sage: B = back_circulant(3)
            sage: B[0, 0] = -1
            sage: print back_circulant(3).vals_in_col(0)
            {0: True, 1: True, 2: True}
        """
        n = self.nrows()
        vals_in_col = {}

        for r in range(n):
            e = self[r, c]
            if e >= 0: vals_in_col[e] = True

        return vals_in_col

def genus(T1, T2):
    """
    Returns the genus of hypermap embedding associated with
    the bitrade (T1, T2). Informally, we compute the
    [tau_1, tau_2, tau_3] permutation representation of the bitrade.
    Each cycle of tau_1, tau_2, and tau_3 gives a rotation scheme for a
    black, white, and star vertex (respectively). The genus then comes
    from Euler's formula. For more details see Carlo Hamalainen: Partitioning
    3-homogeneous latin bitrades. To appear in Geometriae Dedicata,
    available at http://arxiv.org/abs/0710.0938

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: (a, b, c, G) = alternating_group_bitrade_generators(1)
        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: print genus(T1, T2)
        1
        sage: (a, b, c, G) = pq_group_bitrade_generators(3, 7)
        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: print genus(T1, T2)
        3
    """

    cells_map, t1, t2, t3 = tau123(T1, T2)
    return (len(t1.to_cycles()) + len(t2.to_cycles()) + len(t3.to_cycles()) - T1.nr_filled_cells() - 2)/(-2)

def tau123(T1, T2):
    """
    Compute the tau_i representation for a bitrade (T1, T2). See the
    functions tau1, tau2, and tau3 for the mathematical definitions.

    RETURNS:
        (cells_map, t1, t2, t3)

        where cells_map is a map to/from the filled cells of T1, and
        t1, t2, t3 are the tau1, tau2, tau3 permutations.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: (a, b, c, G) = pq_group_bitrade_generators(3, 7)
        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: print T1
        [0 6 4]
        [1 0 5]
        [2 1 6]
        [3 2 0]
        [4 3 1]
        [5 4 2]
        [6 5 3]
        sage: print T2
        [6 4 0]
        [0 5 1]
        [1 6 2]
        [2 0 3]
        [3 1 4]
        [4 2 5]
        [5 3 6]
        sage: (cells_map, t1, t2, t3) = tau123(T1, T2)
        sage: print cells_map
        {1: (0, 0), 2: (0, 1), 3: (0, 2), 4: (1, 0), 5: (1, 1), 6: (1, 2), 7: (2, 0), 8: (2, 1), 9: (2, 2), 10: (3, 0), 11: (3, 1), 12: (3, 2), 13: (4, 0), (2, 1): 8, 15: (4, 2), 16: (5, 0), 17: (5, 1), 18: (5, 2), 19: (6, 0), 20: (6, 1), 21: (6, 2), (5, 1): 17, (4, 0): 13, (1, 2): 6, (3, 0): 10, (5, 0): 16, (2, 2): 9, (4, 1): 14, (1, 1): 5, (3, 2): 12, (0, 0): 1, (6, 0): 19, 14: (4, 1), (4, 2): 15, (1, 0): 4, (0, 1): 2, (6, 1): 20, (3, 1): 11, (2, 0): 7, (6, 2): 21, (5, 2): 18, (0, 2): 3}
        sage: print cells_map_as_square(cells_map, max(T1.nrows(), T1.ncols()))
        [ 1  2  3 -1 -1 -1 -1]
        [ 4  5  6 -1 -1 -1 -1]
        [ 7  8  9 -1 -1 -1 -1]
        [10 11 12 -1 -1 -1 -1]
        [13 14 15 -1 -1 -1 -1]
        [16 17 18 -1 -1 -1 -1]
        [19 20 21 -1 -1 -1 -1]
        sage: t1
        [3, 1, 2, 6, 4, 5, 9, 7, 8, 12, 10, 11, 15, 13, 14, 18, 16, 17, 21, 19, 20]
        sage: t2
        [19, 17, 12, 1, 20, 15, 4, 2, 18, 7, 5, 21, 10, 8, 3, 13, 11, 6, 16, 14, 9]
        sage: print t3
        [5, 9, 13, 8, 12, 16, 11, 15, 19, 14, 18, 1, 17, 21, 4, 20, 3, 7, 2, 6, 10]

        sage: t1.to_cycles()
        [(1, 3, 2),
         (4, 6, 5),
         (7, 9, 8),
         (10, 12, 11),
         (13, 15, 14),
         (16, 18, 17),
         (19, 21, 20)]
        sage: t2.to_cycles()
        [(1, 19, 16, 13, 10, 7, 4),
         (2, 17, 11, 5, 20, 14, 8),
         (3, 12, 21, 9, 18, 6, 15)]
        sage: t3.to_cycles()
        [(1, 5, 12),
         (2, 9, 19),
         (3, 13, 17),
         (4, 8, 15),
         (6, 16, 20),
         (7, 11, 18),
         (10, 14, 21)]

    The product t1*t2*t3 is the identity, i.e. it fixes every point:

        sage: len((t1*t2*t3).fixed_points()) == T1.nr_filled_cells()
        True

    """

    assert is_bitrade(T1, T2)

    cells_map = T1.filled_cells_map()

    t1 = tau1(T1, T2, cells_map)
    t2 = tau2(T1, T2, cells_map)
    t3 = tau3(T1, T2, cells_map)

    return (cells_map, t1, t2, t3)


def isotopism(p):
    """
    Returns a Permutation object that represents an isotopism (for
    rows, columns or symbols of a partial latin square). Since matrices
    in Sage are indexed from 0, this function translates +1 to agree
    with the Permutation class. We also handle PermutationGroupElements:

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: isotopism(5) # identity on 5 points
        [1, 2, 3, 4, 5]

        sage: G = PermutationGroup(['(1,2,3)(4,5)'])
        sage: g = G.gen(0)
        sage: isotopism(g)
        [2, 3, 1, 5, 4]

        sage: isotopism([0,3,2,1]) # 0 goes to 0, 1 goes to 3, etc.
        [1, 4, 3, 2]

        sage: isotopism( (0,1,2) ) # single cycle, presented as a tuple
        [2, 3, 1]

        sage: x = isotopism( ((0,1,2), (3,4)) ) # tuple of cycles
        sage: print x
        [2, 3, 1, 5, 4]
        sage: x.to_cycles()
        [(1, 2, 3), (4, 5)]
    """

    # Identity isotopism on p points:
    if type(p) == Integer:
        return Permutation(range(1, p+1))

    if type(p) == PermutationGroupElement:
        # fixme Ask the Sage mailing list about the tuple/list issue!
        return Permutation(list(p.tuple()))

    if type(p) == list:
        # We expect a list like [0,3,2,1] which means
        # that 0 goes to 0, 1 goes to 3, etc.
        return Permutation(map(lambda x: x+1, p))

    if type(p) == tuple:
        # We have a single cycle:
        if type(p[0]) == Integer:
            return Permutation(tuple(map(lambda x: x+1, p)))

        # We have a tuple of cycles:
        if type(p[0]) == tuple:
            x = isotopism(p[0])

            for i in range(1, len(p)):
                x = x * isotopism(p[i])

            return x

    # Not sure what we got!
    raise NotImplemented

def cells_map_as_square(cells_map, n):
    """
    Returns a LatinSquare with cells numbered from {1, 2, ... } to given
    the dictionary cells_map.

    NOTE:
        The value n should be the maximum of the number of rows and columns of the
        original partial latin square

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: (a, b, c, G) = alternating_group_bitrade_generators(1)
        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: print T1
        [ 0 -1  3  1]
        [-1  1  0  2]
        [ 1  3  2 -1]
        [ 2  0 -1  3]

    There are 12 filled cells in T:

        sage: print cells_map_as_square(T1.filled_cells_map(), max(T1.nrows(), T1.ncols()))
        [ 1 -1  2  3]
        [-1  4  5  6]
        [ 7  8  9 -1]
        [10 11 -1 12]
    """

    assert n > 1

    L = LatinSquare(n, n)

    for r in range(n):
        for c in range(n):
            try:
                L[r, c] = cells_map[ (r,c) ]
            except KeyError:
                # There is no cell (r,c) so skip it
                L[r, c] = -1

    return L

def beta1(rce, T1, T2):
    """
    Find the unique (x, c, e) in T2 such that (r, c, e) is in T1.

    INPUT:
        rce --      tuple (or list) (r, c, e) in T1
        T1, T2 --   latin bitrade

    OUTPUT:
        (x, c, e) in T2.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: print is_bitrade(T1, T2)
        True
        sage: print beta1([0, 0, 0], T1, T2)
        (1, 0, 0)
    """

    r = rce[0]
    c = rce[1]
    e = rce[2]

    assert T1[r, c] == e
    assert e >= 0

    for x in range(T1.nrows()):
        if T2[x, c] == e: return (x, c, e)

    raise PairNotBitrade

def beta2(rce, T1, T2):
    """
    Find the unique (r, x, e) in T2 such that (r, c, e) is in T1.

    INPUT:
        rce --      tuple (or list) (r, c, e) in T1
        T1, T2 --   latin bitrade

    OUTPUT:
        (r, x, e) in T2.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: print is_bitrade(T1, T2)
        True
        sage: print beta2([0, 0, 0], T1, T2)
        (0, 1, 0)
    """

    r = rce[0]
    c = rce[1]
    e = rce[2]

    assert T1[r, c] == e
    assert e >= 0

    for x in range(T1.ncols()):
        if T2[r, x] == e: return (r, x, e)

    raise PairNotBitrade

def beta3(rce, T1, T2):
    """
    Find the unique (r, c, x) in T2 such that (r, c, e) is in T1.

    INPUT:
        rce --      tuple (or list) (r, c, e) in T1
        T1, T2 --   latin bitrade

    OUTPUT:
        (r, c, x) in T2.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: print is_bitrade(T1, T2)
        True
        sage: print beta3([0, 0, 0], T1, T2)
        (0, 0, 4)
    """

    r = rce[0]
    c = rce[1]
    e = rce[2]

    assert T1[r, c] == e
    assert e >= 0

    # fixme this range could be iffy if we
    # work with latin bitrade rectangles...
    for x in range(T1.nrows()):
        if T2[r, c] == x: return (r, c, x)

    raise PairNotBitrade

def tau1(T1, T2, cells_map):
    """
    The definition of tau1 is:

        tau1 : T1 -> T1
        tau1 = beta^(-1)_2 beta_3   (composing left to right)

    where

        beta_i : T2 -> T1

    changes just the i-th coordinate of a triple.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: print is_bitrade(T1, T2)
        True
        sage: (cells_map, t1, t2, t3) = tau123(T1, T2)
        sage: t1 = tau1(T1, T2, cells_map)
        sage: print t1
        [2, 3, 4, 5, 1, 7, 8, 9, 10, 6, 12, 13, 14, 15, 11, 17, 18, 19, 20, 16, 22, 23, 24, 25, 21]
        sage: print t1.to_cycles()
        [(1, 2, 3, 4, 5), (6, 7, 8, 9, 10), (11, 12, 13, 14, 15), (16, 17, 18, 19, 20), (21, 22, 23, 24, 25)]
    """

    # The cells_map has both directions, i.e. integer to
    # cell and cell to integer, so the size of T1 is
    # just half of len(cells_map).
    x = (int(len(cells_map)/2) + 1) * [-1]

    for r in range(T1.nrows()):
        for c in range(T1.ncols()):
            e = T1[r, c]

            if e < 0: continue

            (r2, c2, e2) = beta2( (r,c,e), T1, T2)
            (r3, c3, e3) = beta3( (r2,c2,e2), T2, T1)

            x[ cells_map[(r,c)] ] = cells_map[ (r3,c3) ]

    x.pop(0) # remove the head of the list since we
             # have permutations on 1..(something).

    return Permutation(x)

def tau2(T1, T2, cells_map):
    """
    The definition of tau1 is:

        tau1 : T1 -> T1
        tau1 = beta^(-1)_2 beta_3   (composing left to right)

    where

        beta_i : T2 -> T1

    changes just the i-th coordinate of a triple.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: print is_bitrade(T1, T2)
        True
        sage: (cells_map, t1, t2, t3) = tau123(T1, T2)
        sage: t2 = tau2(T1, T2, cells_map)
        sage: print t2
        [21, 22, 23, 24, 25, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        sage: print t2.to_cycles()
        [(1, 21, 16, 11, 6), (2, 22, 17, 12, 7), (3, 23, 18, 13, 8), (4, 24, 19, 14, 9), (5, 25, 20, 15, 10)]
    """

    # The cells_map has both directions, i.e. integer to
    # cell and cell to integer, so the size of T1 is
    # just half of len(cells_map).
    x = (int(len(cells_map)/2) + 1) * [-1]

    for r in range(T1.nrows()):
        for c in range(T1.ncols()):
            e = T1[r, c]

            if e < 0: continue

            (r2, c2, e2) = beta3( (r,c,e), T1, T2)
            (r3, c3, e3) = beta1( (r2,c2,e2), T2, T1)

            x[ cells_map[(r,c)] ] = cells_map[ (r3,c3) ]

    x.pop(0) # remove the head of the list since we
             # have permutations on 1..(something).

    return Permutation(x)

def tau3(T1, T2, cells_map):
    """
    The definition of tau1 is:

        tau1 : T1 -> T1
        tau1 = beta^(-1)_2 beta_3   (composing left to right)

    where

        beta_i : T2 -> T1

    changes just the i-th coordinate of a triple.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: print is_bitrade(T1, T2)
        True
        sage: (cells_map, t1, t2, t3) = tau123(T1, T2)
        sage: t3 = tau3(T1, T2, cells_map)
        sage: print t3
        [10, 6, 7, 8, 9, 15, 11, 12, 13, 14, 20, 16, 17, 18, 19, 25, 21, 22, 23, 24, 5, 1, 2, 3, 4]
        sage: print t3.to_cycles()
        [(1, 10, 14, 18, 22), (2, 6, 15, 19, 23), (3, 7, 11, 20, 24), (4, 8, 12, 16, 25), (5, 9, 13, 17, 21)]
    """

    # The cells_map has both directions, i.e. integer to
    # cell and cell to integer, so the size of T1 is
    # just half of len(cells_map).
    x = (int(len(cells_map)/2) + 1) * [-1]

    for r in range(T1.nrows()):
        for c in range(T1.ncols()):
            e = T1[r, c]

            if e < 0: continue

            (r2, c2, e2) = beta1( (r,c,e), T1, T2)
            (r3, c3, e3) = beta2( (r2,c2,e2), T2, T1)

            x[ cells_map[(r,c)] ] = cells_map[ (r3,c3) ]

    x.pop(0) # remove the head of the list since we
             # have permutations on 1..(something).

    return Permutation(x)

def back_circulant(n):
    """
    The back-circulant latin square of order n is the
    Cayley table for (Z_n, +), the integers under addition
    modulo n.

    INPUT:
        n --    int; order of the latin square.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: print back_circulant(5)
        [0 1 2 3 4]
        [1 2 3 4 0]
        [2 3 4 0 1]
        [3 4 0 1 2]
        [4 0 1 2 3]
    """

    assert n >= 1

    L = LatinSquare(n, n)

    for r in range(n):
        for c in range(n):
            L[r, c] = (r + c) % n

    return L

def forward_circulant(n):
    """
    The forward-circulant latin square of order n is the
    Cayley table for the operation r + c = (n-c+r) mod n.

    INPUT:
        n --    int; order of the latin square.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: print forward_circulant(5)
        [0 4 3 2 1]
        [1 0 4 3 2]
        [2 1 0 4 3]
        [3 2 1 0 4]
        [4 3 2 1 0]
    """

    assert n >= 1

    L = LatinSquare(n, n)

    for r in range(n):
        for c in range(n):
            L[r, c] = (n-c+r) % n

    return L

def direct_product(L1, L2, L3, L4):
    """
    The 'direct product' of four latin squares L1, L2, L3, L4 of order n
    is the latin square of order 2n consisting of:

    -----------
    | L1 | L2 |
    -----------
    | L3 | L4 |
    -----------

    where the subsquares L2 and L3 have entries offset by n.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: direct_product(back_circulant(4), back_circulant(4), abelian_2group(2), abelian_2group(2))
        [0 1 2 3 4 5 6 7]
        [1 2 3 0 5 6 7 4]
        [2 3 0 1 6 7 4 5]
        [3 0 1 2 7 4 5 6]
        [4 5 6 7 0 1 2 3]
        [5 4 7 6 1 0 3 2]
        [6 7 4 5 2 3 0 1]
        [7 6 5 4 3 2 1 0]
    """

    assert L1.nrows() == L2.nrows() == L3.nrows() == L4.nrows()
    assert L1.ncols() == L2.ncols() == L3.ncols() == L4.ncols()
    assert L1.nrows() == L1.ncols()

    n = L1.nrows()

    D = LatinSquare(2*n, 2*n)

    for r in range(n):
        for c in range(n):
            D[r, c] = L1[r, c]
            D[r, c+n] = L2[r, c] + n
            D[r+n, c] = L3[r, c] + n
            D[r+n, c+n] = L4[r, c]

    return D

def abelian_2group(s):
    """
    Returns the latin square based on the Cayley table for the
    abelian 2-group of order 2^s.

    INPUT:
        s --    int; order of the latin square will be 2^s.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: print abelian_2group(3)
        [0 1 2 3 4 5 6 7]
        [1 0 3 2 5 4 7 6]
        [2 3 0 1 6 7 4 5]
        [3 2 1 0 7 6 5 4]
        [4 5 6 7 0 1 2 3]
        [5 4 7 6 1 0 3 2]
        [6 7 4 5 2 3 0 1]
        [7 6 5 4 3 2 1 0]
    """

    assert s > 0

    if s == 1:
        L = LatinSquare(2, 2)
        L[0, 0] = 0
        L[0, 1] = 1
        L[1, 0] = 1
        L[1, 1] = 0

        return L
    else:
        L_prev = abelian_2group(s-1)
        L = LatinSquare(2**s, 2**s)

        offset = L.nrows()/2

        for r in range(L_prev.nrows()):
            for c in range(L_prev.ncols()):
                L[r, c] = L_prev[r, c]
                L[r+offset, c] = L_prev[r, c] + offset
                L[r, c+offset] = L_prev[r, c] + offset
                L[r+offset, c+offset] = L_prev[r, c]
    return L

def coin():
    """
    Simulates a fair coin (returns True or False)
    using ZZ.random_element(2).

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: x = coin()
        sage: x == 0 or x == 1
        True
    """

    return ZZ.random_element(2) == 0


def next_conjugate(L):
    """
    Permute L[r, c] = e to the conjugate L[c, e] = r

    We assume that L is an n by n matrix and has values
    in the range 0, 1, ..., n-1.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: L = back_circulant(6)
        sage: print L
        [0 1 2 3 4 5]
        [1 2 3 4 5 0]
        [2 3 4 5 0 1]
        [3 4 5 0 1 2]
        [4 5 0 1 2 3]
        [5 0 1 2 3 4]
        sage: print next_conjugate(L)
        [0 1 2 3 4 5]
        [5 0 1 2 3 4]
        [4 5 0 1 2 3]
        [3 4 5 0 1 2]
        [2 3 4 5 0 1]
        [1 2 3 4 5 0]
        sage: L == next_conjugate(next_conjugate(next_conjugate(L)))
        True
    """

    assert L.nrows() == L.ncols()

    n = L.nrows()

    C = LatinSquare(n, n)

    for r in range(n):
        for c in range(n):
            e = L[r, c]
            assert e >= 0 and e < n

            C[c, e] = r

    return C

def row_containing_sym(L, c, x):
    """
    Given an improper latin square L with L[r1, c] = L[r2, c] = x, return
    r1 or r2 with equal probability. This is an internal function and
    should only be used in LatinSquare_generator().

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: L = matrix([(0, 1, 0, 3), (3, 0, 2, 1), (1, 0, 3, 2), (2, 3, 1, 0)])
        sage: print L
        [0 1 0 3]
        [3 0 2 1]
        [1 0 3 2]
        [2 3 1 0]
        sage: c = row_containing_sym(L, 1, 0)
        sage: print c == 1 or c == 2
        True
    """

    r1 = -1
    r2 = -1

    for r in range(L.nrows()):
        if r1 >= 0 and r2 >= 0: break

        if L[r, c] == x and r1 < 0:
            r1 = r
            continue

        if L[r, c] == x and r2 < 0:
            r2 = r
            break

    assert r1 >= 0 and r2 >= 0

    if coin():  return r1
    else:       return r2

def column_containing_sym(L, r, x):
    """
    Given an improper latin square L with L[r, c1] = L[r, c2] = x, return
    c1 or c2 with equal probability. This is an internal function and
    should only be used in LatinSquare_generator().

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: L = matrix([(1, 0, 2, 3), (0, 2, 3, 0), (2, 3, 0, 1), (3, 0, 1, 2)])
        sage: print L
        [1 0 2 3]
        [0 2 3 0]
        [2 3 0 1]
        [3 0 1 2]
        sage: c = column_containing_sym(L, 1, 0)
        sage: print c == 0 or c == 3
        True
    """

    c1 = -1
    c2 = -1

    for c in range(L.ncols()):
        if c1 >= 0 and c2 >= 0: break

        if L[r, c] == x and c1 < 0:
            c1 = c
            continue

        if L[r, c] == x and c2 < 0:
            c2 = c
            break

    assert c1 >= 0 and c2 >= 0

    if coin():  return c1
    else:       return c2

def LatinSquare_generator(L_start, check_assertions = False):
    """
    Generator for a sequence of uniformly distributed latin squares,
    given L_start as the initial latin square. This code implements the
    Markov chain algorithm of Jacobson and Matthews (1996), see below
    for the BibTex entry. This generator will never throw the
    StopIteration exception, so it provides an infinite sequence of
    latin squares.

    EXAMPLES:

    Use the back circulant latin square of order 4 as the initial square
    and print the next two latin squares given by the Markov chain:

        sage: from sage.combinat.matrices.latin import *
        sage: g = LatinSquare_generator(back_circulant(4))
        sage: print g.next().is_latin_square()
        True

    REFERENCE:

    @article{MR1410617,
        AUTHOR = {Jacobson, Mark T. and Matthews, Peter},
         TITLE = {Generating uniformly distributed random {L}atin squares},
       JOURNAL = {J. Combin. Des.},
      FJOURNAL = {Journal of Combinatorial Designs},
        VOLUME = {4},
          YEAR = {1996},
        NUMBER = {6},
         PAGES = {405--437},
          ISSN = {1063-8539},
       MRCLASS = {05B15 (60J10)},
      MRNUMBER = {MR1410617 (98b:05021)},
    MRREVIEWER = {Lars D{\o}vling Andersen},
    }
    """


    if check_assertions: assert L_start.is_latin_square()

    n = L_start.nrows()

    r1 = r2 = c1 = c2 = x = y = z = -1
    proper = True

    L = L_start.copy()

    L_rce = L
    L_cer = LatinSquare(n, n)
    L_erc = LatinSquare(n, n)

    while True:
        if proper:
            if check_assertions: assert L.is_latin_square()

            #################################
            # Update the other two conjugates
            for r in range(n):
                for c in range(n):
                    e = L[r, c]

                    L_cer[c, e] = r
                    L_erc[e, r] = c
            #################################

            yield L

            r1 = ZZ.random_element(n)
            c1 = ZZ.random_element(n)
            x = L[r1, c1]

            y = x
            while y == x:
                y = ZZ.random_element(n)

            # Now find y in row r1 and column c1.
            if check_assertions:
                r2 = 0
                c2 = 0
                while L[r1, c2] != y: c2 += 1
                while L[r2, c1] != y: r2 += 1

                assert L_erc[y, r1] == c2
                assert L_cer[c1, y] == r2

            c2 = L_erc[y, r1]
            r2 = L_cer[c1, y]

            if check_assertions: assert L[r1, c2] == y
            if check_assertions: assert L[r2, c1] == y

            L[r1, c1] = y
            L[r1, c2] = x
            L[r2, c1] = x

            # Now deal with the unknown point.
            # We want to form z + (y - x)
            z = L[r2, c2]

            if z == x:
                L[r2, c2] = y
            else:
                # z and y have positive coefficients
                # x is the improper term with a negative coefficient
                proper = False
        else: # improper square,
            # L[r2, c2] = y + z - x
            # y and z are proper while x is the
            # improper symbol in the cell L[r2, c2].

            r1 = row_containing_sym(L, c2, x)
            c1 = column_containing_sym(L, r2, x)

            if check_assertions: assert L[r1, c2] == x
            if check_assertions: assert L[r2, c1] == x

            # choose one of the proper symbols
            # uniformly at random (we will use whatever
            # lands in variable y).
            if coin():
                y, z = z, y

            # Add/subtract the symbolic difference (y - x)
            L[r2, c2] = z
            L[r1, c2] = y
            L[r2, c1] = y

            if L[r1, c1] == y:
                L[r1, c1] = x
                proper = True
            else: # got another improper square
                z = L[r1, c1]
                x, y = y, x
                r2 = r1
                c2 = c1

                # Now we have L[r2, c2] = z+y-x as
                # usual
                proper = False # for emphasis

def group_to_LatinSquare(G):
    """
    Construct a latin square on the symbols [0, 1, ..., n-1] for
    a group with an n by n Cayley table.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: print group_to_LatinSquare(DihedralGroup(2))
        [0 1 2 3]
        [1 0 3 2]
        [2 3 0 1]
        [3 2 1 0]

        sage: G = gap.Group(PermutationGroupElement((1,2,3)))
        sage: print group_to_LatinSquare(G)
        [0 1 2]
        [1 2 0]
        [2 0 1]
    """

    if isinstance(G, GapElement):
        rows = map(lambda x: list(x), list(gap.MultiplicationTable(G)))
        new_rows = []

        for x in rows:
            new_rows.append(map(lambda x: int(x)-1, x))

        return matrix(new_rows)

    # Otherwise we must have some kind of Sage permutation group object,
    # such as sage.groups.perm_gps.permgroup.PermutationGroup_generic
    # or maybe sage.groups.perm_gps.permgroup_named.
    # We should have a cayley_table() function.

    try:
        C = G.cayley_table()
    except:
        # We got some other kind of group object?
        raise NotImplemented

    L = LatinSquare(C.nrows(), C.ncols())

    k = 0
    entries = {}

    for r in range(C.nrows()):
        for c in range(C.ncols()):
            e_table = C[r, c]

            if entries.has_key(e_table):
                L[r, c] = entries[e_table]
            else:
                entries[e_table] = k
                L[r, c] = k
                k += 1

    return L

def alternating_group_bitrade_generators(m):
    """
    Construct generators a, b, c for the alternating group on 3m+1 points,
    such that a*b*c = 1.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: a, b, c, G = alternating_group_bitrade_generators(1)

        ((1,2,3),
         (1,4,2),
         (2,4,3),
         Permutation Group with generators [(1,2,3), (1,4,2)])

        sage: a*b*c
        ()

        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: print T1
        [ 0 -1  3  1]
        [-1  1  0  2]
        [ 1  3  2 -1]
        [ 2  0 -1  3]
        sage: print T2
        [ 1 -1  0  3]
        [-1  0  2  1]
        [ 2  1  3 -1]
        [ 0  3 -1  2]
    """

    assert m >= 1

    a = tuple(range(1, 2*m+1 + 1))

    b = tuple(range(m+1, 0, -1) + range(2*m+2, 3*m+1 + 1))

    a = PermutationGroupElement(a)
    b = PermutationGroupElement(b)
    c = PermutationGroupElement((a*b)**(-1))

    G = PermutationGroup([a, b])

    return (a, b, c, G)

def pq_group_bitrade_generators(p, q):
    """
    Generators for a group of order pq where p and q are primes
    such that (q % p) == 1.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: pq_group_bitrade_generators(3,7)
        ((2,3,5)(4,7,6), (1,2,3,4,5,6,7), (1,4,2)(3,5,6), Permutation Group with generators [(2,3,5)(4,7,6), (1,2,3,4,5,6,7)])
    """

    assert is_prime(p)
    assert is_prime(q)
    assert (q % p) == 1

    # beta is a primitive root of the
    # congruence x^p = 1 mod q
    F = FiniteField(q)
    fgen = F.multiplicative_generator()
    beta = fgen**((q-1)/p)

    assert beta != 1
    assert (beta**p % q) == 1

    Q = tuple(range(1, q+1))

    P = []
    seenValues = {}
    for i in range(2, q):
        if seenValues.has_key(i):
            continue

        cycle = []
        for k in range(p):
            x = (1 + (i-1)*beta**k) % q
            if x == 0: x = q

            seenValues[x] = True
            cycle.append(x)
        P.append(tuple(cycle))

    G = PermutationGroup([P, Q])
    assert G.order() == p*q
    assert not G.is_abelian()

    a = PermutationGroupElement(P)
    b = PermutationGroupElement(Q)
    c = PermutationGroupElement((a*b)**(-1))

    return (a, b, c, PermutationGroup([P, Q]))

def p3_group_bitrade_generators(p):
    """
    Generators for a group of order p^3 where p is a prime.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: p3_group_bitrade_generators(3)
        ((1,2,3)(4,16,17)(5,20,21)(6,10,14)(7,11,15)(8,24,18)(9,26,23)(12,19,25)(13,22,27),
         (1,4,5)(2,8,9)(3,12,13)(6,18,22)(7,19,23)(10,25,20)(11,16,27)(14,17,26)(15,24,21),
         (1,21,8)(2,23,12)(3,27,4)(5,17,10)(6,13,25)(7,26,16)(9,18,14)(11,22,24)(15,20,19),
         Permutation Group with generators [(1,2,3)(4,16,17)(5,20,21)(6,10,14)(7,11,15)(8,24,18)(9,26,23)(12,19,25)(13,22,27), (1,4,5)(2,8,9)(3,12,13)(6,18,22)(7,19,23)(10,25,20)(11,16,27)(14,17,26)(15,24,21), (1,6,7)(2,10,11)(3,14,15)(4,18,19)(5,22,23)(8,25,16)(9,20,27)(12,17,24)(13,26,21)])
    """

    assert is_prime(p)

    F = gap.new("FreeGroup(3)")

    a = F.gen(1)
    b = F.gen(2)
    c = F.gen(3)

    rels = []
    rels.append( a**p )
    rels.append( b**p )
    rels.append( c**p )
    rels.append( a*b*((b*a*c)**(-1)) )
    rels.append( c*a*((a*c)**(-1)) )
    rels.append( c*b*((b*c)**(-1)) )

    G = F.FactorGroupFpGroupByRels(rels)

    N = gap.NormalClosure(F, gap.Subgroup(F, rels))
    niso = gap.NaturalHomomorphismByNormalSubgroupNC(F, N)

    a = PermutationGroupElement(gap.Image(niso, a))
    b = PermutationGroupElement(gap.Image(niso, b))
    c = PermutationGroupElement(gap.Image(niso, c))

    G = PermutationGroup([a, b, c])

    c = PermutationGroupElement((a*b)**(-1))

    return (a, b, c, G)

def check_bitrade_generators(a, b, c):
    """
    Three group elements a, b, c will generate a bitrade if
    a*b*c = 1 and the subgroups <a>, <b>, <c> intersect (pairwise)
    in just the identity.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: a, b, c, G = p3_group_bitrade_generators(3)
        sage: print check_bitrade_generators(a, b, c)
        True
        sage: print check_bitrade_generators(a, b, G.identity())
        False
    """

    A = PermutationGroup([a])
    B = PermutationGroup([b])
    C = PermutationGroup([c])

    if a*b != c**(-1): return False

    X = gap.Intersection(gap.Intersection(A, B), C)
    return X.Size() == 1

def is_bitrade(T1, T2):
    """
    Combinatorially, a pair (T1, T2) of partial latin squares
    is a bitrade if they are disjoint, have the same shape, and
    have row and column balance. For definitions of each of these
    terms see the relevant function in this file.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: print is_bitrade(T1, T2)
        True
    """

    if not is_disjoint(T1, T2): return False
    if not is_same_shape(T1, T2): return False
    if not is_row_and_col_balanced(T1, T2): return False

    return True

def is_primary_bitrade(a, b, c, G):
    """
    A bitrade generated from elements a, b, c is
    primary if <a, b, c> = G.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: (a, b, c, G) = p3_group_bitrade_generators(5)
        sage: print is_primary_bitrade(a, b, c, G)
        True
    """

    H = PermutationGroup([a, b, c])

    return G == H

def coset_representatives(x, G):
    """
    For some group element x in G, return a list of the canonical
    coset representatives for the subgroup <x> in G.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: a, b, c, G = alternating_group_bitrade_generators(1)
        sage: print coset_representatives(a, G)
        [(), (1,2)(3,4), (1,4)(2,3), (1,3)(2,4)]
    """

    rcosets = gap.RightCosets(G, PermutationGroup([x]))

    cosetReps = []

    for i in range(1, len(rcosets)+1):
        cosetReps.append(gap.Representative(rcosets[i]))

    assert gap.Index(G, PermutationGroup([x])) == len(cosetReps)

    return cosetReps

def entry_label(eLabels, C, x):
    """
    ARGUMENTS:
        eLabels --- list of coset representatives of C in G
        C       --- subgroup G of G
        x       --- element of G

    RETURNS:
        integer i such that C*eLabels[i] == C*x

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: a, b, c, G = alternating_group_bitrade_generators(1)
        sage: print entry_label(coset_representatives(c, G), gap.Group([c]), a)
        1
    """

    # fixme There should be a nicer way to do this using
    # dictionary.

    Cx = gap.RightCoset(C, x)

    for i in range(len(eLabels)):
        y = eLabels[i]

        Cy = gap.RightCoset(C, y)

        if Cx == Cy: return i

    assert False

def bitrade_from_group(a, b, c, G):
    """
    Given group elements a, b, c in G such that abc = 1 and
    the subgroups <a>, <b>, <c> intersect (pairwise) only in the identity,
    construct a bitrade (T1, T2) where rows, columns, and symbols correspond
    to cosets of <a>, <b>, and <c>, respectively.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: a, b, c, G = alternating_group_bitrade_generators(1)
        sage: (T1, T2) = bitrade_from_group(a, b, c, G)
        sage: print T1
        [ 0 -1  3  1]
        [-1  1  0  2]
        [ 1  3  2 -1]
        [ 2  0 -1  3]
        sage: print T2
        [ 1 -1  0  3]
        [-1  0  2  1]
        [ 2  1  3 -1]
        [ 0  3 -1  2]
    """

    T1 = {}
    T2 = {}

    # The rows of the bitrade are labeled using
    # the cosets of <a> in G.
    rLabels = coset_representatives(a, G)
    cLabels = coset_representatives(b, G)
    eLabels = coset_representatives(c, G)

    T1 = LatinSquare(len(rLabels), len(cLabels))
    T2 = LatinSquare(len(rLabels), len(cLabels))

    for i in range(len(rLabels)):
        for j in range(len(cLabels)):
            T1[i, j] = -1
            T2[i, j] = -1

    A = gap.Group([a])
    B = gap.Group([b])
    C = gap.Group([c])

    for i in range(len(rLabels)):
        for j in range(len(cLabels)):
            for k in range(len(eLabels)):
                cosetRow = gap.RightCoset(A, rLabels[i])
                cosetCol = gap.RightCoset(B, cLabels[j])
                cosetEnt = gap.RightCoset(C, eLabels[k])

                g = gap.Intersection(gap.Elements(cosetRow), gap.Elements(cosetCol), gap.Elements(cosetEnt))

                if len(g) == 0: continue

                # At this point we know that a, b, c are the right
                # generators for a group-based bitrade so we can get
                # the *unique* element g in the intersection of
                # cosetRow, cosetCol, and cosetEnt by getting the
                # intersection of just one pair.
                g = gap.Intersection(cosetRow, cosetCol)

                assert len(g) == 1

                g = g[1]
                #eIdx = entry_label(eLabels, C, a**(-1) * g)
                eIdx = entry_label(eLabels, C, a * g)
                # fixme we're doing things with right instead of
                # left cosets now, oops! The published work uses left
                # cosets.

                T1[i, j] = k
                T2[i, j] = eIdx

    return (T1, T2)

def is_disjoint(T1, T2):
    """
    The partial latin squares T1 and T2 are disjoint
    if T1[r, c] != T2[r, c] or T1[r, c] == T2[r, c] == -1 for each
    cell [r, c].

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import is_disjoint, back_circulant, isotopism
        sage: is_disjoint(back_circulant(2), back_circulant(2))
        False

        sage: T1 = back_circulant(5)
        sage: x = isotopism( (0,1,2,3,4) )
        sage: y = isotopism(5) # identity
        sage: z = isotopism(5) # identity
        sage: T2 = T1.apply_isotopism(x, y, z)
        sage: print is_disjoint(T1, T2)
        True
    """

    for i in range(T1.nrows()):
        for j in range(T1.ncols()):
            if T1[i, j] < 0 and T2[i, j] < 0: continue

            if T1[i, j] == T2[i, j]:
                return False

    return True

def is_same_shape(T1, T2):
    """
    Two partial latin squares T1, T2 have the same shape
    if T1[r, c] >= 0 if and only if T2[r, c] >= 0.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: is_same_shape(abelian_2group(2), back_circulant(4))
        True
        sage: is_same_shape(LatinSquare(5), LatinSquare(5))
        True
        sage: is_same_shape(forward_circulant(5), LatinSquare(5))
        False
    """

    for i in range(T1.nrows()):
        for j in range(T1.ncols()):
            if T1[i, j] < 0 and T2[i, j] < 0: continue
            if T1[i, j] >= 0 and T2[i, j] >= 0: continue

            return False

    return True

def is_row_and_col_balanced(T1, T2):
    """
    Partial latin squares T1 and T2 are balanced if the
    symbols appearing in row r of T1 are the same as the symbols
    appearing in row r of T2, for each r, and if the same
    condition holds on columns.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: T1 = matrix([[0,1,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1]])
        sage: T2 = matrix([[0,1,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1]])
        sage: is_row_and_col_balanced(T1, T2)
        True
        sage: T2 = matrix([[0,3,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1], [-1,-1,-1,-1]])
        sage: is_row_and_col_balanced(T1, T2)
        False
    """

    for r in range(T1.nrows()):
        val1 = sets.Set(filter(lambda x: x >= 0, T1.row(r)))
        val2 = sets.Set(filter(lambda x: x >= 0, T2.row(r)))

        if val1 != val2: return False

    for c in range(T1.ncols()):
        val1 = sets.Set(filter(lambda x: x >= 0, T1.column(c)))
        val2 = sets.Set(filter(lambda x: x >= 0, T2.column(c)))

        if val1 != val2: return False

    return True

def dlxcpp_rows_and_map(P):
    """
    Internal function for dlxcpp_find_completions. Given a partial
    latin square P we construct a list of rows of a 0-1
    matrix M such that an exact cover of M corresponds to a
    completion of P to a latin square.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: dlxcpp_rows_and_map(LatinSquare(2))
        ([[0, 4, 8], [1, 5, 8], [2, 4, 9], [3, 5, 9], [0, 6, 10], [1, 7, 10], [2, 6, 11], [3, 7, 11]], {(2, 4, 9): (0, 1, 0), (1, 5, 8): (0, 0, 1), (1, 7, 10): (1, 0, 1), (0, 6, 10): (1, 0, 0), (3, 7, 11): (1, 1, 1), (2, 6, 11): (1, 1, 0), (0, 4, 8): (0, 0, 0), (3, 5, 9): (0, 1, 1)})
    """

    assert P.nrows() == P.ncols()

    n = P.nrows()

    # We will need 3n^2 columns in total:
    #
    # n^2 for the xCy columns
    # n^2 for the xRy columns
    # n^2 for the xy columns

    dlx_rows = []
    cmap = {}

    for r in range(n):
        valsrow = P.vals_in_row(r)

        for c in range(n):
            valscol = P.vals_in_col(c)

            for e in range(n):
                # These should be constants
                c_OFFSET  = e + c*n
                r_OFFSET  = e + r*n + n*n
                xy_OFFSET = 2*n*n + r*n + c

                cmap[(c_OFFSET, r_OFFSET, xy_OFFSET)] = (r,c,e)

                #print "possibility: ", r, c, e, "offsets:", c_OFFSET, r_OFFSET, xy_OFFSET

                #if P[r, c] >= 0: continue

                # We only want the correct value to pop in here
                if P[r, c] >= 0 and P[r, c] != e: continue

                if P[r, c] < 0 and valsrow.has_key(e): continue
                if P[r, c] < 0 and valscol.has_key(e): continue

                dlx_rows.append([c_OFFSET, r_OFFSET, xy_OFFSET])

    return dlx_rows, cmap


def dlxcpp_find_completions(P, nr_to_find = None):
    """
    Returns a list of all latin squares L of the same order as P such that
    P is contained in L. The optional parameter nr_to_find limits the
    number of latin squares that are found.

    EXAMPLES:
        sage: from sage.combinat.matrices.latin import *
        sage: dlxcpp_find_completions(LatinSquare(2))
        [[0 1]
        [1 0], [1 0]
        [0 1]]

        sage: dlxcpp_find_completions(LatinSquare(2), 1)
        [[0 1]
        [1 0]]
    """


    assert P.nrows() == P.ncols()

    n = P.nrows()

    dlx_rows, cmap = dlxcpp_rows_and_map(P)

    SOLUTIONS = {}
    for x in DLXCPP(dlx_rows):
        x.sort()
        SOLUTIONS[tuple(x)] = True

        if nr_to_find is not None and len(SOLUTIONS) >= nr_to_find: break

    comps = []

    for i in SOLUTIONS.keys():
        soln = list(i)

        Q = copy.deepcopy(P)

        for x in soln:
            (r, c, e) = cmap[tuple(dlx_rows[x])]

            if Q[r, c] >= 0:
                assert Q[r, c] == e
            else:
                Q[r, c] = e

        comps.append(Q)

    return comps

