"""
Tableaux

AUTHORS:

    - Mike Hansen (2007): initial version

    - Jason Bandlow (2011): updated to use Parent/Element model, and many
      minor fixes

This file consists of the following major classes:

    Element classes:

        * Tableau
        * SemistandardTableau
        * StandardTableau

    Factory classes:

        * Tableaux
        * SemistandardTableaux
        * StandardTableaux

    Parent classes:

        * Tableaux_all
        * Tableaux_size
        * SemistandardTableaux_all (facade class)
        * SemistandardTableaux_size
        * SemistandardTableaux_size_inf
        * SemistandardTableaux_size_weight
        * SemistandardTableaux_shape
        * SemistandardTableaux_shape_inf
        * SemistandardTableaux_shape_weight
        * StandardTableaux_all (facade class)
        * StandardTableaux_size
        * StandardTableaux_shape

TODO:

    - Move methods that only apply to semistandard tableaux from tableau to
      semistandard tableau

    - Add a class for tableaux of a given shape (eg Tableaux_shape)
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2011 Jason Bandlow <jbandlow@gmail.com>
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
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.infinity import PlusInfinity
from sage.rings.arith import factorial
from sage.rings.integer import Integer
from sage.combinat.combinat import CombinatorialObject
import sage.combinat.skew_tableau
from integer_vector import IntegerVectors
import sage.libs.symmetrica.all as symmetrica
import sage.misc.prandom as random
import copy
import permutation
from sage.misc.flatten import flatten
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.misc.misc import uniq
from sage.misc.sage_unittest import TestSuite
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_cat import Sets
import __builtin__

class Tableau(CombinatorialObject, Element):
    """
    A class to model a tableau.

    INPUT:

    - ``t`` -- a Tableau, a list of lists, or an empty list

    OUTPUT:

    - A Tableau object constructed from ``t``.

    A tableau in Sage is a finite list of lists, whose lengths are weakly
    decreasing, or an empty list, representing the empty tableau.  The entries
    of a tableau can be any sage object.

    Note that Sage uses the English convention for partitions and
    tableaux; the longer rows are displayed on top.

    EXAMPLES::

        sage: t = Tableau([[1,2,3],[4,5]]); t
        [[1, 2, 3], [4, 5]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty print
        1 2 3
        4 5
        sage: t.is_standard()
        True

        sage: Tableau([['a','c','b'],[[],(2,1)]])
        [['a', 'c', 'b'], [[], (2, 1)]]

        sage: Tableau([]) # The empty tableau
        []

    When using code that will generate a lot of tableaux, it is slightly more
    efficient to construct a Tableau from the appropriate Parent object::

        sage: T = Tableaux()
        sage: T([[1, 2, 3], [4, 5]])
        [[1, 2, 3], [4, 5]]

    SEE ALSO:

        - :class:`Tableaux`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`

    TESTS::

        sage: Tableau([[1],[2,3]])
        Traceback (most recent call last):
        ...
        ValueError: A tableau must be a list of lists of weakly decreasing length.
        sage: Tableau([1,2,3])
        Traceback (most recent call last):
        ...
        ValueError: A tableau must be a list of lists.

    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a Tableau is only ever constructed as an
        element_class call of an appropriate parent.

        TESTS::

            sage: t = Tableau([[1,1],[1]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Tableaux
            sage: t.category()
            Category of elements of Tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.Tableaux_all_with_category.element_class'>
        """
        if isinstance(t, Tableau):
            return t

        return Tableaux_all().element_class(Tableaux_all(), t)

    def __init__(self, parent, t):
        r"""
        Initializes a tableau.

        TESTS::

            sage: t = Tableaux()([[1,1],[1]])
            sage: s = Tableaux(3)([[1,1],[1]])
            sage: s==t
            True
            sage: t.parent()
            Tableaux
            sage: s.parent()
            Tableaux of size 3
            sage: r = Tableaux()(s); r.parent()
            Tableaux
            sage: s is t # identical tableaux are distinct objects
            False
        """
        Element.__init__(self, parent)
        CombinatorialObject.__init__(self, t)

        if isinstance(t, Tableau):
            return None

        # CombinatorialObject verifies that t is a list
        # We must verify t is a list of lists
        if not all(isinstance(row, __builtin__.list) for row in t):
            raise ValueError, "A tableau must be a list of lists."

        if not map(len,t) in sage.combinat.partition.Partitions_all():
            raise ValueError, "A tableau must be a list of lists of weakly decreasing length."

    def __setstate__(self, state):
        """
        In order to maintain backwards compatibility and be able to unpickle a
        old pickle from Tableau_class we have to override the default
        __setstate__.

        TESTS::

            sage: loads('x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1+IL\xcaIM,\xe5\n\x81\xd0\xf1\xc99\x89\xc5\xc5\\\x85\x8c\x9a\x8d\x85L\xb5\x85\xcc\x1a\xa1\xac\xf1\x19\x89\xc5\x19\x85,~@VNfqI!kl![l!;\xc4\x9c\xa2\xcc\xbc\xf4b\xbd\xcc\xbc\x92\xd4\xf4\xd4"\xae\xdc\xc4\xec\xd4x\x18\xa7\x90#\x94\xd1\xb05\xa8\x9031\xb14I\x0f\x00\xf6\xae)7')        # indirect doctest for unpickling a Tableau_class element
            [[1]]
            sage: loads(dumps( Tableau([[1]]) ))  # indirect doctest for unpickling a Tableau element
            [[1]]
        """
        if isinstance(state, dict):   # for old pickles from Tableau_class
            self._set_parent(Tableaux())
            self.__dict__ = state
        else:
            self._set_parent(state[0])
            self.__dict__ = state[1]

    def _latex_(self):
        r"""
        Returns a LaTeX version of self.

        EXAMPLES::

            sage: latex(Tableau([[1,1,2],[2,3],[3]]))    # indirect doctest
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{ccc}
            \cline{1-1}\cline{2-2}\cline{3-3}
            \lr{1}&\lr{1}&\lr{2}\\
            \cline{1-1}\cline{2-2}\cline{3-3}
            \lr{2}&\lr{3}\\
            \cline{1-1}\cline{2-2}
            \lr{3}\\
            \cline{1-1}
            \end{array}$}
            }
        """
        if len(self) == 0:
            return "{\emptyset}"
        return self._tex_from_array()

    def _tex_from_array(self):
        r"""
        EXAMPLES::

            sage: print Tableau([[1,2],[3,4]])._tex_from_array()
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{cc}
            \cline{1-1}\cline{2-2}
            \lr{1}&\lr{2}\\
            \cline{1-1}\cline{2-2}
            \lr{3}&\lr{4}\\
            \cline{1-1}\cline{2-2}
            \end{array}$}
            }
        """
        import output
        m = max(len(self), len(self[0]))
        return output.tex_from_array(self)

    def __div__(self, t):
        """
        Returns the skew tableau self/t.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[3,4],[5]])
            sage: t/[1,1]
            [[None, 2, 3], [None, 4], [5]]
            sage: t/[3,1]
            [[None, None, None], [None, 4], [5]]
            sage: t/[2,1,1,1]
            Traceback (most recent call last):
            ...
            ValueError: the shape of the tableau must contain the partition
        """

        #if t is a list, convert to to a partition first
        if isinstance(t, list):
            t = sage.combinat.partition.Partition(t)

        #Check to make sure that tableau shape contains t
        if not self.shape().contains(t):
            raise ValueError, "the shape of the tableau must contain the partition"

        st = copy.deepcopy(self._list)

        for i in range(len(t)):
            for j in range(t[i]):
                st[i][j] = None

        return sage.combinat.skew_tableau.SkewTableau(st)

    def __call__(self, *cell):
        r"""

        INPUT:

        - ``self`` -- a tableau
        - ``cell`` -- a pair of integers, tuple, or list specifying a cell in
          the tableau

        OUTPUT:

        - The value in the corresponding cell.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: t(1,0)
            4
            sage: t((1,0))
            4
            sage: t(3,3)
            Traceback (most recent call last):
            ...
            IndexError: The cell (3,3) is not contained in [[1, 2, 3], [4, 5]]
        """
        if isinstance(cell[0], (int, Integer)):
            i,j = cell[0], cell[1]
        else:
            i,j = cell[0]
        try:
            return self[i][j]
        except IndexError:
            raise IndexError, \
                  "The cell (%d,%d) is not contained in %s"%(i,j,self)


    def shape(self):
        r"""
        Returns the shape of a tableau t.

        EXAMPLES::

            sage: Tableau([[1,2,3],[4,5],[6]]).shape()
            [3, 2, 1]
        """

        return sage.combinat.partition.Partition([len(row) for row in self])

    def size(self):
        """
        Returns the size of the shape of the tableau t.

        EXAMPLES::

            sage: Tableau([[1, 4, 6], [2, 5], [3]]).size()
            6
            sage: Tableau([[1, 3], [2, 4]]).size()
            4
        """
        return sum([len(row) for row in self])

    def corners(self):
        """
        Returns the corners of the tableau t.

        EXAMPLES::

            sage: Tableau([[1, 4, 6], [2, 5], [3]]).corners()
            [(0, 2), (1, 1), (2, 0)]
            sage: Tableau([[1, 3], [2, 4]]).corners()
            [(1, 1)]
        """
        return self.shape().corners()

    def conjugate(self):
        """
        Returns the conjugate of the tableau t.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).conjugate()
            [[1, 3], [2, 4]]
        """
        conj_shape = self.shape().conjugate()

        conj = [[None]*row_length for row_length in conj_shape]

        for i in range(len(conj)):
            for j in range(len(conj[i])):
                conj[i][j] = self[j][i]


        return Tableau(conj)


    def pp(self):
        """
        Returns a pretty print string of the tableau.

        EXAMPLES::

            sage: Tableau([[1,2,3],[3,4],[5]]).pp()
              1  2  3
              3  4
              5
        """
        print '\n'.join([ "".join(map(lambda x: "%3s"%str(x) , row))  for row in self])

    def to_word_by_row(self):
        """
        Returns a word obtained from a row reading of the tableau t.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).to_word_by_row()
            word: 3412
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word_by_row()
            word: 325146
        """
        from sage.combinat.words.word import Word
        w = []
        for row in reversed(self):
            w += row
        return Word(w)

    def to_word_by_column(self):
        """
        Returns the word obtained from a column reading of the tableau t.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).to_word_by_column()
            word: 3142
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word_by_column()
            word: 321546
        """
        from sage.combinat.words.word import Word
        w = []
        conj = self.conjugate()
        for row in conj:
            w += list(reversed(row))
        return Word(w)

    def to_word(self):
        """
        An alias for to_word_by_row.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).to_word()
            word: 3412
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word()
            word: 325146
        """
        return self.to_word_by_row()


    def to_permutation(self):
        """
        Returns a permutation with the entries of self obtained by reading
        self in the reading order.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).to_permutation()
            [3, 4, 1, 2]
        """
        return permutation.Permutation(self.to_word())

    def descents(self):
        """
        Returns a list of the cells (i,j) such that
        self[i][j] > self[i-1][j].

        EXAMPLES::

            sage: Tableau( [[1,4],[2,3]] ).descents()
            [(1, 0)]
            sage: Tableau( [[1,2],[3,4]] ).descents()
            [(1, 0), (1, 1)]
        """
        descents = []
        for i in range(1,len(self)):
            for j in range(len(self[i])):
                if self[i][j] > self[i-1][j]:
                    descents.append((i,j))
        return descents

    def major_index(self):
        """
        Returns the major index of self. The major index is defined to be
        the sum of the number of descents of self and the sum of their
        legs length.

        EXAMPLES

        ::

            sage: Tableau( [[1,4],[2,3]] ).major_index()
            1
            sage: Tableau( [[1,2],[3,4]] ).major_index()
            2
        """
        descents = self.descents()
        p = self.shape()
        return len(descents) + sum([ p.leg_length(*d) for d in descents])

    def attacking_pairs(self):
        """
        Returns a list of the attacking pairs of self. An pair of cells (c,
        d) is said to be attacking if one of the following conditions
        hold:

        1. c and d lie in the same row with c to the west of d

        2. c is in the row immediately to the south of d and c
           lies strictly east of d.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[2,5]])
            sage: t.attacking_pairs()
            [((0, 0), (0, 1)),
             ((0, 0), (0, 2)),
             ((0, 1), (0, 2)),
             ((1, 0), (1, 1)),
             ((1, 1), (0, 0))]
        """
        attacking_pairs = []
        for i in range(len(self)):
            for j in range(len(self[i])):
                #c is in position (i,j)
                #Find the d that satisfy condition 1
                for k in range(j+1,len(self[i])):
                    attacking_pairs.append( ((i,j),(i,k)) )

                #Find the d that satisfy condition 2
                if i == 0:
                    continue
                for k in range(j):
                    attacking_pairs.append( ((i,j),(i-1,k)) )

        return attacking_pairs

    def inversions(self):
        """
        Returns a list of the inversions of self. An inversion is an
        attacking pair (c,d) such that the entry of c in self is greater
        than the entry of d.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[2,5]])
            sage: t.inversions()
            [((1, 1), (0, 0))]
        """
        inversions = []
        for (c,d) in self.attacking_pairs():
            if self.entry(c) > self.entry(d):
                inversions.append( (c,d) )
        return inversions

    def inversion_number(self):
        """
        Returns the inversion number of self.

        The inversion number is defined to be the number of inversion of
        self minus the sum of the arm lengths of the descents of self.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[2,5]])
            sage: t.inversion_number()
            0
        """
        p = self.shape()
        return len(self.inversions()) - sum([ p.arm_length(*cell) for cell in self.descents() ])

    def schuetzenberger_involution(self, n = None):
        """
        Returns the Schuetzenberger involution of the tableau self.
        This method relies on the analogous method on words, which reverts the word
        and then complements all letters within the underlying ordered alphabet.
        If `n` is specified, the underlying alphabet is assumed to be `[1,2,\ldots,n]`.
        If no alphabet is specified, `n` is the maximal letter appearing in self.

        INPUT:

        - ``self`` -- a tableau
        - ``n``    -- an integer specifying the maximal letter in the alphabet (optional)

        OUTPUT:

        - a tableau, the Schuetzenberger involution of self

        EXAMPLES::

           sage: t = Tableau([[1,1,1],[2,2]])
           sage: t.schuetzenberger_involution(3)
           [[2, 2, 3], [3, 3]]

            sage: t = Tableau([[1,2,3],[4,5]])
            sage: t.schuetzenberger_involution()
            [[1, 2, 5], [3, 4]]

            sage: t = Tableau([[1,3,5,7],[2,4,6],[8,9]])
            sage: t.schuetzenberger_involution()
            [[1, 2, 6, 8], [3, 4, 9], [5, 7]]

            sage: t = Tableau([])
            sage: t.schuetzenberger_involution()
            []
        """
        w = self.to_word()
        if w.length() == 0:
            return self
        wi = w.schuetzenberger_involution(n=n)
        t = Tableau([[wi[0]]])
        for k in range(1, w.length()):
            t = t.bump(wi[k])
        return t

    def entries(self):
        """
        Returns a list of all entries of self, in the order obtained
        by reading across the rows.

        EXAMPLES::

            sage: t = Tableau([[1,3], [2]])
            sage: t.entries()
            [1, 3, 2]
        """
        return sum(self, [])

    def entry(self, cell):
        """
        Returns the entry of cell in self. Cell is a tuple (i,j) of
        coordinates.

        EXAMPLES::

            sage: t = Tableau([[1,2],[3,4]])
            sage: t.entry( (0,0) )
            1
            sage: t.entry( (1,1) )
            4
        """
        i,j = cell
        return self[i][j]

    def weight(self):
        """
        Returns the weight of the word corresponding to the tableau ``self``.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).weight()
            [1, 1, 1, 1]
        """
        ed = self.to_word().evaluation_dict()
        entries = ed.keys()
        m = max(entries) + 1 if entries else -1
        return [ed.get(k,0) for k in range(1,m)]

    evaluation = weight

    def is_standard(self):
        """
        Returns True if t is a standard tableau and False otherwise.

        EXAMPLES::

            sage: Tableau([[1, 3], [2, 4]]).is_standard()
            True
            sage: Tableau([[1, 2], [2, 4]]).is_standard()
            False
            sage: Tableau([[2, 3], [2, 4]]).is_standard()
            False
            sage: Tableau([[5, 3], [2, 4]]).is_standard()
            False
        """
        t = self
        #Check to make sure the first position is 1
        fillings = []
        for row in t:
            fillings += row
        fillings.sort()
        if fillings != range(1, t.size()+1):
            return False

        #Check to make sure it is increasing along the rows
        for row in t:
            for i in range(1, len(row)):
                if row[i] <= row[i-1]:
                    return False

        #Check to make sure it is increasing along the columns
        conj = t.conjugate()
        for row in conj:
            for i in range(1, len(row)):
                if row[i] <= row[i-1]:
                    return False

        return True

    def is_rectangular(self):
        """
        Returns True if the tableau t is rectangular and False otherwise.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).is_rectangular()
            True
            sage: Tableau([[1,2,3],[4,5],[6]]).is_rectangular()
            False
        """
        width = len(self[0])
        for row in self:
            if len(row) != width:
                return False
        return True

    def vertical_flip(self):
        """
        Returns the tableau obtained by vertically flipping the tableau t.
        This only works for rectangular tableau.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).vertical_flip()
            [[3, 4], [1, 2]]
        """

        if not self.is_rectangular():
            raise TypeError, "the tableau must be rectangular to use verticl_flip()"

        return Tableau([row for row in reversed(self)])

    def rotate_180(self):
        """
        Returns the tableau obtained by rotating t by 180 degrees.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).rotate_180()
            [[4, 3], [2, 1]]
        """
        if not self.is_rectangular():
            raise TypeError, "the tableau must be rectangular to use verticl_flip()"

        return Tableau([ [l for l in reversed(row)] for row in reversed(self) ])

    def cells(self):
        """
        Returns a list of the coordinates of the cells of self.

        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).cells()
            [(0, 0), (0, 1), (1, 0), (1, 1)]
        """
        s = []
        for i in range(len(self)):
            s += [ (i,j) for j in range(len(self[i])) ]
        return s

    def cells_containing(self, i):
        r"""
        Returns the list of cells in which the letter `i` appears in the tableau `self`.
        The list is ordered with cells appearing from left to right.

        EXAMPLES::

            sage: t = Tableau([[1,1,3],[2,3,5],[4,5]])
            sage: t.cells_containing(5)
            [(2, 1), (1, 2)]
            sage: t.cells_containing(4)
            [(2, 0)]
            sage: t.cells_containing(6)
            []

            sage: t = Tableau([[1,1,2,4],[2,4,4],[4]])
            sage: t.cells_containing(4)
            [(2, 0), (1, 1), (1, 2), (0, 3)]
        """
        list = []
        for r in range(len(self)):
            for c in range(self.shape()[r]-1,-1,-1):
                if self[r][c] == i:
                    list += [(r,c)]
        list.reverse()
        return list

    def k_weight(self, k):
        """
        Returns the k-weight of self.

        EXAMPLES::

            sage: Tableau([[1,2],[2,3]]).k_weight(1)
            [1, 1, 1]
            sage: Tableau([[1,2],[2,3]]).k_weight(2)
            [1, 2, 1]
            sage: t = Tableau([[1,1,1,2,5],[2,3,6],[3],[4]])
            sage: t.k_weight(1)
            [2, 1, 1, 1, 1, 1]
            sage: t.k_weight(2)
            [3, 2, 2, 1, 1, 1]
            sage: t.k_weight(3)
            [3, 1, 2, 1, 1, 1]
            sage: t.k_weight(4)
            [3, 2, 2, 1, 1, 1]
            sage: t.k_weight(5)
            [3, 2, 2, 1, 1, 1]
        """
        res = []
        w = self.weight()
        s = self.cells()

        for l in range(1,len(w)+1):
            new_s = [(i,j) for i,j in s if self[i][j] == l]

            #If there are no elements that meet the condition
            if new_s == []:
                res .append(0)
                continue
            x = uniq([ (i-j)%(k+1) for i,j in new_s ])
            res.append(len(x))

        return res


    def restrict(self, n):
        """
        Returns the restriction of the standard tableau to n.

        EXAMPLES::

            sage: Tableau([[1,2],[3],[4]]).restrict(3)
            [[1, 2], [3]]
            sage: Tableau([[1,2],[3],[4]]).restrict(2)
            [[1, 2]]
            sage: Tableau([[1,1],[2]]).restrict(1)
            [[1, 1]]
        """
        res = [ [y for y in row if y <=n] for row in self]
        return Tableau([row for row in res if row != []])

    def to_chain(self):
        """
        Returns the chain of partitions corresponding to the (semi)standard
        tableau.

        EXAMPLES::

            sage: Tableau([[1,2],[3],[4]]).to_chain()
            [[], [1], [2], [2, 1], [2, 1, 1]]
            sage: Tableau([[1,1],[2]]).to_chain()
            [[], [2], [2, 1]]
            sage: Tableau([[1,1],[3]]).to_chain()
            [[], [2], [2], [2, 1]]
            sage: Tableau([]).to_chain()
            [[]]
        """
        if self == []:
            return [self.shape()]
        m = max(self.to_word())
        return [self.restrict(k).shape() for k in range(m+1)]


    def anti_restrict(self, n):
        """
        Returns the skew tableau formed by removing all of the cells from
        self that are filled with a number less than `n`.

        EXAMPLES::

            sage: t = Tableau([[1,2,3],[4,5]]); t
            [[1, 2, 3], [4, 5]]
            sage: t.anti_restrict(1)
            [[None, 2, 3], [4, 5]]
            sage: t.anti_restrict(2)
            [[None, None, 3], [4, 5]]
            sage: t.anti_restrict(3)
            [[None, None, None], [4, 5]]
            sage: t.anti_restrict(4)
            [[None, None, None], [None, 5]]
        """
        t = list(copy.deepcopy(self))

        for row in xrange(len(t)):
            for col in xrange(len(t[row])):
                if t[row][col] <= n:
                    t[row][col] = None
        return sage.combinat.skew_tableau.SkewTableau( t )


    def up(self):
        """
        An iterator for all the tableaux that can be obtained from self by
        adding a cell.

        EXAMPLES::

            sage: t = Tableau([[1,2]])
            sage: [x for x in t.up()]
            [[[1, 2, 3]], [[1, 2], [3]]]
        """
        #Get a list of all places where we can add a cell
        #to the shape of self

        outside_corners = self.shape().outside_corners()

        n = self.size()

        #Go through and add n+1 to the end of each
        #of the rows
        for row, _ in outside_corners:
            new_t = map(list, self)
            if row != len(self):
                new_t[row] += [n+1]
            else:
                new_t.append([n+1])
            yield Tableau(new_t)

    def up_list(self):
        """
        Returns a list of all the tableaux that can be obtained from self
        by adding a cell.

        EXAMPLES::

            sage: t = Tableau([[1,2]])
            sage: t.up_list()
            [[[1, 2, 3]], [[1, 2], [3]]]
        """
        return list(self.up())

    def down(self):
        """
        An iterator for all the tableaux that can be obtained from self by
        removing a cell. Note that this iterates just over a single
        tableaux. EXAMPLES::

            sage: t = Tableau([[1,2],[3]])
            sage: [x for x in t.down()]
            [[[1, 2]]]
        """
        yield self.restrict( self.size() - 1 )

    def down_list(self):
        """
        Returns a list of all the tableaux that can be obtained from self
        by removing a cell. Note that this is just a single tableaux.

        EXAMPLES::

            sage: t = Tableau([[1,2],[3]])
            sage: t.down_list()
            [[[1, 2]]]
        """
        return list(self.down())

    def to_list(self):
        """
        EXAMPLES::

            sage: t = Tableau([[1,2],[3,4]])
            sage: l = t.to_list(); l
            [[1, 2], [3, 4]]
            sage: l[0][0] = 2
            sage: t
            [[1, 2], [3, 4]]
        """
        return [row[:] for row in self]

    def bump(self, x):
        """
        Schensted's row-bumping (or row-insertion) algorithm.

        EXAMPLES::

            sage: t = Tableau([[1,2],[3]])
            sage: t.bump(1)
            [[1, 1], [2], [3]]
            sage: t
            [[1, 2], [3]]
            sage: t.bump(2)
            [[1, 2, 2], [3]]
            sage: t.bump(3)
            [[1, 2, 3], [3]]
            sage: t
            [[1, 2], [3]]
            sage: t = Tableau([[1,2,2,3],[2,3,5,5],[4,4,6],[5,6]])
            sage: t.bump(2)
            [[1, 2, 2, 2], [2, 3, 3, 5], [4, 4, 5], [5, 6, 6]]
        """
        new_t = self.to_list()
        to_insert = x
        row = 0
        done = False
        while not done:
            #if we are at the end of the tableau
            #add to_insert as the last row
            if row == len(new_t):
                new_t.append([to_insert])
                break

            i = 0
            #try to insert to_insert into row
            while i < len(new_t[row]):
                if to_insert < new_t[row][i]:
                    t = to_insert
                    to_insert = new_t[row][i]
                    new_t[row][i] = t
                    break
                i += 1


            #if we haven't already inserted to_insert
            #append it to the end of row
            if i == len(new_t[row]):
                new_t[row].append(to_insert)
                done = True

            row += 1

        return Tableau(new_t)



    def schensted_insert(self, i, left=False):
        """
        EXAMPLES::

            sage: t = Tableau([[3,5],[7]])
            sage: t.schensted_insert(8)
            [[3, 5, 8], [7]]
            sage: t.schensted_insert(8, left=True)
            [[3, 5], [7], [8]]
        """
        if left:
            return self._left_schensted_insert(i)
        else:
            return self._right_schensted_insert(i)

    def _right_schensted_insert(self, letter):
        """
        EXAMPLES::

            sage: t = Tableau([[3,5],[7]])
            sage: t._right_schensted_insert(8)
            [[3, 5, 8], [7]]
            sage: t._right_schensted_insert(2)
            [[2, 5], [3], [7]]
            sage: t = Tableau([[3,8],[7]])
            sage: t._right_schensted_insert(6)
            [[3, 6], [7, 8]]
        """
        h = self.height()
        if h == 0:
            return Tableau([[letter]])
        h += 1
        rep = self.to_list() + [[]]

        for i in range(h):
            j = len(rep[i]) - 1
            while j >= 0 and rep[i][j] > letter:
                j -= 1
            if j == len(rep[i])-1:
                rep[i].append(letter)
                break
            else:
                new_letter = rep[i][j+1]
                rep[i][j+1] = letter
                letter = new_letter

        return Tableau([ row for row in rep if row != []])

    def _left_schensted_insert(self, letter):
        """
        EXAMPLES::

            sage: t = Tableau([[3,5],[7]])
            sage: t._left_schensted_insert(8)
            [[3, 5], [7], [8]]
            sage: t._left_schensted_insert(6)
            [[3, 5], [6, 7]]
            sage: t._left_schensted_insert(2)
            [[2, 3, 5], [7]]
        """
        h = len(self)
        if h == 0:
            return Tableau([[letter]])
        h1 = h + 1
        rep = self.to_list()
        rep.reverse()

        width = len(rep[h-1])
        heights = self._heights() + [h1]

        for j in range(1, width+2):
            i = heights[j-1]
            while i != h1 and rep[i-1][j-1] >= letter:
                i += 1
            if i == heights[j-1]: #add on top of column j
                if j == 1:
                    rep = [[letter]] + rep
                else:
                    rep[i-2].append(letter)
                break
            elif i == h1 and j == width: #add on right of line i
                if rep[i-2][j-1] < letter:
                    rep[i-2].append(letter)
                else:
                    new_letter = rep[i-2][j-1]
                    rep[i-2][j-1] = letter
                    rep[i-2].append(new_letter)
                break
            else:
                new_letter = rep[i-2][j-1]
                rep[i-2][j-1] = letter
                letter = new_letter

        rep.reverse()
        return Tableau(rep)


    def insert_word(self, w, left=False):
        """
        EXAMPLES::

            sage: t0 = Tableau([])
            sage: w = [1,1,2,3,3,3,3]
            sage: t0.insert_word(w)
            [[1, 1, 2, 3, 3, 3, 3]]
            sage: t0.insert_word(w,left=True)
            [[1, 1, 2, 3, 3, 3, 3]]
            sage: w.reverse()
            sage: t0.insert_word(w)
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t0.insert_word(w,left=True)
            [[1, 1, 3, 3], [2, 3], [3]]
        """
        if left:
            w = [i for i in reversed(w)]
        res = self
        for i in w:
            res = res.schensted_insert(i,left=left)
        return res

    def bump_multiply(left, right):
        """
        Multiply two tableaux using Schensted's bump.

        This product makes the set of tableaux into an associative monoid.
        The empty tableaux is the unit in this monoid.

        Fulton, William. 'Young Tableaux' p11-12

        EXAMPLES::

            sage: t = Tableau([[1,2,2,3],[2,3,5,5],[4,4,6],[5,6]])
            sage: t2 = Tableau([[1,2],[3]])
            sage: t.bump_multiply(t2)
            [[1, 1, 2, 2, 3], [2, 2, 3, 5], [3, 4, 5], [4, 6, 6], [5]]
        """
        if not isinstance(right, Tableau):
            raise TypeError, "right must be a Tableau"

        row = len(right)
        product = copy.deepcopy(left)
        while row > 0:
            row -= 1
            for i in right[row]:
                product = product.bump(i)
        return product

    def slide_multiply(left, right):
        """
        Multiply two tableaux using jeu de taquin.

        This product makes the set of tableaux into an associative monoid.
        The empty tableaux is the unit in this monoid.

        Fulton, William. 'Young Tableaux' p15

        EXAMPLES::

            sage: t = Tableau([[1,2,2,3],[2,3,5,5],[4,4,6],[5,6]])
            sage: t2 = Tableau([[1,2],[3]])
            sage: t.slide_multiply(t2)
            [[1, 1, 2, 2, 3], [2, 2, 3, 5], [3, 4, 5], [4, 6, 6], [5]]
        """
        st = []
        if len(left) == 0:
            return right
        else:
            l = len(left[0])

        for row in range(len(right)):
            st.append([None]*l + right[row])
        for row in range(len(left)):
            st.append(left[row])

        return sage.combinat.skew_tableau.SkewTableau(st).rectify()

    def _slide_up(self, c):
        r"""
        Auxiliary method used for promotion, which removes cell `c` from self,
        slides the letters of self up using jeu de taquin slides, and then fills the empty
        cell at `(0,0)` with the value 0.

        TESTS::

            sage: t = Tableau([[1,1,2],[2,3,5],[4,5]])
            sage: t._slide_up((2,1))
            [[0, 1, 2], [1, 3, 5], [2, 4]]

            sage: t._slide_up((1,2))
            [[0, 1, 2], [1, 2, 3], [4, 5]]

            sage: t = Tableau([[1,1,3],[2,3,5],[4,5]])
            sage: t._slide_up((1,2))
            [[0, 1, 1], [2, 3, 3], [4, 5]]
        """
        new_st = [x[:] for x in self]
        spotl, spotc = c
        while [spotl, spotc] != [0,0]:
            #once moving box is in first column, just move letters up
            if spotc == 0:
                new_st[spotl][spotc] = new_st[spotl-1][spotc]
                spotl -= 1
                continue
            #once moving box is in first row, just move letters up
            elif spotl == 0:
                new_st[spotl][spotc] = new_st[spotl][spotc-1]
                spotc -= 1
                continue
            else:
                #If we get to this stage, we need to compare
                below = new_st[spotl-1][spotc]
                left = new_st[spotl][spotc-1]
                if below >= left:
                    #Swap with the cell below
                    new_st[spotl][spotc] = new_st[spotl-1][spotc]
                    spotl -= 1
                    continue
                else:
                    #Swap with the cell to the left
                    new_st[spotl][spotc] = new_st[spotl][spotc-1]
                    spotc -= 1
                    continue
        #set box in position (0,0) to 0
        new_st[0][0] = 0
        return Tableau(new_st)

    def promotion_inverse(self, n):
        """
        Inverse promotion operator defined on rectangular tableaux using
        jeu de taquin

        EXAMPLES::

            sage: t = Tableau([[1,2],[3,3]])
            sage: t.promotion_inverse(2)
            [[1, 2], [2, 3]]
            sage: t = Tableau([[1,2],[2,3]])
            sage: t.promotion_inverse(2)
            [[1, 1], [2, 3]]
        """
        if not self.is_rectangular():
            raise ValueError, "Tableau is not rectangular"
        s = self.shape()[0]
        l = self.weight()[0]
        word = [i for i in self.to_word() if i>1]
        word = [i-1 for i in word]
        t = Tableau([])
        t = t.insert_word(word)
        t = t.to_list()
        if l<s:
            for i in range(l):
                t[len(t)-1].append(n+1)
        else:
            t.append([n+1 for i in range(s)])
        return Tableau(t)

    def promotion(self, n):
        r"""
        Returns the promotion operator acting on `self` by removing all letters `n+1`
        from tableau `t` (one by one from left to right), applying jeu de taquin to move
        boxes into the empty box, filling 0 into the empty boxes in the first row, and
        finally adding one to each letter.

        EXAMPLES::

            sage: t = Tableau([[1,2],[3,3]])
            sage: t.promotion(2)
            [[1, 1], [2, 3]]

            sage: t = Tableau([[1,1,1],[2,2,3],[3,4,4]])
            sage: t.promotion(3)
            [[1, 1, 2], [2, 2, 3], [3, 4, 4]]

            sage: t = Tableau([[1,2],[2]])
            sage: t.promotion(3)
            [[2, 3], [3]]

            sage: t = Tableau([[1,1,3],[2,2]])
            sage: t.promotion(2)
            [[1, 2, 2], [3, 3]]

            sage: t = Tableau([[1,1,3],[2,3]])
            sage: t.promotion(2)
            [[1, 1, 2], [2, 3]]
        """
        if self.is_rectangular():
            t = self.rotate_180()
            t = [[n+2-i for i in row] for row in t.to_list()]
            t = Tableau(t).promotion_inverse(n)
            t = [[n+2-i for i in row] for row in t.to_list()]
            return Tableau(t).rotate_180()
        p = self
        for c in self.cells_containing(n+1):
            p = p._slide_up(c)
        return Tableau([[i+1 for i in row] for row in p])

    def row_stabilizer(self):
        """
        Return the PermutationGroup corresponding to the row stabilizer of
        self.

        EXAMPLES::

            sage: rs = Tableau([[1,2,3],[4,5]]).row_stabilizer()
            sage: rs.order() == factorial(3)*factorial(2)
            True
            sage: PermutationGroupElement([(1,3,2),(4,5)]) in rs
            True
            sage: PermutationGroupElement([(1,4)]) in rs
            False
            sage: rs = Tableau([[1, 2],[3]]).row_stabilizer()
            sage: PermutationGroupElement([(1,2),(3,)]) in rs
            True
            sage: rs.one().list()
            [1, 2, 3]
            sage: rs = Tableau([[1],[2],[3]]).row_stabilizer()
            sage: rs.order()
            1
        """

        # Ensure that the permutations involve all elements of the
        # tableau, by including the identity permutation on the set [1..k].
        k = max(self.entries())
        gens = [range(1,k+1)]
        for i in range(len(self)):
            for j in range(0, len(self[i])-1):
                gens.append( (self[i][j], self[i][j+1]) )
        return PermutationGroup( gens )

    def column_stabilizer(self):
        """
        Return the PermutationGroup corresponding to the column stabilizer
        of self.

        EXAMPLES::

            sage: cs = Tableau([[1,2,3],[4,5]]).column_stabilizer()
            sage: cs.order() == factorial(2)*factorial(2)
            True
            sage: PermutationGroupElement([(1,3,2),(4,5)]) in cs
            False
            sage: PermutationGroupElement([(1,4)]) in cs
            True
        """

        return self.conjugate().row_stabilizer()

    def height(self):
        """
        Returns the height of the tableau.

        EXAMPLES::

            sage: Tableau([[1,2,3],[4,5]]).height()
            2
            sage: Tableau([[1,2,3]]).height()
            1
            sage: Tableau([]).height()
            0
        """
        return len(self)

    def _heights(self):
        """
        EXAMPLES::

            sage: Tableau([[1,2,3,4],[5,6],[7],[8]])._heights()
            [1, 3, 4, 4]
            sage: Tableau([])._heights()
            []
            sage: Tableau([[1]])._heights()
            [1]
            sage: Tableau([[1,2]])._heights()
            [1, 1]
            sage: Tableau([[1,2],[3],[4]])._heights()
            [1, 3]
        """
        cor = self.corners()
        ncor = len(cor)
        if ncor == 0:
            return []
        k = len(self)
        cor = [ [k-i,j+1]  for i,j in reversed(cor)]

        heights = [1]*(cor[0][1])
        for i in range(1, ncor):
            heights += [ cor[i][0] ]*(cor[i][1]-cor[i-1][1])

        return heights

    def last_letter_lequal(self, tab2):
        """
        Returns True if self is less than or equal to tab2 in the last
        letter ordering.

        EXAMPLES::

            sage: st = StandardTableaux([3,2])
            sage: f = lambda b: 1 if b else 0
            sage: matrix( [ [ f(t1.last_letter_lequal(t2)) for t2 in st] for t1 in st] )
            [1 1 1 1 1]
            [0 1 1 1 1]
            [0 0 1 1 1]
            [0 0 0 1 1]
            [0 0 0 0 1]
        """
        n = self.size()
        if not isinstance(tab2, Tableau):
            try:
                tab2 = Tableau(tab2)
            except StandardError:
                raise TypeError, "tab2 must be a standard tableau"

        if tab2.size() != n:
            raise ValueError, "tab2 must be the same size as self"

        if self == tab2:
            return True

        for j in range(n, 1, -1):
            self_j_pos = None
            for i in range(len(self)):
                if j in self[i]:
                    self_j_pos = i
                    break

            tab2_j_pos = None
            for i in range(len(tab2)):
                if j in tab2[i]:
                    tab2_j_pos = i
                    break

            if self_j_pos < tab2_j_pos:
                return True
            if tab2_j_pos < self_j_pos:
                return False

    def charge(self):
        r"""
        Returns the charge of the reading word of self.  See
        :meth:`sage.combinat.words.finite_word.charge` for more information.

        EXAMPLES::

            sage: Tableau([[1,1],[2,2],[3]]).charge()
            0
            sage: Tableau([[1,1,3],[2,2]]).charge()
            1
            sage: Tableau([[1,1,2],[2],[3]]).charge()
            1
            sage: Tableau([[1,1,2],[2,3]]).charge()
            2
            sage: Tableau([[1,1,2,3],[2]]).charge()
            2
            sage: Tableau([[1,1,2,2],[3]]).charge()
            3
            sage: Tableau([[1,1,2,2,3]]).charge()
            4
        """
        return self.to_word().charge()

    def cocharge(self):
        r"""
        Returns the cocharge of the reading word of self.  See
        :meth:`sage.combinat.words.finite_word.cocharge` for more information.

        EXAMPLES::

            sage: Tableau([[1,1],[2,2],[3]]).cocharge()
            4
            sage: Tableau([[1,1,3],[2,2]]).cocharge()
            3
            sage: Tableau([[1,1,2],[2],[3]]).cocharge()
            3
            sage: Tableau([[1,1,2],[2,3]]).cocharge()
            2
            sage: Tableau([[1,1,2,3],[2]]).cocharge()
            2
            sage: Tableau([[1,1,2,2],[3]]).cocharge()
            1
            sage: Tableau([[1,1,2,2,3]]).cocharge()
            0
        """
        return self.to_word().cocharge()


    ##############
    # katabolism #
    ##############

    def katabolism(self):
        """
        EXAMPLES::

            sage: Tableau([]).katabolism()
            []
            sage: Tableau([[1,2,3,4,5]]).katabolism()
            [[1, 2, 3, 4, 5]]
            sage: Tableau([[1,1,3,3],[2,3],[3]]).katabolism()
            [[1, 1, 2, 3, 3, 3], [3]]
            sage: Tableau([[1, 1, 2, 3, 3, 3], [3]]).katabolism()
            [[1, 1, 2, 3, 3, 3, 3]]
        """
        h = self.height()
        if h == 0:
            return self
        else:
            #Remove the top row and insert it back in
            return Tableau(self[1:]).insert_word(self[0],left=True)

    def katabolism_sequence(self):
        """
        EXAMPLES::

            sage: t = Tableau([[1,2,3,4,5,6,8],[7,9]])
            sage: t.katabolism_sequence()
            [[[1, 2, 3, 4, 5, 6, 8], [7, 9]],
             [[1, 2, 3, 4, 5, 6, 7, 9], [8]],
             [[1, 2, 3, 4, 5, 6, 7, 8], [9]],
             [[1, 2, 3, 4, 5, 6, 7, 8, 9]]]
        """
        h = self.height()
        res = [self]
        while h != 1:
            res.append( res[-1].katabolism() )
            h = res[-1].height()
        return res

    def lambda_katabolism(self, part):
        r"""
        For a partition ``lambda`` and a tableau ``T``, the
        ``lambda``-katabolism of ``T`` is defined by performing the following
        steps.

        1. Truncate the parts of ``lambda`` so that ``lambda`` is contained
        in the shape of ``T``.  Let ``m`` be the length of this partition.

        2. Let ``T_a`` be the first ``m`` rows of ``T``, and ``T_b`` be the
        remaining rows.

        3. Let ``S_a`` be the skew tableau ``T_a / lambda``.

        4. Concatenate the reading words of ``S_a`` and ``T_b``, and insert
        into a tableau.

        EXAMPLES::

            sage: Tableau([[1,1,3],[2,4,5]]).lambda_katabolism([2,1])
            [[3, 5], [4]]
            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: t.lambda_katabolism([])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t.lambda_katabolism([1])
            [[1, 2, 3, 3, 3], [3]]
            sage: t.lambda_katabolism([1,1])
            [[1, 3, 3, 3], [3]]
            sage: t.lambda_katabolism([2,1])
            [[3, 3, 3, 3]]
            sage: t.lambda_katabolism([4,2,1])
            []
            sage: t.lambda_katabolism([5,1])
            [[3, 3]]
            sage: t.lambda_katabolism([4,1])
            [[3, 3]]
        """
        #Reduce the partition if it is too big for the tableau
        part  = [ min(part[i],len(self[i])) for i in range(min(len(self), len(part))) ]
        if self.shape() == part:
            return Tableau([])

        m = len(part)

        w1 = flatten([row for row in reversed(self[m:])])

        w2 = []
        for i,row in enumerate(reversed(self[:m])):
            w2 += row[ part[-1-i] : ]

        return Tableau([]).insert_word(w2+w1)


    def reduced_lambda_katabolism(self, part):
        """
        EXAMPLES::

            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: t.reduced_lambda_katabolism([])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t.reduced_lambda_katabolism([1])
            [[1, 2, 3, 3, 3], [3]]
            sage: t.reduced_lambda_katabolism([1,1])
            [[1, 3, 3, 3], [3]]
            sage: t.reduced_lambda_katabolism([2,1])
            [[3, 3, 3, 3]]
            sage: t.reduced_lambda_katabolism([4,2,1])
            []
            sage: t.reduced_lambda_katabolism([5,1])
            0
            sage: t.reduced_lambda_katabolism([4,1])
            0
        """
        part1 = part

        if self == []:
            return self

        res = self.lambda_katabolism(part)

        if res == []:
            return res

        if res == 0:
            return 0

        a = self[0][0]

        part = [ min(part1[i], len(self[i])) for i in range(min(len(part1),len(self)))]
        tt_part = Tableau([ [a+i]*part[i] for i in range(len(part)) ])
        t_part = Tableau([[self[i][j] for j in range(part[i])] for i in range(len(part))])

        if t_part == tt_part:
            return res
        else:
            return 0

    def katabolism_projector(self, parts):
        """
        EXAMPLES::

            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: t.katabolism_projector([[4,2,1]])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t.katabolism_projector([[1]])
            []
            sage: t.katabolism_projector([[2,1],[1]])
            []
            sage: t.katabolism_projector([[1,1],[4,1]])
            [[1, 1, 3, 3], [2, 3], [3]]
        """
        res = self
        for p in parts:
            res = res.reduced_lambda_katabolism(p)
            if res == 0:
                return 0

        if res == []:
            return self
        else:
            return Tableau([])


    def promotion_operator(self, i):
        """
        EXAMPLES::

            sage: t = Tableau([[1,2],[3]])
            sage: t.promotion_operator(1)
            [[[1, 2], [3], [4]], [[1, 2], [3, 4]], [[1, 2, 4], [3]]]
            sage: t.promotion_operator(2)
            [[[1, 1], [2, 3], [4]],
             [[1, 1, 2], [3], [4]],
             [[1, 1, 4], [2, 3]],
             [[1, 1, 2, 4], [3]]]
            sage: Tableau([[1]]).promotion_operator(2)
            [[[1, 1], [2]], [[1, 1, 2]]]
            sage: Tableau([[1,1],[2]]).promotion_operator(3)
            [[[1, 1, 1], [2, 2], [3]],
             [[1, 1, 1, 2], [2], [3]],
             [[1, 1, 1, 3], [2, 2]],
             [[1, 1, 1, 2, 3], [2]]]

        TESTS::

            sage: Tableau([]).promotion_operator(2)
            [[[1, 1]]]
            sage: Tableau([]).promotion_operator(1)
            [[[1]]]
        """
        chain = self.to_chain()
        part = self.shape()
        weight = self.weight()
        perm = permutation.from_reduced_word(range(1, len(weight)+1))
        l = part.add_horizontal_border_strip(i)
        ltab = [ from_chain( chain + [next] ) for next in l ]
        return [ x.symmetric_group_action_on_values(perm) for x in ltab]


    ##################################
    # actions on tableaux from words #
    ##################################
    def raise_action_from_words(self, f, *args):
        """
        EXAMPLES::

            sage: from sage.combinat.tableau import symmetric_group_action_on_values
            sage: import functools
            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: f = functools.partial(t.raise_action_from_words, symmetric_group_action_on_values)
            sage: f([1,2,3])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: f([3,2,1])
            [[1, 1, 1, 1], [2, 3], [3]]
            sage: f([1,3,2])
            [[1, 1, 2, 2], [2, 2], [3]]
        """
        w = self.to_word()
        w = f(w, *args)
        return from_shape_and_word(self.shape(), w)

    def symmetric_group_action_on_values(self, perm):
        """
        EXAMPLES::

            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: t.symmetric_group_action_on_values([1,2,3])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t.symmetric_group_action_on_values([3,2,1])
            [[1, 1, 1, 1], [2, 3], [3]]
            sage: t.symmetric_group_action_on_values([1,3,2])
            [[1, 1, 2, 2], [2, 2], [3]]
        """
        return self.raise_action_from_words(symmetric_group_action_on_values, perm)

    #########
    # atoms #
    #########
    def socle(self):
        """
        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).socle()
            2
            sage: Tableau([[1,2,3,4]]).socle()
            4
        """
        h = self.height()
        if h == 0:
            return 0
        w1row = self[0]
        i = 0
        while i < len(w1row)-1:
            if w1row[i+1] != w1row[i] + 1:
                break
            i += 1
        return i+1

    def atom(self):
        """
        EXAMPLES::

            sage: Tableau([[1,2],[3,4]]).atom()
            [2, 2]
            sage: Tableau([[1,2,3],[4,5],[6]]).atom()
            [3, 2, 1]
        """
        ll = [ t.socle() for t in self.katabolism_sequence() ]
        lres = ll[:]
        for i in range(1,len(ll)):
            lres[i] = ll[i] - ll[i-1]
        return lres

class SemistandardTableau(Tableau):
    """
    A class to model a semistandard tableau.

    INPUT:

    - ``t`` -- a Tableau, a list of lists, or an empty list

    OUTPUT:

    - A SemistandardTableau object constructed from ``t``.

    A semistandard tableau is a tableau whose entries are positive integers,
    which are weakly increasing in rows and strictly increasing down columns.

    Note that Sage uses the English convention for partitions and
    tableaux; the longer rows are displayed on top.

    EXAMPLES::

        sage: t = SemistandardTableau([[1,2,3],[2,3]]); t
        [[1, 2, 3], [2, 3]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty print
        1 2 3
        2 3
        sage: t = Tableau([[1,2],[2]])
        sage: s = SemistandardTableau(t); s
        [[1, 2], [2]]
        sage: SemistandardTableau([]) # The empty tableau
        []

    When using code that will generate a lot of tableaux, it is slightly more
    efficient to construct a SemistandardTableau from the appropriate Parent object::

        sage: SST = SemistandardTableaux()
        sage: SST([[1, 2, 3], [4, 5]])
        [[1, 2, 3], [4, 5]]

    SEE ALSO:

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`

    TESTS::

        sage: SemistandardTableau([[1,2,3],[1]])
        Traceback (most recent call last):
        ...
        ValueError: Columns are not strictly increasing.

        sage: SemistandardTableau([[1,2,1]])
        Traceback (most recent call last):
        ...
        ValueError: Rows are not weakly increasing.

        sage: SemistandardTableau([[0,1]])
        Traceback (most recent call last):
        ...
        ValueError: Entries must be non-negative integers
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a SemistandardTableau is only ever constructed as an
        element_class call of an appropriate parent.

        TESTS::

            sage: t = SemistandardTableau([[1,1],[2]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Semistandard tableaux
            sage: t.category()
            Category of elements of Semistandard tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.SemistandardTableaux_all_with_category.element_class'>
        """
        if isinstance(t, SemistandardTableau):
            return t

        return SemistandardTableaux_all().element_class(SemistandardTableaux_all(), t)


    def __init__(self, parent, t):
        r"""
        Initializes a semistandard tableau.

        TESTS::

            sage: t = Tableaux()([[1,1],[2]])
            sage: s = SemistandardTableaux(3)([[1,1],[2]])
            sage: s==t
            True
            sage: s.parent()
            Semistandard tableaux of size 3 and maximum entry 3
            sage: r = SemistandardTableaux(3)(t); r.parent()
            Semistandard tableaux of size 3 and maximum entry 3
            sage: isinstance(r, Tableau)
            True
        """
        Tableau.__init__(self, parent, t)

        # Verify row-weak condition and integer entries
        for row in t:
            if not isinstance(row[0], (int, Integer)) or row[0]<=0:
                raise ValueError, "Entries must be non-negative integers"
            for i in range(1, len(row)):
                if not isinstance(row[i], (int, Integer)):
                    raise ValueError, "Entries must be integers"
                if row[i-1] > row[i]:
                    raise ValueError, "Rows are not weakly increasing."


        # Verify col-strict condition
        for i in range(len(t) - 1):
            M = len(t[i+1])
            for j,x in enumerate(t[i]):
                if (j < M and x >= t[i+1][j]):
                    raise ValueError, "Columns are not strictly increasing."

class StandardTableau(SemistandardTableau):
    """
    A class to model a standard tableau.

    INPUT:

    - ``t`` -- a Tableau, a list of lists, or an empty list

    OUTPUT:

    - A StandardTableau object constructed from ``t``.

    A standard tableau is a semistandard tableau whose entries are exactly the
    positive integers from 1 to ``n``, where ``n`` is the size of the tableau.

    EXAMPLES::

        sage: t = StandardTableau([[1,2,3],[4,5]]); t
        [[1, 2, 3], [4, 5]]
        sage: t.shape()
        [3, 2]
        sage: t.pp() # pretty print
        1 2 3
        4 5
        sage: t.is_standard()
        True
        sage: StandardTableau([]) # The empty tableau
        []

    When using code that will generate a lot of tableaux, it is slightly more
    efficient to construct a StandardTableau from the appropriate Parent object::

        sage: ST = StandardTableaux()
        sage: ST([[1, 2, 3], [4, 5]])
        [[1, 2, 3], [4, 5]]

    SEE ALSO:

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`

        sage: StandardTableau([[1,2,3],[4,4]])
        Traceback (most recent call last):
        ...
        ValueError: The entries must consist of exactly the numbers 1..n
        sage: StandardTableau([[1,3,2]])
        Traceback (most recent call last):
        ...
        ValueError: Rows are not weakly increasing.
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a StandardTableau is only ever constructed as an
        element_class call of an appropriate parent.

        TESTS::

            sage: t = StandardTableau([[1,2],[3]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Standard tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.StandardTableaux_all_with_category.element_class'>
        """
        if isinstance(t, StandardTableau):
            return t

        return StandardTableaux_all().element_class(StandardTableaux_all(), t)

    def __init__(self, parent, t):
        r"""
        Initializes a standard tableau.

        TESTS::

            sage: t = Tableaux()([[1,2],[3]])
            sage: s = StandardTableaux(3)([[1,2],[3]])
            sage: s==t
            True
            sage: s.parent()
            Standard tableaux of size 3
            sage: r = StandardTableaux(3)(t); r.parent()
            Standard tableaux of size 3
            sage: isinstance(r, Tableau)
            True
        """
        SemistandardTableau.__init__(self, parent, t)

        # Verify that the entries are exactly the numbers 1..n
        entries = []
        for row in self:
            entries += row
        entries.sort()
        if len(entries)>0 and entries != range(1, max(entries)+1):
            raise ValueError, "The entries must consist of exactly the numbers 1..n"

    def content(self, k):
        """
        Returns the content of ``k`` in a standard tableau. That is, if
        ``k`` appears in row `r` and column `c` of the tableau then we
        return `c-r`.

        EXAMPLES::

            sage: StandardTableau([[1,2],[3,4]]).content(3)
            -1

            sage: StandardTableau([[1,2],[3,4]]).content(6)
            Traceback (most recent call last):
            ...
            ValueError: 6 does not appear in tableau

        """
        for r in range(len(self)):
          try:
            c=self[r].index(k)
            return c-r
          except ValueError:
            pass
        raise ValueError, "%d does not appear in tableau"%k

def from_chain(chain):
    """
    Returns a semistandard tableau from a chain of partitions.

    EXAMPLES::

        sage: from sage.combinat.tableau import from_chain
        sage: from_chain([[], [2], [2, 1], [3, 2, 1]])
        [[1, 1, 3], [2, 3], [3]]
    """
    res = [[0]*chain[-1][i] for i in range(len(chain[-1]))]
    for i in reversed(range(2, len(chain)+1)):
        for j in range(len(chain[i-1])):
            for k in range(chain[i-1][j]):
                res[j][k] = i -1
    return Tableau(res)

def from_shape_and_word(shape, w, order = "French"):
    r"""
    Returns a tableau from a shape and word.

    INPUT:

    - ``shape`` -- a partition
    - ``w`` -- a word whose length equals that of the partition
    - ``order`` -- a string which can take values "French" or "English"; the default is "French"

    OUTPUT:

    A tableau, whose shape is ``shape`` and whose reading word is ``w``.
    If the order is specified to "French", the reading word is to be read starting
    from the top row in French notation (= the bottom row in English notation).
    If the order is specified to "English", the reading word is to be read starting with the
    top row in English notation.

    EXAMPLES::

        sage: from sage.combinat.tableau import from_shape_and_word
        sage: t = Tableau([[1, 3], [2], [4]])
        sage: shape = t.shape(); shape
        [2, 1, 1]
        sage: word = t.to_word(); word
        word: 4213
        sage: from_shape_and_word(shape, word)
        [[1, 3], [2], [4]]
        sage: word = Word(flatten(t))
        sage: from_shape_and_word(shape, word, order = "English")
        [[1, 3], [2], [4]]
    """
    res = []
    j = 0
    if order == "French":
        shape = reversed(shape)
    for l in shape:
        res.append( list(w[j:j+l]) )
        j += l
    if order == "French":
        res.reverse()
    return Tableau(res)

class Tableaux(UniqueRepresentation, Parent):
    """
    A factory class for the various classes of tableaux.

    INPUT:

    - ``n`` (optional) -- a non-negative integer

    OUTPUT:

    - If ``n`` is specified, the class of tableaux of size ``n``. Otherwise,
      the class of all tableaux.

    A tableau in Sage is a finite list of lists, whose lengths are weakly
    decreasing, or an empty list, representing the empty tableau.  The entries
    of a tableau can be any sage object. Because of this, no enumeration
    through the set of Tableaux is possible.

    EXAMPLES::

        sage: T = Tableaux(); T
        Tableaux
        sage: T3 = Tableaux(3); T3
        Tableaux of size 3
        sage: [['a','b']] in T
        True
        sage: [['a','b']] in T3
        False
        sage: t = T3([[1,1,1]]); t
        [[1, 1, 1]]
        sage: t in T
        True
        sage: t.parent()
        Tableaux of size 3
        sage: T([]) # the empty tableau
        []
        sage: T.category()
        Category of sets

    SEE ALSO:

        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`

    TESTS::

        sage: t = Tableaux(3)([[1,2],[3]])
        sage: t.parent()
        Tableaux of size 3
        sage: Tableaux(t)
        Traceback (most recent call last):
        ...
        ValueError: The argument to Tableaux() must be a non-negative integer.
        sage: Tableaux(3)([[1, 1]])
        Traceback (most recent call last):
        ...
        ValueError: [[1, 1]] is not an element of Tableaux of size 3.

        sage: t0 = Tableau([[1]])
        sage: t1 = Tableaux()([[1]])
        sage: t2 = Tableaux()(t1)
        sage: t0 == t1 == t2
        True
        sage: t1 in Tableaux()
        True
        sage: t1 in Tableaux(1)
        True
        sage: t1 in Tableaux(2)
        False

        sage: [[1]] in Tableaux()
        True
        sage: [] in Tableaux(0)
        True
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`Tableaux` for more
        information.

        TESTS::

            sage: Tableaux()
            Tableaux
            sage: Tableaux(3)
            Tableaux of size 3
        """
        if args:
            n = args[0]
        elif 'n' in kwargs:
            n = kwargs[n]
        else:
            n = None

        if n == None:
            return Tableaux_all()
        else:
            if not isinstance(n,(int, Integer)) or n < 0:
                raise ValueError, "The argument to Tableaux() must be a non-negative integer."
            return Tableaux_size(n)

    Element = Tableau

    def _element_constructor_(self, t):
        r"""
        Constructs an object from t as an element of self, if possible. This
        is inherited by all Tableaux, SemistandardTableaux, and
        StandardTableaux classes.

        INPUT:

        - ``t`` -- Data which can be interpreted as a tableau

        OUTPUT:

        - The corresponding tableau object

        TESTS::

            sage: T = Tableaux(3)
            sage: T([[1,2,1]]).parent() is T     # indirect doctest
            True
            sage: T( StandardTableaux(3)([[1, 2, 3]])).parent() is T
            True
            sage: T([[1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 2]] is not an element of Tableaux of size 3.
        """
        if not t in self:
            raise ValueError, "%s is not an element of %s."%(t, self)

        return self.element_class(self, t)

class Tableaux_all(Tableaux):

    def __init__(self):
        r"""
        Initializes the class of all tableaux

        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_all()
            sage: TestSuite(T).run()

        """
        super(Tableaux_all, self).__init__(category=Sets())

    def __contains__(self, x):
        """
        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_all()
            sage: [[1,2],[3,4]] in T
            True
            sage: [[1,2],[3]] in T
            True
            sage: [] in T
            True
            sage: [['a','b']] in T
            True
            sage: Tableau([['a']]) in T
            True

            sage: [1,2,3] in T
            False
            sage: [[1],[1,2]] in T
            False
        """
        if isinstance(x, self.element_class):
            return True
        try:
            x = Tableau(x)
        except ValueError:
            return False
        return True

    def _repr_(self):
        """
        TESTS::

            sage: repr(Tableaux())    # indirect doctest
            'Tableaux'
        """
        return "Tableaux"

    def an_element(self):
        r"""
        Returns a particular element of the class.

        TESTS::

            sage: T = Tableaux()
            sage: T.an_element()
            [[1, 1], [1]]
        """
        return self.element_class(self, [[1, 1], [1]])


    def list(self):
        """
        TESTS::

            sage: Tableaux().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __iter__(self):
        """
        TESTS::

            sage: iter(Tableaux())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError



class Tableaux_size(Tableaux):
    def __init__(self, n):
        r"""
        Initializes the class of tableaux of size n

        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_size(3)
            sage: TestSuite(T).run()

            sage: T = sage.combinat.tableau.Tableaux_size(0)
            sage: TestSuite(T).run()
        """
        super(Tableaux_size, self).__init__(category=Sets())
        self.size = n

    def __contains__(self,x):
        """
        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_size(3)
            sage: [[2,4], [1]] in T
            True

            sage: [[2,4],[1,3]] in T
            False
        """
        if isinstance(x, self.element_class) and sum(map(len,x))==self.size:
            return True
        try:
            x = Tableau(x)
        except ValueError:
            return False
        return sum(map(len,x))==self.size

    def _repr_(self):
        """
        TESTS::

            sage: repr(Tableaux(4))    # indirect doctest
            'Tableaux of size 4'
        """
        return "Tableaux of size %s"%self.size

    def an_element(self):
        r"""
        Returns a particular element of the class.

        TESTS::

            sage: T = sage.combinat.tableau.Tableaux_size(3)
            sage: T.an_element()
            [[1, 1], [1]]
            sage: T = sage.combinat.tableau.Tableaux_size(0)
            sage: T.an_element()
            []
        """
        if self.size==0:
            return self.element_class(self, [])

        if self.size==1:
            return self.element_class(self, [[1]])

        return self.element_class(self, [[1]*(self.size-1),[1]])

    def list(self):
        """
        TESTS::

            sage: Tableaux(3).list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __iter__(self):
        """
        TESTS::

            sage: iter(Tableaux(3))
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError



##########################
# Semi-standard tableaux #
##########################
class SemistandardTableaux(Tableaux):
    """
    A factory class for the various classes of semistandard tableaux.

    INPUT:

    Keyword arguments:

    - ``size`` -- The size of the tableaux
    - ``shape`` -- The shape of the tableaux
    - ``eval`` -- The weight (also called content or weight) of the tableaux
    - `max_entry` -- A maximum entry for the tableaux.  This can be a positive
      integer or infinity (oo). If ``size`` or ``shape`` are specified, `max_entry`
      defaults to be ``size`` or the size of ``shape``.

    Positional arguments:

    - The first argument is interpreted as either ``size`` or ``shape`` according to
      whether it is an integer or a partition
    - The second keyword argument will always be interpreted as ``eval``

    OUTPUT:

    - The appropriate class, after checking basic consistency tests. (For
      example, specifying ``eval`` implies a value for `max_entry`).

    A semistandard tableau is a tableau whose entries are positive integers,
    which are weakly increasing in rows and strictly increasing down columns.
    Note that Sage uses the English convention for partitions and tableaux;
    the longer rows are displayed on top.

    Classes of semistandard tableaux can be iterated over if and only if there is some
    restriction.

    EXAMPLES::

        sage: SST = SemistandardTableaux([2,1]); SST
        Semistandard tableaux of shape [2, 1] and maximum entry 3
        sage: SST.list()
        [[[1, 1], [2]],
         [[1, 1], [3]],
         [[1, 2], [2]],
         [[1, 2], [3]],
         [[1, 3], [2]],
         [[1, 3], [3]],
         [[2, 2], [3]],
         [[2, 3], [3]]]

        sage: SST = SemistandardTableaux(3); SST
        Semistandard tableaux of size 3 and maximum entry 3
        sage: SST.list()
        [[[1, 1, 1]],
         [[1, 1, 2]],
         [[1, 1, 3]],
         [[1, 2, 2]],
         [[1, 2, 3]],
         [[1, 3, 3]],
         [[2, 2, 2]],
         [[2, 2, 3]],
         [[2, 3, 3]],
         [[3, 3, 3]],
         [[1, 1], [2]],
         [[1, 1], [3]],
         [[1, 2], [2]],
         [[1, 2], [3]],
         [[1, 3], [2]],
         [[1, 3], [3]],
         [[2, 2], [3]],
         [[2, 3], [3]],
         [[1], [2], [3]]]

        sage: SST = SemistandardTableaux(3, max_entry=2); SST
        Semistandard tableaux of size 3 and maximum entry 2
        sage: SST.list()
        [[[1, 1, 1]],
         [[1, 1, 2]],
         [[1, 2, 2]],
         [[2, 2, 2]],
         [[1, 1], [2]],
         [[1, 2], [2]]]

        sage: SST = SemistandardTableaux(3, max_entry=oo); SST
        Semistandard tableaux of size 3
        sage: SST[123]
        [[3, 4], [6]]

        sage: SemistandardTableaux(max_entry=2)[11]
        [[1, 1], [2]]

        sage: SemistandardTableaux()[0]
        []

    SEE ALSO:

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`SemistandardTableaux` for more
        information.

        TESTS::

            sage: SemistandardTableaux()
            Semistandard tableaux
            sage: SemistandardTableaux(3)
            Semistandard tableaux of size 3 and maximum entry 3
            sage: SemistandardTableaux(size=3)
            Semistandard tableaux of size 3 and maximum entry 3
            sage: SemistandardTableaux(0)
            Semistandard tableaux of size 0 and maximum entry 0
            sage: SemistandardTableaux([2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux(shape=[2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux([])
            Semistandard tableaux of shape [] and maximum entry 0
            sage: SemistandardTableaux(eval=[2,1])
            Semistandard tableaux of size 3 and weight [2, 1]
            sage: SemistandardTableaux(max_entry=3)
            Semistandard tableaux with maximum entry 3
            sage: SemistandardTableaux(3, [2,1])
            Semistandard tableaux of size 3 and weight [2, 1]
            sage: SemistandardTableaux(3, shape=[2,1])
            Semistandard tableaux of shape [2, 1] and maximum entry 3
            sage: SemistandardTableaux(3, [2,1], shape=[2,1])
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: SemistandardTableaux(3, max_entry=4)
            Semistandard tableaux of size 3 and maximum entry 4
            sage: SemistandardTableaux(3, max_entry=oo)
            Semistandard tableaux of size 3
            sage: SemistandardTableaux([2, 1], max_entry=oo)
            Semistandard tableaux of shape [2, 1]
            sage: SemistandardTableaux([2, 1], [2, 1])
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: mu = Partition([2,1]); SemistandardTableaux(mu, mu)
            Semistandard tableaux of shape [2, 1] and weight [2, 1]
            sage: SemistandardTableaux(3, [2, 1], max_entry=2)
            Semistandard tableaux of size 3 and weight [2, 1]

            sage: SemistandardTableaux(3, shape=[2])
            Traceback (most recent call last):
            ...
            ValueError: size and shape are different sizes

            sage: SemistandardTableaux(3, [2])
            Traceback (most recent call last):
            ...
            ValueError: size and eval are different sizes

            sage: SemistandardTableaux([2],[3])
            Traceback (most recent call last):
            ...
            ValueError: shape and eval are different sizes

            sage: SemistandardTableaux(2,[2], max_entry=4)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must match the weight

            sage: SemistandardTableaux(eval=[2], max_entry=oo)
            Traceback (most recent call last):
            ...
            ValueError: the maximum entry must match the weight

            sage: SemistandardTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: shape must be a partition
        """
        # Process the keyword arguments -- allow for original syntax where
        #   n == size,  p== shape and mu == eval
        n = kwargs.get('n', None)
        size = kwargs.get('size', n)

        p = kwargs.get('p', None)
        shape = kwargs.get('shape', p)

        mu = kwargs.get('eval', None)
        mu = kwargs.get("mu", mu)

        max_entry = kwargs.get('max_entry', None)

        # Process the positional arguments
        if args:
            # The first arg could be either a size or a shape
            if isinstance(args[0], (int, Integer)):
                if size is not None:
                    raise ValueError, "size was specified more than once"
                else:
                    size = args[0]
            else:
                if shape is not None:
                    raise ValueError, "the shape was specified more than once"
                shape = args[0] # we check it's a partition later

        if len(args) == 2:
            # The second non-keyword argument is the weight
            if mu is not None:
                raise ValueError, "the weight was specified more than once"
            else:
                mu = args[1]


        # Consistency checks
        if size is not None:
            if not isinstance(size, (int, Integer)):
                raise ValueError, "size must be an integer"
            elif size < 0:
                raise ValueError, "size must be non-negative"

        if shape is not None:
            # use in (and not isinstance) below so that lists can be used as
            # shorthand
            if not shape in sage.combinat.partition.Partitions():
                raise ValueError, "shape must be a partition"
            shape = sage.combinat.partition.Partition(shape)

        if mu is not None:
            if (not mu in sage.combinat.composition.Compositions()) and\
                    (not mu in sage.combinat.partition.Partitions()):
                raise ValueError, "mu must be a composition"
            mu = sage.combinat.composition.Composition(mu)

        is_inf = max_entry is PlusInfinity()

        if max_entry is not None:
            if not is_inf and not isinstance(max_entry, (int, Integer)):
                raise ValueError, "max_entry must be an integer or PlusInfinity"
            elif max_entry <= 0:
                raise ValueError, "max_entry must be positive"

        if (mu is not None) and (max_entry is not None):
            if max_entry != len(mu):
                raise ValueError, "the maximum entry must match the weight"

        if (size is not None) and (shape is not None):
            if sum(shape) != size:
                # This could return an empty class instead of an error
                raise ValueError, "size and shape are different sizes"

        if (size is not None) and (mu is not None):
            if sum(mu) != size:
                # This could return an empty class instead of an error
                raise ValueError, "size and eval are different sizes"

        # Dispatch appropriately
        if (shape is not None) and (mu is not None):
            if sum(shape) != sum(mu):
                # This could return an empty class instead of an error
                raise ValueError, "shape and eval are different sizes"
            else:
                return SemistandardTableaux_shape_weight(shape, mu)

        if (shape is not None):
            if is_inf:
                return SemistandardTableaux_shape_inf(shape)
            return SemistandardTableaux_shape(shape, max_entry)

        if (mu is not None):
            return SemistandardTableaux_size_weight(sum(mu), mu)

        if (size is not None):
            if is_inf:
                return SemistandardTableaux_size_inf(size)
            return SemistandardTableaux_size(size, max_entry)

        return SemistandardTableaux_all(max_entry)

    Element = SemistandardTableau


    def __getitem__(self, r):
        r"""
        The default implementation of __getitem__ for enumerated sets does not
        allow slices so we override it.

        EXAMPLES::

            sage: StandardTableaux([4,3,3,2])[10:20]     # indirect doctest
            [[[1, 3, 9, 12], [2, 5, 10], [4, 6, 11], [7, 8]],
             [[1, 2, 9, 12], [3, 5, 10], [4, 6, 11], [7, 8]],
             [[1, 3, 9, 12], [2, 4, 10], [5, 6, 11], [7, 8]],
             [[1, 2, 9, 12], [3, 4, 10], [5, 6, 11], [7, 8]],
             [[1, 5, 8, 12], [2, 6, 10], [3, 7, 11], [4, 9]],
             [[1, 4, 8, 12], [2, 6, 10], [3, 7, 11], [5, 9]],
             [[1, 3, 8, 12], [2, 6, 10], [4, 7, 11], [5, 9]],
             [[1, 2, 8, 12], [3, 6, 10], [4, 7, 11], [5, 9]],
             [[1, 4, 8, 12], [2, 5, 10], [3, 7, 11], [6, 9]],
             [[1, 3, 8, 12], [2, 5, 10], [4, 7, 11], [6, 9]]]

            sage: SemistandardTableaux(size=2, max_entry=oo)[5]
            [[2, 3]]

            sage: SemistandardTableaux([2,1], max_entry=oo)[3]
            [[1, 2], [3]]

            sage: SemistandardTableaux(3, max_entry=2)[0:5]    # indirect doctest
            [[[1, 1, 1]],
            [[1, 1, 2]],
            [[1, 2, 2]],
            [[2, 2, 2]],
            [[1, 1], [2]]]

            sage: SemistandardTableaux([2,2], [2, 1, 1])[0]    # indirect doctest
            [[1, 1], [2, 3]]

            sage: SemistandardTableaux([1,1,1], max_entry=4)[0:4]
            [[[1], [2], [3]],
             [[1], [2], [4]],
             [[1], [3], [4]],
             [[2], [3], [4]]]

            sage: SemistandardTableaux(3, [2,1])[1]    # indirect doctest
            [[1, 1], [2]]

            sage: StandardTableaux(3)[:]  # indirect doctest
            [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]

            sage: StandardTableaux([2,2])[1]   # indirect doctest
            [[1, 2], [3, 4]]

        TESTS:

            sage: SemistandardTableaux()[5]
            [[1], [2]]

            sage: SemistandardTableaux(max_entry=2)[5]
            [[2, 2]]

            sage: SemistandardTableaux()[:]
            Traceback (most recent call last):
            ...
            ValueError: infinite set

            sage: SemistandardTableaux(size=2, max_entry=oo)[:]
            Traceback (most recent call last):
            ...
            ValueError: infinite set
        """
        if isinstance(r,(int,Integer)):
            return self.unrank(r)
        elif isinstance(r,slice):
            start=0 if r.start is None else r.start
            stop=r.stop
            if stop is None and not self.is_finite():
                raise ValueError, 'infinite set'
        else:
            raise ValueError, 'r must be an integer or a slice'
        count=0
        tabs=[]
        for t in self:
            if count==stop:
                break
            if count>=start:
                tabs.append(t)
            count+=1

        # this is to cope with empty slices endpoints like [:6] or [:}
        if count==stop or stop is None:
            return tabs
        raise IndexError, 'value out of range'


class SemistandardTableaux_all(DisjointUnionEnumeratedSets, SemistandardTableaux):
    def __init__(self, max_entry=None):
        r"""
        Initializes the class of all semistandard tableaux. Input is not
        checked; please use :class:`SemistandardTableaux` to ensure the
        options are properly parsed.

        TESTS::

            sage: T = sage.combinat.tableau.SemistandardTableaux_all()
            sage: TestSuite(T).run()

            sage: T=sage.combinat.tableau.SemistandardTableaux_all(max_entry=3)
            sage: TestSuite(T).run()
        """
        if max_entry is not PlusInfinity():
            self.max_entry = max_entry
            SST_n = lambda n: SemistandardTableaux_size(n, max_entry)
            DisjointUnionEnumeratedSets.__init__( self,
                    Family(NonNegativeIntegers(), SST_n),
                    facade=True, keepkey = False)

        else:
            self.max_entry = None


    def __contains__(self, t):
        """
        Returns true if ``t`` can be interpreted as a SemistandardTableau

        TESTS::

            sage: T = sage.combinat.tableau.SemistandardTableaux_all()
            sage: [[1,2],[2]] in T
            True
            sage: [] in T
            True
            sage: Tableau([[1]]) in T
            True
            sage: StandardTableau([[1]]) in T
            True

            sage: [[1,2],[1]] in T
            False
            sage: [[1,1],[5]] in T
            True
        """
        if isinstance(t, self.element_class):
            return self.max_entry is None or \
                    t == [] or \
                    max(flatten(t)) <= self.max_entry

        try:
            t = SemistandardTableau(t)
        except ValueError:
            return False

        return self.max_entry is None or \
                t == [] or \
                max(flatten(t)) <= self.max_entry


    def _repr_(self):
        """
        TESTS::

            sage: SemistandardTableaux()    # indirect doctest
            Semistandard tableaux

            sage: SemistandardTableaux(max_entry=3)
            Semistandard tableaux with maximum entry 3
        """
        if self.max_entry is not None:
            return "Semistandard tableaux with maximum entry %s"%str(self.max_entry)
        return "Semistandard tableaux"


    def list(self):
        """
        TESTS::

            sage: SemistandardTableaux().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class SemistandardTableaux_size_inf(SemistandardTableaux):
    def __init__(self, n):
        r"""
        Initializes the class of semistandard tableaux of size ``n`` with no
        maximum entry. Input is not checked; please use
        :class:`SemistandardTableaux` to ensure the options are properly
        parsed.

        TESTS::

            sage: T = sage.combinat.tableau.SemistandardTableaux_size_inf(3)
            sage: TestSuite(T).run()
        """
        super(SemistandardTableaux_size_inf, self).__init__(
              category = InfiniteEnumeratedSets())
        self.size = n


    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux(3, max_entry=oo))    # indirect doctest
            'Semistandard tableaux of size 3'
        """
        return "Semistandard tableaux of size %s"%str(self.size)

    def __contains__(self, t):
        """
        Returns true if ``t`` can be interpreted as an element of the class.

        TESTS::

            sage: T = SemistandardTableaux(3, max_entry=oo)
            sage: [[1,2],[5]] in T
            True
            sage: StandardTableau([[1, 2], [3]]) in T
            True

            sage: [] in T
            False
            sage: Tableau([[1]]) in T
            False
        """
        if isinstance(t, self.element_class):
            return sum(map(len, t)) == self.size

        try:
            t = SemistandardTableau(t)
        except ValueError:
            return False

        return sum(map(len, t)) == self.size

    def __iter__(self):
        """
        EXAMPLES::

            sage: sst = SemistandardTableaux(3, max_entry=oo)
            sage: [sst[t] for t in range(0,5)]
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 2, 2]],
             [[2, 2, 2]],
             [[1, 1], [2]]]
            sage: sst[1000]
            [[2, 12], [7]]
            sage: sst[0].parent() is sst
            True
        """
        # Iterates through with maximum entry as order
        i = 1
        while(True):
            for part in sage.combinat.partition.Partitions(self.size):
                if i != 1:
                    for k in range(1, self.size+1):
                        for c in IntegerVectors(self.size - k, i-1):
                            c.append(k)
                            for sst in SemistandardTableaux_shape_weight(part,
                                    sage.combinat.composition.Composition(c)):
                                yield self.element_class(self, sst)
                else:
                    for sst in SemistandardTableaux_shape_weight(part,
                            sage.combinat.composition.Composition([self.size])):
                        yield self.element_class(self, sst)
            i += 1


    def list(self):
        """
        TESTS::

            sage: SemistandardTableaux(3, max_entry=oo).list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class SemistandardTableaux_shape_inf(SemistandardTableaux):
    def __init__(self, p):
        r"""
        Initializes the class of semistandard tableaux of shape ``p`` and no
        maximum entry. Input is not checked; please use
        :class:`SemistandardTableaux` to ensure the options are properly
        parsed.

        TESTS::

            sage: SST = SemistandardTableaux([2,1], max_entry=oo)
            sage: type(SST)
            <class 'sage.combinat.tableau.SemistandardTableaux_shape_inf_with_category'>
            sage: TestSuite(SST).run()
        """
        super(SemistandardTableaux_shape_inf, self).__init__(
              category = InfiniteEnumeratedSets())

        self.shape = p


    def __contains__(self, x):
        """
        EXAMPLES::

            sage: SST = SemistandardTableaux([2,1], max_entry=oo)
            sage: [[13, 67], [1467]] in SST
            True
            sage: SST = SemistandardTableaux([3,1], max_entry=oo)
            sage: [[13, 67], [1467]] in SST
            False
        """
        if isinstance(x, self.element_class) and map(len,x)==self.shape:
            return True
        try:
            x = SemistandardTableau(x)
        except ValueError:
            return False
        return map(len,x)==self.shape

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux([2,1], max_entry=oo))    # indirect doctest
            'Semistandard tableaux of shape [2, 1]'
        """
        return "Semistandard tableaux of shape %s" %str(self.shape)


    def __iter__(self):
        """
        An iterator for the semistandard partitions of shape p and no maximum entry.
        Iterates through with maximum entry as order.

        EXAMPLES::

            sage: SST = SemistandardTableaux([3, 1], max_entry=oo)
            sage: SST[1000]
            [[1, 1, 10], [6]]
            sage: [ SST[t] for t in range(0, 5) ]
            [[[1, 1, 1], [2]],
             [[1, 1, 2], [2]],
             [[1, 2, 2], [2]],
             [[1, 1, 1], [3]],
             [[1, 1, 2], [3]]]
            sage: SST[0].parent() is SST
            True
        """
        # Iterates through with maximum entry as order
        i = 1
        n = sum(self.shape)
        while(True):
            if i != 1:
                for k in range(1, n+1):
                    for c in IntegerVectors(n - k, i-1):
                        c.append(k)
                        for sst in SemistandardTableaux_shape_weight(self.shape,
                                sage.combinat.composition.Composition(c)):
                            yield self.element_class(self, sst)
            else:
                for sst in SemistandardTableaux_shape_weight(self.shape,
                                sage.combinat.composition.Composition([n])):
                    yield self.element_class(self, sst)
            i += 1


class SemistandardTableaux_size(SemistandardTableaux):
    def __init__(self, n, max_entry=None):
        r"""
        Initializes the class of semistandard tableaux of size ``n``. Input is
        not checked; please use :class:`SemistandardTableaux` to ensure the
        options are properly parsed.

        TESTS::

            sage: SST = SemistandardTableaux(3); SST
            Semistandard tableaux of size 3 and maximum entry 3
            sage: type(SST)
            <class 'sage.combinat.tableau.SemistandardTableaux_size_with_category'>
            sage: TestSuite(SST).run()

            sage: SST = SemistandardTableaux(3, max_entry=6)
            sage: type(SST)
            <class 'sage.combinat.tableau.SemistandardTableaux_size_with_category'>
            sage: TestSuite(SST).run()
        """
        super(SemistandardTableaux_size, self).__init__(
                  category = FiniteEnumeratedSets())

        self.size = n

        if max_entry is None:
            self.max_entry = n
        else:
            self.max_entry = max_entry


    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux(3))    # indirect doctest
            'Semistandard tableaux of size 3 and maximum entry 3'

            sage: repr(SemistandardTableaux(3, max_entry=6))
            'Semistandard tableaux of size 3 and maximum entry 6'
        """
        return "Semistandard tableaux of size %s and maximum entry %s"%(str(self.size), str(self.max_entry))

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[1,2],[3,3]] in SemistandardTableaux(3)
            False
            sage: [[1,2],[3,3]] in SemistandardTableaux(4)
            True
            sage: [[1,2],[3,3]] in SemistandardTableaux(4, max_entry=2)
            False
            sage: SST = SemistandardTableaux(4)
            sage: all([sst in SST for sst in SST])
            True
        """
        if self.size==0:
            return x == []

        if isinstance(x, self.element_class) and \
                sum(map(len,x)) == self.size and \
                max(flatten(x)) <= self.max_entry:
            return True
        try:
            x = SemistandardTableau(x)
        except ValueError:
            return False
        return sum(map(len,x)) == self.size and \
                max(flatten(x)) <= self.max_entry

    def cardinality(self):
        """
        EXAMPLES::

            sage: SemistandardTableaux(3).cardinality()
            19
            sage: SemistandardTableaux(4).cardinality()
            116
            sage: SemistandardTableaux(4, max_entry=2).cardinality()
            9
            sage: SemistandardTableaux(4, max_entry=10).cardinality()
            4225
            sage: ns = range(1, 6)
            sage: ssts = [ SemistandardTableaux(n) for n in ns ]
            sage: all([sst.cardinality() == len(sst.list()) for sst in ssts])
            True
        """
        c = 0
        for part in sage.combinat.partition.Partitions(self.size):
            c += SemistandardTableaux_shape(part, self.max_entry).cardinality()
        return c


    def __iter__(self):
        """
        EXAMPLES::

            sage: [ t for t in SemistandardTableaux(2) ]
            [[[1, 1]], [[1, 2]], [[2, 2]], [[1], [2]]]
            sage: [ t for t in SemistandardTableaux(3) ]
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 1, 3]],
             [[1, 2, 2]],
             [[1, 2, 3]],
             [[1, 3, 3]],
             [[2, 2, 2]],
             [[2, 2, 3]],
             [[2, 3, 3]],
             [[3, 3, 3]],
             [[1, 1], [2]],
             [[1, 1], [3]],
             [[1, 2], [2]],
             [[1, 2], [3]],
             [[1, 3], [2]],
             [[1, 3], [3]],
             [[2, 2], [3]],
             [[2, 3], [3]],
             [[1], [2], [3]]]

            sage: [ t for t in SemistandardTableaux(3, max_entry=2) ]
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 2, 2]],
             [[2, 2, 2]],
             [[1, 1], [2]],
             [[1, 2], [2]]]

            sage: sst = SemistandardTableaux(3)
            sage: sst[0].parent() is sst
            True
        """
        for part in sage.combinat.partition.Partitions(self.size):
            for sst in SemistandardTableaux_shape(part, self.max_entry):
                yield self.element_class(self, sst)


class SemistandardTableaux_shape_weight(SemistandardTableaux):
    def __init__(self, p, mu):
        r"""
        Initializes the class of all semistandard tableaux of shape ``p`` and
        weight ``mu``. Input is not checked; please use
        :class:`SemistandardTableaux` to ensure the options are properly
        parsed.

        TESTS::

            sage: SST = SemistandardTableaux([2,1], [2,1])
            sage: TestSuite(SST).run()
        """
        super(SemistandardTableaux_shape_weight, self).__init__(
              category = FiniteEnumeratedSets())
        self.shape = p
        self.weight = mu
        self.max_entry = len(mu)


    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux([2,1],[2,1]))    # indirect doctest
            'Semistandard tableaux of shape [2, 1] and weight [2, 1]'
        """
        return "Semistandard tableaux of shape %s and weight %s"%(self.shape, self.weight)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: SST = SemistandardTableaux([2,1], [2,1])
            sage: all([sst in SST for sst in SST])
            True
            sage: len(filter(lambda x: x in SST, SemistandardTableaux(3)))
            1
            sage: SST.cardinality()
            1
        """
        if x not in SemistandardTableaux_shape(self.shape, self.max_entry):
            return False
        n = sum(self.shape)

        if n == 0 and len(x) == 0:
            return True

        content = {}
        for row in x:
            for i in row:
                content[i] = content.get(i, 0) + 1
        content_list = [0]*int(max(content))

        for key in content:
            content_list[key-1] = content[key]

        if content_list != self.weight:
            return False

        return True


    def cardinality(self):
        """
        Returns the number of semistandard tableaux of the given shape and
        weight, as computed by kostka_number function of symmetrica.

        EXAMPLES::

            sage: SemistandardTableaux([2,2], [2, 1, 1]).cardinality()
            1
            sage: SemistandardTableaux([2,2,2], [2, 2, 1,1]).cardinality()
            1
            sage: SemistandardTableaux([2,2,2], [2, 2, 2]).cardinality()
            1
            sage: SemistandardTableaux([3,2,1], [2, 2, 2]).cardinality()
            2
        """
        return symmetrica.kostka_number(self.shape,self.weight)

    def __iter__(self):
        """
        TESTS::

            sage: sst = SemistandardTableaux([3,1],[2,1,1])
            sage: [sst[i] for i in range(2)]
            [[[1, 1, 2], [3]], [[1, 1, 3], [2]]]
            sage: sst[0].parent() is sst
            True
        """
        for t in symmetrica.kostka_tab(self.shape, self.weight):
            yield self.element_class(self, t)


    def list(self):
        """
        EXAMPLES::

            sage: SemistandardTableaux([2,2], [2, 1, 1]).list()
            [[[1, 1], [2, 3]]]
            sage: SemistandardTableaux([2,2,2], [2, 2, 1,1]).list()
            [[[1, 1], [2, 2], [3, 4]]]
            sage: SemistandardTableaux([2,2,2], [2, 2, 2]).list()
            [[[1, 1], [2, 2], [3, 3]]]
            sage: SemistandardTableaux([3,2,1], [2, 2, 2]).list()
            [[[1, 1, 2], [2, 3], [3]], [[1, 1, 3], [2, 2], [3]]]
        """
        return symmetrica.kostka_tab(self.shape, self.weight)


class SemistandardTableaux_shape(SemistandardTableaux):
    def __init__(self, p, max_entry=None):
        r"""
        Initializes the class of semistandard tableaux of shape ``p``, with a
        given max_entry. max_entry defaults to the size of ``p``. Input is not
        checked; please use :class:`SemistandardTableaux` to ensure the
        options are properly parsed.

        TESTS::

            sage: SST = SemistandardTableaux([2,1])
            sage: TestSuite(SST).run()

            sage: SST = SemistandardTableaux([2,1], max_entry=5)
            sage: TestSuite(SST).run()
        """
        super(SemistandardTableaux_shape, self).__init__(
              category = FiniteEnumeratedSets())

        self.shape = p

        if max_entry is None:
            self.max_entry = sum(p)
        else:
            self.max_entry = max_entry


    def __iter__(self):
        """
        An iterator for the semistandard partitions of the specified shape.

        EXAMPLES::

            sage: [ t for t in SemistandardTableaux([3]) ]
            [[[1, 1, 1]],
             [[1, 1, 2]],
             [[1, 1, 3]],
             [[1, 2, 2]],
             [[1, 2, 3]],
             [[1, 3, 3]],
             [[2, 2, 2]],
             [[2, 2, 3]],
             [[2, 3, 3]],
             [[3, 3, 3]]]
            sage: [ t for t in SemistandardTableaux([2,1]) ]
            [[[1, 1], [2]],
             [[1, 1], [3]],
             [[1, 2], [2]],
             [[1, 2], [3]],
             [[1, 3], [2]],
             [[1, 3], [3]],
             [[2, 2], [3]],
             [[2, 3], [3]]]
            sage: [ t for t in SemistandardTableaux([1,1,1]) ]
            [[[1], [2], [3]]]

            sage: [ t for t in SemistandardTableaux([1,1,1], max_entry=4) ]
            [[[1], [2], [3]],
             [[1], [2], [4]],
             [[1], [3], [4]],
             [[2], [3], [4]]]

            sage: sst = SemistandardTableaux([3])
            sage: sst[0].parent() is sst
            True
        """
        for c in IntegerVectors(sum(self.shape), self.max_entry):
            for sst in SemistandardTableaux_shape_weight(self.shape,
                    sage.combinat.composition.Composition(c)):
                yield self.element_class(self, sst)


    def __contains__(self, x):
        """
        EXAMPLES::

            sage: SST = SemistandardTableaux([2,1])
            sage: all([sst in SST for sst in SST])
            True
            sage: len(filter(lambda x: x in SST, SemistandardTableaux(3)))
            8
            sage: SST.cardinality()
            8

            sage: SST = SemistandardTableaux([2,1], max_entry=4)
            sage: all([sst in SST for sst in SST])
            True
            sage: SST.cardinality()
            20
        """
        return x in SemistandardTableaux_all(max_entry=self.max_entry)\
                and map(len, x) == self.shape

    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux([2,1]))    # indirect doctest
            'Semistandard tableaux of shape [2, 1] and maximum entry 3'

            sage: repr(SemistandardTableaux([2,1], max_entry=5))
            'Semistandard tableaux of shape [2, 1] and maximum entry 5'
        """
        return "Semistandard tableaux of shape %s and maximum entry %s" %(str(self.shape), str(self.max_entry))

    def cardinality(self):
        """
        EXAMPLES::

            sage: SemistandardTableaux([2,1]).cardinality()
            8
            sage: SemistandardTableaux([2,2,1]).cardinality()
            75
            sage: SymmetricFunctions(QQ).schur()([2,2,1]).expand(5)(1,1,1,1,1) # cross check
            75
            sage: SemistandardTableaux([5]).cardinality()
            126
            sage: SemistandardTableaux([3,2,1]).cardinality()
            896

            sage: SemistandardTableaux([3,2,1], max_entry=7).cardinality()
            2352
        """
        c = 0
        for comp in IntegerVectors(sum(self.shape), self.max_entry):
            c += SemistandardTableaux_shape_weight(self.shape,
                    sage.combinat.composition.Composition(comp)).cardinality()
        return c

class SemistandardTableaux_size_weight(SemistandardTableaux):
    def __init__(self, n, mu):
        r"""
        Initializes the class of semistandard tableaux of size ``n`` and
        weight ``mu``. Input is not checked; please use
        :class:`SemistandardTableaux` to ensure the options are properly
        parsed.

        TESTS::

            sage: SST = SemistandardTableaux(3, [2,1])
            sage: TestSuite(SST).run()
        """
        super(SemistandardTableaux_size_weight, self).__init__(
              category = FiniteEnumeratedSets())
        self.size = n
        self.weight = mu
        self.max_entry = len(mu)


    def _repr_(self):
        """
        TESTS::

            sage: repr(SemistandardTableaux(3, [2,1]))    # indirect doctest
            'Semistandard tableaux of size 3 and weight [2, 1]'
        """
        return "Semistandard tableaux of size %s and weight %s"%(self.size, self.weight)

    def __iter__(self):
        """
        EXAMPLES::

            sage: [ t for t in SemistandardTableaux(3, [2,1]) ]
            [[[1, 1, 2]], [[1, 1], [2]]]
            sage: [ t for t in SemistandardTableaux(4, [2,2]) ]
            [[[1, 1, 2, 2]], [[1, 1, 2], [2]], [[1, 1], [2, 2]]]
            sage: sst = SemistandardTableaux(4, [2,2])
            sage: sst[0].parent() is sst
            True
        """
        for p in sage.combinat.partition.Partitions(self.size):
            for sst in SemistandardTableaux_shape_weight(p, self.weight):
                yield self.element_class(self, sst)


    def cardinality(self):
        """
        EXAMPLES::

            sage: SemistandardTableaux(3, [2,1]).cardinality()
            2
            sage: SemistandardTableaux(4, [2,2]).cardinality()
            3
        """
        c = 0
        for p in sage.combinat.partition.Partitions(self.size):
            c += SemistandardTableaux_shape_weight(p, self.weight).cardinality()
        return c

    def __contains__(self, x):
        """
        TESTS::

            sage: SST = SemistandardTableaux(6, [2,2,2])
            sage: all([sst in SST for sst in SST])
            True
            sage: all([sst in SST for sst in SemistandardTableaux([3,2,1],[2,2,2])])
            True
        """
        return x in SemistandardTableaux_shape_weight(sage.combinat.partition.Partition(map(len,
            x)), self.weight)

########################
# Standard Tableaux    #
########################

class StandardTableaux(SemistandardTableaux):
    """
    A factory for the various classes of standard tableaux.

    INPUT:

    - Either a non-negative integer (possibly specified with the keyword ``n``)
      or a partition.

    OUTPUT:

    - With no argument, the class of all standard tableaux

    - With a non-negative integer argument, ``n``, the class of all standard
      tableaux of size ``n``

    - With a partition argument, the class of all standard tableaux of that
      shape.

    A standard tableau is a semistandard tableaux which contains each of the
    entries from 1 to ``n`` exactly once.

    All classes of standard tableaux are iterable.

    EXAMPLES::

        sage: ST = StandardTableaux(3); ST
        Standard tableaux of size 3
        sage: ST.first()
        [[1, 2, 3]]
        sage: ST.last()
        [[1], [2], [3]]
        sage: ST.cardinality()
        4
        sage: ST.list()
        [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]

    SEE ALSO:

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`SemistandardTableau`
        - :class:`StandardTableau`

    TESTS::

        sage: StandardTableaux()([])
        []
        sage: ST = StandardTableaux([2,2]); ST
        Standard tableaux of shape [2, 2]
        sage: ST.first()
        [[1, 3], [2, 4]]
        sage: ST.last()
        [[1, 2], [3, 4]]
        sage: ST.cardinality()
        2
        sage: ST.list()
        [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`StandardTableaux` for more
        information.

        TESTS::

            sage: StandardTableaux()
            Standard tableaux
            sage: StandardTableaux(3)
            Standard tableaux of size 3
            sage: StandardTableaux([2,1])
            Standard tableaux of shape [2, 1]
            sage: StandardTableaux(0)
            Standard tableaux of size 0

            sage: StandardTableaux(-1)
            Traceback (most recent call last):
            ...
            ValueError: The argument must be a non-negative integer or a partition.
            sage: StandardTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: The argument must be a non-negative integer or a partition.
        """
        if args:
            n = args[0]
        elif 'n' in kwargs:
            n = kwargs[n]
        else:
            n = None

        if n is None:
            return StandardTableaux_all()

        elif n in sage.combinat.partition.Partitions():
            return StandardTableaux_shape(sage.combinat.partition.Partition(n))

        if not isinstance(n,(int, Integer)) or n < 0:
            raise ValueError, "The argument must be a non-negative integer or a partition."

        return StandardTableaux_size(n)

    Element = StandardTableau

class StandardTableaux_all(DisjointUnionEnumeratedSets, StandardTableaux):
    def __init__(self):
        r"""
        Initializes the class of all standard tableaux.

        TESTS::

            sage: ST = StandardTableaux()
            sage: TestSuite(ST).run()
        """
        DisjointUnionEnumeratedSets.__init__( self,
                Family(NonNegativeIntegers(), StandardTableaux_size),
                facade=True, keepkey = False)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[1,1],[2,3]] in StandardTableaux()
            False
            sage: [[1,2],[3,4]] in StandardTableaux()
            True
            sage: [[1,3],[2,4]] in StandardTableaux()
            True
            sage: [] in StandardTableaux()
            True
        """
        if isinstance(x, self.element_class):
            return True

        try:
            x = StandardTableau(x)
        except ValueError:
            return False

        return True

    def _repr_(self):
        """
        TESTS::

            sage: repr(StandardTableaux())    # indirect doctest
            'Standard tableaux'
        """
        return "Standard tableaux"


class StandardTableaux_size(StandardTableaux):
    def __init__(self, n):
        r"""
        Initializes the class of all standard tableaux of size ``n``. Input is
        not checked; please use :class:`StandardTableaux` to ensure the
        options are properly parsed.

        TESTS::

            sage: ST = StandardTableaux(3)
            sage: TestSuite(ST).run()
            sage: ST = StandardTableaux(0)
            sage: TestSuite(ST).run()
        """
        super(StandardTableaux_size, self).__init__(
              category = FiniteEnumeratedSets())
        self.size = n


    def _repr_(self):
        """
        TESTS::

            sage: repr(StandardTableaux(3))    # indirect doctest
            'Standard tableaux of size 3'
        """
        return "Standard tableaux of size %s"%self.size

    def __contains__(self, x):
        """
        TESTS::

            sage: ST3 = StandardTableaux(3)
            sage: all([st in ST3 for st in ST3])
            True
            sage: ST4 = StandardTableaux(4)
            sage: filter(lambda x: x in ST3, ST4)
            []
        """
        if isinstance(x, self.element_class) and sum(map(len, x)) == self.size:
            return True

        try:
            x = StandardTableau(x)
        except ValueError:
            return False

        return sum(map(len, x)) == self.size


    def __iter__(self):
        """
        EXAMPLES::

            sage: [ t for t in StandardTableaux(1) ]
            [[[1]]]
            sage: [ t for t in StandardTableaux(2) ]
            [[[1, 2]], [[1], [2]]]
            sage: [ t for t in StandardTableaux(3) ]
            [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]
            sage: [ t for t in StandardTableaux(4) ]
            [[[1, 2, 3, 4]],
             [[1, 3, 4], [2]],
             [[1, 2, 4], [3]],
             [[1, 2, 3], [4]],
             [[1, 3], [2, 4]],
             [[1, 2], [3, 4]],
             [[1, 4], [2], [3]],
             [[1, 3], [2], [4]],
             [[1, 2], [3], [4]],
             [[1], [2], [3], [4]]]
            sage: ST4 = StandardTableaux(4)
            sage: ST4[0].parent() is ST4
            True
        """
        for p in sage.combinat.partition.Partitions(self.size):
            for st in StandardTableaux(p):
                yield self.element_class(self, st)


    def cardinality(self):
        """
        EXAMPLES::

            sage: StandardTableaux(3).cardinality()
            4
            sage: ns = [1,2,3,4,5,6]
            sage: sts = [StandardTableaux(n) for n in ns]    # indirect doctest
            sage: all([st.cardinality() == len(st.list()) for st in sts])
            True
        """
        c = 0
        for p in sage.combinat.partition.Partitions(self.size):
            c += StandardTableaux(p).cardinality()
        return c


class StandardTableaux_shape(StandardTableaux):
    def __init__(self, p):
        r"""
        Initializes the class of all semistandard tableaux of a given shape.
        Input is not checked; please use :class:`SemistandardTableaux` to
        ensure the options are properly parsed.

        TESTS::

            sage: ST = StandardTableaux([2,1,1])
            sage: TestSuite(ST).run()
        """

        super(StandardTableaux_shape, self).__init__(
              category = FiniteEnumeratedSets())
        self.shape = sage.combinat.partition.Partition(p)



    def __contains__(self, x):
        """
        EXAMPLES::

            sage: ST = StandardTableaux([2,1,1])
            sage: all([st in ST for st in ST])
            True
            sage: len(filter(lambda x: x in ST, StandardTableaux(4)))
            3
            sage: ST.cardinality()
            3
        """
        if isinstance(x, self.element_class) and map(len,x) == self.shape:
            return True

        try:
            x = StandardTableau(x)
        except ValueError:
            return False

        return map(len,x) == self.shape


    def _repr_(self):
        """
        TESTS::

            sage: repr(StandardTableaux([2,1,1]))    # indirect doctest
            'Standard tableaux of shape [2, 1, 1]'
        """
        return "Standard tableaux of shape %s"%str(self.shape)

    def cardinality(self):
        r"""
        Returns the number of standard Young tableaux associated with a
        partition pi

        A formula for the number of Young tableaux associated with a given
        partition. In each cell, write the sum of one plus the number of
        cells horizontally to the right and vertically below the cell (the
        hook length). The number of tableaux is then n! divided by the
        product of all hook lengths.

        For example, consider the partition [3,2,1] of 6 with Ferrers
        Diagram::

            # # #
            # #
            #

        When we fill in the cells with the hook
        lengths, we obtain::

            5 3 1
            3 1
            1

        The hook length formula returns 6!/(5\*3\*1\*3\*1\*1) = 16.

        EXAMPLES::

            sage: StandardTableaux([3,2,1]).cardinality()
            16
            sage: StandardTableaux([2,2]).cardinality()
            2
            sage: StandardTableaux([5]).cardinality()
            1
            sage: StandardTableaux([6,5,5,3]).cardinality()
            6651216

        REFERENCES:

        - http://mathworld.wolfram.com/HookLengthFormula.html
        """
        pi = self.shape

        number = factorial(sum(pi))
        hook = pi.hook_lengths()

        for row in range(len(pi)):
            for col in range(pi[row]):
                #Divide the hook length by the entry
                number /= hook[row][col]

        return Integer(number)

    def __iter__(self):
        r"""
        An iterator for the standard Young tableaux associated to the
        partition pi.

        EXAMPLES::

            sage: [t for t in StandardTableaux([2,2])]
            [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
            sage: [t for t in StandardTableaux([3,2])]
            [[[1, 3, 5], [2, 4]],
             [[1, 2, 5], [3, 4]],
             [[1, 3, 4], [2, 5]],
             [[1, 2, 4], [3, 5]],
             [[1, 2, 3], [4, 5]]]
            sage: st = StandardTableaux([2,1])
            sage: st[0].parent() is st
            True
        """

        pi = self.shape
        #Set the initial tableaux by filling it in going down the columns
        tableau = [[None]*n for n in pi]
        size = sum(pi)
        row = 0
        col = 0
        for i in range(size):
            tableau[row][col] = i+1

            #If we can move down, then do it;
            #otherwise, move to the next column over
            if ( row + 1 < len(pi) and col < pi[row+1]):
                row += 1
            else:
                row = 0
                col += 1

        yield self.element_class(self, tableau)

        if self.cardinality() == 1:
            last_tableau = True
        else:
            last_tableau = False

        while not last_tableau:
            #Convert the tableau to "vector format"
            #tableau_vector[i] is the row that number i
            #is in
            tableau_vector = [None]*size
            for row in range(len(pi)):
                for col in range(pi[row]):
                    tableau_vector[tableau[row][col]-1] = row

            #Locate the smallest integer j such that j is not
            #in the lowest corner of the subtableau T_j formed by
            #1,...,j.  This happens to be first j such that
            #tableau_vector[j]<tableau_vector[j-1].
            #l will correspond to the shape of T_j
            l = [0]*size
            l[0] = 1
            j = 0
            for i in range(1,size):
                l[tableau_vector[i]] += 1
                if ( tableau_vector[i] < tableau_vector[i-1] ):
                    j = i
                    break

            #Find the last nonzero row of l and store it in k
            i = size - 1
            while ( l[i] == 0 ):
                i -= 1
            k = i

            #Find a new row for the letter j (next lowest corner)
            t = l[ 1 + tableau_vector[j] ]
            i = k
            while ( l[i] != t ):
                i -= 1

            #Move the letter j to row i
            tableau_vector[j] = i
            l[i] -= 1

            #Fill in the columns of T_j using 1,...,j-1 in increasing order
            m = 0
            while ( m < j ):
                r = 0
                while ( l[r] != 0 ):
                    tableau_vector[m] = r
                    l[r] -= 1
                    m += 1
                    r += 1

            #Convert the tableau vector back to the regular tableau
            #format
            row_count= [0]*len(pi)
            tableau = [[None]*n for n in pi]

            for i in range(size):
                tableau[tableau_vector[i]][row_count[tableau_vector[i]]] = i+1
                row_count[tableau_vector[i]] += 1

            yield self.element_class(self, tableau)

            #Check to see if we are at the last tableau
            #The last tableau if given by filling in the
            #partition along the rows.  For example, the
            #last partition corresponding to [3,2] is
            #[[1,2,3],
            # [4,5]]
            last_tableau = True
            i = 1
            for row in range(len(pi)):
                for col in range(pi[row]):
                    if tableau[row][col] != i:
                        last_tableau = False
                    i += 1


        return


    def list(self):
        r"""
        Returns a list of the standard Young tableaux of the specified shape.

        EXAMPLES::

            sage: StandardTableaux([2,2]).list()
            [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
            sage: StandardTableaux([5]).list()
            [[[1, 2, 3, 4, 5]]]
            sage: StandardTableaux([3,2,1]).list()
            [[[1, 4, 6], [2, 5], [3]],
             [[1, 3, 6], [2, 5], [4]],
             [[1, 2, 6], [3, 5], [4]],
             [[1, 3, 6], [2, 4], [5]],
             [[1, 2, 6], [3, 4], [5]],
             [[1, 4, 5], [2, 6], [3]],
             [[1, 3, 5], [2, 6], [4]],
             [[1, 2, 5], [3, 6], [4]],
             [[1, 3, 4], [2, 6], [5]],
             [[1, 2, 4], [3, 6], [5]],
             [[1, 2, 3], [4, 6], [5]],
             [[1, 3, 5], [2, 4], [6]],
             [[1, 2, 5], [3, 4], [6]],
             [[1, 3, 4], [2, 5], [6]],
             [[1, 2, 4], [3, 5], [6]],
             [[1, 2, 3], [4, 5], [6]]]
        """
        return [y for y in self]


    def random_element(self):
        """
        Returns a random standard tableau of the given shape using the
        Green-Nijenhuis-Wilf Algorithm.

        EXAMPLES::

            sage: StandardTableaux([2,2]).random_element()
            [[1, 2], [3, 4]]
        """

        p = self.shape

        t = [[None]*n for n in p]


        #Get the cells in the
        cells = []
        for i in range(len(p)):
            for j in range(p[i]):
                cells.append((i,j))

        m = sum(p)
        while m > 0:

            #Choose a cell at random
            cell = random.choice(cells)


            #Find a corner
            inner_corners = p.corners()
            while cell not in inner_corners:
                hooks = []
                for k in range(cell[1], p[cell[0]]):
                    hooks.append((cell[0], k))
                for k in range(cell[0], len(p)):
                    if p[k] > cell[1]:
                        hooks.append((k, cell[1]))

                cell = random.choice(hooks)


            #Assign m to cell
            t[cell[0]][cell[1]] = m

            p = p.remove_cell(cell[0])

            cells.remove(cell)

            m -= 1

        return self.element_class(self, t)


##########################
# Symmetric group action #
##########################
def unmatched_places(w, open, close):
    """
    EXAMPLES::

        sage: from sage.combinat.tableau import unmatched_places
        sage: unmatched_places([2,2,2,1,1,1],2,1)
        ([], [])
        sage: unmatched_places([1,1,1,2,2,2],2,1)
        ([0, 1, 2], [3, 4, 5])
        sage: unmatched_places([], 2, 1)
        ([], [])
        sage: unmatched_places([1,2,4,6,2,1,5,3],2,1)
        ([0], [1])
        sage: unmatched_places([2,2,1,2,4,6,2,1,5,3], 2, 1)
        ([], [0, 3])
        sage: unmatched_places([3,1,1,1,2,1,2], 2, 1)
        ([1, 2, 3], [6])
    """
    lw = len(w)
    places_open = []
    places_close = []
    for i in range(lw):
        letter = w[i]
        if letter == open:
            places_open.append(i)
        elif letter == close:
            if places_open == []:
                places_close.append(i)
            else:
                places_open.pop()
    return places_close, places_open


def symmetric_group_action_on_values(word, perm):
    """
    EXAMPLES::

        sage: from sage.combinat.tableau import symmetric_group_action_on_values
        sage: symmetric_group_action_on_values([1,1,1],[1,3,2])
        [1, 1, 1]
        sage: symmetric_group_action_on_values([1,1,1],[2,1,3])
        [2, 2, 2]
        sage: symmetric_group_action_on_values([1,2,1],[2,1,3])
        [2, 2, 1]
        sage: symmetric_group_action_on_values([2,2,2],[2,1,3])
        [1, 1, 1]
        sage: symmetric_group_action_on_values([2,1,2],[2,1,3])
        [2, 1, 1]
        sage: symmetric_group_action_on_values([2,2,3,1,1,2,2,3],[1,3,2])
        [2, 3, 3, 1, 1, 2, 3, 3]
        sage: symmetric_group_action_on_values([2,1,1],[2,1])
        [2, 1, 2]
        sage: symmetric_group_action_on_values([2,2,1],[2,1])
        [1, 2, 1]
        sage: symmetric_group_action_on_values([1,2,1],[2,1])
        [2, 2, 1]
    """
    w = list(word)
    ts = sage.combinat.permutation.Permutation(perm).reduced_word()
    for j in reversed(range(len(ts))):
        r = ts[j]
        l = r + 1
        places_r, places_l = unmatched_places(w, l, r)

        #Now change the number of l's and r's in the new word
        nbl = len(places_l)
        nbr = len(places_r)
        ma = max(nbl, nbr)
        dif = ma - min(nbl, nbr)
        if ma == nbl:
            for i in range(dif):
                w[places_l[i]] = r
        else:
            for i in range(nbr-dif,ma):
                w[places_r[i]] = l
    return w


# August 2012: Deprecation of internal classes seems to be unnecessarily painful...
from sage.misc.superseded import deprecation

def Tableau_class(*args, **kargs):
    """
    EXAMPLES::

        sage: sage.combinat.tableau.Tableau_class([[3,2]])
        doctest:1: DeprecationWarning: this class is deprecated. Use Tableau_class instead
        See http://trac.sagemath.org/9265 for details.
        [[3, 2]]
    """
    deprecation(9265,'this class is deprecated. Use Tableau_class instead')
    return Tableau(*args, **kargs)

def Tableaux_n(*args, **kargs):
    """
    EXAMPLES::

        sage: sage.combinat.tableau.Tableaux_n(3)
        doctest:1: DeprecationWarning: this class is deprecated. Use Tableaux_size instead
        See http://trac.sagemath.org/9265 for details.
        Tableaux of size 3
    """
    deprecation(9265,'this class is deprecated. Use Tableaux_size instead')
    return Tableaux(*args, **kargs)

def SemistandardTableaux_n(*args, **kargs):
    """
    EXAMPLES::

        sage: sage.combinat.tableau.SemistandardTableaux_n(3)
        doctest:1: DeprecationWarning: this class is deprecated. Use SemistandardTableaux_size instead
        See http://trac.sagemath.org/9265 for details.
        Semistandard tableaux of size 3 and maximum entry 3
    """
    deprecation(9265,'this class is deprecated. Use SemistandardTableaux_size instead')
    return SemistandardTableaux(*args, **kargs)

def SemistandardTableaux_nmu(*args, **kargs):
    """
    EXAMPLES::

        sage: sage.combinat.tableau.SemistandardTableaux_nmu(3,[2,1])
        doctest:1: DeprecationWarning: this class is deprecated. Use SemistandardTableaux_size_weight instead
        See http://trac.sagemath.org/9265 for details.
        Semistandard tableaux of size 3 and weight [2, 1]
    """
    deprecation(9265,'this class is deprecated. Use SemistandardTableaux_size_weight instead')
    return SemistandardTableaux(*args, **kargs)

def SemistandardTableaux_p(*args, **kargs):
    """
    EXAMPLES::

        sage: sage.combinat.tableau.SemistandardTableaux_p([2,1])
        doctest:1: DeprecationWarning: this class is deprecated. Use SemistandardTableaux_shape instead
        See http://trac.sagemath.org/9265 for details.
        Semistandard tableaux of shape [2, 1] and maximum entry 3
    """
    deprecation(9265,'this class is deprecated. Use SemistandardTableaux_shape instead')
    return SemistandardTableaux(*args, **kargs)

def SemistandardTableaux_pmu(*args, **kargs):
    """
    EXAMPLES::

        sage: sage.combinat.tableau.SemistandardTableaux_pmu([2,1],[2,1])
        doctest:1: DeprecationWarning: this class is deprecated. Use SemistandardTableaux_shape_weight instead
        See http://trac.sagemath.org/9265 for details.
        Semistandard tableaux of shape [2, 1] and weight [2, 1]
    """
    deprecation(9265,'this class is deprecated. Use SemistandardTableaux_shape_weight instead')
    return SemistandardTableaux(*args, **kargs)

def StandardTableaux_n(*args, **kargs):
    """
    EXAMPLES::

        sage: sage.combinat.tableau.StandardTableaux_n(2)
        doctest:1: DeprecationWarning: this class is deprecated. Use StandardTableaux_size instead
        See http://trac.sagemath.org/9265 for details.
        Standard tableaux of size 2
    """
    deprecation(9265,'this class is deprecated. Use StandardTableaux_size instead')
    return StandardTableaux(*args, **kargs)

def StandardTableaux_partition(*args, **kargs):
    """
    EXAMPLES::

        sage: sage.combinat.tableau.StandardTableaux_partition([2,1])
        doctest:1: DeprecationWarning: this class is deprecated. Use StandardTableaux_shape instead
        See http://trac.sagemath.org/9265 for details.
        Standard tableaux of shape [2, 1]
    """
    deprecation(9265,'this class is deprecated. Use StandardTableaux_shape instead')
    return StandardTableaux(*args, **kargs)

# October 2012: fixing outdated pickles which use classed being deprecated
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.combinat.tableau', 'Tableau_class',  Tableau)
register_unpickle_override('sage.combinat.tableau', 'Tableaux_n',  Tableaux_size)
register_unpickle_override('sage.combinat.tableau', 'StandardTableaux_n',  StandardTableaux_size)
register_unpickle_override('sage.combinat.tableau', 'StandardTableaux_partition',  StandardTableaux_shape)
register_unpickle_override('sage.combinat.tableau', 'SemistandardTableaux_n',  SemistandardTableaux_size)
register_unpickle_override('sage.combinat.tableau', 'SemistandardTableaux_p',  SemistandardTableaux_shape)
register_unpickle_override('sage.combinat.tableau', 'SemistandardTableaux_nmu',  SemistandardTableaux_size_weight)
register_unpickle_override('sage.combinat.tableau', 'SemistandardTableaux_pmu',  SemistandardTableaux_shape_weight)


