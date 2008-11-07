"""
Tableaux
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
from sage.rings.arith import factorial
from sage.rings.integer import Integer
import sage.combinat.skew_tableau
import partition
from integer_vector import IntegerVectors
import word
import sage.libs.symmetrica.all as symmetrica
import sage.misc.prandom as random
import copy
import permutation
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.misc.misc import uniq
from combinat import CombinatorialClass, CombinatorialObject
import __builtin__

def Tableau(t):
    """
    Returns the tableau object corresponding to t.

    Note that Sage uses the English convention for
    partitions and tableaux.

    EXAMPLES:
        sage: t = Tableau([[1,2,3],[4,5]]); t
        [[1, 2, 3], [4, 5]]
        sage: t.shape()
        [3, 2]
        sage: t.is_standard()
        True
    """
    if isinstance(t, Tableau_class):
        return t
    elif t in Tableaux_all():
        return Tableau_class(t)
    raise ValueError, "invalid tableau"

class Tableau_class(CombinatorialObject):
    def __init__(self, t):
        """
        TESTS:
            sage: t = Tableau([[1,2],[3,4]])
            sage: t == loads(dumps(t))
            True
        """
        CombinatorialObject.__init__(self,t)

    def _latex_(self):
        r"""
        Returns a LaTeX version of self.

        EXAMPLES:
            sage: latex(Tableau([[1,2],[3,4]]))
            {\def\lr#1#2#3{\multicolumn{1}{#1@{\hspace{.6ex}}c@{\hspace{.6ex}}#2}{\raisebox{-.3ex}{$#3$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{cc}
            \cline{1-1}\cline{2-2}%
            \lr{|}{|}{1}&\lr{|}{|}{2}\\ %
            \cline{1-1}\cline{2-2}%
            \lr{|}{|}{3}&\lr{|}{|}{4}\\ %
            \cline{1-1}\cline{2-2}%
            \end{array}$}
            }
        """
        return self._tex_from_array()

    def _tex_from_array(self):
        r"""
        EXAMPLES:
            sage: print Tableau([[1,2],[3,4]])._tex_from_array()
            {\def\lr#1#2#3{\multicolumn{1}{#1@{\hspace{.6ex}}c@{\hspace{.6ex}}#2}{\raisebox{-.3ex}{$#3$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{cc}
            \cline{1-1}\cline{2-2}%
            \lr{|}{|}{1}&\lr{|}{|}{2}\\ %
            \cline{1-1}\cline{2-2}%
            \lr{|}{|}{3}&\lr{|}{|}{4}\\ %
            \cline{1-1}\cline{2-2}%
            \end{array}$}
            }
        """
        import output
        m = max(len(self), len(self[0]))
        return output.tex_from_array(self)

    def __div__(self, t):
        """
        Returns the skew partition self/t.

        EXAMPLES:
            sage: t = Tableau([[1,2,3],[3,4],[5]])
            sage: t/[1,1]
            [[None, 2, 3], [None, 4], [5]]
            sage: t/[3,1]
            [[None, None, None], [None, 4], [5]]
        """

        #if t is a list, convert to to a partition first
        if isinstance(t, list):
            t = partition.Partition(t)

        #Check to make sure that
        if not self.shape().dominates(t):
            raise ValueError, "the partition must dominate t"


        st = copy.deepcopy(self._list)

        for i in range(len(t)):
            for j in range(t[i]):
                st[i][j] = None

        return sage.combinat.skew_tableau.SkewTableau(st)


    def shape(self):
        r"""
        Returns the shape of a tableau t.

        EXAMPLES:
            sage: Tableau([[1,2,3],[4,5],[6]]).shape()
            [3, 2, 1]
        """

        return partition.Partition([len(row) for row in self])

    def size(self):
        """
        Returns the size of the shape of the tableau t.

        EXAMPLES:
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).size()
            6
            sage: Tableau([[1, 3], [2, 4]]).size()
            4

        """
        return sum([len(row) for row in self])

    def corners(self):
        """
        Returns the corners of the tableau t.

        EXAMPLES:
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).corners()
            [[0, 2], [1, 1], [2, 0]]
            sage: Tableau([[1, 3], [2, 4]]).corners()
            [[1, 1]]

        """
        return self.shape().corners()

    def conjugate(self):
        """
        Returns the conjugate of the tableau t.

        EXAMPLES:
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

        EXAMPLES:
            sage: Tableau([[1,2,3],[3,4],[5]]).pp()
              1  2  3
              3  4
              5
        """
        print '\n'.join([ "".join(map(lambda x: "%3s"%str(x) , row))  for row in self])

    def to_word_by_row(self):
        """
        Returns a word obtained from a row reading of the tableau t.

        EXAMPLES:
            sage: Tableau([[1,2],[3,4]]).to_word_by_row()
            [3, 4, 1, 2]
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word_by_row()
            [3, 2, 5, 1, 4, 6]
        """
        w = []
        for row in reversed(self):
            w += row
        return w

    def to_word_by_column(self):
        """
        Returns the word obtained from a column reading of the tableau t.

        EXAMPLES:
            sage: Tableau([[1,2],[3,4]]).to_word_by_column()
            [3, 1, 4, 2]
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word_by_column()
            [3, 2, 1, 5, 4, 6]
        """
        w = []
        conj = self.conjugate()
        for row in conj:
            w += list(reversed(row))
        return w

    def to_word(self):
        """
        An alias for to_word_by_row.

        EXAMPLES:
            sage: Tableau([[1,2],[3,4]]).to_word()
            [3, 4, 1, 2]
            sage: Tableau([[1, 4, 6], [2, 5], [3]]).to_word()
            [3, 2, 5, 1, 4, 6]
        """
        return self.to_word_by_row()


    def to_permutation(self):
        """
        Returns a permutation with the entries of self obtained by
        reading self in the reading order.

        EXAMPLES:
            sage: Tableau([[1,2],[3,4]]).to_permutation()
            [3, 4, 1, 2]
        """
        return permutation.Permutation(self.to_word())

    def descents(self):
        """
        Returns a list of the boxes (i,j) such that
        self[i][j] > self[i-1][j].

        EXAMPLES:
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
        Returns the major index of self.  The major index
        is defined to be the sum of the number of descents
        of self and the sum of their legs.

        EXAMPLES
            sage: Tableau( [[1,4],[2,3]] ).major_index()
            1
            sage: Tableau( [[1,2],[3,4]] ).major_index()
            2
        """
        descents = self.descents()
        p = self.shape()
        return len(descents) + sum([ p.leg(*d) for d in descents])

    def attacking_pairs(self):
        """
        Returns a list of the attacking pairs of self.  An pair of boxes
        (c, d) is said to be attacking if one of the following
        conditions hold:

            1) c and d lie in the same row with c to the west of d
            2) c is in the row immediately to the south of d and c
               lies strictly east of d.

        EXAMPLES:
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
        Returns a list of the inversions of self.  An inversion is an
        attacking pair (c,d) such that the entry of c in self is greater
        than the entry of d.

        EXAMPLES:
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

        The inversion number is defined to be the number of inversion of self
        minus the sum of the arm lengths of the descents of self.

        EXAMPLES:
            sage: t = Tableau([[1,2,3],[2,5]])
            sage: t.inversion_number()
            0
        """
        p = self.shape()
        return len(self.inversions()) - sum([ p.arm(*box) for box in self.descents() ])

    def entry(self, box):
        """
        Returns the entry of box in self.  Box is a tuple (i,j)
        of coordinates.

        EXAMPLES:
            sage: t = Tableau([[1,2],[3,4]])
            sage: t.entry( (0,0) )
            1
            sage: t.entry( (1,1) )
            4
        """
        i,j = box
        return self[i][j]

    def evaluation(self):
        """
        Returns the evaluation of the word from tableau t.

        EXAMPLES:
            sage: Tableau([[1,2],[3,4]]).evaluation()
            [1, 1, 1, 1]
        """

        return word.evaluation(self.to_word())

    weight = evaluation

    def is_standard(self):
        """
        Returns True if t is a standard tableau and False otherwise.

        EXAMPLES:
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

        EXAMPLES:
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

        EXAMPLES:
            sage: Tableau([[1,2],[3,4]]).vertical_flip()
            [[3, 4], [1, 2]]
        """

        if not self.is_rectangular():
            raise TypeError, "the tableau must be rectangular to use verticl_flip()"

        return Tableau([row for row in reversed(self)])

    def rotate_180(self):
        """
        Returns the tableau obtained by rotating t by 180 degrees.

        EXAMPLES:
            sage: Tableau([[1,2],[3,4]]).rotate_180()
            [[4, 3], [2, 1]]
        """
        if not self.is_rectangular():
            raise TypeError, "the tableau must be rectangular to use verticl_flip()"

        return Tableau([ [l for l in reversed(row)] for row in reversed(self) ])

    def boxes(self):
        """
        Returns a list of the coordinates of the boxes of self.

        EXAMPLES:
            sage: Tableau([[1,2],[3,4]]).boxes()
            [(0, 0), (0, 1), (1, 0), (1, 1)]
        """
        s = []
        for i in range(len(self)):
            s += [ (i,j) for j in range(len(self[i])) ]
        return s

    def k_weight(self, k):
        """
        Returns the k-weight of self.

        EXAMPLES:
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
        s = self.boxes()

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

        EXAMPLES:
            sage: Tableau([[1,2],[3],[4]]).restrict(3)
            [[1, 2], [3]]
            sage: Tableau([[1,2],[3],[4]]).restrict(2)
            [[1, 2]]
            sage: Tableau([[1,1],[2]]).restrict(1)
            [[1, 1]]
        """
        res = [ [y for y in row if y <=n] for row in self]
        return Tableau_class([row for row in res if row != []])

    def to_chain(self):
        """
        Returns the chain of partitions corresponding to the
        (semi)standard tableau.

        EXAMPLES:
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
        Returns the skew tableau formed by removing all of the boxes
        from self that are filled with a number less than

        EXAMPLES:
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
        An iterator for all the tableaux that can be obtained from self by adding a box.
        EXAMPLES:
            sage: t = Tableau([[1,2]])
            sage: [x for x in t.up()]
            [[[1, 2, 3]], [[1, 2], [3]]]
        """
        #Get a list of all places where we can add a box
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
        Returns a list of all the tableaux that can be obtained from self by adding a box.

        EXAMPLES:
            sage: t = Tableau([[1,2]])
            sage: t.up_list()
            [[[1, 2, 3]], [[1, 2], [3]]]
        """
        return list(self.up())

    def down(self):
        """
        An iterator for all the tableaux that can be obtained from self by removing a box.  Note that this iterates just over a single tableaux.
        EXAMPLES:
            sage: t = Tableau([[1,2],[3]])
            sage: [x for x in t.down()]
            [[[1, 2]]]
        """
        yield self.restrict( self.size() - 1 )

    def down_list(self):
        """
        Returns a list of all the tableaux that can be obtained from self by removing a box.  Note that this is just a single tableaux.

        EXAMPLES:
            sage: t = Tableau([[1,2],[3]])
            sage: t.down_list()
            [[[1, 2]]]
        """
        return list(self.down())

    def to_list(self):
        """
        EXAMPLES:
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

        EXAMPLES:
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
        EXAMPLES:
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
        EXAMPLES:
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
            return Tableau_class([[letter]])
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

        return Tableau_class([ row for row in rep if row != []])

    def _left_schensted_insert(self, letter):
        """
        EXAMPLES:
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
            return Tableau_class([[letter]])
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
        return Tableau_class(rep)


    def insert_word(self, w, left=False):
        """
        EXAMPLES:
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

        EXAMPLES:
            sage: t = Tableau([[1,2,2,3],[2,3,5,5],[4,4,6],[5,6]])
            sage: t2 = Tableau([[1,2],[3]])
            sage: t.bump_multiply(t2)
            [[1, 1, 2, 2, 3], [2, 2, 3, 5], [3, 4, 5], [4, 6, 6], [5]]

        """
        if not isinstance(right, Tableau_class):
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

        EXAMPLES:
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

    def promotion_inverse(self, n):
	"""
	Inverse promotion operator defined on rectangular tableaux using jeu de taquin

	EXAMPLES:
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
	"""
	Promotion operator defined on rectangular tableaux using jeu de taquin

	EXAMPLES:
	    sage: t = Tableau([[1,2],[3,3]])
	    sage: t.promotion(2)
	    [[1, 1], [2, 3]]
	    sage: t = Tableau([[1,1,1],[2,2,3],[3,4,4]])
	    sage: t.promotion(3)
	    [[1, 1, 2], [2, 2, 3], [3, 4, 4]]
	    sage: t = Tableau([[1,2],[2]])
	    sage: t.promotion(3)
	    Traceback (most recent call last):
	    ...
	    ValueError: Tableau is not rectangular
	"""
	if not self.is_rectangular():
	    raise ValueError, "Tableau is not rectangular"
	t = self.rotate_180()
	t = [[n+2-i for i in row] for row in t.to_list()]
	t = Tableau(t).promotion_inverse(n)
	t = [[n+2-i for i in row] for row in t.to_list()]
	return Tableau(t).rotate_180()


    def row_stabilizer(self):
        """
        Return the PermutationGroup corresponding to the row stabilizer
        of self.

        EXAMPLES:
            sage: rs = Tableau([[1,2,3],[4,5]]).row_stabilizer()
            sage: rs.order() == factorial(3)*factorial(2)
            True
            sage: PermutationGroupElement([(1,3,2),(4,5)]) in rs
            True
            sage: PermutationGroupElement([(1,4)]) in rs
            False
            sage: rs = Tableau([[1],[2],[3]]).row_stabilizer()
            sage: rs.order()
            1
        """

        gens = [ "()" ]
        for i in range(len(self)):
            for j in range(0, len(self[i])-1):
                gens.append( (self[i][j], self[i][j+1]) )
        return PermutationGroup( gens )

    def column_stabilizer(self):
        """
        Return the PermutationGroup corresponding to the column stabilizer
        of self.

        EXAMPLES:
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

        EXAMPLES:
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

        EXAMPLES:
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
        Returns True if self is less than or equal to tab2 in the
        last letter ordering.

        EXAMPLES:
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
        if not isinstance(tab2, Tableau_class):
            try:
                tab2 = Tableau_class(tab2)
            except:
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
        """
        EXAPMLES:
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
        return word.charge([i for i in reversed(self.to_word())])


    def cocharge(self):
        """
        EXAPMLES:
            sage: Tableau([[1,1],[2,2],[3]]).cocharge()
            4
            sage: Tableau([[1,1,3],[2,2]]).cocharge()
            3
            sage: Tableau([[1,1,2],[2],[3]]).cocharge()
            2
            sage: Tableau([[1,1,2],[2,3]]).cocharge()
            2
            sage: Tableau([[1,1,2,3],[2]]).cocharge()
            1
            sage: Tableau([[1,1,2,2],[3]]).cocharge()
            1
            sage: Tableau([[1,1,2,2,3]]).cocharge()
            0

        """
        return word.charge(self.to_word())


    ##############
    # katabolism #
    ##############

    def katabolism(self):
        """
        EXAMPLES:
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
            return Tableau_class(self[1:]).insert_word(self[0],left=True)

    def katabolism_sequence(self):
        """
        EXAMPLES:
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
        """
        EXAMPLES:
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
            return Tableau_class([])

        w1 = Tableau_class( self[len(part):] ).to_word()

        if len(self)-len(part)+11 > 0:
            t2 = Tableau_class(self[:len(part)])
        else:
            t2 = self

        if t2.shape() == part:
            t2 = Tableau_class([])
        else:
            part += [0]*(len(t2)-len(part))
            t2 = [[t2[i][j] for j in range(part[i], len(t2[i]))]  for i in range(len(t2)) ]
            t2 = Tableau_class([ row for row in t2 if row != [] ])

        w2 = t2.to_word()
        return Tableau_class([]).insert_word(w2+w1)


    def reduced_lambda_katabolism(self, part):
        """
        EXAMPLES:
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
        tt_part = Tableau_class([ [a+i]*part[i] for i in range(len(part)) ])
        t_part = Tableau_class([[self[i][j] for j in range(part[i])] for i in range(len(part))])

        if t_part == tt_part:
            return res
        else:
            return 0

    def katabolism_projector(self, parts):
        """
        EXAMPLES:
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
            return Tableau_class([])


    def promotion_operator(self, i):
        """
        EXAMPLES:
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

        TESTS:
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
        EXAMPLES:
            sage: from sage.combinat.word import symmetric_group_action_on_values
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
        EXAMPLES:
            sage: t = Tableau([[1,1,3,3],[2,3],[3]])
            sage: t.symmetric_group_action_on_values([1,2,3])
            [[1, 1, 3, 3], [2, 3], [3]]
            sage: t.symmetric_group_action_on_values([3,2,1])
            [[1, 1, 1, 1], [2, 3], [3]]
            sage: t.symmetric_group_action_on_values([1,3,2])
            [[1, 1, 2, 2], [2, 2], [3]]
        """
        return self.raise_action_from_words(word.symmetric_group_action_on_values, perm)

    #########
    # atoms #
    #########
    def socle(self):
        """
        EXAMPLES:
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
        EXAMPLES:
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


def from_chain(chain):
    """
    Returns a semistandard tableau from a chain of partitions.

    EXAMPLES:
        sage: from sage.combinat.tableau import from_chain
        sage: from_chain([[], [2], [2, 1], [3, 2, 1]])
        [[1, 1, 3], [2, 3], [3]]
    """
    res = [[0]*chain[-1][i] for i in range(len(chain[-1]))]
    for i in reversed(range(2, len(chain)+1)):
        for j in range(len(chain[i-1])):
            for k in range(chain[i-1][j]):
                res[j][k] = i -1
    return Tableau_class(res)

def from_shape_and_word(shape, w):
    """
    Returns a tableau from a shape and word.

    EXAMPLES:
        sage: from sage.combinat.tableau import from_shape_and_word
        sage: t = Tableau([[1, 3], [2], [4]])
        sage: shape = t.shape(); shape
        [2, 1, 1]
        sage: word  = t.to_word(); word
        [4, 2, 1, 3]
        sage: from_shape_and_word(shape, word)
        [[1, 3], [2], [4]]
    """
    res = []
    j = 0
    for i in reversed(range(len(shape))):
        res.append( w[j:j+shape[i]] )
        j += shape[i]
    res.reverse()
    return Tableau_class(res)

def Tableaux(n=None):
    """
    Returns the combinatorial class of tableaux.  If n
    is specified, then it returns the combinatoiral class
    of all tableaux of size n.

    EXAMPLES:
        sage: T = Tableaux(); T
        Tableaux
        sage: [[1,2],[3,4]] in T
        True
        sage: [[1,2],[3]] in T
        True
        sage: [1,2,3] in T
        False

        sage: T = Tableaux(4); T
        Tableaux of size 4
        sage: [[1,2],[3,4]] in T
        True
        sage: [[1,2],[3]] in T
        False
        sage: [1,2,3] in T
        False
    """
    if n == None:
        return Tableaux_all()
    else:
        return Tableaux_n(n)

class Tableaux_all(CombinatorialClass):
    def __init__(self):
        """
        TESTS:
            sage: T = Tableaux()
            sage: T == loads(dumps(T))
            True
        """
        pass

    def __contains__(self, x):
        """
        TESTS:
            sage: T = Tableaux()
            sage: [[1,2],[3,4]] in T
            True
            sage: [[1,2],[3]] in T
            True
            sage: [1,2,3] in T
            False
        """
        if isinstance(x, Tableau_class):
            return True

        if not isinstance(x, __builtin__.list):
            return False

        for row in x:
            if not isinstance(row, __builtin__.list):
                return False

        if map(len, x) not in partition.Partitions_all():
            return False

        return True

    def __repr__(self):
        """
        TESTS:
            sage: repr(Tableaux())
            'Tableaux'
        """
        return "Tableaux"

    def list(self):
        """
        TESTS:
            sage: Tableaux().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def iterator(self):
        """
        TESTS:
            sage: Tableaux().iterator()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class Tableaux_n(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: T = Tableaux(3)
            sage: T == loads(dumps(T))
            True
        """
        self.n = n


    def __repr__(self):
        """
        TESTS:
            sage: repr(Tableaux(4))
            'Tableaux of size 4'
        """
        return "Tableaux of size %s"%self.n

    def __contains__(self,x):
        """
        EXAMPLES:
            sage: [[2,4],[1,3]] in Tableaux(3)
            False
            sage: [[2,4], [1]] in Tableaux(3)
            True
        """
        return x in Tableaux() and sum(map(len, x)) == self.n

    def list(self):
        """
        TESTS:
            sage: Tableaux(3).list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def iterator(self):
        """
        TESTS:
            sage: Tableaux(3).iterator()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

def StandardTableaux(n=None):
    """
    Returns the combinatorial class of standard tableaux.
    If n is specified abd is an integer, then it returns
    the combinatorial class of all standard tableaux of
    size n.  If n is a partition, then it returns the class
    of all standard tableaux of shape n.

    EXAMPLES:
        sage: ST = StandardTableaux(3); ST
        Standard tableaux of size 3
        sage: ST.first()
        [[1, 2, 3]]
        sage: ST.last()
        [[1], [2], [3]]
        sage: ST.count()
        4
        sage: ST.list()
        [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]

        sage: ST = StandardTableaux([2,2]); ST
        Standard tableaux of shape [2, 2]
        sage: ST.first()
        [[1, 3], [2, 4]]
        sage: ST.last()
        [[1, 2], [3, 4]]
        sage: ST.count()
        2
        sage: ST.list()
        [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
    """
    if n == None:
        return StandardTableaux_all()
    elif n in partition.Partitions():
        return StandardTableaux_partition(n)
    else:
        return StandardTableaux_n(n)

class StandardTableaux_all(CombinatorialClass):
    def __init__(self):
        """
        TESTS:
            sage: ST = StandardTableaux()
            sage: ST == loads(dumps(ST))
            True
        """
        pass

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: [[1,1],[2,3]] in StandardTableaux()
            False
            sage: [[1,2],[3,4]] in StandardTableaux()
            True
            sage: [[1,3],[2,4]] in StandardTableaux()
            True
        """
        if x not in Tableaux():
            return False
        else:
            t = Tableau(x)

        #Check to make sure the first position is 1
        fillings = []
        for row in t:
            fillings += row
        fillings.sort()
        if fillings != range(1, max(fillings)+1):
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

    def __repr__(self):
        """
        TESTS:
            sage: repr(StandardTableaux())
            'Standard tableaux'
        """
        return "Standard tableaux"

    def list(self):
        """
        TESTS:
            sage: StandardTableaux().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class StandardTableaux_n(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: ST = StandardTableaux(3)
            sage: ST == loads(dumps(ST))
            True
        """
        self.n = n

    object_class = Tableau_class

    def __repr__(self):
        """
        TESTS:
            sage: repr(StandardTableaux(3))
            'Standard tableaux of size 3'
        """
        return "Standard tableaux of size %s"%self.n

    def __contains__(self, x):
        """
        TESTS:
            sage: ST3 = StandardTableaux(3)
            sage: all([st in ST3 for st in ST3])
            True
            sage: ST4 = StandardTableaux(4)
            sage: filter(lambda x: x in ST3, ST4)
            []
        """
        return x in StandardTableaux() and sum(map(len, x)) == self.n

    def iterator(self):
        """
        EXAMPLES:
            sage: StandardTableaux(1).list()
            [[[1]]]
            sage: StandardTableaux(2).list()
            [[[1, 2]], [[1], [2]]]
            sage: StandardTableaux(3).list()
            [[[1, 2, 3]], [[1, 3], [2]], [[1, 2], [3]], [[1], [2], [3]]]
            sage: StandardTableaux(4).list()
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
        """
        for p in partition.Partitions(self.n):
            for st in StandardTableaux(p):
                yield st

    def count(self):
        """
        EXAMPLES:
            sage: StandardTableaux(3).count()
            4
            sage: ns = [1,2,3,4,5,6]
            sage: sts = [StandardTableaux(n) for n in ns]
            sage: all([st.count() == len(st.list()) for st in sts])
            True
        """
        c = 0
        for p in partition.Partitions(self.n):
            c += StandardTableaux(p).count()
        return c


class StandardTableaux_partition(CombinatorialClass):
    def __init__(self, p):
        """
        TESTS:
            sage: ST = StandardTableaux([2,1,1])
            sage: ST == loads(dumps(ST))
            True
        """
        self.p = partition.Partition(p)

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: ST = StandardTableaux([2,1,1])
            sage: all([st in ST for st in ST])
            True
            sage: len(filter(lambda x: x in ST, StandardTableaux(4)))
            3
            sage: ST.count()
            3

        """
        return x in StandardTableaux() and map(len,x) == self.p

    def __repr__(self):
        """
        TESTS:
            sage: repr(StandardTableaux([2,1,1]))
            'Standard tableaux of shape [2, 1, 1]'
        """
        return "Standard tableaux of shape %s"%str(self.p)

    def count(self):
        r"""
        Returns the number of standard Young tableaux associated with
        a partition pi

        A formula for the number of Young tableaux associated with a given partition.
        In each box, write the sum of one plus the number of boxes horizontally to the right
        and vertically below the box (the hook length).
        The number of tableaux is then n! divided by the product of all hook lengths.

        For example, consider the partition [3,2,1] of 6 with Ferrers Diagram
        * * *
        * *
        *
        When we fill in the boxes with the hook lengths, we obtain
        5 3 1
        3 1
        1
        The hook length formula returns 6!/(5*3*1*3*1*1) = 16.

        EXAMPLES:
            sage: StandardTableaux([3,2,1]).count()
            16
            sage: StandardTableaux([2,2]).count()
            2
            sage: StandardTableaux([5]).count()
            1
            sage: StandardTableaux([6,5,5,3]).count()
            6651216

        REFERENCES:
            http://mathworld.wolfram.com/HookLengthFormula.html
        """
        pi = self.p

        number = factorial(sum(pi))
        hook = pi.hook_lengths()

        for row in range(len(pi)):
            for col in range(pi[row]):
                #Divide the hook length by the entry
                number /= hook[row][col]

        return number

    def iterator(self):
        r"""
        An iterator for the standard Young tableaux associated to the
        partition pi.

        EXAMPLES:
            sage: [p for p in StandardTableaux([2,2])]
            [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
            sage: [p for p in StandardTableaux([3,2])]
            [[[1, 3, 5], [2, 4]],
             [[1, 2, 5], [3, 4]],
             [[1, 3, 4], [2, 5]],
             [[1, 2, 4], [3, 5]],
             [[1, 2, 3], [4, 5]]]

        """

        pi = self.p
        #Set the intial tableaux by filling it in going down the columns
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

        yield Tableau(tableau)

        if self.count() == 1:
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

            yield Tableau(tableau)

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
        Returns a list of the standard Young tableau associated with a
        partition p.

        EXAMPLES:
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
        Returns a random standard tableau of shape p using the
        Green-Nijenhuis-Wilf Algorithm.


        EXAMPLES:
            sage: StandardTableaux([2,2]).random_element()
            [[1, 2], [3, 4]]
        """

        p = self.p

        t = [[None]*n for n in p]


        #Get the cells in the
        cells = []
        for i in range(len(p)):
            for j in range(p[i]):
                cells.append([i,j])

        m = sum(p)
        while m > 0:

            #Choose a cell at random
            cell = random.choice(cells)


            #Find a corner
            inner_corners = p.corners()
            while cell not in inner_corners:
                hooks = []
                for k in range(cell[1], p[cell[0]]):
                    hooks.append([cell[0], k])
                for k in range(cell[0], len(p)):
                    if p[k] > cell[1]:
                        hooks.append([k, cell[1]])

                cell = random.choice(hooks)


            #Assign m to cell
            t[cell[0]][cell[1]] = m

            p = p.remove_box(cell[0])

            cells.remove(cell)

            m -= 1

        return Tableau(t)



##########################
# Semi-standard tableaux #
##########################

def SemistandardTableaux(p=None, mu=None):
    """
    Returns the combinatorial class of semistandard tableaux.

    If p is specified and is a partition, then it returns the
    class of semistandard tableaux of shape p (and max entry
    sum(p))

    If p is specified and is an integer, it returns the class
    of semistandard tableaux of size p.

    If mu is also specified, then it returns the class of
    semistandard tableaux with evaluation/content mu.

    EXAMPLES:
        sage: SST = SemistandardTableaux([2,1]); SST
        Semistandard tableaux of shape [2, 1]
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
        Semistandard tableaux of size 3
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
    """
    if p == None:
        return SemistandardTableaux_all()
    elif p in partition.Partitions():
        if mu == None:
            return SemistandardTableaux_p(p)
        else:
            if sum(p) != sum(mu):
                #Error size mismatch
                raise TypeError, "p and mu must be of the same size"
            else:
                return SemistandardTableaux_pmu(p, mu)
    elif isinstance(p, (int, Integer)):
        if mu == None:
            return SemistandardTableaux_n(p)
        else:
            if p != sum(mu):
                #Error size mismatch
                raise TypeError, "mu must be of size p (= %s)"%p
            else:
                return SemistandardTableaux_nmu(p, mu)
    else:
        raise ValueError

class SemistandardTableaux_all(CombinatorialClass):
    def __init__(self):
        """
        TESTS:
            sage: SST = SemistandardTableaux()
            sage: SST == loads(dumps(SST))
            True
        """

    def __contains__(self, x):
        """
        TESTS:
            sage: [[1,2],[1]] in SemistandardTableaux()
            False
            sage: SST = SemistandardTableaux()
            sage: all([st in SST for st in StandardTableaux(4)])
            True
            sage: [[1,1],[2]] in SemistandardTableaux()
            True
        """
        if x not in Tableaux():
            return False
        else:
            t = Tableau(x)

        #Check to make sure the first position is 1
        for row in t:
            for i in row:
                if not isinstance(i, (int, Integer)):
                    return False

        #Check to make sure it is non-decreasing along the rows
        for row in t:
            for i in range(1, len(row)):
                if row[i] < row[i-1]:
                    return False

        #Check to make sure it is increasing along the columns
        conj = t.conjugate()
        for row in conj:
            for i in range(1, len(row)):
                if row[i] <= row[i-1]:
                    return False

        return True

    def list(self):
        """
        TESTS:
            sage: SemistandardTableaux().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

class SemistandardTableaux_n(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: SST = SemistandardTableaux(3)
            sage: SST == loads(dumps(SST))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS:
            sage: repr(SemistandardTableaux(3))
            'Semistandard tableaux of size 3'
        """
        return "Semistandard tableaux of size %s"%str(self.n)

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: [[1,2],[3,3]] in SemistandardTableaux(3)
            False
            sage: [[1,2],[3,3]] in SemistandardTableaux(4)
            True
            sage: SST = SemistandardTableaux(4)
            sage: all([sst in SST for sst in SST])
            True
        """
        return x in SemistandardTableaux() and sum(map(len, x)) == self.n

    object_class = Tableau_class

    def count(self):
        """
        EXAMPLES:
            sage: SemistandardTableaux(3).count()
            19
            sage: SemistandardTableaux(4).count()
            116
            sage: ns = range(1, 6)
            sage: ssts = [ SemistandardTableaux(n) for n in ns ]
            sage: all([sst.count() == len(sst.list()) for sst in ssts])
            True
        """
        c = 0
        for part in partition.Partitions(self.n):
            c += SemistandardTableaux(part).count()
        return c

    def iterator(self):
        """
        EXAMPLES:
            sage: SemistandardTableaux(2).list()
            [[[1, 1]], [[1, 2]], [[2, 2]], [[1], [2]]]
            sage: SemistandardTableaux(3).list()
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

        """
        for part in partition.Partitions(self.n):
            for sst in SemistandardTableaux(part):
                yield sst

class SemistandardTableaux_pmu(CombinatorialClass):
    def __init__(self, p, mu):
        """
        TESTS:
            sage: SST = SemistandardTableaux([2,1], [2,1])
            sage: SST == loads(dumps(SST))
            True
        """
        self.p = p
        self.mu = mu

    object_class = Tableau_class

    def __repr__(self):
        """
        TESTS:
            sage: repr(SemistandardTableaux([2,1],[2,1]))
            'Semistandard tableaux of shape [2, 1] and evaluation [2, 1]'
        """
        return "Semistandard tableaux of shape %s and evaluation %s"%(self.p, self.mu)

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: SST = SemistandardTableaux([2,1], [2,1])
            sage: all([sst in SST for sst in SST])
            True
            sage: len(filter(lambda x: x in SST, SemistandardTableaux(3)))
            1
            sage: SST.count()
            1
        """
        if not x in SemistandardTableaux(self.p):
            return False
        n = sum(self.p)

        if n == 0 and len(x) == 0:
            return True

        content = {}
        for row in x:
            for i in row:
                content[i] = content.get(i, 0) + 1
        content_list = [0]*int(max(content))

        for key in content:
            content_list[key-1] = content[key]

        if content_list != self.mu:
            return False

        return True


    def count(self):
        """
        EXAMPLES:
            sage: SemistandardTableaux([2,2], [2, 1, 1]).count()
            1
            sage: SemistandardTableaux([2,2,2], [2, 2, 1,1]).count()
            1
            sage: SemistandardTableaux([2,2,2], [2, 2, 2]).count()
            1
            sage: SemistandardTableaux([3,2,1], [2, 2, 2]).count()
            2
        """
        return symmetrica.kostka_number(self.p,self.mu)


    def list(self):
        """
        EXAMPLES:
            sage: SemistandardTableaux([2,2], [2, 1, 1]).list()
            [[[1, 1], [2, 3]]]
            sage: SemistandardTableaux([2,2,2], [2, 2, 1,1]).list()
            [[[1, 1], [2, 2], [3, 4]]]
            sage: SemistandardTableaux([2,2,2], [2, 2, 2]).list()
            [[[1, 1], [2, 2], [3, 3]]]
            sage: SemistandardTableaux([3,2,1], [2, 2, 2]).list()
            [[[1, 1, 2], [2, 3], [3]], [[1, 1, 3], [2, 2], [3]]]
        """
        return symmetrica.kostka_tab(self.p, self.mu)


class SemistandardTableaux_p(CombinatorialClass):
    def __init__(self, p):
        """
        TESTS:
            sage: SST = SemistandardTableaux([2,1])
            sage: SST == loads(dumps(SST))
            True
        """
        self.p = p

    object_class = Tableau_class

    def __contains__(self, x):
        """
        EXAMPLES:
            sage: SST = SemistandardTableaux([2,1])
            sage: all([sst in SST for sst in SST])
            True
            sage: len(filter(lambda x: x in SST, SemistandardTableaux(3)))
            8
            sage: SST.count()
            8
        """
        return x in SemistandardTableaux_all() and map(len, x) == self.p

    def __repr__(self):
        """
        TESTS:
            sage: repr(SemistandardTableaux([2,1]))
            'Semistandard tableaux of shape [2, 1]'
        """
        return "Semistandard tableaux of shape %s" % str(self.p)

    def count(self):
        """
        EXAMPLES:
            sage: SemistandardTableaux([2,1]).count()
            8
            sage: SemistandardTableaux([2,2,1]).count()
            75
            sage: s = SFASchur(QQ)
            sage: s([2,2,1]).expand(5)(1,1,1,1,1)
            75
            sage: SemistandardTableaux([5]).count()
            126
            sage: SemistandardTableaux([3,2,1]).count()
            896
        """
        c = 0
        for comp in IntegerVectors(sum(self.p), sum(self.p)):
            c += SemistandardTableaux_pmu(self.p, comp).count()
        return c

    def iterator(self):
        """
        An iterator for the semistandard partitions of shape p.

        EXAMPLES:
            sage: SemistandardTableaux([3]).list()
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
            sage: SemistandardTableaux([2,1]).list()
            [[[1, 1], [2]],
             [[1, 1], [3]],
             [[1, 2], [2]],
             [[1, 2], [3]],
             [[1, 3], [2]],
             [[1, 3], [3]],
             [[2, 2], [3]],
             [[2, 3], [3]]]
            sage: SemistandardTableaux([1,1,1]).list()
            [[[1], [2], [3]]]
        """
        for c in IntegerVectors(sum(self.p), sum(self.p)):
            for sst in SemistandardTableaux(self.p, c):
                yield sst

class SemistandardTableaux_nmu(CombinatorialClass):
    def __init__(self, n, mu):
        """
        TESTS:
            sage: SST = SemistandardTableaux(3, [2,1])
            sage: SST == loads(dumps(SST))
            True
        """
        self.n = n
        self.mu = mu

    def __repr__(self):
        """
        TESTS:
            sage: repr(SemistandardTableaux(3, [2,1]))
            'Semistandard tableaux of size 3 and evaluation [2, 1]'
        """
        return "Semistandard tableaux of size %s and evaluation %s"%(self.n, self.mu)

    def iterator(self):
        """
        EXAMPLES:
            sage: SemistandardTableaux(3, [2,1]).list()
            [[[1, 1, 2]], [[1, 1], [2]]]
            sage: SemistandardTableaux(4, [2,2]).list()
            [[[1, 1, 2, 2]], [[1, 1, 2], [2]], [[1, 1], [2, 2]]]
        """
        for p in partition.Partitions(self.n):
            for sst in SemistandardTableaux_pmu(p, self.mu):
                yield sst

    def count(self):
        """
        EXAMPLES:
            sage: SemistandardTableaux(3, [2,1]).count()
            2
            sage: SemistandardTableaux(4, [2,2]).count()
            3
        """
        c = 0
        for p in partition.Partitions(self.n):
            c += SemistandardTableaux_pmu(p, self.mu).count()
        return c

    def __contains__(self, x):
        """
        TESTS:
            sage: SST = SemistandardTableaux(6, [2,2,2])
            sage: all([sst in SST for sst in SST])
            True
            sage: all([sst in SST for sst in SemistandardTableaux([3,2,1],[2,2,2])])
            True
        """
        return x in SemistandardTableaux_all() and x in SemistandardTableaux(map(len, x), self.mu)
