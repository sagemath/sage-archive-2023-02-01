r"""
Skew Tableaux
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
#from sage.combinat.partition import Partition
import partition
import sage.combinat.tableau
from sage.combinat.skew_partition import SkewPartition
import partition
import misc
import word
from sage.misc.all import prod
import exceptions
import random
import copy
from combinat import CombinatorialObject, CombinatorialClass
from sage.graphs.graph import DiGraph


def SkewTableau(st=None, expr=None):
    """
    Returns the skew tableau object corresponding to st.

    EXAMPLES:
        sage: st = SkewTableau([[None, 1],[2,3]]); st
        [[None, 1], [2, 3]]
        sage: st.inner_shape()
        [1]
        sage: st.outer_shape()
        [2, 2]

        sage: SkewTableau(expr=[[1,1],[[5],[3,4],[1,2]]])
        [[None, 1, 2], [None, 3, 4], [5]]
    """
    if expr is not None:
        return from_expr(expr)

    for row in st:
        if not isinstance(row, list):
            raise TypeError, "each element of the skew tableau must be a list"
        if row == []:
            raise TypeError, "a skew tableau cannot have an empty list for a row"

    return SkewTableau_class(st)

class SkewTableau_class(CombinatorialObject):
    def __init__(self, t):
        """
        TESTS:
            sage: st = SkewTableau([[None, 1],[2,3]])
            sage: st == loads(dumps(st))
            True
        """
        CombinatorialObject.__init__(self,t)

    def pp(self):
        """
        Returns a pretty print string of the tableau.
        EXAMPLES:
            sage: t = SkewTableau([[None,2,3],[None,4],[5]])
            sage: print t.pp()
              .  2  3
              .  4
              5
        """

        def none_str(x):
            if x is None:
                return "  ."
            else:
                return "%3s"%str(x)

        new_rows = [ "".join(map(none_str , row)) for row in self]
        return '\n'.join(new_rows)

    def outer_shape(self):
        """
        Returns the outer shape of the tableau.

        EXAMPLES:
            sage: SkewTableau([[None,1,2],[None,3],[4]]).outer_shape()
            [3, 2, 1]
        """

        return partition.Partition([len(row) for row in self])


    def inner_shape(self):
        """
        Returns the inner shape of the tableau.

        EXAMPLES:
            sage: SkewTableau([[None,1,2],[None,3],[4]]).inner_shape()
            [1, 1]
        """

        return partition.Partition(filter(lambda x: x != 0, [len(filter(lambda x: x is None, row)) for row in self]))

    def shape(self):
        r"""
        Returns the shape of a tableau t.

        EXAMPLES:
            sage: SkewTableau([[None,1,2],[None,3],[4]]).shape()
            [[3, 2, 1], [1, 1]]
        """

        return SkewPartition([self.outer_shape(), self.inner_shape()])


    def outer_size(self):
        """
        Returns the size of the outer shape of the skew tableau.

        EXAMPLES:
            sage: SkewTableau([[None, 2, 4], [None, 3], [1]]).outer_size()
            6
            sage: SkewTableau([[None, 2], [1, 3]]).outer_size()
            4

        """
        return self.outer_shape().size()



    def inner_size(self):
        """
        Returns the size of the inner shape of the skew tableau.

        EXAMPLES:
            sage: SkewTableau([[None, 2, 4], [None, 3], [1]]).inner_size()
            2
            sage: SkewTableau([[None, 2], [1, 3]]).inner_size()
            1

        """
        return self.inner_shape().size()

    def size(self):
        """
        Returns the number of boxes in the skew tableau.

        EXAMPLES:
            sage: SkewTableau([[None, 2, 4], [None, 3], [1]]).size()
            4
            sage: SkewTableau([[None, 2], [1, 3]]).size()
            3

        """
        return sum([len(filter(lambda x: x is not None,row)) for row in self])



    def conjugate(self):
        """
        Returns the conjugate of the skew tableau.

        EXAMPLES:
            sage: SkewTableau([[None,1],[2,3]]).conjugate()
            [[None, 2], [1, 3]]
        """
        conj_shape = self.outer_shape().conjugate()

        conj = [[None]*row_length for row_length in conj_shape]

        for i in range(len(conj)):
            for j in range(len(conj[i])):
                conj[i][j] = self[j][i]


        return SkewTableau(conj)

    def to_expr_list(self):
        """
        The first list in a result corresponds to the inner partition of the
        skew shape.  The second list is a list of the rows in the skew
        tableau read from the bottom up.

        Provided for compatability with MuPAD-Combinat.  In MuPAD-Combinat, if
        t is a skew tableau, then to_expr_list gives the same result as
        expr(t) would give in MuPAD-Combinat.


        EXAMPLES:
            sage: st = SkewTableau([[None,1,2],[None,3,4],[5]])
            sage: st.to_expr_list()
            [[1, 1], [[5], [3, 4], [1, 2]]]
        """
        is_none = lambda x: x is None
        inner = []
        outer = []
        for row in self:
            inner.append( len( filter(is_none, row) ) )
            outer.append( row[ inner[-1] : ] )
        outer.reverse()
        inner = filter(lambda x: x != 0, inner)
        return [inner, outer]

    def to_word_by_row(self):
        """
        Returns a word obtained from a row reading of the skew tableau.

        EXAMPLES:
            sage: SkewTableau([[None,1],[2,3]]).to_word_by_row()
            [1, 2, 3]
            sage: SkewTableau([[None, 2, 4], [None, 3], [1]]).to_word_by_row()
            [2, 4, 3, 1]
        """
        word = []
        for row in self:
            word += filter(lambda x: x is not None, row)

        return word


    def to_word_by_column(self):
        """
        Returns the word obtained from a column reading of the skew tableau

        EXAMPLES:
            sage: SkewTableau([[None,1],[2,3]]).to_word_by_column()
            [2, 1, 3]
            sage: SkewTableau([[None, 2, 4], [None, 3], [1]]).to_word_by_column()
            [1, 2, 3, 4]

        """
        word = []
        conj = self.conjugate()
        for row in conj:
            word += filter(lambda x: x is not None, row)

        return word

    def to_word_by_reading_order(self):
        """

        EXAMPLES:
            sage: SkewTableau([[None,1],[2,3]]).to_word_by_reading_order()
            [2, 3, 1]
        """
        word = []
        for row in self:
            word = row + word
        return [x for x in word if x is not None]

    def to_word(self):
        """
        An alias for SkewTableau.to_word_by_row().
        """
        return self.to_word_by_row()

    def evaluation(self):
        """
        Returns the evaluation of the word from skew tableau.

        EXAMPLES:
            sage: SkewTableau([[1,2],[3,4]]).evaluation()
            [1, 1, 1, 1]
        """

        return word.evaluation(self.to_word())

    def is_standard(self):
        """
        Returns True if t is a standard skew tableau and False otherwise.

        EXAMPLES:
            sage: SkewTableau([[None, 2], [1, 3]]).is_standard()
            True
            sage: SkewTableau([[None, 2], [2, 4]]).is_standard()
            False
            sage: SkewTableau([[None, 3], [2, 4]]).is_standard()
            False
            sage: SkewTableau([[None, 2], [2, 4]]).is_standard()
            False
        """
        t = self
        #Check to make sure that it is filled with 1...size
        w = self.to_word()
        w.sort()
        if w != range(1, self.size()+1):
            return False



        #Check to make sure it is increasing along the rows
        for row in t:
            for i in range(1, len(row)):
                if row[i-1] is not None and row[i] <= row[i-1]:
                    return False


        #Check to make sure it is increasing along the columns
        conj = t.conjugate()
        for row in conj:
            for i in range(1, len(row)):
                if row[i-1] is not None and row[i] <= row[i-1]:
                    return False

        return True

    def to_tableau(self):
        """
        Returns a tableau with the same filling.  This only works if the
        inner shape of the skew tableau has size zero.

        EXAMPLES:
            sage: SkewTableau([[1,2],[3,4]]).to_tableau()
            [[1, 2], [3, 4]]

        """

        if self.inner_size() != 0:
            raise ValueError, "the inner size of the skew tableau must be 0"
        else:
            return  sage.combinat.tableau.Tableau(self[:])

    def restrict(self, n):
        """
        Returns the restriction of the standard skew tableau to all the numbers less
        than or equal to n.

        EXAMPLES:
            sage: SkewTableau([[None,1],[2],[3]]).restrict(2)
            [[None, 1], [2]]
            sage: SkewTableau([[None,1],[2],[3]]).restrict(1)
            [[None, 1]]
        """
        t = self[:]
        if not self.is_standard():
            raise ValueError, "the skew tableau must be standard to perform the restriction"

        return SkewTableau( filter(lambda z: z != [], map(lambda x: filter(lambda y: y is None or y <= n, x), t)) )

    def to_chain(self):
        """
        Returns the chain of skew partitions corresponding to the standard
        skew tableau.

        EXAMPLES:
            sage: SkewTableau([[None,1],[2],[3]]).to_chain()
            [[[1], [1]], [[2], [1]], [[2, 1], [1]], [[2, 1, 1], [1]]]
        """
        if not self.is_standard():
            raise ValueError, "the skew tableau must be standard to convert to a chain"

        return map(lambda x: self.restrict(x).shape(), range(self.size()+1))

    def slide(self, corner=None):
        """

        Fulton, William. 'Young Tableaux'. p12-13

        EXAMPLES:
            sage: st = SkewTableau([[None, None, None, None,2],[None, None, None, None,6], [None, 2, 4, 4], [2, 3, 6], [5,5]])
            sage: st.slide([2,0])
            [[None, None, None, None, 2], [None, None, None, None, 6], [2, 2, 4, 4], [3, 5, 6], [5]]

        """
        new_st = copy.copy(self[:])
        inner_corners = self.inner_shape().corners()
        outer_corners = self.outer_shape().corners()
        if corner is not None:
            if corner not in inner_corners:
                raise ValueError, "corner must be an inner corner"
        else:
            if len(inner_corners) == 0:
                return self
            else:
                corner = inner_corners[0]

        spot = corner[:]
        while spot not in outer_corners:
            #print spot
            #Check to see if there is nothing to the right
            if spot[1] == len(new_st[spot[0]]) - 1:
                #print "nr"
                #Swap the hole with the box below
                new_st[spot[0]][spot[1]] = new_st[spot[0]+1][spot[1]]
                new_st[spot[0]+1][spot[1]] = None
                spot[0] += 1
                continue

            #Check to see if there is nothing below
            if (spot[0] == len(new_st) - 1) or (len(new_st[spot[0]+1]) <= spot[1]):
                #print "nb"
                #Swap the hole with the box to the right
                new_st[spot[0]][spot[1]] = new_st[spot[0]][spot[1]+1]
                new_st[spot[0]][spot[1]+1] = None
                spot[1] += 1
                continue

            #If we get to this stage, we need to compare
            below = new_st[spot[0]+1][spot[1]]
            right = new_st[spot[0]][spot[1]+1]
            if below <= right:
                #Swap with the box below
                #print "b"
                new_st[spot[0]][spot[1]] = new_st[spot[0]+1][spot[1]]
                new_st[spot[0]+1][spot[1]] = None
                spot[0] += 1
                continue
            else:
                #Swap with the box to the right
                #print "r"
                new_st[spot[0]][spot[1]] = new_st[spot[0]][spot[1]+1]
                new_st[spot[0]][spot[1]+1] = None
                spot[1] += 1
                continue

        #Clean up to remove the "None" at an outside corner
        #Remove the last row if there is nothing left in it
        new_st[spot[0]].pop()
        if len(new_st[spot[0]]) == 0:
            new_st.pop()

        return SkewTableau(new_st)


    def rectify(self):
        """
        Returns a Tableau formed by applying the jeu de taquin process to
        self.

        Fulton, William. 'Young Tableaux'. p15

        EXAMPLES:
        """
        rect = copy.deepcopy(self)
        inner_corners = rect.inner_shape().corners()

        while len(inner_corners) > 0:
            rect = rect.slide()
            inner_corners = rect.inner_shape().corners()

        return rect.to_tableau()




class SkewTableauWithContent(SkewTableau_class):
    def __init__(self, outer, inner, conts=None, maxrows=0):
        self.rows = len(outer)

        if conts is not None:
            self.clen = len(conts)
        else:
            self.clen = 0
        self.clen += self.rows

        self.conts = [0]*self.clen
        if conts is not None:
            for i in range(len(conts)):
                self.conts[i] = conts[i]

        self.outer = outer[:]
        self.inner = inner[:] + [0]*(len(outer)-len(inner))

        if len(self.outer) == 0:
            self.cols = 0
        else:
            self.cols = self.outer[0]


        initial = [ [0]*self.outer[i] for i in range(len(self.outer)) ]
        t = self._setmin(0, self.outer[0]-1, initial=initial)
        SkewTableau.__init__(self,t)

        self.outer_conj = None

        self.maxrows = maxrows
        if (maxrows >= self.clen):
            self.maxrows = 0

        if self.maxrows == 0:
            return

        if self.conts[self.maxrows] != 0:
            raise ValueError, " self.conts[self.maxrows] != 0 "

        self.outer_conj = partition.Partition(self.outer).conjugate()

    def _setmin(self, y, x, initial=None):
        #Compute t which is the "minimum" skew tableau with
        #the given content
        #t = [ [0]*self.outer[i] for i in range(len(self.outer)) ]
        if initial is None:
            t = self.list
        else:
            t = initial

        #this algorithm comes from st_setmin
        #print "Call: (%d, %d)"%(y,x)
        while ( y < self.rows ):
            while x >= self.inner[y]:
                if y == 0 or x < self.inner[y-1]:
                    e = 0
                else:
                    e = t[y-1][x] + 1

                t[y][x] = e
                try:
                    self.conts[e] += 1
                except:
                    raise IndexError, "(%d,%d): %d, %s"%(y,x, e, str(self.conts))

                x -= 1
            y += 1

            if y < self.rows:
                x = self.outer[y] - 1
        return t

    def next(self):
        """
        Returns the next semistandard skew tableau
        """
        #new_st = copy.copy(self)
        #IMPLEMENTATION SPECIFIC
        l = self.list

        for y in reversed(range(self.rows)):
            xlim = self.outer[y]

            for x in range(self.inner[y], xlim):
                if self.maxrows == 0:
                    elim = len(self.conts) -1
                else:
                    elim = self.maxrows + y - self.outer_conj[x]

                if x != xlim - 1:
                    next_cell = l[y][x+1]
                    if next_cell < elim:
                        elim = next_cell

                e = l[y][x]
                self.conts[e] -= 1

                e += 1
                while ( e <= elim and
                        self.conts[e] == self.conts[e-1] ):
                    e += 1

                if e <= elim:
                    l[y][x] = e
                    self.conts[e] += 1

                    #IMPLEMENTATION SPECIFIC
                    #new_st._setmin(y, x-1)
                    self._setmin(y,x-1)
                    return self
                    #return new_st

        return False


def st_iterator(outer, inner, conts=None, maxrows=0):
    st = SkewTableauWithContent(outer,inner,conts=conts,maxrows=maxrows)
    yield st

    next = st.next()
    while next != False:
        yield next
        next = next.next()



def _label_skew(list, sk):
    """
    Returns a filled in a standard skew tableaux given an ordered list of the
    coordinates to filled in.

    EXAMPLES:
        sage: l = [ '0,0', '1,1', '1,0', '0,1' ]
        sage: empty = [[None,None],[None,None]]
        sage: skew_tableau._label_skew(l, empty)
        [[1, 4], [3, 2]]
    """
    i = 1
    skew = copy.deepcopy(sk)
    for coordstring in list:
            coords = coordstring.split(",")
            row = int(coords[0])
            column = int(coords[1])
            skew[row][column] = i
            i += 1
    return skew

def StandardSkewTableaux(skp):
    """
    Returns the combinatorial class of standard skew tableaux of
    shape skp (where skp is a skew partition).

    EXAMPLES:
        sage: StandardSkewTableaux([[3, 2, 1], [1, 1]]).list()
        [[[None, 1, 2], [None, 3], [4]],
         [[None, 1, 2], [None, 4], [3]],
         [[None, 1, 3], [None, 2], [4]],
         [[None, 1, 4], [None, 2], [3]],
         [[None, 1, 3], [None, 4], [2]],
         [[None, 1, 4], [None, 3], [2]],
         [[None, 2, 3], [None, 4], [1]],
         [[None, 2, 4], [None, 3], [1]]]
    """
    return StandardSkewTableaux_skewpartition(SkewPartition(skp))

class StandardSkewTableaux_skewpartition(CombinatorialClass):
    object_class = SkewTableau_class
    def __init__(self, skp):
        """
        TESTS:
            sage: S = StandardSkewTableaux([[3, 2, 1], [1, 1]])
            sage: S == loads(dumps(S))
            True
        """
        self.skp = skp

    def list(self):
        """
        Returns a list for all the standard skew tableaux with shape
        of the skew partition skp. The standard skew tableaux are
        ordered lexicographically by the word obtained from their
        row reading.

        EXAMPLES:
            sage: StandardSkewTableaux([[3, 2, 1], [1, 1]]).list()
            [[[None, 1, 2], [None, 3], [4]],
             [[None, 1, 2], [None, 4], [3]],
             [[None, 1, 3], [None, 2], [4]],
             [[None, 1, 4], [None, 2], [3]],
             [[None, 1, 3], [None, 4], [2]],
             [[None, 1, 4], [None, 3], [2]],
             [[None, 2, 3], [None, 4], [1]],
             [[None, 2, 4], [None, 3], [1]]]

        """
        return [st for st in self]

    def count(self):
        """
        Returns the number of standard skew tableaux with shape of the
        skew partition skp.

        EXAMPLES:
            sage: StandardSkewTableaux([[3, 2, 1], [1, 1]]).count()
            8
        """

        return sum([1 for st in self])

    def iterator(self):
        """
        An iterator for all the standard skew tableau with shape of the
        skew partition skp. The standard skew tableaux are
        ordered lexicographically by the word obtained from their
        row reading.

        EXAMPLES:
            sage: [st for st in StandardSkewTableaux([[3, 2, 1], [1, 1]])]
            [[[None, 1, 2], [None, 3], [4]],
             [[None, 1, 2], [None, 4], [3]],
             [[None, 1, 3], [None, 2], [4]],
             [[None, 1, 4], [None, 2], [3]],
             [[None, 1, 3], [None, 4], [2]],
             [[None, 1, 4], [None, 3], [2]],
             [[None, 2, 3], [None, 4], [1]],
             [[None, 2, 4], [None, 3], [1]]]
        """
        skp = self.skp

        dag = skp.to_dag()
        le_list = list(dag.topological_sort_generator())

        empty = [[None]*row_length for row_length in skp.outer()]

        for le in le_list:
            yield SkewTableau(_label_skew(le, empty))



def from_expr(expr):
    """
    Returns a SkewTableau from a MuPAD-Combinat expr for a
    skew tableau.  The first list in expr is the inner shape
    of the skew tableau.  The second list are the entries
    in the rows of the skew tableau from bottom to top.

    Provided primarily for compatability with MuPAD-Combinat.

    EXAMPLES:
        sage: sage.combinat.skew_tableau.from_expr([[1,1],[[5],[3,4],[1,2]]])
        [[None, 1, 2], [None, 3, 4], [5]]
    """
    skp = []
    outer = expr[1]
    inner = expr[0]+[0]*(len(outer)-len(expr[0]))

    for i in range(len(outer)):
        skp.append( [None]*(inner[i]) + outer[-(i+1)] )

    return SkewTableau(skp)



def from_shape_and_word(shape, word):
    """
    Returns the skew tableau correspnding to the skew partition
    shape and the word obtained from the row reading.

    EXAMPLES:
        sage: t = SkewTableau([[None, 2, 4], [None, 3], [1]])
        sage: shape = t.shape()
        sage: word  = t.to_word()
        sage: skew_tableau.from_shape_and_word(shape, word)
        [[None, 2, 4], [None, 3], [1]]

    """
    st = [ [None]*row_length for row_length in shape[0] ]
    w_count = 0
    for i in range(len(shape[0])):
        for j in range(shape[0][i]):
            if i >= len(shape[1]) or j >= shape[1][i]:
                st[i][j] = word[w_count]
                w_count += 1
    return SkewTableau(st)


#####################
# Under development #
#####################

def SemistandardSkewTableaux(shape, weight):
    return SemistandardSkewTableaux_shapeweight(shape, weight)

class SemistandardSkewTableaux_shapeweight(CombinatorialClass):
    def __init__(self, shape, weight):
        self.shape = shape
        self.weight = weight

    def iterator(self):
        for skt in RibbonTableaux(self.shape, self.weight, 1):
            yield SkewTableau_class(skt)
