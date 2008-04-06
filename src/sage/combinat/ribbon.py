r"""
Ribbons
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

import sage.combinat.misc as misc
import sage.combinat.skew_tableau
import sage.combinat.word as word
import sage.combinat.permutation as permutation
import sage.combinat.partition as partition
import sage.combinat.skew_tableau as skew_tableau
import sage.combinat.skew_partition
import sage.rings.integer
from combinat import CombinatorialClass, CombinatorialObject

SkewTableau   = sage.combinat.skew_tableau.SkewTableau
from_expr     = sage.combinat.skew_tableau.from_expr
SkewPartition = sage.combinat.skew_partition.SkewPartition
Integer       = sage.rings.integer.Integer

def Ribbon(r):
    """
    Returns a ribbon tableau object.

    A ribbon is a skew tableau that does not contain a 2x2 box.
    A ribbon is given by a list of the rows from top to bottom.


    EXAMPLES:
        sage: Ribbon([[2,3],[1,4,5]])
        [[2, 3], [1, 4, 5]]
        sage: Ribbon([[2,3],[1,4,5]]).to_skew_tableau()
        [[None, None, 2, 3], [1, 4, 5]]

    """
    if isinstance(r, list):
        if len(r) == 0:
            return Ribbon_class(r)
        if all(isinstance(i, list) for i in r):
            if all(all(isinstance(j, (int, Integer)) and j>0 for j in i) for i in r):
                return Ribbon_class(r)
    raise TypeError, "r must be a list of positive integers"

class Ribbon_class(CombinatorialObject):
    def ribbon_shape(self):
        """
        Returns the ribbon shape.  The ribbon shape is given
        just by the number of boxes in each row.

        EXAMPLES:
            sage: Ribbon([[2,3],[1,4,5]]).ribbon_shape()
            [2, 3]
        """

        return [len(row) for row in self]

    def height(self):
        """
        Returns the height of the ribbon.

        EXAMPLES:
            sage: Ribbon([[2,3],[1,4,5]]).height()
            2
        """
        return len(self)

    def spin(self):
        """
        Returns the spin of self.

        EXAMPLES:
            sage: Ribbon([[2,3],[1,4,5]]).spin()
            1/2
        """
        return Integer(self.height()-1)/2

    def width(self):
        """
        Returns the width of the ribbon.

        EXAMPLES:
            sage: Ribbon([[2,3],[1,4,5]]).width()
            4
        """
        return 1+sum([len(r)-1 for r in self])

    def size(self):
        """
        Returns the size ( number of boxes ) in the ribbon.

        EXAMPLES:
            sage: Ribbon([[2,3],[1,4,5]]).size()
            5
        """
        return sum([len(r) for r in self])

    def is_standard(self):
        """
        Returns True is the ribbon is standard and False otherwise.
        ribbon are standard if they are filled with the numbers
        1...size and they are increasing along both rows and columns.

        EXAMPLES:
            sage: Ribbon([[2,3],[1,4,5]]).is_standard()
            True
            sage: Ribbon([[2,2],[1,4,5]]).is_standard()
            False
            sage: Ribbon([[4,5],[1,2,3]]).is_standard()
            False
        """

        return self.to_skew_tableau().is_standard()

    def to_skew_tableau(self):
        """
        Returns the skew tableau corresponding to the ribbon
        tableau.

        EXAMPLES:
            sage: Ribbon([[2,3],[1,4,5]]).to_skew_tableau()
            [[None, None, 2, 3], [1, 4, 5]]
        """
        st = []
        space_count = 0
        for row in reversed(self):
            st.append( [None]*space_count + row )
            space_count += len(row) - 1
        st.reverse()
        return sage.combinat.skew_tableau.SkewTableau(st)

    def to_permutation(self):
        """
        Returns the permutation corresponding to the ribbon
        tableau.

        EXAMPLES:
            sage: r = Ribbon([[1], [2,3], [4, 5, 6]])
            sage: r.to_permutation()
            [4, 5, 6, 2, 3, 1]
        """
        return permutation.Permutation(self.to_word())

    def shape(self):
        """
        Returns the skew partition corresponding to the shape of the
        ribbon.

        EXAMPLES:
            sage: Ribbon([[2,3],[1,4,5]]).shape()
            [[4, 3], [2]]
        """
        return self.to_skew_tableau().shape()

    def to_word_by_row(self):
        """
        Returns a word obtained from a row reading of the ribbon.

        EXAMPLES:
            sage: Ribbon([[1],[2,3]]).to_word_by_row()
            [2, 3, 1]
            sage: Ribbon([[2, 4], [3], [1]]).to_word_by_row()
            [1, 3, 2, 4]
        """
        word = []
        for row in self:
            word = row + word
        return word


    def to_word_by_column(self):
        """
        Returns the word obtained from a column reading of the ribbon

        EXAMPLES:
            sage: Ribbon([[1],[2,3]]).to_word_by_column()
            [1, 3, 2]
            sage: Ribbon([[2, 4], [3], [1]]).to_word_by_column()
            [4, 2, 3, 1]

        """
        return self.to_skew_tableau().to_word_by_column()

    def to_word(self):
        """
        An alias for Ribbon.to_word_by_row().

        EXAMPLES:
            sage: Ribbon([[1],[2,3]]).to_word_by_row()
            [2, 3, 1]
            sage: Ribbon([[2, 4], [3], [1]]).to_word_by_row()
            [1, 3, 2, 4]
        """
        return self.to_word_by_row()

    def evaluation(self):
        """
        Returns the evaluation of the word from ribbon.

        EXAMPLES:
            sage: Ribbon([[1,2],[3,4]]).evaluation()
            [1, 1, 1, 1]
        """

        return word.evaluation(self.to_word())



def from_shape_and_word(shape, word):
    """
    Returns the ribbon corresponding to the given
    ribbon shape and word.

    EXAMPLES:
        sage: import sage.combinat.ribbon as ribbon
        sage: ribbon.from_shape_and_word([2,3],[1,2,3,4,5])
        [[1, 2], [3, 4, 5]]
    """
    pos = 0
    r = []
    for l in shape:
        r.append(word[pos:pos+l])
        pos += l
    return Ribbon(r)

def StandardRibbons(shape):
    """
    Returns the combinatorial class of standard ribbon
    tableaux of shape shape.

    EXAMPLES:
        sage: StandardRibbons([2,2])
        Standard ribbon tableaux of shape [2, 2]
        sage: StandardRibbons([2,2]).first()
        [[2, 4], [1, 3]]
        sage: StandardRibbons([2,2]).last()
        [[1, 2], [3, 4]]
        sage: StandardRibbons([2,2]).count()
        5
        sage: StandardRibbons([2,2]).list()
        [[[2, 4], [1, 3]],
         [[2, 3], [1, 4]],
         [[1, 4], [2, 3]],
         [[1, 3], [2, 4]],
         [[1, 2], [3, 4]]]
        sage: StandardRibbons([3,2,2]).count()
        155

    """
    return StandardRibbons_shape(shape)

class StandardRibbons_shape(CombinatorialClass):
    def __init__(self, shape):
        """
        TESTS:
            sage: S = StandardRibbons([2,2])
            sage: S == loads(dumps(S))
            True
        """
        self.shape = shape

    def __repr__(self):
        """
        TESTS:
            sage: repr(StandardRibbons([2,2]))
            'Standard ribbon tableaux of shape [2, 2]'
        """
        return "Standard ribbon tableaux of shape %s"%self.shape


    def first(self):
        """
        Returns the first standard ribbon of
        ribbon shape shape.

        EXAMPLES:
            sage: StandardRibbons([2,2]).first()
            [[2, 4], [1, 3]]
        """
        return from_permutation(permutation.descents_composition_first(self.shape))

    def last(self):
        """
        Returns the first standard ribbon of
        ribbon shape shape.

        EXAMPLES:
            sage: StandardRibbons([2,2]).last()
            [[1, 2], [3, 4]]
        """
        return from_permutation(permutation.descents_composition_last(self.shape))


    def iterator(self):
        """
        An iterator for the standard ribbon of ribbon
        shape shape.

        EXAMPLES:
            sage: [t for t in StandardRibbons([2,2])]
            [[[2, 4], [1, 3]],
             [[2, 3], [1, 4]],
             [[1, 4], [2, 3]],
             [[1, 3], [2, 4]],
             [[1, 2], [3, 4]]]

        """

        for p in permutation.descents_composition_list(self.shape):
            yield from_permutation(p)

def from_permutation(p):
    """
    Returns a standard ribbon of size len(p) from a Permutation p.
    The lengths of each row are given by the distance between the descents
    of the permutation p.

    EXAMPLES:
        sage: import sage.combinat.ribbon as ribbon
        sage: [ribbon.from_permutation(p) for p in Permutations(3)]
        [[[1, 2, 3]], [[2], [1, 3]], [[1, 3], [2]], [[1], [2, 3]], [[1, 2], [3]], [[1], [2], [3]]]

    """
    if p == []:
        return Ribbon([])

    comp = p.descents()

    if comp == []:
        return Ribbon([p[:]])


    #[p[j]$j=compo[i]+1..compo[i+1]] $i=1..nops(compo)-1, [p[j]$j=compo[nops(compo)]+1..nops(p)]
    r = []
    r.append([p[j] for j in range(comp[0]+1)])
    for i in range(len(comp)-1):
        r.append([ p[j] for j in range(comp[i]+1,comp[i+1]+1) ])
    r.append( [ p[j] for j in range(comp[-1]+1, len(p))] )
    r.reverse()
    return Ribbon(r)
