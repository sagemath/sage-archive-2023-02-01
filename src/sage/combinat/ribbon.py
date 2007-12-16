r"""
Ribbon Tableaux
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
from combinat import CombinatorialClass, CombinatorialObject

def Ribbon(r):
    """
    Returns a ribbon tableau object.

    EXAMPLES:
        sage: Ribbon([[2,3],[1,4,5]])
        [[2, 3], [1, 4, 5]]
    """
    return Ribbon_class(r)

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
            [1, 2, 3, 4, 5, 6]
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
            [1, 2, 3]
            sage: Ribbon([[2, 4], [3], [1]]).to_word_by_row()
            [2, 4, 3, 1]
        """
        word = []
        for row in self:
            word += row

        return word


    def to_word_by_column(self):
        """
        Returns the word obtained from a column reading of the ribbon

        EXAMPLES:
            sage: Ribbon([[1],[2,3]]).to_word_by_column()
            [2, 1, 3]
            sage: Ribbon([[2, 4], [3], [1]]).to_word_by_column()
            [2, 3, 1, 4]

        """
        return self.to_skew_tableau().to_word_by_column()

    def to_word(self):
        """
        An alias for Ribbon.to_word_by_row().

        EXAMPLES:
            sage: Ribbon([[1],[2,3]]).to_word_by_row()
            [1, 2, 3]
            sage: Ribbon([[2, 4], [3], [1]]).to_word_by_row()
            [2, 4, 3, 1]
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
        sage: ribbon.from_shape_and_word([2,3],[1,2,3,4,5])
        [[1, 2], [3, 4, 5]]
    """
    pos = 0
    r = []
    for l in shape:
        r.append(word[pos:pos+l])
        pos += l
    return Ribbon(r)

def StandardRibbonTableaux(shape):
    """
    Returns the combinatorial class of standard ribbon
    tableaux of shape shape.

    EXAMPLES:
        sage: StandardRibbonTableaux([2,2])
        Standard ribbon tableaux of shape [2, 2]
        sage: StandardRibbonTableaux([2,2]).first()
        [[1, 3], [2, 4]]
        sage: StandardRibbonTableaux([2,2]).last()
        [[3, 4], [1, 2]]
        sage: StandardRibbonTableaux([2,2]).count()
        5
        sage: StandardRibbonTableaux([2,2]).list()
        [[[1, 3], [2, 4]],
         [[1, 4], [2, 3]],
         [[2, 3], [1, 4]],
         [[2, 4], [1, 3]],
         [[3, 4], [1, 2]]]

        sage: StandardRibbonTableaux([2,2]).list()
        [[[1, 3], [2, 4]],
         [[1, 4], [2, 3]],
         [[2, 3], [1, 4]],
         [[2, 4], [1, 3]],
         [[3, 4], [1, 2]]]
        sage: StandardRibbonTableaux([3,2,2]).count()
        155

    """
    return StandardRibbonTableaux_shape(shape)

class StandardRibbonTableaux_shape(CombinatorialClass):
    def __init__(self, shape):
        """
        TESTS:
            sage: S = StandardRibbonTableaux([2,2])
            sage: S == loads(dumps(S))
            True
        """
        self.shape = shape

    def __repr__(self):
        """
        TESTS:
            sage: repr(StandardRibbonTableaux([2,2]))
            'Standard ribbon tableaux of shape [2, 2]'
        """
        return "Standard ribbon tableaux of shape %s"%self.shape


    def first(self):
        """
        Returns the first standard ribbon of
        ribbon shape shape.

        EXAMPLES:
            sage: StandardRibbonTableaux([2,2]).first()
            [[1, 3], [2, 4]]

        """
        return from_permutation(permutation.descents_composition_first(self.shape))

    def last(self):
        """
        Returns the first standard ribbon of
        ribbon shape shape.

        EXAMPLES:
            sage: StandardRibbonTableaux([2,2]).last()
            [[3, 4], [1, 2]]
        """
        return from_permutation(permutation.descents_composition_last(self.shape))


    def iterator(self):
        """
        An iterator for the standard ribbon of ribbon
        shape shape.

        EXAMPLES:
            sage: [t for t in StandardRibbonTableaux([2,2])]
            [[[1, 3], [2, 4]],
             [[1, 4], [2, 3]],
             [[2, 3], [1, 4]],
             [[2, 4], [1, 3]],
             [[3, 4], [1, 2]]]
        """

        for p in permutation.descents_composition_list(self.shape):
            yield from_permutation(p)

def from_permutation(p):
    """
    Returns a standard ribbon of size len(p) from a Permutation p.
    The lengths of each row are given by the distance between the descents
    of the permutation p.

    EXAMPLES:
        sage: [ribbon.from_permutation(p) for p in Permutations(3)]
        [[[1, 2, 3]],
         [[1, 3], [2]],
         [[2], [1, 3]],
         [[2, 3], [1]],
         [[3], [1, 2]],
         [[3], [2], [1]]]

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

    return Ribbon(r)



#####################
# Under Development #
#####################
class RibbonTableaux_shapeweightlength(CombinatorialClass):
    def __init__(self, shape, weight, length):
        self.shape = shape
        self.weight = weight
        self.length = length

    def list(self):
        pass









def insertion_tableau(skp, perm, evaluation, tableau, length):
    """

    INPUT:
        skp -- skew partitions
        perm, eval -- non-negative integers
        tableau -- skew tableau
        length -- integer

    """
    psave = skp[1]
    partc = skp[1] + [0]*(len(skp[0])-len(skp[1]))

    tableau = tableau.to_expr()[1]
    for k in range(len(tableau)):
        tableau[-(k+1)] += [0]* ( skp[0][k] - partc[k] - len(tableau[-(k+1)]))

    ## We construct a tableau from the southwest corner to the northeast one
    #[ op( revert([[0$(partition[1][k] - partc[k]) ] $k = nops(tableau)+1..nops(partition[1])])) , op(tableau) ];
    tableau =  [[0]*(skp[0][k] - partc[k]) for k in reversed(range(len(tableau), len(skp[0])+1))] + tableau

    tableau = SkewTableau(expr=[skp[1], tableau]).conjugate()
    tableau = tableau.to_expr()[1]

    skp = skp.conjugate().to_list()
    skp[1].append( [0]*(len(skp[0])-len(skp[1])) )

    if len(perm) > len(skp[0]):
        return None

    for k in range(len(perm)):
        if perm[ -(k+1) ] !=0:
            #tableau ... = evaluation
            pass

    return SkewTableau(expr=[psave.conjugate(),tableau]).conjugate()




def list_rec(nexts, current, part, weight, length):
    """

    INPUT:
        nexts, current, part -- skew partitions
        weight -- non-negative integer list
        length -- integer
    """
    if current == [] and nexts == [] and weight == []:
        return [part[1],[]]

    ## Test if the current nodes is not an empty node
    if current == []:
        return None


    ## Test if the current nodes drive us to new solutions
    if next != []:
        res = []
        for i in range(len(current)):
            for j in range(len(nexts[i])):
                res.append( interstion_tableau(part, current[i][1], len(weight), nexts[i][j], length) )
        return res
    else:
        ## The current nodes are at the bottom of the tree
        res = []
        for i in range(len(current)):
            res.append( insertion_tableau(part, current[i][1], len(weight), [[],[]], length) )
        return res


##     //////////////////////////////////////////////////////////////////////////////////////////
##     // Generic function for driving into the graph of partitions coding all ribbons
##     // tableaux of a given shape and weight
##     //////////////////////////////////////////////////////////////////////////////////////////

##     //This function construct the graph of the set of k-ribbon tableaux
##     //of a given skew shape and a given weight.
##     //The first argument is alaways a skew partition.
##     //In the case where the inner partition is empty there is no branch without solutions
##     //In the other cases there is in average a lot of branches without solutions
##     /////////////////////////////////////////////////////////////////////////////////////////



def graph_implementation_rec(skp, weight, lenth, function):

    if sum(weight) == 0:
        weight = []

    partp = Partition(part[1]).conjugate()

    ## Some tests in order to know if the shape and the weight are compatible.
    if weight != [] and weight[-1] <= len(partp):
        perms = Permutations([0]*(len(partp)-weight[-1]) + [length]*(weight[-1]).list())
    else:
        return function([], [], part, weight, length)

    selection = []

    for j in range(perms):
        retire = [(partp[i]+ len(partp) - (i+1) - perms[j][i]) for i in range(partp)]
        retire.sort(reverse=True)
        retire = [ retire[i] - len(partp) + (i+1) for i in range(len(retire))]

        if retire[-1] >= 0 and retire == [i for i in reversed(sorted(retire))]:
            retire = Partition(filter(lambda x: x != 0, retire)).conjugate()


            # Cutting branches if the retired partition has a line strictly included into the inner one
            append = True
            padded_retire = retire + [0]*(len(part[1])-len(retire))
            for k in range(len(part[1])):
                if padded_retire[k] - part[1][k] < 0:
                    append = False
                    break
            if append:
                selection.append([retire, perms[j]])


    #selection contains the list of current nodes


    if len(weight) == 1:
        return function([], selection, part, weight, length)
    else:
        #The recursive calls permit us to construct the list of the sons
        #of all current nodes in selection
        a = [graph_implementation_rec([p[0], part[1]], [weight[i]]*(len(weight)-1), length, function) for p in selection]
        return function(a, selection, part, weight, length)














