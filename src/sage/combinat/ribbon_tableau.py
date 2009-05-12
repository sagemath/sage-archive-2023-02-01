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

from combinat import CombinatorialObject, CombinatorialClass
import skew_tableau, permutation, partition, skew_partition
from sage.rings.all import QQ, ZZ
import sage.calculus.calculus
import functools
from sage.combinat.words.words import Words

def RibbonTableau(rt=None, expr=None):
    """
    Returns a ribbon tableau object.

    EXAMPLES::

        sage: rt = RibbonTableau([[None, 1],[2,3]]); rt
        [[None, 1], [2, 3]]
        sage: rt.inner_shape()
        [1]
        sage: rt.outer_shape()
        [2, 2]

    ::

        sage: RibbonTableau(expr=[[1,1],[[5],[3,4],[1,2]]])
        [[None, 1, 2], [None, 3, 4], [5]]
    """
    if expr is not None:
        return from_expr(expr)

    for row in rt:
        if not isinstance(row, list):
            raise TypeError, "each element of the ribbon tableau must be a list"
        if row == []:
            raise TypeError, "a ribbon tableau cannot have an empty list for a row"

    return RibbonTableau_class(rt)

class RibbonTableau_class(skew_tableau.SkewTableau_class):
    def length(self):
        """
        Returns the length of the ribbons into a ribbon tableau.

        EXAMPLES::

            sage: RibbonTableau([[None, 1],[2,3]]).length()
            1
            sage: RibbonTableau([[1,0],[2,0]]).length()
            2
        """
        if self.to_expr() == [[],[]]:
            return 0

        tableau = self.to_expr()[1]
        l = 0
        t = 0
        for k in range(len(tableau)):
            t += len( [ x for x in tableau[k] if x is not None and x > -1 ] )
            l += len( [ x for x in tableau[k] if x is not None and x > 0  ] )

        if l == 0:
            return t
        else:
            return t/l

    def to_word(self):
        """
        Returns a word obtained from a row reading of self.

        EXAMPLES::

            sage: R = RibbonTableau([[0, 0, 3, 0], [1, 1, 0], [2, 0, 4]])
            sage: R.to_word()
            word: 2041100030
        """
        w = []
        for row in reversed(self):
            w += row
        return Words(alphabet="natural numbers")(w)

    def evaluation(self):
        """
        Returns the evaluation of the ribbon tableau

        EXAMPLES::

            sage: RibbonTableau([[0, 0, 3, 0], [1, 1, 0], [2, 0, 4]]).evaluation()
            [2, 1, 1, 1]
        """
        w = [i for i in self.to_word() if i != 0]
        word = Words(alphabet="positive integers")(w)
        return word.evaluation()

def from_expr(l):
    """
    Returns a RibbonTableau from a MuPAD-Combinat expr for a skew
    tableau. The first list in expr is the inner shape of the skew
    tableau. The second list are the entries in the rows of the skew
    tableau from bottom to top.

    Provided primarily for compatibility with MuPAD-Combinat.

    EXAMPLES::

        sage: import sage.combinat.ribbon_tableau as ribbon_tableau
        sage: sage.combinat.ribbon_tableau.from_expr([[1,1],[[5],[3,4],[1,2]]])
        [[None, 1, 2], [None, 3, 4], [5]]
        sage: type(_)
        <class 'sage.combinat.ribbon_tableau.RibbonTableau_class'>
    """
    return RibbonTableau_class(skew_tableau.from_expr(l)._list)


#####################
# Ribbon Tableaux   #
#####################

def RibbonTableaux(shape, weight, length):
    """
    Returns the combinatorial class of ribbon tableaux of skew shape
    shape and weight weight tiled by ribbons of length length.

    EXAMPLES::

        sage: RibbonTableaux([[2,1],[]],[1,1,1],1)
        Ribbon tableaux of shape [[2, 1], []] and weight [1, 1, 1] with 1-ribbons
    """
    if shape in partition.Partitions():
        shape = partition.Partition(shape)
        shape = skew_partition.SkewPartition([shape, shape.r_core(length)])
    else:
        shape = skew_partition.SkewPartition(shape)

    if shape.size() != length*sum(weight):
        raise ValueError

    weight = [i for i in weight if i != 0]

    return RibbonTableaux_shapeweightlength(shape, weight,length)


class RibbonTableaux_shapeweightlength(CombinatorialClass):
    def __init__(self, shape, weight, length):
        """
        EXAMPLES::

            sage: r = RibbonTableaux([[2,1],[]],[1,1,1],1)
            sage: r == loads(dumps(r))
            True
        """
        self._shape  =  shape
        self._weight =  weight
        self._length =  length
        self._name = "Ribbon tableaux of shape %s and weight %s with %s-ribbons"%(shape, weight, length)


    def __contains__(self, x):
        """
        Note that this just checks to see if x appears in self.list(). This
        should be improved to provide actual checking.

        EXAMPLES::

            sage: r = RibbonTableaux([[2,2],[]],[1,1],2)
            sage: [[0, 0], [1, 2]] in r
            True
            sage: [[1, 0], [2, 0]] in r
            True
            sage: [[0, 1], [2, 0]] in r
            False
        """
        return x in self.list()

    def list(self):
        """
        EXAMPLES::

            sage: RibbonTableaux([[2,1],[]],[1,1,1],1).list()
            [[[1, 3], [2]], [[1, 2], [3]]]
            sage: RibbonTableaux([[2,2],[]],[1,1],2).list()
            [[[0, 0], [1, 2]], [[1, 0], [2, 0]]]
        """
        res = graph_implementation_rec(self._shape, self._weight, self._length, list_rec)
        return [from_expr(x) for x in res]

    def cardinality(self):
        """
        EXAMPLES::

            sage: RibbonTableaux([[2,1],[]],[1,1,1],1).cardinality()
            2
            sage: RibbonTableaux([[2,2],[]],[1,1],2).cardinality()
            2
            sage: RibbonTableaux([[4,3,3],[]],[2,1,1,1],2).cardinality()
            5

        TESTS::

            sage: RibbonTableaux([6,6,6], [4,2], 3).cardinality()
            6
            sage: RibbonTableaux([3,3,3,2,1], [3,1], 3).cardinality()
            1
            sage: RibbonTableaux([3,3,3,2,1], [2,2], 3).cardinality()
            2
            sage: RibbonTableaux([3,3,3,2,1], [2,1,1], 3).cardinality()
            5
            sage: RibbonTableaux([3,3,3,2,1], [1,1,1,1], 3).cardinality()
            12
            sage: RibbonTableaux([5,4,3,2,1],[2,2,1],3).cardinality()
            10

        ::

            sage: RibbonTableaux([8,7,6,5,1,1],[3,2,2,1],3).cardinality()
            85
            sage: RibbonTableaux([5,4,3,2,1,1,1],[2,2,1],3).cardinality()
            10

        ::

            sage: RibbonTableaux([7,7,7,2,1,1],[3,2,0,1,1],3).cardinality()
            25

        Weights with some zeros in the middle and end

        ::

            sage: RibbonTableaux([3,3,3], [0,1,0,2,0], 3).cardinality()
            3
            sage: RibbonTableaux([3,3,3],[1,0,1,0,1,0,0,0],3).cardinality()
            6
        """
        return graph_implementation_rec(self._shape, self._weight, self._length, count_rec)[0]

def insertion_tableau(skp, perm, evaluation, tableau, length):
    """
    INPUT:


    -  ``skp`` - skew partitions

    -  ``perm, evaluation`` - non-negative integers

    -  ``tableau`` - skew tableau

    -  ``length`` - integer


    TESTS::

        sage: from sage.combinat.ribbon_tableau import insertion_tableau
        sage: insertion_tableau([[1], []], [1], 1, [[], []], 1)
        [[], [[1]]]
        sage: insertion_tableau([[2, 1], []], [1, 1], 2, [[], [[1]]], 1)
        [[], [[2], [1, 2]]]
        sage: insertion_tableau([[2, 1], []], [0, 0], 3, [[], [[2], [1, 2]]], 1)
        [[], [[2], [1, 2]]]
        sage: insertion_tableau([[1, 1], []], [1], 2, [[], [[1]]], 1)
        [[], [[2], [1]]]
        sage: insertion_tableau([[2], []], [0, 1], 2, [[], [[1]]], 1)
        [[], [[1, 2]]]
        sage: insertion_tableau([[2, 1], []], [0, 1], 3, [[], [[2], [1]]], 1)
        [[], [[2], [1, 3]]]
        sage: insertion_tableau([[1, 1], []], [2], 1, [[], []], 2)
        [[], [[1], [0]]]
        sage: insertion_tableau([[2], []], [2, 0], 1, [[], []], 2)
        [[], [[1, 0]]]
        sage: insertion_tableau([[2, 2], []], [0, 2], 2, [[], [[1], [0]]], 2)
        [[], [[1, 2], [0, 0]]]
        sage: insertion_tableau([[2, 2], []], [2, 0], 2, [[], [[1, 0]]], 2)
        [[], [[2, 0], [1, 0]]]
        sage: insertion_tableau([[2, 2], [1]], [3, 0], 1, [[], []], 3)
        [[1], [[1, 0], [0]]]
    """
    psave = partition.Partition(skp[1])
    partc = skp[1] + [0]*(len(skp[0])-len(skp[1]))

    tableau = skew_tableau.SkewTableau(expr=tableau).to_expr()[1]

    for k in range(len(tableau)):
        tableau[-(k+1)] += [0]* ( skp[0][k] - partc[k] - len(tableau[-(k+1)]))

    ## We construct a tableau from the southwest corner to the northeast one
    tableau =  [[0]*(skp[0][k] - partc[k]) for k in reversed(range(len(tableau), len(skp[0])))] + tableau

    tableau = skew_tableau.from_expr([skp[1], tableau]).conjugate()
    tableau = tableau.to_expr()[1]

    skp = skew_partition.SkewPartition(skp).conjugate().to_list()
    skp[1].extend( [0]*(len(skp[0])-len(skp[1])) )

    if len(perm) > len(skp[0]):
        return None

    for k in range(len(perm)):
        if perm[ -(k+1) ] !=0:
            tableau[len(tableau)-len(perm)+k][ skp[0][len(perm)-(k+1)] - skp[1][ len(perm)-(k+1) ] - 1 ] = evaluation

    return skew_tableau.SkewTableau(expr=[psave.conjugate(),tableau]).conjugate().to_expr()


def count_rec(nexts, current, part, weight, length):
    """
    INPUT:


    -  ``nexts, current, part`` - skew partitions

    -  ``weight`` - non-negative integer list

    -  ``length`` - integer


    TESTS::

        sage: from sage.combinat.ribbon_tableau import count_rec
        sage: count_rec([], [], [[2, 1, 1], []], [2], 2)
        [0]
        sage: count_rec([[0], [1]], [[[2, 1, 1], [0, 0, 2, 0]], [[4], [2, 0, 0, 0]]], [[4, 1, 1], []], [2, 1], 2)
        [1]
        sage: count_rec([], [[[], [2, 2]]], [[2, 2], []], [2], 2)
        [1]
        sage: count_rec([], [[[], [2, 0, 2, 0]]], [[4], []], [2], 2)
        [1]
        sage: count_rec([[1], [1]], [[[2, 2], [0, 0, 2, 0]], [[4], [2, 0, 0, 0]]], [[4, 2], []], [2, 1], 2)
        [2]
        sage: count_rec([[1], [1], [2]], [[[2, 2, 2], [0, 0, 2, 0]], [[4, 1, 1], [0, 2, 0, 0]], [[4, 2], [2, 0, 0, 0]]], [[4, 2, 2], []], [2, 1, 1], 2)
        [4]
        sage: count_rec([[4], [1]], [[[4, 2, 2], [0, 0, 2, 0]], [[4, 3, 1], [0, 2, 0, 0]]], [[4, 3, 3], []], [2, 1, 1, 1], 2)
        [5]
    """
    if current == []:
        return [0]
    if nexts != []:
        return [sum(sum(j for j in i) for i in nexts)]
    else:
        return [len(current)]

def list_rec(nexts, current, part, weight, length):
    """
    INPUT:


    -  ``nexts, current, part`` - skew partitions

    -  ``weight`` - non-negative integer list

    -  ``length`` - integer


    TESTS::

        sage: from sage.combinat.ribbon_tableau import list_rec
        sage: list_rec([], [[[], [1]]], [[1], []], [1], 1)
        [[[], [[1]]]]
        sage: list_rec([[[[], [[1]]]]], [[[1], [1, 1]]], [[2, 1], []], [1, 2], 1)
        [[[], [[2], [1, 2]]]]
        sage: list_rec([], [[[1], [3, 0]]], [[2, 2], [1]], [1], 3)
        [[[1], [[1, 0], [0]]]]
        sage: list_rec([[[[], [[2]]]]], [[[1], [1, 1]]], [[2, 1], []], [0, 1, 2], 1)
        [[[], [[3], [2, 3]]]]
        sage: list_rec([], [[[], [2]]], [[1, 1], []], [1], 2)
        [[[], [[1], [0]]]]
        sage: list_rec([], [[[], [2, 0]]], [[2], []], [1], 2)
        [[[], [[1, 0]]]]
        sage: list_rec([[[[], [[1], [0]]]], [[[], [[1, 0]]]]], [[[1, 1], [0, 2]], [[2], [2, 0]]], [[2, 2], []], [1, 1], 2)
        [[[], [[1, 2], [0, 0]]], [[], [[2, 0], [1, 0]]]]
        sage: list_rec([], [[[], [2, 2]]], [[2, 2], []], [2], 2)
        [[[], [[1, 1], [0, 0]]]]
        sage: list_rec([], [[[], [1, 1]]], [[2], []], [2], 1)
        [[[], [[1, 1]]]]
        sage: list_rec([[[[], [[1, 1]]]]], [[[2], [1, 1]]], [[2, 2], []], [2, 2], 1)
        [[[], [[2, 2], [1, 1]]]]
    """
    if current == [] and nexts == [] and weight == []:
        return [[part[1],[]]]

    ## Test if the current nodes is not an empty node
    if current == []:
        return []

    ## Test if the current nodes drive us to new solutions
    if nexts != []:
        res = []
        for i in range(len(current)):
            for j in range(len(nexts[i])):
                res.append( insertion_tableau(part, current[i][1], len(weight), nexts[i][j], length) )
        return res
    else:
        ## The current nodes are at the bottom of the tree
        res = []
        for i in range(len(current)):
            res.append( insertion_tableau(part, current[i][1], len(weight), [[],[]], length) )
        return res



#############################
#Spin and Cospin Polynomials#
#############################
def spin_rec(t, nexts, current, part, weight, length):
    """
    Routine used for constructing the spin polynomial.

    INPUT:


    -  ``weight`` - list of non-negative integers

    -  ``length`` - the length of the ribbons we're tiling
       with

    -  ``t`` - the variable


    EXAMPLES::

        sage: from sage.combinat.ribbon_tableau import spin_rec
        sage: sp = SkewPartition
        sage: t = ZZ['t'].gen()
        sage: spin_rec(t, [], [[[], [3, 3]]], sp([[2, 2, 2], []]), [2], 3)
        [t^4]
        sage: spin_rec(t, [[0], [t^4]], [[[2, 1, 1, 1, 1], [0, 3]], [[2, 2, 2], [3, 0]]], sp([[2, 2, 2, 2, 1], []]), [2, 1], 3)
        [t^5]
        sage: spin_rec(t, [], [[[], [3, 3, 0]]], sp([[3, 3], []]), [2], 3)
        [t^2]
        sage: spin_rec(t, [[t^4], [t^3], [t^2]], [[[2, 2, 2], [0, 0, 3]], [[3, 2, 1], [0, 3, 0]], [[3, 3], [3, 0, 0]]], sp([[3, 3, 3], []]), [2, 1], 3)
        [t^6 + t^4 + t^2]
        sage: spin_rec(t, [[t^5], [t^4], [t^6 + t^4 + t^2]], [[[2, 2, 2, 2, 1], [0, 0, 3]], [[3, 3, 1, 1, 1], [0, 3, 0]], [[3, 3, 3], [3, 0, 0]]], sp([[3, 3, 3, 2, 1], []]), [2, 1, 1], 3)
        [2*t^7 + 2*t^5 + t^3]
    """
    from sage.combinat.words.word import Word
    R = ZZ['t']
    if current == []:
        return [R(0)]

    tmp = []
    partp = part[0].conjugate()

    #compute the contribution of the ribbons added at
    #the current node
    for perms in [current[i][1] for i in range(len(current))]:
        perm =  [partp[i] + len(partp)-(i+1)-perms[i] for i in range(len(partp))]
        perm.reverse()
        perm = Word(perm).standard_permutation()
        tmp.append( (weight[-1]*(length-1)-perm.number_of_inversions()) )

    if nexts != []:
        return [sum([sum([t**tmp[i]*nexts[i][j] for j in range(len(nexts[i]))]) for i in range(len(tmp))])]
    else:
        return [sum([t**tmp[i] for i in range(len(tmp))])]


def spin_polynomial_square(part, weight, length):
    """
    Returns the spin polynomial associated with part, weight, and
    length, with the substitution t -> t^2 made.

    EXAMPLES::

        sage: from sage.combinat.ribbon_tableau import spin_polynomial_square
        sage: spin_polynomial_square([6,6,6],[4,2],3)
        t^12 + t^10 + 2*t^8 + t^6 + t^4
        sage: spin_polynomial_square([6,6,6],[4,1,1],3)
        t^12 + 2*t^10 + 3*t^8 + 2*t^6 + t^4
        sage: spin_polynomial_square([3,3,3,2,1], [2,2], 3)
        t^7 + t^5
        sage: spin_polynomial_square([3,3,3,2,1], [2,1,1], 3)
        2*t^7 + 2*t^5 + t^3
        sage: spin_polynomial_square([3,3,3,2,1], [1,1,1,1], 3)
        3*t^7 + 5*t^5 + 3*t^3 + t
        sage: spin_polynomial_square([5,4,3,2,1,1,1], [2,2,1], 3)
        2*t^9 + 6*t^7 + 2*t^5
        sage: spin_polynomial_square([[6]*6, [3,3]], [4,4,2], 3)
        3*t^18 + 5*t^16 + 9*t^14 + 6*t^12 + 3*t^10
    """
    R = ZZ['t']
    t = R.gen()

    if part in partition.Partitions():
        part = skew_partition.SkewPartition([part,partition.Partition_class([])])
    elif part in skew_partition.SkewPartitions():
        part = skew_partition.SkewPartition(part)

    if part == [[],[]] and weight == []:
        return t.parent()(1)

    return R(graph_implementation_rec(part, weight, length, functools.partial(spin_rec,t))[0])

def spin_polynomial(part,weight,length):
    """
    Returns the spin polynomial associated to part, weight, and
    length.

    EXAMPLES::

        sage: from sage.combinat.ribbon_tableau import spin_polynomial
        sage: spin_polynomial([6,6,6],[4,2],3)
        t^6 + t^5 + 2*t^4 + t^3 + t^2
        sage: spin_polynomial([6,6,6],[4,1,1],3)
        t^6 + 2*t^5 + 3*t^4 + 2*t^3 + t^2
        sage: spin_polynomial([3,3,3,2,1], [2,2], 3)
        t^(7/2) + t^(5/2)
        sage: spin_polynomial([3,3,3,2,1], [2,1,1], 3)
        2*t^(7/2) + 2*t^(5/2) + t^(3/2)
        sage: spin_polynomial([3,3,3,2,1], [1,1,1,1], 3)
        3*t^(7/2) + 5*t^(5/2) + 3*t^(3/2) + sqrt(t)
        sage: spin_polynomial([5,4,3,2,1,1,1], [2,2,1], 3)
        2*t^(9/2) + 6*t^(7/2) + 2*t^(5/2)
        sage: spin_polynomial([[6]*6, [3,3]], [4,4,2], 3)
        3*t^9 + 5*t^8 + 9*t^7 + 6*t^6 + 3*t^5
    """
    sp = spin_polynomial_square(part,weight,length)
    t = sage.calculus.calculus.var('t')
    c = sp.coeffs()
    return sum([c[i]*t**(QQ(i)/2) for i in range(len(c))])


def cospin_polynomial(part, weight, length):
    """
    Returns the cospin polynomial associated to part, weight, and
    length.

    EXAMPLES::

        sage: from sage.combinat.ribbon_tableau import cospin_polynomial
        sage: cospin_polynomial([6,6,6],[4,2],3)
        t^4 + t^3 + 2*t^2 + t + 1
        sage: cospin_polynomial([3,3,3,2,1], [3,1], 3)
        1
        sage: cospin_polynomial([3,3,3,2,1], [2,2], 3)
        t + 1
        sage: cospin_polynomial([3,3,3,2,1], [2,1,1], 3)
        t^2 + 2*t + 2
        sage: cospin_polynomial([3,3,3,2,1], [1,1,1,1], 3)
        t^3 + 3*t^2 + 5*t + 3
        sage: cospin_polynomial([5,4,3,2,1,1,1], [2,2,1], 3)
        2*t^2 + 6*t + 2
        sage: cospin_polynomial([[6]*6, [3,3]], [4,4,2], 3)
        3*t^4 + 6*t^3 + 9*t^2 + 5*t + 3
    """
    R = ZZ['t']
    t = R.gen()

    #The power in the spin polynomial are all half integers
    #or all integers.  Manipulation of expressions need to
    #seperate cases
    sp = spin_polynomial_square(part, weight, length)
    if sp == 0:
        return R(0)

    coeffs = [c for c in sp.coeffs() if c != 0]
    d = len(coeffs)-1
    exponents = [d-e for e in range(len(coeffs))]


    return R(sum([ coeffs[i]*t**exponents[i] for i in range(len(coeffs))]))




##     //////////////////////////////////////////////////////////////////////////////////////////
##     // Generic function for driving into the graph of partitions coding all ribbons
##     // tableaux of a given shape and weight
##     //////////////////////////////////////////////////////////////////////////////////////////

##     //This function construct the graph of the set of k-ribbon tableaux
##     //of a given skew shape and a given weight.
##     //The first argument is always a skew partition.
##     //In the case where the inner partition is empty there is no branch without solutions
##     //In the other cases there is in average a lot of branches without solutions
##     /////////////////////////////////////////////////////////////////////////////////////////



def graph_implementation_rec(skp, weight, length, function):
    """
    TESTS::

        sage: from sage.combinat.ribbon_tableau import graph_implementation_rec, list_rec
        sage: graph_implementation_rec(SkewPartition([[1], []]), [1], 1, list_rec)
        [[[], [[1]]]]
        sage: graph_implementation_rec(SkewPartition([[2, 1], []]), [1, 2], 1, list_rec)
        [[[], [[2], [1, 2]]]]
        sage: graph_implementation_rec(SkewPartition([[], []]), [0], 1, list_rec)
        [[[], []]]
    """
    inn = "graph_implementation_rec(%s, %s, %s, %s)"%(skp, weight, length, function)
    #print "!!!", inn
    if sum(weight) == 0:
        weight = []

    partp = skp[0].conjugate()

    ## Some tests in order to know if the shape and the weight are compatible.
    if weight != [] and weight[-1] <= len(partp):
        perms = permutation.Permutations([0]*(len(partp)-weight[-1]) + [length]*(weight[-1])).list()
    else:
        return function([], [], skp, weight, length)

    selection = []

    for j in range(len(perms)):
        retire = [(partp[i]+ len(partp) - (i+1) - perms[j][i]) for i in range(len(partp))]
        retire.sort(reverse=True)
        retire = [ retire[i] - len(partp) + (i+1) for i in range(len(retire))]

        if retire[-1] >= 0 and retire == [i for i in reversed(sorted(retire))]:
            retire = partition.Partition_class(filter(lambda x: x != 0, retire)).conjugate()


            # Cutting branches if the retired partition has a line strictly included into the inner one
            append = True
            padded_retire = retire + [0]*(len(skp[1])-len(retire))
            for k in range(len(skp[1])):
                if padded_retire[k] - skp[1][k] < 0:
                    append = False
                    break
            if append:
                selection.append([retire, perms[j]])


    #selection contains the list of current nodes
    #print "selection", selection

    if len(weight) == 1:
        return function([], selection, skp, weight, length)
    else:
        #The recursive calls permit us to construct the list of the sons
        #of all current nodes in selection
        a = [graph_implementation_rec([p[0], skp[1]], weight[:-1], length, function) for p in selection]
        #print "a", a
        return function(a, selection, skp, weight, length)



##############################################################
def MultiSkewTableau(x):
    """
    Returns a multi skew tableau object which is a tuple of skew
    tableau.

    EXAMPLES::

        sage: s = MultiSkewTableau([ [[None,1],[2,3]], [[1,2],[2]] ])
        sage: s.size()
        6
        sage: s.weight()
        [2, 3, 1]
        sage: s.shape()
        [[[2, 2], [1]], [[2, 1], []]]
    """
    if isinstance(x, MultiSkewTableau_class):
        return x

    return MultiSkewTableau_class([skew_tableau.SkewTableau(i) for i in x] )

class MultiSkewTableau_class(CombinatorialObject):
    def size(self):
        """
        Returns the size of self, which is the sum of the sizes of the skew
        tableaux in self.

        EXAMPLES::

            sage: s = SemistandardSkewTableaux([[2,2],[1]]).list()
            sage: a = MultiSkewTableau([s[0],s[1],s[2]])
            sage: a.size()
            9
        """
        return sum(x.size() for x in self)

    def weight(self):
        """
        Returns the weight of self.

        EXAMPLES::

            sage: s = SemistandardSkewTableaux([[2,2],[1]]).list()
            sage: a = MultiSkewTableau([s[0],s[1],s[2]])
            sage: a.weight()
            [5, 3, 1]
        """
        weights = [x.weight() for x in self]
        m = max([len(x) for x in weights])
        weight = [0]*m
        for w in weights:
            for i in range(len(w)):
                weight[i] += w[i]
        return weight

    def shape(self):
        """
        Returns the shape of self.

        EXAMPLES::

            sage: s = SemistandardSkewTableaux([[2,2],[1]]).list()
            sage: a = MultiSkewTableau([s[0],s[1],s[2]])
            sage: a.shape()
            [[[2, 2], [1]], [[2, 2], [1]], [[2, 2], [1]]]
        """
        return [x.shape() for x in self]

    def inversion_pairs(self):
        """
        Returns a list of the inversion pairs of self.

        EXAMPLES::

            sage: s = MultiSkewTableau([ [[2,3],[5,5]], [[1,1],[3,3]], [[2],[6]] ])
            sage: s.inversion_pairs()
            [((0, (0, 0)), (1, (0, 0))),
             ((0, (1, 0)), (1, (0, 1))),
             ((0, (1, 1)), (1, (0, 0))),
             ((0, (1, 1)), (1, (1, 1))),
             ((0, (1, 1)), (2, (0, 0))),
             ((1, (0, 1)), (2, (0, 0))),
             ((1, (1, 1)), (2, (0, 0)))]
        """
        inv = []
        for k in range(len(self)):
            for b in self[k].boxes():
                inv += self._inversion_pairs_from_position(k,b)
        return inv

    def inversions(self):
        """
        Returns the number of inversion pairs of self.

        EXAMPLES::

            sage: t1 = SkewTableau([[1]])
            sage: t2 = SkewTableau([[2]])
            sage: MultiSkewTableau([t1,t1]).inversions()
            0
            sage: MultiSkewTableau([t1,t2]).inversions()
            0
            sage: MultiSkewTableau([t2,t2]).inversions()
            0
            sage: MultiSkewTableau([t2,t1]).inversions()
            1
            sage: s = MultiSkewTableau([ [[2,3],[5,5]], [[1,1],[3,3]], [[2],[6]] ])
            sage: s.inversions()
            7
        """
        return len(self.inversion_pairs())

    def _inversion_pairs_from_position(self, k, ij):
        """
        Returns the number of inversions at the box position i,j in the kth
        tableaux in self.

        EXAMPLES::

            sage: s = MultiSkewTableau([ [[2,3],[5,5]], [[1,1],[3,3]], [[2],[6]] ])
            sage: s._inversion_pairs_from_position(0, (1,1))
            [((0, (1, 1)), (1, (0, 0))),
             ((0, (1, 1)), (1, (1, 1))),
             ((0, (1, 1)), (2, (0, 0)))]
            sage: s._inversion_pairs_from_position(1, (0,1))
            [((1, (0, 1)), (2, (0, 0)))]
        """
        pk = k
        pi,pj = ij
        c = pi - pj
        value = self[pk][pi][pj]
        pk_boxes = self[pk].boxes_by_content(c)
        same_diagonal  = [ t.boxes_by_content(c) for t in self[pk+1:] ]
        above_diagonal = [ t.boxes_by_content(c+1) for t in self[pk+1:] ]

        res = []
        for i,j in pk_boxes:
            if pi < i and value > self[pk][i][j]:
                res.append( ((pk,(pi,pj)), (pk,(i,j))) )
        for k in range(len(same_diagonal)):
            for i,j in same_diagonal[k]:
                if value > self[pk+k+1][i][j]:
                    res.append( ((pk,(pi,pj)), (pk+k+1,(i,j))) )
        for k in range(len(above_diagonal)):
            for i,j in above_diagonal[k]:
                if value < self[pk+k+1][i][j]:
                    res.append( ((pk,(pi,pj)), (pk+k+1,(i,j))) )
        return res


def SemistandardMultiSkewTableaux(shape, weight):
    """
    Returns the combinatorial class of semistandard multi skew
    tableaux. A multi skew tableau is a k-tuple of skew tableaux of
    givens shape with a specified total weight.

    EXAMPLES::

        sage: s = SemistandardMultiSkewTableaux([ [[2,1],[]], [[2,2],[1]] ], [2,2,2]); s
        Semistandard multi skew tableaux of shape [[[2, 1], []], [[2, 2], [1]]] and weight [2, 2, 2]
        sage: s.list()
        [[[[1, 1], [2]], [[None, 2], [3, 3]]],
         [[[1, 2], [2]], [[None, 1], [3, 3]]],
         [[[1, 3], [2]], [[None, 2], [1, 3]]],
         [[[1, 3], [2]], [[None, 1], [2, 3]]],
         [[[1, 1], [3]], [[None, 2], [2, 3]]],
         [[[1, 2], [3]], [[None, 2], [1, 3]]],
         [[[1, 2], [3]], [[None, 1], [2, 3]]],
         [[[2, 2], [3]], [[None, 1], [1, 3]]],
         [[[1, 3], [3]], [[None, 1], [2, 2]]],
         [[[2, 3], [3]], [[None, 1], [1, 2]]]]
    """
    shape = [skew_partition.SkewPartition(x) for x in shape]
    weight = partition.Partition(weight)

    if sum(weight) != sum(s.size() for s in shape):
        raise ValueError, "the sum of weight must be the sum of the sizes of shape"

    return SemistandardMultiSkewTtableaux_shapeweight(shape, weight)

class SemistandardMultiSkewTtableaux_shapeweight(CombinatorialClass):
    def __init__(self, shape, weight):
        """
        TESTS::

            sage: s = SemistandardMultiSkewTableaux([ [[2,1],[]], [[2,2],[1]] ], [2,2,2])
            sage: s == loads(dumps(s))
            True
        """
        self._shape  = shape
        self._weight = weight
        self._name   = "Semistandard multi skew tableaux of shape %s and weight %s"%(shape, weight)

    def __contains__(self, x):
        """
        TESTS::

            sage: s = SemistandardMultiSkewTableaux([ [[2,1],[]], [[2,2],[1]] ], [2,2,2])
            sage: all(i in s for i in s)
            True
        """
        try:
            x = MultiSkewTableau(x)
        except TypeError:
            return False

        if x.weight() != self._weight:
            return False

        if x.shape() != self._shape:
            return False

        if not all( x[i].is_semistandard() for i in range(len(x)) ):
            return False

        return True

    def list(self):
        """
        EXAMPLES::

            sage: sp = SkewPartitions(3).list()
            sage: SemistandardMultiSkewTableaux([sp[0], sp[-1]],[2,2,2]).list()
            [[[[1], [2], [3]], [[1, 2, 3]]]]

        ::

            sage: a = SkewPartition([[8,7,6,5,1,1],[2,1,1]])
            sage: weight = [3,3,2]
            sage: k = 3
            sage: s = SemistandardMultiSkewTableaux(a.r_quotient(k),weight)
            sage: len(s.list())
            34
            sage: RibbonTableaux(a,weight,k).cardinality()
            34
        """
        parts = self._shape
        mu = self._weight

        #Splitting the partition
        s = [ p.size() for p in parts ]
        parts = [p.to_list() for p in parts]

        #Gluing the partitions
        parttmp = parts[0]
        i = 1
        for i in range(1,len(parts)):
            trans = parttmp[0][0]
            current_part = parts[i]
            current_part[1] += [0]*(len(current_part[0])-len(current_part[1]))
            inner_current = [ trans + j for j in current_part[1] ]
            outer_current = [ trans + j for j in current_part[0] ]
            parttmp = [ outer_current + parttmp[0], inner_current + parttmp[1] ]

        #List the corresponding skew tableaux
        l = [ st.to_word() for st in skew_tableau.SemistandardSkewTableaux(parttmp, mu) ]

        res = []
        for k in range(len(l)):
            pos = 0  #Double check this
            restmp = [ skew_tableau.from_shape_and_word(parts[0], [l[k][j] for j in range(s[0])]) ]
            for i in range(1, len(parts)):
                w = [l[k][j] for j in range(pos+s[i-1], pos+s[i-1]+s[i])]
                restmp.append( skew_tableau.from_shape_and_word(parts[i], w) )
            res.append(MultiSkewTableau_class(restmp))
        return res

