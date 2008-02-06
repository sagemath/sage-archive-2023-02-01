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
import skew_tableau, tableau, word, permutation, partition
from sage.rings.all import QQ, ZZ
import skew_partition
import sage.calculus.calculus
import __builtin__



class RibbonTableau_class(CombinatorialObject):

    def length(self):
        """
        Return the length of the ribbons into a ribbon tableaux.
        """
        if self.to_expr() == [[],[]]:
            return 0

        tableau = self.to_expr()[1]
        l = 0
        t = 0
        for k in range(len(tableau)):
            t += len( [ x for x in tableau[k] if x > -1 ] )
            l += len( [ x for x in tableau[k] if x > 0  ] )

        if l == 0:
            return t
        else:
            return t/l




#####################
# Ribbon Tableaux   #
#####################

def RibbonTableaux(shape, weight, length):
    """
    Returns the combinatorial class of ribbon tableaux of skew shape
    shape and weight weight tiled by ribbons of length length.

    EXAMPLES:
        sage: RibbonTableaux([[2,1],[]],[1,1,1],1)
        Ribbon tableaux of shape [[2, 1], []] and weight [1, 1, 1] with 1-ribbons

    """
    if shape in partition.Partitions():
        shape = partition.Partition(shape)
        shape = skew_partition.SkewPartition([shape, shape.rcore(length)])
    else:
        shape = skew_partition.SkewPartition(shape)

    if shape.size() != length*sum(weight):
        raise ValueError

    weight = [i for i in weight if i != 0]

    return RibbonTableaux_shapeweightlength(shape, weight,length)


class RibbonTableaux_shapeweightlength(CombinatorialClass):
    def __init__(self, shape, weight, length):
        """
        EXAMPLES:
            sage: r = RibbonTableaux([[2,1],[]],[1,1,1],1)
            sage: r == loads(dumps(r))
            True
        """
        self._shape = shape
        self._weight = weight
        self._length = length
        self._name = "Ribbon tableaux of shape %s and weight %s with %s-ribbons"%(shape, weight, length)

    def list(self):
        """
        EXAMPLES:
            sage: RibbonTableaux([[2,1],[]],[1,1,1],1).list()
            [[[], [[2], [1, 3]]], [[], [[3], [1, 2]]]]
        """
        return graph_implementation_rec(self._shape, self._weight, self._length, list_rec)


def insertion_tableau(skp, perm, evaluation, tableau, length):
    """
    INPUT:
        skp -- skew partitions
        perm, evaluation -- non-negative integers
        tableau -- skew tableau
        length -- integer

    TESTS:
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




def list_rec(nexts, current, part, weight, length):
    """
    INPUT:
        nexts, current, part -- skew partitions
        weight -- non-negative integer list
        length -- integer

    TESTS:
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
def spin_rec(nexts, current, part, weight, length, t):
    """
    Routine used for constructing the spin polynomial.

    INPUT:
        nexts, current, part -- skew partitions
        weight -- list of non-negative integers
        length -- the length of the ribbons we're tiling with
        t      -- the variable

    EXAMPLES:
        sage: from sage.combinat.ribbon_tableau import spin_rec
        sage: t = ZZ['t'].gen()
        sage: spin_rec([], [[[], [3, 3]]], [[2, 2, 2], []], [2], 3, t)
        [t^4]
        sage: spin_rec([[0], [t^4]], [[[2, 1, 1, 1, 1], [0, 3]], [[2, 2, 2], [3, 0]]], [[2, 2, 2, 2, 1], []], [2, 1], 3, t)
        [t^5]
        sage: spin_rec([], [[[], [3, 3, 0]]], [[3, 3], []], [2], 3, t)
        [t^2]
        sage: spin_rec([[0], [t^2]], [[[2, 1, 1, 1, 1], [0, 3, 0]], [[3, 3], [3, 0, 0]]], [[3, 3, 1, 1, 1], []], [2, 1], 3, t)
        [t^4]
        sage: spin_rec([], [[[], [3, 3]]], [[2, 2, 2], []], [2], 3, t)
        [t^4]
        sage: spin_rec([], [[[], [3, 3, 0]]], [[3, 2, 1], []], [2], 3, t)
        [t^3]
        sage: spin_rec([], [[[], [3, 3, 0]]], [[3, 3], []], [2], 3, t)
        [t^2]
        sage: spin_rec([[t^4], [t^3], [t^2]], [[[2, 2, 2], [0, 0, 3]], [[3, 2, 1], [0, 3, 0]], [[3, 3], [3, 0, 0]]], [[3, 3, 3], []], [2, 1], 3, t)
        [t^6 + t^4 + t^2]
        sage: spin_rec([[t^5], [t^4], [t^6 + t^4 + t^2]], [[[2, 2, 2, 2, 1], [0, 0, 3]], [[3, 3, 1, 1, 1], [0, 3, 0]], [[3, 3, 3], [3, 0, 0]]], [[3, 3, 3, 2, 1], []], [2, 1, 1], 3, t)
        [2*t^7 + 2*t^5 + t^3]

    """
    R = ZZ['t']
    if current == []:
        return [R(0)]

    res = R(0)
    tmp = []
    partp = partition.Partition_class(part[0]).conjugate()

    #compute the contribution of the ribbons added at
    #the current node
    for perms in [current[i][1] for i in range(len(current))]:
        perm =  [partp[i] + len(partp)-(i+1)-perms[i] for i in range(len(partp))]
        perm.reverse()
        perm = permutation.Permutation( word.standard(perm) )
        tmp.append( (weight[-1]*(length-1)-perm.number_of_inversions()) )

    if nexts != []:
        return [sum([sum([t**tmp[i]*nexts[i][j] for j in range(len(nexts[i]))]) for i in range(len(tmp))])]
    else:
        return [sum([t**tmp[i] for i in range(len(tmp))])]


def spin_polynomial_square(part, weight, length):
    """
    Returns the spin polynomial associated with part, weight, and
    length, with the substitution t -> t^2 made.

    EXAMPLES:
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
        part = [part,[]]

    if part == [[],[]] and weight == []:
        return t.parent()(1)

    def f(*args):
        return spin_rec(*(args+(t,)))

    return R(graph_implementation_rec(part, weight, length, f)[0])

def spin_polynomial(part,weight,length):
    """
    Returns the spin polynomial associated to part, weight, and
    length.

    EXAMPLES:
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
    Returns the cospin polynomial associated to part,
    weight, and length.

    EXAMPLES:
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
##     //The first argument is alaways a skew partition.
##     //In the case where the inner partition is empty there is no branch without solutions
##     //In the other cases there is in average a lot of branches without solutions
##     /////////////////////////////////////////////////////////////////////////////////////////



def graph_implementation_rec(skp, weight, length, function):
    """
    TESTS:
        sage: from sage.combinat.ribbon_tableau import graph_implementation_rec, list_rec
        sage: graph_implementation_rec([[1], []], [1], 1, list_rec)
        [[[], [[1]]]]
        sage: graph_implementation_rec([[2, 1], []], [1, 2], 1, list_rec)
        [[[], [[2], [1, 2]]]]
        sage: graph_implementation_rec([[], []], [0], 1, list_rec)
        [[[], []]]

    """
    inn = "graph_implementation_rec(%s, %s, %s, %s)"%(skp, weight, length, function)
    #print "!!!", inn
    if sum(weight) == 0:
        weight = []

    partp = partition.Partition_class(skp[0]).conjugate()

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


class TupleTableaux(CombinatorialClass):
    pass

def list_ktupletableaux(parts, mu):
    """
    INPUT:
        parts -- list of skew partitions
        mu    -- partition
    """
    #Gluing the partitions
    parttmp = parts[0]
    i = 1
    for i in range(1,len(parts)):
        trans = parttmp[0][0]
        current_part = parts[i].to_list()
        current_part[1] += [0]*(len(current_part[0])-len(current_part[1]))
        inner_current = [ trans + j for j in current_part[1] ]
        outer_current = [ trans + j for j in current_part[0] ]
        parttmp = [ outer_current + parttmp[0], inner_current + parttmp[1] ]

    #List the corresponding skew tableaux
    l = [ st.to_word() for st in SkewTableaux(parttmp, mu) ]

    #Splitting the partition
    s = [ p.size() for p in parts ]

    res = []
    for k in range(len(l)):
        pos = 0  #Double check this
        restmp = [ skew_tableau.from_shape_and_word(parts[0], [l[k][j] for j in range(s[0])]) ]
        for i in range(1, len(parts)):
            w = [l[k][j] for j in range(pos+s[i-1], pos+s[i-1]+s[i])]
            restmp.append( skew_tableau.from_shape_and_word(parts[i], w) )
        res.append(restmp)
    return res

def inversions(tabs):
    """
    Returns the number of inverion pairs of tabs.

    INPUT:
        tabs -- a list of tableaux or skew tableaux

    EXAMPLES:
        sage: from sage.combinat.ribbon_tableau import inversions
        sage: t1 = Tableau([[1]])
        sage: t2 = Tableau([[2]])
        sage: inversions([t1,t1])
        0
        sage: inversions([t1,t2])
        0
        sage: inversions([t2,t2])
        0
        sage: inversions([t2,t1])
        1
    """
    inv = 0
    for k in range(len(tabs)):
        for i in range(len(tabs[k])):
            for j in range(len( tabs[k][-i-1] ) ):
                inv += inversions_from_position(tabs, [k,[i,j]])
    return inv

def inversions_from_position(tabs, pos):
    """
    INPUT:
        tabs -- a list of tableaux or skew tableaux
    """
    pk, (pi, pj) = pos

    #Add checking wether tabs is a list of tableau or
    #skew tableaux
    pass

    #Building the list of inner partitions

    if tabs[0] in tableau.SemistandardTableaux():
        inners = [[0]*len(t) for t in tabs]

        #There is no alignment in the case of straight partitions
        trans = [0]*len(tabs)

    else:
        #tabs = expr(tabs)
        inners = [ t.to_list()[0] for t in tabs ]
        for k in range(len(tabs)):
            pass

        #Computation of the alignment of skew shapes at the
        #southeast / northeast

        skp    = [ t.shape()  for t  in tabs ]
        outers = [ sp.outer() for sp in skp ]
        ma     = max( [ x[0] for x in outers if len(x) > 0 ] )

        #The translation vector
        trans  = [ ma - outers[i][0] for i in range(len(outers)) ]

        #Extracting just the fillings of the skew tableaux
        tabs = [ x[1] for x in tabs ]

    inv = 0

    #Computation of the inversions
    for k in range(len(tabs)):
        for i in range(len(tabs[k])):
            for j in range(len(tabs[k][-i-1])):
                if tabs[pk][-pi-1][pj] < tabs[k][-i-1][j]:
                    #First condition diag(pos) = diag(current) and k < pk
                    #Second condition diag(pos)-1 = diag(current) and k > pk
                    if ((k<pk) and inners[k][i]+trans[k]+j-i == inners[pk][pi]+trans[pk]+pj-pi) or \
                       ((k>pk) and inners[k][i]+trans[k]+j-i == inners[pk][pi]+trans[pk]+pj-pi-1):
                        inv += 1
    return inv

