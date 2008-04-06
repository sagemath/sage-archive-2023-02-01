r"""
Permutations

The Permutations module. Use Permutation? to get information about the
Permutation class, and Permutations? to get information about the
combinatorial class of permutations.
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


from sage.interfaces.all import gap, maxima
from sage.rings.all import QQ, RR, ZZ, polygen, Integer, PolynomialRing, factorial
from sage.rings.arith import binomial
from sage.misc.sage_eval import sage_eval
from sage.libs.all import pari
from sage.matrix.all import matrix
from sage.combinat.tools import transitive_ideal
import sage.combinat.misc as misc
import sage.combinat.subword as subword
import sage.combinat.composition as composition
from sage.combinat.composition import Composition, Compositions, Composition_class
import tableau
import sage.combinat.partition
from permutation_nk import PermutationsNK
import sage.rings.integer
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.misc.prandom import randint, sample
from sage.interfaces.all import gap
from sage.graphs.graph import DiGraph
import itertools
import __builtin__
from combinat import CombinatorialClass, CombinatorialObject, catalan_number
import copy
from necklace import Necklaces
import tableau

permutation_options = {'display':'list', 'mult':'l2r'}

def PermutationOptions(**kwargs):
    """
    Sets the global options for elements of the permutation class.
    The defaults are for permutations to be displayed in list notation
    and the multiplication done from left to right (like in GAP).

    display: 'list'  -- the permutations are displayed in list notation
             'cycle' -- the permutations are displayed in cycle notation
             'singleton' -- the permutations are displayed in cycle notation
                            with singleton cycles shown as well.

    mult: 'l2r' -- the multiplication of permutations is done like composition
                   of functions from left to right. That is, if we think
                   of the permutations p1 and p2 as functions, then
                   (p1*p2)(x) = p2(p1(x)). This is the default in multiplication
                   in GAP.
          'r2l' -- the multiplication of permutations is done right to left
                   so that (p1*p2)(x) = p1(p2(x))

    If no parameters are set, then the function returns a copy of the
    options dictionary.

    Note that these options have no effect on PermutationGroupElements.

    EXAMPLES:
        sage: p213 = Permutation([2,1,3])
        sage: p312 = Permutation([3,1,2])
        sage: PermutationOptions(mult='l2r', display='list')
        sage: po = PermutationOptions()
        sage: po['display']
        'list'
        sage: p213
        [2, 1, 3]
        sage: PermutationOptions(display='cycle')
        sage: p213
        (1,2)
        sage: PermutationOptions(display='singleton')
        sage: p213
        (1,2)(3)
        sage: PermutationOptions(display='list')

        sage: po['mult']
        'l2r'
        sage: p213*p312
        [1, 3, 2]
        sage: PermutationOptions(mult='r2l')
        sage: p213*p312
        [3, 2, 1]
        sage: PermutationOptions(mult='l2r')
    """
    global permutation_options
    if kwargs == {}:
        return copy.copy(permutation_options)

    if 'mult' in kwargs:
        if kwargs['mult'] not in ['l2r', 'r2l']:
            raise ValueError, "mult must be either 'l2r' or 'r2l'"
        else:
            permutation_options['mult'] = kwargs['mult']

    if 'display' in kwargs:
        if kwargs['display'] not in ['list', 'cycle', 'singleton']:
            raise ValueError, "display must be either 'cycle' or 'list'"
        else:
            permutation_options['display'] = kwargs['display']


def Permutation(l):
    """
    Convert l to a Permutation, where l is a list, tuple of integers,
    tuple of tuples of integers, or a string in cycle notation. Returns
    a member of the Permutation class, printed in one-line notation.

    EXAMPLES:
        sage: Permutation([2,1])
        [2, 1]
        sage: Permutation([2, 1, 4, 5, 3])
        [2, 1, 4, 5, 3]
        sage: Permutation('(1,2)')
        [2, 1]
        sage: Permutation('(1,2)(3,4,5)')
        [2, 1, 4, 5, 3]
        sage: Permutation( ((1,2),(3,4,5)) )
        [2, 1, 4, 5, 3]
        sage: Permutation( ((1,2)) )
        [2, 1]
        sage: Permutation( (1,2) )
        [2, 1]
        sage: Permutation( ((1,2),) )
        [2, 1]
        sage: p = Permutation((1, 2, 5)); p
        [2, 5, 3, 4, 1]
        sage: type(p)
        <class 'sage.combinat.permutation.Permutation_class'>


    """
    if isinstance(l, Permutation_class):
        return l
    #if l is a string, then assume it is in cycle notation
    elif isinstance(l, str):
        cycles = l.split(")(")
        cycles[0] = cycles[0][1:]
        cycles[-1] = cycles[-1][:-1]
        cycle_list = []
        for c in cycles:
            cycle_list.append(map(int, c.split(",")))
        return from_cycles(sum([len(c) for c in cycle_list]), cycle_list)
    elif isinstance(l, tuple):
        if len(l) >= 1:
            if isinstance(l[0], tuple):
                n = max( map(max, l) )
                return from_cycles(n, map(list, l))
            else:
                n = max(l)
                l = [list(l)]
                return from_cycles(n, l)
        else:
            raise ValueError, "cannot convert l (= %s) to a Permutation"%l
    else:
        return Permutation_class(l)

class Permutation_class(CombinatorialObject):
    def __hash__(self):
        """
        TESTS:
            sage: d = {}
            sage: p = Permutation([1,2,3])
            sage: d[p] = 1
            sage: d[p]
            1
        """
        if self._hash is None:
            self._hash = str(self).__hash__()
        return self._hash

    def __str__(self):
        """
        TESTS:
            sage: PermutationOptions(display='list')
            sage: p = Permutation([2,1,3])
            sage: str(p)
            '[2, 1, 3]'
            sage: PermutationOptions(display='cycle')
            sage: str(p)
            '(1,2)'
            sage: PermutationOptions(display='singleton')
            sage: str(p)
            '(1,2)(3)'
            sage: PermutationOptions(display='list')
        """
        return repr(self)

    def __repr__(self):
        """
        TESTS:
            sage: PermutationOptions(display='list')
            sage: p = Permutation([2,1,3])
            sage: repr(p)
            '[2, 1, 3]'
            sage: PermutationOptions(display='cycle')
            sage: repr(p)
            '(1,2)'
            sage: PermutationOptions(display='singleton')
            sage: repr(p)
            '(1,2)(3)'
            sage: PermutationOptions(display='list')
        """
        global permutation_options
        display = permutation_options['display']
        if display == 'list':
            return repr(self._list)
        elif display == 'cycle':
            return self.cycle_string()
        elif display == 'singleton':
            return self.cycle_string(singletons=True)
        else:
            raise ValueError, "permutation_options['display'] should be one of 'list', 'cycle', or 'singleton'"

    def cycle_string(self, singletons=False):
        """
        Returns a string of the permutation in cycle notation.

        If singletons=True, it includes 1-cycles in the string.

        EXAMPLES:
            sage: Permutation([1,2,3]).cycle_string()
            '()'
            sage: Permutation([2,1,3]).cycle_string()
            '(1,2)'
            sage: Permutation([2,3,1]).cycle_string()
            '(1,2,3)'
            sage: Permutation([2,1,3]).cycle_string(singletons=True)
            '(1,2)(3)'
        """
        cycles = self.to_cycles(singletons=singletons)
        if cycles == []:
            return "()"
        else:
            return "".join(["("+",".join([str(l) for l in x])+")" for x in cycles])

    def next(self):
        r"""
        Returns the permutation that follows p in lexicographic order.
        If p is the last permutation, then next returns false.

        EXAMPLES:
            sage: p = Permutation([1, 3, 2])
            sage: p.next()
            [2, 1, 3]
            sage: p = Permutation([4,3,2,1])
            sage: p.next()
            False
        """
        p = self[:]
        n = len(self)
        first = -1

        #Starting from the end, find the first o such that
        #p[o] < p[o+1]
        for i in reversed(range(0,n-1)):
            if p[i] < p[i+1]:
                first = i
                break

        #If first is still -1, then we are already at the last permutation
        if first == -1:
            return False

        #Starting from the end, find the first j such that p[j] > p[first]
        j = n - 1
        while p[j] < p[first]:
            j -= 1

        #Swap positions first and j
        (p[j], p[first]) = (p[first], p[j])

        #Reverse the list between first and the end
        first_half = p[:first+1]
        last_half  = p[first+1:]
        last_half.reverse()
        p = first_half + last_half

        return Permutation(p)

    def prev(self):
        r"""
        Returns the permutation that comes directly before p
        in lexicographic order.  If p is the first permutation, then
        it returns False.

        EXAMPLES:
            sage: p = Permutation([1,2,3])
            sage: p.prev()
            False
            sage: p = Permutation([1,3,2])
            sage: p.prev()
            [1, 2, 3]
        """

        p = self[:]
        n = len(self)
        first = -1

        #Starting from the beginning, find the first o such that
        #p[o] > p[o+1]
        for i in range(0, n-1):
            if p[i] > p[i+1]:
                first = i
                break

        #If first is still -1, that is we didn't find any descents,
        #then we are already at the last permutation
        if first == -1:
            return False

        #Starting from the end, find the first j such that p[j] > p[first]
        j = n - 1
        while p[j] > p[first]:
            j -= 1

        #Swap positions first and j
        (p[j], p[first]) = (p[first], p[j])

        #Reverse the list between first+1 and end
        first_half = p[:first+1]
        last_half  = p[first+1:]
        last_half.reverse()
        p = first_half + last_half

        return Permutation(p)


    def to_tableau_by_shape(self, shape):
        """
        Returns a tableau of shape shape with the entries in self.

        EXAMPLES:
            sage: Permutation([3,4,1,2,5]).to_tableau_by_shape([3,2])
            [[1, 2, 5], [3, 4]]
            sage: Permutation([3,4,1,2,5]).to_tableau_by_shape([3,2]).to_permutation()
            [3, 4, 1, 2, 5]
        """
        if sum(shape) != len(self):
            raise ValueError, "the size of the partition must be the length of self"

        t = []
        w = list(self)
        for i in reversed(shape):
            t = [ w[:i] ] + t
            w = w[i:]
        return tableau.Tableau(t)



    def to_cycles(self, singletons=True):
        r"""
        Returns the permutation p as a list of disjoint cycles.

        EXAMPLES:
            sage: Permutation([2,1,3,4]).to_cycles()
            [(1, 2), (3,), (4,)]
            sage: Permutation([2,1,3,4]).to_cycles(singletons=False)
            [(1, 2)]
        """
        p = self[:]
        cycles = []
        toConsider = -1

        #Create the list [1,2,...,len(p)]
        l = [ i+1 for i in range(len(p))]
        cycle = []

        #Go through until we've considered every number between
        #1 and len(p)
        while len(l) > 0:
            #If we are at the end of a cycle
            #then we want to add it to the cycles list
            if toConsider == -1:
                #Add the cycle to the list of cycles
                if singletons:
                    if cycle != []:
                        cycles.append(tuple(cycle))
                else:
                    if len(cycle) > 1:
                        cycles.append(tuple(cycle))


                #Start with the first element in the list
                toConsider = l[0]
                l.remove(toConsider)
                cycle = [ toConsider ]
                cycleFirst = toConsider

            #Figure out where the element under consideration
            #gets mapped to.
            next = p[toConsider - 1]

            #If the next element is the first one in the list
            #then we've reached the end of the cycle
            if next == cycleFirst:
                toConsider = -1
            else:
                cycle.append( next )
                l.remove( next )
                toConsider = next

        #When we're finished, add the last cycle
        if singletons:
            if cycle != []:
                cycles.append(tuple(cycle))
        else:
            if len(cycle) > 1:
                cycles.append(tuple(cycle))

        return cycles

    def to_permutation_group_element(self):
        """
        Returns a PermutationGroupElement equal to self.

        EXAMPLES:
            sage: Permutation([2,1,4,3]).to_permutation_group_element()
            (1,2)(3,4)
            sage: Permutation([1,2,3]).to_permutation_group_element()
            ()
        """
        cycles = self.to_cycles(singletons=False)
        grp = SymmetricGroup(len(self))
        if cycles == []:
            return PermutationGroupElement( '()', parent=grp )
        else:
            return PermutationGroupElement( cycles , parent=grp)

    def signature(p):
        r"""
        Returns the signature of a permutation.

        EXAMPLES:
            sage: Permutation([4, 2, 3, 1, 5]).signature()
            -1

        """
        return (-1)**(len(p)-len(p.to_cycles()))


    def is_even(self):
        r"""
        Returns True if the permutation p is even and false otherwise.

        EXAMPLES:
            sage: Permutation([1,2,3]).is_even()
            True
            sage: Permutation([2,1,3]).is_even()
            False
        """

        if self.signature() == 1:
            return True
        else:
            return False


    def to_matrix(self):
        r"""
        Returns a matrix representing the permutation.

        EXAMPLES:
            sage: Permutation([1,2,3]).to_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: Permutation([1,3,2]).to_matrix()
            [1 0 0]
            [0 0 1]
            [0 1 0]

          Notice that matrix multiplication corresponds to permutation
          multiplication only when the permutation option mult='r2l'

            sage: PermutationOptions(mult='r2l')
            sage: p = Permutation([2,1,3])
            sage: q = Permutation([3,1,2])
            sage: (p*q).to_matrix()
            [0 0 1]
            [0 1 0]
            [1 0 0]
            sage: p.to_matrix()*q.to_matrix()
            [0 0 1]
            [0 1 0]
            [1 0 0]
            sage: PermutationOptions(mult='l2r')
            sage: (p*q).to_matrix()
            [1 0 0]
            [0 0 1]
            [0 1 0]
        """
        p = self[:]
        n = len(p)

        #Build the dictionary of entries since the matrix
        #is extremely sparse
        entries = {}
        for i in range(n):
            entries[(p[i]-1,i)] = 1
        return matrix(n, entries, sparse = True)

    def __mul__(self, rp):
        """
        TESTS:
            sage: p213 = Permutation([2,1,3])
            sage: p312 = Permutation([3,1,2])
            sage: PermutationOptions(mult='l2r')
            sage: p213*p312
            [1, 3, 2]
            sage: PermutationOptions(mult='r2l')
            sage: p213*p312
            [3, 2, 1]
            sage: PermutationOptions(mult='l2r')
        """
        global permutation_options
        if permutation_options['mult'] == 'l2r':
            return self._left_to_right_multiply_on_right(rp)
        else:
            return self._left_to_right_multiply_on_left(rp)

    def __rmul__(self, lp):
        """
        TESTS:
            sage: p213 = Permutation([2,1,3])
            sage: p312 = Permutation([3,1,2])
            sage: PermutationOptions(mult='l2r')
            sage: p213*p312
            [1, 3, 2]
            sage: PermutationOptions(mult='r2l')
            sage: p213*p312
            [3, 2, 1]
            sage: PermutationOptions(mult='l2r')
        """
        global permutation_options
        if permutation_options['mult'] == 'l2r':
            return self._left_to_right_multiply_on_left(lp)
        else:
            return self._left_to_right_multiply_on_right(lp)

    def _left_to_right_multiply_on_left(self,lp):
        """
        EXAMPLES:
            sage: p = Permutation([2,1,3])
            sage: q = Permutation([3,1,2])
            sage: p._left_to_right_multiply_on_left(q)
            [3, 2, 1]
            sage: q._left_to_right_multiply_on_left(p)
            [1, 3, 2]
        """
        #Pad the permutations if they are of
        #different lengths
        new_lp = lp[:] + [i+1 for i in range(len(lp), len(self))]
        new_p1 = self[:] + [i+1 for i in range(len(self), len(lp))]
        return Permutation([ new_p1[i-1] for i in new_lp ])


    def _left_to_right_multiply_on_right(self, rp):
        """
        EXAMPLES:
            sage: p = Permutation([2,1,3])
            sage: q = Permutation([3,1,2])
            sage: p._left_to_right_multiply_on_right(q)
            [1, 3, 2]
            sage: q._left_to_right_multiply_on_right(p)
            [3, 2, 1]
        """
        #Pad the permutations if they are of
        #different lengths
        new_rp = rp[:] + [i+1 for i in range(len(rp), len(self))]
        new_p1 = self[:] + [i+1 for i in range(len(self), len(rp))]
        return Permutation([ new_rp[i-1] for i in new_p1 ])


    ########
    # Rank #
    ########

    def rank(self):
        r"""
        Returns the rank of a permutation in lexicographic
        ordering.

        EXAMPLES:
            sage: Permutation([1,2,3]).rank()
            0
            sage: Permutation([1, 2, 4, 6, 3, 5]).rank()
            10
            sage: perms = Permutations(6).list()
            sage: [p.rank() for p in perms ] == range(factorial(6))
            True
        """
        n = len(self)

        factoradic = self.to_lehmer_code()

        #Compute the index
        rank = 0
        for i in reversed(range(0, n)):
            rank += factoradic[n-1-i]*factorial(i)

        return rank

    ##############
    # Inversions #
    ##############

    def to_inversion_vector(self):
        r"""
        Returns the inversion vector of a permutation p.

        If iv is the inversion vector, then iv[i] is the number of elements
        larger than i that appear to the left of i in the permutation.

        EXAMPLES:
            sage: Permutation([5,9,1,8,2,6,4,7,3]).to_inversion_vector()
            [2, 3, 6, 4, 0, 2, 2, 1, 0]
            sage: Permutation([8,7,2,1,9,4,6,5,10,3]).to_inversion_vector()
            [3, 2, 7, 3, 4, 3, 1, 0, 0, 0]
            sage: Permutation([3,2,4,1,5]).to_inversion_vector()
            [3, 1, 0, 0, 0]

        """
        p = self[:]
        inversion_vector = [0]*len(p)
        for i in range(len(p)):
            j = 0
            while p[j] != i+1:
                if p[j] > i+1:
                    inversion_vector[i] += 1
                j += 1

        return inversion_vector


    def inversions(self):
        r"""
        Returns a list of the inversions of permutation p.

        EXAMPLES:
            sage: Permutation([3,2,4,1,5]).inversions()
            [[0, 1], [0, 3], [1, 3], [2, 3]]


        """
        p = self[:]
        inversion_list = []

        for i in range(len(p)):
            for j in range(i+1,len(p)):
                if  p[i] > p[j]:
                    #inversion_list.append((p[i],p[j]))
                    inversion_list.append([i,j])

        return inversion_list



    def number_of_inversions(self):
        r"""
        Returns the number of inversions in the permutation p.

        An inversion of a permutation is a pair of elements (p[i],p[j])
        with i > j and p[i] < p[j].

        REFERENCES:
            http://mathworld.wolfram.com/PermutationInversion.html

        EXAMPLES:
            sage: Permutation([3,2,4,1,5]).number_of_inversions()
            4
            sage: Permutation([1, 2, 6, 4, 7, 3, 5]).number_of_inversions()
            6

        """

        return sum(self.to_inversion_vector())


    def length(self):
        r"""
        Returns the length of a permutation p.  The length is given by the
        number of inversions of p.

        EXAMPLES:
            sage: Permutation([5, 1, 3, 2, 4]).length()
            5
        """
        return self.number_of_inversions()

    def inverse(self):
        r"""
        Returns the inverse of a permutation

        EXAMPLES:
            sage: Permutation([3,8,5,10,9,4,6,1,7,2]).inverse()
            [8, 10, 1, 6, 3, 7, 9, 2, 5, 4]
            sage: Permutation([2, 4, 1, 5, 3]).inverse()
            [3, 1, 5, 2, 4]


        """
        return Permutation([self.index(i+1)+1 for i in range(len(self))])

    def _icondition(self, i):
        """
        Returns a string which shows the relative positions of i-1,i,i+1
        in self.  Note that i corresponds to a 2 in the string.

        NOTE: An imove can only be applied when the relative positions are
        one of '213', '132', '231', or '312'.  None is returned in the
        other cases to signal that an imove cannot be applied.

        EXAMPLES:
            sage: Permutation([2,1,3])._icondition(2)
            ('213', 1, 0, 2)
            sage: Permutation([1,3,2])._icondition(2)
            ('132', 0, 2, 1)
            sage: Permutation([2,3,1])._icondition(2)
            ('231', 2, 0, 1)
            sage: Permutation([3,1,2])._icondition(2)
            ('312', 1, 2, 0)
            sage: Permutation([1,2,3])._icondition(2)
            (None, 0, 1, 2)
            sage: Permutation([1,3,2,4])._icondition(3)
            ('213', 2, 1, 3)
            sage: Permutation([2,1,3])._icondition(3)
            Traceback (most recent call last):
            ...
            ValueError: i (= 3) must be between 2 and n-1
        """
        if i not in range(2, len(self)):
            raise ValueError, "i (= %s) must be between 2 and n-1"%i
        pos_i   = self.index(i)
        pos_ip1 = self.index(i+1)
        pos_im1 = self.index(i-1)

        if pos_i < pos_im1 and pos_im1 < pos_ip1:
            state = '213'
        elif pos_im1 < pos_ip1 and pos_ip1 < pos_i:
            state =  '132'
        elif pos_i < pos_ip1 and pos_ip1 < pos_im1:
            state =  '231'
        elif pos_ip1 < pos_im1 and pos_im1 < pos_i:
            state = '312'
        else:
            state = None

        return (state, pos_im1, pos_i, pos_ip1)

    def ishift(self, i):
        """
        Returns an the i-shift of self.  If an i-shift of self can't be
        performed, then None is returned.

        An i-shift can be applied when i is not in between i-1 and i+1.
        The i-shift moves i to the other side, and leaves the relative
        positions of i-1 and i+1 in place.

        EXAMPLES:
          Here, 2 is to the left of both 1 and 3.  A 2-shift can be
          applied which moves the 2 to the right and leaves 1 and 3
          in their same relative order.
            sage: Permutation([2,1,3]).ishift(2)
            [1, 3, 2]

          Note that the movement is done in place:
            sage: Permutation([2,4,1,3]).ishift(2)
            [1, 4, 3, 2]

          Since 2 is between 1 and 3 in [1,2,3], an 2-shift cannot be
          applied.
            sage: Permutation([1,2,3]).ishift(2)
            [1, 2, 3]
        """
        state = self._icondition(i)
        if state[0] is None:
            return self

        state, pos_im1, pos_i, pos_ip1 = state
        l = list(self)

        if state == '213':   #goes to 132
            l[pos_i]   = i-1
            l[pos_im1] = i+1
            l[pos_ip1] = i
        elif state == '132': #goes to 213
            l[pos_im1] = i
            l[pos_ip1] = i-1
            l[pos_i]   = i+1
        elif state == '231': #goes to 312
            l[pos_i]   = i+1
            l[pos_ip1] = i-1
            l[pos_im1] = i
        elif state == '312': #goes to 231
            l[pos_ip1] = i
            l[pos_im1] = i+1
            l[pos_i]   = i-1
        else:
            raise ValueError, "invalid state"


        return Permutation_class(l)



    def iswitch(self, i):
        """
        Returns an the i-switch of self.  If an i-switch of self can't be
        performed, then self is returned.


        An i-shift can be applied when i is not in between i-1 and i+1.
        The i-shift moves i to the other side, and switches the relative
        positions of i-1 and i+1 in place.

        EXAMPLES:
          Here, 2 is to the left of both 1 and 3.  A 2-switch can be
          applied which moves the 2 to the right and switches the
          relative order between 1 and 3.
            sage: Permutation([2,1,3]).iswitch(2)
            [3, 1, 2]

          Note that the movement is done in place:
            sage: Permutation([2,4,1,3]).iswitch(2)
            [3, 4, 1, 2]

          Since 2 is between 1 and 3 in [1,2,3], an 2-switch cannot be
          applied.
            sage: Permutation([1,2,3]).iswitch(2)
            [1, 2, 3]

        """
        if i not in range(2, len(self)):
            raise ValueError, "i (= %s) must between 2 and n-1"%i

        state = self._icondition(i)
        if state[0] is None:
            return self

        state, pos_im1, pos_i, pos_ip1 = state
        l = list(self)

        if state == '213':    #goes to 312
            l[pos_i]   = i+1
            l[pos_im1] = i-1
            l[pos_ip1] = i
        elif state == '132':  #goes to 231
            l[pos_im1] = i
            l[pos_ip1] = i+1
            l[pos_i]   = i-1
        elif state == '231':  #goes to 132
            l[pos_i]   = i-1
            l[pos_ip1] = i+1
            l[pos_im1] = i
        elif state == '312':  #goes to 213
            l[pos_ip1] = i
            l[pos_im1] = i-1
            l[pos_i]   = i+1
        else:
            raise ValueError, "invalid state"


        return Permutation_class(l)

    def runs(self):
        r"""
        Returns a list of the runs in the permutation p.

        REFERENCES:
            http://mathworld.wolfram.com/PermutationRun.html

        EXAMPLES:
            sage: Permutation([1,2,3,4]).runs()
            [[1, 2, 3, 4]]
            sage: Permutation([4,3,2,1]).runs()
            [[4], [3], [2], [1]]
            sage: Permutation([2,4,1,3]).runs()
            [[2, 4], [1, 3]]
        """
        p = self[:]
        runs = []
        current_value = p[0]
        current_run = [p[0]]
        for i in range(1, len(p)):
            if p[i] < current_value:
                runs.append(current_run)
                current_run = [p[i]]
            else:
                current_run.append(p[i])

            current_value = p[i]
        runs.append(current_run)

        return runs

    def longest_increasing_subsequence(self):
        r"""
        Returns a list of the longest increasing subsequences of
        the permutation p.

        EXAMPLES:
            sage: Permutation([2,3,4,1]).longest_increasing_subsequence()
            [[2, 3, 4]]
        """
        runs = self.runs()
        lis = []
        max_length = max([len(r) for r in runs])
        for run in runs:
            if len(run) == max_length:
                lis.append(run)
        return lis

    def cycle_type(self):
        r"""
        Returns a partition of len(p) corresponding to the cycle
        type of p.  This is a non-increasing sequence of the
        cycle lengths of p.

        EXAMPLES:
            sage: Permutation([3,1,2,4]).cycle_type()
            [3, 1]
        """
        cycle_type = [len(c) for c in self.to_cycles()]
        cycle_type.sort(reverse=True)
        return sage.combinat.partition.Partition(cycle_type)

    def to_lehmer_code(self):
        r"""
        Returns the Lehmer code of the permutation p.

        EXAMPLES:
            sage: p = Permutation([2,1,3])
            sage: p.to_lehmer_code()
            [1, 0, 0]
            sage: q = Permutation([3,1,2])
            sage: q.to_lehmer_code()
            [2, 0, 0]
        """
        p = self[:]
        n = len(p)
        lehmer = [None]*n
        checked = [0]*n
        ident = range(1,n+1)

        #Compute the factoradic / Lehmer code of the permutation
        for i in range(n):
            j = 0
            pos = 0
            while j < n:
                if checked[j] != 1:
                    if ident[j] == p[i]:
                        checked[j] = 1
                        break
                    pos += 1
                j += 1
            lehmer[i] = pos

        return lehmer


    def to_lehmer_cocode(self):
        r"""
        Returns the Lehmer cocode of p.

        EXAMPLES:
            sage: p = Permutation([2,1,3])
            sage: p.to_lehmer_cocode()
            [0, 1, 0]
            sage: q = Permutation([3,1,2])
            sage: q.to_lehmer_cocode()
            [0, 1, 1]
        """
        p = self[:]
        n = len(p)
        cocode = [0] * n
        for i in range(1, n):
            for j in range(0, i):
                if p[j] > p[i]:
                    cocode[i] += 1
        return cocode



    #################
    # Reduced Words #
    #################

    def reduced_word(self):
        r"""
        Returns the reduced word of a permutation.

        EXAMPLES:
            sage: Permutation([3,5,4,6,2,1]).reduced_word()
            [2, 1, 4, 3, 2, 4, 3, 5, 4, 5]

        """
        code = self.to_lehmer_code()
        reduced_word = []
        for piece in  [ [ i + code[i] - j for j in range(code[i])] for i in range(len(code))]:
            reduced_word += piece

        return reduced_word

    def reduced_words(self):
        r"""
        Returns a list of the reduced words of the permutation p.

        EXAMPLES:
            sage: Permutation([2,1,3]).reduced_words()
            [[1]]
            sage: Permutation([3,1,2]).reduced_words()
            [[2, 1]]
            sage: Permutation([3,2,1]).reduced_words()
            [[1, 2, 1], [2, 1, 2]]
            sage: Permutation([3,2,4,1]).reduced_words()
            [[1, 2, 3, 1], [1, 2, 1, 3], [2, 1, 2, 3]]
        """
        p = self[:]
        n = len(p)
        rws = []
        descents = self.descents()

        if len(descents) == 0:
            return [[]]

        for d in descents:
            pp = p[:d] + [p[d+1]] + [p[d]] + p[d+2:]
            z = lambda x: x + [d+1]
            rws += (map(z, Permutation(pp).reduced_words()))

        return rws



    def reduced_word_lexmin(self):
        r"""
        Returns a lexicographically minimal reduced word of a permutation.

        EXAMPLES:
            sage: Permutation([3,4,2,1]).reduced_word_lexmin()
            [1, 2, 1, 3, 2]
        """
        cocode = self.inverse().to_lehmer_cocode()

        rw = []
        for i in range(len(cocode)):
            piece = [j+1 for j in range(i-cocode[i],i)]
            piece.reverse()
            rw += piece

        return rw


    ################
    # Fixed Points #
    ################

    def fixed_points(self):
        r"""
        Returns a list of the fixed points of the permutation p.

        EXAMPLES:
            sage: Permutation([1,3,2,4]).fixed_points()
            [1, 4]
            sage: Permutation([1,2,3,4]).fixed_points()
            [1, 2, 3, 4]
        """
        fixed_points = []
        for i in range(len(self)):
            if i+1 == self[i]:
                fixed_points.append(i+1)

        return fixed_points

    def number_of_fixed_points(self):
        r"""
        Returns the number of fixed points of the permutation p.

        EXAMPLES:
            sage: Permutation([1,3,2,4]).number_of_fixed_points()
            2
            sage: Permutation([1,2,3,4]).number_of_fixed_points()
            4
        """

        return len(self.fixed_points())


    ############
    # Recoils  #
    ############
    def recoils(self):
        r"""
        Returns the list of the positions of the recoils of the permutation p.

        A recoil of a permutation is an integer i such that i+1 is to the
        left of it.


        EXAMPLES:
            sage: Permutation([1,4,3,2]).recoils()
            [2, 3]
        """
        p = self
        recoils  = []
        for i in range(len(p)):
            if p[i] != len(self) and self.index(p[i]+1) < i:
                recoils.append(i)

        return recoils

    def number_of_recoils(self):
        r"""
        Returns the number of recoils of the permutation p.

        EXAMPLES:
            sage: Permutation([1,4,3,2]).number_of_recoils()
            2

        """
        return len(self.recoils())

    def recoils_composition(self):
        """
        Returns the composition corresponding to recoils of
        the permutation.

        EXAMPLES:
            sage: Permutation([1,3,2,4]).recoils_composition()
            [3]
        """
        d = self.recoils()
        d = [ -1 ] + d
        return [ d[i+1]-d[i] for i in range(len(d)-1)]


    ############
    # Descents #
    ############

    def descents(self, final_descent=False):
        r"""
        Returns the list of the descents of the permutation p.

        A descent of a permutation is an integer i such that p[i]>p[i+1].
        With the final_descent option, the last position of a non empty permutation
        is also considered as a descent.

        EXAMPLES:
            sage: Permutation([1,4,3,2]).descents()
            [1, 2]
            sage: Permutation([1,4,3,2]).descents(final_descent=True)
            [1, 2, 3]
        """
        p = self
        descents = []
        for i in range(len(p)-1):
            if p[i] > p[i+1]:
                descents.append(i)

        if final_descent:
            descents.append(len(p)-1)

        return descents

    def idescents(self, final_descent=False):
        """
        Returns a list of the idescents of self, that is the list of the
        descents of self's inverse.

        With the final_descent option, the last position of a non empty permutation
        is also considered as a descent.

        EXAMPLES:
            sage: Permutation([1,4,3,2]).idescents()
            [1, 2]
            sage: Permutation([1,4,3,2]).idescents(final_descent=True)
            [1, 2, 3]
        """
        return self.inverse().descents(final_descent=final_descent)

    def idescents_signature(self, final_descent=False):
        """
        Each position in self is mapped to -1 if it is an idescent and 1
        if it is not an idescent.

        EXAMPLES:
            sage: Permutation([1,4,3,2]).idescents()
            [1, 2]
            sage: Permutation([1,4,3,2]).idescents_signature()
            [1, -1, -1, 1]
        """
        idescents = self.idescents(final_descent=final_descent)
        d = {True:-1, False:1}
        return [d[i in idescents] for i in range(len(self))]

    def number_of_descents(self, final_descent=False):
        r"""
        Returns the number of descents of the permutation p.

        EXAMPLES:
            sage: Permutation([1,4,3,2]).number_of_descents()
            2
            sage: Permutation([1,4,3,2]).number_of_descents(final_descent=True)
            3
        """
        return len(self.descents(final_descent))

    def number_of_idescents(self, final_descent=False):
        r"""
        Returns the number of descents of the permutation p.

        EXAMPLES:
            sage: Permutation([1,4,3,2]).number_of_idescents()
            2
            sage: Permutation([1,4,3,2]).number_of_idescents(final_descent=True)
            3
        """
        return len(self.idescents(final_descent))


    def descents_composition(self):
        """
        Returns the composition corresponding to the descents of
        the permutation.

        EXAMPLES:
            sage: Permutation([1,3,2,4]).descents_composition()
            [2, 2]
        """
        d = self.descents()
        d = [ -1 ] + d + [len(self)-1]
        return Composition([ d[i+1]-d[i] for i in range(len(d)-1)])

    def descent_polynomial(self):
        r"""
        Returns the descent polynomial of the permutation p.

        The descent polymomial of p is the product of all the
        z[p[i]] where i ranges over the descents of p.

        REFERENCES:
            Garsia and Stanton 1984

        EXAMPLES:
            sage: Permutation([2,1,3]).descent_polynomial()
            z1
            sage: Permutation([4,3,2,1]).descent_polynomial()
            z1*z2^2*z3^3
        """
        p = self
        z = []
        P = PolynomialRing(ZZ, len(p), 'z')
        z = P.gens()
        result = 1
        pol = 1
        for i in range(len(p)-1):
            pol *= z[p[i]-1]
            if p[i] > p[i+1]:
                result *= pol

        return result


    ##############
    # Major Code #
    ##############

    def major_index(self, final_descent=False):
        r"""
        Returns the major index of the permutation p.

        The major index is the sum of the descents of p. Since our permutation
        indices are 0-based, we need to add one the number of descents.

        EXAMPLES:
            sage: Permutation([2,1,3]).major_index()
            1
            sage: Permutation([3,4,1,2]).major_index()
            2
            sage: Permutation([4,3,2,1]).major_index()
            6
        """
        descents = self.descents(final_descent)

        return sum(descents)+len(descents)

    def imajor_index(self, final_descent=False):
        """
        Returns the inverse major index of the permutation self, which is the
        major index of the inverse of self.

        The major index is the sum of the descents of p. Since our permutation
        indices are 0-based, we need to add one the number of descents.

        EXAMPLES:
            sage: Permutation([2,1,3]).imajor_index()
            1
            sage: Permutation([3,4,1,2]).imajor_index()
            2
            sage: Permutation([4,3,2,1]).imajor_index()
            6
        """
        idescents = self.idescents(final_descent)

        return sum(idescents)+len(idescents)

    def to_major_code(self, final_descent=False):
        r"""
        Returns the major code of the permutation p, which is defined
        as the list [m1-m2, m2-m3,..,mn] where mi := maj(pi) is the
        major indices of the permutation math obtained by erasing the
        letters smaller than math in p.

        REFERENCES:
            Carlitz, L. 'q-Bernoulli and Eulerian Numbers' Trans. Amer. Math. Soc. 76 (1954) 332-350
            Skandera, M. 'An Eulerian Partner for Inversions', Sem. Lothar. Combin. 46 (2001) B46d.

        EXAMPLES:
            sage: Permutation([9,3,5,7,2,1,4,6,8]).to_major_code()
            [5, 0, 1, 0, 1, 2, 0, 1, 0]
            sage: Permutation([2,8,4,3,6,7,9,5,1]).to_major_code()
            [8, 3, 3, 1, 4, 0, 1, 0, 0]
        """
        p = self
        major_indices = [0]*(len(p)+1)
        smaller = p[:]
        for i in range(len(p)):
            major_indices[i] = Permutation(smaller).major_index(final_descent)
            #Create the permutation that "erases" all the numbers
            #smaller than i+1
            smaller.remove(1)
            smaller = [i-1 for i in smaller]

        major_code = [ major_indices[i] - major_indices[i+1] for i in range(len(p)) ]
        return major_code

    #########
    # Peaks #
    #########

    def peaks(self):
        r"""
        Returns a list of the peaks of the permutation p.

        A peak of a permutation is an integer i such that
        p[i-1] <= p[i] and p[i] > p[i+1].

        EXAMPLES:
            sage: Permutation([1,3,2,4,5]).peaks()
            [1]
            sage: Permutation([4,1,3,2,6,5]).peaks()
            [2, 4]
        """
        p = self
        peaks = []
        for i in range(1,len(p)-1):
            if p[i-1] <= p[i] and p[i] > p[i+1]:
                peaks.append(i)

        return peaks


    def number_of_peaks(self):
        r"""
        Returns the number of peaks of the permutation p.

        A peak of a permutation is an integer i such that
        p[i-1] <= p[i] and p[i] > p[i+1].

        EXAMPLES:
            sage: Permutation([1,3,2,4,5]).number_of_peaks()
            1
            sage: Permutation([4,1,3,2,6,5]).number_of_peaks()
            2
        """
        return len(self.peaks())

    #############
    # Saliances #
    #############

    def saliances(self):
        r"""
        Returns a list of the saliances of the permutation p.

        A saliance of a permutation p is an integer i such that
        p[i] > p[j] for all j > i.

        EXAMPLES:
            sage: Permutation([2,3,1,5,4]).saliances()
            [3, 4]
            sage: Permutation([5,4,3,2,1]).saliances()
            [0, 1, 2, 3, 4]
        """
        p = self
        saliances = []
        for i in range(len(p)):
            is_saliance = True
            for j in range(i+1, len(p)):
                if p[i] <= p[j]:
                    is_saliance = False
            if is_saliance:
                saliances.append(i)

        return saliances


    def number_of_saliances(self):
        r"""
        Returns the number of saliances of the permutation p.

        EXAMPLES:
            sage: Permutation([2,3,1,5,4]).number_of_saliances()
            2
            sage: Permutation([5,4,3,2,1]).number_of_saliances()
            5
        """
        return len(self.saliances())

    ################
    # Bruhat Order #
    ################
    def bruhat_lequal(self, p2):
        r"""
        Returns True if self is less than p2 in the Bruhat order.

        EXAMPLES:
            sage: Permutation([2,4,3,1]).bruhat_lequal(Permutation([3,4,2,1]))
            True

        """
        p1 = self
        n1 = len(p1)
        n2 = len(p2)

        if n1 == 0:
            return True

        if p1[0] > p2[0] or p1[n1-1] < p2[n1-1]:
            return False

        for i in range(n1):
            c = 0
            for j in range(n1):
                if p2[j] > i+1:
                    c += 1
                if p1[j] > i+1:
                    c -= 1
                if c < 0:
                    return False

        return True

    def weak_excedences(self):
        """
        Returns all the numbers self[i] such that self[i] >= i+1.

        EXAMPLES:
            sage: Permutation([1,4,3,2,5]).weak_excedences()
            [1, 4, 3, 5]
        """
        res = []
        for i in range(len(self)):
            if self[i] >= i + 1:
                res.append(self[i])
        return res


    def bruhat_inversions(self):
        r"""
        Returns the list of inversions of p such that the application of
        this inversion to p decrements its number of inversions.

        Equivalently, it returns the list of pairs (i,j), i<j  such that
        p[i] > p[j] and such that there exists no k between i and j
        satisfying p[i] > p[k].


        EXAMPLES:
            sage: Permutation([5,2,3,4,1]).bruhat_inversions()
            [[0, 1], [0, 2], [0, 3], [1, 4], [2, 4], [3, 4]]
            sage: Permutation([6,1,4,5,2,3]).bruhat_inversions()
            [[0, 1], [0, 2], [0, 3], [2, 4], [2, 5], [3, 4], [3, 5]]
        """
        return __builtin__.list(self.bruhat_inversions_iterator())

    def bruhat_inversions_iterator(self):
        """
        Returns the iterator for the inversions of p such that the
        application of this inversion to p decrements its number of inversions.

        EXAMPLES:
            sage: list(Permutation([5,2,3,4,1]).bruhat_inversions_iterator())
            [[0, 1], [0, 2], [0, 3], [1, 4], [2, 4], [3, 4]]
            sage: list(Permutation([6,1,4,5,2,3]).bruhat_inversions_iterator())
            [[0, 1], [0, 2], [0, 3], [2, 4], [2, 5], [3, 4], [3, 5]]
        """
        p = self
        n = len(p)

        for i in range(n-1):
            for j in range(i+1,n):
                if p[i] > p[j]:
                    ok = True
                    for k in range(i+1, j):
                        if p[i] > p[k] and p[k] > p[j]:
                            ok = False
                            break
                    if ok:
                        yield [i,j]


    def bruhat_succ(self):
        r"""
        Returns a list of the permutations strictly greater than p in the
        Bruhat order such that there is no permutation between one of those
        and p.

        EXAMPLES:
            sage: Permutation([6,1,4,5,2,3]).bruhat_succ()
            [[6, 4, 1, 5, 2, 3],
             [6, 2, 4, 5, 1, 3],
             [6, 1, 5, 4, 2, 3],
             [6, 1, 4, 5, 3, 2]]

        """
        return __builtin__.list(self.bruhat_succ_iterator())

    def bruhat_succ_iterator(self):
        """
        An iterator for the permutations that are strictly greater than p in
        the Bruhat order such that there is no permutation between one of
        those and p.

        EXAMPLES:
            sage: [x for x in Permutation([6,1,4,5,2,3]).bruhat_succ_iterator()]
            [[6, 4, 1, 5, 2, 3],
             [6, 2, 4, 5, 1, 3],
             [6, 1, 5, 4, 2, 3],
             [6, 1, 4, 5, 3, 2]]
        """
        p = self
        n = len(p)

        for z in Permutation(map(lambda x: n+1-x, p)).bruhat_inversions_iterator():
            pp = p[:]
            pp[z[0]] = p[z[1]]
            pp[z[1]] = p[z[0]]
            yield Permutation(pp)



    def bruhat_pred(self):
        r"""
        Returns a list of the permutations strictly smaller than p in the
        Bruhat order such that there is no permutation between one of those
        and p.


        EXAMPLES:
        sage: Permutation([6,1,4,5,2,3]).bruhat_pred()
        [[1, 6, 4, 5, 2, 3],
         [4, 1, 6, 5, 2, 3],
         [5, 1, 4, 6, 2, 3],
         [6, 1, 2, 5, 4, 3],
         [6, 1, 3, 5, 2, 4],
         [6, 1, 4, 2, 5, 3],
         [6, 1, 4, 3, 2, 5]]

        """
        return __builtin__.list(self.bruhat_pred_iterator())

    def bruhat_pred_iterator(self):
        """
        An iterator for the permutations strictly smaller than p in the
        Bruhat order such that there is no permutation between one of those
        and p.


        EXAMPLES:
        sage: [x for x in Permutation([6,1,4,5,2,3]).bruhat_pred_iterator()]
        [[1, 6, 4, 5, 2, 3],
         [4, 1, 6, 5, 2, 3],
         [5, 1, 4, 6, 2, 3],
         [6, 1, 2, 5, 4, 3],
         [6, 1, 3, 5, 2, 4],
         [6, 1, 4, 2, 5, 3],
         [6, 1, 4, 3, 2, 5]]

        """
        p = self
        n = len(p)
        for z in p.bruhat_inversions_iterator():
            pp = p[:]
            pp[z[0]] = p[z[1]]
            pp[z[1]] = p[z[0]]
            yield Permutation(pp)


    def bruhat_smaller(self):
        r"""
        Returns a the combinatorial class of  permutations smaller than or equal
        to p in the Bruhat order.

        EXAMPLES:
            sage: Permutation([4,1,2,3]).bruhat_smaller().list()
            [[1, 2, 3, 4],
             [1, 2, 4, 3],
             [1, 3, 2, 4],
             [1, 4, 2, 3],
             [2, 1, 3, 4],
             [2, 1, 4, 3],
             [3, 1, 2, 4],
             [4, 1, 2, 3]]

        """
        return StandardPermutations_bruhat_smaller(self)


    def bruhat_greater(self):
        r"""
        Returns the combinatorial class of permutations greater than or equal
        to p in the Bruhat order.

        EXAMPLES:
            sage: Permutation([4,1,2,3]).bruhat_greater().list()
            [[4, 1, 2, 3],
             [4, 1, 3, 2],
             [4, 2, 1, 3],
             [4, 2, 3, 1],
             [4, 3, 1, 2],
             [4, 3, 2, 1]]

        """

        return StandardPermutations_bruhat_greater(self)

    ########################
    # Permutohedron  Order #
    ########################

    def permutohedron_lequal(self, p2, side="right"):
        r"""
        Returns True if self is less than p2 in the permutohedron order.

        By default, the computations are done in the right permutohedron.
        If you pass the option side='left', then they will be done in the
        left permutohedron.

        EXAMPLES:
            sage: p = Permutation([3,2,1,4])
            sage: p.permutohedron_lequal(Permutation([4,2,1,3]))
            False
            sage: p.permutohedron_lequal(Permutation([4,2,1,3]), side='left')
            True
        """
        p1 = self
        n1 = len(p1)
        n2 = len(p2)

        l1 = p1.number_of_inversions()
        l2 = p2.number_of_inversions()

        if l1 > l2:
            return False

        if side == "right":
            prod = p1._left_to_right_multiply_on_right(p2.inverse())
        else:
            prod = p1._left_to_right_multiply_on_left(p2.inverse())

        return prod.number_of_inversions() == l2 - l1

    def permutohedron_succ(self, side="right"):
        r"""
        Returns a list of the permutations strictly greater than p in the
        permutohedron order such that there is no permutation between
        one of those and p.

        By default, the computations are done in the right permutohedron.
        If you pass the option side='left', then they will be done in the
        left permutohedron.

        EXAMPLES:
            sage: p = Permutation([4,2,1,3])
            sage: p.permutohedron_succ()
            [[4, 2, 3, 1]]
            sage: p.permutohedron_succ(side='left')
            [[4, 3, 1, 2]]

        """
        p = self
        n = len(p)
        succ = []
        if side == "right":
            rise = lambda perm: filter(lambda i: perm[i] < perm[i+1], range(0,n-1))
            for i in rise(p):
                pp = p[:]
                pp[i] = p[i+1]
                pp[i+1] = p[i]
                succ.append(Permutation(pp))
        else:
            advance = lambda perm: filter(lambda i: perm.index(i) < perm.index(i+1), range(1,n))
            for i in advance(p):
                pp = p[:]
                pp[p.index(i)] = i+1
                pp[p.index(i+1)] = i
                succ.append(Permutation(pp))

        return succ


    def permutohedron_pred(self, side="right"):
        r"""
        Returns a list of the permutations strictly smaller than p in the
        permutohedron order such that there is no permutation between
        one of those and p.

        By default, the computations are done in the right permutohedron.
        If you pass the option side='left', then they will be done in the
        left permutohedron.

        EXAMPLES:
            sage: p = Permutation([4,2,1,3])
            sage: p.permutohedron_pred()
            [[2, 4, 1, 3], [4, 1, 2, 3]]
            sage: p.permutohedron_pred(side='left')
            [[4, 1, 2, 3], [3, 2, 1, 4]]


        """
        p = self
        n = len(p)
        pred = []
        if side == "right":
            for d in p.descents():
                pp = p[:]
                pp[d] = p[d+1]
                pp[d+1] = p[d]
                pred.append(Permutation(pp))
        else:
            recoil = lambda perm: filter(lambda j: perm.index(j) > perm.index(j+1), range(1,n))
            for i in recoil(p):
                pp = p[:]
                pp[p.index(i)] = i+1
                pp[p.index(i+1)] = i
                pred.append(Permutation(pp))
        return pred


    def permutohedron_smaller(self, side="right"):
        r"""
        Returns a list of permutations smaller than or equal to p in the
        permutohedron order.

        By default, the computations are done in the right permutohedron.
        If you pass the option side='left', then they will be done in the
        left permutohedron.

        EXAMPLES:
            sage: Permutation([4,2,1,3]).permutohedron_smaller()
            [[1, 2, 3, 4],
             [1, 2, 4, 3],
             [1, 4, 2, 3],
             [2, 1, 3, 4],
             [2, 1, 4, 3],
             [2, 4, 1, 3],
             [4, 1, 2, 3],
             [4, 2, 1, 3]]

            sage: Permutation([4,2,1,3]).permutohedron_smaller(side='left')
            [[1, 2, 3, 4],
             [1, 3, 2, 4],
             [2, 1, 3, 4],
             [2, 3, 1, 4],
             [3, 1, 2, 4],
             [3, 2, 1, 4],
             [4, 1, 2, 3],
             [4, 2, 1, 3]]


        """

        return transitive_ideal(lambda x: x.permutohedron_pred(side), self)


    def permutohedron_greater(self, side="right"):
        r"""
        Returns a list of permutations greater than or equal to p in the
        permutohedron order.

        By default, the computations are done in the right permutohedron.
        If you pass the option side='left', then they will be done in the
        left permutohedron.

        EXAMPLES:
            sage: Permutation([4,2,1,3]).permutohedron_greater()
            [[4, 2, 1, 3], [4, 2, 3, 1], [4, 3, 2, 1]]
            sage: Permutation([4,2,1,3]).permutohedron_greater(side='left')
            [[4, 2, 1, 3], [4, 3, 1, 2], [4, 3, 2, 1]]


        """

        return transitive_ideal(lambda x: x.permutohedron_succ(side), self)


    ############
    # Patterns #
    ############

    def has_pattern(self, patt):
        r"""
        Returns the boolean answering the question 'Is patt a patter appearing in
        permutation p?'

        EXAMPLES:
            sage: Permutation([3,5,1,4,6,2]).has_pattern([1,3,2])
            True
        """
        p = self
        n = len(p)
        l = len(patt)
        if l > n:
            return False
        for pos in subword.Subwords(range(n),l):
            if to_standard(map(lambda z: p[z] , pos)) == patt:
                return True
        return False

    def avoids(self, patt):
        """
        Returns True if the permutation avoid the pattern patt and False
        otherwise.

        EXAMPLES:
            sage: Permutation([6,2,5,4,3,1]).avoids([4,2,3,1])
            False
            sage: Permutation([6,1,2,5,4,3]).avoids([4,2,3,1])
            True
            sage: Permutation([6,1,2,5,4,3]).avoids([3,4,1,2])
            True
        """
        return not self.has_pattern(patt)

    def pattern_positions(self, patt):
        r"""
        Returns the list of positions where the pattern patt appears
        in p.

        EXAMPLES:
            sage: Permutation([3,5,1,4,6,2]).pattern_positions([1,3,2])
            [[0, 1, 3], [2, 3, 5], [2, 4, 5]]

        """
        p = self

        return __builtin__.list(itertools.ifilter(lambda pos: to_standard(map(lambda z: p[z], pos)) == patt, subword.Subwords(range(len(p)), len(patt)).iterator() ))


    def reverse(self):
        """
        Returns the permutation obtained by reversing the list.

        EXAMPLES:
            sage: Permutation([3,4,1,2]).reverse()
            [2, 1, 4, 3]
            sage: Permutation([1,2,3,4,5]).reverse()
            [5, 4, 3, 2, 1]
        """
        return Permutation_class( [i for i in reversed(self)] )


    def complement(self):
        """
        Returns the complement of the permutation which is obtained
        by replacing each value x in the list with n - x + 1

        EXAMPLES:
            sage: Permutation([1,2,3]).complement()
            [3, 2, 1]
            sage: Permutation([1, 3, 2]).complement()
            [3, 1, 2]
        """
        n = len(self)
        return Permutation_class( map(lambda x: n - x + 1, self) )

    def dict(self):
        """
        Returns a dictionary corresponding to the permutation.

        EXAMPLES:
            sage: p = Permutation([2,1,3])
            sage: d = p.dict()
            sage: d[1]
            2
            sage: d[2]
            1
            sage: d[3]
            3
        """
        d = {}
        for i in range(len(self)):
            d[i+1] = self[i]
        return d

    def action(self, a):
        """
        Return the action of the permutation on a list.

        EXAMPLES:
            sage: p = Permutation([2,1,3])
            sage: a = range(3)
            sage: p.action(a)
            [1, 0, 2]
            sage: b = [1,2,3,4]
            sage: p.action(b)
            Traceback (most recent call last):
            ...
            ValueError: len(a) must equal len(self)
        """
        if len(a) != len(self):
            raise ValueError, "len(a) must equal len(self)"
        return map(lambda i: a[self[i]-1], range(len(a)))

    ######################
    # Robinson-Schensted #
    ######################

    def robinson_schensted(self):
        """
        Returns the pair of standard tableau obtained by running the
        Robinson-Schensted Algorithm on self.

        EXAMPLES:
            sage: p = Permutation([6,2,3,1,7,5,4])
            sage: p.robinson_schensted()
            [[[1, 3, 4], [2, 5], [6, 7]], [[1, 3, 5], [2, 6], [4, 7]]]

        """

        p = [[]]
        q = [[]]

        for i in range(1, len(self)+1):
            #Row insert self[i-1] into p
            row_counter = 0
            r = p[row_counter]
            x = self[i-1]
            while max(r+[0]) > x:
                y = min(filter(lambda z: z > x, r))
                r[r.index(y)] = x
                x = y
                row_counter += 1
                if row_counter == len(p):
                    p.append([])
                r = p[row_counter]
            r.append(x)


            #Insert i into q in the same place as we inserted
            #i into p
            if row_counter == len(q):
                q.append([])
            q[row_counter].append(i)


        return [tableau.Tableau(p),tableau.Tableau(q)]

    def left_tableau(self):
        """
        Returns the right standard tableau after performing the RSK algorithm
        on self.

        EXAMPLES:
            sage: Permutation([1,4,3,2]).left_tableau()
            [[1, 2], [3], [4]]
        """
        return self.robinson_schensted()[0]

    def right_tableau(self):
        """
        Returns the right standard tableau after performing the RSK algorithm
        on self.

        EXAMPLES:
            sage: Permutation([1,4,3,2]).right_tableau()
            [[1, 2], [3], [4]]
        """
        return self.robinson_schensted()[1]

    def remove_extra_fixed_points(self):
        """
        Returns the permutation obtained by removing any
        fixed points at the end of self.

        EXAMPLES:
            sage: Permutation([2,1,3]).remove_extra_fixed_points()
            [2, 1]
            sage: Permutation([1,2,3,4]).remove_extra_fixed_points()
            [1]
        """
        #Strip off all extra fixed points at the end of
        #the permutation.
        i = len(self)-1
        while i >= 1:
            if i != self[i] - 1:
                break
            i -= 1
        return Permutation_class(self[:i+1])



################################################################

def Arrangements(mset, k):
    r"""
    An arrangement of mset is an ordered selection without repetitions
    and is represented by a list that contains only elements from
    mset, but maybe in a different order.

    \code{Arrangements} returns the combinatorial class of arrangements
    of the multiset mset that contain k elements.

    EXAMPLES:
        sage: mset = [1,1,2,3,4,4,5]
        sage: Arrangements(mset,2).list()
        [[1, 1],
         [1, 2],
         [1, 3],
         [1, 4],
         [1, 5],
         [2, 1],
         [2, 3],
         [2, 4],
         [2, 5],
         [3, 1],
         [3, 2],
         [3, 4],
         [3, 5],
         [4, 1],
         [4, 2],
         [4, 3],
         [4, 4],
         [4, 5],
         [5, 1],
         [5, 2],
         [5, 3],
         [5, 4]]
         sage: Arrangements(mset,2).count()
         22
         sage: Arrangements( ["c","a","t"], 2 ).list()
         [['c', 'a'], ['c', 't'], ['a', 'c'], ['a', 't'], ['t', 'c'], ['t', 'a']]
         sage: Arrangements( ["c","a","t"], 3 ).list()
         [['c', 'a', 't'],
          ['c', 't', 'a'],
          ['a', 'c', 't'],
          ['a', 't', 'c'],
          ['t', 'c', 'a'],
          ['t', 'a', 'c']]
    """
    mset = list(mset)
    if map(mset.index, mset) == range(len(mset)):
        return Arrangements_setk(mset, k)
    else:
        return Arrangements_msetk(mset, k)



class Permutations_nk(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS:
            sage: P = Permutations([3,2])
            sage: P == loads(dumps(P))
            True
        """
        self.n = n
        self.k = k

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(3,2))
            'Permutations of {1,...,3} of length 2'
        """
        return "Permutations of {1,...,%s} of length %s"%(self.n, self.k)

    def iterator(self):
        """
        EXAMPLES:
            sage: [p for p in Permutations(3,2)]
            [[1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2]]
            sage: [p for p in Permutations(3,0)]
            [[]]
            sage: [p for p in Permutations(3,4)]
            []
        """
        for x in PermutationsNK(self.n, self.k):
            yield [i+1 for i in x]

    def count(self):
        """
        EXAMPLES:
            sage: Permutations(3,0).count()
            1
            sage: Permutations(3,1).count()
            3
            sage: Permutations(3,2).count()
            6
            sage: Permutations(3,3).count()
            6
            sage: Permutations(3,4).count()
            0
        """
        if self.k <= self.n and self.k >= 0:
            return factorial(self.n)/factorial(self.n-self.k)
        else:
            return 0

    def random(self):
        """
        EXAMPLES:
            sage: Permutations(3,2).random()
            [0, 1]
        """
        return sample(range(self.n), self.k)

class Permutations_mset(CombinatorialClass):
    def __init__(self, mset):
        """
        TESTS:
            sage: S = Permutations(['c','a','c'])
            sage: S == loads(dumps(S))
            True
        """
        self.mset = mset

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(['c','a','c']))
            "Permutations of the multi-set ['c', 'a', 'c']"
        """
        return "Permutations of the multi-set %s"%self.mset

    def iterator(self):
        """
        Algorithm based on:
        http://marknelson.us/2002/03/01/next-permutation/

        EXAMPLES:
            sage: [ p for p in Permutations(['c','t','t'])]
            [['c', 't', 't'], ['t', 'c', 't'], ['t', 't', 'c']]
        """
        mset = self.mset
        n = len(self.mset)
        lmset = __builtin__.list(mset)
        mset_list = map(lambda x: lmset.index(x), lmset)
        mset_list.sort()

        yield [lmset[x] for x in mset_list]

        if n == 1:
            return

        while True:
            one = n - 2
            two = n - 1
            j   = n - 1

            #starting from the end, find the first o such that
            #mset_list[o] < mset_list[o+1]
            while two > 0 and mset_list[one] >= mset_list[two]:
                one -= 1
                two -= 1

            if two == 0:
                return

            #starting from the end, find the first j such that
            #mset_list[j] > mset_list[one]
            while mset_list[j] <= mset_list[one]:
                j -= 1

            #Swap positions one and j
            t = mset_list[one]
            mset_list[one] = mset_list[j]
            mset_list[j] = t


            #Reverse the list between two and last
            i = int((n - two)/2)-1
            #mset_list = mset_list[:two] + [x for x in reversed(mset_list[two:])]
            while i >= 0:
                t = mset_list[ i + two ]
                mset_list[ i + two ] = mset_list[n-1 - i]
                mset_list[n-1 - i] = t
                i -= 1

            #Yield the permutation
            yield [lmset[x] for x in  mset_list]

        def count(self):
            """
            EXAMPLES:
                sage: Permutations([1,2,2]).count()
                3
            """
            lmset = __builtin__.list(mset)
            mset_list = map(lambda x: lmset.index(x), lmset)
            d = {}
            for i in mset_list:
                d[i] = d.get(i, 0) + 1

            c = factorial(len(lmset))
            for i in d:
                if i != 1:
                    c /= factorial(i)

            return c

class Permutations_set(CombinatorialClass):
    def __init__(self, set):
        """
        TESTS:
            sage: S = Permutations(['c','a','t'])
            sage: S == loads(dumps(S))
            True
        """
        self.set = set

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(['c','a','t']))
            "Permutations of the set ['c', 'a', 't']"
        """
        return "Permutations of the set %s"%self.set

    def iterator(self):
        """
        Algorithm based on:
        http://marknelson.us/2002/03/01/next-permutation/

        EXAMPLES:
            sage: [ p for p in Permutations(['c','a','t'])]
            [['c', 'a', 't'],
             ['c', 't', 'a'],
             ['a', 'c', 't'],
             ['a', 't', 'c'],
             ['t', 'c', 'a'],
             ['t', 'a', 'c']]
        """
        set = self.set
        n = len(self.set)
        lset = __builtin__.list(set)
        set_list = map(lambda x: lset.index(x), lset)
        set_list.sort()

        yield [lset[x] for x in set_list]

        if n == 1:
            return

        while True:
            one = n - 2
            two = n - 1
            j   = n - 1

            #starting from the end, find the first o such that
            #set_list[o] < set_list[o+1]
            while two > 0 and set_list[one] >= set_list[two]:
                one -= 1
                two -= 1

            if two == 0:
                return

            #starting from the end, find the first j such that
            #set_list[j] > set_list[one]
            while set_list[j] <= set_list[one]:
                j -= 1

            #Swap positions one and j
            t = set_list[one]
            set_list[one] = set_list[j]
            set_list[j] = t


            #Reverse the list between two and last
            i = int((n - two)/2)-1
            #set_list = set_list[:two] + [x for x in reversed(set_list[two:])]
            while i >= 0:
                t = set_list[ i + two ]
                set_list[ i + two ] = set_list[n-1 - i]
                set_list[n-1 - i] = t
                i -= 1

            #Yield the permutation
            yield [lset[x] for x in set_list]

    def count(self):
        """
        EXAMPLES:
        sage: Permutations([1,2,3]).count()
        6
        """
        return factorial(len(self.set))

    def random(self):
        """
        EXAMPLES:
        sage: Permutations([1,2,3]).random()
        [1, 2, 3]
        """
        return sample(self.set, len(self.set))

class Permutations_msetk(CombinatorialClass):
    def __init__(self, mset, k):
        """
        TESTS:
            sage: P = Permutations([1,2,2],2)
            sage: P == loads(dumps(P))
            True
        """
        self.mset = mset
        self.k = k

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations([1,2,2],2))
            'Permutations of the multi-set [1, 2, 2] of length 2'
        """
        return "Permutations of the multi-set %s of length %s"%(self.mset,self.k)

    def list(self):
        """
        EXAMPLES:
            sage: Permutations([1,2,3],2).list()
            [[1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2]]
            sage: Permutations([1,2,2],2).list()
            [[1, 2], [2, 1], [2, 2]]
        """

        mset = self.mset
        n = len(self.mset)
        lmset = __builtin__.list(mset)
        mset_list = map(lambda x: lmset.index(x), lmset)


        indices = eval(gap.eval('Arrangements(%s,%s)'%(mset_list, self.k)))
        return [[lmset[x] for x in ktuple] for ktuple in indices]


class Permutations_setk(CombinatorialClass):
    def __init__(self, set, k):
        """
        TESTS:
            sage: P = Permutations([1,2,2],2)
            sage: P == loads(dumps(P))
            True
        """
        self.set = set
        self.k = k

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations([1,2,3],2))
            'Permutations of the set [1, 2, 3] of length 2'
        """
        return "Permutations of the set %s of length %s"%(self.set,self.k)

    def iterator(self):
        """
        EXAMPLES:
            sage: [i for i in Permutations([1,2,3],2)]
            [[1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2]]
        """
        for perm in PermutationsNK(len(self.set), self.k):
            yield [self.set[x] for x in perm]

    def random(self):
        """
        EXAMPLES:
            sage: Permutations([1,2,3],2).random()
            [1, 2]
        """
        return sample(self.set, self.k)


class Arrangements_msetk(Permutations_msetk):
    def __repr__(self):
        """
        TESTS:
            sage: repr(Arrangements([1,2,2],2))
            'Arrangements of the multi-set [1, 2, 2] of length 2'
        """
        return "Arrangements of the multi-set %s of length %s"%(self.mset,self.k)

class Arrangements_setk(Permutations_setk):
    def __repr__(self):
        """
        TESTS:
            sage: repr(Arrangements([1,2,3],2))
            'Arrangements of the set [1, 2, 3] of length 2'
        """
        return "Arrangements of the set %s of length %s"%(self.set,self.k)


class StandardPermutations_all(CombinatorialClass):
    def __init__(self):
        """
        TESTS:
            sage: SP = Permutations()
            sage: SP == loads(dumps(SP))
            True
        """
        self.object_class = Permutation_class

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations())
            'Standard permutations'
        """
        return "Standard permutations"

    def __contains__(self,x):
        """
        TESTS:
            sage: [] in Permutations()
            False
            sage: [1] in Permutations()
            True
            sage: [2] in Permutations()
            False
            sage: [1,2] in Permutations()
            True
            sage: [2,1] in Permutations()
            True
            sage: [1,2,2] in Permutations()
            False
            sage: [3,1,5,2] in Permutations()
            False
            sage: [3,4,1,5,2] in Permutations()
            True
        """
        if isinstance(x, Permutation_class):
            return True
        elif isinstance(x, __builtin__.list):
            if len(x) == 0:
                return False
            copy = x[:]
            copy.sort()
            if copy != range(1, len(x)+1):
                return False
            return True
        else:
            return False

    def list(self):
        """
        EXAMPLES:
            sage: Permutations().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class StandardPermutations_n(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: SP = Permutations(3)
            sage: SP == loads(dumps(SP))
            True
        """
        self.n = n
        self.object_class = Permutation_class


    def __contains__(self,x):
        """
        TESTS:
            sage: [1,2] in Permutations(2)
            True
            sage: [1,2] in Permutations(3)
            False
            sage: [3,2,1] in Permutations(3)
            True
        """

        return x in Permutations() and len(x) == self.n

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(3))
            'Standard permutations of 3'
        """
        return "Standard permutations of %s"%self.n

    def iterator(self):
        """
        EXAMPLES:
            sage: [p for p in Permutations(3)]
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        """
        for p in Permutations_set(range(1,self.n+1)):
            yield Permutation_class(p)

    def count(self):
        """
        EXAMPLES:
            sage: Permutations(3).count()
            6
            sage: Permutations(4).count()
            24
        """
        return factorial(self.n)


    def identity(self):
        r"""
        Returns the identity permutation of length n.

        EXAMPLES:
            sage: Permutations(4).identity()
            [1, 2, 3, 4]
        """

        return Permutation_class(range(1,self.n+1))

    def unrank(self, r):
        """
        EXAMPLES:
            sage: SP3 = Permutations(3)
            sage: l = map(SP3.unrank, range(6))
            sage: l == SP3.list()
            True
        """
        if r >= factorial(self.n) or r < 0:
            raise ValueError
        else:
            return from_rank(self.n, r)

    def rank(self, p):
        """
        EXAMPLES:
            sage: SP3 = Permutations(3)
            sage: map(SP3.rank, SP3)
            [0, 1, 2, 3, 4, 5]
        """
        if p in self:
            return Permutation(p).rank()
        else:
            raise ValueError, "x not in self"

    def random(self):
        """
        EXAMPLES:
            sage: Permutations(4).random()
            [1, 3, 2, 4]
        """
        r = randint(0, int(factorial(self.n)-1))
        return self.unrank(r)



#############################
# Constructing Permutations #
#############################
def from_permutation_group_element(pge):
    """
    Returns a Permutation give a PermutationGroupElement pge.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: pge = PermutationGroupElement([(1,2),(3,4)])
        sage: permutation.from_permutation_group_element(pge)
        [2, 1, 4, 3]
    """

    if not isinstance(pge, PermutationGroupElement):
        raise TypeError, "pge (= %s) must be a PermutationGroupElement"%pge

    return Permutation(pge.list())



def from_rank(n, rank):
    r"""
    Returns the permutation with the specified lexicographic
    rank.  The permutation is of the set [1,...,n].

    The permutation is computed without iteratiing through all
    of the permutations with lower rank.  This makes it efficient
    for large permutations.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: Permutation([3, 6, 5, 4, 2, 1]).rank()
        359
        sage: [permutation.from_rank(3, i) for i in range(6)]
        [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        sage: Permutations(6)[10]
        [1, 2, 4, 6, 3, 5]
        sage: permutation.from_rank(6,10)
        [1, 2, 4, 6, 3, 5]

    """

    #Find the factoradic of rank
    factoradic = [None] * n
    for j in range(1,n+1):
        factoradic[n-j] = Integer(rank % j)
        rank = int(rank) / int(j)

    return from_lehmer_code(factoradic)

def from_inversion_vector(iv):
    r"""
    Returns the permutation corresponding to inversion vector iv.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: permutation.from_inversion_vector([3,1,0,0,0])
        [3, 2, 4, 1, 5]
        sage: permutation.from_inversion_vector([2,3,6,4,0,2,2,1,0])
        [5, 9, 1, 8, 2, 6, 4, 7, 3]

    """

    p = [None] * len(iv)
    open_spots = range(len(iv))

    for i in range(len(iv)):
        if iv[i] != 0:
            p[open_spots[iv[i]]] = i+1
            open_spots.remove(open_spots[iv[i]])
        else:
            p[open_spots[0]] = i+1
            open_spots.remove(open_spots[0])
    return Permutation(p)




def from_cycles(n, cycles):
    r"""
    Returns the permutation corresponding to cycles.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: permutation.from_cycles(4, [[1,2]])
        [2, 1, 3, 4]

    """

    p = range(1,n+1)
    for cycle in cycles:
        first = cycle[0]
        for i in range(len(cycle)-1):
            p[cycle[i]-1] = cycle[i+1]
        p[cycle[-1]-1] = first
    return Permutation(p)

def from_lehmer_code(lehmer):
    r"""
    Returns the permutation with Lehmer code lehmer.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: Permutation([2,1,5,4,3]).to_lehmer_code()
        [1, 0, 2, 1, 0]
        sage: permutation.from_lehmer_code(_)
        [2, 1, 5, 4, 3]

    """

    n = len(lehmer)
    perm = [None] * n

    #Convert the factoradic to a permutation
    temp = [None] * n
    for i in range(n):
        lehmer[i] += 1
        temp[i] = lehmer[i]

    perm[n-1] = 1
    for i in reversed(range(n-1)):
        perm[i] = temp[i]
        for j in range(i+1, n):
            if perm[j] >= perm[i]:
                perm[j] += 1

    return Permutation([ p for p in perm ])

def from_reduced_word(rw):
    r"""
    Returns the permutation corresponding to the reduced
    word rw.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: permutation.from_reduced_word([3,2,3,1,2,3,1])
        [3, 4, 2, 1]
        sage: permutation.from_reduced_word([])
        []
    """
    if rw == []:
        return []

    p = [i+1 for i in range(max(rw)+1)]

    for i in rw:
        (p[i-1], p[i]) = (p[i], p[i-1])

    return Permutation(p)


class StandardPermutations_descents(CombinatorialClass):
    def __init__(self, d, n):
        """
        TESTS:
            sage: P = Permutations(descents=([1,0,4,8],12))
            sage: P == loads(dumps(P))
            True
        """
        self.d = d
        self.n = n

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(descents=([1,0,4,8],12)))
            'Standard permutations of 12 with descents [1, 0, 4, 8]'
        """
        return "Standard permutations of %s with descents %s"%(self.n, self.d)

    __object_class = Permutation_class

    def first(self):
        """
        Returns the first permutation with descents d.

        EXAMPLES:
            sage: Permutations(descents=([1,0,4,8],12)).first()
            [3, 2, 1, 4, 6, 5, 7, 8, 10, 9, 11, 12]

        """
        return descents_composition_first(Composition(descents=(self.d,self.n)))


    def last(self):
        """
        Returns the last permutation with descents d.

        EXAMPLES:
            sage: Permutations(descents=([1,0,4,8],12)).last()
            [12, 11, 8, 9, 10, 4, 5, 6, 7, 1, 2, 3]
        """
        return descents_composition_last(Composition(descents=(self.d,self.n)))

    def list(self):
        """
        Returns a list of all the permutations that have the
        descents d.

        EXAMPLES:
             sage: Permutations(descents=([2,4,0],5)).list()
             [[2, 1, 4, 3, 5],
              [2, 1, 5, 3, 4],
              [3, 1, 4, 2, 5],
              [3, 1, 5, 2, 4],
              [4, 1, 3, 2, 5],
              [5, 1, 3, 2, 4],
              [4, 1, 5, 2, 3],
              [5, 1, 4, 2, 3],
              [3, 2, 4, 1, 5],
              [3, 2, 5, 1, 4],
              [4, 2, 3, 1, 5],
              [5, 2, 3, 1, 4],
              [4, 2, 5, 1, 3],
              [5, 2, 4, 1, 3],
              [4, 3, 5, 1, 2],
              [5, 3, 4, 1, 2]]
         """

        return descents_composition_list(Composition(descents=(self.d,self.n)))



def descents_composition_list(dc):
    """
    Returns a list of all the permutations that have a descent
    compositions dc.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: permutation.descents_composition_list([1,2,2])
        [[2, 1, 4, 3, 5],
         [2, 1, 5, 3, 4],
         [3, 1, 4, 2, 5],
         [3, 1, 5, 2, 4],
         [4, 1, 3, 2, 5],
         [5, 1, 3, 2, 4],
         [4, 1, 5, 2, 3],
         [5, 1, 4, 2, 3],
         [3, 2, 4, 1, 5],
         [3, 2, 5, 1, 4],
         [4, 2, 3, 1, 5],
         [5, 2, 3, 1, 4],
         [4, 2, 5, 1, 3],
         [5, 2, 4, 1, 3],
         [4, 3, 5, 1, 2],
         [5, 3, 4, 1, 2]]
    """
    return map(lambda p: p.inverse(), StandardPermutations_recoils(dc).list())

def descents_composition_first(dc):
    r"""
    Computes the smallest element of a descent class having
    a descent decomposition dc.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: permutation.descents_composition_first([1,1,3,4,3])
        [3, 2, 1, 4, 6, 5, 7, 8, 10, 9, 11, 12]
    """

    if not isinstance(dc, Composition_class):
        try:
            dc = Composition(dc)
        except TypeError:
            raise TypeError, "The argument must be of type Composition"

    cpl = [x for x in reversed(dc.conjugate())]
    res = []
    s = 0
    for i in range(len(cpl)):
        res += [s + cpl[i]-j for j in range(cpl[i])]
        s   += cpl[i]

    return Permutation(res)

def descents_composition_last(dc):
    r"""
    Returns the largest element of a descent class having
    a descent decomposition dc.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: permutation.descents_composition_last([1,1,3,4,3])
        [12, 11, 8, 9, 10, 4, 5, 6, 7, 1, 2, 3]

    """
    if not isinstance(dc, Composition_class):
        try:
            dc = Composition(dc)
        except TypeError:
            raise TypeError, "The argument must be of type Composition"
    s = 0
    res = []
    for i in reversed(range(len(dc))):
        res = [j for j in range(s+1,s+dc[i]+1)] + res
        s += dc[i]

    return Permutation(res)


class StandardPermutations_recoilsfiner(CombinatorialClass):
    def __init__(self, recoils):
        """
        TESTS:
            sage: P = Permutations(recoils_finer=[2,2])
            sage: P == loads(dumps(P))
            True
        """
        self.recoils = recoils

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(recoils_finer=[2,2]))
            'Standard permutations whose recoils composition is finer than [2, 2]'
        """
        return "Standard permutations whose recoils composition is finer than %s"%self.recoils

    __object_class = Permutation_class

    def list(self):
        """
        Returns a list of all of the permutations whose
        recoils composition is finer than recoils.

        EXAMPLES:
            sage: Permutations(recoils_finer=[2,2]).list()
            [[1, 2, 3, 4],
             [1, 3, 2, 4],
             [1, 3, 4, 2],
             [3, 1, 2, 4],
             [3, 1, 4, 2],
             [3, 4, 1, 2]]
        """
        recoils = self.recoils
        dag = DiGraph()

        #Add the nodes
        for i in range(1, sum(recoils)+1):
            dag.add_vertex(i)

        #Add the edges to guarantee a finer recoil composition
        pos = 1
        for part in recoils:
            for i in range(part-1):
                dag.add_edge(pos, pos+1)
                pos += 1
            pos += 1

        rcf = []
        for le in dag.topological_sort_generator():
            rcf.append(Permutation(le))
        return rcf


class StandardPermutations_recoilsfatter(CombinatorialClass):
    def __init__(self, recoils):
        """
        TESTS:
            sage: P = Permutations(recoils_fatter=[2,2])
            sage: P == loads(dumps(P))
            True
        """
        self.recoils = recoils

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(recoils_fatter=[2,2]))
            'Standard permutations whose recoils composition is fatter than [2, 2]'
        """
        return "Standard permutations whose recoils composition is fatter than %s"%self.recoils

    __object_class = Permutation_class

    def list(self):
        """
        Returns a list of all of the permutations whose
        recoils composition is fatter than recoils.

        EXAMPLES:
            sage: Permutations(recoils_fatter=[2,2]).list()
            [[1, 3, 2, 4],
             [1, 3, 4, 2],
             [1, 4, 3, 2],
             [3, 1, 2, 4],
             [3, 1, 4, 2],
             [3, 2, 1, 4],
             [3, 2, 4, 1],
             [3, 4, 1, 2],
             [3, 4, 2, 1],
             [4, 1, 3, 2],
             [4, 3, 1, 2],
             [4, 3, 2, 1]]
        """
        recoils = self.recoils
        dag = DiGraph()

        #Add the nodes
        for i in range(1, sum(recoils)+1):
            dag.add_vertex(i)

        #Add the edges to guarantee a fatter recoil composition
        pos = 0
        for i in range(len(recoils)-1):
            pos += recoils[i]
            dag.add_edge(pos+1, pos)


        rcf = []
        for le in dag.topological_sort_generator():
            rcf.append(Permutation(le))
        return rcf

class StandardPermutations_recoils(CombinatorialClass):
    def __init__(self, recoils):
        """
        TESTS:
            sage: P = Permutations(recoils=[2,2])
            sage: P == loads(dumps(P))
            True
        """
        self.recoils = recoils


    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(recoils=[2,2]))
            'Standard permutations whose recoils composition is [2, 2]'
        """
        return "Standard permutations whose recoils composition is %s"%self.recoils

    __object_class = Permutation_class


    def list(self):
        """
        Returns a list of all of the permutations whose
        recoils composition is equal to recoils.

        EXAMPLES:
            sage: Permutations(recoils=[2,2]).list()
            [[1, 3, 2, 4], [1, 3, 4, 2], [3, 1, 2, 4], [3, 1, 4, 2], [3, 4, 1, 2]]
        """

        recoils = self.recoils
        dag = DiGraph()

        #Add all the nodes
        for i in range(1, sum(recoils)+1):
            dag.add_vertex(i)

        #Add the edges which guarantee a finer recoil comp.
        pos = 1
        for part in recoils:
            for i in range(part-1):
                dag.add_edge(pos, pos+1)
                pos += 1
            pos += 1

        #Add the edges which guarantee a fatter recoil comp.
        pos = 0
        for i in range(len(recoils)-1):
            pos += recoils[i]
            dag.add_edge(pos+1, pos)

        rcf = []
        for le in dag.topological_sort_generator():
            rcf.append(Permutation(le))
        return rcf



def from_major_code(mc, final_descent=False):
    r"""
    Returns the permutation corresponding to major code mc.

    REFERENCES:
        Skandera, M. 'An Eulerian Partner for Inversions', Sem. Lothar. Combin. 46 (2001) B46d.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: permutation.from_major_code([5, 0, 1, 0, 1, 2, 0, 1, 0])
        [9, 3, 5, 7, 2, 1, 4, 6, 8]
        sage: permutation.from_major_code([8, 3, 3, 1, 4, 0, 1, 0, 0])
        [2, 8, 4, 3, 6, 7, 9, 5, 1]
        sage: Permutation([2,1,6,4,7,3,5]).to_major_code()
        [3, 2, 0, 2, 2, 0, 0]
        sage: permutation.from_major_code([3, 2, 0, 2, 2, 0, 0])
        [2, 1, 6, 4, 7, 3, 5]

    """
    #define w^(n) to be the one-letter word n
    w = [len(mc)]

    #for i=n-1,..,1 let w^i be the unique word obtained by inserting
    #the letter i into the word w^(i+1) in such a way that
    #maj(w^i)-maj(w^(i+1)) = mc[i]
    for i in reversed(range(1,len(mc))):
        #Lemma 2.2 in Skandera

        #Get the descents of w and place them in reverse order
        d = Permutation(w).descents()
        d.reverse()

        #a is the list of all positions which are not descents
        a = filter(lambda x: x not in d, range(len(w)))

        #k is the number of desecents
        k = len(d)

        #d_k = -1    -- 0 in the lemma, but -1 due to 0-based indexing
        d.append(-1)


        l = mc[i-1]


        indices = d + a
        w.insert(indices[l]+1, i)

    #pi =
    return Permutation(w)


class StandardPermutations_bruhat_smaller(CombinatorialClass):
    def __init__(self, p):
        """
        TESTS:
            sage: P = Permutations(bruhat_smaller=[3,2,1])
            sage: P == loads(dumps(P))
            True
        """
        self.p = p

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(bruhat_smaller=[3,2,1]))
            'Standard permutations that are less than or equal to [3, 2, 1] in the Bruhat order'
        """
        return "Standard permutations that are less than or equal to %s in the Bruhat order"%self.p

    def list(self):
        r"""
        Returns a list of permutations smaller than or equal to p in the
        Bruhat order.

        EXAMPLES:
            sage: Permutations(bruhat_smaller=[4,1,2,3]).list()
            [[1, 2, 3, 4],
             [1, 2, 4, 3],
             [1, 3, 2, 4],
             [1, 4, 2, 3],
             [2, 1, 3, 4],
             [2, 1, 4, 3],
             [3, 1, 2, 4],
             [4, 1, 2, 3]]
        """
        return transitive_ideal(lambda x: x.bruhat_pred(), self.p)



class StandardPermutations_bruhat_greater(CombinatorialClass):
    def __init__(self, p):
        """
        TESTS:
            sage: P = Permutations(bruhat_greater=[3,2,1])
            sage: P == loads(dumps(P))
            True
        """
        self.p = p

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(bruhat_greater=[3,2,1]))
            'Standard permutations that are greater than or equal to [3, 2, 1] in the Bruhat order'
        """
        return "Standard permutations that are greater than or equal to %s in the Bruhat order"%self.p

    def list(self):
        r"""
        Returns a list of permutations greater than or equal to p in the
        Bruhat order.

        EXAMPLES:
            sage: Permutations(bruhat_greater=[4,1,2,3]).list()
            [[4, 1, 2, 3],
             [4, 1, 3, 2],
             [4, 2, 1, 3],
             [4, 2, 3, 1],
             [4, 3, 1, 2],
             [4, 3, 2, 1]]
        """
        return transitive_ideal(lambda x: x.bruhat_succ(), self.p)


################
# Bruhat Order #
################

def bruhat_lequal(p1, p2):
    r"""
    Returns True if p1 is less than p2in the Bruhat order.

    Algorithm from mupad-combinat.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: permutation.bruhat_lequal([2,4,3,1],[3,4,2,1])
        True
    """

    n1 = len(p1)
    n2 = len(p2)

    if n1 == 0:
        return True

    if p1[0] > p2[0] or p1[n1-1] < p2[n1-1]:
        return False

    for i in range(n1):
        c = 0
        for j in range(n1):
            if p2[j] > i+1:
                c += 1
            if p1[j] > i+1:
                c -= 1
            if c < 0:
                return False

    return True



#################
# Permutohedron #
#################

def permutohedron_lequal(p1, p2, side="right"):
    r"""
    Returns True if p1 is less than p2in the permutohedron order.

    By default, the computations are done in the right permutohedron.
    If you pass the option side='left', then they will be done in the
    left permutohedron.


    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: permutation.permutohedron_lequal(Permutation([3,2,1,4]),Permutation([4,2,1,3]))
        False
        sage: permutation.permutohedron_lequal(Permutation([3,2,1,4]),Permutation([4,2,1,3]), side='left')
        True


    """

    n1 = len(p1)
    n2 = len(p2)

    l1 = p1.number_of_inversions()
    l2 = p2.number_of_inversions()

    if l1 > l2:
        return Fal

    if side == "right":
        prod = p1._left_to_right_multiply_on_right(p2.inverse())
    else:
        prod = p1._left_to_right_multiply_on_left(p2.inverse())


    return prod.number_of_inversions() == l2 - l1


############
# Patterns #
############

def to_standard(p):
    r"""
    Returns a standard permutation corresponding to the
    permutation p.

    EXAMPLES:
        sage: import sage.combinat.permutation as permutation
        sage: permutation.to_standard([4,2,7])
        [2, 1, 3]
        sage: permutation.to_standard([1,2,3])
        [1, 2, 3]
    """

    s = p[:]
    biggest = max(p) + 1
    i = 1
    for j in range(len(p)):
        smallest = min(p)
        smallest_index = p.index(smallest)
        s[smallest_index] = i
        i += 1
        p[smallest_index] = biggest

    return Permutation(s)



##########################################################


def CyclicPermutations(mset):
    """
    Returns the combinatorial class of all cyclic permutations of mset
    in cycle notation.  These are the same as necklaces.

    EXAMPLES:
        sage: CyclicPermutations(range(4)).list()
        [[0, 1, 2, 3],
         [0, 1, 3, 2],
         [0, 2, 1, 3],
         [0, 2, 3, 1],
         [0, 3, 1, 2],
         [0, 3, 2, 1]]
        sage: CyclicPermutations([1,1,1]).list()
        [[1, 1, 1]]
    """
    return CyclicPermutations_mset(mset)

class CyclicPermutations_mset(CombinatorialClass):
    def __init__(self, mset):
        """
        TESTS:
            sage: CP = CyclicPermutations(range(4))
            sage: CP == loads(dumps(CP))
            True
        """
        self.mset = mset

    def __repr__(self):
        """
        TESTS:
            sage: repr(CyclicPermutations(range(4)))
            'Cyclic permutations of [0, 1, 2, 3]'
        """
        return "Cyclic permutations of %s"%self.mset

    def list(self, distinct=False):
        """
        EXAMPLES:
            sage: CyclicPermutations(range(4)).list()
            [[0, 1, 2, 3],
             [0, 1, 3, 2],
             [0, 2, 1, 3],
             [0, 2, 3, 1],
             [0, 3, 1, 2],
             [0, 3, 2, 1]]
        """
        return list(self.iterator(distinct=distinct))

    def iterator(self, distinct=False):
        """
        EXAMPLES:
            sage: CyclicPermutations(range(4)).list()
            [[0, 1, 2, 3],
             [0, 1, 3, 2],
             [0, 2, 1, 3],
             [0, 2, 3, 1],
             [0, 3, 1, 2],
             [0, 3, 2, 1]]
             sage: CyclicPermutations([1,1,1]).list()
             [[1, 1, 1]]
             sage: CyclicPermutations([1,1,1]).list(distinct=True)
             [[1, 1, 1], [1, 1, 1]]
        """
        if distinct:
            content = [1]*len(self.mset)
        else:
            content = [0]*len(self.mset)
            index_list = map(self.mset.index, self.mset)
            for i in index_list:
                content[i] += 1

        for necklace in Necklaces(content):
            yield [self.mset[x-1] for x in necklace]

##########################################3

def CyclicPermutationsOfPartition(partition):
    """
    Returns the combinatorial class of all combinations of cyclic
    permutations of each cell of the partition.  This is the same
    as a Cartesian product of necklaces.

    EXAMPLES:
        sage: CyclicPermutationsOfPartition([[1,2,3,4],[5,6,7]]).list()
        [[[1, 2, 3, 4], [5, 6, 7]],
         [[1, 2, 4, 3], [5, 6, 7]],
         [[1, 3, 2, 4], [5, 6, 7]],
         [[1, 3, 4, 2], [5, 6, 7]],
         [[1, 4, 2, 3], [5, 6, 7]],
         [[1, 4, 3, 2], [5, 6, 7]],
         [[1, 2, 3, 4], [5, 7, 6]],
         [[1, 2, 4, 3], [5, 7, 6]],
         [[1, 3, 2, 4], [5, 7, 6]],
         [[1, 3, 4, 2], [5, 7, 6]],
         [[1, 4, 2, 3], [5, 7, 6]],
         [[1, 4, 3, 2], [5, 7, 6]]]

        sage: CyclicPermutationsOfPartition([[1,2,3,4],[4,4,4]]).list()
        [[[1, 2, 3, 4], [4, 4, 4]],
         [[1, 2, 4, 3], [4, 4, 4]],
         [[1, 3, 2, 4], [4, 4, 4]],
         [[1, 3, 4, 2], [4, 4, 4]],
         [[1, 4, 2, 3], [4, 4, 4]],
         [[1, 4, 3, 2], [4, 4, 4]]]

        sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list()
        [[[1, 2, 3], [4, 4, 4]], [[1, 3, 2], [4, 4, 4]]]

        sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list(distinct=True)
        [[[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]],
         [[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]]]
    """
    return CyclicPermutationsOfPartition_partition(partition)

class CyclicPermutationsOfPartition_partition(CombinatorialClass):
    def __init__(self, partition):
        """
        TESTS:
            sage: CP = CyclicPermutationsOfPartition([[1,2,3,4],[5,6,7]])
            sage: CP == loads(dumps(CP))
            True
        """
        self.partition = partition

    def __repr__(self):
        """
        TESTS:
            sage: repr(CyclicPermutationsOfPartition([[1,2,3,4],[5,6,7]]))
            'Cyclic permutations of partition [[1, 2, 3, 4], [5, 6, 7]]'
        """
        return "Cyclic permutations of partition %s"%self.partition

    def iterator(self, distinct=False):
        """
        AUTHOR: Robert Miller

        EXAMPLES:
            sage: CyclicPermutationsOfPartition([[1,2,3,4],[5,6,7]]).list()
            [[[1, 2, 3, 4], [5, 6, 7]],
             [[1, 2, 4, 3], [5, 6, 7]],
             [[1, 3, 2, 4], [5, 6, 7]],
             [[1, 3, 4, 2], [5, 6, 7]],
             [[1, 4, 2, 3], [5, 6, 7]],
             [[1, 4, 3, 2], [5, 6, 7]],
             [[1, 2, 3, 4], [5, 7, 6]],
             [[1, 2, 4, 3], [5, 7, 6]],
             [[1, 3, 2, 4], [5, 7, 6]],
             [[1, 3, 4, 2], [5, 7, 6]],
             [[1, 4, 2, 3], [5, 7, 6]],
             [[1, 4, 3, 2], [5, 7, 6]]]

            sage: CyclicPermutationsOfPartition([[1,2,3,4],[4,4,4]]).list()
            [[[1, 2, 3, 4], [4, 4, 4]],
             [[1, 2, 4, 3], [4, 4, 4]],
             [[1, 3, 2, 4], [4, 4, 4]],
             [[1, 3, 4, 2], [4, 4, 4]],
             [[1, 4, 2, 3], [4, 4, 4]],
             [[1, 4, 3, 2], [4, 4, 4]]]

            sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list()
            [[[1, 2, 3], [4, 4, 4]], [[1, 3, 2], [4, 4, 4]]]

            sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list(distinct=True)
            [[[1, 2, 3], [4, 4, 4]],
             [[1, 3, 2], [4, 4, 4]],
             [[1, 2, 3], [4, 4, 4]],
             [[1, 3, 2], [4, 4, 4]]]
        """

        if len(self.partition) == 1:
            for i in CyclicPermutations_mset(self.partition[0]).iterator(distinct=distinct):
                yield [i]
        else:
            for right in CyclicPermutationsOfPartition_partition(self.partition[1:]).iterator(distinct=distinct):
                for perm in CyclicPermutations_mset(self.partition[0]).iterator(distinct=distinct):
                    yield [perm] + right


    def list(self, distinct=False):
        """
        EXAMPLES:
            sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list()
            [[[1, 2, 3], [4, 4, 4]], [[1, 3, 2], [4, 4, 4]]]
            sage: CyclicPermutationsOfPartition([[1,2,3],[4,4,4]]).list(distinct=True)
            [[[1, 2, 3], [4, 4, 4]],
             [[1, 3, 2], [4, 4, 4]],
             [[1, 2, 3], [4, 4, 4]],
             [[1, 3, 2], [4, 4, 4]]]
        """

        return list(self.iterator(distinct=distinct))



######
#Avoiding


class StandardPermutations_avoiding_12(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: p = Permutations(3, avoiding=[1,2])
            sage: p == loads(dumps(p))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(3, avoiding=[1,2]))
            'Standard permutations of 3 avoiding [1, 2]'
        """
        return "Standard permutations of %s avoiding [1, 2]"%self.n

    def list(self):
        """
        EXAMPLES:
            sage: Permutations(3, avoiding=[1,2]).list()
            [[3, 2, 1]]
        """
        return [Permutation_class(range(self.n, 0, -1))]

class StandardPermutations_avoiding_21(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: p = Permutations(3, avoiding=[2,1])
            sage: p == loads(dumps(p))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(3, avoiding=[2,1]))
            'Standard permutations of 3 avoiding [2, 1]'
        """
        return "Standard permutations of %s avoiding [2, 1]"%self.n

    def list(self):
        """
        EXAMPLES:
            sage: Permutations(3, avoiding=[2,1]).list()
            [[1, 2, 3]]
        """
        return [Permutation_class(range(1, self.n+1))]


class StandardPermutations_avoiding_132(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: p = Permutations(3, avoiding=[1,3,2])
            sage: p == loads(dumps(p))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(3, avoiding=[1,3,2]))
            'Standard permutations of 3 avoiding [1, 3, 2]'
        """
        return "Standard permutations of %s avoiding [1, 3, 2]"%self.n

    def count(self):
        """
        EXAMPLES:
            sage: Permutations(5, avoiding=[1, 3, 2]).count()
            42
            sage: len( Permutations(5, avoiding=[1, 3, 2]).list() )
            42
        """
        return catalan_number(self.n)

    def iterator(self):
        """
        EXAMPLES:
            sage: Permutations(3, avoiding=[1,3,2]).list()
            [[1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
            sage: Permutations(4, avoiding=[1,3,2]).list()
            [[4, 1, 3, 2],
             [4, 2, 1, 3],
             [4, 2, 3, 1],
             [4, 3, 1, 2],
             [4, 3, 2, 1],
             [3, 4, 1, 2],
             [3, 4, 2, 1],
             [2, 3, 4, 1],
             [3, 2, 4, 1],
             [1, 3, 2, 4],
             [2, 1, 3, 4],
             [2, 3, 1, 4],
             [3, 1, 2, 4],
             [3, 2, 1, 4]]

        """
        if self.n == 0:
            return

        elif self.n < 3:
            for p in StandardPermutations_n(self.n):
                yield p
            return

        elif self.n == 3:
            for p in StandardPermutations_n(self.n):
                if p != [1, 2, 3]:
                    yield p
            return



        #Yield all the 132 avoiding permutations to the right.
        for right in StandardPermutations_avoiding_132(self.n - 1):
            yield Permutation_class([self.n] + list(right))

        #yi
        for i in range(1, self.n-1):
            for left in StandardPermutations_avoiding_132(i):
                for right in StandardPermutations_avoiding_132(self.n-i-1):
                    yield Permutation_class( map(lambda x: x+(self.n-i-1), left) + [self.n] + list(right) )


        #Yield all the 132 avoiding permutations to the left
        for left in StandardPermutations_avoiding_132(self.n - 1):
            yield Permutation_class(list(left) + [self.n])


class StandardPermutations_avoiding_123(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: p = Permutations(3, avoiding=[1, 2, 3])
            sage: p == loads(dumps(p))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(3, avoiding=[1, 2, 3]))
            'Standard permutations of 3 avoiding [1, 2, 3]'
        """
        return "Standard permutations of %s avoiding [1, 2, 3]"%self.n

    def count(self):
        """
        EXAMPLES:
            sage: Permutations(5, avoiding=[1, 2, 3]).count()
            42
            sage: len( Permutations(5, avoiding=[1, 2, 3]).list() )
            42
        """
        return catalan_number(self.n)

    def iterator(self):
        """
        EXAMPLES:
            sage: Permutations(3, avoiding=[1, 2, 3]).list()
             [[1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
            sage: Permutations(2, avoiding=[1, 2, 3]).list()
            [[1, 2], [2, 1]]
            sage: Permutations(3, avoiding=[1, 2, 3]).list()
            [[1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        """
        if self.n == 0:
            return

        elif self.n < 3:
            for p in StandardPermutations_n(self.n):
                yield p
            return

        elif self.n == 3:
            for p in StandardPermutations_n(self.n):
                if p != [1, 2, 3]:
                    yield p
            return


        for p in StandardPermutations_avoiding_132(self.n):
            #Convert p to a 123 avoiding permutation by
            m = self.n+1
            minima_pos = []
            minima = []
            for i in range(self.n):
                if p[i] < m:
                    minima_pos.append(i)
                    minima.append(p[i])
                    m = p[i]


            new_p = []
            non_minima = filter(lambda x: x not in minima, range(self.n, 0, -1))
            a = 0
            b = 0
            for i in range(self.n):
                if i in minima_pos:
                    new_p.append( minima[a] )
                    a += 1
                else:
                    new_p.append( non_minima[b] )
                    b += 1

            yield Permutation_class( new_p )


class StandardPermutations_avoiding_321(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: p = Permutations(3, avoiding=[3, 2, 1])
            sage: p == loads(dumps(p))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(3, avoiding=[3, 2, 1]))
            'Standard permutations of 3 avoiding [3, 2, 1]'
        """
        return "Standard permutations of %s avoiding [3, 2, 1]"%self.n

    def count(self):
        """
        EXAMPLES:
            sage: Permutations(5, avoiding=[3, 2, 1]).count()
            42
            sage: len( Permutations(5, avoiding=[3, 2, 1]).list() )
            42
        """
        return catalan_number(self.n)

    def iterator(self):
        """
        EXAMPLES:
            sage: Permutations(3, avoiding=[3, 2, 1]).list() #indirect doctest
            [[2, 3, 1], [3, 1, 2], [1, 3, 2], [2, 1, 3], [1, 2, 3]]
        """
        for p in StandardPermutations_avoiding_123(self.n):
            yield p.reverse()


class StandardPermutations_avoiding_231(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: p = Permutations(3, avoiding=[2, 3, 1])
            sage: p == loads(dumps(p))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(3, avoiding=[2, 3, 1]))
            'Standard permutations of 3 avoiding [2, 3, 1]'
        """
        return "Standard permutations of %s avoiding [2, 3, 1]"%self.n

    def count(self):
        """
        EXAMPLES:
            sage: Permutations(5, avoiding=[2, 3, 1]).count()
            42
            sage: len( Permutations(5, avoiding=[2, 3, 1]).list() )
            42
        """
        return catalan_number(self.n)

    def iterator(self):
        """
        EXAMPLES:
            sage: Permutations(3, avoiding=[2, 3, 1]).list()
            [[2, 3, 1], [3, 1, 2], [1, 3, 2], [2, 1, 3], [1, 2, 3]]

        """
        for p in StandardPermutations_avoiding_132(self.n):
            yield p.reverse()


class StandardPermutations_avoiding_312(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: p = Permutations(3, avoiding=[3, 1, 2])
            sage: p == loads(dumps(p))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(3, avoiding=[3, 1, 2]))
            'Standard permutations of 3 avoiding [3, 1, 2]'
        """
        return "Standard permutations of %s avoiding [3, 1, 2]"%self.n

    def count(self):
        """
        EXAMPLES:
            sage: Permutations(5, avoiding=[3, 1, 2]).count()
            42
            sage: len( Permutations(5, avoiding=[3, 1, 2]).list() )
            42
        """
        return catalan_number(self.n)

    def iterator(self):
        """
        EXAMPLES:
            sage: Permutations(3, avoiding=[3, 1, 2]).list()
            [[3, 1, 2], [2, 3, 1], [2, 1, 3], [1, 3, 2], [1, 2, 3]]

        """
        for p in StandardPermutations_avoiding_132(self.n):
            yield p.complement()


class StandardPermutations_avoiding_213(CombinatorialClass):
    def __init__(self, n):
        """
        TESTS:
            sage: p = Permutations(3, avoiding=[2, 1, 3])
            sage: p == loads(dumps(p))
            True
        """
        self.n = n

    def __repr__(self):
        """
        TESTS:
            sage: repr(Permutations(3, avoiding=[2, 1, 3]))
            'Standard permutations of 3 avoiding [2, 1, 3]'
        """
        return "Standard permutations of %s avoiding [2, 1, 3]"%self.n

    def count(self):
        """
        EXAMPLES:
            sage: Permutations(5, avoiding=[2, 1, 3]).count()
            42
            sage: len( Permutations(5, avoiding=[2, 1, 3]).list() )
            42
        """
        return catalan_number(self.n)

    def iterator(self):
        """
        EXAMPLES:
            sage: Permutations(3, avoiding=[2, 1, 3]).list()
            [[2, 1, 3], [1, 3, 2], [3, 1, 2], [2, 3, 1], [3, 2, 1]]

        """
        for p in StandardPermutations_avoiding_132(self.n):
            yield p.complement().reverse()


class StandardPermutations_avoiding_generic(CombinatorialClass):
    def __init__(self, n, a):
        """
        EXAMPLES:
            sage: P = Permutations(3, avoiding=[[2, 1, 3],[1,2,3]])
            sage: P == loads(dumps(P))
            True
            sage: type(P)
            <class 'sage.combinat.permutation.StandardPermutations_avoiding_generic'>
        """
        self.n = n
        self.a = a

    def __repr__(self):
        """
        EXAMPLES:
            sage: P = Permutations(3, avoiding=[[2, 1, 3],[1,2,3]])
            sage: P.__repr__()
            'Standard permutations of 3 avoiding [[2, 1, 3], [1, 2, 3]]'

        """
        return "Standard permutations of %s avoiding %s"%(self.n, self.a)

    def iterator(self):
        """
        Note that this uses an exteremely inefficient algorithm and should be
        improved.

        EXAMPLES:
            sage: Permutations(3, avoiding=[[2, 1, 3],[1,2,3]]).list()
            [[1, 3, 2], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

        """
        for p in StandardPermutations_n(self.n):
            ls = map(len, self.a)
            found = False
            for l in ls:
                for pos in subword.Subwords(range(self.n),l):
                    if to_standard(map(lambda z: p[z] , pos)) in self.a:
                        found = True
                        break
                if found:
                    break

            if found:
                continue
            else:
                yield p




def Permutations(n=None,k=None, **kwargs):
    """
    Returns a combinatorial class of permutations.

    Permutations(n) returns the class of permutations of n, if n is an
    integer, list, set, or string.

    Permutations(n, k) returns the class of permutations of n (where n
    is any of the above things) of length k; k must be an integer.

    Valid keyword arguments are: 'descents', 'bruhat_smaller',
    'bruhat_greater', 'recoils_finer', 'recoils_fatter', 'recoils', and
    'avoiding'. With the exception of 'avoiding', you cannot specify n
    or k along with a keyword.

    Permutations(descents=list) returns the class of permutations with
    descents in the positions specified by `list'.

    Permutations(bruhat_{smaller,greater}=p) returns the class of
    permutations smaller or greater, respectively, than the given
    permutation in Bruhat order.

    Permutations(recoils=p) returns the class of permutations whose
    recoils composition is p.

    Permutations(recoils_{fatter,finer}=p) returns the class of
    permutations whose recoils composition is fatter or finer,
    respectively, than the given permutation.

    Permutations(n, avoiding=P) returns the class of permutations of n
    avoiding P. Here P may be a single permutation or a list of
    permutations; the returned class will avoid all patterns in P.

    EXAMPLES:

        sage: p = Permutations(3); p
        Standard permutations of 3
        sage: p.list()
        [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

        sage: p = Permutations(3, 2); p
        Permutations of {1,...,3} of length 2
        sage: p.list()
        [[1, 2], [1, 3], [2, 1], [2, 3], [3, 1], [3, 2]]

        sage: p = Permutations(['c', 'a', 't']); p
        Permutations of the set ['c', 'a', 't']
        sage: p.list()
        [['c', 'a', 't'],
         ['c', 't', 'a'],
         ['a', 'c', 't'],
         ['a', 't', 'c'],
         ['t', 'c', 'a'],
         ['t', 'a', 'c']]

        sage: p = Permutations(['c', 'a', 't'], 2); p
        Permutations of the set ['c', 'a', 't'] of length 2
        sage: p.list()
        [['c', 'a'], ['c', 't'], ['a', 'c'], ['a', 't'], ['t', 'c'], ['t', 'a']]


        sage: p = Permutations([1,1,2]); p
        Permutations of the multi-set [1, 1, 2]
        sage: p.list()
        [[1, 1, 2], [1, 2, 1], [2, 1, 1]]


        sage: p = Permutations([1,1,2], 2); p
        Permutations of the multi-set [1, 1, 2] of length 2
        sage: p.list()
        [[1, 1], [1, 2], [2, 1]]

        sage: p = Permutations(descents=[1,3]); p
        Standard permutations of 4 with descents [1, 3]
        sage: p.list()
        [[1, 3, 2, 4], [1, 4, 2, 3], [2, 3, 1, 4], [2, 4, 1, 3], [3, 4, 1, 2]]

        sage: p = Permutations(bruhat_smaller=[1,3,2,4]); p
        Standard permutations that are less than or equal to [1, 3, 2, 4] in the Bruhat order
        sage: p.list()
        [[1, 2, 3, 4], [1, 3, 2, 4]]

        sage: p = Permutations(bruhat_greater=[4,2,3,1]); p
        Standard permutations that are greater than or equal to [4, 2, 3, 1] in the Bruhat order
        sage: p.list()
        [[4, 2, 3, 1], [4, 3, 2, 1]]

        sage: p = Permutations(recoils_finer=[2,1]); p
        Standard permutations whose recoils composition is finer than [2, 1]
        sage: p.list()
        [[1, 2, 3], [1, 3, 2], [3, 1, 2]]

        sage: p = Permutations(recoils_fatter=[2,1]); p
        Standard permutations whose recoils composition is fatter than [2, 1]
        sage: p.list()
        [[1, 3, 2], [3, 1, 2], [3, 2, 1]]

        sage: p = Permutations(recoils=[2,1]); p
        Standard permutations whose recoils composition is [2, 1]
        sage: p.list()
        [[1, 3, 2], [3, 1, 2]]

        sage: p = Permutations(4, avoiding=[1,3,2]); p
        Standard permutations of 4 avoiding [1, 3, 2]
        sage: p.list()
        [[4, 1, 3, 2],
         [4, 2, 1, 3],
         [4, 2, 3, 1],
         [4, 3, 1, 2],
         [4, 3, 2, 1],
         [3, 4, 1, 2],
         [3, 4, 2, 1],
         [2, 3, 4, 1],
         [3, 2, 4, 1],
         [1, 3, 2, 4],
         [2, 1, 3, 4],
         [2, 3, 1, 4],
         [3, 1, 2, 4],
         [3, 2, 1, 4]]

        sage: p = Permutations(5, avoiding=[[3,4,1,2], [4,2,3,1]]); p
        Standard permutations of 5 avoiding [[3, 4, 1, 2], [4, 2, 3, 1]]
        sage: p.count()
        88
        sage: p.random()
        [1, 3, 4, 5, 2]

    """

    valid_args = ['descents', 'bruhat_smaller', 'bruhat_greater',
                  'recoils_finer', 'recoils_fatter', 'recoils', 'avoiding']

    number_of_arguments = 0
    if n is not None:
            number_of_arguments += 1
    else:
        if k is not None:
            number_of_arguments += 1


    #Make sure that exactly one keyword was passed
    for key in kwargs:
        if key not in valid_args:
            raise ValueError, "unknown keyword argument: %s"%key
        if key not in [ 'avoiding' ]:
            number_of_arguments += 1

    if number_of_arguments == 0:
        return StandardPermutations_all()

    if number_of_arguments != 1:
        raise ValueError, "you must specify exactly one argument"

    if n is not None:
        if isinstance(n, (int, Integer)):
            if k is None:
                if 'avoiding' in kwargs:
                    a = kwargs['avoiding']
                    if a in StandardPermutations_all():
                        if a == [1,2]:
                            return StandardPermutations_avoiding_12(n)
                        elif a == [2,1]:
                            return StandardPermutations_avoiding_21(n)
                        elif a == [1,2,3]:
                            return StandardPermutations_avoiding_123(n)
                        elif a == [1,3,2]:
                            return StandardPermutations_avoiding_132(n)
                        elif a == [2,1,3]:
                            return StandardPermutations_avoiding_213(n)
                        elif a == [2,3,1]:
                            return StandardPermutations_avoiding_231(n)
                        elif a == [3,1,2]:
                            return StandardPermutations_avoiding_312(n)
                        elif a == [3,2,1]:
                            return StandardPermutations_avoiding_321(n)
                        else:
                            return StandardPermutations_avoiding_generic(n, [a])
                    elif isinstance(a, __builtin__.list):
                        return StandardPermutations_avoiding_generic(n, a)
                    else:
                        raise ValueError, "do not know how to avoid %s"%a
                else:
                    return StandardPermutations_n(n)
            else:
                return Permutations_nk(n,k)
        else:
            #In this case, we have that n is a list
            if map(n.index, n) == range(len(n)):
                if k is None:
                    return Permutations_set(n)
                else:
                    return Permutations_setk(n,k)
            else:
                if k is None:
                    return Permutations_mset(n)
                else:
                    return Permutations_msetk(n,k)
    elif 'descents' in kwargs:
        if isinstance(kwargs['descents'], tuple):
            return StandardPermutations_descents(*kwargs['descents'])
        else:
            return StandardPermutations_descents(kwargs['descents'], max(kwargs['descents'])+1)
    elif 'bruhat_smaller' in kwargs:
        return StandardPermutations_bruhat_smaller(Permutation(kwargs['bruhat_smaller']))
    elif 'bruhat_greater' in kwargs:
        return StandardPermutations_bruhat_greater(Permutation(kwargs['bruhat_greater']))
    elif 'recoils_finer' in kwargs:
        return StandardPermutations_recoilsfiner(kwargs['recoils_finer'])
    elif 'recoils_fatter' in kwargs:
        return StandardPermutations_recoilsfatter(kwargs['recoils_fatter'])
    elif 'recoils' in kwargs:
        return StandardPermutations_recoils(kwargs['recoils'])
