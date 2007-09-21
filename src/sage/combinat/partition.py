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
from sage.rings.all import QQ, RR, ZZ
from sage.misc.all import prod, sage_eval
from sage.rings.arith import factorial, gcd
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.calculus.calculus import Function_ceil, log, SymbolicVariable, Function_sinh
from sage.calculus.functional import diff
import sage.combinat.generator as generator
import sage.combinat.misc as misc
import sage.combinat.skew_partition
from sage.rings.integer import Integer
from UserList import UserList
import __builtin__
from sage.functions.constants import pi
ceil = Function_ceil()
sinh = Function_sinh()
from combinat import CombinatorialClass, CombinatorialObject
import partitions as partitions_ext
from sage.libs.all import pari

def Partition(l=None, exp=None, core_and_quotient=None):
    """
    Returns a partition object.

    EXAMPLES:
        sage: Partition(exp=[2,1,1])
        [3, 2, 1, 1]
        sage: Partition(core_and_quotient=([2,1], [[2,1],[3],[1,1,1]]))
        [11, 5, 5, 3, 2, 2, 2]
        sage: Partition([3,2,1])
        [3, 2, 1]
    """
    number_of_arguments = 0
    for arg in ['l', 'exp', 'core_and_quotient']:
        if eval(arg) != None:
            number_of_arguments += 1

    if number_of_arguments != 1:
        raise ValueError, "you must specify exactly one argument"

    if l != None:
        if sum(l) == 0:
            l = []
        return Partition_class(l)
    elif exp != None:
        return from_exp(exp)
    else:
        return from_core_and_quotient(*core_and_quotient)


def from_exp(a):
    """
    Returns a partition from its list of multiplicities.

    EXAMPLES:
        sage: partition.from_exp([1,2,1])
        [3, 2, 2, 1]
    """

    p = []
    for i in reversed(range(len(a))):
        p += [i+1]*a[i]

    return Partition(p)



def from_core_and_quotient(core, quotient):
    """
    Returns a partition from its r-core and r-quotient.

    Algorithm from mupad-combinat.

    EXAMPLES:
        sage: partition.from_core_and_quotient([2,1], [[2,1],[3],[1,1,1]])
        [11, 5, 5, 3, 2, 2, 2]
    """
    length = len(quotient)
    k = length*max( [ceil(len(core)/length),len(core)] + [len(q) for q in quotient] ) + length
    v = [ core[i]-(i+1)+1 for i in range(len(core))] + [ -i + 1 for i in range(len(core)+1,k+1) ]
    w = [ filter(lambda x: (x-i) % length == 0, v) for i in range(1, length+1) ]
    new_w = []
    for i in range(length):
        new_w += [ w[i][j] + length*quotient[i][j] for j in range(len(quotient[i]))]
        new_w += [ w[i][j] for j in range(len(quotient[i]), len(w[i])) ]
    new_w.sort()
    new_w.reverse()
    return filter(lambda x: x != 0, [new_w[i-1]+i-1 for i in range(1, len(new_w)+1)])


class Partition_class(CombinatorialObject):
    def ferrers_diagram(self):
        """
        Return the Ferrers diagram of pi.

        INPUT:
            pi -- a partition, given as a list of integers.

        EXAMPLES:
            sage: print Partition([5,5,2,1]).ferrers_diagram()
            *****
            *****
            **
            *
            sage: pi = Partitions(10).list()[11] ## [6,1,1,1,1]
            sage: print pi.ferrers_diagram()
            ******
            *
            *
            *
            *
            sage: pi = Partitions(10).list()[8] ## [6, 3, 1]
            sage: print pi.ferrers_diagram()
            ******
            ***
            *
        """
        return '\n'.join(['*'*p for p in self])


    def __div__(self, p):
        """
        Returns the skew partition self/p.

        EXAMPLES:
            sage: p = Partition([3,2,1])
            sage: p/[1,1]
            [[3, 2, 1], [1, 1]]
            sage: p/[3,2,1]
            [[3, 2, 1], [3, 2, 1]]

        """
        if not self.dominates(Partition_class(p)):
            raise ValueError, "the partition must dominate p"

        return sage.combinat.skew_partition.SkewPartition([self[:], p])

    def power(self,k):
        """
        partition_power( pi, k ) returns the partition corresponding to the
        $k$-th power of a permutation with cycle structure pi
        (thus describes the powermap of symmetric groups).

        Wraps GAP's PowerPartition.

        EXAMPLES:
            sage: p = Partition([5,3])
            sage: p.power(1)
            [5, 3]
            sage: p.power(2)
            [5, 3]
            sage: p.power(3)
            [5, 1, 1, 1]
            sage: p.power(4)
            [5, 3]

         Now let us compare this to the power map on $S_8$:

            sage: G = SymmetricGroup(8)
            sage: g = G([(1,2,3,4,5),(6,7,8)])
            sage: g
            (1,2,3,4,5)(6,7,8)
            sage: g^2
            (1,3,5,2,4)(6,8,7)
            sage: g^3
            (1,4,2,5,3)
            sage: g^4
            (1,5,4,3,2)(6,7,8)

        """
        ans=gap.eval("PowerPartition(%s,%s)"%(self,ZZ(k)))
        return Partition_class(eval(ans))

    def next(self):
        """
        Returns the partition that lexicographically follows
        the partition p.  If p is the last partition, then
        it returns False.

        EXAMPLES:
            sage: Partition([4]).next()
            [3, 1]
            sage: Partition([1,1,1,1]).next()
            False
        """
        p = self
        n = 0
        m = 0
        for i in p:
            n += i
            m += 1

        next_p = p[:] + [1]*(n - len(p))

        #Check to see if we are at the last (all ones) partition
        if p == [1]*n:
            return False

        #
        #If we are not, then run the ZS1 algorithm.
        #

        #Let h be the number of non-one  entries in the
        #partition
        h = 0
        for i in next_p:
            if i != 1:
                h += 1

        if next_p[h-1] == 2:
            m += 1
            next_p[h-1] = 1
            h -= 1
        else:
            r = next_p[h-1] - 1
            t = m - h + 1
            next_p[h-1] = r

            while t >= r :
                h += 1
                next_p[h-1] = r
                t -= r

            if t == 0:
                m = h
            else:
                m = h + 1
                if t > 1:
                    h += 1
                    next_p[h-1] = t

        return Partition_class(next_p[:m])

    def size(self):
        """
        Returns the size of partition p.

        EXAMPLES:
            sage: Partition([2,2]).size()
            4
            sage: Partition([3,2,1]).size()
            6
        """
        return sum(self)

    def sign(self):
        r"""
        partition_sign( pi ) returns the sign of a permutation with cycle structure
        given by the partition pi.

        This function corresponds to a homomorphism from the symmetric group
        $S_n$ into the cyclic group of order 2, whose kernel is exactly the
        alternating group $A_n$. Partitions of sign $1$ are called {\it even partitions}
        while partitions of sign $-1$ are called {\it odd}.

        Wraps GAP's SignPartition.

        EXAMPLES:
            sage: Partition([5,3]).sign()
            1
            sage: Partition([5,2]).sign()
            -1

        {\it Zolotarev's lemma} states that the Legendre symbol
        $ \left(\frac{a}{p}\right)$ for an integer $a \pmod p$ ($p$ a prime number),
        can be computed as sign(p_a), where sign denotes the sign of a permutation
        and p_a the permutation of the residue classes $\pmod p$ induced by
        modular multiplication by $a$, provided $p$ does not divide $a$.

        We verify this in some examples.

            sage: F = GF(11)
            sage: a = F.multiplicative_generator();a
            2
            sage: plist = [int(a*F(x)) for x in range(1,11)]; plist
            [2, 4, 6, 8, 10, 1, 3, 5, 7, 9]

        This corresponds ot the permutation (1, 2, 4, 8, 5, 10, 9, 7, 3, 6)
        (acting the set $\{1,2,...,10\}$) and to the partition [10].

            sage: p = PermutationGroupElement('(1, 2, 4, 8, 5, 10, 9, 7, 3, 6)')
            sage: p.sign()
            -1
            sage: Partition([10]).sign()
            -1
            sage: kronecker_symbol(11,2)
            -1

        Now replace $2$ by $3$:

            sage: plist = [int(F(3*x)) for x in range(1,11)]; plist
            [3, 6, 9, 1, 4, 7, 10, 2, 5, 8]
            sage: range(1,11)
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            sage: p = PermutationGroupElement('(3,4,8,7,9)')
            sage: p.sign()
            1
            sage: kronecker_symbol(3,11)
            1
            sage: Partition([5,1,1,1,1,1]).sign()
            1

        In both cases, Zolotarev holds.

        REFERENCES:
            http://en.wikipedia.org/wiki/Zolotarev's_lemma
        """
        ans=gap.eval("SignPartition(%s)"%(self))
        return sage_eval(ans)


    def up(self):
        r"""
        Returns a generator for partitions that can be obtained from pi
        by adding a box.

        EXAMPLES:
            sage: [p for p in Partition([2,1,1]).up()]
            [[3, 1, 1], [2, 2, 1], [2, 1, 1, 1]]
            sage: [p for p in Partition([3,2]).up()]
            [[4, 2], [3, 3], [3, 2, 1]]
        """
        p = self
        previous = p[0] + 1
        for i, current in enumerate(p):
            if current < previous:
                yield Partition(p[:i] + [ p[i] + 1 ] + p[i+1:])
            previous = current
        else:
            yield Partition(p + [1])

    def up_list(self):
        """
        Returns a list of the partitions that can be formed from the partition
        p by adding a box.

        EXAMPLES:
            sage: Partition([2,1,1]).up_list()
            [[3, 1, 1], [2, 2, 1], [2, 1, 1, 1]]
            sage: Partition([3,2]).up_list()
            [[4, 2], [3, 3], [3, 2, 1]]
        """
        return [p for p in self.up()]

    def down(self):
        r"""
        Returns a generator for partitions that can be obtained from
        p by removing a box.

        EXAMPLES:
            sage: [p for p in Partition([2,1,1]).down()]
            [[1, 1, 1], [2, 1]]
            sage: [p for p in Partition([3,2]).down()]
            [[2, 2], [3, 1]]
            sage: [p for p in Partition([3,2,1]).down()]
            [[2, 2, 1], [3, 1, 1], [3, 2]]
        """
        p = self
        for i in range(len(p)-1):
            if p[i] > p[i+1]:
                yield Partition(p[:i] + [ p[i]-1 ] + p[i+1:])

        last = p[-1]
        if last == 1:
            yield Partition(p[:-1])
        else:
            yield Partition(p[:-1] + [ p[-1] - 1 ])


    def down_list(self):
        """
        Returns a list of the partitions that can be obtained from
        the partition p by removing a box.

        EXAMPLES:
            sage: Partition([2,1,1]).down_list()
            [[1, 1, 1], [2, 1]]
            sage: Partition([3,2]).down_list()
            [[2, 2], [3, 1]]
            sage: Partition([3,2,1]).down_list()
            [[2, 2, 1], [3, 1, 1], [3, 2]]
        """
        return [p for p in self.down()]

    def dominates(self, p2):
        r"""
        Returns True if partition p1 dominates partitions p2. Otherwise,
        it returns False.

        EXAMPLES:
            sage: p = Partition([3,2])
            sage: p.dominates([3,1])
            True
            sage: p.dominates([2,2])
            True
            sage: p.dominates([2,1,1])
            False
            sage: p.dominates([3,3])
            False
        """
        p1 = self
        sum1 = 0
        sum2 = 0
        if len(p2) > len(p1):
            return False
        for i in range(max(len(p1), len(p2))):
            if i < len(p1):
                sum1 += p1[i]
            if i < len(p2):
                sum2 += p2[i]
            if sum2 > sum1:
                return False
            sum1 = 0
            sum2 = 0
        return True

    def conjugate(self):
        """
        conjugate() returns the ``conjugate'' (also called
        ``associated'' in the literature) partition of the partition pi which is
        obtained by transposing the corresponding Ferrers diagram.

        EXAMPLES:
            sage: Partition([2,2]).conjugate()
            [2, 2]
            sage: Partition([6,3,1]).conjugate()
            [3, 2, 2, 1, 1, 1]
            sage: print Partition([6,3,1]).ferrers_diagram()
            ******
            ***
            *
            sage: print Partition([6,3,1]).conjugate().ferrers_diagram()
            ***
            **
            **
            *
            *
            *

        """
        p = self
        if p == []:
            return Partition([])
        else:
            l = len(p)
            conj =  [l]*p[l-1]
            for i in range(2,l+1):
                conj += [l-i+1]*(p[l-i] - p[l-i+1])

            return Partition(conj)


    def associated(self):
        """
        An alias for partition.conjugate(pi).

        EXAMPLES:
            sage: Partition([4,1,1]).associated()
            [3, 1, 1, 1]
            sage: Partition([4,1,1]).conjugate()
            [3, 1, 1, 1]
            sage: Partition([5,4,2,1,1,1]).associated().associated()
            [5, 4, 2, 1, 1, 1]
        """

        return self.conjugate()

    def arm(self, i, j):
        """
        Returns the arm of cell (i,j) in partition p.  The arm of cell (i,j)
        is the number of boxes that appear to the right of cell (i,j).
        Note that i and j are 0-based indices.

        EXAMPLES:
            sage: p = Partition([2,2,1])
            sage: p.arm(0, 0)
            1
            sage: p.arm(0, 1)
            0
            sage: p.arm(2, 0)
            0
            sage: Partition([3,3]).arm(0, 0)
            2
        """
        p = self
        if i < len(p) and j < p[i]:
            return p[i]-(j+1)
        else:
            #Error: invalid coordinates
            pass

    def arm_lengths(self):
        """
        Returns a tableau of shape p where each box is filled its arm.

        EXAMPLES:
            sage: Partition([2,2,1]).arm_lengths()
            [[1, 0], [1, 0], [0]]
            sage: Partition([3,3]).arm_lengths()
            [[2, 1, 0], [2, 1, 0]]
        """
        p = self
        return [[p[i]-(j+1) for j in range(p[i])] for i in range(len(p))]


    def leg(self, i, j):
        """
        Returns the leg of box (i,j) in partition p. The leg of box (i,j)
        is defined to be the number of boxes below it in partition p.


        EXAMPLES:
            sage: p = Partition([2,2,1])
            sage: p.leg(0, 0)
            2
            sage: p.leg(0,1)
            1
            sage: p.leg(2,0)
            0
            sage: Partition([3,3]).leg(0, 0)
            1
        """

        conj = self.conjugate()
        if j < len(conj) and i < conj[j]:
            return conj[j]-(i+1)
        else:
            #Error: invalid coordinates
            pass

    def leg_lengths(self):
        """
        Returns a tableau of shape p with each box filled in with
        its leg.

        EXAMPLES:
            sage: Partition([2,2,1]).leg_lengths()
            [[2, 1], [1, 0], [0]]
            sage: Partition([3,3]).leg_lengths()
            [[1, 1, 1], [0, 0, 0]]
        """
        p = self
        conj = p.conjugate()
        return [[conj[j]-(i+1) for j in range(p[i])] for i in range(len(p))]

    def hook(self, i, j):
        """
        Returns the hook of box (i,j) in the partition p.  The hook of box
        (i,j) is defined to be one more than the sum of number of boxes
        to the right and the number of boxes below.

        EXAMPLES:
            sage: p = Partition([2,2,1])
            sage: p.hook(0, 0)
            4
            sage: p.hook(0, 1)
            2
            sage: p.hook(2, 0)
            1
            sage: Partition([3,3]).hook(0, 0)
            4
        """
        return self.leg(i,j)+self.arm(i,j)+1


    def hook_lengths(self):
        r"""
        Returns a tableau of shape pi with the boxes filled in with the
        hook lengths

        In each box, put the sum of one plus the number of boxes horizontally to the right
        and vertically below the box (the hook length).


        For example, consider the partition [3,2,1] of 6 with Ferrers Diagram
        * * *
        * *
        *
        When we fill in the boxes with the hook lengths, we obtain
        5 3 1
        3 1
        1

        EXAMPLES:
            sage: Partition([2,2,1]).hook_lengths()
            [[4, 2], [3, 1], [1]]
            sage: Partition([3,3]).hook_lengths()
            [[4, 3, 2], [3, 2, 1]]
            sage: Partition([3,2,1]).hook_lengths()
            [[5, 3, 1], [3, 1], [1]]
            sage: Partition([2,2]).hook_lengths()
            [[3, 2], [2, 1]]
            sage: Partition([5]).hook_lengths()
            [[5, 4, 3, 2, 1]]


        REFERENCES:
            http://mathworld.wolfram.com/HookLengthFormula.html
        """
        p = self
        conj = p.conjugate()
        return [[p[i]-(i+1)+conj[j]-(j+1)+1 for j in range(p[i])] for i in range(len(p))]


    def weighted_size(self):
        """
        Returns sum([i*p[i] for i in range(len(p))]).

        EXAMPLES:
            sage: Partition([2,2,1]).weighted_size()
            9
            sage: Partition([3,3]).weighted_size()
            9
        """
        p = self
        return sum([(i+1)*p[i] for i in range(len(p))])


    def to_exp(self, k=0):
        """
        Return a list of the multiplicities of the parts of a partition.

        EXAMPLES:
            sage: Partition([3,2,2,1]).to_exp()
            [0, 1, 2, 1]
        """
        p = self
        if len(p) > 0:
            k = max(k, p[0])
        a = [ 0 ] * (k+1)
        for i in p:
            a[i] += 1
        return a

    def evaluation(self):
        """
        Returns the evaluation of the partition.

        EXAMPLES:
            sage: Partition([4,3,1,1]).evaluation()
            [0, 2, 0, 1, 1]
        """
        return self.to_exp()

    def centralizer_size(self, t=0, q=0):
        """
        Returns the size of the centralizer of any permuation of cycle type p.

        EXAMPLES:
            sage: Partition([2,2,1]).centralizer_size()
            8
            sage: Partition([2,2,2]).centralizer_size()
            48

        REFERENCES:
            Kerber, A. 'Algebraic Combintorics via Finite Group Action', 1.3 p24
        """
        p = self
        a = p.to_exp()
        size = prod([i**a[i]*factorial(a[i]) for i in range(1, len(a))])
        size *= prod( [ (1-q**p[i])/(1-t**p[i]) for i in range(len(p)) ] )

        return size


    def conjugacy_class_size(self):
        """
        Returns the size of the conjugacy class of the symmetric group
        indexed by the partition p.

        EXAMPLES:
            sage: Partition([2,2,2]).conjugacy_class_size()
            15
            sage: Partition([2,2,1]).conjugacy_class_size()
            15
            sage: Partition([2,1,1]).conjugacy_class_size()
            6


        REFERENCES:
            Kerber, A. 'Algebraic Combinatorics via Finite Group Action' 1.3 p24
        """

        return factorial(sum(self))/self.centralizer_size()


    def corners(self):
        """
        Returns a list of the corners of the partitions.  These are the
        positions where we can remove a box.

        EXAMPLES:
            sage: Partition([3,2,1]).corners()
            [[0, 2], [1, 1], [2, 0]]
            sage: Partition([3,3,1]).corners()
            [[1, 2], [2, 0]]

        """
        p = self
        if p == []:
            return []

        lcors = [[0,p[0]-1]]
        nn = len(p)
        if nn == 1:
            return lcors

        lcors_index = 0
        for i in range(1, nn):
            if p[i] == p[i-1]:
                lcors[lcors_index][0] += 1
            else:
                lcors.append([i,p[i]-1])
                lcors_index += 1

        return lcors


    def outside_corners(self):
        """
        Returns a list of the positions where we can add a box so
        that the shape is still a partition.

        EXAMPLES:
            sage: Partition([2,2,1]).outside_corners()
            [[0, 2], [2, 1], [3, 0]]
            sage: Partition([2,2]).outside_corners()
            [[0, 2], [2, 0]]
            sage: Partition([6,3,3,1,1,1]).outside_corners()
            [[0, 6], [1, 3], [3, 1], [6, 0]]

        """
        p = self
        res = [ [0, p[0]] ]
        for i in range(1, len(p)):
            if p[i-1] != p[i]:
                res.append([i,p[i]])
        res.append([len(p), 0])

        return res

    def r_core(self, length):
        """
        Returns the r-core of the partition p.

        EXAMPLES:
            sage: Partition([6,3,2,2]).r_core(3)
            [2, 1, 1]
        """
        p = self
        #Normalize the length
        remainder = len(p) % length
        part = p[:] + [0]*remainder

        #Add the canonical vector to the partition
        part = [part[i-1] + len(part)-i for i in range(1, len(part)+1)]

        for e in range(length):
            k = e
            for i in reversed(range(1,len(part)+1)):
                if part[i-1] % length == e:
                    part[i-1] = k
                    k += length
        part.sort()
        part.reverse()

        #Remove the canonical vector
        part = [part[i-1]-len(part)+i for i in range(1, len(part)+1)]
        #Select the r-core
        return filter(lambda x: x != 0, part)

    def r_quotient(self, length):
        """
        Returns the r-quotient of the partition p.

        EXAMPLES:
            sage: Partition([7,7,5,3,3,3,1]).r_quotient(3) #[[2], [1], [2, 2, 2]]?
            [[1], [2, 2, 2], [2]]
        """
        p = self
        #Normalize the length
        remainder = len(p) % length
        part = p[:] + [0]*remainder


        #Add the canonical vector to the partition
        part = [part[i-1] + len(part)-i for i in range(1, len(part)+1)]

        result = [None]*length

        for e in range(length):
            k = e
            tmp = []
            for i in reversed(range(len(part))):
                if part[i] % length == e:
                    tmp.append((part[i]-k)/length)
                    k += length

            a = filter(lambda x: x != 0, tmp)
            a.reverse()
            result[e] = a

        return result


    def remove_box(self, i, j = None):
        """
        Returns the partiton obtained by removing a box at the end of row i.

        EXAMPLES:
            sage: Partition([2,2]).remove_box(1)
            [2, 1]
            sage: Partition([2,2,1]).remove_box(2)
            [2, 2]
            sage: #Partition([2,2]).remove_box(0)

            sage: Partition([2,2]).remove_box(1,1)
            [2, 1]
            sage: #Partition([2,2]).remove_box(1,0)
        """

        if i >= len(self):
            raise ValueError, "i must be less than the length of the partition"

        if j == None:
            j = self[i] - 1

        if [i,j] not in self.corners():
            raise ValueError, "[%d,%d] is not a corner of the partition" % (i,j)

        if self[i] == 1:
            return Partition(self[:-1])
        else:
            return Partition(self[:i] + [ self[i:i+1][0] - 1 ] + self[i+1:])



    def k_skew(self, k):
        r"""
        Returns the k-skew partition.

        The k-skew diagram of a k-bounded partition is the skew diagram denoted
        $\lambda/^k$ satisfying the conditions:
        (i) row i of $\lambda/^k$ has length $\lambda_i$
        (ii) no cell in $\lambda/^k$ has hook-length exceeding k
        (iii) every square above the diagram of $\lambda/^k$ has hook
        length exceeding k.

        REFERENCES:
            Lapointe, L. and Morse, J. 'Order Ideals in Weak Subposets of Young's Lattice
                and Associated Unimodality Conjectures'

        EXAMPLES:
            sage: p = Partition([4,3,2,2,1,1])
            sage: p.k_skew(4)
            [[9, 5, 3, 2, 1, 1], [5, 2, 1]]
        """

        if len(self) == 0:
            return sage.combinat.skew_partition.SkewPartition([[],[]])

        if self[0] > k:
            raise ValueError, "the partition must be %d-bounded" % k

        #Find the k-skew diagram of the partition formed
        #by removing the first row
        s = Partition(self[1:]).k_skew(k)

        s_inner = s.inner()
        s_outer = s.outer()
        s_conj_rl = s.conjugate().row_lengths()

        #Find the leftmost column with less than
        # or equal to kdiff boxes
        kdiff = k - self[0]

        if s_outer == []:
            spot = 0
        else:
            spot = s_outer[0]

        for i in range(len(s_conj_rl)):
            if s_conj_rl[i] <= kdiff:
                spot = i
                break

        outer = [ self[0] + spot ] + s_outer[:]
        if spot > 0:
            inner = [ spot ] + s_inner[:]
        else:
            inner = s_inner[:]

        return sage.combinat.skew_partition.SkewPartition([outer, inner])


    def k_conjugate(self, k):
        """
        Returns the k-conjugate of the partition.

        The k-conjugate is the partition that is given by the columns
        of the k-skew diagram of the partition.

        EXAMPLES:
            sage: p = Partition([4,3,2,2,1,1])
            sage: p.k_conjugate(4)
            [3, 2, 2, 1, 1, 1, 1, 1, 1]
        """
        return Partition(self.k_skew(k).conjugate().row_lengths())

    def parent(self):
        """
        Returns the combinatorial class of partitions of
        sum(self).

        EXAMPLES:
            sage: Partition([3,2,1]).parent()
            Partitions of the integer 6
        """
        return Partitions(sum(self[:]))

    def arms_legs_coeff(self, i , j):
        """
        EXAMPLES:
            sage: Partition([3,2,1]).arms_legs_coeff(1,1)
            (-t + 1)/(-q + 1)
            sage: Partition([3,2,1]).arms_legs_coeff(0,0)
            (-q^2*t^3 + 1)/(-q^3*t^2 + 1)
        """
        QQqt = PolynomialRing(QQ, ['q', 't'])
        (q,t) = QQqt.gens()
        if i < len(self) and j < self[i]:
            res =  (1-q**self.arm(i,j) * t**(self.leg(i,j)+1))
            res /= (1-q**(self.arm(i,j)+1) * t**self.leg(i,j))
            return res
        else:
            return ZZ(1)

    def macdonald_coeff(self):
        res = 1
        for i in range(len(self)):
            for j in range(self[i]):
                res *= self.arms_legs_coeff(i,j)
        return res

    def jacobi_trudy(self):
        """

        EXAMPLES:
            sage: part = Partition([3,2,1])
            sage: jt = part.jacobi_trudy(); jt
            [h[3] h[1]    0]
            [h[4] h[2]  h[]]
            [h[5] h[3] h[1]]
            sage: s = SFASchur(QQ)
            sage: h = SFAHomogeneous(QQ)
            sage: h( s(part) )
            h[3, 2, 1] - h[3, 3] - h[4, 1, 1] + h[5, 1]
            sage: jt.det()
            h[3, 2, 1] - h[3, 3] - h[4, 1, 1] + h[5, 1]
        """
        return sage.combinat.skew_partition.SkewPartition([ self, [] ]).jacobi_trudy()

##################################################


##################
# Set Partitions #
##################

def partitions_set(S,k=None):
    r"""

    WARNING: Wraps GAP -- hence S must be a list of objects that have
    string representations that can be interpreted by the GAP
    intepreter.  If mset consists of at all complicated SAGE objects,
    this function does *not* do what you expect.  A proper function
    should be written! (TODO!)

    Wraps GAP's PartitionsSet.


    """
    if k==None:
        ans=gap.eval("PartitionsSet(%s)"%S)
    else:
        ans=gap.eval("PartitionsSet(%s,%s)"%(S,k))
    return eval(ans)

def number_of_partitions_set(S,k):
    r"""
    Returns the size of \code{partitions_set(S,k)}.  Wraps GAP's
    NrPartitionsSet.

    The Stirling number of the second kind is the number of partitions
    of a set of size n into k blocks.

    EXAMPLES:
        sage: mset = [1,2,3,4]
        sage: partition.number_of_partitions_set(mset,2)
        7
        sage: stirling_number2(4,2)
        7

    REFERENCES:
        http://en.wikipedia.org/wiki/Partition_of_a_set

    """
    if k==None:
        ans=gap.eval("NrPartitionsSet(%s)"%S)
    else:
        ans=gap.eval("NrPartitionsSet(%s,%s)"%(S,ZZ(k)))
    return ZZ(ans)

def number_of_partitions_list(n,k=None):
    r"""
    Returns the size of partitions_list(n,k).

    Wraps GAP's NrPartitions.

    It is possible to associate with every partition of the integer n
    a conjugacy class of permutations in the symmetric group on n
    points and vice versa.  Therefore p(n) = NrPartitions(n) is the
    number of conjugacy classes of the symmetric group on n points.

    \code{number_of_partitions(n)} is also available in PARI, however
    the speed seems the same until $n$ is in the thousands (in which
    case PARI is faster).

    EXAMPLES:
        sage: partition.number_of_partitions_list(10,2)
        5
        sage: partition.number_of_partitions_list(10)
        42

    A generating function for p(n) is given by the reciprocal of Euler's function:
    \[
    \sum_{n=0}^\infty p(n)x^n = \prod_{k=1}^\infty \left(\frac {1}{1-x^k} \right).
    \]
    SAGE verifies that the first several coefficients do instead agree:

        sage: q = PowerSeriesRing(QQ, 'q', default_prec=9).gen()
        sage: prod([(1-q^k)^(-1) for k in range(1,9)])  ## partial product of
        1 + q + 2*q^2 + 3*q^3 + 5*q^4 + 7*q^5 + 11*q^6 + 15*q^7 + 22*q^8 + O(q^9)
        sage: [partition.number_of_partitions_list(k) for k in range(2,10)]
        [2, 3, 5, 7, 11, 15, 22, 30]

    REFERENCES:
        http://en.wikipedia.org/wiki/Partition_%28number_theory%29

    """
    if k==None:
        ans=gap.eval("NrPartitions(%s)"%(ZZ(n)))
    else:
        ans=gap.eval("NrPartitions(%s,%s)"%(ZZ(n),ZZ(k)))
    return ZZ(ans)


######################
# Ordered Partitions #
######################
def OrderedPartitions(n, k=None):
    return OrderedPartitions_nk(n,k)

class OrderedPartitions_nk(CombinatorialClass):
    def __init__(self, n, k=None):
        self.n = n
        self.k = k

    def __repr__(self):
        string = "Ordered partitions of %s "%self.n
        if self.k != None:
            string += "of length %s"%self.k
        return string

    def list(self):
        n = self.n
        k = self.k
        if self.k==None:
            ans=gap.eval("OrderedPartitions(%s)"%(ZZ(n)))
        else:
            ans=gap.eval("OrderedPartitions(%s,%s)"%(ZZ(n),ZZ(k)))
        result = eval(ans.replace('\n',''))
        result.reverse()
        return result

    def count(self):
        n = self.n
        k = self.k
        if k==None:
            ans=gap.eval("NrOrderedPartitions(%s)"%(n))
        else:
            ans=gap.eval("NrOrderedPartitions(%s,%s)"%(n,k))
        return ZZ(ans)

def ordered_partitions(n,k=None):
    r"""
    An {\it ordered partition of $n$} is an ordered sum
    $$
       n = p_1+p_2 + \cdots + p_k
    $$
    of positive integers and is represented by the list $p = [p_1,p_2,\cdots ,p_k]$.
    If $k$ is omitted then all ordered partitions are returned.

    \code{ordered_partitions(n,k)} returns the set of all (ordered)
    partitions of the positive integer n into sums with k summands.

    Do not call \code{ordered_partitions} with an n much larger than
    15, since the list will simply become too large.

    Wraps GAP's OrderedPartitions.

    The number of ordered partitions $T_n$ of $\{ 1, 2, ..., n \}$ has the
    generating function is
    \[
    \sum_n {T_n \over n!} x^n = {1 \over 2-e^x}.
    \]

    EXAMPLES:
        sage: partition.ordered_partitions(10,2)
        [[1, 9], [2, 8], [3, 7], [4, 6], [5, 5], [6, 4], [7, 3], [8, 2], [9, 1]]

        sage: partition.ordered_partitions(4)
        [[1, 1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 3], [2, 1, 1], [2, 2], [3, 1], [4]]

    REFERENCES:
        http://en.wikipedia.org/wiki/Ordered_partition_of_a_set

    """
    if k==None:
        ans=gap.eval("OrderedPartitions(%s)"%(ZZ(n)))
    else:
        ans=gap.eval("OrderedPartitions(%s,%s)"%(ZZ(n),ZZ(k)))
    result = eval(ans.replace('\n',''))
    return [Partition(p) for p in result]

def number_of_ordered_partitions(n,k=None):
    """
    Returns the size of ordered_partitions(n,k).
    Wraps GAP's NrOrderedPartitions.

    It is possible to associate with every partition of the integer n a conjugacy
    class of permutations in the symmetric group on n points and vice versa.
    Therefore p(n) = NrPartitions(n) is the number of conjugacy classes of the
    symmetric group on n points.


    EXAMPLES:
        sage: partition.number_of_ordered_partitions(10,2)
        9
        sage: partition.number_of_ordered_partitions(15)
        16384
    """
    if k==None:
        ans=gap.eval("NrOrderedPartitions(%s)"%(n))
    else:
        ans=gap.eval("NrOrderedPartitions(%s,%s)"%(n,k))
    return ZZ(ans)

##########################
# Partitions Greatest LE #
##########################
def PartitionsGreatestLE(n,k):
    return PartitionsGreatestLE_nk(n,k)

class PartitionsGreatestLE_nk(CombinatorialClass):
    """
    The combinatorial class of all (unordered) ``restricted'' partitions of the
    integer n having parts less than or equal to the integer k.
    """
    object_class = Partition_class
    def __init__(self, n, k):
        self.n = n
        self.k = k

    def __repr__(self):
        return "Partitions of %s having parts less than or equal to %s"%(self.n, self.k)

    def list(self):
        """
        Returns a list of all (unordered) ``restricted'' partitions of the
        integer n having parts less than or equal to the integer k.

        Wraps GAP's PartitionsGreatestLE.

        EXAMPLES:
            sage: PartitionsGreatestLE(10,2).list()
            [[2, 2, 2, 2, 2],
             [2, 2, 2, 2, 1, 1],
             [2, 2, 2, 1, 1, 1, 1],
             [2, 2, 1, 1, 1, 1, 1, 1],
             [2, 1, 1, 1, 1, 1, 1, 1, 1],
             [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
        """
        result = eval(gap.eval("PartitionsGreatestLE(%s,%s)"%(ZZ(self.n),ZZ(self.k))))
        result.reverse()
        return [Partition(p) for p in result]

##########################
# Partitions Greatest EQ #
##########################
def PartitionsGreatestEQ(n,k):
    return PartitionsGreatestEQ_nk(n,k)

class PartitionsGreatestEQ_nk(CombinatorialClass):
    """
    The combinatorial class of all (unordered) ``restricted'' partitions of the
    integer n having at least one part equal to the integer k.
    """
    object_class = Partition_class
    def __init__(self, n, k):
        self.n = n
        self.k = k

    def __repr__(self):
        return "Partitions of %s having at least one part equal to %s"%(self.n, self.k)

    def list(self):
        """
        Returns a list of all (unordered) ``restricted'' partitions of the
        integer n having at least one part equal to the integer k.

        Wraps GAP's PartitionsGreatestEQ.

        EXAMPLES:
           sage: PartitionsGreatestEQ(10,2).list()
           [[2, 2, 2, 2, 2],
            [2, 2, 2, 2, 1, 1],
            [2, 2, 2, 1, 1, 1, 1],
            [2, 2, 1, 1, 1, 1, 1, 1],
            [2, 1, 1, 1, 1, 1, 1, 1, 1]]
        """
        result = eval(gap.eval("PartitionsGreatestEQ(%s,%s)"%(self.n,self.k)))
        result.reverse()
        return [Partition(p) for p in result]


#########################
# Restricted Partitions #
#########################
def RestrictedPartitions(n, S, k=None):
    return RestrictedPartitions_nsk(n, S, k)

class RestrictedPartitions(CombinatorialClass):
    r"""
    A {\it restricted partition} is, like an ordinary partition, an
    unordered sum $n = p_1+p_2+\ldots+p_k$ of positive integers and is
    represented by the list $p = [p_1,p_2,\ldots,p_k]$, in nonincreasing
    order. The difference is that here the $p_i$ must be elements
    from the set $S$, while for ordinary partitions they may be
    elements from $[1..n]$.

    """
    object_class = Partition_class
    def __init__(self, n, S, k=None):
        self.n = n
        self.S = S
        self.S.sort()
        self.k = k

    def __repr__(self):
        string = "Partitions of %s restricted to the values %s "%(self.n, self.S)
        if self.k != None:
            string += "of length %s" % self.k
        return string

    def list(self):
        r"""
        Returns the set of all restricted partitions of the positive integer
        n into sums with k summands with the summands of the partition coming
        from the set $S$. If k is not given all restricted partitions for all
        k are returned.

        Wraps GAP's RestrictedPartitions.

        EXAMPLES:
            sage: RestrictedPartitions(8,[1,3,5,7]).list()
            [[7, 1],
            [5, 3],
            [5, 1, 1, 1],
            [3, 3, 1, 1],
            [3, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1]]
            sage: RestrictedPartitions(8,[1,3,5,7],2).list()
            [[7, 1], [5, 3]]
        """
        n = self.n
        k = self.k
        S = self.S
        if k==None:
            ans=gap.eval("RestrictedPartitions(%s,%s)"%(n,S))
        else:
            ans=gap.eval("RestrictedPartitions(%s,%s,%s)"%(n,S,k))
        result = eval(ans)
        result.reverse()
        return [Partition(p) for p in result]

    def count(self):
        """
        Returns the size of RestrictedPartitions(n,S,k).
        Wraps GAP's NrRestrictedPartitions.

        EXAMPLES:
            sage: RestrictedPartitions(8,[1,3,5,7]).count()
            6
            sage: RestrictedPartitions(8,[1,3,5,7],2).count()
            2
        """
        n = self.n
        k = self.k
        S = self.S
        if k==None:
            ans=gap.eval("NrRestrictedPartitions(%s,%s)"%(ZZ(n),S))
        else:
            ans=gap.eval("NrRestrictedPartitions(%s,%s,%s)"%(ZZ(n),S,ZZ(k)))
        return ZZ(ans)

####################
# Partition Tuples #
####################
def PartitionTuples(n,k):
    return PartitionTuples_nk(n,k)

class PartitionTuples_nk(CombinatorialClass):
    """
    k-tuples of partitions describe the classes and the characters of
    wreath products of groups with k conjugacy classes with the symmetric
    group $S_n$.

    """
    object_class = Partition_class
    def __init__(self, n, k):
        self.n = n
        self.k = k

    def __repr__(self):
        return "%s-tuples of partitions of %s"%(self.k, self.n)

    def list(self):
        """
        Returns the list of all k-tuples of partitions
        which together form a partition of n.

        Wraps GAP's PartitionTuples.

        EXAMPLES:
            sage: PartitionTuples(3,2).list()
            [[[1, 1, 1], []],
            [[1, 1], [1]],
            [[1], [1, 1]],
            [[], [1, 1, 1]],
            [[2, 1], []],
            [[1], [2]],
            [[2], [1]],
            [[], [2, 1]],
            [[3], []],
            [[], [3]]]
        """
        ans=gap.eval("PartitionTuples(%s,%s)"%(ZZ(self.n),ZZ(self.k)))
        return map(lambda x: map(lambda y: Partition(y), x), eval(ans))

    def count(self):
        r"""
        Returns the number of k-tuples of partitions which together
        form a partition of n.

        Wraps GAP's NrPartitionTuples.

        EXAMPLES:
            sage: PartitionTuples(3,2).count()
            10
            sage: PartitionTuples(8,2).count()
            185

            Now we compare that with the result of the following GAP
            computation:
            \begin{verbatim}
            gap> S8:=Group((1,2,3,4,5,6,7,8),(1,2));
            Group([ (1,2,3,4,5,6,7,8), (1,2) ])
            gap> C2:=Group((1,2));
            Group([ (1,2) ])
            gap> W:=WreathProduct(C2,S8);
            <permutation group of size 10321920 with 10 generators>
            gap> Size(W);
            10321920     ## = 2^8*Factorial(8), which is good:-)
            gap> Size(ConjugacyClasses(W));
            185
            \end{verbatim}
        """

        ans=gap.eval("NrPartitionTuples(%s,%s)"%(ZZ(self.n),ZZ(self.k)))
        return ZZ(ans)

##############
# Partitions #
##############

def Partitions(n=None, **kwargs):
    if n == None:
        return Partitions_all()
    else:
        if len(kwargs) == 0:
            return Partitions_n(n)
        else:
            return Partitions_constraints(n, **kwargs)

class Partitions_all(CombinatorialClass):
    object_class = Partition_class
    def __contains__(self, x):
        if isinstance(x, Partition_class):
            return True
        elif isinstance(x, __builtin__.list):
            for i in range(len(x)):
                if not isinstance(x[i], (int, Integer)):
                    return False
                if x[i] < 0:
                    return False
                if i == 0:
                    prev = x[i]
                    continue
                if x[i] > prev:
                    return False
                prev = x[i]
            return True
        else:
            return False

    def list(self):
        raise NotImplementedError

class Partitions_constraints(CombinatorialClass):
    def __init__(self, n, **kwargs):
        self.n = n
        self.constraints = kwargs

    def __contains__(self, x):
        return x in Partitions_all() and sum(x)==self.n and misc.check_integer_list_constraints(x, singleton=True, **self.constraints)

    def size(self):
        return self.n

    def __repr__(self):
        return "Partitions of the integer %s satisfying constraints"%self.n

    def iterator(self):
        """
        An iterator a list of the partitions of n.

        EXAMPLES:
            sage: [x for x in Partitions(4)]
            [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
            sage: [x for x in Partitions(4, length=2)]
            [[3, 1], [2, 2]]
            sage: [x for x in Partitions(4, min_length=2)]
            [[3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
            sage: [x for x in Partitions(4, max_length=2)]
            [[4], [3, 1], [2, 2]]
            sage: [x for x in Partitions(4, min_length=2, max_length=2)]
            [[3, 1], [2, 2]]
            sage: [x for x in Partitions(4, max_part=2)]
            [[2, 2], [2, 1, 1], [1, 1, 1, 1]]
            sage: [x for x in Partitions(4, min_part=2)]
            [[4], [2, 2]]
            sage: [x for x in Partitions(4, length=3, min_part=0)]
            [[4, 0, 0], [3, 1, 0], [2, 2, 0], [2, 1, 1]]
            sage: [x for x in Partitions(4, min_length=3, min_part=0)]
            [[4, 0, 0], [3, 1, 0], [2, 2, 0], [2, 1, 1], [1, 1, 1, 1]]
            sage: [x for x in Partitions(4, outer=[3,1,1])]
            [[3, 1], [2, 1, 1]]
            sage: [x for x in Partitions(4, outer=['inf', 1, 1])]
            [[4], [3, 1], [2, 1, 1]]
            sage: [x for x in Partitions(4, inner=[1,1,1])]
            [[2, 1, 1], [1, 1, 1, 1]]
            sage: [x for x in Partitions(4, max_slope=-1)]
            [[4], [3, 1]]
            sage: [x for x in Partitions(4, min_slope=-1)]
            [[4], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
            sage: [x for x in Partitions(11, max_slope=-1, min_slope=-3, min_length=2, max_length=4)]
            [[7, 4], [6, 5], [6, 4, 1], [6, 3, 2], [5, 4, 2], [5, 3, 2, 1]]
            sage: [x for x in Partitions(11, max_slope=-1, min_slope=-3, min_length=2, max_length=4, outer=[6,5,2])]
            [[6, 5], [6, 4, 1], [6, 3, 2], [5, 4, 2]]
        """
        n = self.n

        #put the constraints in the namespace
        length = self.constraints.get('length',None)
        min_length = self.constraints.get('min_length',None)
        max_length = self.constraints.get('max_length',None)
        max_part = self.constraints.get('max_part',None)
        min_part = self.constraints.get('min_part',None)
        max_slope = self.constraints.get('max_slope',None)
        min_slope = self.constraints.get('min_slope',None)
        inner = self.constraints.get('inner',None)
        outer = self.constraints.get('outer',None)

        #Preproccess the constraints
        if max_slope == "inf":
            max_slope = n
        if min_slope == "-inf":
            min_slope = -n
        if length != None:
            min_length = length
            max_length = length
        if outer != None:
            for i in range(len(outer)):
                if outer[i] == "inf":
                    outer[i] = n

        l = []
        old_p = Partition_class([self.n])

        #Go through all of the partitions
        while old_p != False:
            meets_constraints = True

            #Handle the case if the user specifies
            #min_part = 0
            if min_part == 0:
                #if length != None and len(p) < length:
                #    p += [0]*(len(p)-length)
                if min_length != None and len(old_p) < min_length:
                    p = old_p + [0]*(min_length-len(old_p))
                else:
                    p = old_p
            else:
                p = old_p

            #if length != None:
            #    if len(p) != length:
            #        meets_constraints = False

            if min_length != None:
                if len(p) < min_length:
                    meets_constraints = False

            if max_length != None:
                if len(p) > max_length:
                    meets_constraints = False

            if max_part != None:
                if max(p) > max_part:
                    meets_constraints = False

            if min_part != None:
                if min(p) < min_part:
                    meets_constraints = False

            if outer != None:
                if len(p) > len(outer):
                    meets_constraints = False
                else:
                    for i in range(len(p)):
                        if p[i] > outer[i]:
                            meets_constraints = False

            if inner != None:
                if len(p) < len(inner):
                    meets_constraints = False
                else:
                    for i in range(len(inner)):
                        if p[i] < inner[i]:
                            meets_constraints = False


            if max_slope != None:
                if max([p[i+1]-p[i] for i in range(0, len(p)-1)]+[-n] ) > max_slope:
                    meets_constraints = False

            if min_slope != None:
                if min([p[i+1]-p[i] for i in range(0, len(p)-1)]+[n] ) < min_slope:
                    meets_constraints = False


            if meets_constraints:
                yield p

            old_p = old_p.next()

class Partitions_n(CombinatorialClass):
    object_class = Partition_class
    def __init__(self, n):
        self.n = n

    def __contains__(self, x):
        return x in Partitions_all() and sum(x)==self.n

    def size(self):
        return self.n

    def __repr__(self):
        return "Partitions of the integer %s"%self.n

    def count(self, algorithm='default'):
        r"""
            algorithm -- (default: 'default')
                'bober' -- use Jonathon Bober's implementation (*very* fast,
                          but new and not well tested yet).
                'gap' -- use GAP (VERY *slow*)
                'pari' -- use PARI.  Speed seems the same as GAP until $n$ is
                          in the thousands, in which case PARI is faster. *But*
                          PARI has a bug, e.g., on 64-bit Linux PARI-2.3.2
                          outputs numbpart(147007)%1000 as 536, but it
                          should be 533!.  So do not use this option.
                'default' -- 'bober' when k is not specified; otherwise
                          use 'gap'.

        IMPLEMENTATION: Wraps GAP's NrPartitions or PARI's numbpart function.

        Use the function \code{partitions(n)} to return a generator over
        all partitions of $n$.

        It is possible to associate with every partition of the integer n
        a conjugacy class of permutations in the symmetric group on n
        points and vice versa.  Therefore p(n) = NrPartitions(n) is the
        number of conjugacy classes of the symmetric group on n points.

        EXAMPLES:
            sage: v = Partitions(5).list(); v
            [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
            sage: len(v)
            7
            sage: Partitions(5).count(algorithm='gap')
            7
            sage: Partitions(5).count(algorithm='pari')
            7
            sage: Partitions(5).count(algorithm='bober')
            7

        The input must be a nonnegative integer or a ValueError is raised.
            sage: Partitions(10).count()
            42
            sage: Partitions(3).count()
            3
            sage: Partitions(10).count()
            42
            sage: Partitions(3).count(algorithm='pari')
            3
            sage: Partitions(10).count(algorithm='pari')
            42
            sage: Partitions(40).count()
            37338
            sage: Partitions(100).count()
            190569292

        A generating function for p(n) is given by the reciprocal of
        Euler's function:

        \[
        \sum_{n=0}^\infty p(n)x^n = \prod_{k=1}^\infty \left(\frac {1}{1-x^k} \right).
        \]

        We use SAGE to verify that the first several coefficients do
        instead agree:

            sage: q = PowerSeriesRing(QQ, 'q', default_prec=9).gen()
            sage: prod([(1-q^k)^(-1) for k in range(1,9)])  ## partial product of
            1 + q + 2*q^2 + 3*q^3 + 5*q^4 + 7*q^5 + 11*q^6 + 15*q^7 + 22*q^8 + O(q^9)
            sage: [Partitions(k) .count()for k in range(2,10)]
            [2, 3, 5, 7, 11, 15, 22, 30]

        REFERENCES:
            http://en.wikipedia.org/wiki/Partition_%28number_theory%29

        """
        n = ZZ(self.n)
        if n < 0:
            raise ValueError, "n (=%s) must be a nonnegative integer"%n
        elif n == 0:
            return ZZ(1)
        if algorithm == 'gap':
            ans=gap.eval("NrPartitions(%s)"%(n))
            return ZZ(ans)
        if algorithm == 'default':
            return partitions_ext.number_of_partitions(n)
        elif algorithm == 'bober':
            return partitions_ext.number_of_partitions(n)
        elif algorithm == 'pari':
            return ZZ(pari(n).numbpart())
        raise ValueError, "unknown algorithm '%s'"%algorithm

    def first(self):
        """
        Returns the lexicographically first partition of a
        positive integer n. This is the partition [n].

        EXAMPLES:
            sage: Partitions(4).first()
            [4]
        """
        return Partition([self.n])

    def last(self):
        """
        Returns the lexicographically last partition of the
        positive integer n.  This is the all-ones partition.

        EXAMPLES:
            sage: Partitions(4).last()
            [1, 1, 1, 1]
        """

        return Partition_class([1]*self.n)


    def iterator(self):
        """
        An iterator a list of the partitions of n.

        EXAMPLES:
            sage: [x for x in Partitions(4)]
            [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
        """
        # base case of the recursion: zero is the sum of the empty tuple
        if self.n == 0:
            yield Partition_class([])
            return

        # modify the partitions of n-1 to form the partitions of n
        for p in Partitions_n(self.n-1):
            if p and (len(p) < 2 or p[-2] > p[-1]):
                yield Partition_class(list(p[:-1]) + [p[-1] + 1])
            yield  Partition_class(p + [1])



def PartitionsInBox(h, w):
    return PartitionsInBox_hw(h, w)

class PartitionsInBox_hw(CombinatorialClass):
    object_class = Partition_class
    def __init__(self, h, w):
        self.h = h
        self.w = w

    def __repr__(self):
        return "Integer Partitions which fit in a %s x %s box" % (self.h, self.w)

    def list(self):
        """
        Returns a list of all the partitions inside a box of height h
        and width w.

        EXAMPLES:
            sage: PartitionsInBox(2,2).list()
            [[], [1], [1, 1], [2], [2, 1], [2, 2]]
        """
        h = self.h
        w = self.w
        if h == 0:
            return [[]]
        else:
            l = [[i] for i in range(0, w+1)]
            add = lambda x: [ x+[i] for i in range(0, x[-1]+1)]
            for i in range(h-1):
                new_list = []
                for element in l:
                    new_list += add(element)
                l = new_list

            return [Partition(filter(lambda x: x!=0, p)) for p in l]
