r"""
Partitions

A partition `p` of a nonnegative integer `n` is a
non-increasing list of positive integers (the *parts* of the
partition) with total sum `n`.

A partition can be depicted by a diagram made of rows of boxes,
where the number of boxes in the `i^{th}` row starting from
the top is the `i^{th}` part of the partition.

The coordinate system related to a partition applies from the top
to the bottom and from left to right. So, the corners of the
partition are [[0,4], [1,2], [2,0]].

EXAMPLES: There are 5 partitions of the integer 4.

::

    sage: Partitions(4).count()
    5
    sage: Partitions(4).list()
    [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

We can use the method .first() to get the 'first' partition of a
number.

::

    sage: Partitions(4).first()
    [4]

Using the method .next(), we can calculute the 'next' partition.
When we are at the last partition, None will be returned.

::

    sage: Partitions(4).next([4])
    [3, 1]
    sage: Partitions(4).next([1,1,1,1]) is None
    True

We can use the .iterator() method to get an object which iterates
over the partitions one by one to save memory. Note that when we do
something like 'for part in Partitions(4)' an iterator is used in
the background.

::

    sage: g = Partitions(4).iterator()
    sage: g.next()
    [4]
    sage: g.next()
    [3, 1]
    sage: g.next()
    [2, 2]
    sage: for p in Partitions(4): print p
    [4]
    [3, 1]
    [2, 2]
    [2, 1, 1]
    [1, 1, 1, 1]

We can add constraints to to the type of partitions we want. For
example, to get all of the partitions of 4 of length 2, we'd do the
following::

    sage: Partitions(4, length=2).list()
    [[3, 1], [2, 2]]

Here is the list of partitions of length at least 2 and the list of
ones with length at most 2::

    sage: Partitions(4, min_length=2).list()
    [[3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    sage: Partitions(4, max_length=2).list()
    [[4], [3, 1], [2, 2]]

The options min_part and max_part can be used to set constraints
on the sizes of all parts. Using max_part, we can select
partitions having only 'small' entries. The following is the list
of the partitions of 4 with parts at most 2::

    sage: Partitions(4, max_part=2).list()
    [[2, 2], [2, 1, 1], [1, 1, 1, 1]]

The min_part options is complementary to max_part and selects
partitions having only 'large' parts. Here is the list of all
partitions of 4 with each part at least 2.

::

    sage: Partitions(4, min_part=2).list()
    [[4], [2, 2]]

The options inner and outer can be used to set part-by-part
constraints. This is the list of partitions of 4 with [3, 1, 1] as
an outer bound::

    sage: Partitions(4, outer=[3,1,1]).list()
    [[3, 1], [2, 1, 1]]

Outer sets max_length to the length of its argument. Moreover, the
parts of outer may be infinite to clear constraints on specific
parts. Here is the list of the partitions of 4 of length at most 3
such that the second and third part are 1 when they exist::

    sage: Partitions(4, outer=[oo,1,1]).list()
    [[4], [3, 1], [2, 1, 1]]

Finally, here are the partitions of 4 with [1,1,1] as an inner
bound. Note that inner sets min_length to the length of its
argument.

::

    sage: Partitions(4, inner=[1,1,1]).list()
    [[2, 1, 1], [1, 1, 1, 1]]

The options min_slope and max_slope can be used to set
constraints on the slope, that is on the difference p[i+1]-p[i] of
two consecutive parts. Here is the list of the strictly decreasing
partitions of 4::

    sage: Partitions(4, max_slope=-1).list()
    [[4], [3, 1]]

The constraints can be combined together in all reasonable ways.
Here are all the partitions of 11 of length between 2 and 4 such
that the difference between two consecutive parts is between -3 and
-1::

    sage: Partitions(11,min_slope=-3,max_slope=-1,min_length=2,max_length=4).list()
    [[7, 4], [6, 5], [6, 4, 1], [6, 3, 2], [5, 4, 2], [5, 3, 2, 1]]

Partition objects can also be created individually with the
Partition function.

::

    sage: Partition([2,1])
    [2, 1]

Once we have a partition object, then there are a variety of
methods that we can use. For example, we can get the conjugate of a
partition. Geometrically, the conjugate of a partition is the
reflection of that partition through its main diagonal. Of course,
this operation is an involution.

::

    sage: Partition([4,1]).conjugate()
    [2, 1, 1, 1]
    sage: Partition([4,1]).conjugate().conjugate()
    [4, 1]

We can go back and forth between the exponential notations of a
partition. The exponential notation can be padded with extra
zeros.

::

    sage: Partition([6,4,4,2,1]).to_exp()
    [1, 1, 0, 2, 0, 1]
    sage: Partition(exp=[1,1,0,2,0,1])
    [6, 4, 4, 2, 1]
    sage: Partition([6,4,4,2,1]).to_exp(5)
    [1, 1, 0, 2, 0, 1]
    sage: Partition([6,4,4,2,1]).to_exp(7)
    [1, 1, 0, 2, 0, 1, 0]
    sage: Partition([6,4,4,2,1]).to_exp(10)
    [1, 1, 0, 2, 0, 1, 0, 0, 0, 0]

We can get the coordinates of the corners of a partition.

::

    sage: Partition([4,3,1]).corners()
    [[0, 3], [1, 2], [2, 0]]

We can compute the r-core and r-quotient of a partition and build
the partition back up from them.

::

    sage: Partition([6,3,2,2]).r_core(3)
    [2, 1, 1]
    sage: Partition([7,7,5,3,3,3,1]).r_quotient(3)
    [[2], [1], [2, 2, 2]]
    sage: p = Partition([11,5,5,3,2,2,2])
    sage: p.r_core(3)
    []
    sage: p.r_quotient(3)
    [[2, 1], [4], [1, 1, 1]]
    sage: Partition(core_and_quotient=([],[[2, 1], [4], [1, 1, 1]]))
    [11, 5, 5, 3, 2, 2, 2]
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

from sage.interfaces.all import gap, gp
from sage.rings.all import QQ, ZZ, infinity, factorial, gcd
from sage.misc.all import prod, sage_eval
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.calculus.all import ceil
import sage.combinat.misc as misc
import sage.combinat.skew_partition
from sage.rings.integer import Integer
import __builtin__
from combinat import CombinatorialClass, CombinatorialObject, cyclic_permutations_iterator
import partitions as partitions_ext
from sage.libs.all import pari
import tableau
import permutation
import sf.sfa
import composition
from integer_vector import IntegerVectors
from cartesian_product import CartesianProduct

def Partition(l=None, exp=None, core_and_quotient=None):
    """
    Returns a partition object.

    Note that Sage uses the English convention for partitions and
    tableaux.

    EXAMPLES::

        sage: Partition(exp=[2,1,1])
        [3, 2, 1, 1]
        sage: Partition(core_and_quotient=([2,1], [[2,1],[3],[1,1,1]]))
        [11, 5, 5, 3, 2, 2, 2]
        sage: Partition([3,2,1])
        [3, 2, 1]
        sage: [2,1] in Partitions()
        True
        sage: [2,1,0] in Partitions()
        False
        sage: Partition([2,1,0])
        [2, 1]
        sage: Partition([1,2,3])
        Traceback (most recent call last):
        ...
        ValueError: [1, 2, 3] is not a valid partition
    """
    number_of_arguments = 0
    for arg in ['l', 'exp', 'core_and_quotient']:
        if locals()[arg] is not None:
            number_of_arguments += 1

    if number_of_arguments != 1:
        raise ValueError, "you must specify exactly one argument"

    if l is not None:
        l = [i for i in l if i != 0]
        if l in Partitions_all():
            return Partition_class(l)
        else:
            raise ValueError, "%s is not a valid partition"%l
    elif exp is not None:
        return from_exp(exp)
    else:
        return from_core_and_quotient(*core_and_quotient)


def from_exp(a):
    """
    Returns a partition from its list of multiplicities.

    EXAMPLES::

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

    EXAMPLES::

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


        -  ``pi`` - a partition, given as a list of integers.


        EXAMPLES::

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

    def pp(self):
        """
        EXAMPLES::

            sage: Partition([5,5,2,1]).pp()
            *****
            *****
            **
            *
        """
        print self.ferrers_diagram()

    def __div__(self, p):
        """
        Returns the skew partition self/p.

        EXAMPLES::

            sage: p = Partition([3,2,1])
            sage: p/[1,1]
            [[3, 2, 1], [1, 1]]
            sage: p/[3,2,1]
            [[3, 2, 1], [3, 2, 1]]
            sage: p/Partition([1,1])
            [[3, 2, 1], [1, 1]]
        """
        if not self.dominates(p):
            raise ValueError, "the partition must dominate p"

        return sage.combinat.skew_partition.SkewPartition([self[:], p])

    def power(self,k):
        """
        partition_power( pi, k ) returns the partition corresponding to
        the `k`-th power of a permutation with cycle structure pi
        (thus describes the powermap of symmetric groups).

        Wraps GAP's PowerPartition.

        EXAMPLES::

            sage: p = Partition([5,3])
            sage: p.power(1)
            [5, 3]
            sage: p.power(2)
            [5, 3]
            sage: p.power(3)
            [5, 1, 1, 1]
            sage: p.power(4)
            [5, 3]
            sage: Partition([3,2,1]).power(3)
            [2, 1, 1, 1, 1]

        Now let us compare this to the power map on `S_8`::

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
        res = []
        for i in self:
            g = gcd(i, k)
            res.extend( [ZZ(i//g)]*int(g) )
        res.sort(reverse=True)
        return Partition_class( res )

    def next(self):
        """
        Returns the partition that lexicographically follows the partition
        p. If p is the last partition, then it returns False.

        EXAMPLES::

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

        EXAMPLES::

            sage: Partition([2,2]).size()
            4
            sage: Partition([3,2,1]).size()
            6
        """
        return sum(self)

    def sign(self):
        r"""
        partition_sign( pi ) returns the sign of a permutation with cycle
        structure given by the partition pi.

        This function corresponds to a homomorphism from the symmetric
        group `S_n` into the cyclic group of order 2, whose kernel
        is exactly the alternating group `A_n`. Partitions of sign
        `1` are called even partitions while partitions of sign
        `-1` are called odd.

        Wraps GAP's SignPartition.

        EXAMPLES::

            sage: Partition([5,3]).sign()
            1
            sage: Partition([5,2]).sign()
            -1

        Zolotarev's lemma states that the Legendre symbol
        `\left(\frac{a}{p}\right)` for an integer
        `a \pmod p` (`p` a prime number), can be computed
        as sign(p_a), where sign denotes the sign of a permutation and
        p_a the permutation of the residue classes `\pmod p`
        induced by modular multiplication by `a`, provided
        `p` does not divide `a`.

        We verify this in some examples.

        ::

            sage: F = GF(11)
            sage: a = F.multiplicative_generator();a
            2
            sage: plist = [int(a*F(x)) for x in range(1,11)]; plist
            [2, 4, 6, 8, 10, 1, 3, 5, 7, 9]

        This corresponds ot the permutation (1, 2, 4, 8, 5, 10, 9, 7, 3, 6)
        (acting the set `\{1,2,...,10\}`) and to the partition
        [10].

        ::

            sage: p = PermutationGroupElement('(1, 2, 4, 8, 5, 10, 9, 7, 3, 6)')
            sage: p.sign()
            -1
            sage: Partition([10]).sign()
            -1
            sage: kronecker_symbol(11,2)
            -1

        Now replace `2` by `3`::

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

        - http://en.wikipedia.org/wiki/Zolotarev's_lemma
        """
        ans=gap.eval("SignPartition(%s)"%(self))
        return sage_eval(ans)

    def up(self):
        r"""
        Returns a generator for partitions that can be obtained from pi by
        adding a box.

        EXAMPLES::

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
        Returns a list of the partitions that can be formed from the
        partition p by adding a box.

        EXAMPLES::

            sage: Partition([2,1,1]).up_list()
            [[3, 1, 1], [2, 2, 1], [2, 1, 1, 1]]
            sage: Partition([3,2]).up_list()
            [[4, 2], [3, 3], [3, 2, 1]]
        """
        return [p for p in self.up()]

    def down(self):
        r"""
        Returns a generator for partitions that can be obtained from p by
        removing a box.

        EXAMPLES::

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
        Returns a list of the partitions that can be obtained from the
        partition p by removing a box.

        EXAMPLES::

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
        Returns True if partition p1 dominates partitions p2. Otherwise, it
        returns False.

        EXAMPLES::

            sage: p = Partition([3,2])
            sage: p.dominates([3,1])
            True
            sage: p.dominates([2,2])
            True
            sage: p.dominates([2,1,1])
            True
            sage: p.dominates([3,3])
            False
            sage: p.dominates([4])
            False
            sage: Partition([4]).dominates(p)
            False
            sage: Partition([]).dominates([1])
            False
            sage: Partition([]).dominates([])
            True
            sage: Partition([1]).dominates([])
            True
        """
        p1 = self
        sum1 = 0
        sum2 = 0
        min_length = min(len(p1), len(p2))
        if min_length == 0:
            return len(p1) >= len(p2)

        for i in range(min_length):
            sum1 += p1[i]
            sum2 += p2[i]
            if sum2 > sum1:
                return False
        return bool(sum(p1) >= sum(p2))

    def boxes(self):
        """
        Return the coordinates of the boxes of self. Coordinates are given
        as (row-index, column-index) and are 0 based.

        EXAMPLES::

            sage: Partition([2,2]).boxes()
            [(0, 0), (0, 1), (1, 0), (1, 1)]
            sage: Partition([3,2]).boxes()
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
        """
        res = []
        for i in range(len(self)):
            for j in range(self[i]):
                res.append( (i,j) )
        return res

    def generalized_pochhammer_symbol(self, a, alpha):
        r"""
        Returns the generalized Pochhammer symbol
        `(a)_{self}^{(\alpha)}`. This is the product over all
        cells (i,j) in p of `a - (i-1) / \alpha + j - 1`.

        EXAMPLES::

            sage: Partition([2,2]).generalized_pochhammer_symbol(2,1)
            12
        """
        res = 1
        for (i,j) in self.boxes():
            res *= (a - (i-1)/alpha+j-1)
        return res

    def conjugate(self):
        """
        conjugate() returns the "conjugate" (also called "associated" in
        the literature) partition of the partition p which is obtained by
        transposing the corresponding Ferrers diagram.

        EXAMPLES::

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
            return Partition_class([])
        else:
            l = len(p)
            conj =  [l]*p[-1]
            for i in xrange(l-1,0,-1):
                conj.extend([i]*(p[i-1] - p[i]))
            return Partition_class(conj)

    def reading_tableau(self):
        """
        EXAMPLES::

            sage: Partition([3,2,1]).reading_tableau()
            [[1, 3, 6], [2, 5], [4]]
        """
        st = tableau.StandardTableaux(self).first()
        word = st.to_word()
        perm = permutation.Permutation(word)
        return perm.robinson_schensted()[1]

    def associated(self):
        """
        An alias for partition.conjugate(p).

        EXAMPLES::

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
        Returns the arm of cell (i,j) in partition p. The arm of cell (i,j)
        is the number of boxes that appear to the right of cell (i,j). Note
        that i and j are 0-based indices. If your coordinates are in the
        form [i,j], use Python's \*-operator.

        EXAMPLES::

            sage: p = Partition([2,2,1])
            sage: p.arm(0, 0)
            1
            sage: p.arm(0, 1)
            0
            sage: p.arm(2, 0)
            0
            sage: Partition([3,3]).arm(0, 0)
            2
            sage: Partition([3,3]).arm(*[0,0])
            2
        """
        p = self
        if i < len(p) and j < p[i]:
            return p[i]-(j+1)
        else:
            #Error: invalid coordinates
            pass

    def arm_lengths(self, flat=False):
        """
        Returns a tableau of shape p where each box is filled its arm. The
        optional boolean parameter flat provides the option of returning a
        flat list.

        EXAMPLES::

            sage: Partition([2,2,1]).arm_lengths()
            [[1, 0], [1, 0], [0]]
            sage: Partition([2,2,1]).arm_lengths(flat=True)
            [1, 0, 1, 0, 0]
            sage: Partition([3,3]).arm_lengths()
            [[2, 1, 0], [2, 1, 0]]
            sage: Partition([3,3]).arm_lengths(flat=True)
            [2, 1, 0, 2, 1, 0]
        """
        p = self
        res = [[p[i]-(j+1) for j in range(p[i])] for i in range(len(p))]
        if flat:
            return sum(res, [])
        else:
            return res


    def leg(self, i, j):
        """
        Returns the leg of box (i,j) in partition p. The leg of box (i,j)
        is defined to be the number of boxes below it in partition p. Note
        that i and j are 0-based. If your coordinates are in the form
        [i,j], use Python's \*-operator.

        EXAMPLES::

            sage: p = Partition([2,2,1])
            sage: p.leg(0, 0)
            2
            sage: p.leg(0,1)
            1
            sage: p.leg(2,0)
            0
            sage: Partition([3,3]).leg(0, 0)
            1
            sage: cell = [0,0]; Partition([3,3]).leg(*cell)
            1
        """

        conj = self.conjugate()
        if j < len(conj) and i < conj[j]:
            return conj[j]-(i+1)
        else:
            #Error: invalid coordinates
            pass

    def leg_lengths(self, flat=False):
        """
        Returns a tableau of shape p with each box filled in with its leg.
        The optional boolean parameter flat provides the option of
        returning a flat list.

        EXAMPLES::

            sage: Partition([2,2,1]).leg_lengths()
            [[2, 1], [1, 0], [0]]
            sage: Partition([2,2,1]).leg_lengths(flat=True)
            [2, 1, 1, 0, 0]
            sage: Partition([3,3]).leg_lengths()
            [[1, 1, 1], [0, 0, 0]]
            sage: Partition([3,3]).leg_lengths(flat=True)
            [1, 1, 1, 0, 0, 0]
        """
        p = self
        conj = p.conjugate()
        res = [[conj[j]-(i+1) for j in range(p[i])] for i in range(len(p))]
        if flat:
            return sum(res, [])
        else:
            return res

    def dominate(self, rows=None):
        """
        Returns a list of the partitions dominated by n. If n is specified,
        then it only returns the ones with = rows rows.

        EXAMPLES::

            sage: Partition([3,2,1]).dominate()
            [[3, 2, 1], [3, 1, 1, 1], [2, 2, 2], [2, 2, 1, 1], [2, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]]
            sage: Partition([3,2,1]).dominate(rows=3)
            [[3, 2, 1], [2, 2, 2]]
        """
        #Naive implementation
        res = [x for x in Partitions_n(self.size()) if self.dominates(x)]
        if rows:
            return [x for x in res if len(x) <= rows]
        else:
            return res

    def contains(self, x):
        """
        Returns True if p2 is a partition whose Ferrers diagram is
        contained in the Ferrers diagram of self

        EXAMPLES::

            sage: p = Partition([3,2,1])
            sage: p.contains([2,1])
            True
            sage: all(p.contains(mu) for mu in Partitions(3))
            True
            sage: all(p.contains(mu) for mu in Partitions(4))
            False
        """
        return len(self) >= len(x) and all(self[i] >= x[i] for i in range(len(x)))

    def hook_product(self, a):
        """
        Returns the Jack hook-product.

        EXAMPLES::

            sage: Partition([3,2,1]).hook_product(x)
            (x + 2)^2*(2*x + 3)
            sage: Partition([2,2]).hook_product(x)
            2*(x + 1)*(x + 2)
        """

        nu = self.conjugate()
        res = 1
        for i in range(len(self)):
            for j in range(self[i]):
                res *= a*(self[i]-j-1)+nu[j]-i
        return res

    def hook_polynomial(self, q, t):
        """
        Returns the two-variable hook polynomial.

        EXAMPLES::

            sage: R.<q,t> = PolynomialRing(QQ)
            sage: a = Partition([2,2]).hook_polynomial(q,t)
            sage: a == (1 - t)*(1 - q*t)*(1 - t^2)*(1 - q*t^2)
            True
            sage: a = Partition([3,2,1]).hook_polynomial(q,t)
            sage: a == (1 - t)^3*(1 - q*t^2)^2*(1 - q^2*t^3)
            True
        """
        nu = self.conjugate()
        res = 1
        for i in range(len(self)):
            for j in range(self[i]):
                res *= 1-q**(self[i]-j-1)*t**(nu[j]-i)
        return res


    def hook(self, i, j):
        """
        Returns the hook of box (i,j) in the partition p. The hook of box
        (i,j) is defined to be one more than the sum of number of boxes to
        the right and the number of boxes below (in the English
        convention). Note that i and j are 0-based. If your coordinates are
        in the form [i,j], use Python's \*-operator.

        EXAMPLES::

            sage: p = Partition([2,2,1])
            sage: p.hook(0, 0)
            4
            sage: p.hook(0, 1)
            2
            sage: p.hook(2, 0)
            1
            sage: Partition([3,3]).hook(0, 0)
            4
            sage: cell = [0,0]; Partition([3,3]).hook(*cell)
            4
        """
        return self.leg(i,j)+self.arm(i,j)+1

    def hooks(self):
        """
        Returns a sorted list of the hook lengths in self.

        EXAMPLES::

            sage: Partition([3,2,1]).hooks()
            [5, 3, 3, 1, 1, 1]
        """
        res = []
        for row in self.hook_lengths():
            res += row
        res.sort(reverse=True)
        return res

    def hook_lengths(self):
        r"""
        Returns a tableau of shape p with the boxes filled in with the
        hook lengths.

        In each box, put the sum of one plus the number of boxes
        horizontally to the right and vertically below the box (the
        hook length).

        For example, consider the partition [3,2,1] of 6 with Ferrers
        Diagram::

            # # #
            # #
            #

        When we fill in the boxes with the hook lengths, we obtain::

            5 3 1
            3 1
            1

        EXAMPLES::

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

        - http://mathworld.wolfram.com/HookLengthFormula.html
        """
        p = self
        conj = p.conjugate()
        return [[p[i]-(i+1)+conj[j]-(j+1)+1 for j in range(p[i])] for i in range(len(p))]

    def upper_hook(self, i, j, alpha):
        r"""
        Returns the upper hook length of the box (i,j) in self. When alpha
        == 1, this is just the normal hook length. As usual, indices are 0
        based.

        The upper hook length of a box (i,j) in a partition
        `\kappa` is defined by

        .. math::

             h_*^\kappa(i,j) = \kappa_j^\prime-i+\alpha(\kappa_i - j+1).



        EXAMPLES::

            sage: p = Partition([2,1])
            sage: p.upper_hook(0,0,1)
            3
            sage: p.hook(0,0)
            3
            sage: [ p.upper_hook(i,j,x) for i,j in p.boxes() ]
            [2*x + 1, x, x]
        """
        p = self
        conj = self.conjugate()
        return conj[j]-(i+1)+alpha*(p[i]-(j+1)+1)

    def upper_hook_lengths(self, alpha):
        r"""
        Returns the upper hook lengths of the partition. When alpha == 1,
        these are just the normal hook lengths.

        The upper hook length of a box (i,j) in a partition
        `\kappa` is defined by

        .. math::

             h_*^\kappa(i,j) = \kappa_j^\prime-i+1+\alpha(\kappa_i - j).



        EXAMPLES::

            sage: Partition([3,2,1]).upper_hook_lengths(x)
            [[3*x + 2, 2*x + 1, x], [2*x + 1, x], [x]]
            sage: Partition([3,2,1]).upper_hook_lengths(1)
            [[5, 3, 1], [3, 1], [1]]
            sage: Partition([3,2,1]).hook_lengths()
            [[5, 3, 1], [3, 1], [1]]
        """
        p = self
        conj = p.conjugate()
        return [[conj[j]-(i+1)+alpha*(p[i]-(j+1)+1) for j in range(p[i])] for i in range(len(p))]

    def lower_hook(self, i, j, alpha):
        r"""
        Returns the lower hook length of the box (i,j) in self. When alpha
        == 1, this is just the normal hook length. Indices are 0-based.

        The lower hook length of a box (i,j) in a partition
        `\kappa` is defined by

        .. math::

             h_*^\kappa(i,j) = \kappa_j^\prime-i+1+\alpha(\kappa_i - j).



        EXAMPLES::

            sage: p = Partition([2,1])
            sage: p.lower_hook(0,0,1)
            3
            sage: p.hook(0,0)
            3
            sage: [ p.lower_hook(i,j,x) for i,j in p.boxes() ]
            [x + 2, 1, 1]
        """
        p = self
        conj = self.conjugate()
        return conj[j]-(i+1)+1+alpha*(p[i]-(j+1))


    def lower_hook_lengths(self, alpha):
        r"""
        Returns the lower hook lengths of the partition. When alpha == 1,
        these are just the normal hook lengths.

        The lower hook length of a box (i,j) in a partition
        `\kappa` is defined by

        .. math::

             h_\kappa^*(i,j) = \kappa_j^\prime-i+\alpha(\kappa_i - j + 1).



        EXAMPLES::

            sage: Partition([3,2,1]).lower_hook_lengths(x)
            [[2*x + 3, x + 2, 1], [x + 2, 1], [1]]
            sage: Partition([3,2,1]).lower_hook_lengths(1)
            [[5, 3, 1], [3, 1], [1]]
            sage: Partition([3,2,1]).hook_lengths()
            [[5, 3, 1], [3, 1], [1]]
        """
        p = self
        conj = p.conjugate()
        return [[conj[j]-(i+1)+1+alpha*(p[i]-(j+1)) for j in range(p[i])] for i in range(len(p))]


    def weighted_size(self):
        """
        Returns sum([i\*p[i] for i in range(len(p))]). It is also the
        sum of the leg of every cell in b, or the sum of
        binomial(s[i],2) for s the conjugate partition of p.

        EXAMPLES::

            sage: Partition([2,2]).weighted_size()
            2
            sage: Partition([3,3,3]).weighted_size()
            9
            sage: Partition([5,2]).weighted_size()
            2
        """
        p = self
        return sum([i*p[i] for i in range(len(p))])


    def length(self):
        """
        Returns the number of parts in self.

        EXAMPLES::

            sage: Partition([3,2]).length()
            2
            sage: Partition([2,2,1]).length()
            3
            sage: Partition([]).length()
            0
        """
        return len(self)

    def to_exp(self, k=0):
        """
        Return a list of the multiplicities of the parts of a partition.
        Use the optional parameter k to get a return list of length at
        least k.

        EXAMPLES::

            sage: Partition([3,2,2,1]).to_exp()
            [1, 2, 1]
            sage: Partition([3,2,2,1]).to_exp(5)
            [1, 2, 1, 0, 0]
        """
        p = self
        if len(p) > 0:
            k = max(k, p[0])
        a = [ 0 ] * (k)
        for i in p:
            a[i-1] += 1
        return a

    def evaluation(self):
        """
        Returns the evaluation of the partition.

        EXAMPLES::

            sage: Partition([4,3,1,1]).evaluation()
            [2, 0, 1, 1]
        """
        return self.to_exp()

    def to_exp_dict(self):
        """
        Returns a dictionary containing the multiplicities of the
        parts of this partition.

        EXAMPLES::

            sage: p = Partition([4,2,2,1])
            sage: d = p.to_exp_dict()
            sage: d[4]
            1
            sage: d[2]
            2
            sage: d[1]
            1
            sage: 5 in d
            False
        """
        d = {}
        for part in self:
            d[part] = d.get(part, 0) + 1
        return d

    def centralizer_size(self, t=0, q=0):
        """
        Returns the size of the centralizer of any permuation of cycle type
        p. If m_i is the multiplicity of i as a part of p, this is given
        by `\prod_i (i^m[i])*(m[i]!)`. Including the optional
        parameters t and q gives the q-t analog which is the former product
        times `\prod_{i=1}^{length(p)} (1 - q^{p[i]}) / (1 - t^{p[i]}).`

        EXAMPLES::

            sage: Partition([2,2,1]).centralizer_size()
            8
            sage: Partition([2,2,2]).centralizer_size()
            48

        REFERENCES:

        - Kerber, A. 'Algebraic Combintorics via Finite Group Action',
          1.3 p24
        """
        p = self
        a = p.to_exp()
        size = prod([(i+1)**a[i]*factorial(a[i]) for i in range(len(a))])
        size *= prod( [ (1-q**p[i])/(1-t**p[i]) for i in range(len(p)) ] )

        return size

    def aut(self):
        r"""
        Returns a factor for the number of permutations with cycle type
        self. self.aut() returns
        `1^{j_1}j_1! \cdots n^{j_n}j_n!` where `j_k`
        is the number of parts in self equal to k.

        The number of permutations having `p` as a cycle type is
        given by

        .. math::

             \frac{n!}{1^{j_1}j_1! \cdots n^{j_n}j_n!} .

        EXAMPLES::

            sage: Partition([2,1]).aut()
            2
        """
        m = self.to_exp()
        return prod([(i+1)**m[i]*factorial(m[i]) for i in range(len(m)) if m[i] > 0])

    def content(self, i, j):
        r"""
        Returns the content statistic of the given cell, which is simply
        defined by j - i. This doesn't technically depend on the partition,
        but is included here because it is often useful in the context of
        partitions.

        EXAMPLES::

            sage: Partition([2,1]).content(0,1)
            1
            sage: p = Partition([3,2])
            sage: sum([p.content(*c) for c in p.boxes()])
            2
        """

        return j - i

    def conjugacy_class_size(self):
        """
        Returns the size of the conjugacy class of the symmetric group
        indexed by the partition p.

        EXAMPLES::

            sage: Partition([2,2,2]).conjugacy_class_size()
            15
            sage: Partition([2,2,1]).conjugacy_class_size()
            15
            sage: Partition([2,1,1]).conjugacy_class_size()
            6

        REFERENCES:

        - Kerber, A. 'Algebraic Combinatorics via Finite Group Action'
          1.3 p24
        """

        return factorial(sum(self))/self.centralizer_size()


    def corners(self):
        """
        Returns a list of the corners of the partitions. These are the
        positions where we can remove a box. Indices are of the form [i,j]
        where i is the row-index and j is the column-index, and are
        0-based.

        EXAMPLES::

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
        Returns a list of the positions where we can add a box so that the
        shape is still a partition. Indices are of the form [i,j] where i
        is the row-index and j is the column-index, and are 0-based.

        EXAMPLES::

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
        Returns the r-core of the partition p. The construction of the
        r-core can be visualized by repeatedly removing border strips of
        size r from p until this is no longer possible. The remaining
        partition is the r-core.

        EXAMPLES::

            sage: Partition([6,3,2,2]).r_core(3)
            [2, 1, 1]
            sage: Partition([]).r_core(3)
            []
            sage: Partition([8,7,7,4,1,1,1,1,1]).r_core(3)
            [2, 1, 1]

        TESTS::

            sage: Partition([3,3,3,2,1]).r_core(3)
            []
            sage: Partition([10,8,7,7]).r_core(4)
            []
            sage: Partition([21,15,15,9,6,6,6,3,3]).r_core(3)
            []
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
        Returns the r-quotient of the partition p. The r-quotient is a list
        of r partitions, constructed in the following way. Label each cell
        in p with its content, modulo r. Let R_i be the set of rows ending
        in a box labelled i, and C_i be the set of columns ending in a box
        labelled i. Then the jth component of the r-quotient of p is the
        partition defined by intersecting R_j with C_j+1.

        EXAMPLES::

            sage: Partition([7,7,5,3,3,3,1]).r_quotient(3)
            [[2], [1], [2, 2, 2]]

        TESTS::

            sage: Partition([8,7,7,4,1,1,1,1,1]).r_quotient(3)
            [[2, 1], [2, 2], [2]]
            sage: Partition([10,8,7,7]).r_quotient(4)
            [[2], [3], [2], [1]]
            sage: Partition([6,3,3]).r_quotient(3)
            [[1], [1], [2]]
            sage: Partition([3,3,3,2,1]).r_quotient(3)
            [[1], [1, 1], [1]]
            sage: Partition([6,6,6,3,3,3]).r_quotient(3)
            [[2, 1], [2, 1], [2, 1]]
            sage: Partition([21,15,15,9,6,6,6,3,3]).r_quotient(3)
            [[5, 2, 1], [5, 2, 1], [7, 3, 2]]
            sage: Partition([21,15,15,9,6,6,3,3]).r_quotient(3)
            [[5, 2], [5, 2, 1], [7, 3, 1]]
            sage: Partition([14,12,11,10,10,10,10,9,6,4,3,3,2,1]).r_quotient(5)
            [[3, 3], [2, 2, 1], [], [3, 3, 3], [1]]
        """
        p = self
        #Normalize the length
        remainder = len(p) % length
        part = p[:] + [0]*(length-remainder)


        #Add the canonical vector to the partition
        part = [part[i-1] + len(part)-i for i in range(1, len(part)+1)]
        result = [None]*length

        #Reducing vector
        for e in range(length):
            k = e
            tmp = []
            for i in reversed(range(len(part))):
                if part[i] % length == e:
                    tmp.append(ZZ((part[i]-k)/length))
                    k += length

            a = [i for i in tmp if i != 0]
            a.reverse()
            result[e] = a

        return result

    def add_box(self, i, j = None):
        r"""
        Returns a partition corresponding to self with a box added in row
        i. i and j are 0-based row and column indices. This does not change
        p.

        Note that if you have coordinates in a list, you can call this
        function with python's \* notation (see the examples below).

        EXAMPLES::

            sage: Partition([3, 2, 1, 1]).add_box(0)
            [4, 2, 1, 1]
            sage: cell = [4, 0]; Partition([3, 2, 1, 1]).add_box(*cell)
            [3, 2, 1, 1, 1]
        """

        if j is None:
            if i >= len(self):
                j = 0
            else:
                j = self[i]

        if [i,j] in self.outside_corners():
            pl = self.to_list()
            if i == len(pl):
                pl.append(1)
            else:
                pl[i] += 1
            return Partition(pl)

        raise ValueError, "[%s, %s] is not an addable box"%(i,j)

    def remove_box(self, i, j = None):
        """
        Returns the partition obtained by removing a box at the end of row
        i.

        EXAMPLES::

            sage: Partition([2,2]).remove_box(1)
            [2, 1]
            sage: Partition([2,2,1]).remove_box(2)
            [2, 2]
            sage: #Partition([2,2]).remove_box(0)

        ::

            sage: Partition([2,2]).remove_box(1,1)
            [2, 1]
            sage: #Partition([2,2]).remove_box(1,0)
        """

        if i >= len(self):
            raise ValueError, "i must be less than the length of the partition"

        if j is None:
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

        The k-skew diagram of a k-bounded partition is the skew diagram
        denoted `\lambda/^k` satisfying the conditions:

        1. row i of `\lambda/^k` has length `\lambda_i`

        2. no cell in `\lambda/^k` has hook-length exceeding k

        3. every square above the diagram of `\lambda/^k` has hook
           length exceeding k.

        REFERENCES:

        - Lapointe, L. and Morse, J. 'Order Ideals in Weak Subposets
          of Young's Lattice and Associated Unimodality Conjectures'

        EXAMPLES::

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

    def to_list(self):
        r"""
        Return self as a list.

        EXAMPLES::

            sage: p = Partition([2,1]).to_list(); p
            [2, 1]
            sage: type(p)
            <type 'list'>

        TESTS::

            sage: p = Partition([2,1])
            sage: pl = p.to_list()
            sage: pl[0] = 0; p
            [2, 1]
        """
        return self._list[:]

    def add_vertical_border_strip(self, k):
        """
        Returns a list of all the partitions that can be obtained by adding
        a vertical border strip of length k to self.

        EXAMPLES::

            sage: Partition([]).add_vertical_border_strip(0)
            [[]]
            sage: Partition([]).add_vertical_border_strip(2)
            [[1, 1]]
            sage: Partition([2,2]).add_vertical_border_strip(2)
            [[3, 3], [3, 2, 1], [2, 2, 1, 1]]
            sage: Partition([3,2,2]).add_vertical_border_strip(2)
            [[4, 3, 2], [4, 2, 2, 1], [3, 3, 3], [3, 3, 2, 1], [3, 2, 2, 1, 1]]
        """
        return [p.conjugate() for p in self.conjugate().add_horizontal_border_strip(k)]

    def add_horizontal_border_strip(self, k):
        """
        Returns a list of all the partitions that can be obtained by adding
        a horizontal border strip of length k to self.

        EXAMPLES::

            sage: Partition([]).add_horizontal_border_strip(0)
            [[]]
            sage: Partition([]).add_horizontal_border_strip(2)
            [[2]]
            sage: Partition([2,2]).add_horizontal_border_strip(2)
            [[2, 2, 2], [3, 2, 1], [4, 2]]
            sage: Partition([3,2,2]).add_horizontal_border_strip(2)
            [[3, 2, 2, 2], [3, 3, 2, 1], [4, 2, 2, 1], [4, 3, 2], [5, 2, 2]]
        """
        conj = self.conjugate().to_list()
        shelf = []
        res = []
        i = 0
        while i < len(conj):
            tmp = 1
            while i+1 < len(conj) and conj[i] == conj[i+1]:
                tmp += 1
                i += 1
            if i == len(conj)-1 and i > 0 and conj[i] != conj[i-1]:
                tmp = 1
            shelf.append(tmp)
            i += 1

        #added the last shelf on the right side of
        #the first line
        shelf.append(k)

        #list all of the positions for cells
        #filling each self from the left to the right
        l = IntegerVectors(k, len(shelf), outer=shelf).list()
        for iv in l:
            tmp = conj + [0]*k
            j = 0
            for t in range(len(iv)):
                while iv[t] > 0:
                    tmp[j] += 1
                    iv[t] -= 1
                    j += 1
                j = sum(shelf[:t+1])
            res.append(Partition_class([i for i in tmp if i != 0]).conjugate())
        return res




    def k_conjugate(self, k):
        """
        Returns the k-conjugate of the partition.

        The k-conjugate is the partition that is given by the columns of
        the k-skew diagram of the partition.

        EXAMPLES::

            sage: p = Partition([4,3,2,2,1,1])
            sage: p.k_conjugate(4)
            [3, 2, 2, 1, 1, 1, 1, 1, 1]
        """
        return Partition(self.k_skew(k).conjugate().row_lengths())

    def parent(self):
        """
        Returns the combinatorial class of partitions of sum(self).

        EXAMPLES::

            sage: Partition([3,2,1]).parent()
            Partitions of the integer 6
        """
        return Partitions(sum(self[:]))

    def arms_legs_coeff(self, i, j):
        """
        This is a statistic on a cell c=[i,j] in the diagram of partition p
        given by

        .. math::

             [ (1 - q^{arm(c)} * t^{leg(c) + 1}) ] /            [ (1 - q^{arm(c) + 1} * t^{leg(c)}) ]



        EXAMPLES::

            sage: Partition([3,2,1]).arms_legs_coeff(1,1)
            (-t + 1)/(-q + 1)
            sage: Partition([3,2,1]).arms_legs_coeff(0,0)
            (-q^2*t^3 + 1)/(-q^3*t^2 + 1)
            sage: Partition([3,2,1]).arms_legs_coeff(*[0,0])
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

    def atom(self):
        """
        Returns a list of the standard tableaux of size self.size() whose
        atom is equal to self.

        EXAMPLES::

            sage: Partition([2,1]).atom()
            [[[1, 2], [3]]]
            sage: Partition([3,2,1]).atom()
            [[[1, 2, 3, 6], [4, 5]], [[1, 2, 3], [4, 5], [6]]]
        """
        res = []
        for tab in tableau.StandardTableaux_n(self.size()):
            if tab.atom() == self:
                res.append(tab)
        return res

    def k_atom(self, k):
        """
        EXAMPLES::

            sage: p = Partition([3,2,1])
            sage: p.k_atom(1)
            []
            sage: p.k_atom(3)
            [[[1, 1, 1], [2, 2], [3]],
             [[1, 1, 1, 2], [2], [3]],
             [[1, 1, 1, 3], [2, 2]],
             [[1, 1, 1, 2, 3], [2]]]
            sage: Partition([3,2,1]).k_atom(4)
            [[[1, 1, 1], [2, 2], [3]], [[1, 1, 1, 3], [2, 2]]]

        TESTS::

            sage: Partition([1]).k_atom(1)
            [[[1]]]
            sage: Partition([1]).k_atom(2)
            [[[1]]]
            sage: Partition([]).k_atom(1)
            [[]]
        """
        res = [ tableau.Tableau_class([]) ]
        for i in range(len(self)):
            res = [ x.promotion_operator( self[-i-1] ) for x in res]
            res = sum(res, [])
            res = [ y.katabolism_projector(Partition_class(self[-i-1:]).k_split(k)) for y in res]
            res = [ i for i in res if i !=0 and i != [] ]
        return res

    def k_split(self, k):
        """
        Returns the k-split of self.

        EXAMPLES::

            sage: Partition([4,3,2,1]).k_split(3)
            []
            sage: Partition([4,3,2,1]).k_split(4)
            [[4], [3, 2], [1]]
            sage: Partition([4,3,2,1]).k_split(5)
            [[4, 3], [2, 1]]
            sage: Partition([4,3,2,1]).k_split(6)
            [[4, 3, 2], [1]]
            sage: Partition([4,3,2,1]).k_split(7)
            [[4, 3, 2, 1]]
            sage: Partition([4,3,2,1]).k_split(8)
            [[4, 3, 2, 1]]
        """
        if self == []:
            return []
        elif k < self[0]:
            return []
        else:
            res = []
            part = list(self)
            while part != [] and part[0]+len(part)-1 >= k:
                p = k - part[0]
                res.append( part[:p+1] )
                part = part[p+1:]
            if part != []:
                res.append(part)
        return res

    def jacobi_trudi(self):
        """
        EXAMPLES::

            sage: part = Partition([3,2,1])
            sage: jt = part.jacobi_trudi(); jt
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
        return sage.combinat.skew_partition.SkewPartition([ self, [] ]).jacobi_trudi()


    def character_polynomial(self):
        r"""
        Returns the character polynomial associated to the partition self.
        The character polynomial `q_\mu` is defined by


        .. math::

           q_\mu(x_1, x_2, \ldots, x_k) = \downarrow \sum_{\alpha \vdash k}\frac{ \chi^\mu_\alpha }{1^{a_1}2^{a_2}\cdots k^{a_k}a_1!a_2!\cdots a_k!} \prod_{i=1}^{k} (ix_i-1)^{a_i}


        where `a_i` is the multiplicity of `i` in
        `\alpha`.

        It is computed in the following manner.

        1. Expand the Schur function `s_\mu` in the power-sum
           basis.

        2. Replace each `p_i` with `ix_i-1`

        3. Apply the umbral operator `\downarrow` to the resulting
           polynomial.

        EXAMPLES::

            sage: Partition([1]).character_polynomial()
            x - 1
            sage: Partition([1,1]).character_polynomial()
            1/2*x0^2 - 3/2*x0 - x1 + 1
            sage: Partition([2,1]).character_polynomial()
            1/3*x0^3 - 2*x0^2 + 8/3*x0 - x2
        """

        #Create the polynomial ring we will use
        k = self.size()
        P = PolynomialRing(QQ, k, 'x')
        x = P.gens()

        #Expand s_mu in the power sum basis
        s = sf.sfa.SFASchur(QQ)
        p = sf.sfa.SFAPower(QQ)
        ps_mu = p(s(self))

        #Replace each p_i by i*x_i-1
        items = ps_mu.monomial_coefficients().items()  #items contains a list of (partition, coeff) pairs
        partition_to_monomial = lambda part: prod([ (i*x[i-1]-1) for i in part ])
        res = [ [partition_to_monomial(mc[0]), mc[1]] for mc in items ]

        #Write things in the monomial basis
        res = [ prod(pair) for pair in res ]
        res = sum( res )

        #Apply the umbral operator and return the result
        return misc.umbral_operation(res)


##################################################

def number_of_partitions_list(n,k=None):
    r"""
    This function will be deprecated in a future version of Sage and
    eventually removed. Use Partitions(n).count() or Partitions(n,
    length=k).count() instead.

    Original docstring follows.

    Returns the size of partitions_list(n,k).

    Wraps GAP's NrPartitions.

    It is possible to associate with every partition of the integer n a
    conjugacy class of permutations in the symmetric group on n points
    and vice versa. Therefore p(n) = NrPartitions(n) is the number of
    conjugacy classes of the symmetric group on n points.

    ``number_of_partitions(n)`` is also available in
    PARI, however the speed seems the same until `n` is in the
    thousands (in which case PARI is faster).

    EXAMPLES::

        sage: partition.number_of_partitions_list(10,2)
        5
        sage: partition.number_of_partitions_list(10)
        42

    A generating function for p(n) is given by the reciprocal of
    Euler's function:

    .. math::

             \sum_{n=0}^\infty p(n)x^n = \prod_{k=1}^\infty \left(\frac {1}{1-x^k} \right).


    Sage verifies that the first several coefficients do instead
    agree::

        sage: q = PowerSeriesRing(QQ, 'q', default_prec=9).gen()
        sage: prod([(1-q^k)^(-1) for k in range(1,9)])  ## partial product of
        1 + q + 2*q^2 + 3*q^3 + 5*q^4 + 7*q^5 + 11*q^6 + 15*q^7 + 22*q^8 + O(q^9)
        sage: [partition.number_of_partitions_list(k) for k in range(2,10)]
        [2, 3, 5, 7, 11, 15, 22, 30]

    REFERENCES:

    - http://en.wikipedia.org/wiki/Partition\_(number\_theory)
    """
    if k is None:
        ans=gap.eval("NrPartitions(%s)"%(ZZ(n)))
    else:
        ans=gap.eval("NrPartitions(%s,%s)"%(ZZ(n),ZZ(k)))
    return ZZ(ans)


######################
# Ordered Partitions #
######################
def OrderedPartitions(n, k=None):
    """
    Returns the combinatoiral class of ordered partitions of n. If k is
    specified, then only the ordered partitions of length k are
    returned.

    .. note::

       It is recommended that you use Compositions instead as
       OrderedPartitions wraps GAP. See also ordered_partitions.

    EXAMPLES::

        sage: OrderedPartitions(3)
        Ordered partitions of 3
        sage: OrderedPartitions(3).list()
        [[3], [2, 1], [1, 2], [1, 1, 1]]
        sage: OrderedPartitions(3,2)
        Ordered partitions of 3 of length 2
        sage: OrderedPartitions(3,2).list()
        [[2, 1], [1, 2]]
    """
    return OrderedPartitions_nk(n,k)

class OrderedPartitions_nk(CombinatorialClass):
    def __init__(self, n, k=None):
        """
        EXAMPLES::

            sage: o = OrderedPartitions(4,2)
            sage: o == loads(dumps(o))
            True
        """
        self.n = n
        self.k = k

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: o = OrderedPartitions(4,2)
            sage: [2,1] in o
            False
            sage: [2,2] in o
            True
            sage: [1,2,1] in o
            False
        """
        return x in composition.Compositions(self.n, length=self.k)

    def __repr__(self):
        """
        EXAMPLES::

            sage: OrderedPartitions(3).__repr__()
            'Ordered partitions of 3'
            sage: OrderedPartitions(3,2).__repr__()
            'Ordered partitions of 3 of length 2'
        """
        string = "Ordered partitions of %s"%self.n
        if self.k is not None:
            string += " of length %s"%self.k
        return string

    def list(self):
        """
        EXAMPLES::

            sage: OrderedPartitions(3).list()
            [[3], [2, 1], [1, 2], [1, 1, 1]]
            sage: OrderedPartitions(3,2).list()
            [[2, 1], [1, 2]]
        """
        n = self.n
        k = self.k
        if self.k is None:
            ans=gap.eval("OrderedPartitions(%s)"%(ZZ(n)))
        else:
            ans=gap.eval("OrderedPartitions(%s,%s)"%(ZZ(n),ZZ(k)))
        result = eval(ans.replace('\n',''))
        result.reverse()
        return result

    def count(self):
        """
        EXAMPLES::

            sage: OrderedPartitions(3).count()
            4
            sage: OrderedPartitions(3,2).count()
            2
        """
        n = self.n
        k = self.k
        if k is None:
            ans=gap.eval("NrOrderedPartitions(%s)"%(n))
        else:
            ans=gap.eval("NrOrderedPartitions(%s,%s)"%(n,k))
        return ZZ(ans)

##########################
# Partitions Greatest LE #
##########################
def PartitionsGreatestLE(n,k):
    """
    Returns the combinatorial class of all (unordered) "restricted"
    partitions of the integer n having parts less than or equal to the
    integer k.

    EXAMPLES::

        sage: PartitionsGreatestLE(10,2)
        Partitions of 10 having parts less than or equal to 2
        sage: PartitionsGreatestLE(10,2).list()
        [[2, 2, 2, 2, 2],
         [2, 2, 2, 2, 1, 1],
         [2, 2, 2, 1, 1, 1, 1],
         [2, 2, 1, 1, 1, 1, 1, 1],
         [2, 1, 1, 1, 1, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
    """
    return PartitionsGreatestLE_nk(n,k)

class PartitionsGreatestLE_nk(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS::

            sage: p = PartitionsGreatestLE(10,2)
            sage: p == loads(dumps(p))
            True
        """
        self.n = n
        self.k = k
        self.object_class = Partition_class
        self._name="Partitions of %s having parts less than or equal to %s"%(self.n, self.k)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [4,3,2,1] in PartitionsGreatestLE(10,2)
            False
            sage: [2,2,2,2,2] in PartitionsGreatestLE(10,2)
            True
        """
        return x in Partitions_n(self.n) and max(x) <= self.k

    def list(self):
        """
        Returns a list of all (unordered) "restricted" partitions of the
        integer n having parts less than or equal to the integer k.

        Wraps GAP's PartitionsGreatestLE.

        EXAMPLES::

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
    """
    Returns combinatorial class of all (unordered) "restricted"
    partitions of the integer n having its greatest part equal to the
    integer k.

    EXAMPLES::

        sage: PartitionsGreatestEQ(10,2)
        Partitions of 10 having greatest part equal to 2
        sage: PartitionsGreatestEQ(10,2).list()
        [[2, 2, 2, 2, 2],
         [2, 2, 2, 2, 1, 1],
         [2, 2, 2, 1, 1, 1, 1],
         [2, 2, 1, 1, 1, 1, 1, 1],
         [2, 1, 1, 1, 1, 1, 1, 1, 1]]
    """
    return PartitionsGreatestEQ_nk(n,k)

class PartitionsGreatestEQ_nk(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS::

            sage: p = PartitionsGreatestEQ(10,2)
            sage: p == loads(dumps(p))
            True
        """
        self.n = n
        self.k = k
        self._name = "Partitions of %s having greatest part equal to %s"%(self.n, self.k)
        self.object_class = Partition_class

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [4,3,2,1] in PartitionsGreatestEQ(10,2)
            False
            sage: [2,2,2,2,2] in PartitionsGreatestEQ(10,2)
            True
            sage: [1]*10 in PartitionsGreatestEQ(10,2)
            False
        """
        return x in Partitions_n(self.n) and max(x) == self.k

    def list(self):
        """
        Wraps GAP's PartitionsGreatestEQ.

        EXAMPLES::

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
    r"""
    A restricted partition is, like an ordinary partition, an unordered
    sum `n = p_1+p_2+\ldots+p_k` of positive integers and is
    represented by the list `p = [p_1,p_2,\ldots,p_k]`, in
    nonincreasing order. The difference is that here the `p_i`
    must be elements from the set `S`, while for ordinary
    partitions they may be elements from `[1..n]`.

    Returns the list of all restricted partitions of the positive
    integer n into sums with k summands with the summands of the
    partition coming from the set S. If k is not given all restricted
    partitions for all k are returned.

    Wraps GAP's RestrictedPartitions.

    EXAMPLES::

        sage: RestrictedPartitions(5,[3,2,1])
        Partitions of 5 restricted to the values [1, 2, 3]
        sage: RestrictedPartitions(5,[3,2,1]).list()
        [[3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
        sage: RestrictedPartitions(5,[3,2,1],4)
        Partitions of 5 restricted to the values [1, 2, 3] of length 4
        sage: RestrictedPartitions(5,[3,2,1],4).list()
        [[2, 1, 1, 1]]
    """
    return RestrictedPartitions_nsk(n, S, k)

class RestrictedPartitions_nsk(CombinatorialClass):
    def __init__(self, n, S, k=None):
        """
        TESTS::

            sage: r = RestrictedPartitions(5,[3,2,1])
            sage: r == loads(dumps(r))
            True
        """
        self.n = n
        self.S = S
        self.S.sort()
        self.k = k
        self.object_class = Partition_class

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [4,1] in RestrictedPartitions(5,[3,2,1])
            False
            sage: [3,2] in RestrictedPartitions(5,[3,2,1])
            True
            sage: [3,2] in RestrictedPartitions(5,[3,2,1],4)
            False
            sage: [2,1,1,1] in RestrictedPartitions(5,[3,2,1],4)
            True
        """
        return x in Partitions_n(self.n) and all(i in self.S for i in x) \
               and (self.k is None or len(x) == self.k)

    def __repr__(self):
        """
        EXAMPLES::

            sage: RestrictedPartitions(5,[3,2,1]).__repr__()
            'Partitions of 5 restricted to the values [1, 2, 3] '
            sage: RestrictedPartitions(5,[3,2,1],4).__repr__()
            'Partitions of 5 restricted to the values [1, 2, 3] of length 4'
        """
        string = "Partitions of %s restricted to the values %s "%(self.n, self.S)
        if self.k is not None:
            string += "of length %s" % self.k
        return string

    def list(self):
        r"""
        Returns the list of all restricted partitions of the positive
        integer n into sums with k summands with the summands of the
        partition coming from the set `S`. If k is not given all
        restricted partitions for all k are returned.

        Wraps GAP's RestrictedPartitions.

        EXAMPLES::

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
        if k is None:
            ans=gap.eval("RestrictedPartitions(%s,%s)"%(n,S))
        else:
            ans=gap.eval("RestrictedPartitions(%s,%s,%s)"%(n,S,k))
        result = eval(ans)
        result.reverse()
        return [Partition(p) for p in result]

    def count(self):
        """
        Returns the size of RestrictedPartitions(n,S,k). Wraps GAP's
        NrRestrictedPartitions.

        EXAMPLES::

            sage: RestrictedPartitions(8,[1,3,5,7]).count()
            6
            sage: RestrictedPartitions(8,[1,3,5,7],2).count()
            2
        """
        n = self.n
        k = self.k
        S = self.S
        if k is None:
            ans=gap.eval("NrRestrictedPartitions(%s,%s)"%(ZZ(n),S))
        else:
            ans=gap.eval("NrRestrictedPartitions(%s,%s,%s)"%(ZZ(n),S,ZZ(k)))
        return ZZ(ans)

####################
# Partition Tuples #
####################
def PartitionTuples(n,k):
    """
    Returns the combinatorial class of k-tuples of partitions of n.
    These are are ordered list of k partitions whose sizes add up to
    n.

    These describe the classes and the characters of wreath products of
    groups with k conjugacy classes with the symmetric group
    `S_n`.

    EXAMPLES::

        sage: PartitionTuples(4,2)
        2-tuples of partitions of 4
        sage: PartitionTuples(3,2).list()
        [[[3], []],
         [[2, 1], []],
         [[1, 1, 1], []],
         [[2], [1]],
         [[1, 1], [1]],
         [[1], [2]],
         [[1], [1, 1]],
         [[], [3]],
         [[], [2, 1]],
         [[], [1, 1, 1]]]
    """
    return PartitionTuples_nk(n,k)

class PartitionTuples_nk(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS::

            sage: p = PartitionTuples(4,2)
            sage: p == loads(dumps(p))
            True
        """
        self.n = n
        self.k = k

    object_class = Partition_class


    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[], [1, 1, 1]] in PartitionTuples(3,2)
            True
            sage: [[1], [1, 1, 1]] in PartitionTuples(3,2)
            False
        """
        p = Partitions_all()
        return len(x) == self.k and all(i in p for i in x)\
               and sum(sum(i) for i in x) == self.n

    def __repr__(self):
        """
        EXAMPLES::

            sage: PartitionTuples(4,2).__repr__()
            '2-tuples of partitions of 4'
        """
        return "%s-tuples of partitions of %s"%(self.k, self.n)

    def iterator(self):
        r"""
        Returns an iterator for all k-tuples of partitions which together
        form a partition of n.

        EXAMPLES::

            sage: PartitionTuples(2,0).list() #indirect doctest
            []
            sage: PartitionTuples(2,1).list() #indirect doctest
            [[[2]], [[1, 1]]]
            sage: PartitionTuples(2,2).list() #indirect doctest
            [[[2], []], [[1, 1], []], [[1], [1]], [[], [2]], [[], [1, 1]]]
            sage: PartitionTuples(3,2).list() #indirect doctest
            [[[3], []], [[2, 1], []], [[1, 1, 1], []], [[2], [1]], [[1, 1], [1]], [[1], [2]], [[1], [1, 1]], [[], [3]], [[], [2, 1]], [[], [1, 1, 1]]]
        """
        p = [Partitions(i) for i in range(self.n+1)]
        for iv in IntegerVectors(self.n,self.k):
            for cp in CartesianProduct(*[p[i] for i in iv]):
                yield cp


    def count(self):
        r"""
        Returns the number of k-tuples of partitions which together form a
        partition of n.

        Wraps GAP's NrPartitionTuples.

        EXAMPLES::

            sage: PartitionTuples(3,2).count()
            10
            sage: PartitionTuples(8,2).count()
            185

        Now we compare that with the result of the following GAP
        computation::

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
        """
        return ZZ(gp.eval('polcoeff(1/eta(x)^%s, %s, x)'%(self.k, self.n)))


##############
# Partitions #
##############

def Partitions(n=None, **kwargs):
    """
    Partitions(n, \*\*kwargs) returns the combinatorial class of
    integer partitions of n, subject to the constraints given by the
    keywords.

    Valid keywords are: starting, ending, min_part, max_part,
    max_length, min_length, length, max_slope, min_slope, inner,
    outer. They have the following meanings:

    - ``starting=p`` specifies that the partitions should all be = p in
      reverse lex order.

    - ``length=k`` specifies that the partitions have
      exactly k parts.

    - ``min_length=k`` specifies that the partitions have
      at least k parts.

    - ``min_part=k`` specifies that all parts of the
      partitions are at least k.

    - ``outer=p`` specifies that the partitions
      be contained inside the partition p.

    - ``min_slope=k`` specifies that the partition have slope at least
      k; the slope is the difference between successive parts.

    The ``max_*`` versions, along with 'inner' and 'ending', work
    analogously.

    EXAMPLES: If no arguments are passed, then the combinatorial class
    of all integer partitions is returned.

    ::

        sage: Partitions()
        Partitions
        sage: [2,1] in Partitions()
        True

    If an integer n is passed, then the combinatorial class of intger
    partitions of n is returned.

    ::

        sage: Partitions(3)
        Partitions of the integer 3
        sage: Partitions(3).list()
        [[3], [2, 1], [1, 1, 1]]

    If starting is passed, then the combinatorial class of partitions
    starting with starting in lexicographic order is returned.

    ::

        sage: Partitions(3, starting=[2,1])
        Partitions of the integer 3 starting with [2, 1]
        sage: Partitions(3, starting=[2,1]).list()
        [[2, 1], [1, 1, 1]]

    If ending is passed, then the combinatorial class of partitions
    starting with ending in lexicographic order is returned.

    ::

        sage: Partitions(3, ending=[2,1])
        Partitions of the integer 3 ending with [2, 1]
        sage: Partitions(3, ending=[2,1]).list()
        [[3], [2, 1]]

    ::

        sage: Partitions(5,min_part=2)
        Partitions of the integer 5 satisfying constraints min_part=2
        sage: Partitions(5,min_part=2).list()
        [[5], [3, 2]]

    ::

        sage: Partitions(3,max_length=2).list()
        [[3], [2, 1]]

    ::

        sage: Partitions(10, min_part=2, length=3).list()
        [[6, 2, 2], [5, 3, 2], [4, 4, 2], [4, 3, 3]]
    """
    if n is None:
        return Partitions_all()
    else:
        if len(kwargs) == 0:
            if isinstance(n, (int,Integer)):
                return Partitions_n(n)
            else:
                raise ValueError, "n must be an integer"
        else:
            if 'starting' in kwargs:
                return Partitions_starting(n, kwargs['starting'])
            elif 'ending' in kwargs:
                return Partitions_ending(n, kwargs['ending'])
            else:
                return Partitions_constraints(n, **kwargs)

class Partitions_all(CombinatorialClass):
    def __init__(self):
        """
        TESTS::

            sage: P = Partitions()
            sage: P == loads(dumps(P))
            True
        """
        pass

    object_class = Partition_class


    def count(self):
        """
        Returns the number of integer partitions.

        EXAMPLES::

            sage: Partitions().count()
            +Infinity
        """
        return infinity

    def __contains__(self, x):
        """
        TESTS::

            sage: P = Partitions()
            sage: Partition([2,1]) in P
            True
            sage: [2,1] in P
            True
            sage: [3,2,1] in P
            True
            sage: [1,2] in P
            False
            sage: [] in P
            True
            sage: [0] in P
            False
        """
        if isinstance(x, Partition_class):
            return True
        elif isinstance(x, __builtin__.list):
            for i in range(len(x)):
                if not isinstance(x[i], (int, Integer)):
                    return False
                if x[i] <= 0:
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

    def __repr__(self):
        """
        TESTS::

            sage: repr(Partitions())
            'Partitions'
        """
        return "Partitions"

    def list(self):
        """
        EXAMPLES::

            sage: Partitions().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

class Partitions_constraints(CombinatorialClass):
    def __init__(self, n, **kwargs):
        """
        TESTS::

            sage: p = Partitions(5, min_part=2)
            sage: p == loads(dumps(p))
            True
        """
        self.n = n
        self.constraints = kwargs

    def __contains__(self, x):
        """
        TESTS::

            sage: p = Partitions(5, min_part=2)
            sage: [2,1] in p
            False
            sage: [2,2,1] in p
            False
            sage: [3,2] in p
            True
        """
        return x in Partitions_all() and sum(x)==self.n and misc.check_integer_list_constraints(x, singleton=True, **self.constraints)

    def __repr__(self):
        """
        TESTS::

            sage: repr( Partitions(5, min_part=2) )
            'Partitions of the integer 5 satisfying constraints min_part=2'
        """
        return "Partitions of the integer %s satisfying constraints %s"%(self.n, ", ".join( ["%s=%s"%(key, self.constraints[key]) for key in sorted(self.constraints.keys())] ))

    def iterator(self):
        """
        An iterator a list of the partitions of n.

        EXAMPLES::

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
            [[2, 1, 1]]
            sage: [x for x in Partitions(4, min_length=3, min_part=0)]
            [[2, 1, 1], [1, 1, 1, 1]]
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
        for p in Partitions_n(self.n)._fast_iterator():
            if misc.check_integer_list_constraints([p], **self.constraints):
                yield Partition_class(p)


class Partitions_n(CombinatorialClass):
    object_class = Partition_class
    def __init__(self, n):
        """
        TESTS::

            sage: p = Partitions(5)
            sage: p == loads(dumps(p))
            True
        """
        self.n = n

    def __contains__(self, x):
        """
        TESTS::

            sage: p = Partitions(5)
            sage: [2,1] in p
            False
            sage: [2,2,1] in p
            True
            sage: [3,2] in p
            True
        """
        return x in Partitions_all() and sum(x)==self.n

    def __repr__(self):
        """
        TESTS::

            sage: repr( Partitions(5) )
            'Partitions of the integer 5'
        """
        return "Partitions of the integer %s"%self.n

    def count(self, algorithm='default'):
        r"""
        INPUT:

        - ``algorithm``

          - (default: ``default``)

          - ``'bober'`` - use Jonathan Bober's implementation (*very*
            fast, but relatively new)

          - ``'gap'`` - use GAP (VERY *slow*)

          - ``'pari'`` - use PARI. Speed seems the same as GAP until
            `n` is in the thousands, in which case PARI is
            faster. *But* PARI has a bug, e.g., on 64-bit Linux
            PARI-2.3.2 outputs numbpart(147007)%1000 as 536, but it
            should be 533!.  So do not use this option.

          - ``'default'`` - ``'bober'`` when k is not specified;
            otherwise use ``'gap'``.

        Use the function ``partitions(n)`` to return a
        generator over all partitions of `n`.

        It is possible to associate with every partition of the integer n a
        conjugacy class of permutations in the symmetric group on n points
        and vice versa. Therefore p(n) = NrPartitions(n) is the number of
        conjugacy classes of the symmetric group on n points.

        EXAMPLES::

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

        ::

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


        .. math::

           \sum_{n=0}^\infty p(n)x^n = \prod_{k=1}^\infty \left(\frac {1}{1-x^k} \right).



        We use Sage to verify that the first several coefficients do
        instead agree::

            sage: q = PowerSeriesRing(QQ, 'q', default_prec=9).gen()
            sage: prod([(1-q^k)^(-1) for k in range(1,9)])  ## partial product of
            1 + q + 2*q^2 + 3*q^3 + 5*q^4 + 7*q^5 + 11*q^6 + 15*q^7 + 22*q^8 + O(q^9)
            sage: [Partitions(k) .count()for k in range(2,10)]
            [2, 3, 5, 7, 11, 15, 22, 30]

        REFERENCES:

        - http://en.wikipedia.org/wiki/Partition\_(number\_theory)
        """
        return number_of_partitions(self.n, algorithm=algorithm)

    def first(self):
        """
        Returns the lexicographically first partition of a positive integer
        n. This is the partition [n].

        EXAMPLES::

            sage: Partitions(4).first()
            [4]
        """
        return Partition([self.n])

    def last(self):
        """
        Returns the lexicographically last partition of the positive
        integer n. This is the all-ones partition.

        EXAMPLES::

            sage: Partitions(4).last()
            [1, 1, 1, 1]
        """

        return Partition_class([1]*self.n)


    def iterator(self):
        """
        An iterator a list of the partitions of n.

        EXAMPLES::

            sage: [x for x in Partitions(4)]
            [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
        """
        for p in self._fast_iterator():
            yield Partition_class(p)


    def _fast_iterator(self):
        """
        A fast iterator for the partitions of n which returns lists and not
        partition types.

        EXAMPLES::

            sage: p = Partitions(4)
            sage: it = p._fast_iterator()
            sage: it.next()
            [4]
            sage: type(_)
            <type 'list'>
        """
        # base case of the recursion: zero is the sum of the empty tuple
        if self.n == 0:
            yield []
            return

        # modify the partitions of n-1 to form the partitions of n
        for p in Partitions_n(self.n-1)._fast_iterator():
            if p and (len(p) < 2 or p[-2] > p[-1]):
                yield p[:-1] + [p[-1] + 1]
            yield p + [1]

class Partitions_starting(CombinatorialClass):
    def __init__(self, n, starting_partition):
        """
        EXAMPLES::

            sage: Partitions(3, starting=[2,1])
            Partitions of the integer 3 starting with [2, 1]
            sage: Partitions(3, starting=[2,1]).list()
            [[2, 1], [1, 1, 1]]

        TESTS::

            sage: p = Partitions(3, starting=[2,1])
            sage: p == loads(dumps(p))
            True
        """
        self.n = n
        self._starting = Partition(starting_partition)

    def __repr__(self):
        """
        EXAMPLES::

            sage: Partitions(3, starting=[2,1]).__repr__()
            'Partitions of the integer 3 starting with [2, 1]'
        """
        return "Partitions of the integer %s starting with %s"%(self.n, self._starting)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: p = Partitions(3, starting=[2,1])
            sage: [1,1] in p
            False
            sage: [2,1] in p
            True
            sage: [1,1,1] in p
            True
            sage: [3] in p
            False
        """
        return x in Partitions_n(self.n) and x <= self._starting

    def first(self):
        """
        EXAMPLES::

            sage: Partitions(3, starting=[2,1]).first()
            [2, 1]
        """
        return self._starting

    def next(self, part):
        """
        EXAMPLES::

            sage: Partitions(3, starting=[2,1]).next(Partition([2,1]))
            [1, 1, 1]
        """
        return part.next()

class Partitions_ending(CombinatorialClass):
    def __init__(self, n, ending_partition):
        """
        EXAMPLES::

            sage: Partitions(4, ending=[1,1,1,1]).list()
            [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
            sage: Partitions(4, ending=[2,2]).list()
            [[4], [3, 1], [2, 2]]
            sage: Partitions(4, ending=[4]).list()
            [[4]]

        TESTS::

            sage: p = Partitions(4, ending=[1,1,1,1])
            sage: p == loads(dumps(p))
            True
        """
        self.n = n
        self._ending = Partition(ending_partition)

    def __repr__(self):
        """
        EXAMPLES::

            sage: Partitions(4, ending=[1,1,1,1]).__repr__()
            'Partitions of the integer 4 ending with [1, 1, 1, 1]'
        """
        return "Partitions of the integer %s ending with %s"%(self.n, self._ending)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: p = Partitions(4, ending=[2,2])
            sage: [4] in p
            True
            sage: [2,1,1] in p
            False
            sage: [2,1] in p
            False
        """
        return x in Partitions_n(self.n) and x >= self._ending

    def first(self):
        """
        EXAMPLES::

            sage: Partitions(4, ending=[1,1,1,1]).first()
            [4]
        """
        return Partition_class([self.n])

    def next(self, part):
        """
        EXAMPLES::

            sage: Partitions(4, ending=[1,1,1,1]).next(Partition([4]))
            [3, 1]
            sage: Partitions(4, ending=[1,1,1,1]).next(Partition([1,1,1,1])) is None
            True
        """
        if part == self._ending:
            return None
        else:
            return part.next()


def PartitionsInBox(h, w):
    """
    Returns the combinatorial class of partitions that fit in a h by w
    box.

    EXAMPLES::

        sage: PartitionsInBox(2,2)
        Integer partitions which fit in a 2 x 2 box
        sage: PartitionsInBox(2,2).list()
        [[], [1], [1, 1], [2], [2, 1], [2, 2]]
    """
    return PartitionsInBox_hw(h, w)

class PartitionsInBox_hw(CombinatorialClass):
    def __init__(self, h, w):
        """
        TESTS::

            sage: p = PartitionsInBox(2,2)
            sage: p == loads(dumps(p))
            True
        """
        self.h = h
        self.w = w
        self._name = "Integer partitions which fit in a %s x %s box" % (self.h, self.w)
        self._object_class = Partition_class

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [] in PartitionsInBox(2,2)
            True
            sage: [2,1] in PartitionsInBox(2,2)
            True
            sage: [3,1] in PartitionsInBox(2,2)
            False
            sage: [2,1,1] in PartitionsInBox(2,2)
            False
        """
        return x in Partitions_all() and len(x) <= self.h \
               and (len(x) == 0 or x[0] <= self.w)

    def list(self):
        """
        Returns a list of all the partitions inside a box of height h and
        width w.

        EXAMPLES::

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


#########################################################################

#### partitions

def partitions_set(S,k=None, use_file=True):
    r"""
    An unordered partition of a set `S` is a set of pairwise
    disjoint nonempty subsets with union `S` and is represented
    by a sorted list of such subsets.

    partitions_set returns the list of all unordered partitions of the
    list `S` of increasing positive integers into k pairwise
    disjoint nonempty sets. If k is omitted then all partitions are
    returned.

    The Bell number `B_n`, named in honor of Eric Temple Bell,
    is the number of different partitions of a set with n elements.

    .. warning::

       Wraps GAP - hence S must be a list of objects that have string
       representations that can be interpreted by the GAP
       intepreter. If mset consists of at all complicated Sage
       objects, this function does *not* do what you expect. See
       SetPartitions in ``combinat/set_partition``.

    .. warning::

       This function is inefficient. The runtime is dominated by
       parsing the output from GAP.

    Wraps GAP's PartitionsSet.

    EXAMPLES::

        sage: S = [1,2,3,4]
        sage: partitions_set(S,2)
        [[[1], [2, 3, 4]],
         [[1, 2], [3, 4]],
         [[1, 2, 3], [4]],
         [[1, 2, 4], [3]],
         [[1, 3], [2, 4]],
         [[1, 3, 4], [2]],
         [[1, 4], [2, 3]]]

    REFERENCES:

    - http://en.wikipedia.org/wiki/Partition_of_a_set
    """
    if k is None:
        ans=gap("PartitionsSet(%s)"%S).str(use_file=use_file)
    else:
        ans=gap("PartitionsSet(%s,%s)"%(S,k)).str(use_file=use_file)
    return eval(ans)

def number_of_partitions_set(S,k):
    r"""
    Returns the size of ``partitions_set(S,k)``. Wraps
    GAP's NrPartitionsSet.

    The Stirling number of the second kind is the number of partitions
    of a set of size n into k blocks.

    EXAMPLES::

        sage: mset = [1,2,3,4]
        sage: number_of_partitions_set(mset,2)
        7
        sage: stirling_number2(4,2)
        7

    REFERENCES

    - http://en.wikipedia.org/wiki/Partition_of_a_set
    """
    if k is None:
        ans=gap.eval("NrPartitionsSet(%s)"%S)
    else:
        ans=gap.eval("NrPartitionsSet(%s,%s)"%(S,ZZ(k)))
    return ZZ(ans)

def partitions_list(n,k=None):
    r"""
    This function will be deprecated in a future version of Sage and
    eventually removed. Use Partitions(n).list() or Partitions(n,
    length=k).list() instead.

    Original docstring follows.

    An unordered partition of `n` is an unordered sum
    `n = p_1+p_2 +\ldots+ p_k` of positive integers and is
    represented by the list `p = [p_1,p_2,\ldots,p_k]`, in
    nonincreasing order, i.e., `p1\geq p_2 ...\geq p_k`.

    INPUT:


    -  ``n, k`` - positive integer


    ``partitions_list(n,k)`` returns the list of all
    (unordered) partitions of the positive integer n into sums with k
    summands. If k is omitted then all partitions are returned.

    Do not call partitions_list with an n much larger than 40, in
    which case there are 37338 partitions, since the list will simply
    become too large.

    Wraps GAP's Partitions.

    EXAMPLES::

        sage: partitions_list(10,2)
        [[5, 5], [6, 4], [7, 3], [8, 2], [9, 1]]
        sage: partitions_list(5)
        [[1, 1, 1, 1, 1], [2, 1, 1, 1], [2, 2, 1], [3, 1, 1], [3, 2], [4, 1], [5]]
    """
    n = ZZ(n)
    if n <= 0:
        raise ValueError, "n (=%s) must be a positive integer"%n
    if k is None:
        ans=gap.eval("Partitions(%s)"%(n))
    else:
        ans=gap.eval("Partitions(%s,%s)"%(n,k))
    return eval(ans.replace('\n',''))

def number_of_partitions(n,k=None, algorithm='default'):
    r"""
    Returns the size of partitions_list(n,k).

    INPUT:


    -  ``n`` - an integer

    -  ``k`` - (default: None); if specified, instead
       returns the cardinality of the set of all (unordered) partitions of
       the positive integer n into sums with k summands.

    -  ``algorithm`` - (default: 'default')

    -  ``'default'`` - If k is not None, then use Gap (very
       slow). If k is None, use Jonathan Bober's highly optimized
       implementation (this is the fastest code in the world for this
       problem).

    -  ``'bober'`` - use Jonathan Bober's implementation

    -  ``'gap'`` - use GAP (VERY *slow*)

    -  ``'pari'`` - use PARI. Speed seems the same as GAP until
       `n` is in the thousands, in which case PARI is
       faster. *But* PARI has a bug, e.g., on 64-bit Linux
       PARI-2.3.2 outputs numbpart(147007)%1000 as 536 when it should
       be 533!. So do not use this option.


    IMPLEMENTATION: Wraps GAP's NrPartitions or PARI's numbpart
    function.

    Use the function ``partitions(n)`` to return a
    generator over all partitions of `n`.

    It is possible to associate with every partition of the integer n a
    conjugacy class of permutations in the symmetric group on n points
    and vice versa. Therefore p(n) = NrPartitions(n) is the number of
    conjugacy classes of the symmetric group on n points.

    EXAMPLES::

        sage: v = list(partitions(5)); v
        [(1, 1, 1, 1, 1), (1, 1, 1, 2), (1, 2, 2), (1, 1, 3), (2, 3), (1, 4), (5,)]
        sage: len(v)
        7
        sage: number_of_partitions(5, algorithm='gap')
        7
        sage: number_of_partitions(5, algorithm='pari')
        7
        sage: number_of_partitions(5, algorithm='bober')
        7

    The input must be a nonnegative integer or a ValueError is raised.

    ::

        sage: number_of_partitions(-5)
        Traceback (most recent call last):
        ...
        ValueError: n (=-5) must be a nonnegative integer

    ::

        sage: number_of_partitions(10,2)
        5
        sage: number_of_partitions(10)
        42
        sage: number_of_partitions(3)
        3
        sage: number_of_partitions(10)
        42
        sage: number_of_partitions(3, algorithm='pari')
        3
        sage: number_of_partitions(10, algorithm='pari')
        42
        sage: number_of_partitions(40)
        37338
        sage: number_of_partitions(100)
        190569292
        sage: number_of_partitions(100000)
        27493510569775696512677516320986352688173429315980054758203125984302147328114964173055050741660736621590157844774296248940493063070200461792764493033510116079342457190155718943509725312466108452006369558934464248716828789832182345009262853831404597021307130674510624419227311238999702284408609370935531629697851569569892196108480158600569421098519

    A generating function for p(n) is given by the reciprocal of
    Euler's function:


    .. math::

             \sum_{n=0}^\infty p(n)x^n = \prod_{k=1}^\infty \left(\frac {1}{1-x^k} \right).



    We use Sage to verify that the first several coefficients do
    instead agree::

        sage: q = PowerSeriesRing(QQ, 'q', default_prec=9).gen()
        sage: prod([(1-q^k)^(-1) for k in range(1,9)])  ## partial product of
        1 + q + 2*q^2 + 3*q^3 + 5*q^4 + 7*q^5 + 11*q^6 + 15*q^7 + 22*q^8 + O(q^9)
        sage: [number_of_partitions(k) for k in range(2,10)]
        [2, 3, 5, 7, 11, 15, 22, 30]

    REFERENCES:

    - http://en.wikipedia.org/wiki/Partition\_(number\_theory)

    TESTS::

        sage: n = 500 + randint(0,500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1500 + randint(0,1500)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000 + randint(0,1000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 100000000 + randint(0,100000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0
        True
        sage: n = 1000000000 + randint(0,1000000000)
        sage: number_of_partitions( n - (n % 385) + 369) % 385 == 0      # takes a long time
        True

    Another consistency test for n up to 500::

        sage: len([n for n in [1..500] if number_of_partitions(n) != number_of_partitions(n,algorithm='pari')])
        0
    """
    n = ZZ(n)
    if n < 0:
        raise ValueError, "n (=%s) must be a nonnegative integer"%n
    elif n == 0:
        return ZZ(1)

    if algorithm == 'default':
        if k is None:
            algorithm = 'bober'
        else:
            algorithm = 'gap'

    if algorithm == 'gap':
        if k is None:
            ans=gap.eval("NrPartitions(%s)"%(ZZ(n)))
        else:
            ans=gap.eval("NrPartitions(%s,%s)"%(ZZ(n),ZZ(k)))
        return ZZ(ans)

    if k is not None:
        raise ValueError, "only the GAP algorithm works if k is specified."

    if algorithm == 'bober':
        return partitions_ext.number_of_partitions(n)

    elif algorithm == 'pari':
        return ZZ(pari(ZZ(n)).numbpart())

    raise ValueError, "unknown algorithm '%s'"%algorithm

def partitions(n):
    r"""
    Generator of all the partitions of the integer `n`.

    INPUT:


    -  ``n`` - int


    To compute the number of partitions of `n` use
    ``number_of_partitions(n)``.

    EXAMPLES::

        sage: partitions(3)
        <generator object at 0x...>
        sage: list(partitions(3))
        [(1, 1, 1), (1, 2), (3,)]

    AUTHORS:

    - Adapted from David Eppstein, Jan Van lent, George Yoshida;
      Python Cookbook 2, Recipe 19.16.
    """
    n = ZZ(n)
    # base case of the recursion: zero is the sum of the empty tuple
    if n == 0:
        yield ( )
        return
    # modify the partitions of n-1 to form the partitions of n
    for p in partitions(n-1):
        yield (1,) + p
        if p and (len(p) < 2 or p[1] > p[0]):
            yield (p[0] + 1,) + p[1:]

def cyclic_permutations_of_partition(partition):
    """
    Returns all combinations of cyclic permutations of each cell of the
    partition.

    AUTHORS:

    - Robert L. Miller

    EXAMPLES::

        sage: from sage.combinat.partition import cyclic_permutations_of_partition
        sage: cyclic_permutations_of_partition([[1,2,3,4],[5,6,7]])
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

    Note that repeated elements are not considered equal::

        sage: cyclic_permutations_of_partition([[1,2,3],[4,4,4]])
        [[[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]],
         [[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]]]
    """
    return list(cyclic_permutations_of_partition_iterator(partition))

def cyclic_permutations_of_partition_iterator(partition):
    """
    Iterates over all combinations of cyclic permutations of each cell
    of the partition.

    AUTHORS:

    - Robert L. Miller

    EXAMPLES::

        sage: from sage.combinat.partition import cyclic_permutations_of_partition
        sage: list(cyclic_permutations_of_partition_iterator([[1,2,3,4],[5,6,7]]))
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

    Note that repeated elements are not considered equal::

        sage: list(cyclic_permutations_of_partition_iterator([[1,2,3],[4,4,4]]))
        [[[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]],
         [[1, 2, 3], [4, 4, 4]],
         [[1, 3, 2], [4, 4, 4]]]
    """
    if len(partition) == 1:
        for i in cyclic_permutations_iterator(partition[0]):
            yield [i]
    else:
        for right in cyclic_permutations_of_partition_iterator(partition[1:]):
            for perm in cyclic_permutations_iterator(partition[0]):
                yield [perm] + right

def ferrers_diagram(pi):
    """
    Return the Ferrers diagram of pi.

    INPUT:


    -  ``pi`` - a partition, given as a list of integers.


    EXAMPLES::

        sage: print ferrers_diagram([5,5,2,1])
        *****
        *****
        **
        *
        sage: pi = partitions_list(10)[30] ## [6,1,1,1,1]
        sage: print ferrers_diagram(pi)
        ******
        *
        *
        *
        *
        sage: pi = partitions_list(10)[33] ## [6, 3, 1]
        sage: print ferrers_diagram(pi)
        ******
        ***
        *
    """
    return '\n'.join(['*'*p for p in pi])


def ordered_partitions(n,k=None):
    r"""
    An ordered partition of `n` is an ordered sum

    .. math::

                n = p_1+p_2 + \cdots + p_k


    of positive integers and is represented by the list
    `p = [p_1,p_2,\cdots ,p_k]`. If `k` is omitted
    then all ordered partitions are returned.

    ``ordered_partitions(n,k)`` returns the list of all
    (ordered) partitions of the positive integer n into sums with k
    summands.

    Do not call ``ordered_partitions`` with an n much
    larger than 15, since the list will simply become too large.

    Wraps GAP's OrderedPartitions.

    The number of ordered partitions `T_n` of
    `\{ 1, 2, ..., n \}` has the generating function is

    .. math::

             \sum_n {T_n \over n!} x^n = {1 \over 2-e^x}.



    EXAMPLES::

        sage: ordered_partitions(10,2)
        [[1, 9], [2, 8], [3, 7], [4, 6], [5, 5], [6, 4], [7, 3], [8, 2], [9, 1]]

    ::

        sage: ordered_partitions(4)
        [[1, 1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 3], [2, 1, 1], [2, 2], [3, 1], [4]]

    REFERENCES:

    - http://en.wikipedia.org/wiki/Ordered_partition_of_a_set
    """
    if k is None:
        ans=gap.eval("OrderedPartitions(%s)"%(ZZ(n)))
    else:
        ans=gap.eval("OrderedPartitions(%s,%s)"%(ZZ(n),ZZ(k)))
    return eval(ans.replace('\n',''))

def number_of_ordered_partitions(n,k=None):
    """
    Returns the size of ordered_partitions(n,k). Wraps GAP's
    NrOrderedPartitions.

    It is possible to associate with every partition of the integer n a
    conjugacy class of permutations in the symmetric group on n points
    and vice versa. Therefore p(n) = NrPartitions(n) is the number of
    conjugacy classes of the symmetric group on n points.

    EXAMPLES::

        sage: number_of_ordered_partitions(10,2)
        9
        sage: number_of_ordered_partitions(15)
        16384
    """
    if k is None:
        ans=gap.eval("NrOrderedPartitions(%s)"%(n))
    else:
        ans=gap.eval("NrOrderedPartitions(%s,%s)"%(n,k))
    return ZZ(ans)

def partitions_greatest(n,k):
    """
    Returns the list of all (unordered) "restricted" partitions of the
    integer n having parts less than or equal to the integer k.

    Wraps GAP's PartitionsGreatestLE.

    EXAMPLES::

        sage: partitions_greatest(10,2)
        [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
         [2, 1, 1, 1, 1, 1, 1, 1, 1],
         [2, 2, 1, 1, 1, 1, 1, 1],
         [2, 2, 2, 1, 1, 1, 1],
         [2, 2, 2, 2, 1, 1],
         [2, 2, 2, 2, 2]]
    """
    return eval(gap.eval("PartitionsGreatestLE(%s,%s)"%(ZZ(n),ZZ(k))))

def partitions_greatest_eq(n,k):
    """
    Returns the list of all (unordered) "restricted" partitions of the
    integer n having at least one part equal to the integer k.

    Wraps GAP's PartitionsGreatestEQ.

    EXAMPLES::

        sage: partitions_greatest_eq(10,2)
        [[2, 1, 1, 1, 1, 1, 1, 1, 1],
         [2, 2, 1, 1, 1, 1, 1, 1],
         [2, 2, 2, 1, 1, 1, 1],
         [2, 2, 2, 2, 1, 1],
         [2, 2, 2, 2, 2]]
    """
    ans = gap.eval("PartitionsGreatestEQ(%s,%s)"%(n,k))
    return eval(ans)

def partitions_restricted(n,S,k=None):
    r"""
    This function will be deprecated in a future version of Sage and
    eventually removed. Use RestrictedPartitions(n, S, k).list()
    instead.

    Original docstring follows.

    A restricted partition is, like an ordinary partition, an unordered
    sum `n = p_1+p_2+\ldots+p_k` of positive integers and is
    represented by the list `p = [p_1,p_2,\ldots,p_k]`, in
    nonincreasing order. The difference is that here the `p_i`
    must be elements from the set `S`, while for ordinary
    partitions they may be elements from `[1..n]`.

    Returns the list of all restricted partitions of the positive
    integer n into sums with k summands with the summands of the
    partition coming from the set `S`. If k is not given all
    restricted partitions for all k are returned.

    Wraps GAP's RestrictedPartitions.

    EXAMPLES::

        sage: partitions_restricted(8,[1,3,5,7])
        [[1, 1, 1, 1, 1, 1, 1, 1],
         [3, 1, 1, 1, 1, 1],
         [3, 3, 1, 1],
         [5, 1, 1, 1],
         [5, 3],
         [7, 1]]
        sage: partitions_restricted(8,[1,3,5,7],2)
        [[5, 3], [7, 1]]
    """
    if k is None:
        ans=gap.eval("RestrictedPartitions(%s,%s)"%(n,S))
    else:
        ans=gap.eval("RestrictedPartitions(%s,%s,%s)"%(n,S,k))
    return eval(ans)

def number_of_partitions_restricted(n,S,k=None):
    """
    This function will be deprecated in a future version of Sage and
    eventually removed. Use RestrictedPartitions(n, S, k).count()
    instead.

    Original docstring follows.

    Returns the size of partitions_restricted(n,S,k). Wraps GAP's
    NrRestrictedPartitions.

    EXAMPLES::

        sage: number_of_partitions_restricted(8,[1,3,5,7])
        6
        sage: number_of_partitions_restricted(8,[1,3,5,7],2)
        2
    """
    if k is None:
        ans=gap.eval("NrRestrictedPartitions(%s,%s)"%(ZZ(n),S))
    else:
        ans=gap.eval("NrRestrictedPartitions(%s,%s,%s)"%(ZZ(n),S,ZZ(k)))
    return ZZ(ans)

def partitions_tuples(n,k):
    """
    partition_tuples( n, k ) returns the list of all k-tuples of
    partitions which together form a partition of n.

    k-tuples of partitions describe the classes and the characters of
    wreath products of groups with k conjugacy classes with the
    symmetric group `S_n`.

    Wraps GAP's PartitionTuples.

    EXAMPLES::

        sage: partitions_tuples(3,2)
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
    ans=gap.eval("PartitionTuples(%s,%s)"%(ZZ(n),ZZ(k)))
    return eval(ans)

def number_of_partitions_tuples(n,k):
    r"""
    number_of_partition_tuples( n, k ) returns the number of
    partition_tuples(n,k).

    Wraps GAP's NrPartitionTuples.

    EXAMPLES::

        sage: number_of_partitions_tuples(3,2)
        10
        sage: number_of_partitions_tuples(8,2)
        185

    Now we compare that with the result of the following GAP
    computation::

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
    """
    ans=gap.eval("NrPartitionTuples(%s,%s)"%(ZZ(n),ZZ(k)))
    return ZZ(ans)

def partition_power(pi,k):
    """
    partition_power( pi, k ) returns the partition corresponding to
    the `k`-th power of a permutation with cycle structure pi
    (thus describes the powermap of symmetric groups).

    Wraps GAP's PowerPartition.

    EXAMPLES::

        sage: partition_power([5,3],1)
        [5, 3]
        sage: partition_power([5,3],2)
        [5, 3]
        sage: partition_power([5,3],3)
        [5, 1, 1, 1]
        sage: partition_power([5,3],4)
        [5, 3]

    Now let us compare this to the power map on `S_8`::

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
    ans=gap.eval("PowerPartition(%s,%s)"%(pi,ZZ(k)))
    return eval(ans)

def partition_sign(pi):
    r"""
    partition_sign( pi ) returns the sign of a permutation with cycle
    structure given by the partition pi.

    This function corresponds to a homomorphism from the symmetric
    group `S_n` into the cyclic group of order 2, whose kernel
    is exactly the alternating group `A_n`. Partitions of sign
    `1` are called even partitions while partitions of sign
    `-1` are called odd.

    Wraps GAP's SignPartition.

    EXAMPLES::

        sage: partition_sign([5,3])
        1
        sage: partition_sign([5,2])
        -1

    Zolotarev's lemma states that the Legendre symbol
    `\left(\frac{a}{p}\right)` for an integer
    `a \pmod p` (`p` a prime number), can be computed
    as sign(p_a), where sign denotes the sign of a permutation and
    p_a the permutation of the residue classes `\pmod p`
    induced by modular multiplication by `a`, provided
    `p` does not divide `a`.

    We verify this in some examples.

    ::

        sage: F = GF(11)
        sage: a = F.multiplicative_generator();a
        2
        sage: plist = [int(a*F(x)) for x in range(1,11)]; plist
        [2, 4, 6, 8, 10, 1, 3, 5, 7, 9]

    This corresponds ot the permutation (1, 2, 4, 8, 5, 10, 9, 7, 3, 6)
    (acting the set `\{1,2,...,10\}`) and to the partition
    [10].

    ::

        sage: p = PermutationGroupElement('(1, 2, 4, 8, 5, 10, 9, 7, 3, 6)')
        sage: p.sign()
        -1
        sage: partition_sign([10])
        -1
        sage: kronecker_symbol(11,2)
        -1

    Now replace `2` by `3`::

        sage: plist = [int(F(3*x)) for x in range(1,11)]; plist
        [3, 6, 9, 1, 4, 7, 10, 2, 5, 8]
        sage: range(1,11)
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        sage: p = PermutationGroupElement('(3,4,8,7,9)')
        sage: p.sign()
        1
        sage: kronecker_symbol(3,11)
        1
        sage: partition_sign([5,1,1,1,1,1])
        1

    In both cases, Zolotarev holds.

    REFERENCES:

    - http://en.wikipedia.org/wiki/Zolotarev's_lemma
    """
    ans=gap.eval("SignPartition(%s)"%(pi))
    return sage_eval(ans)

def partition_associated(pi):
    """
    partition_associated( pi ) returns the "associated" (also called
    "conjugate" in the literature) partition of the partition pi which
    is obtained by transposing the corresponding Ferrers diagram.

    EXAMPLES::

        sage: partition_associated([2,2])
        [2, 2]
        sage: partition_associated([6,3,1])
        [3, 2, 2, 1, 1, 1]
        sage: print ferrers_diagram([6,3,1])
        ******
        ***
        *
        sage: print ferrers_diagram([3,2,2,1,1,1])
        ***
        **
        **
        *
        *
        *
    """
    return list(Partition(pi).conjugate())

