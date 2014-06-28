r"""
Partition/Diagram Algebras
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

from combinat import catalan_number
from combinatorial_algebra import CombinatorialAlgebra, CombinatorialAlgebraElement
from sage.combinat.set_partition import SetPartition, SetPartitions, SetPartitions_set
from sage.sets.set import Set, is_Set
from sage.graphs.graph import Graph
from sage.rings.arith import factorial, binomial
from permutation import Permutations
from sage.rings.all import Integer
from sage.rings.real_mpfr import is_RealNumber
from subset import Subsets
from sage.functions.all import ceil
import functools
import math


def create_set_partition_function(letter, k):
    """
    EXAMPLES::

        sage: from sage.combinat.partition_algebra import create_set_partition_function
        sage: create_set_partition_function('A', 3)
        Set partitions of {1, ..., 3, -1, ..., -3}
    """
    from sage.functions.all import floor
    if isinstance(k, (int, Integer)):
        if k > 0:
            return globals()['SetPartitions' + letter + 'k_k'](k)
    elif is_RealNumber(k):
        if k - math.floor(k) == 0.5:
            return globals()['SetPartitions' + letter + 'khalf_k'](floor(k))

    raise ValueError("k must be an integer or an integer + 1/2")

class SetPartitionsXkElement(SetPartition):
    """
    An element for the classes of ``SetPartitionXk`` where ``X`` is some
    letter.
    """
    def check(self):
        """
        Check to make sure this is a set partition.

        EXAMLPLES::

            sage: A2p5 = SetPartitionsAk(2.5)
            sage: A2p5.first() # random
            {{1, 2, 3, -1, -3, -2}}
        """
        #Check to make sure each element of x is a set
        for s in self:
            assert isinstance(s, (set, frozenset)) or is_Set(s)

#####
#A_k#
#####
SetPartitionsAk = functools.partial(create_set_partition_function,"A")
SetPartitionsAk.__doc__ = (
    """
    Returns the combinatorial class of set partitions of type A_k.

    EXAMPLES::

        sage: A3 = SetPartitionsAk(3); A3
        Set partitions of {1, ..., 3, -1, ..., -3}

        sage: A3.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: A3.last() #random
        {{-1}, {-2}, {3}, {1}, {-3}, {2}}
        sage: A3.random_element()  #random
        {{1, 3, -3, -1}, {2, -2}}

        sage: A3.cardinality()
        203

        sage: A2p5 = SetPartitionsAk(2.5); A2p5
        Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block
        sage: A2p5.cardinality()
        52

        sage: A2p5.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: A2p5.last() #random
        {{-1}, {-2}, {2}, {3, -3}, {1}}
        sage: A2p5.random_element() #random
        {{-1}, {-2}, {3, -3}, {1, 2}}
    """)

class SetPartitionsAk_k(SetPartitions_set):
    def __init__(self, k):
        """
        TESTS::

            sage: A3 = SetPartitionsAk(3); A3
            Set partitions of {1, ..., 3, -1, ..., -3}
            sage: A3 == loads(dumps(A3))
            True
        """
        self.k = k
        SetPartitions_set.__init__(self, frozenset(range(1,k+1) + map(lambda x: -1*x,range(1,k+1))))

    Element = SetPartitionsXkElement

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsAk(3)
            Set partitions of {1, ..., 3, -1, ..., -3}
        """
        return "Set partitions of {1, ..., %s, -1, ..., -%s}"%(self.k, self.k)

class SetPartitionsAkhalf_k(SetPartitions_set):
    def __init__(self, k):
        """
        TESTS::

            sage: A2p5 = SetPartitionsAk(2.5); A2p5
            Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block
            sage: A2p5 == loads(dumps(A2p5))
            True
        """
        self.k = k
        SetPartitions_set.__init__( self, frozenset(range(1,k+2) + map(lambda x: -1*x, range(1,k+1))) )

    Element = SetPartitionsXkElement

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsAk(2.5)
            Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block
        """
        s = self.k+1
        return "Set partitions of {1, ..., %s, -1, ..., -%s} with %s and -%s in the same block"%(s,s,s,s)

    def __contains__(self, x):
        """
        TESTS::

            sage: A2p5 = SetPartitionsAk(2.5)
            sage: all([ sp in A2p5 for sp in A2p5])
            True
            sage: A3 = SetPartitionsAk(3)
            sage: len(filter(lambda x: x in A2p5, A3))
            52
            sage: A2p5.cardinality()
            52
        """
        if x not in SetPartitionsAk_k(self.k+1):
            return False

        for part in x:
            if self.k+1 in part and -self.k-1 not in part:
                return False

        return True

    def __iter__(self):
        """
        TESTS::

            sage: SetPartitionsAk(1.5).list() #random
            [{{1, 2, -2, -1}},
             {{2, -2, -1}, {1}},
             {{2, -2}, {1, -1}},
             {{-1}, {1, 2, -2}},
             {{-1}, {2, -2}, {1}}]

        ::

            sage: ks = [ 1.5, 2.5, 3.5 ]
            sage: aks = map(SetPartitionsAk, ks)
            sage: all([ak.cardinality() == len(ak.list()) for ak in aks])
            True
        """
        kp = Set([-self.k-1])
        for sp in SetPartitions_set.__iter__(self):
            res = []
            for part in sp:
                if self.k+1 in part:
                    res.append( part + kp )
                else:
                    res.append(part)
            yield self.element_class(self, res)

#####
#S_k#
#####
SetPartitionsSk = functools.partial(create_set_partition_function,"S")
SetPartitionsSk.__doc__ = (
    """
    Returns the combinatorial class of set partitions of type S_k.  There
    is a bijection between these set partitions and the permutations
    of 1, ..., k.

    EXAMPLES::

        sage: S3 = SetPartitionsSk(3); S3
        Set partitions of {1, ..., 3, -1, ..., -3} with propagating number 3
        sage: S3.cardinality()
        6

        sage: S3.list()  #random
        [{{2, -2}, {3, -3}, {1, -1}},
         {{1, -1}, {2, -3}, {3, -2}},
         {{2, -1}, {3, -3}, {1, -2}},
         {{1, -2}, {2, -3}, {3, -1}},
         {{1, -3}, {2, -1}, {3, -2}},
         {{1, -3}, {2, -2}, {3, -1}}]
        sage: S3.first() #random
        {{2, -2}, {3, -3}, {1, -1}}
        sage: S3.last() #random
        {{1, -3}, {2, -2}, {3, -1}}
        sage: S3.random_element() #random
        {{1, -3}, {2, -1}, {3, -2}}

        sage: S3p5 = SetPartitionsSk(3.5); S3p5
        Set partitions of {1, ..., 4, -1, ..., -4} with 4 and -4 in the same block and propagating number 4
        sage: S3p5.cardinality()
        6

        sage: S3p5.list() #random
        [{{2, -2}, {3, -3}, {1, -1}, {4, -4}},
         {{2, -3}, {1, -1}, {4, -4}, {3, -2}},
         {{2, -1}, {3, -3}, {1, -2}, {4, -4}},
         {{2, -3}, {1, -2}, {4, -4}, {3, -1}},
         {{1, -3}, {2, -1}, {4, -4}, {3, -2}},
         {{1, -3}, {2, -2}, {4, -4}, {3, -1}}]
        sage: S3p5.first() #random
        {{2, -2}, {3, -3}, {1, -1}, {4, -4}}
        sage: S3p5.last() #random
        {{1, -3}, {2, -2}, {4, -4}, {3, -1}}
        sage: S3p5.random_element() #random
        {{1, -3}, {2, -2}, {4, -4}, {3, -1}}
    """)
class SetPartitionsSk_k(SetPartitionsAk_k):
    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsSk(3)
            Set partitions of {1, ..., 3, -1, ..., -3} with propagating number 3
        """
        return SetPartitionsAk_k._repr_(self) + " with propagating number %s"%self.k

    def __contains__(self, x):
        """
        TESTS::

            sage: A3 = SetPartitionsAk(3)
            sage: S3 = SetPartitionsSk(3)
            sage: all([ sp in S3 for sp in S3])
            True
            sage: S3.cardinality()
            6
            sage: len(filter(lambda x: x in S3, A3))
            6
        """
        if not SetPartitionsAk_k.__contains__(self, x):
            return False

        if propagating_number(x) != self.k:
            return False

        return True

    def cardinality(self):
        """
        Returns k!.

        TESTS::

            sage: SetPartitionsSk(2).cardinality()
            2
            sage: SetPartitionsSk(3).cardinality()
            6
            sage: SetPartitionsSk(4).cardinality()
            24
            sage: SetPartitionsSk(5).cardinality()
            120
        """
        return factorial(self.k)

    def __iter__(self):
        """
        TESTS::

            sage: SetPartitionsSk(3).list() #random
            [{{2, -2}, {3, -3}, {1, -1}},
             {{1, -1}, {2, -3}, {3, -2}},
             {{2, -1}, {3, -3}, {1, -2}},
             {{1, -2}, {2, -3}, {3, -1}},
             {{1, -3}, {2, -1}, {3, -2}},
             {{1, -3}, {2, -2}, {3, -1}}]
            sage: ks = range(1, 6)
            sage: sks = map(SetPartitionsSk, ks)
            sage: all([ sk.cardinality() == len(sk.list()) for sk in sks])
            True
        """
        for p in Permutations(self.k):
            res = []
            for i in range(self.k):
                res.append( Set([ i+1, -p[i] ]) )
            yield self.element_class(self, res)

class SetPartitionsSkhalf_k(SetPartitionsAkhalf_k):
    def __contains__(self, x):
        """
        TESTS::

            sage: S2p5 = SetPartitionsSk(2.5)
            sage: A3 = SetPartitionsAk(3)
            sage: all([sp in S2p5 for sp in S2p5])
            True
            sage: len(filter(lambda x: x in S2p5, A3))
            2
            sage: S2p5.cardinality()
            2
        """
        if not SetPartitionsAkhalf_k.__contains__(self, x):
            return False
        if propagating_number(x) != self.k+1:
            return False
        return True

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsSk(2.5)
            Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and propagating number 3
        """
        s = self.k+1
        return SetPartitionsAkhalf_k._repr_(self) + " and propagating number %s"%s

    def cardinality(self):
        """
        TESTS::

            sage: SetPartitionsSk(2.5).cardinality()
            2
            sage: SetPartitionsSk(3.5).cardinality()
            6
            sage: SetPartitionsSk(4.5).cardinality()
            24

        ::

            sage: ks = [2.5, 3.5, 4.5, 5.5]
            sage: sks = [SetPartitionsSk(k) for k in ks]
            sage: all([ sk.cardinality() == len(sk.list()) for sk in sks])
            True
        """
        return factorial(self.k)

    def __iter__(self):
        """
        TESTS::

            sage: SetPartitionsSk(3.5).list() #random indirect test
            [{{2, -2}, {3, -3}, {1, -1}, {4, -4}},
             {{2, -3}, {1, -1}, {4, -4}, {3, -2}},
             {{2, -1}, {3, -3}, {1, -2}, {4, -4}},
             {{2, -3}, {1, -2}, {4, -4}, {3, -1}},
             {{1, -3}, {2, -1}, {4, -4}, {3, -2}},
             {{1, -3}, {2, -2}, {4, -4}, {3, -1}}]
        """
        for p in Permutations(self.k):
            res = []
            for i in range(self.k):
                res.append( Set([ i+1, -p[i] ]) )

            res.append(Set([self.k+1, -self.k - 1]))
            yield self.element_class(self, res)

#####
#I_k#
#####
SetPartitionsIk = functools.partial(create_set_partition_function,"I")
SetPartitionsIk.__doc__ = (
    """
    Returns the combinatorial class of set partitions of type I_k.  These
    are set partitions with a propagating number of less than k.  Note
    that the identity set partition {{1, -1}, ..., {k, -k}} is not
    in I_k.

    EXAMPLES::

        sage: I3 = SetPartitionsIk(3); I3
        Set partitions of {1, ..., 3, -1, ..., -3} with propagating number < 3
        sage: I3.cardinality()
        197

        sage: I3.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: I3.last() #random
        {{-1}, {-2}, {3}, {1}, {-3}, {2}}
        sage: I3.random_element() #random
        {{-1}, {-3, -2}, {2, 3}, {1}}

        sage: I2p5 = SetPartitionsIk(2.5); I2p5
        Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and propagating number < 3
        sage: I2p5.cardinality()
        50

        sage: I2p5.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: I2p5.last() #random
        {{-1}, {-2}, {2}, {3, -3}, {1}}
        sage: I2p5.random_element() #random
        {{-1}, {-2}, {1, 3, -3}, {2}}

    """)
class SetPartitionsIk_k(SetPartitionsAk_k):
    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsIk(3)
            Set partitions of {1, ..., 3, -1, ..., -3} with propagating number < 3
        """
        return SetPartitionsAk_k._repr_(self) + " with propagating number < %s"%self.k

    def __contains__(self, x):
        """
        TESTS::

            sage: I3 = SetPartitionsIk(3)
            sage: A3 = SetPartitionsAk(3)
            sage: all([ sp in I3 for sp in I3])
            True
            sage: len(filter(lambda x: x in I3, A3))
            197
            sage: I3.cardinality()
            197
        """
        if not SetPartitionsAk_k.__contains__(self, x):
            return False
        if propagating_number(x) >= self.k:
            return False
        return True

    def cardinality(self):
        """
        TESTS::

            sage: SetPartitionsIk(2).cardinality()
            13
        """
        return len(self.list())

    def __iter__(self):
        """
        TESTS::

            sage: SetPartitionsIk(2).list() #random indirect test
                [{{1, 2, -1, -2}},
                 {{2, -1, -2}, {1}},
                 {{2}, {1, -1, -2}},
                 {{-1}, {1, 2, -2}},
                 {{-2}, {1, 2, -1}},
                 {{1, 2}, {-1, -2}},
                 {{2}, {-1, -2}, {1}},
                 {{-1}, {2, -2}, {1}},
                 {{-2}, {2, -1}, {1}},
                 {{-1}, {2}, {1, -2}},
                 {{-2}, {2}, {1, -1}},
                 {{-1}, {-2}, {1, 2}},
                 {{-1}, {-2}, {2}, {1}}]
        """
        for sp in SetPartitionsAk_k.__iter__(self):
            if propagating_number(sp) < self.k:
                yield sp

class SetPartitionsIkhalf_k(SetPartitionsAkhalf_k):
    def __contains__(self, x):
        """
        TESTS::

            sage: I2p5 = SetPartitionsIk(2.5)
            sage: A3 = SetPartitionsAk(3)
            sage: all([ sp in I2p5 for sp in I2p5])
            True
            sage: len(filter(lambda x: x in I2p5, A3))
            50
            sage: I2p5.cardinality()
            50
        """
        if not SetPartitionsAkhalf_k.__contains__(self, x):
            return False
        if propagating_number(x) >= self.k+1:
            return False
        return True

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsIk(2.5)
            Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and propagating number < 3
        """
        return SetPartitionsAkhalf_k._repr_(self) + " and propagating number < %s"%(self.k+1)

    def cardinality(self):
        """
        TESTS::

            sage: SetPartitionsIk(1.5).cardinality()
            4
            sage: SetPartitionsIk(2.5).cardinality()
            50
            sage: SetPartitionsIk(3.5).cardinality()
            871
        """
        return len(self.list())

    def __iter__(self):
        """
        TESTS::

            sage: SetPartitionsIk(1.5).list() #random
            [{{1, 2, -2, -1}},
             {{2, -2, -1}, {1}},
             {{-1}, {1, 2, -2}},
             {{-1}, {2, -2}, {1}}]
        """

        for sp in SetPartitionsAkhalf_k.__iter__(self):
            if propagating_number(sp) < self.k+1:
                yield sp
#####
#B_k#
#####
SetPartitionsBk = functools.partial(create_set_partition_function,"B")
SetPartitionsBk.__doc__ = (
    """
    Returns the combinatorial class of set partitions of type B_k.
    These are the set partitions where every block has size 2.

    EXAMPLES::

        sage: B3 = SetPartitionsBk(3); B3
        Set partitions of {1, ..., 3, -1, ..., -3} with block size 2

        sage: B3.first() #random
        {{2, -2}, {1, -3}, {3, -1}}
        sage: B3.last() #random
        {{1, 2}, {3, -2}, {-3, -1}}
        sage: B3.random_element() #random
        {{2, -1}, {1, -3}, {3, -2}}

        sage: B3.cardinality()
        15

        sage: B2p5 = SetPartitionsBk(2.5); B2p5
        Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and with block size 2

        sage: B2p5.first() #random
        {{2, -1}, {3, -3}, {1, -2}}
        sage: B2p5.last() #random
        {{1, 2}, {3, -3}, {-1, -2}}
        sage: B2p5.random_element() #random
        {{2, -2}, {3, -3}, {1, -1}}

        sage: B2p5.cardinality()
        3
    """)

class SetPartitionsBk_k(SetPartitionsAk_k):
    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsBk(2.5)
            Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and with block size 2
        """
        return SetPartitionsAk_k._repr_(self) + " with block size 2"

    def __contains__(self, x):
        """
        TESTS::

            sage: B3 = SetPartitionsBk(3)
            sage: A3 = SetPartitionsAk(3)
            sage: len(filter(lambda x: x in B3, A3))
            15
            sage: B3.cardinality()
            15
        """
        if not SetPartitionsAk_k.__contains__(self, x):
            return False

        for part in x:
            if len(part) != 2:
                return False

        return True

    def cardinality(self):
        """
        Returns the number of set partitions in B_k where k is an integer.
        This is given by (2k)!! = (2k-1)\*(2k-3)\*...\*5\*3\*1.

        EXAMPLES::

            sage: SetPartitionsBk(3).cardinality()
            15
            sage: SetPartitionsBk(2).cardinality()
            3
            sage: SetPartitionsBk(1).cardinality()
            1
            sage: SetPartitionsBk(4).cardinality()
            105
            sage: SetPartitionsBk(5).cardinality()
            945
        """
        c = 1
        for i in range(1, 2*self.k,2):
            c *= i
        return c

    def __iter__(self):
        """
        TESTS::

            sage: SetPartitionsBk(1).list()
            [{{-1, 1}}]

        ::

            sage: SetPartitionsBk(2).list() #random
            [{{2, -1}, {1, -2}}, {{2, -2}, {1, -1}}, {{1, 2}, {-1, -2}}]
            sage: SetPartitionsBk(3).list() #random
            [{{2, -2}, {1, -3}, {3, -1}},
             {{2, -1}, {1, -3}, {3, -2}},
             {{1, -3}, {2, 3}, {-1, -2}},
             {{3, -1}, {1, -2}, {2, -3}},
             {{3, -2}, {1, -1}, {2, -3}},
             {{1, 3}, {2, -3}, {-1, -2}},
             {{2, -1}, {3, -3}, {1, -2}},
             {{2, -2}, {3, -3}, {1, -1}},
             {{1, 2}, {3, -3}, {-1, -2}},
             {{-3, -2}, {2, 3}, {1, -1}},
             {{1, 3}, {-3, -2}, {2, -1}},
             {{1, 2}, {3, -1}, {-3, -2}},
             {{-3, -1}, {2, 3}, {1, -2}},
             {{1, 3}, {-3, -1}, {2, -2}},
             {{1, 2}, {3, -2}, {-3, -1}}]

        Check to make sure that the number of elements generated is the
        same as what is given by cardinality()

        ::

            sage: bks = [ SetPartitionsBk(i) for i in range(1, 6) ]
            sage: all( [ bk.cardinality() == len(bk.list()) for bk in bks] )
            True
        """
        for sp in SetPartitions(self._set, [2]*(len(self._set)//2)):
            yield self.element_class(self, sp)

class SetPartitionsBkhalf_k(SetPartitionsAkhalf_k):
    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsBk(2.5)
            Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and with block size 2
        """
        return SetPartitionsAkhalf_k._repr_(self) + " and with block size 2"


    def __contains__(self, x):
        """
        TESTS::

            sage: A3 = SetPartitionsAk(3)
            sage: B2p5 = SetPartitionsBk(2.5)
            sage: all([ sp in B2p5 for sp in B2p5 ])
            True
            sage: len(filter(lambda x: x in B2p5, A3))
            3
            sage: B2p5.cardinality()
            3
        """
        if not SetPartitionsAkhalf_k.__contains__(self, x):
            return False
        for part in x:
            if len(part) != 2:
                return False
        return True

    def cardinality(self):
        """
        TESTS::

            sage: B3p5 = SetPartitionsBk(3.5)
            sage: B3p5.cardinality()
            15
        """
        return len(self.list())

    def __iter__(self):
        """
        TESTS::

            sage: B3p5 = SetPartitionsBk(3.5)
            sage: B3p5.cardinality()
            15

        ::

            sage: B3p5.list() #random
            [{{2, -2}, {1, -3}, {4, -4}, {3, -1}},
             {{2, -1}, {1, -3}, {4, -4}, {3, -2}},
             {{1, -3}, {2, 3}, {4, -4}, {-1, -2}},
             {{2, -3}, {1, -2}, {4, -4}, {3, -1}},
             {{2, -3}, {1, -1}, {4, -4}, {3, -2}},
             {{1, 3}, {4, -4}, {2, -3}, {-1, -2}},
             {{2, -1}, {3, -3}, {1, -2}, {4, -4}},
             {{2, -2}, {3, -3}, {1, -1}, {4, -4}},
             {{1, 2}, {3, -3}, {4, -4}, {-1, -2}},
             {{-3, -2}, {2, 3}, {1, -1}, {4, -4}},
             {{1, 3}, {-3, -2}, {2, -1}, {4, -4}},
             {{1, 2}, {-3, -2}, {4, -4}, {3, -1}},
             {{-3, -1}, {2, 3}, {1, -2}, {4, -4}},
             {{1, 3}, {-3, -1}, {2, -2}, {4, -4}},
             {{1, 2}, {-3, -1}, {4, -4}, {3, -2}}]
        """
        set = range(1,self.k+1) + map(lambda x: -1*x, range(1,self.k+1))
        for sp in SetPartitions(set, [2]*(len(set)//2) ):
            yield self.element_class(self, Set(list(sp)) + Set([Set([self.k+1, -self.k -1])]))

#####
#P_k#
#####
SetPartitionsPk = functools.partial(create_set_partition_function,"P")
SetPartitionsPk.__doc__ = (
    """
    Returns the combinatorial class of set partitions of type P_k.
    These are the planar set partitions.

    EXAMPLES::

        sage: P3 = SetPartitionsPk(3); P3
        Set partitions of {1, ..., 3, -1, ..., -3} that are planar
        sage: P3.cardinality()
        132

        sage: P3.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: P3.last() #random
        {{-1}, {-2}, {3}, {1}, {-3}, {2}}
        sage: P3.random_element() #random
        {{1, 2, -1}, {-3}, {3, -2}}

        sage: P2p5 = SetPartitionsPk(2.5); P2p5
        Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and that are planar
        sage: P2p5.cardinality()
        42

        sage: P2p5.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: P2p5.last() #random
        {{-1}, {-2}, {2}, {3, -3}, {1}}
        sage: P2p5.random_element() #random
        {{1, 2, 3, -3}, {-1, -2}}

    """)
class SetPartitionsPk_k(SetPartitionsAk_k):
    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsPk(3)
            Set partitions of {1, ..., 3, -1, ..., -3} that are planar
        """
        return SetPartitionsAk_k._repr_(self) + " that are planar"

    def __contains__(self, x):
        """
        TESTS::

            sage: P3 = SetPartitionsPk(3)
            sage: A3 = SetPartitionsAk(3)
            sage: len(filter(lambda x: x in P3, A3))
            132
            sage: P3.cardinality()
            132
            sage: all([sp in P3 for sp in P3])
            True
        """
        if not SetPartitionsAk_k.__contains__(self, x):
            return False

        if not is_planar(x):
            return False

        return True

    def cardinality(self):
        """
        TESTS::

            sage: SetPartitionsPk(2).cardinality()
            14
            sage: SetPartitionsPk(3).cardinality()
            132
            sage: SetPartitionsPk(4).cardinality()
            1430
        """
        return catalan_number(2*self.k)

    def __iter__(self):
        """
        TESTS::

            sage: SetPartitionsPk(2).list() #random indirect test
            [{{1, 2, -1, -2}},
             {{2, -1, -2}, {1}},
             {{2}, {1, -1, -2}},
             {{-1}, {1, 2, -2}},
             {{-2}, {1, 2, -1}},
             {{2, -2}, {1, -1}},
             {{1, 2}, {-1, -2}},
             {{2}, {-1, -2}, {1}},
             {{-1}, {2, -2}, {1}},
             {{-2}, {2, -1}, {1}},
             {{-1}, {2}, {1, -2}},
             {{-2}, {2}, {1, -1}},
             {{-1}, {-2}, {1, 2}},
             {{-1}, {-2}, {2}, {1}}]
        """
        for sp in SetPartitionsAk_k.__iter__(self):
            if is_planar(sp):
                yield self.element_class(self, sp)

class SetPartitionsPkhalf_k(SetPartitionsAkhalf_k):
    def __contains__(self, x):
        """
        TESTS::

            sage: A3 = SetPartitionsAk(3)
            sage: P2p5 = SetPartitionsPk(2.5)
            sage: all([ sp in P2p5 for sp in P2p5 ])
            True
            sage: len(filter(lambda x: x in P2p5, A3))
            42
            sage: P2p5.cardinality()
            42
        """
        if not SetPartitionsAkhalf_k.__contains__(self, x):
            return False
        if not is_planar(x):
            return False

        return True

    def _repr_(self):
        """
        TESTS::

            sage: repr( SetPartitionsPk(2.5) )
            'Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and that are planar'
        """
        return SetPartitionsAkhalf_k._repr_(self) + " and that are planar"

    def cardinality(self):
        """
        TESTS::

            sage: SetPartitionsPk(2.5).cardinality()
            42
            sage: SetPartitionsPk(1.5).cardinality()
            5
        """
        return len(self.list())

    def __iter__(self):
        """
        TESTS::

            sage: SetPartitionsPk(1.5).list() #random
            [{{1, 2, -2, -1}},
             {{2, -2, -1}, {1}},
             {{2, -2}, {1, -1}},
             {{-1}, {1, 2, -2}},
             {{-1}, {2, -2}, {1}}]
        """
        for sp in SetPartitionsAkhalf_k.__iter__(self):
            if is_planar(sp):
                yield self.element_class(self, sp)


#####
#T_k#
#####
SetPartitionsTk = functools.partial(create_set_partition_function,"T")
SetPartitionsTk.__doc__ = (
    """
    Returns the combinatorial class of set partitions of type T_k.
    These are planar set partitions where every block is of size 2.

    EXAMPLES::

        sage: T3 = SetPartitionsTk(3); T3
        Set partitions of {1, ..., 3, -1, ..., -3} with block size 2 and that are planar
        sage: T3.cardinality()
        5

        sage: T3.first() #random
        {{1, -3}, {2, 3}, {-1, -2}}
        sage: T3.last() #random
        {{1, 2}, {3, -1}, {-3, -2}}
        sage: T3.random_element() #random
        {{1, -3}, {2, 3}, {-1, -2}}

        sage: T2p5 = SetPartitionsTk(2.5); T2p5
        Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and with block size 2 and that are planar
        sage: T2p5.cardinality()
        2

        sage: T2p5.first() #random
        {{2, -2}, {3, -3}, {1, -1}}
        sage: T2p5.last() #random
        {{1, 2}, {3, -3}, {-1, -2}}

    """)
class SetPartitionsTk_k(SetPartitionsBk_k):
    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsTk(3)
            Set partitions of {1, ..., 3, -1, ..., -3} with block size 2 and that are planar
        """
        return SetPartitionsBk_k._repr_(self) + " and that are planar"

    def __contains__(self, x):
        """
        TESTS::

            sage: T3 = SetPartitionsTk(3)
            sage: A3 = SetPartitionsAk(3)
            sage: all([ sp in T3 for sp in T3])
            True
            sage: len(filter(lambda x: x in T3, A3))
            5
            sage: T3.cardinality()
            5
        """
        if not SetPartitionsBk_k.__contains__(self, x):
            return False

        if not is_planar(x):
            return False

        return True

    def cardinality(self):
        """
        TESTS::

            sage: SetPartitionsTk(2).cardinality()
            2
            sage: SetPartitionsTk(3).cardinality()
            5
            sage: SetPartitionsTk(4).cardinality()
            14
            sage: SetPartitionsTk(5).cardinality()
            42
        """
        return catalan_number(self.k)

    def __iter__(self):
        """
        TESTS::

            sage: SetPartitionsTk(3).list() #random
            [{{1, -3}, {2, 3}, {-1, -2}},
             {{2, -2}, {3, -3}, {1, -1}},
             {{1, 2}, {3, -3}, {-1, -2}},
             {{-3, -2}, {2, 3}, {1, -1}},
             {{1, 2}, {3, -1}, {-3, -2}}]
        """
        for sp in SetPartitionsBk_k.__iter__(self):
            if is_planar(sp):
                yield self.element_class(self, sp)

class SetPartitionsTkhalf_k(SetPartitionsBkhalf_k):
    def __contains__(self, x):
        """
        TESTS::

            sage: A3 = SetPartitionsAk(3)
            sage: T2p5 = SetPartitionsTk(2.5)
            sage: all([ sp in T2p5 for sp in T2p5 ])
            True
            sage: len(filter(lambda x: x in T2p5, A3))
            2
            sage: T2p5.cardinality()
            2
        """
        if not SetPartitionsBkhalf_k.__contains__(self, x):
            return False
        if not is_planar(x):
            return False

        return True

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsTk(2.5)
            Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and with block size 2 and that are planar
        """
        return SetPartitionsBkhalf_k._repr_(self) + " and that are planar"

    def cardinality(self):
        """
        TESTS::

            sage: SetPartitionsTk(2.5).cardinality()
            2
            sage: SetPartitionsTk(3.5).cardinality()
            5
            sage: SetPartitionsTk(4.5).cardinality()
            14
        """
        return catalan_number(self.k)

    def __iter__(self):
        """
        TESTS::

            sage: SetPartitionsTk(3.5).list() #random
            [{{1, -3}, {2, 3}, {4, -4}, {-1, -2}},
             {{2, -2}, {3, -3}, {1, -1}, {4, -4}},
             {{1, 2}, {3, -3}, {4, -4}, {-1, -2}},
             {{-3, -2}, {2, 3}, {1, -1}, {4, -4}},
             {{1, 2}, {-3, -2}, {4, -4}, {3, -1}}]
        """
        for sp in SetPartitionsBkhalf_k.__iter__(self):
            if is_planar(sp):
                yield self.element_class(self, sp)



SetPartitionsRk = functools.partial(create_set_partition_function,"R")
SetPartitionsRk.__doc__ = (
    """
    """)
class SetPartitionsRk_k(SetPartitionsAk_k):
    def __init__(self, k):
        """
        TESTS::

            sage: R3 = SetPartitionsRk(3); R3
            Set partitions of {1, ..., 3, -1, ..., -3} with at most 1 positive and negative entry in each block
            sage: R3 == loads(dumps(R3))
            True
        """
        self.k = k
        SetPartitionsAk_k.__init__(self, k)

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsRk(3)
            Set partitions of {1, ..., 3, -1, ..., -3} with at most 1 positive and negative entry in each block
        """
        return SetPartitionsAk_k._repr_(self) + " with at most 1 positive and negative entry in each block"

    def __contains__(self, x):
        """
        TESTS::

            sage: R3 = SetPartitionsRk(3)
            sage: A3 = SetPartitionsAk(3)
            sage: all([ sp in R3 for sp in R3])
            True
            sage: len(filter(lambda x: x in R3, A3))
            34
            sage: R3.cardinality()
            34
        """
        if not SetPartitionsAk_k.__contains__(self, x):
            return False

        for block in x:
            if len(block) > 2:
                return False

            negatives = 0
            positives = 0
            for i in block:
                if i < 0:
                    negatives += 1
                else:
                    positives += 1

                if negatives > 1 or positives > 1:
                    return False

        return True

    def cardinality(self):
        """
        TESTS::

            sage: SetPartitionsRk(2).cardinality()
            7
            sage: SetPartitionsRk(3).cardinality()
            34
            sage: SetPartitionsRk(4).cardinality()
            209
            sage: SetPartitionsRk(5).cardinality()
            1546
        """
        return sum( [ binomial(self.k, l)**2*factorial(l) for l in range(self.k + 1) ] )

    def __iter__(self):
        """
        TESTS::

            sage: len(SetPartitionsRk(3).list() ) == SetPartitionsRk(3).cardinality()
            True
        """
        #The number of blocks with at most two things
        positives = Set(range(1, self.k+1))
        negatives = Set( [ -i for i in positives ] )

        yield self.element_class(self, to_set_partition([],self.k))
        for n in range(1,self.k+1):
            for top in Subsets(positives, n):
                t = list(top)
                for bottom in Subsets(negatives, n):
                    b = list(bottom)
                    for permutation in Permutations(n):
                        l = [ [t[i], b[ permutation[i] - 1 ] ] for i in range(n) ]
                        yield self.element_class(self, to_set_partition(l, k=self.k))

class SetPartitionsRkhalf_k(SetPartitionsAkhalf_k):
    def __contains__(self, x):
        """
        TESTS::

            sage: A3 = SetPartitionsAk(3)
            sage: R2p5 = SetPartitionsRk(2.5)
            sage: all([ sp in R2p5 for sp in R2p5 ])
            True
            sage: len(filter(lambda x: x in R2p5, A3))
            7
            sage: R2p5.cardinality()
            7
        """
        if not SetPartitionsAkhalf_k.__contains__(self, x):
            return False

        for block in x:
            if len(block) > 2:
                return False

            negatives = 0
            positives = 0
            for i in block:
                if i < 0:
                    negatives += 1
                else:
                    positives += 1

                if negatives > 1 or positives > 1:
                    return False


        return True

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsRk(2.5)
            Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and with at most 1 positive and negative entry in each block
        """
        return SetPartitionsAkhalf_k._repr_(self) + " and with at most 1 positive and negative entry in each block"

    def cardinality(self):
        """
        TESTS::

            sage: SetPartitionsRk(2.5).cardinality()
            7
            sage: SetPartitionsRk(3.5).cardinality()
            34
            sage: SetPartitionsRk(4.5).cardinality()
            209
        """
        return sum( [ binomial(self.k, l)**2*factorial(l) for l in range(self.k + 1) ] )

    def __iter__(self):
        """
        TESTS::

            sage: R2p5 = SetPartitionsRk(2.5)
            sage: L = list(R2p5); L #random due to sets
            [{{-2}, {-1}, {3, -3}, {2}, {1}},
             {{-2}, {3, -3}, {2}, {1, -1}},
             {{-1}, {3, -3}, {2}, {1, -2}},
             {{-2}, {2, -1}, {3, -3}, {1}},
             {{-1}, {2, -2}, {3, -3}, {1}},
             {{2, -2}, {3, -3}, {1, -1}},
             {{2, -1}, {3, -3}, {1, -2}}]
            sage: len(L)
            7
        """
        positives = Set(range(1, self.k+1))
        negatives = Set( [ -i for i in positives ] )

        yield self.element_class(self, to_set_partition([[self.k+1, -self.k-1]], self.k+1))
        for n in range(1,self.k+1):
            for top in Subsets(positives, n):
                t = list(top)
                for bottom in Subsets(negatives, n):
                    b = list(bottom)
                    for permutation in Permutations(n):
                        l = [ [t[i], b[ permutation[i] - 1 ] ] for i in range(n) ] + [ [self.k+1, -self.k-1] ]
                        yield self.element_class(self, to_set_partition(l, k=self.k+1))


SetPartitionsPRk = functools.partial(create_set_partition_function,"PR")
SetPartitionsPRk.__doc__ = (
    """
    """)
class SetPartitionsPRk_k(SetPartitionsRk_k):
    def __init__(self, k):
        """
        TESTS::

            sage: PR3 = SetPartitionsPRk(3); PR3
            Set partitions of {1, ..., 3, -1, ..., -3} with at most 1 positive and negative entry in each block and that are planar
            sage: PR3 == loads(dumps(PR3))
            True
        """
        self.k = k
        SetPartitionsRk_k.__init__(self, k)

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsPRk(3)
            Set partitions of {1, ..., 3, -1, ..., -3} with at most 1 positive and negative entry in each block and that are planar
        """
        return SetPartitionsRk_k._repr_(self) + " and that are planar"

    def __contains__(self, x):
        """
        TESTS::

            sage: PR3 = SetPartitionsPRk(3)
            sage: A3 = SetPartitionsAk(3)
            sage: all([ sp in PR3 for sp in PR3])
            True
            sage: len(filter(lambda x: x in PR3, A3))
            20
            sage: PR3.cardinality()
            20
        """
        if not SetPartitionsRk_k.__contains__(self, x):
            return False

        if not is_planar(x):
            return False

        return True

    def cardinality(self):
        """
        TESTS::

            sage: SetPartitionsPRk(2).cardinality()
            6
            sage: SetPartitionsPRk(3).cardinality()
            20
            sage: SetPartitionsPRk(4).cardinality()
            70
            sage: SetPartitionsPRk(5).cardinality()
            252
        """
        return binomial(2*self.k, self.k)

    def __iter__(self):
        """
        TESTS::

            sage: len(SetPartitionsPRk(3).list() ) == SetPartitionsPRk(3).cardinality()
            True
        """
        #The number of blocks with at most two things
        positives = Set(range(1, self.k+1))
        negatives = Set( [ -i for i in positives ] )

        yield self.element_class(self, to_set_partition([], self.k))
        for n in range(1,self.k+1):
            for top in Subsets(positives, n):
                t = sorted(top)
                for bottom in Subsets(negatives, n):
                    b = list(bottom)
                    b.sort(reverse=True)
                    l = [ [t[i], b[ i ] ] for i in range(n) ]
                    yield self.element_class(self, to_set_partition(l, k=self.k))

class SetPartitionsPRkhalf_k(SetPartitionsRkhalf_k):
    def __contains__(self, x):
        """
        TESTS::

            sage: A3 = SetPartitionsAk(3)
            sage: PR2p5 = SetPartitionsPRk(2.5)
            sage: all([ sp in PR2p5 for sp in PR2p5 ])
            True
            sage: len(filter(lambda x: x in PR2p5, A3))
            6
            sage: PR2p5.cardinality()
            6
        """
        if not SetPartitionsRkhalf_k.__contains__(self, x):
            return False

        if not is_planar(x):
            return False

        return True

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitionsPRk(2.5)
            Set partitions of {1, ..., 3, -1, ..., -3} with 3 and -3 in the same block and with at most 1 positive and negative entry in each block and that are planar
        """
        return SetPartitionsRkhalf_k._repr_(self) + " and that are planar"

    def cardinality(self):
        """
        TESTS::

            sage: SetPartitionsPRk(2.5).cardinality()
            6
            sage: SetPartitionsPRk(3.5).cardinality()
            20
            sage: SetPartitionsPRk(4.5).cardinality()
            70
        """
        return binomial(2*self.k, self.k)

    def __iter__(self):
        """
        TESTS::

            sage: L = list(SetPartitionsPRk(2.5)); L
            [{{-3, 3}, {-2}, {-1}, {1}, {2}}, {{-3, 3}, {-2}, {-1, 1}, {2}},
             {{-3, 3}, {-2, 1}, {-1}, {2}}, {{-3, 3}, {-2}, {-1, 2}, {1}},
             {{-3, 3}, {-2, 2}, {-1}, {1}}, {{-3, 3}, {-2, 2}, {-1, 1}}]
            sage: len(L)
            6
        """
        positives = Set(range(1, self.k+1))
        negatives = Set( [ -i for i in positives ] )

        yield self.element_class(self, to_set_partition([[self.k+1, -self.k-1]],k=self.k+1))
        for n in range(1,self.k+1):
            for top in Subsets(positives, n):
                t = sorted(top)
                for bottom in Subsets(negatives, n):
                    b = list(bottom)
                    b.sort(reverse=True)
                    l = [ [t[i], b[ i ] ] for i in range(n) ] + [ [self.k+1, -self.k-1] ]
                    yield self.element_class(self, to_set_partition(l, k=self.k+1))

#########################################################
#Algebras

class PartitionAlgebra_generic(CombinatorialAlgebra):
    def __init__(self, R, cclass, n, k, name=None, prefix=None):
        """
        EXAMPLES::

            sage: from sage.combinat.partition_algebra import *
            sage: s = PartitionAlgebra_sk(QQ, 3, 1)
            sage: s == loads(dumps(s))
            True
        """
        self.k = k
        self.n = n
        self._basis_keys = cclass
        self._name = "Generic partition algebra with k = %s and n = %s and basis %s"%( self.k, self.n, cclass) if name is None else name
        self._one = identity(ceil(self.k))
        self._prefix = "" if prefix is None else prefix
        CombinatorialAlgebra.__init__(self, R)

    def _multiply_basis(self, left, right):
        """
        EXAMPLES::

            sage: from sage.combinat.partition_algebra import *
            sage: s = PartitionAlgebra_sk(QQ, 3, 1)
            sage: t12 = s(Set([Set([1,-2]),Set([2,-1]),Set([3,-3])]))
            sage: t12^2 == s(1) #indirect doctest
            True
        """
        (sp, l) = set_partition_composition(left, right)
        return {sp: self.n**l}
class PartitionAlgebraElement_generic(CombinatorialAlgebraElement):
    pass

class PartitionAlgebraElement_ak(PartitionAlgebraElement_generic):
    pass
class PartitionAlgebra_ak(PartitionAlgebra_generic):
    def __init__(self, R, k, n, name=None):
        """
        EXAMPLES::

            sage: from sage.combinat.partition_algebra import *
            sage: p = PartitionAlgebra_ak(QQ, 3, 1)
            sage: p == loads(dumps(p))
            True
        """
        if name is None:
            name = "Partition algebra A_%s(%s)"%(k, n)
        cclass = SetPartitionsAk(k)
        self._element_class = PartitionAlgebraElement_ak
        PartitionAlgebra_generic.__init__(self, R, cclass, n, k, name=name, prefix="A")

class PartitionAlgebraElement_bk(PartitionAlgebraElement_generic):
    pass
class PartitionAlgebra_bk(PartitionAlgebra_generic):
    def __init__(self, R, k, n, name=None):
        """
        EXAMPLES::

            sage: from sage.combinat.partition_algebra import *
            sage: p = PartitionAlgebra_bk(QQ, 3, 1)
            sage: p == loads(dumps(p))
            True
        """
        if name is None:
            name = "Partition algebra B_%s(%s)"%(k, n)
        cclass = SetPartitionsBk(k)
        self._element_class = PartitionAlgebraElement_bk
        PartitionAlgebra_generic.__init__(self, R, cclass, n, k, name=name, prefix="B")

class PartitionAlgebraElement_sk(PartitionAlgebraElement_generic):
    pass
class PartitionAlgebra_sk(PartitionAlgebra_generic):
    def __init__(self, R, k, n, name=None):
        """
        EXAMPLES::

            sage: from sage.combinat.partition_algebra import *
            sage: p = PartitionAlgebra_sk(QQ, 3, 1)
            sage: p == loads(dumps(p))
            True
        """
        if name is None:
            name = "Partition algebra S_%s(%s)"%(k, n)
        cclass = SetPartitionsSk(k)
        self._element_class = PartitionAlgebraElement_sk
        PartitionAlgebra_generic.__init__(self, R, cclass, n, k, name=name, prefix="S")

class PartitionAlgebraElement_pk(PartitionAlgebraElement_generic):
    pass
class PartitionAlgebra_pk(PartitionAlgebra_generic):
    def __init__(self, R, k, n, name=None):
        """
        EXAMPLES::

            sage: from sage.combinat.partition_algebra import *
            sage: p = PartitionAlgebra_pk(QQ, 3, 1)
            sage: p == loads(dumps(p))
            True
        """
        if name is None:
            name = "Partition algebra P_%s(%s)"%(k, n)
        cclass = SetPartitionsPk(k)
        self._element_class = PartitionAlgebraElement_pk
        PartitionAlgebra_generic.__init__(self, R, cclass, n, k, name=name, prefix="P")

class PartitionAlgebraElement_tk(PartitionAlgebraElement_generic):
    pass
class PartitionAlgebra_tk(PartitionAlgebra_generic):
    def __init__(self, R, k, n, name=None):
        """
        EXAMPLES::

            sage: from sage.combinat.partition_algebra import *
            sage: p = PartitionAlgebra_tk(QQ, 3, 1)
            sage: p == loads(dumps(p))
            True
        """
        if name is None:
            name = "Partition algebra T_%s(%s)"%(k, n)
        cclass = SetPartitionsTk(k)
        self._element_class = PartitionAlgebraElement_tk
        PartitionAlgebra_generic.__init__(self, R, cclass, n, k, name=name, prefix="T")

class PartitionAlgebraElement_rk(PartitionAlgebraElement_generic):
    pass
class PartitionAlgebra_rk(PartitionAlgebra_generic):
    def __init__(self, R, k, n, name=None):
        """
        EXAMPLES::

            sage: from sage.combinat.partition_algebra import *
            sage: p = PartitionAlgebra_rk(QQ, 3, 1)
            sage: p == loads(dumps(p))
            True
        """
        if name is None:
            name = "Partition algebra R_%s(%s)"%(k, n)
        cclass = SetPartitionsRk(k)
        self._element_class = PartitionAlgebraElement_rk
        PartitionAlgebra_generic.__init__(self, R, cclass, n, k, name=name, prefix="R")

class PartitionAlgebraElement_prk(PartitionAlgebraElement_generic):
    pass
class PartitionAlgebra_prk(PartitionAlgebra_generic):
    def __init__(self, R, k, n, name=None):
        """
        EXAMPLES::

            sage: from sage.combinat.partition_algebra import *
            sage: p = PartitionAlgebra_prk(QQ, 3, 1)
            sage: p == loads(dumps(p))
            True
        """
        if name is None:
            name = "Partition algebra PR_%s(%s)"%(k, n)
        cclass = SetPartitionsPRk(k)
        self._element_class = PartitionAlgebraElement_prk
        PartitionAlgebra_generic.__init__(self, R, cclass, n, k, name=name, prefix="PR")

##########################################################

def is_planar(sp):
    """
    Returns True if the diagram corresponding to the set partition is
    planar; otherwise, it returns False.

    EXAMPLES::

        sage: import sage.combinat.partition_algebra as pa
        sage: pa.is_planar( pa.to_set_partition([[1,-2],[2,-1]]))
        False
        sage: pa.is_planar( pa.to_set_partition([[1,-1],[2,-2]]))
        True
    """
    to_consider = map(list, sp)

    #Singletons don't affect planarity
    to_consider = [x for x in to_consider if len(x) > 1]
    n = len(to_consider)

    for i in range(n):
        #Get the positive and negative entries of this
        #part
        ap = [x for x in to_consider[i] if x>0]
        an = [abs(x) for x in to_consider[i] if x<0]
        #print a, ap, an


        #Check if a includes numbers in both the top and bottom rows
        if len(ap) > 0 and len(an) > 0:

            for j in range(n):
                if i == j:
                    continue
                #Get the positive and negative entries of this part
                bp = [x for x in to_consider[j] if x>0]
                bn = [abs(x) for x in to_consider[j] if x<0]

                #Skip the ones that don't involve numbers in both
                #the bottom and top rows
                if len(bn) == 0 or len(bp) == 0:
                    continue

                #Make sure that if min(bp) > max(ap)
                #then min(bn) >  max(an)
                if max(bp) > max(ap):
                    if min(bn) < min(an):
                        return False


        #Go through the bottom and top rows
        for row in [ap, an]:
            if len(row) > 1:
                row.sort()
                for s in range(len(row)-1):
                    if row[s] + 1 == row[s+1]:
                        #No gap, continue on
                        continue
                    else:
                        rng = range(row[s] + 1, row[s+1])

                        #Go through and make sure any parts that
                        #contain numbers in this range are completely
                        #contained in this range
                        for j in range(n):
                            if i == j:
                                continue

                            #Make sure we make the numbers negative again
                            #if we are in the bottom row
                            if row is ap:
                                sr = Set(rng)
                            else:
                                sr = Set(map(lambda x: -1*x, rng))


                            sj = Set(to_consider[j])
                            intersection = sr.intersection(sj)
                            if intersection:
                                if sj != intersection:
                                    return False

    return True


def to_graph(sp):
    """
    Returns a graph representing the set partition sp.

    EXAMPLES::

        sage: import sage.combinat.partition_algebra as pa
        sage: g = pa.to_graph( pa.to_set_partition([[1,-2],[2,-1]])); g
        Graph on 4 vertices

    ::

        sage: g.vertices() #random
        [1, 2, -2, -1]
        sage: g.edges() #random
        [(1, -2, None), (2, -1, None)]
    """
    g = Graph()
    for part in sp:
        part_list = list(part)
        if len(part_list) > 0:
            g.add_vertex(part_list[0])
        for i in range(1, len(part_list)):
            g.add_vertex(part_list[i])
            g.add_edge(part_list[i-1], part_list[i])
    return g

def pair_to_graph(sp1, sp2):
    """
    Return a graph consisting of the disjoint union of the graphs of set
    partitions ``sp1`` and ``sp2`` along with edges joining the bottom
    row (negative numbers) of ``sp1`` to the top row (positive numbers)
    of ``sp2``.

    The vertices of the graph ``sp1`` appear in the result as pairs
    ``(k, 1)``, whereas the vertices of the graph ``sp2`` appear as
    pairs ``(k, 2)``.

    EXAMPLES::

        sage: import sage.combinat.partition_algebra as pa
        sage: sp1 = pa.to_set_partition([[1,-2],[2,-1]])
        sage: sp2 = pa.to_set_partition([[1,-2],[2,-1]])
        sage: g = pa.pair_to_graph( sp1, sp2 ); g
        Graph on 8 vertices

    ::

        sage: g.vertices() #random
        [(1, 2), (-1, 1), (-2, 2), (-1, 2), (-2, 1), (2, 1), (2, 2), (1, 1)]
        sage: g.edges() #random
        [((1, 2), (-1, 1), None),
         ((1, 2), (-2, 2), None),
         ((-1, 1), (2, 1), None),
         ((-1, 2), (2, 2), None),
         ((-2, 1), (1, 1), None),
         ((-2, 1), (2, 2), None)]

    Another example which used to be wrong until :trac:`15958`::

        sage: sp3 = pa.to_set_partition([[1, -1], [2], [-2]])
        sage: sp4 = pa.to_set_partition([[1], [-1], [2], [-2]])
        sage: g = pa.pair_to_graph( sp3, sp4 ); g
        Graph on 8 vertices

        sage: g.vertices()
        [(-2, 1), (-2, 2), (-1, 1), (-1, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
        sage: g.edges()
        [((-2, 1), (2, 2), None), ((-1, 1), (1, 1), None),
         ((-1, 1), (1, 2), None)]
    """
    g = Graph()

    #Add the first set partition to the graph
    for part in sp1:
        part_list = list(part)
        if len(part_list) > 0:
            g.add_vertex( (part_list[0],1) )

            #Add the edge to the second part of the graph
            if part_list[0] < 0:
                g.add_edge( (part_list[0], 1), (abs(part_list[0]), 2)  )

        for i in range(1, len(part_list)):
            g.add_vertex( (part_list[i], 1) )

            #Add the edge to the second part of the graph
            if part_list[i] < 0:
                g.add_edge( (part_list[i], 1), (abs(part_list[i]), 2) )

            #Add the edge between adjacent elements of a part
            g.add_edge( (part_list[i-1], 1), (part_list[i], 1) )

    #Add the second set partition to the graph
    for part in sp2:
        part_list = list(part)
        if len(part_list) > 0:
            g.add_vertex( (part_list[0], 2) )
        for i in range(1, len(part_list)):
            g.add_vertex( (part_list[i], 2) )
            g.add_edge( (part_list[i-1], 2), (part_list[i], 2) )

    return g

def propagating_number(sp):
    """
    Returns the propagating number of the set partition sp. The
    propagating number is the number of blocks with both a positive and
    negative number.

    EXAMPLES::

        sage: import sage.combinat.partition_algebra as pa
        sage: sp1 = pa.to_set_partition([[1,-2],[2,-1]])
        sage: sp2 = pa.to_set_partition([[1,2],[-2,-1]])
        sage: pa.propagating_number(sp1)
        2
        sage: pa.propagating_number(sp2)
        0
    """
    pn = 0
    for part in sp:
        if min(part) < 0  and max(part) > 0:
            pn += 1
    return pn

def to_set_partition(l,k=None):
    """
    Coverts a list of a list of numbers to a set partitions. Each list
    of numbers in the outer list specifies the numbers contained in one
    of the blocks in the set partition.

    If k is specified, then the set partition will be a set partition
    of 1, ..., k, -1, ..., -k. Otherwise, k will default to the minimum
    number needed to contain all of the specified numbers.

    EXAMPLES::

        sage: import sage.combinat.partition_algebra as pa
        sage: pa.to_set_partition([[1,-1],[2,-2]]) == pa.identity(2)
        True
    """
    if k is None:
        if l == []:
            return Set([])
        else:
            k = max( map( lambda x: max( map(abs, x) ), l) )

    to_be_added = Set( range(1, k+1) + map(lambda x: -1*x, range(1, k+1) ) )

    sp = []
    for part in l:
        spart = Set(part)
        to_be_added -= spart
        sp.append(spart)

    for singleton in to_be_added:
        sp.append(Set([singleton]))

    return Set(sp)

def identity(k):
    """
    Returns the identity set partition 1, -1, ..., k, -k

    EXAMPLES::

        sage: import sage.combinat.partition_algebra as pa
        sage: pa.identity(2)
        {{2, -2}, {1, -1}}
    """
    res = []
    for i in range(1, k+1):
        res.append(Set([i, -i]))
    return Set(res)


def set_partition_composition(sp1, sp2):
    """
    Returns a tuple consisting of the composition of the set partitions
    sp1 and sp2 and the number of components removed from the middle
    rows of the graph.

    EXAMPLES::

        sage: import sage.combinat.partition_algebra as pa
        sage: sp1 = pa.to_set_partition([[1,-2],[2,-1]])
        sage: sp2 = pa.to_set_partition([[1,-2],[2,-1]])
        sage: pa.set_partition_composition(sp1, sp2) == (pa.identity(2), 0)
        True
    """
    g = pair_to_graph(sp1, sp2)
    connected_components = g.connected_components()

    res = []
    total_removed = 0
    for cc in connected_components:
        #Remove the vertices that live in the middle two rows
        new_cc = [x for x in cc if not( (x[0]<0 and x[1] == 1) or (x[0]>0 and x[1]==2) )]

        if new_cc == []:
            if len(cc) > 1:
                total_removed += 1
        else:
            res.append( Set(map(lambda x: x[0], new_cc)) )


    return ( Set(res), total_removed )
