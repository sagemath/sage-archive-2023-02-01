r"""
Ordered Set Partitions

An ordered set partition  p of a set s is a partition of s, into subsets called parts and represented as a list of sets. By extension, an ordered set partition of a nonnegative integer n is the set partition of the integers from 1 to n. The number of ordered set partitions of n is called the n-th ordered Bell number.


There is a natural integer composition associated with an ordered set partition, that is the sequence of sizes of all its parts in order.

EXAMPLES:
  There are 13 ordered set partitions of {1,2,3}.

    sage: OrderedSetPartitions(3).count()
    13

  Here is the list of them:

    sage: OrderedSetPartitions(3).list() #random due to the sets
    [[{1}, {2}, {3}],
     [{1}, {3}, {2}],
     [{2}, {1}, {3}],
     [{3}, {1}, {2}],
     [{2}, {3}, {1}],
     [{3}, {2}, {1}],
     [{1}, {2, 3}],
     [{2}, {1, 3}],
     [{3}, {1, 2}],
     [{1, 2}, {3}],
     [{1, 3}, {2}],
     [{2, 3}, {1}],
     [{1, 2, 3}]]

  There are 12 ordered set partitions of {1,2,3,4} whose underlying
  composition is [1,2,1].

    sage: OrderedSetPartitions(4,[1,2,1]).list() #random due to the sets
    [[{1}, {2, 3}, {4}],
     [{1}, {2, 4}, {3}],
     [{1}, {3, 4}, {2}],
     [{2}, {1, 3}, {4}],
     [{2}, {1, 4}, {3}],
     [{3}, {1, 2}, {4}],
     [{4}, {1, 2}, {3}],
     [{3}, {1, 4}, {2}],
     [{4}, {1, 3}, {2}],
     [{2}, {3, 4}, {1}],
     [{3}, {2, 4}, {1}],
     [{4}, {2, 3}, {1}]]

AUTHORS:
    --Mike Hansen
    --MuPAD-Combinat developers (for algorithms and design inspiration)

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

from sage.rings.arith import factorial, binomial
import itertools
import __builtin__
import sage.combinat.composition as composition
import sage.combinat.word as word
import sage.combinat.permutation as permutation
import sage.rings.integer
from sage.combinat.combinat import stirling_number2
from sage.sets.set import Set, is_Set
from combinat import CombinatorialClass, CombinatorialObject
from sage.misc.misc import prod


def OrderedSetPartitions(s, c=None):
    """
    Returns the combinatorial class of ordered set partitions
    of s.

    EXAMPLES:
        sage: OS = OrderedSetPartitions([1,2,3,4]); OS
        Ordered set partitions of {1, 2, 3, 4}
        sage: OS.count()
        75
        sage: OS.first()
        [{1}, {2}, {3}, {4}]
        sage: OS.last()
        [{1, 2, 3, 4}]
        sage: OS.random() #random
        [{1}, {3}, {2, 4}]

        sage: OS = OrderedSetPartitions([1,2,3,4], [2,2]); OS
        Ordered set partitions of {1, 2, 3, 4} into parts of size [2, 2]
        sage: OS.count()
        6
        sage: OS.first()
        [{1, 2}, {3, 4}]
        sage: OS.last()
        [{3, 4}, {1, 2}]
        sage: OS.list()
        [[{1, 2}, {3, 4}],
         [{1, 3}, {2, 4}],
         [{1, 4}, {2, 3}],
         [{2, 3}, {1, 4}],
         [{2, 4}, {1, 3}],
         [{3, 4}, {1, 2}]]

         sage: OS = OrderedSetPartitions("cat"); OS
         Ordered set partitions of ['c', 'a', 't']
         sage: OS.list()
         [[{'a'}, {'c'}, {'t'}],
          [{'a'}, {'t'}, {'c'}],
          [{'c'}, {'a'}, {'t'}],
          [{'t'}, {'a'}, {'c'}],
          [{'c'}, {'t'}, {'a'}],
          [{'t'}, {'c'}, {'a'}],
          [{'a'}, {'c', 't'}],
          [{'c'}, {'a', 't'}],
          [{'t'}, {'a', 'c'}],
          [{'a', 'c'}, {'t'}],
          [{'a', 't'}, {'c'}],
          [{'c', 't'}, {'a'}],
          [{'a', 'c', 't'}]]
    """
    if isinstance(s, (int, sage.rings.integer.Integer)):
        if s < 0:
            raise ValueError, "s must be non-negative"
        set = Set(range(1, s+1))
    elif isinstance(s, str):
        set = [x for x in s]
    else:
        set = Set(s)

    if c is None:
        return OrderedSetPartitions_s(set)
    else:
        if isinstance(c, (int, sage.rings.integer.Integer)):
            return OrderedSetPartitions_sn(set, c)
        elif c not in composition.Compositions(len(set)):
            raise ValueError, "c must be a composition of %s"%len(set)
        else:
            return OrderedSetPartitions_scomp(set,c)


class OrderedSetPartitions_s(CombinatorialClass):
    def __init__(self, s):
        """
        TESTS:
            sage: OS = OrderedSetPartitions([1,2,3,4])
            sage: OS == loads(dumps(OS))
            True
        """
        self.s = s

    def __repr__(self):
        """
        TESTS:
            sage: repr(OrderedSetPartitions([1,2,3,4]))
            'Ordered set partitions of {1, 2, 3, 4}'
        """
        return "Ordered set partitions of %s"%self.s

    def __contains__(self, x):
        """
        TESTS:
            sage: OS = OrderedSetPartitions([1,2,3,4])
            sage: all([sp in OS for sp in OS])
            True
        """
        #x must be a list
        if not isinstance(x, list):
            return False

        #The total number of elements in the list
        #should be the same as the number is self.s
        if sum(map(len, x)) != len(self.s):
            return False

        #Check to make sure each element of the list
        #is a set
        u = Set([])
        for s in x:
            if not is_Set(s):
                return False
            u = u.union(s)

        #Make sure that the union of all the
        #sets is the original set
        if u != Set(self.s):
            return False

        return True

    def count(self):
        """
        EXAMPLES:
            sage: OrderedSetPartitions(0).count()
            1
            sage: OrderedSetPartitions(1).count()
            1
            sage: OrderedSetPartitions(2).count()
            3
            sage: OrderedSetPartitions(3).count()
            13
            sage: OrderedSetPartitions([1,2,3]).count()
            13
            sage: OrderedSetPartitions(4).count()
            75
            sage: OrderedSetPartitions(5).count()
            541
        """
        set = self.s
        return sum([factorial(k)*stirling_number2(len(set),k) for k in range(len(set)+1)])



    def iterator(self):
        """
        EXAMPLES:
            sage: OrderedSetPartitions([1,2,3]).list()
            [[{1}, {2}, {3}],
             [{1}, {3}, {2}],
             [{2}, {1}, {3}],
             [{3}, {1}, {2}],
             [{2}, {3}, {1}],
             [{3}, {2}, {1}],
             [{1}, {2, 3}],
             [{2}, {1, 3}],
             [{3}, {1, 2}],
             [{1, 2}, {3}],
             [{1, 3}, {2}],
             [{2, 3}, {1}],
             [{1, 2, 3}]]
        """
        for x in composition.Compositions(len(self.s)):
            for z in OrderedSetPartitions(self.s,x):
                yield z


class OrderedSetPartitions_sn(CombinatorialClass):
    def __init__(self, s, n):
        """
        TESTS:
            sage: OS = OrderedSetPartitions([1,2,3,4], 2)
            sage: OS == loads(dumps(OS))
            True
        """
        self.s = s
        self.n = n

    def __contains__(self, x):
        """
        TESTS:
            sage: OS = OrderedSetPartitions([1,2,3,4], 2)
            sage: all([sp in OS for sp in OS])
            True
            sage: OS.count()
            14
            sage: len(filter(lambda x: x in OS, OrderedSetPartitions([1,2,3,4])))
            14
        """
        return x in OrderedSetPartitions_s(self.s) and len(x) == self.n

    def __repr__(self):
        """
        TESTS:
            sage: repr(OrderedSetPartitions([1,2,3,4], 2))
            'Ordered set partitions of {1, 2, 3, 4} into 2 parts'
        """
        return "Ordered set partitions of %s into %s parts"%(self.s,self.n)

    def count(self):
        """

        EXAMPLES:
            sage: OrderedSetPartitions(4,2).count()
            14
            sage: OrderedSetPartitions(4,1).count()
            1

        """
        set = self.s
        n   = self.n
        return factorial(n)*stirling_number2(len(set),n)

    def iterator(self):
        """
        EXAMPLES:
            sage: OrderedSetPartitions([1,2,3,4], 2).list()
            [[{1}, {2, 3, 4}],
             [{2}, {1, 3, 4}],
             [{3}, {1, 2, 4}],
             [{4}, {1, 2, 3}],
             [{1, 2}, {3, 4}],
             [{1, 3}, {2, 4}],
             [{1, 4}, {2, 3}],
             [{2, 3}, {1, 4}],
             [{2, 4}, {1, 3}],
             [{3, 4}, {1, 2}],
             [{1, 2, 3}, {4}],
             [{1, 2, 4}, {3}],
             [{1, 3, 4}, {2}],
             [{2, 3, 4}, {1}]]
        """
        for x in composition.Compositions(len(self.s),length=self.n):
            for z in OrderedSetPartitions_scomp(self.s,x):
                yield z

class OrderedSetPartitions_scomp(CombinatorialClass):
    def __init__(self, s, comp):
        """
        TESTS:
            sage: OS = OrderedSetPartitions([1,2,3,4], [2,1,1])
            sage: OS == loads(dumps(OS))
            True
        """
        self.s = s
        self.c = composition.Composition(comp)

    def __repr__(self):
        """
        TESTS:
            sage: repr(OrderedSetPartitions([1,2,3,4], [2,1,1]))
            'Ordered set partitions of {1, 2, 3, 4} into parts of size [2, 1, 1]'
        """
        return "Ordered set partitions of %s into parts of size %s"%(self.s,self.c)

    def __contains__(self, x):
        """
        TESTS:
            sage: OS = OrderedSetPartitions([1,2,3,4], [2,1,1])
            sage: all([ sp in OS for sp in OS])
            True
            sage: OS.count()
            12
            sage: len(filter(lambda x: x in OS, OrderedSetPartitions([1,2,3,4])))
            12
        """
        return x in OrderedSetPartitions_s(self.s) and map(len, x) == self.c

    def count(self):
        """
        EXAMPLES:
            sage: OrderedSetPartitions(5,[2,3]).count()
            10
            sage: OrderedSetPartitions(0, []).count()
            1
            sage: OrderedSetPartitions(0, [0]).count()
            1
            sage: OrderedSetPartitions(0, [0,0]).count()
            1
            sage: OrderedSetPartitions(5, [2,0,3]).count()
            10
        """
        return factorial(len(self.s))/prod([factorial(i) for i in self.c])

    def iterator(self):
        """
        TESTS:
            sage: OrderedSetPartitions([1,2,3,4], [2,1,1]).list()
            [[{1, 2}, {3}, {4}],
             [{1, 2}, {4}, {3}],
             [{1, 3}, {2}, {4}],
             [{1, 4}, {2}, {3}],
             [{1, 3}, {4}, {2}],
             [{1, 4}, {3}, {2}],
             [{2, 3}, {1}, {4}],
             [{2, 4}, {1}, {3}],
             [{3, 4}, {1}, {2}],
             [{2, 3}, {4}, {1}],
             [{2, 4}, {3}, {1}],
             [{3, 4}, {2}, {1}]]
        """
        comp = self.c
        lset = [x for x in self.s]
        l = len(self.c)
        dcomp = [-1] + comp.descents(final_descent=True)

        p = []
        for j in range(l):
            p += [j+1]*comp[j]

        for x in permutation.Permutations(p):
            res = permutation.Permutation_class(range(1,len(lset)))*word.standard(x).inverse()
            res =[lset[x-1] for x in res]
            yield [ Set( res[dcomp[i]+1:dcomp[i+1]+1] ) for i in range(l)]




