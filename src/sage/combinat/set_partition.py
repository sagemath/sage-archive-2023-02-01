r"""
Set Partitions

A set partition s of a set set is a partition of set, into subsets
called parts and represented as a set of sets. By extension, a set
partition of a nonnegative integer n is the set partition of the
integers from 1 to n. The number of set partitions of n is called
the n-th Bell number.

There is a natural integer partition associated with a set
partition, that is the non-decreasing sequence of sizes of all its
parts.

There is a classical lattice associated with all set partitions of
n. The infimum of two set partitions is the set partition obtained
by intersecting all the parts of both set partitions. The supremum
is obtained by transitive closure of the relation i related to j
iff they are in the same part in at least one of the set
partitions.

EXAMPLES: There are 5 set partitions of the set 1,2,3.

::

    sage: SetPartitions(3).count()
    5

Here is the list of them

::

    sage: SetPartitions(3).list() #random due to the sets
    [{{1, 2, 3}}, {{2, 3}, {1}}, {{1, 3}, {2}}, {{1, 2}, {3}}, {{2}, {3}, {1}}]

There are 6 set partitions of 1,2,3,4 whose underlying partition is
[2, 1, 1]::

    sage: SetPartitions(4, [2,1,1]).list() #random due to the sets
    [{{3, 4}, {2}, {1}},
     {{2, 4}, {3}, {1}},
     {{4}, {2, 3}, {1}},
     {{1, 4}, {2}, {3}},
     {{1, 3}, {4}, {2}},
     {{1, 2}, {4}, {3}}]

AUTHORS:

- Mike Hansen

- MuPAD-Combinat developers (for algorithms and design inspiration).
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

from sage.sets.set import Set, EnumeratedSet, is_Set, Set_object_enumerated
import sage.combinat.partition as partition
import sage.rings.integer
import __builtin__
import itertools
from cartesian_product import CartesianProduct
import sage.combinat.subset as subset
import sage.combinat.set_partition_ordered as set_partition_ordered
import copy
from combinat import CombinatorialClass, CombinatorialObject, bell_number, stirling_number2
from subword import Subwords
from misc import IterableFunctionCall


def SetPartitions(s, part=None):
    """
    An unordered partition of a set `S` is a set of pairwise
    disjoint nonempty subsets with union `S` and is represented
    by a sorted list of such subsets.

    SetPartitions(s) returns the class of all set partitions of the set
    s, which can be a set or a string; if a string, each character is
    considered an element.

    SetPartitions(n), where n is an integer, returns the class of all
    set partitions of the set [1, 2,..., n].

    You may specify a second argument k. If k is an integer,
    SetPartitions returns the class of set partitions into k parts; if
    it is an integer partition, SetPartitions returns the class of set
    partitions whose block sizes correspond to that integer partition.

    The Bell number `B_n`, named in honor of Eric Temple Bell,
    is the number of different partitions of a set with n elements.

    EXAMPLES::

        sage: S = [1,2,3,4]
        sage: SetPartitions(S,2)
        Set partitions of [1, 2, 3, 4] with 2 parts
        sage: SetPartitions([1,2,3,4], [3,1]).list()
        [{{2, 3, 4}, {1}}, {{1, 3, 4}, {2}}, {{3}, {1, 2, 4}}, {{4}, {1, 2, 3}}]
        sage: SetPartitions(7, [3,3,1]).count()
        70

    In strings, repeated letters are considered distinct::

        sage: SetPartitions('aabcd').count()
        52
        sage: SetPartitions('abcde').count()
        52

    REFERENCES:

    - http://en.wikipedia.org/wiki/Partition_of_a_set
    """
    if isinstance(s, (int, sage.rings.integer.Integer)):
        set = range(1, s+1)
    elif isinstance(s, str):
        set = [x for x in s]
    else:
        set = s

    if part is not None:
        if isinstance(part, (int, sage.rings.integer.Integer)):
            if len(set) < part:
                raise ValueError, "part must be <= len(set)"
            else:
                return SetPartitions_setn(set,part)
        else:
            if part not in partition.Partitions(len(set)):
                raise ValueError, "part must be a partition of %s"%len(set)
            else:
                return SetPartitions_setparts(set, [partition.Partition(part)])
    else:
        return SetPartitions_set(set)



class SetPartitions_setparts(CombinatorialClass):
    object_class = Set_object_enumerated
    def __init__(self, set, parts):
        """
        TESTS::

            sage: S = SetPartitions(4, [2,2])
            sage: S == loads(dumps(S))
            True
        """
        self.set = set
        self.parts = parts

    def __repr__(self):
        """
        TESTS::

            sage: repr(SetPartitions(4, [2,2]))
            'Set partitions of [1, 2, 3, 4] with sizes in [[2, 2]]'
        """
        return "Set partitions of %s with sizes in %s"%(self.set, self.parts)

    def __contains__(self, x):
        """
        TESTS::

            sage: S = SetPartitions(4, [2,2])
            sage: all([sp in S for sp in S])
            True
        """
        #x must be a set
        if not is_Set(x):
            return False

        #Make sure that the number of elements match up
        if sum(map(len, x)) != len(self.set):
            return False

        #Check to make sure each element of x
        #is a set
        u = Set([])
        for s in x:
            if not is_Set(s):
                return False
            u = u.union(s)

        #Make sure that the union of all the
        #sets is the original set
        if u != Set(self.set):
            return False

        return True

    def count(self):
        """
        Returns the number of set partitions of set. This number is given
        by the n-th Bell number where n is the number of elements in the
        set.

        If a partition or partition length is specified, then count will
        generate all of the set partitions.

        EXAMPLES::

            sage: SetPartitions([1,2,3,4]).count()
            15
            sage: SetPartitions(3).count()
            5
            sage: SetPartitions(3,2).count()
            3
            sage: SetPartitions([]).count()
            1
        """
        return len(self.list())


    def _iterator_part(self, part):
        """
        Returns an iterator for the set partitions with block sizes
        corresponding to the partition part.

        INPUT:


        -  ``part`` - a Partition object


        EXAMPLES::

            sage: S = SetPartitions(3)
            sage: it = S._iterator_part(Partition([1,1,1]))
            sage: list(sorted(map(list, it.next())))
            [[1], [2], [3]]
            sage: S21 = SetPartitions(3,Partition([2,1]))
            sage: len(list(S._iterator_part(Partition([2,1])))) == S21.count()
            True
        """
        set = self.set

        nonzero = []
        expo = [0]+part.to_exp()

        for i in range(len(expo)):
            if expo[i] != 0:
                nonzero.append([i, expo[i]])

        taillesblocs = map(lambda x: (x[0])*(x[1]), nonzero)

        blocs = set_partition_ordered.OrderedSetPartitions(copy.copy(set), taillesblocs).list()

        for b in blocs:
            lb = [ IterableFunctionCall(_listbloc, nonzero[i][0], nonzero[i][1], b[i]) for i in range(len(nonzero)) ]
            for x in itertools.imap(lambda x: _union(x), CartesianProduct( *lb )):
                yield x



    def iterator(self):
        """
        An iterator for all the set partitions of the set.

        EXAMPLES::

            sage: SetPartitions(3).list()
            [{{1, 2, 3}}, {{2, 3}, {1}}, {{1, 3}, {2}}, {{1, 2}, {3}}, {{2}, {3}, {1}}]
        """
        for p in self.parts:
            for sp in self._iterator_part(p):
                yield sp


class SetPartitions_setn(SetPartitions_setparts):
    def __init__(self, set, n):
        """
        TESTS::

            sage: S = SetPartitions(5, 3)
            sage: S == loads(dumps(S))
            True
        """
        self.n = n
        SetPartitions_setparts.__init__(self, set, partition.Partitions(len(set), length=n).list())

    def __repr__(self):
        """
        TESTS::

            sage: repr(SetPartitions(5, 3))
            'Set partitions of [1, 2, 3, 4, 5] with 3 parts'
        """
        return "Set partitions of %s with %s parts"%(self.set,self.n)

    def count(self):
        """
        The Stirling number of the second kind is the number of partitions
        of a set of size n into k blocks.

        EXAMPLES::

            sage: SetPartitions(5, 3).count()
            25
            sage: stirling_number2(5,3)
            25
        """
        return stirling_number2(len(self.set), self.n)

class SetPartitions_set(SetPartitions_setparts):
    def __init__(self, set):
        """
        TESTS::

            sage: S = SetPartitions([1,2,3])
            sage: S == loads(dumps(S))
            True
        """
        SetPartitions_setparts.__init__(self, set, partition.Partitions(len(set)))

    def __repr__(self):
        """
        TESTS::

            sage: repr( SetPartitions([1,2,3]) )
            'Set partitions of [1, 2, 3]'
        """
        return "Set partitions of %s"%(self.set)

    def count(self):
        """
        EXAMPLES::

            sage: SetPartitions(4).count()
            15
            sage: bell_number(4)
            15
        """
        return bell_number(len(self.set))



def _listbloc(n, nbrepets, listint=None):
    """
    listbloc decomposes a set of n\*nbrepets integers (the list
    listint) in nbrepets parts.

    It is used in the algorithm to generate all set partitions.

    Not to be called by the user.

    EXAMPLES::

        sage: list(sage.combinat.set_partition._listbloc(2,1))
        [{{1, 2}}]
        sage: l = [Set([Set([3, 4]), Set([1, 2])]), Set([Set([2, 4]), Set([1, 3])]), Set([Set([2, 3]), Set([1, 4])])]
        sage: list(sage.combinat.set_partition._listbloc(2,2,[1,2,3,4])) == l
        True
    """
    if isinstance(listint, (int, sage.rings.integer.Integer)) or listint is None:
        listint = Set(range(1,n+1))


    if nbrepets == 1:
        yield Set([listint])
        return

    l = __builtin__.list(listint)
    l.sort()
    smallest = Set(l[:1])
    new_listint = Set(l[1:])

    f = lambda u, v: u.union(_set_union([smallest,v]))

    for ssens in subset.Subsets(new_listint, n-1):
        for z in _listbloc(n, nbrepets-1, new_listint-ssens):
            yield f(z,ssens)



def _union(s):
    """
    TESTS::

        sage: s = Set([ Set([1,2]), Set([3,4]) ])
        sage: sage.combinat.set_partition._union(s)
        {1, 2, 3, 4}
    """
    result = Set([])
    for ss in s:
        result = result.union(ss)
    return result

def _set_union(s):
    """
    TESTS::

        sage: s = Set([ Set([1,2]), Set([3,4]) ])
        sage: sage.combinat.set_partition._set_union(s)
        {{1, 2, 3, 4}}
    """
    result = Set([])
    for ss in s:
        result = result.union(ss)
    return EnumeratedSet([result])

def inf(s,t):
    """
    Returns the infimum of the two set partitions s and t.

    EXAMPLES::

        sage: sp1 = Set([Set([2,3,4]),Set([1])])
        sage: sp2 = Set([Set([1,3]), Set([2,4])])
        sage: s = Set([ Set([2,4]), Set([3]), Set([1])]) #{{2, 4}, {3}, {1}}
        sage: sage.combinat.set_partition.inf(sp1, sp2) == s
        True
    """
    temp = [ss.intersection(ts) for ss in s for ts in t]
    temp = filter(lambda x: x != Set([]), temp)
    return EnumeratedSet(temp)

def sup(s,t):
    """
    Returns the supremum of the two set partitions s and t.

    EXAMPLES::

        sage: sp1 = Set([Set([2,3,4]),Set([1])])
        sage: sp2 = Set([Set([1,3]), Set([2,4])])
        sage: s = Set([ Set([1,2,3,4]) ])
        sage: sage.combinat.set_partition.sup(sp1, sp2) == s
        True
    """
    res = s
    for p in t:
        inters = Set(filter(lambda x: x.intersection(p) != Set([]), __builtin__.list(res)))
        res = res.difference(inters).union(_set_union(inters))

    return res


def standard_form(sp):
    """
    Returns the set partition as a list of lists.

    EXAMPLES::

        sage: map(sage.combinat.set_partition.standard_form, SetPartitions(4, [2,2]))
        [[[3, 4], [1, 2]], [[2, 4], [1, 3]], [[2, 3], [1, 4]]]
    """

    return [__builtin__.list(x) for x in sp]


def less(s, t):
    """
    Returns True if s t otherwise it returns False.

    EXAMPLES::

        sage: z = SetPartitions(3).list()
        sage: sage.combinat.set_partition.less(z[0], z[1])
        False
        sage: sage.combinat.set_partition.less(z[4], z[1])
        True
        sage: sage.combinat.set_partition.less(z[4], z[0])
        True
        sage: sage.combinat.set_partition.less(z[3], z[0])
        True
        sage: sage.combinat.set_partition.less(z[2], z[0])
        True
        sage: sage.combinat.set_partition.less(z[1], z[0])
        True
        sage: sage.combinat.set_partition.less(z[0], z[0])
        False
    """

    if _union(s) != _union(t):
        raise ValueError, "cannont compare partitions of different sets"

    if s == t:
        return False

    for p in s:
        f = lambda z: z.intersection(p) != Set([])
        if len(filter(f, __builtin__.list(t)) ) != 1:
            return False

    return True


