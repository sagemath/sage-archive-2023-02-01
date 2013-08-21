r"""
Set Partitions

AUTHORS:

- Mike Hansen

- MuPAD-Combinat developers (for algorithms and design inspiration).

- Travis Scrimshaw (2013-02-28): Removed ``CombinatorialClass`` and added
  entry point through :class:`SetPartition`.
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

from sage.sets.set import Set, is_Set, Set_object_enumerated

import itertools
import copy

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.infinity import infinity
from sage.rings.integer import Integer

from sage.combinat.cartesian_product import CartesianProduct
from sage.combinat.misc import IterableFunctionCall
from sage.combinat.combinatorial_map import combinatorial_map
import sage.combinat.subset as subset
from sage.combinat.partition import Partition, Partitions
from sage.combinat.set_partition_ordered import OrderedSetPartitions
from sage.combinat.combinat import bell_number, stirling_number2
from sage.combinat.permutation import Permutation

class SetPartition(ClonableArray):
    """
    A partition of a set.

    A set partition `p` of a set `S` is a partition of `S` into subsets
    called parts and represented as a set of sets. By extension, a set
    partition of a nonnegative integer `n` is the set partition of the
    integers from 1 to `n`. The number of set partitions of `n` is called
    the `n`-th Bell number.

    There is a natural integer partition associated with a set partition,
    that is the non-decreasing sequence of sizes of all its parts.

    There is a classical lattice associated with all set partitions of
    `n`. The infimum of two set partitions is the set partition obtained
    by intersecting all the parts of both set partitions. The supremum
    is obtained by transitive closure of the relation `i` related to `j`
    if and only if they are in the same part in at least one of the set
    partitions.

    EXAMPLES: There are 5 set partitions of the set `\{1,2,3\}`.

    ::

        sage: SetPartitions(3).cardinality()
        5

    Here is the list of them

    ::

        sage: SetPartitions(3).list()
        [{{1, 2, 3}},
         {{1}, {2, 3}},
         {{1, 3}, {2}},
         {{1, 2}, {3}},
         {{1}, {2}, {3}}]

    There are 6 set partitions of `\{1,2,3,4\}` whose underlying partition is
    `[2, 1, 1]`::

        sage: SetPartitions(4, [2,1,1]).list()
        [{{1}, {2}, {3, 4}},
         {{1}, {2, 4}, {3}},
         {{1}, {2, 3}, {4}},
         {{1, 4}, {2}, {3}},
         {{1, 3}, {2}, {4}},
         {{1, 2}, {3}, {4}}]

    Since :trac:`14140`, we can create a set partition directly by
    :class:`SetPartition` which creates the parent object by taking the
    union of the partitions passed in. However it is recommended and
    (marginally) faster to create the parent first and then create the set
    partition from that. ::

        sage: s = SetPartition([[1,3],[2,4]]); s
        {{1, 3}, {2, 4}}
        sage: s.parent()
        Set partitions of {1, 2, 3, 4}
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, parts):
        """
        Create a set partition from ``parts`` with the appropriate parent.

        EXAMPLES::

            sage: s = SetPartition([[1,3],[2,4]]); s
            {{1, 3}, {2, 4}}
            sage: s.parent()
            Set partitions of {1, 2, 3, 4}
        """
        P = SetPartitions( _union(Set(map(Set, parts))) )
        return P.element_class(P, parts)

    def __init__(self, parent, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: s = S([[1,3],[2,4]])
            sage: TestSuite(s).run()
            sage: SetPartition([])
            {}
        """
        ClonableArray.__init__(self, parent, sorted(map(Set, s), key=min))

    def check(self):
        """
        Check that we are a valid ordered set partition.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: s = OS([[1, 3], [2, 4]])
            sage: s.check()
        """
        assert self in self.parent()

    def __mul__(self, other):
        r"""
        The product of two set partitions (resp. compositions) is defined
        as the intersection (resp. ordered) of sets appearing in the
        partitions (resp. compositions).

        EXAMPLES::

            sage: x = SetPartition([ [1,2], [3,5,4] ])
            sage: y = SetPartition(( (3,1,2), (5,4) ))
            sage: x * y
            {{1, 2}, {3}, {4, 5}}
        """
        new_composition = []
        for B in self:
           for C in other:
                BintC = B.intersection(C)
                if BintC:
                    new_composition.append(BintC)
        return self.__class__(self.parent(), new_composition)

    def _repr_(self):
        """
        Return string representation of ``self``.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: S([[1,3],[2,4]])
            {{1, 3}, {2, 4}}
        """
        return '{' + repr(list(self))[1:-1] + '}'

    def _latex_(self):
        r"""
        Return a `\LaTeX` string representation.

        EXAMPLES::

            sage: x = SetPartition([[1,2], [3,5,4]])
            sage: latex(x)
            \{\{1, 2\}, \{3, 4, 5\}\}
        """
        return repr(self).replace("{",r"\{").replace("}",r"\}")

    cardinality = ClonableArray.__len__

    def inf(self, t):
        """
        Return the infimum of ``self`` and ``t`` in the classical set partition
        lattice.

        The infimum of two set partitions is the set partition obtained
        by intersecting all the parts of both set partitions.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: sp1 = S([[2,3,4], [1]])
            sage: sp2 = S([[1,3], [2,4]])
            sage: s = S([[2,4], [3], [1]])
            sage: sp1.inf(sp2) == s
            True
        """
        temp = [ss.intersection(ts) for ss in self for ts in t]
        temp = filter(lambda x: x != Set([]), temp)
        return self.__class__(self.parent(), temp)

    def sup(self,t):
        """
        Return the supremum of ``self`` and ``t`` in the classical set
        partition lattice.

        The supremum is obtained by transitive closure of the relation `i`
        related to `j` if and only if they are in the same part in at least
        one of the set partitions.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: sp1 = S([[2,3,4], [1]])
            sage: sp2 = S([[1,3], [2,4]])
            sage: s = S([[1,2,3,4]])
            sage: sp1.sup(sp2) == s
            True
        """
        res = Set(list(self))
        for p in t:
            inters = Set(filter(lambda x: x.intersection(p) != Set([]), list(res)))
            res = res.difference(inters).union(_set_union(inters))
        return self.__class__(self.parent(), res)

    @combinatorial_map(name='shape')
    def shape(self):
        r"""
        Return the integer partition whose parts are the sizes of the sets
        in ``self``.

        EXAMPLES::

            sage: S = SetPartitions(5)
            sage: x = S([[1,2], [3,5,4]])
            sage: x.shape()
            [3, 2]
            sage: y = S([[2], [3,1], [5,4]])
            sage: y.shape()
            [2, 2, 1]
        """
        from sage.combinat.partition import Partition
        return Partition(sorted(map(len, self), reverse=True))

    # we define aliases for shape()
    shape_partition = shape
    to_partition = shape

    @combinatorial_map(name='to permutation')
    def to_permutation(self):
        """
        Convert ``self`` to a permutation by considering the partitions as
        cycles.

        EXAMPLES::

            sage: s = SetPartition([[1,3],[2,4]])
            sage: s.to_permutation()
            [3, 4, 1, 2]
        """
        return Permutation(tuple( map(tuple, self.standard_form()) ))

    def standard_form(self):
        """
        Return ``self`` as a list of lists.

        EXAMPLES::

            sage: [x.standard_form() for x in SetPartitions(4, [2,2])]
            [[[1, 2], [3, 4]], [[1, 3], [2, 4]], [[1, 4], [2, 3]]]
        """
        return list(map(list, self))

    def apply_permutation(self, p):
        r"""
        Apply ``p`` to the underlying set of ``self``.

        INPUT:

        - ``p`` -- A permutation

        EXAMPLES::

            sage: x = SetPartition([[1,2], [3,5,4]])
            sage: p = Permutation([2,1,4,5,3])
            sage: x.apply_permutation(p)
            {{1, 2}, {3, 4, 5}}
            sage: q = Permutation([3,2,1,5,4])
            sage: x.apply_permutation(q)
            {{1, 4, 5}, {2, 3}}
        """
        return self.__class__(self.parent(), [Set(map(p, B)) for B in self])

    def is_noncrossing(self):
        r"""
        Check if ``self`` is non crossing.

        EXAMPLES::

            sage: x = SetPartition([[1,2],[3,4]])
            sage: x.is_noncrossing()
            True
            sage: x = SetPartition([[1,3],[2,4]])
            sage: x.is_noncrossing()
            False

        AUTHOR: Florent Hivert
        """
        l = list(self)
        mins = map(min, l)
        maxs = map(max, l)

        for i in range(1, len(l)):
            for j in range(i):
                poss = [mins[i], maxs[i], mins[j], maxs[j]]
                possort = sorted(poss)
                cont = [possort.index(mins[i]), possort.index(maxs[i])]
                if cont == [0,2] or cont == [1,3]:
                    return False
                if (cont == [0,3] and
                    any(mins[j] < x < maxs[j] for x in l[i])):
                    return False
                if (cont == [1,2] and
                    any(mins[i] < x < maxs[i] for x in l[j])):
                    return False
        return True

class SetPartitions(Parent, UniqueRepresentation):
    r"""
    An unordered partition of a set `S` is a set of pairwise
    disjoint nonempty subsets with union `S` and is represented
    by a sorted list of such subsets.

    ``SetPartitions(s)`` returns the class of all set partitions of the set
    ``s``, which can be a set or a string; if a string, each character is
    considered an element.

    ``SetPartitions(n)``, where n is an integer, returns the class of all
    set partitions of the set `\{1, 2, \ldots, n\}`.

    You may specify a second argument `k`. If `k` is an integer,
    :class:`SetPartitions` returns the class of set partitions into `k` parts;
    if it is an integer partition, :class:`SetPartitions` returns the class of
    set partitions whose block sizes correspond to that integer partition.

    The Bell number `B_n`, named in honor of Eric Temple Bell,
    is the number of different partitions of a set with `n` elements.

    EXAMPLES::

        sage: S = [1,2,3,4]
        sage: SetPartitions(S,2)
        Set partitions of {1, 2, 3, 4} with 2 parts
        sage: SetPartitions([1,2,3,4], [3,1]).list()
        [{{1}, {2, 3, 4}}, {{1, 3, 4}, {2}}, {{1, 2, 4}, {3}}, {{1, 2, 3}, {4}}]
        sage: SetPartitions(7, [3,3,1]).cardinality()
        70

    In strings, repeated letters are not considered distinct as of
    :trac:`14140`::

        sage: SetPartitions('abcde').cardinality()
        52
        sage: SetPartitions('aabcd').cardinality()
        15

    REFERENCES:

    - :wikipedia:`Partition_of_a_set`
    """
    @staticmethod
    def __classcall_private__(cls, s, part=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: T = SetPartitions([1,2,3,4])
            sage: S is T
            True
        """
        if isinstance(s, (int, Integer)):
            s = frozenset(range(1, s+1))
        else:
            try:
                if s.cardinality() == infinity:
                    raise ValueError("The set must be finite")
            except AttributeError:
                pass
            s = frozenset(s)

        if part is not None:
            if isinstance(part, (int, Integer)):
                if len(s) < part:
                    raise ValueError("part must be <= len(set)")
                else:
                    return SetPartitions_setn(s, part)
            else:
                if part not in Partitions(len(s)):
                    raise ValueError("part must be a partition of %s"%len(s))
                else:
                    return SetPartitions_setparts(s, Partition(part))
        else:
            return SetPartitions_set(s)

    def __init__(self, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SetPartitions(3)
            sage: TestSuite(S).run()
            sage: SetPartitions(0).list()
            [{}]
            sage: SetPartitions([]).list()
            [{}]
        """
        self._set = s
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __contains__(self, x):
        """
        TESTS::

            sage: S = SetPartitions(4, [2,2])
            sage: all([sp in S for sp in S])
            True
        """
        #x must be a set, if it's not, make it into one
        if not (isinstance(x, (SetPartition, set, frozenset)) or is_Set(x)):
            return False

        #Make sure that the number of elements match up
        if sum(map(len, x)) != len(self._set):
            return False

        #Check to make sure each element of x is a set
        u = Set([])
        for s in x:
            if not (isinstance(s, (set, frozenset)) or is_Set(s)):
                return False
            u = u.union(Set(s))

        #Make sure that the union of all the
        #sets is the original set
        if u != Set(self._set):
            return False

        return True

    def _element_constructor_(self, s):
        """
        Construct an element of ``self`` from ``s``.

        INPUT:

        - ``s`` -- A set of sets

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: S([[1,3],[2,4]])
            {{1, 3}, {2, 4}}
            sage: S = SetPartitions([])
            sage: S([])
            {}
        """
        if isinstance(s, SetPartition):
            if s.parent() == self:
                return s
            raise ValueError("Cannot convert %s into an element of %s"%(s, self))
        return self.element_class(self, s)

    Element = SetPartition

    def _iterator_part(self, part):
        """
        Return an iterator for the set partitions with block sizes
        corresponding to the partition part.

        INPUT:

        -  ``part`` -- a :class:`Partition` object

        EXAMPLES::

            sage: S = SetPartitions(3)
            sage: it = S._iterator_part(Partition([1,1,1]))
            sage: list(sorted(map(list, it.next())))
            [[1], [2], [3]]
            sage: S21 = SetPartitions(3,Partition([2,1]))
            sage: len(list(S._iterator_part(Partition([2,1])))) == S21.cardinality()
            True
        """
        nonzero = []
        expo = [0]+part.to_exp()

        for i in range(len(expo)):
            if expo[i] != 0:
                nonzero.append([i, expo[i]])

        taillesblocs = map(lambda x: (x[0])*(x[1]), nonzero)

        blocs = OrderedSetPartitions(self._set, taillesblocs)

        for b in blocs:
            lb = [IterableFunctionCall(_listbloc, nonzero[i][0], nonzero[i][1], b[i]) for i in range(len(nonzero))]
            for x in itertools.imap(lambda x: _union(x), CartesianProduct( *lb )):
                yield x

    def is_less_than(self, s, t):
        """
        Check if `s < t` in the ordering on set partitions.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: s = S([[1,3],[2,4]])
            sage: t = S([[1],[2],[3],[4]])
            sage: S.is_less_than(t, s)
            True
        """
        if s.parent()._set != t.parent()._set:
            raise ValueError("cannot compare partitions of different sets")

        if s == t:
            return False

        for p in s:
            f = lambda z: z.intersection(p) != Set([])
            if len(filter(f, list(t)) ) != 1:
                return False
        return True

    lt = is_less_than

class SetPartitions_setparts(SetPartitions):
    r"""
    Class of all set partitions with fixed partition sizes corresponding to
    and integer partition `\lambda`.
    """
    @staticmethod
    def __classcall_private__(cls, s, parts):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = SetPartitions(4, [2,2])
            sage: T = SetPartitions([1,2,3,4], Partition([2,2]))
            sage: S is T
            True
        """
        if isinstance(s, (int, Integer)):
            s = xrange(1, s+1)
        return super(SetPartitions_setparts, cls).__classcall__(cls, frozenset(s), Partition(parts))

    def __init__(self, s, parts):
        """
        TESTS::

            sage: S = SetPartitions(4, [2,2])
            sage: S == loads(dumps(S))
            True
        """
        SetPartitions.__init__(self, s)
        self.parts = parts

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitions(4, [2,2])
            Set partitions of {1, 2, 3, 4} with sizes in [2, 2]
        """
        return "Set partitions of %s with sizes in %s"%(Set(self._set), self.parts)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: SetPartitions(3, [2,1]).cardinality()
            3
        """
        return len(self.list())

    def __iter__(self):
        """
        An iterator for all the set partitions of the set.

        EXAMPLES::

            sage: SetPartitions(3, [2,1]).list()
            [{{1}, {2, 3}}, {{1, 3}, {2}}, {{1, 2}, {3}}]
        """
        for sp in self._iterator_part(self.parts):
            yield self.element_class(self, sp)

class SetPartitions_setn(SetPartitions):
    @staticmethod
    def __classcall_private__(cls, s, n):
        """
        Normalize ``s`` to ensure a unique representation.

        EXAMPLES::

            sage: S1 = SetPartitions(set([2,1,4]), 2)
            sage: S2 = SetPartitions([4,1,2], 2)
            sage: S3 = SetPartitions((1,2,4), 2)
            sage: S1 is S2, S1 is S3
            (True, True)
        """
        return super(SetPartitions_setn, cls).__classcall__(cls, frozenset(s), n)

    def __init__(self, s, n):
        """
        TESTS::

            sage: S = SetPartitions(5, 3)
            sage: TestSuite(S).run()
        """
        self.n = n
        SetPartitions.__init__(self, s)

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitions(5, 3)
            Set partitions of {1, 2, 3, 4, 5} with 3 parts
        """
        return "Set partitions of %s with %s parts"%(Set(self._set), self.n)

    def cardinality(self):
        """
        The Stirling number of the second kind is the number of partitions
        of a set of size `n` into `k` blocks.

        EXAMPLES::

            sage: SetPartitions(5, 3).cardinality()
            25
            sage: stirling_number2(5,3)
            25
        """
        return stirling_number2(len(self._set), self.n)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: SetPartitions(3).list()
            [{{1, 2, 3}}, {{1}, {2, 3}}, {{1, 3}, {2}}, {{1, 2}, {3}}, {{1}, {2}, {3}}]
        """
        for p in Partitions(len(self._set), length=self.n):
            for sp in self._iterator_part(p):
                yield self.element_class(self, sp)

class SetPartitions_set(SetPartitions):
    """
    Set partitions of a fixed set `S`.
    """
    @staticmethod
    def __classcall_private__(cls, s):
        """
        Normalize ``s`` to ensure a unique representation.

        EXAMPLES::

            sage: S1 = SetPartitions(set([2,1,4]))
            sage: S2 = SetPartitions([4,1,2])
            sage: S3 = SetPartitions((1,2,4))
            sage: S1 is S2, S1 is S3
            (True, True)
        """
        return super(SetPartitions_set, cls).__classcall__(cls, frozenset(s))

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitions([1,2,3])
            Set partitions of {1, 2, 3}
        """
        return "Set partitions of %s"%(Set(self._set))

    def cardinality(self):
        """
        Return the number of set partitions of the set `S`.

        The cardinality is given by the `n`-th Bell number where `n` is the
        number of elements in the set `S`.

        EXAMPLES::

            sage: SetPartitions([1,2,3,4]).cardinality()
            15
            sage: SetPartitions(3).cardinality()
            5
            sage: SetPartitions(3,2).cardinality()
            3
            sage: SetPartitions([]).cardinality()
            1
        """
        return bell_number(len(self._set))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: SetPartitions(3).list()
            [{{1, 2, 3}}, {{1}, {2, 3}}, {{1, 3}, {2}}, {{1, 2}, {3}}, {{1}, {2}, {3}}]
        """
        for p in Partitions(len(self._set)):
            for sp in self._iterator_part(p):
                yield self.element_class(self, sp)

def _listbloc(n, nbrepets, listint=None):
    r"""
    Decompose a set of `n \times n` ``brepets`` integers (the list
    ``listint``) in ``nbrepets`` parts.

    It is used in the algorithm to generate all set partitions.

    .. WARNING::

        Internal function that is not to be called by the user.

    EXAMPLES::

        sage: list(sage.combinat.set_partition._listbloc(2,1))
        [{{1, 2}}]
        sage: l = [Set([Set([3, 4]), Set([1, 2])]), Set([Set([2, 4]), Set([1, 3])]), Set([Set([2, 3]), Set([1, 4])])]
        sage: list(sage.combinat.set_partition._listbloc(2,2,[1,2,3,4])) == l
        True
    """
    if isinstance(listint, (int, Integer)) or listint is None:
        listint = Set(range(1,n+1))

    if nbrepets == 1:
        yield Set([listint])
        return

    l = list(listint)
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
    return Set([result])

def inf(s,t):
    """
    Deprecated in :trac:`14140`. Use :meth:`SetPartition.inf()` instead.

    EXAMPLES::

        sage: sp1 = Set([Set([2,3,4]),Set([1])])
        sage: sp2 = Set([Set([1,3]), Set([2,4])])
        sage: s = Set([ Set([2,4]), Set([3]), Set([1])]) #{{2, 4}, {3}, {1}}
        sage: sage.combinat.set_partition.inf(sp1, sp2) == s
        doctest:1: DeprecationWarning: inf(s, t) is deprecated. Use s.inf(t) instead.
        See http://trac.sagemath.org/14140 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(14140, 'inf(s, t) is deprecated. Use s.inf(t) instead.')
    temp = [ss.intersection(ts) for ss in s for ts in t]
    temp = filter(lambda x: x != Set([]), temp)
    return Set(temp)

def sup(s,t):
    """
    Deprecated in :trac:`14140`. Use :meth:`SetPartition.sup()` instead.

    EXAMPLES::

        sage: sp1 = Set([Set([2,3,4]),Set([1])])
        sage: sp2 = Set([Set([1,3]), Set([2,4])])
        sage: s = Set([ Set([1,2,3,4]) ])
        sage: sage.combinat.set_partition.sup(sp1, sp2) == s
        doctest:1: DeprecationWarning: sup(s, t) is deprecated. Use s.sup(t) instead.
        See http://trac.sagemath.org/14140 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(14140, 'sup(s, t) is deprecated. Use s.sup(t) instead.')
    res = s
    for p in t:
        inters = Set(filter(lambda x: x.intersection(p) != Set([]), list(res)))
        res = res.difference(inters).union(_set_union(inters))
    return res

def standard_form(sp):
    """
    Deprecated in :trac:`14140`. Use :meth:`SetPartition.standard_form()`
    instead.

    EXAMPLES::

        sage: map(sage.combinat.set_partition.standard_form, SetPartitions(4, [2,2]))
        doctest:1: DeprecationWarning: standard_form(sp) is deprecated. Use sp.standard_form() instead.
        See http://trac.sagemath.org/14140 for details.
        [[[1, 2], [3, 4]], [[1, 3], [2, 4]], [[1, 4], [2, 3]]]
    """
    from sage.misc.superseded import deprecation
    deprecation(14140, 'standard_form(sp) is deprecated. Use sp.standard_form() instead.')
    return [list(x) for x in sp]

def less(s, t):
    """
    Deprecated in :trac:`14140`. Use :meth:`SetPartitions.is_less_than()`
    instead.

    EXAMPLES::

        sage: z = SetPartitions(3).list()
        sage: sage.combinat.set_partition.less(z[0], z[1])
        doctest:1: DeprecationWarning: less(s, t) is deprecated. Use SetPartitions.is_less_tan(s, t) instead.
        See http://trac.sagemath.org/14140 for details.
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(14140, 'less(s, t) is deprecated. Use SetPartitions.is_less_tan(s, t) instead.')
    if _union(s) != _union(t):
        raise ValueError("cannot compare partitions of different sets")
    if s == t:
        return False
    for p in s:
        f = lambda z: z.intersection(p) != Set([])
        if len(filter(f, list(t)) ) != 1:
            return False
    return True

def cyclic_permutations_of_set_partition(set_part):
    """
    Returns all combinations of cyclic permutations of each cell of the
    set partition.

    AUTHORS:

    - Robert L. Miller

    EXAMPLES::

        sage: from sage.combinat.set_partition import cyclic_permutations_of_set_partition
        sage: cyclic_permutations_of_set_partition([[1,2,3,4],[5,6,7]])
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
    """
    return list(cyclic_permutations_of_set_partition_iterator(set_part))


def cyclic_permutations_of_set_partition_iterator(set_part):
    """
    Iterates over all combinations of cyclic permutations of each cell
    of the set partition.

    AUTHORS:

    - Robert L. Miller

    EXAMPLES::

        sage: from sage.combinat.set_partition import cyclic_permutations_of_set_partition_iterator
        sage: list(cyclic_permutations_of_set_partition_iterator([[1,2,3,4],[5,6,7]]))
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
    """
    from sage.combinat.permutation import CyclicPermutations
    if len(set_part) == 1:
        for i in CyclicPermutations(set_part[0]):
            yield [i]
    else:
        for right in cyclic_permutations_of_set_partition_iterator(set_part[1:]):
            for perm in CyclicPermutations(set_part[0]):
                yield [perm] + right



