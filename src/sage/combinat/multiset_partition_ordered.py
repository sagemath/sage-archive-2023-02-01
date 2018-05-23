r"""
Ordered Multiset Partitions

An ordered multiset partition `c` of a multiset `X` is a list of nonempty
subsets of `X` (not multisets), called the blocks of `c`, whose multi-union
is `X`.

This module provides tools for manipulating ordered multiset partitions.

.. NOTE::

    An (ordered) multiset partition connotates an (ordered) collection
    of *multisets* (not sets), whose multi-union is a given multiset.
    Perhaps a better name would be "strict" ordered multiset partition.
    Terminology is taken from:

    - [HRWxx] - arXiv:1509.07058
    - [HRSxx] - arXiv:1609.07575
    - [BCHOPSYxx] - arXiv:1707.08709v2

EXAMPLES::

    sage: from sage.combinat.multiset_partition_ordered import *
    sage: OrderedMultisetPartition([[5, 3], [1, 3]])
    [{3,5}, {1,3}]
    sage: list(OrderedMultisetPartitions(4))
    [[{4}], [{1,3}], [{3}, {1}], [{1,2}, {1}], [{2}, {2}], [{2}, {1}, {1}],
     [{1}, {3}], [{1}, {1,2}], [{1}, {2}, {1}], [{1}, {1}, {2}], [{1}, {1}, {1}, {1}]]
    sage: OrderedMultisetPartitions([1,4], 3).list()
    [[{1,4}, {1}], [{1,4}, {4}], [{1}, {1,4}], [{4}, {1,4}], [{1}, {1}, {1}],
     [{1}, {1}, {4}], [{1}, {4}, {1}], [{1}, {4}, {4}], [{4}, {1}, {1}],
     [{4}, {1}, {4}], [{4}, {4}, {1}], [{4}, {4}, {4}]]

AUTHORS:

- Aaron Lauve (2018)
"""
#*****************************************************************************
#       Copyright (C) 2018 Aaron Lauve       <lauve at math.luc.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import, division
from six.moves import range
from six import add_metaclass

from functools import reduce
from itertools import chain

from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.cartesian_product import cartesian_product
from sage.categories.sets_cat import EmptySetError
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.tensor import tensor
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.flatten import flatten
from sage.misc.all import prod
from sage.sets.set import Set, Set_object
from sage.sets.family import Family
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.functions.other import binomial
from sage.calculus.var import var

from sage.combinat.subset import Subsets
from sage.combinat.combinat import CombinatorialElement
from sage.combinat.composition import Composition, Compositions
from sage.combinat.permutation import Permutations_mset
from sage.combinat.partition import RegularPartitions_n
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.shuffle import ShuffleProduct, ShuffleProduct_overlapping
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.combinat.root_system.cartan_type import CartanType

@add_metaclass(InheritComparisonClasscallMetaclass)
class OrderedMultisetPartition(ClonableArray):
    r"""
    Ordered Multiset Partition

    An ordered multiset partition `c` of a multiset `X` is a list of nonempty
    subsets (not multisets), called the blocks of `c`, whose multi-union is `X`.

    EXAMPLES:

    The simplest way to create a ordered multiset partition is by specifying its
    entries as a list, tuple::

        sage: OrderedMultisetPartition([[3],[2,1]])
        [{3}, {1,2}]
        sage: OrderedMultisetPartition(((3,), (1,2)))
        [{3}, {1,2}]
        sage: OrderedMultisetPartition([set([i]) for i in range(2,5)])
        [{2}, {3}, {4}]

    You can also create a ordered multiset partition `c` from a list of positive integers
    and from a list of nonnegative integers. In the former case, each integer is given
    its own block of `c`. In the latter case, zeros separate the blocks of `c`::

        sage: OrderedMultisetPartition([i for i in range(2,5)])
        [{2}, {3}, {4}]
        sage: OrderedMultisetPartition([1, 1, 0, 3, 5, 0, 2, 1])
        Traceback (most recent call last):
        ...
        ValueError: Ordered multiset partitions do not have repeated entries
        within blocks ([[1, 1], [3, 5], [2, 1]] received)
        sage: OrderedMultisetPartition([1, 0, 1, 3, 5, 0, 2, 1])
        [{1}, {1,3,5}, {1,2}]

    TESTS::

        sage: C = OrderedMultisetPartition([[3,1],[2]])
        sage: TestSuite(C).run()
    """
    @staticmethod
    def __classcall_private__(cls, co):
        """
        Create an ordered multiset partition (i.e., a list of sets) from the passed
        arguments with the appropriate parent.

        EXAMPLES::

            sage: OrderedMultisetPartition([[3], [2,1]])
            [{3}, {1,2}]
            sage: OrderedMultisetPartition([2, 3, 4, 5])
            [{2}, {3}, {4}, {5}]
            sage: OrderedMultisetPartition([1, 1, 0, 3, 5, 0, 2, 1])
            Traceback (most recent call last):
            ...
            ValueError: Ordered multiset partitions do not have repeated entries
             within blocks ([[1, 1], [3, 5], [2, 1]] received)
            sage: OrderedMultisetPartition([1, 0, 1, 3, 5, 0, 2, 1])
            [{1}, {1,3,5}, {1,2}]
            sage: OrderedMultisetPartition([2, 3, 4, 5])
            [{2}, {3}, {4}, {5}]
            sage: OrderedMultisetPartition([1, 0, 1, 3, 5, 0, 2, 1])
            [{1}, {1,3,5}, {1,2}]
        """
        if not co:
            P = OrderedMultisetPartitions([])
            return P.element_class(P, [])
        if isinstance(co[0], (list, tuple, set, frozenset, Set_object)):
            # standard input
            X = _concatenate(co)
            P = OrderedMultisetPartitions(_get_weight(X))
            return P.element_class(P, co)
        else:
            # user shortcuts
            weight = _get_weight(c for c in co if c not in {0, '0'})
            P = OrderedMultisetPartitions_X(tuple(weight.iteritems()))
            return P.element_class(P, P.from_list(co))

    def __init__(self, *args):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: OMPs = OrderedMultisetPartitions()
            sage: co = OMPs([[1,3], [1,'a']])
            sage: OrderedMultisetPartition([[1,3], [1,'a']]) == co
            True

        TESTS::

            sage: OMP = OrderedMultisetPartition
            sage: co = OMP([[1,3], [1,'a']])
            sage: TestSuite(co).run()
            sage: co = OMP('abc0ab0a')
            sage: TestSuite(co).run()
            sage: OMP([]) == OMP('')
            True
            sage: OMP([1,2,3,1]) == OMP([[1],[2],[3],[1]])
            True
            sage: OMP([1,2,3,0,1]) == OMP([[1,2,3],[1]])
            True
        """
        if len(args) == 1:
            parent = OrderedMultisetPartitions()
        else:
            parent = args[0]
        co = args[-1]
        if co and not isinstance(co[0], (list, tuple, set, frozenset, Set_object)):
            co = parent.from_list(co)
        ClonableArray.__init__(self, parent, [Set(list(k)) for k in co])
        self._multiset = _get_multiset(co)
        self._weight = _get_weight(self._multiset)
        self._order = sum([len(block) for block in self])
        if all((a in ZZ and a > 0) for a in self._multiset):
            self._n = ZZ(sum(self._multiset))
        else:
            self._n = None

    def check(self):
        """
        Check that we are a valid ordered multiset partition.

        .. TODO:: make tests work.

        EXAMPLES::

            sage: OMP4 = OrderedMultisetPartitions(4)
            sage: co = OMP4([[1], [1,2]])
            sage: co.check()

            sage: OMP = OrderedMultisetPartitions()
            sage: co = OMP([[1], [1], ['a']])
            sage: co.check()

        TESTS::

            sage: co = OMP4([[1, 3], [1, 4]], check=False)
            sage: co.check()
            Traceback (most recent call last):
            ...
            AssertionError: [{1,3}, {1,4}] is not an element
             of Ordered Multiset Partitions of integer 4

            sage: co = OMP([[1, 1], [1, 4]])
            Traceback (most recent call last):
            ...
            AssertionError: cannot convert [[1, 1], [1, 4]] into an element
             of Ordered Multiset Partitions
        """
        if self not in self.parent():
            raise ValueError("{} not an element of {}".format(self, self.parent()))

    def _repr_(self):
        return self._repr_tight()

    def _repr_normal(self):
        # TODO: simplify if/once ``_repr_`` method for ``Set`` sorts its elements.
        string_parts = map(lambda k: str(sorted(k)), self)
        string_parts = ", ".join(string_parts).replace("[","{").replace("]","}")
        return "[" + string_parts + "]"

    def _repr_tight(self):
        repr = self._repr_normal()
        # eliminate spacing within blocks
        return repr.replace(", ", ",").replace("},{", "}, {")

    def __hash__(self):
        """
        Return the hash of ``self``.

        The parent is not included as part of the hash.

        EXAMPLES::

            sage: OMP = OrderedMultisetPartitions(4)
            sage: A = OMP([[1], [1, 2]])
            sage: B = OMP([{1}, {1, 2}])
            sage: hash(A) == hash(B)
            True
        """
        return sum(hash(x) for x in self)

    def __eq__(self, y):
        """
        Check equality of ``self`` and ``y``.

        The parent is not included as part of the equality check.

        EXAMPLES::

            sage: OMP = OrderedMultisetPartitions(4)
            sage: A = OMP([[1], [1, 2]])
            sage: B = OMP([{1}, {3}])
            sage: A == B
            False
            sage: C = OMP([[1,2], [1]])
            sage: A == C
            False
            sage: D = OMP.from_zero_list([1,0,1,2])
            sage: A == D
            True
        """
        if not isinstance(y, OrderedMultisetPartition):
            return False
        return list(self) == list(y)

    def __ne__(self, y):
        """
        Check lack of equality of ``self`` and ``y``.

        The parent is not included as part of the equality check.
        """
        return not (self == y)

    def __add__(self, other):
        """
        Return the concatenation of two ordered multiset partitions.

        EXAMPLES::

            sage: OrderedMultisetPartition([1, 1, 3]) + OrderedMultisetPartition([4, 1, 0, 2])
            [{1}, {1}, {3}, {1,4}, {2}]

        TESTS::

            sage: OrderedMultisetPartition([]) + OrderedMultisetPartition([]) == OrderedMultisetPartition([])
            True
            sage: OrderedMultisetPartition([1,0,1,0,3]) + OrderedMultisetPartition([1,4,0,2]) == OrderedMultisetPartition([{1}, {1}, {3}, {1,4}, {2}])
            True
        """
        co = list(self) + list(other)
        X = _concatenate(co)
        return OrderedMultisetPartitions(_get_weight(X))(co)

    @combinatorial_map(order=2, name='reversal')
    def reversal(self):
        r"""
        Return the reverse ordered multiset partition of ``self``.

        The reverse of a ordered multiset partition `(B_1, B_2, \ldots, B_k)`
        is defined as the ordered multiset partition `(B_k, B_{k-1}, \ldots, B_1)`.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([1, 0, 1, 3, 0, 2, 3, 4]); C
            [{1}, {1,3}, {2,3,4}]
            sage: C.reversal()
            [{2,3,4}, {1,3}, {1}]
        """
        return self.parent()(list(reversed(self)))

    def shape_from_cardinality(self):
        r"""
        Return a composition that records the cardinality of each block of ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([3, 4, 1, 0, 2, 0, 1, 2, 3, 7]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.shape_from_cardinality()
            [3, 1, 4]
            sage: OrderedMultisetPartition([]).shape_from_cardinality() == Composition([])
            True
        """
        return Composition([len(k) for k in self])

    def shape_from_size(self):
        r"""
        Return a composition that records the sum of entries of each block of ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([3, 4, 1, 0, 2, 0, 1, 2, 3, 7]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.shape_from_size()
            [8, 2, 13]

        TESTS::

            sage: OrderedMultisetPartition([]).shape_from_size() == Composition([])
            True
            sage: D = OrderedMultisetPartition([['a', 'b'], ['a']]); D
            [{'a','b'}, {'a'}]
            sage: D.shape_from_size() == None
            True
        """
        if self._n is not None:
            return Composition([sum(k) for k in self])

    def letters(self):
        r"""
        Return the set of distinct elements occurring within the blocks of ``self``.

        .. SEEALSO::

            :meth:`Words.letters`
        """
        return _union_of_sets(list(self))

    def multiset(self, as_dict=False):
        r"""
        Return the multiset corresponding to ``self`` as a list or as a dictionary.

        Counts frequency of each letter in ``self.letters()``

        EXAMPLES::

            sage: C = OrderedMultisetPartition([3, 4, 1, 0, 2, 0, 1, 2, 3, 7]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.multiset()
            (1, 1, 2, 2, 3, 3, 4, 7)
            sage: C.multiset(as_dict=True)
            {1:2, 2:2, 3:3, 4:1, 7:1}
            sage: OrderedMultisetPartition([]).multiset() == ()
            True
        """
        if as_dict:
            return self._weight
        else:
            return self._multiset

    def max_letter(self):
        r"""
        Return the maximum letter appearing in ``self.letters()``.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([3, 4, 1, 0, 2, 0, 1, 2, 3, 7])
            sage: C.max_letter()
            7
            sage: D = OrderedMultisetPartition('abc0ab0a0bcf0cd')
            sage: D.max_letter()
            'f'
        """
        if not self.letters():
            return None
        else:
            return max(self.letters())

    def size(self):
        """
        Return the size of ``self``, that is the sum of all integers in all blocks.
        Else, return ``None``.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([3, 4, 1, 0, 2, 0, 1, 2, 3, 7]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.size()
            23
            sage: C.size() == sum(k for k in C.shape_from_size())
            True
            sage: OrderedMultisetPartition([[7,1],[3]]).size()
            11
            sage: OrderedMultisetPartition([]).size() == 0
            True
        """
        return self._n

    def order(self):
        """
        Return the total number of integers in all blocks.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([3, 4, 1, 0, 2, 0, 1, 2, 3, 7]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.order()
            8
            sage: C.order() == sum(C.weight().values())
            True
            sage: C.order() == sum(k for k in C.shape_from_cardinality())
            True
            sage: OrderedMultisetPartition([[7,1],[3]]).order()
            3
        """
        return self._order

    def length(self):
        """
        Return the number of blocks.

        EXAMPLES::

            sage: OrderedMultisetPartition([[7,1],[3]]).length()
            2
        """
        return len(self)

    def weight(self, as_weak_comp=False):
        """
        Return a dictionary, with keys being the letters in ``self.letters()``
        and values being their (positive) frequency.

        Alternatively, if ``as_weak_comp`` is ``True``, count the number of instances
        `n_i` for each distinct positive integer `i` across all blocks of ``self``.
        Return as a list `[n_1, n_2, n_3, ..., n_k]`, where `k` is the max letter
        appearing in ``self.letters()``.

        EXAMPLES::

            sage: OrderedMultisetPartition([[6,1],[1,3],[1,3,6]]).weight(as_weak_comp=True)
            [3, 0, 2, 0, 0, 2]
            sage: OrderedMultisetPartition([[6,1],[1,3],[1,3,6]]).weight()
            {1:3, 3:2, 6:2}
            sage: OrderedMultisetPartition([]).weight() == {}
            True
        """
        w = self._weight
        if as_weak_comp:
            w = [w.get(i, 0) for i in range(1,self.max_letter()+1)]
        return w

    def deconcatenate(self, k=2):
        """
        Choose `k-1` cut points to deconcatenate the ordered list ``self`` into
        k (possibly) empty lists. (A "balls and boxes" operation.)

        This is not to be confused with ``self.split()``,
        which splits each block of ``self`` before making k-tuples of ordered multiset partitions.

        Returns a set of k-tuples of ordered multiset partitions.

        EXAMPLES::

            sage: OrderedMultisetPartition([[7,1],[3]]).deconcatenate()
            [([], [{1,7}, {3}]), ([{1,7}], [{3}]), ([{1,7}, {3}], [])]
            sage: OrderedMultisetPartition('bc0a').deconcatenate()
            [([], [{'b','c'}, {'a'}]), ([{'b','c'}], [{'a'}]), ([{'b','c'}, {'a'}], [])]
            sage: OrderedMultisetPartition('abc0').deconcatenate(3)
            [([], [], [{'a','b','c'}]),
             ([], [{'a','b','c'}], []),
             ([{'a','b','c'}], [], [])]

        TESTS::

            sage: C = OrderedMultisetPartition('abcde'); C
            [{'a'}, {'b'}, {'c'}, {'d'}, {'e'}]
            sage: all( C.deconcatenate(k).cardinality()
            ....:      == binomial(C.length() + k-1, k-1)
            ....:      for k in range(1, 5) )
            True
        """
        P = OrderedMultisetPartitions(alphabet=self.letters(),max_length=self.length())
        out = []
        for c in IntegerListsLex(self.length(), length=k):
            ps = [sum(c[:i]) for i in range(k+1)]
            out.append(tuple([P(self[ps[i]:ps[i+1]]) for i in range(len(ps)-1)]))
        return Set(out)

    def split(self, k=2):
        """
        Return a set of `k`-tuples of ordered multiset partitions, each member of which is
        formed by first splitting each block of ``self`` into `k` (possibly empty) sets.

        First: form a `k`-splitting `(b_i^1,\ldots, b_i^k)` for each block `B_i`
        of ``self``; that is, `b_i^1 \cup \ldots \cup b_i^k = B_i`.
        Next, the `k`-tuple `(A_1, \ldots, A_k)` is formed by first creating (for all `j`),
        the list `[b_1^j, b_2^j, \ldots]`, then excising any empty sets `b_i^j` appearing
        here before creating the ordered multiset partition `A_j`.

        This is not to be confused with ``self.deconcatenate()``, which serves as
        the preimage of the ``+`` operation on ordered multiset partitions.

        EXAMPLES::

            sage: sorted(OrderedMultisetPartition([[1,2],[4]]).split())
            [([], [{1,2}, {4}]), ([{2}, {4}], [{1}]), ([{1,2}], [{4}]),
             ([{4}], [{1,2}]), ([{1}], [{2}, {4}]), ([{1}, {4}], [{2}]),
             ([{2}], [{1}, {4}]), ([{1,2}, {4}], [])]
            sage: sorted(OrderedMultisetPartition([[1,2]]).split(3))
            [([], [], [{1,2}]), ([], [{2}], [{1}]), ([], [{1}], [{2}]),
             ([], [{1,2}], []), ([{2}], [{1}], []), ([{1}], [], [{2}]),
             ([{1}], [{2}], []), ([{2}], [], [{1}]), ([{1,2}], [], [])]

        TESTS::

            sage: C = OrderedMultisetPartition([1,2,0,4,5,6]); C
            [{1,2}, {4,5,6}]
            sage: C.split().cardinality() == 2**len(C[0]) * 2**len(C[1])
            True
            sage: C.split(3).cardinality() == (1+2)**len(C[0]) * (1+2)**len(C[1])
            True
            sage: C = OrderedMultisetPartition([])
            sage: C.split(3) == Set([(C, C, C)])
            True
        """
        P = OrderedMultisetPartitions(alphabet=self.letters(),max_length=self.length())

        # corner case
        if not self:
            return Set([tuple([self]*k)])
        else:
            out = set()
            tmp = cartesian_product([_split_block(block, k) for block in self])
            for t in tmp:
                out.add(tuple([P([k for k in c if len(k)>0]) for c in zip(*t)]))
            return Set(out)

    def finer(self, strong=False):
        """
        Return the set of ordered multiset partitions that are finer than ``self``.

        An ordered multiset partition `A` is finer than another `B` if, reading left-to-right,
        every block of `B` is the union of some consecutive blocks of `A`.

        If optional argument ``strong`` is set to ``True``, then return only those `A`
        where blocks are deconcatenation of blocks of `B`.


        EXAMPLES::

            sage: C = OrderedMultisetPartition([[3,2]]).finer()
            sage: C.cardinality()
            3
            sage: C.list()
            [[{2}, {3}], [{2,3}], [{3}, {2}]]
            sage: OrderedMultisetPartition([]).finer()
            {[]}
            sage: O = OrderedMultisetPartitions([1, 1, 'a', 'b'])
            sage: o = O.an_element()
            sage: o.finer()
            {[{'b'}, {'a'}, {1}, {1}], [{'a'}, {'b'}, {1}, {1}], [{'a','b'}, {1}, {1}]}
            sage: o.finer() & o.fatter() == Set([o])
            True
        """
        P = OrderedMultisetPartitions(self._multiset)

        if not self:
            return Set([self])
        else:
            tmp = cartesian_product([_refine_block(block, strong) for block in self])
            return Set([P(_concatenate(map(list,c))) for c in tmp])

    def is_finer(self, co):
        """
        Return ``True`` if the ordered multiset partition ``self`` is finer than the
        composition ``co``; otherwise, return ``False``.

        EXAMPLES::

            sage: OrderedMultisetPartition([[4],[1],[2]]).is_finer([[1,4],[2]])
            True
            sage: OrderedMultisetPartition([[1],[4],[2]]).is_finer([[1,4],[2]])
            True
            sage: OrderedMultisetPartition([[1,4],[1],[1]]).is_finer([[1,4],[2]])
            False
        """
        X = _concatenate(co)
        if self.weight() != OrderedMultisetPartitions(_get_weight(X))(co).weight():
            return False

        # trim common prefix and suffix to make the search-space smaller
        co1 = map(Set,self)
        co2 = map(Set,co)
        while co1[0] == co2[0]:
            co1 = co1[1:]; co2 = co2[1:]
        while co1[-1] == co2[-1]:
            co1 = co1[:-1]; co2 = co2[:-1]

        co1 = OrderedMultisetPartition(co1)
        co2 = OrderedMultisetPartition(co2)
        return co1 in co2.finer()

    def fatten(self, grouping):
        r"""
        Return the ordered multiset partition fatter than ``self``, obtained by grouping
        together consecutive parts according to ``grouping`` (whenever this does not violate the strictness condition).

        INPUT:

        - ``grouping`` -- a composition (or list) whose sum is the length of ``self``

        EXAMPLES:

        Let us start with the composition::

            sage: C = OrderedMultisetPartition([4,1,5,0,2,0,7,1]); C
            [{1,4,5}, {2}, {1,7}]

        With ``grouping`` equal to `(1, 1, 1)`, `C` is left unchanged::

            sage: C.fatten([1,1,1])
            [{1,4,5}, {2}, {1,7}]

        With ``grouping`` equal to `(2,1)` or `(1,2)`, a union of consecutive parts is achieved::

            sage: C.fatten([2,1])
            [{1,2,4,5}, {1,7}]
            sage: C.fatten([1,2])
            [{1,4,5}, {1,2,7}]

        However, the ``grouping`` `(3)` will throw an error, as `1` cannot appear twice in any block of ``C``::

            sage: C.fatten(Composition([3]))
            Traceback (most recent call last):
            ...
            ValueError: [{1,4,5,2,1,7}] is not a valid ordered multiset partition
        """
        if sum(list(grouping)) != self.length():
            raise ValueError("%s is not a composition of ``self.length()`` (=%s)"%(grouping, self.length()))

        valid = True
        result = []
        for i in range(len(grouping)):
            result_i = self[sum(grouping[:i]) : sum(grouping[:i+1])]
            # check that grouping[i] is allowed, i.e., `|A\cup B| = |A| + |B|`
            strict_size = sum(map(len,result_i))
            size = len(_union_of_sets(result_i))
            if size < strict_size:
                valid = False
            result.append(_concatenate(result_i))
        if not valid:
            str_rep = '['
            for i in range(len(grouping)):
                st = ",".join(str(k) for k in result[i])
                str_rep += "{" + st+ "}"
            str_rep = str_rep.replace("}{", "}, {") + "]"
            raise ValueError("%s is not a valid ordered multiset partition"%(str_rep))
        else:
            return OrderedMultisetPartitions(self._multiset)(result)

    def fatter(self):
        """
        Return the set of ordered multiset partitions which are fatter than ``self``.

        An ordered multiset partition `A` is fatter than another `B` if, reading left-to-right,
        every block of `A` is the union of some consecutive blocks of `B`.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([{1,4,5}, {2}, {1,7}]).fatter()
            sage: C.cardinality()
            3
            sage: list(C)
            [[{1,4,5}, {2}, {1,7}], [{1,4,5}, {1,2,7}], [{1,2,4,5}, {1,7}]]
            sage: list(OrderedMultisetPartition([['a','b'],['c'],['a']]).fatter())
            [[{'a','b'}, {'a','c'}], [{'a','b','c'}, {'a'}], [{'a','b'}, {'c'}, {'a'}]]

        Some extreme cases::

            sage: list(OrderedMultisetPartition('abc0').fatter())
            [[{'a','b','c'}]]
            sage: list(OrderedMultisetPartition([]).fatter())
            [[]]
            sage: OrderedMultisetPartition([1,2,3,4]).fatter().issubset(OrderedMultisetPartition([[1,2,3,4]]).finer())
            True
        """
        out = set()
        for c in Compositions(self.length()):
            try:
                out.add(self.fatten(c))
            except ValueError:
                pass
        return Set(out)


    def minimaj(self):
        r"""
        Return the minimaj statistic on ordered multiset partitions.

        We define `minimaj` via an example:

        1. Sort the block in ``self`` as prescribed by ``self.minimaj_ordering()``,
           keeping track of the original separation into blocks.
           - in:   [{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}]
           - out:  [ 5,7,1 /  2,4 /  5,6 /  4,6,8 /  3,1 /  1,2,3 ]

        2. Record the indices where descents in this word occur.
           - word:      [5, 7, 1 / 2, 4 / 5, 6 / 4, 6, 8 / 3, 1 / 1, 2, 3]
           - indices:    1  2  3   4  5   6  7   8  9 10  11 12  13 14 15
           - descents:  {   2,               7,       10, 11             }

        3. Compute the sum of the descents
           - minimaj = 2 + 7 + 10 + 11 = 30

        See:
        [HRWxx] - arXiv:1509.07058

        EXAMPLES::

            sage: C = OrderedMultisetPartition([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}])
            sage: C, C.minimaj_word()
            ([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}],
             [5, 7, 1, 2, 4, 5, 6, 4, 6, 8, 3, 1, 1, 2, 3])
            sage: C.minimaj()
            30
            sage: C = OrderedMultisetPartition([{2,4}, {1,2,3}, {1,6,8}, {2,3}])
            sage: C, C.minimaj_word()
            ([{2,4}, {1,2,3}, {1,6,8}, {2,3}], [2, 4, 1, 2, 3, 6, 8, 1, 2, 3])
            sage: C.minimaj()
            9
            sage: OrderedMultisetPartition([]).minimaj()
            0
            sage: C = OrderedMultisetPartition('bd0abc0b')
            sage: C, C.minimaj_word()
            ([{'b','d'}, {'a','b','c'}, {'b'}], ['d', 'b', 'c', 'a', 'b', 'b'])
            sage: C.minimaj()
            4
        """
        D = _descents(self.minimaj_word())
        return sum(D) + len(D)

    def minimaj_word(self):
        r"""
        Return an ordering of ``self._multiset`` derived by concatenating the blocks of
        ``self`` after first ordering them in a prescribed order, which we define
        via an example below.

        Sort the blocks `[B_1,...,B_k]` of ``self`` from right to left via:

        1. Sort the last block `B_k` in increasing order, call it the word `W_k`

        2. If blocks `B_{i+1}, \ldots, B_k` have been converted to words
           `W_{i+1}, \ldots, W_k`, use the letters in `B_i` to make the unique
           word `W_i` that has a factorization `W_i=(u,v)` satisfying:
           - letters of `u` and `v` appear in increasing order, with `v` possibly empty
           - letters in `vu` appear in increasing order
           - ``v[-1]`` is the largest letter `a \in B_i` satisfying ``a <= W_{i+1}[0]``

        EXAMPLES::

            sage: C = OrderedMultisetPartition([2,1,0,1,2,3,0,1,2,0,3,0,1]); C
            [{1,2}, {1,2,3}, {1,2}, {3}, {1}]
            sage: C.minimaj_word()
            [1, 2, 2, 3, 1, 1, 2, 3, 1]
            sage: C = OrderedMultisetPartition([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}]); C
            [{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}]
            sage: C.minimaj_word()
            [5, 7, 1, 2, 4, 5, 6, 4, 6, 8, 3, 1, 1, 2, 3]
            sage: OrderedMultisetPartition([]).minimaj_word()
            []
        """
        return [_ for block in self.minimaj_blocks() for _ in block]

    def minimaj_blocks(self):
        r"""
        Return the ordering on blocks of ``self`` leading to ``self.minimaj_word()``.
        """
        if len(self) == 0:
            return []

        C = [sorted(self[-1])]
        for i in range(1,len(self)):
            lower = []; upper = []
            for j in self[-1-i]:
                if j <= C[0][0]:
                    lower.append(j)
                else:
                    upper.append(j)
            C = [sorted(upper)+sorted(lower)] + C
        return tuple(map(tuple,C))

    def to_tableau(self):
        r"""
        Return a sequence of row words corresponding to (skew-)tableaux.

        Output is the minimaj bijection `\varphi` of [BCHOPSYxx]_
        applied to ``self``.

        .. TODO::

            Implement option for mapping to sequence of (skew-)tableaux?

        EXAMPLES::

            sage: mu = ((1,2,4),(4,5),(3,),(4,6,1),(2,3,1),(1,),(2,5))
            sage: to_tableau(mu)
            [[5, 1], [3, 1], [6], [5, 4, 2], [1, 4, 3, 4, 2, 1, 2]]
        """
        if not self:
            return []
        bb = self.minimaj_blocks()
        b = [block[0] for block in bb]
        beginning = _partial_sum(self.shape_from_cardinality())
        w = _concatenate(bb)
        D = [0] + _descents(w) + [len(w)]
        pieces = [b]
        for i in range(len(D)-1):
            p = [w[j] for j in range(D[i]+1,D[i+1]+1) if j not in beginning]
            pieces = [p[::-1]] + pieces
        return pieces

    def major_index(self):
        r"""
        A statistic on ordered multiset partitions, which we define via an example.

        1. Sort each block in the list ``self`` in descending order to create
           a word `w`, keeping track of the original separation into blocks:

           - in:  [{3,4,5}, {2,3,4}, {1}, {4,5}]
           - out: [ 5,4,3 /  4,3,2 /  1 /  5,4 ]

        2. Create a sequence `v = (v_0, v_1, v_2, \ldots)`  of length
           ``self.order()+1``, built recursively by:

           - `v_0 = 0`
           - `v_j = v_{j-1} + \delta(j)`, where `\delta(j) = 1` if `j` is
             the index of an end of a block, and zero otherwise.

             * in:    [ 5,4,3 /  4,3,2 /  1 /  5,4]
             * out: (0, 0,0,1,   1,1,2,   3,   3,4)

        3. Compute `\sum_j v_j`, restricted to descent positions in `w`, i.e.,
           sum over those `j` with `w_j > w_{j+1}`:

           - in:  w:   [5, 4, 3, 4, 3, 2, 1, 5, 4]
                  v: (0 0, 0, 1, 1, 1, 2, 3, 3, 4)
           - maj :=     0 +0    +1 +1 +2    +3     = 7

        See:

        [Wil16] -  An extension of MacMahon's equidistribution theorem to ordered multiset partitions

        EXAMPLES::

            sage: C = OrderedMultisetPartition([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}])
            sage: C.major_index()
            27
            sage: C = OrderedMultisetPartition([{3,4,5}, {2,3,4}, {1}, {4,5}])
            sage: C.major_index()
            7
        """
        ew = [enumerate(sorted(k)) for k in self]
        w = []
        v = [0]
        for eblock in ew:
            for (i,wj) in sorted(eblock, reverse=True):
                vj = v[-1]
                if i == 0:
                    vj += 1
                v.append(vj)
                w.append(wj)
        maj = [v[j+1] for j in range(len(w)-1) if w[j] > w[j+1]]
        return sum(maj)

    def _multiset_partition(self):
        """
        Return the (strict) (unordered) multiset partition, realized as a list of sets
        (built by sorting ``self`` by block size and least element in each block).

        EXAMPLES::

            sage: C = OrderedMultisetPartition([2,1,0,1,2,3,0,1,2,0,3,0,1]); C
            [{1,2}, {1,2,3}, {1,2}, {3}, {1}]
            sage: C._multiset_partition()
            [{1}, {3}, {1,2}, {1,2}, {1,2,3}]
        """
        return sorted(self, key=lambda x: (len(x),min(x)))


    def shuffle_product(self, other, overlap=False):
        r"""
        Return the shuffles (with multiplicity) of parts of ``self`` with parts of ``other``.
        In case overlap is True, it also takes the union of adjacent disjoint parts (one from ``self`` and the other from ``other``)

        .. SEEALSO::

            :meth:`Composition.shuffle`

        EXAMPLES::

            sage: A = OrderedMultisetPartition([2,1,3,0,1,2]); A
            [{1,2,3}, {1,2}]
            sage: B = OrderedMultisetPartition([[3,4]]); B
            [{3,4}]
            sage: C = OrderedMultisetPartition([[4,5]]); C
            [{4,5}]
            sage: list(A.shuffle_product(B))
            [[{1,2,3}, {1,2}, {3,4}], [{3,4}, {1,2,3}, {1,2}], [{1,2,3}, {3,4}, {1,2}]]
            sage: list(A.shuffle_product(B, overlap=True))
            [[{1,2,3}, {1,2}, {3,4}], [{1,2,3}, {3,4}, {1,2}],
             [{3,4}, {1,2,3}, {1,2}], [{1,2,3}, {1,2,3,4}]]

            sage: list(A.shuffle_product(C, overlap=True))
            [[{1,2,3}, {1,2}, {4,5}], [{1,2,3}, {4,5}, {1,2}], [{4,5}, {1,2,3}, {1,2}],
             [{1,2,3,4,5}, {1,2}], [{1,2,3}, {1,2,4,5}]]
        """
        other = OrderedMultisetPartition(other)
        if overlap is True:
            for term in ShuffleProduct_overlapping(self, other):
                if list(term._multiset) == sorted(self._multiset + other._multiset):
                    yield term
        else:
            for term in ShuffleProduct(self, other):
                yield term

    # A sorting scheme.
    # An alternative to that served up by ``OrderedMultisetPartitions().__iter__()``
    def soll_gt(self, other):
        """
        Return ``True`` if the composition ``self`` is greater than the
        composition ``other`` with respect to the soll-ordering; otherwise,
        return ``False``.

        The soll-ordering is a total order on the set of all ordered multiset partitions.
        ("soll" is short for "size, order, length, lexicographic".)

        An ordered multiset partition `I` is greater than an ordered multiset partition
        `J` if and only if one of the following conditions holds:

        - The size of `I` is greater than size of `J`. (Recall size of `I` is ``None``
          if `I` is not numeric.)

        - The sizes of `I` and `J` are equal, and the order of `I` is
          greater than the order of `J`.

        - The sizes and orders of `I` and `J` coincide, but the length of `I`
          is greater than the length of `J`.

        - The sizes/orders/lengths of `I` and `J` coincide, but `I` is lexicographically
          greater than `J`. (Here, lexicographic order is taken with respect to blocks,
          not entries within the blocks.)

        EXAMPLES::

            sage: C = OrderedMultisetPartition([2,1,0,1,2,3,0,1,2,0,3,0,1])
            sage: D1 = OrderedMultisetPartition([2,4,0,2,3])
            sage: D2 = OrderedMultisetPartition([{1,2}, {1,2,3}, {1,2}, {1}, {1,2}])
            sage: D3 = OrderedMultisetPartition([1,0,1,2,3,0,1,2,3,0,1,2])
            sage: D4 = OrderedMultisetPartition([1,2,0,2,3,0,1,3,0,1,2,0,1])
            sage: C.soll_gt(D1) #size
            True
            sage: C.soll_gt(D2) #order
            False
            sage: C.soll_gt(D3) #length
            True
            sage: C.soll_gt(D4) #lex
            True

        TESTS::

            sage: C = OrderedMultisetPartition('ba0abc0ab0c0a')
            sage: D2 = OrderedMultisetPartition('ab0abc0ab0a0ab')
            sage: D3 = OrderedMultisetPartition('a0abc0abc0ab')
            sage: D4 = OrderedMultisetPartition('ab0bc0ac0ab0a')
            sage: C.soll_gt(D2) #order
            False
            sage: C.soll_gt(D3) #length
            True
            sage: C.soll_gt(D4) #lex
            False
        """
        co1 = self
        co2 = OrderedMultisetPartition(other)
        # check size
        if (co1).size() > (co2).size():
            return True
        elif (co1).size() < (co2).size():
            return False
        # check order
        if (co1).order() > (co2).order():
            return True
        elif (co1).order() < (co2).order():
            return False
        # check length
        if (co1).length() > (co2).length():
            return True
        if (co1).length() < (co2).length():
            return False
        # check lex.
        if co1._n and co2._n:
            # give preference to blocks with larger sum.
            co1 = [Word([sum(k)]+sorted(k)) for k in co1]
            co2 = [Word([sum(k)]+sorted(k)) for k in co2]
        else:
            co1 = [Word(sorted(k)) for k in co1]
            co2 = [Word(sorted(k)) for k in co2]
        return co1 > co2

##############################################################

class OrderedMultisetPartitions(UniqueRepresentation, Parent):
    r"""
    Ordered Multiset Partitions

    An ordered multiset partition `c` of a multiset `X` is a list of nonempty subsets
    (not multisets), called the blocks of `c`, whose multi-union is `X`.

    Alternatively:
    An ordered multiset partition `c` of a nonnegative integer `n` is a list of
    subsets of positive integers (called blocks) with total sum of all integers
    appearing in all blocks equal to `n`.

    .. TODO::

    - Fix ``check`` error.
      See, e.g., below, which should be the empty list.
      sage: OrderedMultisetPartitions(4, weight=[0,0,0,0,1]).list()

    INPUT:

    Expects one or two arguments, with different behaviors resulting:
    - One Argument:
        + `X` -- a dictionary (representing a multiset for `c`),
                  or an integer (representing the size of `c`)
    - Two Arguments:
        + `alph` -- a list (representing allowable letters within `c`),
                    or a positive integer (representing the maximal allowable letter)
        + `ord`  -- a nonnegative integer (the total number of letters within `c`)

    Optional keyword arguments are as follows:
    (See corresponding methods in see :class:`OrderedMultisetPartition` for more details.)

    - ``weight=X``     (iterable or dictionary `X`) specifies the multiset for `c`
    - ``size=n``       (integer `n`) specifies sum of integers across all blocks in the partition
    - ``alphabet=S``   (iterable `S`) specifies allowable elements for the blocks of `c`
    - ``length=k``     (integer `k`) specifies the number of blocks in the partition
    - ``min_length=k`` (integer `k`) specifies minimum number of blocks in the partition
    - ``max_length=k`` (integer `k`) specifies maximum number of blocks in the partition
    - ``order=n``      (integer `n`) specifies the number of integers in the partition (=sum of cardinalities of its blocks)
    - ``min_order=n``  (integer `n`) specifies minimum number of integers in the partition (=sum of cardinalities of its blocks)
    - ``max_order=n``  (integer `n`) specifies maximum number of integers in the partition (=sum of cardinalities of its blocks)

    EXAMPLES:

    There are 5 ordered multiset partitions of multiset `{{1, 1, 4}}`::

        sage: OrderedMultisetPartitions([1,1,4]).cardinality()
        5

    Here is the list of them::

        sage: OrderedMultisetPartitions([1,1,4]).list()
        [[{1}, {1}, {4}], [{1}, {1,4}], [{1}, {4}, {1}], [{1,4}, {1}], [{4}, {1}, {1}]]

    There are 5 ordered multiset partitions of integer 3::

        sage: OrderedMultisetPartitions(3).cardinality()
        5

    Here is the list of them::

        sage: OrderedMultisetPartitions(3).list()
        [[{3}], [{1,2}], [{2}, {1}], [{1}, {2}], [{1}, {1}, {1}]]

    You can use the ``.first()`` method to get the 'first' ordered multiset partition of
    a number::

        sage: OrderedMultisetPartitions(3).first()
        [{3}]

    You can also calculate the 'next' ordered multiset partition given the current
    one::

        sage: P = OrderedMultisetPartitions(3)
        sage: P.next(P([[2],[1]]))
        [{1}, {2}]

    If no arguments are specified, this returns the combinatorial class of
    all ordered multiset partitions::

        sage: OrderedMultisetPartitions()
        Ordered Multiset Partitions
        sage: [] in OrderedMultisetPartitions()
        True
        sage: [[2,3], [1]] in OrderedMultisetPartitions()
        True
        sage: [['a','b'], ['a']] in OrderedMultisetPartitions()
        True
        sage: [[-2,3], [3]] in OrderedMultisetPartitions()
        True
        sage: [[2], [3,3]] in OrderedMultisetPartitions()
        False

    If `n` is specified, it returns the class of ordered multiset partitions of `n`::

        sage: OrderedMultisetPartitions(3)
        Ordered Multiset Partitions of integer 3
        sage: list(OrderedMultisetPartitions(3))
        [[{3}], [{1,2}], [{2}, {1}], [{1}, {2}], [{1}, {1}, {1}]]
        sage: OrderedMultisetPartitions(3).cardinality()
        5

    If, alternatively, two arguments `alpha` and `ord` are specified, it returns the class of
    ordered multiset partitions of order `ord` over the alphabet `alph`::

        sage: OrderedMultisetPartitions([1, 4], 2)
        Ordered Multiset Partitions of order 2 over alphabet {1, 4}
        sage: list(OrderedMultisetPartitions([1, 4], 2))
        [[{1,4}], [{1}, {1}], [{1}, {4}], [{4}, {1}], [{4}, {4}]]
        sage: OrderedMultisetPartitions([1, 4], 2).cardinality()
        5

    The following examples show how to test whether or not an object
    is an ordered multiset partition::

        sage: [[3,2],[2]] in OrderedMultisetPartitions()
        True
        sage: [[3,2],[2]] in OrderedMultisetPartitions(7)
        True
        sage: [[3,2],[2]] in OrderedMultisetPartitions([2,2,3])
        True
        sage: [[3,2],[2]] in OrderedMultisetPartitions(5)
        False

    The options ``length``, ``min_length``, and ``max_length`` can be used
    to set length constraints on the ordered multiset partitions. For example, the
    ordered multiset partitions of 4 of length equal to, at least, and at most 2 are
    given by::

        sage: OrderedMultisetPartitions(4, length=2).list()
        [[{3}, {1}], [{1,2}, {1}], [{2}, {2}], [{1}, {3}], [{1}, {1,2}]]
        sage: OrderedMultisetPartitions(4, min_length=3).list()
        [[{2}, {1}, {1}], [{1}, {2}, {1}], [{1}, {1}, {2}], [{1}, {1}, {1}, {1}]]
        sage: OrderedMultisetPartitions(4, max_length=2).list()
        [[{4}], [{1,3}], [{3}, {1}], [{1,2}, {1}], [{2}, {2}], [{1}, {3}], [{1}, {1,2}]]

    The option ``alphabet`` constrains which integers appear across all blocks of
    the ordered multiset partition. For example, the ordered multiset partitions of 4
    are listed for different choices of alphabet below. Note that ``alphabet``
    is allowed to be an integer or an iterable::

        sage: OMPs = OrderedMultisetPartitions
        sage: OMPs(4, alphabet=3).list()
        [[{1,3}], [{3}, {1}], [{1,2}, {1}], [{2}, {2}], [{2}, {1}, {1}], [{1}, {3}],
         [{1}, {1,2}], [{1}, {2}, {1}], [{1}, {1}, {2}], [{1}, {1}, {1}, {1}]]
        sage: OMPs(4, alphabet=3) == OMPs(4, alphabet=[1,2,3])
        True
        sage: OMPs(4, alphabet=[3]).list()
        []
        sage: OMPs(4, alphabet=[1,3]).list()
        [[{1,3}], [{3}, {1}], [{1}, {3}], [{1}, {1}, {1}, {1}]]
        sage: OMPs(4, alphabet=[2]).list()
        [[{2}, {2}]]
        sage: OMPs(4, alphabet=[1,2]).list()
        [[{1,2}, {1}], [{2}, {2}], [{2}, {1}, {1}], [{1}, {1,2}],
         [{1}, {2}, {1}], [{1}, {1}, {2}], [{1}, {1}, {1}, {1}]]
        sage: OMPs(4, alphabet=4).list() == OMPs(4).list()
        True

    The option ``weight`` can also be used to constrain the ordered multiset partitions.
    It is a refinement of ``alphabet`` in that it specifies the number of times each integer
    appears. In principle, it is a dictionary, but weak compositions are also allowed.
    For example, the ordered multiset partitions of 4 are listed by weight below::

        sage: OrderedMultisetPartitions(4, weight=[0,0,0,1])
        Ordered Multiset Partitions of integer 4 with constraint: weight={4: 1}
        sage: OrderedMultisetPartitions(4, weight=[0,0,0,1]).list()
        [[{4}]]
        sage: OrderedMultisetPartitions(4, weight=[1,0,1]).list()
        [[{1}, {3}], [{3}, {1}], [{1,3}]]
        sage: OrderedMultisetPartitions(4, weight=[0,2]).list()
        [[{2}, {2}]]
        sage: OrderedMultisetPartitions(4, weight=[0,1,1]).list()
        []
        sage: OrderedMultisetPartitions(4, weight=[2,1]).list()
        [[{1}, {1}, {2}], [{1}, {2}, {1}], [{1}, {1,2}], [{2}, {1}, {1}], [{1,2}, {1}]]
        sage: O1 = OrderedMultisetPartitions(weight=[2,0,1])
        sage: O2 = OrderedMultisetPartitions(weight={1:2, 3:1})
        sage: O1 == O2
        True
        sage: OrderedMultisetPartitions(4, weight=[4]).list()
        [[{1}, {1}, {1}, {1}]]

    The ``order`` options are coarse versions of ``weight`` when the multiset is not
    specified. The ordered multiset partitions of integer 4 are listed by order below::

        sage: OrderedMultisetPartitions(4, order=1).list()
        [[{4}]]
        sage: OrderedMultisetPartitions(4, order=2).list()
        [[{1,3}], [{3}, {1}], [{2}, {2}], [{1}, {3}]]
        sage: OrderedMultisetPartitions(4, order=3).list()
        [[{1,2}, {1}], [{2}, {1}, {1}], [{1}, {1,2}], [{1}, {2}, {1}], [{1}, {1}, {2}]]
        sage: OrderedMultisetPartitions(4, order=4).list()
        [[{1}, {1}, {1}, {1}]]
        sage: OrderedMultisetPartitions(4, max_order=2).list()
        [[{4}], [{1,3}], [{3}, {1}], [{2}, {2}], [{1}, {3}]]

    An explicit use of keyword ``order`` is overwritten if it is also implicitly
    provided by the user via the required arguments::

        sage: OrderedMultisetPartitions([1,1,4], order=2).list() # order should be three
        [[{1, 4}, {1}], [{1}, {1, 4}], [{1}, {1}, {4}], [{1}, {4}, {1}], [{4}, {1}, {1}]]
        sage: OrderedMultisetPartitions([1,4], 3, order=2).list() # order should be three
        [[{1,4}, {1}], [{1}, {1,4}], [{1}, {1}, {1}], [{1}, {1}, {4}],
         [{1}, {4}, {1}], [{1}, {4}, {4}], [{4}, {1}, {1}], [{4}, {1}, {4}],
         [{4}, {4}, {1}], [{4}, {4}, {4}]]

    TESTS::

        sage: C = OrderedMultisetPartitions(8, length=3); C.cardinality()
        72
        sage: C == loads(dumps(C))
        True
    """
    @staticmethod
    def __classcall_private__(self, *args, **constraints):
        """
        Return the correct parent based upon the input.

        EXAMPLES::

            sage: OrderedMultisetPartitions()
            Ordered Multiset Partitions
            sage: OrderedMultisetPartitions(4)
            Ordered Multiset Partitions of integer 4
            sage: OrderedMultisetPartitions(4, max_order=2)
            Ordered Multiset Partitions of integer 4 with constraint: max_order=2

            sage: OrderedMultisetPartitions({1:2, 3:1})
            Ordered Multiset Partitions of multiset {{1, 1, 3}}
            sage: OrderedMultisetPartitions({1:2, 3:1}) == OrderedMultisetPartitions([1,1,3])
            True
            sage: OrderedMultisetPartitions({'a':2, 'c':1}, length=2)
            Ordered Multiset Partitions of multiset {{a, a, c}} with constraint: length=2
            sage: OrderedMultisetPartitions({'a':2, 'c':1}, length=4).list()
            []

            sage: OrderedMultisetPartitions(4, 3)
            Ordered Multiset Partitions of order 3 over alphabet {1, 2, 3, 4}
            sage: OrderedMultisetPartitions(['a', 'd'], 3)
            Ordered Multiset Partitions of order 3 over alphabet {a, d}
            sage: OrderedMultisetPartitions([2,4], 3, min_length=2)
            Ordered Multiset Partitions of order 3 over alphabet {2, 4}
             with constraint: min_length=2
        """
        constraints = dict(constraints)
        if "weight" in constraints:
            # Should be a 'dictionary' of letter-frequencies, but accept a weak composition
            w = constraints["weight"]
            if not isinstance(w, dict):
                # make sure we didn't receive ``some_dict.iteritems()``
                if len(w) > 0 and isinstance(w[0], (list, tuple)):
                    w = dict(w)
                else:
                    w = {i+1:w[i] for i in range(len(w)) if w[i] > 0}
            if not all((a in ZZ and a > 0) for a in w.values()):
                raise ValueError("%s must be a dictionary of letter-frequencies or a weak composition"%w)
            else:
                constraints["weight"] = tuple(w.iteritems())

        if "alphabet" in constraints:
            A = constraints["alphabet"]
            if A in ZZ:
                A = range(1,A+1)
            constraints["alphabet"] = Set(A)

        if len(args) == 2: # treat as `alphabet` & `order`
            alph = args[0]; ord = args[1]
            if alph in ZZ:
                alph = range(1,alph+1)
            if (alph and len(Set(alph)) == len(alph)) and (ord in ZZ and ord >= 0):
                constraints.pop("alphabet", None)
                constraints.pop("order", None)
                if constraints == {}:
                    return OrderedMultisetPartitions_A(Set(alph), ord)
                else:
                    return OrderedMultisetPartitions_A_constraints(Set(alph), ord, **constraints)
            elif Set(alph) == Set([]) and ord == 0:
                return OrderedMultisetPartitions_A_constraints(Set(alph), ord, **constraints)
            else:
                raise ValueError("alphabet=%s must be a nonempty set and order=%s must be a nonnegative integer"%(alph,ord))

        elif len(args) == 1: # treat as `size` or `multiset`
            X = args[0]
            if isinstance(X, (list, tuple)):
                tmp = {}
                for i in X:
                    tmp[i] = tmp.get(i, 0) + 1
                X = tmp
            if isinstance(X, dict):
                constraints.pop("size", None)
                constraints.pop("weight", None)
                constraints.pop("alphabet", None)
                constraints.pop("order", None)
                X_items = tuple(X.iteritems())
                if constraints == {}:
                    return OrderedMultisetPartitions_X(X_items)
                else:
                    return OrderedMultisetPartitions_X_constraints(X_items, **constraints)

            elif X in ZZ and X >= 0:
                constraints.pop("size", None)
                if constraints == {}:
                    return OrderedMultisetPartitions_n(X)
                else:
                    return OrderedMultisetPartitions_n_constraints(X, **constraints)

            else:
                raise ValueError("%s must be a nonnegative integer or a list or dictionary representing a multiset"%X)

        else: # generic parent
            return OrderedMultisetPartitions_all_constraints(**constraints)

    def __init__(self, is_finite=None, **constraints):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: OS = OrderedMultisetPartitions(4)
            sage: TestSuite(OS).run()
        """

        constraints = dict(constraints)

        # standardize values for certain keywords
        if "alphabet" in constraints:
            if constraints["alphabet"] in ZZ:
                constraints["alphabet"] = Set(range(1, constraints["alphabet"]+1))
            else:
                constraints["alphabet"] = Set(constraints["alphabet"])

        if "weight" in constraints:
            X = dict(constraints["weight"])
            constraints["weight"] = X
            #constraints["alphabet"] = Set(X.keys())
            #constraints["order"] = sum(X.values())
            #if all(a in ZZ for a in X):
            #    constraints["size"] = sum(a*X[a] for a in X)
            constraints.pop("alphabet",None)
            constraints.pop("min_order",None)
            constraints.pop("order",None)
            constraints.pop("max_order",None)
            constraints.pop("size",None)

        if "length" in constraints:
            constraints.pop("min_length",None)
            constraints.pop("max_length",None)
        if constraints.get("min_length",-1) == constraints.get("max_length",infinity):
            constraints["length"] = constraints.pop("min_length")
            constraints.pop("max_length",None)

        if "order" in constraints:
           constraints.pop("min_order",None)
           constraints.pop("max_order",None)
        if constraints.get("min_order",-1) == constraints.get("max_order",infinity):
            constraints["order"] = constraints.pop("min_order")
            constraints.pop("max_order",None)

        # pop keys with empty values, push numeric values to ``int`` type.
        self.constraints = {}
        for (key,val) in constraints.iteritems():
            if val:
                self.constraints[key] = val

        if is_finite or _is_finite(self.constraints):
            Parent.__init__(self, category=FiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Ordered Multiset Partitions"

    def _constraint_repr_(self, cdict=None):
        if not cdict:
            cdict = self.constraints
        constr = ""
        ss = ['%s=%s'%(key, val) for (key,val) in cdict.iteritems()]
        if len(ss) > 1:
            constr = " with constraints: " + ", ".join(ss)
        elif len(ss) == 1:
            constr = " with constraint: " + ", ".join(ss)
        return constr

    def _element_constructor_(self, lst):
        """
        Construct an element of ``self`` from ``lst``.

        .. TODO:: correctly implement ``check`` here and within OrderedMultisetPartition.

        EXAMPLES::

            sage: P = OrderedMultisetPartitions()
            sage: P([[3],[3,1]]) # indirect doctest
            [{3}, {1,3}]
            sage: P1 = OrderedMultisetPartitions(7, alphabet=3)
            sage: A1 = P1([[3],[3,1]]); A1
            [{3}, {1,3}]
            sage: P2 = OrderedMultisetPartitions(alphabet=3)
            sage: A2 = P2([[3],[3,1]]); A2
            [{3}, {1,3}]
            sage: A1 == A2
            True
            sage: P = OrderedMultisetPartitions(3)
            sage: P([[3],[3,1]])
            Traceback (most recent call last):
            ...
            AssertionError: cannot convert [[3], [3, 1]] into an element of
             Ordered Multiset Partitions of integer 3
        """
        ## construction for user shorthands
        if not lst:
            omp = []
        elif not isinstance(lst[0], (list, tuple, set, frozenset, Set_object)):
            omp = self.from_list(list(lst))

        ## construction for standard input
        else:
            omp = map(list, lst)

        if omp in self:
            return self.element_class(self, map(Set, omp))
        else:
            raise AssertionError("cannot convert %s into an element of %s"%(lst, self))

    Element = OrderedMultisetPartition

    def __contains__(self, x):
        """
        TESTS::

            sage: OMP = OrderedMultisetPartitions
            sage: [[2,1], [1,3]] in OMP()
            True
            sage: [[2,1], [1,3]] in OMP(7)
            True
            sage: [[2,2], [1,3]] in OMP()
            False
            sage: [] in OMP() and [] in OMP(0)
            True
            sage: [] in OMP(2)
            False
            sage: [[2, 1]] in OMP(3, length=2)
            False
            sage: [[2, -1]] in OMP()
            True
        """
        if not isinstance(x, (OrderedMultisetPartition, list, tuple)):
            return False
        else:
            return self._has_valid_blocks(x)

    def _has_valid_blocks(self, x):
        """
        Blocks should be nonempty sets/lists/tuples of distinct elements.
        """
        for block in x:
            if not isinstance(block, (list, tuple, set, Set_object)):
                return False
            if not block or len(Set(block)) != len(block):
                return False
        return True

    def _satisfies_constraints(self, x):
        X = _concatenate(x)
        co = OrderedMultisetPartitions(_get_weight(X))(x)
        def pass_test(co, (key,tst)):
            if key == 'size':
                return co.size() == tst
            if key == 'length':
                return co.length() == tst
            if key == 'min_length':
                return co.length() >= tst
            if key == 'max_length':
                return co.length() <= tst
            if key == 'weight':
                return co.weight() == dict(tst)
            if key == 'alphabet':
                if tst in ZZ:
                    tst = range(1,tst+1)
                return set(tst).issuperset(set(co.letters()))
            if key == 'order':
                return co.order() == tst
            if key == 'min_order':
                return co.order() >= tst
            if key == 'max_order':
                return co.order() <= tst

        return all(pass_test(co, (key,val)) for (key,val) in self.constraints.iteritems() if val)

    def from_list(self, lst):
        """
        Return an ordered multiset partition `c` into singleton blocks, whose singletons
        are the elements `lst`.

        INPUT:

        - ``lst`` -- an iterable

        EXAMPLES::

            sage: OrderedMultisetPartitions().from_list([1,4,8])
            [{1}, {4}, {8}]
            sage: OrderedMultisetPartitions().from_list('abaa')
            [{'a'}, {'b'}, {'a'}, {'a'}]
        """
        if all(a in ZZ for a in lst) and any(a < 0 for a in lst):
            raise ValueError("Something is wrong: `from_list` does not expect to see negative integers; received {}.".format(str(lst)))
        if 0 in list(lst) or '0' in list(lst):
            return self._from_zero_list(lst)
        else:
            d = [Set([x]) for x in lst]
            return self.element_class(self, d)

    def _from_zero_list(self, lst_with_zeros):
        r"""
        Return an ordered multiset partition from a list of nonnegative integers.
        Blocks are separated by zeros. Consecutive zeros are ignored.

        EXAMPLES::

            sage: OrderedMultisetPartitions().from_list([1,2,0,2,0,1,3,4,0,0,4])
            [{1,2}, {2}, {1,3,4}, {4}]
            sage: OrderedMultisetPartitions()._from_zero_list([1,2,4])
            [{1,2,4}]
            sage: OrderedMultisetPartitions()._from_zero_list('abc0a0b')
            [{'a','b','c'}, {'a'}, {'b'}]
        """
        from_zero_lst = list(lst_with_zeros)
        if from_zero_lst[-1] not in {0,'0'}:
            from_zero_lst += [0]
        co = []; block=[]
        for a in from_zero_lst:
            if a in {0,'0'}:
                if block:
                    co.append(block)
                    block = []
            else:
                block.append(a)
        if self._has_valid_blocks(co):
            return self.element_class(self, map(Set, co))
        else:
            raise ValueError("ordered multiset partitions do not have repeated entries within blocks (%s received)"%str(co))

    def _iterator_weight(self, weight):
        """
        Return an iterator for the ordered multiset partitions with weight given by
        the dictionary (or weak composition) ``weight``.

        EXAMPLES::

            sage: OMP = OrderedMultisetPartitions()
            sage: it = OMP._iterator_weight(weight={'a':2, 'b':1})
            sage: list(next(it))
            [[{'a'}, {'a'}, {'b'}], [{'a'}, {'a','b'}], [{'a'}, {'b'}, {'a'}],
             [{'a','b'}, {'a'}], [{'b'}, {'a'}, {'a'}]]
            sage: OMP = OrderedMultisetPartitions()
            sage: it = OMP._iterator_weight(weight=[3,0,1])
            sage: list(next(it))
            [[{1}, {1}, {1}, {3}], [{1}, {1}, {1,3}], [{1}, {1}, {3}, {1}],
             [{1}, {1,3}, {1}], [{1}, {3}, {1}, {1}],
             [{1,3}, {1}, {1}], [{3}, {1}, {1}, {1}]]
            sage: OMPw = OrderedMultisetPartitions(weight=[2,1])
            sage: len(list(OMP._iterator_weight(weight=[2,1]))) == OMPw.cardinality()
            True
        """
        if isinstance(weight, (list, tuple)):
            weight = {k+1: val for k,val in enumerate(weight) if val > 0}
        if isinstance(weight, dict):
            multiset = tuple([k for k in sorted(weight) for _ in range(weight[k])])
        P = OrderedMultisetPartitions(weight)

        if multiset == ():
            yield self.element_class(self, [])
        else:
            # We build ordered multiset partitions of `X` by permutation + deconcatenation
            for alpha in Permutations_mset(multiset):
                co = _break_at_descents(alpha, weak=True)
                for A in P(co).finer(strong=True):
                    yield A

    def an_element(self):
        """
        Return an element of ``self``.

        Rudimentary. Picks the first valid element served up by ``self.__iter__()``.
        """
        try:
            iteration = self.__iter__()
            while True:
                co = next(iteration)
                if co in self:
                    return co
        except StopIteration:
            raise EmptySetError("%s is the empty set"%self)

###############

class OrderedMultisetPartitions_all_constraints(OrderedMultisetPartitions):
    """
    Class of all ordered multiset partitions (with or without constraints).

    TESTS::

        sage: C = OrderedMultisetPartitions(); repr(C)
        'Ordered Multiset Partitions'
        sage: C == loads(dumps(C))
        True
        sage: TestSuite(C).run()

        sage: C = OrderedMultisetPartitions(weight=[2,0,1], length=2); repr(C)
        'Ordered Multiset Partitions with constraints: length=2, weight={1: 2, 3: 1}'
        sage: C == loads(dumps(C))
        True
        sage: TestSuite(C).run()

        sage: D1 = OrderedMultisetPartitions(weight={1:2, 3:1}, min_length=2, max_length=2)
        sage: D2 = OrderedMultisetPartitions({1:2, 3:1}, min_length=2, max_length=2)
        sage: D3 = OrderedMultisetPartitions(5, weight={1:2, 3:1}, length=2)
        sage: D4 = OrderedMultisetPartitions([1,3], 3, weight={1:2, 3:1}, length=2)
        sage: D5 = OrderedMultisetPartitions([1,3], 3, size=5, length=2)
        sage: all(C == eval('D'+str(i)) for i in range(1,6))
        True
        sage: E = OrderedMultisetPartitions({1:2, 3:1}, size=5, min_length=2)
        sage: C == E
        False
    """
    def __init__(self, **constraints):
        """
        Initialize ``self``.
        """
        OrderedMultisetPartitions.__init__(self, None, **constraints)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Ordered Multiset Partitions" + self._constraint_repr_()

    def subset(self, *args):
        """
        Return a subset of all ordered multiset partitions.

        Expects one or two arguments, with different subsets resulting:

        - One Argument:

          * ``X`` -- a dictionary (representing a multiset for `c`),
            or an integer (representing the size of `c`)

        - Two Arguments:

          * ``alph`` -- a list (representing allowable letters within `c`),
            or a positive integer (representing the maximal allowable letter)
          * ``ord``  -- a nonnegative integer (the total number of letters
            within `c`)

        EXAMPLES::

            sage: C = OrderedMultisetPartitions()
            sage: C.subset(3)
            Ordered Multiset Partitions of integer 3
            sage: C.subset(3) == OrderedMultisetPartitions(3)
            True
            sage: C.subset([1,1,4])
            Ordered Multiset Partitions of multiset {{1, 1, 4}}
            sage: C.subset([1,4], 2)
            Ordered Multiset Partitions of order 2 over alphabet {1, 4}
        """
        if not args:
            return self
        return OrderedMultisetPartitions(*args, **self.constraints)

    def __iter__(self):
        """
        Iterate over ordered multiset partitions. If no constraints, then
        iteration is passed to :class:`OrderedMultisetPartitions_n`
        for all sizes `n`.

        EXAMPLES::

            sage: OrderedMultisetPartitions(3).list()
            [[{3}], [{1,2}], [{2}, {1}], [{1}, {2}], [{1}, {1}, {1}]]
            sage: OrderedMultisetPartitions(0).list()
            [[]]
            sage: C = OrderedMultisetPartitions()
            sage: it = C.__iter__()
            sage: [next(it) for i in range(16)]
            [[], [{1}], [{2}], [{1}, {1}], [{3}], [{1,2}], [{2}, {1}],
             [{1}, {2}], [{1}, {1}, {1}], [{4}], [{1,3}], [{3}, {1}],
             [{1,2}, {1}], [{2}, {2}], [{2}, {1}, {1}], [{1}, {3}]]

        TESTS::

            sage: OrderedMultisetPartitions(alphabet=[1,3], max_length=2).list()
            [[], [{1}], [{3}], [{1,3}], [{1}, {1}], [{1}, {3}],
             [{3}, {1}], [{3}, {3}], [{1,3}, {1}], [{1,3}, {3}],
             [{1}, {1,3}], [{3}, {1,3}], [{1,3}, {1,3}]]
            sage: C = OrderedMultisetPartitions(min_length=2, max_order=2)
            sage: it = C.__iter__()
            sage: [next(it) for i in range(15)]
            [[{1}, {1}], [{2}, {1}], [{1}, {2}], [{3}, {1}], [{2}, {2}],
             [{1}, {3}], [{4}, {1}], [{3}, {2}], [{2}, {3}], [{1}, {4}],
             [{5}, {1}], [{4}, {2}], [{3}, {3}], [{2}, {4}], [{1}, {5}]]
        """
        constraints = self.constraints
        if self.category() == FiniteEnumeratedSets():
            if "weight" in constraints:
                iterator = self._iterator_weight(constraints["weight"])
            elif "size" in constraints:
                n = constraints["size"]
                iterator = OrderedMultisetPartitions_n(n).__iter__()
            elif "alphabet" in constraints:
                A = constraints["alphabet"]
                n = len(A)
                if "order" in constraints:
                    ord = constraints["order"]
                    iterator = OrderedMultisetPartitions_A(A, ord).__iter__()
                elif "length" in constraints:
                    k = constraints["length"]
                    iterator = chain(
                        *(OrderedMultisetPartitions_A(A, ord).__iter__(k, k) \
                            for ord in range(k, n*k + 1)))
                else:
                    maxk = constraints.get("max_length", constraints.get("max_order", infinity))
                    mink = constraints.get("min_length", 0)
                    max_ord = min(n * maxk, constraints.get("max_order", infinity))
                    iterator = chain(
                        *(OrderedMultisetPartitions_A(A, ord).__iter__(mink, maxk) \
                            for ord in range(mink, max_ord+1)))
            for co in iterator:
                if self._satisfies_constraints(co):
                    yield self.element_class(self, co)
        else:
            # iterate over partitions of multisets of positive integers
            # or letters over an alphabet
            if "alphabet" in self.constraints:
                A = constraints["alphabet"]
                ell = 0
                while True:
                    for co in OrderedMultisetPartitions_A(A, ell):
                        if self._satisfies_constraints(co):
                            yield self.element_class(self, co)
                    ell += 1
            else:
                n = 0
                while True:
                    for co in OrderedMultisetPartitions_n(n):
                        if self._satisfies_constraints(co):
                            yield self.element_class(self, co)
                    n += 1

###############

class OrderedMultisetPartitions_n(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of a fixed integer `n`.
    """
    def __init__(self, n):
        """
        TESTS::

            sage: C = OrderedMultisetPartitions(Integer(3))
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()
            sage: C2 = OrderedMultisetPartitions(int(3))
            sage: C is C2
            True
            sage: C3 = OrderedMultisetPartitions(7/2)
            Traceback (most recent call last):
            ...
            ValueError:  7/2 must be a nonnegative integer or a list or
             a dictionary representing a multiset
        """
        OrderedMultisetPartitions.__init__(self, True)
        self._n = n

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(OrderedMultisetPartitions(3))
            'Ordered Multiset Partitions of integer 3'
        """
        return "Ordered Multiset Partitions of integer %s" % self._n

    def _has_valid_blocks(self, x):
        """
        Check that blocks of ``x`` are valid.

        Blocks should be nonempty sets/lists/tuples of distinct positive
        integers. Also the sum of all integers should be ``self._n``.
        """
        no_repeats = OrderedMultisetPartitions._has_valid_blocks(self, x)
        nonnegative = all((i in ZZ and i > 0) for block in x for i in block)
        valid_sum = sum(map(sum, x)) == self._n
        return no_repeats and nonnegative and valid_sum

    def cardinality(self):
        """
        Return the number of ordered multiset partitions of `n`.
        """
        # Dispense with the complex computation for small orders.
        orders = {0:1, 1:1, 2:2, 3:5, 4:11, 5:25}
        if self._n <= 5:
            return orders[self._n]

        # We view an ordered multiset partition as a list of 2-regular integer partitions.
        #
        # The 2-regular partitions have a nice generating function (see OEIS:A000009).
        # Below, we take (products of) coefficients of polynomials to compute cardinality.
        t = var('t')
        partspoly = prod([1+t**k for k in range(1,self._n+1)]).coefficients()
        def partspoly_coeff(d): return partspoly[d][0]
        deg = 0
        for alpha in Compositions(self._n):
            deg += prod([partspoly_coeff(d) for d in alpha])
        return deg

    def an_element(self):
        r"""
        Return a typical element of ``OrderedMultisetPartition_n``.
        """
        #output will have at most three blocks, each of size 1, 2, or 3.
        alpha = Compositions(self._n, max_part=self._n//3+1).an_element()
        out = []
        for a in alpha:
            if a in {1, 2, 4}:
                out.append([a])
            else:
                if a % 2:
                    out.append([a//2+1, a//2])
                else:
                    out.append([a//2, a//2-1, 1])
        return self.element_class(self, map(Set, out))

    def random_element(self):
        r"""
        Return a random ``OrderedMultisetPartition_n`` with uniform probability.

        .. TODO::
        Is it really uniform probability?

        This method generates a random composition and then
        then creates new blocks after positions that are not ascents.

        EXAMPLES::

            sage: OrderedMultisetPartitions(5).random_element() # random
            [{1,2}, {1}, {1}]
            sage: OrderedMultisetPartitions(5).random_element() # random
            [{2}, {1,2}]
        """
        C = Compositions(self._n).random_element()
        co = _break_at_descents(C)
        return self.element_class(self, map(Set, co))

    def __iter__(self):
        """
        Iterate over the ordered multiset partitions of integer `n`.

        The degree `n` part of ordered multiset partition contains all sequences of
        subsets of `\NN_+` whose total sum adds up to `n`.

        TESTS::

            sage: OrderedMultisetPartitions(3).list()
            [[{3}], [{1,2}], [{2}, {1}], [{1}, {2}], [{1}, {1}, {1}]]
            sage: OrderedMultisetPartitions(0).list()
            [[]]
        """
        # iteration scheme:
        # - start from an integer composition ``alpha``.
        # - for each ``a`` in ``alpha``, pick distinct integers that sum to ``a``
        for alpha in Compositions(self._n, length=None):
            for p in cartesian_product([RegularPartitions_n(a, 2) for a in alpha]):
                yield self.element_class(self, [Set(list(k)) for k in p])

class OrderedMultisetPartitions_n_constraints(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of a fixed integer `n` satisfying constraints.

    .. TODO:: add ``an_element`` method?
    """
    def __init__(self, n, **constraints):
        """
        Mimic class ``OrderedMultisetPartitions_n`` to initialize.
        """
        OrderedMultisetPartitions.__init__(self, True, **constraints)
        self._n = n

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        cdict = dict(self.constraints)
        cdict.pop("size", None)
        base_repr = "Ordered Multiset Partitions of integer %s" % self._n
        return base_repr + self._constraint_repr_(cdict)

    def __iter__(self):
        """
        Iterate over the ordered multiset partitions of integer ``self._n`` satisfying
        the constraints in ``self.constraints``.
        """
        if "weight" in self.constraints:
            iterator = self._iterator_weight(self.constraints["weight"])
        else:
            iterator = OrderedMultisetPartitions_n(self._n)
        for co in iterator:
            if self._satisfies_constraints(co):
                yield self.element_class(self, co)

    def _has_valid_blocks(self, x):
        """
        Check that blocks of ``x`` are valid and satisfy constraints.
        """
        valid = OrderedMultisetPartitions_n(self._n)._has_valid_blocks(x)
        return valid and self._satisfies_constraints(x)

###############

class OrderedMultisetPartitions_X(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of a fixed multiset `X`.
    """
    def __init__(self, X):
        """
        TESTS::

            sage: C = OrderedMultisetPartitions([1,1,4])
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()
            sage: C2 = OrderedMultisetPartitions({1:2, 4:1})
            sage: C is C2
            True
        """
        OrderedMultisetPartitions.__init__(self, True)
        self._X = X
        self._Xlist = tuple([k for (k,v) in sorted(X) for _ in range(v)])

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(OrderedMultisetPartitions([1,1,4]))
            'Ordered Multiset Partitions of multiset {{1, 1, 4}}'
        """
        ms_rep = "{{" + ", ".join(map(str, self._Xlist)) + "}}"
        return "Ordered Multiset Partitions" + " of multiset %s"%ms_rep

    def _has_valid_blocks(self, x):
        """
        Check that blocks of ``x`` are valid.

        Blocks should be nonempty sets/lists/tuples whose union is the given multiset.
        """
        no_repeats = OrderedMultisetPartitions._has_valid_blocks(self, x)
        valid_partition = _get_multiset(x) == self._Xlist
        return no_repeats and valid_partition

    def cardinality(self):
        """
        Return the number of ordered partitions of multiset ``X``.
        """
        if self._Xlist == ():
            return 0

        # We build ordered multiset partitions of `X` by permutation + deconcatenation
        # Is there a balls-and-boxes formula for this?

        deg = 0
        for alpha in Permutations_mset(self._Xlist):
            fattest = _break_at_descents(alpha)
            deg += prod(2**(len(k)-1) for k in fattest)
        return deg

    def an_element(self):
        r"""
        Return a typical ``OrderedMultisetPartition_X``.
        """
        if not self._Xlist:
            return self.element_class(self, [])
        alpha = Permutations_mset(self._Xlist).an_element()
        co = _break_at_descents(alpha)

        # construct "an element" by breaking the first fat block of `co` in two
        elt = []
        for i in range(len(co)):
            if len(co[i])==1:
                elt.append(co[i])
            else:
                break
        elt.append(co[i][:len(co[i])//2 + 1])
        elt.append(co[i][len(co[i])//2 + 1:])
        elt.extend(co[i+1:])
        return self.element_class(self, map(Set, elt))

    def random_element(self):
        r"""
        Return a random ``OrderedMultisetPartition`` with uniform probability.

        .. TODO::

            Is it really uniform probability?

        This method:

        - generates a random permutation of the multiset, then
        - creates new blocks after positions that are not ascents to
          build ``fat``, then
        - takes a random element of ``fat.finer()``.

        EXAMPLES::

            sage: OrderedMultisetPartitions([1,1,3]).random_element() # random
            [{1}, {1,3}]
            sage: OrderedMultisetPartitions([1,1,3]).random_element() # random
            [{3}, {1}, {1}]
        """
        if not self._Xlist:
            return self.element_class(self, [])

        alpha = Permutations_mset(self._Xlist).random_element()
        co = _break_at_descents(alpha)
        finer = self.element_class(self, map(Set,co)).finer()
        return finer.random_element()

    def __iter__(self):
        """
        Iterate over the ordered multiset partitions of `X`.

        TESTS::

            sage: OrderedMultisetPartitions([1,1,3]).list()
            [[{1}, {1}, {3}], [{1}, {1,3}], [{1}, {3}, {1}], [{1,3}, {1}], [{3}, {1}, {1}]]
            sage: OrderedMultisetPartitions({}).list()
            [[]]
        """
        for A in self._iterator_weight(dict(self._X)):
            yield A

class OrderedMultisetPartitions_X_constraints(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of a fixed multiset `X`
    satisfying constraints.

    .. TODO:: add ``an_element`` method?
    """
    def __init__(self, X, **constraints):
        """
        Mimic class ``OrderedMultisetPartitions_X`` to initialize.
        """
        OrderedMultisetPartitions.__init__(self, True, **constraints)
        self._X = X
        self._Xlist = tuple(k for (k,v) in sorted(X) for _ in range(v))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        cdict = dict(self.constraints)
        cdict.pop("weight", None)
        ms_rep = "{{" + ", ".join(map(str, self._Xlist)) + "}}"
        base_repr = "Ordered Multiset Partitions" + " of multiset %s"%ms_rep
        return base_repr + self._constraint_repr_(cdict)

    def _has_valid_blocks(self, x):
        """
        Check that blocks of ``x`` are valid and satisfy constraints.
        """
        valid = OrderedMultisetPartitions_X(self._X)._has_valid_blocks(x)
        return valid and self._satisfies_constraints(x)

    def __iter__(self):
        """
        Iterate over the ordered multiset partitions of multiset ``self._X`` satisfying
        the constraints in ``self.constraints``.
        """
        for co in OrderedMultisetPartitions_X(self._X):
            if self._satisfies_constraints(co):
                yield self.element_class(self, co)


###############

class OrderedMultisetPartitions_A(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of specified order over a fixed alphabet.
    """
    def __init__(self, A, order):
        """
        TESTS::

            sage: C = OrderedMultisetPartitions(3, 2)
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()
            sage: C2 = OrderedMultisetPartitions([1,2,3], 2)
            sage: C is C2
            True
            sage: list(OrderedMultisetPartitions([1,2,3], 2))
            [[{1,2}], [{1,3}], [{2,3}], [{1}, {1}], [{1}, {2}], [{1}, {3}], [{2}, {1}],
             [{2}, {2}], [{2}, {3}], [{3}, {1}], [{3}, {2}], [{3}, {3}]]
            sage: list(OrderedMultisetPartitions([1,2,3], 2, length=1))
            [[{1,2}], [{1,3}], [{2,3}]]
        """
        OrderedMultisetPartitions.__init__(self, True)
        self._alphabet = A
        self._order = order

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(OrderedMultisetPartitions(3, 2))
            'Ordered Multiset Partitions of order 2 over alphabet {1, 2, 3}'
            sage: repr(OrderedMultisetPartitions([1,3], 2))
            'Ordered Multiset Partitions of order 2 over alphabet {1, 3}'
        """
        A_rep = "Ordered Multiset Partitions of order " + str(self._order)
        A_rep += " over alphabet {%s}"%(", ".join(map(str, sorted(self._alphabet))))
        return A_rep

    def _has_valid_blocks(self, x):
        """
        Check that blocks of ``x`` are valid.

        Blocks should be nonempty sets/lists/tuples of order ``self._order``
        all of whose elements are taken from ``self._alphabet``.
        """
        no_repeats = OrderedMultisetPartitions._has_valid_blocks(self, x)
        valid_order = sum(map(len,x)) == self._order
        valid_letters = self._alphabet.issuperset(Set(_concatenate(x)))
        return no_repeats and valid_order and valid_letters

    def an_element(self):
        r"""
        Return a typical ``OrderedMultisetPartition_A``.
        """
        alpha = Compositions(self._order).an_element()
        co = [Subsets(self._alphabet, a).an_element() for a in alpha]
        return self.element_class(self, map(Set, co))

    def random_element(self):
        r"""
        Return a random ``OrderedMultisetPartition_A`` with uniform probability.

        .. TODO::

            Is it really uniform probability?

        This method:

        - generates a random permutation of the multiset, then
        - creates new blocks after positions that are not ascents
          to build ``fat``, then
        - takes a random element of ``fat.finer()``.

        EXAMPLES::

            sage: OrderedMultisetPartitions([1,4], 3).random_element() # random
            [{4}, {1,4}]
            sage: OrderedMultisetPartitions([1,3], 4).random_element() # random
            [{1,3}, {1}, {3}]
        """
        if not self._alphabet:
            return self.element_class(self, [])

        alpha = Compositions(self._order, max_part=len(self._alphabet)).random_element()
        co = [Subsets(self._alphabet, a).random_element() for a in alpha]
        return self.element_class(self, map(Set, co))

    def __iter__(self, min=None, max=None):
        """
        Iterate over the ordered multiset partitions over alphabet `A`.

        TESTS::

            sage: OrderedMultisetPartitions([1,4], 3).list()
            [[{1,4}, {1}], [{1,4}, {4}], [{1}, {1,4}], [{4}, {1,4}], [{1}, {1}, {1}],
             [{1}, {1}, {4}], [{1}, {4}, {1}], [{1}, {4}, {4}], [{4}, {1}, {1}],
             [{4}, {1}, {4}], [{4}, {4}, {1}], [{4}, {4}, {4}]]
            sage: OrderedMultisetPartitions([], 0).list()
            [[]]
        """
        if self._order == 0:
            yield self.element_class(self,[])
        else:
            # iteration scheme:
            # start from an integer composition ``alpha`` of ``self._order``.
            # for each ``a`` in ``alpha``, pick ``a`` letters from ``alphabet``
            if min:
                min_length = min
            else:
                min_length = self._order // len(self._alphabet)
            if max:
                max_length = max
            else:
                max_length = self._order

            deg = 0
            for k in range(min_length, max_length+1):
                for alpha in Compositions(self._order, length=k,
                                          max_part=len(self._alphabet)):
                    for co in cartesian_product([Subsets(self._alphabet, a) for a in alpha]):
                        yield self.element_class(self, co)

    def cardinality(self):
        """
        Return the number of ordered partitions of order ``self._order`` on
        alphabet ``self._alphabet``.

        TESTS::

            sage: len(OrderedMultisetPartitions([0,42], 3).list())
            12
        """
        if self._order == 0:
            return 0

        # iteration scheme:
        # - start from an integer composition ``alpha`` of ``self._order``.
        # - for each ``a`` in ``alpha``, pick ``a`` letters from ``alphabet``
        min_length = self._order // len(self._alphabet)
        max_length = self._order

        deg = 0
        for k in range(min_length, max_length+1):
            for alpha in Compositions(self._order, length=k, max_part=len(self._alphabet)):
                deg += prod(binomial(len(self._alphabet), a) for a in alpha)
        return deg

class OrderedMultisetPartitions_A_constraints(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of specified order over a fixed alphabet
    satisfying constraints.
    """
    def __init__(self, A, order, **constraints):
        """
        Mimic class ``OrderedMultisetPartitions_A`` to initialize.
        """
        OrderedMultisetPartitions.__init__(self, True, **constraints)
        self._alphabet = A
        self._order = order

    def _repr_(self):
        r"""
        Return a string representation of ``self``.
        """
        cdict = dict(self.constraints)
        cdict.pop("alphabet", None)
        cdict.pop("order", None)
        base_repr = "Ordered Multiset Partitions of order " + str(self._order)
        base_repr += " over alphabet {%s}"%(", ".join(map(str, sorted(self._alphabet))))
        return base_repr + self._constraint_repr_(cdict)

    def _has_valid_blocks(self, x):
        """
        Check that blocks of ``x`` are valid and satisfy constraints.
        """
        valid = OrderedMultisetPartitions_A(self._alphabet, self._order)._has_valid_blocks(x)
        return valid and self._satisfies_constraints(x)

    def an_element(self):
        r"""
        Return a typical ``OrderedMultisetPartition_A``.
        """
        keys = self.constraints.keys()
        n = len(self._alphabet)
        ell = self._order
        if list(keys) == ["length"]:
            kmin = kmax = k = self.constraints["length"]
            if (ell < kmin) or (n * kmax < ell):
                raise EmptySetError("%s is the empty set"%self)
            alpha = Compositions(ell, length=k, max_part=n).an_element()
            co = [Subsets(self._alphabet, a).an_element() for a in alpha]
            #assume ``co`` satisfies all constraints
            return self.element_class(self, map(Set, co))
        else:
            try:
                return next(self.__iter__())
            except StopIteration:
                raise EmptySetError("%s is the empty set"%self)

    def __iter__(self):
        """
        Iterate over the ordered multiset partitions of order ``self._order``over the
        alphabet ``self._alphabet`` satisfying the constraints in ``self.constraints``.

        TESTS::

            sage: len(OrderedMultisetPartitions([0,'a','b'], 3, length=2).list())
            18
            sage: len(OrderedMultisetPartitions([0,'a','b'], 3, max_length=2).list())
            19
            sage: len(OrderedMultisetPartitions([0,'a','b'], 3, length=4).list())
            0
        """
        if "weight" in self.constraints:
            iterator = self._iterator_weight(self.constraints["weight"])
        else:
            if len(self._alphabet) == 0:
                iterator = OrderedMultisetPartitions_n(0).__iter__()
            else:
                minval = self._order // len(self._alphabet)
                minval = self.constraints.get("min_length", minval)
                minval = self.constraints.get("length", minval)

                maxval = self.constraints.get("max_length", self._order)
                maxval = self.constraints.get("length", maxval)

                # FIXME: This has a strong code smell. You probably want to have
                #   a separate iterator function.
                iterator = OrderedMultisetPartitions_A(self._alphabet, self._order).__iter__(minval, maxval)

        for co in iterator:
            if self._satisfies_constraints(co):
                yield self.element_class(self, co)

###############

def _get_multiset(co):
    """
    Construct the multiset (as a sorted tuple) suggested by the lists of lists ``co``.
    """
    return tuple(sorted(_concatenate(co)))

def _get_weight(lst):
    """
    Construct the multiset (as a dictionary) suggested by the multiset-as-list ``lst``.
    """
    out = {}
    for k in lst:
        out[k] = out.get(k,0) + 1
    return out

def _union_of_sets(list_of_sets):
    """
    Return the union of a list of iterables as a Set object.
    """
    return reduce(lambda a,b: Set(a)|Set(b), list_of_sets, Set([]))

def _concatenate(list_of_iters):
    """
    Return the concatenation of a list of iterables as a tuple.
    """
    #if not list_of_iters:
    #    return []
    #return reduce(lambda a,b: a+b, list_of_iters)
    return tuple(_ for block in list_of_iters for _ in block)

def _is_finite(constraints):
    r"""
    Return `True` if the ``constraints`` dictionary describes a finite number
    of ordered multiset partitions.
    """
    if Set(["size", "weight"]).intersection(Set(constraints.keys())):
        return True
    return ("alphabet" in constraints) and \
                    Set(["length", "max_length", "order", "max_order"]).intersection(Set(constraints.keys()))

def _is_initial_segment(lst):
    """
    Return True if ``lst`` is an interval in `ZZ` of the form `[0, 1, \ldots, n]`.
    """
    return list(range(max(lst)+1)) == lst

def _partial_sum(lst):
    """
    Return partial sums of elements in ``lst``.

    EXAMPLES::

        sage: list = [1,3,5]
        sage: partial_sum(list)
        [0, 1, 4, 9]
    """
    result = [0]
    for i in range(len(lst)):
        result.append(result[-1]+lst[i])
    return result

def _descents(w):
    r"""
    Return descent positions in the word ``w``.
    """
    return [j for j in range(len(w)-1) if w[j] > w[j+1]]

def _break_at_descents(alpha, weak=True):
    if not alpha:
        return []

    Blocks = []
    block = [alpha[0]]
    for i in range(1,len(alpha)):
        if (alpha[i-1] > alpha[i]) or (alpha[i-1] == alpha[i] and weak):
            Blocks.append(block)
            block = [alpha[i]]
        else:
            block.append(alpha[i])
    if block:
        Blocks.append(block)
    return Blocks

def _refine_block(S, strong=False):
    r"""
    Return the set of all possible refinements of a set `S`.

    A refinement of `S` is a tuple of nonempty subsets whose union is `S`.

    If optional argument ``strong`` is set to ``True``, then only those
    refinements that are deconcatenations of the list `sorted(S)` are returned.
    """
    X = list(S)
    n = len(X)
    out = []
    if not strong:
        WordSet = IntegerListsLex(min_part=0, max_part=n-1, length=n)
    else:
        WordSet = IntegerListsLex(min_part=0, max_part=n-1, length=n, min_slope=0)

    for w in WordSet:
        if _is_initial_segment(sorted(set(w))):
            a = [set([]) for _ in range(max(w)+1)]
            for pos in range(n):
                a[w[pos]].add(X[pos])
            out.append(a)
    return Set([tuple(map(Set,x)) for x in out])

def _split_block(S, k=2):
    r"""
    Return the set of all possible splittings of a set `S` into `k` parts.

    A splitting of `S` is a tuple of (possibly empty) subsets whose union is `S`.
    """
    X = list(S)
    n = len(X)
    out = []
    WordSet = IntegerListsLex(min_part=0, max_part=k-1, length=n)
    for w in WordSet:
        a = [set([]) for _ in range(k)]
        for pos in range(n):
            a[w[pos]].add(X[pos])
        out.append(a)
    return Set([tuple(map(Set,x)) for x in out])

def _from_tableau(T):
    """
        Return an ordered multiset partition from a sequence of row words
        corresponding to (skew-)tableaux.

        Output is the minimaj bijection `\varphi^{-1}` of [BCHOPSYxx]_ applied to ``T``.

    .. TODO:: implement option for mapping from sequence of actual (skew-)tableaux?

    EXAMPLES::

        sage: mu = OrderedMultisetPartitions(14).an_element(); mu
        [{2,3}, {2,3}, {4}]
        sage: mu.to_tableau()
        [[3, 2], [], [3, 2, 4]]
        sage: _from_tableau(mu.to_tableau()) == mu
        True
    """
    X = _concatenate(T)
    mu = [(i,) for i in T[-1]]
    breaks = [0] + _descents(T[-1]) + [len(mu)-1]
    T = [T[i][::-1] for i in range(len(T)-1)][::-1]
    for f in range(len(breaks)-1):
        for j in range(breaks[f],breaks[f+1]+1):
            mu[j] += tuple(i for i in T[f] if (mu[j][0] < i or j == breaks[f])
                                               and (j == breaks[f+1] or i <= mu[j+1][0]))
    return OrderedMultisetPartitions(X)(mu)

###############

class MinimajCrystal(UniqueRepresentation, Parent):
    r"""
    Crystal of ordered multiset partitions with `ell` letters from alphabet
    `\{1,2,...,n\}` divided into `k` blocks.

    EXAMPLES::

        sage: list(crystals.Minimaj(2,3,2))
        [((2, 1), (1,)), ((2,), (1, 2)), ((1,), (1, 2)), ((1, 2), (2,))]

    TESTS::

        sage: list(crystals.Minimaj(3,2,4)) # more blocks than letters
        []
        sage: list(crystals.Minimaj(3,6,2))
        [((1, 2), (2, 1), (1, 2))]
        sage: list(crystals.Minimaj(2,5,2)) # blocks too fat for alphabet
        []
    """
    def __init__(self, k, ell, n):
        Parent.__init__(self, category=ClassicalCrystals())
        self.n = n
        self.k = k
        self.ell = ell
        self._cartan_type = CartanType(['A',n-1])
        self._full_module = OrderedMultisetPartitions(n, ell, length=k)
        if n == 1:
            self.module_generators = [self(b) for b in self._full_module]
        else:
            self.module_generators = []
            for mu in self._full_module:
                b = self(mu)
                if Set([b.e(i) for i in range(1,n)]) == Set([None]):
                    self.module_generators.append(b)

    def _repr_(self):
        """
        EXAMPLES::

            sage: B = crystals.Minimaj(2,4,3); B
            Minimaj Crystal of type A_3 of words of length 4 into 2 blocks
        """
        return "Minimaj Crystal of type A_%s of words of length %s into %s blocks"%(self.n-1, self.ell, self.k)

    def an_element(self):
        """
        Return a typical element of ``self``.
        """
        return self(OrderedMultisetPartitions(self.n, self.ell, length=self.k).an_element())

    def __contains__(self, x):
        return x in self._full_module

    def from_tableau(self, T):
        r"""
        Return the minimaj bijection `\varphi^{-1}` of [BCHOPSYxx]_ applied to ``T``.

        EXAMPLES::

            sage: B = crystals.Minimaj(3,6,3)
            sage: b = B.an_element(); b
            ((1, 2, 3), (3, 1), (2,))
            sage: t = b.to_tableau(); t
            [[1], [3, 2], [1, 3, 2]]

        TESTS::

            sage: B = crystals.Minimaj(3,6,3)
            sage: all(mu == B.from_tableau(mu.to_tableau()) for mu in B)
            True
            sage: t = B.an_element().to_tableau()
            sage: B1 = crystals.Minimaj(2,6,3)
            sage: B1.from_tableau(t)
            Traceback (most recent call last):
            ...
            AssertionError: [{1,2,3}, {1,3}, {2}] is not an element of Minimaj Crystal
             of type A_2 of words of length 6 into 2 blocks
        """
        mu = _from_tableau(T)
        if mu in self:
            return self(mu)
        else:
            raise AssertionError("%s is not an element of %s"%(mu, self))

    def val(self, q='q'):
        r"""
        Return `Val` polynomial corresponding to ``self``.

        EXAMPLES:

        Verifying Example 4.5 from [BCHOPSYxx]_::

            sage: B = crystals.Minimaj(2, 4, 3) # for `Val_{4,1}^{(3)}`
            sage: B.val()
            (q^2+q+1)*s[2, 1, 1] + q*s[2, 2]
        """
        H = self.module_generators
        Sym = SymmetricFunctions(QQ[q])
        q = Sym.base_ring().gens()[0]
        s = Sym.schur()
        return sum((q**(t.minimaj()) * s[sorted(t.weight().values(), reverse=True)]
                   for t in H), Sym.zero())

    class Element(OrderedMultisetPartition):
        def _repr_(self):
            return repr(self.minimaj_blocks())

        def e(self, i):
            r"""
            Return `e_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.Minimaj(2,3,4)
                sage: b = B.an_element(); b
                ((3,), (2, 3))
                sage: [b.e(i) for i in range(1,4)]
                [((3,), (1, 3)), ((3, 2), (2,)), None]
            """
            t = self.to_tableau()
            B = CrystalOfLetters(['A', self.parent().n-1])
            T = tensor([B]*self.parent().ell)
            blocks = [len(h) for h in t]
            partial = _partial_sum(blocks)
            b = T(*[B(a) for a in flatten(t)])
            if b.e(i) is None:
                return None
            b = b.e(i)
            b = [[ZZ(b[a].value) for a in range(partial[j], partial[j+1])]
                 for j in range(len(partial)-1)]
            return self.parent().from_tableau(b)

        def f(self,i):
            r"""
            Return `f_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.Minimaj(2,3,4)
                sage: b = B.an_element(); b
                ((3,), (2, 3))
                sage: [b.f(i) for i in range(1,4)]
                [None, None, ((4,), (2, 3))]
            """
            t = self.to_tableau()
            B = CrystalOfLetters(['A',self.parent().n-1])
            T = tensor([B]*self.parent().ell)
            blocks = [len(h) for h in t]
            partial = _partial_sum(blocks)
            b = T(*[B(a) for a in flatten(t)])
            if b.f(i) is None:
                return None
            b = b.f(i)
            b = [[ZZ(b[a].value) for a in range(partial[j],partial[j+1])]
                 for j in range(len(partial)-1)]
            return self.parent().from_tableau(b)

