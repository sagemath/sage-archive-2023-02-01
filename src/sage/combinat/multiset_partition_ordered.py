r"""
Ordered Multiset Partitions and the Minimaj Crystal

This module provides element and parent classes for ordered multiset partitions.
It also implements the minimaj crystal of Benkart, et al. (See :class:`MinimajCrystal`.)

AUTHORS:

- Aaron Lauve (2018): initial implementation. First draft of minimaj crystal code
  provided by Anne Schilling.

REFERENCES:

- [BCHOPSY2017]_
- [HRW2015]_
- [HRS2016]_
- [LM2018]_

EXAMPLES:

An ordered multiset partition of the multiset `\{\{1, 3, 3, 5\}\}`::

    sage: OrderedMultisetPartition([[5, 3], [1, 3]])
    [{3,5}, {1,3}]

Ordered multiset partitions of the multiset `\{\{1, 3, 3\}\}`::

    sage: OrderedMultisetPartitions([1,1,3]).list()
    [[{1}, {1}, {3}], [{1}, {1,3}], [{1}, {3}, {1}], [{1,3}, {1}], [{3}, {1}, {1}]]

Ordered multiset partitions of the integer 4::

    sage: OrderedMultisetPartitions(4).list()
    [[{4}], [{1,3}], [{3}, {1}], [{1,2}, {1}], [{2}, {2}], [{2}, {1}, {1}],
     [{1}, {3}], [{1}, {1,2}], [{1}, {2}, {1}], [{1}, {1}, {2}], [{1}, {1}, {1}, {1}]]

Ordered multiset partitions on the alphabet `\{1, 4\}` of order 3::

    sage: OrderedMultisetPartitions([1,4], 3).list()
    [[{1,4}, {1}], [{1,4}, {4}], [{1}, {1,4}], [{4}, {1,4}], [{1}, {1}, {1}],
     [{1}, {1}, {4}], [{1}, {4}, {1}], [{1}, {4}, {4}], [{4}, {1}, {1}],
     [{4}, {1}, {4}], [{4}, {4}, {1}], [{4}, {4}, {4}]]

Crystal of ordered multiset partitions on the alphabet `\{1,2,3\}` with 4 letters
divided into 2 blocks::

    sage: crystals.Minimaj(3, 4, 2).list()
    [((2, 3, 1), (1,)), ((2, 3), (1, 2)), ((2, 3), (1, 3)), ((2, 1), (1, 2)),
     ((3, 1), (1, 2)), ((3, 1, 2), (2,)), ((3, 1), (1, 3)), ((3, 1), (2, 3)),
     ((3, 2), (2, 3)), ((2, 1), (1, 3)), ((2,), (1, 2, 3)), ((3,), (1, 2, 3)),
     ((1,), (1, 2, 3)), ((1, 2), (2, 3)), ((1, 2, 3), (3,))]
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
from sage.structure.element_wrapper import ElementWrapper
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.all import prod
from sage.sets.set import Set, Set_object
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.functions.other import binomial
from sage.calculus.var import var

from sage.combinat.subset import Subsets_sk
from sage.combinat.composition import Composition, Compositions
from sage.combinat.permutation import Permutations_mset
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.shuffle import ShuffleProduct, ShuffleProduct_overlapping
from sage.combinat.crystals.letters import CrystalOfLetters as Letters
from sage.combinat.root_system.cartan_type import CartanType

@add_metaclass(InheritComparisonClasscallMetaclass)
class OrderedMultisetPartition(ClonableArray):
    r"""
    Ordered Multiset Partition

    An *ordered multiset partition* `c` of a multiset `X` is a list
    `[c_1, \ldots, c_r]` of nonempty subsets of `X` (note: not
    sub-multisets), called the *blocks* of `c`, whose multi-union is `X`.

    EXAMPLES:

    The simplest way to create an ordered multiset partition is by specifying
    its blocks as a list, tuple::

        sage: OrderedMultisetPartition([[3],[2,1]])
        [{3}, {1,2}]
        sage: OrderedMultisetPartition(((3,), (1,2)))
        [{3}, {1,2}]
        sage: OrderedMultisetPartition([set([i]) for i in range(2,5)])
        [{2}, {3}, {4}]

    REFERENCES:

    - [HRW2015]_
    - [HRS2016]_
    - [LM2018]_
    """
    @staticmethod
    def __classcall_private__(cls, co):
        """
        Create an ordered multiset partition (i.e., a list of sets) from the passed
        arguments with the appropriate parent.

        EXAMPLES::

            sage: OrderedMultisetPartition([[3], [2,1]])
            [{3}, {1,2}]
            sage: c = OrderedMultisetPartitions()([{2}, {3}, {4}, {5}]); c
            [{2}, {3}, {4}, {5}]
            sage: d = OrderedMultisetPartitions((1,1,1,2,3,5))([{1}, {5, 1, 3}, {2, 1}]); d
            [{1}, {1,3,5}, {1,2}]

        TESTS::

            sage: c.parent() == OrderedMultisetPartitions([2,3,4,5])
            False
            sage: d.parent() == OrderedMultisetPartitions([1,1,1,2,3,5])
            True
            sage: repr(OrderedMultisetPartition([]).parent())
            'Ordered Multiset Partitions of multiset {{}}'
        """
        if not co:
            P = OrderedMultisetPartitions([])
            return P.element_class(P, [])
        else:
            X = _concatenate(co)
            P = OrderedMultisetPartitions(_get_weight(X))
            return P.element_class(P, co)

    def __init__(self, *args):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: c = OrderedMultisetPartitions(7)([[1,3], [1,2]])
            sage: OrderedMultisetPartition([[1,3], [1,2]]) == c
            True
            sage: c.weight()
            {1: 2, 2: 1, 3: 1}

        TESTS::

            sage: OMP = OrderedMultisetPartition
            sage: c0 = OMP([])
            sage: OMP([[]]) == c0
            True
            sage: TestSuite(c0).run()

            sage: d = OMP([[1, 3], [1, 'a', 'b']])
            sage: TestSuite(d).run()

            sage: OMPs = OrderedMultisetPartitions()
            sage: d = OMPs([['a','b','c'],['a','b'],['a']])
            sage: TestSuite(d).run()

            sage: c.size() == 7
            True
            sage: d.size() == None
            True
        """
        if len(args) == 1:
            parent = OrderedMultisetPartitions()
        else:
            parent = args[0]
        # Delte empty blocks
        co = [block for block in args[-1] if block]
        if _has_nonempty_sets(co):
            ClonableArray.__init__(self, parent, [frozenset(list(k)) for k in co])
            self._multiset = _get_multiset(co)
            self._weight = _get_weight(self._multiset)
            self._order = sum([len(block) for block in self])
            if all((a in ZZ and a > 0) for a in self._multiset):
                self._n = ZZ(sum(self._multiset))
            else:
                self._n = None
        else:
            raise ValueError("cannot view %s as an ordered partition of %s"%(co, parent._Xtup))


    def check(self):
        """
        Check that we are a valid ordered multiset partition.

        EXAMPLES::

            sage: c = OrderedMultisetPartitions(4)([[1], [1,2]])
            sage: c.check()

            sage: OMPs = OrderedMultisetPartitions()
            sage: c = OMPs([[1], [1], ['a']])
            sage: c.check()

        TESTS::

            sage: c = OMPs([[1, 1], [1, 4]])
            Traceback (most recent call last):
            ...
            ValueError: cannot convert [[1, 1], [1, 4]] into an element
             of Ordered Multiset Partitions
        """
        if self not in self.parent():
            raise ValueError("{} not an element of {}".format(self, self.parent()))

    def _repr_(self):
        """
        Return a string representation of ``self.``

        EXAMPLES::

            sage: A = OrderedMultisetPartition([[4],[1,2,4],[2,3], [1]])
            sage: A._repr_()
            '[{4}, {1,2,4}, {2,3}, {1}]'
        """
        return self._repr_tight()

    def _repr_normal(self):
        r"""
        Viewing ``self`` as a list `[A_1, \ldots, A_r]` of sets,
        return the standard Sage string representation of `[A_1, \ldots, A_r]`.

        EXAMPLES::

            sage: OrderedMultisetPartition([[4,1,3], [3,2,5]])._repr_normal()
            '[{1, 3, 4}, {2, 3, 5}]'
            sage: OrderedMultisetPartition([[4,1,3,11], [3,'a',5]])._repr_normal()
            "[{1, 11, 3, 4}, {3, 5, 'a'}]"
        """
        # TODO: simplify if/once ``_repr_`` method for ``Set`` sorts its elements.
        if self._n:
            string_parts = map(lambda k: str(sorted(k)), self)
        else:
            string_parts = map(lambda k: str(sorted(k, key=str)), self)
        string_parts = ", ".join(string_parts).replace("[","{").replace("]","}")
        return "[" + string_parts + "]"

    def _repr_tight(self):
        r"""
        Starting from the standard Sage string representation of ``self``
        as a list `[A_1, \ldots, A_r]` of sets, return the shorter string
        gotten by deleting spaces within ``repr(A_i)``.

        EXAMPLES::

            sage: A = OrderedMultisetPartition([[4], [1,2,4], [2,3], [1]])
            sage: A._repr_normal()
            '[{4}, {1, 2, 4}, {2, 3}, {1}]'
            sage: A._repr_tight()
            '[{4}, {1,2,4}, {2,3}, {1}]'
        """
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

        TESTS::

            sage: OMP_n = OrderedMultisetPartitions(4)
            sage: OMP_X = OrderedMultisetPartitions([1,1,2])
            sage: OMP_A = OrderedMultisetPartitions(2, 3)
            sage: mu = [[1], [1, 2]]
            sage: OMP_n(mu) == OMP_X(mu) == OMP_A(mu)
            True
            sage: OMP_n(mu) == mu
            False
            sage: OMP_n(mu) == OMP_n([{1}, {3}])
            False
            sage: OMP_n(mu) == OMP_X([[1], [1,2]])
            True
        """
        if not isinstance(y, OrderedMultisetPartition):
            return False
        return list(self) == list(y)

    def __ne__(self, y):
        """
        Check lack of equality of ``self`` and ``y``.

        The parent is not included as part of the equality check.

        TESTS::

            sage: OMP = OrderedMultisetPartitions(4)
            sage: mu = [[1], [1, 2]]
            sage: OMP(mu).__ne__(mu)
            True
            sage: nu = [[1], [2], [1]]
            sage: OMP(mu).__ne__(OMP(nu))
            True
        """
        return not (self == y)

    def __add__(self, other):
        """
        Return the concatenation of two ordered multiset partitions.

        This operation represents the product in Hopf algebra of ordered multiset
        partitions in its natural basis [LM2018]_.

        EXAMPLES::

            sage: OMP = OrderedMultisetPartition
            sage: OMP([[1],[1],[1,3]]) + OMP([[4,1],[2]])
            [{1}, {1}, {1,3}, {1,4}, {2}]

        TESTS::

            sage: OMP([]) + OMP([]) == OMP([])
            True
            sage: OMP([[1],[1],[1,3]]) + OMP([]) == OMP([[1],[1],[1,3]])
            True
        """
        co = list(self) + list(other)
        X = _concatenate(co)
        return OrderedMultisetPartitions(_get_weight(X))(co)

    @combinatorial_map(order=2, name='reversal')
    def reversal(self):
        r"""
        Return the reverse ordered multiset partition of ``self``.

        The reverse of an ordered multiset partition `(B_1, B_2, \ldots, B_k)`
        is defined as the ordered multiset partition `(B_k, B_{k-1}, \ldots, B_1)`.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([[1], [1, 3], [2, 3, 4]]); C
            [{1}, {1,3}, {2,3,4}]
            sage: C.reversal()
            [{2,3,4}, {1,3}, {1}]
        """
        return self.parent()(list(reversed(self)))

    def shape_from_cardinality(self):
        """
        Return a composition that records the cardinality of each block of ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([[3, 4, 1], [2], [1, 2, 3, 7]]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.shape_from_cardinality()
            [3, 1, 4]
            sage: OrderedMultisetPartition([]).shape_from_cardinality() == Composition([])
            True
        """
        return Composition([len(k) for k in self])

    def shape_from_size(self):
        """
        Return a composition that records the sum of entries of each block of ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([[3, 4, 1], [2], [1, 2, 3, 7]]); C
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
        """
        Return the set of distinct elements occurring within the blocks of ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([[3, 4, 1], [2], [1, 2, 3, 7]]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.letters()
            frozenset({1, 2, 3, 4, 7})
        """
        return _union_of_sets(list(self))

    def multiset(self, as_dict=False):
        """
        Return the multiset corresponding to ``self`` as a tuple or as a dictionary.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([[3, 4, 1], [2], [1, 2, 3, 7]]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.multiset()
            (1, 1, 2, 2, 3, 3, 4, 7)
            sage: C.multiset(as_dict=True)
            {1: 2, 2: 2, 3: 2, 4: 1, 7: 1}
            sage: OrderedMultisetPartition([]).multiset() == ()
            True
        """
        if as_dict:
            return self._weight
        else:
            return self._multiset

    def max_letter(self):
        """
        Return the maximum letter appearing in ``self.letters()``.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([[3, 4, 1], [2], [1, 2, 3, 7]])
            sage: C.max_letter()
            7
            sage: D = OrderedMultisetPartition([['a','b','c'],['a','b'],['a'],['b','c','f'],['c','d']])
            sage: D.max_letter()
            'f'
            sage: C = OrderedMultisetPartition([])
            sage: C.max_letter()
        """
        if not self.letters():
            return None
        else:
            return max(self.letters())

    def size(self):
        """
        Return the size of ``self`` (that is, the sum of all integers in all blocks)
        if ``self`` is a list of subsets of positive integers.

        Else, return ``None``.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([[3, 4, 1], [2], [1, 2, 3, 7]]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.size()
            23
            sage: C.size() == sum(k for k in C.shape_from_size())
            True
            sage: OrderedMultisetPartition([[7,1],[3]]).size()
            11

        TESTS::

            sage: OrderedMultisetPartition([]).size() == 0
            True
            sage: OrderedMultisetPartition([['a','b'],['a','b','c']]).size() is None
            True
        """
        return self._n

    def order(self):
        """
        Return the total number of elements in all blocks.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([[3, 4, 1], [2], [1, 2, 3, 7]]); C
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
        r"""
        Return a dictionary, with keys being the letters in ``self.letters()``
        and values being their (positive) frequency.

        Alternatively, if ``as_weak_comp`` is ``True``, count the number of instances
        `n_i` for each distinct positive integer `i` across all blocks of ``self``.
        Return as a list `[n_1, n_2, n_3, ..., n_k]`, where `k` is the max letter
        appearing in ``self.letters()``.

        EXAMPLES::

            sage: c = OrderedMultisetPartition([[6,1],[1,3],[1,3,6]])
            sage: c.weight()
            {1: 3, 3: 2, 6: 2}
            sage: c.weight(as_weak_comp=True)
            [3, 0, 2, 0, 0, 2]

        TESTS::

            sage: OrderedMultisetPartition([]).weight() == {}
            True

            sage: c = OrderedMultisetPartition([['a','b'],['a','b','c'],['b'],['b'],['c']])
            sage: c.weight()
            {'a': 2, 'b': 4, 'c': 2}
            sage: c.weight(as_weak_comp=True)
            Traceback (most recent call last):
            ...
            ValueError: {'a': 2, 'c': 2, 'b': 4} is not a numeric multiset
        """
        w = self._weight
        if as_weak_comp:
            if all(v in ZZ for v in w):
                w = [w.get(i, 0) for i in range(1, self.max_letter()+1)]
            else:
                raise ValueError("%s is not a numeric multiset"%w)
        return w

    def deconcatenate(self, k=2):
        r"""
        Return the set of `k`-deconcatenations of ``self``.

        A `k`-tuple `(C_1, \ldots, C_k)` of ordered multiset partitions represents
        a `k`-deconcatenation of an ordered multiset partition `C` if
        `C_1 + \cdots + C_k = C`.

        .. NOTE::

            This is not to be confused with ``self.split_blocks()``, which splits each block
            of ``self`` before making `k`-tuples of ordered multiset partitions.

        EXAMPLES::

            sage: sorted(OrderedMultisetPartition([[7,1],[3,4,5]]).deconcatenate())
            [([], [{1,7}, {3,4,5}]), ([{1,7}], [{3,4,5}]), ([{1,7}, {3,4,5}], [])]
            sage: OrderedMultisetPartition([['b','c'],['a']]).deconcatenate()
            {([], [{'b','c'}, {'a'}]), ([{'b','c'}, {'a'}], []), ([{'b','c'}], [{'a'}])}
            sage: OrderedMultisetPartition([['a','b','c']]).deconcatenate(3)
            {([{'a','b','c'}], [], []),
             ([], [{'a','b','c'}], []),
             ([], [], [{'a','b','c'}])}

        TESTS::

            sage: C = OrderedMultisetPartition([['a'],['b'],['c'],['d'],['e']]); C
            [{'a'}, {'b'}, {'c'}, {'d'}, {'e'}]
            sage: all( C.deconcatenate(k).cardinality()
            ....:      == binomial(C.length() + k-1, k-1)
            ....:      for k in range(1, 5) )
            True
        """
        P = OrderedMultisetPartitions(alphabet=self.letters(), max_length=self.length())
        out = []
        for c in IntegerListsLex(self.length(), length=k):
            ps = [sum(c[:i]) for i in range(k+1)]
            out.append(tuple([P(self[ps[i]:ps[i+1]]) for i in range(len(ps)-1)]))
        return Set(out)

    def split_blocks(self, k=2):
        r"""
        Return a dictionary representing the `k`-splittings of ``self``.

        A `k`-tuple `(A^1, \ldots, A^k)` of ordered multiset partitions represents
        a `k`-splitting of an ordered multiset partition `A = [b_1, \ldots, b_r]`
        if one can express each block `b_i` as an (ordered) disjoint union of sets
        `b_i = b^1_i \sqcup \cdots \sqcup b^k_i` (some possibly empty) so that
        each `A^j` is the ordered multiset partition corresponding to the list
        `[b^j_1, b^j_2, \ldots, b^j_r]`, excising empty sets appearing therein.

        This operation represents the coproduct in Hopf algebra of ordered multiset
        partitions in its natural basis [LM2018]_.

        EXAMPLES::

            sage: sorted(OrderedMultisetPartition([[1,2],[3,4]]).split_blocks())
            [([], [{1,2}, {3,4}]), ([{3,4}], [{1,2}]),
             ([{2}, {4}], [{1}, {3}]), ([{2}, {3,4}], [{1}]),
             ([{1}, {4}], [{2}, {3}]), ([{3}], [{1,2}, {4}]),
             ([{2}], [{1}, {3,4}]), ([{1}], [{2}, {3,4}]),
             ([{1}, {3}], [{2}, {4}]), ([{1}, {3,4}], [{2}]),
             ([{2}, {3}], [{1}, {4}]), ([{1,2}], [{3,4}]),
             ([{1,2}, {4}], [{3}]), ([{1,2}, {3,4}], []),
             ([{4}], [{1,2}, {3}]), ([{1,2}, {3}], [{4}])]
            sage: OrderedMultisetPartition([[1,2]]).split_blocks(3)
            {([], [], [{1,2}]): 1, ([], [{1}], [{2}]): 1, ([], [{2}], [{1}]): 1,
             ([], [{1,2}], []): 1, ([{2}], [], [{1}]): 1, ([{1}], [], [{2}]): 1,
             ([{1}], [{2}], []): 1, ([{2}], [{1}], []): 1, ([{1,2}], [], []): 1}
            sage: OrderedMultisetPartition([[4],[4]]).split_blocks()
            {([], [{4}, {4}]): 1, ([{4}], [{4}]): 2, ([{4}, {4}], []): 1}

        TESTS::

            sage: C = OrderedMultisetPartition([[1,2],[4,5,6]]); C
            [{1,2}, {4,5,6}]
            sage: sum(C.split_blocks().values()) == 2**len(C[0]) * 2**len(C[1])
            True
            sage: sum(C.split_blocks(3).values()) == (1+2)**len(C[0]) * (1+2)**len(C[1])
            True
            sage: C = OrderedMultisetPartition([])
            sage: C.split_blocks(3) == {(C, C, C): 1}
            True
        """
        P = OrderedMultisetPartitions(alphabet=self.letters(), max_length=self.length())

        # corner case
        if not self:
            return {tuple([self]*k): 1}
        else:
            out = {}
            tmp = cartesian_product([_split_block(block, k) for block in self])
            for t in tmp:
                tt = tuple([P([l for l in c if len(l)>0]) for c in zip(*t)])
                out[tt] = out.get(tt,0) + 1
            return out

    def finer(self, strong=False):
        """
        Return the set of ordered multiset partitions that are finer than ``self``.

        An ordered multiset partition `A` is finer than another `B` if,
        reading left-to-right, every block of `B` is the union of some consecutive
        blocks of `A`.

        If optional argument ``strong`` is set to ``True``, then return only those
        `A` whose blocks are deconcatenations of blocks of `B`. (Here, we view
        blocks of `B` as sorted lists instead of sets.)

        EXAMPLES::

            sage: C = OrderedMultisetPartition([[3,2]]).finer()
            sage: C.cardinality()
            3
            sage: C.list()
            [[{2}, {3}], [{2,3}], [{3}, {2}]]
            sage: OrderedMultisetPartition([]).finer()
            {[]}
            sage: O = OrderedMultisetPartitions([1, 1, 'a', 'b'])
            sage: o = O([{1}, {'a', 'b'}, {1}])
            sage: o.finer()
            {[{1}, {'b'}, {'a'}, {1}], [{1}, {'a'}, {'b'}, {1}], [{1}, {'a','b'}, {1}]}
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
        """
        Return the ordered multiset partition fatter than ``self``, obtained by
        grouping together consecutive parts according to ``grouping`` (whenever
        this does not violate the strictness condition).

        INPUT:

        - ``grouping`` -- a composition (or list) whose sum is the length of ``self``

        EXAMPLES:

        Let us start with the composition::

            sage: C = OrderedMultisetPartition([[4,1,5], [2], [7,1]]); C
            [{1,4,5}, {2}, {1,7}]

        With ``grouping`` equal to `(1, 1, 1)`, `C` is left unchanged::

            sage: C.fatten([1,1,1])
            [{1,4,5}, {2}, {1,7}]

        With ``grouping`` equal to `(2,1)` or `(1,2)`, a union of consecutive parts
        is achieved::

            sage: C.fatten([2,1])
            [{1,2,4,5}, {1,7}]
            sage: C.fatten([1,2])
            [{1,4,5}, {1,2,7}]

        However, the ``grouping`` `(3)` will throw an error, as `1` cannot appear twice
        in any block of ``C``::

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

        An ordered multiset partition `A` is fatter than another `B` if, reading
        left-to-right, every block of `A` is the union of some consecutive blocks
        of `B`.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([{1,4,5}, {2}, {1,7}]).fatter()
            sage: C.cardinality()
            3
            sage: list(C)
            [[{1,4,5}, {2}, {1,7}], [{1,4,5}, {1,2,7}], [{1,2,4,5}, {1,7}]]
            sage: list(OrderedMultisetPartition([['a','b'],['c'],['a']]).fatter())
            [[{'a','b'}, {'a','c'}], [{'a','b','c'}, {'a'}], [{'a','b'}, {'c'}, {'a'}]]

        Some extreme cases::

            sage: list(OrderedMultisetPartition([['a','b','c']]).fatter())
            [[{'a','b','c'}]]
            sage: list(OrderedMultisetPartition([]).fatter())
            [[]]
            sage: OrderedMultisetPartition([[1], [2], [3], [4]]).fatter().issubset(OrderedMultisetPartition([[1,2,3,4]]).finer())
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
        """
        Return the minimaj statistic on ordered multiset partitions.

        We define `minimaj` via an example:

        1. Sort the block in ``self`` as prescribed by ``self.minimaj_word()``,
           keeping track of the original separation into blocks.
           - in:   [{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}]
           - out:  ( 5,7,1 /  2,4 /  5,6 /  4,6,8 /  3,1 /  1,2,3 )

        2. Record the indices where descents in this word occur.
           - word:      (5, 7, 1 / 2, 4 / 5, 6 / 4, 6, 8 / 3, 1 / 1, 2, 3)
           - indices:    1  2  3   4  5   6  7   8  9 10  11 12  13 14 15
           - descents:  {   2,               7,       10, 11             }

        3. Compute the sum of the descents
           - minimaj = 2 + 7 + 10 + 11 = 30

        REFERENCES:

        - [HRW2015]_

        EXAMPLES::

            sage: C = OrderedMultisetPartition([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}])
            sage: C, C.minimaj_word()
            ([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}],
             (5, 7, 1, 2, 4, 5, 6, 4, 6, 8, 3, 1, 1, 2, 3))
            sage: C.minimaj()
            30
            sage: C = OrderedMultisetPartition([{2,4}, {1,2,3}, {1,6,8}, {2,3}])
            sage: C, C.minimaj_word()
            ([{2,4}, {1,2,3}, {1,6,8}, {2,3}], (2, 4, 1, 2, 3, 6, 8, 1, 2, 3))
            sage: C.minimaj()
            9
            sage: OrderedMultisetPartition([]).minimaj()
            0
            sage: C = OrderedMultisetPartition([['b','d'],['a','b','c'],['b']])
            sage: C, C.minimaj_word()
            ([{'b','d'}, {'a','b','c'}, {'b'}], ('d', 'b', 'c', 'a', 'b', 'b'))
            sage: C.minimaj()
            4
        """
        D = _descents(self.minimaj_word())
        return sum(D) + len(D)

    def minimaj_word(self):
        """
        Return an ordering of ``self._multiset`` derived from the minimaj ordering
        on blocks of ``self``.

        .. SEEALSO::

            :meth:`OrderedMultisetPartition.minimaj_blocks()`.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([[2,1], [1,2,3], [1,2], [3], [1]]); C
            [{1,2}, {1,2,3}, {1,2}, {3}, {1}]
            sage: C.minimaj_blocks()
            ((1, 2), (2, 3, 1), (1, 2), (3,), (1,))
            sage: C.minimaj_word()
            (1, 2, 2, 3, 1, 1, 2, 3, 1)
        """
        return _concatenate(self.minimaj_blocks())

    def minimaj_blocks(self):
        r"""
        Return the minimaj ordering on blocks of ``self``.

        We define the ordering via the example below.

        Sort the blocks `[B_1,...,B_k]` of ``self`` from right to left via:

        1. Sort the last block `B_k` in increasing order, call it the word `W_k`

        2. If blocks `B_{i+1}, \ldots, B_k` have been converted to words
           `W_{i+1}, \ldots, W_k`, use the letters in `B_i` to make the unique
           word `W_i` that has a factorization `W_i=(u,v)` satisfying:

           - letters of `u` and `v` appear in increasing order, with `v` possibly empty;
           - letters in `vu` appear in increasing order;
           - ``v[-1]`` is the largest letter `a \in B_i` satisfying ``a <= W_{i+1}[0]``.

        EXAMPLES::

            sage: OrderedMultisetPartition([[1,5,7], [2,4], [5,6], [4,6,8], [1,3], [1,2,3]])
            [{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}]
            sage: _.minimaj_blocks()
            ((5, 7, 1), (2, 4), (5, 6), (4, 6, 8), (3, 1), (1, 2, 3))
            sage: OrderedMultisetPartition([]).minimaj_blocks()
            ()
        """
        if len(self) == 0:
            return ()

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
        Return a sequence of lists corresponding to row words of (skew-)tableaux.

        OUTPUT:

        The minimaj bijection `\varphi` of Benkart, et al. applied to ``self``.

        .. TODO::

            Implement option for mapping to sequence of (skew-)tableaux?

        REFERENCES:

        - [BCHOPSY2017]_

        EXAMPLES::

            sage: co = ((1,2,4),(4,5),(3,),(4,6,1),(2,3,1),(1,),(2,5))
            sage: OrderedMultisetPartition(co).to_tableau()
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
        Return the major index of ``self``.

        The major index is a statistic on ordered multiset partitions, which
        we define here via an example.

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

        REFERENCES:

        - [HRW2015]_

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

    def shuffle_product(self, other, overlap=False):
        r"""
        Return the shuffles (with multiplicity) of blocks of ``self``
        with blocks of ``other``.

        In case optional argument ``overlap`` is ``True``, instead return the allowable
        overlapping shuffles. An overlapping shuffle `C` is allowable if, whenever one
        of its blocks `c` comes from the union `c = a \cup b` of a block of ``self``
        and a block of ``other``, then this union is disjoint.

        .. SEEALSO::

            :meth:`Composition.shuffle_product()`

        EXAMPLES::

            sage: A = OrderedMultisetPartition([[2,1,3], [1,2]]); A
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
        P = OrderedMultisetPartitions(self._multiset + other._multiset)
        if not overlap:
            for term in ShuffleProduct(self, other, element_constructor=P):
                yield term
        else:
            A = map(tuple, self); B = map(tuple, other)
            for term in ShuffleProduct_overlapping(A, B):
                if len(_concatenate(map(frozenset, term))) == len(P._Xtup):
                    yield P(term)

##############################################################

class OrderedMultisetPartitions(UniqueRepresentation, Parent):
    r"""
    Ordered Multiset Partitions

    An *ordered multiset partition* `c` of a multiset `X` is a list of nonempty subsets
    (not multisets), called the *blocks* of `c`, whose multi-union is `X`.

    The number of blocks of `c` is called its *length*. The *order* of `c` is the
    cardinality of the multiset `X`. If, additionally, `X` is a multiset of positive
    integers, then the *size* of `c` is the sum of all elements of `X`.

    The user may wish to focus on ordered multiset partitions of a given size, or
    over a given alphabet. Hence, this class allows a variety of arguments as input.

    INPUT:

    Expects one or two arguments, with different behaviors resulting:
    - One Argument:
        + `X` -- a dictionary or list or tuple
                 (representing a multiset for `c`),
                 or an integer (representing the size of `c`)
    - Two Arguments:
        + `alph` -- a list (representing allowable letters within blocks of `c`),
                    or a positive integer (representing the maximal allowable letter)
        + `ord`  -- a nonnegative integer (the total number of letters within `c`)

    Optional keyword arguments are as follows:
    (See corresponding methods in see :class:`OrderedMultisetPartition` for more details.)

    - ``weight=X``     (list or dictionary `X`) specifies the multiset for `c`
    - ``size=n``       (integer `n`) specifies the size of `c`
    - ``alphabet=S``   (iterable `S`) specifies allowable elements for the blocks of `c`
    - ``length=k``     (integer `k`) specifies the number of blocks in the partition
    - ``min_length=k`` (integer `k`) specifies minimum number of blocks in the partition
    - ``max_length=k`` (integer `k`) specifies maximum number of blocks in the partition
    - ``order=n``      (integer `n`) specifies the cardinality of the multiset that `c` partitions
    - ``min_order=n``  (integer `n`) specifies minimum number of elements in the partition
    - ``max_order=n``  (integer `n`) specifies maximum number of elements in the partition

    EXAMPLES:

    Passing one argument to ``OrderedMultisetPartitions``:

    There are 5 ordered multiset partitions of the multiset `\{\{1, 1, 4\}\}`::

        sage: OrderedMultisetPartitions([1,1,4]).cardinality()
        5

    Here is the list of them::

        sage: OrderedMultisetPartitions([1,1,4]).list()
        [[{1}, {1}, {4}], [{1}, {1,4}], [{1}, {4}, {1}], [{1,4}, {1}], [{4}, {1}, {1}]]

    By chance, there are also 5 ordered multiset partitions of the integer 3::

        sage: OrderedMultisetPartitions(3).cardinality()
        5

    Here is the list of them::

        sage: OrderedMultisetPartitions(3).list()
        [[{3}], [{1,2}], [{2}, {1}], [{1}, {2}], [{1}, {1}, {1}]]

    Passing two argument to ``OrderedMultisetPartitions``:

    There are also 5 ordered multiset partitions of order 2 over the alphabet `\{1, 4\}`::

        sage: OrderedMultisetPartitions([1, 4], 2)
        Ordered Multiset Partitions of order 2 over alphabet {1, 4}
        sage: OrderedMultisetPartitions([1, 4], 2).cardinality()
        5

    Here is the list of them:

        sage: OrderedMultisetPartitions([1, 4], 2).list()
        [[{1,4}], [{1}, {1}], [{1}, {4}], [{4}, {1}], [{4}, {4}]]

    If no arguments are passed to ``OrderedMultisetPartitions``, then the code returns
    the combinatorial class of all ordered multiset partitions::

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

    Specifying optional arguments:

    - The options ``length``, ``min_length``, and ``max_length`` can be used
      to set length constraints on the ordered multiset partitions. For example, the
      ordered multiset partitions of 4 of length equal to, at least, and at most 2 are
      given by::

        sage: OrderedMultisetPartitions(4, length=2).list()
        [[{3}, {1}], [{1,2}, {1}], [{2}, {2}], [{1}, {3}], [{1}, {1,2}]]
        sage: OrderedMultisetPartitions(4, min_length=3).list()
        [[{2}, {1}, {1}], [{1}, {2}, {1}], [{1}, {1}, {2}], [{1}, {1}, {1}, {1}]]
        sage: OrderedMultisetPartitions(4, max_length=2).list()
        [[{4}], [{1,3}], [{3}, {1}], [{1,2}, {1}], [{2}, {2}], [{1}, {3}], [{1}, {1,2}]]

    - The option ``alphabet`` constrains which integers appear across all blocks of
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

    - The option ``weight`` specifies which multiset `X` is to be considered, if it was
      not passed as one of the required arguments for ``OrderedMultisetPartitions``.
      In principle, it is a dictionary, but weak compositions are also allowed.
      For example, the ordered multiset partitions of 4 are listed by weight below::

        sage: OrderedMultisetPartitions(4, weight=[0,0,0,1])
        Ordered Multiset Partitions of integer 4 with constraint: weight={4: 1}
        sage: OrderedMultisetPartitions(4, weight=[0,0,0,1]).list()
        [[{4}]]
        sage: OrderedMultisetPartitions(4, weight=[1,0,1]).list()
        [[{1,3}], [{1}, {3}], [{3}, {1}]]
        sage: OrderedMultisetPartitions(4, weight=[0,2]).list()
        [[{2}, {2}]]
        sage: OrderedMultisetPartitions(4, weight=[0,1,1]).list()
        []
        sage: OrderedMultisetPartitions(4, weight=[2,1]).list()
        [[{1}, {1,2}], [{1}, {1}, {2}], [{1,2}, {1}], [{1}, {2}, {1}], [{2}, {1}, {1}]]
        sage: O1 = OrderedMultisetPartitions(weight=[2,0,1])
        sage: O2 = OrderedMultisetPartitions(weight={1:2, 3:1})
        sage: O1 == O2
        True
        sage: OrderedMultisetPartitions(4, weight=[4]).list()
        [[{1}, {1}, {1}, {1}]]

      This option is ignored if a multiset `X` is passed as a required argument::

        sage: OrderedMultisetPartitions([1,4], weight={1:3, 2:1}).list()
        [[{1}, {4}], [{1,4}], [{4}, {1}]]

    - The (max/min) ``order`` options place constraints on ordered multiset partitions
      when the multiset `X` is not given as required argument or via the ``weight``
      keyword. The ordered multiset partitions of integer 4 are listed by order below::

        sage: OrderedMultisetPartitions(4, order=1).list()
        [[{4}]]
        sage: OrderedMultisetPartitions(4, order=2).list()
        [[{1,3}], [{3}, {1}], [{2}, {2}], [{1}, {3}]]
        sage: OrderedMultisetPartitions(4, order=3).list()
        [[{1,2}, {1}], [{2}, {1}, {1}], [{1}, {1,2}], [{1}, {2}, {1}], [{1}, {1}, {2}]]
        sage: OrderedMultisetPartitions(4, order=4).list()
        [[{1}, {1}, {1}, {1}]]

      And here is a use of ``max_order``, giving the ordered multiset partitions of
      integer 4 with order 1 or 2::

        sage: OrderedMultisetPartitions(4, max_order=2).list()
        [[{4}], [{1,3}], [{3}, {1}], [{2}, {2}], [{1}, {3}]]

      An explicit use of keyword ``order`` is ignored if the order is also implicitly
      provided by the user via other required (or keyword argument ``weight``). In
      each example below, the order must be 3, so the call ``order=2`` is ignored::

        sage: OrderedMultisetPartitions([1,1,4], order=2).list()
        [[{1}, {1}, {4}], [{1}, {1,4}], [{1}, {4}, {1}], [{1,4}, {1}], [{4}, {1}, {1}]]
        sage: OrderedMultisetPartitions([1,4], 3, order=2).list()
        [[{1,4}, {1}], [{1,4}, {4}], [{1}, {1,4}], [{4}, {1,4}], [{1}, {1}, {1}],
         [{1}, {1}, {4}], [{1}, {4}, {1}], [{1}, {4}, {4}], [{4}, {1}, {1}],
         [{4}, {1}, {4}], [{4}, {4}, {1}], [{4}, {4}, {4}]]
        sage: OrderedMultisetPartitions(6, weight=[2,0,0,1], order=2).list()
        [[{1}, {1}, {4}], [{1}, {1,4}], [{1}, {4}, {1}], [{1,4}, {1}], [{4}, {1}, {1}]]

    TESTS::

        sage: C = OrderedMultisetPartitions(8, length=3); C.cardinality()
        72
        sage: C == loads(dumps(C))
        True
    """
    @staticmethod
    def __classcall_private__(self, *args, **constraints):
        """
        Return the correct parent based upon the input:

            * :class:`OrderedMultisetPartitions_all_constraints`
            * :class:`OrderedMultisetPartitions_X`
            * :class:`OrderedMultisetPartitions_X_constraints`
            * :class:`OrderedMultisetPartitions_n`
            * :class:`OrderedMultisetPartitions_n_constraints`
            * :class:`OrderedMultisetPartitions_A`
            * :class:`OrderedMultisetPartitions_A_constraints`

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
            constraints["alphabet"] = frozenset(A)

        if len(args) == 2: # treat as `alphabet` & `order`
            alph = args[0]; ord = args[1]
            if alph in ZZ:
                alph = range(1,alph+1)
            if (alph and len(Set(alph)) == len(alph)) and (ord in ZZ and ord >= 0):
                constraints.pop("alphabet", None)
                constraints.pop("order", None)
                if constraints == {}:
                    return OrderedMultisetPartitions_A(frozenset(alph), ord)
                else:
                    return OrderedMultisetPartitions_A_constraints(frozenset(alph), ord, **constraints)
            elif frozenset(alph) == frozenset() and ord == 0:
                return OrderedMultisetPartitions_A_constraints(frozenset(alph), ord, **constraints)
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

        .. TODO::

            Fix code to yield outputs indicated below.

        TESTS::

            sage: c = {"length":4, "max_order":6, "alphabet":[2,4,5,6]}
            sage: OrderedMultisetPartitions(**c).constraints
            {'alphabet': frozenset({2, 4, 5, 6}), 'length': 4, 'max_order': 6}
            sage: OrderedMultisetPartitions(17, **c).constraints
            {'alphabet': frozenset({2, 4, 5, 6}), 'length': 4, 'max_order': 6}
            sage: OrderedMultisetPartitions(17, **c).full_constraints
            {'alphabet': frozenset({2, 4, 5, 6}), 'length': 4, 'max_order': 6, 'size': 17}

            sage: c = {"length":4, "min_length":5, "max_order":6, "order":5, "alphabet":4}
            sage: OrderedMultisetPartitions(**c).constraints
            {'alphabet': frozenset({1, 2, 3, 4}), 'length': 4, 'order': 5}
            sage: OrderedMultisetPartitions(4, 5, **c).constraints
            {'length': 4}
            sage: OrderedMultisetPartitions(4, 5, **c).full_constraints
            {'alphabet': frozenset({1, 2, 3, 4}), 'length': 4, 'order': 5}

            sage: c = {"weight":[2,2,0,3], "min_length":5, "max_order":6, "order":5, "alphabet":4}
            sage: OrderedMultisetPartitions(**c).constraints
            {'min_length': 5, 'weight': {1: 2, 2: 2, 4: 3}}
            sage: OrderedMultisetPartitions([1,1,2,2,4,4,4], **c).constraints
            {'min_length': 5}
            sage: OrderedMultisetPartitions([1,1,2,2,4,4,4], **c).full_constraints
            {'min_length': 5, 'weight': {1: 2, 2: 2, 4: 3}}
        """
        constraints = dict(constraints)

        # standardize values for certain keywords
        if "alphabet" in constraints:
            if constraints["alphabet"] in ZZ:
                constraints["alphabet"] = frozenset(range(1, constraints["alphabet"]+1))
            else:
                constraints["alphabet"] = frozenset(constraints["alphabet"])

        if "weight" in constraints:
            X = dict(constraints["weight"])
            constraints["weight"] = X
            constraints.pop("alphabet",None)
            constraints.pop("min_order",None)
            constraints.pop("order",None)
            constraints.pop("max_order",None)
            constraints.pop("size",None)

        if "length" in constraints:
            constraints.pop("min_length",None)
            constraints.pop("max_length",None)
        min_k = constraints.get("min_length",0)
        max_k = constraints.get("max_length",infinity)
        assert min_k <= max_k, "min_length=%s <= max_length=%s"%(min_k, max_k)
        if min_k == max_k:
            constraints["length"] = \
                constraints.pop("min_length", constraints.pop("max_length"))

        if "order" in constraints:
           constraints.pop("min_order",None)
           constraints.pop("max_order",None)
        min_ord = constraints.get("min_order",0)
        max_ord = constraints.get("max_order",infinity)
        assert min_ord <= max_ord, "min_order=%s <= max_order=%s"%(min_ord, max_ord)
        if min_ord == max_ord:
            constraints["order"] = \
                constraints.pop("min_order", constraints.pop("max_order"))

        # pop keys with empty values, with the exception of 'size' or 'order'
        self.constraints = {}
        for (key,val) in constraints.iteritems():
            if val:
                self.constraints[key] = val
            elif key in ("size", "order", "length") and val is not None:
                self.constraints[key] = val

        self.full_constraints = dict(self.constraints)
        if hasattr(self, "_X"):
            self.full_constraints["weight"] = dict(self._X)
            self.constraints.pop("weight", None)
        if hasattr(self, "_n"):
            self.full_constraints["size"] = self._n
            self.constraints.pop("size", None)
        if hasattr(self, "_alphabet"):
            self.full_constraints["alphabet"] = self._alphabet
            self.constraints.pop("alphabet", None)
            self.full_constraints["order"] = self._order
            self.constraints.pop("order", None)

        if is_finite or _is_finite(constraints):
            Parent.__init__(self, category=FiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: OrderedMultisetPartitions()._repr_()
            'Ordered Multiset Partitions'
        """
        return "Ordered Multiset Partitions"

    def _constraint_repr_(self, cdict=None):
        """
        Return a string representation of all constraints
        appearing within ``self.constraints``.

        A helper method for ``self._repr_()``.

        EXAMPLES::

            sage: OMPs = OrderedMultisetPartitions()
            sage: c = {"length":4, "max_order":6, "alphabet":frozenset([2,4,5,6])}
            sage: OMPs._constraint_repr_(c)
            ' with constraints: alphabet={2, 4, 5, 6}, length=4, max_order=6'
            sage: c = {"size":14}
            sage: OMPs._constraint_repr_(c)
            ' with constraint: size=14'
        """
        if not cdict:
            cdict = dict(self.constraints)
        if "alphabet" in cdict:
            cdict["alphabet"] = Set(cdict["alphabet"])
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

        EXAMPLES::

            sage: P = OrderedMultisetPartitions()
            sage: A = P([[3],[3,1]]) ; A # indirect doctest
            [{3}, {1,3}]
            sage: P1 = OrderedMultisetPartitions(7, alphabet=3)
            sage: A1 = P1([[3],[3,1]]); A1
            [{3}, {1,3}]
            sage: P2 = OrderedMultisetPartitions(alphabet=3)
            sage: A2 = P2([[3],[3,1]]); A2
            [{3}, {1,3}]
            sage: A == A1 == A2
            True
            sage: P = OrderedMultisetPartitions(3)
            sage: P([[3],[3,1]])
            Traceback (most recent call last):
            ...
            ValueError: cannot convert [[3], [3, 1]] into an element of
             Ordered Multiset Partitions of integer 3
        """
        if not lst:
            omp = []
        else:
            omp = map(list, lst)

        if omp in self:
            return self.element_class(self, map(frozenset, omp))
        else:
            raise ValueError("cannot convert %s into an element of %s"%(lst, self))

    Element = OrderedMultisetPartition

    def __contains__(self, x):
        """
        TESTS::

            sage: [[2,1], [1,3]] in OrderedMultisetPartitions()
            True
            sage: [[2,1], [1,3]] in OrderedMultisetPartitions(7)
            True
            sage: [[2,2], [1,3]] in OrderedMultisetPartitions()
            False
            sage: [] in OrderedMultisetPartitions()
            True
            sage: [] in OrderedMultisetPartitions(0)
            True
            sage: [] in OrderedMultisetPartitions(2)
            False
            sage: [[2, 1]] in OrderedMultisetPartitions(3, length=2)
            False
            sage: [[2, -1]] in OrderedMultisetPartitions()
            True
        """
        if not isinstance(x, (OrderedMultisetPartition, list, tuple)):
            return False
        else:
            return _has_nonempty_sets(x) and self._satisfies_constraints(x)

    def _satisfies_constraints(self, x):
        """
        Check whether or not ``x`` satisfies all of the constraints
        appearing within ``self.full_constraints`` (Boolean output).

        NOTES::

            This test will cause an infinite recursion with
            ``self._element_constructor()`` if the ``__contains__``
            method in ``OrderedMultisetPartitions_X`` is removed.

        TESTS::

            sage: c = {"length":3, "max_order":5, "alphabet":[1,2,4], "size":12}
            sage: OMPs = OrderedMultisetPartitions(**c)
            sage: OMPs._satisfies_constraints([{2,4}, {1}, {1,4}])
            True
            sage: failures = {((2,4), (2,4)), ((1,2,4), (1,), (1,4)), \
                              ((2,4), (3,), (3,)), ((2,4), (1,), (2,4))}
            sage: any(OMPs._satisfies_constraints(x) for x in failures)
            False
            sage: c = {"max_length":4, "order":6, "weight":{1:2, 2:1, 4:2}}
            sage: OMPs = OrderedMultisetPartitions(**c) # order constraint is ignored
            sage: OMPs._satisfies_constraints([{2,4}, {1}, {1,4}])
            True
            sage: failures = {((2,), (4,), (1,), (1,), (4,)), ((1,), (1,), (2,4), (2,4))}
            sage: any(OMPs._satisfies_constraints(x) for x in failures)
            False
        """
        X = _concatenate(x)
        P = OrderedMultisetPartitions_X(tuple(_get_weight(X).iteritems()))
        x = P.element_class(P, [frozenset(block) for block in x])
        def pass_test(co, key, tst):
            # define simple tests for each possible constraint
            if key == 'size':
                return co.size() == tst
            if key == 'length':
                return co.length() == tst
            if key == 'min_length':
                return co.length() >= tst
            if key == 'max_length':
                return co.length() <= tst
            if key == 'weight':
                return co.weight() == tst
            if key == 'alphabet':
                return frozenset(co.letters()).issubset(tst)
            if key == 'order':
                return co.order() == tst
            if key == 'min_order':
                return co.order() >= tst
            if key == 'max_order':
                return co.order() <= tst

        return all(pass_test(x, key, tst) for (key, tst) in self.full_constraints.iteritems() if tst)

    def _from_list(self, lst):
        """
        Return an ordered multiset partition of singleton blocks, whose singletons
        are the elements ``lst``.

        INPUT:

        - ``lst`` -- an iterable

        EXAMPLES::

            sage: OMPs = OrderedMultisetPartitions()
            sage: OMPs._from_list([1,4,0,8])
            [{1,4}, {8}]
            sage: OMPs._from_list([1,4,8])
            [{1}, {4}, {8}]
            sage: OMPs._from_list([1,4,8,0]) == OrderedMultisetPartition([[1,4,8]])
            True
            sage: OMPs._from_list('abaa')
            [{'a'}, {'b'}, {'a'}, {'a'}]
            sage: OMPs._from_list('ab0a0a')
            [{'a','b'}, {'a'}, {'a'}]

        TESTS::

            sage: OMPs._from_list([1,0,2,3,1]) == OrderedMultisetPartition([[1], [2,3,1]])
            True
            sage: OMPs._from_list([1,2,'3',0,1]) == OrderedMultisetPartition([{1,2,'3'}, [1]])
            True
        """
        if all(a in ZZ for a in lst) and any(a < 0 for a in lst):
            raise ValueError("Something is wrong: `_from_list` does not expect to see negative integers; received {}.".format(str(lst)))
        if 0 in list(lst) or '0' in list(lst):
            return self._from_list_with_zeros(lst)
        else:
            d = [frozenset([x]) for x in lst]
            c = self.element_class(self, d)
            # give a better parent, if `self` is generic
            if isinstance(self, OrderedMultisetPartitions_all_constraints):
                P = OrderedMultisetPartitions(_get_weight(lst))
                return P.element_class(P, c)
            else:
                return self.element_class(self, c)

    def _from_list_with_zeros(self, lst_with_zeros):
        r"""
        Return an ordered multiset partition from a list of nonnegative integers.
        Blocks are separated by zeros. Consecutive zeros are ignored.

        EXAMPLES::

            sage: OrderedMultisetPartitions()._from_list([1,2,4])
            [{1}, {2}, {4}]
            sage: OrderedMultisetPartitions()._from_list_with_zeros([1,2,4])
            [{1,2,4}]
            sage: OrderedMultisetPartitions()._from_list_with_zeros([1,0,2,0,0,4])
            [{1}, {2}, {4}]
            sage: OrderedMultisetPartitions()._from_list_with_zeros('abc00a0b')
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
        if co in self:
            c = self.element_class(self, map(frozenset, co))
            # give a better parent, if `self` is generic
            if isinstance(self, OrderedMultisetPartitions_all_constraints):
                P = OrderedMultisetPartitions(c.weight())
                return P.element_class(P, c)
            else:
                return c
        else:
            raise ValueError("ordered multiset partitions do not have repeated entries within blocks (%s received)"%str(co))

    def an_element(self):
        """
        Return an element of ``self``.

        Rudimentary. Picks the first valid element served up by ``self.__iter__()``.

        EXAMPLES::

            sage: OMP = OrderedMultisetPartitions
            sage: OMP().an_element()
            []
            sage: OMP(length=4).an_element()
            [{1}, {1}, {1}, {1}]
            sage: OMP(length=4, min_order=6).an_element()
            [{1,2}, {1,2}, {1}, {1}]
            sage: OMP(length=4, min_order=6, alphabet=[1,2,'a']).an_element()
            [{1,2,'a'}, {'a'}, {'a'}, {'a'}]
        """
        try:
            iteration = self.__iter__()
            while True:
                co = next(iteration)
                if co in self:
                    return co
        except StopIteration:
            raise EmptySetError("%s is the empty set"%self)

    def __iter__(self):
        """
        Iterate over ordered multiset partitions.

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
            sage: OrderedMultisetPartitions(alphabet=[1,3], min_length=2).list()
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot list an infinite set
        """
        # Look for evidence of ``FiniteEnumeratedSets()`` among constraints.
        # ``_base_iterator`` ignores most constraints in ``self.full_constraints``.
        iterator = _base_iterator(self.full_constraints)
        if iterator:
            for co in iterator:
                if self._satisfies_constraints(co):
                    yield self.element_class(self, co)
        else:
            # iterate over blocks of letters over an alphabet
            if "alphabet" in self.constraints:
                A = self.constraints["alphabet"]
                # establish a cutoff order `max_ell`
                max = self.constraints.get("max_length", infinity)
                max = self.constraints.get("length", max)
                max = max * len(A)
                max = self.constraints.get("max_order", max)
                max_ell = self.constraints.get("order", max)
                ell = 0
                while True and ell <= max_ell:
                    for co in _iterator_order(A, ell):
                        if self._satisfies_constraints(co):
                            yield self.element_class(self, co)
                    ell += 1
            # or iterate over partitions of multisets of positive integers
            else:
                n = 0
                while True:
                    for co in _iterator_size(n):
                        if self._satisfies_constraints(co):
                            yield self.element_class(self, co)
                    n += 1

###############

class OrderedMultisetPartitions_all_constraints(OrderedMultisetPartitions):
    """
    Class of all ordered multiset partitions (with or without constraints).

    EXAMPLES::

        sage: C = OrderedMultisetPartitions(); C
        Ordered Multiset Partitions
        sage: [[1],[1,'a']] in C
        True

        sage: OrderedMultisetPartitions(weight=[2,0,1], length=2)
        Ordered Multiset Partitions with constraints: length=2, weight={1: 2, 3: 1}
    """
    def __init__(self, **constraints):
        """
        Initialize ``self``.

        TESTS::

            sage: C = OrderedMultisetPartitions()
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()

            sage: C = OrderedMultisetPartitions(weight=[2,0,1], length=2)
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()

            sage: D1 = OrderedMultisetPartitions(weight={1:2, 3:1}, min_length=2, max_length=2)
            sage: D2 = OrderedMultisetPartitions({1:2, 3:1}, min_length=2, max_length=2)
            sage: D3 = OrderedMultisetPartitions(5, weight={1:2, 3:1}, length=2)
            sage: D4 = OrderedMultisetPartitions([1,3], 3, weight={1:2, 3:1}, length=2)
            sage: D5 = OrderedMultisetPartitions([1,3], 3, size=5, length=2)
            sage: any(C == eval('D'+str(i)) for i in range(1,6))
            False
            sage: all(Set(C) == Set(eval('D'+str(i))) for i in range(1,6))
            True
            sage: E = OrderedMultisetPartitions({1:2, 3:1}, size=5, min_length=2)
            sage: Set(C) == Set(E)
            False
        """
        OrderedMultisetPartitions.__init__(self, None, **constraints)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: OrderedMultisetPartitions(min_length=3, max_order=5)._repr_()
            'Ordered Multiset Partitions with constraints: min_length=3, max_order=5'
            sage: OrderedMultisetPartitions(min_length=3, max_order=5, alphabet=[1,'a'])._repr_()
            "Ordered Multiset Partitions with constraints: alphabet={'a', 1}, max_order=5, min_length=3"
        """
        return "Ordered Multiset Partitions" + self._constraint_repr_()

    def subset(self, size):
        """
        Return a subset of all ordered multiset partitions.

        INPUT:

        - ``size`` -- an integer representing a slice of all ordered
                      multiset partitions.

        The slice alluded to above is taken with respect to length, or
        to order or to size, depending on the constraints of  ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartitions(weight={2:2, 3:1, 5:1})
            sage: C.subset(3)
            Ordered Multiset Partitions of multiset {{2, 2, 3, 5}} with constraint: length=3
            sage: C = OrderedMultisetPartitions(alphabet=[2,3,5])
            sage: C.subset(3)
            Ordered Multiset Partitions of order 3 over alphabet {2, 3, 5}
            sage: C = OrderedMultisetPartitions()
            sage: C.subset(3)
            Ordered Multiset Partitions of integer 3
            sage: C.subset(3) == OrderedMultisetPartitions(3)
            True
        """
        fc = self.full_constraints
        if "weight" in fc:
            return OrderedMultisetPartitions(fc["weight"], length=size, **self.constraints)
        elif "alphabet" in fc:
            return OrderedMultisetPartitions(fc["alphabet"], size, **self.constraints)
        else:
            return OrderedMultisetPartitions(size, **self.constraints)

###############

class OrderedMultisetPartitions_n(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of a fixed integer `n`.
    """
    def __init__(self, n):
        """

        TESTS::

            sage: C = OrderedMultisetPartitions(Integer(4))
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()
            sage: C2 = OrderedMultisetPartitions(int(4))
            sage: C is C2
            True
            sage: C3 = OrderedMultisetPartitions(7/2)
            Traceback (most recent call last):
            ...
            ValueError:  7/2 must be a nonnegative integer or a list or
             dictionary representing a multiset
        """
        self._n = n
        OrderedMultisetPartitions.__init__(self, True)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(OrderedMultisetPartitions(3))
            'Ordered Multiset Partitions of integer 3'
        """
        return "Ordered Multiset Partitions of integer %s" % self._n

    def cardinality(self):
        """
        Return the number of ordered multiset partitions of integer ``self._n``.

        TESTS::

            sage: len(OrderedMultisetPartitions(10).list())
            1500
            sage: OrderedMultisetPartitions(10).cardinality()
            1500
        """
        # Dispense with the complex computation for small orders.
        orders = {0:1, 1:1, 2:2, 3:5, 4:11, 5:25}
        if self._n <= 5:
            return ZZ(orders[self._n])

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
        return ZZ(deg)

    def an_element(self):
        """
        Return a typical element of ``OrderedMultisetPartition_n``.

        EXAMPLES::

            sage: OrderedMultisetPartitions(13).an_element()
            [{2,3}, {2,3}, {1,2}]
            sage: OrderedMultisetPartitions(14).an_element()
            [{2,3}, {2,3}, {4}]
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
        return self.element_class(self, map(frozenset, out))

    def random_element(self):
        """
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
        return self.element_class(self, map(frozenset, co))

    def __iter__(self):
        """
        Iterate over ``self``.

        TESTS::

            sage: O = OrderedMultisetPartitions(6)
            sage: it = O.__iter__()
            sage: [next(it) for _ in range(10)]
            [[{6}], [{2,4}], [{1,5}], [{1,2,3}],
             [{5}, {1}], [{2,3}, {1}], [{1,4}, {1}],
             [{4}, {2}], [{1,3}, {2}], [{4}, {1}, {1}]]
        """
        return _iterator_size(self._n)

class OrderedMultisetPartitions_n_constraints(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of a fixed integer `n` satisfying constraints.
    """
    def __init__(self, n, **constraints):
        """
        Mimic class ``OrderedMultisetPartitions_n`` to initialize.

        TESTS::

            sage: C = OrderedMultisetPartitions(6, length=3)
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()

            sage: C = OrderedMultisetPartitions(6, weight=[3,0,1], length=3)
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()
        """
        self._n = n
        OrderedMultisetPartitions.__init__(self, True, size=n, **constraints)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = OrderedMultisetPartitions(14, length=4, max_order=6, alphabet={2,4,5,6})
            sage: O._repr_()
            'Ordered Multiset Partitions of integer 14 with constraints:
             alphabet={2, 4, 5, 6}, length=4, max_order=6'
        """
        cdict = dict(self.constraints)
        cdict.pop("size", None)
        base_repr = "Ordered Multiset Partitions of integer %s" % self._n
        return base_repr + self._constraint_repr_(cdict)

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
        self._X = X
        # sort the multiset
        if all((k in ZZ and k > 0) for (k,v) in X):
            self._Xtup = tuple([k for (k,v) in sorted(X) for _ in range(v)])
        else:
            self._Xtup = tuple([k for (k,v) in sorted(X, key=str) for _ in range(v)])
        OrderedMultisetPartitions.__init__(self, True)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(OrderedMultisetPartitions([1,1,4]))
            'Ordered Multiset Partitions of multiset {{1, 1, 4}}'
        """
        ms_rep = "{{" + ", ".join(map(str, self._Xtup)) + "}}"
        return "Ordered Multiset Partitions" + " of multiset %s"%ms_rep

    def __contains__(self, x):
        """
        TESTS::

            sage: from sage.combinat.multiset_partition_ordered import OrderedMultisetPartitions_X as OMPX
            sage: [[2,1], [1,3]] in OMPX(((1,2), (2,1), (3,1)))
            True
            sage: co = OrderedMultisetPartition([[2,1], [1,3]])
            sage: co in OMPX(((1,2), (2,1), (3,1)))
            True
            sage: [[2,1], [2,3]] in OMPX(((1,2), (2,1), (3,1)))
            False
            sage: [] in OMPX(())
            True
            sage: [[2, -1], [2,'a']] in OMPX(((2,2), (-1,1), ('a',1)))
            True
        """
        if not isinstance(x, (OrderedMultisetPartition, list, tuple)):
            return False
        else:
            x_Xtup = sorted(_concatenate(x), key=str)
            self_Xtup = sorted(self._Xtup, key=str)
            return _has_nonempty_sets(x) and x_Xtup == self_Xtup

    def cardinality(self):
        """
        Return the number of ordered partitions of multiset ``X``.

        TESTS::

            sage: len(OrderedMultisetPartitions([2,2,2,3,4,5]).list())
            535
            sage: OrderedMultisetPartitions([2,2,2,3,4,5]).cardinality()
            535
        """
        if self._Xtup == ():
            return ZZ(0)

        # We build ordered multiset partitions of `X` by permutation + deconcatenation
        # Is there a balls-and-boxes formula for this?

        deg = 0
        for alpha in Permutations_mset(self._Xtup):
            fattest = _break_at_descents(alpha)
            deg += prod(2**(len(k)-1) for k in fattest)
        return ZZ(deg)

    def an_element(self):
        """
        Return a typical ``OrderedMultisetPartition_X``.

        EXAMPLES::

            sage: OrderedMultisetPartitions([2,2,2,3,4,5]).an_element()
            [{2}, {2}, {2,3,4}, {5}]
            sage: OrderedMultisetPartitions([2,2,2,3,4,4,5]).an_element()
            [{2}, {2}, {2,3}, {4}, {4,5}]
        """
        if not self._Xtup:
            return self.element_class(self, [])
        alpha = Permutations_mset(self._Xtup).an_element()
        co = _break_at_descents(alpha)

        # construct "an element" by breaking the first fat block of `co` in two
        elt = []
        for i in range(len(co)):
            if len(co[i]) == 1:
                elt.append(co[i])
            else:
                break
        elt.append(co[i][:len(co[i])//2 + 1])
        elt.append(co[i][len(co[i])//2 + 1:])
        elt.extend(co[i+1:])
        return self.element_class(self, map(frozenset, elt))

    def random_element(self):
        """
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
        if not self._Xtup:
            return self.element_class(self, [])

        alpha = Permutations_mset(self._Xtup).random_element()
        co = _break_at_descents(alpha)
        finer = self.element_class(self, map(frozenset,co)).finer()
        return finer.random_element()

    def __iter__(self):
        """
        Iterate over ``self``.

        TESTS::

            sage: O = OrderedMultisetPartitions([1, 1, 'a'])
            sage: it = O.__iter__()
            sage: sorted([next(it) for _ in range(O.cardinality())], key=str)
            [[{'a'}, {1}, {1}],
             [{1,'a'}, {1}],
             [{1}, {'a'}, {1}],
             [{1}, {1,'a'}],
             [{1}, {1}, {'a'}]]
            sage: O = OrderedMultisetPartitions([1, 1, 2])
            sage: it = O.__iter__()
            sage: [next(it) for _ in range(O.cardinality())]
            [[{1}, {1,2}], [{1}, {1}, {2}], [{1,2}, {1}],
             [{1}, {2}, {1}], [{2}, {1}, {1}]]
        """
        return _iterator_weight(weight=dict(self._X))

class OrderedMultisetPartitions_X_constraints(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of a fixed multiset `X`
    satisfying constraints.
    """
    def __init__(self, X, **constraints):
        """
        Mimic class ``OrderedMultisetPartitions_X`` to initialize.

        TESTS::

            sage: C = OrderedMultisetPartitions([1,1,2,4], length=3)
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()

            sage: C = OrderedMultisetPartitions([1,1,2,4], weight=[3,0,1], max_length=3) # weight ignored
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()
        """
        self._X = X
        self._Xtup = tuple(k for (k,v) in sorted(X) for _ in range(v))
        OrderedMultisetPartitions.__init__(self, True, weight=X, **constraints)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = OrderedMultisetPartitions([2,2,2,3,4,4,5], min_length=4, max_length=5)
            sage: O._repr_()
            'Ordered Multiset Partitions of multiset {{2, 2, 2, 3, 4, 4, 5}}
             with constraints: min_length=4, max_length=5'
        """
        cdict = dict(self.constraints)
        cdict.pop("weight", None)
        ms_rep = "{{" + ", ".join(map(str, self._Xtup)) + "}}"
        base_repr = "Ordered Multiset Partitions" + " of multiset %s"%ms_rep
        return base_repr + self._constraint_repr_(cdict)

###############

class OrderedMultisetPartitions_A(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of specified order `d`
    over a fixed alphabet `A`.
    """
    def __init__(self, A, d):
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
        """
        self._alphabet = A
        self._order = d
        OrderedMultisetPartitions.__init__(self, True)

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

    def an_element(self):
        """
        Return a typical ``OrderedMultisetPartition_A``.

        EXAMPLES::

            sage: OrderedMultisetPartitions([2,3,4,5], 3).an_element()
            [{4}, {4}, {4}]
        """
        alpha = Compositions(self._order).an_element()
        co = [Subsets_sk(self._alphabet, a).an_element() for a in alpha]
        return self.element_class(self, map(frozenset, co))

    def random_element(self):
        """
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
        co = [Subsets_sk(self._alphabet, a).random_element() for a in alpha]
        return self.element_class(self, map(frozenset, co))

    def __iter__(self):
        """
        Iterate over ``self``.

        TESTS::

            sage: O = OrderedMultisetPartitions([1, 'a'], 3)
            sage: it = O.__iter__()
            sage: [next(it) for _ in range(O.cardinality())]
            [[{1,'a'}, {'a'}], [{1,'a'}, {1}], [{'a'}, {1,'a'}],
             [{1}, {1,'a'}], [{'a'}, {'a'}, {'a'}], [{'a'}, {'a'}, {1}],
             [{'a'}, {1}, {'a'}], [{'a'}, {1}, {1}], [{1}, {'a'}, {'a'}],
             [{1}, {'a'}, {1}], [{1}, {1}, {'a'}], [{1}, {1}, {1}]]
        """
        return _iterator_order(self._alphabet, self._order)

    def cardinality(self):
        """
        Return the number of ordered partitions of order ``self._order`` on
        alphabet ``self._alphabet``.

        TESTS::

            sage: len(OrderedMultisetPartitions([1, 'a'], 3).list())
            12
            sage: OrderedMultisetPartitions([1, 'a'], 3).cardinality()
            12
        """
        if self._order == 0:
            return ZZ(0)

        # iteration scheme:
        # - start from an integer composition ``alpha`` of ``self._order``.
        # - for each ``a`` in ``alpha``, pick ``a`` letters from ``alphabet``
        min_length = self._order // len(self._alphabet)
        max_length = self._order

        deg = 0
        for k in range(min_length, max_length+1):
            for alpha in Compositions(self._order, length=k, max_part=len(self._alphabet)):
                deg += prod(binomial(len(self._alphabet), a) for a in alpha)
        return ZZ(deg)

class OrderedMultisetPartitions_A_constraints(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of specified order `d`
    over a fixed alphabet `A` satisfying constraints.
    """
    def __init__(self, A, d, **constraints):
        """
        Mimic class ``OrderedMultisetPartitions_A`` to initialize.

        EXAMPLES::

            sage: list(OrderedMultisetPartitions(3, 2, length=3)) # should be empty
            []
            sage: list(OrderedMultisetPartitions([1,2,4], 2, length=1))
            [[{1,2}], [{1,4}], [{2,4}]]

        TESTS::

            sage: C = OrderedMultisetPartitions(3, 2, length=3)
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()
            sage: C = OrderedMultisetPartitions([1,2,4], 4, min_length=3)
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()
        """
        self._alphabet = A
        self._order = d
        OrderedMultisetPartitions.__init__(self, True, alphabet=A, order=d, **constraints)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = OrderedMultisetPartitions([2,3,4,5], 4, length=3)
            sage: O._repr_()
            'Ordered Multiset Partitions of order 4 over alphabet {2, 3, 4, 5}
            with constraint: length=3'
            sage: O = OrderedMultisetPartitions([2,3,4,5], 4, max_length=4, size=10)
            sage: O._repr_()
            'Ordered Multiset Partitions of order 4 over alphabet {2, 3, 4, 5}
            with constraints: max_length=4, size=10'
        """
        cdict = dict(self.constraints)
        cdict.pop("alphabet", None)
        cdict.pop("order", None)
        base_repr = "Ordered Multiset Partitions of order " + str(self._order)
        base_repr += " over alphabet {%s}"%(", ".join(map(str, sorted(self._alphabet))))
        return base_repr + self._constraint_repr_(cdict)

    def an_element(self):
        """
        Return a typical element of ``self``.

        If ``length`` is the only constraint, pick something interesting.
        Else, return the first element satisfying the given constraints.

        EXAMPLES::

            sage: from sage.combinat.multiset_partition_ordered import OrderedMultisetPartitions_A_constraints
            sage: S = frozenset([2,3,4,5])
            sage: O1 = OrderedMultisetPartitions_A_constraints(S, 4, length=2)
            sage: O1.an_element()
            [{2,4,5}, {4}]
            sage: O2 = OrderedMultisetPartitions_A_constraints(S, 4, max_length=4)
            sage: O2.an_element()
            [{2,3,4,5}]
        """
        keys = self.constraints.keys()
        n = len(self._alphabet)
        ell = self._order
        if list(keys) == ["length"]:
            kmin = kmax = k = self.constraints["length"]
            if (ell < kmin) or (n * kmax < ell):
                raise EmptySetError("%s is the empty set"%self)
            alpha = Compositions(ell, length=k, max_part=n).an_element()
            co = [Subsets_sk(self._alphabet, a).an_element() for a in alpha]
            #assume ``co`` satisfies all constraints
            return self.element_class(self, map(frozenset, co))
        else:
            try:
                return next(self.__iter__())
            except StopIteration:
                raise EmptySetError("%s is the empty set"%self)

###############

def _get_multiset(co):
    """
    Construct the multiset (as a sorted tuple) suggested by the lists of lists ``co``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _get_multiset
        sage: L = ((1,), (1, 6), (6, 7), (1,), (1, 3))
        sage: _get_multiset(L)
        (1, 1, 1, 1, 3, 6, 6, 7)
    """
    return tuple(sorted(_concatenate(co)))

def _get_weight(lst):
    """
    Construct the multiset (as a dictionary) suggested by the multiset-as-list ``lst``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _get_weight
        sage: L = (1, 1, 1, 3, 1, 6, 6, 7)
        sage: _get_weight(L)
        {1: 4, 3: 1, 6: 2, 7: 1}
    """
    out = {}
    for k in lst:
        out[k] = out.get(k,0) + 1
    return out

def _has_nonempty_sets(x):
    """
    Blocks should be nonempty sets/lists/tuples of distinct elements.

    TESTS::

        sage: from sage.combinat.multiset_partition_ordered import _has_nonempty_sets
        sage: _has_nonempty_sets([[2,4], {1}, (1,4)])
        True
        sage: _has_nonempty_sets([[2,4], {}, (1,4)])
        False
        sage: _has_nonempty_sets([(2,4), (1,1), (1,4)])
        False
    """
    for block in x:
        if not isinstance(block, (list, tuple, set, frozenset, Set_object)):
            return False
        if not tuple(block) or Set(block).cardinality() != len(tuple(block)):
            return False
    return True

def _union_of_sets(list_of_sets):
    """
    Return the union of a list of iterables as a frozenset.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _union_of_sets
        sage: L = ([1,2,3], Set([1,5,6]), [], range(5,8))
        sage: _union_of_sets(L)
        frozenset({1, 2, 3, 5, 6, 7})
    """
    return reduce(lambda a,b: frozenset(a)|frozenset(b), list_of_sets, frozenset())

def _concatenate(list_of_iters):
    """
    Return the concatenation of a list of iterables as a tuple.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _concatenate
        sage: L = ([1,2,3], Set([4,5,6]), [], range(7,11))
        sage: _concatenate(L)
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    """
    #if not list_of_iters:
    #    return []
    #return reduce(lambda a,b: a+b, list_of_iters)
    return tuple(_ for block in list_of_iters for _ in block)

def _is_finite(constraints):
    """
    Return ``True`` if the dictionary ``constraints`` corresponds to
    a finite collection of ordered multiset partitions.

    If either ``weight`` or ``size`` is among the constraints, then
    the constraints represent a finite collection of ordered multiset
    partitions. If both are absent, one needs ``alphabet`` to be
    present (plus a bound on length or order) in order to have a
    finite collection of ordered multiset partitions.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _is_finite
        sage: W = {"weight": {1:3, 2:3, 4:1}, "length": 5}
        sage: S = {"size": 44, "min_length": 5}
        sage: AO = {"alphabet": range(44), "max_order": 5}
        sage: all(_is_finite(constr) for constr in (W, S, AO))
        True
        sage: AL = {"alphabet": range(44), "min_order": 5}
        sage: _is_finite(AL)
        False
    """
    if "weight" in constraints or "size" in constraints:
        return True
    elif "alphabet" in constraints:
        # Assume the alphabet is finite
        Bounds = Set(["length", "max_length", "order", "max_order"])
        return Bounds.intersection(Set(constraints)) != Set()

def _base_iterator(constraints):
    """
    Return a base iterator for ordered multiset partitions or ``None``.

    If the keys within ``constraints`` dictionary correspond to a finite set
    of ordered multiset partitions, return an iterator. Else, return ``None``.

    EXAMPLES:

    If key ``weight`` is present, ignore all other constraints
    (passes to ``_iterator_weight``)::
    ::

        sage: from sage.combinat.multiset_partition_ordered import _base_iterator
        sage: constraints = {"weight": {1:3, 2:3, 4:1}, "length": 5}
        sage: it = _base_iterator(constraints)
        sage: [next(it) for _ in range(4)] # note the partitions of length 6 and 7
        [[{1}, {1}, {1,2}, {2}, {2}, {4}], [{1}, {1}, {1}, {2}, {2}, {2}, {4}],
         [{1}, {1}, {1,2}, {2}, {2,4}], [{1}, {1}, {1}, {2}, {2}, {2,4}]]

    If key ``size`` is present, pass to ``_iterator_size``, which then
    takes into account whether or not keys ``length`` and ``alphabet``
    are among the constraints::

        sage: constraints = {"size": 5}
        sage: it = _base_iterator(constraints)
        sage: [next(it) for _ in range(8)]
        [[{5}], [{2,3}], [{1,4}], [{4}, {1}], [{1,3}, {1}],
         [{3}, {2}], [{1,2}, {2}], [{3}, {1}, {1}]]

        sage: constraintsL = {"size": 6, "length":2}
        sage: it = _base_iterator(constraintsL)
        sage: [next(it) for _ in range(8)]
        [[{5}, {1}], [{2,3}, {1}], [{1,4}, {1}], [{4}, {2}],
         [{1,3}, {2}], [{3}, {3}], [{3}, {1,2}], [{1,2}, {3}]]

        sage: constraintsA = {"size": 6, "alphabet":frozenset([2, 3])}
        sage: it = _base_iterator(constraintsA)
        sage: list(it)
        [[{3}, {3}], [{2}, {2}, {2}]]

    If key ``alphabet`` is present, the slice may still be infinite, in
    which case ``None`` is returned. Else, use to ``_iterator_order``::

        sage: constraints = {"alphabet": frozenset([3, 4]), "min_length":2}
        sage: _base_iterator(constraints)  is None
        True
        sage: constraints = {"alphabet": frozenset([3, 4]), "max_length":2}
        sage: it = _base_iterator(constraints)
        sage: list(it)
        [[], [{3}], [{4}], [{3,4}], [{3}, {3}], [{3}, {4}],
         [{4}, {3}], [{4}, {4}], [{3,4}, {3}], [{3,4}, {4}],
         [{3}, {3,4}], [{4}, {3,4}], [{3,4}, {3,4}]]
    """
    if "weight" in constraints:
        return _iterator_weight(constraints["weight"])
    elif "size" in constraints:
        return _iterator_size(constraints["size"], \
            constraints.get("length",None), constraints.get("alphabet",None))
    elif "alphabet" in constraints:
        A = constraints["alphabet"]
        # assumes `alphabet` is finite
        min_k = constraints.get("min_length", 0)
        max_k = constraints.get("max_length", infinity)
        min_ord = constraints.get("min_order", 0)
        max_ord = constraints.get("max_order", max_k * len(A))
        max_k = min(max_k, max_ord)
        if "length" in constraints:
            min_k = max_k = constraints["length"]
            min_ord = max(min_ord, min_k)
            max_ord = min(max_ord, len(A) * max_k)
        if "order" in constraints:
            min_ord = max_ord = constraints["order"]
            max_k = min(max_k, max_ord)
            if min_ord:
                min_k =max(1, min_k, min_ord // len(A))
        if infinity not in (max_k, max_ord):
            return chain(*(_iterator_order(A, ord, range(min_k, max_k+1)) \
                        for ord in range(min_ord, max_ord+1)))
    # else
    return None

def _iterator_weight(weight):
    """
    An iterator for the ordered multiset partitions with weight given by
    the dictionary (or weak composition) ``weight``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _iterator_weight
        sage: list(_iterator_weight({'a':2, 'b':1}))
        [[{'a'}, {'a'}, {'b'}], [{'a'}, {'a','b'}], [{'a'}, {'b'}, {'a'}],
         [{'a','b'}, {'a'}], [{'b'}, {'a'}, {'a'}]]
        sage: list(_iterator_weight([3,0,1]))
        [[{1}, {1}, {1}, {3}], [{1}, {1}, {1,3}], [{1}, {1}, {3}, {1}],
         [{1}, {1,3}, {1}], [{1}, {3}, {1}, {1}],
         [{1,3}, {1}, {1}], [{3}, {1}, {1}, {1}]]
    """
    if isinstance(weight, (list, tuple)):
        weight = {k+1: val for k,val in enumerate(weight) if val > 0}
    if isinstance(weight, dict):
        multiset = tuple([k for k in sorted(weight) for _ in range(weight[k])])
    P = OrderedMultisetPartitions_X(tuple(weight.iteritems()))

    if not multiset:
        yield P([])
    else:
        # We build ordered multiset partitions of `X` by permutation + deconcatenation
        for alpha in Permutations_mset(multiset):
            co = _break_at_descents(alpha, weak=True)
            for A in P(co).finer(strong=True):
                yield A

def _iterator_size(size, length=None, alphabet=None):
    r"""
    An iterator for the ordered multiset partitions of integer `n`.

    The degree `n` part of ordered multiset partitions contains all sequences of
    subsets of `\NN_+` whose total sum adds up to `n`.

    If optional argument ``alphabet`` is given, it should be a ``Set`` object.
    Then only yield those `c` with all letters taken from ``alphabet``.

    TESTS::

        sage: from sage.combinat.multiset_partition_ordered import _iterator_size
        sage: list(_iterator_size(3))
        [[{3}], [{1,2}], [{2}, {1}], [{1}, {2}], [{1}, {1}, {1}]]
        sage: list(_iterator_size(5, alphabet={1,3}))
        [[{1,3}, {1}], [{3}, {1}, {1}], [{1}, {1,3}], [{1}, {3}, {1}],
         [{1}, {1}, {3}], [{1}, {1}, {1}, {1}, {1}]]
    """
    # iteration scheme:
    # - start from an integer composition ``alpha``.
    # - for each ``a`` in ``alpha``, pick distinct integers that sum to ``a``
    P = OrderedMultisetPartitions_n(size)
    if alphabet:
        min_p = min(alphabet)
        max_p = max(alphabet)
        for alpha in Compositions(size, length=length):
            for p in cartesian_product([IntegerListsLex(a, min_slope=1, \
                    min_part=min_p, max_part=min(a, max_p)) for a in alpha]):
                if frozenset(_concatenate(p)).issubset(frozenset(alphabet)):
                    yield P([frozenset(list(k)) for k in p])
    else:
        for alpha in Compositions(size, length=length):
            for p in cartesian_product([IntegerListsLex(a, min_slope=1, \
                    min_part=1) for a in alpha]):
                yield P([frozenset(list(k)) for k in p])

def _iterator_order(A, d, lengths=None):
    """
    An iterator for the ordered multiset partitions of order `d` over alphabet `A`.

    If optional argument ``lengths`` is given, it should be a list of integers.
    Then only yield ordered multiset partitions with length in ``lengths``.

    TESTS::

        sage: from sage.combinat.multiset_partition_ordered import _iterator_order
        sage: list(_iterator_order({1,4}, 3))
        [[{1,4}, {1}], [{1,4}, {4}], [{1}, {1,4}], [{4}, {1,4}], [{1}, {1}, {1}],
         [{1}, {1}, {4}], [{1}, {4}, {1}], [{1}, {4}, {4}], [{4}, {1}, {1}],
         [{4}, {1}, {4}], [{4}, {4}, {1}], [{4}, {4}, {4}]]
        sage: list(_iterator_order([1,4], 3, [3]))
        [[{1}, {1}, {1}], [{1}, {1}, {4}], [{1}, {4}, {1}], [{1}, {4}, {4}],
         [{4}, {1}, {1}], [{4}, {1}, {4}], [{4}, {4}, {1}], [{4}, {4}, {4}]]
        sage: list(_iterator_order([1,2,4], 3, [1,2]))[:10]
        [[{1,2,4}],  [{1,2}, {1}], [{1,2}, {2}], [{1,2}, {4}], [{1,4}, {1}],
         [{1,4}, {2}], [{1,4}, {4}], [{2,4}, {1}], [{2,4}, {2}], [{2,4}, {4}]]
        sage: list(_iterator_order([1,4], 3, [4]))
        []
        sage: list(_iterator_order([1,4], 0, [3]))
        []
        sage: list(_iterator_order([1,4], 0, [0,3]))
        [[]]
        sage: list(_iterator_order([1,4], 0))
        [[]]
    """
    A = frozenset(A)
    P = OrderedMultisetPartitions_A(A, d)

    # iteration scheme:
    # start from an integer composition ``alpha`` of ``d``.
    # for each ``a`` in ``alpha``, pick ``a`` letters from ``A``
    n = len(A)
    if not lengths:
        if d:
            lengths = range(max(1, d // n), d+1)
        else:
            lengths = (0,)

    for k in lengths:
        if not k and not d:
            yield P([])
        else:
            for alpha in Compositions(d, length=k, max_part=n):
                for co in cartesian_product([Subsets_sk(A, a) for a in alpha]):
                    yield P(co)

def _partial_sum(lst):
    """
    Return partial sums of elements in ``lst``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _partial_sum
        sage: lst = [1,3,5]
        sage: _partial_sum(lst)
        [0, 1, 4, 9]
    """
    result = [0]
    for i in range(len(lst)):
        result.append(result[-1]+lst[i])
    return result

def _descents(w):
    r"""
    Return descent positions in the word ``w``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _descents
        sage: _descents([1, 2, 3, 2, 2, 1, 4, 3]) == [2, 4, 6]
        True
        sage: _descents([])
        []
    """
    return [j for j in range(len(w)-1) if w[j] > w[j+1]]

def _break_at_descents(alpha, weak=True):
    r"""
    Return the deconcatenation of the composition ``alpha`` at its
    set of descents.

    OUTPUT:

    A list `[a_1, \ldots, a_r]` of nonempty lists whose concatenation
    is ``list(alpha)`` with the property that ``alpha[i] >= alpha[i+1]``
    if and only if positions `i` and `i+1` correspond to different
    lists. (Specifically, ``alpha[i]`` is the last letter of some
    `a_j` and ``alpha[i+1]`` is the first letter of `a_{j+1}`.)

    If the optional argument ``weak`` is ``False``, then only make
    breaks when ``alpha[i] > alpha[i+1]``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _break_at_descents
        sage: _break_at_descents([1, 2, 3, 2, 2, 1, 4, 3])
        [[1, 2, 3], [2], [2], [1, 4], [3]]
        sage: _break_at_descents([1, 2, 3, 2, 2, 1, 4, 3], weak=False)
        [[1, 2, 3], [2, 2], [1, 4], [3]]
        sage: _break_at_descents([])
        []
    """
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
    refinements that are deconcatenations of the list ``sorted(S)`` are returned.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _refine_block
        sage: _refine_block([1, 2, 3])
        {({3}, {1, 2}), ({2}, {3}, {1}), ({1, 2, 3},), ({1, 3}, {2}),
         ({1}, {2, 3}), ({3}, {2}, {1}), ({1}, {2}, {3}), ({3}, {1}, {2}),
         ({2, 3}, {1}), ({2}, {1}, {3}), ({1, 2}, {3}), ({2}, {1, 3}),
         ({1}, {3}, {2})}
        sage: _refine_block([1, 2, 3], strong=True)
        {({1}, {2, 3}), ({1, 2}, {3}), ({1}, {2}, {3}), ({1, 2, 3},)}

    TESTS::

        sage: len(_refine_block([1, 2, 3, 4])) == 1 + binomial(4,1)*2 + binomial(4,2) + binomial(4,2)*factorial(3) + factorial(4)
        True
        sage: len(_refine_block([1, 2, 3, 4], strong=True)) == 2**3
        True
        sage: _refine_block([])
        Traceback (most recent call last):
        ...
        ValueError: S (=[]) must be nonempty
    """
    if not S:
        raise ValueError("S (=%s) must be nonempty"%S)

    if all(s in ZZ for s in S):
        X = sorted(S)
    else:
        X = sorted(S, key=str)
    n = len(X)
    out = []
    if not strong:
        WordSet = IntegerListsLex(min_part=0, max_part=n-1, length=n)
    else:
        WordSet = IntegerListsLex(min_part=0, max_part=n-1, length=n, min_slope=0)

    for w in WordSet:
        if _is_initial_segment(sorted(set(w))):
            a = [set() for _ in range(max(w)+1)]
            for pos in range(n):
                a[w[pos]].add(X[pos])
            out.append(a)
    return Set([tuple(map(Set,x)) for x in out])

def _is_initial_segment(lst):
    r"""
    Return True if ``lst`` is an interval in `\ZZ` of the form `[0, 1, \ldots, n]`.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _is_initial_segment
        sage: _is_initial_segment([1, 2, 3])
        False
        sage: _is_initial_segment([0, 1, 2, 3])
        True
        sage: _is_initial_segment([0])
        True
    """
    return list(range(max(lst)+1)) == lst

def _split_block(S, k=2):
    """
    Return the set of all possible splittings of a set `S` into `k` parts.

    A splitting of `S` is a tuple of (possibly empty) subsets whose union is `S`.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _split_block
        sage: S = [1, 2, 3]
        sage: _split_block(S, 1)
        {({1, 2, 3},)}
        sage: _split_block(S, 2)
        {({3}, {1, 2}), ({}, {1, 2, 3}), ({1, 3}, {2}), ({1}, {2, 3}),
         ({2, 3}, {1}), ({1, 2}, {3}), ({1, 2, 3}, {}), ({2}, {1, 3})}
        sage: _split_block({2, 4}, 3)
        {({}, {2, 4}, {}), ({}, {2}, {4}), ({4}, {2}, {}),
         ({2, 4}, {}, {}), ({4}, {}, {2}), ({}, {4}, {2}),
         ({2}, {4}, {}), ({2}, {}, {4}), ({}, {}, {2, 4})}
    """
    if all(s in ZZ for s in S):
        X = sorted(S)
    else:
        X = sorted(S, key=str)
    n = len(X)
    out = []
    WordSet = IntegerListsLex(min_part=0, max_part=k-1, length=n)
    for w in WordSet:
        a = [set([]) for _ in range(k)]
        for pos in range(n):
            a[w[pos]].add(X[pos])
        out.append(a)
    return Set([tuple(map(Set,x)) for x in out])

def _to_minimaj_blocks(T):
    r"""
    Return a tuple of tuples, representing an ordered multiset partition in
    the minimaj ordering on blocks

    INPUT:

    - ``T`` -- a sequence of row words corresponding to (skew-)tableaux.

    OUTPUT:

    The minimaj bijection `\varphi^{-1}` of [BCHOPSY2017]_ applied to ``T``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_ordered import _to_minimaj_blocks
        sage: co = OrderedMultisetPartitions(14).an_element(); co
        [{2,3}, {2,3}, {4}]
        sage: co.to_tableau()
        [[3, 2], [], [3, 2, 4]]
        sage: _to_minimaj_blocks(co.to_tableau()) == co.minimaj_blocks()
        True
    """
    mu = [(i,) for i in T[-1]]
    breaks = [0] + _descents(T[-1]) + [len(mu)-1]
    T = [T[i][::-1] for i in range(len(T)-1)][::-1]
    for f in range(len(breaks)-1):
        for j in range(breaks[f],breaks[f+1]+1):
            mu[j] += tuple(i for i in T[f] if (mu[j][0] < i or j == breaks[f])
                                               and (j == breaks[f+1] or i <= mu[j+1][0]))
    return tuple(mu)


###############

class MinimajCrystal(UniqueRepresentation, Parent):
    r"""
    Crystal of ordered multiset partitions with `ell` letters from alphabet
    `\{1,2,...,n\}` divided into `k` blocks.

    Elements are represented in the minimaj ordering of blocks as in Benkart, et al.

    .. NOTE::

        Elements are not stored internally as ordered multiset partitions,
        but as certain (pairs of) words stemming from the minimaj bijection
        `\varphi` of Benkart, et al.. See :class:`crystals.Minimaj.Element`
        for further details.

    REFERENCES:

    - [BCHOPSY2017]_

    AUTHORS:

    - Anne Schilling (2018): initial draft
    - Aaron Lauve (2018): changed to use ``Letters`` crystal for elements

    EXAMPLES::

        sage: list(crystals.Minimaj(2,3,2))
        [((2, 1), (1,)), ((2,), (1, 2)), ((1,), (1, 2)), ((1, 2), (2,))]

        sage: b = crystals.Minimaj(3, 5, 2).an_element(); b
        ((2, 3, 1), (1, 3))
        sage: b.e(2)
        ((2, 3, 1), (1, 2))
    """
    def __init__(self, n, ell, k):
        """
        Initialize ``self``.

        TESTS::

            sage: B = crystals.Minimaj(2,3,2)
            sage: B == loads(dumps(B))
            True
            sage: B = crystals.Minimaj(3, 5, 2)
            sage: TestSuite(B).run()

            sage: list(crystals.Minimaj(2,6,3))
            [((1, 2), (2, 1), (1, 2))]
            sage: list(crystals.Minimaj(2,5,2)) # blocks too fat for alphabet
            []
            sage: list(crystals.Minimaj(4,2,3)) # more blocks than letters
            Traceback (most recent call last):
            ...
            ValueError: n (=4), ell (=2), and k (=3) must all be positive integers
        """
        Parent.__init__(self, category=ClassicalCrystals())
        self.n = n
        self.ell = ell
        self.k = k
        if not all([n in ZZ, ell in ZZ, k in ZZ]):
            raise TypeError("n (=%s), ell (=%s), and k (=%s) must all be positive integers"%(n, ell, k))
        if not all([n > 0, ell >= k, k > 0]):
            raise ValueError("n (=%s), ell (=%s), and k (=%s) must all be positive integers"%(n, ell, k))
        self._cartan_type = CartanType(['A',n-1])
        B = Letters(['A', n-1])
        T = tensor([B]*ell)
        self._BT = (B, T)
        self._OMPs = OrderedMultisetPartitions(n, ell, length=k)
        self.module_generators = []
        for co in self._OMPs:
            t = co.to_tableau()
            word = T(*[B(a) for a in _concatenate(t)])
            blocks = [len(h) for h in t]
            breaks = tuple(_partial_sum(blocks))
            mu = self.element_class(self, (word, breaks))
            self.module_generators.append(mu)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: B = crystals.Minimaj(3,4,2); B
            Minimaj Crystal of type A_2 of words of length 4 into 2 blocks
        """
        return "Minimaj Crystal of type A_%s of words of length %s into %s blocks"%(self.n-1, self.ell, self.k)

    def an_element(self):
        """
        Return a typical element of ``self``.

        EXAMPLES::

            sage: B = crystals.Minimaj(4,5,3)
            sage: B.an_element()
            ((4, 1, 3), (3,), (3,))
            sage: B = crystals.Minimaj(2,2,1)
            sage: B.an_element()
            ((1, 2),)
            sage: B = crystals.Minimaj(1,2,1)
            sage: B.an_element()
            Traceback (most recent call last):
            ...
            EmptySetError: Ordered Multiset Partitions of order 2 over alphabet {1}
             with constraint: length=1 is the empty set
        """
        t = self._OMPs.an_element().to_tableau()
        breaks = tuple(_partial_sum([len(h) for h in t]))
        B,T = self._BT
        return self.element_class(self, (T(*[B(a) for a in _concatenate(t)]), breaks))

    def _element_constructor_(self, x):
        """
        Build an element of Minimaj from the ordered multiset partition ``x``.

        EXAMPLES::

            sage: B1 = crystals.Minimaj(4,5,3); b = B1.an_element(); b
            ((4, 1, 3), (3,), (3,))
            sage: B1._element_constructor_(list(b))
            ((4, 1, 3), (3,), (3,))
            sage: B1._element_constructor_([[1,2,3], [2], [2]])
            ((3, 1, 2), (2,), (2,))
            sage: B2 = crystals.Minimaj(5,5,3)
            sage: B2._element_constructor_(b)
            ((4, 1, 3), (3,), (3,))
        """
        # Allow ``x`` to be either of:
        # - an ordered multiset partition in ``self._OMPs``;
        # - an element of another Minimaj crystal with
        #   + same `ell` and `k`, and
        #   + all letters smaller or equal to ``self._n``.
        x = list(x)
        if x in self:
            t = self._OMPs(x).to_tableau()
            breaks = tuple(_partial_sum([len(h) for h in t]))
            B,T = self._BT
            return self.element_class(self, (T(*[B(a) for a in _concatenate(t)]), breaks))
        else:
            raise ValueError("cannot convert %s into an element of %s"%(x, self))

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is an element of ``self`` or an ordered multiset partition.

        EXAMPLES::

            sage: B1 = crystals.Minimaj(4,5,3); b1 = B1.an_element(); b1
            ((4, 1, 3), (3,), (3,))
            sage: B2 = crystals.Minimaj(5,5,3); b2 = B2.an_element(); b2
            ((4, 5, 1), (3,), (3,))
            sage: b2a = B2(((1,2,3), (2,), (3,))); b2a
            ((3, 1, 2), (2,), (3,))
            sage: b1 in B2
            True
            sage: b2 in B1
            False
            sage: b2a in B1
            True
        """
        if isinstance(x, MinimajCrystal.Element):
            if x.parent() == self:
                return True
            else:
                return list(x) in self._OMPs
        else:
            return x in self._OMPs

    def from_tableau(self, t):
        r"""
        Return the bijection `\varphi^{-1}` of [BCHOPSY2017]_ applied to ``t``.

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
            sage: B1 = crystals.Minimaj(3,6,2)
            sage: B1.from_tableau(t)
            Traceback (most recent call last):
            ...
            ValueError: ((1, 2, 3), (3, 1), (2,)) is not an element of
             Minimaj Crystal of type A_2 of words of length 6 into 2 blocks
        """
        mu = _to_minimaj_blocks(t)
        if mu in self:
            return self(mu)
        else:
            raise ValueError("%s is not an element of %s"%(mu, self))

    def val(self, q='q'):
        """
        Return `Val` polynomial corresponding to ``self``.

        EXAMPLES:

        Verifying Example 4.5 from [BCHOPSY2017]_::

            sage: B = crystals.Minimaj(3, 4, 2) # for `Val_{4,1}^{(3)}`
            sage: B.val()
            (q^2+q+1)*s[2, 1, 1] + q*s[2, 2]
        """
        H = [self._OMPs(list(b)) for b in self.highest_weight_vectors()]
        Sym = SymmetricFunctions(ZZ[q])
        q = Sym.base_ring().gens()[0]
        s = Sym.schur()
        return sum((q**(t.minimaj()) * s[sorted(t.weight().values(), reverse=True)]
                   for t in H), Sym.zero())

    class Element(ElementWrapper):
        r"""
        An element of a Minimaj crystal.

        .. NOTE::

            Minimaj elements `b` are stored internally as pairs `(w, breaks)`, where:

            - `w` is a word of length ``self.parent().ell`` over the letters
                `1` up to ``self.parent().n``;
            - `breaks` is a list of de-concatenation points to turn `w` into a list
                of row words of (skew-)tableaux that represent `b` under the minimaj
                bijection `\varphi` of [BCHOPSY2017]_.

            The pair `(w, breaks)` may be recovered via ``b.value``.
        """
        def _repr_(self):
            """
            Return the string representation of ``self``.

            EXAMPLES::

                sage: b = crystals.Minimaj(4,5,3).an_element(); b._repr_()
                '((4, 1, 3), (3,), (3,))'
            """
            return repr(self._minimaj_blocks_from_word_pair())

        def __iter__(self):
            """
            Iterate over ``self._minimaj_blocks_from_word_pair()``.

            EXAMPLES::

                sage: b = crystals.Minimaj(4,5,3).an_element(); b
                ((4, 1, 3), (3,), (3,))
                sage: b.value
                ([3, 1, 4, 3, 3], (0, 2, 2, 5))
                sage: list(b)
                [(4, 1, 3), (3,), (3,)]
            """
            return self._minimaj_blocks_from_word_pair().__iter__()

        def _minimaj_blocks_from_word_pair(self):
            """
            Return the tuple of tuples (in the minimaj ordering on blocks
            of ordered multiset partitions) corresponding to ``self``.

            EXAMPLES::

                sage: B = crystals.Minimaj(4,5,3)
                sage: b = B.an_element(); b.value
                ([3, 1, 4, 3, 3], (0, 2, 2, 5))
                sage: b._minimaj_blocks_from_word_pair()
                ((4, 1, 3), (3,), (3,))
            """
            return _to_minimaj_blocks(self.to_tableau())

        def to_tableau(self):
            r"""
            Return the image of the ordered multiset partition ``self`` under
            the minimaj bijection `\varphi` of [BCHOPSY2017]_.

            EXAMPLES::

                sage: B = crystals.Minimaj(4,5,3)
                sage: b = B.an_element(); b
                ((4, 1, 3), (3,), (3,))
                sage: b.to_tableau()
                [[3, 1], [], [4, 3, 3]]
            """
            w, breaks = self.value
            return [[ZZ(w[a].value) for a in range(breaks[j], breaks[j+1])]
                        for j in range(len(breaks)-1)]

        def e(self, i):
            r"""
            Return `e_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.Minimaj(4,3,2)
                sage: b = B.an_element(); b
                ((2, 3), (3,))
                sage: [b.e(i) for i in range(1,4)]
                [((1, 3), (3,)), ((2,), (2, 3)), None]
            """
            P = self.parent()
            w, breaks = self.value
            if w.e(i) is None:
                return None
            w = w.e(i)
            return P.element_class(P, (w, breaks))

        def f(self,i):
            r"""
            Return `f_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.Minimaj(4,3,2)
                sage: b = B.an_element(); b
                ((2, 3), (3,))
                sage: [b.f(i) for i in range(1,4)]
                [None, None, ((2, 3), (4,))]
            """
            P = self.parent()
            w, breaks = self.value
            if w.f(i) is None:
                return None
            w = w.f(i)
            return P.element_class(P, (w, breaks))
