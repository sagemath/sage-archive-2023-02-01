r"""
Ordered multiset partitions

An ordered multiset partition `c` of a multiset `X` is a list of subsets of `X` (not multisets),
called the blocks of `c`, whose multi-union is `X`.

This module provides tools for manipulating ordered multiset partitions.

REMARK:

    An (ordered) multiset partition connotates an (ordered) collection of *multisets* (not sets),
    whose multi-union is a given multiset.
    Perhaps a better name would be "strict" ordered multiset partition.
    Terminology is taken from:
      - [HRWxx] - arXiv:1509.07058
      - [HRSxx] - arXiv:1609.07575
      - [BCHOPSYxx] - arXiv:1707.08709v2

EXAMPLES::

    sage: OrderedMultisetPartition([[5, 3], [1, 3]])
    [{3,5}, {1,3}]
    sage: list(OrderedMultisetPartitions(4))
    [[{1}, {1}, {1}, {1}], [{1}, {1}, {2}], [{1}, {2}, {1}], [{1}, {1,2}], [{1}, {3}],
    [{2}, {1}, {1}], [{2}, {2}], [{1,2}, {1}], [{3}, {1}], [{1,4}], [{4}]]

AUTHORS:

- Aaron Lauve, based on ``composition.py``
"""
#*****************************************************************************
#       Copyright (C) 2018 Aaron Lauve       <lauve at math.luc.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#              http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.cartesian_product import cartesian_product
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.sets.set import Set
from sage.sets.family import Family
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.rings.integer_ring import ZZ

from sage.combinat.combinat import CombinatorialElement
from sage.combinat.composition import Compositions
from sage.combinat.partition import RegularPartitions_n
from sage.combinat.words.words import Words
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.shuffle import ShuffleProduct, ShuffleProduct_overlapping
from functools import reduce
from sage.calculus.var import var

class OrderedMultisetPartition(CombinatorialElement):
    r"""
    Ordered multiset partitions

    An ordered multiset partition `c` of a multiset `X` is a list of subsets (not multisets),
    called the blocks of `c`, whose multi-union is `X`.

    EXAMPLES:

    The simplest way to create a ordered multiset partition is by specifying its
    entries as a list, tuple (or other iterable)::

        sage: OrderedMultisetPartition([[3],[2,1]])
        [{3}, {1,2}]
        sage: OrderedMultisetPartition(((3,), (1,2)))
        [{3}, {1,2}]
        sage: OrderedMultisetPartition(set([i]) for i in range(2,5))
        [{2}, {3}, {4}]

    You can also create a ordered multiset partition `c` from a list of positive integers
    and from a list of nonnegative integers. In the former case, each integer is given
    its own block of `c`. In the latter case, zeros separate the blocks of `c`::

        sage: OrderedMultisetPartition(from_list=[i for i in range(2,5)])
        [{2}, {3}, {4}]
        sage: OrderedMultisetPartition(from_zero_list=[1, 1, 0, 3, 5, 0, 2, 1])
        Error: {1, 1} is not a valid block
        sage: OrderedMultisetPartition(from_zero_list=[1, 0, 1, 3, 5, 0, 2, 1])
        [{1}, {1,3,5}, {1,2}]

    EXAMPLES::

        sage: C = OrderedMultisetPartition([[3,1],[2]])
        sage: TestSuite(C).run()
    """
    @staticmethod
    def __classcall_private__(cls, co=None, from_list=None, from_zero_list=None):
        """
        This constructs a list of lists from optional arguments, delegating the
        construction of a :class:`OrderedMultisetPartition` to the ``element_class()``
        call of the appropriate parent.

        EXAMPLES::

            sage: OrderedMultisetPartition([[3], [2,1]])
            [{3}, {1,2}]
            sage: OrderedMultisetPartition(from_list=[2, 3, 4, 5])
            [{2}, {3}, {4}, {5}]
            sage: OrderedMultisetPartition(from_zero_list=[1, 1, 0, 3, 5, 0, 2, 1])
            Error: {1, 1} is not a valid block
            sage: OrderedMultisetPartition(from_zero_list=[1, 0, 1, 3, 5, 0, 2, 1])
            [{1}, {1,3,5}, {1,2}]
            sage: OrderedMultisetPartition([2, 3, 4, 5])
            [{2}, {3}, {4}, {5}]
            sage: OrderedMultisetPartition([1, 0, 1, 3, 5, 0, 2, 1])
            [{1}, {1,3,5}, {1,2}]
        """
        if from_list is not None:
            return OrderedMultisetPartitions().from_list(from_list)
        elif from_zero_list is not None:
            return OrderedMultisetPartitions().from_zero_list(from_zero_list)
        else:
            return OrderedMultisetPartitions()(co)

    def _repr_(self):
        # each block of `self` should be stored as a sorted list
        string_parts = map(lambda block: str(block).replace(" ","").replace("[","{").replace("]","}"), self)
        return "[" + ", ".join(string_parts) + "]"

    def __getitem__(self, i):
        """
        Return the concatenation of two ordered multiset partitions.

        EXAMPLES::
        """
        return list(self)[i]

    def __add__(self, other):
        """
        Return the concatenation of two ordered multiset partitions.

        EXAMPLES::

            sage: OrderedMultisetPartition(from_list=[1, 1, 3]) + OrderedMultisetPartition(from_zero_list=[4, 1, 0, 2])
            [{1}, {1}, {3}, {1,4}, {2}]

        TESTS::

            sage: OrderedMultisetPartition([]) + OrderedMultisetPartition([]) == OrderedMultisetPartition([])
            True
            sage: OrderedMultisetPartition([1,0,1,0,3]) + OrderedMultisetPartition([1,4,0,3]) == OrderedMultisetPartition([{1}, {1}, {3}, {1,4}, {2}])
            True
        """
        return OrderedMultisetPartitions()(list(self)+list(other))

    @combinatorial_map(order=2, name='reversal')
    def reversed(self):
        r"""
        Return the reverse ordered multiset partition of ``self``.

        The reverse of a ordered multiset partition `(B_1, B_2, \ldots, B_k)`
        is defined as the ordered multiset partition `(B_k, B_{k-1}, \ldots, B_1)`.

        EXAMPLES::

            sage: C = OrderedMultisetPartition(from_zero_list=[1, 0, 1, 3, 0, 2, 3, 4]); C
            [{1}, {1,3}, {2,3,4}]
            sage: C.reversed()
            [{2,3,4}, {1,3}, {1}]
        """
        return self.parent()(reversed(self))

    def shape_from_cardinality(self):
        r"""
        Return a composition that records the cardinality of each block of ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartition(from_zero_list=[3, 4, 1, 0, 2, 0, 1, 2, 3, 7]); C
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
            sage: OrderedMultisetPartition([]).shape_from_size() == Composition([])
            True
        """
        return Composition([sum(k) for k in self])

    def letters(self):
        r"""
        Return a sorted list of distinct integers occurring within the blocks of ``self``.

        .. SEEALSO::

            :meth:`Words.letters`
        """
        s = union_of_sets(list(self))
        return sorted(s)

    def multiset(self, as_dict=False):
        r"""
        Return the multiset corresponding to ``self`` as a list or as a dictionary.

        Counts frequency of each letter in ``self.letters()``

        EXAMPLES::

            sage: C = OrderedMultisetPartition([3, 4, 1, 0, 2, 0, 1, 2, 3, 7]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.multiset()
            [1, 1, 2, 2, 3, 3, 4, 7]
            sage: C.multiset(as_dict=True)
            {1:2, 2:2, 3:3, 4:1, 7:1}
            sage: OrderedMultisetPartition([]).multiset() == []
            True
        """
        w = self.weight()
        if as_dict:
            return {i+1:w[i] for i in range(len(w))}
        else:
            return [(i+1) for _ in range(w[i]) for i in range(len(w))]

    def max_letter(self):
        r"""
        Return the maximum integer appearing in ``self.letters()``.
        """
        if len(self) == 0:
            return 0
        else:
            return max(self.letters())

    def size(self):
        """
        Return the size of ``self``, that is the sum of all integers in all blocks.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([3, 4, 1, 0, 2, 0, 1, 2, 3, 7]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.size()
            23
            sage: C.size() == sum(C.shape_from_size())
            True
            sage: OrderedMultisetPartition([[7,1],[3]]).size()
            11
            sage: OrderedMultisetPartition([]).size() == 0
            True
        """
        return sum([sum(block) for block in self])
        ## alternatively...
        #w = self.weight()
        #return sum([w[i]*(i+1) for i in range(len(w))])

    def order(self):
        """
        Return the total number of integers in all blocks.

        EXAMPLES::

            sage: C = OrderedMultisetPartition([3, 4, 1, 0, 2, 0, 1, 2, 3, 7]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.order()
            8
            sage: C.order() == sum(C.weight())
            True
            sage: C.order() == sum(C.shape_from_cardinality())
            True
            sage: OrderedMultisetPartition([[7,1],[3]]).order()
            3
        """
        return sum([len(block) for block in self])

    def length(self):
        """
        Return the number of blocks.

        EXAMPLES::

            sage: OrderedMultisetPartition([[7,1],[3]]).length()
            2
        """
        return len(self)

    def weight(self, as_dict=False):
        """
        Count the number of instances `n_i` for each distinct positive integer `i`
        across all blocks of ``self``. Return as a list `[n_1, n_2, n_3, ..., n_k]`,
        where `k` is the max letter in ``self.letters()``.

        Alternatively, return as a dictionary, with keys being only the letters
        in ``self.letters()`` and values being their (positive) with positive frequency.

        EXAMPLES::

            sage: OrderedMultisetPartition([[6,1],[1,3],[1,3,6]]).weight()
            [3, 0, 2, 0, 0, 2]
            sage: OrderedMultisetPartition([[6,1],[1,3],[1,3,6]]).weight(as_dict=True)
            {1:3, 3:2, 6:2}
            sage: OrderedMultisetPartition([]).weight() == []
            True
        """
        a = [0 for _ in range(self.max_letter())]
        for block in self:
            for i in block:
                a[i-1] += 1
        if as_dict:
            return {i+1:a[i] for i in range(len(a)) if a[i] !=0}
        else:
            return list(a)

    def deconcatenate(self, k=2):
        """
        Choose `k-1` cut points to deconcatenate the ordered list ``self`` into
        k (possibly) empty lists. (A "balls and boxes" operation.)

        This is not to be confused with ``self.split()``,
        which splits each block of ``self`` before making k-tuples of ordered multiset partitions.

        Returns a set of k-tuples of ordered multiset partitions.

        EXAMPLES::

            sage: OrderedMultisetPartition([[7,1],[3]]).deconcatenate()
            {([], [{1,7}, {3}]), ([{1,7}], [{3}]), ([{1,7}, {3}], [])}
            sage: OrderedMultisetPartition([[1,3,7]]).deconcatenate(3)
            {([], [], [{1,3,7}]), ([], [{1,3,7}], []), ([{1,3,7}], [], [])}

        TESTS::

            sage: C = OrderedMultisetPartition([1,2,3,4,5]); C
            [{1}, {2}, {3}, {4}, {5}]
            sage: all( C.deconcatenate(k).cardinality()
            ....:      == binomial(C.length() + k, k)
            ....:      for k in range(1, 5) )
            True
        """
        P = OrderedMultisetPartitions()
        out = Set([])
        l = self.length()
        for c in Compositions(l):
            ps = [0]+c.partial_sums()
            out.append(tuple([P(self[ps[i]:ps[i+1]]) for i in range(len(ps)-1)]))
        return out

    def split(self, k=2):
        """
        Split each block of ``self`` into a list of k (possibly empty) sets, then arrange
        the corresponding sets, block-by-block, to produce k ordered multiset partitions.

        This is not to be confused with ``self.deconcatenate()``,
        which may be viewed as the preimage of the `+` operation on ordered multiset partitions.

        Returns a set of k-tuples of ordered multiset partitions.

        EXAMPLES::

            sage: OrderedMultisetPartition([[1,2],[4]]).split()
            {([{1,2}, {4}], []), ([], [{1,2}, {4}]), ([{4}], [{1,2}]), ([{1}, {4}], [{2}]),
            ([{1}], [{2}, {4}]), ([{2}], [{1}, {4}]), ([{2}, {4}], [{1}]), ([{1,2}], [{4}])}
            sage: OrderedMultisetPartition([[1,2]]).split(3)
            {([{1}], [], [{2}]), ([], [{1,2}], []), ([{2}], [], [{1}]), ([], [{2}], [{1}]),
            ([{2}], [{1}], []), ([], [{1}], [{2}]), ([{1}], [{2}], []), ([{1,2}], [], []), ([], [], [{1,2}])}

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
        def P(cartprod):
            return OrderedMultisetPartitions()([k for k in cartprod if len(k)>0])
        # corner case
        if not self:
            return Set([tuple([self]*k)])
        else:
            out = set()
            tmp = cartesian_product([split_block(block, k) for block in self])
            for t in tmp:
                out.add(tuple([P(c) for c in zip(*t)]))
            return Set(out)

    def finer(self):
        """
        Return the set of ordered multiset partitions that are finer than ``self``.

        An ordered multiset partition `A` is finer than another `B` if, reading left-to-right,
        every block of `B` is the union of some consecutive blocks of `A`.

        EXAMPLES::

            sage: C =  OrderedMultisetPartition([[3,2]]).finer()
            sage: C.cardinality()
            3
            sage: C.list()
            [[{2,3}], [{3}, {2}], [{2}, {3}]]
            sage: OrderedMultisetPartition([]).finer()
            {[]}
        """
        if not self:
            return Set([self])
        else:
            tmp = cartesian_product([refine_block(block) for block in self])
            print list(tmp)
            return Set([OrderedMultisetPartitions()(concatenate(map(list,c))) for c in tmp])

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
        if self.weight() != OrderedMultisetPartitions()(co).weight():
            return False

        # trim common prefix and suffix
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

        - ``grouping`` -- a composition (or list or tuple) whose sum is the length of ``self``

        EXAMPLES:

        Let us start with the composition::

            sage: C = OrderedMultisetPartition(from_zero_list=[4,1,5,0,2,0,7,1]); C
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
            ValueError: [{1,2,4,5,1,7}] is not a valid ordered multiset partition

        """
        if sum(grouping) != self.length():
            raise ValueError("%s is a composition of `self.length()` (=%s)"%(grouping, self.length()))

        valid = True
        result = [None] * len(grouping)
        j = 0
        for i in range(len(grouping)):
            result[i] = self[j:j+grouping[i]]
            # check that grouping[i] is allowed, i.e., `|A\cup B| = |A| + |B|`
            strict_size = sum(map(len,result[i]))
            size = len(union_of_sets(result[i]))
            if size < strict_size:
                valid = False
            j += grouping[i]

        if not valid:
            str_rep = '['
            for i in range(len(grouping)):
                st = ",".join(map(lambda k: str(Set(k)), result[i])).replace('},{',',').replace(' ', '')
                str_rep += st
            str_rep = str_rep.replace("}{", "}, {") + "]"
            raise ValueError("%s is not a valid ordered multiset partition"%(str_rep))
        else:
            out = []
            for i in range(len(grouping)):
                out.append(union_of_sets(result[i]))
            return OrderedMultisetPartition(out)

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
            [[{1,4,5}, {1,2,7}], [{1,2,4,5}, {1,7}], [{1,4,5}, {2}, {1,7}]]

        Some extreme cases::

            sage: list(OrderedMultisetPartition([[5]]).fatter())
            [[{5}]]
            sage: list(OrderedMultisetPartition([]).fatter())
            [[]]
            sage: OrderedMultisetPartition(from_list=[1,2,3,4]).fatter().issubset(OrderedMultisetPartition([[1,2,3,4]]).finer())
            True
        """
        out = set([])
        for c in Compositions(self.length()):
            try:
                out.add(self.fatten(c))
            except ValueError:
                pass
        return Set(out)


    def minimaj(self):
        r"""
        A statistic on ordered multiset partitions, which we define via an example.

        1. Sort each block in the list ``self`` as in ``self.minimaj_order()``
           - in:   [{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}]
           - out:  [[5,7,1], [2,4], [5,6], [4,6,8], [3,1], [1,2,3]]

        2. Record the indices where descents in this word occur.
           - word:      [5, 7, 1 / 2, 4 / 5, 6 / 4, 6, 8 / 3, 1 / 1, 2, 3]
           - indices:    1  2  3   4  5   6  7   8  9 10  11 12  13 14 15
           - descents:  {   2,               7,       10, 11}

        3. Compute the sum of the descents
           - minimaj = 2 + 7 + 10 + 11 = 30

        See:
        [HRWxx] - arXiv:1509.07058

        EXAMPLES::

            sage: C = OrderedMultisetPartition([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}])
            sage: C, C.to_minimaj_composition()
            ([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}],
             [5, 7, 1, 2, 4, 5, 6, 4, 6, 8, 3, 1, 1, 2, 3])
            sage: C.minimaj()
            ([5, 7, 1, 2, 4, 5, 6, 4, 6, 8, 3, 1, 1, 2, 3], 30)
            sage: C = OrderedMultisetPartition([{2,4}, {1,2,3}, {1,6,8}, {2,3}])
            sage: C, C.to_minimaj_composition()
            ([{2,4}, {1,2,3}, {1,6,8}, {2,3}], [2, 4, 1, 2, 3, 6, 8, 1, 2, 3])
            sage: C.minimaj()
            9
            sage: OrderedMultisetPartition([]).minimaj()
            0
        """
        C = self.to_minimaj_composition()
        D = [j+1 for j in range(len(C)-1) if C[j]>C[j+1]]
        return sum(D)

    def to_minimaj_composition(self):
        r"""
        Returns a composition derived by concatenating the blocks of ``self`` after first ordering them
        in a prescribed order, which we define via an example.

        Sort the blocks `[B_1,...,B_k]` of ``self`` from right to left via:
            
        1. Sort the last block `B_k` in increasing order, call it the word `W_k`
        
        2. If blocks `B_{i+1}, \ldots, B_k` have been converted to words `W_{i+1}, \ldots, W_k`,
           use the letters in `B_i` to make the unique word `W_i` that has a
           factorization `W_i=(u,v)` satisfying:
            - letters of `u` and `v` appear in increasing order, with `v` possibly empty
            - letters in `vu` appear in increasing order
            - ``v[-1]`` is the largest letter `a \in B_i` satisfying ``a <= W_{i+1}[0]``

        EXAMPLES::

            sage: C = OrderedMultisetPartition(from_zero_list=[2,1,0,1,2,3,0,1,2,0,3,0,1]); C
            [{1,2}, {1,2,3}, {1,2}, {3}, {1}]
            sage: C.to_minimaj_composition()
            [1, 2, 2, 3, 1, 1, 2, 3, 1]
            sage: C = OrderedMultisetPartition([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}]); C
            [{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}]
            sage: C.to_minimaj_composition()
            [5, 7, 1, 2, 4, 5, 6, 4, 6, 8, 3, 1, 1, 2, 3]
            sage: OrderedMultisetPartition([]).to_minimaj_composition()
            []
        """
        B = [sorted(block) for block in self]
        if len(B) == 0:
            return Composition([])
        C = B[-1]
        for i in range(1,len(B)):
            Bi = B[-1-i]
            j0 = -1
            for j in range(len(Bi)):
                if Bi[j]<=C[0]:
                    j0 = j
                else:
                    break
            Bi = Bi[j0+1:]+Bi[:j0+1]
            C = Bi + C
        return Composition(C)

    def major_index(self):
        r"""
        A statistic on ordered multiset partitions, which we define via an example.

        1. Sort each block in the list ``self`` in descending order to create a word `w`,
           keeping track of the original separation into blocks:
           - in:  [{3,4,5}, {2,3,4}, {1}, {4,5}]
           - out: [ 5,4,3 /  4,3,2 /  1 /  5,4] 

        2. Create a sequence `v = (v_0, v_1, v_2, \ldots)`  of length ``self.order()+1``, 
           built recursively by:
           - `v_0 = 0`
           - `v_j = v_{j-1} + \delta(j)`, where `\delta(j) = 1` if `j` is the index of 
             an end of a block, and zero otherwise.
              + in:    [ 5,4,3 /  4,3,2 /  1 /  5,4] 
              + out: (0, 0,0,1,   1,1,2,   3,   3,4)

        3. Compute the `\sum_j v_j` restricted to descent positions: those `j` with `w_j > w_{j+1}`:
           - in:    [5, 4, 3, 4, 3, 2, 1, 5, 4]
                  (0 0, 0, 1, 1, 1, 2, 3, 3, 4)
           - maj :=  0 +0    +1 +1 +2    +3     = 7

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
        for block in ew:
            for (i,wi) in sorted(block, reverse=True):
                vj = v[-1]
                if i == 0:
                    vj += 1
                v.append(vj)
                w.append(wi)
        maj = [v[j+1] for j in range(len(w)-1) if w[j]>w[j+1]]
        return sum(maj)

    def _multiset_partition(self):
        """
        Return the (strict) (unordered) multiset partition obtained by sorting ``self`` by least element in each block,
        realized as a list of sets.

        EXAMPLES::

            sage: C = OrderedMultisetPartition(from_zero_list=[2,1,0,1,2,3,0,1,2,0,3,0,1]); C
            [{1,2}, {1,2,3}, {1,2}, {3}, {1}]
            sage: C._multiset_partition()
            [{1}, {1,2}, {1,2}, {1,2,3}, {3}]
        """
        return sorted(map(Set, self))


    def shuffle_product(self, other, overlap=False):
        r"""
        Return the shuffles (with multiplicity) of parts of ``self`` with parts of ``other``.
        In case overlap is True, it also takes the union of adjacent disjoint parts (one from ``self`` and the other from ``other``)

        .. SEEALSO::

            :meth:`Composition.shuffle`

        EXAMPLES::

            sage: A = OrderedMultisetPartition(from_zero_list=[2,1,3,0,1,2]); A
            [{1,2,3}, {1,2}]
            sage: B = OrderedMultisetPartition([[3,4]]); B
            [{3,4}]
            sage: C = OrderedMultisetPartition([[4,5]]); C
            [{4,5}]
            sage: A.shuffle_product(B)
            [[{1,2,3}, {1,2}, {3,4}], [{1,2,3}, {3,4}, {1,2}], [{3,4}, {1,2,3}, {1,2}]]
            sage: A.shuffle_product(B, overlap=True)
            [[{1,2,3}, {1,2}, {3,4}], [{1,2,3}, {3,4}, {1,2}], [{1,2,3}, {1,2,3,4}], [{3,4}, {1,2,3}, {1,2}]]
            sage: A.shuffle_product(C, overlap=True)
            [[{1,2,3}, {1,2}, {4,5}], [{1,2,3}, {4,5}, {1,2}], [{1,2,3}, {1,2,4,5}], [{4,5}, {1,2,3}, {1,2}], [{1,2,3,4,5}, {1,2}]]
        """
        other = OrderedMultisetPartition(other)
        if overlap is True:
            for term in ShuffleProduct_overlapping(self, other):
                if term.size() == self.size() + other.size():
                    yield term
        else:
            for term in ShuffleProduct(self, other):
                yield term

    # provided as an alternative sorting scheme to
    # that served up by `OrderedMultisetPartitions().__iter__()`
    def sll_gt(self, co2):
        """
        Return ``True`` if the composition ``self`` is greater than the
        composition ``co2`` with respect to the sll-ordering; otherwise,
        return ``False``.

        The sll-ordering is a total order on the set of all ordered multiset partitions
        defined as follows: A ordered multiset partition `I` is greater than a
        ordered multiset partition `J` if and only if one of the following conditions
        holds:

        - The size of `I` is greater than the size of `J`.

        - The size of `I` equals the size of `J`, but the length of `I`
          is greater than the length of `J`. (Here, `length` is measured by total
          number of integers appearing across all blocks of the ordered multiset partition.)

        - The size of `I` equals the size of `J`, and the length of `I`
          equals the length of `J`, but `I` is lexicographically
          greater than `J`. (Here, letters are subsets, so for the lexicographic ordering,
          we say letters `a` > `b` whenever `sum(a)` > `sum(b)` or `sorted(a)` > `sorted(b)`.)

        ("sll-ordering" is short for "size, length, lexicographic
        ordering".)

        EXAMPLES::

            sage: C = OrderedMultisetPartition(from_zero_list=[2,1,0,1,2,3,0,1,2,0,3,0,1])
            sage: D1 = OrderedMultisetPartition(from_zero_list=[2,4,0,2,3])
            sage: D2 = OrderedMultisetPartition(from_zero_list=[2,4,0,2,4,0,2,3])
            sage: D3 = OrderedMultisetPartition(from_zero_list=[2,4,0,2,4,0,1,3])
            sage: D4 = OrderedMultisetPartition(from_zero_list=[1,2,0,2,4,0,2,0,1,0,1,2,0,1])
            sage: print(repr(C) + ' >? ' + repr(D1)); C.sll_gt(D1)  #size
            [{1,2}, {1,2,3}, {1,2}, {3}, {1}] >? [{2,4}, {2,3}]
            True
            sage: print(repr(C) + ' >? ' + repr(D2)); C.sll_gt(D2)  #size
            [{1,2}, {1,2,3}, {1,2}, {3}, {1}] >? [{2,4}, {2,4}, {2,3}]
            False
            sage: print(repr(C) + ' >? ' + repr(D3)); C.sll_gt(D3)  #length
            [{1,2}, {1,2,3}, {1,2}, {3}, {1}] >? [{2,4}, {2,4}, {1,3}]
            True
            sage: print(repr(C) + ' >? ' + repr(D4)); C.sll_gt(D4)  #lex
            [{1,2}, {1,2,3}, {1,2}, {3}, {1}] >? [{1,2}, {2,4}, {2}, {1}, {1,2}, {1}]
            False
        """
        co1 = self
        # check size
        if (co1).size() > (co2).size():
            return True
        elif (co1).size() < (co2).size():
            return False
        # check length
        if (co1).length() > (co2).length():
            return True
        if (co1).length() < (co2).length():
            return False
        # check lex.
        for i in range(min(len(co1), len(co2))):
            if sum(co1[i]) > sum(co2[i]):
                return True
            elif sum(co1[i]) < sum(co2[i]):
                return False
            if sorted(co1[i]) > sorted(co2[i]):
                return True
            elif sorted(co1[i]) < sorted(co2[i]):
                return False
        return False

##############################################################

class OrderedMultisetPartitions(UniqueRepresentation, Parent):
    r"""
    Set of ordered multiset partitions.

    An ordered multiset partition `c` of a multiset `X` is a list of subsets (not multisets),
    called the blocks of `c`, whose multi-union is `X`.

    Alternatively:
    An ordered multiset partition `c` of a nonnegative integer `n` is a list of
    subsets of positive integers (called blocks) with total sum of all integers
    appearing in all blocks equal to `n`.

    Valid keywords are:
    (See corresponding methods in see :class`OrderedMultisetPartition` for more details.)
    - ``length=k``     (integer `k`) specifies the number of blocks in the partition
    - ``min_length=k`` (integer `k`) specifies minimum number of blocks in the partition
    - ``max_length=k`` (integer `k`) specifies maximum number of blocks in the partition
    - ``alphabet=S``  (iterable `S`) specifies allowable positive integers for use in the partition
    - ``weight=C``     (weak composition `C`) specifies frequency of occurrences of each positive integer
    - ``order=n``      (integer `n`) specifies the number of integers in the partition (=sum of cardinalities of its blocks)
    - ``min_order=n``  (integer `n`) specifies minimum number of integers in the partition (=sum of cardinalities of its blocks)
    - ``max_order=n``  (integer `n`) specifies maximum number of integers in the partition (=sum of cardinalities of its blocks)

    EXAMPLES:

    There are 5 ordered multiset partitions of 3::

        sage: OrderedMultisetPartitions(3).cardinality()
        5

    Here is the list of them::

        sage: OrderedMultisetPartitions(3).list()
        [[{1}, {1}, {1}], [{1,2}], [{2}, {1}], [{1}, {2}], [{3}]]

    You can use the ``.first()`` method to get the 'first' ordered multiset partition of
    a number::

        sage: OrderedMultisetPartitions(3).first()
        [{1}, {1}, {1}]

    You can also calculate the 'next' ordered multiset partition given the current
    one::

        sage: OrderedMultisetPartitions(3).next([[2],[1]])
        [{3}]

    If `n` is not specified, this returns the combinatorial class of
    all (non-negative) ordered multiset partitions::

        sage: OrderedMultisetPartitions()
        Ordered Multiset Partitions of non-negative integers
        sage: [] in OrderedMultisetPartitions()
        True
        sage: [[2,3],[1]] in OrderedMultisetPartitions()
        True
        sage: [[-2,3],[3]] in OrderedMultisetPartitions()
        False
        sage: [[2],[3,3]] in OrderedMultisetPartitions()
        False

    If `n` is specified, it returns the class of ordered multiset partitions of `n`::

        sage: OrderedMultisetPartitions(3)
        Ordered Multiset Partitions of 3
        sage: list(OrderedMultisetPartitions(3))
        [[{1}, {1}, {1}], [{1,2}], [{2}, {1}], [{1}, {2}], [{3}]]
        sage: OrderedMultisetPartitions(3).cardinality()
        5

    The following examples show how to test whether or not an object
    is an ordered multiset partition::

        sage: [[3,2],[2]] in OrderedMultisetPartitions()
        True
        sage: [[3,2],[2]] in OrderedMultisetPartitions(7)
        True
        sage: [[3,2],[2]] in OrderedMultisetPartitions(5)
        False

    The options ``length``, ``min_length``, and ``max_length`` can be used
    to set length constraints on the ordered multiset partitions. For example, the
    ordered multiset partitions of 4 of length equal to, at least, and at most 2 are
    given by::

        sage: OrderedMultisetPartitions(4, length=2).list()
        [[{1}, {3}], [{1}, {1,2}], [{2}, {2}], [{3}, {1}], [{1,2}, {1}]]
        sage: OrderedMultisetPartitions(4, min_length=3).list()
        [[{1}, {1}, {1}, {1}], [{1}, {1}, {2}], [{1}, {2}, {1}], [{2}, {1}, {1}]]
        sage: OrderedMultisetPartitions(4, max_length=2).list()
        [[{1}, {3}], [{1}, {1,2}], [{2}, {2}], [{3}, {1}], [{1,2}, {1}], [{4}], [{1,3}]]

    The option ``alphabet`` constrains which integers appear across all blocks of 
    the ordered multiset partition. For example, the ordered multiset partitions of 4 
    are listed for different choices of alphabet below. Note that ``alphabet`` 
    is allowed to be an integer or an iterable::

        sage: OMPs = OrderedMultisetPartitions
        sage: OMPs(4, alphabet=3).list()
        [[{1}, {1}, {1}, {1}], [{1}, {1}, {2}], [{1}, {2}, {1}], [{1}, {3}], [{1}, {1,2}],
         [{2}, {1}, {1}], [{2}, {2}], [{3}, {1}], [{1,2}, {1}], [{1,3}]]
        sage: OMPs(4, alphabet=3) == OMPs(4, alphabet=[1,2,3])
        True
        sage: OMPs(4, alphabet=[3]).list()
        []
        sage: OMPs(4, alphabet=[1,3]).list()
        [[{1}, {1}, {1}, {1}], [{1}, {3}], [{3}, {1}], [{1,3}]]
        sage: OMPs(4, alphabet=[2]).list()
        [[{2}, {2}]]
        sage: OMPs(4, alphabet=[1,2]).list()
        [[{1}, {1}, {1}, {1}], [{1}, {1}, {2}], [{1}, {2}, {1}], [{1}, {1,2}],
         [{2}, {1}, {1}], [{2}, {2}], [{1,2}, {1}]]
        sage: OMPs(4, alphabet=4).list() == OMPs(4).list()
        True

    The option ``weight`` can also be used to constrain the ordered multiset partitions. 
    It is a refinement of ``alphabet`` in that it specifies the number of times each integer
    appears. For example, the ordered multiset partitions of 4 are listed by weight below::

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
        sage: OrderedMultisetPartitions(4, weight=[4]).list()
        [[{1}, {1}, {1}, {1}]]

    The ``order`` options are coarse versions of ``weight``. For example, the ordered 
    multiset partitions of 4 are listed by order below::

        sage: OrderedMultisetPartitions(4, order=1).list()
        [[{4}]]
        sage: OrderedMultisetPartitions(4, order=2).list()
        [[{1}, {3}], [{2}, {2}], [{3}, {1}], [{1,3}]]
        sage: OrderedMultisetPartitions(4, order=3).list()
        [[{1}, {1}, {2}], [{1}, {2}, {1}], [{1}, {1,2}], [{2}, {1}, {1}], [{1,2}, {1}]]
        sage: OrderedMultisetPartitions(4, order=4).list()
        [[{1}, {1}, {1}, {1}]]
        sage: OrderedMultisetPartitions(4, max_order=2).list()
        [[{1}, {3}], [{2}, {2}], [{3}, {1}], [{4}], [{1,3}]]

    TESTS::

        sage: C = OrderedMultisetPartitions(8, length=3); C.cardinality()
        72
        sage: C == loads(dumps(C))
        True
    """
    @staticmethod
    def __classcall_private__(self, n=None, **kwargs):
        """
        Return the correct parent based upon the input.

        EXAMPLES::

            sage: C = OrderedMultisetPartitions(3)
            sage: C2 = OrderedMultisetPartitions(int(3))
            sage: C is C2
            True
        """
        if n is None:
            is_infinite = True
            if len(kwargs) != 0:
                # TODO: rework? The corresponding slices would be infinite, but otherwise
                #       there's no problem with, e.g., calling `OrderedMultisetPartitions(min_length=3)`
                raise ValueError("Incorrect number of arguments")
            return OrderedMultisetPartitions_all()
        else:
            is_infinite = False
            if len(kwargs) == 0:
                if n in ZZ:
                    return OrderedMultisetPartitions_n(n)
                else:
                    raise ValueError("n must be an integer")
            else:
                if len(kwargs) == 1:
                    ss = ''
                else:
                    ss = 's'
                name = "Ordered Multiset Partitions of the integer %s satisfying constraint%s %s"%(n, ss, ", ".join( ["%s=%s"%(key, kwargs[key]) for key in sorted(kwargs)] ))
                def rename(): return name
                F = Family(apply_strict_multiset_constraints(OrderedMultisetPartitions_n(n), kwargs))
                F._repr_ = rename
                F.constraints = dict(kwargs)

                return F

    def __init__(self, is_infinite=False):
        """
        Initialize ``self``.

        TESTS::

            sage: C = OrderedMultisetPartitions()
            sage: TestSuite(C).run()
            sage: C = OrderedMultisetPartitions(5, weight=[2,0,1])
            sage: C._repr_ == 'Ordered Multiset Partitions of the integer 5 satisfying constraint weight=[2, 0, 1]'
            True
            sage: C.cardinality() == 5
            sage: C = OrderedMultisetPartitions(5, weight=[2,0,1], min_length=4)
            sage: C.list() == []
            True
            sage: C._repr_ == 'Ordered Multiset Partitions of the integer 5 satisfying constraints min_length=4, weight=[2, 0, 1]'
            True
            sage: C = OrderedMultisetPartitions()
            sage: TestSuite(C).run()
        """
        if is_infinite:
            Parent.__init__(self, category=InfiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=FiniteEnumeratedSets())

    Element = OrderedMultisetPartition

    def _element_constructor_(self, lst):
        """
        Construct an element with ``self`` as parent.

        EXAMPLES::

            sage: P = OrderedMultisetPartitions()
            sage: P([[3],[3,1]]) # indirect doctest
            [{3}, {1,3}]
        """
        if len(lst) == 0:
            return self.element_class(self, [])
        if lst[0] in ZZ:
            if 0 in lst:
                return OrderedMultisetPartitions().from_zero_list(list(lst))
            else:
                return OrderedMultisetPartitions().from_list(list(lst))

        lst = map(list,lst)
        if not self._has_valid_blocks(lst):
            raise ValueError("%s not in %s"%(lst, self))
        elt = self.element_class(self, map(Set, lst))
        if elt not in self:
            raise ValueError("%s not in %s"%(elt, self))
        return elt

    def _has_valid_blocks(self, x):
        """
        Blocks should be nonempty sets/lists/tuples of distinct elements.

        TODO: In case parent is not `OrderedMultisetPartitions_n`, then allow entries
        of blocks to be more exotic objects than positive integers.
        """
        for block in x:
            if (not isinstance(block, (list, tuple, set, Set))):
                return False
            if len(block) == 0 or (len(Set(block)) != len(list(block))):
                return False
            if not all([(i in ZZ and i > 0) for i in block]):
                return False
        return True

    def __contains__(self, x):
        """
        TODO: make these work.

            sage: [[2, -1, 'a']] in OrderedMultisetPartitions()
            True
            sage: [[2, -1]] in OrderedMultisetPartitions_n()
            False

        TESTS::

            sage: [[2,1], [1,3]] in OrderedMultisetPartitions()
            True
            sage: [[2,1], [1,3]] in OrderedMultisetPartitions(7)
            True
            sage: [[2,2], [1,3]] in OrderedMultisetPartitions()
            False
            sage: [] in OrderedMultisetPartitions()
            True
            sage: [[2, -1]] in OrderedMultisetPartitions()
            False
        """
        if isinstance(x, OrderedMultisetPartition):
            return True
        elif isinstance(x, (list,tuple)):
            return self._has_valid_blocks(x)
        else:
            return False

    def from_list(self, lst):
        """
        Return an ordered multiset partition from a list of integers,
        each meant to represent a singleton block.

        INPUT:

        - ``lst`` -- an iterable

        OUTPUT:

        - The composition of ``sum(lst)`` whose blocks are the singletons listed in
          ``lst``.

        EXAMPLES::

            sage: OrderedMultisetPartitions().from_list([1,4,8])
            [{1}, {4}, {8}]
            sage: OrderedMultisetPartitions().from_list([1,4,1,8])
            [{1}, {4}, {1}, {8}]
        """
        if isinstance(lst, (list, tuple)) and (0 not in lst) and (lst[0] in ZZ):
            d = [[x] for x in lst]
            return self.element_class(self, d)
        else:
            raise ValueError("Something is wrong: `from_list` expects a list of positive integers, received {}.".format(str(lst)))

    def from_zero_list(self, lst_with_zeros):
        r"""
        Return an ordered multiset partition from a list of nonnegative integers.
        Blocks are separated by zeros. Consecutive zeros are ignored.

        EXAMPLES::

            sage: OrderedMultisetPartitions().from_zero_list([1,2,0,2,0,1,3,4,0,0,4])
            [{1,2}, {2}, {1,3,4}, {4}]
            sage: OrderedMultisetPartitions().from_zero_list([1,2,4])
            [{1,2,4}]
        """
        from_zero_lst = list(lst_with_zeros)
        if not from_zero_lst[-1] == 0:
            from_zero_lst += [0]
        lst = []; block=[]
        for a in from_zero_lst:
            if a == 0 and block !=[]:
                lst.append([b for b in block if b != 0])
                block = []
            else:
                block.append(a)
        if self._has_valid_blocks(lst):
            return self.element_class(self, map(Set, lst))
        else:
            raise ValueError("Ordered Multiset Partitions do not have repeated entries within blocks; received %s."%str(lst))

# Hack to pass ordered multiset partitions with constraints.
def apply_strict_multiset_constraints(n_or_OMPs, kwargs):
    """
    Filter output of ``OrderedMultisetPartitions_n``,
    subject to passed constraints dictionary ``kwargs``.

    Allowable keywords are: 
    length, min_length, max_length, alphabet, weight, order, min_order, max_order.
    """
    if n_or_OMPs in ZZ:
        OMPs = OrderedMultisetPartitions_n(n_or_OMPs)
    else:
        OMPs = n_or_OMPs

    def pass_test(co, (key,tst)):
        if key == 'length':
            return co.length() == tst
        if key == 'min_length':
            return co.length() >= tst
        if key == 'max_length':
            return co.length() <= tst
        if key == 'alphabet':
            if tst in ZZ:
                tst = range(1,tst+1)
            return set(tst).issuperset(set(co.letters()))
        if key == 'weight':
            t = len(tst)
            return co.weight()[:t] == list(tst)
        if key == 'order':
            return co.order() == tst
        if key == 'min_order':
            return co.order() >= tst
        if key == 'max_order':
            return co.order() <= tst
    def passes_tests(co):
        return all([pass_test(co, (key,tst)) for (key,tst) in kwargs.iteritems()])

    for co in OMPs:
        if passes_tests(co):
            yield co

class OrderedMultisetPartitions_all(OrderedMultisetPartitions):
    """
    Class of all ordered multiset partitions.
    """
    def __init__(self):
        """
        Initialize ``self``.

        TESTS::

            sage: C = OrderedMultisetPartitions()
            sage: TestSuite(C).run()
        """
        OrderedMultisetPartitions.__init__(self, True)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(OrderedMultisetPartitions())
            'Ordered Multiset Partitions of non-negative integers'
        """
        return "Ordered Multiset Partitions of non-negative integers"

    def subset(self, size=None):
        """
        Return the set of ordered multiset partitions of the given size.

        EXAMPLES::

            sage: C = OrderedMultisetPartitions()
            sage: C.subset(4)
            Ordered Multiset Partitions of 4
            sage: C.subset(size=3)
            Ordered Multiset Partitions of 3
        """
        if size is None:
            return self
        return OrderedMultisetPartitions(size)

    def __iter__(self):
        """
        Iterate over all ordered multiset partitions.

        TESTS::

            sage: C = OrderedMultisetPartitions()
            sage: it = C.__iter__()
            sage: [next(it) for i in range(12)]
            [[], [{1}], [{1}, {1}], [{2}],
             [{1}, {1}, {1}], [{1}, {2}], [{2}, {1}], [{3}], [{1,2}],
             [{1}, {1}, {1}, {1}], [{1}, {1}, {2}], [{1}, {2}, {1}]]
        """
        n = 0
        while True:
            for c in OrderedMultisetPartitions(n):
                yield self.element_class(self, list(c))  # should this just be `c` ?
            n += 1

class OrderedMultisetPartitions_n(OrderedMultisetPartitions):
    """
    Class of ordered multiset partitions of a fixed `n`.
    """
    @staticmethod
    def __classcall_private__(cls, n):
        """
        Standardize input to ensure a unique representation.

        EXAMPLES::

            sage: C = OrderedMultisetPartitions(5)
            sage: C2 = OrderedMultisetPartitions(int(5))
            sage: C3 = OrderedMultisetPartitions(ZZ(5))
            sage: C is C2
            True
            sage: C is C3
            True
        """
        return super(OrderedMultisetPartitions_n, cls).__classcall__(cls, Integer(n))

    def __init__(self, n):
        """
        TESTS::

            sage: C = OrderedMultisetPartitions(3)
            sage: C == loads(dumps(C))
            True
            sage: TestSuite(C).run()
        """
        self.n = n
        OrderedMultisetPartitions.__init__(self, False)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(OrderedMultisetPartitions(3))
            'Ordered Multiset Partitions of 3'
        """
        return "Ordered Multiset Partitions of %s"%self.n

    def __contains__(self, x):
        """
        TESTS::

            sage: [[2,1,3]] in OrderedMultisetPartitions(6)
            True
            sage: [[2,1],[2]] in OrderedMultisetPartitions(6)
            False
            sage: [] in OrderedMultisetPartitions(0)
            True
            sage: [0] in OrderedMultisetPartitions(0)
            False
            sage: [[]] in OrderedMultisetPartitions(0)
            False
        """
        return x in OrderedMultisetPartitions() and sum(map(sum,x)) == self.n

    def cardinality(self):
        """
        Return the number of ordered multiset partitions of `n`.
        """
        # Dispense with the complex computation for small orders.
        orders = {0:1, 1:1, 2:2, 3:5, 4:11, 5:25}
        if self.n <= 5:
            return orders[self.n]

        # We view an ordered multiset partition as a list of 2-regular integer partitions.
        #
        # The 2-regular partitions have a nice generating function (see OEIS:A000009).
        # Below, we take (products of) coefficients of polynomials to compute cardinality.
        t = var('t')
        def prod(lst): return reduce(lambda a,b: a*b, lst, 1)
        partspoly = prod([1+t**k for k in range(1,self.n+1)]).coefficients()
        def partspoly_coeff(d): return partspoly[d][0]
        deg = 0
        for alpha in Compositions(self.n):
            deg += prod([partspoly_coeff(d) for d in alpha])
        #orders[n] = deg
        return deg

    def an_element(self):
        r"""
        Return a typical ``OrderedMultisetPartition``.
        """
        
    def random_element(self):
        r"""
        Return a random ``OrderedMultisetPartition`` with uniform probability.

        ::TODO::
        Is it really uniform probability?

        This method generates a random composition and then
        then creates new blocks after positions that are not ascents.

        EXAMPLES::

            sage: OrderedMultisetPartitions(5).random_element() # random
            [{1}, {1}, {1}, {1}, {1}]
            sage: OrderedMultisetPartitions(5).random_element() # random
            [{1}, {1,3}]
        """
        C = Compositions(self.n).random_element()
        C = list(C)+[0]
        out = []; block = []
        for i in range(len(C)-1):
            block.append(C[i])
            if C[i] >= C[i+1]:
                out.append(block)
                block = []
        return self.element_class(self, map(Set, out))

    def __iter__(self):
        """
        Iterate over the ordered multiset partitions of `n`.

        The degree d part of ordered multiset partition has all sequences of subsets of NN_+ whose total sum adds up to d.

        TESTS::

            sage: OrderedMultisetPartitions(3).list()
            [[{1}, {1}, {1}], [{1}, {2}], [{2}, {1}], [{3}], [{1,2}]]
            sage: OrderedMultisetPartitions(0).list()
            [[]]
        """
        for alpha in Compositions(self.n):
            for p in cartesian_product([RegularPartitions_n(a, 2) for a in alpha]):
                c = map(sorted, p)
                yield self.element_class(self, map(Set, c))



###############

def is_initial_segment(lst):
    return list(range(max(lst)+1)) == lst

def refine_block(X):
    r"""
    Return the set of all possible refinements of a set `X`.

    A refinement of `X` is a tuple of nonempty subsets whose union is `X`.
    """
    XX = list(X)
    n = len(XX)
    out = []
    for w in Words(range(n),n):
        if is_initial_segment(sorted(w.letters())):
            a = [set([]) for _ in range(max(w)+1)]
            for pos in range(n):
                a[w[pos]].add(XX[pos])
            out.append(a)
    return Set([tuple(map(Set,x)) for x in out])

def split_block(X, k=2):
    r"""
    Return the set of all possible splittings of a set `X`.

    A splitting of `X` is a tuple of (possibly empty) subsets whose union is `X`.
    """
    XX = list(X)
    n = len(XX)
    out = []
    for w in Words(range(k),n):
        a = [set([]) for _ in range(k)]
        for pos in range(n):
            a[w[pos]].add(XX[pos])
        out.append(a)
    return Set([tuple(map(Set,x)) for x in out])


def union_of_sets(list_of_sets):
    """
    Return the union of a list of iterables as a Set object.
    """
    return reduce(lambda a,b: Set(a)|Set(b), list_of_sets, Set([]))
def concatenate(list_of_iters):
    return reduce(lambda a,b: a+b, list_of_iters, [])

## a related function (needs an adjective like "strict")
## not used anywhere above
def unordered_multiset_partitions_n(d):
    tmp = OrderedMultisetPartitions(d)
    bas = set([OrderedMultisetPartition(sorted(k)) for k in tmp])
    return bas