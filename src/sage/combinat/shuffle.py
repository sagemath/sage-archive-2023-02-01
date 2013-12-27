# -*- coding: utf-8 -*-
"""
Shuffle product of iterable objects.

The shuffle product of two sequences of lengths `m` and `n` is a
sum over the `\\binomial{m+n}{n}` ways of interleaving the two sequences.

That could be defined inductively by:

MATH::

    (a_n)_{n \geqslant 0} \Cup (b_m)_{m \geqslant 0} = a_0 \cdot \left((a_n)_{n \geqslant 1} \Cup (b_m)_{m \geqslant 0}
    + b_0 \cdot \left((a_n)_{n \geqslant 0} \Cup (b_m)_{m \geqslant 1}\right)

with `(a_n)` and `(b_m)` two non-empty sequences and if one of them is empty
then the product is equals to the other.

The shuffle product has been introduced by S. Eilenberg & S. Mac Lane in 1953 _[EilLan53].

References:
-----------

.. [EilLan53] On the groups `H(\pi, n)`, I,
    Samuel Eilenberg and
    Saunders Mac Lane,
    1953.

Author:
-------

 - Jean-Baptiste Priez
"""
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.graphs.digraph import DiGraph
from sage.graphs.linearextensions import LinearExtensions
from sage.misc.misc import forall
from sage.misc.cachefunc import cached_function
from sage.rings.arith import binomial
from sage.structure.sage_object import SageObject
import collections, itertools

@cached_function
def Sh(n, m):
    """
    TESTS::

        sage: from sage.combinat.shuffle import Sh
        sage: list(Sh(3,2))
        [[1, 2, 3, 4, 5],
         [1, 2, 4, 3, 5],
         [1, 2, 4, 5, 3],
         [1, 4, 2, 3, 5],
         [1, 4, 2, 5, 3],
         [1, 4, 5, 2, 3],
         [4, 1, 2, 3, 5],
         [4, 1, 2, 5, 3],
         [4, 1, 5, 2, 3],
         [4, 5, 1, 2, 3]]
    """
    # [1,...,n].shifted_shuffle([1,...,m])
    # <=> [1,...,n,n+1,...,n+m].right_permutohedron_interval(
    #     [n+1,...,n+m,1,...,n])
    # <=>
    d = DiGraph()
    d.add_vertices(xrange(1, n + m + 1))
    d.add_edges(itertools.chain(
            itertools.combinations(xrange(1, n + 1), 2),
            itertools.combinations(xrange(n + 1, n + m + 1), 2)
    ))
    return LinearExtensions(d)

class SetShuffleProduct(SageObject):
    """
    The union of all possible shuffle product of two sets of iterable.

    TESTS::

        sage: from sage.combinat.shuffle import SetShuffleProduct
        sage: TestSuite(SetShuffleProduct).run()

    """

    def __init__(self, l1, l2, elem_constructor=None):
        """
        TESTS::

            sage: from sage.combinat.shuffle import SetShuffleProduct
            sage: SetShuffleProduct({(1,2,3), (2,3,4)}, {(5,)})
            sage: list(SetShuffleProduct({(1,2,3), (2,3,4)}, {(5,)}))
            [[2, 3, 4, 5],
             [2, 3, 5, 4],
             [2, 5, 3, 4],
             [5, 2, 3, 4],
             [1, 2, 3, 5],
             [1, 2, 5, 3],
             [1, 5, 2, 3],
             [5, 1, 2, 3]]

        """
        assert(isinstance(l1, collections.Iterable) and
               isinstance(l2, collections.Iterable)
        )
        assert(forall(l1, lambda elem: isinstance(elem, collections.Iterable)))
        assert(forall(l2, lambda elem: isinstance(elem, collections.Iterable)))
        self._l1 = list(l1)
        self._l2 = list(l2)

        if elem_constructor:
            self._element_constructor = elem_constructor
        else:
            try:
                e = iter(l1).next()
                if hasattr(e, "parent") and hasattr(e.parent(), "_element_constructor"):
                    self._set_constructor = e.parent()._element_constructor_
                else:
                    self._set_constructor = list
            except StopIteration:
                self._element_constructor = list

    def _repr_(self):
        """
        TESTS::

            sage: SetShuffleProduct([[1,2],[3,4]], [[1,4]])
            Shuffle set product of : [[1, 2], [3, 4]] and [[1, 4]]

        """
        return "Shuffle set product of : %s and %s" %(self._element_constructor(self._l1),
                                                      self._element_constructor(self._l2))

    def _ascii_art_(self):
        """
        TESTS::

            sage: SetShuffleProduct([[BinaryTree()], [BinaryTree([]), BinaryTree([[],[]])]],
            ....: [[1,4]])
            Shuffle set product of:
            [       [          ] ]
            [       [ o,   o   ] ]
            [       [     / \  ] ]     [ [      ] ]
            [ [  ], [    o   o ] ] and [ [ 1, 4 ] ]

        """
        from sage.misc.ascii_art import ascii_art, ascii_art_list
        return ascii_art("Set shuffle product of:") * \
            (ascii_art_list(self._l1) + ascii_art(" and ") +
             ascii_art_list(self._l2))

    def __iter__(self):
        """
        TESTS::

            sage: from sage.combinat.shuffle import SetShuffleProduct
            sage: list(SetShuffleProduct([[],[]], [[]]))
            [[], []]
            sage: list(SetShuffleProduct([[1,2],[3]], [[4]]))
            [[1, 2, 4], [1, 4, 2], [4, 1, 2], [3, 4], [4, 3]]
            sage: list(SetShuffleProduct([[1,2],[3,4]], [[1,4]], elem_constructor=set))
            [set([1, 2, 4]),
             set([1, 2, 4]),
             set([1, 2, 4]),
             set([1, 2, 4]),
             set([1, 2, 4]),
             set([1, 2, 4]),
             set([1, 3, 4]),
             set([1, 3, 4]),
             set([1, 3, 4]),
             set([1, 3, 4]),
             set([1, 3, 4]),
             set([1, 3, 4])]
        """
        return itertools.chain(
            *(itertools.imap(
                lambda (el1, el2): iter(ShuffleProduct(el1, el2, self._element_constructor)),
                itertools.product(self._l1, self._l2)
            ))
        )

    def cardinality(self):
        """
        TESTS::

            sage: SetShuffleProduct([[1,2],[3,4]], [[1,4]], elem_constructor=set).cardinality()
            12
        """
        def comp_binom((el1, el2)):
            ll1 = len(list(el1)); ll2 = len(list(el2))
            return binomial(ll1 + ll2, ll2)

        return sum(itertools.imap(
            comp_binom,
            itertools.product(self._l1, self._l2)
        ))

class ShuffleProduct(SageObject):
    """
    Shuffle product of two iterable.

    """

    def __init__(self, l1, l2, elem_constructor=None):
        """
        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: ShuffleProduct([1,2,3],[4,5])
            Shuffle product of : [1, 2, 3] and [4, 5]
            sage: list(ShuffleProduct(Word("aa"), Word("bbb"), Word))
            [word: aabbb,
             word: ababb,
             word: abbab,
             word: abbba,
             word: baabb,
             word: babab,
             word: babba,
             word: bbaab,
             word: bbaba,
             word: bbbaa]

        """
        assert(isinstance(l1, collections.Iterable) and
               isinstance(l2, collections.Iterable)
        )
        self._l1 = list(l1)
        self._l2 = list(l2)

        if not elem_constructor:
            if hasattr(l1, "parent") and hasattr(l1.parent(), "_element_constructor"):
                self._element_constructor = l1.parent()._element_constructor_
            else:
                self._element_constructor = list
        else:
            self._element_constructor = elem_constructor

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: ShuffleProduct([1,2,3],[4,5])
            Shuffle product of : [1, 2, 3] and [4, 5]
            sage: B = BinaryTree
            sage: ShuffleProduct([B(), B([[],[]])], [])
            Shuffle product of : [., [[., .], [., .]]] and []
        """
        return "Shuffle product of : %s and %s" % (self._l1, self._l2)

    def _ascii_art_(self):
        """
        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: ascii_art(ShuffleProduct([1,2,3],[4,5]))
            Shuffle product of:
            [         ]     [      ]
            [ 1, 2, 3 ] and [ 4, 5 ]
            sage: B = BinaryTree
            sage: ascii_art(ShuffleProduct([B([]), B([[],[]])],
            ....:   [B([[[],[]],[[],None]])]))
            Shuffle product of:
                             [     __o__   ]
            [          ]     [    /     \  ]
            [ o,   o   ]     [   o       o ]
            [     / \  ]     [  / \     /  ]
            [    o   o ] and [ o   o   o   ]
        """
        from sage.misc.ascii_art import ascii_art, ascii_art_list
        return ascii_art("Shuffle product of:") * \
            (ascii_art_list(self._l1) + ascii_art(" and ") +
             ascii_art_list(self._l2))

    def __iter__(self):
        """
        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: list(ShuffleProduct([1,2,3],[4,5]))
            [[1, 2, 3, 4, 5],
             [1, 2, 4, 3, 5],
             [1, 2, 4, 5, 3],
             [1, 4, 2, 3, 5],
             [1, 4, 2, 5, 3],
             [1, 4, 5, 2, 3],
             [4, 1, 2, 3, 5],
             [4, 1, 2, 5, 3],
             [4, 1, 5, 2, 3],
             [4, 5, 1, 2, 3]]
            sage: B = BinaryTree
            sage: ascii_art(list(ShuffleProduct([B([]), B([[],[]])],
            ....:   [B([[[],[]],[[],None]])])))
            [ [                       ]  [                       ]
            [ [ o,   o        __o__   ]  [ o,     __o__      o   ]  [     __o__    o    o
            [ [     / \      /     \  ]  [       /     \    / \  ]  [    /     \       / \
            [ [    o   o,   o       o ]  [      o       o  o   o ]  [   o       o     o   o
            [ [            / \     /  ]  [     / \     /         ]  [  / \     /
            [ [           o   o   o   ], [    o   o   o  ,       ], [ o   o   o  ,  ,
            <BLANKLINE>
               ]
             ] ]
             ] ]
             ] ]
             ] ]
             ] ]

        """
        n1 = len(self._l1)
        n2 = len(self._l2)
        OC = self._element_constructor
        for sig in Sh(n1, n2):
            yield OC([self._l1[i-1] if i<= n1 else self._l2[i-n1-1] for i in sig])


    def __contains__(self, li):
        """
        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: sh = ShuffleProduct([1,2,3],[4,5])
            sage: list(range(1,6)) in sh
            True
            sage: list(range(1,7)) in sh
            False
            sage: [3,4,5,1,2] in sh
            False
            sage: [1,4,2,5,3] in sh
            True
            sage: [1,4,2,2,5,3] in sh
            False
        """
        l1 = self._l1
        l2 = self._l2
        len_l1 = len(l1)
        len_l2 = len(l2)
        i_l2 = i_l1 = 0
        li = list(li)
        for i, el in enumerate(li):
            if l1[i_l1] == el:
                i_l1 += 1
            elif l2[i_l2] == el:
                i_l2 += 1
            else:
                return False
            if i_l1 == len_l1:
                return li[i + 1:] == l2[i_l2:]
            if i_l2 == len_l2:
                return li[i + 1:] == l1[i_l1:]
        return (i_l1 + 1 == len(l1)) and (i_l2 + 1 == len(l2))

    def cardinality(self):
        """
        TESTS::

            sage: from sage.combinat.shuffle import ShuffleProduct
            sage: ShuffleProduct([3,1,2], [4,2,1,3]).cardinality()
            35
            sage: ShuffleProduct([3,1,2,5,6,4], [4,2,1,3]).cardinality() == binomial(10,4)
            True
        """
        ll1 = len(self._l1)
        ll2 = len(self._l2)
        return binomial(ll1 + ll2, ll1)


#from sage.combinat.map_reduce import SearchForestMapReduce

def parallel_shuffle(l1, l2):
    """
    FIXME:: use SearchForestParallelIterator !
    """
    def children(node):
        (l1, l2, acc) = node
        if len(l1) == 0:
            if len(l2) == 0:
                return
            yield ([], [], acc + l2)
        elif len(l2) == 0:
            yield ([], [], acc + l1)
        else:
            yield (l1[1:], l2, acc + [l1[0]])
            yield (l1, l2[1:], acc + [l2[0]])

    def post_proc(node):
        if len(node[0]) > 0 or len(node[1]) > 0:
            return []
        return [node[2]]

    return None
    """SearchForestMapReduce(
            roots=[(l1, l2, [])],
            children=children,
            map_function=post_proc,
            reduce_init=[]
    ).run()"""