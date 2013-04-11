"""
Binary trees

This module deals with binary trees as mathematical (in particular immmutable)
objects.

.. NOTE :: If you need the data-structure for example to represent sets or hash
           tables with AVL trees, you should have a look at
           :mod:`sage.misc.sagex_ds`.

**AUTHORS:**

- Florent Hivert (2010-2011): initial implementation.
"""
#*****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.list_clone import ClonableArray
from sage.combinat.abstract_tree import (AbstractClonableTree,
                                         AbstractLabelledClonableTree)
from sage.combinat.ordered_tree import LabelledOrderedTrees
from sage.rings.integer import Integer
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute

class BinaryTree(AbstractClonableTree, ClonableArray):
    """
    The class of binary trees

    Binary trees here mean ordered (a.k.a. plane) binary trees,
    meaning that the children of each node are ordered.

    INPUT:

    - ``children`` -- ``None`` (default) or a list, tuple or iterable of
      length 2 of binary trees or convertible objects. This corresponds to
      the standard recursive definition of a binary tree as either a leaf
      or a pair of binary trees. Syntactic sugar allows leaving out all
      but the outermost calls of the ``BinaryTree()`` constructor, so that,
      e. g., ``BinaryTree([BinaryTree(None),BinaryTree(None)])`` can be
      simplified to ``BinaryTree([None,None])``. It is also allowed to
      abbreviate ``[None, None]`` by ``[]``.

    - ``check`` -- (default to ``True``) whether check for binary should be
      performed or not.

    EXAMPLES::

        sage: BinaryTree()
        .
        sage: BinaryTree(None)
        .
        sage: BinaryTree([])
        [., .]
        sage: BinaryTree([None, None])
        [., .]
        sage: BinaryTree([None, []])
        [., [., .]]
        sage: BinaryTree([[], None])
        [[., .], .]
        sage: BinaryTree("[[], .]")
        [[., .], .]

        sage: BinaryTree([[], None, []])
        Traceback (most recent call last):
        ...
        ValueError: this is not a binary tree

    TESTS::

        sage: t1 = BinaryTree([[None, [[],[[], None]]],[[],[]]])
        sage: t2 = BinaryTree([[[],[]],[]])
        sage: with t1.clone() as t1c:
        ...       t1c[1,1,1] = t2
        sage: t1 == t1c
        False
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        Ensure that binary trees created by the enumerated sets and directly
        are the same and that they are instances of :class:`BinaryTree`

        TESTS::

            sage: from sage.combinat.binary_tree import BinaryTrees_all
            sage: issubclass(BinaryTrees_all().element_class, BinaryTree)
            True
            sage: t0 = BinaryTree([[],[[], None]])
            sage: t0.parent()
            Binary trees
            sage: type(t0)
            <class 'sage.combinat.binary_tree.BinaryTrees_all_with_category.element_class'>

            sage: t1 = BinaryTrees()([[],[[], None]])
            sage: t1.parent() is t0.parent()
            True
            sage: type(t1) is type(t0)
            True

            sage: t1 = BinaryTrees(4)([[],[[], None]])
            sage: t1.parent() is t0.parent()
            True
            sage: type(t1) is type(t0)
            True
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        """
        The automatic parent of the element of this class

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

            sage: BinaryTree._auto_parent
            Binary trees
            sage: BinaryTree([None, None]).parent()
            Binary trees
         """
        return BinaryTrees_all()

    def __init__(self, parent, children = None, check = True):
        """
        TESTS::

            sage: BinaryTree([None, None]).parent()
            Binary trees
            sage: BinaryTree("[., [., [., .]]]")
            [., [., [., .]]]
            sage: BinaryTree("[.,.,.]")
            Traceback (most recent call last):
            ...
            ValueError: this is not a binary tree
            sage: all(BinaryTree(repr(bt)) == bt for i in range(6) for bt in BinaryTrees(i))
            True
        """
        if (type(children) is str):  # if the input is the repr of a binary tree
            children = children.replace(".","None")
            from ast import literal_eval
            children = literal_eval(children)
        if children is None:
            children = []
        elif (children == [] or children == ()
              or isinstance(children, (Integer, int))):
            children = [None, None]
        if (children.__class__ is self.__class__ and
            children.parent() == parent):
            children = list(children)
        else:
            children = [self.__class__(parent, x) for x in children]
        ClonableArray.__init__(self, parent, children, check=check)

    def check(self):
        """
        Checks that ``self`` is a binary tree

        EXAMPLES::

            sage: BinaryTree([[], []])     # indirect doctest
            [[., .], [., .]]
            sage: BinaryTree([[], [], []]) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: this is not a binary tree
            sage: BinaryTree([[]])         # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: this is not a binary tree
        """
        if not (not self or len(self) == 2):
            raise ValueError("this is not a binary tree")

    def _repr_(self):
        """
        TESTS::

            sage: t1 = BinaryTree([[], None]); t1  # indirect doctest
            [[., .], .]
            sage: BinaryTree([[None, t1], None])   # indirect doctest
            [[., [[., .], .]], .]
        """
        if not self:
            return "."
        else:
            return super(BinaryTree, self)._repr_()

    def is_empty(self):
        """
        Returns whether ``self`` is  empty.

        EXAMPLES::

            sage: BinaryTree().is_empty()
            True
            sage: BinaryTree([]).is_empty()
            False
            sage: BinaryTree([[], None]).is_empty()
            False
        """
        return not self

    def graph(self):
        """
        Convert ``self`` to a digraph

        EXAMPLE::

            sage: t1 = BinaryTree([[], None])
            sage: t1.graph()
            Digraph on 5 vertices

            sage: t1 = BinaryTree([[], [[], None]])
            sage: t1.graph()
            Digraph on 9 vertices
            sage: t1.graph().edges()
            [(0, 1, None), (0, 4, None), (1, 2, None), (1, 3, None), (4, 5, None), (4, 8, None), (5, 6, None), (5, 7, None)]
        """
        from sage.graphs.graph import DiGraph
        res = DiGraph()
        def rec(tr, idx):
            if not tr:
                return
            else:
                nbl = 2*tr[0].node_number() + 1
                res.add_edges([[idx,idx+1], [idx,idx+1+nbl]])
                rec(tr[0], idx + 1)
                rec(tr[1], idx + nbl + 1)
        rec(self, 0)
        return res

    def canonical_labelling(self,shift=1):
        """
        Returns a labelled version of ``self``.

        The actual canonical labelling is currently unspecified. However, it
        is guaranteed to have labels in `1...n` where `n` is the number of
        nodes of the tree. Moreover, two (unlabelled) trees compare as equal if
        and only if their canonical labelled trees compare as equal.

        EXAMPLES::

            sage: BinaryTree().canonical_labelling()
            .
            sage: BinaryTree([]).canonical_labelling()
            1[., .]
            sage: BinaryTree([[[], [[], None]], [[], []]]).canonical_labelling()
            5[2[1[., .], 4[3[., .], .]], 7[6[., .], 8[., .]]]
        """
        LTR = self.parent().labelled_trees()
        if self:
            sz0 = self[0].node_number()
            return LTR([self[0].canonical_labelling(shift),
                        self[1].canonical_labelling(shift+1+sz0)],
                       label=shift+sz0)
        else:
            return LTR(None)

    def show(self):
        """
        TESTS::

            sage: t1 = BinaryTree([[], [[], None]])
            sage: t1.show()
        """
        self.graph().show(layout='tree', tree_root=0, tree_orientation="down")

    def make_node(self, child_list = [None, None]):
        """
        Modify ``self`` so that it becomes a node with children ``childlist``

        INPUT:

        - ``child_list`` -- a pair of binary trees (or objects convertible to)

        .. NOTE:: ``self`` must be in a mutable state.

        .. SEEALSO::
            :meth:`make_leaf <sage.combinat.binary_tree.BinaryTree.make_leaf>`

        EXAMPLES::

            sage: t = BinaryTree()
            sage: t.make_node([None, None])
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with t.clone() as t1:
            ...    t1.make_node([None, None])
            sage: t, t1
            (., [., .])
            sage: with t.clone() as t:
            ...    t.make_node([BinaryTree(), BinaryTree(), BinaryTree([])])
            Traceback (most recent call last):
            ...
            ValueError: the list must have length 2
            sage: with t1.clone() as t2:
            ...    t2.make_node([t1, t1])
            sage: with t2.clone() as t3:
            ...    t3.make_node([t1, t2])
            sage: t1, t2, t3
            ([., .], [[., .], [., .]], [[., .], [[., .], [., .]]])
        """
        self._require_mutable()
        child_lst = [self.__class__(self.parent(), x) for x in child_list]
        if not(len(child_lst) == 2):
            raise ValueError("the list must have length 2")
        self.__init__(self.parent(), child_lst, check=False)

    def make_leaf(self):
        """
        Modify ``self`` so that it became a leaf

        .. NOTE:: ``self`` must be in a mutable state.

        .. SEEALSO::
            :meth:`make_node <sage.combinat.binary_tree.BinaryTree.make_node>`

        EXAMPLES::

            sage: t = BinaryTree([None, None])
            sage: t.make_leaf()
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: with t.clone() as t1:
            ...    t1.make_leaf()
            sage: t, t1
            ([., .], .)
        """
        self._require_mutable()
        self.__init__(self.parent(), None)

    def _to_dyck_word_rec(self):
        r"""
        EXAMPLES::

            sage: BinaryTree()._to_dyck_word_rec()
            []
            sage: BinaryTree([])._to_dyck_word_rec()
            [1, 0]
            sage: BinaryTree([[[], [[], None]], [[], []]])._to_dyck_word_rec()
            [1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0]
        """
        if self:
            return ([1]+self[0]._to_dyck_word_rec()+
                    [0]+self[1]._to_dyck_word_rec())
        else:
            return []

    def to_dyck_word(self):
        r"""
        Return the Dyck word associated to ``self``

        The bijection is defined recursively as follows:

        - a leaf is associated to the empty Dyck Word

        - a tree with chidren `l,r` is associated to the Dyck word
          `1 T(l) 0 T(r)` where `T(l)` and `T(r)` are the Dyck words
          associated to `l` and `r`.

        EXAMPLES::

            sage: BinaryTree().to_dyck_word()
            []
            sage: BinaryTree([]).to_dyck_word()
            [1, 0]
            sage: BinaryTree([[[], [[], None]], [[], []]]).to_dyck_word()
            [1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0]
        """
        from sage.combinat.dyck_word import DyckWord
        return DyckWord(self._to_dyck_word_rec())

    def canopee(self):
        """
        Returns the canopee of ``self``

        The *canopee* of a non empty binary tree `T` with `n` internal nodes is
        the list `l` of `0` and `1` of length `n-1` obtained by going along the
        leaves of `T` from left to right except the two extremal ones, writing
        `0` if the leaf is a right leaf and `1` if the leaf is a left leaf.

        EXAMPLES::

            sage: BinaryTree([]).canopee()
            []
            sage: BinaryTree([None, []]).canopee()
            [1]
            sage: BinaryTree([[], None]).canopee()
            [0]
            sage: BinaryTree([[], []]).canopee()
            [0, 1]
            sage: BinaryTree([[[], [[], None]], [[], []]]).canopee()
            [0, 1, 0, 0, 1, 0, 1]

        The number of pairs `(t_1, t_2)` of binary trees of size `n` such that
        the canopee of `t_1` is the complementary of the canopee of `t_2` is
        also the number of Baxter permutations (see [DG]_, see
        also :oeis:`A001181`). We check this in small cases::

            sage: [len([(u,v) for u in BinaryTrees(n) for v in BinaryTrees(n)
            ...       if map(lambda x:1-x, u.canopee()) == v.canopee()])
            ...    for n in range(1, 5)]
            [1, 2, 6, 22]

        Here is a less trivial implementation of this::

            sage: from sage.sets.finite_set_map_cy import fibers
            sage: from sage.misc.all import attrcall
            sage: def baxter(n):
            ...      f = fibers(lambda t: tuple(t.canopee()),
            ...                   BinaryTrees(n))
            ...      return sum(len(f[i])*len(f[tuple(1-x for x in i)])
            ...                 for i in f)
            sage: [baxter(n) for n in range(1, 7)]
            [1, 2, 6, 22, 92, 422]

        TESTS::

            sage: t = BinaryTree().canopee()
            Traceback (most recent call last):
            ...
            ValueError: canopee is only defined for non empty binary trees

        REFERENCES:

            .. [DG] S. Dulucq and O, Guibert. Mots de piles, tableaux
               standards et permutations de Baxter, proceedings of
               Formal Power Series and Algebraic Combinatorics, 1994.
        """
        if not self:
            raise ValueError("canopee is only defined for non empty binary trees")
        res = []
        def add_leaf_rec(tr):
            for i in range(2):
                if tr[i]:
                    add_leaf_rec(tr[i])
                else:
                    res.append(1-i)
        add_leaf_rec(self)
        return res[1:-1]


from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.classcall_metaclass import ClasscallMetaclass

from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.misc.cachefunc import cached_method


# Abstract class to serve as a Factory no instance are created.
class BinaryTrees(UniqueRepresentation, Parent):
    """
    Factory for binary trees.

    INPUT:

    - ``size`` -- (optional) an integer

    OUPUT:

    - the set of all binary trees (of the given ``size`` if specified)

    EXAMPLES::

        sage: BinaryTrees()
        Binary trees

        sage: BinaryTrees(2)
        Binary trees of size 2

    .. NOTE:: this in a factory class whose constructor returns instances of
              subclasses.

    .. NOTE:: the fact that OrderedTrees is a class instead of a simple callable
              is an implementation detail. It could be changed in the future
              and one should not rely on it.
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        """
        TESTS::

            sage: from sage.combinat.binary_tree import BinaryTrees_all, BinaryTrees_size
            sage: isinstance(BinaryTrees(2), BinaryTrees)
            True
            sage: isinstance(BinaryTrees(), BinaryTrees)
            True
            sage: BinaryTrees(2) is BinaryTrees_size(2)
            True
            sage: BinaryTrees(5).cardinality()
            42
            sage: BinaryTrees() is BinaryTrees_all()
            True
        """
        if n is None:
            return BinaryTrees_all()
        else:
            if not (isinstance(n, (Integer, int)) and n >= 0):
                raise ValueError("n must be a non negative integer")
            return BinaryTrees_size(Integer(n))

    @cached_method
    def leaf(self):
        """
        Return a left tree with ``self`` as parent.

        EXAMPLES::

            sage: BinaryTrees().leaf()
            .

        TEST::

            sage: (BinaryTrees().leaf() is
            ...    sage.combinat.binary_tree.BinaryTrees_all().leaf())
            True
        """
        return self(None)

#################################################################
# Enumerated set of all binary trees
#################################################################
class BinaryTrees_all(DisjointUnionEnumeratedSets, BinaryTrees):

    def __init__(self):
        """
        TESTS::

            sage: from sage.combinat.binary_tree import BinaryTrees_all
            sage: B = BinaryTrees_all()
            sage: B.cardinality()
            +Infinity

            sage: it = iter(B)
            sage: (it.next(), it.next(), it.next(), it.next(), it.next())
            (., [., .], [., [., .]], [[., .], .], [., [., [., .]]])
            sage: it.next().parent()
            Binary trees
            sage: B([])
            [., .]

            sage: B is BinaryTrees_all()
            True
            sage: TestSuite(B).run()
            """
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), BinaryTrees_size),
            facade=True, keepkey = False)

    def _repr_(self):
        """
        TEST::

            sage: BinaryTrees()   # indirect doctest
            Binary trees
        """
        return "Binary trees"

    def __contains__(self, x):
        """
        TESTS::

            sage: S = BinaryTrees()
            sage: 1 in S
            False
            sage: S([]) in S
            True
        """
        return isinstance(x, self.element_class)

    def __call__(self, x=None, *args, **keywords):
        """
        Ensure that ``None`` instead of ``0`` is passed by default.

        TESTS::

            sage: B = BinaryTrees()
            sage: B()
            .
        """
        return super(BinaryTrees, self).__call__(x, *args, **keywords)

    def unlabelled_trees(self):
        """
        Returns the set of unlabelled trees associated to ``self``

        EXAMPLES::

            sage: BinaryTrees().unlabelled_trees()
            Binary trees
        """
        return self

    def labelled_trees(self):
        """
        Returns the set of labelled trees associated to ``self``

        EXAMPLES::

            sage: BinaryTrees().labelled_trees()
            Labelled binary trees
        """
        return LabelledBinaryTrees()

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: B = BinaryTrees()
            sage: B._element_constructor_([])
            [., .]
            sage: B([[],[]]) # indirect doctest
            [[., .], [., .]]
            sage: B(None)    # indirect doctest
            .
        """
        return self.element_class(self, *args, **keywords)

    Element = BinaryTree

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from combinat import catalan_number
#################################################################
# Enumerated set of binary trees of a given size
#################################################################
class BinaryTrees_size(BinaryTrees):
    """
    The enumerated sets of binary trees of given size

    TESTS::

        sage: from sage.combinat.binary_tree import BinaryTrees_size
        sage: for i in range(6): TestSuite(BinaryTrees_size(i)).run()
    """
    def __init__(self, size):
        """
        TESTS::

            sage: S = BinaryTrees(3)
            sage: S == loads(dumps(S))
            True

            sage: S is BinaryTrees(3)
            True
        """
        super(BinaryTrees_size, self).__init__(category = FiniteEnumeratedSets())
        self._size = size

    def _repr_(self):
        """
        TESTS::

            sage: BinaryTrees(3)   # indirect doctest
            Binary trees of size 3
        """
        return "Binary trees of size %s"%(self._size)

    def __contains__(self, x):
        """
        TESTS::

            sage: S = BinaryTrees(3)
            sage: 1 in S
            False
            sage: S([[],[]]) in S
            True
        """
        return isinstance(x, self.element_class) and x.node_number() == self._size

    def _an_element_(self):
        """
        TESTS::

            sage: BinaryTrees(0).an_element()  # indirect doctest
            .
        """
        return self.first()

    def cardinality(self):
        """
        The cardinality of ``self``

        This is a Catalan number.

        TESTS::

            sage: BinaryTrees(0).cardinality()
            1
            sage: BinaryTrees(5).cardinality()
            42
        """
        return catalan_number(self._size)

    def __iter__(self):
        """
        A basic generator.

        .. TODO:: could be optimized.

        TESTS::

            sage: BinaryTrees(0).list()
            [.]
            sage: BinaryTrees(1).list()
            [[., .]]
            sage: BinaryTrees(3).list()
            [[., [., [., .]]], [., [[., .], .]], [[., .], [., .]], [[., [., .]], .], [[[., .], .], .]]
        """
        if self._size == 0:
            yield self._element_constructor_()
        else:
            for i in range(0, self._size):
                for lft in self.__class__(i):
                    for rgt in self.__class__(self._size-1-i):
                        yield self._element_constructor_([lft, rgt])

    @lazy_attribute
    def _parent_for(self):
        """
        The parent of the element generated by ``self``

        TESTS::

            sage: S = BinaryTrees(3)
            sage: S._parent_for
            Binary trees
        """
        return BinaryTrees_all()

    @lazy_attribute
    def element_class(self):
        """
        TESTS::

            sage: S = BinaryTrees(3)
            sage: S.element_class
            <class 'sage.combinat.binary_tree.BinaryTrees_all_with_category.element_class'>
            sage: S.first().__class__ == BinaryTrees().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: S = BinaryTrees(0)
            sage: S([])   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: wrong number of nodes
            sage: S(None)   # indirect doctest
            .

            sage: S = BinaryTrees(1)   # indirect doctest
            sage: S([])
            [., .]
        """
        res = self.element_class(self._parent_for, *args, **keywords)
        if res.node_number() != self._size:
            raise ValueError("wrong number of nodes")
        return res



class LabelledBinaryTree(AbstractLabelledClonableTree, BinaryTree):
    """
    The class of labelled binary tree

    EXAMPLE::

        sage: LBT = LabelledBinaryTree
        sage: t1 = LBT([[LBT([], label=2), None], None], label=4); t1
        4[None[2[., .], .], .]
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        Ensure that trees created by the sets and directly are the same and
        that they are instances of :class:`LabelledTree`

        TESTS::

            sage: issubclass(LabelledBinaryTrees().element_class, LabelledBinaryTree)
            True
            sage: t0 = LabelledBinaryTree([[],[[], None]], label = 3)
            sage: t0.parent()
            Labelled binary trees
            sage: type(t0)
            <class 'sage.combinat.binary_tree.LabelledBinaryTrees_with_category.element_class'>
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        """
        The automatic parent of the element of this class

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

            sage: LabelledBinaryTree._auto_parent
            Labelled binary trees
            sage: LabelledBinaryTree([], label = 3).parent()
            Labelled binary trees
         """
        return LabelledBinaryTrees()

    def _repr_(self):
        """
        TESTS::

            sage: LBT = LabelledBinaryTree
            sage: t1 = LBT([[LBT([], label=2), None], None], label=4); t1
            4[None[2[., .], .], .]
            sage: LBT([[],[[], None]], label = 3)   # indirect doctest
            3[None[., .], None[None[., .], .]]
        """
        if not self:
            if self._label is not None:
                return repr(self._label)
            else:
                return "."
        else:
            return "%s%s"%(self._label, self[:])

    def binary_search_insert(self, letter):
        """
        Insert a letter in a binary search tree

        INPUT:

        - ``letter`` -- any object comparable with the label of ``self``

        .. NOTE:: ``self`` is supposed to be a binary search tree. No check is
                  performed.

        EXAMPLES::

            sage: LBT = LabelledBinaryTree
            sage: LBT(None).binary_search_insert(3)
            3[., .]
            sage: LBT([], label = 1).binary_search_insert(3)
            1[., 3[., .]]
            sage: LBT([], label = 3).binary_search_insert(1)
            3[1[., .], .]
            sage: res = LBT(None)
            sage: for i in [3,1,5,2,4,6]:
            ...       res = res.binary_search_insert(i)
            sage: res
            3[1[., 2[., .]], 5[4[., .], 6[., .]]]
        """
        LT = self.parent()._element_constructor_
        if not self:
            return LT([], label = letter)
        else:
            if letter <= self.label():
                fils = self[0].binary_search_insert(letter)
                return LT([fils, self[1]], label=self.label())
            else:
                fils = self[1].binary_search_insert(letter)
                return LT([self[0], fils], label=self.label())

    _UnLabelled = BinaryTree


class LabelledBinaryTrees(LabelledOrderedTrees):
    """
    This is a parent stub to serve as a factory class for trees with various
    labels constraints
    """
    def _repr_(self):
        """
        TESTS::

            sage: LabelledBinaryTrees()   # indirect doctest
            Labelled binary trees
        """
        return "Labelled binary trees"

    def _an_element_(self):
        """
        Returns a labelled binary tree

        EXAMPLE::

            sage: LabelledBinaryTrees().an_element()   # indirect doctest
            toto[42[3[., .], 3[., .]], 5[None[., .], None[., .]]]
        """
        LT = self._element_constructor_
        t  = LT([], label = 3)
        t1 = LT([t,t], label = 42)
        t2  = LT([[], []], label = 5)
        return LT([t1,t2], label = "toto")

    def unlabelled_trees(self):
        """
        Returns the set of unlabelled trees associated to ``self``

        EXAMPLES::

            sage: LabelledBinaryTrees().unlabelled_trees()
            Binary trees

        This is used to compute the shape::

            sage: t = LabelledBinaryTrees().an_element().shape(); t
            [[[., .], [., .]], [[., .], [., .]]]
            sage: t.parent()
            Binary trees

        TESTS::

            sage: t = LabelledBinaryTrees().an_element()
            sage: t.canonical_labelling()
            4[2[1[., .], 3[., .]], 6[5[., .], 7[., .]]]
        """
        return BinaryTrees_all()

    def labelled_trees(self):
        """
        Returns the set of labelled trees associated to ``self``

        EXAMPLES::

            sage: LabelledBinaryTrees().labelled_trees()
            Labelled binary trees
        """
        return self

    Element = LabelledBinaryTree



################################################################
# Interface attempt with species...
#
# Kept here for further reference when species will be improved
################################################################
# from sage.combinat.species.library import (
#     CombinatorialSpecies, EmptySetSpecies, SingletonSpecies)

# BT = CombinatorialSpecies()
# F =  EmptySetSpecies()
# N =  SingletonSpecies()
# BT.define(F+N*(BT*BT))
# # b3 = BT.isotypes(range(3))
# # tr = b3.list()[1]

# def BTsp_to_bintrees(bt):
#     """
#     sage: from sage.combinat.binary_tree import BT, BTsp_to_bintrees
#     sage: BTsp_to_bintrees(BT.isotypes(range(5))[0])
#     [., [., [., [., [., .]]]]]
#     sage: def spls(size):
#     ...    return map(BTsp_to_bintrees, BT.isotypes(range(size)).list())
#     sage: spls(3)
#     [[., [., [., .]]], [., [[., .], .]], [[., .], [., .]], [[., [., .]], .], [[[., .], .], .]]
#     sage: all(spls(i) == BinaryTrees(i).list() for i in range(5))
#     True
#     """
#     if len(bt) == 0:
#         return BinaryTree()
#     elif len(bt) == 2 and len(bt[1]) == 2:
#         return BinaryTree(map(BTsp_to_bintrees, list(bt[1])))
#     else:
#         raise ValueError
