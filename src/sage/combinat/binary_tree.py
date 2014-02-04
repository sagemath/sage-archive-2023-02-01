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
from sage.combinat.combinatorial_map import combinatorial_map

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

    def _ascii_art_( self ):
        r"""
        TESTS::

            sage: ascii_art(BinaryTree())
            <BLANKLINE>
            sage: ascii_art(BinaryTree([]))
            o
            sage: for bt in BinaryTrees(3):
            ....:     print ascii_art(bt)
            o
             \
              o
               \
                o
            o
             \
              o
             /
            o
              o
             / \
            o   o
              o
             /
            o
             \
              o
                o
               /
              o
             /
            o
            sage: ascii_art(BinaryTree([None,[]]))
            o
             \
              o
            sage: ascii_art(BinaryTree([None,[None,[]]]))
            o
             \
              o
               \
                o
            sage: ascii_art(BinaryTree([None,[[],None]]))
            o
             \
              o
             /
            o
            sage: ascii_art(BinaryTree([None,[[[],[]],[]]]))
               o
                \
                _o_
               /   \
              o     o
             / \
            o   o
            sage: ascii_art(BinaryTree([None,[[None,[[],[]]],None]]))
            o
             \
              o
             /
            o
             \
              o
             / \
            o   o
            sage: ascii_art(BinaryTree([[],None]))
              o
             /
            o
            sage: ascii_art(BinaryTree([[[[],None], None],None]))
                  o
                 /
                o
               /
              o
             /
            o
            sage: ascii_art(BinaryTree([[[],[]],None]))
                o
               /
              o
             / \
            o   o
            sage: ascii_art(BinaryTree([[[None,[]],[[[],None],None]], None]))
                   o
                  /
              ___o___
             /       \
            o         o
             \       /
              o     o
                   /
                  o
            sage: ascii_art(BinaryTree([[None,[[],[]]],None]))
              o
             /
            o
             \
              o
             / \
            o   o
            sage: ascii_art(BinaryTree([[],[]]))
              o
             / \
            o   o
            sage: ascii_art(BinaryTree([[],[[],None]]))
              _o_
             /   \
            o     o
                 /
                o
            sage: ascii_art(BinaryTree([[None,[]],[[[],None],None]]))
              ___o___
             /       \
            o         o
             \       /
              o     o
                   /
                  o
            sage: ascii_art(BinaryTree([[[],[]],[[],None]]))
                __o__
               /     \
              o       o
             / \     /
            o   o   o
            sage: ascii_art(BinaryTree([[[],[]],[[],[]]]))
                __o__
               /     \
              o       o
             / \     / \
            o   o   o   o
            sage: ascii_art(BinaryTree([[[[],[]],[[],[]]],[]]))
                    ___o___
                   /       \
                __o__       o
               /     \
              o       o
             / \     / \
            o   o   o   o
            sage: ascii_art(BinaryTree([[],[[[[],[]],[[],[]]],[]]]))
              _____o______
             /            \
            o           ___o___
                       /       \
                    __o__       o
                   /     \
                  o       o
                 / \     / \
                o   o   o   o
        """
        node_to_str = lambda bt: str(bt.label()) if hasattr(bt, "label") else "o"

        if self.is_empty():
            from sage.misc.ascii_art import empty_ascii_art
            return empty_ascii_art

        from sage.misc.ascii_art import AsciiArt
        if self[0].is_empty() and self[1].is_empty():
            bt_repr = AsciiArt( [node_to_str(self)] )
            bt_repr._root = 1
            return bt_repr
        if self[0].is_empty():
            node = node_to_str(self)
            rr_tree = self[1]._ascii_art_()
            if rr_tree._root > 2:
                f_line = " " ** Integer( rr_tree._root - 3 ) + node
                s_line = " " ** Integer( len( node ) + rr_tree._root - 3 ) + "\\"
                t_repr = AsciiArt( [f_line, s_line] ) * rr_tree
                t_repr._root = rr_tree._root - 2
            else:
                f_line = node
                s_line = " " + "\\"
                t_line = " " ** Integer( len( node ) + 1 )
                t_repr = AsciiArt( [f_line, s_line] ) * ( AsciiArt( [t_line] ) + rr_tree )
                t_repr._root = rr_tree._root
            t_repr._baseline = t_repr._h - 1
            return t_repr
        if self[1].is_empty():
            node = node_to_str(self)
            lr_tree = self[0]._ascii_art_()
            f_line = " " ** Integer( lr_tree._root + 1 ) + node
            s_line = " " ** Integer( lr_tree._root ) + "/"
            t_repr = AsciiArt( [f_line, s_line] ) * lr_tree
            t_repr._root = lr_tree._root + 2
            t_repr._baseline = t_repr._h - 1
            return t_repr
        node = node_to_str(self)
        lr_tree = self[0]._ascii_art_()
        rr_tree = self[1]._ascii_art_()
        nb_ = lr_tree._l - lr_tree._root + rr_tree._root - 1
        nb_L = int( nb_ / 2 )
        nb_R = nb_L + ( 1 if nb_ % 2 == 1 else 0 )
        f_line = " " ** Integer( lr_tree._root + 1 ) + "_" ** Integer( nb_L ) + node
        f_line += "_" ** Integer( nb_R )
        s_line = " " ** Integer( lr_tree._root ) + "/" + " " ** Integer( len( node ) + rr_tree._root - 1 + ( lr_tree._l - lr_tree._root ) ) + "\\"
        t_repr = AsciiArt( [f_line, s_line] ) * ( lr_tree + AsciiArt( [" " ** Integer( len( node ) + 2 )] ) + rr_tree )
        t_repr._root = lr_tree._root + nb_L + 2
        t_repr._baseline = t_repr._h - 1
        return t_repr

    def is_empty(self):
        """
        Return whether ``self`` is  empty.

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
        Return a labelled version of ``self``.

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

    def _to_dyck_word_rec(self, usemap="1L0R"):
        r"""
        EXAMPLES::

            sage: BinaryTree()._to_dyck_word_rec()
            []
            sage: BinaryTree([])._to_dyck_word_rec()
            [1, 0]
            sage: BinaryTree([[[], [[], None]], [[], []]])._to_dyck_word_rec()
            [1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0]
            sage: BinaryTree([[None,[]],None])._to_dyck_word_rec("L1R0")
            [1, 1, 0, 0, 1, 0]
            sage: BinaryTree([[[], [[], None]], [[], []]])._to_dyck_word_rec("L1R0")
            [1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0]
        """
        if self:
            w = []
            for l in usemap:
                if l == "L": w += self[0]._to_dyck_word_rec(usemap)
                elif l == "R": w+=self[1]._to_dyck_word_rec(usemap)
                elif l == "1": w+=[1]
                elif l == "0": w+=[0]
            return w
        else:
            return []

    @combinatorial_map(name = "to the Tamari corresponding Dyck path")
    def to_dyck_word_tamari(self):
        r"""
        Return the Dyck word associated with ``self`` in consistency with
        the Tamari order on dyck words and binary trees.

        The bijection is defined recursively as follows:

        - a leaf is associated with an empty Dyck word

        - a tree with children `l,r` is associated with the Dyck word
          `T(l) 1 T(r) 0`

        EXAMPLES::

            sage: BinaryTree().to_dyck_word_tamari()
            []
            sage: BinaryTree([]).to_dyck_word_tamari()
            [1, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word_tamari()
            [1, 1, 0, 0, 1, 0]
            sage: BinaryTree([[[], [[], None]], [[], []]]).to_dyck_word_tamari()
            [1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0]
        """
        from sage.combinat.dyck_word import DyckWord
        return self.to_dyck_word("L1R0")
        
    def tamari_interval(self, other):
        r"""
        Return the Tamari interval between ``self`` and ``other`` as a
        :class:`TamariIntervalPoset`.

        A "Tamari interval" is an interval in the Tamari order
        (:meth:`tamari_greater`).
        
        INPUT:

        - ``other`` -- a binary tree greater or equal to ``self``
          in the Tamari order

        EXAMPLES::

            sage: bt = BinaryTree([[None, [[], None]], None])
            sage: ip = bt.tamari_interval(BinaryTree([None, [[None, []], None]])); ip
            The tamari interval of size 4 induced by relations [(2, 4), (3, 4), (3, 1), (2, 1)]
            sage: ip.lower_binary_tree()
            [[., [[., .], .]], .]
            sage: ip.upper_binary_tree()
            [., [[., [., .]], .]]
            sage: ip.interval_cardinality()
            4
            sage: ip.length_of_maximal_chain()
            3
            sage: list(ip.binary_trees())
            [[., [[., [., .]], .]],
             [[., [., [., .]]], .],
             [., [[[., .], .], .]],
             [[., [[., .], .]], .]]
            sage: bt.tamari_interval(BinaryTree([[None,[]],[]]))
            Traceback (most recent call last):
            ...
            ValueError: The two binary trees are not comparable on the Tamari lattice.

        TESTS:

        Setting ``other`` equal to ``bt`` gives an interval consisting of
        just one element::

            sage: ip = bt.tamari_interval(bt)
            sage: ip
            The tamari interval of size 4 induced by relations [(1, 4), (2, 3), (3, 4), (3, 1), (2, 1)]
            sage: list(ip.binary_trees())
            [[[., [[., .], .]], .]]

        Empty trees work well::

            sage: bt = BinaryTree()
            sage: ip = bt.tamari_interval(bt)
            sage: ip
            The tamari interval of size 0 induced by relations []
            sage: list(ip.binary_trees())
            [.]

        """
        from sage.combinat.interval_posets import TamariIntervalPosets
        return TamariIntervalPosets.from_binary_trees(self, other)

    @combinatorial_map(name="to Dyck paths: up step, left tree, down step, right tree")
    def to_dyck_word(self, usemap="1L0R"):
        r"""
        INPUT:

        - ``usemap`` -- a string, either ``1L0R``, ``1R0L``, ``L1R0``, ``R1L0``

        Return the Dyck word associated with ``self`` using the given map.

        The bijection is defined recursively as follows:

        - a leaf is associated to the empty Dyck Word

        - a tree with children `l,r` is associated with the Dyck word
          described by ``usemap`` where `L` and `R` are respectively the
          Dyck words associated with the `l` and `r`.

        EXAMPLES::

            sage: BinaryTree().to_dyck_word()
            []
            sage: BinaryTree([]).to_dyck_word()
            [1, 0]
            sage: BinaryTree([[[], [[], None]], [[], []]]).to_dyck_word()
            [1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word()
            [1, 1, 0, 1, 0, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word("1R0L")
            [1, 0, 1, 1, 0, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word("L1R0")
            [1, 1, 0, 0, 1, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word("R1L0")
            [1, 1, 0, 1, 0, 0]
            sage: BinaryTree([[None,[]],None]).to_dyck_word("R10L")
            Traceback (most recent call last):
            ...
            ValueError: R10L is not a correct map

        TESTS::

            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt == bt.to_dyck_word().to_binary_tree()
            True
            sage: bt == bt.to_dyck_word("1R0L").to_binary_tree("1R0L")
            True
            sage: bt == bt.to_dyck_word("L1R0").to_binary_tree("L1R0")
            True
            sage: bt == bt.to_dyck_word("R1L0").to_binary_tree("R1L0")
            True
        """
        from sage.combinat.dyck_word import DyckWord
        if usemap not in ["1L0R", "1R0L", "L1R0", "R1L0"]:
            raise ValueError, "%s is not a correct map"%(usemap)
        return DyckWord(self._to_dyck_word_rec(usemap))

    def _to_ordered_tree(self, bijection="left", root=None):
        r"""
        Internal recursive method to obtain an ordered tree from a binary
        tree.

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt._to_ordered_tree()
            [[], [[]]]
            sage: bt._to_ordered_tree(bijection="right")
            [[[]], []]
            sage: bt._to_ordered_tree(bijection="none")
            Traceback (most recent call last):
            ...
            ValueError: the bijection argument should be either left or right
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt._to_ordered_tree()
            [[], [[], []], [[], [[]]]]
            sage: bt._to_ordered_tree(bijection="right")
            [[[[]], [[]]], [[]], []]
        """
        close_root = False
        if(root is None):
            from sage.combinat.ordered_tree import OrderedTree
            root = OrderedTree().clone()
            close_root = True
        if(self):
            left, right = self[0],self[1]
            if(bijection == "left"):
                root = left._to_ordered_tree(bijection=bijection,root=root)
                root.append(right._to_ordered_tree(bijection=bijection,root=None))
            elif(bijection =="right"):
                root.append(left._to_ordered_tree(bijection=bijection, root=None))
                root = right._to_ordered_tree(bijection=bijection,root=root)
            else:
                raise ValueError("the bijection argument should be either left or right")
        if(close_root):
            root.set_immutable()
        return root

    @combinatorial_map(name="To ordered tree, left child = left brother")
    def to_ordered_tree_left_branch(self):
        r"""
        Return an ordered tree of size `n+1` by the following recursive rule:

        - if `x` is the left child of `y`, `x` becomes the left brother
          of `y`
        - if `x` is the right child of `y`, `x` becomes the last child
          of `y`

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt.to_ordered_tree_left_branch()
            [[], [[]]]
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt.to_ordered_tree_left_branch()
            [[], [[], []], [[], [[]]]]
        """
        return self._to_ordered_tree()

    @combinatorial_map(name="To ordered tree, right child = right brother")
    def to_ordered_tree_right_branch(self):
        r"""
        Return an ordered tree of size n+1 by the following recursive rule:

        - if `x` is the right child of `y`, `x` becomes the right brother
          of `y`
        - if `x` is the left child of `y`, `x` becomes the first child
          of `y`

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt.to_ordered_tree_right_branch()
            [[[]], []]
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt.to_ordered_tree_right_branch()
            [[[[]], [[]]], [[]], []]
        """
        return self._to_ordered_tree(bijection="right")

    def _postfix_word(self, left_first = True, start = 1):
        r"""
        Internal recursive method to obtain a postfix canonical read of the
        binary tree.

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt._postfix_word()
            [1, 3, 2]
            sage: bt._postfix_word(left_first=False)
            [3, 1, 2]
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt._postfix_word()
            [1, 3, 4, 2, 6, 8, 7, 5]
            sage: bt._postfix_word(left_first=False)
            [8, 6, 7, 3, 4, 1, 2, 5]
        """
        if not self:
            return []
        left = self[0]._postfix_word(left_first, start)
        label = start + self[0].node_number()
        right = self[1]._postfix_word(left_first, start = label +1)
        if left_first:
            left.extend(right)
            left.append(label)
            return left
        else:
            right.extend(left)
            right.append(label)
            return right

    @combinatorial_map(name="To 312 avoiding permutation")
    def to_312_avoiding_permutation(self):
        r"""
        Return a 312-avoiding permutation corresponding to the binary tree.

        The linear extensions of a binary tree form an interval of the weak
        order called the sylester class of the tree. This permutation is
        the minimal element of this sylvester class.

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt.to_312_avoiding_permutation()
            [1, 3, 2]
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt.to_312_avoiding_permutation()
            [1, 3, 4, 2, 6, 8, 7, 5]

        TESTS::

            sage: bt = BinaryTree([[],[]])
            sage: bt == bt.to_312_avoiding_permutation().binary_search_tree_shape(left_to_right=False)
            True
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt == bt.to_312_avoiding_permutation().binary_search_tree_shape(left_to_right=False)
            True
        """
        from sage.combinat.permutation import Permutation
        return Permutation(self._postfix_word())

    @combinatorial_map(name="To complete tree")
    def as_ordered_tree(self,with_leaves=True):
        r"""
        Return the same tree seen as an ordered tree. By default, leaves
        are transformed into actual nodes.

        EXAMPLES::

            sage: bt = BinaryTree([]); bt
            [., .]
            sage: bt.as_ordered_tree()
            [[], []]
            sage: bt.as_ordered_tree(with_leaves = False)
            []
            sage: bt = bt.canonical_labelling(); bt
            1[., .]
            sage: bt.as_ordered_tree()
            1[None[], None[]]
            sage: bt.as_ordered_tree(with_leaves=False)
            1[]
        """
        if with_leaves:
            children = [child.as_ordered_tree(with_leaves) for child in self]
        else:
            children = [child.as_ordered_tree(with_leaves) for child in self if not child.is_empty()]
        if self in LabelledBinaryTrees():
            from sage.combinat.ordered_tree import LabelledOrderedTree
            return LabelledOrderedTree(children, label = self.label())
        else:
            from sage.combinat.ordered_tree import OrderedTree
            return OrderedTree(children)

    @combinatorial_map(name="To graph")
    def to_undirected_graph(self, with_leaves = False):
        r"""
        Return the undirected graph obtained from the tree nodes and edges.
        Leafs are ignored by default but can set ``with_leaves`` to ``True``
        to obtain the graph of the complete tree.

        EXAMPLES::

            sage: bt = BinaryTree([])
            sage: bt.to_undirected_graph()
            Graph on 1 vertex
            sage: bt.to_undirected_graph(with_leaves=True)
            Graph on 3 vertices

        If the tree is labelled, we use its labelling to label the graph.
        Otherwise, we use the graph canonical labelling which means that
        two different trees can have the same graph.

        EXAMPLES::

            sage: bt = BinaryTree([[],[None,[]]])
            sage: bt.canonical_labelling()
            2[1[., .], 3[., 4[., .]]]
            sage: bt.canonical_labelling().to_undirected_graph().edges()
            [(1, 2, None), (2, 3, None), (3, 4, None)]
            sage: bt.to_undirected_graph().edges()
            [(0, 3, None), (1, 2, None), (2, 3, None)]
            sage: bt.canonical_labelling().to_undirected_graph() == bt.to_undirected_graph()
            False
            sage: BinaryTree([[],[]]).to_undirected_graph() == BinaryTree([[[],None],None]).to_undirected_graph()
            True
        """
        return self.as_ordered_tree(with_leaves).to_undirected_graph()

    @combinatorial_map(name="To poset")
    def to_poset(self, with_leaves = False, root_to_leaf=False):
        r"""
        Return the poset obtained by interpreting the tree as a hasse
        diagram.

        The default orientation is from leaves to root but you can
        pass ``root_to_leaf=True`` to obtain the inverse orientation.

        Leafs are ignored by default but can set ``with_leaves`` to ``True``
        to obtain the poset of the complete tree.

        INPUT:

        - ``with_leaves`` -- boolean, true if leaves should be added to the poset
        - ``root_to_leaf`` -- boolean, true if the poset orientation should
          be from root to leaves. It is false by default.

        EXAMPLES::

            sage: bt = BinaryTree([])
            sage: bt.to_poset()
            Finite poset containing 1 elements
            sage: bt.to_poset(with_leaves=True)
            Finite poset containing 3 elements
            sage: bt.to_poset(with_leaves=True).cover_relations()
            [[0, 2], [1, 2]]
            sage: bt = BinaryTree([])
            sage: bt.to_poset(with_leaves=True,root_to_leaf=True).cover_relations()
            [[0, 1], [0, 2]]

        If the tree is labelled, we use its labelling to label the poset.
        Otherwise, we use the poset canonical labelling::

            sage: bt = BinaryTree([[],[None,[]]]).canonical_labelling()
            sage: bt
            2[1[., .], 3[., 4[., .]]]
            sage: bt.to_poset().cover_relations()
            [[4, 3], [3, 2], [1, 2]]
        """
        return self.as_ordered_tree(with_leaves).to_poset(root_to_leaf)

    @combinatorial_map(name="To 132 avoiding permutation")
    def to_132_avoiding_permutation(self):
        r"""
        Return a 132-avoiding permutation corresponding to the binary tree.

        The linear extensions of a binary tree form an interval of the weak
        order called the sylester class of the tree. This permutation is
        the maximal element of this sylvester class.

        EXAMPLES::

            sage: bt = BinaryTree([[],[]])
            sage: bt.to_132_avoiding_permutation()
            [3, 1, 2]
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt.to_132_avoiding_permutation()
            [8, 6, 7, 3, 4, 1, 2, 5]

        TESTS::

            sage: bt = BinaryTree([[],[]])
            sage: bt == bt.to_132_avoiding_permutation().binary_search_tree_shape(left_to_right=False)
            True
            sage: bt = BinaryTree([[[], [[], None]], [[], []]])
            sage: bt == bt.to_132_avoiding_permutation().binary_search_tree_shape(left_to_right=False)
            True
        """
        from sage.combinat.permutation import Permutation
        return Permutation(self._postfix_word(left_first=False))

    @combinatorial_map(order = 2, name="Left-right symmetry")
    def left_right_symmetry(self):
        r"""
        Return the left-right symmetrized tree of ``self``.

        EXAMPLES::

            sage: BinaryTree().left_right_symmetry()
            .
            sage: BinaryTree([]).left_right_symmetry()
            [., .]
            sage: BinaryTree([[],None]).left_right_symmetry()
            [., [., .]]
            sage: BinaryTree([[None, []],None]).left_right_symmetry()
            [., [[., .], .]]
        """
        if not self:
            return BinaryTree()
        tree = [self[1].left_right_symmetry(),self[0].left_right_symmetry()]
        if(not self in LabelledBinaryTrees()):
            return BinaryTree(tree)
        return LabelledBinaryTree(tree, label = self.label())

    @combinatorial_map(order=2, name="Left border symmetry")
    def left_border_symmetry(self):
        r"""
        Return the tree where a symmetry has been applied recursively on
        all left borders. If a tree is made of three trees `[T_1, T_2,
        T_3]` on its left border, it becomes `[T_3', T_2', T_1']` where
        same symmetry has been applied to `T_1, T_2, T_3`.

        EXAMPLES::

            sage: BinaryTree().left_border_symmetry()
            .
            sage: BinaryTree([]).left_border_symmetry()
            [., .]
            sage: BinaryTree([[None,[]],None]).left_border_symmetry()
            [[., .], [., .]]
            sage: BinaryTree([[None,[None,[]]],None]).left_border_symmetry()
            [[., .], [., [., .]]]
            sage: bt = BinaryTree([[None,[None,[]]],None]).canonical_labelling()
            sage: bt
            4[1[., 2[., 3[., .]]], .]
            sage: bt.left_border_symmetry()
            1[4[., .], 2[., 3[., .]]]
        """
        if not self:
            return BinaryTree()
        border = []
        labelled = self in LabelledBinaryTrees()
        labels = []
        t = self
        while(t):
            border.append(t[1].left_border_symmetry())
            if labelled: labels.append(t.label())
            t = t[0]
        tree = BinaryTree()
        for r in border:
            if labelled:
                tree = LabelledBinaryTree([tree,r],label=labels.pop(0))
            else:
                tree = BinaryTree([tree,r])
        return tree

    def canopee(self):
        """
        Return the canopee of ``self``

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
        Return the set of unlabelled trees associated to ``self``

        EXAMPLES::

            sage: BinaryTrees().unlabelled_trees()
            Binary trees
        """
        return self

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``

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
        Return a labelled binary tree

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
        Return the set of unlabelled trees associated to ``self``

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
        Return the set of labelled trees associated to ``self``

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
