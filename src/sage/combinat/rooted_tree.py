r"""
Rooted (Unordered) Trees

AUTHORS:

- Florent Hivert (2011): initial revision
"""
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.sets_cat import Sets
from sage.combinat.ranker import on_fly
from sage.combinat.abstract_tree import (AbstractClonableTree,
                                         AbstractLabelledClonableTree)
from sage.misc.cachefunc import cached_function
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.list_clone import NormalizedClonableList
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


@cached_function
def number_of_rooted_trees(n):
    r"""
    Return the number of rooted trees with `n` nodes.

    Compute the number `a(n)` of rooted trees with `n` nodes using the
    recursive formula ([SL000081]_):

    .. math::

        a(n+1) = 1/n \sum_{k=1}^{n} ( \sum_{d|k} d a(d) ) a(n-k+1)

    EXAMPLES::

        sage: from sage.combinat.rooted_tree import number_of_rooted_trees
        sage: [number_of_rooted_trees(i) for i in range(10)]
        [0, 1, 1, 2, 4, 9, 20, 48, 115, 286]

    REFERENCES:

    .. [SL000081] Sloane's :oeis:`A000081`
    """
    if n == 0:
        return Integer(0)
    if n == 1:
        return Integer(1)
    n = Integer(n)
    return sum(sum(d * number_of_rooted_trees(d) for d in k.divisors()) *
               number_of_rooted_trees(n - k)
               for k in ZZ.range(1, n)) // (n - 1)


class RootedTree(AbstractClonableTree, NormalizedClonableList):
    r"""
    The class for unordered rooted trees.

    The *unordered rooted trees* are an inductive datatype defined
    as follows: An unordered rooted tree is either a leaf
    (carrying no information), or a multiset of unordered rooted
    trees. In the latter case, the trees that belong to this
    multiset are said to be the *children* of the tree.

    The *labelled rooted trees* (:class:`LabelledRootedTree`)
    form a subclass of this class.

    One can create a tree from any list (or more generally iterable)
    of trees or objects convertible to a tree.

    EXAMPLES::

        sage: RootedTree([])
        []
        sage: RootedTree([[], [[]]])
        [[], [[]]]
        sage: RootedTree([[[]], []])
        [[], [[]]]

    One can also enter a rooted tree by using a simple numerical
    encoding of rooted trees::

        sage: from sage.combinat.abstract_tree import from_hexacode
        sage: RT = RootedTrees()
        sage: from_hexacode('32001010', RT)
        [[[]], [[]], [[], []]]
    """
    # Standard auto-parent trick
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        Ensure that rooted trees created by the enumerated sets and directly
        are the same and that they are instances of :class:`RootedTree`

        TESTS::

            sage: from sage.combinat.rooted_tree import RootedTrees_all
            sage: issubclass(RootedTrees_all().element_class, RootedTree)
            True
            sage: t0 = RootedTree([[],[[]]])
            sage: t0.parent()
            Rooted trees
            sage: type(t0)
            <class 'sage.combinat.rooted_tree.RootedTrees_all_with_category.element_class'>

            sage: t1 = RootedTrees()([[],[[]]])
            sage: t1.parent() is t0.parent()
            True
            sage: type(t1) is type(t0)
            True

            sage: t1 = RootedTrees(4)([[],[[]]])
            sage: t1.parent() is t0.parent()
            True
            sage: type(t1) is type(t0)
            True
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        """
        The automatic parent of the elements of this class

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

            sage: RootedTree._auto_parent
            Rooted trees
            sage: RootedTree([]).parent()
            Rooted trees
         """
        return RootedTrees_all()

    def __init__(self, parent=None, children=[], check=True):
        """
        TESTS::

            sage: RT4 = RootedTrees(4)
            sage: t1 = RT4([[],[[]]])
            sage: TestSuite(t1).run()

        Some bad inputs are refused::

            sage: RT4(69)
            Traceback (most recent call last):
            ...
            TypeError: input (69) is not a valid tree
        """
        try:
            children = list(children)
        except TypeError:
            raise TypeError("input ({}) is not a valid tree".format(children))
        tst = (children.__class__ is self.__class__
               and children.parent() == parent)
        if not tst:
            children = [self.__class__(parent, x) for x in children]
        NormalizedClonableList.__init__(self, parent, children, check=check)

    _unordered_ranker = on_fly()

    def normalize(self):
        r"""
        Normalize ``self``.

        This function is at the core of the implementation of rooted
        (unordered) trees. The underlying structure is provided by
        ordered rooted trees. Every rooted tree is represented by a
        normalized element in the set of its planar embeddings.

        There should be no need to call ``normalize`` directly as it
        is called automatically upon creation and cloning or
        modification (by ``NormalizedClonableList``).

        The normalization has a recursive definition. It means first
        that every sub-tree is itself normalized, and also that
        sub-trees are sorted. Here the sort is performed according to
        the rank function.

        EXAMPLES::

            sage: RT = RootedTree
            sage: RT([[],[[]]]) == RT([[[]],[]])  # indirect doctest
            True
            sage: rt1 = RT([[],[[]]])
            sage: rt2 = RT([[[]],[]])
            sage: rt1 is rt2
            False
            sage: rt1._get_list() is rt2._get_list()
            True
        """
        rank, unrank = self._unordered_ranker
        self._require_mutable()
        for st in self:
            assert st.is_immutable(), "Subtree {} is not normalized".format(st)
        self._get_list().sort(key=rank)

        # ensure unique representation
        self.set_immutable()
        res = unrank(rank(self))
        if self is not res:
            self._set_list(res._get_list())

    def is_empty(self):
        r"""
        Return if ``self`` is the empty tree.

        For rooted trees, returns always ``False``

        .. NOTE::

            This is not the same as ``bool(t)`` which returns whether
            ``t`` has some child or not.

        EXAMPLES::

            sage: t = RootedTrees(4)([[],[[]]])
            sage: t.is_empty()
            False
            sage: bool(t)
            True
        """
        return False

    def graft_list(self, other):
        """
        Return the list of trees obtained by grafting ``other`` on ``self``.

        Here grafting means that one takes the disjoint union of
        ``self`` and ``other``, and add one edge from the root of
        other to a vertex of ``self``. The root of the resulting
        tree is the root of ``self``.

        This is useful for free pre-Lie algebras.

        EXAMPLES::

            sage: RT = RootedTree
            sage: x = RT([])
            sage: y = RT([x, x])
            sage: x.graft_list(x)
            [[[]]]
            sage: y.graft_list(x)
            [[[], [], []], [[], [[]]], [[], [[]]]]
        """
        resu = []
        with self.clone() as t:
            t.append(other)
            resu += [t]
        for i, sub in enumerate(self):
            for new_sub in sub.graft_list(other):
                with self.clone() as t:
                    t[i] = new_sub
                    resu += [t]
        return resu

    def graft_on_root(self, other):
        """
        Return the tree obtained by grafting ``other`` on the root of ``self``.

        Here grafting means that one takes the disjoint union of
        ``self`` and ``other``, and add one edge from the root of
        other to the root of ``self``. The root of the resulting
        tree is the root of ``self``.

        This is useful for free Nap algebras.

        EXAMPLES::

            sage: RT = RootedTree
            sage: x = RT([])
            sage: y = RT([x, x])
            sage: x.graft_on_root(x)
            [[]]
            sage: y.graft_on_root(x)
            [[], [], []]
            sage: x.graft_on_root(y)
            [[[], []]]
        """
        with self.clone() as t:
            t.append(other)
            resu = t
        return resu


class RootedTrees(UniqueRepresentation, Parent):
    """
    Factory class for rooted trees

    INPUT:

    - ``size`` -- (optional) an integer

    OUPUT:

    - the set of all rooted trees (of the given ``size`` if specified)

    EXAMPLES::

        sage: RootedTrees()
        Rooted trees

        sage: RootedTrees(2)
        Rooted trees with 2 nodes
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        """
        TESTS::

            sage: from sage.combinat.rooted_tree import (RootedTrees_all,
            ....:    RootedTrees_size)
            sage: RootedTrees(2) is RootedTrees_size(2)
            True
            sage: RootedTrees(5).cardinality()
            9
            sage: RootedTrees() is RootedTrees_all()
            True

        TESTS::

            sage: RootedTrees(0)
            Traceback (most recent call last):
            ...
            ValueError: n must be a positive integer
        """
        if n is None:
            return RootedTrees_all()
        else:
            if not (isinstance(n, (Integer, int)) and n >= 1):
                raise ValueError("n must be a positive integer")
            return RootedTrees_size(Integer(n))


class RootedTrees_all(DisjointUnionEnumeratedSets, RootedTrees):

    def __init__(self):
        """
        TESTS::

           sage: from sage.combinat.rooted_tree import RootedTrees_all

            sage: sum(x**len(t) for t in
            ....:     set(RootedTree(t) for t in OrderedTrees(6)))
            x^5 + x^4 + 3*x^3 + 6*x^2 + 9*x
            sage: sum(x**len(t) for t in RootedTrees(6))
            x^5 + x^4 + 3*x^3 + 6*x^2 + 9*x

            sage: TestSuite(RootedTrees()).run()
        """
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), RootedTrees_size),
            facade=True, keepkey=False)

    def _repr_(self):
        r"""
        TESTS::

            sage: RootedTrees()   # indirect doctest
            Rooted trees
        """
        return "Rooted trees"

    def __contains__(self, x):
        """
        TESTS::

            sage: S = RootedTrees()
            sage: 1 in S
            False
            sage: S([]) in S
            True
        """
        return isinstance(x, self.element_class)

    def unlabelled_trees(self):
        """
        Return the set of unlabelled trees associated to ``self``.

        EXAMPLES::

            sage: RootedTrees().unlabelled_trees()
            Rooted trees
        """
        return self

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``.

        EXAMPLES::

            sage: RootedTrees().labelled_trees()
            Labelled rooted trees

        As a consequence::

            sage: lb = RootedTrees()([[],[[], []]]).canonical_labelling()
            sage: lb
            1[2[], 3[4[], 5[]]]
            sage: lb.__class__
            <class 'sage.combinat.rooted_tree.LabelledRootedTrees_with_category.element_class'>
            sage: lb.parent()
            Labelled rooted trees
        """
        return LabelledRootedTrees()

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: B = RootedTrees()
            sage: B._element_constructor_([])
            []
            sage: B([[],[]]) # indirect doctest
            [[], []]
        """
        return self.element_class(self, *args, **keywords)

    Element = RootedTree


class RootedTrees_size(RootedTrees):
    """
    The enumerated set of rooted trees with a given number of nodes

    TESTS::

        sage: from sage.combinat.rooted_tree import RootedTrees_size
        sage: for i in range(1, 6): TestSuite(RootedTrees_size(i)).run()
    """
    def __init__(self, n):
        """
        TESTS::

            sage: for i in range(1, 6):
            ....:     TestSuite(RootedTrees(i)).run()
        """
        super(RootedTrees_size, self).__init__(category=FiniteEnumeratedSets())
        self._n = n

    def _repr_(self):
        r"""
        TESTS::

            sage: RootedTrees(4)    # indirect doctest
            Rooted trees with 4 nodes
        """
        return "Rooted trees with {} nodes".format(self._n)

    def __contains__(self, x):
        """
        TESTS::

            sage: S = RootedTrees(3)
            sage: 1 in S
            False
            sage: S([[],[]]) in S
            True
        """
        return isinstance(x, self.element_class) and x.node_number() == self._n

    def _an_element_(self):
        """
        TESTS::

            sage: RootedTrees(4).an_element()  # indirect doctest
            [[[[]]]]
        """
        return self.first()

    def __iter__(self):
        """
        An iterator for ``self``

        This generates the rooted trees of given size. The algorithm
        first picks a partition for the sizes of subtrees, then picks
        appropriate tuples of smaller trees.

        EXAMPLES::

            sage: from sage.combinat.rooted_tree import *
            sage: RootedTrees(1).list()
            [[]]
            sage: RootedTrees(2).list()
            [[[]]]
            sage: RootedTrees(3).list()
            [[[[]]], [[], []]]
        """
        if self._n == 1:
            yield self._element_constructor_([])
        else:
            from sage.combinat.partition import Partitions
            from sage.combinat.multichoose_nk import MultichooseNK
            from sage.combinat.cartesian_product import CartesianProduct
            for part in Partitions(self._n - 1):
                mults = part.to_exp_dict()
                choices = []
                for p, mp in mults.items():
                    lp = self.__class__(p).list()
                    new_choice = MultichooseNK(len(lp), mp).map(
                        lambda l: [lp[i] for i in l]).list()
                    choices.append(new_choice)
                for c in CartesianProduct(*choices):
                    yield self._element_constructor_(sum(c, []))

    def check_element(self, el, check=True):
        r"""
        Check that a given tree actually belongs to ``self``.

        This just checks the number of vertices.

        EXAMPLES::

            sage: RT3 = RootedTrees(3)
            sage: RT3([[],[]])     # indirect doctest
            [[], []]
            sage: RT3([[],[],[]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: wrong number of nodes
        """
        if el.node_number() != self._n:
            raise ValueError("wrong number of nodes")

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: RootedTrees(1).cardinality()
            1
            sage: RootedTrees(3).cardinality()
            2
        """
        return number_of_rooted_trees(self._n)

    @lazy_attribute
    def _parent_for(self):
        """
        The parent of the elements generated by ``self``.

        TESTS::

            sage: S = RootedTrees(3)
            sage: S._parent_for
            Rooted trees
        """
        return RootedTrees_all()

    @lazy_attribute
    def element_class(self):
        """
        TESTS::

            sage: S = RootedTrees(3)
            sage: S.element_class
            <class 'sage.combinat.rooted_tree.RootedTrees_all_with_category.element_class'>
            sage: S.first().__class__ == RootedTrees().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: S = RootedTrees(2)
            sage: S([])   # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: wrong number of nodes
            sage: S([[]])   # indirect doctest
            [[]]

            sage: S = RootedTrees(1)   # indirect doctest
            sage: S([])
            []
        """
        res = self.element_class(self._parent_for, *args, **keywords)
        if res.node_number() != self._n:
            raise ValueError("wrong number of nodes")
        return res


class LabelledRootedTree(AbstractLabelledClonableTree, RootedTree):
    """
    Labelled rooted trees.

    A labelled rooted tree is a rooted tree with a label
    attached at each node.

    More formally:
    The *labelled rooted trees* are an inductive datatype defined
    as follows: A labelled rooted tree is either a leaf endowed
    with a label (which can be any object, including ``None``),
    or a multiset of labelled rooted trees, again endowed with a
    label (which may and may not be of the same type as the labels
    of the trees from the multiset). In the latter case, the trees
    that belong to this multiset are said to be the *children* of
    the tree.

    INPUT:

    - ``children`` -- a list or tuple or more generally any iterable
      of trees or objects convertible to trees

    - ``label`` -- any Sage object (default is ``None``)

    EXAMPLES::

        sage: x = LabelledRootedTree([], label = 3); x
        3[]
        sage: LabelledRootedTree([x, x, x], label = 2)
        2[3[], 3[], 3[]]
        sage: LabelledRootedTree((x, x, x), label = 2)
        2[3[], 3[], 3[]]
        sage: LabelledRootedTree([[],[[], []]], label = 3)
        3[None[], None[None[], None[]]]

    Children are reordered in a session dependant order::

        sage: y = LabelledRootedTree([], label = 5); x
        3[]
        sage: xyy2 = LabelledRootedTree((x, y, y), label = 2); xyy2  #random
        2[3[], 5[], 5[]]
        sage: yxy2 = LabelledRootedTree((y, x, y), label = 2); yxy2  #random
        2[3[], 5[], 5[]]
        sage: xyy2 == yxy2
        True

    TESTS::

        sage: xyy2._get_list() is yxy2._get_list()
        True
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        """
        Ensure that trees created by the sets and directly are the same and
        that they are instances of :class:`LabelledRootedTree`

        TESTS::

            sage: issubclass(LabelledRootedTrees().element_class, LabelledRootedTree)
            True
            sage: t0 = LabelledRootedTree([[],[[], []]], label = 3)
            sage: t0.parent()
            Labelled rooted trees
            sage: type(t0)
            <class 'sage.combinat.rooted_tree.LabelledRootedTrees_with_category.element_class'>
        """
        return cls._auto_parent.element_class(cls._auto_parent, *args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        """
        The automatic parent of the element of this class

        When calling the constructor of an element of this class, one needs a
        parent. This class attribute specifies which parent is used.

        EXAMPLES::

            sage: LabelledRootedTree._auto_parent
            Labelled rooted trees
            sage: LabelledRootedTree([], label = 3).parent()
            Labelled rooted trees
         """
        return LabelledRootedTrees()

    _UnLabelled = RootedTree


class LabelledRootedTrees(UniqueRepresentation, Parent):
    """
    This is a parent stub to serve as a factory class for labelled rooted trees

    EXAMPLES::

        sage: LRT = LabelledRootedTrees(); LRT
        Labelled rooted trees
        sage: x = LRT([], label = 3); x
        3[]
        sage: x.parent() is LRT
        True
        sage: y = LRT([x, x, x], label = 2); y
        2[3[], 3[], 3[]]
        sage: y.parent() is LRT
        True

    .. TODO::

        add the possibility to restrict the labels to a fixed set.
    """
    def __init__(self, category=None):
        """
        TESTS::

            sage: TestSuite(LabelledRootedTrees()).run()
        """
        if category is None:
            category = Sets()
        Parent.__init__(self, category=category)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        TESTS::

            sage: LabelledRootedTrees()   # indirect doctest
            Labelled rooted trees
        """
        return "Labelled rooted trees"

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLE::

            sage: LabelledRootedTrees().cardinality()
            +Infinity
        """
        return Infinity

    def _an_element_(self):
        """
        Return a labelled tree.

        EXAMPLE::

            sage: LabelledRootedTrees().an_element()   # indirect doctest
            alpha[3[], 42[3[], 3[]], 5[None[]]]
        """
        LT = self._element_constructor_
        t = LT([], label=3)
        t1 = LT([t, t], label=42)
        t2 = LT([[]], label=5)
        return LT([t, t1, t2], label="alpha")

    def unlabelled_trees(self):
        """
        Return the set of unlabelled trees associated to ``self``.

        EXAMPLES::

            sage: LabelledRootedTrees().unlabelled_trees()
            Rooted trees
        """
        return RootedTrees_all()

    def labelled_trees(self):
        """
        Return the set of labelled trees associated to ``self``.

        EXAMPLES::

            sage: LabelledRootedTrees().labelled_trees()
            Labelled rooted trees
        """
        return self

    Element = LabelledRootedTree
