r"""
Rooted (Cayley, Unordered) Trees

AUTHORS:

- Florent Hivert (2011): initial revision
"""

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_cat import Sets
from sage.combinat.ranker import on_fly
from sage.combinat.abstract_tree import (
    AbstractClonableTree, AbstractLabelledClonableTree)
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
import sage.structure.set_factories as factories


@cached_function
def number_of_rooted_trees(n):
    """
    Number of rooted trees with `n` nodes

    Compute the number `a(n)` of rooted trees with `n` nodes using the
    recursive formula ([SL000081]_):

    .. math::

        a(n+1) = (1/n) * sum_{k=1..n} ( sum_{d|k} d*a(d) ) * a(n-k+1)

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
    The class for unordered rooted trees

    One can create a tree from any list (or more generally iterable)
    of trees or objects convertible to a tree.

    EXAMPLES::

        sage: RootedTree([])
        []
        sage: RootedTree([[], [[]]])
        [[], [[]]]
        sage: RootedTree([[[]], []])
        [[], [[]]]
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

            sage: t1 = RootedTrees(4)([[],[[]]])
            sage: TestSuite(t1).run()
        """
        tst = (children.__class__ is self.__class__
               and children.parent() == parent)
        if tst:
            children = list(children)
        else:
            children = [self.__class__(parent, x) for x in children]
        NormalizedClonableList.__init__(self, parent, children, check=check)

    _cayley_ranker = on_fly()

    def normalize(self):
        r"""
        Normalize ``self``

        There should be no need to call ``normalize`` directly as it
        is called automatically upon creation and cloning or
        modification.

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
        rank, unrank = self._cayley_ranker
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
        Return if ``self`` is the empty tree

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
        Return the list of trees obtained by grafting ``other`` on ``self``

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
        Return the tree obtained by grafting ``other`` on the root of ``self``

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


class RootedTreesFactory(factories.SetFactory):
    """
    Factory class for rooted trees
    """

    def __call__(self, n=(), policy=None):
        """
        TESTS::

            sage: from sage.combinat.rooted_tree import RootedTrees_all, RootedTrees_size
            sage: RootedTrees(2) is RootedTrees_size(2)
            True
            sage: RootedTrees(5).cardinality()
            9
            sage: RootedTrees() is RootedTrees_all()
            True
        """
        if policy is None:
            policy = self._default_policy

        if n == ():
            return RootedTrees_all(policy)
        elif isinstance(n, (Integer, int)) and n >= 0:
            return RootedTrees_size(n, policy)
        msg = "Do not know how to compute {}({})".format(self, n)
        raise NotImplementedError(msg)

        # try:
        #     from sage.combinat.partition import Partition
        #     lst = Partition(n)
        # except ValueError:
        #     raise NotImplementedError("Don't know how to compute %%(%s)"%(
        #         self, n))
        # else:
        #     return RootedTrees_part(lst, policy)

    def add_constraints(self, cons, (n, st)):
        """
        EXAMPLES::

            sage: RootedTrees.add_constraints((), ((3,), {}))
            (3,)
            sage: RootedTrees.add_constraints((3,), ((), {}))
            (3,)
        """
        return cons + n

    @lazy_attribute
    def _default_policy(self):
        r"""
        TESTS::

            sage: from sage.combinat.rooted_tree import RootedTreesFactory
            sage: RT = RootedTreesFactory()
            sage: RT._default_policy
            Set factory policy for <class 'sage.combinat.rooted_tree.RootedTree'> with parent Rooted trees[=Rooted tree set factory(())]
        """
        return factories.TopMostParentPolicy(self, (), RootedTree)

    def __repr__(self):
        """
        TESTS::

            sage: from sage.combinat.rooted_tree import RootedTreesFactory
            sage: RootedTreesFactory()
            Rooted tree set factory
        """
        return "Rooted tree set factory"

RootedTrees = RootedTreesFactory()


class RootedTrees_all(factories.ParentWithSetFactory,
                      DisjointUnionEnumeratedSets):
    """
    TESTS::

        sage: sum(x^len(t) for t in
        ...       set(RootedTree(t) for t in OrderedTrees(6)))
        x^5 + x^4 + 3*x^3 + 6*x^2 + 9*x
        sage: sum(x^len(t) for t in RootedTrees(6))
        x^5 + x^4 + 3*x^3 + 6*x^2 + 9*x
    """
    @staticmethod
    def __classcall_private__(cls, policy=RootedTrees._default_policy):
        """
        Input normalization

        TESTS::

            sage: from sage.combinat.rooted_tree import RootedTrees_all
            sage: RootedTrees() is RootedTrees_all()
            True
        """
        return super(RootedTrees_all, cls).__classcall__(cls, policy)

    def __init__(self, policy=RootedTrees._default_policy):
        """
        TESTS::

            sage: TestSuite(RootedTrees()).run()
        """
        factories.ParentWithSetFactory.__init__(self, (), policy,
                                                category=InfiniteEnumeratedSets())
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), self._of_size),
            facade=True, keepkey=False,
            category=self.category())

    def _repr_(self):
        r"""
        TESTS::

            sage: RootedTrees()   # indirect doctest
            Rooted trees
        """
        return "Rooted trees"

    def _of_size(self, size):
        r"""
        The sub-enumerated set of trees of a given size

        Passed to :class:DisjointUnionEnumeratedSets

        TESTS::

            sage: RootedTrees()._of_size(4)
            Rooted trees with 4 nodes
        """
        return RootedTrees_size(size, policy=self.facade_policy())

    def check_element(self, el, check):
        r"""
        Check that a given tree actually belongs to ``self``

        See :class:`sage.structure.set_factories.ParentWithSetFactory`

        TESTS::

            sage: RT = RootedTrees()
            sage: RT([[],[]])     # indirect doctest
            [[], []]
        """
        pass

    def unlabelled_trees(self):
        """
        Returns the set of unlabelled trees associated to ``self``

        EXAMPLES::

            sage: RootedTrees().unlabelled_trees()
            Rooted trees
        """
        return self

    def labelled_trees(self):
        """
        Returns the set of labelled trees associated to ``self``

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


class RootedTrees_size(factories.ParentWithSetFactory, UniqueRepresentation):
    """
    The enumerated set of rooted trees with a given number of nodes

    EXAMPLES::

    """
    @staticmethod
    def __classcall_private__(cls, n, policy=RootedTrees._default_policy):
        """
        Input normalization

        TESTS::

            sage: from sage.combinat.rooted_tree import RootedTrees_size
            sage: RootedTrees(4) is RootedTrees_size(4)
            True
            sage: RootedTrees_size(4) is RootedTrees_size(int(4))
            True
        """
        return super(RootedTrees_size, cls).__classcall__(
            cls, Integer(n), policy)

    def __init__(self, n, policy):
        """
        TESTS::

            sage: for i in range(0, 6):
            ....:     TestSuite(RootedTrees(i)).run()
        """
        self._n = n
        factories.ParentWithSetFactory.__init__(self, (n,), policy,
                                                category=FiniteEnumeratedSets())

    def _repr_(self):
        r"""
        TESTS::

            sage: RootedTrees(4)    # indirect doctest
            Rooted trees with 4 nodes
        """
        return "Rooted trees with {} nodes".format(self._n)

    def __iter__(self):
        """
        An iterator for ``self``

        EXAMPLES::

            sage: from sage.combinat.rooted_tree import *
            sage: RootedTrees(0).list()
            []
            sage: RootedTrees(1).list()
            [[]]
            sage: RootedTrees(2).list()
            [[[]]]
            sage: RootedTrees(3).list()
            [[[[]]], [[], []]]
        """
        if self._n == 0:
            pass
        elif self._n == 1:
            yield self._element_constructor_([])
        else:
            from sage.combinat.partition import Partitions
            from sage.combinat.multichoose_nk import MultichooseNK
            from sage.combinat.cartesian_product import CartesianProduct
            for part in Partitions(self._n - 1):
                mults = part.to_exp_dict()
                choices = []
                for p, mp in mults.items():
                    lp = self.__class__(p, self.policy()).list()
                    new_choice = MultichooseNK(len(lp), mp).map(
                        lambda l: [lp[i] for i in l]).list()
                    choices.append(new_choice)
                for c in CartesianProduct(*choices):
                    yield self._element_constructor_(sum(c, []))

    def check_element(self, el, check=True):
        r"""
        Check that a given tree actually belongs to ``self``

        See :class:`sage.structure.set_factories.ParentWithSetFactory`

        EXAMPLES::

            sage: RT3 = RootedTrees(3)
            sage: RT3([[],[]])     # indirect doctest
            [[], []]
            sage: RT3([[],[],[]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Wrong number of nodes
        """
        if el.node_number() != self._n:
            raise ValueError("Wrong number of nodes")

    def cardinality(self):
        r"""
        Return the cardinality of ``self``

        EXAMPLES::

            sage: RootedTrees(1).cardinality()
            1
            sage: RootedTrees(3).cardinality()
            2
            sage: RootedTrees(0).cardinality()
            0
        """
        return number_of_rooted_trees(self._n)


class LabelledRootedTree(AbstractLabelledClonableTree, RootedTree):
    """
    Labelled rooted trees

    A labelled rooted tree is a rooted tree with a label attached at each node

    INPUT:

    - ``children`` -- a list or tuple or more generally any iterable
                      of trees or object convertible to trees
    - ``label`` -- any Sage object default to ``None``

    EXAMPLES::

        sage: x = LabelledRootedTree([], label = 3); x
        3[]
        sage: LabelledRootedTree([x, x, x], label = 2)
        2[3[], 3[], 3[]]
        sage: LabelledRootedTree((x, x, x), label = 2)
        2[3[], 3[], 3[]]
        sage: LabelledRootedTree([[],[[], []]], label = 3)
        3[None[], None[None[], None[]]]

    Children are reordered in a session dependant order:

        sage: y = LabelledRootedTree([], label = 5); x
        3[]
        sage: xyy2 = LabelledRootedTree((x, y, y), label = 2); xyy2 #random
        2[3[], 5[], 5[]]
        sage: yxy2 = LabelledRootedTree((y, x, y), label = 2); yxy2 #random
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
    This is a parent stub to serve as a factory class for trees with various
    labels constraints

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
        TESTS::

            sage: LabelledRootedTrees()   # indirect doctest
            Labelled rooted trees
        """
        return "Labelled rooted trees"

    def cardinality(self):
        """
        Returns the cardinality of `self`

        EXAMPLE::

            sage: LabelledRootedTrees().cardinality()
            +Infinity
        """
        return Infinity

    def _an_element_(self):
        """
        Returns a labelled tree

        EXAMPLE::

            sage: LabelledRootedTrees().an_element()   # indirect doctest
            toto[3[], 42[3[], 3[]], 5[None[]]]
        """
        LT = self._element_constructor_
        t = LT([], label=3)
        t1 = LT([t, t], label=42)
        t2 = LT([[]], label=5)
        return LT([t, t1, t2], label="toto")

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: T = LabelledRootedTrees()
            sage: T([], label=2)     # indirect doctest
            2[]
        """
        return self.element_class(self, *args, **keywords)

    def unlabelled_trees(self):
        """
        Returns the set of unlabelled trees associated to ``self``

        EXAMPLES::

            sage: LabelledRootedTrees().unlabelled_trees()
            Rooted trees
        """
        return RootedTrees_all()

    def labelled_trees(self):
        """
        Returns the set of labelled trees associated to ``self``

        EXAMPLES::

            sage: LabelledRootedTrees().labelled_trees()
            Labelled rooted trees
        """
        return self

    Element = LabelledRootedTree
