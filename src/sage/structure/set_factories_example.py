r"""
An example of set factory
=========================

The goal of this module is to exemplify the use of set factories. Note
that the code is intentionally kept minimal; many things and in
particular several iterators could be written in a more efficient way.

.. SEEALSO::

    :mod:`.set_factories` for an introduction to set
    factories, their specifications, and examples of their use and
    implementation based on this module.
"""
#*****************************************************************************
#  Copyright (C) 2012 Florent Hivert <florent.hivert at lri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.set_factories import (
    SetFactory, ParentWithSetFactory, TopMostParentPolicy)
from sage.sets.all import DisjointUnionEnumeratedSets
from sage.sets.family import LazyFamily
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.rings.integer import Integer
from sage.misc.lazy_attribute import lazy_attribute

MAX = 5


class XYPairsFactory(SetFactory):
    r"""
    An example of set factory, for sets of pairs of integers.

    .. SEEALSO::

        :mod:`.set_factories` for an introduction to set factories.
    """
    def __call__(self, x=None, y=None, policy=None):
        r"""
        Construct the subset from constraints.

        Consider the set `S` of couples `(x,y)` with `x` and `y` in
        `I:=\{0,1,2,3,4\}`. Returns the subsets of element of `S` satisfying
        some constraints.

        INPUT:

        - ``x=a`` -- where ``a`` is an integer (default to ``None``).
        - ``y=b`` -- where ``b`` is an integer (default to ``None``).
        - ``policy`` -- the policy passed to the created set.

        .. SEEALSO::

            :class:`.set_factories.SetFactoryPolicy`

        EXAMPLES:

        Let us first create the set factory::

            sage: from sage.structure.set_factories_example import XYPairsFactory
            sage: XYPairs = XYPairsFactory()

        One can then use the set factory to construct a set::

            sage: P = XYPairs(); P.list()
            [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (0, 3), (1, 3), (2, 3), (3, 3), (4, 3), (0, 4), (1, 4), (2, 4), (3, 4), (4, 4)]

        .. NOTE::

            This function is actually the ``__call__`` method of
            :class:`XYPairsFactory`.
        """
        if policy is None:
            policy = self._default_policy

        if isinstance(x, (Integer, int)):
            if isinstance(y, (Integer, int)):
                return SingletonPair(x, y, policy)
            return PairsX_(x, policy)
        elif isinstance(y, (Integer, int)):
                return Pairs_Y(y, policy)
        return AllPairs(policy)

    def add_constraints(self, cons, (args, opts)):
        r"""
        Add constraints to the set ``cons`` as per
        :meth:`SetFactory.add_constraints<.set_factories.SetFactory.add_constraints>`.

        This is a crude implementation for the sake of the demonstration which
        should not be taken as an example.

        EXAMPLES::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs.add_constraints((3,None), ((2,), {}))
            Traceback (most recent call last):
            ...
            ValueError: Duplicate value for constraints 'x': was 3 now 2
            sage: XYPairs.add_constraints((), ((2,), {}))
            (2, None)
            sage: XYPairs.add_constraints((), ((2,), {'y':3}))
            (2, 3)
        """
        res = list(cons)
        res += [None] * (2 - len(res))

        def set_args(argss):
            for i, v in enumerate(argss):
                if res[i] is not None and v is not None:
                    raise ValueError("Duplicate value for constraints '{}': "
                                     "was {} now {}".format(['x', 'y'][i],
                                                            res[i], v))
                if v is not None:
                    res[i] = v
        set_args(args)

        def parse_args(x=None, y=None):
            set_args((x, y))
        parse_args(**opts)
        if res == (None, None):
            return ()
        return tuple(res)

    @lazy_attribute
    def _default_policy(self):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairsFactory
            sage: XYPairs = XYPairsFactory()
            sage: XYPairs._default_policy
            Set factory policy for <class 'sage.structure.set_factories_example.XYPair'> with parent AllPairs[=Factory for XY pairs(())]
        """
        return TopMostParentPolicy(self, (), XYPair)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs   # indirect doctest
            Factory for XY pairs
        """
        return "Factory for XY pairs"

XYPairs = XYPairsFactory()
XYPairs.__doc__ = XYPairsFactory.__call__.__doc__


class XYPair(ElementWrapper):
    r"""
    A class for Elements `(x,y)` with `x` and `y` in `\{0,1,2,3,4\}`.

    EXAMPLES::

        sage: from sage.structure.set_factories_example import *
        sage: p = XYPair(Parent(), (0,1)); p
        (0, 1)
        sage: p = XYPair(Parent(), (0,8))
        Traceback (most recent call last):
        ...
        ValueError: numbers must be in range(5)
    """
    def __init__(self, parent, value, check=True):
        """
        TESTS::

            sage: from sage.structure.set_factories_example import *
            sage: P = XYPairs(); p = P.list()[0]
            sage: TestSuite(p).run()
        """
        if check:
            if not isinstance(value, tuple):
                raise ValueError("Value {} must be a tuple".format(value))
            if len(value) != 2:
                raise ValueError("Value must be of length 2")
            if not all(int(x) in range(MAX) for x in value):
                raise ValueError("numbers must be in range({})".format(MAX))
        ElementWrapper.__init__(self, parent, value)


class AllPairs(ParentWithSetFactory, DisjointUnionEnumeratedSets):
    r"""
    This parent shows how one can use set factories together with
    :class:`DisjointUnionEnumeratedSets`.

    TESTS::

        sage: from sage.structure.set_factories_example import XYPairs
        sage: P = XYPairs(); P.list()
        [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (0, 2), (1, 2), (2, 2), (3, 2), (4, 2), (0, 3), (1, 3), (2, 3), (3, 3), (4, 3), (0, 4), (1, 4), (2, 4), (3, 4), (4, 4)]
    """
    def __init__(self, policy):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: TestSuite(XYPairs()).run()
        """
        ParentWithSetFactory.__init__(self, (), policy,
                                      category=FiniteEnumeratedSets())
        DisjointUnionEnumeratedSets.__init__(self,
                                             LazyFamily(range(MAX),
                                                        self._single_pair),
                                             facade=True, keepkey=False,
                                             category=self.category())

    def _single_pair(self, letter):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs()._single_pair(1)
            {(a, 1) | a in range(5)}
        """
        return Pairs_Y(letter, policy=self.facade_policy())

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs()   # indirect doctest
            AllPairs
        """
        return "AllPairs"

    def check_element(self, el, check):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: P = XYPairs()
            sage: P.check_element(P.an_element(), True)
            sage: XYPairs()((7, 0))  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: numbers must be in range(5)
        """
        pass


class PairsX_(ParentWithSetFactory, UniqueRepresentation):
    r"""
    The set of pairs `(x, 0), (x, 1), ..., (x, 4)`.

    TESTS::

        sage: from sage.structure.set_factories_example import XYPairs
        sage: P = XYPairs(0); P.list()
        [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4)]
    """
    def __init__(self, x, policy):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: TestSuite(XYPairs(0)).run()
        """
        self._x = x
        ParentWithSetFactory.__init__(self, (x, None), policy,
                                      category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs(x=1)
            {(1, b) | b in range(5)}
        """
        return "{(%s, b) | b in range(%s)}" % (self._x, MAX)

    def an_element(self):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: P = XYPairs(x=0); P.an_element()
            (0, 0)
        """
        return self._element_constructor_((self._x, 0), check=False)

    def check_element(self, el, check):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: P = XYPairs(x=1)
            sage: P.check_element(P.an_element(), True)
            sage: XYPairs(x=1)((0, 0))  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Wrong first coordinate
        """
        (x, y) = el.value
        if x != self._x:
            raise ValueError("Wrong first coordinate")

    def __iter__(self):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: list(XYPairs(x=1))
            [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4)]
        """
        for i in range(MAX):
            yield self._element_constructor_((self._x, i), check=False)


class Pairs_Y(ParentWithSetFactory, DisjointUnionEnumeratedSets):
    r"""
    The set of pairs `(0, y), (1, y), ..., (4, y)`.

    .. WARNING::

        Put a nice warning _single_pair

    TESTS::

        sage: from sage.structure.set_factories_example import XYPairs
        sage: P = XYPairs(y=1); P.list()
        [(0, 1), (1, 1), (2, 1), (3, 1), (4, 1)]
    """
    def __init__(self, y, policy):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: TestSuite(XYPairs(y=1)).run()
        """
        self._y = y
        ParentWithSetFactory.__init__(self, (None, y), policy,
                                      category=FiniteEnumeratedSets())
        DisjointUnionEnumeratedSets.__init__(
            self, LazyFamily(range(MAX), self._single_pair),
            facade=True, keepkey=False,
            category=self.category())  # TODO remove and fix disjoint union.

    def _repr_(self):
        """
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs(y=1)
            {(a, 1) | a in range(5)}
        """
        return "{(a, %s) | a in range(%s)}" % (self._y, MAX)

    def an_element(self):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs(y=1).an_element()
            (0, 1)
        """
        return self._element_constructor_((0, self._y), check=False)

    def _single_pair(self, letter):
        r"""
        Comment that and put link to documentation caveat....

        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs(y=1)._single_pair(0)
            {(0, 1)}
        """
        return SingletonPair(letter, self._y, policy=self.facade_policy())

    def check_element(self, el, check):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: P = XYPairs(y=1)
            sage: P.check_element(P.an_element(), True)
            sage: XYPairs(y=1)((1, 0))  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: Wrong second coordinate
        """
        (x, y) = el.value
        if y != self._y:
            raise ValueError("Wrong second coordinate")


class SingletonPair(ParentWithSetFactory, UniqueRepresentation):
    r"""
    TESTS::

        sage: from sage.structure.set_factories_example import XYPairs
        sage: P = XYPairs(0,1); P.list()
        [(0, 1)]
    """
    def __init__(self, x, y, policy):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: TestSuite(XYPairs(0,1)).run()
        """
        self._xy = (x, y)
        ParentWithSetFactory.__init__(self, (x, y), policy,
                                      category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs(x=2, y=1)
            {(2, 1)}
        """
        return "{%s}" % (self._xy,)

    def check_element(self, el, check):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: XYPairs(0,1).check_element(XYPairs()((0,1)), True)
            sage: XYPairs(0,1).check_element(XYPairs()((1,0)), True)
            Traceback (most recent call last):
            ...
            ValueError: Wrong coordinate
            sage: XYPairs(0,1)((1,1))
            Traceback (most recent call last):
            ...
            ValueError: Wrong coordinate
        """
        if el.value != self._xy:
            raise ValueError("Wrong coordinate")

    def __iter__(self):
        r"""
        TESTS::

            sage: from sage.structure.set_factories_example import XYPairs
            sage: list(XYPairs(0,1))
            [(0, 1)]
        """
        yield self._element_constructor_(self._xy, check=False)
