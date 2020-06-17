r"""
Backtracking

This library contains generic tools for constructing large sets whose
elements can be enumerated by exploring a search space with a (lazy)
tree or graph structure.

- :class:`GenericBacktracker`: Depth first search through a tree
  described by a ``children`` function, with branch pruning, etc.

Deprecated classes (use :func:`RecursivelyEnumeratedSet` instead):

- :class:`TransitiveIdeal`: Depth first search through a
  graph described by a ``neighbours`` relation.

- :class:`TransitiveIdealGraded`: Breadth first search
  through a graph described by a ``neighbours`` relation.

Deprecation details:

- ``TransitiveIdeal(succ, seeds)`` keeps the same behavior as before
  :trac:`6637` and is now the same as ``RecursivelyEnumeratedSet(seeds,
  succ, structure=None, enumeration='naive')``.

- ``TransitiveIdealGraded(succ, seeds, max_depth)`` keeps the same behavior
  as before :trac:`6637` and is now the same as
  ``RecursivelyEnumeratedSet(seeds, succ, structure=None,
  enumeration='breadth', max_depth=max_depth)``.

.. todo::

    - Deprecate ``TransitiveIdeal`` and ``TransitiveIdealGraded``.

    - Once the deprecation has been there for enough time: delete
      ``TransitiveIdeal`` and ``TransitiveIdealGraded``.

"""
#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
#                     2009 Nicolas M. Thiery <nthiery at users.sf.net>
#                     2010 Nicolas Borie <nicolas.borie at math.u-psud.fr>
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
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.monoids import Monoids
from sage.structure.parent import Parent
from sage.categories.commutative_additive_semigroups import (
        CommutativeAdditiveSemigroups)
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.sets.recursively_enumerated_set import (
    RecursivelyEnumeratedSet_generic, RecursivelyEnumeratedSet_forest)

class GenericBacktracker(object):
    r"""
    A generic backtrack tool for exploring a search space organized as a tree,
    with branch pruning, etc.

    See also :class:`RecursivelyEnumeratedSet_forest` and :class:`TransitiveIdeal` for
    handling simple special cases.
    """
    def __init__(self, initial_data, initial_state):
        r"""
        EXAMPLES::

            sage: from sage.combinat.backtrack import GenericBacktracker
            sage: p = GenericBacktracker([], 1)
            sage: loads(dumps(p))
            <sage.combinat.backtrack.GenericBacktracker object at 0x...>
        """
        self._initial_data = initial_data
        self._initial_state = initial_state

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.permutation import PatternAvoider
            sage: p = PatternAvoider(Permutations(4), [[1,3,2]])
            sage: len(list(p))
            14
        """
        #Initialize the stack of generators with the initial data.
        #The generator in stack[i] is a generator for the i^th level
        #of the search tree.
        stack = []
        stack.append(self._rec(self._initial_data, self._initial_state))

        done = False
        while not done:
            #Try to get the next object in this level
            try:
                obj, state, yld = next(stack[-1])
            except StopIteration:
                #If there are no more, go back up the tree
                #We also need to check if we've exhausted all
                #possibilities
                stack.pop()
                done = len(stack) == 0
                continue

            #If the return state is None, then obj is a leaf
            #of the search tree.  If yld is True, then obj
            #should be yielded.
            if yld is True:
                yield obj
            if state is not None:
                stack.append( self._rec(obj, state) )

class PositiveIntegerSemigroup(UniqueRepresentation, RecursivelyEnumeratedSet_forest):
    r"""
    The commutative additive semigroup of positive integers.

    This class provides an example of algebraic structure which
    inherits from :class:`RecursivelyEnumeratedSet_forest`. It builds the positive
    integers a la Peano, and endows it with its natural commutative
    additive semigroup structure.

    EXAMPLES::

        sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
        sage: PP = PositiveIntegerSemigroup()
        sage: PP.category()
        Join of Category of monoids and Category of commutative additive semigroups and Category of infinite enumerated sets and Category of facade sets
        sage: PP.cardinality()
        +Infinity
        sage: PP.one()
        1
        sage: PP.an_element()
        1
        sage: some_elements = list(PP.some_elements()); some_elements
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100]

    TESTS::

        sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
        sage: PP = PositiveIntegerSemigroup()

    We factor out the long test from the ``TestSuite``::

        sage: TestSuite(PP).run(skip='_test_enumerated_set_contains')
        sage: PP._test_enumerated_set_contains()  # long time
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
        """
        RecursivelyEnumeratedSet_forest.__init__(self, facade = ZZ, category=(InfiniteEnumeratedSets(), CommutativeAdditiveSemigroups(), Monoids()))

    def roots(self):
        r"""
        Return the single root of ``self``.

        EXAMPLES::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
            sage: list(PP.roots())
            [1]
        """
        return [ZZ(1)]

    def children(self, x):
        r"""
        Return the single child ``x+1`` of the integer ``x``

        EXAMPLES::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
            sage: list(PP.children(1))
            [2]
            sage: list(PP.children(42))
            [43]
        """
        return [ZZ(x+1)]

    def one(self):
        r"""
        Return the unit of ``self``.

        EXAMPLES::

            sage: from sage.combinat.backtrack import PositiveIntegerSemigroup
            sage: PP = PositiveIntegerSemigroup()
            sage: PP.one()
            1
        """
        return self.first()

class TransitiveIdeal(RecursivelyEnumeratedSet_generic):
    r"""
    Generic tool for constructing ideals of a relation.

    INPUT:

    - ``relation`` -- a function (or callable) returning a list (or iterable)
    - ``generators`` -- a list (or iterable)

    Returns the set `S` of elements that can be obtained by repeated
    application of ``relation`` on the elements of ``generators``.

    Consider ``relation`` as modeling a directed graph (possibly with
    loops, cycles, or circuits). Then `S` is the ideal generated by
    ``generators`` under this relation.

    Enumerating the elements of `S` is achieved by depth first search
    through the graph. The time complexity is `O(n+m)` where `n` is
    the size of the ideal, and `m` the number of edges in the
    relation. The memory complexity is the depth, that is the maximal
    distance between a generator and an element of `S`.

    See also :class:`RecursivelyEnumeratedSet_forest` and :class:`TransitiveIdealGraded`.

    EXAMPLES::

        sage: from sage.combinat.backtrack import TransitiveIdeal
        sage: [i for i in TransitiveIdeal(lambda i: [i+1] if i<10 else [], [0])]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+1,3)], [0])]
        [0, 1, 2]
        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+2,3)], [0])]
        [0, 2, 1]
        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+2,10)], [0])]
        [0, 2, 4, 6, 8]
        sage: sorted(i for i in TransitiveIdeal(lambda i: [mod(i+3,10),mod(i+5,10)], [0]))
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: sorted(i for i in TransitiveIdeal(lambda i: [mod(i+4,10),mod(i+6,10)], [0]))
        [0, 2, 4, 6, 8]
        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+3,9)], [0,1])]
        [0, 1, 3, 4, 6, 7]

        sage: [p for p in TransitiveIdeal(lambda x:[x],[Permutation([3,1,2,4]), Permutation([2,1,3,4])])]
        [[2, 1, 3, 4], [3, 1, 2, 4]]

    We now illustrate that the enumeration is done lazily, by depth first
    search::

        sage: C = TransitiveIdeal(lambda x: [x-1, x+1], (-10, 0, 10))
        sage: f = C.__iter__()
        sage: [ next(f) for i in range(6) ]
        [0, 1, 2, 3, 4, 5]

    We compute all the permutations of 3::

        sage: sorted(p for p in TransitiveIdeal(attrcall("permutohedron_succ"), [Permutation([1,2,3])]))
        [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]

    We compute all the permutations which are larger than [3,1,2,4],
    [2,1,3,4] in the right permutohedron::

        sage: sorted(p for p in TransitiveIdeal(attrcall("permutohedron_succ"),
        ....:     [Permutation([3,1,2,4]), Permutation([2,1,3,4])]))
        [[2, 1, 3, 4], [2, 1, 4, 3], [2, 3, 1, 4], [2, 3, 4, 1],
         [2, 4, 1, 3], [2, 4, 3, 1], [3, 1, 2, 4], [3, 1, 4, 2],
         [3, 2, 1, 4], [3, 2, 4, 1], [3, 4, 1, 2], [3, 4, 2, 1],
         [4, 2, 1, 3], [4, 2, 3, 1], [4, 3, 1, 2], [4, 3, 2, 1]]

    Using TransitiveIdeal people have been using the ``__contains__``
    method provided from the ``__iter__`` method. We need to make sure that
    this continues to work::

        sage: T = TransitiveIdeal(lambda a:[a+7,a+5], [0])
        sage: 12 in T
        True

    """
    def __init__(self, succ, generators):
        r"""
        TESTS::

            sage: from sage.combinat.backtrack import TransitiveIdeal
            sage: C = TransitiveIdeal(factor, (1, 2, 3))
            sage: C._succ
            <function factor at ...>
            sage: C._generators
            (1, 2, 3)
            sage: loads(dumps(C))   # should test for equality with C, but equality is not implemented
        """
        RecursivelyEnumeratedSet_generic.__init__(self, seeds=generators, successors=succ, enumeration='naive')
        self._generators = self._seeds
        self._succ = self.successors

    def __iter__(self):
        r"""
        Return an iterator on the elements of ``self``.

        TESTS::

            sage: from sage.combinat.backtrack import TransitiveIdeal
            sage: C = TransitiveIdeal(lambda x: [1,2], ())
            sage: list(C) # indirect doctest
            []

            sage: C = TransitiveIdeal(lambda x: [1,2], (1,))
            sage: list(C) # indirect doctest
            [1, 2]

            sage: C = TransitiveIdeal(lambda x: [], (1,2))
            sage: list(C) # indirect doctest
            [1, 2]

        """
        return self.naive_search_iterator()

class TransitiveIdealGraded(RecursivelyEnumeratedSet_generic):
    r"""
    Generic tool for constructing ideals of a relation.

    INPUT:

    - ``relation`` -- a function (or callable) returning a list (or iterable)

    - ``generators`` -- a list (or iterable)

    - ``max_depth`` -- (Default: infinity) Specifies the maximal depth to
      which elements are computed

    Return the set `S` of elements that can be obtained by repeated
    application of ``relation`` on the elements of ``generators``.

    Consider ``relation`` as modeling a directed graph (possibly with
    loops, cycles, or circuits). Then `S` is the ideal generated by
    ``generators`` under this relation.

    Enumerating the elements of `S` is achieved by breadth first search
    through the graph; hence elements are enumerated by increasing
    distance from the generators. The time complexity is `O(n+m)`
    where `n` is the size of the ideal, and `m` the number of edges in
    the relation. The memory complexity is the depth, that is the
    maximal distance between a generator and an element of `S`.

    See also :class:`RecursivelyEnumeratedSet_forest` and :class:`TransitiveIdeal`.

    EXAMPLES::

        sage: from sage.combinat.backtrack import TransitiveIdealGraded
        sage: [i for i in TransitiveIdealGraded(lambda i: [i+1] if i<10 else [], [0])]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    We now illustrate that the enumeration is done lazily, by breadth first search::

        sage: C = TransitiveIdealGraded(lambda x: [x-1, x+1], (-10, 0, 10))
        sage: f = C.__iter__()

    The elements at distance 0 from the generators::

        sage: sorted([ next(f) for i in range(3) ])
        [-10, 0, 10]

    The elements at distance 1 from the generators::

        sage: sorted([ next(f) for i in range(6) ])
        [-11, -9, -1, 1, 9, 11]

    The elements at distance 2 from the generators::

        sage: sorted([ next(f) for i in range(6) ])
        [-12, -8, -2, 2, 8, 12]

    The enumeration order between elements at the same distance is not specified.

    We compute all the permutations which are larger than [3,1,2,4] or
    [2,1,3,4] in the permutohedron::

        sage: sorted(p for p in TransitiveIdealGraded(attrcall("permutohedron_succ"),
        ....:     [Permutation([3,1,2,4]), Permutation([2,1,3,4])]))
        [[2, 1, 3, 4], [2, 1, 4, 3], [2, 3, 1, 4], [2, 3, 4, 1],
         [2, 4, 1, 3], [2, 4, 3, 1], [3, 1, 2, 4], [3, 1, 4, 2],
         [3, 2, 1, 4], [3, 2, 4, 1], [3, 4, 1, 2], [3, 4, 2, 1],
         [4, 2, 1, 3], [4, 2, 3, 1], [4, 3, 1, 2], [4, 3, 2, 1]]
    """
    def __init__(self, succ, generators, max_depth=float("inf")):
        r"""
        TESTS::

            sage: from sage.combinat.backtrack import TransitiveIdealGraded
            sage: C = TransitiveIdealGraded(factor, (1, 2, 3))
            sage: C._succ
            <function factor at ...>
            sage: C._generators
            (1, 2, 3)
            sage: loads(dumps(C))   # should test for equality with C, but equality is not implemented
        """
        RecursivelyEnumeratedSet_generic.__init__(self, seeds=generators, successors=succ, enumeration='breadth', max_depth=max_depth)
        self._generators = self._seeds
        self._succ = self.successors

    def __iter__(self):
        r"""
        Return an iterator on the elements of ``self``.

        TESTS::

            sage: from sage.combinat.backtrack import TransitiveIdealGraded
            sage: C = TransitiveIdealGraded(lambda x: [1,2], ())
            sage: list(C) # indirect doctest
            []

            sage: C = TransitiveIdealGraded(lambda x: [1,2], (1,))
            sage: list(C) # indirect doctest
            [1, 2]

            sage: C = TransitiveIdealGraded(lambda x: [], (1,2))
            sage: list(C) # indirect doctest
            [1, 2]

        ::

            sage: fn = lambda i: [i+1] if i<10 else []
            sage: C = TransitiveIdealGraded(fn, [0], max_depth=1)
            sage: list(C)
            [0, 1]
        """
        return self.breadth_first_search_iterator()

