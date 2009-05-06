"""
Backtracking

This library contains generic tools for constructing large sets whose
elements can be enumerated by exploring a search space with a (lazy)
tree or graph structure.

 - SearchForest:
   Depth first search through a tree descrived by a `child` function
 - GenericBacktracker:
   Depth first search through a tree descrived by a `child` function, with branch pruning, ...
 - TransitiveIdeal:
   Depth first search through a graph described by a `neighbours` relation
 - TransitiveIdealGraded:
   Breath first search through a graph described by a `neighbours` relation

Todo: find a good and consistent naming scheme!!! Do we want to
emphasize the underlying graph/tree structure? The branch&bound
aspect? The transitive closure of a relation point of view?

Todo: do we want TransitiveIdeal(relation, generators) or TransitiveIdeal(generators, relation)?
The code needs to be standardized once the choice is done.

"""
#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
#                     2009 Nicolas M. Thiery <nthiery at users.sf.net>
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
class GenericBacktracker(object):
    """
    A generic backtrack tool for exploring a search space organized as
    a tree, with branch pruning, ...

    See also ``SearchForest`` and ``TransitiveIdeal`` for handling
    simple special cases.
    """

    def __init__(self, initial_data, initial_state):
        """
        EXAMPLES::

            sage: from sage.combinat.backtrack import GenericBacktracker
            sage: p = GenericBacktracker([], 1)
            sage: loads(dumps(p))
            <sage.combinat.backtrack.GenericBacktracker object at 0x...>
        """
        self._initial_data = initial_data
        self._initial_state = initial_state

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.permutation import PatternAvoider
            sage: p = PatternAvoider(4, [[1,3,2]])
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
                obj, state, yld = stack[-1].next()
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


def search_forest_iterator(roots, childs):
    """
    INPUT:

     - ``roots``: a list (or iterable)

     - ``childs``: a function returning a list (or iterable)

    Returns an iterator on the nodes of the forest having the given
    roots, and where ``child(x)`` returns the childs of the node ``x``
    of the forest.

    EXAMPLES::

        sage: from sage.combinat.backtrack import search_forest_iterator
        sage: list(search_forest_iterator([[]], lambda l: [l+[0], l+[1]] if len(l) < 3 else []))
        [[], [0], [0, 0], [0, 0, 0], [0, 0, 1], [0, 1], [0, 1, 0], [0, 1, 1], [1], [1, 0], [1, 0, 0], [1, 0, 1], [1, 1], [1, 1, 0], [1, 1, 1]]

    """

    #Invariant: stack[i] contains an iterator for the siblings of the i-th node of the current branch
    stack = [iter(roots)]
    while len(stack) > 0:
        # Try to get the next node at this depth
        try:
            node = stack[-1].next()
        except StopIteration:
            #If there are no more, go back up the tree
            # We also need to check if we've exhausted all
            # possibilities
            stack.pop()
            continue

        yield node
        stack.append( iter(childs(node)) )

from sage.combinat.combinat import CombinatorialClass
class SearchForest(CombinatorialClass):
    """
    INPUT::

     - ``roots``: a list (or iterable)

     - ``childs``: a function returning a list (or iterable)

    Returns the set of nodes of the forest having the given roots, and
    where ``child(x)`` returns the childs of the node ``x`` of the forest.

    See also ``GenericBacktracker``, ``TransitiveIdeal``, and ``TransitiveIdealGraded``.

    EXAMPLES::

        sage: list(SearchForest([[]], lambda l: [l+[0], l+[1]] if len(l) < 3 else []))
        [[], [0], [0, 0], [0, 0, 0], [0, 0, 1], [0, 1], [0, 1, 0], [0, 1, 1], [1], [1, 0], [1, 0, 0], [1, 0, 1], [1, 1], [1, 1, 0], [1, 1, 1]]
    """
    def __init__(self, roots, childs):
        """
        TESTS::

            sage: C = SearchForest((1,), lambda x: [x+1])
            sage: C._roots
            (1,)
            sage: C._childs
            <function <lambda> at ...>

        """
        self._roots = roots
        self._childs = childs

    def _repr_(self):
        """
        TESTS::
            sage: SearchForest((1,), lambda x: [x+1])	# Todo: improve!
            An enumerated set
        """
        return "An enumerated set"

    def __iter__(self):
        """
        Returns an iterator on the elements of self.

        EXAMPLES::

            sage: def succ(l):
            ...        return [l+[0], l+[1]]
            ...
            sage: C = SearchForest(([],), succ)
            sage: f = C.__iter__()
            sage: f.next()
            []
            sage: f.next()
            [0]
            sage: f.next()
            [0, 0]

            sage: import __main__
            sage: __main__.succ = succ # just because succ has been defined interactively
            sage: loads(dumps(C))
            An enumerated set
        """
        return search_forest_iterator(self._roots, self._childs)

class TransitiveIdeal():
    """
    Generic tool for constructing ideals of a relation.

    INPUT::

    - ``relation``: a function (or callable) returning a list (or iterable)

    - ``generators``: a list (or iterable)

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

    See also ``SearchForest` and ``TransitiveIdealGraded``.

    EXAMPLES::

        sage: [i for i in TransitiveIdeal(lambda i: [i+1] if i<10 else [], [0])]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+1,3)], [0])]
        [0, 1, 2]
        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+2,3)], [0])]
        [0, 2, 1]
        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+2,10)], [0])]
        [0, 2, 4, 6, 8]
        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+3,10),mod(i+5,10)], [0])]
        [0, 3, 8, 1, 4, 5, 6, 7, 9, 2]
        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+4,10),mod(i+6,10)], [0])]
        [0, 4, 8, 2, 6]
        sage: [i for i in TransitiveIdeal(lambda i: [mod(i+3,9)], [0,1])]
        [0, 1, 3, 4, 6, 7]

        sage: [p for p in TransitiveIdeal(lambda x:[x],[Permutation([3,1,2,4]), Permutation([2,1,3,4])])]
        [[2, 1, 3, 4], [3, 1, 2, 4]]

    We now illustrates that the enumeration is done lazily, by depth
    first search::

        sage: C = TransitiveIdeal(lambda x: [x-1, x+1], (-10, 0, 10))
        sage: f = C.__iter__()
        sage: [ f.next() for i in range(6) ]
        [0, 1, 2, 3, 4, 5]

    We compute all the permutations of 3::

        sage: [p for p in TransitiveIdeal(attrcall("permutohedron_succ"), [Permutation([1,2,3])])]
        [[1, 2, 3], [2, 1, 3], [1, 3, 2], [2, 3, 1], [3, 2, 1], [3, 1, 2]]

    We compute all the permutations which are larger than [3,1,2,4],
    [2,1,3,4] in the right permutohedron::

        sage: [p for p in TransitiveIdeal(attrcall("permutohedron_succ"), [Permutation([3,1,2,4]), Permutation([2,1,3,4])])]
        [[2, 1, 3, 4], [2, 1, 4, 3], [2, 4, 1, 3], [4, 2, 1, 3], [4, 2, 3, 1], [4, 3, 2, 1], [3, 1, 2, 4], [2, 4, 3, 1], [3, 2, 1, 4], [2, 3, 1, 4], [2, 3, 4, 1], [3, 2, 4, 1], [3, 1, 4, 2], [3, 4, 2, 1], [3, 4, 1, 2], [4, 3, 1, 2]]

    """
    def __init__(self, succ, generators):
        """
        TESTS::

            sage: C = TransitiveIdeal(factor, (1, 2, 3))
            sage: C._succ
            <function factor at ...>
            sage: C._generators
            (1, 2, 3)
            sage: loads(dumps(C))   # should test for equality with C, but equality is not implemented
        """
        self._succ = succ
        self._generators = generators

    def __iter__(self):
        """
        Returns an iterator on the elements of self.

        TESTS::

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
        known = set(self._generators)
        todo = known.copy()
        while len(todo) > 0:
            x = todo.pop()
            yield x
            for y in self._succ(x):
                if y == None or y in known:
                    continue
                todo.add(y)
                known.add(y)
        return


class TransitiveIdealGraded(TransitiveIdeal):
    """
    Generic tool for constructing ideals of a relation.

    INPUT::

    - ``relation``: a function (or callable) returning a list (or iterable)

    - ``generators``: a list (or iterable)

    Returns the set `S` of elements that can be obtained by repeated
    application of ``relation`` on the elements of ``generators``.

    Consider ``relation`` as modeling a directed graph (possibly with
    loops, cycles, or circuits). Then `S` is the ideal generated by
    ``generators`` under this relation.

    Enumerating the elements of `S` is achieved by breath first search
    through the graph; hence elements are enumerated by increasing
    distance from the generators. The time complexity is `O(n+m)`
    where `n` is the size of the ideal, and `m` the number of edges in
    the relation. The memory complexity is the depth, that is the
    maximal distance between a generator and an element of `S`.

    See also ``SearchForest` and ``TransitiveIdeal``.

    EXAMPLES::

        sage: [i for i in TransitiveIdealGraded(lambda i: [i+1] if i<10 else [], [0])]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    We now illustrates that the enumeration is done lazily, by breath first search::

        sage: C = TransitiveIdealGraded(lambda x: [x-1, x+1], (-10, 0, 10))
        sage: f = C.__iter__()

    The elements at distance 0 from the generators::

        sage: sorted([ f.next() for i in range(3) ])
        [-10, 0, 10]

    The elements at distance 1 from the generators::

        sage: sorted([ f.next() for i in range(6) ])
        [-11, -9, -1, 1, 9, 11]

    The elements at distance 2 from the generators::

        sage: sorted([ f.next() for i in range(6) ])
        [-12, -8, -2, 2, 8, 12]

    The enumeration order between elements at the same distance is not specified.

    We compute all the permutations which are larger than [3,1,2,4] or
    [2,1,3,4] in the permutohedron::

          sage: [p for p in TransitiveIdealGraded(attrcall("permutohedron_succ"), [Permutation([3,1,2,4]), Permutation([2,1,3,4])])]
          [[3, 1, 2, 4], [2, 1, 3, 4], [2, 1, 4, 3], [3, 2, 1, 4], [2, 3, 1, 4], [3, 1, 4, 2], [2, 3, 4, 1], [3, 4, 1, 2], [3, 2, 4, 1], [2, 4, 1, 3], [2, 4, 3, 1], [4, 3, 1, 2], [4, 2, 1, 3], [3, 4, 2, 1], [4, 2, 3, 1], [4, 3, 2, 1]]

    """
    def __init__(self, succ, generators):
        """
        TESTS::

            sage: C = TransitiveIdealGraded(factor, (1, 2, 3))
            sage: C._succ
            <function factor at ...>
            sage: C._generators
            (1, 2, 3)
            sage: loads(dumps(C))   # should test for equality with C, but equality is not implemented
        """
        self._succ = succ
        self._generators = generators

    def __iter__(self):
        """
        Returns an iterator on the elements of self.

        TESTS::

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
        current_level = self._generators
        known = set(current_level)
        while len(current_level) > 0:
            next_level = set()
            for x in current_level:
                yield x
                for y in self._succ(x):
                    if y == None or y in known:
                        continue
                    next_level.add(y)
                    known.add(y)
            current_level = next_level
        return
