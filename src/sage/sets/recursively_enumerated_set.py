r"""
Recursively enumerated set

Recursively enumerated set described by some seeds and a successor
function. The successor function may have some structure (symmetric,
graded, forest) or not. Many kinds of iterators are provided: depth first
search, breadth first search or elements of given depth.

See http://en.wikipedia.org/wiki/Recursively_enumerable_set

AUTHORS:

- Sebastien Labbe, April 2014, at Sage Days 57, Cernay-la-ville

EXAMPLES:

We construct a recursively enumerated set with symmetric structure and
depth first search for default enumeration algorithm::

    sage: succ = lambda a: [(a[0]-1,a[1]), (a[0],a[1]-1), (a[0]+1,a[1]), (a[0],a[1]+1)]
    sage: seeds = [(0,0)]
    sage: C = RecursivelyEnumeratedSet(seeds, succ, structure='symmetric', algorithm='depth')
    sage: C
    A recursively enumerated set with a symmetric structure (depth first search)

In this case, depth first search is the default algorithm for iteration::

    sage: it_depth = iter(C)
    sage: [next(it_depth) for _ in range(10)]
    [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9)]

Breadth first search::

    sage: it_breadth = C.breadth_first_search_iterator()
    sage: [next(it_breadth) for _ in range(10)]
    [(0, 0), (0, 1), (0, -1), (1, 0), (-1, 0), (-1, 1), (-2, 0), (0, 2), (2, 0), (-1, -1)]

Level (elements of given depth) iterator::

    sage: it_level = C.level_iterator()
    sage: for _ in range(3): sorted(next(it_level))
    [(0, 0)]
    [(-1, 0), (0, -1), (0, 1), (1, 0)]
    [(-2, 0), (-1, -1), (-1, 1), (0, -2), (0, 2), (1, -1), (1, 1), (2, 0)]

"""
#*****************************************************************************
#       Copyright (C) 2014 Sebastien Labbe <slabqc at gmail.com>
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
from sage.structure.parent import Parent
from sage.categories.enumerated_sets import EnumeratedSets
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
class RecursivelyEnumeratedSet(Parent):
    r"""
    Return a recursively enumerated set.

    A set S is called recursively enumerable if there is an algorithm that
    enumerates the members of S. We consider here the recursively
    enumerated set that are described by some ``seeds`` and a successor
    function ``succ``. 

    Let `U` be a set and ``succ`` `:U\to 2^U` be a successor function
    associating to each element of `U` a subset of `U`. Let ``seeds`` be a
    subset of `U`. Let `S\subseteq U` be the set of elements of `U` that
    can be reached from a seed by applying recursively the ``succ`` function.
    This class provides different kinds of iterators (breadth first, depth
    first, elements of given depth, etc.) for the elements of `S`.

    See http://en.wikipedia.org/wiki/Recursively_enumerable_set

    INPUT:

    - ``seeds`` -- list (or iterable) of hashable objects
    - ``succ`` -- function (or callable) returning a list (or iterable)
    - ``structure`` -- string (optional, default: ``None``), structure of the
      set, possible values are:

      - ``None`` -- nothing is known about the structure of the set.
      - ``'forest'`` -- if the ``succ`` function generates a *forest*, that
        is, each element can be reached uniquely from a seed.
      - ``'graded'`` -- if the ``succ`` function is *graded*, that is, all
        paths from a seed to a given element have equal length.
      - ``'symmetric'`` -- if the relation is *symmetric*, that is, 
        ``y in succ(x)`` if and only if ``x in succ(y)``

    - ``algorithm`` -- ``'depth'``, ``'breadth'``, ``'naive'`` or ``None``
      (optional, default: ``None``)
    - ``max_depth`` -- integer (optional, default: ``float("inf")``), only
      for breadth first search
    - ``post_process`` -- (optional, default: ``None``), for forest only
    - ``facade`` -- (optional, default: ``None``)
    - ``category`` -- (optional, default: ``None``)

    EXAMPLES:

    A recursive set with no other information::

        sage: f = lambda a: [a+3, a+5]
        sage: C = RecursivelyEnumeratedSet([0], f)
        sage: C
        A recursively enumerated set (naive search)
        sage: it = iter(C)
        sage: [next(it) for _ in range(10)]
        [0, 3, 8, 11, 5, 6, 9, 10, 12, 13]

    A recursive set with a forest structure::

        sage: f = lambda a: [2*a,2*a+1]
        sage: C = RecursivelyEnumeratedSet([1], f, structure='forest')
        sage: C
        An enumerated set with a forest structure
        sage: it = C.depth_first_search_iterator()
        sage: [next(it) for _ in range(7)]
        [1, 2, 4, 8, 16, 32, 64]
        sage: it = C.breadth_first_search_iterator()
        sage: [next(it) for _ in range(7)]
        [1, 2, 3, 4, 5, 6, 7]

    A recursive set given by a symmetric relation::

        sage: f = lambda a: [a-1,a+1]
        sage: C = RecursivelyEnumeratedSet([10, 15], f, structure='symmetric')
        sage: C
        A recursively enumerated set with a symmetric structure (breadth first search)
        sage: it = iter(C)
        sage: [next(it) for _ in range(7)]
        [10, 15, 16, 9, 11, 14, 8]

    A recursive set given by a graded relation::

        sage: f = lambda a: [a+1, a+I]
        sage: C = RecursivelyEnumeratedSet([0], f, structure='graded')
        sage: C
        A recursively enumerated set with a graded structure (breadth first search)
        sage: it = iter(C)
        sage: [next(it) for _ in range(7)]
        [0, 1, I, I + 1, 2, 2*I, I + 2]

    .. WARNING:: 

        If you do not set the good structure, you might obtain bad results,
        like elements generated twice::

            sage: f = lambda a: [a-1,a+1]
            sage: C = RecursivelyEnumeratedSet([0], f, structure='graded')
            sage: it = iter(C)
            sage: [next(it) for _ in range(7)]
            [0, 1, -1, 0, 2, -2, 1]

    TESTS::

        sage: C = RecursivelyEnumeratedSet((1, 2, 3), factor)
        sage: C._succ
        <function factor at ...>
        sage: C._seeds
        (1, 2, 3)
        sage: loads(dumps(C))
        A recursively enumerated set (naive search)
    """
    __metaclass__ = ClasscallMetaclass
    @staticmethod
    def __classcall_private__(cls, seeds, succ, structure=None,
            algorithm=None, max_depth=float("inf"), post_process=None,
            facade=None, category=None):
        r"""
        EXAMPLES::

            sage: f = lambda a:[a+1]

        Different structure for the sets::

            sage: RecursivelyEnumeratedSet([0], f, structure=None)
            A recursively enumerated set (naive search)
            sage: RecursivelyEnumeratedSet([0], f, structure='graded')
            A recursively enumerated set with a graded structure (breadth first search)
            sage: RecursivelyEnumeratedSet([0], f, structure='symmetric')
            A recursively enumerated set with a symmetric structure (breadth first search)
            sage: RecursivelyEnumeratedSet([0], f, structure='forest')
            An enumerated set with a forest structure

        Different default enumeration algorithms::

            sage: RecursivelyEnumeratedSet([0], f, algorithm='breadth')
            A recursively enumerated set (breadth first search)
            sage: RecursivelyEnumeratedSet([0], f, algorithm='naive')
            A recursively enumerated set (naive search)
            sage: RecursivelyEnumeratedSet([0], f, algorithm='depth')
            A recursively enumerated set (depth first search)

        """
        if structure is None:
            if algorithm is None: algorithm = 'naive'
            return typecall(RecursivelyEnumeratedSet, seeds, succ,
                    algorithm, max_depth, facade=facade, category=category)
        elif structure == 'symmetric':
            if algorithm is None: algorithm = 'breadth'
            return RecursivelyEnumeratedSet_symmetric(seeds, succ,
                    algorithm, max_depth, facade=facade, category=category)
        elif structure == 'forest':
            if algorithm is None: algorithm = 'depth'
            from sage.combinat.backtrack import SearchForest
            return SearchForest(roots=seeds, children=succ, algorithm=algorithm,
                 post_process=post_process, facade=facade, category=category)
        elif structure == 'graded':
            if algorithm is None: algorithm = 'breadth'
            return RecursivelyEnumeratedSet_graded(seeds, succ, algorithm,
                    max_depth, facade=facade, category=category)
        else:
            raise ValueError("Unknown value for structure (=%s)" % structure)
    def __init__(self, seeds, succ, 
                 algorithm='depth', max_depth=float("inf"), 
                 facade = None, category=None):
        r"""
        TESTS::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: C
            A recursively enumerated set (naive search)
        """
        self._seeds = seeds
        self._succ = succ
        assert algorithm in ['naive', 'depth', 'breadth'], "unknown algorithm(=%s)" % self._algorithm
        self._algorithm = algorithm
        self._max_depth = max_depth
        Parent.__init__(self, facade=facade, category=EnumeratedSets().or_subcategory(category))

    # Disable __len__ from Parent (#12955)
    # Because Python assumes __len__ is fast and we can't have a fast
    # default implmentation
    __len__ = None

    def __iter__(self):
        r"""
        Return an iterator on the elements of ``self``. 
        
        The enumeration is done depth first or breadth first depending on
        the value of ``self._algorithm``.

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: it_naive = iter(RecursivelyEnumeratedSet([0], f, algorithm='naive'))
            sage: it_depth = iter(RecursivelyEnumeratedSet([0], f, algorithm='depth'))
            sage: it_breadth = iter(RecursivelyEnumeratedSet([0], f, algorithm='breadth'))
            sage: [next(it_naive) for _ in range(10)]
            [0, 3, 8, 11, 5, 6, 9, 10, 12, 13]
            sage: [next(it_depth) for _ in range(10)]
            [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]
            sage: [next(it_breadth) for _ in range(10)]
            [0, 3, 5, 8, 10, 6, 9, 11, 13, 15]

        """
        if self._algorithm == 'naive':
            return self.naive_search_iterator()
        elif self._algorithm == 'breadth':
            return self.breadth_first_search_iterator(max_depth=self._max_depth)
        elif self._algorithm == 'depth':
            return self.depth_first_search_iterator()
        else:
            raise ValueError("unknown value for algorithm(=%s)" % self._algorithm)
    def _repr_(self):
        r"""
        TESTS::

            sage: RecursivelyEnumeratedSet([1], lambda x: [x+1, x+I], structure=None)
            A recursively enumerated set (naive search)

        ::

            sage: RecursivelyEnumeratedSet([1], lambda x: [x+1, x+I], structure='graded')
            A recursively enumerated set with a graded structure (breadth first search)

        ::

            sage: RecursivelyEnumeratedSet([1], lambda x: [x-1, x+1], structure='symmetric')
            A recursively enumerated set with a symmetric structure (breadth first search)
        """
        L = ["A recursively enumerated set"]
        classname = self.__class__.__name__
        if classname.startswith('RecursivelyEnumeratedSet_graded'):
            L.append("with a graded structure")
        elif classname.startswith('RecursivelyEnumeratedSet_symmetric'):
            L.append("with a symmetric structure")
        elif classname.startswith('RecursivelyEnumeratedSet'):
            pass
        else:
            pass
        if self._algorithm in ['depth', 'breadth']:
            L.append("(%s first search)" % self._algorithm)
        else:
            L.append("(%s search)" % self._algorithm)
        return " ".join(L)

    def _breadth_first_search_iterator_from_level_iterator(self, max_depth=None):
        r"""
        Returns an iterator on the elements of self (breadth first).

        INPUT:

        - ``max_depth`` -- (Default: ``None``) Specifies the maximal depth
          to which elements are computed. If None, the value of
          ``self._max_depth`` is used.

        .. NOTE:: 
        
            It should be slower than the other one since it must generates
            the whole level before yielding the first element of each
            level. It is used for test only.

        EXAMPLES::

            sage: f = lambda a: [(a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='graded')
            sage: it = C._breadth_first_search_iterator_from_level_iterator(max_depth=3)
            sage: list(it)
            [(0, 0), (0, 1), (1, 0), (2, 0), (1, 1), (0, 2)]

        This iterator is used by default for symmetric structure::

            sage: f = lambda a: [a-1,a+1]
            sage: S = RecursivelyEnumeratedSet([10], f, structure='symmetric')
            sage: it = iter(S)
            sage: [next(it) for _ in range(7)]
            [10, 9, 11, 8, 12, 13, 7]
        """
        if max_depth is None:
            max_depth = self._max_depth
        it = self.level_iterator()
        i = 0
        while i < max_depth:
            level = next(it)
            for a in level:
                yield a
            i += 1

    def level_iterator(self):
        r"""
        Returns an iterator over the levels of self.

        A level is a set of elements of the same depth.

        It is currently implemented only for herited classes.

        OUTPUT:

            an iterator of sets

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: it = C.level_iterator()    # todo: not implemented
        """
        raise NotImplementedError

    def breadth_first_search_iterator(self, max_depth=None):
        r"""
        Returns an iterator on the elements of self (breadth first).

        This code remembers every elements generated.

        INPUT:

        - ``max_depth`` -- (Default: ``None``) Specifies the
            maximal depth to which elements are computed.
            If None, the value of ``self._max_depth`` is used.

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: it = C.breadth_first_search_iterator()
            sage: [next(it) for _ in range(10)]
            [0, 3, 5, 8, 10, 6, 9, 11, 13, 15]
        """
        if max_depth is None:
            max_depth = self._max_depth
        current_level = self._seeds
        known = set(current_level)
        depth = 0
        while len(current_level) > 0 and depth <= max_depth:
            next_level = set()
            for x in current_level:
                yield x
                for y in self._succ(x):
                    if y == None or y in known:
                        continue
                    next_level.add(y)
                    known.add(y)
            current_level = next_level
            depth += 1

    def _breadth_first_search_iterator_using_queue(self):
        r"""
        Returns an iterator on the elements of self (breadth first).

        This code remembers every elements generated and uses python
        queues. It is 3 times slower than the other one.

        See http://en.wikipedia.org/wiki/Breadth-first_search

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: it = C._breadth_first_search_iterator_using_queue()
            sage: [next(it) for _ in range(10)]
            [0, 3, 5, 6, 8, 10, 9, 11, 13, 15]

        """
        known = set(self._seeds)
        from Queue import Queue
        q = Queue()
        for s in self._seeds:
            q.put(s)
        while q.not_empty:
            x = q.get()
            yield x
            for y in self._succ(x):
                if y is None or y in known:
                    continue
                q.put(y)
                known.add(y)

    def naive_search_iterator(self):
        r"""
        Returns an iterator on the elements of self (in no particular order).

        This code remembers every elements generated.

        TESTS:

        We compute all the permutations of 3::

            sage: seeds = [Permutation([1,2,3])]
            sage: succ = attrcall("permutohedron_succ")
            sage: R = RecursivelyEnumeratedSet(seeds, succ)
            sage: list(R.naive_search_iterator())
            [[1, 2, 3], [2, 1, 3], [1, 3, 2], [2, 3, 1], [3, 2, 1], [3, 1, 2]]

        """
        known = set(self._seeds)
        todo = known.copy()
        while len(todo) > 0:
            x = todo.pop()
            yield x
            for y in self._succ(x):
                if y == None or y in known:
                    continue
                todo.add(y)
                known.add(y)

    def depth_first_search_iterator(self):
        r"""
        Returns an iterator on the elements of self (depth first).

        This code remembers every elements generated.

        See http://en.wikipedia.org/wiki/Depth-first_search

        EXAMPLES::

            sage: f = lambda a: [a+3, a+5]
            sage: C = RecursivelyEnumeratedSet([0], f)
            sage: it = C.depth_first_search_iterator()
            sage: [next(it) for _ in range(10)]
            [0, 5, 10, 15, 20, 25, 30, 35, 40, 45]

        """
        stack = []    
        for s in self._seeds:
            stack.append(s)
        known = set()
        while stack:
            x = stack.pop()
            if x is None or x in known:
                continue
            yield x
            known.add(x)
            for y in self._succ(x):
                stack.append(y)

class RecursivelyEnumeratedSet_symmetric(RecursivelyEnumeratedSet):
    r"""
    Generic tool for constructing ideals of a symmetric relation.

    INPUT:

    - ``seeds`` -- list (or iterable) of hashable objects
    - ``succ`` -- function (or callable) returning a list (or iterable)
    - ``algorithm`` -- ``'depth'``, ``'breadth'`` or ``None`` (default: ``None``)
    - ``max_depth`` -- integer (default: ``float("inf")``)

    EXAMPLES::

        sage: f = lambda a: [a-1,a+1]
        sage: C = RecursivelyEnumeratedSet([0], f, structure='symmetric')
        sage: C
        A recursively enumerated set with a symmetric structure (breadth first search)
        sage: it = iter(C)
        sage: [next(it) for _ in range(7)]
        [0, 1, -1, 2, -2, 3, -3]

    TESTS:

    Do not use lambda functions for saving purposes::

        sage: f = lambda a: [a-1,a+1]
        sage: C = RecursivelyEnumeratedSet([0], f, structure='symmetric')
        sage: loads(dumps(C))
        Traceback (most recent call last):
        ...
        PicklingError: Can't pickle <type 'function'>: attribute lookup __builtin__.function failed

    This works in the command line but apparently not as a doctest::

        sage: def f(a): return [a-1,a+1]
        sage: C = RecursivelyEnumeratedSet([0], f, structure='symmetric')
        sage: loads(dumps(C))
        Traceback (most recent call last):
        ...
        PicklingError: Can't pickle <type 'function'>: attribute lookup __builtin__.function failed

    """
    breadth_first_search_iterator = RecursivelyEnumeratedSet._breadth_first_search_iterator_from_level_iterator
    def elements_of_depth_iterator(self, depth=0):
        r"""
        Returns an iterator over the elements of ``self`` of given depth.
        An element of depth `n` can be obtained applying `n` times the
        children function from a root.

        EXAMPLES::

            sage: f = lambda a: [a-1, a+1]
            sage: S = RecursivelyEnumeratedSet([5, 10], f, structure='symmetric')
            sage: it = S.elements_of_depth_iterator(2)    #todo: not implemented
            sage: sorted(it)                              #todo: not implemented
            [3, 7, 8, 12]
        """
        raise NotImplementedError

    def level_iterator(self):
        r"""
        Returns an iterator over the levels of self.

        A level is a set of elements of the same depth.

        The algorithm remembers only the last two level generated since the
        structure is symmetric.

        OUTPUT:

            an iterator of sets

        EXAMPLES::

            sage: f = lambda a: [a-1, a+1]
            sage: S = RecursivelyEnumeratedSet([10], f, structure='symmetric')
            sage: it = S.level_iterator()
            sage: [sorted(next(it)) for _ in range(5)]
            [[10], [9, 11], [8, 12], [7, 13], [6, 14]]

        Starting with two generators::

            sage: f = lambda a: [a-1, a+1]
            sage: S = RecursivelyEnumeratedSet([5, 10], f, structure='symmetric')
            sage: it = S.level_iterator()
            sage: [sorted(next(it)) for _ in range(5)]
            [[5, 10], [4, 6, 9, 11], [3, 7, 8, 12], [2, 13], [1, 14]]

        Gaussian integers::

            sage: f = lambda a: [a+1, a+I]
            sage: S = RecursivelyEnumeratedSet([0], f, structure='symmetric')
            sage: it = S.level_iterator()
            sage: [sorted(next(it)) for _ in range(7)]
            [[0],
             [I, 1],
             [2*I, I + 1, 2],
             [3*I, 2*I + 1, I + 2, 3],
             [4*I, 3*I + 1, 2*I + 2, I + 3, 4],
             [5*I, 4*I + 1, 3*I + 2, 2*I + 3, I + 4, 5],
             [6*I, 5*I + 1, 4*I + 2, 3*I + 3, 2*I + 4, I + 5, 6]]
        """
        A = set()
        B = set(self._seeds)
        while len(B) > 0:
            yield B
            A,B = B, self._get_next_level(A, B)

    def _get_next_level(self, A, B):
        r"""
        Return the set of elements of depth n+1.

        INPUT:

        - ``A`` -- set, the set of elements of depth n-1
        - ``B`` -- set, the set of elements of depth n

        OUTPUT:

        - ``C`` -- set, the set of elements of depth n+1

        EXAMPLES::

            sage: f = lambda a: [a-1, a+1]
            sage: S = RecursivelyEnumeratedSet([5, 10], f, structure='symmetric')
            sage: sorted(S._get_next_level([2,8], [3,7]))
            [4, 6]
            sage: sorted(S._get_next_level([3,7], [2,8]))
            [1, 9]
        """
        C = set()
        for x in B:
            for y in self._succ(x):
                if (y is None or y in A or y in B):
                    continue
                C.add(y)
        return C

class RecursivelyEnumeratedSet_graded(RecursivelyEnumeratedSet):
    r"""
    Generic tool for constructing ideals of a graded relation.

    INPUT:

    - ``seeds`` -- list (or iterable) of hashable objects
    - ``succ`` -- function (or callable) returning a list (or iterable)
    - ``algorithm`` -- ``'depth'``, ``'breadth'`` or ``None`` (default: ``None``)
    - ``max_depth`` -- integer (default: ``float("inf")``)

    EXAMPLES::

        sage: f = lambda a: [(a[0]+1,a[1]), (a[0],a[1]+1)]
        sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='graded', max_depth=3)
        sage: C
        A recursively enumerated set with a graded structure (breadth first search)
        sage: sorted(C)
        [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (3, 0)]

    """
    def breadth_first_search_iterator(self, max_depth=None):
        r"""
        Return an iterator on the elements of self (breadth first).

        This iterator make use of the graded structure by remembering only
        the elements of the current depth.

        INPUT:

        - ``max_depth`` -- (Default: ``None``) Specifies the
          maximal depth to which elements are computed.
          If None, the value of ``self._max_depth`` is used.

        EXAMPLES::

            sage: f = lambda a: [(a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='graded')
            sage: it = C.breadth_first_search_iterator(max_depth=3)
            sage: list(it)
            [(0, 0), (0, 1), (1, 0), (2, 0), (1, 1), (0, 2), (3, 0), (1, 2), (0, 3), (2, 1)]
        """
        if max_depth is None:
            max_depth = self._max_depth
        current_level = self._seeds
        depth = 0
        while len(current_level) > 0 and depth <= max_depth:
            next_level = set()
            for x in current_level:
                yield x
                for y in self._succ(x):
                    if y == None or y in next_level:
                        continue
                    next_level.add(y)
            current_level = next_level
            depth += 1
    def level_iterator(self):
        r"""
        Returns an iterator over the levels of self.

        A level is a set of elements of the same depth.

        The algorithm remembers only the current level generated since the
        structure is graded.

        OUTPUT:

            an iterator of sets

        EXAMPLES::

            sage: f = lambda a: [(a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='graded', max_depth=3)
            sage: it = C.level_iterator()
            sage: for _ in range(4): sorted(next(it))
            [(0, 0)]
            [(0, 1), (1, 0)]
            [(0, 2), (1, 1), (2, 0)]
            [(0, 3), (1, 2), (2, 1), (3, 0)]
        """
        B = self._seeds
        while len(B) > 0:
            yield B
            B = self._get_next_level(B)

    def _get_next_level(self, B):
        r"""
        Return the set of elements of depth n+1.

        INPUT:

        - ``B`` -- set, the set of elements of depth n

        OUTPUT:

        - ``C`` -- set, the set of elements of depth n+1

        EXAMPLES::

            sage: f = lambda a: [(a[0]+1,a[1]), (a[0],a[1]+1)]
            sage: C = RecursivelyEnumeratedSet([(0,0)], f, structure='graded')
            sage: sorted(C._get_next_level(C._seeds))
            [(0, 1), (1, 0)]
            sage: sorted(C._get_next_level(_))
            [(0, 2), (1, 1), (2, 0)]
            sage: sorted(C._get_next_level(_))
            [(0, 3), (1, 2), (2, 1), (3, 0)]
        """
        C = set()
        for x in B:
            for y in self._succ(x):
                if (y is None or y in B):
                    continue
                C.add(y)
        return C

