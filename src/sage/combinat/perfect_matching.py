r"""
Perfect matchings

A perfect matching of a set `S` is a partition into 2-element sets. If `S` is
the set `\{1,...,n\}`, it is equivalent to fixpoint-free involutions. These
simple combinatorial objects appear in different domains such as combinatoric
of orthogonal polynomials and of the hyperoctaedral groups (see [MV]_, [McD]_
and also [CM]_):

AUTHOR:

    - Valentin Feray, 2010 : initial version
    - Martin Rubey, 2017: inherit from SetPartition

EXAMPLES:

    Create a perfect matching::

        sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')]);m
        [('a', 'e'), ('c', 'b'), ('d', 'f')]

    Count its crossings, if the ground set is totally ordered::

        sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
        [(1, 3), (8, 2), (4, 7), (5, 6)]
        sage: n.number_of_crossings()
        1

    List the perfect matchings of a given ground set::

        sage: PerfectMatchings(4).list()
        [[(1, 2), (3, 4)], [(1, 3), (2, 4)], [(1, 4), (2, 3)]]

REFERENCES:

    .. [MV] combinatorics of orthogonal polynomials (A. de Medicis et
       X.Viennot, Moments des q-polynomes de Laguerre et la bijection de
       Foata-Zeilberger, Adv. Appl. Math., 15 (1994), 262-304)

    .. [McD] combinatorics of hyperoctahedral group, double coset algebra and
       zonal polynomials (I. G. Macdonald, Symmetric functions and Hall
       polynomials, Oxford University Press, second edition, 1995, chapter
       VII).

    .. [CM] Benoit Collins, Sho Matsumoto, On some properties of
       orthogonal Weingarten functions, :arxiv:`0903.5143`.
"""
#*****************************************************************************
#       Copyright (C) 2010 Valentin Feray <feray@labri.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
# python3
from __future__ import division, print_function
from six.moves import range
from six import add_metaclass

from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.combinat.permutation import Permutation
from sage.sets.set import Set, is_Set
from sage.misc.misc_c import prod
from sage.matrix.constructor import matrix
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.set_partition import SetPartition, SetPartitions
from sage.combinat.partition import Partition
from sage.rings.infinity import infinity


class PerfectMatching(SetPartition):
    r"""
    Class of perfect matching.

    An instance of the class can be created from a list of pairs or from a
    fixed point-free involution as follows::

        sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')]);m
        [('a', 'e'), ('c', 'b'), ('d', 'f')]
        sage: n = PerfectMatching([3,8,1,7,6,5,4,2]);n
        [(1, 3), (8, 2), (4, 7), (5, 6)]
        sage: isinstance(m,PerfectMatching)
        True

    The parent, which is the set of perfect matchings of the ground set, is
    automatically created::

        sage: n.parent()
        Set of perfect matchings of {1, 2, 3, 4, 5, 6, 7, 8}

    If the ground set is ordered, one can, for example, ask if the matching is
    non crossing::

        sage: PerfectMatching([(1, 4), (2, 3), (5, 6)]).is_non_crossing()
        True

    TESTS::

        sage: m = PerfectMatching([]); m
        []
        sage: m.parent()
        Set of perfect matchings of {}
    """
    @staticmethod
    def __classcall_private__(cls, parts, check=True):
        """
        Create a perfect matching from ``parts`` with the appropriate parent.

        This function tries to recognize the input (it can be either a list or
        a tuple of pairs, or a fix-point free involution given as a list or as
        a permutation), constructs the parent (enumerated set of
        PerfectMatchings of the ground set) and calls the __init__ function to
        construct our object.

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')]);m
            [('a', 'e'), ('c', 'b'), ('d', 'f')]
            sage: isinstance(m, PerfectMatching)
            True
            sage: n = PerfectMatching([3, 8, 1, 7, 6, 5, 4, 2]);n
            [(1, 3), (8, 2), (4, 7), (5, 6)]
            sage: n.parent()
            Set of perfect matchings of {1, 2, 3, 4, 5, 6, 7, 8}
            sage: PerfectMatching([(1, 4), (2, 3), (5, 6)]).is_non_crossing()
            True

        The function checks that the given list or permutation is a valid perfect
        matching (i.e. a list of pairs with pairwise disjoint elements  or a
        fixpoint-free involution) and raises a ValueError otherwise::

            sage: PerfectMatching([(1, 2, 3), (4, 5)])
            Traceback (most recent call last):
            ...
            ValueError: [(1, 2, 3), (4, 5)] is not a valid perfect matching: all elements of the list must be pairs

        If you know your datas are in a good format, use directly
        ``PerfectMatchings(objects)(data)``.

        TESTS::

             sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')])
             sage: TestSuite(m).run()
             sage: m = PerfectMatching([])
             sage: TestSuite(m).run()
             sage: PerfectMatching(6)
             Traceback (most recent call last):
             ...
             ValueError: cannot convert p (= 6) to a PerfectMatching
             sage: PerfectMatching([(1,2,3)])
             Traceback (most recent call last):
             ...
             ValueError: [(1, 2, 3)] is not a valid perfect matching:
             all elements of the list must be pairs
             sage: PerfectMatching([(1,1)])
             Traceback (most recent call last):
             ...
             ValueError: [(1, 1)] is not a valid perfect matching:
             there are some repetitions
             sage: PerfectMatching(Permutation([4,2,1,3]))
             Traceback (most recent call last):
             ...
             ValueError: The permutation p (= [4, 2, 1, 3]) is not a fixed point free involution
        """
        if ((isinstance(parts, list) and
             all((isinstance(x, (int, Integer)) for x in parts)))
            or isinstance(parts, Permutation)):
            s = Permutation(parts)
            if not all(e == 2 for e in s.cycle_type()):
                raise ValueError("The permutation p (= %s) is not a "
                                 "fixed point free involution" % s)
            parts = s.to_cycles()

        base_set = reduce(lambda x,y: x.union(y), map(Set, parts), Set([]))
        P = PerfectMatchings(base_set)
        return P.element_class(P, parts, check=check)

    def _repr_(self):
        r"""
        Return a string representation of the matching ``self``.

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')]);m
            [('a', 'e'), ('c', 'b'), ('d', 'f')]
            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]);n
            [(1, 3), (8, 2), (4, 7), (5, 6)]
        """
        return str([(a,b) for a,b in self])

    def _latex_(self):
        r"""
        A latex representation of ``self`` using the tikzpicture package.

        EXAMPLES::

            sage: P = PerfectMatching([(1,3),(2,5),(4,6)])
            sage: latex(P)  # optional - dot2tex; random
            \begin{tikzpicture}
            ...
            \end{tikzpicture}

        TESTS:

        Above we added ``random`` since warnings might be displayed
        once. The second time, there should be no warnings::

            sage: print(P._latex_())  # optional - dot2tex
            \begin{tikzpicture}
            ...
            \end{tikzpicture}

        ..TODO::

            This should probably call the latex method of
            :class:`SetPartition` with appropriate defaults.
        """
        G = self.to_graph()
        G.set_pos(G.layout_circular())
        G.set_latex_options(
            vertex_size=0.4,
            edge_thickness=0.04,
        )
        return G._latex_()

    def partner(self, x):
        r"""
        Returns the element in the same pair than ``x`` in the matching ``self``.

        EXAMPLES::

            sage: m = PerfectMatching([(-3, 1), (2, 4), (-2, 7)]); m.partner(4)
            2
            sage: n = PerfectMatching([('c','b'),('d','f'),('e','a')])
            sage: n.partner('c')
            'b'
        """
        for a, b in self:
            if a == x:
                return b
            if b == x:
                return a
        raise ValueError("%s in not an element of the %s" % (x, self))

    def conjugate_by_permutation(self, p):
        r"""
        This is deprecated. Use :meth:`apply_permutation` instead.

        Returns the conjugate of the perfect matching ``self`` by the
        permutation ``p`` of the ground set.

        EXAMPLES::

            sage: m = PerfectMatching([(1,4),(2,6),(3,5)])
            sage: m.conjugate_by_permutation(Permutation([4,1,5,6,3,2]))
            doctest:...: DeprecationWarning: conjugate_by_permutation is deprecated; use apply_permutation instead
            See http://trac.sagemath.org/23982 for details.

            [(1, 2), (3, 5), (4, 6)]

        TESTS::

            sage: PerfectMatching([]).conjugate_by_permutation(Permutation([]))
            []

        """
        from sage.misc.superseded import deprecation
        deprecation(23982, "conjugate_by_permutation is deprecated; use apply_permutation instead")

        return self.apply_permutation(p)

    def loops_iterator(self, other=None):
        r"""
        INPUT:

             - ``other`` -- a perfect matching of the same set of ``self``.
               (if the second argument is empty, the method :meth:`an_element` is
               called on the parent of the first)

        OUTPUT:

            If we draw the two perfect matchings simultaneously as edges of a
            graph, the graph obtained is a union of cycles of even lengths.
            The function returns an iterator for these cycles (each cycle is
            given as a list).

        EXAMPLES::

            sage: o = PerfectMatching([(1, 7), (2, 4), (3, 8), (5, 6)])
            sage: p = PerfectMatching([(1, 6), (2, 7), (3, 4), (5, 8)])
            sage: it = o.loops_iterator(p)
            sage: next(it)
            [1, 7, 2, 4, 3, 8, 5, 6]
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration
        """
        if other is None:
            other = self.parent().an_element()
        elif self.parent() != other.parent():
            s = "%s is not a matching of the ground set of %s" % (other, self)
            raise ValueError(s)
        remain = self.base_set().set()
        while len(remain) > 0:
            a = remain.pop()
            b = self.partner(a)
            remain.remove(b)
            loop = [a, b]
            c = other.partner(b)
            while c != a:
                b = self.partner(c)
                remain.remove(c)
                loop.append(c)
                remain.remove(b)
                loop.append(b)
                c = other.partner(b)
            yield loop

    def loops(self, other=None):
        r"""
        INPUT:

             - ``other`` -- a perfect matching of the same set of ``self``.
               (if the second argument is empty, the method :meth:`an_element` is
               called on the parent of the first)

        OUTPUT:

            If we draw the two perfect matchings simultaneously as edges of a
            graph, the graph obtained is a union of cycles of even lengths.
            The function returns the list of these cycles (each cycle is given
            as a list).

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')])
            sage: n = PerfectMatching([('a','b'),('d','f'),('e','c')])
            sage: m.loops(n)
            [['a', 'e', 'c', 'b'], ['d', 'f']]

            sage: o = PerfectMatching([(1, 7), (2, 4), (3, 8), (5, 6)])
            sage: p = PerfectMatching([(1, 6), (2, 7), (3, 4), (5, 8)])
            sage: o.loops(p)
            [[1, 7, 2, 4, 3, 8, 5, 6]]
        """
        return list(self.loops_iterator(other))

    def loop_type(self, other=None):
        r"""
        INPUT:

             - ``other`` -- a perfect matching of the same set of ``self``.
               (if the second argument is empty, the method :meth:`an_element` is
               called on the parent of the first)

        OUTPUT:

            If we draw the two perfect matchings simultaneously as edges of a
            graph, the graph obtained is a union of cycles of even
            lengths. The function returns the ordered list of the semi-length
            of these cycles (considered as a partition)

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')])
            sage: n = PerfectMatching([('a','b'),('d','f'),('e','c')])
            sage: m.loop_type(n)
            [2, 1]

        TESTS::

            sage: m = PerfectMatching([]); m.loop_type()
            []
        """
        return Partition(reversed(
                sorted([len(l)//2 for l in self.loops_iterator(other)])))

    def number_of_loops(self, other=None):
        r"""
        INPUT:

            - ``other`` -- a perfect matching of the same set of ``self``.
              (if the second argument is empty, the method :meth:`an_element` is
              called on the parent of the first)

        OUTPUT:

            If we draw the two perfect matchings simultaneously as edges of a
            graph, the graph obtained is a union of cycles of even lengths.
            The function returns their numbers.

        EXAMPLES::

            sage: m = PerfectMatching([('a','e'),('b','c'),('d','f')])
            sage: n = PerfectMatching([('a','b'),('d','f'),('e','c')])
            sage: m.number_of_loops(n)
            2
        """
        c = Integer(0)
        one = Integer(1)
        for _ in self.loops_iterator(other):
            c += one
        return c

    def crossings_iterator(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns an iterator over the pairs of
            crossing lines (as a line correspond to a pair, the iterator
            produces pairs of pairs).

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (8, 2), (4, 7), (5, 6)]
            sage: it = n.crossings_iterator();
            sage: next(it)
            ((1, 3), (8, 2))
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration
        """
        x = list(self)
        while x:
            (i, j) = x.pop(0)
            for (a, b) in x:
                # if (i<a<j<b) or (i<b<j<a) or (j<a<i<b) or (j<b<i<a) or (
                #        a<i<b<j) or (a<j<b<i) or (b<i<a<j) or (b<j<a<i):
                labij = sorted([a, b, i, j])
                posij = sorted([labij.index(i), labij.index(j)])
                if posij == [0, 2] or posij == [1, 3]:
                    yield ((i, j), (a, b))

    def crossings(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns the list of the pairs of
            crossing lines (as a line correspond to a pair, it returns a list
            of pairs of pairs).

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (8, 2), (4, 7), (5, 6)]
            sage: n.crossings()
            [((1, 3), (8, 2))]

        TESTS::

            sage: m = PerfectMatching([]); m.crossings()
            []
        """
        return list(self.crossings_iterator())

    def number_of_crossings(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns the number the pairs of crossing
            lines.

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (8, 2), (4, 7), (5, 6)]
            sage: n.number_of_crossings()
            1
        """
        c = Integer(0)
        one = Integer(1)
        for _ in self.crossings_iterator():
            c += one
        return c

    def is_non_crossing(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns ``True`` if the picture obtained
            this way has no crossings.

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (8, 2), (4, 7), (5, 6)]
            sage: n.is_non_crossing()
            False
            sage: PerfectMatching([(1, 4), (2, 3), (5, 6)]).is_non_crossing()
            True
        """
        it = self.crossings_iterator()
        try:
            next(it)
        except StopIteration:
            return True
        else:
            return False

    def nestings_iterator(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns an iterator over the pairs of
            nesting lines (as a line correspond to a pair, the iterator
            produces pairs of pairs).

        EXAMPLES::

            sage: n = PerfectMatching([(1, 6), (2, 7), (3, 5), (4, 8)])
            sage: it = n.nestings_iterator();
            sage: next(it)
            ((1, 6), (3, 5))
            sage: next(it)
            ((2, 7), (3, 5))
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration
        """
        x = list(self)
        while x:
            (i, j) = x.pop(0)
            for (a, b) in x:
                # if (i<a<j<b) or (i<b<j<a) or (j<a<i<b) or (j<b<i<a) or (
                #        a<i<b<j) or (a<j<b<i) or (b<i<a<j) or (b<j<a<i):
                labij = sorted([a, b, i, j])
                posij = sorted([labij.index(i), labij.index(j)])
                if posij == [0, 3] or posij == [1, 2]:
                    yield ((i, j), (a, b))

    def nestings(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns the list of the pairs of
            nesting lines (as a line correspond to a pair, it returns a list
            of pairs of pairs).

        EXAMPLES::

            sage: m = PerfectMatching([(1, 6), (2, 7), (3, 5), (4, 8)])
            sage: m.nestings()
            [((1, 6), (3, 5)), ((2, 7), (3, 5))]

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (8, 2), (4, 7), (5, 6)]
            sage: n.nestings()
            [((8, 2), (4, 7)), ((8, 2), (5, 6)), ((4, 7), (5, 6))]

        TESTS::

            sage: m = PerfectMatching([]); m.nestings()
            []
        """
        return list(self.nestings_iterator())

    def number_of_nestings(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns the number the pairs of nesting
            lines.

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (8, 2), (4, 7), (5, 6)]
            sage: n.number_of_nestings()
            3
        """
        c = Integer(0)
        one = Integer(1)
        for _ in self.nestings_iterator():
            c += one
        return c

    def is_non_nesting(self):
        r"""
        INPUT:

            A perfect matching on a *totally ordered* ground set.

        OUTPUT:

            We place the element of a ground set and draw the perfect matching
            by linking the elements of the same pair in the upper
            half-plane. This function returns ``True`` if the picture obtained
            this way has no nestings.

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (8, 2), (4, 7), (5, 6)]
            sage: n.is_non_nesting()
            False
            sage: PerfectMatching([(1, 3), (2, 5), (4, 6)]).is_non_nesting()
            True
        """
        it = self.nestings_iterator()
        try:
            next(it)
        except StopIteration:
            return True
        else:
            return False

    def Weingarten_function(self, d, other=None):
        r"""
        Returns the Weingarten function of two pairings.

        This function is the value of some integrals over the orthogonal
        groups `O_N`.  With the convention of [CM]_, the method returns
        `Wg^{O(d)}(other,self)`.

        EXAMPLES::

            sage: var('N')
            N
            sage: m = PerfectMatching([(1,3),(2,4)])
            sage: n = PerfectMatching([(1,2),(3,4)])
            sage: factor(m.Weingarten_function(N,n))
            -1/((N + 2)*(N - 1)*N)
        """
        if other is None:
            other = self.parent().an_element()
        W = self.parent().Weingarten_matrix(d)
        return W[other.rank()][self.rank()]

    def to_graph(self):
        r"""
        Returns the graph corresponding to the perfect matching.

        OUTPUT:

            The realization of ``self`` as a graph.

        EXAMPLES::

            sage: PerfectMatching([[1,3], [4,2]]).to_graph().edges(labels=False)
            [(1, 3), (2, 4)]
            sage: PerfectMatching([[1,4], [3,2]]).to_graph().edges(labels=False)
            [(1, 4), (2, 3)]
            sage: PerfectMatching([]).to_graph().edges(labels=False)
            []
        """
        from sage.graphs.graph import Graph
        G = Graph()
        for a, b in self:
            G.add_edge((a, b))
        return G

    def to_non_crossing_set_partition(self):
        r"""
        Returns the noncrossing set partition (on half as many elements)
        corresponding to the perfect matching if the perfect matching is
        noncrossing, and otherwise gives an error.

        OUTPUT:

            The realization of ``self`` as a noncrossing set partition.

        EXAMPLES::

            sage: PerfectMatching([[1,3], [4,2]]).to_non_crossing_set_partition()
            Traceback (most recent call last):
            ...
            ValueError: matching must be non-crossing
            sage: PerfectMatching([[1,4], [3,2]]).to_non_crossing_set_partition()
            {{1, 2}}
            sage: PerfectMatching([]).to_non_crossing_set_partition()
            {}
        """
        from sage.combinat.set_partition import SetPartition
        if not self.is_non_crossing():
            raise ValueError("matching must be non-crossing")
        else:
            perm = self.to_permutation()
            perm2 = Permutation([perm[2 * i] // 2
                                 for i in range(len(perm) // 2)])
        return SetPartition(perm2.cycle_tuples())

class PerfectMatchings(SetPartitions):
    r"""
    Class of perfect matchings of a ground set. At the creation, the set
    can be given as any iterable object. If the argument is an integer `n`, it
    will be transformed into `[1 .. n]`::

        sage: M = PerfectMatchings(6);M
        Set of perfect matchings of {1, 2, 3, 4, 5, 6}
        sage: PerfectMatchings([-1, -3, 1, 2])
        Set of perfect matchings of {1, 2, -3, -1}

    One can ask for the list, the cardinality or an element of a set of
    perfect matching::

        sage: PerfectMatchings(4).list()
        [[(1, 2), (3, 4)], [(1, 3), (2, 4)], [(1, 4), (2, 3)]]
        sage: PerfectMatchings(8).cardinality()
        105
        sage: M = PerfectMatchings(('a', 'e', 'b', 'f', 'c', 'd'))
        sage: M.an_element()
        [('a', 'b'), ('c', 'd'), ('e', 'f')]
        sage: all(PerfectMatchings(i).an_element() in PerfectMatchings(i)
        ....:      for i in range(2,11,2))
        True

    TESTS::

        sage: PerfectMatchings(5).list()
        []
        sage: TestSuite(PerfectMatchings(6)).run()
        sage: TestSuite(PerfectMatchings([])).run()
    """

    @staticmethod
    def __classcall_private__(cls, s):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = PerfectMatchings(4)
            sage: T = PerfectMatchings([1,2,3,4])
            sage: S is T
            True
        """
        if isinstance(s, (int, Integer)):
            s = frozenset(range(1, s+1))
        else:
            try:
                if s.cardinality() == infinity:
                    raise ValueError("The set must be finite")
            except AttributeError:
                pass
            s = frozenset(s)
        return super(PerfectMatchings, cls).__classcall__(cls, s)

    def _element_constructor_(self, s, check=True):
        """
        Construct an element of ``self`` from ``s``.

        INPUT:

        - ``s`` -- a set of sets

        EXAMPLES::

            sage: S = PerfectMatchings(4)
            sage: elt = S([[1,3],[2,4]]); elt
            [(1, 3), (2, 4)]
            sage: S = PerfectMatchings([])
            sage: S([])
            []
        """
        if isinstance(s, PerfectMatching):
            if isinstance(s.parent(), PerfectMatchings):
                return self.element_class(self, list(s), check=check)
            raise ValueError("cannot convert %s into an element of %s"%(s, self))

        return self.element_class(self, s, check=check)

    def __init__(self, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: M = PerfectMatchings(6)
            sage: TestSuite(M).run()
            sage: PerfectMatchings(0).list()
            [[]]
            sage: PerfectMatchings([]).list()
            [[]]
        """
        self._set = s
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a description of ``self``.

        TESTS::

            sage: PerfectMatchings([-1, -3, 1, 2])
            Set of perfect matchings of {1, 2, -3, -1}
        """
        return "Set of perfect matchings of %s"%(Set(self._set))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: PerfectMatchings(4).list()
            [[(1, 2), (3, 4)], [(1, 3), (2, 4)], [(1, 4), (2, 3)]]
        """
        def iter_aux(s):
            n = len(s)
            if n == 0:
                yield []
            elif n == 1:
                pass
            else:
                a = s[0]
                for i in range(1, n):
                    b = s[i]
                    for p in iter_aux(s[1:i] + s[i+1:]):
                        yield [(a, b)]+p


        for p in iter_aux(list(self._set)):
            yield self.element_class(self, p)

    def __contains__(self, x):
        """
        Test if ``x`` is an element of ``self``.


        EXAMPLES::

            sage: m = PerfectMatching([(1,2),(4,3)])
            sage: m in PerfectMatchings(4)
            True
            sage: m in PerfectMatchings((0, 1, 2, 3))
            False
            sage: all(m in PerfectMatchings(6) for m in PerfectMatchings(6))
            True

        Note that the class of ``x`` does not need to be ``PerfectMatching``:
        if the data defines a perfect matching of the good set, the function
        returns ``True``::

            sage: [(1, 4), (2, 3)] in PerfectMatchings(4)
            True
            sage: [(1, 3, 6), (2, 4), (5,)] in PerfectMatchings(6)
            False
            sage: [('a', 'b'), ('a', 'c')] in PerfectMatchings(
            ....:      ('a', 'b', 'c', 'd'))
            False

        TESTS::

            sage: SA = PerfectMatchings([1,2,3,7])
            sage: Set([Set([1,2]),Set([3,7])]) in SA
            True
            sage: Set([Set([1,2]),Set([2,3])]) in SA
            False
            sage: Set([]) in SA
            False
        """
        if not all(len(p) == 2 for p in x):
            return False

        base_set = reduce(lambda x,y: x.union(y), map(Set, x), Set([]))
        return len(base_set) == 2*len(x) and base_set == Set(self._set)

    def base_set(self):
        """
        Return the base set of ``self``.

        EXAMPLES::

            sage: PerfectMatchings(3).base_set()
            {1, 2, 3}
        """
        return Set(self._set)

    def base_set_cardinality(self):
        """
        Return the cardinality of the base set of ``self``.

        EXAMPLES::

            sage: PerfectMatchings(3).base_set_cardinality()
            3
        """
        return len(self._set)

    def cardinality(self):
        """
        Return the cardinality of the set of perfect matchings ``self``.

        This is `1*3*5*...*(2n-1)`, where `2n` is the size of the ground set.

        EXAMPLES::

            sage: PerfectMatchings(8).cardinality()
            105
            sage: PerfectMatchings([1,2,3,4]).cardinality()
            3
            sage: PerfectMatchings(3).cardinality()
            0
            sage: PerfectMatchings([]).cardinality()
            1

        """
        n = len(self._set)
        if n % 2 == 1:
            return 0
        else:
            return Integer(prod(i for i in range(n) if i % 2 == 1))

    def random_element(self):
        r"""
        Returns a random element of self.

        ..TODO::

            This really belongs to :class:`SetPartition`!

        EXAMPLES::

            sage: M = PerfectMatchings(('a', 'e', 'b', 'f', 'c', 'd'))
            sage: M.an_element()
            [('a', 'b'), ('c', 'd'), ('e', 'f')]
            sage: all(PerfectMatchings(2*i).an_element() in PerfectMatchings(2*i)
            ....:      for i in range(2,11,2))
            True

        TESTS::

            sage: p = PerfectMatchings(13).random_element()
            Traceback (most recent call last):
            ...
            ValueError: there is no perfect matching on an odd number of elements

        """
        n = len(self._set)

        if n % 2 == 1:
            raise ValueError("there is no perfect matching on an odd number of elements")

        k = n//2

        from sage.combinat.permutation import Permutations
        p = Permutations(n).random_element()
        l = list(self._set)
        return self([(l[p[2*i]-1], l[p[2*i+1]-1]) for i in range(k)])

    an_element = random_element

    @cached_method
    def Weingarten_matrix(self, N):
        r"""
        Returns the Weingarten matrix corresponding to the set of
        PerfectMatchings ``self``.

        It is a useful theoretical tool to compute polynomial integral
        over the orthogonal group `O_N` (see [CM]_).

        EXAMPLES::

            sage: M = PerfectMatchings(4).Weingarten_matrix(var('N'))
            sage: N*(N-1)*(N+2)*M.apply_map(factor)
            [N + 1    -1    -1]
            [   -1 N + 1    -1]
            [   -1    -1 N + 1]
        """
        G = matrix([[N**(p1.number_of_loops(p2)) for p1 in self]
                    for p2 in self])
        return G**(-1)

    Element = PerfectMatching
