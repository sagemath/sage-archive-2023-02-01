from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.magmas import Magmas
from sage.categories.monoids import Monoids
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.cachefunc import cached_method
from paths import QuiverPath
from quiver import Quiver

class FreeSmallCategory(UniqueRepresentation, Parent):
    """
    The free small category that is associated with a quiver.

    .. NOTE::

        The vertices of the quiver correspond to the objects of the category,
        which thus form a set. The directed paths in the quiver correspond
        to morphisms between objects of the category.

    EXAMPLES::

        sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
        sage: F = Q.free_small_category()
        sage: F
        Free small category of Quiver on 3 vertices
        sage: F.variable_names()
        ('e_1', 'e_2', 'e_3', 'a', 'b', 'c', 'd')
        sage: F.gens()
        (e_1, e_2, e_3, a, b, c, d)
        sage: F.category()
        Join of Category of finite enumerated sets and Category of magmas
        sage: TestSuite(F).run()

    If there is only a single vertex, the free small category is a monoid. If
    the underlying quiver has cycles or loops, then the free small category is
    of course only an infinite enumerated set::

        sage: Q = Quiver({1:{1:['a','b', 'c', 'd']}})
        sage: F = Q.free_small_category()
        sage: F.category()
        Join of Category of infinite enumerated sets and Category of monoids
        sage: TestSuite(F).run()

    """
    @staticmethod
    def __classcall__(cls, *args, **kwds):
        """
        Preprocessing: Convert input into a quiver before using the cache.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: from sage.quivers.free_small_category import FreeSmallCategory
            sage: F = FreeSmallCategory({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: F is FreeSmallCategory(Q) # indirect doctest
            True
            sage: F is loads(dumps(F))
            True

        """
        if args:
            Q = args[0]
            if isinstance(Q,Quiver):
                return super(FreeSmallCategory, cls).__classcall__(cls,Q)
        Q = Quiver(*args,**kwds)
        return super(FreeSmallCategory, cls).__classcall__(cls,Q)

    Element = QuiverPath
    def __init__(self, Q):
        """
        INPUT:

        - a :class:`~sage.quivers.quiver.Quiver`.

        NOTE:

        There is a pre-processing. The class' constructor will also
        accept data which a quiver can be created from.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: from sage.quivers.free_small_category import FreeSmallCategory
            sage: F = FreeSmallCategory({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: F is FreeSmallCategory(Q) # indirect doctest
            True
            sage: F
            Free small category of Quiver on 3 vertices

        """
        self._quiver = Q
        nvert = len(Q.vertices())
        if nvert == 0:
            raise TypeError("The quiver must not be empty")
        if Q.is_acyclic_quiver():
            cat = FiniteEnumeratedSets()
        else:
            cat = InfiniteEnumeratedSets()
        names=['e_{0}'.format(v) for v in Q.vertices()] + [e[2] for e in Q.edges()]
        if nvert == 1:
            Parent.__init__(self, names=names, category=cat.join([cat,Monoids()]))
        else:
            Parent.__init__(self, names=names, category=cat.join([cat,Magmas()]))

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: Q.free_small_category()
            Free small category of Quiver on 3 vertices

        """
        return "Free small category of %s"%self._quiver

    def _coerce_map_from_(self, other):
        """
        A coercion exists from A to B exists, if the underlying quiver
        of A is a sub-quiver of the underlying quiver of B (preserving
        names).

        EXAMPLES::

            sage: Q1 = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: Q2 = Quiver({1:{2:['a'], 3:['c']}})
            sage: Q3 = Quiver({1:{2:['a','x'], 3:['c']}, 3:{1:['d']}})
            sage: F1 = Q1.free_small_category()
            sage: F2 = Q2.free_small_category()
            sage: F3 = Q3.free_small_category()
            sage: F1.has_coerce_map_from(F2)   # indirect doctest
            True
            sage: F1.has_coerce_map_from(F3)
            False
            sage: d = F1([(3,1,'d')]); d
            d
            sage: c = F2([(1,3,'c')]); c
            c
            sage: c.parent() is F1
            False
            sage: c in F1    # indirect doctest
            True
            sage: d*c        # indirect doctest
            d*c
            sage: (d*c).parent() is F1
            True
            sage: c3 = F3([(1,3,'c')]); c3
            c
            sage: c3 in F1   # indirect doctest
            False
            sage: d*c3
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Free small category of Quiver on 3 vertices' and 'Free small category of Quiver on 3 vertices'

        """
        if not isinstance(other, FreeSmallCategory):
            return
        return other.quiver().is_subgraph(self._quiver, induced=False)

    @cached_method
    def zero(self):
        """
        A free algebra has a unique "invalid" element, that behaves like a zero.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: F = Q.free_small_category()
            sage: F.zero()
            invalid path

        """
        return self.element_class(None, parent=self)

    @cached_method
    def arrows(self):
        """
        The elements corresponding to edges of the underlying quiver.

        EXAMPLES::

            sage: F = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}}).free_small_category()
            sage: F.arrows()
            (a, b, c, d)
        """
        return tuple(self([e]) for e in self._quiver.edges())

    @cached_method
    def idempotents(self):
        """
        The idempotents corresponding to the vertices of the underlying quiver.

        EXAMPLES::

            sage: F = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}}).free_small_category()
            sage: F.idempotents()
            (e_1, e_2, e_3)

        """
        return tuple(self((v,v)) for v in self._quiver.vertices())

    def ngens(self):
        """
        The number of generators (:meth:`arrows` and :meth:`idempotents`)

        EXAMPLES::

            sage: F = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}}).free_small_category()
            sage: F.ngens()
            7

        """
        Q = self._quiver
        return Q.num_verts() + Q.num_edges()

    @cached_method
    def gen(self, i):
        """
        Generator number i

        INPUT:

        i, and integer

        OUTPUT:

        An idempotent, if i is smaller than the number of vertices,
        or an arrow otherwise.

        EXAMPLES::

            sage: F = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}}).free_small_category()
            sage: F.1         # indirect doctest
            e_2
            sage: F.idempotents()[1]
            e_2
            sage: F.5
            c
            sage: F.gens()[5]
            c

        """
        Q = self._quiver
        nv = Q.num_verts()
        if i < nv:
            v = Q.vertices()[i]
            return self((v,v))
        e = Q.edges()[i-nv]
        return self([e])

    def gens(self):
        """
        The tuple of generators.

        NOTE:

        This coincides with the sum of the output of :meth:`idempotents`
        and :meth:`arrows`.

        EXAMPLES::

            sage: F = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}}).free_small_category()
            sage: F.gens()
            (e_1, e_2, e_3, a, b, c, d)
            sage: F.gens() == F.idempotents() + F.arrows()
            True

        """
        return self.idempotents() + self.arrows()

    def is_finite(self):
        """
        This free small category is finite if and only if the underlying
        quiver is acyclic.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: Q.free_small_category().is_finite()
            True
            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: Q.free_small_category().is_finite()
            False

        """
        return self._quiver.is_acyclic_quiver()

    def __len__(self):
        """
        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: F = Q.free_small_category()
            sage: len(F)
            9
            sage: list(F)
            [e_1, e_2, e_3, b, a, c, d, b*d, a*d]
            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: F = Q.free_small_category()
            sage: len(F)
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite list

        """
        return len(self.all_paths())

    def cardinality(self):
        """
        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: F = Q.free_small_category()
            sage: F.cardinality()
            9
            sage: A = F.algebra(QQ)
            sage: list(A.basis())
            [e_1, e_2, e_3, b, a, c, d, b*d, a*d]
            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: F = Q.free_small_category()
            sage: F.cardinality()
            +Infinity
            sage: A = F.algebra(QQ)
            sage: list(A.basis())
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite list

        """
        from sage.all import ZZ
        if self._quiver.is_acyclic_quiver():
            return ZZ(len(self))
        from sage.rings.infinity import Infinity
        return Infinity

    def __iter__(self):
        """
        Iterate over the elements of self, i.e., over quiver paths

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: F = Q.free_small_category()
            sage: list(F)
            [e_1, e_2, e_3, b, a, c, d, b*d, a*d]

        The elements are sorted by length. Of course, the list of elements
        is infinite for quivers with cycles. ::

            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: F = Q.free_small_category()
            sage: F.is_finite()
            False
            sage: list(F)
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite list

         However, one can iterate::

            sage: counter = 0
            sage: for p in F:
            ....:     counter += 1
            ....:     print p
            ....:     if counter==20:
            ....:         break
            e_1
            e_2
            e_3
            b
            a
            d
            c
            b*d
            a*d
            d*c
            c*b
            c*a
            b*d*c
            a*d*c
            d*c*b
            d*c*a
            c*b*d
            c*a*d
            b*d*c*b
            b*d*c*a

        """
        d = 0
        length_d_available = True
        while length_d_available:
            length_d_available = False
            for v in self._quiver.vertices():
                for w in self.iter_paths_by_length_and_startpoint(d,v):
                    length_d_available = True
                    yield w
            d += 1

    def iter_paths_by_length_and_startpoint(self, d, v):
        """
        An iterator over quiver paths with a fixed length and start point.

        INPUT:

        - d, an integer, path length
        - v, a vertex, start point of the paths

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: F = Q.free_small_category()
            sage: F.is_finite()
            False
            sage: list(F.iter_paths_by_length_and_startpoint(4,1))
            [b*d*c*b, b*d*c*a, a*d*c*b, a*d*c*a]
            sage: list(F.iter_paths_by_length_and_startpoint(5,1))
            [b*d*c*b*d, b*d*c*a*d, a*d*c*b*d, a*d*c*a*d]
            sage: list(F.iter_paths_by_length_and_startpoint(5,2))
            [d*c*b*d*c, d*c*a*d*c]

        TEST::

             sage: Q = Quiver({1:{1:['a','b', 'c', 'd']}})
             sage: F = Q.free_small_category()
             sage: list(F.iter_paths_by_length_and_startpoint(2,1))
             [d*d,
              d*c,
              d*b,
              d*a,
              c*d,
              c*c,
              c*b,
              c*a,
              b*d,
              b*c,
              b*b,
              b*a,
              a*d,
              a*c,
              a*b,
              a*a]
             sage: len(list(F.iter_paths_by_length_and_startpoint(2,1)))
             16

        """
        # iterate over length d paths starting at vertex v
        if not d>=0:
            raise ValueError("Path length must be a non-negative integer")
        if v not in self._quiver:
            raise ValueError("The starting point %s is not a vertex of the underlying quiver"%v)
        if not d:
            yield self([(v,v)],check=False)
        else:
            for w in self.iter_paths_by_length_and_startpoint(d-1, v):
                for a in self._quiver._backend.iterator_out_edges([w.terminal_vertex()],True):
                    yield self(list(w)+[a],check=False)

    def iter_paths_by_length_and_endpoint(self, d, v):
        """
        An iterator over quiver paths with a fixed length and end point.

        INPUT:

        - d, an integer, path length
        - v, a vertex, end point of the paths

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: F = Q.free_small_category()
            sage: F.is_finite()
            False
            sage: list(F.iter_paths_by_length_and_endpoint(4,1))
            [c*b*d*c, c*a*d*c]
            sage: list(F.iter_paths_by_length_and_endpoint(5,1))
            [d*c*b*d*c, d*c*a*d*c]
            sage: list(F.iter_paths_by_length_and_endpoint(5,2))
            [c*b*d*c*b, c*a*d*c*b, c*b*d*c*a, c*a*d*c*a]

        """
        # iterate over length d paths ending at vertex v
        if not d>=0:
            raise ValueError("Path length must be a non-negative integer")
        if v not in self._quiver:
            raise ValueError("The starting point %s is not a vertex of the underlying quiver"%v)
        if not d:
            yield self([(v,v)],check=False)
        else:
            for w in self.iter_paths_by_length_and_endpoint(d-1, v):
                for a in self._quiver._backend.iterator_in_edges([w.initial_vertex()],True):
                    yield self([a]+list(w),check=False)

    def quiver(self):
        """
        The underlying quiver of this free small category.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: F = Q.free_small_category()
            sage: F.quiver() is Q
            True
        """
        return self._quiver

    def reverse(self):
        """
        The free small category of the reverse quiver.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: F = Q.free_small_category()
            sage: F.reverse() is Q.reverse().free_small_category()
            True

        """
        return self._quiver.reverse().free_small_category()

    def algebra(self, R):
        """
        Path algebra of the underlying quiver.

        INPUT:

        R, a commutative ring

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b']}, 2:{3:['d']}, 3:{1:['c']}})
            sage: F = Q.free_small_category()
            sage: F.algebra(GF(3)) is Q.algebra(GF(3))
            True

        """
        return self._quiver.algebra(R)

    def all_paths(self, start=None, end=None):
        """
        List of all paths between a pair of vertices (start, end).

        INPUT:

        - ``start`` - integer or None (default: None), the initial vertex of
          the paths in the output.  If None is given then the initial vertex
          is arbitrary.
        - ``end`` - integer or None (default: None), the terminal vertex of
          the paths in the output.  If None is given then the terminal vertex
          is arbitrary.

        OUTPUT:

        - list of paths, excluding the invalid path

        .. NOTE::

            If there are multiple edges between two vertices, the method
            :meth:`sage.graphs.digraph.all_paths` will not differentiate
            between them. But this method, which is not for digraphs but for
            the free associative magma associated with it, will.

        EXAMPLES::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: F = Q.free_small_category()
            sage: F.all_paths(1, 3)
            [b*d, a*d, c]

        If start=end then we expect only the trivial path at that vertex::

            sage: F.all_paths(1, 1)
            [e_1]

        The empty list is returned if there are no paths between the given vertices::

            sage: F.all_paths(3, 1)
            []

        If end=None then all edge paths beginning at start are returned, including the
        trivial path::

            sage: F.all_paths(2)
            [e_2, d]

        If start=None then all edge paths ending at end are returned, including the
        trivial path.  Note that the two edges from vertex 1 to vertex 2 count as two
        different edge paths::

            sage: F.all_paths(None, 2)
            [b, a, e_2]
            sage: F.all_paths(end=2)
            [b, a, e_2]

        If start=end=None then all edge paths are returned, including trivial paths::

            sage: F.all_paths()
            [e_1, b, a, b*d, a*d, c, e_2, d, e_3]

        The vertex given must be a vertex of the quiver::

            sage: F.all_paths(1, 4)
            Traceback (most recent call last):
            ...
            ValueError: The end vertex 4 is not a vertex of the quiver.

        If the underlying quiver is cyclic, a ``NotImplementedError`` is raised::

            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: F = Q.free_small_category()
            sage: F.all_paths()
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite list

        """
        if not self.is_finite():
            raise NotImplementedError("infinite list")
        # Check that given arguments are vertices
        if start is not None and start not in self._quiver:
            raise ValueError("The start vertex " + str(start) + " is not a vertex of the quiver.")
        if end is not None and end not in self._quiver:
            raise ValueError("The end vertex " + str(end) + " is not a vertex of the quiver.")

        # Handle quivers with cycles
        Q = self._quiver
        if not Q.is_acyclic_quiver():
            raise ValueError("The underlying quiver has cycles, thus, there may be an infinite of directed paths")

        # Handle start=None
        if start == None:
            results = []
            for v in Q:
                results += self.all_paths(v, end)
            return results

        # Handle end=None
        if end == None:
            results = []
            for v in Q:
                results += self.all_paths(start, v)
            return results

        # Handle the trivial case
        if start == end:
            return [self([(start, end)],check=False)]

        # This function will recursively convert a path given in terms of
        # vertices to a list of QuiverPaths.
        def _v_to_e(path):
            if len(path) == 1:
                return [self([(path[0], path[0])],check=False)]
            paths = []
            for a in Q.edge_label(path[0], path[1]):
                for b in _v_to_e(path[1:]):
                    paths.append(self([(path[0], path[1], a)]+list(b),check=False))
            return paths

        # For each vertex path we append the resulting edge paths
        result = []
        for path in Q.all_paths(start, end):
            result += _v_to_e(path)

        # The result is all paths from start to end
        return result
