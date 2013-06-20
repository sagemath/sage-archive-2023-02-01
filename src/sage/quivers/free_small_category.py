from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.magmas import Magmas
from sage.categories.monoids import Monoids
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

        sage: from sage.quivers.quiver import Quiver
        sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
        sage: F = Q.free_small_category()
        sage: F
        Free small category of Quiver on 3 vertices
        sage: F.category()
        Category of magmas
        sage: TestSuite(F).run()

    If there is only a single vertex, the free algebra is a monoid. ::

        sage: Q = Quiver({1:{1:['a','b', 'c', 'd']}})
        sage: F = Q.free_small_category()
        sage: F.category()
        Category of monoids
        sage: TestSuite(F).run()

    """
    @staticmethod
    def __classcall__(cls, *args, **kwds):
        if args:
            Q = args[0]
            if isinstance(Q,Quiver):
                return super(FreeSmallCategory, cls).__classcall__(cls,Q)
        Q = Quiver(*args,**kwds)
        return super(FreeSmallCategory, cls).__classcall__(cls,Q)

    Element = QuiverPath
    def __init__(self, Q):
        self._quiver = Q
        nvert = len(Q.vertices())
        if nvert == 0:
            raise TypeError("The quiver must not be empty")
        if nvert == 1:
            Parent.__init__(self, category=Monoids())
        else:
            Parent.__init__(self, category=Magmas())

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: from sage.quivers.quiver import Quiver
            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 2:{3:['d']}})
            sage: Q.free_small_category()
            Free small category of Quiver on 3 vertices

        """
        return "Free small category of %s"%self._quiver

    def _coerce_map_from_(self, other):
        if not isinstance(other, FreeSmallCategory):
            return
        return other.quiver().is_subgraph(self._quiver, induced=False)

    @cached_method
    def zero(self):
        """
        A free algebra has a unique "invalid" element, that behaves like a zero.

        EXAMPLES::

            sage: from sage.quivers.quiver import Quiver
            sage: Q = Quiver({1:{2:['a','b'], 3:['c']}, 3:{1:['d']}})
            sage: F = Q.free_small_category()
            sage: F.zero()
            invalid path

        """
        return self.element_class(None, parent=self)

    def is_finite(self):
        """
        This free small category is finite if and only if the underlying
        quiver is acyclic.

        EXAMPLES::

            sage: from sage.quivers.quiver import Quiver
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

            sage: from sage.quivers.quiver import Quiver
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

    def __iter__(self):
        """
        Iterate over the elements of self, i.e., over quiver paths

        EXAMPLES::

            sage: from sage.quivers.quiver import Quiver
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
        return self._quiver

    def reverse(self):
        return self._quiver.reverse().free_small_category()

    def algebra(self, R):
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

            sage: from sage.quivers.quiver import Quiver
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


#    def __contains__(self, other):
#        """
#        Implements the ``in`` keyword.
#
#        If other is a QuiverPath then ``other in self`` returns True if and only if
#        each edge of other is an edge of self.  In the case of a trivial path at a
#        vertex the vertex must be a vertex of the graph.  The invalid path belongs to
#        every quiver.  If other is not a QuiverPath then the call is passed up to
#        ``DiGraph.__contains__()`` which returns True if and only if other is a vertex.
#
#        TESTS::
#
#            sage: from sage.modules.quiver_module import Quiver, QuiverPath
#            sage: Q = Quiver({1: {2:['a']}, 2:{3:['b']}})
#            sage: x = QuiverPath([(1, 2, 'a'), (2, 3, 'b')])
#            sage: y = x*(3, 4, 'c')
#            sage: x in Q
#            True
#            sage: y in Q
#            False
#            sage: 1 in Q
#            True
#            sage: 4 in Q
#            False
#            sage: QuiverPath((1, 1)) in Q
#            True
#            sage: QuiverPath((4, 4)) in Q
#            False
#            sage: QuiverPath([(1, 1), (2, 2)]) in Q
#            True
#        """
#
#        # Pass up if not a QuiverPath
#        if not isinstance(other, QuiverPath):
#            return super(Quiver, self).__contains__(other)
#
#        # Empty paths are in every quiver
#        if not other._path:
#            return True
#
#        # If its a trivial path just check the vertex
#        if other._path[0][0] == other._path[0][1]:
#            return super(Quiver, self).__contains__(other._path[0][0])
#
#        # Otherwise check each edge
#        for e in other:
#            if (not super(Quiver, self).__contains__(e[0]) or
#                    not super(Quiver, self).__contains__(e[1]) or
#                    e[2] not in self.edge_label(e[0], e[1])):
#                return False
#        return True
