r"""
Base class for polyhedra over `\ZZ`
"""

########################################################################
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################



from sage.rings.all import ZZ, QQ, gcd
from sage.misc.all import cached_method
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from constructor import Polyhedron
from base import Polyhedron_base



#########################################################################
class Polyhedron_ZZ(Polyhedron_base):
    """
    Base class for Polyhedra over `\ZZ`

    TESTS::

        sage: p = Polyhedron([(0,0)], base_ring=ZZ);  p
        A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex
        sage: TestSuite(p).run(skip='_test_pickling')
    """
    def _is_zero(self, x):
        """
        Test whether ``x`` is zero.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=ZZ)
            sage: p._is_zero(0)
            True
            sage: p._is_zero(1/100000)
            False
        """
        return x==0

    def _is_nonneg(self, x):
        """
        Test whether ``x`` is nonnegative.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=ZZ)
            sage: p._is_nonneg(1)
            True
            sage: p._is_nonneg(-1/100000)
            False
        """
        return x>=0

    def _is_positive(self, x):
        """
        Test whether ``x`` is positive.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], base_ring=ZZ)
            sage: p._is_positive(1)
            True
            sage: p._is_positive(0)
            False
        """
        return x>0

    _base_ring = ZZ

    def is_lattice_polytope(self):
        r"""
        Return whether the polyhedron is a lattice polytope.

        OUTPUT:

        ``True`` if the polyhedron is compact and has only integral
        vertices, ``False`` otherwise.

        EXAMPLES::

            sage: polytopes.cross_polytope(3).is_lattice_polytope()
            True
            sage: polytopes.regular_polygon(5).is_lattice_polytope()
            False
        """
        return True

    @cached_method
    def polar(self):
        """
        Return the polar (dual) polytope.

        The polytope must have the IP-property (see
        :meth:`has_IP_property`), that is, the origin must be an
        interior point. In particular, it must be full-dimensional.

        OUTPUT:

        The polytope whose vertices are the coefficient vectors of the
        inequalities of ``self`` with inhomogeneous term normalized to
        unity.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(1,0,0),(0,1,0),(0,0,1),(-1,-1,-1)], base_ring=ZZ)
            sage: p.polar()
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: type(_)
            <class 'sage.geometry.polyhedron.backend_ppl.Polyhedra_ZZ_ppl_with_category.element_class'>
            sage: p.polar().base_ring()
            Integer Ring
        """
        if not self.has_IP_property():
            raise ValueError('The polytope must have the IP property.')

        vertices = [ ieq.A()/ieq.b() for
                     ieq in self.inequality_generator() ]
        if all( all(v_i in ZZ for v_i in v) for v in vertices):
            return Polyhedron(vertices=vertices, base_ring=ZZ)
        else:
            return Polyhedron(vertices=vertices, base_ring=QQ)

    @cached_method
    def is_reflexive(self):
        """
        EXAMPLES::

            sage: p = Polyhedron(vertices=[(1,0,0),(0,1,0),(0,0,1),(-1,-1,-1)], base_ring=ZZ)
            sage: p.is_reflexive()
            True
        """
        return self.polar().is_lattice_polytope()

    @cached_method
    def has_IP_property(self):
        """
        Test whether the polyhedron has the IP property.

        The IP (interior point) property means that

        * ``self`` is compact (a polytope).

        * ``self`` contains the origin as an interior point.

        This implies that

        * ``self`` is full-dimensional.

        * The dual polyhedron is again a polytope (that is, a compact
          polyhedron), though not necessarily a lattice polytope.

        EXAMPLES::

            sage: Polyhedron([(1,1),(1,0),(0,1)], base_ring=ZZ).has_IP_property()
            False
            sage: Polyhedron([(0,0),(1,0),(0,1)], base_ring=ZZ).has_IP_property()
            False
            sage: Polyhedron([(-1,-1),(1,0),(0,1)], base_ring=ZZ).has_IP_property()
            True

        REFERENCES:

        ..  [PALP]
            Maximilian Kreuzer, Harald Skarke:
            "PALP: A Package for Analyzing Lattice Polytopes
            with Applications to Toric Geometry"
            Comput.Phys.Commun. 157 (2004) 87-106
            :arxiv:`math/0204356`
        """
        return self.is_compact() and self.interior_contains(self.ambient_space().zero())

    def fibration_generator(self, dim):
        """
        Generate the lattice polytope fibrations.

        For the purposes of this function, a lattice polytope fiber is
        a sub-lattice polytope. Projecting the plane spanned by the
        subpolytope to a point yields another lattice polytope, the
        base of the fibration.

        INPUT:

        - ``dim`` -- integer. The dimension of the lattice polytope
          fiber.

        OUTPUT:

        A generator yielding the distinct lattice polytope fibers of
        given dimension.

        EXAMPLES::

            sage: P = Polyhedron(toric_varieties.P4_11169().fan().rays(), base_ring=ZZ)
            sage: list( P.fibration_generator(2) )
            [A 2-dimensional polyhedron in ZZ^4 defined as the convex hull of 3 vertices]
        """
        from sage.combinat.combination import Combinations
        if not self.is_compact():
            raise ValueError('Only polytopes (compact polyhedra) are allowed.')

        nonzero_points = [p for p in self.integral_points() if not p.is_zero()]
        origin = [[0]*self.ambient_dim()]
        fibers = set()
        parent = self.parent()

        for points in Combinations(nonzero_points, dim):
                plane = parent.element_class(parent, [origin,[],points], None)
                if plane.dim() != dim:
                    continue
                fiber = self.intersection(plane)
                if fiber.base_ring() is not ZZ:
                    continue
                fiber_vertices = tuple(sorted(tuple(v) for v in fiber.vertex_generator()))
                if fiber_vertices not in fibers:
                    yield fiber
                    fibers.update([fiber_vertices])
                plane.delete()

    def find_translation(self, translated_polyhedron):
        """
        Return the translation vector to ``translated_polyhedron``.

        INPUT:

        - ``translated_polyhedron`` -- a polyhedron.

        OUTPUT:

        A `\ZZ`-vector that translates ``self`` to
        ``translated_polyhedron``. A ``ValueError`` is raised if
        ``translated_polyhedron`` is not a translation of ``self``,
        this can be used to check that two polyhedra are not
        translates of each other.

        EXAMPLES::

            sage: X = polytopes.n_cube(3)
            sage: X.find_translation(X + vector([2,3,5]))
            (2, 3, 5)
            sage: X.find_translation(2*X)
            Traceback (most recent call last):
            ...
            ValueError: polyhedron is not a translation of self
        """
        no_translation_exception = ValueError('polyhedron is not a translation of self')
        if ( set(self.rays()) != set(translated_polyhedron.rays()) or
             set(self.lines()) != set(translated_polyhedron.lines()) or
             self.n_vertices() != translated_polyhedron.n_vertices() ):
            raise no_translation_exception
        sorted_vertices = sorted(map(vector, self.vertices()))
        sorted_translated_vertices = sorted(map(vector, translated_polyhedron.vertices()))
        v = sorted_translated_vertices[0] - sorted_vertices[0]
        if any(vertex+v != translated_vertex
               for vertex, translated_vertex in zip(sorted_vertices, sorted_translated_vertices)):
            raise no_translation_exception
        return v

    def _subpoly_parallel_facets(self):
        """
        Generator for all lattice sub-polyhedra with parallel facets.

        In a sub-polyhedron `Y\subset X` not all edges of `Y` need to
        be parallel to `X`. This method iterates over all
        sub-polyhedra where they are parallel, up to an overall
        translation of the sub-polyhedron. Degenerate sub-polyhedra of
        dimension strictly smaller are included.

        OUTPUT:

        A generator yielding `\ZZ`-polyhedra. By construction, each
        facet of the returned polyhedron is parallel to one of the
        facets of ``self``.

        EXAMPLES::

            sage: X = Polyhedron(vertices=[(0,0), (0,1), (1,0), (1,1)])
            sage: X._subpoly_parallel_facets()
            <generator object _subpoly_parallel_facets at 0x...>
            sage: for p in X._subpoly_parallel_facets():
            ...       print p.Vrepresentation()
            (A vertex at (0, 0),)
            (A vertex at (0, -1), A vertex at (0, 0))
            (A vertex at (-1, 0), A vertex at (0, 0))
            (A vertex at (-1, -1), A vertex at (-1, 0), A vertex at (0, -1), A vertex at (0, 0))

        TESTS::

            sage: X = Polyhedron(vertices=[(0,), (3,)])
            sage: [ p.vertices() for p in X._subpoly_parallel_facets() ]
            [(A vertex at (0),),
             (A vertex at (-1), A vertex at (0)),
             (A vertex at (-2), A vertex at (0)),
             (A vertex at (-3), A vertex at (0))]
            sage: list( Polyhedron(vertices=[[0,0]])._subpoly_parallel_facets() )
            [A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex]
            sage: list( Polyhedron()._subpoly_parallel_facets() )
            [The empty polyhedron in ZZ^0]
        """
        if self.dim()>2 or not self.is_compact():
            raise NotImplementedError('only implemented for bounded polygons')
        from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d
        vertices = cyclic_sort_vertices_2d(self.vertices())
        n = len(vertices)
        if n==1:  # single point
            yield self
            return
        edge_vectors = []
        for i in range(0,n):
            v = vertices[(i+1) % n].vector() - vertices[i].vector()
            d = gcd(list(v))
            v_prim = (v/d).change_ring(ZZ)
            edge_vectors.append([ v_prim*i for i in range(d+1) ])
        origin = self.ambient_space().zero()
        parent = self.parent()
        from sage.combinat.cartesian_product import CartesianProduct
        for edges in CartesianProduct(*edge_vectors):
            v = []
            point = origin
            for e in edges:
                point += e
                v.append(point)
            if point!=origin:   # does not close up, not a subpolygon
                continue
            yield parent([v, [], []], None)

    @cached_method
    def Minkowski_decompositions(self):
        """
        Return all Minkowski sums that add up to the polyhedron.

        OUTPUT:

        A tuple consisting of pairs `(X,Y)` of `\ZZ`-polyhedra that
        add up to ``self``. All pairs up to exchange of the summands
        are returned, that is, `(Y,X)` is not included if `(X,Y)`
        already is.

        EXAMPLES::

            sage: square = Polyhedron(vertices=[(0,0),(1,0),(0,1),(1,1)])
            sage: square.Minkowski_decompositions()
            ((A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices),
             (A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices,
              A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices))

        Example from http://cgi.di.uoa.gr/~amantzaf/geo/ ::

            sage: Q = Polyhedron(vertices=[(4,0), (6,0), (0,3), (4,3)])
            sage: R = Polyhedron(vertices=[(0,0), (5,0), (8,4), (3,2)])
            sage: (Q+R).Minkowski_decompositions()
            ((A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 7 vertices),
             (A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices),
             (A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 7 vertices),
             (A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 5 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices),
             (A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 7 vertices),
             (A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 5 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices),
             (A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 7 vertices),
             (A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices,
              A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 6 vertices))

           sage: [ len(square.dilation(i).Minkowski_decompositions())
           ...     for i in range(6) ]
           [1, 2, 5, 8, 13, 18]
           sage: [ ceil((i^2+2*i-1)/2)+1 for i in range(10) ]
           [1, 2, 5, 8, 13, 18, 25, 32, 41, 50]
        """
        if self.dim()>2 or not self.is_compact():
            raise NotImplementedError('only implemented for bounded polygons')
        summands = []
        def is_known_summand(poly):
            for summand in summands:
                try:
                    poly.find_translation(summand)
                    return True
                except ValueError:
                    pass
        decompositions = []
        for X in self._subpoly_parallel_facets():
            if is_known_summand(X):
                continue
            Y = self - X
            if X+Y != self:
                continue
            decompositions.append((X,Y))
            summands += [X, Y]
        return tuple(decompositions)

