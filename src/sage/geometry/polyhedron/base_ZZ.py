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



from sage.rings.all import ZZ, QQ
from sage.misc.all import cached_method
from sage.matrix.constructor import matrix

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
            http://arxiv.org/abs/math/0204356
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

