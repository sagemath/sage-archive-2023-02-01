"""
The Python backend

EXAMPLES::

    sage: p0 = (0, 0)
    sage: p1 = (1, 0)
    sage: p2 = (1/2, AA(3).sqrt()/2)
    sage: equilateral_triangle = Polyhedron([p0, p1, p2])
"""

from sage.rings.all import ZZ, QQ
from sage.rings.integer import LCM_list
from sage.matrix.constructor import matrix
from base import Polyhedron_base


class Polyhedron_field(Polyhedron_base):
    """
    Polyhedra over all field supported by Sage

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(AA(2).sqrt(),0),(0,AA(3).sqrt())], 
        ....:                rays=[(1,1)], lines=[], backend='field')
        sage: TestSuite(p).run()

    TESTS::

        sage: p = Polyhedron([(0,sqrt(2))], base_ring=AA);  p
        A 0-dimensional polyhedron in AA^2 defined as the convex hull of 1 vertex
        sage: TestSuite(p).run()
    """
    def _is_zero(self, x):
        """
        Test whether ``x`` is zero.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)
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

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)
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

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)
            sage: p._is_positive(1)
            True
            sage: p._is_positive(0)
            False
        """
        return x>0

    def _init_from_Vrepresentation(self, vertices, rays, lines, minimize=True, verbose=False):
        """
        Construct polyhedron from V-representation data.

        INPUT:

        - ``vertices`` -- list of point. Each point can be specified
           as any iterable container of
           :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``rays`` -- list of rays. Each ray can be specified as any
          iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``lines`` -- list of lines. Each line can be specified as
          any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``verbose`` -- boolean (default: ``False``). Whether to print
          verbose output for debugging purposes.

        EXAMPLES::

            sage: p = Polyhedron(backend='field')
            sage: from sage.geometry.polyhedron.backend_field import Polyhedron_field
            sage: Polyhedron_field._init_from_Vrepresentation(p, [], [], [])
        """
        from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep, Vrep2Hrep
        H = Vrep2Hrep(self.base_ring(), self.ambient_dim(), vertices, rays, lines)
        V = Hrep2Vrep(self.base_ring(), self.ambient_dim(), 
                      H.inequalities, H.equations)
        self._init_Vrepresentation_backend(V)
        self._init_Hrepresentation_backend(H)

    def _init_from_Hrepresentation(self, ieqs, eqns, minimize=True, verbose=False):
        """
        Construct polyhedron from H-representation data.

        INPUT:

        - ``ieqs`` -- list of inequalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``eqns`` -- list of equalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``verbose`` -- boolean (default: ``False``). Whether to print
          verbose output for debugging purposes.

        EXAMPLES::

            sage: p = Polyhedron(backend='field')
            sage: from sage.geometry.polyhedron.backend_field import Polyhedron_field
            sage: Polyhedron_field._init_from_Hrepresentation(p, [], [])
        """
        from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep, Vrep2Hrep
        assert len(eqns)==0   #TODO
        V = Hrep2Vrep(self.base_ring(), self.ambient_dim(), ieqs)
        H = Vrep2Hrep(self.base_ring(), self.ambient_dim(), 
                      V.vertices, V.rays, V.lines)
        self._init_Vrepresentation_backend(V)
        self._init_Hrepresentation_backend(H)

    def _init_Vrepresentation_backend(self, Vrep):
        """
        Create the Vrepresentation objects from the ppl polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/sqrt(2)),(sqrt(2),0),(4,sqrt(5)/6)],
            ...                  base_ring=AA, backend='field')  # indirect doctest
            sage: p.Hrepresentation()
            (An inequality (-0.1582178750233332?, 1.097777812326429?) x + 0.2237538646678492? >= 0,
             An inequality (-0.1419794359520263?, -1.698172434277148?) x + 1.200789243901438? >= 0,
             An inequality (0.3001973109753594?, 0.600394621950719?) x - 0.4245431085692869? >= 0)
            sage: p.Vrepresentation()
            (A vertex at (0.?e-15, 0.707106781186548?),
             A vertex at (1.414213562373095?, 0),
             A vertex at (4.000000000000000?, 0.372677996249965?))
        """
        self._Vrepresentation = []
        parent = self.parent()
        for v in Vrep.vertices:
            parent._make_Vertex(self, v)
        for r in Vrep.rays:
            parent._make_Ray(self, r)
        for l in Vrep.lines:
            parent._make_Line(self, l)
        self._Vrepresentation = tuple(self._Vrepresentation)

    def _init_Hrepresentation_backend(self, Hrep):
        """
        Create the Hrepresentation objects from the ppl polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/sqrt(2)),(sqrt(2),0),(4,sqrt(5)/6)],
            ...                  base_ring=AA, backend='field')  # indirect doctest
            sage: p.Hrepresentation()
            (An inequality (-0.1582178750233332?, 1.097777812326429?) x + 0.2237538646678492? >= 0,
             An inequality (-0.1419794359520263?, -1.698172434277148?) x + 1.200789243901438? >= 0,
             An inequality (0.3001973109753594?, 0.600394621950719?) x - 0.4245431085692869? >= 0)
            sage: p.Vrepresentation()
            (A vertex at (0.?e-15, 0.707106781186548?),
             A vertex at (1.414213562373095?, 0),
             A vertex at (4.000000000000000?, 0.372677996249965?))
        """
        self._Hrepresentation = []
        parent = self.parent()
        for ieq in Hrep.inequalities:
            parent._make_Inequality(self, ieq)
        for eqn in Hrep.equations:
            parent._make_Equation(self, eqn)
        self._Hrepresentation = tuple(self._Hrepresentation)

    def _init_empty_polyhedron(self):
        """
        Initializes an empty polyhedron.

        TESTS::

            sage: empty = Polyhedron(backend='field', base_ring=AA); empty
            The empty polyhedron in AA^0
            sage: empty.Vrepresentation()
            ()
            sage: empty.Hrepresentation()
            (An equation -1 == 0,)
            sage: Polyhedron(vertices = [], backend='field')
            The empty polyhedron in QQ^0
            sage: Polyhedron(backend='field')._init_empty_polyhedron()
        """
        super(Polyhedron_field, self)._init_empty_polyhedron()


