"""
The Normaliz backend for polyhedral computations
"""
from __future__ import absolute_import

from sage.rings.all import ZZ, QQ
from sage.rings.integer import LCM_list
from sage.misc.functional import denominator
from sage.matrix.constructor import matrix

from .base import Polyhedron_base
from .base_QQ import Polyhedron_QQ
from .base_ZZ import Polyhedron_ZZ


#########################################################################
class Polyhedron_normaliz(Polyhedron_base):
    """
    Polyhedra with normaliz

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[], backend='normaliz')
        sage: TestSuite(p).run(skip='_test_pickling')

    Two ways to get the full space::

        sage: Polyhedron(eqns=[[0, 0, 0]], backend='normaliz')
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines
        sage: Polyhedron(ieqs=[[0, 0, 0]], backend='normaliz')
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines

    A lower-dimensional affine cone; we test that there are no mysterious
    inequalities coming in from the homogenization::

        sage: P = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)], backend='normaliz')
        sage: P.n_inequalities()
        1
        sage: P.equations()
        (An equation (1, 0) x - 1 == 0,)

    TESTS:

    Tests copied from various methods in base.py::

        sage: p = Polyhedron(vertices = [[1,0,0],[0,1,0],[0,0,1]], backend='normaliz')
        sage: p.n_equations()
        1
        sage: p.n_inequalities()
        3

        sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in range(6)], backend='normaliz')
        sage: p.n_facets()
        8

        sage: p = Polyhedron(vertices = [[1,0],[0,1],[1,1]], rays=[[1,1]], backend='normaliz')
        sage: p.n_vertices()
        2

        sage: p = Polyhedron(vertices = [[1,0],[0,1]], rays=[[1,1]], backend='normaliz')
        sage: p.n_rays()
        1

        sage: p = Polyhedron(vertices = [[0,0]], rays=[[0,1],[0,-1]], backend='normaliz')
        sage: p.n_lines()
        1

    """

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

            sage: p = Polyhedron(backend='normaliz')
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz
            sage: Polyhedron_normaliz._init_from_Vrepresentation(p, [], [], [])
        """
        import PyNormaliz
        if vertices is None:
            vertices = []
        nmz_vertices = []
        for v in vertices:
            d = LCM_list([denominator(v_i) for v_i in v])
            dv = [ d*v_i for v_i in v ]
            nmz_vertices.append(dv + [d])
        if rays is None:
            rays = []
        nmz_rays = []
        for r in rays:
            d = LCM_list([denominator(r_i) for r_i in r])
            dr = [ d*r_i for r_i in r ]
            nmz_rays.append(dr)
        if lines is None: lines = []
        nmz_lines = []
        for l in lines:
            d = LCM_list([denominator(l_i) for l_i in l])
            dl = [ d*l_i for l_i in l ]
            nmz_lines.append(dl)
        if not nmz_vertices and not nmz_rays and not nmz_lines:
            # Special case to avoid:
            #   error: Some error in the normaliz input data detected:
            #   All input matrices empty!
            self._init_empty_polyhedron()
        else:
            data = ["vertices", nmz_vertices,
                    "cone", nmz_rays,
                    "subspace", nmz_lines]
            self._normaliz_cone = PyNormaliz.NmzCone(data)
            assert self._normaliz_cone, "NmzCone({}) did not return a cone".format(data)
            self._init_Vrepresentation_from_normaliz(minimize)
            self._init_Hrepresentation_from_normaliz(minimize)

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

            sage: p = Polyhedron(backend='normaliz')
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz
            sage: Polyhedron_normaliz._init_from_Hrepresentation(p, [], [])
        """
        import PyNormaliz
        if ieqs is None: ieqs = []
        nmz_ieqs = []
        for ieq in ieqs:
            d = LCM_list([denominator(ieq_i) for ieq_i in ieq])
            dieq = [ ZZ(d*ieq_i) for ieq_i in ieq ]
            b = dieq[0]
            A = dieq[1:]
            nmz_ieqs.append(A + [b])
        if not nmz_ieqs:
            # If normaliz gets an empty list of inequalities, it adds
            # nonnegativities. So let's add a tautological inequality to work
            # around this.
            nmz_ieqs.append([0]*self.ambient_dim() + [0])
        if eqns is None: eqns = []
        nmz_eqns = []
        for eqn in eqns:
            d = LCM_list([denominator(eqn_i) for eqn_i in eqn])
            deqn = [ ZZ(d*eqn_i) for eqn_i in eqn ]
            b = deqn[0]
            A = deqn[1:]
            nmz_eqns.append(A + [b])
        data = ["inhom_equations", nmz_eqns,
                "inhom_inequalities", nmz_ieqs]
        self._normaliz_cone = PyNormaliz.NmzCone(data)
        assert self._normaliz_cone, "NmzCone({}) did not return a cone".format(data)
        self._init_Vrepresentation_from_normaliz(minimize)
        self._init_Hrepresentation_from_normaliz(minimize)

    def _init_Vrepresentation_from_normaliz(self, minimize):
        """
        Create the Vrepresentation objects from the normaliz polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],
            ...                  backend='normaliz')  # indirect doctest
            sage: set(p.Hrepresentation())
            {An inequality (1, 4) x - 2 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (-5, 12) x + 10 >= 0}
            sage: set(p.Vrepresentation())
            {A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6)}

        """
        import PyNormaliz
        self._Vrepresentation = []
        parent = self.parent()
        base_ring = self.base_ring()
        cone = self._normaliz_cone
        for g in PyNormaliz.NmzResult(cone, "VerticesOfPolyhedron"):
            d = g[-1]
            if d == 1:
                parent._make_Vertex(self, g[:-1])
            else:
                parent._make_Vertex(self, [base_ring(x)/d for x in g[:-1]])
        for g in PyNormaliz.NmzResult(cone, "ExtremeRays"):
            parent._make_Ray(self, g[:-1])
        for g in PyNormaliz.NmzResult(cone, "MaximalSubspace"):
            parent._make_Line(self, g[:-1])
        self._Vrepresentation = tuple(self._Vrepresentation)

    def _init_Hrepresentation_from_normaliz(self, minimize):
        """
        Create the Hrepresentation objects from the normaliz polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],
            ...                  backend='normaliz')  # indirect doctest
            sage: set(p.Hrepresentation())
            {An inequality (1, 4) x - 2 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (-5, 12) x + 10 >= 0}
            sage: set(p.Vrepresentation())
            {A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6)}

        """
        import PyNormaliz
        self._Hrepresentation = []
        base_ring = self.base_ring()
        cone = self._normaliz_cone
        parent = self.parent()
        for g in PyNormaliz.NmzResult(cone, "SupportHyperplanes"):
            if all(x==0 for x in g[:-1]):
                # Ignore vertical inequality
                pass
            else:
                parent._make_Inequality(self, (g[-1],) + tuple(g[:-1]))
        for g in PyNormaliz.NmzResult(cone, "Equations"):
            parent._make_Equation(self, (g[-1],) + tuple(g[:-1]))
        self._Hrepresentation = tuple(self._Hrepresentation)

    def _init_empty_polyhedron(self):
        """
        Initializes an empty polyhedron.

        TESTS::

            sage: empty = Polyhedron(backend='normaliz'); empty
            The empty polyhedron in ZZ^0
            sage: empty.Vrepresentation()
            ()
            sage: empty.Hrepresentation()
            (An equation -1 == 0,)
            sage: Polyhedron(vertices = [], backend='normaliz')
            The empty polyhedron in ZZ^0
            sage: Polyhedron(backend='normaliz')._init_empty_polyhedron()
        """
        super(Polyhedron_normaliz, self)._init_empty_polyhedron()
        # Can't seem to set up an empty _normaliz_cone.
        # For example, PyNormaliz.NmzCone(['vertices', []]) gives
        # error: Some error in the normaliz input data detected: All input matrices empty!
        self._normaliz_cone = None


#########################################################################
class Polyhedron_QQ_normaliz(Polyhedron_normaliz, Polyhedron_QQ):
    """
    Polyhedra over `\QQ` with normaliz

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[],
        ...                  backend='normaliz', base_ring=QQ)
        sage: TestSuite(p).run(skip='_test_pickling')
    """
    pass


#########################################################################
class Polyhedron_ZZ_normaliz(Polyhedron_normaliz, Polyhedron_ZZ):
    """
    Polyhedra over `\ZZ` with normaliz

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[])
        ...                  backend='normaliz', base_ring=ZZ)
        sage: TestSuite(p).run(skip='_test_pickling')
    """
    pass
