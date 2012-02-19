"""
The PPL (Parma Polyhedra Library) backend for polyhedral computations
"""

from sage.rings.all import ZZ, QQ
from sage.rings.arith import lcm
from sage.misc.functional import denominator
from sage.matrix.constructor import matrix
from sage.libs.ppl import (
    C_Polyhedron, Constraint_System, Generator_System,
    Variable, Linear_Expression,
    line, ray, point )

from base_QQ import Polyhedron_QQ
from representation import (
    PolyhedronRepresentation,
    Hrepresentation,
    Inequality, Equation,
    Vrepresentation,
    Vertex, Ray, Line )



#########################################################################
class Polyhedron_QQ_ppl(Polyhedron_QQ):
    """
    Polyhedra over `\QQ` with ppl

    INPUT:

    - ``ambient_dim`` -- integer. The dimension of the ambient space.

    - ``Vrep`` -- a list ``[vertices, rays, lines]``.

    - ``Hrep`` -- a list ``[ieqs, eqns]``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[], backend='ppl')
        sage: TestSuite(p).run()
    """

    def _init_from_Vrepresentation(self, ambient_dim, vertices, rays, lines, minimize=True):
        """
        Construct polyhedron from V-representation data.

        INPUT:

        - ``ambient_dim`` -- integer. The dimension of the ambient space.

        - ``vertices`` -- list of point. Each point can be specified
           as any iterable container of
           :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``rays`` -- list of rays. Each ray can be specified as any
          iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``lines`` -- list of lines. Each line can be specified as
          any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        EXAMPLES::

            sage: p = Polyhedron(backend='ppl')
            sage: from sage.geometry.polyhedron.backend_ppl import Polyhedron_QQ_ppl
            sage: Polyhedron_QQ_ppl._init_from_Vrepresentation(p, 2, [], [], [])
        """
        gs = Generator_System()
        if vertices is None: vertices = []
        for v in vertices:
            d = lcm([denominator(v_i) for v_i in v])
            dv = [ ZZ(d*v_i) for v_i in v ]
            gs.insert(point(Linear_Expression(dv, 0), d))
        if rays is None: rays = []
        for r in rays:
            d = lcm([denominator(r_i) for r_i in r])
            dr = [ ZZ(d*r_i) for r_i in r ]
            gs.insert(ray(Linear_Expression(dr, 0)))
        if lines is None: lines = []
        for l in lines:
            d = lcm([denominator(l_i) for l_i in l])
            dl = [ ZZ(d*l_i) for l_i in l ]
            gs.insert(line(Linear_Expression(dl, 0)))
        self._ppl_polyhedron = C_Polyhedron(gs)
        self._init_Vrepresentation_from_ppl(minimize)
        self._init_Hrepresentation_from_ppl(minimize)


    def _init_from_Hrepresentation(self, ambient_dim, ieqs, eqns, minimize=True):
        """
        Construct polyhedron from H-representation data.

        INPUT:

        - ``ambient_dim`` -- integer. The dimension of the ambient space.

        - ``ieqs`` -- list of inequalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``eqns`` -- list of equalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        EXAMPLES::

            sage: p = Polyhedron(backend='ppl')
            sage: from sage.geometry.polyhedron.backend_ppl import Polyhedron_QQ_ppl
            sage: Polyhedron_QQ_ppl._init_from_Hrepresentation(p, 2, [], [])
        """
        cs = Constraint_System()
        if ieqs is None: ieqs = []
        for ieq in ieqs:
            d = lcm([denominator(ieq_i) for ieq_i in ieq])
            dieq = [ ZZ(d*ieq_i) for ieq_i in ieq ]
            b = dieq[0]
            A = dieq[1:]
            cs.insert(Linear_Expression(A, b) >= 0)
        if eqns is None: eqns = []
        for eqn in eqns:
            d = lcm([denominator(eqn_i) for eqn_i in eqn])
            deqn = [ ZZ(d*eqn_i) for eqn_i in eqn ]
            b = deqn[0]
            A = deqn[1:]
            cs.insert(Linear_Expression(A, b) == 0)
        self._ppl_polyhedron = C_Polyhedron(cs)
        self._init_Vrepresentation_from_ppl(minimize)
        self._init_Hrepresentation_from_ppl(minimize)


    def _init_Vrepresentation_from_ppl(self, minimize):
        """
        Create the Vrepresentation objects from the ppl polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],
            ...                  backend='ppl')  # indirect doctest
            sage: p.Hrepresentation()
            (An inequality (1, 4) x - 2 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (-5, 12) x + 10 >= 0)
            sage: p._ppl_polyhedron.minimized_constraints()
            Constraint_System {x0+4*x1-2>=0, x0-12*x1+6>=0, -5*x0+12*x1+10>=0}
            sage: p.Vrepresentation()
            (A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6))
            sage: p._ppl_polyhedron.minimized_generators()
            Generator_System {point(0/2, 1/2), point(2/1, 0/1), point(24/6, 5/6)}
        """
        self._Vrepresentation = []
        gs = self._ppl_polyhedron.minimized_generators()
        for g in gs:
            if g.is_point():
                d = g.divisor()
                Vertex(self, [x/d for x in g.coefficients()])
            elif g.is_ray():
                Ray(self, g.coefficients())
            elif g.is_line():
                Line(self, g.coefficients())
            else:
                assert False
        self._Vrepresentation = tuple(self._Vrepresentation)


    def _init_Hrepresentation_from_ppl(self, minimize):
        """
        Create the Vrepresentation objects from the ppl polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],
            ...                  backend='ppl')  # indirect doctest
            sage: p.Hrepresentation()
            (An inequality (1, 4) x - 2 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (-5, 12) x + 10 >= 0)
            sage: p._ppl_polyhedron.minimized_constraints()
            Constraint_System {x0+4*x1-2>=0, x0-12*x1+6>=0, -5*x0+12*x1+10>=0}
            sage: p.Vrepresentation()
            (A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6))
            sage: p._ppl_polyhedron.minimized_generators()
            Generator_System {point(0/2, 1/2), point(2/1, 0/1), point(24/6, 5/6)}
        """
        self._Hrepresentation = []
        cs = self._ppl_polyhedron.minimized_constraints()
        for c in cs:
            if c.is_inequality():
                Inequality(self, (c.inhomogeneous_term(),) + c.coefficients())
            elif c.is_equality():
                Equation(self, (c.inhomogeneous_term(),) + c.coefficients())
        self._Hrepresentation = tuple(self._Hrepresentation)


    def _init_empty_polyhedron(self, ambient_dim):
        """
        Initializes an empty polyhedron.

        INPUT:

        - ``ambient_dim`` -- integer. The dimension of the ambient space.

        TESTS::

            sage: empty = Polyhedron(backend='ppl'); empty
            The empty polyhedron in QQ^0
            sage: empty.Vrepresentation()
            ()
            sage: empty.Hrepresentation()
            (An equation -1 == 0,)
            sage: Polyhedron(vertices = [], backend='ppl')
            The empty polyhedron in QQ^0
            sage: Polyhedron(backend='ppl')._init_empty_polyhedron(0)
        """
        super(Polyhedron_QQ_ppl, self)._init_empty_polyhedron(ambient_dim)
        self._ppl_polyhedron = C_Polyhedron(ambient_dim, 'empty')


