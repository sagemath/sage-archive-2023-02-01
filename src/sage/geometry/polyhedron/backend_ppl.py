"""
The PPL (Parma Polyhedra Library) backend for polyhedral computations
"""

from sage.structure.element import Element
from sage.rings.all import ZZ
from sage.rings.integer import Integer
from sage.arith.functions import LCM_list
from sage.misc.functional import denominator
from .base import Polyhedron_base
from .base_QQ import Polyhedron_QQ
from .base_ZZ import Polyhedron_ZZ

from sage.misc.lazy_import import lazy_import
from sage.features import PythonModule
lazy_import('ppl', ['C_Polyhedron', 'Generator_System', 'Constraint_System',
                    'Linear_Expression', 'line', 'ray', 'point'],
                    feature=PythonModule("ppl", spkg="pplpy"))


#########################################################################
class Polyhedron_ppl(Polyhedron_base):
    """
    Polyhedra with ppl

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[], backend='ppl')
        sage: TestSuite(p).run()
    """

    def __init__(self, parent, Vrep, Hrep, ppl_polyhedron=None, **kwds):
        """
        Initializes the polyhedron.

        See :class:`Polyhedron_ppl` for a description of the input
        data.

        TESTS::

            sage: p = Polyhedron()
            sage: TestSuite(p).run()
            sage: p = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)])
            sage: TestSuite(p).run()
            sage: q = polytopes.cube()
            sage: p = q.parent().element_class(q.parent(), None, None, q._ppl_polyhedron)
            sage: TestSuite(p).run()
        """
        if ppl_polyhedron:
            if Hrep is not None or Vrep is not None:
                raise ValueError("only one of Vrep, Hrep, or ppl_polyhedron can be different from None")
            Element.__init__(self, parent=parent)
            minimize = True if 'minimize' in kwds and kwds['minimize'] else False
            self._init_from_ppl_polyhedron(ppl_polyhedron, minimize)
        else:
            Polyhedron_base.__init__(self, parent, Vrep, Hrep, **kwds)

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

            sage: p = Polyhedron(backend='ppl')
            sage: from sage.geometry.polyhedron.backend_ppl import Polyhedron_ppl
            sage: Polyhedron_ppl._init_from_Vrepresentation(p, [], [], [])
        """
        gs = Generator_System()
        if vertices is None: vertices = []
        for v in vertices:
            gs.insert(self._convert_generator_to_ppl(v, 2))
        if rays is None: rays = []
        for r in rays:
            gs.insert(self._convert_generator_to_ppl(r, 3))
        if lines is None: lines = []
        for l in lines:
            gs.insert(self._convert_generator_to_ppl(l, 4))
        if gs.empty():
            ppl_polyhedron = C_Polyhedron(self.ambient_dim(), 'empty')
        else:
            ppl_polyhedron = C_Polyhedron(gs)
        self._init_from_ppl_polyhedron(ppl_polyhedron, minimize)

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

            sage: p = Polyhedron(backend='ppl')
            sage: from sage.geometry.polyhedron.backend_ppl import Polyhedron_ppl
            sage: Polyhedron_ppl._init_from_Hrepresentation(p, [], [])
        """
        cs = Constraint_System()
        if ieqs is None: ieqs = []
        for ieq in ieqs:
            cs.insert(self._convert_constraint_to_ppl(ieq, 0))
        if eqns is None: eqns = []
        for eqn in eqns:
            cs.insert(self._convert_constraint_to_ppl(eqn, 1))
        if cs.empty():
            ppl_polyhedron = C_Polyhedron(self.ambient_dim(), 'universe')
        else:
            ppl_polyhedron = C_Polyhedron(cs)
        self._init_from_ppl_polyhedron(ppl_polyhedron, minimize)

    def _init_from_ppl_polyhedron(self, ppl_polyhedron, minimize=True):
        """
        Create the V-/Hrepresentation objects from the ppl polyhedron.

        TESTS::

            sage: p = Polyhedron(backend='ppl')
            sage: from sage.geometry.polyhedron.backend_ppl import Polyhedron_ppl
            sage: Polyhedron_ppl._init_from_Hrepresentation(p, [], [])  # indirect doctest
        """
        self._ppl_polyhedron = ppl_polyhedron
        self._init_Vrepresentation_from_ppl(minimize)
        self._init_Hrepresentation_from_ppl(minimize)

    def _init_Vrepresentation_from_ppl(self, minimize):
        """
        Create the Vrepresentation objects from the ppl polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],
            ....:                backend='ppl')  # indirect doctest
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
        parent = self.parent()
        for g in gs:
            coefficients = [Integer(mpz) for mpz in g.coefficients()]
            if g.is_point():
                d = Integer(g.divisor())
                if d.is_one():
                    parent._make_Vertex(self, coefficients)
                else:
                    parent._make_Vertex(self, [x/d for x in coefficients])
            elif g.is_ray():
                parent._make_Ray(self, coefficients)
            elif g.is_line():
                parent._make_Line(self, coefficients)
            else:
                assert False
        self._Vrepresentation = tuple(self._Vrepresentation)

    def _init_Hrepresentation_from_ppl(self, minimize):
        """
        Create the Hrepresentation objects from the ppl polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],
            ....:                backend='ppl')  # indirect doctest
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
        parent = self.parent()
        for c in cs:
            if c.is_inequality():
                parent._make_Inequality(self, (c.inhomogeneous_term(),) + c.coefficients())
            elif c.is_equality():
                parent._make_Equation(self, (c.inhomogeneous_term(),) + c.coefficients())
        self._Hrepresentation = tuple(self._Hrepresentation)

    def _init_empty_polyhedron(self):
        """
        Initializes an empty polyhedron.

        TESTS::

            sage: empty = Polyhedron(backend='ppl'); empty
            The empty polyhedron in ZZ^0
            sage: empty.Vrepresentation()
            ()
            sage: empty.Hrepresentation()
            (An equation -1 == 0,)
            sage: Polyhedron(vertices = [], backend='ppl')
            The empty polyhedron in ZZ^0
            sage: Polyhedron(backend='ppl')._init_empty_polyhedron()
        """
        super(Polyhedron_ppl, self)._init_empty_polyhedron()
        self._ppl_polyhedron = C_Polyhedron(self.ambient_dim(), 'empty')

    @staticmethod
    def _convert_generator_to_ppl(v, typ):
        r"""
        Convert a generator to ``ppl``.

        INPUT:

        - ``v`` -- a vertex, ray, or line.

        - ``typ`` -- integer; 2 -- vertex; 3 -- ray; 4 -- line

        EXAMPLES::

            sage: P = Polyhedron()
            sage: P._convert_generator_to_ppl([1, 1/2, 3], 2)
            point(2/2, 1/2, 6/2)
            sage: P._convert_generator_to_ppl([1, 1/2, 3], 3)
            ray(2, 1, 6)
            sage: P._convert_generator_to_ppl([1, 1/2, 3], 4)
            line(2, 1, 6)
        """
        if typ == 2:
            ob = point
        elif typ == 3:
            ob = ray
        else:
            ob = line

        d = LCM_list([denominator(v_i) for v_i in v])
        if d.is_one():
            return ob(Linear_Expression(v, 0))
        else:
            dv = [ d*v_i for v_i in v ]
            if typ == 2:
                return ob(Linear_Expression(dv, 0), d)
            else:
                return ob(Linear_Expression(dv, 0))

    @staticmethod
    def _convert_constraint_to_ppl(c, typ):
        r"""
        Convert a constraint to ``ppl``.

        INPUT:

        - ``c`` -- an inequality or equation.

        - ``typ`` -- integer; 0 -- inequality; 3 -- equation

        EXAMPLES::

            sage: P = Polyhedron()
            sage: P._convert_constraint_to_ppl([1, 1/2, 3], 0)
            x0+6*x1+2>=0
            sage: P._convert_constraint_to_ppl([1, 1/2, 3], 1)
            x0+6*x1+2==0
        """
        d = LCM_list([denominator(c_i) for c_i in c])
        dc = [ ZZ(d*c_i) for c_i in c ]
        b = dc[0]
        A = dc[1:]
        if typ == 0:
            return Linear_Expression(A, b) >= 0
        else:
            return Linear_Expression(A, b) == 0




#########################################################################
class Polyhedron_QQ_ppl(Polyhedron_ppl, Polyhedron_QQ):
    r"""
    Polyhedra over `\QQ` with ppl

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[],
        ....:                backend='ppl', base_ring=QQ)
        sage: TestSuite(p).run(skip='_test_pickling')
    """
    pass


#########################################################################
class Polyhedron_ZZ_ppl(Polyhedron_ppl, Polyhedron_ZZ):
    r"""
    Polyhedra over `\ZZ` with ppl

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)], lines=[],
        ....:                backend='ppl', base_ring=ZZ)
        sage: TestSuite(p).run(skip='_test_pickling')
    """
    pass
