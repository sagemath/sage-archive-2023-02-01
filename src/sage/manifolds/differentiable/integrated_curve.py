r"""
Integrated Curves in Manifolds

Given a differentiable manifold `M`, an *integrated curve* curve in `M`
is a differentiable curve constructed as a numerical solution to a
system of second order differential equations.

Integrated curves are implemented by :class:`IntegratedCurve`, which the
classes :class:`IntegratedAutoparallelCurve` and
:class:`IntegratedGeodesic` inherit.

AUTHORS:

- Karim Van Aelst (2017): initial version

"""

#***********************************************************************
#       Copyright (C) 2017 Karim Van Aelst <karim.van-aelst@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#***********************************************************************

from sage.symbolic.expression import Expression
from sage.rings.infinity import Infinity
from sage.calculus.desolvers import desolve_system_rk4
from sage.calculus.desolvers import desolve_odeint
from sage.manifolds.chart import Chart
from sage.manifolds.differentiable.curve import DifferentiableCurve
from sage.manifolds.differentiable.tangent_vector import TangentVector
from sage.calculus.interpolation import Spline
from sage.misc.decorators import options
from sage.misc.functional import numerical_approx

class IntegratedCurve(DifferentiableCurve):
    r"""
    Given a chart with coordinates denoted :MATH:`(x_{1}, ..., x_{n})`,
    an instance of :class:`IntegratedCurve` is a curve
    :MATH:`t \mapsto (x_{1}(t), ..., x_{n}(t))` constructed as a
    numerical solution to a system of second order differential
    equations satisfied by the coordinate curves
    :MATH:`t \mapsto x_{i}(t)`.

    INPUT:

    - ``parent`` --
      :class:`~sage.manifolds.differentiable.manifold_homset.IntegratedCurveSet`
      the set of curves `\mathrm{Hom}(I, M)` to which the curve belongs
    - ``equations_rhs`` -- list of the right-hand sides of the equations
      on the velocities only (the term *velocity* referring to the
      derivatives :MATH:`d x_{i} / dt` of the coordinate curves)
    - ``velocities`` -- list of the symbolic expressions used in
      ``equations_rhs`` to denote the velocities
    - ``curve_parameter`` -- symbolic expression used in
      ``equations_rhs`` to denote the parameter of the curve (denoted
      :MATH:`t` in the descriptions above)
    - ``initial_tangent_vector`` --
      :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
      initial tangent vector of the curve
    - ``chart`` -- (default: ``None``) chart on the manifold in
      which the equations are given; if ``None`` the default chart
      of the manifold is assumed
    - ``name`` -- (default: ``None``) string; symbol given to the curve
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the curve; if none is provided, ``name`` will be used
    - ``is_isomorphism`` -- (default: ``False``) determines whether the
      constructed object is a diffeomorphism; if set to ``True``,
      then `M` must have dimension one
    - ``is_identity`` -- (default: ``False``) determines whether the
      constructed object is the identity map; if set to ``True``,
      then `M` must coincide with the domain of the curve

    EXAMPLE:

    Motion of a charged particle in an axial magnetic field linearly
    increasing in time and exponentially decreasing in space:

    .. MATH::

        \mathbf{B}(t,\mathbf{x}) = \frac{B_{0}t}{T} \exp \left(
        -\frac{ x_{1}^{2} + x_{2}^{2} }{ L^{2} } \right) \mathbf{e_{3}}.

    Equations of motion are:

    .. MATH::

        \ddot{x}_{1}(t) &= \frac{qB(t,\mathbf{x}(t))}{m}
                           \dot{x}_{2}(t)                             \\
        \ddot{x}_{2}(t) &= - \frac{qB(t, \mathbf{x}(t))}{m}
                           \dot{x}_{1}(t)                             \\
        \ddot{x}_{3}(t) &= 0

    Start with declaring a chart on a 3-dimensional manifold and the
    symbolic expressions denoting the velocities and the various
    parameters::

        sage: M = Manifold(3, 'M', start_index=1)
        sage: X.<x1,x2,x3> = M.chart()
        sage: var('t B_0 m q L T')
        (t, B_0, m, q, L, T)
        sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
        sage: D = X.symbolic_velocities(); D
        [Dx1, Dx2, Dx3]
        sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]

    Set the initial conditions::

        sage: p = M.point((0,0,0), name='p')
        sage: Tp = M.tangent_space(p)
        sage: v = Tp((1,0,1))

    Declare an integrated curve and display information relative to it::

        sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
        ....:                                         verbose=False) ; c
        Integrated curve c in the 3-dimensional differentiable
         manifold M
        sage: sys = c.system()
        Curve c in the 3-dimensional differentiable manifold M
         integrated over the Real interval (0, 5) as a solution to the
         following system, written w.r.t.
         Chart (M, (x1, x2, x3)):
        <BLANKLINE>
        Initial point: Point p on the 3-dimensional differentiable
         manifold M with coordinates [0, 0, 0] w.r.t.
         Chart (M, (x1, x2, x3))
        Initial tangent vector: Tangent vector at Point p on
         the 3-dimensional differentiable manifold M with
         components [1, 0, 1] w.r.t. Chart (M, (x1, x2, x3))
        <BLANKLINE>
        d(x1)/dt = Dx1
        d(x2)/dt = Dx2
        d(x3)/dt = Dx3
        d(Dx1)/dt = B_0*Dx2*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
        d(Dx2)/dt = -B_0*Dx1*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
        d(Dx3)/dt = 0
        <BLANKLINE>

    Generate a solution of the system and an interpolation of this
    solution::

        sage: sol = c.solve(step=0.2,
        ....:         parameters_values={B_0:1, m:1, q:1, L:10, T:1},
        ....:         solution_key='carac time 1')
        Performing 4th order Runge-Kutta integration with Maxima by
         default...
        <BLANKLINE>
        Numerical integration completed.
        Checking all points are in the chart domain...
        <BLANKLINE>
        All points are in the chart domain.
        The resulting list of points was associated with the key
         'carac time 1' (if this key already referred to a former
         numerical solution, such a solution was erased).
        sage: interp = c.interpolate(solution_key='carac time 1',
        ....:                              interpolation_key='interp 1')
        Performing cubic spline interpolation by default...
        Interpolation completed and associated with the key 'interp 1'
         (if this key already referred to a former interpolation,
         such an interpolation was erased).

    Such an interpolation is required to evaluate the curve and the
    vector tangent to the curve for any value of the curve parameter::

        sage: c(1.9)
        Evaluating point coordinates from the interpolation associated
         with the key 'interp 1' by default...
        [1.3776707219621374, -0.9000776970132945, 1.9]
        sage: v = c.tangent_vector_eval_at(4.3)
        Evaluating tangent vector components from the interpolation
         associated with the key 'interp 1' by default...
        sage: v
        Tangent vector at Point on the 3-dimensional differentiable
         manifold M
        sage: v.display()
        -0.9303968397216424 d/dx1 - 0.3408080563014475 d/dx2 +
         1.0000000000000004 d/dx3

    Plotting a numerical solution (with or without its tangent vector
    field) also requires the solution to be interpolated at least once::

        sage: c_plot_2d_1 = c.plot_integrated(ambient_coords=[x1, x2],
        ....:               interpolation_key='interp 1', thickness=2.5,
        ....:               display_tangent=True, plot_points=200,
        ....:               plot_points_tangent=10, scale=0.5,
        ....:               color='blue', color_tangent='red')
        A tiny final offset equal to the value of 'end_point_offset[1]'
         (= 0.001) was introduced in order to safely compute the last
         tangent vector from the interpolation.
        sage: c_plot_2d_1.show()

    .. PLOT::

        M = Manifold(3, 'M')
        X = M.chart('x1 x2 x3')
        var('x1 x2 x3 t B_0 m q L T')
        B = B_0*t/T*exp(-(x1**2 + x2**2)/L**2)
        D = X.symbolic_velocities()
        eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
        p = M.point((0,0,0), name='p')
        Tp = M.tangent_space(p)
        v = Tp((1,0,1))
        c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c')
        c.solve(step=0.2,
                parameters_values={B_0:1, m:1, q:1, L:10, T:1},
                solution_key='carac time 1')
        c.interpolate(solution_key='carac time 1',
                      interpolation_key='interp 1')
        c_plot_2d_1 = c.plot_integrated(ambient_coords=[x1, x2],
                        interpolation_key='interp 1', thickness=2.5,
                        display_tangent=True, plot_points=200,
                        plot_points_tangent=10, scale=0.5, color='blue',
                        color_tangent='red')
        sphinx_plot(c_plot_2d_1)

    An instance of :class:`IntegratedCurve` may store several numerical
    solutions and interpolations::

        sage: sol = c.solve(step=0.2,
        ....:         parameters_values={B_0:1, m:1, q:1, L:10, T:100},
        ....:         solution_key='carac time 100', verbose=False)
        sage: interp = c.interpolate(solution_key='carac time 100',
        ....:             interpolation_key='interp 100', verbose=False)
        sage: c_plot_3d_100=c.plot_integrated(interpolation_key='interp 100',
        ....:                   thickness=2.5, display_tangent=True,
        ....:                   plot_points=200, plot_points_tangent=10,
        ....:                   scale=0.5, color='green',
        ....:                   color_tangent='orange', verbose=False)
        sage: c_plot_3d_1=c.plot_integrated(interpolation_key='interp 1',
        ....:                   thickness=2.5, display_tangent=True,
        ....:                   plot_points=200, plot_points_tangent=10,
        ....:                   scale=0.5, color='blue',
        ....:                   color_tangent='red', verbose=False)
        sage: graph = c_plot_3d_1 + c_plot_3d_100
        sage: graph.show()

    .. PLOT::

        M = Manifold(3, 'M')
        X = M.chart('x1 x2 x3')
        var('x1 x2 x3 t B_0 m q L T')
        B = B_0*t/T*exp(-(x1**2 + x2**2)/L**2)
        D = X.symbolic_velocities()
        eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
        p = M.point((0,0,0), name='p')
        Tp = M.tangent_space(p)
        v = Tp((1,0,1))
        c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c')
        c.solve(step=0.2,
                parameters_values={B_0:1, m:1, q:1, L:10, T:1},
                solution_key='carac time 1')
        c.interpolate(solution_key='carac time 1',
                      interpolation_key='interp 1')
        c.solve(step=0.2,
                parameters_values={B_0:1, m:1, q:1, L:10, T:100},
                solution_key='carac time 100')
        c.interpolate(solution_key='carac time 100',
                      interpolation_key='interp 100')
        c_plot_3d_1 = c.plot_integrated(interpolation_key='interp 1',
                        thickness=2.5, display_tangent=True,
                        plot_points=200, plot_points_tangent=10,
                        scale=0.5, color='blue', color_tangent='red',
                        verbose=False)
        c_plot_3d_100 = c.plot_integrated(interpolation_key='interp 100',
                            thickness=2.5, display_tangent=True,
                            plot_points=200, plot_points_tangent=10,
                            scale=0.5, color='green',
                            color_tangent='orange', verbose=False)
        graph = c_plot_3d_1 + c_plot_3d_100
        sphinx_plot(graph)

    """

    def __init__(self, parent, equations_rhs, velocities,
                 curve_parameter, initial_tangent_vector, chart=None,
                 name=None, latex_name=None,
                 is_isomorphism=False, is_identity=False, verbose=True):
        r"""
        Construct a numerical curve.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns + [x1], D, (t, 0, 5), v,
            ....:                                              name='c')
            Traceback (most recent call last):
            ...
            ValueError: Number of equations should equal codomain
             dimension.
            sage: c = M.integrated_curve(eqns, D + [x1], (t, 0, 5), v,
            ....:                                              name='c')
            Traceback (most recent call last):
            ...
            ValueError: Number of velocities should equal codomain
             dimension.
            sage: c = M.integrated_curve(eqns, D,(t,-oo,5), v, name='c')
            Traceback (most recent call last):
            ...
            ValueError: Both boundaries of the interval defining the
             domain of a Homset of integrated curves need to be finite.
            sage: c = M.integrated_curve(eqns, D, (t,0,5), x1, name='c')
            Traceback (most recent call last):
            ...
            TypeError: x1 should be a tangent vector.
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c',
            ....:                                     verbose=False) ; c
            Integrated curve c in the 3-dimensional differentiable
             manifold M
            sage: #TestSuite(c).run()

        """

        from sage.symbolic.ring import SR

        # start with parent class method to initialize the four last
        # arguments:
        DifferentiableCurve.__init__(self, parent, name=name,
                   latex_name=latex_name, is_isomorphism=is_isomorphism,
                   is_identity=is_identity) # (coord_expression=None)

        # check argument 'parent': 't_min' and 't_max' below are only
        # allowed to be either expressions of finite real values:
        domain = self.domain()
        t_min = domain.lower_bound()
        t_max = domain.upper_bound()
        if t_min == -Infinity or t_max == +Infinity:
            raise ValueError("Both boundaries of the interval " +
                             "need to be finite.")

        codomain = self.codomain()

        if self._is_identity: # init of DifferentiableCurve should have
        # set this attribute to 'True' if argument is_identity of this
        # __init__ was set to 'True'.
        # In this case, init of DifferentiableCurve has checked that the
        # domain and codomain coincide.
        # Then, at this stage, they necessarily are a certain finite
        # real interval.
        # Arguments 'equations_rhs', 'initial_tangent_vector' and
        # 'chart' are then modified whatever their values were in order
        # to construct an integrated version of the identity map on the
        # finite interval.
            equations_rhs = [0]
            chart = codomain.default_chart()
            p = codomain.point([t_min + (t_max-t_min)/10**(6)])#slightly
            # shifting the initial point inside the interval is required
            # since the interval is open
            Tp = codomain.tangent_space(p)
            initial_tangent_vector = Tp([1 - 10**(-6)]) # taking a
            # derivative slightly less than 1 is required for the last
            # point to be inside the open interval as well (this
            # obviously damages the approximation of the identity
            # provided by 'self')

        # check argument 'equations_rhs':
        dim = codomain.dim()
        if len(equations_rhs) != dim:
            raise ValueError("Number of equations should equal " +
                             "codomain dimension.")

        # check the chart:
        if chart is not None:
            if chart not in codomain.atlas():
                raise ValueError("{} should be a chart ".format(chart) +
                                 "on the {}".format(codomain))
        else:
            chart = codomain.default_chart()

        # check argument 'velocities':
        if len(velocities) != dim:
            raise ValueError("Number of velocities should equal " +
                             "codomain dimension.")
        # in particular, check that no velocity coincides with a
        # coordinate:
        for vel in velocities:
            if vel in chart[:]:
                str_error = "{} should not be used as a ".format(vel)
                str_error += "velocity since it also denotes "
                str_error += "a coordinate."
                raise ValueError(str_error)

        # check argument 'curve_parameter':
        if not isinstance(curve_parameter, Expression):
            raise TypeError("{} should be ".format(curve_parameter) +
                             "a symbolic expression.")
        # in particular, check that it does not coincide with a
        # coordinate or a velocity:
        coords_vels = list(chart[:]) + list(velocities)
        if curve_parameter in coords_vels:
            str_error = "{} should not be used ".format(curve_parameter)
            str_error += "as the curve parameter since it also denotes "
            str_error += "a coordinate or a velocity."
            raise ValueError(str_error)
        # the various algorithms called in 'solve' method are in charge
        # of raising errors about possibly remaining problems regarding
        # 'curve_parameter'

        # check argument 'initial_tangent_vector':
        if not isinstance(initial_tangent_vector, TangentVector):
            raise TypeError("{} ".format(initial_tangent_vector) +
                            "should be a tangent vector.")
        # in particular, check that its base point sits in the domain
        # of the chart (if its coordinates are explicitly given):
        initial_pt = initial_tangent_vector.parent().base_point()
        # line above retrieves the initial point as the base point of
        # the tangent space to which the initial tangent vector belongs
        initial_pt_coords = initial_pt.coordinates(chart)
        i0 = chart.manifold().start_index()
        for i in range(dim):
            if not isinstance(initial_pt_coords[i], Expression):
                coord_value = initial_pt_coords[i]
                coord_min = chart.coord_bounds(i+i0)[0][0]
                coord_max = chart.coord_bounds(i+i0)[1][0]
                if coord_value <= coord_min or coord_value >= coord_max:
                    raise ValueError("Initial point should be in the " +
                                     "domain of the chart.")

        # prepare attribute '_parameters':
        announced_variables = set(coords_vels + [curve_parameter])
        parameters = set()
        # extract all the variables appearing in the equations:
        for eqn in equations_rhs:
            if isinstance(eqn, Expression): # some right hand sides
            # might merely be real numbers and not expressions, so that
            # they do not contain any variable, and method 'variables'
            # could not be called on them
                parameters = parameters.union(eqn.variables())
        # remove the Expressions that should not be treated as
        # parameters (i.e. the coordinate functions, the velocities and
        # the curve parameter):
        parameters = parameters.difference(announced_variables)
        # extract all the variables appearing in the boundaries:
        if isinstance(t_min, Expression):
            parameters = parameters.union(t_min.variables())
        if isinstance(t_max, Expression):
            parameters = parameters.union(t_max.variables())
        # extract all the variables appearing in the initial point
        # coordinates:
        for coord in initial_pt_coords:
            if isinstance(coord,Expression):
                parameters = parameters.union(coord.variables())
        # extract all the variables appearing in the initial tangent
        # vector components:
        initial_coord_basis = chart.frame().at(initial_pt)
        initial_tgt_vec_comps=initial_tangent_vector[initial_coord_basis,:]
        for comp in initial_tgt_vec_comps:
            if isinstance(comp, Expression):
                parameters = parameters.union(comp.variables())

        # check at this stage that no parameter coincides with a
        # coordinate, a velocity, or the curve parameter; this would
        # mean that an Expression used to denote either a bound, a
        # coordinate of the initial point or a component of the initial
        # tangent vector coincides with a coordinate, a velocity or the
        # curve parameter (which would make no sense):
        if len(parameters) != 0:
            for param in parameters:
                if param in announced_variables:
                    str_error = "{} should not be used ".format(param)
                    str_error += "as a parameter since it also denotes "
                    str_error += "a coordinate, a velocity or the "
                    str_error += "curve parameter."
                    raise ValueError(str_error)

        # define all attributes
        self._equations_rhs = list(equations_rhs) # converts to list
        # since might not already be a list (which is later required)
        self._velocities = list(velocities) # converts to list
        # since might not already be a list (which is later required)
        self._curve_parameter = curve_parameter
        self._initial_tangent_vector = initial_tangent_vector
        self._chart = chart
        self._parameters = parameters
        self._ode_solver = None # if needed, becomes an instance of
        # 'ode_solver', which performs most of the numerical integrations
        # offered by method 'solve'
        self._solutions = {} # dictionary containing all numerically
        # computed lists of points of the curve, the keys being chosen
        # by the user when calling method 'solve'
        self._interpolations = {} # dictionary containing lists of
        # interpolation objects, each interpolation object implementing
        # the interpolation of one of the numerical coordinate curves,
        # and the keys being chosen by the user when calling
        # method 'interpolate'

        if verbose:
            print("The curve was correctly set.")
            if len(self._parameters) != 0:
                print("Parameters appearing in the differential " +
                      "system defining the curve are " +
                      "{}.".format(self._parameters))
            else:
                print("No parameter appears in the differential " +
                      "system defining the curve.")

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                                     verbose=False) ; c
            Integrated curve in the 3-dimensional differentiable
             manifold M
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v,
            ....:                           name='c', verbose=False) ; c
            Integrated curve c in the 3-dimensional differentiable
             manifold M

        """

        description = "Integrated curve "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description

    def __reduce__(self):
        r"""
        Reduction function for the pickle protocole.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c',
            ....:                                         verbose=False)
            sage: c.__reduce__()
            (<class 'sage.manifolds.differentiable.integrated_curve.IntegratedCurveSet_with_category.element_class'>,
             (Set of Morphisms from Real interval (0, 5) to
              3-dimensional differentiable manifold M in Category of
              homsets of subobjects of sets and topological spaces which
              actually are integrated curves,
              [B_0*Dx2*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m),
               -B_0*Dx1*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m),
               0],
              [Dx1, Dx2, Dx3],
              t,
              Tangent vector at Point p on the 3-dimensional differentiable manifold M,
              Chart (M, (x1, x2, x3)),
              'c',
              'c',
              False,
              False))

        Test of pickling::

            sage: loads(dumps(c))
            Integrated curve c in the 3-dimensional differentiable
             manifold M

        """

        return (type(self), (self.parent(), self._equations_rhs,
                self._velocities, self._curve_parameter,
                self._initial_tangent_vector, self._chart,
                self._name, self._latex_name, self._is_isomorphism,
                self._is_identity))

    def system(self, verbose=True):
        r"""
        Provide a detailed description of the system defining the curve
        and returns the system defining it: chart, equations and initial
        conditions.

        INPUT:

        - ``verbose`` -- (default: ``True``) prints a detailed
          description of the curve

        OUTPUT:

        - list containing the attributes :attr:`equations_rhs`,
          :attr:`initial_tangent_vector` and :attr:`chart`

        EXAMPLE:

        System defining an integrated curve::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c',
            ....:                                         verbose=False)
            sage: sys = c.system()
            Curve c in the 3-dimensional differentiable manifold M
             integrated over the Real interval (0, 5) as a solution to
             the following system, written w.r.t.
             Chart (M, (x1, x2, x3)):
            <BLANKLINE>
            Initial point: Point p on the 3-dimensional differentiable
             manifold M with coordinates [0, 0, 0] w.r.t.
             Chart (M, (x1, x2, x3))
            Initial tangent vector: Tangent vector at Point p on the
             3-dimensional differentiable manifold M with
             components [1, 0, 1] w.r.t. Chart (M, (x1, x2, x3))
            <BLANKLINE>
            d(x1)/dt = Dx1
            d(x2)/dt = Dx2
            d(x3)/dt = Dx3
            d(Dx1)/dt = B_0*Dx2*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
            d(Dx2)/dt = -B_0*Dx1*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
            d(Dx3)/dt = 0
            <BLANKLINE>
            sage: sys_mute = c.system(verbose=False)
            sage: sys_mute == sys
            True

        """

        v0 = self._initial_tangent_vector
        chart = self._chart

        if verbose:
            initial_tgt_space = v0.parent()
            initial_pt = initial_tgt_space.base_point() # retrieves
            # the initial point as the base point of the tangent space
            # to which initial tangent vector belongs
            initial_pt_coords = list(initial_pt.coordinates(chart))
            # previous line converts to list since would otherwise be a
            # tuple ; will raise error if coordinates in chart are not
            # known

            initial_coord_basis = chart.frame().at(initial_pt)
            initial_tgt_vec_comps = v0[initial_coord_basis,:] # will
            # raise error if components in coordinate basis are not
            # known

            description = "Curve "
            if self._name is not None:
                description += self._name + " "
            description += "in the {} ".format(self.codomain())
            description += "integrated over the "
            description += "{} ".format(self.domain())
            description += "as a solution to the following system, "
            description += "written w.r.t. "
            description += "{}:\n\n".format(chart)

            description += "Initial point: {} ".format(initial_pt)
            description += "with coordinates "
            description += "{} ".format(initial_pt_coords)
            description += "w.r.t. {}\n".format(chart)

            description += "Initial tangent vector: {} ".format(v0)
            description += "with components "
            description +="{}".format(initial_tgt_vec_comps)
            description += " w.r.t. {}\n\n".format(chart)

            for coord_func,velocity in zip(chart[:],self._velocities):
                description += "d({})/d{} = {}\n".format(coord_func,
                                                  self._curve_parameter,
                                                  velocity)

            for velocity,eqn in zip(self._velocities,self._equations_rhs):
                description += "d({})/d{} = {}\n".format(velocity,
                                                  self._curve_parameter,
                                                  eqn)

            print(description)

        return [self._equations_rhs, v0, chart]

    def solve_analytical(self, verbose=True):
        r"""
        Solve analytically the differential system defining the curve
        using Maxima via Sage solver ``desolve_system``.
        In case of success, the analytical expressions are added to the
        dictionnary of expressions representing the curve.
        Pay attention to the fact that ``desolve_system`` only considers
        initial conditions given at an initial parameter value equal to
        zero, although the parameter range may not contain zero.
        Yet, assuming that it does, values of the coordinates functions
        at such zero initial parameter value are denoted by the name of
        the coordinate function followed by the string "_0".

        OUTPUT:

        - list of the analytical expressions of the coordinate functions
          (when the differential system could be solved analytically),
          or boolean 'FALSE' (in case the differential system could not
          be solved analytically)

        EXAMPLE:

        Analytical expression of the trajectory of a charged particle in
        a uniform, stationnary magnetic field::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q] = var('t B_0 m q')
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B_0/m*D[1], -q*B_0/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c',
            ....:                                         verbose=False)
            sage: sys = c.system()
            Curve c in the 3-dimensional differentiable manifold M
             integrated over the Real interval (0, 5) as a solution to
             the following system, written w.r.t.
             Chart (M, (x1, x2, x3)):
            <BLANKLINE>
            Initial point: Point p on the 3-dimensional differentiable
             manifold M with coordinates [0, 0, 0] w.r.t.
             Chart (M, (x1, x2, x3))
            Initial tangent vector: Tangent vector at Point p on the
             3-dimensional differentiable manifold M with components
             [1, 0, 1] w.r.t. Chart (M, (x1, x2, x3))
            <BLANKLINE>
            d(x1)/dt = Dx1
            d(x2)/dt = Dx2
            d(x3)/dt = Dx3
            d(Dx1)/dt = B_0*Dx2*q/m
            d(Dx2)/dt = -B_0*Dx1*q/m
            d(Dx3)/dt = 0
            <BLANKLINE>
            sage: sol = c.solve_analytical(verbose=False)
            sage: c.expr()
            ((B_0*q*x1_0 - Dx2_0*m*cos(B_0*q*t/m) +
               Dx1_0*m*sin(B_0*q*t/m) + Dx2_0*m)/(B_0*q),
             (B_0*q*x2_0 + Dx1_0*m*cos(B_0*q*t/m) +
              Dx2_0*m*sin(B_0*q*t/m) - Dx1_0*m)/(B_0*q),
             Dx3_0*t + x3_0)

        """

        from sage.calculus.var import function
        from sage.calculus.functional import diff
        from sage.calculus.desolvers import desolve_system
        from sage.symbolic.assumptions import assume, forget
        from sage.symbolic.ring import var

        dim = self.codomain().dim()
        i0 = self.codomain().start_index()
        des = self._velocities + self._equations_rhs
        par = self._curve_parameter

        if len(self._parameters) != 0:
            for param in self._parameters:
                assume(param != 0)

        y = []
        for i in range(2*dim):
            name = "y{}".format(i+i0)
            y += [function(name)(par)]

        for i in range(dim):
            vel = self._velocities[i]
            des[i] = des[i].substitute({vel:y[dim+i]})
            des[i] = diff(y[i],par) == des[i]
            for j in range(dim):
                coord = self._chart[:][j] # important to use '[:]' on
                # 'chart' to avoid problems due to non zero starting
                # index (i0)
                veloc = self._velocities[j]
                des[dim+i]=des[dim+i].substitute({coord:y[j]})
                des[dim+i]=des[dim+i].substitute({veloc:y[dim+j]})
            des[dim+i] = diff(y[dim+i],par) == des[dim+i]

        dvars = y
        ics = [0]
        y_ics_first_half = []
        y_ics_second_half = []
        for i in range(dim):
            coord = self._chart[:][i] # important to use '[:]'
            # on 'chart' to avoid problems due to non zero
            # starting index (i0)
            veloc = self._velocities[i]
            str_var_coord = "{}_0".format(coord)
            str_var_veloc = "{}_0".format(veloc)
            y_coord_0 = var(str_var_coord)
            y_veloc_0 = var(str_var_veloc)
            y_ics_first_half += [y_coord_0]
            y_ics_second_half += [y_veloc_0]
        ics += y_ics_first_half + y_ics_second_half

        try:
            sol = desolve_system(des, dvars, ivar=self._curve_parameter,
                                                                ics=ics)
        except NotImplementedError:
            coords_sol_expr = False
            if verbose:
                print("The system could not be solved analytically.")
        else:
            coords_sol_expr = []
            for relation in sol[0:dim]:
                expr = relation.rhs().simplify_full()
                coords_sol_expr += [expr]
            self.add_expr(self.domain().default_chart(), self._chart,
                                                        coords_sol_expr)

        if len(self._parameters) != 0:
            for param in self._parameters:
                forget(param != 0)

        return tuple(coords_sol_expr)

    def solve(self, step=None, method=None, solution_key=None,
              parameters_values=None, verbose=True):
        r"""
        Integrate the curve numerically over the domain of integration.

        INPUT:

        - ``step`` -- (default: ``0.1``) step of integration
        - ``method`` -- (default: ``None``) numerical scheme to use for
          the integration of the curve; algorithms available are:

          * 'rk4_maxima' - 4th order classical Runge-Kutta, which makes
            use of Maxima's dynamics package via Sage solver
            ``desolve_system_rk4``
          * 'ode_int' - makes use of ``odeint`` from scipy.integrate
            module via Sage solver ``desolve_odeint``

        and those provided by ``GSL`` via Sage class
        :class:`~sage.calculus.ode.ode_solver`:

          * 'rk2' - embedded Runge-Kutta (2,3)
          * 'rk4' - 4th order classical Runge-Kutta
          * 'rkf45' - Runge-Kutta-Felhberg (4,5)
          * 'rkck' - embedded Runge-Kutta-Cash-Karp (4,5)
          * 'rk8pd' - Runge-Kutta prince-dormand (8,9)
          * 'rk2imp' - implicit 2nd order Runge-Kutta at Gaussian points
          * 'rk4imp' - implicit 4th order Runge-Kutta at Gaussian points
          * 'gear1' - M=1 implicit Gear
          * 'gear2' - M=2 implicit Gear
          * 'bsimp' - implicit Burlisch-Stoer (requires Jacobian)

        - ``solution_key`` -- (default: ``None``) key which the
          resulting numerical solution will be associated to ; a default
          value is given if none is provided
        - ``parameters_values`` -- (default: ``None``) list of numerical
          values of the parameters present in the system defining the
          curve, to be substituted in the equations before integration
        - ``verbose`` -- (default: ``True``) prints information about
          the computation in progress

        OUTPUT:

        - list of the numerical points of the solution computed

        EXAMPLE:

        Computing a numerical solution::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c',
            ....:                                         verbose=False)
            sage: sol = c.solve(parameters_values={m:1, q:1, L:10, T:1})
            Traceback (most recent call last):
            ...
            ValueError: Numerical values should be provided for each of
             the parameters set([B_0, m, q, L, T]).
            sage: sol = c.solve(method='my method',
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            Traceback (most recent call last):
            ...
            ValueError: No available method of integration referred to
             as 'my method'.
            sage: sol = c.solve(
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            Performing 4th order Runge-Kutta integration with Maxima by
             default...
            Resulting list of points will be associated with the key
             'rk4_maxima' by default.
            <BLANKLINE>
            Numerical integration completed.
            Checking all points are in the chart domain...
            <BLANKLINE>
            All points are in the chart domain.
            The resulting list of points was associated with the key
             'rk4_maxima' (if this key already referred to a former
             numerical solution, such a solution was erased).
            sage: sol_mute = c.solve(verbose=False,
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            sage: sol_mute == sol
            True

        """

        from sage.symbolic.ring import SR

        if method is None:
            method = 'rk4_maxima'
            if verbose:
                print("Performing 4th order Runge-Kutta integration " +
                      "with Maxima by default...")

        if solution_key is None:
            solution_key = method
            if verbose:
                print("Resulting list of points will be associated " +
                      "with the key '{}' ".format(solution_key) +
                      "by default.")

        t_min = self.domain().lower_bound()
        t_max = self.domain().upper_bound()

        eqns_rhs = self._equations_rhs

        v0 = self._initial_tangent_vector
        chart = self._chart

        initial_tgt_space = v0.parent()
        initial_pt = initial_tgt_space.base_point() # retrieves
        # the initial point as the base point of the tangent space
        # to which the initial tangent vector belongs
        initial_pt_coords = list(initial_pt.coordinates(chart))
        # previous line converts to list since would otherwise be a
        # tuple (yet might need to be added to [t_min] later); will
        # raise error if coordinates in chart cannot be obtained

        initial_coord_basis = chart.frame().at(initial_pt)
        initial_tgt_vec_comps = list(v0[initial_coord_basis,:])#idem

        dim = self.codomain().dim()

        if len(self._parameters) != 0:
            if parameters_values is None or len(parameters_values)!=len(self._parameters):
                raise ValueError("Numerical values should be " +
                                 "provided for each of the " +
                                 "parameters "
                                 "{}.".format(self._parameters))
            for key in parameters_values.keys():
                parameters_values[key]=numerical_approx(parameters_values[key])#gets
                # numerical values in case some parameters values
                # contain expressions such as pi; will raise error if
                # any element of parameters_values is not numerical

            if isinstance(t_min, Expression):
                t_min = parameters_values[t_min]
                if t_min == -Infinity or t_min == +Infinity:
                    raise ValueError("Both boundaries of the " +
                                      "interval need to be finite.")

            if isinstance(t_max, Expression):
                t_max = parameters_values[t_max]
                if t_max == -Infinity or t_max == +Infinity:
                    raise ValueError("Both boundaries of the " +
                                     "interval need to be finite.")

            for i in range(dim):
                if isinstance(eqns_rhs[i], Expression): # some right
                # hand sides might merely be real numbers and not
                # expressions, so that they do not contain any variable,
                # and method 'variables' could not be called on them
                    eqns_rhs[i]=eqns_rhs[i].substitute(parameters_values)

            for i in range(dim):
                if isinstance(initial_pt_coords[i],Expression):
                    AUX = initial_pt_coords[i]
                    AUX = AUX.substitute(parameters_values)
                    initial_pt_coords[i] = AUX
                if isinstance(initial_tgt_vec_comps[i],Expression):
                    AUX2 = initial_tgt_vec_comps[i]
                    AUX2 = AUX2.substitute(parameters_values)
                    initial_tgt_vec_comps[i] = AUX2
                # 'AUX' and 'AUX2' only used for the lines of
                # source code to be shorter

        t_min = numerical_approx(t_min)
        t_max = numerical_approx(t_max)

        for i in range(dim):
            if not isinstance(eqns_rhs[i], Expression): # in case of a
            # right hand side that is not an Expression (and then is a
            # number), it is needed to be converted to an Expression
            # since some solvers called below require only expressions
                eqns_rhs[i] = SR(eqns_rhs[i])

        if step is None:
            step = (t_max - t_min)/100

        initial_pt_coords = [numerical_approx(coord) for coord
                             in initial_pt_coords]
        initial_tgt_vec_comps = [numerical_approx(comp) for comp
                                 in initial_tgt_vec_comps]
        # the last two instructions retrieve numerical values even
        # if no parameters had to be substituted, in case some
        # coordinates or components contain expressions such as pi,
        # or are not RealNumber, since variable 'ics' of
        # 'desolve_system_rk4' used below needs to be a list of
        # RealNumber

        if not chart.valid_coordinates(*initial_pt_coords):
            raise ValueError("Initial point should be in the " +
                             "domain of the chart.")

        ode_solver_methods = ["rk2","rk4","rkf45","rkck","rk8pd"]
        ode_solver_methods+= ["rk2imp","rk4imp","gear1","gear2","bsimp"]

        if method == 'rk4_maxima':
            des = self._velocities + eqns_rhs
            dvars = list(chart[:]) + self._velocities
            ics = [t_min] + initial_pt_coords + initial_tgt_vec_comps

            sol = desolve_system_rk4(des, dvars,
                                     ivar=self._curve_parameter,
                                     ics=ics,
                                     end_points=[t_min, t_max],
                                     step=step)
        elif method == "ode_int":
            des = self._velocities + eqns_rhs
            ics = initial_pt_coords + initial_tgt_vec_comps
            times = srange(t_min, t_max, step, include_endpoint=True)
            dvars = list(chart[:]) + self._velocities

            sol0 = desolve_odeint(des, ics, times, dvars,
                                             ivar=self._curve_parameter)

            # rewrite the solution to prepare for the extraction (which
            # removes information about the velocities), and convert
            # elements of type 'numpy.float64' to standard type 'float'
            sol = []
            for t, coords_array in zip(times, sol0):
                coords_values = [float(coord_value) for coord_value
                                 in coords_array ]
                sol += [ [t] + coords_values ]
        elif method in ode_solver_methods:
            T = self._ode_solver

            if T is None:
                def system(t,y):
                    syst = self._velocities + eqns_rhs
                    par = self._curve_parameter
                    for i in range(dim):
                        vel = self._velocities[i]
                        syst[i] = syst[i].substitute({vel:y[dim+i]})
                        syst[dim+i] = syst[dim+i].substitute({par:t})
                        for j in range(dim):
                            coord = chart[:][j] # important to use '[:]'
                            # on 'chart' to avoid problems due to non
                            # zero starting index (i0)
                            veloc = self._velocities[j]
                            syst[dim+i]=syst[dim+i].substitute({coord:y[j]})
                            syst[dim+i]=syst[dim+i].substitute({veloc:y[dim+j]})
                    return syst
                from sage.calculus.ode import ode_solver
                T = ode_solver(function=system)

            T.algorithm = method
            y_0 = initial_pt_coords + initial_tgt_vec_comps
            t_span = srange(t_min, t_max, step, include_endpoint=True)

            if method == "bsimp":
                # this method requires the expression of the Jacobian
                # matrix of the application defining the right-hand side
                # of the system to be provided

                if T.jacobian is None:
                    def jacobian(t,y):
                        jac = []
                        par = self._curve_parameter
                        for i in range(dim):
                            new_row = [0 for j in range(2*dim)]
                            new_row[dim + i] = 1
                            jac += [new_row]

                        for i in range(dim):
                            semi_row_coords = []
                            semi_row_vels = []
                            for j in range(dim):
                                coord = chart[:][j] # important to use
                                # '[:]' on 'chart' to avoid problems due
                                # to non zero starting index (i0)
                                vel = self._velocities[j]
                                AUX = eqns_rhs[i].derivative(coord)
                                AUX2 = eqns_rhs[i].derivative(vel)
                                AUX = AUX.substitute({par:t})
                                AUX2 = AUX2.substitute({par:t})
                                for k in range(dim):
                                    coordin = chart[:][k] # important to
                                    # use '[:]' on 'chart' to avoid
                                    # problems due to non zero starting
                                    # index (i0)
                                    veloc = self._velocities[k]
                                    AUX = AUX.substitute({coordin:y[k]})
                                    AUX = AUX.substitute({veloc:y[dim+k]})
                                    AUX2 = AUX2.substitute({coordin:y[k]})
                                    AUX2 = AUX2.substitute({veloc:y[dim+k]})
                                semi_row_coords += [AUX]
                                semi_row_vels += [AUX2]
                            jac += [semi_row_coords + semi_row_vels]

                        last_semi_row_coords = [0 for j in range(dim)]
                        last_semi_row_vels = []
                        for j in range(dim):
                            AUX3 = eqns_rhs[j].derivative(par)
                            AUX3 = AUX3.substitute({par:t})
                            for m in range(dim):
                                coordin = chart[:][m] # important to use
                                # '[:]' on 'chart' to avoid problems due
                                # to non zero starting index (i0)
                                veloc = self._velocities[m]
                                AUX3 = AUX3.substitute({coordin:y[m]})
                                AUX3 = AUX3.substitute({veloc:y[dim+m]})
                            last_semi_row_vels += [AUX3]
                        jac += [last_semi_row_coords + last_semi_row_vels]
                        # 'AUX', 'AUX2' and 'AUX3' only used for the lines
                        # of source code to be shorter
                        return jac
                    T.jacobian = jacobian

                T.ode_solve(jacobian=jacobian, y_0=y_0, t_span=t_span)
            else:
                T.ode_solve(y_0=y_0, t_span=t_span)

            sol0 = T.solution
            sol = []
            for point in sol0:
                sol += [[point[0]] + point[1]]
            # above loop rewrites the solution in the same form than
            # that provided by other methods ('rk4_maxima' and
            # 'ode_int'), in order to extract the time and corresponding
            # coordinate values a few lines below, in the same way for
            # all methods

        else:
            raise ValueError("No available method of integration " +
                             "referred to as '{}'.".format(method))

        # eventually, extract the time and corresponding coordinate
        # values from each point of the solution computed (thus removing
        # information about the values of the velocities ; should the
        # latter be conserved ? They could turn useful in method
        # 'tangent_vector_eval_at', and in 'plot' when plotting the
        # tangent vectors.)
        coords_sol = [point[0:dim+1] for point in sol]

        if verbose:
            print("\nNumerical integration completed.\n" +
                  "Checking all points are in the chart domain...")

        N = len(coords_sol)
        n = 0
        while n < N and chart.valid_coordinates(*coords_sol[n][1:dim+1]):
            n += 1

        if n < N:
            raise ValueError("The {}th point ".format(n) +
                             "of the numerical solution (obtained at " +
                             "time {}) is out ".format(t_min + n*step) +
                             "of the chart domain. A curve with a " +
                             "smaller maximal value of the curve " +
                             "parameter, or a smaller initial tangent "+
                             "vector might be considered.")
        else:
            self._solutions[solution_key] = coords_sol
            if verbose:
                print("\nAll points are in the chart domain.\n" +
                      "The resulting list of points was associated " +
                      "with the key '{}' ".format(solution_key) +
                      "(if this key already referred to a former " +
                      "numerical solution, such a solution was erased).")
            return self._solutions[solution_key]

    def solution(self, solution_key=None, verbose=True):
        r"""
        Return the solution (list of points) associated with the given
        key.

        INPUT:

        - ``solution_key`` -- (default: ``None``) key which the
          requested numerical solution is associated to ; a default
          value is chosen if none is provided
        - ``verbose`` -- (default: ``True``) prints information about
          the solution returned

        OUTPUT:

        - list of the numerical points of the solution requested

        EXAMPLE:

        Requesting a numerical solution previously computed::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c',
            ....:                                         verbose=False)
            sage: sol = c.solve(verbose=False,
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1},
            ....:        solution_key='sol_T1')
            sage: sol_bis = c.solution()
            Returning the numerical solution associated with the key
             'sol_T1' by default...
            sage: sol_bis == sol
            True
            sage: sol_ter = c.solution(solution_key='sol_T1')
            sage: sol_ter == sol
            True
            sage: sol_mute = c.solution(verbose=False)
            sage: sol_mute == sol
            True

        """

        if solution_key is None:
            if 'rk4_maxima' in self._solutions.keys():
                solution_key = 'rk4_maxima'
            else:
                solution_key = self._solutions.keys()[0] # will raise
                # error if self._solutions empty
            if verbose:
                print("Returning the numerical solution associated " +
                      "with the key '{}' ".format(solution_key) +
                      "by default...")
        elif solution_key not in self._solutions.keys():
            raise ValueError("No existing key " +
                             "'{}' ".format(solution_key) +
                             "referring to any numerical solution.")

        return self._solutions[solution_key]

    def interpolate(self, solution_key=None, method=None,
                    interpolation_key=None, verbose=True):
        r"""
        Interpolate the chosen numerical solution using the given
        interpolation method.

        INPUT:

        - ``solution_key`` -- (default: ``None``) key which the
          numerical solution to interpolate is associated to ; a default
          value is chosen if none is provided
        - ``method`` -- (default: ``None``) interpolation scheme to use;
          algorithms available are

          * 'cubic spline', which makes use of ``GSL`` via Sage class
            :class:`~sage.calculus.interpolation.Spline`

        - ``interpolation_key`` -- (default: ``None``) key which the
          resulting interpolation will be associated to ; a default
          value is given if none is provided
        - ``verbose`` -- (default: ``True``) prints information about
          the interpolation in progress

        OUTPUT:

        - built interpolation object

        EXAMPLE:

        Interpolating a numerical solution previously computed::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c',
            ....:                                         verbose=False)
            sage: sol = c.solve(method='rk4_maxima',
            ....:        solution_key='sol_T1',
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1},
            ....:        verbose=False)
            sage: interp = c.interpolate(solution_key='my solution')
            Traceback (most recent call last):
            ...
            ValueError: No existing key 'my solution' referring to any
             numerical solution.
            sage: interp = c.interpolate(solution_key='sol_T1',
            ....:                        method='my method')
            Traceback (most recent call last):
            ...
            ValueError: No available method of interpolation referred to
             as 'my method'.
            sage: interp = c.interpolate(method='cubic spline',
            ....:                        solution_key='sol_T1',
            ....:                        interpolation_key='interp_T1')
            Interpolation completed and associated with the key
             'interp_T1' (if this key already referred to a former
             interpolation, such an interpolation was erased).
            sage: interp = c.interpolate()
            Interpolating the numerical solution associated with the
             key 'sol_T1' by default...
            Performing cubic spline interpolation by default...
            Resulting interpolation will be associated with the key
             'cubic spline-interp-sol_T1' by default.
            Interpolation completed and associated with the key
             'cubic spline-interp-sol_T1' (if this key already referred
             to a former interpolation, such an interpolation was
             erased).

        """

        if solution_key is None:
            if 'rk4_maxima' in self._solutions.keys():
                solution_key = 'rk4_maxima'
            else:
                solution_key = self._solutions.keys()[0] # will raise
                # error if self._solutions empty
            if verbose:
                print("Interpolating the numerical solution " +
                      "associated with the key " +
                      "'{}' ".format(solution_key) +
                      "by default...")
        elif solution_key not in self._solutions.keys():
            raise ValueError("No existing key " +
                             "'{}' ".format(solution_key) +
                             "referring to any numerical solution.")

        if method is None:
            method = 'cubic spline'
            if verbose:
                print("Performing cubic spline interpolation by "
                      "default...")

        if interpolation_key is None:
            interpolation_key = "{}-interp-".format(method)
            interpolation_key += "{}".format(solution_key)
            if verbose:
                print("Resulting interpolation will be associated " +
                      "with the key '{}' ".format(interpolation_key) +
                      "by default.")

        if method=='cubic spline':
            self._interpolations[interpolation_key] = []
            dim = self.codomain().dim()
            for i in range(dim):
                coordinate_curve = []
                for point in self._solutions[solution_key]:
                    coordinate_curve += [[point[0], point[i+1]]]
                self._interpolations[interpolation_key]+=[Spline(coordinate_curve)]
        else:
            raise ValueError("No available method of interpolation " +
                             "referred to as '{}'.".format(method))

        if verbose:
            print("Interpolation completed and associated with the " +
                  "key '{}' ".format(interpolation_key) +
                  "(if this key already referred to a former " +
                  "interpolation, such an interpolation was erased).")

        return self._interpolations[interpolation_key]

    def interpolation(self, interpolation_key=None, verbose=True):
        r"""
        Return the interpolation object associated with the given key.

        INPUT:

        - ``interpolation_key`` -- (default: ``None``) key which the
          requested interpolation is associated to ; a default
          value is chosen if none is provided
        - ``verbose`` -- (default: ``True``) prints information about
          the interpolation object returned

        OUTPUT:

        - requested interpolation object

        EXAMPLE:

        Requesting an interpolation object previously computed::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c',
            ....:                                         verbose=False)
            sage: sol = c.solve(method='rk4_maxima',
            ....:        solution_key='sol_T1',
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1},
            ....:        verbose=False)
            sage: interp = c.interpolate(method='cubic spline',
            ....:                  solution_key='sol_T1', verbose=False,
            ....:                  interpolation_key='interp_T1')
            sage: c.interpolation(interpolation_key='my interp')
            Traceback (most recent call last):
            ...
            ValueError: No existing key 'my interp' referring to any
             interpolation.
            sage: default_interp = c.interpolation()
            Returning the interpolation associated with the key
             'interp_T1' by default...
            sage: default_interp == interp
            True
            sage: interp_mute = c.interpolation(verbose=False)
            sage: interp_mute == interp
            True

        """

        if interpolation_key==None:
            if 'cubic spline' in self._interpolations.keys():
                interpolation_key = 'cubic spline'
            else:
                interpolation_key = self._interpolations.keys()[0]#will
                # raise error if self._interpolations empty
            if verbose:
                print("Returning the interpolation associated with " +
                      "the key '{}' ".format(interpolation_key) +
                      "by default...")
        elif interpolation_key not in self._interpolations.keys():
            raise ValueError("No existing key " +
                             "'{}' ".format(interpolation_key) +
                             "referring to any interpolation.")

        return self._interpolations[interpolation_key]

    def __call__(self, t, interpolation_key=None,
                 verbose=True):
        r"""
        Return the image of the curve for the given value of the curve
        parameter, using the chosen interpolation.

        INPUT:

        - ``t'' -- curve parameter value at which the curve is evaluated
        - ``interpolation_key`` -- (default: ``None``) key which the
          interpolation requested to compute the point is associated to ;
          a default value is chosen if none is provided
        - ``verbose`` -- (default: ``True``) prints information about
          the interpolation used

        OUTPUT:

        - list of the numerical coordinates of the point


        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c',
            ....:                                         verbose=False)
            sage: sol = c.solve(method='rk4_maxima',
            ....:        solution_key='sol_T1',
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1},
            ....:        verbose=False)
            sage: interp = c.interpolate(method='cubic spline',
            ....:                  solution_key='sol_T1', verbose=False,
            ....:                  interpolation_key='interp_T1')
            sage: c(1.1, interpolation_key='my interp')
            Traceback (most recent call last):
            ...
            ValueError: No existing key 'my interp' referring to any
             interpolation.
            sage: c(1.1)
            Evaluating point coordinates from the interpolation
             associated with the key 'interp_T1' by default...
            [1.060743343394347, -0.2153835404373033, 1.1]
            sage: pt = c(1.1, verbose=False); pt
            [1.060743343394347, -0.2153835404373033, 1.1]

        """

        if interpolation_key==None:
            if 'cubic spline' in self._interpolations.keys():
                interpolation_key = 'cubic spline'
            else:
                interpolation_key = self._interpolations.keys()[0]#will
                # raise error if self._interpolations empty
            if verbose:
                print("Evaluating point coordinates from the " +
                  "interpolation associated with the key " +
                  "'{}' by default...".format(interpolation_key))
        elif interpolation_key not in self._interpolations.keys():
            raise ValueError("No existing key " +
                             "'{}' ".format(interpolation_key) +
                             "referring to any interpolation.")

        interpolation = self._interpolations[interpolation_key]

        if isinstance(interpolation[0], Spline): #partial test, in case
        # future interpolation objects do not contain lists of instances
        # of the Spline class
            interpolated_coordinates=[coordinate_curve_spline(t)
                           for coordinate_curve_spline in interpolation]
            return interpolated_coordinates

        raise TypeError("Unexpected type of interpolation object.")

    def tangent_vector_eval_at(self, t,
                               interpolation_key=None, verbose=True):
        r"""
        Return the vector tangent to the curve at the given curve
        parameter with components evaluated from the given
        interpolation.

        INPUT:

        - ``t`` -- curve parameter value at which the tangent vector is
          evaluated
        - ``interpolation_key`` -- (default: ``None``) key which the
          interpolation requested to compute the tangent vector is
          associated to ; a default value is chosen if none is provided
        - ``verbose`` -- (default: ``True``) prints information about
          the interpolation used

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
          tangent vector with numerical components

        EXAMPLE:

        Evaluating a vector tangent to the curve::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c',
            ....:                                         verbose=False)
            sage: sol = c.solve(method='rk4_maxima',
            ....:        solution_key='sol_T1',
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1},
            ....:        verbose=False)
            sage: interp = c.interpolate(method='cubic spline',
            ....:                  solution_key='sol_T1', verbose=False,
            ....:                  interpolation_key='interp_T1')
            sage: tg_vec = c.tangent_vector_eval_at(1.22,
            ....:                         interpolation_key='my interp')
            Traceback (most recent call last):
            ...
            ValueError: No existing key 'my interp' referring to any
             interpolation.
            sage: tg_vec = c.tangent_vector_eval_at(1.22); tg_vec
            Evaluating tangent vector components from the interpolation
             associated with the key 'interp_T1' by default...
            Tangent vector at Point on the 3-dimensional differentiable
             manifold M
            sage: tg_vec[:]
            [0.7392639473853356, -0.6734182305341726, 1.0000000000000007]
            sage: tg_vec_mute = c.tangent_vector_eval_at(1.22,
            ....:                         verbose=False,
            ....:                         interpolation_key='interp_T1')
            sage: tg_vec_mute == tg_vec
            True

        """

        if interpolation_key==None:
            if 'cubic spline' in self._interpolations.keys():
                interpolation_key = 'cubic spline'
            else:
                interpolation_key = self._interpolations.keys()[0]#will
                # raise error if self._interpolations empty
            if verbose:
                print("Evaluating tangent vector components from the " +
                      "interpolation associated with the key " +
                      "'{}' by default...".format(interpolation_key))
        elif interpolation_key not in self._interpolations.keys():
            raise ValueError("No existing key " +
                             "'{}' ".format(interpolation_key) +
                             "referring to any interpolation.")

        interpolation = self._interpolations[interpolation_key]

        if isinstance(interpolation[0], Spline): #partial test, in case
        # future interpolation objects do not contain lists of instances
        # of the Spline class
            interpolated_coordinates=[coordinate_curve_spline(t)
                           for coordinate_curve_spline in interpolation]
            M = self.codomain()
            p = M.point(interpolated_coordinates, chart=self._chart,
                        name=None)
            Tp = M.tangent_space(p)

            evaluated_tgt_vec_comp=[coordinate_curve_spline.derivative(t)
                    for coordinate_curve_spline in interpolation] # by
            # default, order=1 in method 'derivative' of a class Spline
            basis = self._chart.frame().at(p)
            v = Tp(evaluated_tgt_vec_comp, basis=basis)
            return v

        raise TypeError("Unexpected type of interpolation object.")

    @options(thickness=1, width_tangent=1, plot_points=75,
             aspect_ratio='automatic', plot_points_tangent=10, scale=1)
    def plot_integrated(self, chart=None, ambient_coords=None,
             mapping=None, prange=None, interpolation_key=None,
             include_end_point=(True, True),
             end_point_offset=(0.001, 0.001), verbose=True, color='red',
             style='-', label_axes=True, display_tangent=False,
             color_tangent='blue', **kwds):
        r"""
        Plot the 2D or 3D projection of the curve onto the space of the
        chosen two or three ambient coordinates, based on the
        interpolation of a numerical solution previously computed.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve.plot`
            for complete information about the input.

        ADDITIONAL INPUT:

        - ``interpolation_key`` -- (default: ``None``) key associated to
          the interpolation object used for the plot ; a default value
          is chosen if none is provided
        - ``verbose`` -- (default: ``True``) prints information about
          the interpolation object used and the plotting in progress

        EXAMPLE:

        Trajectory of a particle of unit mass and unit charge in an
        unit, axial, uniform, stationnary magnetic field::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: var('t')
            t
            sage: D = X.symbolic_velocities()
            sage: eqns = [D[1], -D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 6), v,
            ....:                               name='c', verbose=False)
            sage: sol = c.solve(verbose=False)
            sage: interp = c.interpolate(verbose=False)
            sage: c_plot_2d = c.plot_integrated(ambient_coords=[x1, x2],
            ....:                 thickness=2.5,
            ....:                 display_tangent=True, plot_points=200,
            ....:                 plot_points_tangent=10, scale=0.5,
            ....:                 color='blue', color_tangent='red')
            Plotting from the interpolation associated with the key
             'cubic spline-interp-rk4_maxima' by default...
            A tiny final offset equal to the value of
             'end_point_offset[1]' (= 0.001) was introduced in order to
             safely compute the last tangent vector from the
             interpolation.
            sage: c_plot_2d.show()

        .. PLOT::

            M = Manifold(3, 'M')
            X = M.chart('x1 x2 x3')
            var('x1 x2 x3 t')
            D = X.symbolic_velocities()
            eqns = [D[1], -D[0], 0]
            p = M.point((0,0,0), name='p')
            Tp = M.tangent_space(p)
            v = Tp((1,0,1))
            c = M.integrated_curve(eqns, D, (t, 0, 6), v, name='c',
            ....:                                         verbose=False)
            c.solve(verbose=False)
            c.interpolate(verbose=False)
            c_plot_2d_1 = c.plot_integrated(ambient_coords=[x1, x2],
                            thickness=2.5,
                            display_tangent=True, plot_points=200,
                            plot_points_tangent=10, scale=0.5,
                            color='blue', color_tangent='red')
            sphinx_plot(c_plot_2d_1)

        """

        from sage.manifolds.chart import RealChart

        #
        # Get the @options from kwds
        #
        thickness = kwds.pop('thickness')
        plot_points = kwds.pop('plot_points')
        aspect_ratio = kwds.pop('aspect_ratio')

        #
        # Interpolation to use
        #
        if interpolation_key==None:
            if 'cubic spline' in self._interpolations.keys():
                interpolation_key = 'cubic spline'
            else:
                interpolation_key = self._interpolations.keys()[0]#will
                # raise error if self._interpolations empty
            if verbose:
                print("Plotting from the interpolation associated " +
                  "with the key '{}' ".format(interpolation_key) +
                  "by default...")
        elif interpolation_key not in self._interpolations.keys():
            raise ValueError("No existing key " +
                             "'{}' ".format(interpolation_key) +
                             "referring to any interpolation.")

        interpolation = self._interpolations[interpolation_key]

        #
        # The mapping, if present, and the chart w.r.t. which the curve
        # is plotted
        #
        if mapping is None:
            i0 = self.codomain().start_index()
            if chart is None:
                chart = self._chart
            else:
                if not isinstance(chart, RealChart):
                    raise TypeError("{} is not a real " +
                                    "chart".format(chart))
                mapping = self.codomain().identity_map()
        else:
            i0 = mapping.codomain().start_index()
            if chart is None:
                chart = mapping.codomain().default_chart()
            elif not isinstance(chart, RealChart):
                raise TypeError("{} is not a real chart".format(chart))

        #
        # Coordinates of the above chart w.r.t. which the curve is
        # plotted
        #
        if ambient_coords is None:
            ambient_coords = chart[:]  # all chart coordinates are used
        n_pc = len(ambient_coords)
        if n_pc != 2 and n_pc !=3:
            raise ValueError("the number of coordinates involved in " +
                             "the plot must be either 2 or 3, " +
                             "not {}".format(n_pc))

        # From now on, 'pc' will denote coordinates in terms of which
        # the curve is plotted (i.e. the "ambient coordinates"), while
        # 'coord' will denote coordinates on self.domain(); of course,
        # when, for instance, the mapping is the identity map, these may
        # be the same.

        # indices of plot coordinates
        ind_pc = [chart[:].index(pc)+i0 for pc in ambient_coords] # will
        # raise an error if ambient_coords are not associated with chart

        #
        # Maximal parameter range for the plot of the chosen
        # interpolation
        #
        param_min = interpolation[0][0][0]
        param_max = interpolation[0][-1][0]

        if prange is None:
            prange = (param_min, param_max)
        elif not isinstance(prange, (tuple, list)):
            raise TypeError("{} is neither ".format(prange) +
                            "a tuple nor a list")
        elif len(prange) != 2:
            raise ValueError("the argument prange must be a " +
                             "tuple/list of 2 elements")
        else:
            p = prange #'p' declared only for the line below to be shorter
            if p[0]<param_min or p[0]>param_max or p[1]<param_min or p[1]>param_max:
                raise ValueError("Parameter range should be a " +
                                "subinterval of the curve domain " +
                                "({}).".format(self.domain()))

        tmin = numerical_approx(prange[0])
        tmax = numerical_approx(prange[1])

        if not include_end_point[0]:
            tmin += numerical_approx(end_point_offset[0])

        if not include_end_point[1]:
            tmax -= numerical_approx(end_point_offset[1])

        if mapping is None:
            if isinstance(interpolation[0], Spline): #partial test,
            # in case future interpolation objects do not contain lists
            # of instances of the Spline class

                #
                # List of points for the plot curve
                #
                plot_curve = []
                dt = (tmax - tmin) / (plot_points - 1)
                t = tmin

                for k in range(plot_points):
                    plot_curve.append([interpolation[j-i0](t) for j in ind_pc])
                    t += dt

                if display_tangent:
                    from sage.plot.graphics import Graphics
                    from sage.plot.arrow import arrow2d
                    from sage.plot.plot3d.shapes import arrow3d

                    scale = kwds.pop('scale')
                    plot_points_tangent=kwds.pop('plot_points_tangent')
                    width_tangent = kwds.pop('width_tangent')

                    if not tmax < param_max:
                        tmax=numerical_approx(tmax- end_point_offset[1])
                        if verbose:
                            print("A tiny final offset equal to " +
                                  "the value of " +
                                  "'end_point_offset[1]' " +
                                  "(= {}) ".format(end_point_offset[1])+
                                  "was introduced " +
                                  "in order to safely compute the " +
                                  "last tangent vector from the " +
                                  "interpolation.")

                    plot_vectors = Graphics()
                    dt = (tmax - tmin) / (plot_points_tangent - 1)
                    t = tmin

                    for k in range(plot_points_tangent):
                        # interpolated ambient coordinates:
                        xp = [interpolation[j-i0](t) for j in ind_pc]

                        # tangent vector ambiant components evaluated
                        # from the interpolation:
                        vec = [coordinate_curve_spline.derivative(t)
                           for coordinate_curve_spline in interpolation]

                        coord_tail = xp
                        coord_head = [xp[j] + scale*vec[j]
                                      for j in range(len(xp))]

                        if coord_head != coord_tail:
                            if n_pc == 2:
                                plot_vectors+=arrow2d(tailpoint=coord_tail,
                                                   headpoint=coord_head,
                                                   color=color_tangent,
                                                   width=width_tangent)
                            else:
                                plot_vectors+=arrow3d(coord_tail,
                                                    coord_head,
                                                    color=color_tangent,
                                                    width=width_tangent)
                        t += dt

                    return plot_vectors+DifferentiableCurve._graphics(self,
                                             plot_curve, ambient_coords,
                                             thickness=thickness,
                                             aspect_ratio=aspect_ratio,
                                             color=color,
                                             style=style,
                                             label_axes=label_axes)

                return DifferentiableCurve._graphics(self, plot_curve,
                                 ambient_coords, thickness=thickness,
                                 aspect_ratio=aspect_ratio, color=color,
                                 style=style, label_axes=label_axes)

            raise TypeError("Unexpected type of interpolation object.")
        else:
            #
            # The coordinate expressions of the mapping and the
            # coordinates involved
            #
            for chart_pair in mapping._coord_expression.keys():
                subs=(chart_pair[0]._subcharts,chart_pair[1]._subcharts)
                # 'subs' declared only for the line below to be shorter
                if self._chart in subs[0] and chart in subs[1]:
                    transf = {}
                    required_coords = set()
                    for pc in ambient_coords:
                        j = chart[:].index(pc)
                        AUX = mapping._coord_expression[chart_pair]
                        # 'AUX' used only for the lines of source code
                        # to be shorter
                        transf[pc] = AUX.expr()[j]
                        AUX2 = transf[pc].variables() # idem
                        required_coords=required_coords.union(AUX2)
                    break
            else:
                raise ValueError("No expression has been found for " +
                                 "{} in terms of {}".format(self,chart))

            if isinstance(interpolation[0], Spline): # partial test, in
            # case future interpolation objects do not contain lists of
            # instances of the Spline class

                #
                # List of points for the plot curve
                #
                plot_curve = []
                dt = (tmax - tmin) / (plot_points - 1)
                t = tmin
                required_coords_values = {}

                for k in range(plot_points):
                    for coord in required_coords:
                        i = self._chart[:].index(coord)
                        required_coords_values[coord]=interpolation[i](t)

                    xp = []
                    for j in ind_pc:
                        pc = chart[j]
                        AUX = transf[pc]
                        AUX = AUX.substitute(required_coords_values)
                        # 'AUX' only used for the lines of source code
                        #  to be shorter
                        xp+=[numerical_approx(AUX)]

                    plot_curve.append(xp)
                    t += dt

                if display_tangent:
                    from sage.plot.graphics import Graphics
                    from sage.plot.arrow import arrow2d
                    from sage.plot.plot3d.shapes import arrow3d

                    scale = kwds.pop('scale')
                    plot_points_tangent=kwds.pop('plot_points_tangent')
                    width_tangent = kwds.pop('width_tangent')

                    if not tmax < param_max:
                        tmax=numerical_approx(tmax- end_point_offset[1])
                        if verbose:
                            print("A tiny final offset equal to " +
                                  "the value of " +
                                  "'end_point_offset[1]' " +
                                  "(= {}) ".format(end_point_offset[1])+
                                  "was introduced " +
                                  "in order to safely compute the " +
                                  "last tangent vector from the " +
                                  "interpolation.")

                    plot_vectors = Graphics()
                    dt = (tmax - tmin) / (plot_points_tangent - 1)
                    t = tmin
                    Dcoord_Dt = {}

                    Dpc_Dcoord = {}
                    for pc in ambient_coords:
                        Dpc_Dcoord[pc] = {}
                        for coord in transf[pc].variables():
                            Dpc_Dcoord[pc][coord]=transf[pc].derivative(coord)

                    for k in range(plot_points_tangent):
                        for coord in required_coords:
                            i = self._chart[:].index(coord)
                            AUX = interpolation[i] # 'AUX' only used
                            # for the lines below to be shorter
                            required_coords_values[coord] = AUX(t)
                            Dcoord_Dt[coord] = AUX.derivative(t)

                        xp = []
                        pushed_vec = []
                        for j in ind_pc:
                            pc = chart[j]
                            AUX = transf[pc]
                            AUX = AUX.substitute(required_coords_values)
                            # 'AUX' only used for the lines of code to
                            # be shorter
                            xp+=[numerical_approx(AUX)]

                            pushed_comp = 0
                            for coord in transf[pc].variables():
                                D = Dpc_Dcoord[pc][coord]
                                D = D.substitute(required_coords_values)
                                D=numerical_approx(D)
                                pushed_comp += Dcoord_Dt[coord] * D

                            pushed_vec += [pushed_comp]

                        coord_tail = xp
                        coord_head =[xp[j] + scale*pushed_vec[j]
                                     for j in range(len(xp))]

                        if coord_head != coord_tail:
                            if n_pc == 2:
                                plot_vectors+=arrow2d(tailpoint=coord_tail,
                                                   headpoint=coord_head,
                                                   color=color_tangent,
                                                   width=width_tangent)
                            else:
                                plot_vectors+=arrow3d(coord_tail,
                                                    coord_head,
                                                    color=color_tangent,
                                                    width=width_tangent)
                        t += dt

                    return plot_vectors+DifferentiableCurve._graphics(self,
                                             plot_curve, ambient_coords,
                                             thickness=thickness,
                                             aspect_ratio=aspect_ratio,
                                             color=color,
                                             style=style,
                                             label_axes=label_axes)

                return DifferentiableCurve._graphics(self, plot_curve,
                                 ambient_coords, thickness=thickness,
                                 aspect_ratio=aspect_ratio, color=color,
                                 style=style, label_axes=label_axes)

            raise TypeError("Unexpected type of interpolation object.")

class IntegratedAutoparallelCurve(IntegratedCurve):
    r"""
    Numerical autoparallel curve on the manifold with
    respect to a given affine connection.

    INPUT:

    - ``parent`` --
      :class:`~sage.manifolds.differentiable.manifold_homset.DifferentiableCurveSet`
      the set of curves `\mathrm{Hom}(I, M)` to which the curve belongs
    - ``affine_connection`` --
      :class:`~sage.manifolds.differentiable.affine_connection.AffineConnection`
      affine connection with respect to which the curve is
      autoparallel
    - ``curve_parameter`` -- symbolic expression to be used as the
      parameter of the curve (the equations defining an instance of
      IntegratedAutoparallelCurve are such that ``t`` will actually be
      an affine parameter of the curve)
    - ``initial_tangent_vector`` --
      :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
      initial tangent vector of the curve
    - ``chart`` -- (default: ``None``) chart on the manifold in terms of
      which the equations are expressed; if ``None`` the default chart
      of the manifold is assumed
    - ``name`` -- (default: ``None``) string; symbol given to the curve
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the curve; if none is provided, ``name`` will be used
    - ``is_isomorphism`` -- (default: ``False``) determines whether the
      constructed object is a diffeomorphism; if set to ``True``,
      then `M` must have dimension one
    - ``is_identity`` -- (default: ``False``) determines whether the
      constructed object is the identity map; if set to ``True``,
      then `M` must coincide with the domain of the curve

    EXAMPLE:

    Autoparallel curves associated with the Mercator projection of the
    unit 2-sphere :MATH:`\mathbb{S}^{2}`.

    .. SEEALSO::

        https://idontgetoutmuch.wordpress.com/2016/11/24/mercator-a-connection-with-torsion/
        for more details about Mercator projection

    On the Mercator projection, the lines of longitude all appear
    vertical and then all parallel w.r.t each others.
    Likewise, all the lines of latitude appear horizontal and parallel
    w.r.t each others.
    These curves may be recovered as autoparallel curves of a certain
    connection :MATH:`\nabla` to be made explicit.

    Start with declaring the standard polar coordinates
    :MATH:`(\theta, \phi)` on :MATH:`\mathbb{S}^{2}` and the
    corresponding coordinate frame :MATH:`(e_{\theta}, e_{\phi})`::

        sage: S2 = Manifold(2, 'S^2', start_index=1)
        sage: polar.<th,ph>=S2.chart()
        sage: epolar = polar.frame()

    Normalizing :MATH:`e_{\phi}` provides an orthonormal basis::

        sage: ch_basis = S2.automorphism_field()
        sage: ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        sage: epolar_ON = epolar.new_frame(ch_basis,'epolar_ON')

    Denote :MATH:`(\hat{e}_{\theta}, \hat{e}_{\phi})` such an
    orthonormal frame field.
    In any point, the vector field :MATH:`\hat{e}_{\theta}` is
    normalized and tangent to the line of longitude through the point.
    Likewise, :MATH:`\hat{e}_{\phi}` is normalized and tangent to the
    line of latitude.

    Now, set an affine connection with respect to which such fields are
    parallely transported in all directions, that is:
    :MATH:`\nabla \hat{e}_{\theta} = \nabla \hat{e}_{\phi} = 0`.
    This is equivalent to setting all the connection coefficients to
    zero w.r.t. this frame::

        sage: nab = S2.affine_connection('nab')
        sage: nab.set_coef(epolar_ON)[:]
        [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]

    This connection is such that two vectors are parallel if their
    angles to a given meridian are the same.
    Check that this connection is compatible with the Euclidean
    metric tensor :MATH:`g` induced on :MATH:`\mathbb{S}^{2}`::

        sage: g = S2.metric('g')
        sage: g[1,1], g[2,2] = 1, (sin(th))^2
        sage: nab(g)[:]
        [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]

    Yet, this connection is not the Levi-Civita connection, which
    implies that it has non-vanishing torsion::

        sage: nab.torsion()[:]
        [[[0, 0], [0, 0]], [[0, cos(th)/sin(th)], [-cos(th)/sin(th), 0]]]

    Set generic initial conditions for the autoparallel curves to
    compute::

        sage: [th0, ph0, v_th0, v_ph0] = var('th0 ph0 v_th0 v_ph0')
        sage: p = S2.point((th0, ph0), name='p')
        sage: Tp = S2.tangent_space(p)
        sage: v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))

    Note here that the components ``(v\_th0, v\_ph0)`` of the initial
    tangent vector ``v`` refer to the basis
    ``epolar\_ON`` :MATH:` = (\hat{e}_{\theta}, \hat{e}_{\phi})`
    and not the coordinate basis
    ``epolar`` :MATH:` = (e_{\theta}, e_{\phi})`.
    This is merely to help picture the aspect of the tangent vector in
    the usual embedding of :MATH:`\mathbb{S}^{2}` in
    :MATH:`\mathbb{R}^{3}` thanks to using an orthonormal frame,
    since providing the components w.r.t. the coordinate basis would
    require mutliplying the second component (i.e. in :MATH:`\phi`) in
    order to picture the vector in the same way.
    This subtlety will need to be taken into account later when the
    numerical curve will be compared to the analytical solution.

    Now, declare the corresponding integrated autoparallel curve and display
    the differential system it satisfies::

        sage: [t, tmin, tmax] = var('t tmin tmax')
        sage: c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax),
        ....:                   v, chart=polar, name='c', verbose=False)
        sage: sys = c.system()
        Autoparallel curve c in the 2-dimensional differentiable
         manifold S^2 equipped with Affine connection nab on the
         2-dimensional differentiable manifold S^2, and integrated over
         the Real interval (tmin, tmax) as a solution to the following
         equations, written w.r.t. Chart (S^2, (th, ph)):
        <BLANKLINE>
        Initial point: Point p on the 2-dimensional differentiable
         manifold S^2 with coordinates [th0, ph0] w.r.t.
         Chart (S^2, (th, ph))
        Initial tangent vector: Tangent vector at Point p on the
         2-dimensional differentiable manifold S^2 with
         components [v_th0, v_ph0/sin(th0)] w.r.t. Chart (S^2, (th, ph))
        <BLANKLINE>
        d(th)/dt = Dth
        d(ph)/dt = Dph
        d(Dth)/dt = 0
        d(Dph)/dt = -Dph*Dth*cos(th)/sin(th)
        <BLANKLINE>

    Set a dictionnary providing the parameter range and the initial
    conditions for a line of latitude and a line of longitude::

        sage: dict_params={'latit':{tmin:0,tmax:3,th0:pi/4,ph0:0.1,v_th0:0,v_ph0:1},
        ....:   'longi':{tmin:0,tmax:3,th0:0.1,ph0:0.1,v_th0:1,v_ph0:0}}

    Declare the Mercator coordinates :MATH:`(\xi, \zeta)` and the
    corresponding coordinate change from the polar coordinates::

        sage: mercator.<xi,ze>=S2.chart(r'xi:(-oo,oo):\xi ze:(0,2*pi):\zeta')
        sage: polar.transition_map(mercator, (log(tan(th/2)), ph))
        Change of coordinates from Chart (S^2, (th, ph)) to Chart
         (S^2, (xi, ze))

    Ask for the identity map in terms of these charts in order to add
    this coordinate change to its dictionnary of expressions. This is
    required to plot the curve w.r.t the Mercator chart::

        sage: identity = S2.identity_map()
        sage: identity.coord_functions(polar, mercator)
        Coordinate functions (log(sin(1/2*th)/cos(1/2*th)), ph) on the
         Chart (S^2, (th, ph))

    Solve, interpolate and prepare the plot for the solutions
    corresponding to the two initial conditions previously set::

        sage: graph2D_mercator = Graphics()
        sage: for key in dict_params.keys():
        ....:     sol = c.solve(solution_key='sol-'+key,
        ....:         parameters_values=dict_params[key], verbose=False)
        ....:     interp = c.interpolate(solution_key='sol-'+key,
        ....:            interpolation_key='interp-'+key, verbose=False)
        ....:     graph2D_mercator+=c.plot_integrated(interpolation_key='interp-'+key,
        ....:                            chart=mercator, thickness=2,
        ....:                            verbose=False)

    Prepare a grid of Mercator coordinates lines, and plot the curves
    over it::

        sage: graph2D_mercator_coords=mercator.plot(chart=mercator,
        ....:                           number_values=8,color='yellow')
        sage: (graph2D_mercator + graph2D_mercator_coords).show()

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON')
        nab = S2.affine_connection('nab')
        nab.set_coef(epolar_ON)[:]
        [t,tmin,tmax,th0,ph0,v_th0,v_ph0]=var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                                  chart=polar, name='c')
        dict_params={'latit':{tmin:0,tmax:3,th0:pi/4,ph0:0.1,v_th0:0,v_ph0:1},
                'longi':{tmin:0,tmax:3,th0:0.1,ph0:0.1,v_th0:1,v_ph0:0}}
        mercator = S2.chart(r'xi:(-oo,oo):\xi ze:(0,2*pi):\zeta')
        [xi,ze] = var('xi ze')
        polar.transition_map(mercator, (log(tan(th/2)), ph))
        identity = S2.identity_map()
        identity.coord_functions(polar, mercator)
        graph2D_mercator = Graphics()
        for key in dict_params.keys():
            c.solve(solution_key='sol-'+key,
                      parameters_values=dict_params[key], verbose=False)
            c.interpolate(solution_key='sol-'+key,
                         interpolation_key='interp-'+key, verbose=False)
            graph2D_mercator += c.plot_integrated(interpolation_key='interp-'+key,
                             chart=mercator, thickness=2, verbose=False)
        graph2D_mercator_coords = mercator.plot(chart=mercator,
                                        number_values=8, color='yellow')
        sphinx_plot(graph2D_mercator + graph2D_mercator_coords)

    The resulting curves are horizontal and vertical as expected.
    It is easier to check that these are latitude and longitude lines
    respectively when plotting them on :MATH:`\mathbb{S}^{2}`.
    To do so, use :MATH:`\mathbb{R}^{3}` as the codomain of the standard
    map embedding :MATH:`(\mathbb{S}^{2}, (\theta, \phi))` in the
    3-dimensional Euclidean space::

        sage: R3 = Manifold(3, 'R3', start_index=1)
        sage: cart.<X,Y,Z> = R3.chart()
        sage: euclid_embedding = S2.diff_map(R3,
        ....:  {(polar, cart):[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})

    Plot the resulting curves on the grid of polar coordinates lines on
    :MATH:`\mathbb{S}^{2}`::

        sage: graph3D_embedded_curves = Graphics()
        sage: for key in dict_params.keys():
        ....:     graph3D_embedded_curves+=c.plot_integrated(interpolation_key='interp-'+key,
        ....:         mapping=euclid_embedding, thickness=5,
        ....:         display_tangent=True, scale=0.4, width_tangent=0.5,
        ....:         verbose=False)
        sage: graph3D_embedded_polar_coords = polar.plot(chart=cart,
        ....:                          mapping=euclid_embedding,
        ....:                          number_values=15, color='yellow')
        sage: graph=graph3D_embedded_curves+graph3D_embedded_polar_coords
        sage: graph.show()

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON')
        nab = S2.affine_connection('nab')
        nab.set_coef(epolar_ON)[:]
        [t,tmin,tmax,th0,ph0,v_th0,v_ph0]=var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                                  chart=polar, name='c')
        dict_params={'latit':{tmin:0,tmax:3,th0:pi/4,ph0:0.1,v_th0:0,v_ph0:1},
                'longi':{tmin:0,tmax:3,th0:0.1,ph0:0.1,v_th0:1,v_ph0:0}}
        R3 = Manifold(3, 'R3', start_index=1)
        cart = R3.chart('X Y Z')
        [X,Y,Z] = var('X Y Z')
        euclid_embedding = S2.diff_map(R3,
              {(polar, cart):[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})
        graph3D_embedded_curves = Graphics()
        for key in dict_params.keys():
            c.solve(solution_key='sol-'+key,
                      parameters_values=dict_params[key], verbose=False)
            c.interpolate(solution_key='sol-'+key,
                         interpolation_key='interp-'+key, verbose=False)
            graph3D_embedded_curves+=c.plot_integrated(interpolation_key='interp-'+key,
                mapping=euclid_embedding, thickness=5,
                display_tangent=True, scale=0.4, width_tangent=0.5,
                verbose=False)
        graph3D_embedded_polar_coords = polar.plot(chart=cart,
                                       mapping=euclid_embedding,
                                       number_values=15, color='yellow')
        graph = graph3D_embedded_curves+graph3D_embedded_polar_coords
        sphinx_plot(graph)

    Finally, one may plot a general autoparallel curve w.r.t
    :MATH:`\nabla` that is neither a line of latitude or longitude.
    The vectors tangent to such a curve make an angle different from 0
    or :MATH:`\pi/2` with the lines of latitude and longitude.
    Then, compute a curve such that both components of its initial
    tangent vectors are non zero::

        sage: sol = c.solve(solution_key='sol-angle',
        ....:  parameters_values={tmin:0,tmax:2,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:8},
        ....:  verbose=False)
        sage: interp = c.interpolate(solution_key='sol-angle',
        ....:           interpolation_key='interp-angle', verbose=False)

    Plot the resulting curve in the Mercator plane.
    This generates a straight line, as expected::

        sage: graph2D_mercator_angle_curve=c.plot_integrated(interpolation_key='interp-angle',
        ....:         chart=mercator, thickness=1, display_tangent=True,
        ....:         scale=0.2, width_tangent=0.2, verbose=False)
        sage: graph2D_mercator_angle_curve.show()

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON')
        nab = S2.affine_connection('nab')
        nab.set_coef(epolar_ON)[:]
        [t,tmin,tmax,th0,ph0,v_th0,v_ph0]=var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                                  chart=polar, name='c')
        mercator = S2.chart(r'xi:(-oo,oo):\xi ze:(0,2*pi):\zeta')
        [xi,ze] = var('xi ze')
        polar.transition_map(mercator, (log(tan(th/2)), ph))
        identity = S2.identity_map()
        identity.coord_functions(polar, mercator)
        sol = c.solve(solution_key='sol-angle',
         parameters_values={tmin:0,tmax:2,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:8},
         verbose=False)
        interp = c.interpolate(solution_key='sol-angle',
                        interpolation_key='interp-angle', verbose=False)
        graph2D_mercator_angle_curve=c.plot_integrated(interpolation_key='interp-angle',
                      chart=mercator, thickness=1, display_tangent=True,
                      scale=0.2, width_tangent=0.2, verbose=False)
        sphinx_plot(graph2D_mercator_angle_curve)

    One may eventually plot such a curve on :MATH:`\mathbb{S}^{2}`::

        sage: graph3D_embedded_angle_curve=c.plot_integrated(interpolation_key='interp-angle',
        ....:        mapping=euclid_embedding, thickness=5,
        ....:        display_tangent=True, scale=0.1, width_tangent=0.5,
        ....:        verbose=False)
        sage: graph=graph3D_embedded_angle_curve+graph3D_embedded_polar_coords
        sage: graph.show()

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON)
        nab = S2.affine_connection('nab')
        nab.set_coef(epolar_ON)[:]
        [t,tmin,tmax,th0,ph0,v_th0,v_ph0]=var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                                  chart=polar, name='c')
        R3 = Manifold(3, 'R3', start_index=1)
        cart = R3.chart('X Y Z')
        [X,Y,Z] = var('X Y Z')
        euclid_embedding = S2.diff_map(R3,
              {(polar, cart):[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})
        sol = c.solve(solution_key='sol-angle',
         parameters_values={tmin:0,tmax:2,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:8},
         verbose=False)
        interp = c.interpolate(solution_key='sol-angle',
                        interpolation_key='interp-angle', verbose=False)
        graph3D_embedded_angle_curve=c.plot_integrated(interpolation_key='interp-angle',
            mapping=euclid_embedding, thickness=5, display_tangent=True,
            scale=0.1, width_tangent=0.5, verbose=False)
        graph3D_embedded_polar_coords = polar.plot(chart=cart,
             mapping=euclid_embedding, number_values=15, color='yellow')
        graph=graph3D_embedded_angle_curve+graph3D_embedded_polar_coords
        sphinx_plot(graph)

    All the curves presented are loxodromes, and the differential system
    defining them (displayed above) may be solved analytically,
    providing the following expressions:

    .. MATH::

        \theta(t) &= \theta_{0} + \dot{\theta}_{0} (t - t_{0})      \\
        \phi(t) &= \phi_{0} - \frac{1}{\tan \alpha} \left(
        \ln \tan \frac{\theta_{0} + \dot{\theta}_{0} (t - t_{0})}{2} -
        \ln \tan \frac{\theta_{0}}{2} \right)

    where :MATH:`\alpha` is the angle between the curve and any latitude
    line it crosses; then, one finds
    :MATH:`\tan \alpha = - \dot{\theta}_{0}/(\dot{\phi}_{0} \sin \theta_{0})`)
    (then :MATH:`\tan \alpha \leq 0` when the initial tangent vector
    points towards the southeast).

    In order to use these expressions to compare with the result
    provided by the numerical integration, remember that the components
    ``(v\_th0, v\_ph0)`` of the initial
    tangent vector ``v`` refer to the basis
    ``epolar\_ON`` :MATH:`= (\hat{e}_{\theta}, \hat{e}_{\phi})` and not the
    coordinate basis
    ``epolar`` :MATH:`= (e_{\theta}, e_{\phi})`.
    Therefore, the following relations hold:
    ``v\_ph0`` = :MATH:`\dot{\phi}_{0} \sin \theta_{0}` (and not merely
    :MATH:`\dot{\phi}_{0}`), while ``v\_th0`` clearly is
    :MATH:`\dot{\theta}_{0}`.

    With this in mind, plot an analytical curve to compare with a
    numerical solution::

        sage: graph2D_mercator_angle_curve=c.plot_integrated(interpolation_key='interp-angle',
        ....:                chart=mercator, thickness=1, verbose=False)
        sage: expr_ph = ph0+v_ph0/v_th0*(ln(tan((v_th0*t+th0)/2))-ln(tan(th0/2)))
        sage: c_loxo = S2.curve({polar:[th0+v_th0*t, expr_ph]}, (t,0,2),
        ....:                                             name='c_loxo')

    Ask for the expression of the loxodrome in terms of the Mercator
    chart in order to add it to its dictionnary of expressions.
    It is a particularly long expression, and there is no particular
    need to diplay it, which is why it may simply be affected to an
    arbitrary variable ``expr_mercator``, which will never be used
    again.
    But adding the expression to the dictionnary is required to plot the
    curve w.r.t the Mercator chart::

        sage: expr_mercator = c_loxo.expression(chart2=mercator)

    Plot the curves (for clarity, set a 2 degrees shift in the initial
    value of :MATH:`\theta_{0}` so that the curves do not overlap)::

        sage: graph2D_mercator_loxo = c_loxo.plot(chart=mercator,
        ....:  parameters={th0:pi/4+2*pi/180, ph0:0.1, v_th0:1, v_ph0:8},
        ....:  thickness=1, color='blue')
        sage: (graph2D_mercator_angle_curve+graph2D_mercator_loxo).show()

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON')
        nab = S2.affine_connection('nab')
        nab.set_coef(epolar_ON)[:]
        [t, tmin, tmax, th0, ph0] = var('t tmin tmax th0 ph0')
        [v_th0, v_ph0, alpha] = var('v_th0 v_ph0 alpha')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                                  chart=polar, name='c')
        mercator = S2.chart(r'xi:(-oo,oo):\xi ze:(0,2*pi):\zeta')
        [xi,ze] = var('xi ze')
        polar.transition_map(mercator, (log(tan(th/2)), ph))
        identity = S2.identity_map()
        identity.coord_functions(polar, mercator)
        sol = c.solve(solution_key='sol-angle',
         parameters_values={tmin:0,tmax:2,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:8},
         verbose=False)
        interp = c.interpolate(solution_key='sol-angle',
                        interpolation_key='interp-angle', verbose=False)
        graph2D_mercator_angle_curve=c.plot_integrated(interpolation_key='interp-angle',
                             chart=mercator, thickness=1, verbose=False)
        expr_ph = ph0+v_ph0/v_th0*(ln(tan((v_th0*t+th0)/2))-ln(tan(th0/2)))
        c_loxo = S2.curve({polar:[th0+v_th0*t, expr_ph]}, (t,0,2),
                                                          name='c_loxo')
        c_loxo.expression(chart2=mercator)
        graph2D_mercator_loxo = c_loxo.plot(chart=mercator,
              parameters={th0:pi/4+2*pi/180, ph0:0.1, v_th0:1, v_ph0:8},
              thickness=1, color='blue')
        sphinx_plot(graph2D_mercator_angle_curve+graph2D_mercator_loxo)

    Both curves do have the same aspect.
    One may eventually compare these curves on :MATH:`\mathbb{S}^{2}`::

        sage: graph3D_embedded_angle_curve=c.plot_integrated(interpolation_key='interp-angle',
        ....:    mapping=euclid_embedding, thickness=3, verbose=False)
        sage: graph3D_embedded_loxo = c_loxo.plot(mapping=euclid_embedding,
        ....:  parameters={th0:pi/4+2*pi/180, ph0:0.1, v_th0:1, v_ph0:8},
        ....:  thickness=3, color = 'blue')
        sage: graph=graph3D_embedded_angle_curve + graph3D_embedded_loxo
        sage: graph += graph3D_embedded_polar_coords
        sage: graph.show()

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON')
        nab = S2.affine_connection('nab')
        nab.set_coef(epolar_ON)[:]
        [t, tmin, tmax, th0, ph0] = var('t tmin tmax th0 ph0')
        [v_th0, v_ph0, alpha] = var('v_th0 v_ph0 alpha')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                                  chart=polar, name='c')
        R3 = Manifold(3, 'R3', start_index=1)
        cart = R3.chart('X Y Z')
        [X,Y,Z] = var('X Y Z')
        euclid_embedding = S2.diff_map(R3,
              {(polar, cart):[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})
        sol = c.solve(solution_key='sol-angle',
         parameters_values={tmin:0,tmax:2,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:8},
         verbose=False)
        interp = c.interpolate(solution_key='sol-angle',
                        interpolation_key='interp-angle', verbose=False)
        graph3D_embedded_angle_curve=c.plot_integrated(interpolation_key='interp-angle',
                   mapping=euclid_embedding, thickness=3, verbose=False)
        expr_ph = ph0+v_ph0/v_th0*(ln(tan((v_th0*t+th0)/2))-ln(tan(th0/2)))
        c_loxo = S2.curve({polar:[th0+v_th0*t, expr_ph]}, (t,0,2),
                                                          name='c_loxo')
        graph3D_embedded_loxo = c_loxo.plot(mapping=euclid_embedding,
              parameters={th0:pi/4+2*pi/180, ph0:0.1, v_th0:1, v_ph0:8},
              thickness=3, color='blue')
        graph3D_embedded_polar_coords = polar.plot(chart=cart,
             mapping=euclid_embedding, number_values=15, color='yellow')
        graph = graph3D_embedded_angle_curve + graph3D_embedded_loxo
        graph += graph3D_embedded_polar_coords
        sphinx_plot(graph)

    """

    def __init__(self, parent, affine_connection, curve_parameter,
                 initial_tangent_vector, chart=None, name=None,
                 latex_name=None, is_isomorphism=False,
                 is_identity=False, verbose=True):
        r"""Construct an autoparallel curve with respect to the given
        affine connection with the given initial tangent vector.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, A, B] = var('t A B')
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[X.frame(),0,0,1],nab[X.frame(),2,1,2]=A*x1^2,B*x2*x3
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_autoparallel_curve(nab, (t, 0, 5), v,
            ....:                           name='c', verbose=False) ; c
            Integrated autoparallel curve c in the 3-dimensional
             differentiable manifold M
            sage: TestSuite(c).run()

        """

        # setting the chart to gain access to the coordinate functions
        if chart is None:
            chart = parent.codomain().default_chart()

        coordinate_functions = chart[:]
        velocities = chart.symbolic_velocities()

        if is_identity:
            equations_rhs = None
        else:
            dim = parent.codomain().dim()
            i0 = parent.codomain().start_index()
            equations_rhs = []

            gamma = affine_connection.coef()

            for alpha in range(dim):
                rhs = 0
                for mu in range(dim):
                    for nu in range(dim):
                        AUX = velocities[mu] * velocities[nu]
                        rhs-= gamma[alpha+i0, mu+i0, nu+i0].expr() * AUX
                        # 'AUX' only used for the line above to be shorter
                equations_rhs += [rhs.simplify_full()]

        IntegratedCurve.__init__(self, parent, equations_rhs,
                                 velocities, curve_parameter,
                                 initial_tangent_vector, chart=chart,
                                 name=name, latex_name=latex_name,
                                 is_isomorphism=is_isomorphism,
                                 is_identity=is_identity,
                                 verbose=verbose)

        self._affine_connection = affine_connection

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, A, B] = var('t A B')
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[X.frame(),0,0,1],nab[X.frame(),2,1,2]=A*x1^2,B*x2*x3
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_autoparallel_curve(nab, (t,0,5), v,
            ....:                                     verbose=False) ; c
            Integrated autoparallel curve in the 3-dimensional
             differentiable manifold M
            sage: c = M.integrated_autoparallel_curve(nab, (t, 0, 5), v,
            ....:                           name='c', verbose=False) ; c
            Integrated autoparallel curve c in the 3-dimensional
             differentiable manifold M

        """

        description = "Integrated autoparallel curve "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description

    def __reduce__(self):
        r"""
        Reduction function for the pickle protocole.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, A, B] = var('t A B')
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[X.frame(),0,0,1],nab[X.frame(),2,1,2]=A*x1^2,B*x2*x3
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_autoparallel_curve(nab, (t, 0, 5), v,
            ....:                               name='c', verbose=False)
            sage: c.__reduce__()
            (<class 'sage.manifolds.differentiable.integrated_curve.IntegratedAutoparallelCurveSet_with_category.element_class'>,
             (Set of Morphisms from Real interval (0, 5) to
              3-dimensional differentiable manifold M in Category of
              homsets of subobjects of sets and topological spaces which
              actually are integrated autoparallel curves w.r.t a
              certain affine connection,
              Affine connection nabla on the 3-dimensional
              differentiable manifold M,
              t,
              Tangent vector at Point p on the 3-dimensional differentiable manifold M,
              Chart (M, (x1, x2, x3)),
              'c',
              'c',
              False,
              False))

        Test of pickling::

            sage: loads(dumps(c))
            Integrated autoparallel curve c in the 3-dimensional differentiable manifold M

        """

        return (type(self), (self.parent(), self._affine_connection,
                self._curve_parameter, self._initial_tangent_vector,
                self._chart, self._name, self._latex_name,
                self._is_isomorphism, self._is_identity))

    def system(self, verbose=True):
        r"""
        Provide a detailed description of the system defining the
        autoparallel curve and returns the system defining it: chart,
        equations and initial conditions.

        INPUT:

        - ``verbose`` -- (default: ``True``) prints a detailed
          description of the curve

        OUTPUT:

        - list containing the attributes :attr:`equations_rhs`,
          :attr:`initial_tangent_vector` and :attr:`chart`

        EXAMPLE:

        System defining an autoparallel curve::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, A, B] = var('t A B')
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[X.frame(),0,0,1],nab[X.frame(),2,1,2]=A*x1^2,B*x2*x3
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_autoparallel_curve(nab, (t, 0, 5), v,
            ....:                                         verbose=False)
            sage: sys = c.system()
            Autoparallel curve in the 3-dimensional differentiable
             manifold M equipped with Affine connection nabla on the
             3-dimensional differentiable manifold M, and integrated
             over the Real interval (0, 5) as a solution to the
             following equations, written w.r.t.
             Chart (M, (x1, x2, x3)):
            <BLANKLINE>
            Initial point: Point p on the 3-dimensional differentiable
             manifold M with coordinates [0, 0, 0] w.r.t.
             Chart (M, (x1, x2, x3))
            Initial tangent vector: Tangent vector at Point p on the
             3-dimensional differentiable manifold M with
             components [1, 0, 1] w.r.t. Chart (M, (x1, x2, x3))
            <BLANKLINE>
            d(x1)/dt = Dx1
            d(x2)/dt = Dx2
            d(x3)/dt = Dx3
            d(Dx1)/dt = -A*Dx1*Dx2*x1^2
            d(Dx2)/dt = 0
            d(Dx3)/dt = -B*Dx2*Dx3*x2*x3
            <BLANKLINE>
            sage: sys_bis = c.system(verbose=False)
            sage: sys_bis == sys
            True

        """

        v0 = self._initial_tangent_vector
        chart = self._chart

        if verbose:
            initial_tgt_space = v0.parent()
            initial_pt = initial_tgt_space.base_point() # retrieves
            # the initial point as the base point of the tangent space
            # to which initial tangent vector belongs
            initial_pt_coords = list(initial_pt.coordinates(chart))
            # previous line converts to list since would otherwise be a
            # tuple ; will raise error if coordinates in chart are not
            # known

            initial_coord_basis = chart.frame().at(initial_pt)
            initial_tgt_vec_comps = v0[initial_coord_basis,:] # will
            # raise error if components in coordinate basis are not
            # known

            description = "Autoparallel curve "
            if self._name is not None:
                description += self._name + " "
            description += "in the {} ".format(self.codomain())
            description += "equipped with "
            description += "{}, ".format(self._affine_connection)
            description += "and integrated over the "
            description += "{} ".format(self.domain())
            description += "as a solution to the following equations, "
            description += "written w.r.t. "
            description += "{}:\n\n".format(chart)

            description += "Initial point: {} ".format(initial_pt)
            description += "with coordinates "
            description += "{} ".format(initial_pt_coords)
            description += "w.r.t. {}\n".format(chart)

            description += "Initial tangent vector: {} ".format(v0)
            description += "with components "
            description +="{}".format(initial_tgt_vec_comps)
            description += " w.r.t. {}\n\n".format(chart)

            for coord_func,velocity in zip(chart[:],self._velocities):
                description += "d({})/d{} = {}\n".format(coord_func,
                                                  self._curve_parameter,
                                                  velocity)

            for velocity,eqn in zip(self._velocities,self._equations_rhs):
                description += "d({})/d{} = {}\n".format(velocity,
                                                  self._curve_parameter,
                                                  eqn)

            print(description)

        return [self._equations_rhs, v0, chart]

class IntegratedGeodesic(IntegratedAutoparallelCurve):
    r"""
    Numerical geodesic on the manifold with respect to a
    given metric.

    INPUT:

    - ``parent`` --
      :class:`~sage.manifolds.differentiable.manifold_homset.DifferentiableCurveSet`
      the set of curves `\mathrm{Hom}(I, M)` to which the curve belongs
    - ``metric`` --
      :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
      metric with respect to which the curve is a geodesic
    - ``curve_parameter`` -- symbolic expression to be used as the
      parameter of the curve (the equations defining an instance of
      IntegratedGeodesic are such that ``t`` will actually be an affine
      parameter of the curve);
    - ``initial_tangent_vector`` --
      :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
      initial tangent vector of the curve
    - ``chart`` -- (default: ``None``) chart on the manifold in terms of
      which the equations are expressed; if ``None`` the default chart
      of the manifold is assumed
    - ``name`` -- (default: ``None``) string; symbol given to the curve
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
      the curve; if none is provided, ``name`` will be used
    - ``is_isomorphism`` -- (default: ``False``) determines whether the
      constructed object is a diffeomorphism; if set to ``True``,
      then `M` must have dimension one
    - ``is_identity`` -- (default: ``False``) determines whether the
      constructed object is the identity map; if set to ``True``,
      then `M` must coincide with the domain of the curve

    EXAMPLE:

    Geodesics of the unit 2-sphere :MATH:`\mathbb{S}^{2}`.
    Start with declaring the standard polar coordinates
    :MATH:`(\theta, \phi)` on :MATH:`\mathbb{S}^{2}` and the
    corresponding coordinate frame :MATH:`(e_{\theta}, e_{\phi})`::

        sage: S2 = Manifold(2, 'S^2', start_index=1)
        sage: polar.<th,ph>=S2.chart('th ph')
        sage: epolar = polar.frame()

    Set the Euclidean metric tensor :MATH:`g` induced on
    :MATH:`\mathbb{S}^{2}`::

        sage: g = S2.metric('g')
        sage: g[1,1], g[2,2] = 1, (sin(th))^2

    Set generic initial conditions for the geodesics to compute::

        sage: [th0, ph0, v_th0, v_ph0] = var('th0 ph0 v_th0 v_ph0')
        sage: p = S2.point((th0, ph0), name='p')
        sage: Tp = S2.tangent_space(p)
        sage: v = Tp((v_th0, v_ph0), basis=epolar.at(p))

    Declare the corresponding integrated geodesic and display the
    differential system it satisfies::

        sage: [t, tmin, tmax] = var('t tmin tmax')
        sage: c = S2.integrated_geodesic(g, (t, tmin, tmax), v,
        ....:                      chart=polar, name='c', verbose=False)
        sage: sys = c.system()
        Geodesic c in the 2-dimensional differentiable manifold S^2
         equipped with Riemannian metric g on the 2-dimensional
         differentiable manifold S^2, and integrated over the Real
         interval (tmin, tmax) as a solution to the following geodesic
         equations, written w.r.t. Chart (S^2, (th, ph)):
        <BLANKLINE>
        Initial point: Point p on the 2-dimensional differentiable
        manifold S^2 with coordinates [th0, ph0] w.r.t.
        Chart (S^2, (th, ph))
        Initial tangent vector: Tangent vector at Point p on the
        2-dimensional differentiable manifold S^2 with
        components [v_th0, v_ph0] w.r.t. Chart (S^2, (th, ph))
        <BLANKLINE>
        d(th)/dt = Dth
        d(ph)/dt = Dph
        d(Dth)/dt = Dph^2*cos(th)*sin(th)
        d(Dph)/dt = -2*Dph*Dth*cos(th)/sin(th)
        <BLANKLINE>

    Set a dictionnary providing the parameter range and the initial
    conditions for various geodesics::

        sage: dict_params={'equat':{tmin:0,tmax:3,th0:pi/2,ph0:0.1,v_th0:0,v_ph0:1},
        ....:   'longi':{tmin:0,tmax:3,th0:0.1,ph0:0.1,v_th0:1,v_ph0:0},
        ....:   'angle':{tmin:0,tmax:3,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:1}}

    Use :MATH:`\mathbb{R}^{3}` as the codomain of the standard map
    embedding :MATH:`(\mathbb{S}^{2}, (\theta, \phi))` in the
    3-dimensional Euclidean space::

        sage: R3 = Manifold(3, 'R3', start_index=1)
        sage: cart.<X,Y,Z> = R3.chart()
        sage: euclid_embedding = S2.diff_map(R3,
        ....:  {(polar, cart):[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})

    Solve, interpolate and prepare the plot for the solutions
    corresponding to the three initial conditions previously set::

        sage: graph3D_embedded_geods = Graphics()
        sage: for key in dict_params.keys():
        ....:     sol = c.solve(solution_key='sol-'+key,
        ....:         parameters_values=dict_params[key], verbose=False)
        ....:     interp = c.interpolate(solution_key='sol-'+key,
        ....:            interpolation_key='interp-'+key, verbose=False)
        ....:     graph3D_embedded_geods+=c.plot_integrated(interpolation_key='interp-'+key,
        ....:                      mapping=euclid_embedding, thickness=5,
        ....:                      display_tangent=True, scale=0.3,
        ....:                      width_tangent=0.5, verbose=False)

    Plot the resulting geodesics on the grid of polar coordinates lines
    on :MATH:`\mathbb{S}^{2}` and check that these are great circles::

        sage: graph3D_embedded_polar_coords = polar.plot(chart=cart,
        ....:                          mapping=euclid_embedding,
        ....:                          number_values=15, color='yellow')
        sage: graph=graph3D_embedded_geods+graph3D_embedded_polar_coords
        sage: graph.show()

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph')
        [th, ph] = var('th ph')
        epolar = polar.frame()
        g = S2.metric('g')
        g[1,1], g[2,2] = 1, (sin(th))**2
        [t,tmin,tmax,th0,ph0,v_th0,v_ph0]=var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar.at(p))
        c = S2.integrated_geodesic(g, (t, tmin, tmax), v, chart=polar,
                                                               name='c')
        dict_params={'equat':{tmin:0,tmax:3,th0:pi/2,ph0:0.1,v_th0:0,v_ph0:1},
               'longi':{tmin:0,tmax:3,th0:0.1,ph0:0.1,v_th0:1,v_ph0:0},
               'angle':{tmin:0,tmax:3,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:1}}
        R3 = Manifold(3, 'R3', start_index=1)
        cart = R3.chart('X Y Z')
        [X,Y,Z] = var('X Y Z')
        euclid_embedding = S2.diff_map(R3,
              {(polar, cart):[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})
        graph3D_embedded_geods = Graphics()
        for key in dict_params.keys():
            sol = c.solve(solution_key='sol-'+key,
                      parameters_values=dict_params[key], verbose=False)
            interp = c.interpolate(solution_key='sol-'+key,
                         interpolation_key='interp-'+key, verbose=False)
            graph3D_embedded_geods+=c.plot_integrated(interpolation_key='interp-'+key,
                                  mapping=euclid_embedding, thickness=5,
                                  display_tangent=True, scale=0.3,
                                  width_tangent=0.5, verbose=False)
        graph3D_embedded_polar_coords = polar.plot(chart=cart,
                                       mapping=euclid_embedding,
                                       number_values=15, color='yellow')
        graph = graph3D_embedded_geods + graph3D_embedded_polar_coords
        sphinx_plot(graph)

    """

    def __init__(self, parent, metric, curve_parameter,
                 initial_tangent_vector, chart=None, name=None,
                 latex_name=None, is_isomorphism=False,
                 is_identity=False, verbose=True):

        r"""Construct a geodesic curve with respect to the given metric
        with the given initial tangent vector.

        TESTS::

            sage: S2 = Manifold(2, 'S^2')
            sage: X.<theta,phi> = S2.chart()
            sage: [t, A] = var('t A')
            sage: g = S2.metric('g')
            sage: g[0,0] = A
            sage: g[1,0] = 0
            sage: g[1,1] = A*sin(theta)^2
            sage: p = S2.point((pi/2,0), name='p')
            sage: Tp = S2.tangent_space(p)
            sage: v = Tp((1/sqrt(2),1/sqrt(2)))
            sage: c = S2.integrated_geodesic(g, (t,0,pi), v, name='c',
            ....:                                     verbose=False) ; c
            Integrated geodesic c in the 2-dimensional differentiable
             manifold S^2
            sage: TestSuite(c).run()

        """

        if is_identity:
            affine_connection = None
        else:
            affine_connection = metric.connection()

        IntegratedAutoparallelCurve.__init__(self, parent,
                                   affine_connection, curve_parameter,
                                   initial_tangent_vector, chart=chart,
                                   name=name, latex_name=latex_name,
                                   is_isomorphism=is_isomorphism,
                                   is_identity=is_identity,
                                   verbose=verbose)

        self._metric = metric

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: S2 = Manifold(2, 'S^2')
            sage: X.<theta,phi> = S2.chart()
            sage: [t, A] = var('t A')
            sage: g = S2.metric('g')
            sage: g[0,0] = A
            sage: g[1,0] = 0
            sage: g[1,1] = A*sin(theta)^2
            sage: p = S2.point((pi/2,0), name='p')
            sage: Tp = S2.tangent_space(p)
            sage: v = Tp((1/sqrt(2),1/sqrt(2)))
            sage: c = S2.integrated_geodesic(g, (t, 0, pi), v,
            ....:                                     verbose=False) ; c
            Integrated geodesic in the 2-dimensional differentiable
             manifold S^2
            sage: c = S2.integrated_geodesic(g, (t,0,pi), v,
            ....:                           name='c', verbose=False) ; c
            Integrated geodesic c in the 2-dimensional differentiable
             manifold S^2

        """

        description = "Integrated geodesic "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description

    def __reduce__(self):
        r"""
        Reduction function for the pickle protocole.

        TESTS::

            sage: S2 = Manifold(2, 'S^2')
            sage: X.<theta,phi> = S2.chart()
            sage: [t, A] = var('t A')
            sage: g = S2.metric('g')
            sage: g[0,0] = A
            sage: g[1,0] = 0
            sage: g[1,1] = A*sin(theta)^2
            sage: p = S2.point((pi/2,0), name='p')
            sage: Tp = S2.tangent_space(p)
            sage: v = Tp((1/sqrt(2),1/sqrt(2)))
            sage: c = S2.integrated_geodesic(g, (t, 0, pi), v, name='c',
            ....:                                         verbose=False)
            sage: c.__reduce__()
            (<class 'sage.manifolds.differentiable.integrated_curve.IntegratedGeodesicSet_with_category.element_class'>,
             (Set of Morphisms from Real interval (0, pi) to
              2-dimensional differentiable manifold S^2 in Category of
              homsets of subobjects of sets and topological spaces which
              actually are integrated geodesics w.r.t a certain metric,
              Riemannian metric g on the 2-dimensional differentiable
              manifold S^2,
              t,
              Tangent vector at Point p on the 2-dimensional differentiable manifold S^2,
              Chart (S^2, (theta, phi)),
              'c',
              'c',
              False,
              False))

        Test of pickling::

            sage: loads(dumps(c))
            Integrated geodesic c in the 2-dimensional differentiable manifold S^2

        """

        return (type(self), (self.parent(), self._metric,
                self._curve_parameter, self._initial_tangent_vector,
                self._chart, self._name, self._latex_name,
                self._is_isomorphism, self._is_identity))

    def system(self, verbose=True):
        r"""
        Return the system defining the geodesic : chart, equations and
        initial conditions

        INPUT:

        - ``verbose`` -- (default: ``True``) prints a detailed
          description of the curve

        OUTPUT:

        - list containing the attributes :attr:`equations_rhs`,
          :attr:`initial_tangent_vector` and :attr:`chart`

        EXAMPLE:

        System defining a geodesic::

            sage: S2 = Manifold(2, 'S^2')
            sage: X.<theta,phi> = S2.chart()
            sage: [t, A] = var('t A')
            sage: g = S2.metric('g')
            sage: g[0,0] = A
            sage: g[1,0] = 0
            sage: g[1,1] = A*sin(theta)^2
            sage: p = S2.point((pi/2,0), name='p')
            sage: Tp = S2.tangent_space(p)
            sage: v = Tp((1/sqrt(2),1/sqrt(2)))
            sage: c = S2.integrated_geodesic(g, (t, 0, pi), v, name='c',
            ....:                                         verbose=False)
            sage: sys = c.system()
            Geodesic c in the 2-dimensional differentiable manifold S^2
             equipped with Riemannian metric g on the 2-dimensional
             differentiable manifold S^2, and integrated over the Real
             interval (0, pi) as a solution to the following geodesic
             equations, written w.r.t. Chart (S^2, (theta, phi)):
            <BLANKLINE>
            Initial point: Point p on the 2-dimensional differentiable
             manifold S^2 with coordinates [1/2*pi, 0] w.r.t.
             Chart (S^2, (theta, phi))
            Initial tangent vector: Tangent vector at Point p on the
             2-dimensional differentiable manifold S^2 with
             components [1/2*sqrt(2), 1/2*sqrt(2)] w.r.t.
             Chart (S^2, (theta, phi))
            <BLANKLINE>
            d(theta)/dt = Dtheta
            d(phi)/dt = Dphi
            d(Dtheta)/dt = Dphi^2*cos(theta)*sin(theta)
            d(Dphi)/dt = -2*Dphi*Dtheta*cos(theta)/sin(theta)
            <BLANKLINE>
            sage: sys_bis = c.system(verbose=False)
            sage: sys_bis == sys
            True

        """

        v0 = self._initial_tangent_vector
        chart = self._chart

        if verbose:
            initial_tgt_space = v0.parent()
            initial_pt = initial_tgt_space.base_point()#retrieves
            # the initial point as the base point of the tangent space
            # to which initial tangent vector belongs
            initial_pt_coords = list(initial_pt.coordinates(chart))
            # previous line converts to list since would otherwise be a
            # tuple ; will raise error if coordinates in chart are not
            # known

            initial_coord_basis = chart.frame().at(initial_pt)
            initial_tgt_vec_comps = v0[initial_coord_basis,:] # will
            # raise error if components in coordinate basis are not
            # known

            description = "Geodesic "
            if self._name is not None:
                description += self._name + " "
            description += "in the {} ".format(self.codomain())
            description += "equipped with "
            description += "{}, ".format(self._metric)
            description += "and integrated over the "
            description += "{} ".format(self.domain())
            description += "as a solution to the following "
            description += "geodesic equations, written w.r.t. "
            description += "{}:\n\n".format(chart)

            description += "Initial point: {} ".format(initial_pt)
            description += "with coordinates "
            description += "{} ".format(initial_pt_coords)
            description += "w.r.t. {}\n".format(chart)

            description += "Initial tangent vector: {} ".format(v0)
            description += "with components "
            description +="{}".format(initial_tgt_vec_comps)
            description += " w.r.t. {}\n\n".format(chart)

            for coord_func,velocity in zip(chart[:],self._velocities):
                description += "d({})/d{} = {}\n".format(coord_func,
                                                  self._curve_parameter,
                                                  velocity)

            for velocity,eqn in zip(self._velocities,self._equations_rhs):
                description += "d({})/d{} = {}\n".format(velocity,
                                                  self._curve_parameter,
                                                  eqn)

            print(description)

        return [self._equations_rhs, v0, chart]
